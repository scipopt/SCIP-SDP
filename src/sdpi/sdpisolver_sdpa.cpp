/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG  *//* shows all added nonzero entries */
/* #define SCIP_DEBUG_PRINTTOFILE */ /* prints each problem inserted into SDPA to the file sdpa.dat-s and the starting point to sdpa.ini-s */

/**@file   sdpisolver_sdpa.cpp
 * @brief  interface for SDPA
 * @author Tristan Gally
 * @author Ambros Gleixner
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"

/* turn off warnings for sdpa (doesn't seem to work) */
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include "sdpa_call.h"                       /* SDPA callable library interface */
#pragma GCC diagnostic warning "-Wstrict-prototypes"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */


/* Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMRY. */
#define BMS_CALL(x)   do                                                                                     \
                      {                                                                                      \
                         if( NULL == (x) )                                                                   \
                         {                                                                                   \
                            SCIPerrorMessage("No memory in function call.\n");                               \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/* This will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem. */
#define CHECK_IF_SOLVED(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            SCIPABORT();                                                                     \
                            return SCIP_LPERROR;                                                               \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )


/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SDPA*                 sdpa;               /**< solver-object */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in SDPA (nvars minus the number of variables with lb = ub) */
   int*                  inputtosdpamapper;  /**< entry i gives the index of input variable i in sdpa (starting from 1) or
                                               *  -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                               *  this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  sdpatoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in sdpa (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real*            fixedvarsobj;       /**< entry i gives the objective value of the i-th fixed variable */
   int                   nvarbounds;         /**< number of variable bounds given to sdpa, length of sdpavarboundpos */
   int*                  varboundpos;        /**< maps position of variable bounds in the variable bound part of the LP-block in sdpa to the sdpa-indices
                                               *  of the corresponding variables, -n means lower bound of variable n, +n means upper bound */
#ifndef NDEBUG
   SCIP_Bool             infeasible;         /**< true if the problem is infeasible during insertion/presolving (if constraints without active variables present) */
#endif
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< this is used to check if the SDP-Constraint is feasible */
   SCIP_Real             objlimit;           /**< objective limit for SDP solver */
   int                   threads;            /**< number of threads */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             penalty;            /**< was the problem last solved using a penalty formulation */
};


/*
 * Local Functions
 */

#ifndef NDEBUG
/** Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,          /**< pointer to an SDP interface solver structure */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( lb < ub + sdpisolver->epsilon );

   return ( REALABS(ub-lb) <= sdpisolver->epsilon );
}
#else
#define isFixed(sdpisolver,lb,ub) (REALABS(ub-lb) <= sdpisolver->epsilon)
#endif


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver, version number can only be printed to file/stdout */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "SDPA";
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Primal-dual Interior Point Solver for Semidefinite Programming developed by Katsuki Fujisawa et al. (sdpa.sourceforge.net)";
}

/** Does the solver have a way to solve a penalty formulation on its own or must one be provided? */
SCIP_Bool SCIPsdpiSolverKnowsPenalty(
   void
   )
{
   return TRUE;
}

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver           /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->sdpa;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;

   /* this will be properly initialized then calling solve */
   (*sdpisolver)->sdpa = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtosdpamapper = NULL;
   (*sdpisolver)->sdpatoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobj = NULL;
   (*sdpisolver)->nvarbounds = 0;
   (*sdpisolver)->varboundpos = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;
#ifndef NDEBUG
   (*sdpisolver)->infeasible = FALSE;
#endif

   (*sdpisolver)->epsilon = 1e-3;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->threads = 1;
   (*sdpisolver)->sdpinfo = FALSE;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   if (((*sdpisolver)->sdpa) != NULL)
   {
      /* free SDPA object using destructor and free memory via blockmemshell */
      //(*sdpisolver)->sdpa->~SDPA();
      delete (*sdpisolver)->sdpa;
      //(*sdpi)->sdpa.terminate(); //TODO which one to use?
      //BMSfreeMemory(&((*sdpisolver)->sdpa));
   }

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->nvars);

   if ( (*sdpisolver)->nvars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->inputtosdpamapper, (*sdpisolver)->nvars);

   if ( (*sdpisolver)->nactivevars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->sdpatoinputmapper, (*sdpisolver)->nactivevars);

   if ( (*sdpisolver)->nvars >= (*sdpisolver)->nactivevars )
   {
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsobj, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);
   }

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/**@} */


/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** loads and solves an SDP
 *
 *  For the non-constant SDP- and the LP-part, the original arrays before fixings should be given, for the constant
 *  SDP-part the arrays AFTER fixings should be given. In addition, an array needs to be given, that for every block and
 *  every row/col index within that block either has value -1, meaning that this index should be deleted, or a
 *  non-negative integer stating the number of indices before it that are to be deleated, meaning that this index will
 *  be decreased by that number, in addition to that the total number of deleted indices for each block should be given.
 *  Optionally an array start may be given with a starting point for the solver (if this is NULL then the solver should
 *  start from scratch).
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< this returns the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this, i.e.,
                                              *   the value to decrease this index by, this array should have memory allocated in the size
                                              *   sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int					    noldlpcons,		   /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*			          lprownactivevars,   /**< number of active variables for each lp constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start               /**< NULL or a starting point for the solver, this should have length nvars */
   )//TODO: start needs to include X,y,Z for SDPA
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
               nremovedinds, nlpcons, noldlpcons, lprhs, lprownactivevars, lpnnonz, lprow, lpcol, lpval, start);
}

/** loads and solves an SDP using a penalty formulation
 *
 *  The penalty formulation of the SDP is:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathds{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u.\f}
 *  Alternatively withObj can be set to false to set \f$ b \f$ to 0 and only check for feasibility (if the optimal objective value is
 *  bigger than 0 the problem is infeasible, otherwise it's feasible).
 *  For the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the constant SDP-part the arrays AFTER fixings
 *  should be given. In addition, an array needs to be given, that for every block and every row/col index within that block either has value
 *  -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices before it that are to be deleated,
 *  meaning that this index will be decreased by that number. Moreover, the total number of deleted indices for each block should be given.
 *  An optional starting point for the solver may be given; if it is NULL, the solver will start from scratch.
 *
 *  @warning This only works for some solvers, check with SCIPsdpiKnowsPenalty first, otherwise this returns an error (in which case you should form
 *  the penalty formulation yourself and pass it via LoadAndSolve).
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   SCIP_Real             gamma,              /**< the penalty parameter above, needs to be >= 0 */
   SCIP_Bool             withObj,            /**< if this is false, the objective is set to 0 */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< this returns the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                              *   the value to decrease this index by, this array should have memory allocated in the size
                                              *   sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*                  lprownactivevars,   /**< number of active variables for each lp constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start               /**< NULL or a starting point for the solver, this should have length nvars */
) //TODO: start needs to include X,y,Z for SDPA
{
   SCIP_Real* sdpavarbounds;
   int i;
   int k;
   int block;
   int nfixedvars;
   SCIP_Bool checkinput; /* should the input be checked for consistency in SDPA ? */
   int lpconsind;
   int lastrow;

   assert( sdpisolver != NULL );

   checkinput = FALSE;

#ifndef SCIP_DEBUG
   sdpisolver->infeasible = FALSE;
   checkinput = TRUE;
#endif

   /* set the penalty flag accordingly */
   sdpisolver->penalty = (gamma == 0.0) ? FALSE : TRUE;

   /* allocate memory for inputtosdpamapper, sdpatoinputmapper and the fixed variable information, for the latter this will later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtosdpamapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->sdpatoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsobj), sdpisolver->nvars - sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) );

   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;

   /* find the fixed variables */
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         sdpisolver->fixedvarsobj[nfixedvars] = obj[i];
         sdpisolver->fixedvarsval[nfixedvars] = lb[i]; /* if lb=ub, than this is the value the variable will have in every solution */
         nfixedvars++;
         sdpisolver->inputtosdpamapper[i] = -nfixedvars;
      }
      else
      {
         sdpisolver->sdpatoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->nactivevars++;
         sdpisolver->inputtosdpamapper[i] = sdpisolver->nactivevars; /* sdpa starts counting at 1, so we do this after increasing nactivevars */
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* shrink the fixedvars and sdpatoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsobj), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->sdpatoinputmapper), nvars, sdpisolver->nactivevars) );

   /* insert data */
   SCIPdebugMessage("Inserting Data into SDPA for SDP (%d) \n", ++sdpisolver->sdpcounter);

   if ( sdpisolver->sdpa != 0 )
   {
      /* if the SDPA solver has already been created, clear the current problem instance */
      //sdpisolver->sdpa->terminate();
      delete sdpisolver->sdpa;
      sdpisolver->sdpa = new SDPA();
   }
   else
   {
      /* we use this construction to allocate the memory for the SDPA object also via the blockmemshell */
      /* sdpisolver->sdpa = static_cast<SDPA*>(BMSallocMemoryCPP(sizeof(SDPA))); */
      sdpisolver->sdpa = new SDPA();
   }
   assert( sdpisolver->sdpa != 0 );

   /* compute number of variable bounds and save them in sdpavarbounds */
   sdpisolver->nvarbounds = 0;
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &sdpavarbounds, 2 * sdpisolver->nactivevars) );
   if ( sdpisolver->varboundpos == NULL )
   {
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->varboundpos), 2 * sdpisolver->nvars) );
   }
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      assert( 0 <= sdpisolver->sdpatoinputmapper[i] && sdpisolver->sdpatoinputmapper[i] < nvars );
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->sdpatoinputmapper[i]]))
      {
         sdpavarbounds[sdpisolver->nvarbounds] = lb[sdpisolver->sdpatoinputmapper[i]];
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = -(i + 1); /* negative sign means lower bound, i + 1 because sdpa starts counting from one */
         (sdpisolver->nvarbounds)++;
      }
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->sdpatoinputmapper[i]]))
      {
         sdpavarbounds[sdpisolver->nvarbounds] = ub[sdpisolver->sdpatoinputmapper[i]];
         sdpisolver->varboundpos[sdpisolver->nvarbounds] = +(i + 1); /* positive sign means upper bound, i + 1 because sdpa starts counting from one */
         (sdpisolver->nvarbounds)++;
      }
   }

   /* if we use a penalty formulation, we need the constraint r >= 0 */
   if ( gamma != 0.0 )
      sdpisolver->nvarbounds++;

   if ( sdpisolver->sdpinfo )
      sdpisolver->sdpa->setDisplay(stdout);
   else
      sdpisolver->sdpa->setDisplay(0);

#ifdef SCIP_MORE_DEBUG
   FILE* fpOut = fopen("output.tmp", "w");
   if ( ! fpOut )
      exit(-1);
   sdpisolver->sdpa->setResultFile(fpOut);
#endif

   /* if we want to use a starting point we have to tell SDPA to allocate memory for it */
   if ( start != NULL )
      sdpisolver->sdpa->setInitPoint(true);

   /* initialize blockstruct */
   if ( gamma == 0.0 ) /* we initialize this with an exact 0.0 in Solve without penalty */
      sdpisolver->sdpa->inputConstraintNumber(sdpisolver->nactivevars);
   else
      sdpisolver->sdpa->inputConstraintNumber(sdpisolver->nactivevars + 1); /* the additional variable is r which is multiplied with the identity matrix */

   /* if there are any lp-cons/variable-bounds, we get an extra block for those, lastrow - nshifts is the number of lp constraints added */
   sdpisolver->sdpa->inputBlockNumber((nlpcons + sdpisolver->nvarbounds > 0) ? nsdpblocks + 1 : nsdpblocks);

   /* block+1 because SDPA starts counting at 1 */
   for (block = 0; block < nsdpblocks; block++)
   {
      sdpisolver->sdpa->inputBlockSize(block + 1, sdpblocksizes[block] - nremovedinds[block]);
      sdpisolver->sdpa->inputBlockType(block + 1, SDPA::SDP);
   }
   if ( nlpcons + sdpisolver->nvarbounds > 0 )
   {
      /* the last block is the lp block, the size has a negative sign */
      sdpisolver->sdpa->inputBlockSize(nsdpblocks + 1, -(nlpcons + sdpisolver->nvarbounds));
      sdpisolver->sdpa->inputBlockType(nsdpblocks + 1, SDPA::LP);
   }
   sdpisolver->sdpa->initializeUpperTriangleSpace();

   /* set objective values */
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withObj )
      {
         /* insert objective value, SDPA counts from 1 to n instead of 0 to n-1 */
         sdpisolver->sdpa->inputCVec(i + 1, obj[sdpisolver->sdpatoinputmapper[i]]);
      }
      if (gamma != 0.0)
         sdpisolver->sdpa->inputCVec(sdpisolver->nactivevars + 1, gamma); /* set the objective of the additional var to gamma */
   }

   /* start inserting the non-constant SDP-Constraint-Matrices */
   if ( sdpnnonz > 0 )
   {
      int v;
      int blockvar;

      assert( nsdpblocks > 0 );
      assert( sdpnblockvarnonz != NULL );
      assert( sdpnblockvars != NULL );
      assert( sdpcol != NULL );
      assert( sdprow != NULL );
      assert( sdpval != NULL );
      assert( sdpvar != NULL );
      assert( indchanges != NULL );
      assert( nremovedinds != NULL );

      for (block = 0; block < nsdpblocks; block++)
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   -> building block %d (%d)\n", block + 1, sdpisolver->sdpcounter);
#endif
         for (i = 0; i < sdpisolver->nactivevars; i++)
         {
            /* we iterate over all non-fixed variables, so add them to sdpa for this block/var combination */
            v = sdpisolver->sdpatoinputmapper[i];

#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("      -> adding coefficient matrix for variable %d which becomes variable %d in SDPA (%d)\n", v, i, sdpisolver->sdpcounter);
#endif

            /* find the position of variable v in this block */
            blockvar = -1;
            for (k = 0; k < sdpnblockvars[block]; k++)
            {
               if (v == sdpvar[block][k])
               {
                  blockvar = k;
                  break;
               }
            }

            if ( blockvar > -1)  /* the variable exists in this block */
            {
               for (k = 0; k < sdpnblockvarnonz[block][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( indchanges[block][sdprow[block][blockvar][k]] > -1 && indchanges[block][sdpcol[block][blockvar][k]] > -1 );

                  assert( indchanges[block][sdprow[block][blockvar][k]] <= sdprow[block][blockvar][k]);
                  assert( indchanges[block][sdpcol[block][blockvar][k]] <= sdpcol[block][blockvar][k]);

                  assert( 0 <= sdprow[block][blockvar][k] && sdprow[block][blockvar][k] < sdpblocksizes[block] );
                  assert( 0 <= sdpcol[block][blockvar][k] && sdpcol[block][blockvar][k] < sdpblocksizes[block] );

                  /* rows and columns start with one in SDPA, so we have to add 1 to the indices */
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMessage("         -> adding nonzero %g at (%d,%d) (%d)\n",
                     sdpval[block][blockvar][k],
                     sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]] + 1,
                     sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]] + 1,
                     sdpisolver->sdpcounter);
#endif

                  sdpisolver->sdpa->inputElement(i + 1, block + 1,
                     sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]] + 1,
                     sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]] + 1,
                     sdpval[block][blockvar][k], checkinput);
               }
            }
         }
         /* insert the identity matrix if we are using a penalty formulation */
         if ( gamma != 0.0 )
         {
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("      -> adding coefficient matrix for variable penalty variable r in SDPA (%d)\n", i, v, sdpisolver->sdpcounter);
#endif
            for (i = 0; i < sdpblocksizes[i] - nremovedinds[i]; i++)
            {
#ifdef SCIP_MORE_DEBUG
                  SCIPdebugMessage("         -> adding nonzero 1.0 at (%d,%d) (%d)\n", i + 1, i + 1, sdpisolver->sdpcounter);
#endif

                  sdpisolver->sdpa->inputElement(sdpisolver->nactivevars + 1, block + 1, i + 1, i + 1, 1.0, checkinput);
            }
         }
      }
   }

   /* start inserting the constant matrix */
   if ( sdpconstnnonz > 0 )
   {
      assert( nsdpblocks > 0 );
      assert( sdpconstnblocknonz!= NULL );
      assert( sdpconstcol != NULL );
      assert( sdpconstrow != NULL );
      assert( sdpconstval != NULL );

      for (block = 0; block < nsdpblocks; block++)
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("   -> building block %d (%d)\n", block + 1, sdpisolver->sdpcounter);
#endif
         for (k = 0; k < sdpconstnblocknonz[block]; k++)
         {
            /* rows and cols with active nonzeros should not be removed */
            assert( indchanges[block][sdpconstrow[block][k]] > -1 && indchanges[block][sdpconstcol[block][k]] > -1 );

            assert( indchanges[block][sdpconstrow[block][k]] <= sdpconstrow[block][k]);
            assert( indchanges[block][sdpconstcol[block][k]] <= sdpconstcol[block][k]);

            assert (0 <= sdpconstrow[block][k] && sdpconstrow[block][k] < sdpblocksizes[block]);
            assert (0 <= sdpconstcol[block][k] && sdpconstcol[block][k] < sdpblocksizes[block]);

            /* rows and columns start with one in SDPA, so we have to add 1 to the indices, the constant matrix is given as variable 0 */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding constant nonzero %g at (%d,%d) (%d)\n", sdpconstval[block][k],
               sdpconstcol[block][k] - indchanges[block][sdpconstcol[block][k]] + 1,
               sdpconstrow[block][k] - indchanges[block][sdpconstrow[block][k]] + 1,
               sdpisolver->sdpcounter);
#endif
            sdpisolver->sdpa->inputElement(0, block + 1,
               sdpconstcol[block][k] - indchanges[block][sdpconstcol[block][k]] + 1,
               sdpconstrow[block][k] - indchanges[block][sdpconstrow[block][k]] + 1,
               sdpconstval[block][k], checkinput);
         }
      }
   }

   /* inserting LP nonzeros */
   lpconsind = 0;
   lastrow = -1; /* this together means, that we start numbering the rows at one, like sdpa wants it */
   for (i = 0; i < lpnnonz; i++)
   {
      assert( 0 <= lprow[i] && lprow[i] < noldlpcons );
      assert( 0 <= lpcol[i] && lpcol[i] < nvars );
      assert( REALABS(lpval[i]) > sdpisolver->epsilon );

      /* if the variable is active and the constraint is more than a bound, we add it */
      if ( sdpisolver->inputtosdpamapper[lpcol[i]] > 0 )
      {
       /* as this is an active variable, there should be at least one in the constraint */
         assert( lprownactivevars[lprow[i]] > 0 );
         if ( lprownactivevars[lprow[i]] > 1 )
         {
            if ( lprow[i] > lastrow )  /* we update the lpcons-counter */
            {
               lpconsind++;
               lastrow = lprow[i];
               /* if we use a penalty formulation, add the r * Identity entry */
               if ( gamma != 0.0 )
               {
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding nonzero 1.0 at (%d,%d) for penalty variable r in SDPA (%d)\n",
               lpconsind, lpconsind, sdpisolver->sdpcounter);
#endif
                  /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with blocks starting at 1, as are rows), the
                   * r-variable is variable nactivevars + 1 */
                  sdpisolver->sdpa->inputElement(sdpisolver->nactivevars + 1, nsdpblocks + 1, lpconsind, lpconsind, 1.0, checkinput);
               }
            }
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding nonzero %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
               lpval[i], lpconsind, lpconsind, lpcol[i], sdpisolver->inputtosdpamapper[lpcol[i]], sdpisolver->sdpcounter);
#endif
            /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with blocks starting at 1, as are rows) */
            sdpisolver->sdpa->inputElement(sdpisolver->inputtosdpamapper[lpcol[i]], nsdpblocks + 1,
               lpconsind, lpconsind, lpval[i], checkinput);
         }
      }
   }
   assert( lpconsind == nlpcons ); /* this is equal, because we number from one to nlpcons in sdpa */

   /* inserting LP right-hand-sides for active constraints*/
   for (i = 0; i < nlpcons; i++)
   {
      if ( REALABS(lprhs[i]) > sdpisolver->epsilon )
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("         -> adding rhs %g at (%d,%d) (%d)\n", lprhs[i], i+1, i+1, sdpisolver->sdpcounter);
#endif
         /* LP constraints are added as diagonal entries of the last block, right-hand-side is added as variable zero */
         sdpisolver->sdpa->inputElement(0, nsdpblocks + 1, i + 1, i + 1, lprhs[i], checkinput);
      }
   }

   /* print each LP-constraint as one formatted constraint in addition to the single entries inserted into SDPA */
#ifdef SCIP_MORE_DEBUG
   lpconsind = 0;
   lastrow = -1; /* this together means that we start numbering the rows at one, like sdpa wants it */
   for (i = 0; i < lpnnonz; i++)
   {
      assert( 0 <= lprow[i] && lprow[i] < noldlpcons );
      assert( 0 <= lpcol[i] && lpcol[i] < nvars );
      assert( REALABS(lpval[i]) > sdpisolver->epsilon );

      /* if the variable is active and the constraint is more than a bound, we add it */
      if ( sdpisolver->inputtosdpamapper[lpcol[i]] > 0 )
      {
         printf("+ %f <x%d> ", lpval[i], lpcol[i]);
       /* as this is an active variable, there should be at least one in the constraint */
         assert( lprownactivevars[lprow[i]] > 0 );
         if ( lprownactivevars[lprow[i]] > 1 )
         {
            if ( lprow[i] > lastrow )  /* we update the lpcons-counter */
            {
               printf(" >= %f\n", lprhs[lpconsind]);
               lpconsind++;
               lastrow = lprow[i];
            }
         }
      }
   }
   printf(" >= %f\n", lprhs[lpconsind]);
   assert( lpconsind == nlpcons ); /* this is equal, because we number from one to nlpcons in sdpa */
#endif

   /* insert variable bounds, these are also added as LP-constraints and therefore diagonal entries of the LP block
    * if we work with the penalty formulation, we get an extra entry for r >= 0, but this we will add afterwards */
   for (i = 0; i < ((gamma == 0.0) ? sdpisolver->nvarbounds : sdpisolver->nvarbounds - 1); i++)
   {
      assert( 0 < abs(sdpisolver->varboundpos[i]) && abs(sdpisolver->varboundpos[i] <= sdpisolver->nactivevars) ); /* the indices are already those for SDPA */

      /* for lower bound */
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* add it as an lp-constraint for this variable (- because we saved -n for the lower bound), at the position
          * (nactivelpcons + 1) + varbound-index, because we have >= the variable has coefficient +1 */
         sdpisolver->sdpa->inputElement(-sdpisolver->varboundpos[i], nsdpblocks + 1, nlpcons + 1 + i, nlpcons + 1 + i, 1.0, checkinput);

         if ( REALABS(sdpavarbounds[i]) > sdpisolver->epsilon )
         {
            /* the bound is added as the rhs and therefore variable zero */
            sdpisolver->sdpa->inputElement(0, nsdpblocks + 1, nlpcons + 1 + i, nlpcons + 1 + i, sdpavarbounds[i], checkinput);
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding lower bound %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  sdpavarbounds[i], nlpcons + 1 + i, nlpcons + 1 + i, sdpisolver->sdpatoinputmapper[-sdpisolver->varboundpos[i] - 1],
                  -sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
         }
         else
         {
            /* as the bound is zero, we don't need to add a right hand side */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding lower bound 0 at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  nlpcons + 1 + i, nlpcons + 1 + i, sdpisolver->sdpatoinputmapper[-sdpisolver->varboundpos[i] - 1],
               -sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
         }
      }
      else
      {
         /* this is an upper bound */

         /* add it as an lp-constraint for this variable, at the position nactivelpcons + varbound-index, because we have >= but we
          * want <= for the upper bound, we have to multiply by -1 and therefore the variable has coefficient -1 */
         sdpisolver->sdpa->inputElement(sdpisolver->varboundpos[i], nsdpblocks + 1, nlpcons + 1 + i, nlpcons + 1 + i, -1.0, checkinput);

         if ( REALABS(sdpavarbounds[i]) > sdpisolver->epsilon )
         {
            /* the bound is added as the rhs and therefore variable zero, we multiply by -1 for <= */
            sdpisolver->sdpa->inputElement(0, nsdpblocks + 1, nlpcons + 1 + i, nlpcons + 1 + i, -sdpavarbounds[i], checkinput);
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding upper bound %g at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  sdpavarbounds[i], nlpcons + 1 + i, nlpcons + 1 + i, sdpisolver->sdpatoinputmapper[sdpisolver->varboundpos[i] - 1],
                  sdpisolver->varboundpos[i], sdpisolver->sdpcounter);
#endif
         }
         else
         {
            /* as the bound is zero, we don't need to add a right hand side */
#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("         -> adding upper bound 0 at (%d,%d) for variable %d which became variable %d in SDPA (%d)\n",
                  0, nlpcons + 1 + i, nlpcons + 1 + i, sdpisolver->sdpatoinputmapper[sdpisolver->varboundpos[i] - 1],
                  sdpisolver->varboundpos[i]);
#endif
         }
      }
   }

   if ( gamma != 0.0 )
   {
      /* we add the variable bound r >= 0 */
      sdpisolver->sdpa->inputElement(sdpisolver->nactivevars + 1, nsdpblocks + 1, nlpcons + 1 + i, nlpcons + 1 + i, 1.0, checkinput);
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("         -> adding lower bound r >= 0 at (%d,%d)  in SDPA (%d)\n", nlpcons + 1 + i, nlpcons + 1 + i, sdpisolver->sdpcounter);
#endif
   }

   /* free the arrays used for counting and saving variable bounds and LP-right-hand-sides */
   BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpavarbounds, 2 * sdpisolver->nactivevars);

   /* transform the matrices to a more efficient form */
   sdpisolver->sdpa->initializeUpperTriangle();

#ifdef SCIP_DEBUG_PRINTTOFILE
   /* if necessary, dump input data and initial point */
   sdpisolver->sdpa->writeInputSparse(const_cast<char*>("sdpa.dat-s"), const_cast<char*>("%+8.3e"));
   sdpisolver->sdpa->writeInitSparse(const_cast<char*>("sdpa.ini-s"), const_cast<char*>("%+8.3e"));
#endif

   sdpisolver->sdpa->initializeSolve();

   /* set the starting solution */
   if (start != NULL)
   {
      //TODO: needs to be changed to y, Z and penalty formulation
      for (i = 1; i <= sdpisolver->nactivevars; i++) /* we iterate over the variables in sdpa */
         sdpisolver->sdpa->inputInitXVec(i, start[sdpisolver->sdpatoinputmapper[i] - 1]);
   }

   /* initialize settings */
   if ( gamma == 0.0 )
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_UNSTABLE_BUT_FAST);
   else
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW); /* if we already had problems with this problem, there is no reason to try fast */
   sdpisolver->sdpa->setParameterLowerBound(-1e20);
   sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->epsilon);
   sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->feastol);
   /* set the objective limit */
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
   else
      sdpisolver->sdpa->setParameterUpperBound(1e20);
#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif

   /* sdpisolver->sdpa->setParameterLambdaStar(1e5); */

   /* set number of threads */
   char str[1024];
   snprintf(str, 1024, "OMP_NUM_THREADS=%d", sdpisolver->threads);
   int status = putenv(str);
   if ( status )
   {
      SCIPdebugMessage("Setting the number of threads failed (%d, %d).\n", status, errno);
      return SCIP_LPERROR;
   }

   SCIPdebugMessage("Calling SDPA solve (SDP: %d, threads: %lld)\n", sdpisolver->sdpcounter, sdpisolver->sdpa->getNumThreads());
   sdpisolver->sdpa->solve();
   sdpisolver->solved = TRUE;

   /* check whether problem has been stably solved, if it wasn't and we didn't yet run the stable parametersettings (for the penalty formulation we do so), try
    * again with more stable parameters */
   if ( (! SCIPsdpiSolverIsAcceptable(sdpisolver)) && gamma == 0.0 )
   {
      SCIPdebugMessage("Numerical troubles -- solving SDP %d again ...\n", sdpisolver->sdpcounter);
      printf("Numerical troubles -- solving SDP %d again ...\n", sdpisolver->sdpcounter);

      /* initialize settings */
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_DEFAULT);
      sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->epsilon);
      sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->feastol);
      /* set the objective limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
      else
         sdpisolver->sdpa->setParameterUpperBound(1e20);
#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif
      sdpisolver->sdpa->solve();
      sdpisolver->solved = TRUE;
      /* if we still didn't converge, set the parameters even more conservativly */
      if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) )
      {
         SCIPdebugMessage("Numerical troubles -- solving SDP %d again^2 ...\n", sdpisolver->sdpcounter);
         printf("Numerical troubles -- solving SDP %d again^2 ...\n", sdpisolver->sdpcounter);

         /* initialize settings */
         sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);
         sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->epsilon);
         sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->feastol);
         /* set the objective limit */
         if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
            sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
         else
            sdpisolver->sdpa->setParameterUpperBound(1e20);
#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif
         sdpisolver->sdpa->solve();
         sdpisolver->solved = TRUE;
      }
   }

#ifdef SCIP_DEBUG
   /* print the phase value , i.e. whether solving was successfull */
   char phase_string[15];
   sdpisolver->sdpa->getPhaseString((char*)phase_string);
   SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are the same as in our formulation)\n", phase_string);
#endif

#ifdef SCIP_MORE_DEBUG
   (void) fclose(fpOut);
#endif

   return SCIP_OKAY;
}

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/* general comment: what we call dual is called primal in SDPA, so return the opposite of what sounds right */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   return sdpisolver->solved;
}

/** returns true if the solver could determine whether or not the problem is feasible
 *
 *  So it returns true if the solver knows that the problem is feasible/infeasible/unbounded, it returns false if the
 *  solver doesn't know anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be
 *  used
 */
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
      return FALSE;

   return TRUE;
}

/** gets information about primal and dual feasibility of the current SDP solution
 *  only call this after SCIPsdpiSolverFeasibilityKnown returned true */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      *primalfeasible = FALSE;
      *dualfeasible = FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   switch ( phasetype )
   {
   case SDPA::pdOPT:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pdFEAS:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pFEAS_dINF:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pINF_dFEAS:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;
   case SDPA::pUNBD:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      SCIPdebugMessage("SDPA stopped because dual objective became smaller than lower bound\n");
      break;
   case SDPA::dUNBD:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      SCIPdebugMessage("SDPA stopped because primal objective became bigger than upper bound\n");
      break;
   default: /* contains noInfo, pFeas, dFeas, pdInf */
      SCIPerrorMessage("SDPA doesn't know if primal and dual solutions are feasible\n");
      SCIPABORT();
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *
 *  This does not necessarily mean, that the solver knows and can return the primal ray.
 *  This is not implemented for all Solvers, always returns false (and a debug message) if it isn't.
 */
SCIP_Bool SCIPsdpiSolverExistsPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in SDPA!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 *
 *  This is not implemented for all Solvers, always returns false (and a debug message) if it isn't.
 */
SCIP_Bool SCIPsdpiSolverHasPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in SDPA!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is infeasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is feasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS || phasetype == SDPA::pdOPT || phasetype == SDPA::dFEAS  || phasetype == SDPA::pdFEAS )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point)
 *
 *  This does not necessarily mean, that the solver knows and can return the dual ray.
 *  This is not implemented for all Solvers, will always return false (and a debug message) if it isn't.
 */
SCIP_Bool SCIPsdpiSolverExistsDualRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in SDPA!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 *
 *  This is not implemented for all Solvers, will always return false (and a debug message) if it isn't.
 */
SCIP_Bool SCIPsdpiSolverHasDualRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in SDPA!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is infeasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is feasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF || phasetype == SDPA::pdOPT || phasetype == SDPA::pFEAS  || phasetype == SDPA::pdFEAS )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::dUNBD )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdFEAS )
   {
      if ( sdpisolver->sdpa->getParameterMaxIteration() == sdpisolver->sdpa->getIteration() )
         return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in SDPA!\n");
   return SCIP_LPERROR;
}

/** returns the internal solution status of the solver, which has the following meaning:
 * -1: solver wasn't started
 *  0: converged
 *  1: infeasible start
 *  2: numerical problems
 *  3: objective limit reached
 *  4: iteration limit reached
 *  5: time limit reached
 *  6: user termination
 *  7: other */
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return -1;
   }
#endif

   if ( sdpisolver->sdpa == NULL )
      return -1;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT || phasetype == SDPA::pFEAS_dINF || phasetype == SDPA::pINF_dFEAS )
      return 0;
   if ( phasetype == SDPA::pdINF )
      return 1;
   if ( phasetype == SDPA::dUNBD)
      return 3;
   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdFEAS )
      return 4;
   else /* should include pUNBD */
      return 7;
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 *  which is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED( sdpisolver );
#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   phasetype = sdpisolver->sdpa->getPhaseValue();

   /* we are happy if we converged, or we reached the objective limit (dUNBD) or we could show that our (dual, primal for SDPA) problem is
    * infeasible [except for numerics] */
   if ( SCIPsdpiSolverIsConverged(sdpisolver) || phasetype == SDPA::dUNBD || phasetype == SDPA::pINF_dFEAS)
      return TRUE;
   else
   {
      double pobj;
      double dobj;
      double gap;

      /* if it didn't converge check the optimality gap */
      pobj = sdpisolver->sdpa->getDualObj();
      dobj = sdpisolver->sdpa->getPrimalObj();

      gap = REALABS(pobj - dobj);

      if ( (gap < sdpisolver->epsilon) || ((gap / (0.5 * (REALABS(pobj) + REALABS(dobj)))) < sdpisolver->epsilon) ) /* this is the duality gap used in SDPA */
         return TRUE;
   }
   return FALSE;

/* TODO: also check for primal feasibility, as this is also needed for optimality */
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");

   /* todo: change settings to stable */
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   int v;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   *objval = sdpisolver->sdpa->getPrimalObj();

#ifndef NDEBUG
   SCIP_Real primalval = sdpisolver->sdpa->getDualObj();
   SCIP_Real gap = (REALABS(*objval - primalval) / (0.5 * (REALABS(primalval) + REALABS(*objval)))); /* duality gap used in SDPA */
   if ( gap > sdpisolver->epsilon )
      SCIPdebugMessage("Attention: got objective value (before adding values of fixed variables) of %f in SCIPsdpiSolverGetObjval, "
            "but primal objective is %f with duality gap %f!\n", *objval, primalval, gap );
#endif

   /* todo: compute this value when setting up fixed variables */
   /* as we didn't add the fixed (lb = ub) variables to sdpa, we have to add their contributions to the objective by hand */
   for (v = 0; v < sdpisolver->nvars; v++)
   {
      if (sdpisolver->inputtosdpamapper[v] < 0)
         *objval += sdpisolver->fixedvarsobj[-sdpisolver->inputtosdpamapper[v] - 1] * sdpisolver->fixedvarsval[-sdpisolver->inputtosdpamapper[v] - 1];
   }

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs
 *
 *  If dualsollength isn't equal to the number of variables this will return the needed length and a debug message.
 */
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   )
{
   SCIP_Real* sdpasol;
   int v;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   assert( dualsollength != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   if ( objval != NULL )
   {
      *objval = sdpisolver->sdpa->getPrimalObj();

#ifndef NDEBUG
      SCIP_Real primalval = sdpisolver->sdpa->getDualObj();
      SCIP_Real gap = (REALABS(*objval - primalval) / (0.5 * (REALABS(primalval) + REALABS(*objval)))); /* duality gap used in SDPA */
      if ( gap > sdpisolver->epsilon )
      {
         SCIPdebugMessage("Attention: got objective value (before adding values of fixed variables) of %f in SCIPsdpiSolverGetSol, "
            "but primal objective is %f with duality gap %f!\n", *objval, primalval, gap );
      }
#endif

      /* todo: compute this value when setting up fixed variables */
      /* as we didn't add the fixed (lb = ub) variables to sdpa, we have to add their contributions to the objective by hand */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if (sdpisolver->inputtosdpamapper[v] < 0)
            *objval += sdpisolver->fixedvarsobj[-sdpisolver->inputtosdpamapper[v] - 1] * sdpisolver->fixedvarsval[-sdpisolver->inputtosdpamapper[v] - 1];
      }
   }

   if ( *dualsollength > 0 )
   {
      assert( dualsol != NULL );
      if ( *dualsollength < sdpisolver->nvars )
      {
         SCIPdebugMessage("The given array in SCIPsdpiSolverGetSol only had length %d, but %d was needed", *dualsollength, sdpisolver->nvars);
         *dualsollength = sdpisolver->nvars;

         return SCIP_OKAY;
      }

      /* get the solution from sdpa */
      assert( sdpisolver->nactivevars == sdpisolver->sdpa->getConstraintNumber() );
      sdpasol = sdpisolver->sdpa->getResultXVec();
      /* insert the entries into dualsol, for non-fixed vars we copy those from sdpa, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if (sdpisolver->inputtosdpamapper[v] > 0)
         {
            /* minus one because the inputtosdpamapper gives the sdpa indices which start at one, but the array starts at 0 */
            dualsol[v] = sdpasol[sdpisolver->inputtosdpamapper[v] - 1];
         }
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            assert( -sdpisolver->inputtosdpamapper[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtosdpamapper[v]) - 1];
         }
      }
   }
   return SCIP_OKAY;
}

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown. Note: if a variable is either fixed or unbounded in the dual
 *  problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_RETCODE SCIPsdpiSolverGetPrimalBoundVars(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            lbvars,             /**< returns the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< returns the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{
   int i;
   SCIP_Real* X; /* block of primal solution matrix corresponding to the LP-part */
   int lpblockind;
   int nlpcons;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );
   CHECK_IF_SOLVED( sdpisolver );

   /* check if the arrays are long enough */
   if ( *arraylength < sdpisolver->nvars )
   {
      *arraylength = sdpisolver->nvars;
      SCIPdebugMessage("Insufficient length of array in SCIPsdpiSolverGetPrimalBoundVars (gave %d, needed %d)\n", *arraylength, sdpisolver->nvars);
      return SCIP_OKAY;
   }

   /* initialize the return-arrays with zero */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      lbvars[i] = 0.0;
      ubvars[i] = 0.0;
   }

   /* if no variable bounds were added, we return the zero vector (we do this separately, because in this case there might be no LP-block) */
   if ( sdpisolver->nvarbounds == 0 )
   {
      SCIPdebugMessage("Asked for PrimalBoundVars, but there were no variable bounds in sdpa, returning zero vector !");
      return SCIP_OKAY;
   }

   /* get the block of primal solution matrix corresponding to the LP-part from sdpa */
   lpblockind = sdpisolver->sdpa->getBlockNumber(); /* the LP block is the last one and sdpa counts from one */
   assert( sdpisolver->sdpa->getBlockType(lpblockind) == SDPA::LP );
   nlpcons = sdpisolver->sdpa->getBlockSize(lpblockind);
   assert( nlpcons >= 0 );

   X = sdpisolver->sdpa->getResultYMat(lpblockind);

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars || (sdpisolver->nvarbounds <= 2 * sdpisolver->nvars + 1 && sdpisolver->penalty ) );
   /* if we solved a penalty formulation, the last variable bound belongs to the penalty variable, which we aren't interested in here */
   for (i = 0; i < ((sdpisolver->penalty) ? sdpisolver->nvarbounds - 1 : sdpisolver->nvarbounds); i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* this is a lower bound */
         lbvars[sdpisolver->sdpatoinputmapper[-1 * sdpisolver->varboundpos[i] -1]] = X[nlpcons - sdpisolver->nvarbounds + i]; /* the first nlpcons entries correspond to lp-constraints */
      }
      else
      {
         /* this is an upper bound */
         ubvars[sdpisolver->sdpatoinputmapper[+1 * sdpisolver->varboundpos[i] - 1]] = X[nlpcons - sdpisolver->nvarbounds + i]; /* the first nlpcons entries correspond to lp-constraints */
      }
   }

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   assert( iterations != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      *iterations = 0;
      return SCIP_OKAY;
   }
#endif

   *iterations = sdpisolver->sdpa->getIteration();
   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP solver */
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   return 1E+20; /* default infinity from SCIP */
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
}

/** returns highest penalty parameter to be used */
SCIP_Real SCIPsdpiSolverMaxPenParam(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   return 1E+10;  /* not yet implemented */
}

/** checks if given value is greater or equal to the highest penalty parameter to be used */
SCIP_Bool SCIPsdpiSolverIsGEMaxPenParam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             val                 /**< value to be compared to maximum penalty parameter */
   )
{
   return ((val <= -SCIPsdpiSolverMaxPenParam(sdpisolver)) || (val >= SCIPsdpiSolverMaxPenParam(sdpisolver)));
}

/** gets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   assert( sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_EPSILON:
      *dval = sdpisolver->epsilon;
      break;
   case SCIP_SDPPAR_FEASTOL:
      *dval = sdpisolver->feastol;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverSetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_EPSILON:
      sdpisolver->epsilon = dval;
      SCIPdebugMessage("Setting sdpisolver epsilon to %f.\n", dval);
      break;
   case SCIP_SDPPAR_FEASTOL:
      sdpisolver->feastol = dval;
      SCIPdebugMessage("Setting sdpisolver feastol to %f.\n", dval);
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %f.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_THREADS:
      *ival = sdpisolver->threads;
      SCIPdebugMessage("Getting sdpisolver number of threads: %d.\n", *ival);
      break;
   case SCIP_SDPPAR_SDPINFO:
      *ival = sdpisolver->sdpinfo;
      SCIPdebugMessage("Getting sdpisolver information output (%d).\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_THREADS:
      sdpisolver->threads = ival;
      SCIPdebugMessage("Setting sdpisolver number of threads to %d.\n", ival);
      break;
   case SCIP_SDPPAR_SDPINFO:
      sdpisolver->sdpinfo = ival;
      SCIPdebugMessage("Setting sdpisolver information output (%d).\n", ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   assert( fname != NULL );

   sdpisolver->sdpa->writeInputSparse(const_cast<char*>(fname), const_cast<char*>("%8.3f"));

   return SCIP_OKAY;
}

/**@} */
