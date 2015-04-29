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

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/

/**@file   sdpisolver_dsdp.c
 * @brief  interface for DSDP
 * @author Tristan Gally
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"

/* turn off warning for DSDSP */
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */
#pragma GCC diagnostic warning "-Wstrict-prototypes"


#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */

/* Calls a DSDP-Function and transforms the return-code to a SCIP_LPERROR if needed. */
#define DSDP_CALL(x)  do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            SCIPABORT();                                                                     \
                            return SCIP_LPERROR;                                                               \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/* Same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced. */
#define DSDP_CALLM(x) do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            SCIPABORT();                                                                     \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

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
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP to add variable bounds to */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in DSDP (nvars minus the number of variables with lb = ub) */
   int*                  inputtodsdpmapper;  /**< entry i gives the index of input variable i in dsdp (starting from 1) or
                                               *  -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                               *  this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  dsdptoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in dsdp (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
#ifndef NDEBUG
   SCIP_Bool             infeasible;         /**< true if the problem is infeasible during insertion/presolving (if constraints without active variables present) */
#endif
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< this is used to check if the SDP-Constraint is feasible */
   SCIP_Real             objlimit;           /**< objective limit for SDP solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
};


/*
 * Local Functions
 */

/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static
int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}

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

   return ( REALABS(ub-lb) <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (REALABS(ub-lb) <= sdpisolver->epsilon)
#endif

/** sort the given row, col and val arrays first by non-decreasing col-indices, than for those with identical col-indices by non-increasing row-indices */
static
void sortColRow(
   int*                  row,                /* row indices */
   int*                  col,                /* column indices */
   SCIP_Real*            val,                /* values */
   int                   length              /* length of the given arrays */
   )
{
   int firstentry;
   int nextentry = 0;

   /* first sort by col indices */
   SCIPsortIntIntReal(col, row, val, length);

   /* for those with identical col-indices now sort by non-decreasing row-index, first find all entries with the same col-index */
   while (nextentry < length)
   {
      firstentry = nextentry; /* the next col starts where the last one ended */

      while (nextentry < length && col[nextentry] == col[firstentry]) /* as long as the row still matches, increase nextentry */
         ++nextentry;

      /* now sort all entries between firstentry and nextentry-1 by their row-indices */
      SCIPsortIntReal(row + firstentry, val + firstentry, nextentry - firstentry);
   }
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver, getting version doesn't seem to be supported by DSDP */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "DSDP";
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Dual-Scaling Interior Point Solver for Semidefinite Programming developed by Steve Benson, Yinyu Ye, and Xiong Zhang (http://www.mcs.anl.gov/hs/software/DSDP/)";
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
   return (void*) sdpisolver->dsdp;
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

   /* the following four variables will be properly initialized only immediatly prior to solving because DSDP and the
    * SDPCone need information about the number of variables and sdpblocks during creation */
   (*sdpisolver)->dsdp = NULL;
   (*sdpisolver)->sdpcone = NULL;
   (*sdpisolver)->lpcone = NULL;
   (*sdpisolver)->bcone = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtodsdpmapper = NULL;
   (*sdpisolver)->dsdptoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;
#ifndef NDEBUG
   (*sdpisolver)->infeasible = FALSE;
#endif

   (*sdpisolver)->epsilon = 1e-3;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
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

   if ( (*sdpisolver)->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy((*sdpisolver)->dsdp) );
   }

   if ( (*sdpisolver)->nvars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->inputtodsdpmapper, (*sdpisolver)->nvars);

   if ( (*sdpisolver)->nactivevars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->dsdptoinputmapper, (*sdpisolver)->nactivevars);

   if ( (*sdpisolver)->nvars >= (*sdpisolver)->nactivevars )
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);

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
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int                   noldlpcons,         /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*                  rownactivevars,     /**< number of active variables for each lp constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start               /**< NULL or a starting point for the solver, this should have length nvars */
   )
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars,
           sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
           indchanges, nremovedinds, blockindchanges, nremovedblocks, nlpcons, noldlpcons, lprhs, rownactivevars, lpnnonz, lprow, lpcol, lpval, start, NULL);
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
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   int	                noldlpcons,		   /**< number of LP-constraints including those with less than two active nonzeros */
   SCIP_Real*            lprhs,              /**< right hand sides of active LP rows after fixings (may be NULL if nlpcons = 0) */
   int*	                rownactivevars,	   /**< number of active variables for each lp constraint */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   SCIP_Bool*            feasorig            /**< is the solution to the penalty-formulation feasible for the original problem? (may be NULL if gamma = 0) */
)
{
   int* dsdpconstind = NULL;  /* indices for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpconstval = NULL; /* non-zero values for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdpind = NULL;       /* indices for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpval = NULL;    /* non-zero values for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplpbegcol = NULL;  /* starting-indices for all columns in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplprow = NULL;     /* row indices in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdplpval = NULL;  /* nonzeroes in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int i;
   int j;
   int ind;
   int block;
   int startind;
   int nfixedvars;
   int dsdpnlpnonz;

#ifdef SCIP_DEBUG
   DSDPTerminationReason reason; /* this will later be used to check if DSDP converged */
#endif

   SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", ++sdpisolver->sdpcounter);

   assert( sdpisolver != NULL );
   assert( gamma > - sdpisolver->epsilon );

#ifndef NDEBUG
   sdpisolver->infeasible = FALSE;
#endif

   /* allocate memory for inputtodsdpmapper, dsdptoinputmapper and the fixed variable information, for the latter this will later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtodsdpmapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) );

   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;

   /* find the fixed variables */
   sdpisolver->fixedvarsobjcontr = 0.0;
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         nfixedvars++;
         sdpisolver->inputtodsdpmapper[i] = -nfixedvars;
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* this is the value this variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars - 1] = lb[i]; /* if lb=ub, than this is the value the variable will have in every solution */
         SCIPdebugMessage("Fixing variable %d locally to %f for SDP %d in DSDP\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
         sdpisolver->dsdptoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->nactivevars++;
         sdpisolver->inputtodsdpmapper[i] = sdpisolver->nactivevars; /* dsdp starts counting at 1, so we do this after increasing nactivevars */
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in DSDP\n", i, sdpisolver->inputtodsdpmapper[i], sdpisolver->sdpcounter);
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* shrink the fixedvars and dsdptoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), nvars, sdpisolver->nactivevars) );

   /* insert data */

   if ( sdpisolver->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy(sdpisolver->dsdp) ); /* if there already exists a DSDP-instance, destroy the old one */
   }

   DSDP_CALLM( DSDPCreate(sdpisolver->nactivevars, &(sdpisolver->dsdp)) );
   DSDP_CALLM( DSDPCreateSDPCone(sdpisolver->dsdp, nsdpblocks - nremovedblocks, &(sdpisolver->sdpcone)) );
   DSDP_CALLM( DSDPCreateLPCone(sdpisolver->dsdp, &(sdpisolver->lpcone)) );
   DSDP_CALLM( DSDPCreateBCone(sdpisolver->dsdp, &(sdpisolver->bcone)) );

#ifdef SCIP_MORE_DEBUG
         printf("setting objective values for SDP %d:\n", sdpisolver->sdpcounter);
#endif
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withObj )
      {
         /* insert objective value, DSDP counts from 1 to n instead of 0 to n-1, *(-1) because DSDP maximizes instead of minimizing */
         DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, i+1, -1.0 * obj[sdpisolver->dsdptoinputmapper[i]]) );
#ifdef SCIP_MORE_DEBUG
         printf("var %d (was var %d): %f, ", i+1, sdpisolver->dsdptoinputmapper[i], obj[sdpisolver->dsdptoinputmapper[i]]);
#endif
      }
      else
      {
         DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, i+1, 0.0) );
      }
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->dsdptoinputmapper[i]]) )
      {
         /* insert lower bound, DSDP counts from 1 to n instead of 0 to n-1 */
         DSDP_CALL( BConeSetLowerBound(sdpisolver->bcone, i+1, lb[sdpisolver->dsdptoinputmapper[i]]) );
      }
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->dsdptoinputmapper[i]]) )
      {
         /* insert upper bound, DSDP counts from 1 to n instead of 0 to n-1 */
         DSDP_CALL(BConeSetUpperBound(sdpisolver->bcone, i+1, ub[sdpisolver->dsdptoinputmapper[i]]));
      }
   }

#ifdef SCIP_MORE_DEBUG
   printf("\n");
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpisolver->bcone);
#endif

   /* set blocksizes */
   for (block = 0; block < nsdpblocks; ++block)
   {
      /* only insert blocksizes for the blocks we didn't remove */
      if ( blockindchanges[block] > -1 )
      {
         /* (blocks are counted from 0 to m-1) */
         DSDP_CALL( SDPConeSetBlockSize(sdpisolver->sdpcone, block- blockindchanges[block], sdpblocksizes[block] - nremovedinds[block]) );
      }
   }

   /* start inserting the non-constant SDP-Constraint-Matrices */
   if ( sdpnnonz > 0 )
   {
      int v;
      int k;
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

      /* allocate memory */
      /* This needs to be one long array, because DSDP uses it for solving so all nonzeros have to be in it and it may not be freed before the
       * problem is solved. The distinct blocks/variables (for the i,j-parts) are then given by dsdpind + startind, which gives a pointer to the
       * first array-element belonging to this block and then the number of elements in this block is given to DSDP for iterating over it. */

      /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz) );
      /* values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz) );

      ind = 0; /* this will be used for iterating over the nonzeroes */

      for (block = 0; block < nsdpblocks; block++)
      {
         for (i = 0; i < sdpisolver->nactivevars; i++)
         {
            /* we iterate over all non-fixed variables, so add them to the dsdp arrays for this block/var combination */
            v = sdpisolver->dsdptoinputmapper[i];

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

            startind = ind;

            if (blockvar > -1) /* the variable exists in this block */
            {
               for (k = 0; k < sdpnblockvarnonz[block][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( indchanges[block][sdprow[block][blockvar][k]] > -1 && indchanges[block][sdpcol[block][blockvar][k]] > -1 );

                  /* substract the number of removed indices before the row and col to get the indices after fixings */
                  dsdpind[ind] = compLowerTriangPos(sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]],
                                                    sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]]);
                  dsdpval[ind] = -1 * sdpval[block][blockvar][k];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
                  ind++;
               }

               /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
               SCIPsortIntReal(dsdpind + startind, dsdpval + startind, sdpnblockvarnonz[block][blockvar]);

               assert( blockindchanges[block] > -1 ); /* we shouldn't insert into blocks we removed */

               /* i + 1 because DSDP starts counting the variables at 1, adding startind shifts the arrays to the first
                * nonzero belonging to this block and this variable */
               DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], i + 1, sdpblocksizes[block] - nremovedinds[block],
                     1.0, 0, dsdpind + startind,dsdpval + startind, sdpnblockvarnonz[block][blockvar]));
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

      /* allocate memory */

      /* DSDP uses these for solving, so they may not be freed before the problem is solved. */

      /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz) );
      /* values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz) );

      ind = 0;

      for (block = 0; block < nsdpblocks; block++)
      {
         startind = ind; /* starting index of this block in the dsdpconst arrays */

         if ( sdpconstnblocknonz[block] > 0 )
         {
            /* insert the constant-nonzeros */
            for (i = 0; i < sdpconstnblocknonz[block]; i++)
            {
               /* rows and cols with nonzeros should not be removed */
               assert( indchanges[block][sdpconstrow[block][i]] > -1 && indchanges[block][sdpconstcol[block][i]] > -1 );

               /* substract the number of deleted indices before this to get the index after variable fixings */
               dsdpconstind[ind] = compLowerTriangPos(sdpconstrow[block][i] - indchanges[block][sdpconstrow[block][i]],
                                                      sdpconstcol[block][i] - indchanges[block][sdpconstcol[block][i]]);
               dsdpconstval[ind] = -1 * sdpconstval[block][i]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j) should be positive semidefinite */
               ind++;
            }

            /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
            SCIPsortIntReal(dsdpconstind + startind, dsdpconstval + startind, sdpconstnblocknonz[block]);

            assert( blockindchanges[block] > -1 ); /* we shouldn't insert into a block we removed */

            /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block by adding
             * startind, ind - startind gives the number of elements for this block */
            DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], 0, sdpblocksizes[block] - nremovedinds[block],
                                               1.0, 0, dsdpconstind + startind, dsdpconstval + startind, ind - startind));
         }
      }
   }

#ifdef SCIP_MORE_DEBUG
   SDPConeView2(sdpisolver->sdpcone);
#endif

   /* start inserting the LP constraints */
   if ( nlpcons > 0 || lpnnonz > 0 || ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
   {
	   int nextcol;
	   int* rowmapper;
	   int pos;

      assert( noldlpcons > 0 );
      assert( lprhs != NULL );
      assert( lpcol != NULL );
      assert( lprow != NULL );
      assert( lpval != NULL );

      /* memory allocation */

      /* these arrays are needed in DSDP during solving, so they may only be freed afterwards */
      /* dsdplpbegcol[i] gives the number of nonzeros in column 0 (right hand side) till i-1 (i going from 1 till m, with extra entries 0 (always 0) and m+1 (always lpcons + lpnnonz)) */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2) );
      /* dsdplprow saves the row indices of the LP for DSDP */
      /* worst-case length is lpnnonz + nlpcons, because right hand sides are also included in the vector, this will be
       * shortened after the exact length after fixings is known, in case we have an objective limit, this is increased by one entry for the right-hand-side
       * and at most nvars entries for the nonzeros */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpcons + lpnnonz) );
      else
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpcons + 1) + lpnnonz + nvars) );

      /* values given to DSDP */
      /* worst-case length is lpnnonz + nlpcons, because right hand sides are also included in the vector, this will be
       * shortened after the exact length after fixings is known, in case we have an objective limit, this is increased by one entry for the right-hand-side
       * and at most nvars entries for the nonzeros */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz) );
      else
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpcons + 1) + lpnnonz + nvars) );

      /* allocate memory to save which lpconstraints are mapped to which index */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &rowmapper, noldlpcons) );

      /* compute the rowmapper */
      pos = 0;
      for (i = 0; i < noldlpcons; i++)
      {
    	  if ( rownactivevars[i] >= 2 )
    	  {
    		  rowmapper[i] = pos;
    		  pos++;
    	  }
    	  else
    		  rowmapper[i] = -1;
      }

      /* add all right-hand-sides that are greater than zero */
      dsdpnlpnonz = 0;
      for (i = 0; i < nlpcons; i++)
      {
    	  if (REALABS(lprhs[i]) > sdpisolver->epsilon )
    	  {
    		  dsdplprow[dsdpnlpnonz] = i;
    		  dsdplpval[dsdpnlpnonz] = -lprhs[i]; /* we multiply by -1 because DSDP wants <= instead of >= */
    		  dsdpnlpnonz++;
    	  }
      }

      /* add the right-hand-side for the objective bound */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if (REALABS(sdpisolver->objlimit) > sdpisolver->epsilon )
         {
            dsdplprow[dsdpnlpnonz] = nlpcons; /* this is the last lp-constraint, as DSDP counts from 0 to nlpcons-1, this is number nlpcons */
            dsdplpval[dsdpnlpnonz] = sdpisolver->objlimit; /* as we want <= upper bound, this is the correct type of inequality for DSDP */
            dsdpnlpnonz++;
         }
      }

      /* now add the nonzeros */

      /* for this we have to sort the nonzeros by col first and then by row, as this is the sorting DSDP wants */
      sortColRow(lprow, lpcol, lpval, lpnnonz);

      /* iterate over all nonzeros to add the active ones to the dsdp arrays and compute dsdplpbegcol */
      nextcol = 0;
      dsdplpbegcol[0] = 0;
      for (i = 0; i < lpnnonz; i++)
      {
         /* if a new variable starts, set the corresponding dsdplpbegcol-entry */
         if (lpcol[i] >= nextcol)
         {
            /* set the dsdplpbegcol entries, as there might be active variables which appear only in the sdp but not the lp-part, we also have to set
             * the starting values for all variables in between to the same value (as we also set the entry for the found variable, this for-queue
             * will always have at least one index in the index set) */
            for (j = nextcol; j <= lpcol[i]; j++)
            {
               if (sdpisolver->inputtodsdpmapper[j] >= 0)
               {
                  assert( ! (isFixed(sdpisolver, lb[j], ub[j])) );
                  dsdplpbegcol[sdpisolver->inputtodsdpmapper[j]] = dsdpnlpnonz;

                  /* add the entry to the objlimit-lp-constraint for the last variables */
                  if ( (! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit)) && (REALABS( obj[j] ) > sdpisolver->epsilon))
                  {
                     dsdplprow[dsdpnlpnonz] = nlpcons;
                     dsdplpval[dsdpnlpnonz] = obj[j];
                     dsdpnlpnonz++;
                  }
               }
            }
            nextcol = j; /* this also equals lpcol[i]+1 */
         }
         /* add the nonzero, if it isn't fixed and the row isn't to be deleted (because it is only a variable bound) */
         if ( ! isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]) && rownactivevars[lprow[i]] > 1 )
         {
            /* the index is adjusted for deleted lp rows, also rows are numbered 0,...,nlpcons-1 in DSDP, as they are
             * here, nlpcons is added to the index as the first nlpcons entries correspond to the right hand sides */
            assert( rowmapper[lprow[i]] > -1 );
            dsdplprow[dsdpnlpnonz] = rowmapper[lprow[i]];
            dsdplpval[dsdpnlpnonz] = -lpval[i]; /* - because dsdp wants <= instead of >= constraints */
            dsdpnlpnonz++;
         }
#ifndef SCIP_NDEBUG
         /* if this is an active nonzero for the row, it should have at least one active var */
         else
            assert( isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]) || rownactivevars[lprow[i]] == 1 );
#endif
      }

      /* set the begcol array for all remaining variables (that aren't part of the LP-part), and also set the objlimit-constraint-entries */
      for (j = nextcol; j < nvars; j++)
      {
         if (sdpisolver->inputtodsdpmapper[j] >= 0)
         {
            assert( ! (isFixed(sdpisolver, lb[j], ub[j])) );
            dsdplpbegcol[sdpisolver->inputtodsdpmapper[j]] = dsdpnlpnonz;
            /* add the entry to the objlimit-lp-constraint for the last variables */
            if ( (! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit)) && (REALABS( obj[j] ) > sdpisolver->epsilon))
            {
               dsdplprow[dsdpnlpnonz] = nlpcons;
               dsdplpval[dsdpnlpnonz] = obj[j];
               dsdpnlpnonz++;
            }
         }
      }

      dsdplpbegcol[sdpisolver->nactivevars + 1] = dsdpnlpnonz;

      /* free the memory for the rowshifts-array */
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &rowmapper, noldlpcons);

      /* shrink the dsdplp-arrays */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpcons + lpnnonz, dsdpnlpnonz) );
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz, dsdpnlpnonz) );
      }
      else
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, (nlpcons + 1) + lpnnonz + nvars, dsdpnlpnonz) );
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, (nlpcons + 1) + lpnnonz + nvars, dsdpnlpnonz) );
      }

      /* add the arrays to dsdp */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         DSDP_CALL( LPConeSetData(sdpisolver->lpcone, nlpcons, dsdplpbegcol, dsdplprow, dsdplpval) );
      }
      else
      {
         DSDP_CALL( LPConeSetData(sdpisolver->lpcone, nlpcons + 1, dsdplpbegcol, dsdplprow, dsdplpval) );
      }
#ifdef SCIP_MORE_DEBUG
      LPConeView(sdpisolver->lpcone);
#endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpisolver->sdpcounter);

   DSDP_CALL( DSDPSetGapTolerance(sdpisolver->dsdp, sdpisolver->epsilon) );  /* set DSDP's tolerance for duality gap */
   DSDP_CALL( DSDPSetRTolerance(sdpisolver->dsdp, sdpisolver->feastol) );    /* set DSDP's tolerance for the SDP-constraints */
   if ( sdpisolver-> sdpinfo )
   {
      DSDP_CALL( DSDPSetStandardMonitor(sdpisolver->dsdp, 1) );   /* output DSDP information after every 1 iteration */
   }


   /* set the penalty parameter */
   if ( gamma >= sdpisolver->epsilon ) /* in sdpisolverSolve this is called with an exact 0 */
   {
      DSDP_CALL( DSDPSetPenaltyParameter(sdpisolver->dsdp, gamma) );
      DSDP_CALL( DSDPUsePenalty(sdpisolver->dsdp, 1) );
   }

   /* set the starting solution */
   if (start != NULL)
   {
      for (i = 1; i <= sdpisolver->nactivevars; i++) /* we iterate over the variables in DSDP */
      {
         DSDP_CALL( DSDPSetY0(sdpisolver->dsdp, i, start[sdpisolver->dsdptoinputmapper[i]]) );
      }
   }

   /* start the solving process */
   DSDP_CALLM( DSDPSetup(sdpisolver->dsdp) );
   DSDP_CALL( DSDPSolve(sdpisolver->dsdp) );

   DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) ); /* computes X and determines feasibility and unboundedness of the solution */
   sdpisolver->solved = TRUE;

   /*these arrays were used to give information to DSDP and were needed during solving and for computing X, so they may only be freed now*/
   if ( sdpconstnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz);
   }

   if ( sdpnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz);
   }

   if (nlpcons > 0 || lpnnonz > 0)
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, dsdpnlpnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, dsdpnlpnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2);
   }

#ifdef SCIP_DEBUG
   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, &reason));

   switch ( reason )
   {
   case DSDP_CONVERGED:
      SCIPdebugMessage("DSDP converged!\n");
      break;

   case DSDP_INFEASIBLE_START:
      SCIPdebugMessage("DSDP started with an infeasible point!\n");
      break;

   case DSDP_SMALL_STEPS:
      SCIPdebugMessage("Short step lengths created by numerical difficulties prevented progress in DSDP!\n");
      break;

   case DSDP_INDEFINITE_SCHUR_MATRIX:
      SCIPdebugMessage("Schur Matrix in DSDP was indefinite but should have been positive semidefinite!\n");
      break;

   case DSDP_MAX_IT:
      SCIPdebugMessage("DSDP reached maximum number of iterations!\n");
      break;

   case DSDP_NUMERICAL_ERROR:
      SCIPdebugMessage("A numerical error occured in DSDP!\n");
      break;

   case DSDP_UPPERBOUND:
      SCIPdebugMessage("Dual objective value in DSDP reached upper bound.\n");
      break;

   case DSDP_USER_TERMINATION:
      SCIPdebugMessage("DSDP didn't stop solving, did you?\n");
      break;

   case CONTINUE_ITERATING:
      SCIPdebugMessage("DSDP wants to continue iterating but somehow was stopped!\n");
      break;

   default:
      SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
      break;
   }
#endif

   if ( gamma >= sdpisolver->epsilon )
   {
      double rval;

      DSDP_CALL( DSDPGetR(sdpisolver->dsdp, &rval) );
      *feasorig = (rval < sdpisolver->feastol );
   }

   return SCIP_OKAY;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

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
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   if ( pdfeasible == DSDP_PDUNKNOWN )
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
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
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

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   switch ( pdfeasible )
   {
   case DSDP_PDFEASIBLE:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;

   case DSDP_UNBOUNDED:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;

   case DSDP_INFEASIBLE:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;

   default: /* should only include DSDP_PDUNKNOWN */
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
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
   SCIPdebugMessage("Not implemented in DSDP!\n");
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
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return FALSE;

   return TRUE;
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
   SCIPdebugMessage("Not implemented in DSDP!\n");
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
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible));
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   DSDP_CALL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   if ( reason == DSDP_CONVERGED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Method not implemented for DSDP, as objective limit is given as an ordinary LP-constraint ! \n");

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, &reason));

   if ( reason == DSDP_MAX_IT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return SCIP_LPERROR;
}

/** returns the internal solution status of the solver
 *  -1: solver wasn't started
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
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return -1;
   }
#endif

   if ( sdpisolver->dsdp == NULL )
      return -1;

   DSDP_CALL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   switch ( reason )
   {
   case DSDP_CONVERGED:
      return 0;

   case DSDP_INFEASIBLE_START:
      return 1;

   case DSDP_SMALL_STEPS:
      return 2;

   case DSDP_INDEFINITE_SCHUR_MATRIX:
      return 2;

   case DSDP_MAX_IT:
      return 4;

   case DSDP_NUMERICAL_ERROR:
      return 2;

   case DSDP_UPPERBOUND:
      return 3;

   case DSDP_USER_TERMINATION:
      return 6;

   default:
      return 7;
   }
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   return (SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverIsPrimalFeasible(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 *  which is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   if ( SCIPsdpiSolverIsConverged(sdpisolver) )
      return TRUE;

   return FALSE;
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert( sdpisolver != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

#ifndef NDEBUG
   if ( sdpisolver->infeasible )
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
   *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

   /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective by hand */
   *objval += sdpisolver->fixedvarsobjcontr;

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
   double* dsdpsol;
   int v;

   assert( sdpisolver != NULL );
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
      DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
      *objval *= -1; /* DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

      /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective by hand */
      *objval += sdpisolver->fixedvarsobjcontr;
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

      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, sdpisolver->nactivevars) );
      DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, sdpisolver->nactivevars) ); /* last entry needs to be the number of variables, will return an error otherwise */

      /* insert the entries into dualsol, for non-fixed vars we copy those from dsdp, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if (sdpisolver->inputtodsdpmapper[v] > -1)
         {
            /* minus one because the inputtodsdpmapper gives the dsdp indices which start at one, but the array starts at 0 */
            dualsol[v] = dsdpsol[sdpisolver->inputtodsdpmapper[v] - 1];
         }
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtodsdpmapper[v]) - 1];
         }
      }
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, sdpisolver->nactivevars);
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
   double* lbvarsdsdp;
   double* ubvarsdsdp;
   int i;

   assert( sdpisolver != NULL );
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

   /* allocate memory for the arrays given to DSDP */
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &lbvarsdsdp, sdpisolver->nactivevars) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &ubvarsdsdp, sdpisolver->nactivevars) );

   /* get the values for the active variables from DSDP */
   DSDP_CALL( BConeCopyX(sdpisolver->bcone, lbvarsdsdp, ubvarsdsdp, sdpisolver->nactivevars) );

   /* copy them to the right spots of lbvars & ubvars */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      if (sdpisolver->inputtodsdpmapper[i] < 0)
      {
         /* if the variable was fixed, it didn't exist in the relaxation, so we set the value to 0
          * (as DSDP already uses this value for unbounded vars) */
         lbvars[i] = 0;
         ubvars[i] = 0;
      }
      else
      {
         lbvars[i] = lbvarsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
         ubvars[i] = ubvarsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
      }
   }

   /* free allocated memory */
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &ubvarsdsdp, sdpisolver->nactivevars);
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &lbvarsdsdp, sdpisolver->nactivevars);

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpisolver != NULL );
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

   DSDP_CALL( DSDPGetIts(sdpisolver->dsdp, iterations) );
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
   return 1E+10;  /* DSDP will start with penalty param 10^10 if called normally */
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
      /* DSDP only allows to set a dual bound, but as we want to solve the dual problem in DSDP, we would need to set a primal bound, which doesn't exist in
       * DSDP, so we can't do anything in this case. */
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
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/**@} */
