/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define SCIP_DEBUG
#define SCIP_MORE_DEBUG

/**@file   sdpi_dsdp.c
 * @brief  interface for dsdp
 * @author Tristan Gally
 */

#include <assert.h>

#include "sdpi/sdpi.h"

#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */

/* calls a DSDP-Function and transforms the return-code to a SCIP_ERROR if needed */
#define DSDP_CALL(x)  do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/* same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced */
#define DSDP_CALLM(x) do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_NOMEMORY;                                                            \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/* Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMRY */
#define BMS_CALL(x)   do                                                                                      \
                       {                                                                                      \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/* this will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem */
#define CHECK_IF_SOLVED(sdpisolver)  do                                                                       \
                        {                                                                                     \
                           if (!(sdpisolver->solved))                                                         \
                           {                                                                                  \
                              SCIPerrorMessage("Tried to access solution information ahead of solving! \n");  \
                              SCIPABORT();                                                                    \
                              return SCIP_ERROR;                                                              \
                           }                                                                                  \
                        }                                                                                     \
                        while( FALSE )

struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP to add variable bounds to */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in DSDP (nvars minus the number of variables with lb=ub) */
   int*                  inputtodsdpmapper;  /**< entry i gives the index of input variable i in dsdp (starting from 1) or
                                               *  -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                               *  this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  dsdptoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in dsdp (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real*            fixedvarsobj;       /**< entry i gives the objective value of the i-th fixed variable */
#ifndef NDEBUG
   SCIP_Bool             infeasible;         /**< this is set to true if the problem is found infeasible during insertion/presolving (if there are
                                               *  LP constraints without active variables l */
#endif
   SCIP_Bool             solved;             /**< was the SDP solved since the problem was last changed */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< this is used for checking if primal and dual objective are equal */
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
/**
 * Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
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

/**
 * sort the given row, col and val arrays first by non-decreasing col-indices, than for those by identical col-indices by non-increasing row-indices
 */
static
void sortColRow(
   int*                  row,                /* row indices */
   int*                  col,                /* column indices */
   SCIP_Real*            val,                /* values */
   int                   length              /* length of the given arrays */
   )
{
   int firstentry;
   int nextentry;

   /* first sort by col indices */
   SCIPsortIntIntReal(col, row, val, length);

   /* for those with identical col-indices now sort by non-decreasing row-index, first find all entries with the same col-index */
   nextentry = 0;
   while (nextentry < length)
   {
      firstentry = nextentry; /* the next col starts where the last one ended*/

      while (nextentry < length && col[nextentry] == col[firstentry]) /* as long as the row still matches, increase nextentry */
      {
         nextentry++;
      }

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

/** Does the solver have a way to solve a penalty formulation on its own or must one be provided */
const SCIP_Bool SCIPsdpiSolverKnowsPenalty(
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
   DSDP newdsdp;
   SDPCone newsdpcone;
   LPCone newlpcone;
   BCone newbcone;

   assert ( sdpisolver != NULL );
   assert ( blkmem != NULL );

   /* these will be properly initialized only immediatly prior to solving because DSDP and the SDPCone need information about the number
    * of variables and sdpblocks during creation */
   newdsdp = NULL;
   newsdpcone = NULL;
   newlpcone = NULL;
   newbcone = NULL;

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL(BMSallocBlockMemory(blkmem, sdpisolver));

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->dsdp = newdsdp;
   (*sdpisolver)->sdpcone = newsdpcone;
   (*sdpisolver)->lpcone = newlpcone;
   (*sdpisolver)->bcone = newbcone;
   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtodsdpmapper = NULL;
   (*sdpisolver)->dsdptoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobj = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;
#ifndef NDEBUG
   (*sdpisolver)->infeasible = FALSE;
#endif

   (*sdpisolver)->epsilon = 1e-6;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert ( sdpisolver != NULL );
   assert ( *sdpisolver != NULL );

   if (((*sdpisolver)->dsdp) != NULL)
   {
   DSDP_CALL(DSDPDestroy((*sdpisolver)->dsdp));
   }

   if ((*sdpisolver)->nvars > 0)
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->inputtodsdpmapper, (*sdpisolver)->nvars);

   if ((*sdpisolver)->nactivevars > 0)
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->dsdptoinputmapper, (*sdpisolver)->nvars);

   if ((*sdpisolver)->nvars > (*sdpisolver)->nactivevars)
   {
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsobj, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);
   }

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/**@} */


/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** inserts the SDP, then solves the SDP, for the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the
 *  constant SDP-part the arrays AFTER fixings should be given, in addition to that an array needs to be given, that for every block and every row/col
 *  index within that block either has value -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices
 *  before it that are to be deleated, meaning that this index will be decreased by that number, in addition to that the total number of deleted
 *  indices for each block should be given
 *
 *  attention: the given lp arrays will be sorted in their original position */
EXTERN
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks AFTER FIXINGS*/
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                *  number of entries  of sdpconst row/col/val [i] AFTER FIXINGS*/
   int**                 sdpconstrow,        /**< pointers to row-indices for each block  AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                               *  in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< this returns the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                               *  the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                               *  the value to decrease this index by, this array should have memory allocated in the size
                                               *  sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, will get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, will get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries, will get sorted (may be NULL if lpnnonz = 0) */
   )
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars,
           sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
           indchanges, nremovedinds, nlpcons, lprhs, lpnnonz, lprow, lpcol, lpval);
}

/** inserts the SDP, then solves the following penalty formulation of the SDP:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u}
 *  \f
 *  alternatively withObj can be to false to set \f b \f to false and only check for feasibility (if the optimal
 *  objective value is bigger than 0 the problem is infeasible, otherwise it's feasible)
 *  For the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the constant SDP-part the arrays AFTER fixings
 *  should be given, in addition to that an array needs to be given, that for every block and every row/col index within that block either has value
 *  -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices before it that are to be deleated,
 *  meaning that this index will be decreased by that number, in addition to that the total number of deleted indices for each block should be given.
 *
 *  attention: this only works for some solvers, check with SCIPsdpiKnowsPenalty first, otherwise this returns an error (in that case you should form
 *  the penalty formulation yourself and pass it via LoadAndSolve
 *
 *  attention: the given lp arrays will be sorted in their original position */
EXTERN
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   SCIP_Real             penaltyParam,       /**< the penalty parameter Gamma above, needs to be >= 0 */
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
                                                *  number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                               *  in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< this returns the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                               *  the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                               *  the value to decrease this index by, this array should have memory allocated in the size
                                               *  sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, will be sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, will be sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries, will be sorted (may be NULL if lpnnonz = 0) */
)
{
   int* dsdpconstind;         /* indices for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpconstval;      /* non-zero values for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdpind;              /* indices for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpval;          /* non-zero values for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplpbegcol;         /* starting-indices for all columns in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplprow;         /* row indices in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdplpval;         /* nonzeroes in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int dsdplparraylength;
   int i;
   int j;
   int ind;
   int block;
   int startind;
   int** fixedind; /* these two arrays will save the fixed nonzeros for each block for later adding them to the constant arrays */
   int** fixedval; /* "----------------------------------------------------------------------------------------" */
   int* nfixednonz;
   int totalnfixednonz; /* sum of fixed nonzeros over all blocks */
   int nfixedvars;

#ifdef SCIP_DEBUG
   DSDPTerminationReason* reason; /* this will later be used to check if DSDP converged */
#endif

   assert ( sdpisolver != NULL );
   assert ( penaltyParam >= 0.0 );

   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;
#ifndef NDEBUG
   sdpisolver->infeasible = FALSE;
#endif

   /* allocate memory for inputtosdpmapper and the fixed variable information, for the latter this will later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtodsdpmapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsobj), sdpisolver->nvars - sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) );

   /* find the fixed variables */
   for (i = 0; i < nvars; i++)
   {
      if (isFixed(sdpisolver, lb[i], ub[i]))
      {
         nfixedvars++;
         sdpisolver->inputtodsdpmapper[i] = -nfixedvars;
         sdpisolver->fixedvarsobj[nfixedvars - 1] = obj[i];
         sdpisolver->fixedvarsval[nfixedvars - 1] = lb[i]; /* if lb=ub, than this is the value the variable will have in every solution */
      }
      else
      {
         sdpisolver->dsdptoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->nactivevars++;
         sdpisolver->inputtodsdpmapper[i] = sdpisolver->nactivevars;
      }
   }

   assert ( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* shrink the fixedvars arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsobj), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );

   /* if there are fixed variables, prepare arrays to save fixed nonzeros for later adding them to the constant nonzeros */
   if (sdpisolver->nactivevars < sdpisolver->nvars)
   {
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &fixedind, nsdpblocks));
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &fixedval, nsdpblocks));

      for (i = 0; i < nsdpblocks; i++)
      {
         nfixednonz[block] = 0;
         BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &(fixedind[block]), sdpnnonz));
         BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &(fixedval[block]), sdpnnonz));
      }
   }

   /* insert data */

   SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", sdpisolver->sdpcounter);

   if (sdpisolver->dsdp != NULL)
   {
      DSDP_CALL(DSDPDestroy(sdpisolver->dsdp)); /* if there already exists a DSDP-instance, destroy the old one */
   }

   DSDP_CALLM(DSDPCreate(sdpisolver->nactivevars, &(sdpisolver->dsdp)));
   DSDP_CALLM(DSDPCreateSDPCone(sdpisolver->dsdp, nsdpblocks, &(sdpisolver->sdpcone)));
   DSDP_CALLM(DSDPCreateLPCone(sdpisolver->dsdp, &(sdpisolver->lpcone)));
   DSDP_CALLM(DSDPCreateBCone(sdpisolver->dsdp, &(sdpisolver->bcone)));

   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if (withObj)
      {
         DSDP_CALL(DSDPSetDualObjective(sdpisolver->dsdp, i+1, -1.0 * obj[sdpisolver->dsdptoinputmapper[i]])); /* insert objective value, DSDP counts
                                                           * from 1 to n instead of 0 to n-1, *(-1) because DSDP maximizes instead of minimizing */
      }
      else
      {
         DSDP_CALL(DSDPSetDualObjective(sdpisolver->dsdp, i+1, 0.0));
      }
      if (!SCIPsdpiSolverIsInfinity(sdpisolver, lb[sdpisolver->dsdptoinputmapper[i]]))
      {
         DSDP_CALL(BConeSetLowerBound(sdpisolver->bcone, i+1, lb[sdpisolver->dsdptoinputmapper[i]])); /* insert lower bound, DSDP counts from 1 to n
                                                                                                       * instead of 0 to n-1 */
      }
      if (!SCIPsdpiSolverIsInfinity(sdpisolver, ub[sdpisolver->dsdptoinputmapper[i]]))
      {
         DSDP_CALL(BConeSetUpperBound(sdpisolver->bcone, i+1, ub[sdpisolver->dsdptoinputmapper[i]])); /* insert upper bound, DSDP counts from 1 to n
                                                                                                       * instead of 0 to n-1 */
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpisolver->bcone);
#endif

   /* set blocksizes */
   for(i = 0; i < nsdpblocks; i++)
   {
      DSDP_CALL(SDPConeSetBlockSize(sdpisolver->sdpcone, i, sdpblocksizes[i] - nremovedinds[i])); /* (blocks are counted from 0 to m-1) */
   }


   /*start inserting the non-constant SDP-Constraint-Matrices */
   if(sdpnnonz > 0)
   {
      int var;
      int k;

      assert ( nsdpblocks > 0 );
      assert ( sdpnblockvarnonz != NULL );
      assert ( sdpnblockvars != NULL );
      assert ( sdpcol != NULL );
      assert ( sdprow != NULL );
      assert ( sdpval != NULL );
      assert ( sdpvar != NULL );
      assert ( indchanges != NULL );
      assert ( nremovedinds != NULL );

      /*allocate memory */
      /*This needs to be one long array, because DSDP uses it for solving so all nonzeros have to be in it and it may not be freed before the problem is solved. The distinct blocks/variables
       *(for the i,j-parts) are then given by dsdpind + startind, which gives a pointer to the first array-element belonging to this block and then the number of
       *elements in this block is given to DSDP for iterating over it */

      /*indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz));
      /*values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz));

      ind = 0; /* this will be used for iterating over the nonzeroes */

      for(block = 0; block < nsdpblocks; block++)
      {
         for(var = 0; var < sdpnblockvars[block]; var++)
         {
            if (sdpisolver->inputtodsdpmapper[sdpvar[block][var]] > -1)
            {
               /* this variable isn't fixed, so add it to the dsdp arrays for this block/var combination (otherwise we can ignore it, as the constant
                * matrix after fixings has already been computed */
               startind = ind;

               for (k = 0; k < sdpnblockvarnonz[block][var]; k++)
               {
                  assert ( indchanges[block][sdprow[block][var][k]] > -1 && indchanges[block][sdpcol[block][var][k]] > -1 ); /* rows and cols with
                                                                                                         * active nonzeros should not be removed */
                  dsdpind[ind] = compLowerTriangPos(sdprow[block][var][k] - indchanges[block][sdprow[block][var][k]],
                                                    sdpcol[block][var][k] - indchanges[block][sdpcol[block][var][k]]); /* substract the number of
                                                    * removed indices before the row and col to get the indices after fixings */
                  dsdpval[ind] = -1 * sdpval[block][var][k];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
                  ind++;
               }

               /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
               SCIPsortIntReal(dsdpind + startind, dsdpval + startind, sdpnblockvarnonz[block][var]);

               DSDP_CALL(SDPConeSetASparseVecMat(sdpisolver->sdpcone, block, sdpisolver->inputtodsdpmapper[sdpvar[block][var]], sdpblocksizes[block], 1.0, 0, dsdpind + startind,
                  dsdpval + startind, sdpnblockvarnonz[block][var])); /* the inputtosdpsolver already adds the +1 shift needed for dsdp, adding startind shifts
                                                                       * the arrays to the first nonzero belonging to this block and this variable */
            }
         }
      }

   }

   /* start inserting the constant matrix */
   if ( sdpconstnnonz > 0 )
   {

#ifndef NDEBUG
      if (sdpconstnnonz > 0)
      {
         assert ( nsdpblocks > 0 );
         assert ( sdpconstnblocknonz!= NULL );
         assert ( sdpconstcol != NULL );
         assert ( sdpconstrow != NULL );
         assert ( sdpconstval != NULL );
      }
#endif

      /*allocate memory*/

      /* DSDP uses these for solving, so they may not be freed before the problem is solved. */

      /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz));
      /* values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz));

      ind = 0;

      for(block = 0; block < nsdpblocks; block++)
      {
         startind = ind; /* starting index of this block in the dsdpconst arrays */

         /* insert the constant-nonzeros */
         for(i = 0; i < sdpconstnblocknonz[block]; i++)
         {
            assert ( indchanges[block][sdpconstrow[block][i]] > -1 && indchanges[block][sdpconstcol[block][i]] > -1 ); /* rows and cols with
                                                                                                               *  nonzeros should not be removed */
            dsdpconstind[ind] = compLowerTriangPos(sdpconstrow[block][i] - indchanges[block][sdpconstrow[block][i]],
                                                   sdpconstcol[block][i] - indchanges[block][sdpconstcol[block][i]]); /* substract the number of
                                                   * deleted indices before this to get the index after variable fixings */
            dsdpconstval[ind] = -1 * sdpconstval[block][i]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j)
                                                             * should be positive semidefinite */
            ind++;
         }

         /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
         SCIPsortIntReal(dsdpconstind + startind, dsdpconstval + startind, sdpconstnblocknonz[block]);

      }

      DSDP_CALL(SDPConeSetASparseVecMat(sdpisolver->sdpcone, block, 0, sdpblocksizes[block], 1.0, 0, dsdpconstind + startind,
            dsdpconstval + startind, ind - startind));   /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block
                                                          * by adding startind, ind - startind gives the number of elements for this block */
   }

#ifdef SCIP_MORE_DEBUG
      SDPConeView2(sdpisolver->sdpcone);
#endif


   /*start inserting the LP constraints */
   if(nlpcons > 0 || lpnnonz > 0)
   {
      int* rowshifts;
      int nshifts;
      int lastrow;
      int lastcol;
      int nremovedlpcons;
      SCIP_Real rowactive;

      assert ( nlpcons > 0 );
      assert ( lprhs != NULL );
      assert ( lpcol != NULL );
      assert ( lprow != NULL );
      assert ( lpval != NULL );

      /* memory allocation */

      /* these arrays are needed in DSDP during solving, so they may only be freed afterwards */
      /* dsdplpbegcol[i] gives the number of nonzeros in column 0 (right hand side) till i-1 (i going from 1 till m, with extra entries 0 (always 0) and m+1 (always lpcons + lpnnonz)) */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2));
      /* dsdplprow saves the row indices of the LP for DSDP */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpcons + lpnnonz)); /* worst-case length is lpnnonz + nlpcons, because
                                     * right hand sides are also included in the vector, this will be shortened after the exact length after fixings is known */
      /* values given to DSDP */
      BMS_CALL(BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz)); /* worst-case length is lpnnonz + nlpcons, because
                                     * right hand sides are also included in the vector, this will be shortened after the exact length after fixings is known */

      /* array to save which rows can be deleted as they are empty, in that case rowshifts[rowind] = -1, otherwise it gives the
       * number of rowindices deleted before it, therefore the index should be decreased by that number */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &rowshifts, nlpcons) );

      /* first sort the arrays by rows to check for empty constraints */
      SCIPsortIntIntReal(lprow, lpcol, lpval, lpnnonz);

      /* iterate over all nonzeros to see if any rows can be deleted, and add the rhs for all rows that weren't removed */
      nshifts = 0;
      lastrow = -1;
      rowactive = TRUE; /* this is just a technical start point to not need another check in the if below, we don't want to do anything for row -1 */

      for (i = 0; i < lpnnonz; i++)
      {
         assert ( 0 <= lpcol[i] && lpcol[i] < nvars );
         assert ( 0 <= lprow[i] && lprow[i] < nlpcons );

         if (lprow[i] > lastrow)
         {
            /* the next row starts */

            /* first we check if the last row (which is now finished) can be deleted */
            if ( ! rowactive)
            {
#ifndef NDEBUG
               /* if there were no active nonzeros, we don't need the lp row, we only check if the rhs is positive (then we have 0 >= rhs > 0 which
                * renders the SDP infeasible) [but actually we test rhs < 0 because dsdp works with <= in the lp part, so we already multiplied
                * everything with -1], this is only done in Debug mode as this should normally be taken care of by SCIP itself */
               if (dsdplpval[lastrow - nshifts] < -(sdpisolver->epsilon))
               {
                  sdpisolver->infeasible = TRUE;
                  return SCIP_OKAY; /* we know the problem is infeasible, so we don't have to solve it */
               }
#endif

               rowshifts[lprow[i]] = -1;
               nshifts++;
            }

            /* we check if there were any rows in between without any nonzeros (if this is just the next row, the for-queue is empty) */
            for (j = lastrow + 1; j < lprow[i]; j++)
            {
#ifndef NDEBUG
               /* as these rows didn't have any nonzeros, they can surely be eliminated, we just check again if rhs > 0 [< 0 for dsdp] */
               if (lprhs[j] < -(sdpisolver->epsilon))
               {
                  sdpisolver->infeasible = TRUE;
                  return SCIP_OKAY; /* we know the problem is infeasible, so we don't have to solve it */
               }
#endif
               rowshifts[j] = -1;
               nshifts++;
            }

            /* now we initialize the new row */
            lastrow = lprow[i];
            if (isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]))
            {
               /* we don't know yet if it is active, so we set the bool to false and compute the rhs, then continue checking the rest of the nonzeros
                * in the next iterations */
               rowactive = FALSE;
               dsdplpval[lastrow - nshifts] = -lprhs[lastrow] + lpval[i]; /* the rhs is multiplied by -1 as dsdp wants <= instead of >=, the + for
                                                                              * lpval comes because of this and rhs - lhs, so - * - = + */
            }
            else
            {
               /* we have found an active nonzero, so this row can't be deleted */
               rowshifts[lprow[i]] = nshifts;
               dsdplpval[lastrow - nshifts] = -lprhs[lastrow]; /* the rhs is multiplied by -1 as dsdp wants <= instead of >= */
               rowactive = TRUE;
            }
         }
         else
         {
            /* as the row index didn't increase, we have another nonzero for the row we are currently handling */
            if (isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]))
            {
               /* we just add the fixed value to the rhs, this doesn't change anything about the activity status */
               dsdplpval[lastrow - nshifts] += lpval[i]; /* + = - * - because of rhs - lhs and <= instead of >= */
            }
            else
            {
               /* we found an active variable, so this will definitly stay active */
               rowactive = TRUE;
               rowshifts[lprow[i]] = nshifts;
            }
         }
      }

      nremovedlpcons = nshifts;


      /* now add the nonzeros */

      /* for this we have to sort the nonzeros by col first and then by row, as this is the sorting DSDP wants */
      sortColRow(lprow, lpcol, lpval, lpnnonz);

      /* iterate over all nonzeros to add the active ones to the dsdp arrays and compute dsdplpbegcol */
      lastcol = -1;
      dsdplpbegcol[0] = 0;
      dsdplpbegcol[1] = nlpcons - nshifts; /* the number of LP-constraints that will be given to dsdp */
      for (i = 0; i < lpnnonz; i++)
      {
         if (sdpisolver->inputtodsdpmapper[lpcol[i]] > lastcol)
         {
            /* set the dsdplpbegcol entries, as there might be active variables which appear only in the sdp but not the lp-part, we also have to set
             * the starting values for all variable in between to the same value (as we also set the entry for the found variable, this for-queue
             * will always have at least one index in the index set */
            for (j = lastcol + 1; j <= sdpisolver->inputtodsdpmapper[lpcol[i]]; j++)
               dsdplpbegcol[j] = nlpcons + i - nshifts; /* the nonzeros only start after the rhs, the are shifted nshift positions to the left */

            /* add the nonzero, if it isn't fixed, otherwise note the needed shift for the rest */
            if (isFixed(sdpisolver, lb[lpcol[i]], ub[lpcol[i]]))
               nshifts++;
            else
            {
               dsdplprow[i - nshifts] = lprow[i] - rowshifts[lprow[i]]; /* the index is adjusted for deleted lp rows, also rows are numbered
                                                                         * 0,...,nlpcons-1 in DSDP, as they are here */
               dsdplpval[i - nshifts] = -lpval[i]; /* - because dsdp wants <= instead of >= constraints */
            }
         }
      }
      dsdplpbegcol[nvars + 1] = nlpcons + lpnnonz - nshifts; /* the length of the dsdplp arrays */

      /* free the memory for the rowshifts-array */
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &rowshifts, nlpcons);

      /* shrink the dsdplp-arrays */
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, nlpcons + lpnnonz, nlpcons + lpnnonz - nshifts) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz, nlpcons + lpnnonz - nshifts) );
      dsdplparraylength = nlpcons + lpnnonz - nshifts;

      /* add the arrays to dsdp */
      DSDP_CALL( LPConeSetData(sdpisolver->lpcone, nlpcons - nremovedlpcons, dsdplpbegcol, dsdplprow, dsdplpval) );
      #ifdef SCIP_MORE_DEBUG
      LPConeView(sdpisolver->lpcone);
      #endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpisolver->sdpcounter);

   DSDP_CALL(DSDPSetGapTolerance(sdpisolver->dsdp, 1e-3));  /* set DSDP's tolerance for duality gap */
   DSDP_CALL(DSDPSetRTolerance(sdpisolver->dsdp, 1e-6));    /* set DSDP's tolerance for the SDP-constraints */


   /* set the penalty parameter */
   if (penaltyParam != 0.0) /* in sdpisolverSolve this is called with an exact 0 */
   {
      DSDPSetPenaltyParameter(sdpisolver->dsdp, penaltyParam);
      DSDPUsePenalty(sdpisolver->dsdp, 1);
   }

   DSDP_CALLM(DSDPSetup(sdpisolver->dsdp));
   DSDP_CALL(DSDPSolve(sdpisolver->dsdp));

   DSDP_CALL(DSDPComputeX(sdpisolver->dsdp)); /*computes X and determines feasibility and unboundedness of the solution */
   sdpisolver->solved = TRUE;

   /*these arrays were used to give information to DSDP and were needed during solving and for computing X, so they may only be freed now*/
   if ( sdpconstnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz);
   }

   if ( sdpnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz);
   }

   if(nlpcons > 0 || lpnnonz > 0)
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, sdpisolver->nactivevars + 2);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplprow, dsdplparraylength);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, dsdplparraylength);
   }

#ifdef SCIP_DEBUG
   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &reason));
   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   switch ( *reason ) /* TODO: perhaps also check for feasibility and call the penalty-method here in that case */
   {
      case DSDP_CONVERGED:
         SCIPdebugMessage("DSDP converged!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_INFEASIBLE_START:
         SCIPdebugMessage("DSDP started with an infeasible point!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_SMALL_STEPS:
         SCIPdebugMessage("Short step lengths created by numerical difficulties prevented progress in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_INDEFINITE_SCHUR_MATRIX:
         SCIPdebugMessage("Schur Matrix in DSDP was indefinite but should have been positive semidefinite!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_MAX_IT:
         SCIPdebugMessage("DSDP reached maximum number of iterations!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_NUMERICAL_ERROR:
         SCIPdebugMessage("A numerical error occured in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_UPPERBOUND:
         SCIPdebugMessage("Dual objective value in DSDP reached upper bound.\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_USER_TERMINATION:
         SCIPdebugMessage("DSDP didn't stop solving, did you?\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case CONTINUE_ITERATING:
         SCIPdebugMessage("DSDP wants to continue iterating but somehow was stopped!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      default:
         SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;
   }

#endif

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
   assert ( sdpisolver != NULL );
   return sdpisolver->solved;
}

/** returns true if the solver could determine whether or not the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver doesn't know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be used */
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   assert ( primalfeasible != NULL );
   assert ( dualfeasible != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      *primalfeasible = FALSE;
      *dualfeasible = FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));

   switch ( *pdfeasible)
   {
      case DSDP_PDFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      case DSDP_UNBOUNDED:
         *primalfeasible = FALSE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      case DSDP_INFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = FALSE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      default: /* should only include DSDP_PDUNKNOWN */
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
         SCIPABORT();
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
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
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
 */
SCIP_Bool SCIPsdpiSolverHasPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
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
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
SCIP_Bool SCIPsdpiSolverHasDualRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible));
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &reason));

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_CONVERGED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &reason));

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_UPPERBOUND)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return FALSE;
   }
#endif

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &reason));

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_MAX_IT)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return SCIP_ERROR;
}

/** returns the internal solution status of the solver */
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return -1;
   }
#endif

   if (sdpisolver->dsdp == NULL)
      return -1;

   BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &reason));

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   switch ( *reason)
   {
      case DSDP_CONVERGED:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 0;
      }
      case DSDP_INFEASIBLE_START:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 1;
      }
      case DSDP_SMALL_STEPS:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_INDEFINITE_SCHUR_MATRIX:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_MAX_IT:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 4;
      }
      case DSDP_NUMERICAL_ERROR:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_UPPERBOUND:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 3;
      }
      case DSDP_USER_TERMINATION:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 6;
      }
      default:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 7;
      }
   }
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   return (SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverIsPrimalFeasible(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 * that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving.");
      return TRUE;
   }
#endif

   if (SCIPsdpiSolverIsConverged(sdpisolver))
   {
      return TRUE;
   }
   else
   {
      double* pobj;
      double* dobj;
      double gap;

      printf("Numerical Trouble in DSDP!\n");

      /* if it didn't converge check the optimality gap */
      BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &pobj));
      BMS_CALL(BMSallocBlockMemory(sdpisolver->blkmem, &dobj));

      DSDP_CALL(DSDPGetPObjective(sdpisolver->dsdp, pobj));
      DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, dobj));

      gap = abs(*pobj - *dobj);

      if ((gap < sdpisolver->epsilon) || ((gap / (0.5 * (abs(*pobj) + abs(*dobj)))) < sdpisolver->epsilon)) /* this is the duality gap used in SDPA */
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &pobj);
         BMSfreeBlockMemory(sdpisolver->blkmem, &dobj);
         return TRUE;
      }
      else
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &pobj);
         BMSfreeBlockMemory(sdpisolver->blkmem, &dobj);
         return FALSE;
      }
   }
/* TODO: also check for primal feasibility, as this is also needed for optimality */
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   int v;

   assert ( sdpisolver != NULL );
   assert ( objval != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, objval));
   *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

   /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective by hand */
   for (v = 0; v < sdpisolver->nvars - sdpisolver->nactivevars; v++)
   {
      if (sdpisolver->inputtodsdpmapper[v] == -1)
         *objval += sdpisolver->fixedvarsobj[v] * sdpisolver->fixedvarsval[v];
   }

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return , the needed length and
 *  a debug message */
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                               *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   )
{
   SCIP_Real* dsdpsol;
   int v;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   if ( objval != NULL )
   {
      DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, objval));
      *objval *= -1; /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

      /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective by hand */
      for (v = 0; v < sdpisolver->nvars - sdpisolver->nactivevars; v++)
      {
         if (sdpisolver->inputtodsdpmapper[v] == -1)
            *objval += sdpisolver->fixedvarsobj[v] * sdpisolver->fixedvarsval[v];
      }
   }

   if (*dualsollength > 0)
   {
      assert(dualsol != NULL);
      if (*dualsollength < sdpisolver->nvars)
      {
         SCIPdebugMessage("The given array in SCIPsdpiSolverGetSol only had length %d, but %d was needed", *dualsollength, sdpisolver->nvars);
         *dualsollength = sdpisolver->nvars;
      }

      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, sdpisolver->nactivevars) );
      DSDP_CALL(DSDPGetY(sdpisolver->dsdp, dsdpsol, sdpisolver->nactivevars)); /*last entry needs to be the number of variables, will return an error otherwise */

      /* insert the entries into dualsol, for non-fixed vars we copy those from dsdp, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if (sdpisolver->inputtodsdpmapper[v] > -1)
            dualsol[v] = dsdpsol[sdpisolver->inputtodsdpmapper[v] - 1]; /* minus one because the inputtosdpmapper gives the dsdp indices which start at one,
                                                                      * but the array starts at 0 */
         else
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtodsdpmapper[v]) - 1]; /* this is the value that was saved when inserting, as this variable has lb=ub */
      }
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, sdpisolver->nactivevars);
   }
   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

#ifndef NDEBUG
   if (sdpisolver->infeasible)
   {
      SCIPdebugMessage("Problem wasn't given to solver as dual infeasibility was detected during insertion/presolving, so no solution exists.");
      return SCIP_OKAY;
   }
#endif

   DSDP_CALL(DSDPGetIts(sdpisolver->dsdp, iterations));
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
   SCIP_Real            val                 /**< value to be checked for infinity */
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
   SCIP_Real            val                 /**< value to be compared to maximum penalty parameter */
   )
{
   return ((val <= -SCIPsdpiSolverMaxPenParam(sdpisolver)) || (val >= SCIPsdpiSolverMaxPenParam(sdpisolver)));
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
   return SCIP_ERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */
