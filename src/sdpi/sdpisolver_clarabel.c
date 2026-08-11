/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2025 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2025 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpisolver_clarabel.c
 * @brief  interface for Clarabel
 * @author Marc Pfetsch
 *
 * For solving the following problem (the penalty formulation)
 *   \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r\\
 *      \mbox{s.t.} & & \sum_{i \in I} A_i^{(k)} y_i - A_0^{(k)} + I r \succeq 0 \quad \forall \ k \in K, \\
 *      & & \sum_{i \in I} d_{ij} y_i \geq c_j \quad \forall \ j \in J, \\
 *      & & \ell_i \leq y_i \leq u_i \quad \forall \ i \in I,\ r \geq 0,
 *   \f}
 * we use the following problem for Clarabel:
 *   \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r\\
 *      \mbox{s.t.} & & - \sum_{i \in I} A_i^{(k)} y_i - I r + Z^{(k)} = - A_0^{(k)} \quad \forall \ k \in K, \\
 *      & & - \sum_{i \in I} d_{ij} y_i + z_j = - c_j \quad \forall \ j \in J, \\
 *      & & y_i + s_i = u_i \quad \forall i \in I,\\
 *      & & - y_i + t_i = - \ell_i \quad \forall \ i \in I,\\
 *      & & Z \succeq 0,\\
 *      & & z, r, s, t \geq 0.
 *   \f}
 *
 * Some notes:
 * - Clarabel expects the problem to be in the form \f$Ax + s = b\f$ and \f$s\f$ has to be contained in a cone.
 *   Thus, for each constraint there is a cone to which the corresponding slack variable belongs.
 * - For the semidefinite cone, Clarabel only uses the upper triangle and scales the off-diagonal entries by \f$\sqrt{2}\f$.
 * - Therefore, when writing down \f$- \sum_{i \in I} A_i^{(k)} y_i - I r + Z^{(k)} = - A_0^{(k)}\f$, we also scale the
 *   off-diagonal entries by \f$\sqrt{2}\f$. Similarly, the solutions in the dual problem have to be scaled.
 *
 * @todo Use SCIP message handler.
 * @todo Possibly use clarabelsettings.tol_gap_rel.
 * @todo Check for objective cutoff.
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"
#include "sdpi/sdpiclock.h"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for SCIPsnprintf() */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */
#include "scip/pub_message.h"                /* for debug and error message */

#include <clarabel.h>

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e10     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e6      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */


/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   SCIP_Real             opttime;            /**< time spend in optimziation */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in the dual problem in Clarabel (nvars minus the number of fixed variables) */
   int                   nconss;             /**< total number of constraints in Clarabel SDP */
   int                   nsdpconss;          /**< number of constraints in Clarabel SDP corresponding to SDP part */
   int                   nlpconss;           /**< number of constraints in Clarabel SDP corresponding to LP part */
   int                   maxnvars;           /**< size of the arrays inputtocblmap, cbltoinputmap, fixedvarsval, and objcoefs */
   int*                  inputtocblmap;      /**< entry i gives the index of input variable i in Clarabel (starting from 0) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  cbltoinputmap;      /**< entry i gives the original index of the i-th variable in Clarabel (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   int                   nvarbounds;         /**< number of variable bounds given to Clarabel, length of varboundpos */
   int*                  varboundpos;        /**< maps position of primal variable corresponding to variable bound to the positions
                                              *   of the corresponding variables, -n means lower bound of variable n, +n means upper bound;
                                              *   entry i gives variable bound corresponding to the primal variable in the i-th position
                                              *   of the boundblock */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< tolerance used for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance for the SDP-solver */
   SCIP_Real             objlimit;           /**< objective limit for SDP solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             usepresolving;      /**< Should presolving be used? */
   SCIP_Bool             usescaling;         /**< Should the SDP-solver use scaling? */
   SCIP_Bool             penalty;            /**< Was the problem last solved using a penalty formulation? */
   SCIP_Bool             feasorig;           /**< Was the last problem solved with a penalty formulation and with original objective coefficents
                                              *   and the solution was feasible for the original problem? */
   SCIP_Bool             rbound;             /**< Was the penalty parameter bounded during the last solve call? */
   SCIP_Bool             timelimit;          /**< Was the solver stopped because of the time limit? */
   int                   nthreads;           /**< number of threads the SDP solver should use (-1 = number of cores) */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
   SCIP_Bool             scaleobj;           /**< whether the objective should be scaled */
   SCIP_Real             objscalefactor;     /**< objective scaling factor */
   ClarabelSolverStatus  solstat;            /**< Clarabel solving status */
   ClarabelDefaultSolver* clarabelsolver;    /**< Clarabel solver */
};


/*
 * Local Methods
 */

/* defines for lexicographic sorting macros */
#define SORTTPL_NAMEEXT     UIntPtrUIntPtrReal
#define SORTTPL_KEYTYPE     uintptr_t
#define SORTTPL_FIELD1TYPE  SCIP_Real
#include "../scipsdp/sorttpllex.c" /*lint !e451*/


/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY. */
#define BMS_CALL(x)   do                                                                                     \
                      {                                                                                      \
                         if( NULL == (x) )                                                                   \
                         {                                                                                   \
                            SCIPerrorMessage("No memory in function call.\n");                               \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )
/** This will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem. */
#define CHECK_IF_SOLVED(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** This is the same as CHECK_IF_SOLVED, but will be called for methods returning a bool instead of a SCIP_RETURNCODE */
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                      \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            SCIPABORT();                                                                     \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )


#ifndef NDEBUG
/** Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( sdpisolver != NULL );
   assert( lb < ub + sdpisolver->feastol );

   return (ub-lb <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (ub-lb <= sdpisolver->epsilon)
#endif

/** calculate memory size for dynamically allocated arrays */
static
int calcGrowSize(
   int                   initsize,           /**< initial size of array */
   int                   num                 /**< minimum number of entries to store */
   )
{
   int oldsize;
   int size;

   assert( initsize >= 0 );
   assert( num >= 0 );

   /* calculate the size with loop, such that the resulting numbers are always the same (-> block memory) */
   initsize = MAX(initsize, SCIP_DEFAULT_MEM_ARRAYGROWINIT);
   size = initsize;
   oldsize = size - 1;

   /* second condition checks against overflow */
   while ( size < num && size > oldsize )
   {
      oldsize = size;
      size = (int)(SCIP_DEFAULT_MEM_ARRAYGROWFAC * size + initsize);
   }

   /* if an overflow happened, set the correct value */
   if ( size <= oldsize )
      size = num;

   assert( size >= initsize );
   assert( size >= num );

   return size;
}

/** ensure size of mapping data */
static
SCIP_RETCODE ensureMappingDataMemory(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver structure */
   int                   nvars               /**< number of variables */
   )
{
   int newsize;

   assert( sdpisolver != NULL );

   if ( nvars > sdpisolver->maxnvars )
   {
      newsize = calcGrowSize(sdpisolver->maxnvars, nvars);

      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->varboundpos), 2 * sdpisolver->maxnvars, 2 * newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtocblmap), sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->cbltoinputmap), sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), sdpisolver->maxnvars, newsize) );
      sdpisolver->maxnvars = newsize;
   }

   return SCIP_OKAY;
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

char solvername[SCIP_MAXSTRLEN];

/** gets name and version (if available) of SDP-solver */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "Clarabel";
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Interior-point solver for conic problems by P. Goulart and Y. Chen (https://clarabel.org/)";
}

/** gets pointer to SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->clarabelsolver;
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return 8;
}

/** Should primal solution values be saved for warmstarting purposes? */
SCIP_Bool SCIPsdpiSolverDoesWarmstartNeedPrimal(
   void
   )
{
   return FALSE;
}

/**@} */


/*
 * SDP Solver Interface Creation and Destruction Methods
 */

/**@name SDP Solver InterfaceCreation and Destruction Methods */
/**@{ */

/** creates an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   )
{
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );
   assert( bufmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiSolverCreate.\n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

   (*sdpisolver)->opttime = 0.0;
   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->nconss = 0;
   (*sdpisolver)->nsdpconss = 0;
   (*sdpisolver)->nlpconss = 0;
   (*sdpisolver)->inputtocblmap = NULL;
   (*sdpisolver)->cbltoinputmap = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->objcoefs = NULL;
   (*sdpisolver)->maxnvars = 0;
   (*sdpisolver)->nvarbounds = 0;
   (*sdpisolver)->varboundpos = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-4;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->usepresolving = TRUE;
   (*sdpisolver)->usescaling = TRUE;
   (*sdpisolver)->nthreads = -1;
   (*sdpisolver)->solstat = ClarabelUnsolved;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->scaleobj = FALSE;
   (*sdpisolver)->objscalefactor = 1.0;

   (*sdpisolver)->solstat = ClarabelUnsolved;
   (*sdpisolver)->clarabelsolver = NULL;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtocblmap, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->cbltoinputmap, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->maxnvars);

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   SCIPdebugMessage("Resetting counter of SDP-Interface from %d to 0.\n", sdpisolver->sdpcounter);
   sdpisolver->sdpcounter = 0;

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
 *  @note starting point needs to be given with original indices (before any local presolving), last block should be the LP block with indices
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independent of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective coefficients of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   const int*            sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int* const*           sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int* const*           sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real* const*     sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int* const*           sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int* const*           sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int** const*          sdprow,             /**< pointer to the row-indices for each block and variable */
   int** const*          sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real** const*    sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int* const*           indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   const int*            nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   const int*            blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number LP-constraints */
   const int*            lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   const SCIP_Real*      lplhs,              /**< left-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   const SCIP_Real*      lprhs,              /**< right-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   const int*            lpbeg,              /**< start index of each row in ind- and val-array */
   const int*            lpind,              /**< column-index for each entry in lpval-array */
   const SCIP_Real*      lpval,              /**< values of LP-constraint matrix entries */
   const SCIP_Real*      starty,             /**< NULL or dual vector y as starting point for the solver, this should have length nvars */
   const int*            startZnblocknonz,   /**< dual matrix Z = sum Ai yi as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int* const*           startZrow,          /**< dual matrix Z = sum Ai yi as starting point for the solver: row indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int* const*           startZcol,          /**< dual matrix Z = sum Ai yi as starting point for the solver: column indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   SCIP_Real* const*     startZval,          /**< dual matrix Z = sum Ai yi as starting point for the solver: values for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   const int*            startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int* const*           startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   int* const*           startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_Real* const*     startXval,          /**< primal matrix X as starting point for the solver: values for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP) */
   SDPI_CLOCK*           usedsdpitime        /**< clock to measure how much time has been used for the current solve */
   )
{
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, FALSE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
      sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
      nremovedinds, blockindchanges, nremovedblocks, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval,
      starty, startZnblocknonz, startZrow, startZcol, startZval, startXnblocknonz, startXrow, startXcol, startXval, startsettings,
      timelimit, usedsdpitime, NULL, NULL);
}

/** loads and solves an SDP using a penalty formulation
 *
 *  The penalty formulation of the SDP is:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy + r \cdot \mathbb{I} \geq d \\
 *      & & l \leq y \leq u \\
 *      & & r \geq 0.\f}
 *  Alternatively withobj can be set to false to set b to 0 and only check for feasibility (if the optimal objective value is
 *  bigger than 0 the problem is infeasible, otherwise it's feasible), and rbound can be set to false to remove the non-negativity condition on r.
 *  For the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the constant SDP-part the arrays AFTER fixings
 *  should be given. In addition, an array needs to be given, that for every block and every row/col index within that block either has value
 *  -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices before it that are to be deleated,
 *  meaning that this index will be decreased by that number. Moreover, the total number of deleted indices for each block should be given.
 *  An optional starting point for the solver may be given; if it is NULL, the solver will start from scratch.
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 *  @note starting point needs to be given with original indices (before any local presolving), last block should be the LP block with indices
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independent of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */ /*lint --e{715}*/
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< Should r be non-negative? */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective coefficients of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   const int*            sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int* const*           sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int* const*           sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real* const*     sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int* const*           sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int* const*           sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int** const*          sdprow,             /**< pointer to the row-indices for each block and variable */
   int** const*          sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real** const*    sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int* const*           indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   const int*            nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   const int*            blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of LP-constraints */
   const int*            lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   const SCIP_Real*      lplhs,              /**< left-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   const SCIP_Real*      lprhs,              /**< right-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   const int*            lpbeg,              /**< start index of each row in ind- and val-array */
   const int*            lpind,              /**< column-index for each entry in lpval-array */
   const SCIP_Real*      lpval,              /**< values of LP-constraint matrix entries */
   const SCIP_Real*      starty,             /**< NULL or dual vector y as starting point for the solver, this should have length nvars */
   const int*            startZnblocknonz,   /**< dual matrix Z = sum Ai yi as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int* const*           startZrow,          /**< dual matrix Z = sum Ai yi as starting point for the solver: row indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int* const*           startZcol,          /**< dual matrix Z = sum Ai yi as starting point for the solver: column indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   SCIP_Real* const*     startZval,          /**< dual matrix Z = sum Ai yi as starting point for the solver: values for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   const int*            startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int* const*           startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   int* const*           startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_Real* const*     startXval,          /**< primal matrix X as starting point for the solver: values for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP) */
   SDPI_CLOCK*           usedsdpitime,       /**< clock to measure how much time has been used for the current solve */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                              *   (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                              *   this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
   )
{/*lint --e{715}*/
   SCIP_Real maxabsobjcoef = 0.0;
   SCIP_Real solvertimelimit;
   SCIP_Real val;
   SCIP_Real sqrt2;
   int blockvar;
   int nfixedvars;
   int sumblocksizes;
   int b;
   int i;
   int j;
   int v;
   int k;
   int ind;
   int row;
   int col;
   int cnt;
   uintptr_t pos;

   /* Clarabel specific variables: */
   ClarabelDefaultSettings clarabelsettings;
   ClarabelCscMatrix clarabelmatrix;
   ClarabelSupportedConeT* clarabelcones;
   ClarabelDefaultInfo clarabelinfo;
   ClarabelCscMatrix clarabelquadraticmatrix;
   SCIP_Real* clarabelvals;
   SCIP_Real* clarabelbounds;
   SCIP_Real* clarabelobj;
   SCIP_Real* clarabelrhs;
   int* clarabelblocksizes;
   int* clarabelblockstart;
   uintptr_t* clarabelrows;
   uintptr_t* clarabelcols;
   uintptr_t* clarabelcolptr;
   uintptr_t* clarabelquadraticcolptr = NULL;
   int clarabellpnnonz;
   int clarabelsdpnnonz;
   int clarabelnnonz;
   int clarabelnvars;
   int clarabelnsdpconss;
   int clarabelnlpconss;
   int clarabelnconss;
   int clarabelncones;
   int clarabelstart;

   assert( sdpisolver != NULL );
   assert( penaltyparam > -1 * sdpisolver->epsilon );
   assert( penaltyparam < sdpisolver->epsilon || ( feasorig != NULL ) );
   assert( nvars > 0 );
   assert( obj != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nsdpblocks >= 0 );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( nsdpblocks == 0 || sdpnblockvars != NULL );
   assert( sdpconstnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpnblockvarnonz != NULL );
   assert( nsdpblocks == 0 || sdpvar != NULL );
   assert( nsdpblocks == 0 || sdprow != NULL );
   assert( nsdpblocks == 0 || sdpcol != NULL );
   assert( nsdpblocks == 0 || sdpval != NULL );
   assert( nsdpblocks == 0 || indchanges != NULL );
   assert( nsdpblocks == 0 || nremovedinds != NULL );
   assert( nsdpblocks == 0 || blockindchanges != NULL );
   assert( 0 <= nremovedblocks && nremovedblocks <= nsdpblocks );
   assert( nlpcons >= 0 );
   assert( nlpcons == 0 || lpindchanges != NULL );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( nlpcons == 0 || lpbeg != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lpind != NULL );
   assert( nlpcons == 0 || lpval != NULL );

   /* compute the time limit to set for the solver */
   solvertimelimit = timelimit;
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
      solvertimelimit -= SDPIclockGetTime(usedsdpitime);

   sdpisolver->niterations = 0;
   sdpisolver->nsdpcalls = 0;
   sdpisolver->objscalefactor = 1.0;

   /* check the timelimit */
   if ( solvertimelimit <= 0.0 )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->solved = FALSE;
      return SCIP_OKAY;
   }
   sdpisolver->timelimit = FALSE;
   sdpisolver->feasorig = FALSE;
   sqrt2 = sqrt(2);

   /* set the penalty and rbound flags accordingly */
   if ( penaltyparam < sdpisolver->epsilon )
   {
      sdpisolver->penalty = FALSE;
      sdpisolver->rbound = FALSE;
   }
   else
   {
      sdpisolver->penalty = TRUE;
      sdpisolver->rbound = rbound;
   }

   /* Increase SDP-counter if we do not use the penalty formulation because for the outside this is still the same SDP. */
   if ( sdpisolver->penalty )
   {
      SCIPdebugMessage("Inserting Data again into Clarabel for penalty formulation of SDP (%d).\n", sdpisolver->sdpcounter);
   }
   else
   {
      ++sdpisolver->sdpcounter;
      SCIPdebugMessage("Inserting data into Clarabel for SDP (%d).\n", sdpisolver->sdpcounter);
   }

   /* ensure memory for varboundpos, inputtocblmap, cbltoinputmap and the fixed and active variable information */
   SCIP_CALL( ensureMappingDataMemory(sdpisolver, nvars) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelbounds, 2 * nvars) ); /*lint !e647*/

   /* find fixed variables */
   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   sdpisolver->nvarbounds = 0;
   sdpisolver->fixedvarsobjcontr = 0.0;
   nfixedvars = 0;
   for (i = 0; i < nvars; ++i)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* value that this fixed variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars] = lb[i];    /* if lb=ub, then this is the value the variable will have in every solution */
         ++nfixedvars;
         sdpisolver->inputtocblmap[i] = -nfixedvars;
         SCIPdebugMessage("Fixing variable %d locally to %g for SDP %d in Clarabel.\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in Clarabel.\n", i, sdpisolver->nactivevars, sdpisolver->sdpcounter);
#endif

         sdpisolver->cbltoinputmap[sdpisolver->nactivevars] = i;
         sdpisolver->inputtocblmap[i] = sdpisolver->nactivevars;
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];

         if ( withobj )
         {
            if ( REALABS(obj[i]) > maxabsobjcoef )
               maxabsobjcoef = REALABS(obj[i]);
         }

         if ( lb[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            clarabelbounds[sdpisolver->nvarbounds] = - lb[i]; /* -1 because the standard is <= */
            sdpisolver->varboundpos[sdpisolver->nvarbounds++] = -(sdpisolver->nactivevars + 1); /* negative sign means lower bound */
         }

         if ( ub[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            clarabelbounds[sdpisolver->nvarbounds] = ub[i];
            sdpisolver->varboundpos[sdpisolver->nvarbounds++] = +(sdpisolver->nactivevars + 1); /* positive sign means upper bound */
         }

         ++sdpisolver->nactivevars;
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* adjust maxabsobjcoef in penalty formulation */
   if ( sdpisolver->penalty )
   {
      if ( penaltyparam > maxabsobjcoef )
         maxabsobjcoef = penaltyparam;
   }

   /* determine total number of sides in LP-constraints and nonzeros */
   clarabellpnnonz = 0;
   clarabelnlpconss = 0;
   if ( nlpcons > 0 )
   {
      for (i = 0; i < nlpcons; ++i)
      {
         int delta = 0;  /* count how many nonzeros are added depending on lhs/rhs */

         if ( lpindchanges[i] < 0 )
            continue;

         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            ++clarabelnlpconss;
            ++delta;
         }

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            ++clarabelnlpconss;
            ++delta;
         }

         /* the following overestimates the true number a bit, because of fixed variables */
         if ( i == nlpcons - 1 )
            clarabellpnnonz += delta * (lpnnonz - lpbeg[i]);
         else
            clarabellpnnonz += delta * (lpbeg[i+1] - lpbeg[i]);
      }
   }
   assert( clarabellpnnonz <= 2 * lpnnonz ); /* factor 2 comes from left- and right-hand-sides */

   /* compute SDP sizes for Clarabel model */
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelblocksizes, nsdpblocks - nremovedblocks) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelblockstart, nsdpblocks) );

   /* compute number of constraints for the Clarabel matrix in the following */
   clarabelnsdpconss = 0;
   clarabelsdpnnonz = 0;
   sumblocksizes = 0;
   for (b = 0; b < nsdpblocks; ++b)
   {
      clarabelblockstart[b] = clarabelnsdpconss;
      if ( blockindchanges[b] > -1 )
      {
         int size;
         int dim;

         assert( 0 <= blockindchanges[b] && blockindchanges[b] <= b && (b - blockindchanges[b]) <= (nsdpblocks - nremovedblocks) );
         dim = sdpblocksizes[b] - nremovedinds[b];
         size = dim * (dim + 1) / 2;

         clarabelblocksizes[b - blockindchanges[b]] = size; /*lint !e679*/
         clarabelnsdpconss += size;
         sumblocksizes += dim;

         for (blockvar = 0; blockvar < sdpnblockvars[b]; ++blockvar)
         {
            if ( sdpisolver->inputtocblmap[sdpvar[b][blockvar]] >= 0 )
               clarabelsdpnnonz += sdpnblockvarnonz[b][blockvar];
         }

         if ( sdpisolver->penalty )
            clarabelsdpnnonz += dim;
      }
   }

   /* determine number of variables */
   clarabelnvars = sdpisolver->nactivevars;   /* original variables are preserved */
   if ( sdpisolver->penalty )
      ++clarabelnvars;  /* for penalty variable r */

   /* determine number of linear constraints (equations) */
   clarabelnconss = clarabelnsdpconss;
   clarabelnconss += clarabelnlpconss;
   clarabelnconss += sdpisolver->nvarbounds;
   if ( sdpisolver->rbound )
   {
      assert( sdpisolver->penalty );
      ++clarabelnconss;
   }

   sdpisolver->nconss = clarabelnconss;
   sdpisolver->nsdpconss = clarabelnsdpconss;
   sdpisolver->nlpconss = clarabelnlpconss;

   /* create rhs vector */
   BMS_CALL( BMSallocClearBufferMemoryArray(sdpisolver->bufmem, &clarabelrhs, clarabelnconss) );

   /* set rhs values for the SDP constraints */
   if ( sdpconstnnonz > 0 )
   {
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] >= 0 )
         {
            /* if some indices in the block were removed, we need to change indices accordingly */
            for (k = 0; k < sdpconstnblocknonz[b]; k++)
            {
               /* rows and cols with active nonzeros should not be removed */
               assert( -1 < indchanges[b][sdpconstrow[b][k]] && indchanges[b][sdpconstrow[b][k]] <= sdpconstrow[b][k] );
               assert( -1 < indchanges[b][sdpconstcol[b][k]] && indchanges[b][sdpconstcol[b][k]] <= sdpconstcol[b][k] );

               assert( 0 <= sdpconstrow[b][k] && sdpconstrow[b][k] <= sdpblocksizes[b] );
               assert( 0 <= sdpconstcol[b][k] && sdpconstcol[b][k] <= sdpblocksizes[b] );

               row = sdpconstrow[b][k] - indchanges[b][sdpconstrow[b][k]];
               col = sdpconstcol[b][k] - indchanges[b][sdpconstcol[b][k]];
               ind = clarabelblockstart[b] + row * (row + 1)/2 + col;

               assert( 0 <= ind && ind < clarabelnconss );

               if ( row == col )
                  clarabelrhs[ind] = - sdpconstval[b][k];  /* -1 because we move to the right hand side */
               else
                  clarabelrhs[ind] = - sdpconstval[b][k] * sqrt2;  /* -1 because we move to the right hand side */
            }
         }
      }
   }

   /* set rhs values for the LPs */
   ind = 0;
   clarabelstart = clarabelnsdpconss;
   for (i = 0; i < nlpcons; i++)
   {
      if ( lpindchanges[i] < 0 )
         continue;

      /* left hand side */
      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
      {
         clarabelrhs[clarabelstart + ind] = - lplhs[i];  /* -1 because of <= conversion */
         assert( clarabelstart + ind < clarabelnconss );
         ++ind;
      }

      /* right hand side */
      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         clarabelrhs[ind + clarabelstart] = lprhs[i];
         assert( clarabelstart + ind < clarabelnconss );
         ++ind;
      }
   }
   assert( ind == clarabelnlpconss );

   /* add rhs values for variable bounds */
   clarabelstart += clarabelnlpconss;
   assert( sdpisolver->nvarbounds + clarabelstart == clarabelnconss - sdpisolver->rbound );
   for (i = 0; i < sdpisolver->nvarbounds; i++)
      clarabelrhs[clarabelstart + i] = clarabelbounds[i];
   /* rhs for penalty variable (if rbound is true) is automatically 0 */

   /* possibly scale objective */
   if ( sdpisolver->scaleobj )
   {
      if ( REALABS(maxabsobjcoef) > 1.0 )
      {
         assert( withobj );
         sdpisolver->objscalefactor = maxabsobjcoef;
         SCIPdebugMessage("Scaling objective by %g.\n", 1.0 / sdpisolver->objscalefactor);
         maxabsobjcoef = 1.0;  /* this is now 1 because of scaling */
      }
   }
   else
      assert( sdpisolver->objscalefactor == 1.0 );

   /* create Clarabel objective */
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelobj, clarabelnvars) );
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      SCIP_Real objcoef = 0.0;
      if ( withobj )
      {
         objcoef = sdpisolver->objcoefs[i];
         assert( objcoef == obj[sdpisolver->cbltoinputmap[i]] );
         objcoef /= sdpisolver->objscalefactor;
      }
      clarabelobj[i] = objcoef;
   }

   /* the penalty constraint has objective coefficient Gamma */
   if ( sdpisolver->penalty )
      clarabelobj[sdpisolver->nactivevars] = penaltyparam / sdpisolver->objscalefactor;
   assert( clarabelnvars == sdpisolver->nactivevars + sdpisolver->penalty ? 1 : 0 );

   /* start inserting the constraints: we first generate the matrix in row format and convert later to column format */

   /* compute number of nonzeros in matrix */
   clarabelnnonz = clarabelsdpnnonz;
   clarabelnnonz += clarabellpnnonz;          /* add LP nonzeros */
   clarabelnnonz += sdpisolver->nvarbounds;   /* add variable bounds */
   if ( sdpisolver->penalty )
   {
      clarabelnnonz += sumblocksizes;         /* add one entry for each block dimension for identity matrix */
      if ( sdpisolver->rbound )
         ++clarabelnnonz;                     /* add one entry for nonnegativity bound */
   }

   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelrows, clarabelnnonz) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelcols, clarabelnnonz) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelvals, clarabelnnonz) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &clarabelcolptr, clarabelnvars + 1) );

   /* construct the matrix */
   cnt = 0;
   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] >= 0 )
      {
         for (blockvar = 0; blockvar < sdpnblockvars[b]; blockvar++)
         {
            v = sdpisolver->inputtocblmap[sdpvar[b][blockvar]];

            /* check if the variable is active */
            if ( v > -1 )
            {
               assert( v < sdpisolver->nactivevars );

               assert( sdpnblockvarnonz[b][blockvar] <= sdpnnonz );
               for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( 0 <= indchanges[b][sdprow[b][blockvar][k]] && indchanges[b][sdprow[b][blockvar][k]] <= sdprow[b][blockvar][k] );
                  assert( 0 <= indchanges[b][sdpcol[b][blockvar][k]] && indchanges[b][sdpcol[b][blockvar][k]] <= sdpcol[b][blockvar][k] );

                  assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                  assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                  row = sdprow[b][blockvar][k] - indchanges[b][sdprow[b][blockvar][k]];
                  col = sdpcol[b][blockvar][k] - indchanges[b][sdpcol[b][blockvar][k]];
                  ind = clarabelblockstart[b] + row * (row + 1)/2 + col;

                  clarabelrows[cnt] = ind;
                  clarabelcols[cnt] = v;
                  if ( row == col )
                     clarabelvals[cnt] = -sdpval[b][blockvar][k];
                  else
                     clarabelvals[cnt] = -sdpval[b][blockvar][k] * sqrt2;
                  assert( cnt < clarabelnnonz );
                  ++cnt;
               }
            }
         }

         /* add the identity matrix corresponding to the penalty variable */
         if ( sdpisolver->penalty )
         {
            int dim;

            dim = sdpblocksizes[b] - nremovedinds[b];
            for (k = 0; k < dim; ++k)
            {
               clarabelrows[cnt] = clarabelblockstart[b] + k * (k + 1) / 2;
               clarabelcols[cnt] = clarabelnvars - 1;
               clarabelvals[cnt] = -1.0;
               assert( cnt < clarabelnnonz );
               ++cnt;
            }
         }
      }
   }
   assert( cnt <= clarabelsdpnnonz );
   assert( cnt <= clarabelnnonz );

   /* add the LP constraints */
   if ( lpnnonz > 0 )
   {
      ind = 0;
      for (i = 0; i < nlpcons; ++i)
      {
         int nextbeg;
         int delta = 0;

         if ( lpindchanges[i] < 0 )
            continue;

         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
            ++delta;

         assert( 0 <= lpbeg[i] && lpbeg[i] < lpnnonz );
         if ( i == nlpcons - 1 )
            nextbeg = lpnnonz;
         else
            nextbeg = lpbeg[i+1];

         for (j = lpbeg[i]; j < nextbeg; ++j)
         {
            assert( 0 <= lpind[j] && lpind[j] < nvars );
            v = sdpisolver->inputtocblmap[lpind[j]];
            if ( v >= 0 )
            {
               assert( v < sdpisolver->nactivevars );
               if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
               {
                  clarabelrows[cnt] = clarabelnsdpconss + ind;
                  clarabelcols[cnt] = v;
                  clarabelvals[cnt] = - lpval[j];  /* -1 because of <= conversion */
                  assert( clarabelcols[cnt] < clarabelnvars );
                  assert( cnt < clarabelnnonz );
                  ++cnt;
               }

               if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
               {
                  clarabelrows[cnt] = clarabelnsdpconss + ind + delta; /* delta = 1 if lhs != 0 and 0 otherwise */
                  clarabelcols[cnt] = v;
                  clarabelvals[cnt] = lpval[j];
                  assert( clarabelcols[cnt] < clarabelnvars );
                  assert( cnt < clarabelnnonz );
                  ++cnt;
               }
            }
         }
         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
            ++ind;
         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
            ++ind;
      }
      assert( ind == clarabelnlpconss );
   }

   /* finally add the entries corresponding to varbounds */
   ind = 0;
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* lower bound */
         v = - sdpisolver->varboundpos[i] - 1;
         assert( 0 <= v && v < sdpisolver->nactivevars );
         clarabelrows[cnt] = clarabelnsdpconss + clarabelnlpconss + ind;
         clarabelcols[cnt] = v;
         clarabelvals[cnt] = -1.0; /* -1 because of <= conversion */
         assert( clarabelcols[cnt] < clarabelnvars );
         assert( cnt < clarabelnnonz );
         ++cnt;
         ++ind;
      }
      else
      {
         /* upper bound */
         v = sdpisolver->varboundpos[i] - 1;
         assert( 0 <= v && v < sdpisolver->nactivevars );
         clarabelrows[cnt] = clarabelnsdpconss + clarabelnlpconss + ind;
         clarabelcols[cnt] = v;
         clarabelvals[cnt] = 1.0;
         assert( clarabelcols[cnt] < clarabelnvars );
         assert( cnt < clarabelnnonz );
         ++cnt;
         ++ind;
      }
   }
   assert( ind == sdpisolver->nvarbounds );
   assert( cnt <= clarabelnnonz );

   if ( sdpisolver->rbound )
   {
      assert( sdpisolver->penalty );
      clarabelrows[cnt] = clarabelnconss - 1;
      clarabelcols[cnt] = clarabelnvars - 1;
      clarabelvals[cnt] = -1.0; /* -1 because of <= conversion */
      assert( cnt < clarabelnnonz );
      ++cnt;
   }

   /* resort according to columns and then rows */
   SCIPlexSortUIntPtrUIntPtrReal(clarabelcols, clarabelrows, clarabelvals, cnt);

   /* compute beginnings of new columns */
   pos = 0;
   for (k = 0; k < clarabelnvars; ++k)
   {
      clarabelcolptr[k] = pos;
      while ( pos < cnt && clarabelcols[pos] == k )
         ++pos;
   }
   assert( pos == cnt );
   clarabelcolptr[clarabelnvars] = cnt;

   /* possibly free previous solver */
   if ( sdpisolver->clarabelsolver != NULL )
   {
      clarabel_DefaultSolver_free(sdpisolver->clarabelsolver);
      sdpisolver->clarabelsolver = NULL;
   }

   solvertimelimit = timelimit;
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
      solvertimelimit -= SDPIclockGetTime(usedsdpitime);

   if ( solvertimelimit <= 0.0 )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->solved = FALSE;
   }
   else
   {
      SCIP_Real feastol;
      SCIP_Real gaptol;

      /* construct matrix */
      clarabel_CscMatrix_init(&clarabelmatrix, clarabelnconss, clarabelnvars, clarabelcolptr, clarabelrows, clarabelvals);

      /* get default settings */
      clarabelsettings = clarabel_DefaultSettings_default();

#ifdef SCIP_DEBUG
      /* turn on output */
      sdpisolver->sdpinfo = TRUE;
#endif

      /* set output */
      if ( sdpisolver->sdpinfo )
         clarabelsettings.verbose = true;
      else
         clarabelsettings.verbose = false;

      /* set number of threads */
      if ( sdpisolver->nthreads > 0 )
         clarabelsettings.max_threads = sdpisolver->nthreads;
      else
         clarabelsettings.max_threads = 1;

      /* set time limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
         clarabelsettings.time_limit = solvertimelimit;

      /* we want absolute tolerances */
      clarabelsettings.tol_feas = sdpisolver->sdpsolverfeastol;
      clarabelsettings.tol_gap_abs = sdpisolver->gaptol;
      clarabelsettings.tol_infeas_abs = sdpisolver->feastol;

      clarabelsettings.tol_infeas_rel = sdpisolver->feastol / (1.0 + maxabsobjcoef);
      SCIPdebugMessage("tol_infeas_abs = %g, tol_infeas_rel = %g\n", clarabelsettings.tol_infeas_abs, clarabelsettings.tol_infeas_rel);

      /* Clarabel does not seem to have an objective cutoff. */

      /* turn presolving on/off */
      if ( sdpisolver->usepresolving )
      {
         clarabelsettings.presolve_enable = true;
         SCIPdebugMessage("Turning presolving on.\n");
      }
      else
      {
         clarabelsettings.presolve_enable = false;
         SCIPdebugMessage("Turning presolving off.\n");
      }

      /* turn scaling on/off */
      if ( sdpisolver->usescaling )
         clarabelsettings.equilibrate_enable = true;
      else
         clarabelsettings.equilibrate_enable = false;

      /* set up cones: one cone for each SDP block (SDP) and one for the LP (nonnegative) */
      clarabelncones = nsdpblocks - nremovedblocks + 1;
      BMS_CALL( BMSallocClearBufferMemoryArray(sdpisolver->bufmem, &clarabelcones, clarabelncones) );
      ind = 0;
      for (b = 0; b < nsdpblocks; ++b)
      {
         if ( blockindchanges[b] >= 0 )
            clarabelcones[ind++] = ClarabelPSDTriangleConeT(sdpblocksizes[b] - nremovedinds[b]);
      }
      if ( sdpisolver->rbound )
         clarabelcones[ind++] = ClarabelNonnegativeConeT(clarabelnlpconss + sdpisolver->nvarbounds + 1);
      else
         clarabelcones[ind++] = ClarabelNonnegativeConeT(clarabelnlpconss + sdpisolver->nvarbounds);
      assert( ind == nsdpblocks - nremovedblocks + 1 );

      /* set up empty quadratic matrix */
      BMSallocClearBufferMemoryArray(sdpisolver->bufmem, &clarabelquadraticcolptr, clarabelnvars + 1);
      clarabel_CscMatrix_init(&clarabelquadraticmatrix, clarabelnvars, clarabelnvars, clarabelquadraticcolptr, NULL, NULL);

      /* create solver */
      sdpisolver->clarabelsolver = clarabel_DefaultSolver_new(&clarabelquadraticmatrix, clarabelobj, &clarabelmatrix, clarabelrhs, clarabelncones, clarabelcones, &clarabelsettings);
      assert( sdpisolver->clarabelsolver != NULL );

#ifdef SCIP_MORE_DEBUG
      clarabel_DefaultSolver_save_to_file(sdpisolver->clarabelsolver, "debug.clarabel");
#endif

      /* actually solve SDP */
      clarabel_DefaultSolver_solve(sdpisolver->clarabelsolver);
      clarabelinfo = clarabel_DefaultSolver_info(sdpisolver->clarabelsolver);

      /* update number of SDP-iterations and -calls */
      ++sdpisolver->nsdpcalls;
      sdpisolver->niterations = clarabelinfo.iterations;
      sdpisolver->solstat = clarabelinfo.status;
      sdpisolver->opttime = clarabelinfo.solve_time;

      SCIPdebugMessage("Solved problem using Clarabel, status %d.\n", sdpisolver->solstat);

      sdpisolver->solved = TRUE;

      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelquadraticcolptr);
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelcones);
   }

   /* free memory */
   BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &clarabelquadraticcolptr);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelcolptr);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelvals);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelcols);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelrows);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelobj);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelrhs);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelblockstart);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelblocksizes);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &clarabelbounds);

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

/** returns true if the solver could determine whether the problem is feasible
 *
 *  So it returns true if the solver knows that the problem is feasible/infeasible/unbounded, it returns false if the
 *  solver does not know anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be
 *  used.
 */
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelUnsolved:
   case ClarabelAlmostSolved:
   case ClarabelAlmostPrimalInfeasible:
   case ClarabelAlmostDualInfeasible:
   case ClarabelMaxIterations:
   case ClarabelMaxTime:
   case ClarabelNumericalError:
   case ClarabelInsufficientProgress:
   case ClarabelCallbackTerminated:
      return FALSE;

   case ClarabelDualInfeasible:
   {
      ClarabelDefaultSolution clarabelsol;
      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_prim < sdpisolver->epsilon )
         return TRUE;
      return FALSE;
   }
   case ClarabelSolved:
   case ClarabelPrimalInfeasible:
      return TRUE;
   default:
      SCIPdebugMessage("Unknown return code in SCIPsdpiSolverFeasibilityKnown\n");
      return FALSE;
   }
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   assert( sdpisolver != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelSolved:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case ClarabelPrimalInfeasible:
   {
      ClarabelDefaultSolution clarabelsol;

      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_dual < sdpisolver->epsilon )
         *primalfeasible = TRUE;
      else
         *primalfeasible = FALSE;
      *dualfeasible = FALSE;
      break;
   }
   case ClarabelDualInfeasible:
   {
      ClarabelDefaultSolution clarabelsol;

      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_prim < sdpisolver->epsilon )
         *dualfeasible = TRUE;
      else
         *dualfeasible = FALSE;
      *primalfeasible = FALSE;
      break;
   }
   default:
      SCIPdebugMessage("Clarabel does not know about feasibility of solutions.\n");
      return SCIP_LPERROR;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   if ( sdpisolver->solstat == ClarabelPrimalInfeasible )
   {
      ClarabelDefaultSolution clarabelsol;
      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_dual < sdpisolver->epsilon )
         return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelDualInfeasible:
      return TRUE;
   case ClarabelSolved:
      break;
   default:
      SCIPdebugMessage("Clarabel does not know about feasibility of solutions.\n");
      break;
   }
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelSolved:
      return TRUE;
   case ClarabelDualInfeasible:
      break;
   case ClarabelPrimalInfeasible:
   {
      ClarabelDefaultSolution clarabelsol;
      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_dual < sdpisolver->epsilon )
         return TRUE;
   }
   default:
      SCIPdebugMessage("Clarabel does not know about feasibility of solutions.\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   if ( sdpisolver->solstat == ClarabelDualInfeasible )
   {
      ClarabelDefaultSolution clarabelsol;
      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_prim < sdpisolver->epsilon )
         return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelPrimalInfeasible:
      return TRUE;
   case ClarabelSolved:
      break;
   default:
      SCIPdebugMessage("Clarabel does not know about feasibility of solutions.\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelSolved:
      return TRUE;
   case ClarabelPrimalInfeasible:
      break;
   case ClarabelDualInfeasible:
   {
      ClarabelDefaultSolution clarabelsol;
      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      if ( clarabelsol.r_prim < sdpisolver->epsilon )
         return TRUE;
      break;
   }
   default:
      SCIPdebugMessage("Clarabel does not know about feasibility of solutions.\n");
      break;
   }/*lint !e788*/
   return FALSE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   CHECK_IF_SOLVED_BOOL( sdpisolver );

   switch ( sdpisolver->solstat )
   {
   case ClarabelSolved:
   case ClarabelPrimalInfeasible:
   case ClarabelDualInfeasible:
      return TRUE;

   case ClarabelAlmostSolved:
   case ClarabelAlmostPrimalInfeasible:
   case ClarabelAlmostDualInfeasible:
   case ClarabelMaxIterations:
   case ClarabelMaxTime:
   case ClarabelNumericalError:
   case ClarabelInsufficientProgress:
   case ClarabelCallbackTerminated:
      return FALSE;

   default:
      SCIPdebugMessage("Clarabel does not know about convergence.\n");
   }
   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   /* Clarabel does not have an objective limit */
   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   if ( sdpisolver->solstat == ClarabelMaxIterations )
      return TRUE;
   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( sdpisolver->solstat == ClarabelMaxTime )
      return TRUE;
   return FALSE;
}

/** returns the internal solution status of the solver, which has the following meaning:<br>
 * -1: solver was not started<br>
 *  0: converged<br>
 *  1: infeasible start<br>
 *  2: numerical problems<br>
 *  3: objective limit reached<br>
 *  4: iteration limit reached<br>
 *  5: time limit reached<br>
 *  6: user termination<br>
 *  7: other
 */
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( ! sdpisolver->solved )
      return -1;

   if ( sdpisolver->timelimit )
      return 5;

   switch ( sdpisolver->solstat )
   {
   case ClarabelUnsolved:
      return -1;

   case ClarabelSolved:
   case ClarabelPrimalInfeasible:
   case ClarabelDualInfeasible:
      return 0;

   case ClarabelAlmostSolved:
   case ClarabelAlmostPrimalInfeasible:
   case ClarabelAlmostDualInfeasible:
   case ClarabelNumericalError:
   case ClarabelInsufficientProgress:
      return 2;

   case ClarabelMaxIterations:
      return 4;

   case ClarabelMaxTime:
      return 5;

   case ClarabelCallbackTerminated:
   default:
      return 7;
   }
   return -1;
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   CHECK_IF_SOLVED_BOOL( sdpisolver );

   if ( sdpisolver->solstat != ClarabelSolved )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework
 */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   return SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverFeasibilityKnown(sdpisolver);
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{
   ClarabelDefaultSolution clarabelsol;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( objval != NULL );

   /* check for unboundedness */
   if ( SCIPsdpiSolverIsDualUnbounded(sdpisolver) || SCIPsdpiSolverIsPrimalInfeasible(sdpisolver) )
   {
      *objval = -SCIPsdpiSolverInfinity(sdpisolver);
      return SCIP_OKAY;
   }

   clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
   assert( clarabelsol.x_length == sdpisolver->nactivevars + sdpisolver->penalty ? 1 : 0 );

   if ( sdpisolver->penalty && ! sdpisolver->feasorig )
   {
      /* in this case we cannot really trust the solution given by Clarabel, since changes in the value of r much less than epsilon can
       * cause huge changes in the objective, so using the objective value given by Clarabel is numerically more stable */
      *objval *= clarabelsol.obj_val * sdpisolver->objscalefactor;
   }
   else
   {
      int v;

      /* since the objective value given by Clarabel sometimes differs slightly from the correct value for the given solution,
       * we get the solution from Clarabel and compute the correct objective value */
      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
         *objval += clarabelsol.x[v] * sdpisolver->objcoefs[v];
   }

   /* as we didn't add the fixed (lb = ub) variables to Clarabel, we have to add their contributions to the objective as well */
   *objval += sdpisolver->fixedvarsobjcontr;

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs */
SCIP_RETCODE SCIPsdpiSolverGetDualSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval,             /**< pointer to store the objective value (or NULL) */
   SCIP_Real*            dualsol             /**< array of length nvars to store the dual solution vector (or NULL) */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( dualsol != NULL )
   {
      ClarabelDefaultSolution clarabelsol;
      int v;

      clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
      assert( clarabelsol.x_length == sdpisolver->nactivevars + sdpisolver->penalty ? 1 : 0 );

      /* insert the entries into dualsol, for non-fixed vars we copy those from Clarabel, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtocblmap[v] >= 0 )
            dualsol[v] = clarabelsol.x[sdpisolver->inputtocblmap[v]];  /* use dual solution */
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            assert( -sdpisolver->inputtocblmap[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtocblmap[v]) - 1]; /*lint !e679*/ /* -1 because we wanted strictly negative vals */
         }
      }

      /* if both solution and objective should be returned, we can use the solution to compute the objective */
      if ( objval != NULL )
      {
         if ( sdpisolver->penalty && ! sdpisolver->feasorig )
         {
            /* reverse scaling */
            *objval = clarabelsol.obj_val * sdpisolver->objscalefactor;
         }
         else
         {
            /* since the objective value given by Clarabel sometimes differs slightly from the correct value for the given solution,
             * we get the solution from Clarabel and compute the correct objective value */
            *objval = 0.0;
            for (v = 0; v < sdpisolver->nactivevars; v++)
               *objval += clarabelsol.x[v] * sdpisolver->objcoefs[v];
         }

         /* as we didn't add the fixed (lb = ub) variables to Clarabel, we have to add their contributions to the objective as well */
         *objval += sdpisolver->fixedvarsobjcontr;
      }
   }
   else if ( objval != NULL )
   {
      SCIP_CALL( SCIPsdpiSolverGetObjval(sdpisolver, objval) );
   }

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X for the preoptimal solution */
SCIP_RETCODE SCIPsdpiSolverGetPreoptimalPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz */
   int*                  startXnblocknonz    /**< array to return number of nonzeros for row/col/val-arrays in each block
                                              *   or first entry equal to -1 if no primal solution is available */
   )
{ /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet.\n");

   return SCIP_NOTIMPLEMENTED;
}

/** gets preoptimal dual solution vector and primal matrix for warmstarting purposes
 *
 *  @note The last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity.
 */
SCIP_RETCODE SCIPsdpiSolverGetPreoptimalSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success,            /**< Could a preoptimal solution be returned? */
   SCIP_Real*            dualsol,            /**< array to return the dual solution vector (or NULL) */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks (+ 1)) or -1 if no primal matrix should be returned */
   int*                  startXnblocknonz,   /**< input: size of row/col/val-arrays in each block (or NULL if nblocks = -1)
                                              *   output: number of nonzeros in each block or first entry -1 if no primal solution is available */
   int**                 startXrow,          /**< array for returning row indices of X (or NULL if nblocks = -1) */
   int**                 startXcol,          /**< array for returning column indices of X (or NULL if nblocks = -1) */
   SCIP_Real**           startXval           /**< array for returning values of X (or NULL if nblocks = -1) */
   )
{/*lint !e1784*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_NOTIMPLEMENTED;
}/*lint !e715*/

/** gets the solution corresponding to the lower and upper variable-bounds in the primal problem
 *
 *  The arrays need to have size nvars.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_RETCODE SCIPsdpiSolverGetPrimalBoundVars(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            lbvals,             /**< array to store the values of the variables corresponding to lower bounds in the primal problems */
   SCIP_Real*            ubvals              /**< array to store the values of the variables corresponding to upper bounds in the primal problems */
   )
{
   ClarabelDefaultSolution clarabelsol;
   int ind;
   int i;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( lbvals != NULL );
   assert( ubvals != NULL );

   /* initialize the return-arrays with zero */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      lbvals[i] = 0.0;
      ubvals[i] = 0.0;
   }

   if ( sdpisolver->nvarbounds <= 0 )
      return SCIP_OKAY;

   /* get solution from Clarabel */
   clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
   assert( clarabelsol.z_length == sdpisolver->nconss );

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars );

   ind = sdpisolver->nsdpconss + sdpisolver->nlpconss;
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      /* this is a lower bound */
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* the last nvarbounds entries correspond to the varbounds; we need to unscale these values */
         lbvals[sdpisolver->cbltoinputmap[- sdpisolver->varboundpos[i] -1]] = clarabelsol.z[ind];
         ++ind;
      }
      else
      {  /* this is an upper bound */
         assert( sdpisolver->varboundpos[i] > 0 );

         /* the last nvarbounds entries correspond to the varbounds; we need to unscale these values */
         ubvals[sdpisolver->cbltoinputmap[sdpisolver->varboundpos[i] - 1]] = clarabelsol.z[ind];
         ++ind;
      }
   }
   assert( ind <= sdpisolver->nconss );

   return SCIP_OKAY;
}

/** gets the primal solution corresponding to the LP row sides */
SCIP_RETCODE SCIPsdpiSolverGetPrimalLPSides(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int                   nlpcons,            /**< number of LP rows */
   int*                  lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   SCIP_Real*            lplhs,              /**< lhs of LP rows */
   SCIP_Real*            lprhs,              /**< rhs of LP rows */
   SCIP_Real*            lhsvals,            /**< array to store the values of the variables corresponding to LP lhs */
   SCIP_Real*            rhsvals             /**< array to store the values of the variables corresponding to LP rhs */
   )
{
   ClarabelDefaultSolution clarabelsol;
   int ind;
   int i;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( lplhs != NULL );
   assert( lprhs != NULL );
   assert( lhsvals != NULL );
   assert( rhsvals != NULL );

   if ( nlpcons <= 0 )
      return SCIP_OKAY;

   /* get solution from Clarabel */
   clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
   assert( clarabelsol.z_length == sdpisolver->nconss );

   /* loop through LP rows */
   ind = sdpisolver->nsdpconss;
   for (i = 0; i < nlpcons; i++)
   {
      if ( lpindchanges[i] < 0 )
      {
         lhsvals[i] = 0.0;
         rhsvals[i] = 0.0;
         continue;
      }

      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
      {
         lhsvals[i] = clarabelsol.z[ind];
         ++ind;
      }
      else
         lhsvals[i] = 0.0;

      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         rhsvals[i] = clarabelsol.z[ind];
         ++ind;
      }
      else
         rhsvals[i] = 0.0;
   }
   assert( ind <= sdpisolver->nsdpconss + sdpisolver->nlpconss );

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X (including lp block) */
SCIP_RETCODE SCIPsdpiSolverGetPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet.\n");
   return SCIP_NOTIMPLEMENTED;
}

/** returns the primal matrix X
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */
SCIP_RETCODE SCIPsdpiSolverGetPrimalMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block
                                              *   output: number of nonzeros in each block */
   int**                 startXrow,          /**< pointer to store row indices of X */
   int**                 startXcol,          /**< pointer to store column indices of X */
   SCIP_Real**           startXval           /**< pointer to store values of X */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet.\n");
   return SCIP_NOTIMPLEMENTED;
}

/** returns the primal solution matrix (without LP rows) */
SCIP_RETCODE SCIPsdpiSolverGetPrimalSolutionMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nsdpblocks,         /**< number of blocks */
   int*                  sdpblocksizes,      /**< sizes of the blocks */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< pointer to store index change for each block, system is the same as for indchanges */
   SCIP_Real**           primalmatrices      /**< pointer to store values of the primal matrices */
   )
{
   ClarabelDefaultSolution clarabelsol;
   SCIP_Real sqrt2;
   int blockstart = 0;
   int b;

   assert( sdpisolver != NULL );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( indchanges != NULL );
   assert( nremovedinds != NULL );
   assert( blockindchanges != NULL );
   assert( primalmatrices != NULL );

   if ( nsdpblocks <= 0 )
      return SCIP_OKAY;

   /* get solution from Clarabel */
   clarabelsol = clarabel_DefaultSolver_solution(sdpisolver->clarabelsolver);
   assert( clarabelsol.z_length == sdpisolver->nconss );

   sqrt2 = sqrt(2);

   /* loop over all SDP blocks */
   for (b = 0; b < nsdpblocks; b++)
   {
      int blocksize;
      int j;

      assert( primalmatrices[b] != NULL );

      blocksize = sdpblocksizes[b];

      /* initialize solution matrix with 0s */
      for (j = 0; j < blocksize * blocksize; ++j)
         primalmatrices[b][j] = 0.0;

      /* treat blocks that were not removed */
      if ( blockindchanges[b] >= 0 )
      {
         SCIP_Real val;
         int redsize;
         int row;
         int col;
         int ind;
         int i;

         redsize = blocksize - nremovedinds[b];

         /* fill in matrix */
         for (i = 0; i < blocksize; ++i)
         {
            if ( indchanges[b][i] >= 0 )
            {
               row = i - indchanges[b][i];
               assert( 0 <= row && row < redsize );

               for (j = 0; j <= i; ++j)
               {
                  if ( indchanges[b][j] >= 0 )
                  {
                     col = j - indchanges[b][j];
                     assert( 0 <= col && col < redsize );

                     ind = blockstart + row * (row + 1) / 2 + col;
                     assert( ind < sdpisolver->nsdpconss );
                     val = clarabelsol.z[ind];

                     if ( row != col )
                        val /= sqrt2;

                     primalmatrices[b][i * blocksize + j] = val;
                     primalmatrices[b][j * blocksize + i] = val;
                  }
               }
            }
         }
         blockstart += redsize * (redsize + 1) / 2;
      }
   }
   return SCIP_OKAY;
}

/** return the maximum absolute value of the optimal primal matrix */
SCIP_Real SCIPsdpiSolverGetMaxPrimalEntry(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet.\n");
   return SCIP_INVALID;
}

/** gets the time for the last SDP optimization call of solver */
SCIP_RETCODE SCIPsdpiSolverGetTime(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real*            opttime             /**< pointer to store the time for optimization of the solver */
   )
{
   assert( sdpisolver != NULL );
   assert( opttime != NULL );

   *opttime = sdpisolver->opttime;

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

   *iterations = sdpisolver->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{/*lint --e{715,1784}*/
   assert( sdpisolver != NULL );
   assert( calls != NULL );

   *calls = sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP interface solver structure */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP solver */
   )
{
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );

   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) )
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   else if ( sdpisolver->penalty )
      *usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;
   else
      *usedsetting = SCIP_SDPSOLVERSETTING_FAST;

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{/*lint --e{715}*/
   return 1.0e20;
}

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
}

/** gets floating point parameter of SDP-Solver */
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
   case SCIP_SDPPAR_GAPTOL:
      *dval = sdpisolver->gaptol;
      break;
   case SCIP_SDPPAR_FEASTOL:
      *dval = sdpisolver->feastol;
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      *dval = sdpisolver->sdpsolverfeastol;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP-Solver */
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
      SCIPdebugMessage("Setting sdpisolver epsilon to %g.\n", dval);
      break;
   case SCIP_SDPPAR_GAPTOL:
      sdpisolver->gaptol = dval;
      SCIPdebugMessage("Setting sdpisolver gaptol to %g.\n", dval);
      break;
   case SCIP_SDPPAR_FEASTOL:
      sdpisolver->feastol = dval;
      SCIPdebugMessage("Setting sdpisolver feastol to %g.\n", dval);
      break;
   case SCIP_SDPPAR_SDPSOLVERFEASTOL:
      sdpisolver->sdpsolverfeastol = dval;
      SCIPdebugMessage("Setting sdpisolver sdpsolverfeastol to %g.\n", dval);
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %g.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** gets integer parameter of SDP-Solver */
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
      *ival = (int) sdpisolver->sdpinfo;
      SCIPdebugMessage("Getting sdpisolver information output (%d).\n", *ival);
      break;
   case SCIP_SDPPAR_NTHREADS:
      *ival = sdpisolver->nthreads;
      SCIPdebugMessage("Getting sdpisolver number of threads: %d.\n", *ival);
      break;
   case SCIP_SDPPAR_USEPRESOLVING:
      *ival = (int) sdpisolver->usepresolving;
      SCIPdebugMessage("Getting usepresolving (%d).\n", *ival);
      break;
   case SCIP_SDPPAR_USESCALING:
      *ival = (int) sdpisolver->usescaling;
      SCIPdebugMessage("Getting usescaling (%d).\n", *ival);
      break;
   case SCIP_SDPPAR_SCALEOBJ:
      *ival = (int) sdpisolver->scaleobj;
      SCIPdebugMessage("Getting scaleobj (%d).\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpisolver != NULL );

   switch( type )
   {
   case SCIP_SDPPAR_NTHREADS:
      sdpisolver->nthreads = ival;
      SCIPdebugMessage("Setting sdpisolver number of threads to %d.\n", ival);
      break;
   case SCIP_SDPPAR_SDPINFO:
      assert( 0 <= ival && ival <= 1 );
      sdpisolver->sdpinfo = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting sdpisolver information output (%d).\n", ival);
      break;
   case SCIP_SDPPAR_USEPRESOLVING:
      assert( 0 <= ival && ival <= 1 );
      sdpisolver->usepresolving = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting usepresolving (%d).\n", ival);
      break;
   case SCIP_SDPPAR_USESCALING:
      assert( 0 <= ival && ival <= 1 );
      sdpisolver->usescaling = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting usescaling (%d).\n", ival);
      break;
   case SCIP_SDPPAR_SCALEOBJ:
      assert( 0 <= ival && ival <= 1 );
      sdpisolver->scaleobj = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting scaleobj (%d).\n", ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }/*lint !e788*/

   return SCIP_OKAY;
}

/** compute and set lambdastar (only used for SDPA) */
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Lambdastar parameter not used by Clarabel"); /* this parameter is only used by SDPA */

   return SCIP_OKAY;
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{/*lint --e{1784}*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( penaltyparam != NULL );

   compval = PENALTYPARAM_FACTOR * maxcoeff;

   if ( compval < MIN_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %g.\n", MIN_PENALTYPARAM);
      *penaltyparam = MIN_PENALTYPARAM;
   }
   else if ( compval > MAX_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %g.\n", MAX_PENALTYPARAM);
      *penaltyparam = MAX_PENALTYPARAM;
   }
   else
   {
      SCIPdebugMessage("Setting penaltyparameter to %g.\n", compval);
      *penaltyparam = compval;
   }

   return SCIP_OKAY;
}

/** compute and set the maximum penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputeMaxPenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   )
{/*lint --e{1784}*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( maxpenaltyparam != NULL );

   compval = penaltyparam * MAXPENALTYPARAM_FACTOR;

   if ( compval < MAX_MAXPENALTYPARAM )
   {
      *maxpenaltyparam = compval;
      SCIPdebugMessage("Setting maximum penaltyparameter to %g.\n", compval);
   }
   else
   {
      *maxpenaltyparam = MAX_MAXPENALTYPARAM;
      SCIPdebugMessage("Setting penaltyparameter to %g.\n", MAX_MAXPENALTYPARAM);
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
   assert( sdpisolver != NULL );
   assert( fname != NULL );

#ifdef FEATURE_SERDE
   if ( sdpisolver->clarabelsolver != NULL )
   {
      clarabel_DefaultSolver_free(sdpisolver->clarabelsolver);
      sdpisolver->clarabelsolver = NULL;
   }
   sdpisolver->clarabelsolver = clarabel_DefaultSolver_load_from_file(fname);
   return SCIP_OKAY;
#else
   SCIPerrorMessage("Clarabel not compiled with FEATURE_SERDE.\n");
   return SCIP_LPERROR;
#endif
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   assert( sdpisolver != NULL );
   assert( fname != NULL );

#if FEATURE_SERDE
   clarabel_DefaultSolver_save_to_file(sdpisolver->clarabelsolver, fname);
   return SCIP_OKAY;
#else
   SCIPerrorMessage("Clarabel not compiled with FEATURE_SERDE.\n");
   return SCIP_LPERROR;
#endif
}

/**@} */
