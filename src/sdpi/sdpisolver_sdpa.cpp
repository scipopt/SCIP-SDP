/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define SCIP_DEBUG */
/* #define SCIP_MORE_DEBUG */
/* #define SCIP_DEBUG_PRINTTOFILE */ /* prints each problem inserted into SDPA to the file sdpa.dat-s and the starting point to sdpa.ini-s */

/**@file   sdpisolver_sdpa.cpp
 * @brief  interface for SDPA
 * @author Tristan Gally
 * @author Ambros Gleixner
 * @author Marc Pfetsch
 *
 *
 * SDPA solve the dual problem as SCIP-SDP, but with a slightly changed notation:
 *   \f{eqnarray*}{
 *      \min & & c^T y \\
 *      \mbox{s.t.} & & X^{(k)} = \sum_{i \in I} A_i^{(k)} x_i - A_0^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & \sum_{i \in I} d_{ij} x_i \geq b_j \quad \forall \ j \in J, \\
 *      & & \ell_i \leq y_i \leq u_i \quad \forall \ i \in I
 *   \f}
 * The primal problem is then
 *   \f{eqnarray*}{
 *      \max & & \sum_{k \in K} A_0^{(k)} \bullet Y^{(k)} + \sum_{j \in J} b_j y_j - \sum_{i \in I_u} u_i v_i + \sum_{i \in I_\ell} \ell_i w_i \\
 *      \mbox{s.t.} & & \sum_{k \in K} A_i^{(k)} \bullet Y^{(k)} + \sum_{j \in J} d_{ij} y_j - 1_{\{u_i < \infty\}} v_i + 1_{\{\ell_i > -\infty\}} w_i = c_i \quad \forall \ i \in I,\\
 *      & & Y^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & y_j \geq 0 \quad \forall \ j \in J,\\
 *      & & v_i \geq 0 \quad \forall \ i \in I_u,\\
 *      & & w_i \geq 0 \quad \forall \ i \in I_\ell,
 *   \f}
 * In the code below, we use \f$ X\f$ for the SDPA primal variable \f$ Y\f$ and we use \f$ Z\f$ for the SDPA dual variable \f$ X\f$.
 */

#include <assert.h>
#ifdef OPENBLAS
#include <cblas.h>
#endif

#include "sdpi/sdpisolver.h"

/* turn off warnings for sdpa (doesn't seem to work) */
#pragma GCC diagnostic ignored "-Wcast-qual"
#pragma GCC diagnostic ignored "-Wpedantic"
#include "sdpa_call.h"                       /* SDPA callable library interface */
#pragma GCC diagnostic warning "-Wcast-qual"
#pragma GCC diagnostic warning "-Wpedantic"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */
#include "scip/pub_message.h"                /* for debug and error message */

/* turn off lint warnings for whole file: */
/*lint --e{1784}*/


/* local defines */
#define GAPTOLCHANGE                1        /**< change gaptol by this factor when switching from fast to default and from default to stable settings */
#define FEASTOLCHANGE               1        /**< change feastol by this factor when switching from fast to default and from default to stable settings */
#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */

#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1E-9     /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance */

#define MIN_LAMBDASTAR              1e-6     /**< if lambda star is to be computed, this is the minimum value it will take */
#define MAX_LAMBDASTAR              1e8      /**< if lambda star is to be computed, this is the maximum value it will take */
#define LAMBDASTAR_FACTOR           1e0      /**< if lambda star is to be computed, the biggest guess of the SDP blocks is multiplied by this value */
#define LAMBDASTAR_TWOPOINTS        TRUE     /**< if lambda star is to be computed, should we use only a low and a high value or instead a continuous interval */
#define LAMBDASTAR_THRESHOLD        1e1      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then we distinguish between low and high using this */
#define LAMBDASTAR_LOW              1.5      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then this is the value for below the threshold */
#define LAMBDASTAR_HIGH             1e5      /**< if lambda star is to be computed and LAMBDASTAR_TWOPOINTS=TRUE, then this is the value for above the threshold */
#define LAMBDASTAR_DEFAULT          1e2      /**< default value of lambda star taken from SDPA fast settings */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e12     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e1      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */


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
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                 \
                      {                                                                                      \
                         if (!(sdpisolver->solved))                                                          \
                         {                                                                                   \
                            SCIPerrorMessage("Tried to access solution information for SDP %d ahead of solving!\n", sdpisolver->sdpcounter);  \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )


/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   SDPA*                 sdpa;               /**< solver-object */
   int                   nvars;              /**< number of input variables */
   int                   maxnvars;           /**< size of the arrays inputtomosekmapper, mosektoinputmapper, fixedvarsval, and objcoefs */
   int                   nactivevars;        /**< number of variables present in SDPA (nvars minus the number of variables with lb = ub) */
   int*                  inputtosdpamapper;  /**< entry i gives the index of input variable i in sdpa (starting from 1) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  sdpatoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in sdpa (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   int                   nvarbounds;         /**< number of variable bounds given to sdpa, length of sdpavarboundpos */
   int*                  varboundpos;        /**< maps position of variable bounds in the variable bound part of the LP-block in sdpa to the sdpa-indices
                                              *   of the corresponding variables, -n means lower bound of variable n, +n means upper bound */
   int*                  inputtovbmapper;    /**< maps lower and upper bounds of input variables to positions in varboundarray: entry 2i gives
                                              *   position of lower bound, entry 2i+1 gives position of upper bound */
   int                   nlpcons;            /**< number of LP constraints */
   int                   maxnlpcons;         /**< maximal number of LP constraints */
   int                   nlpineqs;           /**< number of inequalities in SDPA */
   int                   nsdpblocks;         /**< number of given sdp blocks; length of inputtoblockmapper */
   int                   maxnsdpblocks;      /**< maximal number of given sdp blocks */
   int*                  maxsdpblocksizes;   /**< array of maximal block sizes */
   int*                  rowmapper;          /**< entry 2i gives SDPA-index of left side of ranged row i, 2i+1 index of rhs */
   int*                  rowtoinputmapper;   /**< if rowtoinputmapper[i] = 2j, then i-th SDPA row corresponds to lhs of ranged row j, rhs if 2j + 1*/
   int*                  inputtoblockmapper; /**< entry i gives sdpa-index of i-th sdp-block (or -1 if removed); length is nsdpblocks */
   int**                 blockindmapper;     /**< entry b,i gives original index of row/col i in (sdpa-)block b; length is sdpa-nsdpblocks (LP block not included) times sdpa-blocksize[b] */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   SCIP_Bool             timelimit;          /**< Was the SDP not given to the solver because the time limit was already reached? */
   int                   sdpcounter;         /**< used for debug messages */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
   SCIP_Real             epsilon;            /**< tolerance for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance given to SDP-solver */
   SCIP_Real             objlimit;           /**< objective limit for SDP-solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP-solver output information to the screen? */
   SCIP_Bool             penalty;            /**< was the problem last solved using a penalty formulation */
   SCIP_Bool             feasorig;           /**< was the last problem solved with a penalty formulation and with original objective coefficents
                                              *   and the solution was feasible for the original problem? */
   SCIP_Bool             rbound;             /**< was the penalty parameter bounded during the last solve call */
   SCIP_SDPSOLVERSETTING usedsetting;        /**< setting used to solve the last SDP */
   SCIP_Real             lambdastar;         /**< lambda star parameter to give to SDPA for initial point */
   SCIP_Real*            preoptimalsol;      /**< first feasible solution with gap less or equal preoptimalgap */
   SCIP_Real**           preoptimalsolx;     /**< dense blockwise primal matrix for first feasible solution with gap less or equal preoptimalgap */
   SCIP_Real*            preoptimalsolxlp;   /**< LP part of primal solution for first feasible solution with gap less or equal preoptimalgap */
   SCIP_Bool             preoptimalsolexists;/**< saved feasible solution with gap less or equal preoptimalgap */
   SCIP_Real             preoptimalgap;      /**< gap at which a preoptimal solution should be saved for warmstarting purposes */
   int                   nthreads;           /**< number of threads the SDP solver should use (-1 = number of cores) */
};


/*
 * Local Functions
 */

#ifndef NDEBUG
/** Test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( sdpisolver != NULL );
   assert( lb < ub + sdpisolver->feastol );

   if ( ub - lb <= sdpisolver->epsilon )
      return TRUE;
   return FALSE;
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
   int                   nvars,              /**< number of variables */
   int                   nsdpblocks,         /**< number of SDP blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP constraints */
   SCIP_Bool             usepreoptimalsol    /**< whether preoptimalsol(x) should be set up */
   )
{
   int newsize;
   int oldmaxnsdpblocks;
   int oldmaxnlpcons;
   int oldmaxnvars;
   int blocksize;
   int b;

   assert( sdpisolver != NULL );

   oldmaxnlpcons = sdpisolver->maxnlpcons;
   oldmaxnvars = sdpisolver->maxnvars;

   /* treat arrays for variables */
   if ( nvars > sdpisolver->maxnvars )
   {
      newsize = calcGrowSize(sdpisolver->maxnvars, nvars);

      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->inputtosdpamapper, sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->sdpatoinputmapper, sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->fixedvarsval, sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->objcoefs, sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->inputtovbmapper, 2 * sdpisolver->maxnvars, 2 * newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->varboundpos, 2 * sdpisolver->maxnvars, 2 * newsize) );

      if ( usepreoptimalsol )
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->preoptimalsol, sdpisolver->maxnvars, newsize) );
      }

      sdpisolver->maxnvars = newsize;
   }

   /* treat arrays for SDP constraints */
   if ( nsdpblocks > sdpisolver->maxnsdpblocks )
   {
      oldmaxnsdpblocks = sdpisolver->maxnsdpblocks;
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->inputtoblockmapper, sdpisolver->maxnsdpblocks, nsdpblocks) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->blockindmapper, sdpisolver->maxnsdpblocks, nsdpblocks) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->maxsdpblocksizes, sdpisolver->maxnsdpblocks, nsdpblocks) );
      if ( usepreoptimalsol )
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->preoptimalsolx, sdpisolver->maxnsdpblocks, nsdpblocks) );
      }
      sdpisolver->maxnsdpblocks = nsdpblocks;

      for (b = 0; b < oldmaxnsdpblocks; ++b)
      {
         blocksize = sdpblocksizes[b];
         if ( blocksize > sdpisolver->maxsdpblocksizes[b] )
         {
            BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->blockindmapper[b], sdpisolver->maxsdpblocksizes[b], blocksize);

            if ( usepreoptimalsol )
            {
               BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->preoptimalsolx[b], sdpisolver->maxsdpblocksizes[b] * sdpisolver->maxsdpblocksizes[b], blocksize * blocksize);
            }
            sdpisolver->maxsdpblocksizes[b] = blocksize;
         }
      }

      for (b = oldmaxnsdpblocks; b < nsdpblocks; ++b)
      {
         blocksize = sdpblocksizes[b];
         BMSallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->blockindmapper[b], blocksize);
         if ( usepreoptimalsol )
         {
            BMSallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->preoptimalsolx[b], blocksize * blocksize);
         }
         sdpisolver->maxsdpblocksizes[b] = blocksize;
      }
   }

   /* treat arrays for LP constraints */
   if ( nlpcons > sdpisolver->maxnlpcons )
   {
      newsize = calcGrowSize(sdpisolver->maxnlpcons, nlpcons);
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->rowmapper, 2 * sdpisolver->maxnlpcons, 2 * newsize) ); /*lint !e647*/
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->rowtoinputmapper, 2 * sdpisolver->maxnlpcons, 2 * newsize) );
      sdpisolver->maxnlpcons = newsize;
   }

   /* special case: treat preoptimalsolxlp */
   if ( usepreoptimalsol && ( nlpcons > oldmaxnlpcons || nvars > oldmaxnvars ) )
   {
      assert( nlpcons <= sdpisolver->maxnlpcons );
      assert( nvars <= sdpisolver->maxnvars );

      BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->preoptimalsolxlp, 2 * (oldmaxnlpcons + oldmaxnvars), 2 * (sdpisolver->maxnlpcons + sdpisolver->maxnvars));
   }

   return SCIP_OKAY;
}


/** If the problem is feasible for SDPA but not within our feasibility tolerance, adjust feasibility tolerance in
 *  SDPA and resolve until feasibility in SDPA and feasibility with regards to our tolerance match
 */
static
SCIP_RETCODE checkFeastolAndResolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< penalty parameter Gamma */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indices will be modified by these, see indchanges */
   int                   nlpcons,            /**< number LP-constraints */
   int*                  lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lpbeg,              /**< start index of each row in ind- and val-array */
   int*                  lpind,              /**< column-index for each entry in lpval-array */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries */
   SCIP_Real*            feastol,            /**< current feasibility tolerance for SDPA */
   SCIP_Real*            gaptol              /**< current gap tolerance for SDPA */
   )
{
#ifdef SCIP_DEBUG
   char phase_string[15];
#endif

   assert( feastol != NULL );
   assert( gaptol != NULL );

   if ( penaltyparam >= sdpisolver->epsilon )
      return SCIP_OKAY;

   while ( SCIPsdpiSolverIsAcceptable(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver) && *feastol >= INFEASMINFEASTOL && *gaptol >= INFEASMINFEASTOL )
   {
      SCIP_Real* solvector;
      SCIP_Real primalobj;
      SCIP_Real dualobj;
      SCIP_Bool infeasible;
      SCIP_Bool solveagain = FALSE;

      /* get current solution */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solvector, nvars) ); /*lint !e530*/
      SCIP_CALL( SCIPsdpiSolverGetDualSol(sdpisolver, NULL, solvector) );

      /* check the solution for feasibility with regards to our tolerance */
      SCIP_CALL( SCIPsdpSolcheckerCheck(sdpisolver->bufmem, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
            indchanges, nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval,
            solvector, sdpisolver->feastol, sdpisolver->epsilon, &infeasible) );
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solvector);

      if ( infeasible )
      {
         *feastol *= INFEASFEASTOLCHANGE;
         if ( *feastol >= INFEASMINFEASTOL )
         {
            SCIPdebugMessage("Solution feasible for SDPA but outside feasibility tolerance, changing SDPA feasibility tolerance to %g.\n", *feastol);
            sdpisolver->sdpa->setParameterEpsilonDash(*feastol);
            solveagain = TRUE;
         }
      }

      /* check whether duality gap is small enough */
      dualobj = sdpisolver->sdpa->getPrimalObj();
      primalobj = sdpisolver->sdpa->getDualObj();
      if ( REALABS(dualobj - primalobj) >= sdpisolver->gaptol )
      {
         infeasible = TRUE;
         *gaptol *= INFEASFEASTOLCHANGE;
         if ( *gaptol >= INFEASMINFEASTOL )
         {
            SCIPdebugMessage("Solution feasible, but duality gap is too large, changing SDPA gap tolerance to %g.\n", *gaptol);
            sdpisolver->sdpa->setParameterEpsilonStar(*gaptol);
            solveagain = TRUE;
         }
      }

      if ( solveagain )
      {
#ifdef SCIP_MORE_DEBUG
         sdpisolver->sdpa->printParameters(stdout);
#endif
         sdpisolver->sdpa->setInitPoint(false);
         sdpisolver->sdpa->resetParameters();
         sdpisolver->sdpa->solve();

         /* update number of SDP-iterations and -calls */
         sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
         ++sdpisolver->nsdpcalls;

#ifdef SCIP_DEBUG
         /* print the phase value , i.e. whether solving was successfull */
         sdpisolver->sdpa->getPhaseString((char*)phase_string);
         SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif
      }
      else
      {
         if ( infeasible )
         {
            sdpisolver->solved = FALSE;
            SCIPmessagePrintInfo(sdpisolver->messagehdlr, "SDPA failed to reach required feasibility tolerance!\n");
         }
         break;
      }
   }

   return SCIP_OKAY;
}

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name and version (if available) of SDP-solver */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{/*lint !e1784*/
   return "SDPA"; /* version number can only be printed to file/stdout */
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{/*lint !e1784*/
   return "Primal-dual Interior Point Solver for SDPs developed by K. Fujisawa et al. (sdpa.sourceforge.net)";
}

/** gets pointer to SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->sdpa;
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return 2;
}

/** Should primal solution values be saved for warmstarting purposes? */
SCIP_Bool SCIPsdpiSolverDoesWarmstartNeedPrimal(
   void
   )
{
   return TRUE;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );
   assert( bufmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

   /* this will be properly initialized then calling solve */
   (*sdpisolver)->sdpa = NULL;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->maxnvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->nsdpblocks = 0;
   (*sdpisolver)->maxnsdpblocks = 0;
   (*sdpisolver)->inputtosdpamapper = NULL;
   (*sdpisolver)->sdpatoinputmapper = NULL;
   (*sdpisolver)->inputtoblockmapper = NULL;
   (*sdpisolver)->blockindmapper = NULL;
   (*sdpisolver)->maxsdpblocksizes = NULL;
   (*sdpisolver)->maxnlpcons = 0;
   (*sdpisolver)->rowmapper = NULL;
   (*sdpisolver)->rowtoinputmapper = NULL;
   (*sdpisolver)->fixedvarsval = NULL;
   (*sdpisolver)->fixedvarsobjcontr = 0.0;
   (*sdpisolver)->objcoefs = NULL;
   (*sdpisolver)->nvarbounds = 0;
   (*sdpisolver)->varboundpos = NULL;
   (*sdpisolver)->inputtovbmapper = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->sdpcounter = 0;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->nsdpcalls = 0;
   (*sdpisolver)->nlpcons = 0;
   (*sdpisolver)->nlpineqs = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-6;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   (*sdpisolver)->lambdastar = LAMBDASTAR_DEFAULT;

   (*sdpisolver)->preoptimalsol = NULL;
   (*sdpisolver)->preoptimalsolx = NULL;
   (*sdpisolver)->preoptimalsolxlp = NULL;
   (*sdpisolver)->preoptimalsolexists = FALSE;
   (*sdpisolver)->preoptimalgap = -1.0;

   (*sdpisolver)->nthreads = -1;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP-solver interface */
   )
{/*lint !e1784*/
   int b;
   int blocksize;

   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   for (b = 0; b < (*sdpisolver)->maxnsdpblocks; b++)
   {
      blocksize = (*sdpisolver)->maxsdpblocksizes[b];
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &((*sdpisolver)->blockindmapper[b]), blocksize);

      if ( (*sdpisolver)->preoptimalsolx != NULL )
      {
         BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &((*sdpisolver)->preoptimalsolx[b]), blocksize * blocksize);
      }
   }
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->preoptimalsolxlp, 2 * ((*sdpisolver)->maxnlpcons + (*sdpisolver)->maxnvars));

   /* free LP block of preoptimal solution X array */
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->preoptimalsolx, (*sdpisolver)->maxnsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->blockindmapper, (*sdpisolver)->maxnsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->maxsdpblocksizes, (*sdpisolver)->maxnsdpblocks);

   /* free SDPA object using destructor and free memory via blockmemshell */
   if ( (*sdpisolver)->sdpa != NULL)
      delete (*sdpisolver)->sdpa;

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->preoptimalsol, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtoblockmapper, (*sdpisolver)->maxnsdpblocks);

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->rowtoinputmapper, 2 * (*sdpisolver)->maxnlpcons);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->rowmapper, 2 * (*sdpisolver)->maxnlpcons);

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtovbmapper, 2 * (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtosdpamapper, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->sdpatoinputmapper, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP-Counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{/*lint !e1784*/
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
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independant of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint-matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number LP-constraints */
   int*                  lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lpbeg,              /**< start index of each row in ind- and val-array */
   int*                  lpind,              /**< column-index for each entry in lpval-array */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries */
   SCIP_Real*            starty,             /**< NULL or dual vector y as starting point for the solver, this should have length nvars */
   int*                  startZnblocknonz,   /**< dual matrix Z = sum Ai yi as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startZrow,          /**< dual matrix Z = sum Ai yi as starting point for the solver: row indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int**                 startZcol,          /**< dual matrix Z = sum Ai yi as starting point for the solver: column indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   SCIP_Real**           startZval,          /**< dual matrix Z = sum Ai yi as starting point for the solver: values for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int*                  startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   int**                 startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_Real**           startXval,          /**< primal matrix X as starting point for the solver: values for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   SDPI_CLOCK*           usedsdpitime        /**< clock to measure how much time has been used for the current solve */
   )
{/*lint !e1784*/
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
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independant of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< should r be non-negative ? */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint-matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number LP-constraints */
   int*                  lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP-constraints after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lpbeg,              /**< start index of each row in ind- and val-array */
   int*                  lpind,              /**< column-index for each entry in lpval-array */
   SCIP_Real*            lpval,              /**< values of LP-constraint matrix entries */
   SCIP_Real*            starty,             /**< NULL or dual vector y as starting point for the solver, this should have length nvars */
   int*                  startZnblocknonz,   /**< dual matrix Z = sum Ai yi as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startZrow,          /**< dual matrix Z = sum Ai yi as starting point for the solver: row indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int**                 startZcol,          /**< dual matrix Z = sum Ai yi as starting point for the solver: column indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   SCIP_Real**           startZval,          /**< dual matrix Z = sum Ai yi as starting point for the solver: values for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int*                  startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   int**                 startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_Real**           startXval,          /**< primal matrix X as starting point for the solver: values for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP and MOSEK, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   SDPI_CLOCK*           usedsdpitime,       /**< clock to measure how much time has been used for the current solve */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                              *   (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                              *   this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
   )
{/*lint !e1784*/
   SCIP_Real feastol;
   SCIP_Real gaptol;
   SCIP_Real solvertimelimit;
   SCIP_Bool usepreoptimalsol;
   int nfixedvars;
   bool checkinput; /* should the input be checked for consistency in SDPA? */
   int nsdpablocks;
   int nlpineqs;
   int ind;
   int i;
   int k;
   int b;

#ifdef SCIP_DEBUG
   char phase_string[15];
#endif

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
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( nlpcons == 0 || lpbeg != NULL );
   assert( nlpcons == 0 || lpindchanges != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lpind != NULL );
   assert( nlpcons == 0 || lpval != NULL );
   assert( startZnblocknonz == NULL || startZrow != NULL );
   assert( startZnblocknonz == NULL || startZcol != NULL );
   assert( startZnblocknonz == NULL || startZval != NULL );
   assert( startXnblocknonz == NULL || startXrow != NULL );
   assert( startXnblocknonz == NULL || startXcol != NULL );
   assert( startXnblocknonz == NULL || startXval != NULL );

   sdpisolver->niterations = 0;
   sdpisolver->nsdpcalls = 0;
   sdpisolver->feasorig = FALSE;

   /* compute the time limit to set for the solver */
   solvertimelimit = timelimit;
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
      solvertimelimit -= SDPIclockGetTime(usedsdpitime);

   if ( solvertimelimit <= 0.0 )
   {
      sdpisolver->solved = FALSE;
      sdpisolver->timelimit = TRUE;
      return SCIP_OKAY;
   }
   else
      sdpisolver->timelimit = FALSE;

#ifndef NDEBUG
   checkinput = true;
#else
   checkinput = false;
#endif

   sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP) */
   if ( penaltyparam < sdpisolver->epsilon )
      SCIPdebugMessage("Inserting Data into SDPA for SDP (%d) \n", ++sdpisolver->sdpcounter);
   else
      SCIPdebugMessage("Inserting Data again into SDPA for SDP (%d) \n", sdpisolver->sdpcounter);

   /* free old SDPA solver */
   if ( sdpisolver->sdpa != 0 )
   {
      /* if the SDPA solver has already been created, clear the current problem instance */
      delete sdpisolver->sdpa;
   }
   sdpisolver->sdpa = new SDPA();
   assert( sdpisolver->sdpa != 0 );

   /* set number of threads (do this early, since this might affect factorization) */
   if ( sdpisolver->nthreads > 0 )
   {
#ifdef OPENBLAS
      openblas_set_num_threads(sdpisolver->nthreads);
#endif
      sdpisolver->sdpa->setNumThreads(sdpisolver->nthreads);
   }

   /* set the penalty and rbound flags accordingly */
   sdpisolver->penalty = (penaltyparam < sdpisolver->epsilon) ? FALSE : TRUE;
   sdpisolver->rbound = rbound;

   if ( sdpisolver->preoptimalgap >= 0.0 )
      usepreoptimalsol = TRUE;
   else
      usepreoptimalsol = FALSE;

   /* ensure memory for working space */
   SCIP_CALL( ensureMappingDataMemory(sdpisolver, nvars, nsdpblocks, sdpblocksizes, nlpcons, usepreoptimalsol) );

   assert( sdpisolver->nvars <= sdpisolver->maxnvars );
   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   sdpisolver->nvarbounds = 0;
   sdpisolver->preoptimalsolexists = FALSE;
   sdpisolver->nsdpblocks = nsdpblocks;
   sdpisolver->nlpcons = nlpcons;
   nsdpablocks = nsdpblocks - nremovedblocks;
   nfixedvars = 0;

   /* find fixed variables and set up mapping functions */
   sdpisolver->fixedvarsobjcontr = 0.0;
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* this is the value this fixed variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars] = lb[i]; /* if lb=ub, then this is the value the variable will have in every solution */
         nfixedvars++;
         sdpisolver->inputtosdpamapper[i] = -nfixedvars;
         SCIPdebugMessage("Fixing variable %d locally to %g in SDPA.\n", i, lb[i]);
         sdpisolver->inputtovbmapper[2 * i] = -1;
         sdpisolver->inputtovbmapper[2 * i + 1] = -1;
      }
      else
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d in SDPA.\n", i, sdpisolver->nactivevars + 1);
#endif
         sdpisolver->sdpatoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->inputtosdpamapper[i] = sdpisolver->nactivevars + 1; /* sdpa starts counting at 1 */
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];

         if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[i]) )
         {
            sdpisolver->varboundpos[sdpisolver->nvarbounds] = -(sdpisolver->nactivevars + 1); /* negative sign means lower bound, +1 because of SDPA */
            sdpisolver->inputtovbmapper[2 * i] = sdpisolver->nvarbounds++;
         }
         else
            sdpisolver->inputtovbmapper[2 * i] = -1;

         if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[i]) )
         {
            sdpisolver->varboundpos[sdpisolver->nvarbounds] = +(sdpisolver->nactivevars + 1); /* positive sign means upper bound, +1 because of SDPA */
            sdpisolver->inputtovbmapper[2 * i + 1] = sdpisolver->nvarbounds++;
         }
         else
            sdpisolver->inputtovbmapper[2 * i + 1] = -1;

         sdpisolver->nactivevars++;
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nactivevars );

   /* if we use a penalty formulation, we need the constraint r >= 0 */
   if ( penaltyparam >= sdpisolver->epsilon && rbound )
      sdpisolver->nvarbounds++;

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* determine total number of sides in LP-constraints */
   nlpineqs = 0;
   for (i = 0; i < nlpcons; i++)
   {
      if ( lpindchanges[i] < 0 )
      {
         sdpisolver->rowmapper[2 * i] = -1;
         sdpisolver->rowmapper[2 * i + 1] = -1;
         continue;
      }

      assert( lplhs[i] < SCIPsdpiSolverInfinity(sdpisolver) );
      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
      {
         sdpisolver->rowmapper[2 * i] = nlpineqs + 1;
         sdpisolver->rowtoinputmapper[nlpineqs] = 2 * i;
         ++nlpineqs;
      }
      else
         sdpisolver->rowmapper[2 * i] = -1;

      assert( lprhs[i] > -SCIPsdpiSolverInfinity(sdpisolver) );
      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         sdpisolver->rowmapper[2 * i + 1] = nlpineqs + 1;
         sdpisolver->rowtoinputmapper[nlpineqs] = 2 * i + 1;
         ++nlpineqs;
      }
      else
         sdpisolver->rowmapper[2 * i + 1] = -1;
   }
   assert( nlpineqs <= 2 * nlpcons );
   sdpisolver->nlpineqs = nlpineqs;

   /* compute blockmapper and blockindmapper */
   ind = 0;
   assert( nsdpblocks <= sdpisolver->maxnsdpblocks );
   for (b = 0; b < nsdpblocks; ++b)
   {
      if ( blockindchanges[b] >= 0 )
      {
         int blockind = 0;

         sdpisolver->inputtoblockmapper[b] = ind + 1; /* +1 because SDPA start counting at one */

         for (i = 0; i < sdpblocksizes[b]; i++)
         {
            if ( indchanges[b][i] >= 0 )
            {
               assert( ind <= nsdpblocks );
               sdpisolver->blockindmapper[ind][blockind] = i;
               assert( i - indchanges[b][i] == blockind );
               blockind++;
            }
         }
         assert( blockind == sdpblocksizes[b] - nremovedinds[b] );
         ind++;
      }
      else
         sdpisolver->inputtoblockmapper[b] = -1;
   }

   /* initialize blockstruct */
   if ( penaltyparam < sdpisolver->epsilon ) /* we initialize this with an exact 0.0 in solve without penalty */
      sdpisolver->sdpa->inputConstraintNumber((SDPA_INT) sdpisolver->nactivevars);
   else
      sdpisolver->sdpa->inputConstraintNumber((SDPA_INT) sdpisolver->nactivevars + 1); /* the additional variable is r for the penalty formulation */

   /* if there are any LP-cons/variable-bounds, we get an extra block for those */
   if ( nlpineqs + sdpisolver->nvarbounds > 0 )
      sdpisolver->sdpa->inputBlockNumber((SDPA_INT) (nsdpablocks + 1));
   else
      sdpisolver->sdpa->inputBlockNumber((SDPA_INT) (nsdpablocks));

   /* initialize SDP blocks */
   for (b = 0; b < nsdpblocks; ++b)
   {
      if ( blockindchanges[b] >= 0 )
      {
         /* +1 because of SDPA counting */
         SCIPdebugMessage("adding block %d to SDPA as block %d with size %d.\n", b, b - blockindchanges[b] + 1, sdpblocksizes[b] - nremovedinds[b]);
         sdpisolver->sdpa->inputBlockSize((SDPA_INT) (b - blockindchanges[b] + 1), (SDPA_INT) (sdpblocksizes[b] - nremovedinds[b]));
         sdpisolver->sdpa->inputBlockType((SDPA_INT) (b - blockindchanges[b] + 1), SDPA::SDP);
      }
   }

   /* initialize LP block */
   if ( nlpineqs + sdpisolver->nvarbounds > 0 )
   {
      /* the last block is the LP block, the size has a negative sign */
      SCIPdebugMessage("adding LP block to SDPA as block %d with %d LP inequalities and %d variable bounds.\n", nsdpablocks + 1, nlpineqs, sdpisolver->nvarbounds);
      sdpisolver->sdpa->inputBlockSize((SDPA_INT) (nsdpablocks + 1), (SDPA_INT) -(nlpineqs + sdpisolver->nvarbounds));
      sdpisolver->sdpa->inputBlockType((SDPA_INT) (nsdpablocks + 1), SDPA::LP);
   }

   /* initialize space */
   sdpisolver->sdpa->initializeUpperTriangleSpace();

   /* set objective values */
   if ( withobj )
   {
      for (i = 0; i < sdpisolver->nactivevars; i++)
      {
         /* insert objective value, SDPA counts from 1 to n instead of 0 to n-1 */
         sdpisolver->sdpa->inputCVec((SDPA_INT) (i + 1), obj[sdpisolver->sdpatoinputmapper[i]]);
      }
   }

   if ( penaltyparam >= sdpisolver->epsilon )
      sdpisolver->sdpa->inputCVec((SDPA_INT) (sdpisolver->nactivevars + 1), penaltyparam); /* set the objective of the additional var to penaltyparam */

   /* inserting non-constant SDP constraint matrices */
   if ( sdpnnonz > 0 )
   {
      int blockvar;
      int v;

      assert( nsdpblocks > 0 );
      assert( sdpnblockvarnonz != NULL );
      assert( sdpnblockvars != NULL );
      assert( sdpcol != NULL );
      assert( sdprow != NULL );
      assert( sdpval != NULL );
      assert( sdpvar != NULL );
      assert( indchanges != NULL );
      assert( nremovedinds != NULL );

      for (b = 0; b < nsdpblocks; ++b)
      {
         /* if the block has no entries, we skip it */
         if ( blockindchanges[b] < 0 )
            continue;

         assert( b - blockindchanges[b] + 1 <= nsdpablocks );

         /* iterate over all variables in this block */
         for (blockvar = 0; blockvar < sdpnblockvars[b]; blockvar++)
         {
            v = sdpisolver->inputtosdpamapper[sdpvar[b][blockvar]];

            /* check if the variable is active */
            if ( v > 0 )
            {
               assert( v <= sdpisolver->nactivevars );

               for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( 0 <= indchanges[b][sdprow[b][blockvar][k]] && indchanges[b][sdprow[b][blockvar][k]] <= sdprow[b][blockvar][k] );
                  assert( 0 <= indchanges[b][sdpcol[b][blockvar][k]] && indchanges[b][sdpcol[b][blockvar][k]] <= sdpcol[b][blockvar][k] );

                  assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                  assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                  /* rows and columns start with 1 in SDPA, so we have to add 1 to the indices */
                  sdpisolver->sdpa->inputElement((SDPA_INT) v, (SDPA_INT) (b - blockindchanges[b] + 1),
                     (SDPA_INT) (sdpcol[b][blockvar][k] - indchanges[b][sdpcol[b][blockvar][k]] + 1),
                     (SDPA_INT) (sdprow[b][blockvar][k] - indchanges[b][sdprow[b][blockvar][k]] + 1),
                     sdpval[b][blockvar][k], checkinput);
               }
            }
         }

         /* insert the identity matrix if we are using a penalty formulation */
         if ( penaltyparam >= sdpisolver->epsilon )
         {
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               sdpisolver->sdpa->inputElement((SDPA_INT) (sdpisolver->nactivevars + 1), (SDPA_INT) (b - blockindchanges[b] + 1),
                  (SDPA_INT) (i + 1), (SDPA_INT) (i + 1), 1.0, checkinput);
            }
         }
      }
   }

   /* inserting constant matrix */
   if ( sdpconstnnonz > 0 )
   {
      assert( nsdpblocks > 0 );
      assert( sdpconstnblocknonz!= NULL );
      assert( sdpconstcol != NULL );
      assert( sdpconstrow != NULL );
      assert( sdpconstval != NULL );

      for (b = 0; b < nsdpblocks; ++b)
      {
         if ( blockindchanges[b] < 0 )
            continue;

         for (k = 0; k < sdpconstnblocknonz[b]; k++)
         {
            /* rows and cols with active nonzeros should not be removed */
            assert( -1 < indchanges[b][sdpconstrow[b][k]] && indchanges[b][sdpconstrow[b][k]] <= sdpconstrow[b][k] );
            assert( -1 < indchanges[b][sdpconstcol[b][k]] && indchanges[b][sdpconstcol[b][k]] <= sdpconstcol[b][k] );

            assert( 0 <= sdpconstrow[b][k] && sdpconstrow[b][k] < sdpblocksizes[b] );
            assert( 0 <= sdpconstcol[b][k] && sdpconstcol[b][k] < sdpblocksizes[b] );

            /* rows and columns start with 1 in SDPA, so we have to add 1 to the indices, the constant matrix is given as variable 0 */
            sdpisolver->sdpa->inputElement((SDPA_INT) 0, (SDPA_INT) (b - blockindchanges[b] + 1),
               (SDPA_INT) (sdpconstcol[b][k] - indchanges[b][sdpconstcol[b][k]] + 1),
               (SDPA_INT) (sdpconstrow[b][k] - indchanges[b][sdpconstrow[b][k]] + 1),
               sdpconstval[b][k], checkinput);
         }
      }
   }

   /* inserting LP nonzeros */
   if ( lpnnonz > 0 )
   {
      SDPA_INT lhsrhscnt = 1;   /* 1 because SDPA starts with 1 */

      for (i = 0; i < nlpcons; ++i)
      {
         SCIP_Bool lhspresent = FALSE;
         SCIP_Bool rhspresent = FALSE;
         int nextbeg;
         int v;
         int j;

         if ( lpindchanges[i] < 0 )
            continue;

         /* determine presence of new lhs/rhs */
         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
            lhspresent = TRUE;

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
            rhspresent = TRUE;

         assert( lhspresent || rhspresent );
         assert( 0 <= lpbeg[i] && lpbeg[i] < lpnnonz );

         if ( i == nlpcons - 1 )
            nextbeg = lpnnonz;
         else
            nextbeg = lpbeg[i+1];

         for (j = lpbeg[i]; j < nextbeg; ++j)
         {
            assert( 0 <= lpind[j] && lpind[j] < nvars );
            assert( REALABS(lpval[j]) > sdpisolver->epsilon );

            /* if the variable is active */
            v = sdpisolver->inputtosdpamapper[lpind[j]];
            if ( v > 0 )
            {
               /* add the LP-nonzero to the lhs-inequality if it exists: */
               if ( lhspresent )
               {
                  /* LP nonzeros are added as diagonal entries of the last block */
                  sdpisolver->sdpa->inputElement((SDPA_INT) v, (SDPA_INT) (nsdpablocks + 1), lhsrhscnt, lhsrhscnt, lpval[j], checkinput);
               }

               /* add the LP-nonzero to the rhs-inequality if it exists: */
               if ( rhspresent )
               {
                  /* LP nonzeros are added as diagonal entries of the last block; negative value because this is a <=-constraint, while SDPA works with >= */
                  sdpisolver->sdpa->inputElement((SDPA_INT) v, (SDPA_INT) (nsdpablocks + 1), lhsrhscnt + lhspresent, lhsrhscnt + lhspresent, - lpval[j], checkinput);
               }
            }
         }

         /* move counter forward */
         if ( lhspresent )
         {
            /* LP constraints are added as diagonal entries of the last block, left-hand-side is added as variable zero */
            sdpisolver->sdpa->inputElement(0LL, (SDPA_INT) (nsdpablocks + 1), lhsrhscnt, lhsrhscnt, lplhs[i], checkinput);

            /* if we use a penalty formulation, add the r * Identity entry */
            if ( penaltyparam >= sdpisolver->epsilon )
            {
               /* LP nonzeros are added as diagonal entries of the last block, the r-variable is variable nactivevars + 1 */
               sdpisolver->sdpa->inputElement((SDPA_INT) (sdpisolver->nactivevars + 1), (SDPA_INT) (nsdpablocks + 1), lhsrhscnt, lhsrhscnt, 1.0, checkinput);
            }
            ++lhsrhscnt;
         }

         if ( rhspresent )
         {
            /* LP constraints are added as diagonal entries of the last block, right-hand side is added as variable zero with negative value */
            sdpisolver->sdpa->inputElement(0LL, (SDPA_INT) (nsdpablocks + 1), lhsrhscnt, lhsrhscnt, - lprhs[i], checkinput);

            /* if we use a penalty formulation, add the r * Identity entry */
            if ( penaltyparam >= sdpisolver->epsilon )
            {
               /* LP nonzeros are added as diagonal entries of the last block (coming after the last SDP-block, with
                * blocks starting at 1, as are rows), the r-variable is variable nactivevars + 1 */
               sdpisolver->sdpa->inputElement((SDPA_INT) (sdpisolver->nactivevars + 1), (SDPA_INT) (nsdpablocks + 1), lhsrhscnt, lhsrhscnt, 1.0, checkinput);
            }
            ++lhsrhscnt;
         }
      }
   }

   /* inserting variable bounds */
   if ( sdpisolver->nvarbounds > 0 )
   {
      SDPA_INT varboundcnt;
      int v;

      varboundcnt = (SDPA_INT) nlpineqs + 1;  /* +1 because of SDPA */

      for (i = 0; i < nvars; i++)
      {
         v = sdpisolver->inputtosdpamapper[i];
         if ( v > 0 )
         {
            assert( 0 < v && v <= sdpisolver->nactivevars );
            if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, lb[i]) )
            {
               /* add bound as LP constraint for this variable */
               sdpisolver->sdpa->inputElement((SDPA_INT) v, (SDPA_INT) (nsdpablocks + 1),  varboundcnt, varboundcnt, 1.0, checkinput);

               if ( REALABS(lb[i]) > sdpisolver->epsilon )
               {
                  /* the bound is added as the rhs and therefore variable zero */
                  sdpisolver->sdpa->inputElement(0LL, (SDPA_INT) (nsdpablocks + 1), varboundcnt, varboundcnt, lb[i], checkinput);
               }
               ++varboundcnt;
            }

            if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, ub[i]) )
            {
               /* add bound as LP constraint for this variable; - because of <= */
               sdpisolver->sdpa->inputElement((SDPA_INT) v, (SDPA_INT) (nsdpablocks + 1), varboundcnt, varboundcnt, -1.0, checkinput);

               if ( REALABS(ub[i]) > sdpisolver->epsilon )
               {
                  /* the bound is added as the rhs and therefore variable zero, we multiply by -1 for <= */
                  sdpisolver->sdpa->inputElement(0LL, (SDPA_INT) (nsdpablocks + 1), varboundcnt, varboundcnt, -ub[i], checkinput);
               }
               ++varboundcnt;
            }
         }
      }

      if ( penaltyparam >= sdpisolver->epsilon && rbound )
      {
         /* we add the variable bound r >= 0 */
         sdpisolver->sdpa->inputElement((SDPA_INT) (sdpisolver->nactivevars + 1), (SDPA_INT) (nsdpablocks + 1), varboundcnt, varboundcnt, 1.0, checkinput);
         ++varboundcnt;
      }
      assert( varboundcnt == nlpineqs + 1 + sdpisolver->nvarbounds );
   }


   /* initialize settings (this needs to be done before inserting the problem as the initial point depends on the settings) */
   if ( penaltyparam >= sdpisolver->epsilon || startsettings == SCIP_SDPSOLVERSETTING_STABLE || startsettings == SCIP_SDPSOLVERSETTING_PENALTY )
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW); /* if we already had problems with this problem, there is no reason to try fast */
      /* as we want to solve with stable settings, we also update epsilon and the feasibility tolerance, as we skip the default settings, we multpy twice */
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      feastol = FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol;
      gaptol = sdpisolver->gaptol;
      SCIPdebugMessage("Start solving process with stable settings ...\n");
   }
   else if ( startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST )
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_UNSTABLE_BUT_FAST);
      sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->sdpsolverfeastol);
      feastol = sdpisolver->sdpsolverfeastol;
      gaptol = sdpisolver->gaptol;
      SCIPdebugMessage("Start solving process with fast settings\n");
   }
   else if ( startsettings == SCIP_SDPSOLVERSETTING_MEDIUM )
   {
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_DEFAULT);
      /* as we want to solve with stable settings, we also update epsilon and the feasibility tolerance, as we skip the default settings, we multiply once */
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      feastol = FEASTOLCHANGE * sdpisolver->sdpsolverfeastol;
      gaptol = sdpisolver->gaptol;
      SCIPdebugMessage("Start solving process with medium settings ...\n");
   }
   else
   {
      SCIPdebugMessage("Unknown setting for start-settings: %d!\n", startsettings);  \
      return SCIP_LPERROR;
   }

   sdpisolver->sdpa->setParameterLowerBound(-1e20);

   /* set the objective limit */
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
   else
      sdpisolver->sdpa->setParameterUpperBound(1e8);

   /* increase Lambda Star, this seems to help the numerics */
   sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);
   //sdpisolver->sdpa->setParameterBetaStar(0.01);
   //sdpisolver->sdpa->setParameterBetaBar(0.02);

#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpa->printParameters(stdout);
#endif

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
   if ( starty != NULL && startZnblocknonz != NULL && startXnblocknonz != NULL && penaltyparam < sdpisolver->epsilon )
      sdpisolver->sdpa->setInitPoint(true);
   else
      sdpisolver->sdpa->setInitPoint(false);

   /* set the starting solution */
   if ( starty != NULL && startZnblocknonz != NULL && startXnblocknonz != NULL && penaltyparam < sdpisolver->epsilon )
   {
      int sdpaind;
      int varboundind;

      /* set dual vector */
      for (i = 1; i <= sdpisolver->nactivevars; i++)
         sdpisolver->sdpa->inputInitXVec((SDPA_INT) i, starty[sdpisolver->sdpatoinputmapper[i - 1]]);/*lint !e747*/

      /* set SDP-blocks of dual matrix */
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] >= 0 )
         {
            for (i = 0; i < startZnblocknonz[b]; i++)
            {
               if ( indchanges[b][startZrow[b][i]] > -1 && indchanges[b][startZcol[b][i]] > -1 ) /* the row/col may have been fixed to zero in the meantime */
               {
                  assert( ! SCIPsdpiSolverIsInfinity(sdpisolver, startZval[b][i]) );
                  sdpisolver->sdpa->inputInitXMat((SDPA_INT) b + 1 - blockindchanges[b], (SDPA_INT) startZrow[b][i] + 1 - indchanges[b][startZrow[b][i]],
                     (SDPA_INT) startZcol[b][i] + 1 - indchanges[b][startZcol[b][i]], startZval[b][i]);
               }
            }
         }
      }

      /* set LP-block of dual matrix */
      for (i = 0; i < startZnblocknonz[nsdpblocks]; i++)
      {
         assert( startZrow[nsdpblocks][i] == startZcol[nsdpblocks][i] );
         assert( 0 <= startZrow[nsdpblocks][i] && startZrow[nsdpblocks][i] < 2 * nlpcons + 2 * nvars );

         if ( startZrow[nsdpblocks][i] < 2 * nlpcons ) /* linear constraint */
         {
            int rowidx;

            rowidx = startZrow[nsdpblocks][i];
            assert( 0 <= rowidx && rowidx < 2 * nlpcons );

            sdpaind = sdpisolver->rowmapper[rowidx];
            if ( sdpaind >= 0 )
            {
               assert( ! SCIPsdpiSolverIsInfinity(sdpisolver, startZval[nsdpblocks][i]) );
               sdpisolver->sdpa->inputInitXMat((SDPA_INT) nsdpblocks + 1, sdpaind, sdpaind, startZval[nsdpblocks][i]);
            }
         }
         else /* varbound */
         {
            int varidx;

            varidx = startZrow[nsdpblocks][i] - 2 * nlpcons;
            assert( 0 <= varidx && varidx < 2 * nvars );
            varboundind = sdpisolver->inputtovbmapper[varidx];
            if ( varboundind > -1 ) /* rhs */
            {
               assert( ! SCIPsdpiSolverIsInfinity(sdpisolver, startZval[nsdpblocks][i]) );
               sdpisolver->sdpa->inputInitXMat((SDPA_INT) nsdpblocks + 1, nlpineqs + varboundind + 1,
                  nlpineqs + varboundind + 1, startZval[nsdpblocks][i]);
            }
         }
      }

      /* set SDP-blocks of primal matrix */
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] >= 0 )
         {
            for (i = 0; i < startXnblocknonz[b];  i++)
            {
               if ( indchanges[b][startXrow[b][i]] > -1 && indchanges[b][startXcol[b][i]] > -1 ) /* the row/col may have been fixed to zero in the meantime */
               {
                  sdpisolver->sdpa->inputInitYMat((SDPA_INT) b + 1 - blockindchanges[b], (SDPA_INT) startXrow[b][i] + 1 - indchanges[b][startXrow[b][i]],
                     (SDPA_INT) startXcol[b][i] + 1 - indchanges[b][startXcol[b][i]], startXval[b][i]);
               }
            }
         }
      }

      /* set LP-block of primal matrix */
      for (i = 0; i < startXnblocknonz[nsdpblocks]; i++)
      {
         assert( startXrow[nsdpblocks][i] == startXcol[nsdpblocks][i] );

         if ( startXrow[nsdpblocks][i] < 2 * nlpcons ) /* linear constraint */
         {
            int rowidx;

            rowidx = startXrow[nsdpblocks][i];
            assert( 0 <= rowidx && rowidx < 2 * nlpcons );

            sdpaind = sdpisolver->rowmapper[rowidx];
            if ( sdpaind >= 0 )
               sdpisolver->sdpa->inputInitYMat((SDPA_INT) nsdpblocks + 1, sdpaind, sdpaind, startXval[nsdpblocks][i]);
         }
         else /* varbound */
         {
            int varidx;

            varidx = startXrow[nsdpblocks][i] - 2 * nlpcons;
            assert( 0 <= varidx && varidx < 2 * nvars );
            varboundind = sdpisolver->inputtovbmapper[varidx];
            if ( varboundind > -1 ) /*lhs */
            {
               sdpisolver->sdpa->inputInitYMat((SDPA_INT) nsdpblocks + 1, nlpineqs + varboundind + 1,
                     nlpineqs + varboundind + 1, startXval[nsdpblocks][i]);
            }
         }
      }
   }
   else if ( penaltyparam >= sdpisolver->epsilon )
      SCIPdebugMessage("Skipping insertion of starting point for penalty formulation.\n");

   /* transform the matrices to a more efficient form */
   sdpisolver->sdpa->initializeUpperTriangle(checkinput);
   sdpisolver->sdpa->initializeSolve();

#ifdef SCIP_DEBUG_PRINTTOFILE
   /* if necessary, dump input data and initial point */
   sdpisolver->sdpa->writeInputSparse(const_cast<char*>("sdpa.dat-s"), const_cast<char*>("%+8.3e"));
   sdpisolver->sdpa->writeInitSparse(const_cast<char*>("sdpa.ini-s"), const_cast<char*>("%+8.3e"));
#endif

   SCIPdebugMessage("Calling SDPA solve (SDP: %d) ...\n", sdpisolver->sdpcounter);

   /* check if we want to save a preoptimal solution for warmstarting purposes */
   if ( usepreoptimalsol && (startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST) )
   {
      SCIP_Real* sdpasol;
      int sdpablocksize;

      /* first solve up to gaptol */
      sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->preoptimalgap);
      sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->preoptimalgap);
      sdpisolver->sdpa->solve();

      /* save preoptimal solution (only if solving succeeded) */
      if ( SCIPsdpiSolverIsAcceptable(sdpisolver) )
      {
         SCIPdebugMessage("Saving preoptimal solution for warmstarting purposes.\n");
         sdpasol =  sdpisolver->sdpa->getResultXVec();

         /* copy dual solution vector */
         for (i = 0; i < sdpisolver->nactivevars; i++)
            sdpisolver->preoptimalsol[i] = sdpasol[i];

         /* copy primal matrix */
         for (b = 0; b < sdpisolver->sdpa->getBlockNumber(); b++)
         {
            sdpasol = sdpisolver->sdpa->getResultYMat(b + 1);
            sdpablocksize = (int) sdpisolver->sdpa->getBlockSize(b + 1);
            assert( sdpablocksize <= 2 * (sdpisolver->maxnlpcons + sdpisolver->maxnvars) );
            if ( sdpisolver->sdpa->getBlockType(b + 1) == SDPA::LP )
            {
               /* only the diagonal entries need to be saved */
               for (i = 0; i < sdpablocksize; i++)
                  sdpisolver->preoptimalsolxlp[i] = sdpasol[i];
            }
            else
            {
               /* TODO: think about saving this in sparse format or at least only the upper triangular matrix */
               for (i = 0; i < sdpablocksize * sdpablocksize; i++)
                  sdpisolver->preoptimalsolx[b][i] = sdpasol[i];
            }
         }

         sdpisolver->preoptimalsolexists = TRUE;
      }
      else
      {
         SCIPdebugMessage("Solving to gaptol failed! No preoptimal solution available.\n");
      }

      /* add number of iterations */
      sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();

      /* copy current iterate and resolve with real tolerance */
      sdpisolver->sdpa->setInitPoint(true);
      sdpisolver->sdpa->copyCurrentToInit();
      sdpisolver->sdpa->setParameterEpsilonStar(sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(sdpisolver->sdpsolverfeastol);
      sdpisolver->sdpa->solve();
   }
   else
      sdpisolver->sdpa->solve();

   sdpisolver->solved = TRUE;

   /* update number of SDP-iterations and -calls */
   sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
   ++sdpisolver->nsdpcalls;

#ifdef SCIP_DEBUG
   /* print the phase value , i.e. whether solving was successfull */
   sdpisolver->sdpa->getPhaseString((char*)phase_string);
   SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif

   /* remember settings */
   if ( SCIPsdpiSolverIsAcceptable(sdpisolver) && penaltyparam < sdpisolver->epsilon )
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_FAST;
   else if ( penaltyparam >= sdpisolver->epsilon )
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;

   /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
    * reports feasibility, resolve it with adjusted tolerance */
   SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
         sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
         nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval, &feastol, &gaptol) );

   /* check whether problem has been stably solved, if it wasn't and we didn't yet use the default parametersettings (for the penalty formulation we do so), try
    * again with more stable parameters */
   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) && penaltyparam < sdpisolver->epsilon &&
      (startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST) )
   {
      SCIPdebugMessage("Numerical troubles -- solving SDP %d again ...\n", sdpisolver->sdpcounter);

      /* initialize settings */
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_DEFAULT);
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      sdpisolver->sdpa->setParameterLowerBound(-1e20);

      /* set the objective limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
      else
         sdpisolver->sdpa->setParameterUpperBound(1e8);

      /* increase Lambda Star, this seems to help the numerics */
      sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);

#ifdef SCIP_MORE_DEBUG
      sdpisolver->sdpa->printParameters(stdout);
#endif
      sdpisolver->sdpa->setInitPoint(false);
      sdpisolver->sdpa->resetParameters();
      sdpisolver->sdpa->solve();
      sdpisolver->solved = TRUE;

      /* update number of SDP-iterations and -calls */
      sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
      ++sdpisolver->nsdpcalls;

      /* remember setting */
      if ( SCIPsdpiSolverIsAcceptable(sdpisolver) )
         sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_MEDIUM;

#ifdef SCIP_DEBUG
      /* print the phase value , i.e. whether solving was successfull */
      sdpisolver->sdpa->getPhaseString((char*)phase_string);
      SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in contrast to our formulation)\n", phase_string);
#endif

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
      SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
            nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval, &feastol, &gaptol) );
   }

   /* if we still didn't converge, and did not yet use the stable settings, set the parameters even more conservativly */
   if ( (! SCIPsdpiSolverIsAcceptable(sdpisolver)) && penaltyparam < sdpisolver->epsilon &&
      (startsettings == SCIP_SDPSOLVERSETTING_UNSOLVED || startsettings == SCIP_SDPSOLVERSETTING_FAST || startsettings == SCIP_SDPSOLVERSETTING_MEDIUM) )
   {
      SCIPdebugMessage("Numerical troubles -- solving SDP %d again^2 ...\n", sdpisolver->sdpcounter);

      /* initialize settings */
      sdpisolver->sdpa->setParameterType(SDPA::PARAMETER_STABLE_BUT_SLOW);
      sdpisolver->sdpa->setParameterEpsilonStar(GAPTOLCHANGE * GAPTOLCHANGE * sdpisolver->gaptol);
      sdpisolver->sdpa->setParameterEpsilonDash(FEASTOLCHANGE * FEASTOLCHANGE * sdpisolver->sdpsolverfeastol);
      sdpisolver->sdpa->setParameterLowerBound(-1e20);

      /* set the objective limit */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         sdpisolver->sdpa->setParameterUpperBound(sdpisolver->objlimit);
      else
         sdpisolver->sdpa->setParameterUpperBound(1e8);

      /* increase Lambda Star, this seems to help the numerics */
      sdpisolver->sdpa->setParameterLambdaStar(sdpisolver->lambdastar);

#ifdef SCIP_MORE_DEBUG
      sdpisolver->sdpa->printParameters(stdout);
#endif
      sdpisolver->sdpa->setInitPoint(false);
      sdpisolver->sdpa->resetParameters();
      sdpisolver->sdpa->solve();
      sdpisolver->solved = TRUE;

      /* update number of SDP-iterations and -calls */
      sdpisolver->niterations += (int) sdpisolver->sdpa->getIteration();
      ++sdpisolver->nsdpcalls;

      /* remember setting */
      if ( SCIPsdpiSolverIsAcceptable(sdpisolver) )
         sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_STABLE;

#ifdef SCIP_DEBUG
      /* print the phase value , i.e. whether solving was successfull */
      sdpisolver->sdpa->getPhaseString((char*)phase_string);
      SCIPdebugMessage("SDPA solving finished with status %s (primal and dual here are switched in constrast to our formulation)\n", phase_string);
#endif

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
      SCIP_CALL( checkFeastolAndResolve(sdpisolver, penaltyparam, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, indchanges,
            nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval, &feastol, &gaptol) );
   }

#ifdef SCIP_MORE_DEBUG
   (void) fclose(fpOut);
#endif

   /* if we solved a penalty formulation, check if the solution is feasible for the original problem (which is the case iff r < feastol) */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      SCIP_Real* sdpasol;
      SCIP_Real* X;
      int nblockssdpa;
      int nrow;
      SCIP_Real trace = 0.0;

      /* in the second case we have r as an additional variable */
      assert( (sdpisolver->nactivevars + 1 == sdpisolver->sdpa->getConstraintNumber()) );  /*lint !e776*/

      sdpasol = sdpisolver->sdpa->getResultXVec();

      /* we get r as the last variable in SDPA */
      assert( feasorig != NULL );
      if ( sdpasol[sdpisolver->nactivevars] < sdpisolver->feastol )
         *feasorig = TRUE;
      else
         *feasorig = FALSE;

      /* only set sdpisolver->feasorig to true if we solved with objective, because only in this case we want to compute
       * the objective value by hand since it is numerically more stable then the result returned by SDPA */
      if ( withobj )
         sdpisolver->feasorig = *feasorig;

      /* if r > 0 or we are in debug mode, also check the primal bound */
#ifdef NDEBUG
      if ( ! *feasorig && penaltybound != NULL )
      {
#endif
         SCIPdebugMessage("Solution not feasible in original problem, r = %f\n", sdpasol[sdpisolver->nactivevars]);

         /* compute Tr(X) */

         /* iterate over all blocks (SDPA starts counting at one and includes the LP block) */
         nblockssdpa = (int) sdpisolver->sdpa->getBlockNumber();
         for (b = 1; b <= nblockssdpa; b++)
         {
            /* get the block from SDPA */
            X = sdpisolver->sdpa->getResultYMat((SDPA_INT) b);/*lint !e747*/
            nrow = (int) sdpisolver->sdpa->getBlockSize((SDPA_INT) b);/*lint !e747*/
            assert( nrow >= 0 );

            /* if it is the LP-block, we omit the variable bounds as the penalty variable is not added to them */
            if ( sdpisolver->sdpa->getBlockType((SDPA_INT) b) == SDPA::LP )/*lint !e747*/
            {
               /* iterate over all diagonal entries (until we reach the varbound part), adding them to the trace */
               for (i = 0; i < nrow - sdpisolver->nvarbounds; i++)
                  trace += X[i]; /* get entry (i+1,i+1) for the diagonal matrix X */
            }
            else
            {
               /* iterate over all diagonal entries and add them to the trace */
               for (i = 0; i < nrow; i++)
                  trace += X[i + i*nrow];/*lint !e679*/ /* get entry (i+1,i+1) in X */
            }
         }

         /* if the relative gap is smaller than the tolerance, we return equality */
         if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )/*lint !e414*/
         {
            if ( penaltybound != NULL )
               *penaltybound = TRUE;

            SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible\n",
               trace, penaltyparam);
         }
         else if ( penaltybound != NULL )
            *penaltybound = FALSE;

#ifdef NDEBUG
      }
#endif
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
      return FALSE;

   return TRUE;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   switch ( phasetype )/*lint --e{788}*/
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
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;
   case SDPA::pINF_dFEAS:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;
   case SDPA::pUNBD:
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      SCIPdebugMessage("SDPA stopped because dual objective became smaller than lower bound\n");
      break;
   case SDPA::dUNBD:
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      SCIPdebugMessage("SDPA stopped because primal objective became bigger than upper bound\n");
      break;
   default: /* contains noInfo, pFeas, dFeas, pdInf */
      SCIPerrorMessage("SDPA doesn't know if primal and dual solutions are feasible\n");
      return SCIP_LPERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
*/
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
*/
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is infeasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
*/
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
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

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
*/
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is unbounded");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
*/
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if dual problem is infeasible.\n");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pFEAS_dINF )
      return TRUE;
   else if ( phasetype == SDPA::pUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because primal objective became bigger than upper bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::noINFO || phasetype == SDPA::dFEAS || phasetype == SDPA::pdINF )
   {
      SCIPdebugMessage("SDPA doesn't know if primal problem is feasible");
      return FALSE;
   }
   else if ( phasetype ==  SDPA::pINF_dFEAS || phasetype == SDPA::pdOPT || phasetype == SDPA::dFEAS  || phasetype == SDPA::pdFEAS )
      return TRUE;
   else if ( phasetype == SDPA::dUNBD )
   {
      SCIPdebugMessage("SDPA was stopped because dual objective became smaller than lower bound");
      return TRUE;
   }

   return FALSE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pUNBD )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL);
   CHECK_IF_SOLVED_BOOL( sdpisolver );

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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   return sdpisolver->timelimit;
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( sdpisolver->sdpa == NULL || (! sdpisolver->solved) )
      return -1;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT || phasetype == SDPA::pFEAS_dINF || phasetype == SDPA::pINF_dFEAS )
      return 0;
   if ( phasetype == SDPA::pdINF )
      return 1;
   if ( phasetype == SDPA::pUNBD )
      return 3;
   if ( phasetype == SDPA::noINFO || phasetype == SDPA::pFEAS || phasetype == SDPA::dFEAS || phasetype == SDPA::pdFEAS )
      return 4;
   else /* should include dUNBD */
      return 7;
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( ! sdpisolver->solved )
      return FALSE;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   if ( phasetype == SDPA::pdOPT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework
 */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{/*lint !e1784*/
   SDPA::PhaseType phasetype;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   phasetype = sdpisolver->sdpa->getPhaseValue();

   /* we are happy if we converged, or we reached the objective limit (pUNBD) or we could show that our problem is
    * infeasible [except for numerics], or unbounded */
   if ( SCIPsdpiSolverIsConverged(sdpisolver) || phasetype == SDPA::pUNBD || phasetype == SDPA::pINF_dFEAS || phasetype == SDPA::pFEAS_dINF )
      return TRUE;

   return FALSE;
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{/*lint !e1784*/
   SCIPdebugMessage("Not implemented yet\n");

   return SCIP_LPERROR;
}/*lint !e715*/

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( sdpisolver->penalty && ! sdpisolver->feasorig )
   {
      *objval = sdpisolver->sdpa->getPrimalObj();

#ifndef NDEBUG
      SCIP_Real primalval = sdpisolver->sdpa->getDualObj();
      SCIP_Real gap = REALABS(*objval - primalval) / (0.5 * (REALABS(primalval) + REALABS(*objval))); /* duality gap used in SDPA */
      if ( gap > sdpisolver->gaptol )
      {
         SCIPdebugMessage("Attention: got objective value (before adding values of fixed variables) of %g in SCIPsdpiSolverGetSol, "
            "but primal objective is %g with duality gap %g!\n", *objval, primalval, gap );
      }
#endif
   }
   else
   {
      SCIP_Real* sdpasol;
      int v;

      /* since the objective value given by SDPA sometimes differs slightly from the correct value for the given solution,
       * we get the solution from SDPA and compute the correct objective value */
      assert( (sdpisolver->penalty && sdpisolver->nactivevars + 1 == sdpisolver->sdpa->getConstraintNumber()) || /*lint !e776*/
         sdpisolver->nactivevars == sdpisolver->sdpa->getConstraintNumber() ); /* in the second case we have r as an additional variable */

      sdpasol = sdpisolver->sdpa->getResultXVec();

      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
         *objval += sdpasol[v] * sdpisolver->objcoefs[v];
   }

   /* add contribution of fixed variables */
   *objval += sdpisolver->fixedvarsobjcontr;

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs */
SCIP_RETCODE SCIPsdpiSolverGetDualSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval,             /**< pointer to store the objective value (or NULL) */
   SCIP_Real*            dualsol             /**< array of length nvars to store the dual solution vector (or NULL) */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( objval != NULL )
   {
      SCIP_CALL( SCIPsdpiSolverGetObjval(sdpisolver, objval) );
   }

   if ( dualsol != NULL )
   {
      SCIP_Real* sdpasol;
      int v;

      /* get the solution from sdpa */
      assert( (sdpisolver->penalty && sdpisolver->nactivevars + 1 == sdpisolver->sdpa->getConstraintNumber()) || /*lint !e776*/
         sdpisolver->nactivevars == sdpisolver->sdpa->getConstraintNumber() ); /* in the second case we have r as an additional variable */

      sdpasol = sdpisolver->sdpa->getResultXVec();

      /* insert the entries into dualsol, for non-fixed vars we copy those from SDPA */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtosdpamapper[v] > 0 )
         {
            /* minus one because the inputtosdpamapper starts at 1, but the array starts at 0 */
            dualsol[v] = sdpasol[sdpisolver->inputtosdpamapper[v] - 1];
         }
         else
         {
            /* the fixed value was saved at the beginning */
            assert( - sdpisolver->inputtosdpamapper[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[- sdpisolver->inputtosdpamapper[v] - 1]; /*lint !e679*/
         }
      }
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
{
   int b;
   int sdpablock;
   int sdpaind;
   int i;
   int c;
   int r;
   int blocksize;

   assert( sdpisolver != NULL );
   assert( nblocks > 0 );
   assert( startXnblocknonz != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( ! sdpisolver->preoptimalsolexists )
   {
      SCIPdebugMessage("Failed to retrieve preoptimal solution for warmstarting purposes.\n");
      startXnblocknonz[0] = -1;
      return SCIP_OKAY;
   }

   assert( nblocks <= sdpisolver->nsdpblocks + 1 );
   assert( nblocks == sdpisolver->sdpa->getBlockNumber() );

   /* iterate over all SDP-blocks, get the solution and count the nonzeros */
   for (b = 0; b < sdpisolver->nsdpblocks; b++)
   {
      sdpablock = sdpisolver->inputtoblockmapper[b];

      startXnblocknonz[b] = 0;

      if ( sdpablock != -1 )
      {
         blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);

         /* iterate once over the upper triangular part of the matrix (saving the corresponding entries in the lower triangular part for the SDPI) */
         for (r = 0; r < blocksize; r++)
         {
            for (c = r; c < blocksize; c++)
            {
               sdpaind = r + blocksize * c;

               if ( REALABS(sdpisolver->preoptimalsolx[sdpablock - 1][sdpaind]) > sdpisolver->epsilon )
                  startXnblocknonz[b]++;
            }
         }
      }
   }

   /* compute the entry for the LP-block */
   startXnblocknonz[sdpisolver->nsdpblocks] = 0;
   if ( sdpisolver->sdpa->getBlockType(sdpisolver->nsdpblocks) == SDPA::LP )
   {
      blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpisolver->nsdpblocks);

      for (i = 0; i < blocksize; i++)
      {
         if ( REALABS(sdpisolver->preoptimalsolxlp[i]) > sdpisolver->epsilon )
            startXnblocknonz[sdpisolver->nsdpblocks]++;
      }
   }

   return SCIP_OKAY;
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
{
   int v;
   int b;
   int r;
   int c;
   int sdpaind;
   int sdpablock;
   int blocksize;

   assert( sdpisolver != NULL );
   assert( success != NULL );
   assert( dualsol != NULL );
   assert( startXnblocknonz != NULL || nblocks == -1 );
   assert( startXrow != NULL || nblocks == -1 );
   assert( startXcol != NULL || nblocks == -1 );
   assert( startXval != NULL || nblocks == -1 );

   if ( ! sdpisolver->preoptimalsolexists )
   {
      SCIPdebugMessage("Failed to retrieve preoptimal solution for warmstarting purposes. \n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (v = 0; v < sdpisolver->nvars; v++)
   {
      sdpaind = sdpisolver->inputtosdpamapper[v];
      if ( sdpaind > 0)
      {
         /* minus one because the array starts at 0 */
         dualsol[v] = sdpisolver->preoptimalsol[sdpaind - 1];
      }
      else
      {
         /* this is the value that was saved when inserting, as this variable has lb=ub */
         dualsol[v] = sdpisolver->fixedvarsval[- sdpaind - 1]; /*lint !e679*/
      }
   }

   if ( nblocks > -1 )
   {
      int blocknnonz = 0;

      assert( startXnblocknonz != NULL );
      assert( startXrow != NULL );
      assert( startXcol != NULL );
      assert( startXval != NULL );
      assert( nblocks <= sdpisolver->nsdpblocks + 1 );
      assert( nblocks == sdpisolver->sdpa->getBlockNumber() );

      /* iterate over all SDP-blocks and get the solution */
      for (b = 0; b < sdpisolver->nsdpblocks; b++)
      {
         assert( startXrow[b] != NULL );
         assert( startXcol[b] != NULL );
         assert( startXval[b] != NULL );

         sdpablock = sdpisolver->inputtoblockmapper[b];

         if ( sdpablock >= 0 )
         {
            /* since we reset the preoptimalsolution for every solve, the blocksize should have stayed the same */
            blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);

            /* iterate once over the upper triangular part of the matrix (saving the corresponding entries in the lower triangular part for the SDPI) */
            for (r = 0; r < blocksize; r++)
            {
               for (c = r; c < blocksize; c++)
               {
                  sdpaind = r + blocksize * c;

                  if ( REALABS(sdpisolver->preoptimalsolx[sdpablock - 1][sdpaind]) > sdpisolver->epsilon )
                  {
                     assert( blocknnonz < startXnblocknonz[b] );
                     startXrow[b][blocknnonz] = sdpisolver->blockindmapper[b][c];
                     startXcol[b][blocknnonz] = sdpisolver->blockindmapper[b][r];
                     startXval[b][blocknnonz] = sdpisolver->preoptimalsolx[sdpablock - 1][sdpaind]; /* -1 because sdpa starts counting at 1 */
                     blocknnonz++;
                  }
               }
            }

            startXnblocknonz[b] = blocknnonz;
         }
         else
            startXnblocknonz[b] = 0;
      }

      /* compute entries for the LP-block */
      if ( sdpisolver->sdpa->getBlockType(sdpisolver->nsdpblocks) == SDPA::LP )
      {
         int i;

         /* since we reset the preoptimalsolution for every solve, the blocksize should have stayed the same */
         blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpisolver->nsdpblocks);
         assert( blocksize == sdpisolver->nvarbounds + sdpisolver->nlpineqs );

         /* iterate over LP constraints */
         assert( b == sdpisolver->nsdpblocks );
         blocknnonz = 0;
         for (i = 0; i < sdpisolver->nlpineqs; i++)
         {
            if ( REALABS(sdpisolver->preoptimalsolxlp[i]) > sdpisolver->epsilon )
            {
               assert( blocknnonz < startXnblocknonz[b] );

               startXrow[b][blocknnonz] = sdpisolver->rowtoinputmapper[i];
               startXcol[b][blocknnonz] = sdpisolver->rowtoinputmapper[i];
               startXval[b][blocknnonz] = sdpisolver->preoptimalsolxlp[i];
               assert( 0 <= startXrow[b][blocknnonz] && startXrow[b][blocknnonz] < 2 * sdpisolver->nlpcons );
               blocknnonz++;
               assert( blocknnonz <= startXnblocknonz[b] );
            }
         }

         /* iterate over varbounds */
         for (i = 0; i < sdpisolver->nvarbounds; i++)
         {
            SCIP_Real val;

            val = sdpisolver->preoptimalsolxlp[sdpisolver->nlpineqs + i];

            if ( REALABS(val) > sdpisolver->epsilon )
            {
               int inputpos;
               int vbpos;

               vbpos = sdpisolver->varboundpos[i];
               if ( vbpos  > 0 ) /* rhs */
                  inputpos = 2 * sdpisolver->nlpcons + 2 * sdpisolver->sdpatoinputmapper[vbpos - 1] + 1;
               else
                  inputpos = 2 * sdpisolver->nlpcons + 2 * sdpisolver->sdpatoinputmapper[- vbpos - 1];

               startXrow[b][blocknnonz] = inputpos;
               startXcol[b][blocknnonz] = inputpos;
               startXval[b][blocknnonz] = val;
               assert( 2 * sdpisolver->nlpcons <= inputpos && inputpos < 2 * sdpisolver->nlpcons + 2 * sdpisolver->nvars );
               blocknnonz++;
               assert( blocknnonz <= startXnblocknonz[b] );
            }
         }
         startXnblocknonz[b] = blocknnonz;
      }
      else
         startXnblocknonz[b] = 0;
   }

   *success = TRUE;

   return SCIP_OKAY;
}

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
{  /*lint !e1784*/
   SCIP_Real* X; /* block of primal solution matrix corresponding to the LP-part */
   SDPA_INT lpblockind;
   int penoffset = 0;
   int nlpcons;
   int sdpapos;
   int pos;
   int i;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( lbvals != NULL );
   assert( ubvals != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   /* initialize the return-arrays with zero */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      lbvals[i] = 0.0;
      ubvals[i] = 0.0;
   }

   /* if no variable bounds were added, we return the zero vector (we do this separately, because in this case there might be no LP-block) */
   if ( sdpisolver->nvarbounds == 0 )
   {
      SCIPdebugMessage("Asked for PrimalBoundVars, but there were no variable bounds in SDPA, returning zero vector!\n");
      return SCIP_OKAY;
   }

   /* get the block of primal solution matrix corresponding to the LP-part from SDPA */
   lpblockind = (int) sdpisolver->sdpa->getBlockNumber(); /* the LP block is the last one and SDPA counts from one */
   if ( sdpisolver->sdpa->getBlockType(lpblockind) != SDPA::LP )
   {
      /* if the last block is not an LP-block, no variable bounds existed */
      return SCIP_OKAY;
   }
   nlpcons = (int) sdpisolver->sdpa->getBlockSize(lpblockind);
   assert( nlpcons >= 0 );

   X = sdpisolver->sdpa->getResultYMat(lpblockind);

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars || (sdpisolver->nvarbounds <= 2 * sdpisolver->nvars + 1 && sdpisolver->penalty ) );

   /* if we solved a penalty formulation, the last variable bound belongs to the penalty variable, which we aren't interested in here */
   if ( sdpisolver->penalty )
      penoffset = 1;

   for (i = 0; i < sdpisolver->nvarbounds - penoffset; i++)
   {
      /* if it is a lower bound */
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         sdpapos = - sdpisolver->varboundpos[i] - 1;
         assert( 0 <= sdpapos && sdpapos < sdpisolver->nactivevars );
         pos = sdpisolver->sdpatoinputmapper[sdpapos];
         assert( 0 <= pos && pos < sdpisolver->nvars );
         /* the last nvarbounds entries correspond to the varbounds */
         lbvals[pos] = X[nlpcons - sdpisolver->nvarbounds + i];
      }
      else /* if it is an upper bound */
      {
         sdpapos = sdpisolver->varboundpos[i] - 1;
         assert( 0 <= sdpapos && sdpapos < sdpisolver->nactivevars );
         pos = sdpisolver->sdpatoinputmapper[sdpapos];
         assert( 0 <= pos && pos < sdpisolver->nvars );
         /* the last nvarbounds entries correspond to the varbounds */
         ubvals[pos] = X[nlpcons - sdpisolver->nvarbounds + i]; /*lint !e679, !e834 */
      }
   }

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
{  /*lint --e{715}*/
   SCIP_Real* X; /* block of primal solution matrix corresponding to the LP-part */
   SDPA_INT lpblockind;
   int pos = 0;
   int i;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( lpindchanges != NULL );
   assert( lplhs != NULL );
   assert( lprhs != NULL );
   assert( lhsvals != NULL );
   assert( rhsvals != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( nlpcons <= 0 )
      return SCIP_OKAY;

   /* get the block of primal solution matrix corresponding to the LP-part from SDPA */
   lpblockind = (int) sdpisolver->sdpa->getBlockNumber(); /* the LP block is the last one and SDPA counts from one */
   if ( sdpisolver->sdpa->getBlockType(lpblockind) != SDPA::LP )
   {
      /* if the last block is not an LP-block, no variable bounds existed */
      return SCIP_OKAY;
   }

   X = sdpisolver->sdpa->getResultYMat(lpblockind);

   /* loop through LP rows */
   for (i = 0; i < nlpcons; i++)
   {
      if ( lpindchanges[i] < 0 )
      {
         lhsvals[i] = 0.0;
         rhsvals[i] = 0.0;
         continue;
      }

      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         lhsvals[i] = X[pos++];
      else
         lhsvals[i] = 0.0;

      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         rhsvals[i] = X[pos++];
      else
         rhsvals[i] = 0.0;
   }

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X (including lp block) */
SCIP_RETCODE SCIPsdpiSolverGetPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   )
{
   SCIP_Real* X;
   int b;
   int sdpablock;
   int sdpaind;
   int i;
   int r;
   int c;
   int blocksize;

   assert( sdpisolver != NULL );
   assert( nblocks > 0 );
   assert( startXnblocknonz != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( nblocks != sdpisolver->nsdpblocks + 1 )
   {
      SCIPerrorMessage("SCIPsdpiSolverGetPrimalNonzeros expected nblocks = %d but got %d\n", sdpisolver->nsdpblocks + 1, nblocks);
      return SCIP_LPERROR;
   }

   /* iterate over all SDP-blocks, get the solution and count the nonzeros */
   for (b = 0; b < sdpisolver->nsdpblocks; b++)
   {
      sdpablock = sdpisolver->inputtoblockmapper[b];

      startXnblocknonz[b] = 0;

      if ( sdpablock != -1 )
      {
         X = sdpisolver->sdpa->getResultYMat(sdpablock);
         blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);

         /* iterate once over the upper triangular part of the matrix (saving the corresponding entries in the lower triangular part for the SDPI) */
         for (r = 0; r < blocksize; r++)
         {
            for (c = r; c < blocksize; c++)
            {
               sdpaind = r + blocksize * c;

               if ( REALABS(X[sdpaind]) > sdpisolver->epsilon )
                  startXnblocknonz[b]++;
            }
         }
      }
   }

   /* compute the entry for the LP-block */
   startXnblocknonz[nblocks - 1] = 0;
   sdpablock = (int) sdpisolver->sdpa->getBlockNumber();

   if ( sdpisolver->sdpa->getBlockType(sdpablock) == SDPA::LP )
   {
      X = sdpisolver->sdpa->getResultYMat(sdpablock);
      blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);

      for (i = 0; i < blocksize; i++)
      {
         if ( REALABS(X[i]) > sdpisolver->epsilon )
            startXnblocknonz[b]++;
      }
   }

   return SCIP_OKAY;
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
{
   SCIP_Real* X;
   int b;
   int r;
   int c;
   int sdpaind;
   int sdpablock;
   int blocksize;
   int blocknnonz;

   assert( sdpisolver != NULL );
   assert( nblocks > 0 );
   assert( startXnblocknonz != NULL );
   assert( startXrow != NULL );
   assert( startXcol != NULL );
   assert( startXval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( nblocks != sdpisolver->nsdpblocks + 1 )
   {
      SCIPerrorMessage("SCIPsdpiSolverGetPrimalNonzeros expected nblocks = %d but got %d\n", sdpisolver->nsdpblocks + 1, nblocks);
      return SCIP_LPERROR;
   }

   /* iterate over all SDP-blocks and get the solution */
   for (b = 0; b < sdpisolver->nsdpblocks; b++)
   {
      sdpablock = sdpisolver->inputtoblockmapper[b];

      blocknnonz = 0;

      if ( sdpablock != -1 )
      {
         X = sdpisolver->sdpa->getResultYMat(sdpablock);
         blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);
         blocknnonz = 0;

         /* iterate once over the upper triangular part of the matrix (saving the corresponding entries in the lower triangular part for the SDPI) */
         for (r = 0; r < blocksize; r++)
         {
            for (c = r; c < blocksize; c++)
            {
               sdpaind = r + blocksize * c;

               if ( REALABS(X[sdpaind]) > sdpisolver->epsilon )
               {
                  startXrow[b][blocknnonz] = sdpisolver->blockindmapper[b][c];
                  startXcol[b][blocknnonz] = sdpisolver->blockindmapper[b][r];
                  startXval[b][blocknnonz] = X[sdpaind];
                  blocknnonz++;
                  assert( blocknnonz <= startXnblocknonz[b] );
               }
            }
         }

         startXnblocknonz[b] = blocknnonz;
      }
      else
         startXnblocknonz[b] = 0;
   }

   /* compute entries for the LP-block */
   if ( sdpisolver->sdpa->getBlockType(sdpisolver->nsdpblocks) == SDPA::LP )
   {
      int i;

      X = sdpisolver->sdpa->getResultYMat(sdpisolver->nsdpblocks);
      blocksize = (int) sdpisolver->sdpa->getBlockSize(sdpisolver->nsdpblocks);
      assert( blocksize == sdpisolver->nlpineqs + sdpisolver->nvarbounds );

      /* iterate over LP constraints */
      assert( b == sdpisolver->nsdpblocks );
      blocknnonz = 0;
      for (i = 0; i < sdpisolver->nlpineqs; i++)
      {
         if ( REALABS(X[i]) > sdpisolver->epsilon )
         {
            startXrow[b][blocknnonz] = sdpisolver->rowtoinputmapper[i];
            startXcol[b][blocknnonz] = sdpisolver->rowtoinputmapper[i];
            startXval[b][blocknnonz] = X[i];
            assert( 0 <= startXrow[b][blocknnonz] && startXrow[b][blocknnonz] < 2 * sdpisolver->nlpcons );
            blocknnonz++;
            assert( blocknnonz <= startXnblocknonz[b] );
         }
      }

      /* iterate over varbounds */
      for (i = 0; i < sdpisolver->nvarbounds; i++)
      {
         SCIP_Real val;

         val = X[sdpisolver->nlpineqs + i];

         if ( REALABS(val) > sdpisolver->epsilon )
         {
            int inputpos;
            int vbpos;

            vbpos = sdpisolver->varboundpos[i];
            if ( vbpos > 0 ) /* rhs */
               inputpos = 2 * sdpisolver->nlpcons + 2 * sdpisolver->sdpatoinputmapper[vbpos - 1] + 1;
            else
               inputpos = 2 * sdpisolver->nlpcons + 2 * sdpisolver->sdpatoinputmapper[- vbpos - 1];

            startXrow[b][blocknnonz] = inputpos;
            startXcol[b][blocknnonz] = inputpos;
            startXval[b][blocknnonz] = val;
            assert( 2 * sdpisolver->nlpcons <= inputpos && inputpos < 2 * sdpisolver->nlpcons + 2 * sdpisolver->nvars );
            blocknnonz++;
            assert( blocknnonz <= startXnblocknonz[b] );
         }
      }
      startXnblocknonz[b] = blocknnonz;
   }
   else
      startXnblocknonz[b] = 0;

   return SCIP_OKAY;
}

/** returns the primal solution matrix (without LP rows) */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverGetPrimalSolutionMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nsdpblocks,         /**< number of blocks */
   int*                  sdpblocksizes,      /**< sizes of the blocks */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< pointer to store index change for each block, system is the same as for indchanges */
   SCIP_Real**           primalmatrices      /**< pointer to store values of the primal matrix */
   )
{  /*lint --e{715}*/
   int b;

   assert( sdpisolver != NULL );
   assert( nsdpblocks != sdpisolver->nsdpblocks + 1 );

   /* loop over all SDP blocks */
   for (b = 0; b < nsdpblocks; b++)
   {
      SCIP_Real* X;
      int blocksize;
      int redsize;
      int idx = 0;
      int sdpablock;
      int i;
      int j;

      assert( primalmatrices[b] != NULL );

      blocksize = sdpblocksizes[b];

      /* initialize solution matrix with 0s */
      for (j = 0; j < blocksize * blocksize; ++j)
         primalmatrices[b][j] = 0.0;

      /* treat blocks that were not removed */
      if ( blockindchanges[b] >= 0 )
      {
         sdpablock = sdpisolver->inputtoblockmapper[b];
         assert( sdpablock != -1 );
         X = sdpisolver->sdpa->getResultYMat(sdpablock);
         redsize = (int) sdpisolver->sdpa->getBlockSize(sdpablock);

         /* fill in matrix */
         for (i = 0; i < blocksize; ++i)
         {
            if ( indchanges[b][i] >= 0 )
            {
               assert( 0 <= i - indchanges[b][i] && i - indchanges[b][i] < redsize );
               for (j = 0; j < blocksize; ++j)
               {
                  if ( indchanges[b][j] >= 0 )
                  {
                     assert( 0 <= j - indchanges[b][j] && j - indchanges[b][j] < redsize );
                     assert( 0 <= idx && idx < redsize * redsize );
                     primalmatrices[b][i * blocksize + j] = X[idx++];
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** return the maximum absolute value of the optimal primal matrix */
SCIP_Real SCIPsdpiSolverGetMaxPrimalEntry(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   SCIP_Real* X;
   SCIP_Real maxentry;
   int nblocks;
   int blocksize;
   int arraylength;
   int b;
   int i;

   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );

   nblocks = (int) sdpisolver->sdpa->getBlockNumber();
   maxentry = 0.0;

   /* iterate over all blocks */
   for (b = 0; b < nblocks; b++)
   {
      /* get the length of the X-array, if it is an LP block, the array has length n, otherwise n*(n+1)/2 */
      blocksize = (int) sdpisolver->sdpa->getBlockSize((SDPA_INT) b + 1);
      arraylength = sdpisolver->sdpa->getBlockType((SDPA_INT) b + 1) == SDPA::LP ? blocksize : (blocksize * (blocksize + 1) / 2);

      /* iterate over all entries in this block */
      X = sdpisolver->sdpa->getResultYMat((SDPA_INT) b + 1);/*lint !e747*/
      for (i = 0; i < arraylength; i++)
      {
         if ( X[i] > maxentry )
            maxentry = X[i];
      }
   }

   return maxentry;
}

/** gets the time for the last SDP optimization call of solver */
SCIP_RETCODE SCIPsdpiSolverGetTime(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real*            opttime             /**< pointer to store the time for optimization of the solver */
   )
{
   assert( sdpisolver != NULL );
   assert( opttime != NULL );

   /* not implemented */
   *opttime = 0.0;

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( iterations != NULL );

   *iterations = sdpisolver->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( sdpisolver->sdpa != NULL );
   assert( calls != NULL );

   *calls = sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   *usedsetting = sdpisolver->usedsetting;

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */ /*lint -e{715}*/
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{  /*lint --e{715}*/
   return 1E+20; /* default infinity from SCIP */
}

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{/*lint !e1784*/
   if ( val <= -SCIPsdpiSolverInfinity(sdpisolver) || val >= SCIPsdpiSolverInfinity(sdpisolver) )
      return TRUE;
   return FALSE;
}

/** gets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   assert( sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )/*lint --e{788}*/
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
   case SCIP_SDPPAR_PENALTYPARAM:
      *dval = 0.0;
      SCIPdebugMessage("Parameter SCIP_SDPPAR_PENALTYPARAM not used by SDPA"); /* this parameter is only used by DSDP */
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      *dval = sdpisolver->lambdastar;
      break;
   case SCIP_SDPPAR_WARMSTARTPOGAP:
      *dval = sdpisolver->preoptimalgap;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
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
   case SCIP_SDPPAR_PENALTYPARAM:
      SCIPdebugMessage("Parameter SCIP_SDPPAR_PENALTYPARAM not used by SDPA\n."); /* this parameter is only used by DSDP */
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %g.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIPdebugMessage("Setting sdpisolver lambdastar parameter to %g.\n", dval);
      sdpisolver->lambdastar = dval;
      break;
   case SCIP_SDPPAR_WARMSTARTPOGAP:
      SCIPdebugMessage("Setting sdpisolver preoptgap to %g.\n", dval);
      sdpisolver->preoptimalgap = dval;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** gets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< parameter value */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
      *ival = (int) sdpisolver->sdpinfo;
      SCIPdebugMessage("Getting sdpisolver information output (%d).\n", *ival);
      break;
   case SCIP_SDPPAR_NTHREADS:
      *ival = sdpisolver->nthreads;
      SCIPdebugMessage("Getting sdpisolver number of threads: %d.\n", *ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** sets integer parameter of SDP-Solver */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{/*lint !e1784*/
   assert( sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_SDPINFO:
      sdpisolver->sdpinfo = (SCIP_Bool) ival;
      SCIPdebugMessage("Setting sdpisolver information output (%d).\n", ival);
      break;
   case SCIP_SDPPAR_NTHREADS:
      sdpisolver->nthreads = ival;
      SCIPdebugMessage("Setting sdpisolver number of threads to %d.\n", ival);
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

   return SCIP_OKAY;
}

/** compute and set lambdastar (only used for SDPA) */
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{/*lint !e1784*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );

   /* we set the value to min{max{MIN_LAMBDASTAR, LAMBDASTAR_FACTOR * MAX_GUESS}, MAX_LAMBDASTAR}, where MAX_GUESS is the maximum of the guesses
    * of the SDP-Blocks, if the define LAMBDASTAR_TWOPOINTS is set, we instead choose either LAMBDASTAR_LOW or HIGH depending on LAMBASTAR_THRESHOLD */

#ifdef LAMBDASTAR_TWOPOINTS
   if ( maxguess < LAMBDASTAR_THRESHOLD )
      compval = LAMBDASTAR_LOW;
   else
      compval = LAMBDASTAR_HIGH;
#else
      compval = LAMBDASTAR_FACTOR * maxguess;
#endif

   if ( compval < MIN_LAMBDASTAR )
   {
      sdpisolver->lambdastar = MIN_LAMBDASTAR;
      SCIPdebugMessage("Setting lambdastar to %g.\n", MIN_LAMBDASTAR);
   }
   else if ( compval > MAX_LAMBDASTAR )
   {
      sdpisolver->lambdastar = MAX_LAMBDASTAR;
      SCIPdebugMessage("Setting lambdastar to %g.\n", MAX_LAMBDASTAR);
   }
   else
   {
      sdpisolver->lambdastar = compval;
      SCIPdebugMessage("Setting lambdastar to %g.\n", compval);
   }

   return SCIP_OKAY;
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{/*lint !e1784*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( penaltyparam != NULL );

   compval = PENALTYPARAM_FACTOR * maxcoeff;

   if ( compval < MIN_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penalty parameter to %g.\n", MIN_PENALTYPARAM);
      *penaltyparam = MIN_PENALTYPARAM;
   }
   else if ( compval > MAX_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penalty parameter to %g.\n", MAX_PENALTYPARAM);
      *penaltyparam = MAX_PENALTYPARAM;
   }
   else
   {
      SCIPdebugMessage("Setting penalty parameter to %g.\n", compval);
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
{/*lint !e1784*/
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( maxpenaltyparam != NULL );

   compval = penaltyparam * MAXPENALTYPARAM_FACTOR;

   if ( compval < MAX_MAXPENALTYPARAM )
   {
      *maxpenaltyparam = compval;
      SCIPdebugMessage("Setting maximum penalty parameter to %g.\n", compval);
   }
   else
   {
      *maxpenaltyparam = MAX_MAXPENALTYPARAM;
      SCIPdebugMessage("Setting penalty parameter to %g.\n", MAX_MAXPENALTYPARAM);
   }

   return SCIP_OKAY;
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{/*lint !e1784*/
   assert( fname != NULL );

   sdpisolver->sdpa->writeInputSparse(const_cast<char*>(fname), const_cast<char*>("%8.3f"));

   return SCIP_OKAY;
}

/**@} */
