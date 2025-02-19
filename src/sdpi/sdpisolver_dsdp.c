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

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/

/**@file   sdpisolver_dsdp.c
 * @brief  interface for DSDP
 * @author Tristan Gally
 * @author Marc Pfetsch
 *
 * Note: The memory used to pass data to DSDP needs to be available during the whole solution process and for accessing
 * solutions. We therefore allocate the memory once, reallocated when needed and free it only at the end.
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"

/* turn off warning for DSDSP */
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */
#pragma GCC diagnostic warning "-Wstrict-prototypes"

#ifdef OPENBLAS
#include <cblas.h>
#endif

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */
#include "scip/pub_message.h"                /* for debug and error message */


#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e12     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e4      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */
#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1e-15    /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance */


/** Calls a DSDP-Function and transforms the return-code to a SCIP_LPERROR if needed. */
#define DSDP_CALL(x)  do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as DSDP_CALL, but used for functions returning a boolean. */
#define DSDP_CALL_BOOL(x)  do                                                                                \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as DSDP_CALL, but used for functions returning an int. */
#define DSDP_CALL_INT(x)  do                                                                                 \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return _dsdperrorcode_;                                                          \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced. */
#define DSDP_CALLM(x) do                                                                                     \
                      {                                                                                      \
                         int _dsdperrorcode_;                                                                \
                         if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                         {                                                                                   \
                            SCIPerrorMessage("DSDP-Error <%d> in function call.\n", _dsdperrorcode_);        \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

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

/** Calls a gettimeofday and transforms the return-code to a SCIP_ERROR if needed. */
#define TIMEOFDAY_CALL(x)  do                                                                                \
                      {                                                                                      \
                         int _errorcode_;                                                                    \
                         if ( (_errorcode_ = (x)) != 0 )                                                     \
                         {                                                                                   \
                            SCIPerrorMessage("Error in gettimeofday! \n");                                   \
                            return SCIP_ERROR;                                                               \
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

/** This is the same as CHECK_IF_SOLVED, but will be called for methods returning a bool instead of a SCIP_RETURNCODE. */
#define CHECK_IF_SOLVED_BOOL(sdpisolver)  do                                                                      \
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
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP to add variable bounds to */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in DSDP (nvars minus the number of variables with lb = ub) */
   int*                  inputtodsdpmapper;  /**< entry i gives the index of input variable i in dsdp (starting from 1) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  dsdptoinputmapper;  /**< entry i gives the original index of the (i+1)-th variable in dsdp (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   SCIP_Bool             solved;             /**< Was the SDP solved since the problem was last changed? */
   int                   sdpcounter;         /**< used for debug messages */
   SCIP_Real             epsilon;            /**< tolerance for absolute checks */
   SCIP_Real             gaptol;             /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< feasibility tolerance that should be achieved */
   SCIP_Real             sdpsolverfeastol;   /**< feasibility tolerance given to the SDP-solver */
   SCIP_Real             penaltyparam;       /**< the penalty parameter Gamma used for the penalty formulation if the SDP-solver didn't converge */
   SCIP_Real             objlimit;           /**< objective limit for SDP-solver */
   SCIP_Bool             sdpinfo;            /**< Should the SDP-solver output information to the screen? */
   int                   nthreads;           /**< number of threads the SDP solver should use (-1 = number of cores) */
   SCIP_Bool             penalty;            /**< Did the last solve use a penalty formulation? */
   SCIP_Bool             penaltyworbound;    /**< Was a penalty formulation solved without bounding r? */
   SCIP_Bool             feasorig;           /**< was the last problem solved with a penalty formulation and with original objective coefficents
                                              *   and the solution was feasible for the original problem? */
   SCIP_SDPSOLVERSETTING usedsetting;        /**< setting used to solve the last SDP */
   SCIP_Bool             timelimit;          /**< was the solver stopped because of the time limit? */
   SCIP_Bool             timelimitinitial;   /**< was the problem not even given to the solver because of the time limit? */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   SCIP_Real             opttime;            /**< time spend in optimziation */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
   SCIP_Real*            preoptimalsol;      /**< first feasible solution with gap less or equal preoptimalgap */
   SCIP_Bool             preoptimalsolexists;/**< saved feasible solution with gap less or equal preoptimalgap */
   SCIP_Real             preoptimalgap;      /**< gap at which a preoptimal solution should be saved for warmstarting purposes */

   /* tempory data used for passing data to DSDP */
   int                   dsdpconstsize;      /**< size of dsdpconstind and dsdpconstval */
   int                   dsdpsize;           /**< size of dsdpind and dsdpval */
   int                   dsdplpsize;         /**< size of dsdprow, dsdpcol, dsdplpval */
   int                   dsdplpncols;        /**< size of dsdplpbeg */
   int*                  dsdpconstind;       /**< indices for constant SDP-constraint-matrices */
   SCIP_Real*            dsdpconstval;       /**< non-zero values for constant SDP-constraint-matrices */
   int*                  dsdpind;            /**< indices for SDP-constraint-matrices */
   SCIP_Real*            dsdpval;            /**< non-zero values for SDP-constraint-matrices */
   int*                  dsdplpbeg;          /**< indices at which the columns of the LP constraints start */
   int*                  dsdplprow;          /**< row indices of the columns in the LP part */
   int*                  dsdplpcol;          /**< column indices of the columns (needed for transposing the matrix) */
   SCIP_Real*            dsdplpval;          /**< matrix values of the LP part */
};

/** data used for checking the timelimit */
typedef struct Timings
{
   SDPI_CLOCK*           usedsdpitime;       /**< clock to measure how much time has been used for the current solve */
   SCIP_Real             timelimit;          /**< timelimit in seconds */
   SCIP_Bool             stopped;            /**< was the solver stopped because of the time limit? */
} Timings;


/*
 * Local Functions
 */

/** for given row and column (i,j) computes the position in the lower triangular part
 *  numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
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
/** test if a lower bound lb is not smaller than an upper bound ub, meaning that lb > ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( lb < ub + sdpisolver->feastol );

   return (ub-lb <= sdpisolver->epsilon);
}
#else
#define isFixed(sdpisolver,lb,ub) (ub-lb <= sdpisolver->epsilon)
#endif

/** sort the given row, col and val arrays first by non-decreasing col-indices, than for those with identical col-indices by non-increasing row-indices */
static
void sortColRow(
   int*                  row,                /**< row indices */
   int*                  col,                /**< column indices */
   SCIP_Real*            val,                /**< values */
   int                   length              /**< length of the given arrays */
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

/** check the time limit after each iteration in DSDP */
static
int checkTimeLimitDSDP(
   DSDP                  dsdp,               /**< DSDP-pointer */
   void*                 ctx                 /**< pointer to data of iteration monitor */
   )
{
   Timings* timings;

   assert( dsdp != NULL );
   assert( ctx != NULL );

   timings = (Timings*) ctx;

   if ( SDPIclockGetTime(timings->usedsdpitime) > timings->timelimit )
   {
      DSDP_CALL( DSDPSetConvergenceFlag(dsdp, DSDP_USER_TERMINATION) );/*lint !e641 */
      timings->stopped = TRUE;
      SCIPdebugMessage("Time limit reached! Stopping DSDP.\n");
   }

   return 0;
}

/** check gap and set preoptimal solution if small enough */
static
int checkGapSetPreoptimalSol(
   DSDP                  dsdp,               /**< DSDP-pointer */
   void*                 ctx                 /**< pointer to data of iteration monitor */
   )
{
   SCIP_Real absgap;
   SCIP_Real relgap;
   SCIP_Real pobj;
   SCIP_Real dobj;
   SCIP_Real r;

   /* we only need to set the preoptimal solution once and we do not save it if the penalty formulation was used (in that case we won't warmstart
    * anyways and, most importantly, if solving without bound on r, we added another variable, so the memory would not be enough)
    */
   if ( ((SCIP_SDPISOLVER*) ctx)->preoptimalsolexists || ((SCIP_SDPISOLVER*) ctx)->penalty || ((SCIP_SDPISOLVER*) ctx)->penaltyworbound )
      return 0;

   DSDP_CALL_INT( DSDPGetPPObjective(dsdp,&pobj) );
   DSDP_CALL_INT( DSDPGetDDObjective(dsdp,&dobj) );
   DSDP_CALL_INT( DSDPGetDualityGap(dsdp,&absgap) );

   relgap = absgap / (1.0 + (REALABS(dobj)/2) + (REALABS(pobj)/2) ); /* compare dsdpconverge.c */

   /* check feasibility through penalty variable r */
   DSDP_CALL_INT( DSDPGetR(dsdp,&r) );

   if ( r < ((SCIP_SDPISOLVER*) ctx)->feastol && relgap < ((SCIP_SDPISOLVER*) ctx)->preoptimalgap )
   {
      DSDP_CALL_INT( DSDPGetY(dsdp, ((SCIP_SDPISOLVER*) ctx)->preoptimalsol, ((SCIP_SDPISOLVER*) ctx)->nactivevars) );
      ((SCIP_SDPISOLVER*) ctx)->preoptimalsolexists = TRUE;
      SCIPdebugMessage("penalty variable %f, gap %f -> saving preoptimal solution\n", r, relgap);
   }

   return 0;
}


/** ensure size of the SDP data */
static
SCIP_RETCODE ensureSDPDataSize(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   sdpconstnnonz,      /**< required space for the constant matrix */
   int                   sdpnnonz            /**< required space for the matrices */
   )
{
   if ( sdpconstnnonz > sdpisolver->dsdpconstsize )
   {
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpconstind, sdpisolver->dsdpconstsize, sdpconstnnonz) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpconstval, sdpisolver->dsdpconstsize, sdpconstnnonz) );
      sdpisolver->dsdpconstsize = sdpconstnnonz;
   }

   if ( sdpnnonz > sdpisolver->dsdpsize )
   {
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpind, sdpisolver->dsdpsize, sdpnnonz) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpval, sdpisolver->dsdpsize, sdpnnonz) );
      sdpisolver->dsdpsize = sdpnnonz;
   }

   return SCIP_OKAY;
}

/** ensure size of the LP data */
static
SCIP_RETCODE ensureLPDataSize(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   ncols,              /**< required number of columns */
   int                   lpnnonz             /**< required space for LP constraints */
   )
{
   if ( ncols > sdpisolver->dsdplpncols )
   {
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpbeg, sdpisolver->dsdplpncols, ncols) );
      sdpisolver->dsdplpncols = ncols;
   }

   if ( lpnnonz > sdpisolver->dsdplpsize )
   {
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplprow, sdpisolver->dsdplpsize, lpnnonz) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpcol, sdpisolver->dsdplpsize, lpnnonz) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpval, sdpisolver->dsdplpsize, lpnnonz) );
      sdpisolver->dsdplpsize = lpnnonz;
   }

   return SCIP_OKAY;
}


/** free SDP and LP data */
static
SCIP_RETCODE freeDataSize(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   if ( sdpisolver->dsdpconstsize > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpconstval, sdpisolver->dsdpconstsize);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpconstind, sdpisolver->dsdpconstsize);
      sdpisolver->dsdpconstsize = 0;
   }

   if ( sdpisolver->dsdpsize > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpval, sdpisolver->dsdpsize);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdpind, sdpisolver->dsdpsize);
      sdpisolver->dsdpsize = 0;
   }

   if ( sdpisolver->dsdplpncols > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpbeg, sdpisolver->dsdplpncols);
      sdpisolver->dsdplpncols = 0;
   }

   if ( sdpisolver->dsdplpsize > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpval, sdpisolver->dsdplpsize);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplprow, sdpisolver->dsdplpsize);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sdpisolver->dsdplpcol, sdpisolver->dsdplpsize);
      sdpisolver->dsdplpsize = 0;
   }

   return SCIP_OKAY;
}

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name and version (if available) of SDP-solver*/
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "DSDP"; /* getting the version is not supported in DSDP */
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Dual-Scaling Interior Point SDP-Solver by S. Benson, Y. Ye, and X. Zhang (http://www.mcs.anl.gov/hs/software/DSDP/)";
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
{
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->dsdp;
}

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   )
{
   return 10;
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
{
   assert( sdpisolver != NULL );
   assert( blkmem != NULL );
   assert( bufmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

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
   (*sdpisolver)->objcoefs = NULL;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->timelimitinitial = FALSE;
   (*sdpisolver)->penalty = FALSE;
   (*sdpisolver)->penaltyworbound = FALSE;
   (*sdpisolver)->feasorig = FALSE;
   (*sdpisolver)->sdpcounter = 0;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->opttime = 0.0;
   (*sdpisolver)->nsdpcalls = 0;

   (*sdpisolver)->epsilon = 1e-9;
   (*sdpisolver)->gaptol = 1e-6;
   (*sdpisolver)->feastol = 1e-6;
   (*sdpisolver)->sdpsolverfeastol = 1e-6;
   (*sdpisolver)->penaltyparam = 1e5;
   (*sdpisolver)->objlimit = SCIPsdpiSolverInfinity(*sdpisolver);
   (*sdpisolver)->sdpinfo = FALSE;
   (*sdpisolver)->nthreads = -1;
   (*sdpisolver)->usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   (*sdpisolver)->preoptimalsolexists = FALSE;
   (*sdpisolver)->preoptimalgap = -1.0;

   (*sdpisolver)->dsdpconstsize = 0;
   (*sdpisolver)->dsdpsize = 0;
   (*sdpisolver)->dsdplpsize = 0;
   (*sdpisolver)->dsdplpncols = 0;
   (*sdpisolver)->dsdpconstind = NULL;
   (*sdpisolver)->dsdpconstval = NULL;
   (*sdpisolver)->dsdpind = NULL;
   (*sdpisolver)->dsdpval = NULL;
   (*sdpisolver)->dsdplpbeg = NULL;
   (*sdpisolver)->dsdplprow = NULL;
   (*sdpisolver)->dsdplpcol = NULL;
   (*sdpisolver)->dsdplpval = NULL;

   return SCIP_OKAY;
}

/** deletes an SDP solver interface */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );
   assert( *sdpisolver != NULL );

   SCIPdebugMessage("Freeing SDPISolver\n");

   SCIP_CALL( freeDataSize(*sdpisolver) );

   if ( (*sdpisolver)->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy((*sdpisolver)->dsdp) );
   }

   if ( (*sdpisolver)->nvars > 0 )
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->inputtodsdpmapper, (*sdpisolver)->nvars);

   if ( (*sdpisolver)->nactivevars > 0 )
   {
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->preoptimalsol, (*sdpisolver)->nactivevars);
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->dsdptoinputmapper, (*sdpisolver)->nactivevars);
      BMSfreeBlockMemoryArray((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->nactivevars);
   }

   if ( (*sdpisolver)->nvars >= (*sdpisolver)->nactivevars )
      BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->nvars - (*sdpisolver)->nactivevars);

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/** increases the SDP-Counter */
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );

   sdpisolver->sdpcounter++;

   return SCIP_OKAY;
}

/** reset the SDP-Counter to zero */
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
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
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independant of some lhs/rhs being infinity (the starting point
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
   return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars,
      sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
      indchanges, nremovedinds, blockindchanges, nremovedblocks, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind,
      lpval, starty, startZnblocknonz, startZrow, startZcol, startZval, startXnblocknonz, startXrow, startXcol, startXval, startsettings,
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
 */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< should r be non-negative ? */
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
{  /*lint --e{413,715}*/
   int maxlpnnonz = 0;
   int i;
   int j;
   int ind;
   int block;
   int startind;
   int nfixedvars;
   int dsdpnlpnonz = 0;
   int nrnonz = 0;
   int oldnactivevars;
   SCIP_Real feastol;
   SCIP_Real gaptol;
   SCIP_Real solvertimelimit;
   SCIP_Real oldsdpitime;
   Timings timings;

#ifdef SCIP_DEBUG
   DSDPTerminationReason reason; /* this will later be used to check if DSDP converged */
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

   if ( solvertimelimit <= 0.0 )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->timelimitinitial = TRUE;
      sdpisolver->solved = FALSE;
      return SCIP_OKAY;
   }
   else
      sdpisolver->timelimitinitial = FALSE;

   sdpisolver->feasorig = FALSE;
   sdpisolver->penalty = penaltyparam > sdpisolver->epsilon;

   /* start the timing */
   timings.usedsdpitime = usedsdpitime;
   timings.timelimit = timelimit;
   timings.stopped = FALSE;

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP), also remember settings for statistics */
   if ( penaltyparam < sdpisolver->epsilon )
   {
      SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", ++sdpisolver->sdpcounter);
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_FAST;
   }
   else
   {
      SCIPdebugMessage("Inserting Data again into DSDP for SDP (%d) \n", sdpisolver->sdpcounter);
      sdpisolver->usedsetting = SCIP_SDPSOLVERSETTING_PENALTY;
   }

   /* allocate memory for inputtomosekmapper, mosektoinputmapper and the fixed and active variable information, for the latter this will
    * later be shrinked if the needed size is known */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtodsdpmapper), sdpisolver->nvars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), sdpisolver->nactivevars, nvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), sdpisolver->nvars - sdpisolver->nactivevars, nvars) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), sdpisolver->nactivevars, nvars) ); /*lint !e776*/

   sdpisolver->nvars = nvars;
   oldnactivevars = sdpisolver->nactivevars;
   sdpisolver->nactivevars = 0;
   nfixedvars = 0;
   sdpisolver->niterations = 0;
   sdpisolver->nsdpcalls = 0;

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
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];
         sdpisolver->nactivevars++;
         sdpisolver->inputtodsdpmapper[i] = sdpisolver->nactivevars; /* dsdp starts counting at 1, so we do this after increasing nactivevars */
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in DSDP\n", i, sdpisolver->inputtodsdpmapper[i], sdpisolver->sdpcounter);
#endif
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      SCIPdebugMessage("Variable %d is the slack variable for the explicit penalty formulation.\n", sdpisolver->nactivevars + 1);
   }

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* shrink the fixedvars, objcoefs and mosektoinputmapper arrays to the right size */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->objcoefs), nvars, sdpisolver->nactivevars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->fixedvarsval), nvars, nfixedvars) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->dsdptoinputmapper), nvars, sdpisolver->nactivevars) );

   /* adjust length of preoptimal solution array */
   if ( sdpisolver->nactivevars != oldnactivevars )
   {
      if ( oldnactivevars == 0 )
      {
         BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->preoptimalsol), sdpisolver->nactivevars) );
      }
      else
      {
         BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->preoptimalsol), oldnactivevars, sdpisolver->nactivevars) );
      }
   }
   sdpisolver->preoptimalsolexists = FALSE;

   /* insert data */

   if ( sdpisolver->dsdp != NULL )
   {
      DSDP_CALL( DSDPDestroy(sdpisolver->dsdp) ); /* if there already exists a DSDP-instance, destroy the old one */
   }

   /* in case we don't want to bound r, we can't use the penalty formulation in DSDP and have to give r explicitly */
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      DSDP_CALLM( DSDPCreate(sdpisolver->nactivevars + 1, &(sdpisolver->dsdp)) );
      sdpisolver->penaltyworbound = TRUE;
   }
   else
   {
      DSDP_CALLM( DSDPCreate(sdpisolver->nactivevars, &(sdpisolver->dsdp)) );
      sdpisolver->penaltyworbound = FALSE;
   }
   DSDP_CALLM( DSDPCreateSDPCone(sdpisolver->dsdp, nsdpblocks - nremovedblocks, &(sdpisolver->sdpcone)) );
   if ( nlpcons > 0 )
   {
      DSDP_CALLM( DSDPCreateLPCone(sdpisolver->dsdp, &(sdpisolver->lpcone)) );
   }
   DSDP_CALLM( DSDPCreateBCone(sdpisolver->dsdp, &(sdpisolver->bcone)) );

   /* determine sizes */
   nrnonz = 0;
   if ( penaltyparam > sdpisolver->epsilon && ! rbound )
   {
      /* we need to compute the total number of nonzeros for the slack variable r, which equals the total number of diagonal entries */
      for (block = 0; block < nsdpblocks; block++)
         nrnonz += sdpblocksizes[block] - nremovedinds[block];
      assert( nrnonz >= 0 );
   }

   SCIP_CALL( ensureSDPDataSize(sdpisolver, sdpconstnnonz, sdpnnonz + nrnonz) );

#ifdef SCIP_MORE_DEBUG
   SCIPmessagePrintInfo(sdpisolver->messagehdlr, "setting objective values for SDP %d:\n", sdpisolver->sdpcounter);
#endif

   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withobj )
      {
         /* insert objective value, DSDP counts from 1 to n instead of 0 to n-1, *(-1) because DSDP maximizes instead of minimizing */
         DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, i+1, -1.0 * obj[sdpisolver->dsdptoinputmapper[i]]) );
#ifdef SCIP_MORE_DEBUG
         SCIPmessagePrintInfo(sdpisolver->messagehdlr, "var %d (was var %d): %f, ", i+1, sdpisolver->dsdptoinputmapper[i], obj[sdpisolver->dsdptoinputmapper[i]]);
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

   /* insert the objective value for r if solving without rbound, it is variable nactivevars + 1 and the objective is multiplied by -1 as we maximize */
   if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
   {
      DSDP_CALL( DSDPSetDualObjective(sdpisolver->dsdp, sdpisolver->nactivevars + 1, -1.0 * penaltyparam) );
#ifdef SCIP_MORE_DEBUG
      SCIPmessagePrintInfo(sdpisolver->messagehdlr, "slack variable r: %f, ", penaltyparam);
#endif
   }

#ifdef SCIP_MORE_DEBUG
   SCIPmessagePrintInfo(sdpisolver->messagehdlr, "\n");
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
               if ( v == sdpvar[block][k] )
               {
                  blockvar = k;
                  break;
               }
            }

            startind = ind;

            if ( blockvar > -1 ) /* the variable exists in this block */
            {
               for (k = 0; k < sdpnblockvarnonz[block][blockvar]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( indchanges[block][sdprow[block][blockvar][k]] > -1 && indchanges[block][sdpcol[block][blockvar][k]] > -1 );

                  /* substract the number of removed indices before the row and col to get the indices after fixings */
                  sdpisolver->dsdpind[ind] = compLowerTriangPos(sdprow[block][blockvar][k] - indchanges[block][sdprow[block][blockvar][k]],
                     sdpcol[block][blockvar][k] - indchanges[block][sdpcol[block][blockvar][k]]);
                  sdpisolver->dsdpval[ind] = - sdpval[block][blockvar][k];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
                  ind++;
               }

               /* sort the arrays for this matrix (by non decreasing indices) as this might help the solving time of DSDP */
               SCIPsortIntReal(sdpisolver->dsdpind + startind, sdpisolver->dsdpval + startind, sdpnblockvarnonz[block][blockvar]);

               assert( blockindchanges[block] > -1 ); /* we shouldn't insert into blocks we removed */

               /* i + 1 because DSDP starts counting the variables at 1, adding startind shifts the arrays to the first
                * nonzero belonging to this block and this variable */
               DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], i + 1, sdpblocksizes[block] - nremovedinds[block],
                     1.0, 0, sdpisolver->dsdpind + startind, sdpisolver->dsdpval + startind, sdpnblockvarnonz[block][blockvar]));
            }
         }
      }

      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
      {
         startind = ind;
         /* add r * Identity for each block */
         for (block = 0; block < nsdpblocks; block++)
         {
            if ( blockindchanges[block] > -1 )
            {
               for (i = 0; i < sdpblocksizes[block] - nremovedinds[block]; i++)
               {
                  sdpisolver->dsdpind[ind] = compLowerTriangPos(i, i);
                  sdpisolver->dsdpval[ind] = -1.0; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0 + r*I) should be positive semidefinite */
                  ind++;
               }
               DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], sdpisolver->nactivevars + 1,
                     sdpblocksizes[block] - nremovedinds[block], 1.0, 0, sdpisolver->dsdpind + ind - (sdpblocksizes[block] - nremovedinds[block]) ,
                     sdpisolver->dsdpval + ind - (sdpblocksizes[block] - nremovedinds[block]), sdpblocksizes[block] - nremovedinds[block]) ); /*lint !e679*/
            }
         }
         assert( ind - startind == nrnonz );
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
               sdpisolver->dsdpconstind[ind] = compLowerTriangPos(sdpconstrow[block][i] - indchanges[block][sdpconstrow[block][i]],
                  sdpconstcol[block][i] - indchanges[block][sdpconstcol[block][i]]);
               sdpisolver->dsdpconstval[ind] = - sdpconstval[block][i]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j) should be positive semidefinite */
               ind++;
            }

            /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
            SCIPsortIntReal(sdpisolver->dsdpconstind + startind, sdpisolver->dsdpconstval + startind, sdpconstnblocknonz[block]);

            assert( blockindchanges[block] > -1 ); /* we shouldn't insert into a block we removed */

            /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block by adding
             * startind, ind - startind gives the number of elements for this block */
            DSDP_CALL( SDPConeSetASparseVecMat(sdpisolver->sdpcone, block - blockindchanges[block], 0, sdpblocksizes[block] - nremovedinds[block],
                  1.0, 0, sdpisolver->dsdpconstind + startind, sdpisolver->dsdpconstval + startind, ind - startind));
         }
      }
   }

#ifdef SCIP_MORE_DEBUG
   SDPConeView2(sdpisolver->sdpcone);
#endif

   /* start inserting the LP constraints */
   if ( nlpcons > 0 || lpnnonz > 0 || ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
   {
      int pos;
      int lhsrhscnt = 0;
      int colidx;
      int nlpineqs = 0;

      assert( lpbeg != NULL );
      assert( lpind != NULL );
      assert( lpval != NULL );

      /* compute the number of inequalities and nonzeros (consider splitting of ranged rows) */
      for (i = 0; i < nlpcons; i++)
      {
         int nextbeg;

         if ( lpindchanges[i] < 0 )
            continue;

         if ( i == nlpcons - 1 )
            nextbeg = lpnnonz;
         else
            nextbeg = lpbeg[i + 1];

         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            ++nlpineqs;
            maxlpnnonz += nextbeg - lpbeg[i] + 1;  /* +1 for lhs */
         }

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            ++nlpineqs;
            maxlpnnonz += nextbeg - lpbeg[i] + 1;  /* +1 for rhs */
         }
      }
      assert( nlpineqs <= 2 * nlpcons );
      assert( maxlpnnonz <= 2 * lpnnonz + 2 * nlpcons );

      /* determine upper bound on number of nonnzeros */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
         maxlpnnonz += sdpisolver->nactivevars + 1;   /* + 1 for objective limit value */

      if ( penaltyparam > sdpisolver->epsilon && (! rbound) )
         maxlpnnonz += nlpineqs;

      SCIP_CALL( ensureLPDataSize(sdpisolver, sdpisolver->nactivevars + 3, maxlpnnonz) );

      /* iterate over LP rows */
      for (i = 0; i < nlpcons; ++i)
      {
         SCIP_Bool lhspresent = FALSE;
         SCIP_Bool rhspresent = FALSE;
         int nextbeg;
         int v;

         if ( lpindchanges[i] < 0 )
            continue;

         /* determine presence of new lhs/rhs */
         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
            lhspresent = TRUE;

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
            rhspresent = TRUE;

         assert( 0 <= lpbeg[i] && lpbeg[i] < lpnnonz );
         assert( lhspresent || rhspresent );

         if ( i == nlpcons - 1 )
            nextbeg = lpnnonz;
         else
            nextbeg = lpbeg[i + 1];

         for (j = lpbeg[i]; j < nextbeg; ++j)
         {
            assert( 0 <= lpind[j] && lpind[j] < nvars );
            assert( REALABS(lpval[j]) > sdpisolver->epsilon );

            v = sdpisolver->inputtodsdpmapper[lpind[j]];
            if ( v >= 0 )
            {
               assert( ! isFixed(sdpisolver, lb[lpind[j]], ub[lpind[j]]) );
               assert( 1 <= v && v <= sdpisolver->nactivevars );

               if ( lhspresent )
               {
                  sdpisolver->dsdplpcol[dsdpnlpnonz] = v - 1;   /* minus 1 because inputtodsdpmapper starts at 1 */
                  sdpisolver->dsdplprow[dsdpnlpnonz] = lhsrhscnt;
                  sdpisolver->dsdplpval[dsdpnlpnonz] = -lpval[j]; /* - because dsdp wants <= instead of >= constraints */
                  ++dsdpnlpnonz;
               }

               if ( rhspresent )
               {
                  sdpisolver->dsdplpcol[dsdpnlpnonz] = v - 1;   /* minus 1 because inputtodsdpmapper starts at 1 */
                  sdpisolver->dsdplprow[dsdpnlpnonz] = lhsrhscnt + lhspresent;
                  sdpisolver->dsdplpval[dsdpnlpnonz] = lpval[j];
                  ++dsdpnlpnonz;
               }
            }
         }

         /* move counter forward */
         if ( lhspresent )
            ++lhsrhscnt;
         if ( rhspresent )
            ++lhsrhscnt;
      }
      assert( lhsrhscnt == nlpineqs );
      assert( dsdpnlpnonz <= maxlpnnonz );

      /* possibly add row for objective cutoff */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         for (j = 0; j < nvars; ++j)
         {
            int v;

            v = sdpisolver->inputtodsdpmapper[j];
            if ( v >= 0 )
            {
               assert( ! isFixed(sdpisolver, lb[j], ub[j]) );
               assert( 1 <= v && v <= sdpisolver->nactivevars );

               if ( REALABS(obj[j]) > sdpisolver->epsilon )
               {
                  sdpisolver->dsdplpcol[dsdpnlpnonz] = v - 1;   /* minus 1 because inputtodsdpmapper starts at 1 */
                  sdpisolver->dsdplprow[dsdpnlpnonz] = nlpineqs;
                  sdpisolver->dsdplpval[dsdpnlpnonz] = obj[j];
                  ++dsdpnlpnonz;
               }
            }
         }
      }
      assert( dsdpnlpnonz <= maxlpnnonz );

      /* transpose matrix by sorting first by column and then by rows */
      sortColRow(sdpisolver->dsdplprow, sdpisolver->dsdplpcol, sdpisolver->dsdplpval, dsdpnlpnonz);

      /* determine next column index */
      colidx = sdpisolver->nactivevars;

      /* add r * Identity if using a penalty formulation without a bound on r */
      if ( penaltyparam > sdpisolver->epsilon && ! rbound )
      {
         for (i = 0; i < nlpineqs; i++)
         {
            sdpisolver->dsdplpcol[dsdpnlpnonz] = colidx;
            sdpisolver->dsdplprow[dsdpnlpnonz] = i;
            sdpisolver->dsdplpval[dsdpnlpnonz] = -1.0; /* for >=-inequalities we would add a +1, but then we have to multiply these with -1 for DSDP */
            ++dsdpnlpnonz;
         }
         ++colidx;
      }
      assert( dsdpnlpnonz <= maxlpnnonz );

      /* add all left- and right-hand-sides as last column */
      pos = 0;
      for (i = 0; i < nlpcons; i++)
      {
         if ( lpindchanges[i] < 0 )
            continue;

         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            if ( REALABS(lplhs[i]) > sdpisolver->epsilon )
            {
               sdpisolver->dsdplpcol[dsdpnlpnonz] = colidx;
               sdpisolver->dsdplprow[dsdpnlpnonz] = pos;
               sdpisolver->dsdplpval[dsdpnlpnonz] = -lplhs[i]; /* we multiply by -1 because DSDP wants <= instead of >= */
               ++dsdpnlpnonz;
            }
            ++pos;
         }

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            if ( REALABS(lprhs[i]) > sdpisolver->epsilon )
            {
               sdpisolver->dsdplpcol[dsdpnlpnonz] = colidx;
               sdpisolver->dsdplprow[dsdpnlpnonz] = pos;
               sdpisolver->dsdplpval[dsdpnlpnonz] = lprhs[i];
               ++dsdpnlpnonz;
            }
            ++pos;
         }
      }
      assert( pos == nlpineqs );
      assert( dsdpnlpnonz <= maxlpnnonz );

      /* add the right-hand-side for the objective bound */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         if ( REALABS(sdpisolver->objlimit) > sdpisolver->epsilon )
         {
            sdpisolver->dsdplpcol[dsdpnlpnonz] = colidx;
            sdpisolver->dsdplprow[dsdpnlpnonz] = nlpineqs;
            sdpisolver->dsdplpval[dsdpnlpnonz] = sdpisolver->objlimit; /* as we want <= upper bound, this is the correct type of inequality for DSDP */
            dsdpnlpnonz++;
         }
      }
      assert( dsdpnlpnonz <= maxlpnnonz );
      ++colidx;

      /* set up begcol array */
      pos = 0;
      for (j = 0; j < colidx; ++j)
      {
         while ( pos < dsdpnlpnonz && j > sdpisolver->dsdplpcol[pos] )
            ++pos;
         sdpisolver->dsdplpbeg[j] = pos;
      }
      sdpisolver->dsdplpbeg[colidx] = dsdpnlpnonz;
      assert( colidx <= sdpisolver->nactivevars + 3 );

#ifndef NDEBUG
      for (j = 0; j < colidx; ++j)
      {
         int l;

         for (l = sdpisolver->dsdplpbeg[j]; l < sdpisolver->dsdplpbeg[j+1]; ++l)
         {
            assert( sdpisolver->dsdplpcol[l] == j );
            assert( 0 <= sdpisolver->dsdplprow[l] && sdpisolver->dsdplprow[l] < nlpineqs + 1 );
            assert( REALABS(sdpisolver->dsdplpval[l]) > sdpisolver->epsilon );
         }
      }
#endif

      /* add the arrays to dsdp (in this case we need no additional if for the penalty version without bounds, as we already added the extra var,
       * so DSDP knows, that there is an additional entry in dsdplpbeg which then gives the higher number of nonzeros) */
      if ( SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         DSDP_CALL( LPConeSetData2(sdpisolver->lpcone, nlpineqs, sdpisolver->dsdplpbeg, sdpisolver->dsdplprow, sdpisolver->dsdplpval) );
      }
      else
      {
         DSDP_CALL( LPConeSetData2(sdpisolver->lpcone, nlpineqs + 1, sdpisolver->dsdplpbeg, sdpisolver->dsdplprow, sdpisolver->dsdplpval) );
      }
#ifdef SCIP_MORE_DEBUG
      LPConeView2(sdpisolver->lpcone);
#endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpisolver->sdpcounter);

   DSDP_CALL( DSDPSetGapTolerance(sdpisolver->dsdp, sdpisolver->gaptol) );  /* set DSDP's tolerance for duality gap */
   DSDP_CALL( DSDPSetRTolerance(sdpisolver->dsdp, sdpisolver->sdpsolverfeastol) );    /* set DSDP's tolerance for the SDP-constraints */
#ifndef SCIP_DEBUG
   if ( sdpisolver->sdpinfo )
#endif
   {
      DSDP_CALL( DSDPSetStandardMonitor(sdpisolver->dsdp, 1) );   /* output DSDP information after every 1 iteration */
   }

   /* set the penalty parameter (only if rbound = TRUE, otherwise we had to add everything ourselves) */
   if ( penaltyparam >= sdpisolver->epsilon && rbound ) /* in sdpisolverSolve this is called with an exact 0 */
   {
      DSDP_CALL( DSDPSetPenaltyParameter(sdpisolver->dsdp, penaltyparam) );
      DSDP_CALL( DSDPUsePenalty(sdpisolver->dsdp, 1) );
   }
   else
   {
      /* set the penalty parameter to the default value */
      DSDP_CALL( DSDPSetPenaltyParameter(sdpisolver->dsdp, sdpisolver->penaltyparam) );
   }

   /* set the starting solution */
   if ( starty != NULL )
   {
      for (i = 0; i < sdpisolver->nactivevars; i++) /* we iterate over the variables in DSDP */
      {
         DSDP_CALL( DSDPSetY0(sdpisolver->dsdp, i + 1, starty[sdpisolver->dsdptoinputmapper[i]]) ); /* i+1 since DSDP uses indices 1 to n */
      }
   }

#ifdef OPENBLAS
   if ( sdpisolver->nthreads > 0 )
      openblas_set_num_threads(sdpisolver->nthreads);
#endif

   /* start the solving process */
   DSDP_CALLM( DSDPSetup(sdpisolver->dsdp) );

   /* if there is a timelimit, set the corresponding callback */
   if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, timelimit) )
   {
      DSDP_CALLM( DSDPSetMonitor(sdpisolver->dsdp, checkTimeLimitDSDP, (void*) &timings) );
   }

   /* if preoptimal solutions should be saved for warmstarting purposes, set the corresponding callback */
   if ( sdpisolver->preoptimalgap >= 0.0 )
   {
      DSDP_CALL( DSDPSetMonitor(sdpisolver->dsdp, checkGapSetPreoptimalSol, (void*) sdpisolver) );
   }
   oldsdpitime = SDPIclockGetTime(usedsdpitime);
   DSDP_CALL( DSDPSolve(sdpisolver->dsdp) );
   sdpisolver->opttime = SDPIclockGetTime(usedsdpitime) - oldsdpitime;

   DSDP_CALL( DSDPGetIts(sdpisolver->dsdp, &(sdpisolver->niterations)) );
   sdpisolver->nsdpcalls++;

   /* check if solving was stopped because of the time limit */
   if ( timings.stopped )
   {
      sdpisolver->timelimit = TRUE;
      sdpisolver->solved = FALSE;
   }
   else
   {
      sdpisolver->timelimit = FALSE;
      DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) ); /* computes X and determines feasibility and unboundedness of the solution */
      sdpisolver->solved = TRUE;
   }

   /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
    * reports feasibility, resolve it with adjusted tolerance */
   feastol = sdpisolver->sdpsolverfeastol;
   gaptol = sdpisolver->gaptol;

   while ( SCIPsdpiSolverIsAcceptable(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver) && penaltyparam < sdpisolver->epsilon && feastol >= INFEASMINFEASTOL )
   {
      SCIP_Real* solvector;
      SCIP_Real primalobj;
      SCIP_Real dualobj;
      SCIP_Bool infeasible;
      SCIP_Bool solveagain = FALSE;
      int newiterations;

      /* get current solution */
      BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solvector, nvars) );
      SCIP_CALL( SCIPsdpiSolverGetDualSol(sdpisolver, NULL, solvector) );

      /* check the solution for feasibility with regards to our tolerance */
      SCIP_CALL( SCIPsdpSolcheckerCheck(sdpisolver->bufmem, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
            indchanges, nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval,
            solvector, sdpisolver->feastol, sdpisolver->epsilon, &infeasible) );
      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solvector);

      if ( infeasible )
      {
         feastol *= INFEASFEASTOLCHANGE;
         if ( feastol >= INFEASMINFEASTOL )
         {
            SCIPdebugMessage("Solution feasible for DSDP but outside feasibility tolerance, changing DSDP feasibility tolerance to %g.\n", feastol);
            DSDP_CALL( DSDPSetRTolerance(sdpisolver->dsdp, feastol) );    /* set DSDP's tolerance for the SDP-constraints */
            solveagain = TRUE;
         }
      }

      /* Check whether duality gap is small enough: We use (PP) and (DD), since these are the problems that are actually solved. */
      DSDP_CALL( DSDPGetDDObjective(sdpisolver->dsdp, &dualobj) );
      DSDP_CALL( DSDPGetPPObjective(sdpisolver->dsdp, &primalobj) );
      if ( REALABS(dualobj - primalobj) >= sdpisolver->gaptol )
      {
         infeasible = TRUE;
         gaptol *= INFEASFEASTOLCHANGE;
         if ( gaptol >= INFEASMINFEASTOL )
         {
            SCIPdebugMessage("Solution feasible, but duality gap %g is too large, changing DSDP gap tolerance to %g.\n", REALABS(dualobj - primalobj), gaptol);
            DSDP_CALL( DSDPSetGapTolerance(sdpisolver->dsdp, gaptol) );  /* set DSDP's tolerance for duality gap */
            solveagain = TRUE;
         }
      }

      if ( solveagain )
      {
         oldsdpitime = SDPIclockGetTime(usedsdpitime);
         DSDP_CALL( DSDPSolve(sdpisolver->dsdp) );
         sdpisolver->opttime += SDPIclockGetTime(usedsdpitime) - oldsdpitime;

         /* update number of SDP-iterations and -calls */
         sdpisolver->nsdpcalls++;
         DSDP_CALL( DSDPGetIts(sdpisolver->dsdp, &newiterations) );
         sdpisolver->niterations += newiterations;

         /* check if solving was stopped because of the time limit */
         if ( timings.stopped )
         {
            sdpisolver->timelimit = TRUE;
            sdpisolver->solved = FALSE;
         }
         else
         {
            sdpisolver->timelimit = FALSE;
            DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) ); /* computes X and determines feasibility and unboundedness of the solution */
            sdpisolver->solved = TRUE;
         }
      }
      else
      {
         if ( infeasible )
         {
            sdpisolver->solved = FALSE;
            SCIPmessagePrintInfo(sdpisolver->messagehdlr, "DSDP failed to reach required feasibility tolerance (feastol: %g, gaptol: %g)!\n", feastol, gaptol);
         }
         break;
      }
   }

#ifdef SCIP_DEBUG
   DSDP_CALL( DSDPStopReason(sdpisolver->dsdp, &reason) );

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

   if ( penaltyparam >= sdpisolver->epsilon && sdpisolver->solved )
   {
      if ( rbound )
      {
         /* in this case we used the penalty-formulation of DSDP, so we can check their value of r */
         SCIP_Real rval;
         SCIP_Real trace;

         DSDP_CALL( DSDPGetR(sdpisolver->dsdp, &rval) );

         *feasorig = (rval < sdpisolver->feastol );

         /* only set sdpisolver->feasorig to true if we solved with objective, because only in this case we want to compute
          * the objective value by hand since it is numerically more stable then the result returned by DSDP */
         if ( withobj )
            sdpisolver->feasorig = *feasorig;

         /* if r > 0 or we are in debug mode, also check the primal bound */
         if ( ! *feasorig )
         {
            if ( penaltybound != NULL )
            {
               SCIPdebugMessage("Solution not feasible in original problem, r = %g.\n", rval);

               /* get the trace of X to compare it with the penalty parameter */
               DSDP_CALL( DSDPGetTraceX(sdpisolver->dsdp, &trace) );

#if 0 /* DSDP doesn't seem to adhere to its own feasiblity tolerance */
               assert( trace < penaltyparam + sdpisolver->feastol ); /* solution should be primal feasible */
#endif

               /* if the relative gap is smaller than the tolerance, we return equality */
               if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )
               {
                  *penaltybound = TRUE;
                  SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible\n",
                     trace, penaltyparam);
               }
               else
                  *penaltybound = FALSE;
            }
         }
      }
      else
      {
         SCIP_Real* dsdpsol;
         SCIP_Real trace;

         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &dsdpsol, sdpisolver->nactivevars + 1) ); /*lint !e776*/
         /* last entry of DSDPGetY needs to be the number of variables, will return an error otherwise */
         DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, sdpisolver->nactivevars + 1) );

         *feasorig = (dsdpsol[sdpisolver->nactivevars] < sdpisolver->feastol); /* r is the last variable in DSDP, so the last entry gives us the value */
         if ( ! *feasorig )
         {
            if ( penaltybound != NULL )
            {
               SCIPdebugMessage("Solution not feasible in original problem, r = %g.\n", dsdpsol[sdpisolver->nactivevars]);

               /* get the trace of X to compare it with the penalty parameter */
               DSDP_CALL( DSDPGetTraceX(sdpisolver->dsdp, &trace) );

#if 0 /* DSDP doesn't seem to adhere to its own feasiblity tolerance */
               assert( trace < penaltyparam + sdpisolver->feastol ); /* solution should be primal feasible */
#endif

               /* if the relative gap is smaller than the tolerance, we return equality */
               if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )
               {
                  *penaltybound = TRUE;
                  SCIPdebugMessage("Tr(X) = %f == %f = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible.\n",
                        trace, penaltyparam);
               }
               else
                  *penaltybound = FALSE;
            }
         }
          BMSfreeBufferMemoryArray(sdpisolver->bufmem, &dsdpsol);
      }
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
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   if ( pdfeasible == DSDP_PDUNKNOWN )
      return FALSE;

   return TRUE;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   assert( primalfeasible != NULL );
   assert( dualfeasible != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   DSDP_CALL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   switch ( pdfeasible )
   {
   case DSDP_PDFEASIBLE:
      /* primal and dual feasible */
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;

   case DSDP_UNBOUNDED:
      /* dual unbounded and primal infeasible */
      *primalfeasible = FALSE;
      *dualfeasible = TRUE;
      break;

   case DSDP_INFEASIBLE:
      /* dual infeasible and primal unbounded */
      *primalfeasible = TRUE;
      *dualfeasible = FALSE;
      break;

   default: /* should only include DSDP_PDUNKNOWN */
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
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
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_LPERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );
   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_UNBOUNDED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL(DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible));

   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible) );

   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
      return FALSE;

   return TRUE;
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return FALSE;

   if ( ! sdpisolver->solved )
      return FALSE;

   DSDP_CALL_BOOL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   if ( reason == DSDP_CONVERGED )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the objective limit was reached */ /*lint -e{715}*/
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Method not implemented for DSDP, as objective limit is given as an ordinary LP-constraint, so in case the objective limit was "
      "exceeded, the problem will be reported as infeasible!\n");

   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPTerminationReason reason;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL( DSDPStopReason(sdpisolver->dsdp, &reason) );

   if ( reason == DSDP_MAX_IT )
      return TRUE;

   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
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
{
   DSDPTerminationReason reason;
   int dsdpreturn;

   assert( sdpisolver != NULL );

   if ( sdpisolver->dsdp == NULL || (! sdpisolver->solved) )
      return -1;

   if ( sdpisolver->timelimit )
      return 5;

   dsdpreturn = DSDPStopReason(sdpisolver->dsdp, &reason);

   if (dsdpreturn != 0)
   {
      SCIPerrorMessage("DSDP-Error <%d> in function call.\n", dsdpreturn);
      return 7;
   }

   switch ( reason )/*lint --e{788}*/
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
      return 5;

   default:
      return 7;
   }
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   assert( sdpisolver != NULL );

   return (SCIPsdpiSolverIsConverged(sdpisolver) && SCIPsdpiSolverIsPrimalFeasible(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework
 */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP-solver interface */
   )
{
   DSDPSolutionType pdfeasible;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   DSDP_CALL_BOOL(DSDPGetSolutionType(sdpisolver->dsdp, &pdfeasible));

   if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsConverged(sdpisolver);
}

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{
   assert( sdpisolver != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   if ( sdpisolver->penalty && ! sdpisolver->feasorig )
   {
      /* in this case we cannot really trust the solution given by DSDP, since changes in the value of r much less than epsilon can
       * cause huge changes in the objective, so using the objective value given by DSDP is numerically more stable */
      DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
      *objval = - (*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   }
   else
   {
      SCIP_Real* dsdpsol;
      int dsdpnvars;
      int v;

      /* since the objective value given by DSDP sometimes differs slightly from the correct value for the given solution,
       * we get the solution from DSDP and compute the correct objective value */
      dsdpnvars = sdpisolver->penaltyworbound ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars; /* in the first case we added r as an explicit var */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars) );
      DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, dsdpnvars) ); /* last entry needs to be the number of variables, will return an error otherwise */

      /* use the solution to compute the correct objective value */
      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
         *objval += sdpisolver->objcoefs[v] * dsdpsol[v];

      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars);/*lint !e737 */
   }

   /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective as well */
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
      SCIP_Real* dsdpsol;
      int dsdpnvars;
      int v;

      dsdpnvars = sdpisolver->penaltyworbound ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars; /* in the first case we added r as an explicit var */
      BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars) );
      DSDP_CALL( DSDPGetY(sdpisolver->dsdp, dsdpsol, dsdpnvars) ); /* last entry needs to be the number of variables, will return an error otherwise */

      /* insert the entries into dualsol, for non-fixed vars we copy those from DSDP */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtodsdpmapper[v] > -1 )
         {
            /* minus one because the inputtosdpamapper starts at 1, but the array starts at 0 */
            dualsol[v] = dsdpsol[sdpisolver->inputtodsdpmapper[v] - 1];
         }
         else
         {
            /* the fixed value was saved at the beginning */
            dualsol[v] = sdpisolver->fixedvarsval[- sdpisolver->inputtodsdpmapper[v] - 1]; /*lint !e679*/
         }
      }

      if ( objval != NULL )
      {
         if ( sdpisolver->penalty && ! sdpisolver->feasorig )
         {
            /* in this case we cannot really trust the solution given by DSDP, since changes in the value of r much less than epsilon can
             * cause huge changes in the objective, so using the objective value given by DSDP is numerically more stable */
            DSDP_CALL( DSDPGetDObjective(sdpisolver->dsdp, objval) );
            *objval = - (*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
         }
         else
         {
            /* use the solution to compute the correct objective value */
            *objval = 0.0;
            for (v = 0; v < sdpisolver->nactivevars; v++)
               *objval += sdpisolver->objcoefs[v] * dsdpsol[v];
         }

         /* as we didn't add the fixed (lb = ub) variables to dsdp, we have to add their contributions to the objective as well */
         *objval += sdpisolver->fixedvarsobjcontr;
      }

      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpsol, dsdpnvars);/*lint !e737 */
   }
   else if ( objval != NULL )
   {
      SCIP_CALL( SCIPsdpiSolverGetObjval(sdpisolver, objval) );
   }

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X for the preoptimal solution */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverGetPreoptimalPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz */
   int*                  startXnblocknonz    /**< array to return number of nonzeros for row/col/val-arrays in each block
                                              *   or first entry equal to -1 if no primal solution is available */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");

   return SCIP_PLUGINNOTFOUND;
}

/** gets preoptimal dual solution vector and primal matrix for warmstarting purposes
 *
 *  @note The last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity.
 */ /*lint -e{715}*/
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
{  /*lint --e{715}*/
   int v;

   assert( sdpisolver != NULL );
   assert( success != NULL );
   assert( dualsol != NULL );

   /* we do not want to return a primal matrix, since it is not used by DSDP for warmstarting purposes */
   assert( nblocks == -1 );

   if ( ! sdpisolver->preoptimalsolexists )
   {
      SCIPdebugMessage("Failed to retrieve preoptimal solution for warmstarting purposes. \n");
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (v = 0; v < sdpisolver->nvars; v++)
   {
      if (sdpisolver->inputtodsdpmapper[v] > -1)
      {
         /* minus one because the inputtodsdpmapper gives the dsdp indices which start at one, but the array starts at 0 */
         dualsol[v] = sdpisolver->preoptimalsol[sdpisolver->inputtodsdpmapper[v] - 1];
      }
      else
      {
         /* this is the value that was saved when inserting, as this variable has lb=ub */
         dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtodsdpmapper[v]) - 1]; /*lint !e679*/
      }
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
{
   SCIP_Real* lbvalsdsdp;
   SCIP_Real* ubvalsdsdp;
   int i;

   assert( sdpisolver != NULL );
   assert( lbvals != NULL );
   assert( ubvals != NULL );
   CHECK_IF_SOLVED( sdpisolver );

   /* allocate memory for the arrays given to DSDP */
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &lbvalsdsdp, sdpisolver->nactivevars) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpisolver->blkmem, &ubvalsdsdp, sdpisolver->nactivevars) );

   /* get the values for the active variables from DSDP */
   DSDP_CALL( BConeCopyX(sdpisolver->bcone, lbvalsdsdp, ubvalsdsdp, sdpisolver->nactivevars) );

   /* copy them to the right spots of lbvars & ubvars */
   for (i = 0; i < sdpisolver->nvars; i++)
   {
      if ( sdpisolver->inputtodsdpmapper[i] < 0 )
      {
         /* if the variable was fixed, it didn't exist in the relaxation, so we set the value to 0
          * (as DSDP already uses this value for unbounded vars) */
         lbvals[i] = 0;
         ubvals[i] = 0;
      }
      else
      {
         lbvals[i] = lbvalsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
         ubvals[i] = ubvalsdsdp[sdpisolver->inputtodsdpmapper[i] - 1];
      }
   }

   /* free allocated memory */
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &ubvalsdsdp, sdpisolver->nactivevars);
   BMSfreeBlockMemoryArrayNull(sdpisolver->blkmem, &lbvalsdsdp, sdpisolver->nactivevars);

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
   SCIP_Real* primalvals;
   int nprimalvals;
   int ind = 0;
   int i;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( lplhs != NULL );
   assert( lprhs != NULL );
   assert( lhsvals != NULL );
   assert( rhsvals != NULL );

   if ( nlpcons <= 0 )
      return SCIP_OKAY;

   /* get primal solution for LP part from DSDP */
   DSDP_CALL( LPConeGetXArray(sdpisolver->lpcone, &primalvals, &nprimalvals) );

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
      {
         lhsvals[i] = primalvals[ind];
         ++ind;
      }
      else
         lhsvals[i] = 0.0;

      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         rhsvals[i] = primalvals[ind];
         ++ind;
      }
      else
         rhsvals[i] = 0.0;

      assert( ind <= nprimalvals );
   }
   assert( ind == nprimalvals );

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X (including lp block) */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverGetPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** returns the primal matrix X
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverGetPrimalMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block
                                              *   output: number of nonzeros in each block */
   int**                 startXrow,          /**< pointer to store row indices of X */
   int**                 startXcol,          /**< pointer to store column indices of X */
   SCIP_Real**           startXval           /**< pointer to store values of X */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
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
   SCIP_Real**           primalmatrices      /**< pointer to store values of the primal matrix */
   )
{  /*lint --e{715}*/
   int b;

   assert( sdpisolver != NULL );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( indchanges != NULL );
   assert( nremovedinds != NULL );
   assert( blockindchanges != NULL );
   assert( primalmatrices != NULL );

   DSDP_CALL( DSDPComputeX(sdpisolver->dsdp) );

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
         SCIP_Real* X;   /* the upper triangular entries of matrix X */
         SCIP_Real val;
         int idx = 0;
         int n;
         int i;

#ifndef NDEBUG
         int redsize;
         redsize = blocksize - nremovedinds[b];
#endif
         DSDP_CALL( SDPConeGetXArray(sdpisolver->sdpcone, b - blockindchanges[b], &X, &n) );
         assert( n == redsize * (redsize + 1)/2 );

         /* DSDP_CALL( SDPConeViewX(sdpisolver->sdpcone, b - blockindchanges[b], blocksize, X, n) ); */

         /* fill in matrix */
         for (j = 0; j < blocksize; ++j)
         {
            if ( indchanges[b][j] >= 0 )
            {
               assert( 0 <= j - indchanges[b][j] && j - indchanges[b][j] < redsize );

               for (i = 0; i <= j; ++i)
               {
                  if ( indchanges[b][i] >= 0 )
                  {
                     assert( 0 <= i - indchanges[b][i] && i - indchanges[b][i] < redsize );
                     assert( 0 <= idx && idx < redsize * (redsize + 1)/2 );
                     val = X[idx++];
                     primalmatrices[b][i * blocksize + j] = val;
                     primalmatrices[b][j * blocksize + i] = val;
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** return the maximum absolute value of the optimal primal matrix */ /*lint -e{715}*/
SCIP_Real SCIPsdpiSolverGetMaxPrimalEntry(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
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
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpisolver != NULL );
   assert( iterations != NULL );

   if ( sdpisolver->timelimitinitial )
      *iterations = 0;
   else
      *iterations = sdpisolver->niterations;

   return SCIP_OKAY;
}

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   )
{
   assert( sdpisolver != NULL );
   assert( calls != NULL );

   if ( sdpisolver->timelimitinitial )
      *calls = 0;
   else
      *calls = sdpisolver->nsdpcalls;

   return SCIP_OKAY;
}

/** gets the settings used by the SDP solver for the last solve call */
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   )
{
   assert( sdpisolver != NULL );
   assert( usedsetting != NULL );

   if ( ! SCIPsdpiSolverIsAcceptable(sdpisolver) )
      *usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   else
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
{
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
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
      *dval = sdpisolver->penaltyparam;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      *dval = sdpisolver->objlimit;
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
{
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
      sdpisolver->penaltyparam = dval;
      SCIPdebugMessage("Setting sdpisolver penaltyparameter to %g.\n", dval);
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIPdebugMessage("Setting sdpisolver objlimit to %g.\n", dval);
      sdpisolver->objlimit = dval;
      break;
   case SCIP_SDPPAR_LAMBDASTAR:
      SCIPdebugMessage("Parameter SCIP_SDPPAR_LAMBDASTAR not used by DSDP.\n"); /* this parameter is only used by SDPA */
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
{
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
{
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

/** compute and set lambdastar (only used for SDPA) */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Lambdastar parameter not used by DSDP.\n"); /* this parameter is only used by SDPA */

   return SCIP_OKAY;
}

/** compute and set the penalty parameter */
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   )
{
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( penaltyparam != NULL );

   compval = PENALTYPARAM_FACTOR * maxcoeff;

   if ( compval < MIN_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MIN_PENALTYPARAM);
      sdpisolver->penaltyparam = MIN_PENALTYPARAM;
      *penaltyparam = MIN_PENALTYPARAM;
   }
   else if ( compval > MAX_PENALTYPARAM )
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_PENALTYPARAM);
      sdpisolver->penaltyparam = MAX_PENALTYPARAM;
      *penaltyparam = MAX_PENALTYPARAM;
   }
   else
   {
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", compval);
      sdpisolver->penaltyparam = compval;
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
{
   SCIP_Real compval;

   assert( sdpisolver != NULL );
   assert( maxpenaltyparam != NULL );

   compval = penaltyparam * MAXPENALTYPARAM_FACTOR;

   if ( compval < MAX_MAXPENALTYPARAM )
   {
      *maxpenaltyparam = compval;
      SCIPdebugMessage("Setting maximum penaltyparameter to %f.\n", compval);
   }
   else
   {
      *maxpenaltyparam = MAX_MAXPENALTYPARAM;
      SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_MAXPENALTYPARAM);
   }

   /* if the maximum penalty parameter is smaller than the initial penalty paramater, we decrease the initial one correspondingly */
   if ( sdpisolver->penaltyparam > *maxpenaltyparam )
   {
      SCIPdebugMessage("Decreasing penaltyparameter of %f to maximum penalty paramater of %f.\n", sdpisolver->penaltyparam, *maxpenaltyparam);
      sdpisolver->penaltyparam = *maxpenaltyparam;
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

/** writes SDP to a file */ /*lint -e{715}*/
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/**@} */
