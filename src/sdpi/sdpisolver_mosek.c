/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2023 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2023 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/
/*#define SCIP_DEBUG_PRINTTOFILE  *//* prints each problem inserted into MOSEK to the file mosek.task */
/* #define SCIP_PRINT_SOLU        /\* prints each solution (primal and dual) reported from MOSEK to the screen *\/ */

/**@file   sdpisolver_mosek.c
 * @brief  interface for MOSEK
 * @author Tristan Gally
 * @author Marc Pfetsch
 *
 * As MOSEK solves the primal instead of the dual problem, for solving the problem
 *
 *   \f{eqnarray*}{
 *      \min & & b^T y \\
 *      \mbox{s.t.} & & \sum_{i \in I} A_i^{(k)} y_i - A_0^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & \sum_{i \in I} d_{ij} y_i \geq c_j \quad \forall \ j \in J, \\
 *      & & \ell_i \leq y_i \leq u_i \quad \forall \ i \in I
 *   \f}
 *
 * we insert the problem
 *
 *   \f{eqnarray*}{
 *      \max & & \sum_{k \in K} A_0^{(k)} \bullet X^{(k)} + \sum_{j \in J} c_j x_j - \sum_{i \in I_u} u_i v_i + \sum_{i \in I_\ell} \ell_i w_i \\
 *      \mbox{s.t.} & & \sum_{k \in K} A_i^{(k)} \bullet X^{(k)} + \sum_{j \in J} d_{ij} x_j - 1_{\{u_i < \infty\}} v_i + 1_{\{\ell_i > -\infty\}} w_i = b_i \quad \forall \ i \in I,\\
 *      & & X^{(k)} \succeq 0 \quad \forall \ k \in K, \\
 *      & & x_j \geq 0 \quad \forall \ j \in J,\\
 *      & & v_i \geq 0 \quad \forall \ i \in I_u,\\
 *      & & w_i \geq 0 \quad \forall \ i \in I_\ell,
 *   \f}
 *
 * where \f$I_\ell := \{i \in I: \ell_i > -\infty\}\f$ and \f$I_u := \{i \in I: u < \infty\}\f$.
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"
#include "sdpi/sdpiclock.h"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for SCIPsnprintf() */
#include "mosek.h"                           /* for MOSEK routines */
#include "sdpi/sdpsolchecker.h"              /* to check solution with regards to feasibility tolerance */
#include "scip/pub_message.h"                /* for debug and error message */
#include "tinycthread/tinycthread.h"         /* for thread local environments */

/* @todo Use MSK_putexitfunc to catch errors.
 * @todo Think about what to do with near optimality etc. (If MOSEK cannot compute a solution that has the prescribed accuracy, then it will
 *       multiply the termination tolerances with MSK_DPAR_INTPNT_CO_TOL_NEAR_REL. If the solution then satisfies the termination criteria, then
 *       the solution is denoted near optimal, near feasible and so forth.)
 */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e10     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR         1e6      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */
#define PENALTYBOUNDTOL             1E-3     /**< if the relative gap between Tr(X) and penaltyparam for a primal solution of the penaltyformulation
                                              *   is bigger than this value, it will be reported to the sdpi */
#define INFEASFEASTOLCHANGE         0.1      /**< change feastol by this factor if the solution was found to be infeasible with regards to feastol */
#define INFEASMINFEASTOL            1E-9      /**< minimum value for feasibility tolerance when encountering problems with regards to tolerance;
                                              *   @todo Think about doing this for absolute feastol */
#define CONVERT_ABSOLUTE_TOLERANCES TRUE     /**< should absolute tolerances be converted to relative tolerances for MOSEK */
#define CHECKVIOLATIONS             0        /**< Should the violations reported from MOSEK should be checked? */
#if MSK_VERSION_MAJOR >= 9
#define NEAR_REL_TOLERANCE           1.0     /**< MOSEK will multiply all tolerances with this factor after stalling */
#endif

#ifdef SCIP_THREADSAFE
   #if defined(_Thread_local)
      /* Use thread local environment in order to not create a new environment for each new LP. */
      static _Thread_local MSKenv_t reusemskenv = NULL;
      static _Thread_local int numsdp = 0;
      #define SCIP_REUSEENV
   #endif
#else
   /* Global Mosek environment in order to not create a new environment for each new LP. This is not thread safe. */
   static MSKenv_t reusemskenv = NULL;
   static int numsdp = 0;
   #define SCIP_REUSEENV
#endif


/** data used for SDP interface */
struct SCIP_SDPiSolver
{
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler for printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   BMS_BUFMEM*           bufmem;             /**< buffer memory */
   MSKenv_t              mskenv;             /**< MOSEK environement */
   MSKtask_t             msktask;            /**< MOSEK task */
   SCIP_Real             opttime;            /**< time spend in optimziation */
   int                   nvars;              /**< number of input variables */
   int                   nactivevars;        /**< number of variables present in the dual problem in MOSEK (nvars minus the number of variables with lb = ub) */
   int                   maxnvars;           /**< size of the arrays inputtomosekmapper, mosektoinputmapper, fixedvarsval, and objcoefs */
   int*                  inputtomosekmapper; /**< entry i gives the index of input variable i in MOSEK (starting from 0) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int*                  mosektoinputmapper; /**< entry i gives the original index of the i-th variable in MOSEK (indices go from 0 to nactivevars-1) */
   SCIP_Real*            fixedvarsval;       /**< entry i gives the lower and upper bound of the i-th fixed variable */
   SCIP_Real             fixedvarsobjcontr;  /**< total contribution to the objective of all fixed variables, computed as sum obj * val */
   SCIP_Real*            objcoefs;           /**< objective coefficients of all active variables */
   int                   nvarbounds;         /**< number of variable bounds given to MOSEK, length of varboundpos */
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
   SCIP_Bool             penalty;            /**< was the problem last solved using a penalty formulation */
   SCIP_Bool             feasorig;           /**< was the last problem solved with a penalty formulation and with original objective coefficents
                                              *   and the solution was feasible for the original problem? */
   SCIP_Bool             rbound;             /**< was the penalty parameter bounded during the last solve call */
   MSKrescodee           terminationcode;    /**< reason for termination of the last call to the MOSEK-optimizer */
   MSKsolstae            solstat;            /**< solution status of last call to MOSEK-optimizer */
   SCIP_Bool             timelimit;          /**< was the solver stopped because of the time limit? */
   int                   nthreads;           /**< number of threads the SDP solver should use (-1 = number of cores) */
   int                   niterations;        /**< number of SDP-iterations since the last solve call */
   int                   nsdpcalls;          /**< number of SDP-calls since the last solve call */
   SCIP_Bool             scaleobj;           /**< whether the objective should be scaled */
   SCIP_Real             objscalefactor;     /**< objective scaling factor */
#ifdef SCIP_REUSEENV
   int*                  numsdp;             /**< pointer to count on number of tasks in environment */
   MSKenv_t*             reusemskenv;        /**< pointer to reused Mosek environment */
#endif
};


/*
 * Local Methods
 */

/** Calls a MOSEK function and transforms the return-code to a SCIP_LPERROR if needed. */
#define MOSEK_CALL(x)  do                                                                                    \
                      {                                                                                      \
                         MSKrescodee _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK ) \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", (int)_mosekerrorcode_); \
                            return SCIP_LPERROR;                                                             \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as MOSEK_CALL, but used for functions returning a boolean. */
#define MOSEK_CALL_BOOL(x)  do                                                                               \
                      {                                                                                      \
                         MSKrescodee _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK ) \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", (int)_mosekerrorcode_); \
                            return FALSE;                                                                    \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

/** Same as MOSEK_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced. */
#define MOSEK_CALLM(x) do                                                                                    \
                      {                                                                                      \
                         MSKrescodee _mosekerrorcode_;                                                               \
                         if ( (_mosekerrorcode_ = (x)) != MSK_RES_OK ) \
                         {                                                                                   \
                            SCIPerrorMessage("MOSEK-Error <%d> in function call.\n", (int)_mosekerrorcode_); \
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

/** print string using message handler of SCIP */
static
void MSKAPI printstr(
#if MSK_VERSION_MAJOR >= 9
   MSKuserhandle_t       handler,            /**< error handler */
   const char*           str                 /**< string to print */
#else
   void*                 handler,            /**< error handler */
   MSKCONST char         str[]               /**< string to print */
#endif
   )
{  /*lint --e{715}*/
#if 0
   char errstr[32];
   snprintf(errstr, 32, "MOSEK Error %d", MSK_RES_ERR_DUP_NAME);
   if ( strncmp(errstr, str, strlen(errstr)) == 0 )
      return;
#endif

   SCIPmessagePrintInfo((SCIP_MESSAGEHDLR *) handler, "MOSEK: %s", str);
}


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
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->inputtomosekmapper), sdpisolver->maxnvars, newsize) );
      BMS_CALL( BMSreallocBlockMemoryArray(sdpisolver->blkmem, &(sdpisolver->mosektoinputmapper), sdpisolver->maxnvars, newsize) );
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
   MSKrescodee rescodee;
   int majorver = 0;
   int minorver = 0;
#if MSK_VERSION_MAJOR < 9
   int build = 0;
#endif
   int revision = 0;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif

#if MSK_VERSION_MAJOR < 9
   rescodee = MSK_getversion(&majorver, &minorver, &build, &revision);/*lint !e123*/
#else
   rescodee = MSK_getversion(&majorver, &minorver, &revision);/*lint !e123*/
#endif

   if ( rescodee != MSK_RES_OK )
      return "MOSEK";

#ifndef NDEBUG
#if MSK_VERSION_MAJOR < 9
   snprintfreturn = SCIPsnprintf(solvername, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d.%d", majorver, minorver, build, revision);/*lint !e123*/
#else
   snprintfreturn = SCIPsnprintf(solvername, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d", majorver, minorver, revision);/*lint !e123*/
#endif
   assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
#if MSK_VERSION_MAJOR < 9
   (void) SCIPsnprintf(solvername, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d.%d", majorver, minorver, build, revision);
#else
   (void) SCIPsnprintf(solvername, SCIP_MAXSTRLEN, "MOSEK %d.%d.%d", majorver, minorver, revision);
#endif
#endif

   return solvername;
}

/** gets description of SDP-solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Homogeneous, self-dual interior-point solver for semidefinite programming developed by MOSEK ApS (http://www.mosek.com)";
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
   return (void*) NULL;
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
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
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

   SCIPdebugMessage("Calling SCIPsdpiCreate \n");

   BMS_CALL( BMSallocBlockMemory(blkmem, sdpisolver) );

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->bufmem = bufmem;

#ifdef SCIP_REUSEENV
   if ( reusemskenv == NULL )
   {
      assert( numsdp == 0 );
      MOSEK_CALL( MSK_makeenv(&reusemskenv, NULL) );
   }
   (*sdpisolver)->mskenv = reusemskenv;
   ++numsdp;

   /* remember address of numsdp and reusemskenv, in case they are thread-local and SCIPsdpisolverFree is called from different thread */
   (*sdpisolver)->numsdp = &numsdp;
   (*sdpisolver)->reusemskenv = &reusemskenv;
#else
   MOSEK_CALL( MSK_makeenv(&((*sdpisolver)->mskenv), NULL) );/*lint !e641*/ /* the NULL-argument is a debug file, but setting this will spam the whole folder */
#endif

   /* this will be properly initialized then calling solve */
   (*sdpisolver)->msktask = NULL;
   (*sdpisolver)->opttime = 0.0;

   (*sdpisolver)->nvars = 0;
   (*sdpisolver)->nactivevars = 0;
   (*sdpisolver)->inputtomosekmapper = NULL;
   (*sdpisolver)->mosektoinputmapper = NULL;
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
   (*sdpisolver)->terminationcode = MSK_RES_OK;
   (*sdpisolver)->solstat = MSK_SOL_STA_UNKNOWN;
   (*sdpisolver)->timelimit = FALSE;
   (*sdpisolver)->niterations = 0;
   (*sdpisolver)->scaleobj = FALSE;
   (*sdpisolver)->objscalefactor = 1.0;

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

   if ( (*sdpisolver)->msktask != NULL )
   {
      MOSEK_CALL( MSK_deletetask(&((*sdpisolver)->msktask)) );/*lint !e641*/
   }

#ifdef SCIP_REUSEENV
   /* decrement the numsdp that belongs to the thread where SCIPsdpisolverCreate was called */
   assert( *(*sdpisolver)->numsdp > 0 );
   --(*(*sdpisolver)->numsdp);

   /* if numsdp reached zero, then also free the Mosek environment (that belongs to the thread where SCIPsdpisolverCreate was called) */
   if ( *(*sdpisolver)->numsdp == 0 )
   {
      /* free reused environment */
      MOSEK_CALL( MSK_deleteenv((*sdpisolver)->reusemskenv) );
      *(*sdpisolver)->reusemskenv = NULL;
   }
#else
   MOSEK_CALL( MSK_deleteenv(&((*sdpisolver)->mskenv)) );/*lint !e641*/
#endif

   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->varboundpos, 2 * (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->inputtomosekmapper, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->mosektoinputmapper, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->fixedvarsval, (*sdpisolver)->maxnvars);
   BMSfreeBlockMemoryArrayNull((*sdpisolver)->blkmem, &(*sdpisolver)->objcoefs, (*sdpisolver)->maxnvars);

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

/** reset the SDP-Counter to zero */
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
 */ /*lint --e{715}*/
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
{/*lint --e{715}*/
   int b;
   int i;
   int j;
   int blockvar;
   int v;
   int k;
   MSKint64t mosekindex;
   int ind;
   SCIP_Real* mosekvarbounds;
   int nfixedvars;
   int nlpvars;
   int* mosekblocksizes;
   SCIP_Real one = 1.0; /* MOSEK always wants a pointer to factors for a sum of matrices, we always use a single matrix with factor one */
   int* mosekrow;
   int* mosekcol;
   SCIP_Real* mosekval;
   int row;
   SCIP_Real val;
   SCIP_Real solvertimelimit;
#ifdef SCIP_MORE_DEBUG
   char name[SCIP_MAXSTRLEN];
#endif
#if CONVERT_ABSOLUTE_TOLERANCES
   SCIP_Real maxrhscoef = 0.0; /* MOSEK uses a relative feasibility tolerance, the largest rhs-coefficient is needed for converting the absolute tolerance */
#endif
   SCIP_Real maxabsobjcoef = 0.0;

   assert( sdpisolver != NULL );
   assert( sdpisolver->mskenv != NULL );
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

   /* create an empty task (second and third argument are guesses for maximum number of constraints and variables), if there already is one, delete it */
   if ( sdpisolver->msktask != NULL )
   {
      MOSEK_CALL( MSK_deletetask(&sdpisolver->msktask) );/*lint !e641*/
   }
   if ( penaltyparam < sdpisolver->epsilon )
   {
      MOSEK_CALLM( MSK_maketask(sdpisolver->mskenv, nvars, nsdpblocks - nremovedblocks + nlpcons + 2 * nvars, &sdpisolver->msktask) );/*lint !e641*/
   }
   else
   {
      MOSEK_CALLM( MSK_maketask(sdpisolver->mskenv, nvars + 1, nsdpblocks - nremovedblocks + nlpcons + 2 * nvars, &sdpisolver->msktask) );/*lint !e641*/
   }

#if MSK_VERSION_MAJOR >= 9
   MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_NEAR_REL, NEAR_REL_TOLERANCE) );
#endif

   /* redirect output to SCIP message handler */
#ifdef SCIP_MORE_DEBUG
   sdpisolver->sdpinfo = TRUE;
#endif

   if ( sdpisolver->sdpinfo )
   {
      MOSEK_CALL( MSK_linkfunctotaskstream(sdpisolver->msktask, MSK_STREAM_LOG, (MSKuserhandle_t) sdpisolver->messagehdlr, printstr) );/*lint !e641*/
   }

   /* set number of threads */
   if ( sdpisolver->nthreads > 0 )
   {
      MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_NUM_THREADS, sdpisolver->nthreads) );/*lint !e641*/
   }

   /* only increase the counter if we don't use the penalty formulation to stay in line with the numbers in the general interface (where this is still the
    * same SDP) */
   if ( penaltyparam < sdpisolver->epsilon )
   {
      SCIPdebugMessage("Inserting data into MOSEK for SDP (%d) \n", ++sdpisolver->sdpcounter);
   }
   else
   {
      SCIPdebugMessage("Inserting Data again into MOSEK for penalty formulation of SDP (%d) \n", sdpisolver->sdpcounter);
   }

   /* set the penalty and rbound flags accordingly */
   sdpisolver->penalty = (penaltyparam < sdpisolver->epsilon) ? FALSE : TRUE;
   sdpisolver->rbound = rbound;

   /* ensure memory for varboundpos, inputtomosekmapper, mosektoinputmapper and the fixed and active variable information */
   SCIP_CALL( ensureMappingDataMemory(sdpisolver, nvars) );
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekvarbounds, 2 * nvars) ); /*lint !e647*/

   /* find fixed variables */
   sdpisolver->nvars = nvars;
   sdpisolver->nactivevars = 0;
   sdpisolver->nvarbounds = 0;
   nfixedvars = 0;
   sdpisolver->fixedvarsobjcontr = 0.0;
   for (i = 0; i < nvars; i++)
   {
      if ( isFixed(sdpisolver, lb[i], ub[i]) )
      {
         sdpisolver->fixedvarsobjcontr += obj[i] * lb[i]; /* this is the value this fixed variable contributes to the objective */
         sdpisolver->fixedvarsval[nfixedvars] = lb[i]; /* if lb=ub, then this is the value the variable will have in every solution */
         nfixedvars++;
         sdpisolver->inputtomosekmapper[i] = -nfixedvars;
         SCIPdebugMessage("Fixing variable %d locally to %g for SDP %d in MOSEK.\n", i, lb[i], sdpisolver->sdpcounter);
      }
      else
      {
#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Variable %d becomes variable %d for SDP %d in MOSEK.\n", i, sdpisolver->nactivevars, sdpisolver->sdpcounter);
#endif

         sdpisolver->mosektoinputmapper[sdpisolver->nactivevars] = i;
         sdpisolver->inputtomosekmapper[i] = sdpisolver->nactivevars;
         sdpisolver->objcoefs[sdpisolver->nactivevars] = obj[i];

         if ( withobj )
         {
            /* the objective moves to the rhs */
            if ( REALABS(obj[i]) > maxabsobjcoef )
               maxabsobjcoef = REALABS(obj[i]);
         }

         if ( lb[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            mosekvarbounds[sdpisolver->nvarbounds] = lb[i];
            sdpisolver->varboundpos[sdpisolver->nvarbounds++] = -(sdpisolver->nactivevars + 1); /* negative sign means lower bound */
         }

         if ( ub[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            mosekvarbounds[sdpisolver->nvarbounds] = - ub[i]; /* we give the upper bounds a negative sign for the objective */
            sdpisolver->varboundpos[sdpisolver->nvarbounds++] = +(sdpisolver->nactivevars + 1); /* positive sign means upper bound */
         }

         sdpisolver->nactivevars++;
      }
   }
   assert( sdpisolver->nactivevars + nfixedvars == sdpisolver->nvars );

   /* adjust maxabsobjcoef in penalty formulation */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      if ( penaltyparam > maxabsobjcoef )
         maxabsobjcoef = penaltyparam;
   }

   /* if we want to solve without objective, we reset fixedvarsobjcontr */
   if ( ! withobj )
      sdpisolver->fixedvarsobjcontr = 0.0;

   /* determine total number of sides in LP-constraints */
   nlpvars = 0;
   if ( nlpcons > 0 )
   {
      for (i = 0; i < nlpcons; i++)
      {
         if ( lpindchanges[i] < 0 )
            continue;

         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
#if CONVERT_ABSOLUTE_TOLERANCES
            /* update largest rhs-entry */
            if ( REALABS(lplhs[i]) > maxrhscoef )
               maxrhscoef = REALABS(lplhs[i]);
#endif
            ++nlpvars;
         }

         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
#if CONVERT_ABSOLUTE_TOLERANCES
            /* update largest rhs-entry */
            if ( REALABS(lprhs[i]) > maxrhscoef )
               maxrhscoef = REALABS(lprhs[i]);
#endif
            ++nlpvars;
         }
      }
      assert( nlpvars <= 2 * nlpcons ); /* factor 2 comes from left- and right-hand-sides */
   }

   /* create matrix variables */
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekblocksizes, nsdpblocks - nremovedblocks) ); /*lint !e679 !e776*/

   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] > -1 )
      {
         assert( 0 <= blockindchanges[b] && blockindchanges[b] <= b && (b - blockindchanges[b]) <= (nsdpblocks - nremovedblocks) );
         mosekblocksizes[b - blockindchanges[b]] = sdpblocksizes[b] - nremovedinds[b];/*lint !e679*/
      }
   }
   MOSEK_CALLM( MSK_appendbarvars(sdpisolver->msktask, nsdpblocks - nremovedblocks, mosekblocksizes) );/*lint !e641*/

   /* create scalar variables (since we solve the primal problem, these are not the active variables but the dual variables for the
    * lp constraints and variable bounds) */
   MOSEK_CALLM( MSK_appendvars(sdpisolver->msktask, nlpvars + sdpisolver->nvarbounds) );/*lint !e641*/

   /* the variables for the LP constraints and variable bounds are non-negative */
   MOSEK_CALLM( MSK_putvarboundsliceconst(sdpisolver->msktask, 0, nlpvars + sdpisolver->nvarbounds, MSK_BK_LO, 0.0, MSK_INFINITY) );/*li

   /* append empty constraints (since we solve the primal problem, we get one constraint for each active variable) */
   MOSEK_CALLM( MSK_appendcons(sdpisolver->msktask, (penaltyparam < sdpisolver->epsilon) ? sdpisolver->nactivevars : sdpisolver->nactivevars + 1) );/*lint !e641*/

   /* set objective values for the matrix variables */
   if ( sdpconstnnonz > 0 )
   {
      i = 0;

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] >= 0 )
         {
            /* if some indices in the block were removed, we need to change indices accordingly */
            if ( nremovedinds[b] > 0 )
            {
               int* moseksdpconstrow;
               int* moseksdpconstcol;

               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstrow, sdpconstnblocknonz[b]) );
               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstcol, sdpconstnblocknonz[b]) );

               for (k = 0; k < sdpconstnblocknonz[b]; k++)
               {
                  /* rows and cols with active nonzeros should not be removed */
                  assert( -1 < indchanges[b][sdpconstrow[b][k]] && indchanges[b][sdpconstrow[b][k]] <= sdpconstrow[b][k] );
                  assert( -1 < indchanges[b][sdpconstcol[b][k]] && indchanges[b][sdpconstcol[b][k]] <= sdpconstcol[b][k] );

                  assert( 0 <= sdpconstrow[b][k] && sdpconstrow[b][k] <= sdpblocksizes[b] );
                  assert( 0 <= sdpconstcol[b][k] && sdpconstcol[b][k] <= sdpblocksizes[b] );

                  moseksdpconstrow[k] = sdpconstrow[b][k] - indchanges[b][sdpconstrow[b][k]];
                  moseksdpconstcol[k] = sdpconstcol[b][k] - indchanges[b][sdpconstcol[b][k]];

#if CONVERT_ABSOLUTE_TOLERANCES
                  /* update largest rhs-entry */
                  if ( REALABS(sdpconstval[b][k]) > maxrhscoef )
                     maxrhscoef = REALABS(sdpconstval[b][k]);
#endif
               }

               MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], sdpconstnblocknonz[b],
                     moseksdpconstrow, moseksdpconstcol, sdpconstval[b], &mosekindex) );/*lint !e641, !e679, !e747*/

               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstcol);
               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksdpconstrow);
            }
            else
            {
               MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], sdpconstnblocknonz[b],
                     sdpconstrow[b], sdpconstcol[b], sdpconstval[b], &mosekindex) );/*lint !e641, !e679, !e747*/
            }
            MOSEK_CALL( MSK_putbarcj(sdpisolver->msktask, i, 1, &mosekindex, &one) );/*lint !e641, !e747*/
            i++;
         }
      }
   }

   /* set objective values for the scalar variables */
   /* TODO: check for equality constraints */
   /* first for those corresponding to LP constraints in the dual */
   ind = 0;
   for (i = 0; i < nlpcons; i++)
   {
      if ( lpindchanges[i] < 0 )
         continue;

      /* left hand side */
      if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
      {
         MOSEK_CALL( MSK_putcj(sdpisolver->msktask, ind, lplhs[i]) );/*lint !e641*/
#ifdef SCIP_MORE_DEBUG
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lhs#%d", i);
         MOSEK_CALL( MSK_putvarname(sdpisolver->msktask, ind, name) );
#endif
         ++ind;
      }

      /* right hand side */
      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         MOSEK_CALL( MSK_putcj(sdpisolver->msktask, ind, - lprhs[i]) );/*lint !e641, !e644*/
#ifdef SCIP_MORE_DEBUG
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "rhs#%d", i - 1);
         MOSEK_CALL( MSK_putvarname(sdpisolver->msktask, ind, name) );
#endif
         ++ind;
      }
   }
   assert( ind == nlpvars );

   /* finally for those corresponding to variable bounds in the dual */
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      MOSEK_CALL( MSK_putcj(sdpisolver->msktask, nlpvars + i, mosekvarbounds[i]) );/*lint !e641*/ /* for the ub's we already added a negative sign in mosekvarbounds */

#if 0
      /* We currently do not include variable bounds in maxrhscoef, because it does not seem to be beneficial
       * overall. The bounds are very relevant for cardinality least square instances in which all variables are binary,
       * except for one continuous variable representing the objective value. The objective value can be
       * large. Enlarging maxrhscoef will not particularly help in this context, since the objective values are measured
       * relatively and the bounds are filtered out later anyway. */
#if CONVERT_ABSOLUTE_TOLERANCES
      if ( REALABS(mosekvarbounds[i]) > maxrhscoef )
         maxrhscoef = REALABS(mosekvarbounds[i]);
#endif
#endif

#ifdef SCIP_MORE_DEBUG
      if ( sdpisolver->varboundpos[i] < 0 ) /* lower bound */
      {
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lb#%d", sdpisolver->mosektoinputmapper[-1 * sdpisolver->varboundpos[i] - 1]);
         MOSEK_CALL( MSK_putvarname(sdpisolver->msktask, nlpvars + i, name) );
      }
      else /* upper bound */
      {
         assert( sdpisolver->varboundpos[i] > 0 ); /* we should not have value 0 so that we can clearly differentiate between positive and negative */
         /* give the variable a meaningful name for debug output */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ub#%d", sdpisolver->mosektoinputmapper[sdpisolver->varboundpos[i] - 1]);
         MOSEK_CALL( MSK_putvarname(sdpisolver->msktask, nlpvars + i, name) );
      }
#endif
   }

   /* set objective sense (since we want to minimize in the dual, we maximize in the primal) */
   MOSEK_CALL( MSK_putobjsense(sdpisolver->msktask, MSK_OBJECTIVE_SENSE_MAXIMIZE) );/*lint !e641*/

   /* start inserting the constraints */

   /* first add the matrices A_i */
   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] > -1 )
      {
         /* prepare memory */
         if ( nremovedinds[b] > 0 )
         {
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, sdpnnonz) );
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekcol, sdpnnonz) );
         }
         else
         {
            mosekrow = NULL;
            mosekcol = NULL;
         }

         for (blockvar = 0; blockvar < sdpnblockvars[b]; blockvar++)
         {
            v = sdpisolver->inputtomosekmapper[sdpvar[b][blockvar]];

            /* check if the variable is active */
            if ( v > -1 )
            {
               assert( v < sdpisolver->nactivevars );

               /* if there are removed indices, we have to adjust the column and row indices accordingly */
               if ( nremovedinds[b] > 0 )
               {
                  assert( sdpnblockvarnonz[b][blockvar] <= sdpnnonz );
                  for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
                  {
                     /* rows and cols with active nonzeros should not be removed */
                     assert( 0 <= indchanges[b][sdprow[b][blockvar][k]] && indchanges[b][sdprow[b][blockvar][k]] <= sdprow[b][blockvar][k] );
                     assert( 0 <= indchanges[b][sdpcol[b][blockvar][k]] && indchanges[b][sdpcol[b][blockvar][k]] <= sdpcol[b][blockvar][k] );

                     assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                     assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                     mosekrow[k] = sdprow[b][blockvar][k] - indchanges[b][sdprow[b][blockvar][k]];
                     mosekcol[k] = sdpcol[b][blockvar][k] - indchanges[b][sdpcol[b][blockvar][k]];
                  }
                  assert( k == sdpnblockvarnonz[b][blockvar] );

                  MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) k,
                        mosekrow, mosekcol, sdpval[b][blockvar], &mosekindex) );/*lint !e641, !e679*/
               }
               else
               {
                  MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) sdpnblockvarnonz[b][blockvar],
                        sdprow[b][blockvar], sdpcol[b][blockvar], sdpval[b][blockvar], &mosekindex) );/*lint !e641, !e679*/
               }

               MOSEK_CALL( MSK_putbaraij(sdpisolver->msktask, v, b - blockindchanges[b], (long long) 1, &mosekindex, &one) );/*lint !e641*/
            }
         }
         BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekcol);
         BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekrow);
      }
   }

   /* add the identity matrix corresponding to the penalty variable */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      int* identityindices;
      SCIP_Real* identityvalues;

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &identityindices, mosekblocksizes[b - blockindchanges[b]]) );
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &identityvalues, mosekblocksizes[b - blockindchanges[b]]) );

            for (i = 0; i < mosekblocksizes[b - blockindchanges[b]]; i++)
            {
               identityindices[i] = i;
               identityvalues[i] = 1.0;
            }
            MOSEK_CALL( MSK_appendsparsesymmat(sdpisolver->msktask, mosekblocksizes[b - blockindchanges[b]], (long long) mosekblocksizes[b - blockindchanges[b]],
                  identityindices, identityindices, identityvalues, &mosekindex) );/*lint !e641, !e679*/
            MOSEK_CALL( MSK_putbaraij(sdpisolver->msktask, sdpisolver->nactivevars, b - blockindchanges[b], (long long) 1, &mosekindex, &one) );/*lint !e641, !e679*/

            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &identityvalues);
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &identityindices);
         }
      }
   }

   /* add the entries corresponding to the lp-constraints in the dual problem */
   if ( lpnnonz > 0 )
   {
      int varcnt = 0;

      if ( penaltyparam < sdpisolver->epsilon )
      {
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, lpnnonz) );
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekval, lpnnonz) );
      }
      else
      {
         /* one extra entry for the penalty-constraint Trace = Gamma */
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekrow, lpnnonz + 1) );/*lint !e776*/
         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &mosekval, lpnnonz + 1) );/*lint !e776*/
      }

      for (i = 0; i < nlpcons; ++i)
      {
         int mosekind = 0;
         int nextbeg;

         if ( lpindchanges[i] < 0 )
            continue;

         assert( 0 <= lpbeg[i] && lpbeg[i] < lpnnonz );
         if ( i == nlpcons - 1 )
            nextbeg = lpnnonz;
         else
            nextbeg = lpbeg[i+1];

         for (j = lpbeg[i]; j < nextbeg; ++j)
         {
            assert( 0 <= lpind[j] && lpind[j] < nvars );
            v = sdpisolver->inputtomosekmapper[lpind[j]];
            if ( v >= 0 )
            {
               assert( v < sdpisolver->nactivevars );
               mosekrow[mosekind] = v;
               mosekval[mosekind++] = lpval[j];
            }
         }
         /* all rows with at most one nonzero should have been sorted out, except when checking the Slater condition */
         assert( mosekind >= 1 );

         /* add the additional entry for the penalty constraint Trace = Gamma */
         if ( penaltyparam >= sdpisolver->epsilon )
         {
            mosekrow[mosekind] = sdpisolver->nactivevars;
            mosekval[mosekind++] = 1.0;
         }
         assert( mosekind <= lpnnonz + 1 );

         /* treat left hand side */
         if ( lplhs[i] > - SCIPsdpiSolverInfinity(sdpisolver) )
         {
            MOSEK_CALL( MSK_putacol(sdpisolver->msktask, varcnt++, mosekind, mosekrow, mosekval) );/*lint !e641*/
         }

         /* treat right hand side */
         if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
         {
            /* multiply column by -1 */
            for (j = 0; j < (penaltyparam < sdpisolver->epsilon ? mosekind : mosekind - 1); j++)/*lint !e644*/
               mosekval[j] *= -1.0;

            MOSEK_CALL( MSK_putacol(sdpisolver->msktask, varcnt++, mosekind, mosekrow, mosekval) );/*lint !e641*/
         }
      }
      assert( varcnt == nlpvars );

      BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekval);
      BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &mosekrow);
   }

   /* finally add the entries corresponding to varbounds in the dual problem, we get exactly one entry per variable */
   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* lower bound */
         row = - sdpisolver->varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
         val = 1.0;
      }
      else
      {
         /* upper bound */
         assert( sdpisolver->varboundpos[i] > 0 ); /* we should not have a zero as we wanted a clear differentiation between positive and negative */
         row = sdpisolver->varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
         val = -1.0;
      }
      MOSEK_CALL( MSK_putacol(sdpisolver->msktask, nlpvars + i, 1, &row, &val) );/*lint !e641*/
   }

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

   /* make all constraints equality constraints with right-hand side b_i (or 0 if we solve without objective) */
   for (i = 0; i < sdpisolver->nactivevars; i++)
   {
      if ( withobj )
      {
         SCIP_Real objcoef;

         objcoef = sdpisolver->objcoefs[i];
         assert( objcoef == obj[sdpisolver->mosektoinputmapper[i]] );
         objcoef /= sdpisolver->objscalefactor;
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, i, MSK_BK_FX, objcoef, objcoef) );/*lint !e641*/
      }
      else
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, i, MSK_BK_FX, 0.0, 0.0) );/*lint !e641*/
      }

#ifdef SCIP_MORE_DEBUG
      /* give the constraint a meaningful name for debug output */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "var#%d", sdpisolver->mosektoinputmapper[i]);
      MOSEK_CALL( MSK_putconname(sdpisolver->msktask, i, name) );
#endif
   }

   /* the penalty constraint has right-hand side Gamma, it is a <=-inequality if r was bounded and an equality constraint otherwise */
   if ( penaltyparam >= sdpisolver->epsilon )
   {
      SCIP_Real p;

      p = penaltyparam / sdpisolver->objscalefactor;
      if ( rbound )
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, sdpisolver->nactivevars, MSK_BK_UP, - MSK_INFINITY, p) );/*lint !e641*/
      }
      else
      {
         MOSEK_CALL( MSK_putconbound(sdpisolver->msktask, sdpisolver->nactivevars, MSK_BK_FX, p, p) );/*lint !e641*/
      }
#ifdef SCIP_MORE_DEBUG
      /* give the constraint a meaningful name for debug output */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "penalty");
      MOSEK_CALL( MSK_putconname(sdpisolver->msktask, i, name) );
#endif
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

      /* set epsilon and feasibility tolerance (note that the dual in MOSEK is the problem we are interested in, so this is where we use feastol,
       * since MOSEK works with relative tolerance, we adjust our absolute tolerance accordingly, so that any solution satisfying the relative
       * tolerance in MOSEK satisfies our absolute tolerance) */
#if CONVERT_ABSOLUTE_TOLERANCES
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, sdpisolver->gaptol) );
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, sdpisolver->sdpsolverfeastol / (1.0 + maxrhscoef)) );
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, sdpisolver->sdpsolverfeastol / MAX(1.0, maxabsobjcoef)) );
      SCIPdebugMessage("Setting tolerances for MOSEK: feastol = %.12g (maxrhscoef = %.12g); gaptol = %.12g; infeastol = %.12g (maxobjcoef = %.12g)\n",
         sdpisolver->sdpsolverfeastol / (1.0 + maxrhscoef), maxrhscoef, sdpisolver->gaptol, sdpisolver->sdpsolverfeastol / (1.0 + maxabsobjcoef), maxabsobjcoef);
#else
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_PFEAS, sdpisolver->gaptol) );
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, sdpisolver->sdpsolverfeastol) );
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, sdpisolver->sdpsolverfeastol) );
#endif
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_MU_RED, sdpisolver->gaptol) );
      MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, sdpisolver->gaptol) );

      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
      {
         MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_OPTIMIZER_MAX_TIME, solvertimelimit) );/*lint !e641*/
      }

      /* set objective cutoff */
      if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->objlimit) )
      {
         MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_UPPER_OBJ_CUT, sdpisolver->objlimit / sdpisolver->objscalefactor) );
      }

      /* turn presolving on/off */
      if ( sdpisolver->usepresolving )
      {
         SCIPdebugMessage("Turning presolving on.\n");
         MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_ON) );
      }
      else
      {
         SCIPdebugMessage("Turning presolving off.\n");
         MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF) );
      }

      /* turn scaling on/off */
      if ( sdpisolver->usescaling )
      {
         MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_FREE) );
      }
      else
      {
         MOSEK_CALL( MSK_putintparam(sdpisolver->msktask, MSK_IPAR_INTPNT_SCALING, MSK_SCALING_NONE) );
      }

      /* print parameters if asked to */
#ifdef SCIP_PRINT_PARAMETERS
      MOSEK_CALL( MSK_printparam(sdpisolver->msktask) );
#endif

      /* write to file if asked to */
#ifdef SCIP_DEBUG_PRINTTOFILE
      SCIP_CALL( SCIPsdpiSolverWriteSDP(sdpisolver, "mosek.task") );
#endif

      /* solve the problem */
      MOSEK_CALL( MSK_optimizetrm(sdpisolver->msktask, &sdpisolver->terminationcode) );/*lint !e641*/
      MOSEK_CALL( MSK_getdouinf(sdpisolver->msktask, MSK_DINF_OPTIMIZER_TIME, &sdpisolver->opttime) );
      MOSEK_CALL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &sdpisolver->solstat) );/*lint !e641*/

      if ( sdpisolver->sdpinfo )
      {
         MOSEK_CALL( MSK_optimizersummary(sdpisolver->msktask, MSK_STREAM_LOG) );/*lint !e641*/
         MOSEK_CALL( MSK_analyzesolution(sdpisolver->msktask, MSK_STREAM_LOG, MSK_SOL_ITR) );/*lint !e641*/
      }

      SCIPdebugMessage("Solved problem using MOSEK, return code %d.\n", sdpisolver->terminationcode);

      sdpisolver->solved = TRUE;

      /* update number of SDP-iterations and -calls */
      ++sdpisolver->nsdpcalls;
      MOSEK_CALL( MSK_getnaintinf(sdpisolver->msktask, "MSK_IINF_INTPNT_ITER", &(sdpisolver->niterations)) );/*lint !e641*/

      /* possibly repair status */
      if ( sdpisolver->terminationcode == MSK_RES_TRM_STALL || (sdpisolver->solstat == MSK_SOL_STA_UNKNOWN && sdpisolver->terminationcode != MSK_RES_TRM_MAX_TIME) )
      {
         SCIP_Real pobj;
         SCIP_Real pviolcon;
         SCIP_Real pviolvar;
         SCIP_Real pviolbarvar;
         SCIP_Real dobj;
         SCIP_Real dviolcon;
         SCIP_Real dviolvar;
         SCIP_Real dviolbarvar;

         MOSEK_CALL( MSK_getsolutioninfo(sdpisolver->msktask, MSK_SOL_ITR, &pobj, &pviolcon, &pviolvar, &pviolbarvar, NULL, NULL,
               &dobj, &dviolcon, &dviolvar, &dviolbarvar, NULL) );

         SCIPdebugMessage("Absolute primal violations: constraints: %g, variables: %g, SDP: %g.\n", pviolcon, pviolvar, pviolbarvar);
         SCIPdebugMessage("Absolute dual violations: constraints: %g, variables: %g, SDP: %g.\n", dviolcon, dviolvar, dviolbarvar); 
         if ( pviolcon <= sdpisolver->feastol && pviolvar <= sdpisolver->feastol && pviolbarvar <= sdpisolver->feastol
            && dviolcon <= sdpisolver->feastol && dviolvar <= sdpisolver->feastol && dviolbarvar <= sdpisolver->feastol )
         {
            if ( REALABS(dobj - pobj) <= sdpisolver->gaptol )
            {
               sdpisolver->terminationcode = MSK_RES_OK;
               sdpisolver->solstat = MSK_SOL_STA_OPTIMAL;
               SCIPdebugMessage("Detected stalling - repairing termination code and solution status to 'optimal'.\n");
            }
         }
      }

      /* if the problem has been stably solved but did not reach the required feasibility tolerance, even though the solver
       * reports feasibility, resolve it with adjusted tolerance */
#if CONVERT_ABSOLUTE_TOLERANCES
      feastol = sdpisolver->sdpsolverfeastol / (1 + maxrhscoef);
#else
      feastol = sdpisolver->sdpsolverfeastol;
#endif
      gaptol = sdpisolver->gaptol;

      while ( SCIPsdpiSolverIsAcceptable(sdpisolver) && SCIPsdpiSolverIsDualFeasible(sdpisolver) && penaltyparam < sdpisolver->epsilon && feastol >= INFEASMINFEASTOL )
      {
         SCIP_Real* solvector;
         SCIP_Bool infeasible;
         SCIP_Bool solveagain = FALSE;
         SCIP_Real opttime;
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

#if CHECKVIOLATIONS == 1
         {
            /* The following code obtains the primal and dual solution reported from MOSEK and checks the violations of
             * the variable bounds, and linear as well as SDP constraints. Since MOSEK seems to use a different LAPACK
             * routine to compute the eigenvalues of a symmetric matrix, the violations of the SDP constraints in the
             * primal and dual reported from MOSEK and computed by SCIP-SDP differ slightly. Apart from that, the primal
             * violations (variables and linear constraints) SCIP-SDP computes the same violations as MOSEK. For the
             * dual problem however, MOSEK uses slack variables to obtain equalities, whereas the original (dual)
             * problem in SCIP-SDP may have inequalities. Thus, MOSEK reports violations in the dual with respect to the
             * slack variables and the corresponding equalities, whereas SCIP-SDP computes violations with respect to
             * the original inequalities. Thus, those violations differ. Moreover, the dual variable violations reported
             * from MOSEK seem to be the violations of the constraints in the dual problem, and the dual constraint
             * violations from MOSEK seem to be the violation of the variables in the dual problem, since the dual of
             * the primal constraints are the dual variables.
             */
            SCIP_Real* solvectorprimal;
            SCIP_Real** solmatrices;
            SCIP_Real maxabsviolbndsd;
            SCIP_Real sumabsviolbndsd;
            SCIP_Real maxabsviolconsd;
            SCIP_Real sumabsviolconsd;
            SCIP_Real maxabsviolsdpd;
            SCIP_Real sumabsviolsdpd;
            SCIP_Real maxabsviolbndsp;
            SCIP_Real sumabsviolbndsp;
            SCIP_Real maxabsviolconsp;
            SCIP_Real sumabsviolconsp;
            SCIP_Real maxabsviolsdpp;
            SCIP_Real sumabsviolsdpp;
            SCIP_Bool checkinfeas;
            SCIP_Real pobj;
            SCIP_Real pviolcon;
            SCIP_Real pviolvar;
            SCIP_Real pviolbarvar;
            SCIP_Real dobj;
            SCIP_Real dviolcon;
            SCIP_Real dviolvar;
            SCIP_Real dviolbarvar;
            int nprimalvars;
            int nprimalmatrixvars;
            int blocksize;

            /* get violations from Mosek */
            MOSEK_CALL( MSK_getsolutioninfo(sdpisolver->msktask, MSK_SOL_ITR, &pobj, &pviolcon, &pviolvar, &pviolbarvar, NULL, NULL,
               &dobj, &dviolcon, &dviolvar, &dviolbarvar, NULL) );

            nprimalvars = nlpvars + sdpisolver->nvarbounds;
            nprimalmatrixvars = nsdpblocks - nremovedblocks;

            /* get current solution */
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solvectorprimal, nprimalvars) );
            BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solmatrices, nprimalmatrixvars) );
            for (i = 0; i < nprimalmatrixvars; i++)
            {
               blocksize = mosekblocksizes[i];
               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &solmatrices[i], blocksize) );
            }

            MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, solvectorprimal) );
            for (i = 0; i < nprimalmatrixvars; i++)
               MOSEK_CALL( MSK_getbarxj(sdpisolver->msktask, MSK_SOL_ITR, i, solmatrices[i]) );

#ifdef SCIP_PRINT_SOLU
            /* print primal and dual solution reported from MOSEK */
            SCIPdebugMessage("Dual solution reported from MOSEK transformed to our problem:\n");
            for (i = 0; i < nvars; i++)
               SCIPdebugMessage("y[%d] = %.15g\n", i, solvector[i]);

            SCIPdebugMessage("Primal solution reported from MOSEK:\n");
            for (i = 0; i < nprimalvars; i++)
               SCIPdebugMessage("x[%d] = %.15g\n", i, solvectorprimal[i]);

            for (i = 0; i < nprimalmatrixvars; i++)
            {
               blocksize = mosekblocksizes[i];
               for (j = 0; j <  blocksize * (blocksize + 1) / 2; j++)
                  SCIPdebugMessage("X_%d[%d] = %.15g\n", i, j, solmatrices[i][j]);
            }
#endif

            /* check violations reported from Mosek */
            SCIP_CALL( SCIPsdpSolcheckerCheckAndGetViolDual(sdpisolver->bufmem, nvars, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz,
                  sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval,
                  indchanges, nremovedinds, blockindchanges, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind, lpval,
                  solvector, sdpisolver->feastol, sdpisolver->epsilon, &maxabsviolbndsd, &sumabsviolbndsd, &maxabsviolconsd, &sumabsviolconsd,
                  &maxabsviolsdpd, &sumabsviolsdpd, &checkinfeas) );
            if ( checkinfeas )
               assert( infeasible );

            SCIP_CALL( SCIPsdpSolcheckerCheckAndGetViolPrimal(sdpisolver->bufmem, nvars, obj, lb, ub, sdpisolver->inputtomosekmapper, nsdpblocks,
                  sdpblocksizes, sdpnblockvars, sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz,
                  sdpvar, sdprow, sdpcol, sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nlpcons, lpindchanges, lplhs, lprhs, lpnnonz, lpbeg, lpind,
                  lpval, solvectorprimal, solmatrices, sdpisolver->feastol, sdpisolver->epsilon, &maxabsviolbndsp, &sumabsviolbndsp, &maxabsviolconsp,
                  &sumabsviolconsp, &maxabsviolsdpp, &sumabsviolsdpp, &checkinfeas) );
            if ( checkinfeas )
               assert( infeasible );

            /* dviolcon is the maximal violation of the dual variables, since they are the dual of the primal constraints */
            SCIPdebugMessage("Maximal violations for the dual problem: vars: %g (Mosek), %g (SCIP-SDP); cons: %g (Mosek), %g (SCIP-SDP); SDP: %g (Mosek), %g (SCIP-SDP)\n",
               dviolcon, maxabsviolbndsd, dviolvar, maxabsviolconsd, dviolbarvar, maxabsviolsdpd);
            SCIPdebugMessage("Maximal violations for the primal problem: vars: %g (Mosek), %g (SCIP-SDP); cons: %g (Mosek), %g (SCIP-SDP); SDP: %g (Mosek), %g (SCIP-SDP)\n",
               pviolvar, maxabsviolbndsp, pviolcon, maxabsviolconsp, pviolbarvar, maxabsviolsdpp);
            SCIPdebugMessage("Sum of violations for the dual problem: vars: %g (SCIP-SDP); cons: %g (SCIP-SDP); SDP: %g (SCIP-SDP)\n",
               sumabsviolbndsd, sumabsviolconsd, sumabsviolsdpd);
            SCIPdebugMessage("Sum of violations for the primal problem: vars: %g (SCIP-SDP); cons: %g (SCIP-SDP); SDP: %g (SCIP-SDP)\n",
               sumabsviolbndsp, sumabsviolconsp, sumabsviolsdpp);

            for (i = 0; i < nprimalmatrixvars; i++)
               BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solmatrices[i]);
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solmatrices);
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solvectorprimal);
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &solvector);
         }
#endif

         if ( infeasible )
         {
            feastol *= INFEASFEASTOLCHANGE;
            if ( feastol >= INFEASMINFEASTOL )
            {
               SCIPdebugMessage("Solution feasible for Mosek but outside feasibility tolerance, changing Mosek feasibility tolerance to %g.\n", feastol);
               MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_DFEAS, feastol) );
               MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_INTPNT_CO_TOL_INFEAS, feastol) );
               solveagain = TRUE;
            }
         }

         if ( solveagain )
         {
            /* set the time limit */
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
               if ( ! SCIPsdpiSolverIsInfinity(sdpisolver, solvertimelimit) )
               {
                  MOSEK_CALL( MSK_putdouparam(sdpisolver->msktask, MSK_DPAR_OPTIMIZER_MAX_TIME, solvertimelimit) );/*lint !e641*/
               }

               /* solve the problem */
               MOSEK_CALL( MSK_optimizetrm(sdpisolver->msktask, &(sdpisolver->terminationcode)) );/*lint !e641*/
               MOSEK_CALL( MSK_getsolsta(sdpisolver->msktask, MSK_SOL_ITR, &sdpisolver->solstat) );/*lint !e641*/
               MOSEK_CALL( MSK_getdouinf(sdpisolver->msktask, MSK_DINF_OPTIMIZER_TIME, &opttime) );
               sdpisolver->opttime += opttime;

               if ( sdpisolver->sdpinfo )
               {
                  MOSEK_CALL( MSK_optimizersummary(sdpisolver->msktask, MSK_STREAM_LOG) );/*lint !e641*/
                  MOSEK_CALL( MSK_analyzesolution(sdpisolver->msktask, MSK_STREAM_LOG, MSK_SOL_ITR) );/*lint !e641*/
               }

               /* update number of SDP-iterations and -calls */
               ++sdpisolver->nsdpcalls;
               MOSEK_CALL( MSK_getnaintinf(sdpisolver->msktask, "MSK_IINF_INTPNT_ITER", &newiterations) );/*lint !e641*/
               sdpisolver->niterations += newiterations;

               /* possibly repair status */
               if ( sdpisolver->terminationcode == MSK_RES_TRM_STALL || sdpisolver->solstat == MSK_SOL_STA_UNKNOWN )
               {
                  SCIP_Real pobj;
                  SCIP_Real pviolcon;
                  SCIP_Real pviolvar;
                  SCIP_Real pviolbarvar;
                  SCIP_Real dobj;
                  SCIP_Real dviolcon;
                  SCIP_Real dviolvar;
                  SCIP_Real dviolbarvar;

                  MOSEK_CALL( MSK_getsolutioninfo(sdpisolver->msktask, MSK_SOL_ITR, &pobj, &pviolcon, &pviolvar, &pviolbarvar, NULL, NULL,
                        &dobj, &dviolcon, &dviolvar, &dviolbarvar, NULL) );

                  SCIPdebugMessage("Absolute primal violations: constraints: %g, variables: %g, SDP: %g.\n", pviolcon, pviolvar, pviolbarvar);
                  SCIPdebugMessage("Absolute dual violations: constraints: %g, variables: %g, SDP: %g.\n", dviolcon, dviolvar, dviolbarvar);
                  if ( pviolcon <= sdpisolver->feastol && pviolvar <= sdpisolver->feastol && pviolbarvar <= sdpisolver->feastol
                     && dviolcon <= sdpisolver->feastol && dviolvar <= sdpisolver->feastol && dviolbarvar <= sdpisolver->feastol )
                  {
                     if ( REALABS(dobj - pobj) <= sdpisolver->gaptol )
                     {
                        sdpisolver->terminationcode = MSK_RES_OK;
                        sdpisolver->solstat = MSK_SOL_STA_OPTIMAL;
                        SCIPdebugMessage("Detected stalling - repairing termination code and solution status to 'optimal'.\n");
                     }
                  }
               }
            }
         }
         else
         {
            if ( infeasible )
            {
               sdpisolver->solved = FALSE;
               SCIPmessagePrintInfo(sdpisolver->messagehdlr, "MOSEK failed to reach required feasibility tolerance (feastol: %g, gaptol: %g)!\n", feastol, gaptol);
            }
            break;
         }
      }

      if ( sdpisolver->solved )
      {
         /* if using a penalty formulation, check if the solution is feasible for the original problem
          * we should always count it as infeasible if the penalty problem was unbounded */
         if ( penaltyparam >= sdpisolver->epsilon && sdpisolver->solstat == MSK_SOL_STA_PRIM_INFEAS_CER )
         {
            assert( feasorig != NULL );
            *feasorig = FALSE;
            SCIPdebugMessage("Penalty Problem unbounded!\n");
         }
         else if ( penaltyparam >= sdpisolver->epsilon && ! sdpisolver->timelimit && sdpisolver->terminationcode != MSK_RES_TRM_MAX_TIME )
         {
            SCIP_Real* moseksol;
            SCIP_Real trace = 0.0;
            SCIP_Real* x;

            assert( feasorig != NULL );

            /* get the r variable in the dual problem */
            BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->nactivevars + 1);/*lint !e776*/

            MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

            *feasorig = (moseksol[sdpisolver->nactivevars] < sdpisolver->feastol); /*lint !e413*/

            /* only set sdpisolver->feasorig to true if we solved with objective, because only in this case we want to compute
             * the objective value by hand since it is numerically more stable then the result returned by MOSEK */
            if ( withobj )
               sdpisolver->feasorig = *feasorig;

            /* if r > 0 also check the primal bound */
            if ( ! *feasorig && penaltybound != NULL )
            {
               SCIPdebugMessage("Solution not feasible in original problem, r = %g.\n", moseksol[sdpisolver->nactivevars]);

               /* compute Tr(X) */

               /* start with the diagonal entries of the primal semidefinite variables */
               for (b = 0; b < nsdpblocks; b++)
               {
                  if ( blockindchanges[b] > -1 )
                  {
                     SCIP_Real* X; /* the upper triangular entries of matrix X */
                     int size;

                     size = sdpblocksizes[b] - nremovedinds[b];

                     BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &X, size * (size + 1) / 2) );
                     MOSEK_CALL( MSK_getbarxj(sdpisolver->msktask, MSK_SOL_ITR, b - blockindchanges[b], X) );/*lint !e641*/

                     /* iterate over all diagonal entries */
                     for (i = 0; i < size; i++)
                     {
                        /* get index in the lower triangular part */
                        ind = i * (i + 3) / 2;/*lint !e776*/ /*  i*(i+1)/2 + i  */
                        assert( ind < size * (size + 1) / 2 );
                        trace += X[ind];
                     }

                     BMSfreeBufferMemoryArray(sdpisolver->bufmem, &X);
                  }
               }

               /* add primal lp-variables */
               BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &x, nlpvars + sdpisolver->nvarbounds) );

               MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, x) );/*lint !e641*/

               for (i = 0; i < nlpvars; i++)
                  trace += x[i];

               BMSfreeBufferMemoryArrayNull(sdpisolver->bufmem, &x);

               /* if the relative gap is smaller than the tolerance, we return equality */
               if ( (penaltyparam - trace) / penaltyparam < PENALTYBOUNDTOL )/*lint !e414*/
               {
                  assert( penaltybound != NULL );
                  *penaltybound = TRUE;
                  SCIPdebugMessage("Tr(X) = %g == %g = Gamma, penalty formulation not exact, Gamma should be increased or problem is infeasible.\n",
                     trace, penaltyparam);
               }
               else
                  *penaltybound = FALSE;
            }
            BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
         }
      }
   }

   /* free memory */
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &mosekblocksizes);
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &mosekvarbounds);

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
   case MSK_SOL_STA_UNKNOWN:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
      return FALSE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   default:
      SCIPdebugMessage("Unknown return code in SCIPsdpiSolverFeasibilityKnown\n"); /* TODO: add illposed_cer */
      return FALSE;
   }/*lint !e788*/
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
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      *primalfeasible = TRUE;
      *dualfeasible = TRUE;
      break;
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      *primalfeasible = FALSE;
      *dualfeasible = FALSE;
      break;
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      *primalfeasible = FALSE;
      *dualfeasible = FALSE;
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
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

   switch ( sdpisolver->solstat )
   {
   case MSK_SOL_STA_DUAL_INFEAS_CER:
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
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
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
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
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      return TRUE;
   case MSK_SOL_STA_DUAL_INFEAS_CER:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
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

   switch ( sdpisolver->solstat )
   {
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
      break;
   }/*lint !e788*/
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
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      return TRUE;
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
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
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      return TRUE;
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
      break;
   default:
      SCIPdebugMessage("MOSEK does not know about feasibility of solutions\n");
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

   /* check if Mosek stalled when it was already acceptable */
   if ( sdpisolver->terminationcode == MSK_RES_TRM_STALL )
   {
      SCIP_Real pobj;
      SCIP_Real dobj;
      SCIP_Real gapnormalization;

      /* check the solution status */
      switch ( sdpisolver->solstat )
      {
      case MSK_SOL_STA_UNKNOWN:
      case MSK_SOL_STA_PRIM_FEAS:
      case MSK_SOL_STA_DUAL_FEAS:
         return FALSE;
      case MSK_SOL_STA_OPTIMAL:
      case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      case MSK_SOL_STA_PRIM_INFEAS_CER:
      case MSK_SOL_STA_DUAL_INFEAS_CER:
         /* check duality gap */
         MOSEK_CALL_BOOL( MSK_getdualobj(sdpisolver->msktask, MSK_SOL_ITR, &dobj) );
         MOSEK_CALL_BOOL( MSK_getprimalobj(sdpisolver->msktask, MSK_SOL_ITR, &pobj) );
         /* for the relative gap we divide by max(1.0, min(pobj, dobj)), as this is also done in Mosek */
         gapnormalization = dobj > pobj ? (pobj > 1.0 ? pobj : 1.0) : (dobj > 1.0 ? dobj : 1.0);
         if ( REALABS((pobj-dobj) / gapnormalization) < sdpisolver->gaptol )
            return TRUE;
         else
            return FALSE;
      default:
         return FALSE;
      }
   }

   return sdpisolver->terminationcode == MSK_RES_OK;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   return sdpisolver->terminationcode == MSK_RES_TRM_OBJECTIVE_RANGE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED_BOOL( sdpisolver );

   return sdpisolver->terminationcode == MSK_RES_TRM_MAX_ITERATIONS;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );

   if ( sdpisolver->timelimit )
      return TRUE;

   if ( ! sdpisolver->solved )
      return FALSE;

   return sdpisolver->terminationcode == MSK_RES_TRM_MAX_TIME;
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

   switch ( sdpisolver->terminationcode )
   {
   case MSK_RES_OK:
      return 0;
   case MSK_RES_TRM_MAX_NUM_SETBACKS:
   case MSK_RES_TRM_NUMERICAL_PROBLEM:
   case MSK_RES_TRM_STALL:
      return 2;
   case MSK_RES_TRM_OBJECTIVE_RANGE:
      return 3;
   case MSK_RES_TRM_MAX_ITERATIONS:
      return 4;
   case MSK_RES_TRM_MAX_TIME:
      return 5;
   default:
      return 7;
   }/*lint !e788*/
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

   if ( sdpisolver->terminationcode != MSK_RES_OK )
      return FALSE;

   if ( sdpisolver->solstat != MSK_SOL_STA_OPTIMAL )
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
   SCIP_Real* moseksol;

   assert( sdpisolver != NULL );
   CHECK_IF_SOLVED( sdpisolver );
   assert( objval != NULL );

   /* check for unboundedness */
   if ( SCIPsdpiSolverIsDualUnbounded(sdpisolver) || SCIPsdpiSolverIsPrimalInfeasible(sdpisolver) )
   {
      *objval = -SCIPsdpiSolverInfinity(sdpisolver);
      return SCIP_OKAY;
   }

   if ( sdpisolver->penalty && ! sdpisolver->feasorig )
   {
      /* in this case we cannot really trust the solution given by MOSEK, since changes in the value of r much less than epsilon can
       * cause huge changes in the objective, so using the objective value given by MOSEK is numerically more stable */
      MOSEK_CALL( MSK_getdualobj(sdpisolver->msktask, MSK_SOL_ITR, objval) );

      /* reverse scaling */
      *objval *= sdpisolver->objscalefactor;
   }
   else
   {
      int v;

      /* since the objective value given by MOSEK sometimes differs slightly from the correct value for the given solution,
       * we get the solution from MOSEK and compute the correct objective value */
      BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->penalty ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars);
      MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

      *objval = 0.0;
      for (v = 0; v < sdpisolver->nactivevars; v++)
         *objval += moseksol[v] * sdpisolver->objcoefs[v];

      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
   }

   /* as we didn't add the fixed (lb = ub) variables to MOSEK, we have to add their contributions to the objective as well */
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
      SCIP_Real* moseksol;
      int v;

      BMSallocBufferMemoryArray(sdpisolver->bufmem, &moseksol, sdpisolver->penalty ? sdpisolver->nactivevars + 1 : sdpisolver->nactivevars);

      MOSEK_CALL( MSK_gety(sdpisolver->msktask, MSK_SOL_ITR, moseksol) );/*lint !e641*/

      /* insert the entries into dualsol, for non-fixed vars we copy those from MOSEK, the rest are the saved entries from inserting (they equal lb=ub) */
      for (v = 0; v < sdpisolver->nvars; v++)
      {
         if ( sdpisolver->inputtomosekmapper[v] >= 0 )
            dualsol[v] = moseksol[sdpisolver->inputtomosekmapper[v]];
         else
         {
            /* this is the value that was saved when inserting, as this variable has lb=ub */
            assert( -sdpisolver->inputtomosekmapper[v] <= sdpisolver->nvars - sdpisolver->nactivevars );
            dualsol[v] = sdpisolver->fixedvarsval[(-1 * sdpisolver->inputtomosekmapper[v]) - 1]; /*lint !e679*/ /* -1 because we wanted strictly negative vals */
         }
      }

      /* if both solution and objective should be printed, we can use the solution to compute the objective */
      if ( objval != NULL )
      {
         if ( sdpisolver->penalty && ! sdpisolver->feasorig )
         {
            /* in this case we cannot really trust the solution given by MOSEK, since changes in the value of r much less than epsilon can
             * cause huge changes in the objective, so using the objective value given by MOSEK is numerically more stable */
            MOSEK_CALL( MSK_getdualobj(sdpisolver->msktask, MSK_SOL_ITR, objval) );

            /* reverse scaling */
            *objval *= sdpisolver->objscalefactor;
         }
         else
         {
            /* since the objective value given by MOSEK sometimes differs slightly from the correct value for the given solution,
             * we get the solution from MOSEK and compute the correct objective value */
            *objval = 0.0;
            for (v = 0; v < sdpisolver->nactivevars; v++)
               *objval += moseksol[v] * sdpisolver->objcoefs[v];
         }

         /* as we didn't add the fixed (lb = ub) variables to MOSEK, we have to add their contributions to the objective as well */
         *objval += sdpisolver->fixedvarsobjcontr;
      }

      BMSfreeBufferMemoryArray(sdpisolver->bufmem, &moseksol);
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
   SCIPdebugMessage("Not implemented yet\n");

   return SCIP_PLUGINNOTFOUND;
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
   return SCIP_LPERROR;
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
   SCIP_Real* primalvals;
   int nprimalvars;
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

   /* get primal solution from MOSEK */
   MOSEK_CALL( MSK_getnumvar(sdpisolver->msktask, &nprimalvars) );/*lint !e641*/
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &primalvals, nprimalvars) );
   MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, primalvals) );/*lint !e641*/

   /* iterate over all variable bounds and insert the corresponding primal variables in the right positions of the return-arrays */
   assert( sdpisolver->nvarbounds <= 2 * sdpisolver->nvars );

   for (i = 0; i < sdpisolver->nvarbounds; i++)
   {
      /* this is a lower bound */
      if ( sdpisolver->varboundpos[i] < 0 )
      {
         /* the last nvarbounds entries correspond to the varbounds; we need to unscale these values */
         lbvals[sdpisolver->mosektoinputmapper[- sdpisolver->varboundpos[i] -1]] = primalvals[nprimalvars - sdpisolver->nvarbounds + i] * sdpisolver->objscalefactor;
      }
      else
      {  /* this is an upper bound */
         assert( sdpisolver->varboundpos[i] > 0 );

         /* the last nvarbounds entries correspond to the varbounds; we need to unscale these values */
         ubvals[sdpisolver->mosektoinputmapper[sdpisolver->varboundpos[i] - 1]] = primalvals[nprimalvars - sdpisolver->nvarbounds + i] * sdpisolver->objscalefactor;
      }
   }

   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &primalvals);

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
   SCIP_Real* primalvars;
   int nprimalvars;
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

   /* get primal solution from Mosek */
   MOSEK_CALL( MSK_getnumvar(sdpisolver->msktask, &nprimalvars) );/*lint !e641*/
   BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &primalvars, nprimalvars) );
   MOSEK_CALL( MSK_getxx(sdpisolver->msktask, MSK_SOL_ITR, primalvars) );/*lint !e641*/

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
         lhsvals[i] = primalvars[ind] * sdpisolver->objscalefactor;
         ++ind;
      }
      else
         lhsvals[i] = 0.0;

      if ( lprhs[i] < SCIPsdpiSolverInfinity(sdpisolver) )
      {
         rhsvals[i] = primalvars[ind] * sdpisolver->objscalefactor;
         ++ind;
      }
      else
         rhsvals[i] = 0.0;

      assert( ind <= nprimalvars );
   }
   BMSfreeBufferMemoryArray(sdpisolver->bufmem, &primalvars);

   return SCIP_OKAY;
}

/** return number of nonzeros for each block of the primal solution matrix X (including lp block) */
SCIP_RETCODE SCIPsdpiSolverGetPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
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
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** returns the primal solution matrix (without LP rows)
 *
 *  The solution of Mosek is given as an array in lower triangular column-wise stacked form. The position of entry (i,j)
 *  (with i >= j) in this array is given by \sum_{k=0}^{j-1} (n - k) + (i - j), because:
 *  - The number of entries before column j is given by the first sum (from first column: 0, n, n + (n-1), ...).
 *  - Then we add the number of entries in column j, which is i - j.
 */
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
   int b;

   assert( sdpisolver != NULL );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( indchanges != NULL );
   assert( nremovedinds != NULL );
   assert( blockindchanges != NULL );
   assert( primalmatrices != NULL );

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
         int redsize;
         int row;
         int col;
         int idx;
         int i;

         redsize = blocksize - nremovedinds[b];

         BMS_CALL( BMSallocBufferMemoryArray(sdpisolver->bufmem, &X, redsize * (redsize + 1) / 2) );
         MOSEK_CALL( MSK_getbarxj(sdpisolver->msktask, MSK_SOL_ITR, b - blockindchanges[b], X) );/*lint !e641*/

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

                     if ( row >= col )
                        idx = (redsize * col) - (col - 1) * col/2 + (row - col);
                     else
                        idx = (redsize * row) - (row - 1) * row/2 + (col - row);
                     assert( 0 <= idx && idx < redsize * (redsize + 1)/2 );

                     val = X[idx];
                     primalmatrices[b][i * blocksize + j] = val;
                     primalmatrices[b][j * blocksize + i] = val;
                  }
               }
            }
         }

         BMSfreeBufferMemoryArray(sdpisolver->bufmem, &X);
      }
   }

   return SCIP_OKAY;
}

/** return the maximum absolute value of the optimal primal matrix */
SCIP_Real SCIPsdpiSolverGetMaxPrimalEntry(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
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
   return 1.0e16;
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
   SCIPdebugMessage("Lambdastar parameter not used by MOSEK"); /* this parameter is only used by SDPA */

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
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   assert( sdpisolver != NULL );
   assert( fname != NULL );

   MOSEK_CALL( MSK_writedata(sdpisolver->msktask, fname) );/*lint !e641*/

   return SCIP_OKAY;
}

/**@} */
