/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2017 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2017 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_sdp.c
 * @ingroup RELAXATORS
 * @brief  SDP-relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG   *//* displays complete solution for each relaxation */
/* #define SCIP_EVEN_MORE_DEBUG  *//* shows number of deleted empty cols/rows for every relaxation and variable status &
 * bounds as well as all constraints in the beginning */
/* #define SCIP_PRINT_WARMSTART  *//* print initial point given for warmstarts */
#define SLATERSOLVED_ABSOLUTE /* uncomment this to return the absolute number of nodes for, e.g., solved fast with slater in addition to percentages */

#include "relax_sdp.h"

#include "assert.h"                     /*lint !e451*/
#include "string.h"                     /* for strcmp */
#include <sys/time.h>                   /* for timeofday */

#include "SdpVarmapper.h"
#include "SdpVarfixer.h"
#include "sdpi/sdpi.h"
#include "sdpi/lapack.h"
#include "scipsdp/cons_sdp.h"
#include "scipsdp/cons_savesdpsol.h"
#include "scipsdp/cons_savedsdpsettings.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define RELAX_NAME                  "SDP"
#define RELAX_DESC                  "SDP-relaxator"
#define RELAX_PRIORITY              1
#define RELAX_FREQ                  1

/* default values for parameters: */
#define DEFAULT_PENALTYPARAM        -1.0     /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#define DEFAULT_LAMBDASTAR          -1.0     /**< the parameter lambda star used by SDPA to set the initial point */
#define DEFAULT_MAXPENALTYPARAM     -1.0     /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#define DEFAULT_WARMSTARTIPFACTOR   0.50     /**< factor for interior point in convexcombination of IP and parent solution, if warmstarts are enabled */
#define DEFAULT_WARMSTARTPRIMALTYPE 3        /**< how to warmstart the primal problem? 1: scaled identity, 2: elementwise reciprocal, 3: saved primal sol */
#define DEFAULT_WARMSTARTIPTYPE     1        /**< which interior point to use for convex combination for warmstarts? 1: scaled identity, 2: analytic center */
#define DEFAULT_WARMSTARTPROJECT    2        /**< how to update dual matrix for new bounds? 1: use old bounds, 2: use new bounds, 3: use new bounds and project on psd cone, 4: use new bounds and solve rounding problem */
#define DEFAULT_WARMSTARTPROJMINEV  -1       /**< minimum eigenvector to allow when projecting onto the positive (semi-)definite cone */
#define DEFAULT_WARMSTARTPROJPDSAME TRUE     /**< Should one shared minimum eigenvalue be computed for primal and dual problem instead of different ones if warmstartpmevpar = -1 ? */
#define DEFAULT_WARMSTART_PREOPTIMAL_SOL FALSE /**< Should a preoptimal solution (with higher epsilon) instead of the optimal solution be used for warmstarts (currently only implemented fo DSDP) */
#define DEFAULT_WARMSTARTPREOPTGAP  1e-2     /**< If warmstartpreoptimalsol is TRUE, this is the gap where the preoptimal solution is saved (currently only implemented fo DSDP) */
#define DEFAULT_WARMSTARTROUNDONLYINF FALSE  /**< Only use solution of roundingproblem to detect infeasibility (only has an effect for warmstartproject = 4) */
#define DEFAULT_SLATERCHECK         0        /**< Should the Slater condition be checked ? */
#define DEFAULT_OBJLIMIT            FALSE    /**< Should an objective limit be given to the SDP-Solver ? */
#define DEFAULT_RESOLVE             TRUE     /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
#define DEFAULT_TIGHTENVB           TRUE     /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
#define DEFAULT_SDPINFO             FALSE    /**< Should the SDP solver output information to the screen? */
#define DEFAULT_WARMSTART           FALSE    /**< Should the SDP solver try to use warmstarts? */
#define DEFAULT_DISPLAYSTAT         FALSE    /**< Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ? */
#define DEFAULT_SETTINGSRESETFREQ   -1       /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
#define DEFAULT_SETTINGSRESETOFS    0        /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */
#define DEFAULT_SDPSOLVERTHREADS    -1       /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */

#define WARMSTART_MINVAL            0.01     /**< if we get a value less than this when warmstarting (currently only for the linear part when combining with analytic center), the value is set to this */
#define WARMSTART_PROJ_MINRHSOBJ    1        /**< minimum value for rhs/obj when computing minimum eigenvalue for warmstart-projection */
#define WARMSTART_PROJ_FACTOR       0.1      /**< factor to multiply maximum rhs/obj/coef with when computing minimum eigenvalue for warmstart-projection */
#define WARMSTART_PROJ_FACTOR_LHS   10       /**< factor to multiply maximum SDP coefficient with before applying WARMSTART_PROJ_FACTOr (to account for summation of lhs entries) */
#define WARMSTART_PROJ_FACTOR_PRIMAL 0.1     /**< factor to multiply maximum obj with when computing minimum eigenvalue for warmstart-projection in the primal */
#define WARMSTART_PROJ_FACTOR_DUAL  0.1      /**< factor to multiply maximum rhs with when computing minimum eigenvalue for warmstart-projection in the dual */
#define WARMSTART_PREOPT_MIN_Z_LPVAL 0.01    /**< minimal (diagonal) entry for LP block of dual matrix for preoptimal warmstarts */

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

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_SDPI*            sdpi;               /**< general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SCIP_LPI*             lpi;                /**< LP interface; used for rounding problems */
   SdpVarmapper*         varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */

   SCIP_Real             objval;             /**< objective value of the last SDP-relaxation */
   SCIP_Bool             origsolved;         /**< solved original problem to optimality (not only a penalty or probing formulation) */
   SCIP_Bool             probingsolved;      /**< was the last probing SDP solved successfully? */
   long int              lastsdpnode;        /**< number of the SCIP node the current SDP-solution belongs to */
   SCIP_Bool             feasible;           /**< was the last solved SDP feasible */

   SCIP_Real             sdpsolvergaptol;    /**< the stopping criterion for the duality gap the sdpsolver should use */
   SCIP_Real             sdpsolverfeastol;   /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
   SCIP_Real             penaltyparam;       /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   SCIP_Real             maxpenaltyparam;    /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   SCIP_Real             lambdastar;         /**< the parameter lambda star used by SDPA to set the initial point */
   SCIP_Real             computedlambdastar; /**< computed value for lambda star parameter used by SDPA to set the initial point */
   int                   npenaltyincr;       /**< maximum number of times the penalty parameter will be increased if penalty formulation failed */
   int                   slatercheck;        /**< Should the Slater condition for the dual problem be check ahead of solving every SDP ? */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             displaystat;        /**< Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ? */
   SCIP_Bool             objlimit;           /**< Should an objective limit be given to the SDP solver? */
   SCIP_Bool             resolve;            /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
   SCIP_Bool             tightenvb;          /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
   int                   settingsresetfreq;  /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
   int                   settingsresetofs;   /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */
   int                   sdpsolverthreads;   /**< number of threads the SDP solver should use, currently only supported for MOSEK (-1 = number of cores) */

   int                   sdpcalls;           /**< number of solved SDPs (used to compute average SDP iterations), different settings tried are counted as multiple calls */
   int                   sdpinterfacecalls;  /**< number of times the SDP interfaces was called (used to compute slater statistics) */
   int                   sdpiterations;      /**< saves the total number of sdp-iterations */
   int                   solvedfast;         /**< number of SDPs solved with fast settings */
   int                   solvedmedium;       /**< number of SDPs solved with medium settings */
   int                   solvedstable;       /**< number of SDPs solved with stable settings */
   int                   solvedpenalty;      /**< number of SDPs solved using penalty formulation */
   int                   unsolved;           /**< number of SDPs that could not be solved even using a penalty formulation */
   int                   stablewslater;      /**< number of instances solved with fastest settings where primal and dual slater held */
   int                   unstablewslater;    /**< number of instances solved with stable settings where primal and dual slater held */
   int                   penaltywslater;     /**< number of instances solved with penalty formulation where primal and dual slater held */
   int                   boundedwslater;     /**< number of instances we could compute a bound for via the penalty approach where primal and dual slater held */
   int                   unsolvedwslater;    /**< number of instances that could not be solved where primal and dual slater held */
   int                   stablenoslater;     /**< number of instances solved with fastest setting where either primal or dual slater did not hold */
   int                   unstablenoslater;   /**< number of instances solved with stable settings where either primal or dual slater did not hold */
   int                   penaltynoslater;    /**< number of instances solved with penalty formulation where either primal or dual slater did not hold */
   int                   boundednoslater;    /**< number of instances we could compute a bound for via the penalty approach where either primal or dual slater did not hold */
   int                   unsolvednoslater;   /**< number of instances that could not be solved where either primal or dual slater did not hold */
   int                   nslaterholds;       /**< number of SDPs for which primal and dual slater condition held */
   int                   nnoslater;          /**< number of SDPs for which either primal or dual slater condition did not hold (including those where we could not check the other) */
   int                   nslatercheckfailed; /**< number of SDPs for which we failed to check the slater condition (this only includes SDPs where both checks failed or
                                              *   one checked returned slater holds and the other failed but not those where the first check already returned that it does not hold */
   int                   npslaterholds;      /**< number of SDPs for which primal slater condition held */
   int                   npnoslater;         /**< number of SDPs for which primal slater condition did not hold */
   int                   npslatercheckfailed;/**< number of SDPs for which we failed to check the dual slater condition */
   int                   ndslaterholds;      /**< number of SDPs for which dual slater condition held */
   int                   ndnoslater;         /**< number of SDPs for which dual slater condition did not hold */
   int                   ndslatercheckfailed;/**< number of SDPs for which we failed to check the dual slater condition */
   int                   nslaterinfeasible;  /**< number of SDPs for which we detected infeasibility during the Slater check */
   int                   stableinfeasible;   /**< number of instances solved with fastest settings where the dual slater check showed that the problem is infeasible */
   int                   unstableinfeasible; /**< number of instances solved with stable settings where the dual slater check showed that the problem is infeasible */
   int                   penaltyinfeasible;  /**< number of instances solved with penalty formulation where the dual slater check showed that the problem is infeasible */
   int                   boundedinfeasible;  /**< number of instances we could compute a bound for via the penalty approach where the dual slater check showed that the problem is infeasible */
   int                   unsolvedinfeasible; /**< number of instances that could not be solved where the dual slater check showed that the problem is infeasible */
   int                   roundingprobinf;    /**< number of instances that where detected infeasible through the primal rounding problem */
   int                   primalroundfails;   /**< number of instances where the primal rounding problem failed */
   int                   dualroundfails;     /**< number of instances where the dual rounding problem failed */
   int                   roundstartsuccess;  /**< number of instances that could be warmstarted using the solution of the rounding problems */
   int                   roundingoptimal;    /**< number of instances where the optimal solution was found by the rounding problem */
   int                   roundingcutoff;     /**< number of instances that could be cut off through bounding by the rounding problem */
   SCIP_Real             roundingprobtime;   /**< total time spent in rouding problems for warmstarting/infeasibility detection */

   SCIP_Bool             warmstart;          /**< Should the SDP solver try to use warmstarts? */
   SCIP_Real             warmstartipfactor;  /**< factor for interior point in convexcombination of IP and parent solution, if warmstarts are enabled */
   int                   warmstartprimaltype;/**< how to warmstart the primal problem? 1: scaled identity/analytic center, 2: elementwise reciprocal, 3: saved primal sol
                                               *  TODO: should probably remove elementwise reciprocal, since this doesn't work from a theoretical point of view*/
   int                   warmstartproject;   /**< how to update dual matrix for new bounds? 1: use old bounds, 2: use new bounds, 3: use new bounds and project on psd cone, 4: use new bounds and solve rounding problem */
   SCIP_Real             warmstartpmevprimalpar; /**< SCIP parameter for min eigenvalue when projecting primal onto positive definite cone; -1 for automatic computation */
   SCIP_Real             warmstartpmevdualpar; /**< SCIP parameter for min eigenvalue when projecting dual onto positive definite cone; -1 for automatic computation */
   SCIP_Real             warmstartprojminevprimal; /**< minimum eigenvalue to allow when projecting onto the positive (semi-)definite cone in the primal */
   SCIP_Real             warmstartprojminevdual; /**< minimum eigenvalue to allow when projecting onto the positive (semi-)definite cone in the dual */
   SCIP_Bool             warmstartprojpdsame;/**< Should one shared minimum eigenvalue respectively maximum entry be computed for primal and dual problem instead of different ones for primal and dual and each block for projection or convex combination ? */
   int                   warmstartiptype;    /**< which interior point to use for convex combination for warmstarts? 1: scaled identity, 2: analytic center */
   SCIP_Bool             warmstartpreoptsol; /**< Should a preoptimal solution (with higher epsilon) instead of the optimal solution be used for warmstarts (currently only implemented fo DSDP) */
   SCIP_Real             warmstartpreoptgap; /**< In case a preoptimal solution should be used for warmstarts, this gives the gap where the solution should be saved (currently only implemented fo DSDP) */
   SCIP_Bool             warmstartroundonlyinf; /**< Only use solution of roundingproblem to detect infeasibility (only has an effect for warmstartproject = 4) */
   int                   nblocks;            /**< number of blocks INCLUDING lp-block */
   SCIP_Bool             ipXexists;          /**< has an interior point for primal matrix X been successfully computed */
   int*                  ipXnblocknonz;      /**< interior point for primal matrix X for convex combination for warmstarts: number of nonzeros for each block
                                              *   if computation of analytic center failed, first entry will be -1 */
   int**                 ipXrow;             /**< interior point for primal matrix X for convex combination for warmstarts: row indices */
   int**                 ipXcol;             /**< interior point for primal matrix X for convex combination for warmstarts: column indices */
   SCIP_Real**           ipXval;             /**< interior point for primal matrix X for convex combination for warmstarts: values */
   SCIP_Bool             ipZexists;          /**< has an interior point for dual matrix Z (and corresponding vector y) been successfully computed */
   SCIP_SOL*             ipy;                /**< interior point for dual vector y for convex combination for warmstarts */
   int*                  ipZnblocknonz;      /**< interior point for dual matrix Z for convex combination for warmstarts: number of nonzeros for each block
                                               *   if computation of analytic center failed, first entry will be -1 */
   int**                 ipZrow;             /**< interior point for dual matrix Z for convex combination for warmstarts: row indices */
   int**                 ipZcol;             /**< interior point for dual matrix Z for convex combination for warmstarts: column indices */
   SCIP_Real**           ipZval;             /**< interior point for dual matrix Z for convex combination for warmstarts: values */
};

/** expand sparse matrix to full matrix format needed by LAPACK */
static
SCIP_RETCODE expandSparseMatrix(
   int                   nnonz,              /**< number of nonzeros and length of row/col/val arrays */
   int                   blocksize,          /**< size of matrix (and squareroot of memory allocated for fullmat) */
   int*                  row,                /**< row indices */
   int*                  col,                /**< column indices */
   SCIP_Real*            val,                /**< values */
   SCIP_Real*            fullmat             /**< pointer to store full matrix */
   )
{
   int i;
   int matrixsize;

   assert( nnonz >= 0 );
   assert( row != NULL );
   assert( col != NULL );
   assert( val != NULL );
   assert( fullmat != NULL );

   matrixsize = blocksize * blocksize;

   /* initialize matrix with zeros */
   for (i = 0; i < matrixsize; i++)
      fullmat[i] = 0.0;

   for (i = 0; i < nnonz; i++)
   {
      assert( row[i] * blocksize + col[i] <= matrixsize );
      fullmat[row[i] * blocksize + col[i]] = val[i];
      assert( col[i] * blocksize + row[i] <= matrixsize );
      fullmat[col[i] * blocksize + row[i]] = val[i];
   }

   return SCIP_OKAY;
}

/** multiplies all entries in the i-th column by scale[i] */
static
SCIP_RETCODE scaleTransposedMatrix(
   int                   blocksize,          /* number of rows and columns */
   SCIP_Real*            matrix,             /* matrix entries given as blocksize^2 array */
   SCIP_Real*            scale               /* array of length blocksize to multiply the columns of matrix with */
   )
{
   int r;
   int c;

   assert( blocksize >= 0 );
   assert( matrix != NULL );
   assert( scale != NULL );

   for (r = 0; r < blocksize; r++)
   {
      for (c = 0; c < blocksize; c++)
      {
         matrix[r * blocksize + c] *= scale[c];
      }
   }

   return SCIP_OKAY;
}

/** inserts all the SDP data into the corresponding SDP Interface */
static
SCIP_RETCODE putSdpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper,          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Bool             primalobj,          /**< should the primal objective coefficients (constant part of the SDP constraint) be used ? */
   SCIP_Bool             boundprimal         /**< should the primal problem be bounded (Tr(X)<=1) through a penalty term in the dual ? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_CONS** conss;
   SCIP_VAR** blockvars;
   SCIP_VAR** vars;
   SCIP_Real*** val;
   SCIP_Real** constval;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real param;
   int*** row;
   int*** col;
   int** nblockvarnonz;
   int** constrow;
   int** constcol;
   int** sdpvar;
   int* sdpblocksizes;
   int* nblockvars;
   int* nconstblocknonz;
   int constnnonzcounter;
   int blocknnonz;
   int sdpconstnnonz;
   int sdpnnonz;
   int nsdpblocks;
   int constlength;
   int nvars;
   int nvarspen;
   int nconss;
   int ind;
   int i;
   int j;

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolvergaptol", &param) );

   SCIPdebugMessage("Putting SDP Data in general SDP interface!\n");

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   nvarspen = boundprimal ? nvars + 1 : nvars; /* if the primal should be bounded, an additional penalty variable is added to the dual */

   /* prepare arrays of objective values and bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvarspen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvarspen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvarspen) );

   for (i = 0; i < nvars; i++)
   {
      obj[i] = SCIPvarGetObj(vars[i]);
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
   }
   if ( boundprimal )
   {
      obj[nvars] = 1.0; /* this objective coefficient together with lb = 0 and the identity matrix leads to constraint Tr(X) <= 1 */
      lb[nvars] = 0.0;
      ub[nvars] = SCIPinfinity(scip);
   }

   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);

   /* count the number of sdpblocks and compute the number of nonzeros */
   nsdpblocks = 0;
   sdpnnonz = 0;
   sdpconstnnonz = 0;

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

#ifdef SCIP_EVEN_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      if ( strcmp(conshdlrname, "SDP") == 0 )
      {
         nsdpblocks++;

         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[i], &blocknnonz, &constnnonzcounter) );
         sdpnnonz += blocknnonz;
         sdpconstnnonz += constnnonzcounter;
      }
   }

   /* create the sdp- and sdpconst-arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpblocksizes, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nblockvarnonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nconstblocknonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &col, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constcol, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constrow, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &constval, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nblockvars, nsdpblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpvar, nsdpblocks) );

   for (i = 0; i < nsdpblocks; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(nblockvarnonz[i]), nvarspen) );
      SCIP_CALL( SCIPallocBufferArray(scip, &col[i], nvarspen) );
      SCIP_CALL( SCIPallocBufferArray(scip, &row[i], nvarspen) );
      SCIP_CALL( SCIPallocBufferArray(scip, &val[i], nvarspen) );
   }

   /* get the SDP-data */
   ind = 0; /* index of the current sdp block in the complete sdp */
   SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nvars) );

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(conshdlrname, "SDP") == 0 )
      {
         assert( ind < nsdpblocks );

         /* allocate memory for the constant nonzeros */
         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[i], NULL, &constlength) );
         nconstblocknonz[ind] = constlength;
         SCIP_CALL( SCIPallocBufferArray(scip, &(constcol[ind]), constlength) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(constrow[ind]), constlength) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(constval[ind]), constlength) );

         /* get the data */
         SCIP_CALL( SCIPconsSdpGetData(scip, conss[i], &nblockvars[ind], &blocknnonz, &sdpblocksizes[ind], &nvars, nblockvarnonz[ind], col[ind],
            row[ind], val[ind], blockvars, &nconstblocknonz[ind], constcol[ind], constrow[ind], constval[ind]) );

         /* nvars and nconstblocknonz[ind] would have been overwritten if the space in the given arrays hadn't been sufficient */
         assert( nvars == SCIPgetNVars(scip) );
         assert( nconstblocknonz[ind] <= constlength );

         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar[ind]), boundprimal ? nblockvars[ind] + 1 : nblockvars[ind]) );

         /* get global variable indices */
         for (j = 0; j < nblockvars[ind]; j++)
            sdpvar[ind][j] = SCIPsdpVarmapperGetSdpIndex(varmapper, blockvars[j]);

         if ( boundprimal )
         {
            /* penalty variable is added as final variable to bound the primal */
            sdpvar[ind][nblockvars[ind]] = SCIPsdpVarmapperGetNVars(varmapper);
            nblockvarnonz[ind][nblockvars[ind]] = sdpblocksizes[ind];

            /* add identity matrix times penalty variable */
            SCIP_CALL( SCIPallocBufferArray(scip, &col[ind][nblockvars[ind]], sdpblocksizes[ind]) );
            SCIP_CALL( SCIPallocBufferArray(scip, &row[ind][nblockvars[ind]], sdpblocksizes[ind]) );
            SCIP_CALL( SCIPallocBufferArray(scip, &val[ind][nblockvars[ind]], sdpblocksizes[ind]) );

            for (j = 0; j < sdpblocksizes[ind]; j++)
            {
               col[ind][nblockvars[ind]][j] = j;
               row[ind][nblockvars[ind]][j] = j;
               val[ind][nblockvars[ind]][j] = 1.0;
            }
            nblockvars[ind]++;
         }

         ind++;
      }
   }

   /* free the memory that is no longer needed */
   SCIPfreeBufferArray(scip, &blockvars);

   /* load data into SDPI */
   if ( primalobj )
   {
      SCIP_CALL( SCIPsdpiLoadSDP(sdpi, nvarspen,  obj, lb, ub, nsdpblocks, sdpblocksizes, nblockvars, sdpconstnnonz, nconstblocknonz, constrow,
                               constcol, constval, sdpnnonz, nblockvarnonz, sdpvar, row, col, val, 0,
                               NULL, NULL, 0, NULL, NULL, NULL) ); /* insert the SDP part, add an empty LP part */
   }
   else
   {
      /* overwrite nconstblocknonz */
      for (i = 0; i < nsdpblocks; i++)
         nconstblocknonz[i] = 0;

      SCIP_CALL( SCIPsdpiLoadSDP(sdpi, nvarspen,  obj, lb, ub, nsdpblocks, sdpblocksizes, nblockvars, 0, nconstblocknonz, NULL,
                               NULL, NULL, sdpnnonz, nblockvarnonz, sdpvar, row, col,  val, 0,
                               NULL, NULL, 0, NULL, NULL, NULL) ); /* insert the SDP part, add an empty LP part */
   }


   /* free the remaining memory */
   for (i = 0; i < nsdpblocks; i++)
   {
      if ( boundprimal )
      {
         SCIPfreeBufferArrayNull(scip, &val[i][nblockvars[i] - 1]);
         SCIPfreeBufferArrayNull(scip, &row[i][nblockvars[i] - 1]);
         SCIPfreeBufferArrayNull(scip, &col[i][nblockvars[i] - 1]);
      }
      SCIPfreeBufferArrayNull(scip, &(sdpvar[i]));
      SCIPfreeBufferArrayNull(scip, &val[i]);
      SCIPfreeBufferArrayNull(scip, &row[i]);
      SCIPfreeBufferArrayNull(scip, &col[i]);
      SCIPfreeBufferArrayNull(scip, &(nblockvarnonz[i]));
      SCIPfreeBufferArrayNull(scip, &(constval[i]));
      SCIPfreeBufferArrayNull(scip, &(constrow[i]));
      SCIPfreeBufferArrayNull(scip, &(constcol[i]));
   }

   SCIPfreeBufferArrayNull(scip, &sdpvar);
   SCIPfreeBufferArrayNull(scip, &nblockvars);
   SCIPfreeBufferArrayNull(scip, &constval);
   SCIPfreeBufferArrayNull(scip, &constrow);
   SCIPfreeBufferArrayNull(scip, &constcol);
   SCIPfreeBufferArrayNull(scip, &val);
   SCIPfreeBufferArrayNull(scip, &row);
   SCIPfreeBufferArrayNull(scip, &col);
   SCIPfreeBufferArrayNull(scip, &nconstblocknonz);
   SCIPfreeBufferArrayNull(scip, &nblockvarnonz);
   SCIPfreeBufferArrayNull(scip, &sdpblocksizes);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);
   SCIPfreeBufferArray(scip, &obj);

   return SCIP_OKAY;
}

/** inserts all the LP data (including bounds and objective) into the corresponding SDP Interface */
static
SCIP_RETCODE putLpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper,          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Bool             primalobj,          /**< should the primal objective coefficients (lhs/rhs of LP-constraints) be used ? */
   SCIP_Bool             dualobj             /**< should the dual objective coefficients be used ? */
   )
{
   SCIP_VAR** vars;
   SCIP_COL** rowcols;
   SCIP_ROW** rows;
   SCIP_Bool tightenvb;
   SCIP_Real* rowvals;
   SCIP_Real* lhs;
   SCIP_Real* rhs;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   SCIP_Real* val;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   int* inds;
   int* objinds;
   int* rowind;
   int* colind;
   int nrowssdpi;
   int nrows;
   int rownnonz;
   int nvars;
   int nconss;
   int scipnnonz;
   int nnonz;
   int i;
   int j;

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/tightenvb", &tightenvb) );

   SCIPdebugMessage("inserting %d LPRows into the interface.\n", nrows);

   /* compute the total number of LP nonzeroes in SCIP */
   scipnnonz = 0;
   for (i = 0; i < nrows; i++)
   {
      assert( rows[i] != NULL );
      scipnnonz += SCIProwGetNNonz(rows[i]);
   }

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, scipnnonz) );

   /* insert the nonzeroes */
   nnonz = 0; /* this is recomputed for the sdpi, because of the possible duplication of non-zeroes for lhs and rhs */
   nconss = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      SCIP_ROW* row;
      SCIP_Bool tightened = FALSE;
      SCIP_Real tightenedval = 0.0;
      SCIP_Bool swapped = FALSE;

      row = rows[i];
      assert( row != 0 );
      rownnonz = SCIProwGetNNonz(row);

      rowvals = SCIProwGetVals(row);
      rowcols = SCIProwGetCols(row);
      sciplhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      sciprhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

      /* check whether we have a variable bound and can strenghten the big-M */
      if ( tightenvb && rownnonz == 2 && (SCIPisZero(scip, sciplhs) || SCIPisZero(scip, sciprhs) ) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Real val1;
         SCIP_Real val2;

         val1 = rowvals[0];
         val2 = rowvals[1];

         assert( rowcols[0] != NULL );
         assert( rowcols[1] != NULL );
         var1 = SCIPcolGetVar(rowcols[0]);
         var2 = SCIPcolGetVar(rowcols[1]);
         assert( var1 != NULL );
         assert( var2 != NULL );

         /* check that variables are not locally fixed */
         if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
         {
            /* one coefficient must be 1 and the other negative */
            if ( (SCIPisEQ(scip, val1, 1.0) || SCIPisEQ(scip, val2, 1.0)) && ( SCIPisNegative(scip, val1) || SCIPisNegative(scip, val2) ) )
            {
               /* We want x - a z <= 0 or x - a z >= 0, where var1 = x and var2 = z; possibly swap variables otherwise */
               if ( ! SCIPisEQ(scip, val1, 1.0) || ! SCIPisNegative(scip, val2) )
               {
                  SCIPswapPointers((void**) &var1, (void**) &var2);

                  val2 = val1;
                  swapped = TRUE;
               }

               /* var2 needs to be binary */
               if ( SCIPvarIsBinary(var2) )
               {
                  if ( SCIPisZero(scip, sciprhs) )
                  {
                     if ( SCIPisLT(scip, SCIPvarGetUbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetUbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }

                  if ( SCIPisZero(scip, sciplhs) )
                  {
                     if ( SCIPisGT(scip, SCIPvarGetLbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetLbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }
               }
            }
         }
      }

      for (j = 0; j < rownnonz; j++)
      {
         /* if the Big-M was tightened, we use the new value (the position where this new value is used is dependant on wheter we needed to swap) */
         if ( tightened && ( (swapped && (j == 0)) || ((! swapped) && (j == 1)) ) ) /* use the tightened value */
         {
            if ( SCIPisFeasGT(scip, REALABS(tightenedval), 0.0) )
            {
               assert( SCIPcolGetVar(rowcols[j]) != 0 );
               colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
               rowind[nnonz] = nconss;
               val[nnonz] = tightenedval;
               nnonz++;
            }
         }
         else if ( SCIPisFeasGT(scip, REALABS(rowvals[j]), 0.0))
         {
            assert( SCIPcolGetVar(rowcols[j]) != 0 );
            colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = nconss;
            val[nnonz] = rowvals[j];
            nnonz++;
         }
      }
      lhs[nconss] = primalobj ? sciplhs : (SCIPisInfinity(scip, -sciplhs) ? -sciplhs : 0.0);
      rhs[nconss] = primalobj ? sciprhs : (SCIPisInfinity(scip, sciprhs) ? sciprhs : 0.0);
      nconss++;
   }

   /* delete the old LP-block from the sdpi */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &nrowssdpi) );
   if ( nrowssdpi > 0 )
   {
      SCIP_CALL( SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1) );
   }

   /* add the LP-block to the sdpi */
   SCIP_CALL( SCIPsdpiAddLPRows(sdpi, nconss, lhs, rhs, nnonz, (const int*)rowind, (const int*)colind, val) );

   /* free the remaining arrays */
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &colind);
   SCIPfreeBufferArray(scip, &rowind);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

   /* update bounds */

   /* get the variables */
   vars = SCIPgetVars(scip);
   assert( vars != NULL );

   /* prepare arrays of bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );

   /* get new bounds and objective coefficients */
   for (i = 0; i < nvars; i++)
   {
      assert( vars[i] != NULL );
      lb[i] = primalobj ? SCIPvarGetLbLocal(vars[i]) : (SCIPisInfinity(scip, -SCIPvarGetLbLocal(vars[i])) ? SCIPvarGetLbLocal(vars[i]) : 0.0);
      ub[i] = primalobj ? SCIPvarGetUbLocal(vars[i]) : (SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[i])) ? SCIPvarGetUbLocal(vars[i]) : 0.0);
      inds[i] = i; /* we want to change all bounds, so all indices are included in inds */
      obj[i] = dualobj ? SCIPvarGetObj(vars[i]) : 0.0;
      objinds[i] = i;
   }

   /* inform interface */
   SCIP_CALL( SCIPsdpiChgBounds(sdpi, nvars, inds, lb, ub) );
   SCIP_CALL( SCIPsdpiChgObj(sdpi, nvars, objinds, obj) );

   /* free the bounds-arrays */
   SCIPfreeBufferArray(scip, &objinds);
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}

/** calculate relaxation and process the relaxation results */
static
SCIP_RETCODE calcRelax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< data of the relaxator */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound          /**< pointer to store lowerbound */
   )
{
   char saveconsname[SCIP_MAXSTRLEN];
   SCIP_SDPSOLVERSETTING startsetting;
   SCIP_SDPSOLVERSETTING usedsetting;
   SCIP_CONS* savedsetting;
   SCIP_CONS** conss;
   SCIP_VAR** vars;
   SCIP_SDPI* sdpi;
   SCIP_Bool rootnode;
   SCIP_Bool enforceslater;
   SCIP_Real timelimit;
   SCIP_Real objforscip;
   SCIP_Real* solforscip;
   SCIP_SDPSLATERSETTING slatersetting;
   SCIP_SDPSLATER primalslater;
   SCIP_SDPSLATER dualslater;
   int naddediters;
   int naddedsdpcalls;
   int nblocks;
   int nvars;
   int b;
   int i;
   int v;

   SCIPdebugMessage("calcRelax called\n");

   assert( scip != NULL );
   assert( relaxdata != NULL );
   assert( result != NULL );
   assert( lowerbound != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars (scip);

   sdpi = relaxdata->sdpi;
   assert( sdpi != NULL );

   if ( relaxdata->objlimit )
   {
      /* set the objective limit */
      assert( SCIPgetUpperbound(scip) > -SCIPsdpiInfinity(sdpi) );
      SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPgetUpperbound(scip)) );
   }
   /* if this is the root node and we cannot solve the problem, we want to check for the Slater condition independent from the SCIP parameter */
   rootnode = ! SCIPnodeGetParent(SCIPgetCurrentNode(scip));

   /* find settings to use for this relaxation */
   if ( rootnode || (SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == relaxdata->settingsresetofs) ||
      ( relaxdata->settingsresetfreq > 0 && ((SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) - relaxdata->settingsresetofs) % relaxdata->settingsresetfreq == 0)) ||
      (strcmp(SCIPsdpiGetSolverName(), "DSDP") == 0) || (strstr(SCIPsdpiGetSolverName(), "Mosek") != NULL))
   {
      startsetting = SCIP_SDPSOLVERSETTING_UNSOLVED; /* in the root node we have no information, at each multiple of resetfreq we reset */
   }
   else
   {
      SCIP_CONSHDLR* conshdlr;
      int parentconsind;

      /* get constraint handler */
      conshdlr = SCIPfindConshdlr(scip, "Savedsdpsettings");
      if ( conshdlr == NULL )
      {
         SCIPerrorMessage("Savedsdpsettings constraint handler not found!\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* get startsettings of parent node, usually it will be the last active constraint of the corresponding constraint handler, so we iterate from
       * the end of the list until we find the correct one */
      conss = SCIPconshdlrGetConss(conshdlr);
      parentconsind = SCIPconshdlrGetNActiveConss(conshdlr) - 1;
      (void) SCIPsnprintf(saveconsname, SCIP_MAXSTRLEN, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(scip))));

      while ( parentconsind >= 0 && strcmp(saveconsname, SCIPconsGetName(conss[parentconsind])) )
         parentconsind--;
      if ( parentconsind >= 0 )
         startsetting = SCIPconsSavedsdpsettingsGetSettings(scip, conss[parentconsind]);
      else
      {
         SCIPdebugMessage("Startsetting from parent node not found, restarting with fastest settings!\n");
         startsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
      }
   }

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if ( timelimit <= 0.0 )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }
   }

   /* if no dual bound is known (we are in the root node and not only repropagating), we will have to abort, so we want
    * to check the Slater condition in this case */
   enforceslater = SCIPisInfinity(scip, -1 * SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip)));

   /* solve the problem (using warmstarts if parameter is true and we are not in the root node and all neccessary data is available) */
   if ( ( ! SCIPnodeGetParent(SCIPgetCurrentNode(scip))) || ( ! relaxdata->warmstart ) || ((relaxdata->warmstartiptype == 2) &&
         SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) && ((SCIPsdpiDoesWarmstartNeedPrimal() && ! relaxdata->ipXexists) || (! relaxdata->ipZexists))) )
   {
      SCIP_CALL(SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, startsetting, enforceslater, timelimit));
   }
   else if ( relaxdata->warmstart && (relaxdata->warmstartprimaltype != 2) && (relaxdata->warmstartiptype == 2) && SCIPisEQ(scip, relaxdata->warmstartipfactor, 1.0) )
   {
      SCIP_Real* ipy;

      /* if we warmstart with the analytic center, we give a pointer to the arrays in relaxdata (just the sol needs to be transformed to a vector first) */

      SCIP_CALL( SCIPallocBufferArray(scip, &ipy, nvars) );
      for (v = 0; v < nvars; v++)
         ipy[v] = SCIPgetSolVal(scip, relaxdata->ipy, SCIPsdpVarmapperGetSCIPvar(relaxdata->varmapper, v));

#ifdef SCIP_PRINT_WARMSTART
      SCIPdebugMessage("warmstart using the following analytic centers:\n");
      for (v = 0; v < nvars; v++)
         SCIPdebugMessage("y[%d] = %f\n", v, ipy[v]);
      if ( SCIPsdpiDoesWarmstartNeedPrimal() )
      {
         for (b = 0; b < relaxdata->nblocks; b++)
         {
            SCIPdebugMessage("dual block %d\n", b);
            for (i = 0; i < relaxdata->ipZnblocknonz[b]; i++)
            {
               SCIPdebugMessage("Z(%d,%d)=%f\n", relaxdata->ipZrow[b][i], relaxdata->ipZcol[b][i], relaxdata->ipZval[b][i]);
            }
         }
         for (b = 0; b < relaxdata->nblocks; b++)
         {
            SCIPdebugMessage("primal block %d\n", b);
            for (i = 0; i < relaxdata->ipXnblocknonz[b]; i++)
            {
               SCIPdebugMessage("X(%d,%d)=%f\n", relaxdata->ipXrow[b][i], relaxdata->ipXcol[b][i], relaxdata->ipXval[b][i]);
            }
         }
      }
#endif

      SCIP_CALL(SCIPsdpiSolve(sdpi, ipy, relaxdata->ipZnblocknonz, relaxdata->ipZrow, relaxdata->ipZcol, relaxdata->ipZval, relaxdata->ipXnblocknonz,
                        relaxdata->ipXrow, relaxdata->ipXcol, relaxdata->ipXval, startsetting, enforceslater, timelimit));

      SCIPfreeBufferArray(scip, &ipy);
   }
   else
   {
      SCIP_CONSHDLR* conshdlr;
      SCIP_SOL* dualsol;
      SCIP_Real* starty = NULL;
      int* startZnblocknonz = NULL;
      int** startZrow = NULL;
      int** startZcol = NULL;
      SCIP_Real** startZval = NULL;
      int* startXnblocknonz = NULL;
      int** startXrow = NULL;
      int** startXcol = NULL;
      SCIP_Real** startXval = NULL;
      int parentconsind;
      SCIP_Longint parentnodenumber;
      SCIP_VAR* var;

      /* find starting solution as optimal solution of parent node */

      /* get constraint handler */
      conshdlr = SCIPfindConshdlr(scip, "Savesdpsol");
      if ( conshdlr == NULL )
      {
         SCIPerrorMessage("Savesdpsol constraint handler not found\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* get saveconstraint of parent node, usually it will be the last active constraint of the corresponding constraint handler, so we iterate from
       * the end of the list until we find the correct one */
      conss = SCIPconshdlrGetConss(conshdlr);
      parentconsind = SCIPconshdlrGetNActiveConss(conshdlr) - 1;
      parentnodenumber = SCIPnodeGetNumber(SCIPnodeGetParent(SCIPgetCurrentNode(scip)));

      while ( parentconsind >= 0 && SCIPconsSavesdpsolGetNodeIndex(scip, conss[parentconsind]) != parentnodenumber)
         parentconsind--;

      /* if there are no savesdpsol constraints (e.g. because the parent node couldn't be solved successfully), solve
       * without warmstart
       */
      if ( parentconsind < 0 )
      {
         SCIPdebugMessage("Starting SDP-Solving from scratch since no warmstart information available for node %lld\n", parentnodenumber);
         SCIP_CALL(SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, startsetting, enforceslater, timelimit));
      }
      else
      {
         SCIPdebugMessage("Using warmstartinformation from node %lld\n", parentnodenumber);

         /* get solution */
         dualsol = SCIPconsSavesdpsolGetDualVector(scip, conss[parentconsind]);

         /* allocate memory */
         SCIP_CALL( SCIPallocBufferArray(scip, &starty, nvars) );

         /* transform solution to vector for SDPI and check if it is still feasible for the variable bounds, otherwise round it */
         for (v = 0; v < nvars; v++)
         {
            var = SCIPsdpVarmapperGetSCIPvar(relaxdata->varmapper, v);
            starty[v] = SCIPgetSolVal(scip, dualsol, var);
            /* correct solution to new bounds unless warmstartproject == 1 */
            if (SCIPisLT(scip, starty[v], SCIPvarGetLbLocal(var)) && (relaxdata->warmstartproject == 2 || relaxdata->warmstartproject == 3 || relaxdata->warmstartproject == 4))
            {
               starty[v] = SCIPvarGetLbLocal(var);
               /* update solution (used to compute dual matrix) according to new bounds */
               SCIP_CALL( SCIPsetSolVal(scip, dualsol, var, SCIPvarGetLbLocal(var)) );
            }
            else if (SCIPisGT(scip, starty[v], SCIPvarGetUbLocal(var)) && (relaxdata->warmstartproject == 2 || relaxdata->warmstartproject == 3 || relaxdata->warmstartproject == 4))
            {
               starty[v] = SCIPvarGetUbLocal(var);
               /* update solution (used to compute dual matrix) according to new bounds */
               SCIP_CALL( SCIPsetSolVal(scip, dualsol, var, SCIPvarGetUbLocal(var)) );
            }

            /* if we take a convex combination, adjust y accordingly (if we use rounding problems, we recompute y later anyways) */
            if ( SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) && relaxdata->warmstartproject != 4 )
            {
               if ( relaxdata->warmstartiptype == 1 )
               {
                  /* we take a convex combination with 0, so we just scale */
                  starty[v] *= 1 - relaxdata->warmstartipfactor;
               }
               else if ( relaxdata->warmstartiptype == 2 )
               {
                  /* take the convex combination with the saved analytic center */
                  starty[v] = (1 - relaxdata->warmstartipfactor) * starty[v] + relaxdata->warmstartipfactor * SCIPgetSolVal(scip, relaxdata->ipy, var);
               }
            }
         }

         /* if the SDP-solver needs the primal solution (and the dual matrix) in addition to the dual vector... */
         if ( SCIPsdpiDoesWarmstartNeedPrimal() )
         {
            SCIP_CONSHDLR* sdpconshdlr;
            SCIP_CONS** sdpblocks;
            SCIP_COL** rowcols;
            SCIP_ROW** rows;
            int blocksize;
            int nrows;
            int rownnonz;
            int r;
            SCIP_Real maxprimalentry = 0.0;
            SCIP_Real maxdualentry;
            SCIP_Real identitydiagonal = 0.0;
            SCIP_Real rowval;
            SCIP_Real* rowvals;
            SCIP_Bool* diagentryexists;

            sdpconshdlr = SCIPfindConshdlr(scip, "SDP");
            nblocks = SCIPconshdlrGetNConss(sdpconshdlr);
            sdpblocks = SCIPconshdlrGetConss(sdpconshdlr);
            SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

            SCIP_CALL( SCIPallocBufferArray(scip, &startZnblocknonz, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startZrow, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startZcol, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startZval, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startXnblocknonz, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startXrow, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startXcol, nblocks + 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &startXval, nblocks + 1) );

            /* compute the scaling factor for the dual identity matrix (for numerical stability, this should be at least 1) */
            if ( relaxdata->warmstartiptype == 1 )
            {
               maxprimalentry = SCIPconsSavesdpsolGetMaxPrimalEntry(scip, conss[0]);
               if ( SCIPisLT(scip, maxprimalentry, 1.0) )
                  maxprimalentry = 1.0;
            }

            /* iterate over all blocks and fill X and Z */
            for (b = 0; b < nblocks; b++)
            {
               blocksize = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);

               if ( relaxdata->warmstartproject == 3 || relaxdata->warmstartproject == 4 )
               {
                  /* since we later take the projection onto the psd cone, we cannot a priori determine the size, so we take the maximum possible */
                  startZnblocknonz[b] = blocksize * (blocksize + 1) / 2;
               }
               else
               {
                  startZnblocknonz[b] = SCIPconsSdpComputeUbSparseSdpMatrixLength(scip, sdpblocks[b]);

                  /* since we take a convex combination with either the identity matrix or the analytic center, we have to allocate memory for that as well */
                  if ( SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) )
                  {
                     if ( relaxdata->warmstartiptype == 1 )
                        startZnblocknonz[b] += blocksize;
                     else if ( relaxdata->warmstartiptype == 2 )
                        startZnblocknonz[b] += relaxdata->ipZnblocknonz[b];
                  }
               }

               SCIP_CALL( SCIPallocBufferArray(scip, &startZrow[b], startZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZcol[b], startZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZval[b], startZnblocknonz[b]) );

               /* compute Z matrix */
               SCIP_CALL( SCIPconsSdpComputeSparseSdpMatrix(scip, sdpblocks[b], dualsol, &(startZnblocknonz[b]), startZrow[b], startZcol[b], startZval[b]) );

               /* compute projection onto psd cone (computed as U * diag(lambda_i_+) * U^T where U consists of the eigenvectors of the matrix) */
               if ( relaxdata->warmstartproject == 3 )
               {
                  SCIP_Real* fullZmatrix;
                  SCIP_Real* eigenvalues;
                  SCIP_Real* eigenvectors;
                  SCIP_Real* scaledeigenvectors;
                  SCIP_Real matrixsize;
                  SCIP_Real epsilon;
                  int c;
                  int matrixpos;

                  matrixsize = blocksize * blocksize;

                  SCIP_CALL( SCIPallocBufferArray(scip, &fullZmatrix, matrixsize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &eigenvectors, matrixsize) );

                  SCIP_CALL( expandSparseMatrix(startZnblocknonz[b], blocksize, startZrow[b], startZcol[b], startZval[b], fullZmatrix) );

                  SCIP_CALL( SCIPlapackComputeEigenvectorDecomposition(SCIPbuffer(scip), blocksize, fullZmatrix, eigenvalues, eigenvectors) );

                  /* duplicate memory of eigenvectors to compute diag(lambda_i_+) * U^T */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &scaledeigenvectors, eigenvectors, matrixsize) );

                  /* set all negative eigenvalues to zero (using the property that LAPACK returns them in ascending order) */
                  i = 0;
                  while (i < blocksize && SCIPisLT(scip, eigenvalues[i], relaxdata->warmstartprojminevdual) )
                  {
                     eigenvalues[i] = relaxdata->warmstartprojminevdual;
                     i++;
                  }

                  /* compute diag(lambda_i_+) * U^T */
                  SCIP_CALL( scaleTransposedMatrix(blocksize, scaledeigenvectors, eigenvalues) );

                  /* compute U * [diag(lambda_i_+) * U^T] (note that transposes are switched because LAPACK uses column-first-format) */
                  SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksize, blocksize, eigenvectors, TRUE, blocksize, blocksize, scaledeigenvectors,
                        FALSE, fullZmatrix) );

                  /* extract sparse matrix from projection */
                  startZnblocknonz[b] = 0;
                  epsilon = SCIPepsilon(scip);
                  for (r = 0; r < blocksize; r++)
                  {
                     for (c = r; c < blocksize; c++)
                     {
                        matrixpos = r * blocksize + c;
                        if ( REALABS(fullZmatrix[matrixpos]) > epsilon )
                        {
                           startZrow[b][startZnblocknonz[b]] = r;
                           startZcol[b][startZnblocknonz[b]] = c;
                           startZval[b][startZnblocknonz[b]] = fullZmatrix[matrixpos];
                           startZnblocknonz[b]++;
                        }
                     }
                  }

                  /* free memory */
                  SCIPfreeBufferArray(scip, &scaledeigenvectors);
                  SCIPfreeBufferArray(scip, &eigenvectors);
                  SCIPfreeBufferArray(scip, &eigenvalues);
                  SCIPfreeBufferArray(scip, &fullZmatrix);
               }

               if ( relaxdata->warmstartprimaltype == 1 )
               {
                  /* we set X to maxprimalentry times the identity matrix */
                  if ( relaxdata->warmstartiptype == 1 )
                  {
                     startXnblocknonz[b] = blocksize;
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b]) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b]) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b]) );
                     for (i = 0; i < startXnblocknonz[b]; i++)
                     {
                        startXrow[b][i] = i;
                        startXcol[b][i] = i;
                        startXval[b][i] = maxprimalentry;
                     }
                  }
               }
               else if ( relaxdata->warmstartprimaltype == 2 )
               {
                  startXnblocknonz[b] = startZnblocknonz[b];
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b]) );
                  for (i = 0; i < startZnblocknonz[b]; i++)
                  {
                     startXrow[b][i] = startZrow[b][i];
                     startXcol[b][i] = startZcol[b][i];
                     startXval[b][i] = 1 / startZval[b][i];
                  }
               }
               else if ( relaxdata->warmstartprimaltype != 3 && relaxdata->warmstartproject != 4 )
               {
                  SCIPerrorMessage("Unknown value %d for warmstartprimaltype.\n", relaxdata->warmstartprimaltype);
                  SCIPABORT();
               }
            }

            if ( relaxdata->warmstartproject != 4 )
            {
               /** fill LP-block */
               SCIP_CALL( SCIPallocBufferArray(scip, &startZrow[b], 2 * nrows + 2 * nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZcol[b], 2 * nrows + 2 * nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZval[b], 2 * nrows + 2 * nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], 2 * nrows + 2 * nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], 2 * nrows + 2 * nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], 2 * nrows + 2 * nvars) );

               /* to get a positive definite matrix, all the entries need to be strictly positive */
               startZnblocknonz[b] = 2 * nrows + 2 * nvars;
               startXnblocknonz[b] = 2 * nrows + 2 * nvars;

               for (r = 0; r < nrows; r++)
               {
                  /* compute row value for current solution */
                  rowval = 0.0;
                  rownnonz = SCIProwGetNNonz(rows[r]);
                  rowvals = SCIProwGetVals(rows[r]);
                  rowcols = SCIProwGetCols(rows[r]);
                  for (i = 0; i < rownnonz; i++)
                     rowval += SCIPgetSolVal(scip, dualsol, SCIPcolGetVar(rowcols[i])) * rowvals[i];

                  startZrow[b][2*r] = 2*r;
                  startZcol[b][2*r] = 2*r;
                  startZval[b][2*r] = rowval - (SCIProwGetLhs(rows[r]) - SCIProwGetConstant(rows[r]));

                  if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[b][2*r], relaxdata->warmstartprojminevdual) )
                     startZval[b][2*r] = relaxdata->warmstartprojminevdual;
                  /* we only take the convex combination if the value is less than one, since the maxblockentry is equal to the value
                   * otherwise, so taking the convex combination doesn't change anything in that case
                   */
                  else if ( relaxdata->warmstartiptype == 1 && SCIPisLT(scip, startZval[b][2*r], 1.0) )
                  {
                     /* since we want the value to be strictly positive, if the original entry is negative we just set it to warmstartipfactor */
                     if ( SCIPisLT(scip, startZval[b][2*r], 0.0) )
                        startZval[b][2*r] = relaxdata->warmstartipfactor;
                     else
                        startZval[b][2*r] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*r] + relaxdata->warmstartipfactor;
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     startZval[b][2*r] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*r] + relaxdata->warmstartipfactor * relaxdata->ipZval[b][2*r];

                     /* if this is non-positive, we shift it to a strictly positive value */
                     if ( SCIPisLT(scip, startZval[b][2*r], WARMSTART_MINVAL) )
                        startZval[b][2*r] = WARMSTART_MINVAL;
                  }

                  if ( relaxdata->warmstartpreoptsol && startZval[b][2*r] < WARMSTART_PREOPT_MIN_Z_LPVAL )
                     startZval[b][2*r] = WARMSTART_PREOPT_MIN_Z_LPVAL;

                  startZrow[b][2*r + 1] = 2*r + 1;
                  startZcol[b][2*r + 1] = 2*r + 1;
                  startZval[b][2*r + 1] = SCIProwGetRhs(rows[r]) - SCIProwGetConstant(rows[r]) - rowval;

                  if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[b][2*r + 1], relaxdata->warmstartprojminevdual) )
                     startZval[b][2*r + 1] = relaxdata->warmstartprojminevdual;
                  else if ( relaxdata->warmstartiptype == 1 && SCIPisLT(scip, startZval[b][2*r + 1], 1.0) )
                  {
                     /* since we want the value to be strictly positive, if the original entry is negative we just set it to warmstartipfactor */
                     if ( SCIPisLT(scip, startZval[b][2*r + 1], 0.0) )
                        startZval[b][2*r + 1] = relaxdata->warmstartipfactor;
                     else
                        startZval[b][2*r + 1] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*r + 1] + relaxdata->warmstartipfactor;
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     startZval[b][2*r + 1] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*r + 1] + relaxdata->warmstartipfactor * relaxdata->ipZval[b][2*r + 1];

                     /* if this is non-positive, we shift it to a strictly positive value */
                     if ( SCIPisLT(scip, startZval[b][2*r + 1], WARMSTART_MINVAL) )
                        startZval[b][2*r + 1] = WARMSTART_MINVAL;
                  }

                  if ( relaxdata->warmstartpreoptsol && startZval[b][2*r + 1] < WARMSTART_PREOPT_MIN_Z_LPVAL )
                     startZval[b][2*r + 1] = WARMSTART_PREOPT_MIN_Z_LPVAL;

                  if ( relaxdata->warmstartprimaltype == 1 && relaxdata->warmstartiptype == 1 )
                  {
                     startXrow[b][2*r] = 2*r;
                     startXcol[b][2*r] = 2*r;
                     startXval[b][2*r] = maxprimalentry;
                     startXrow[b][2*r + 1] = 2*r + 1;
                     startXcol[b][2*r + 1] = 2*r + 1;
                     startXval[b][2*r + 1] = maxprimalentry;
                  }
                  else if ( relaxdata->warmstartprimaltype == 2 )
                  {
                     startXrow[b][2*r] = startZrow[b][2*r];
                     startXcol[b][2*r] = startZcol[b][2*r];
                     startXval[b][2*r] = 1 / startZval[b][2*r];
                     startXrow[b][2*r + 1] = startZrow[b][2*r + 1];
                     startXcol[b][2*r + 1] = startZcol[b][2*r + 1];
                     startXval[b][2*r + 1] = 1 / startZval[b][2*r + 1];
                  }
                  else if ( relaxdata->warmstartprimaltype != 3 && relaxdata->warmstartiptype == 1 )
                  {
                     SCIPerrorMessage("Unknown value %d for warmstartprimaltype.\n", relaxdata->warmstartprimaltype);
                     SCIPABORT();
                  }
               }

               for (v = 0; v < nvars; v++)
               {
                  startZrow[b][2*nrows + 2*v] = 2*nrows + 2*v;
                  startZcol[b][2*nrows + 2*v] = 2*nrows + 2*v;
                  startZval[b][2*nrows + 2*v] = SCIPgetSolVal(scip, dualsol, vars[v]) - SCIPvarGetLbLocal(vars[v]);
                  if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[b][2*nrows + 2*v], relaxdata->warmstartprojminevdual) )
                     startZval[b][2*nrows + 2*v] = relaxdata->warmstartprojminevdual;
                  else if ( relaxdata->warmstartiptype == 1 && SCIPisLT(scip, startZval[b][2*nrows + 2*v], 1.0) )
                  {
                     /* since we want the value to be strictly positive, if the original entry is negative we just set it to warmstartipfactor */
                     if ( SCIPisLT(scip, startZval[b][2*nrows + 2*v], 0.0) )
                        startZval[b][2*nrows + 2*v] = relaxdata->warmstartipfactor;
                     else
                        startZval[b][2*nrows + 2*v] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*nrows + 2*v] + relaxdata->warmstartipfactor;
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     startZval[b][2*nrows + 2*v] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*nrows + 2*v] + relaxdata->warmstartipfactor * relaxdata->ipZval[b][2*nrows + 2*v];

                     /* if this is non-positive, we shift it to a strictly positive value */
                     if ( SCIPisLT(scip, startZval[b][2*nrows + 2*v], WARMSTART_MINVAL) )
                        startZval[b][2*nrows + 2*v] = WARMSTART_MINVAL;
                  }

                  if ( relaxdata->warmstartpreoptsol && startZval[b][2*nrows + 2*v] < WARMSTART_PREOPT_MIN_Z_LPVAL )
                     startZval[b][2*nrows + 2*v] = WARMSTART_PREOPT_MIN_Z_LPVAL;

                  startZrow[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                  startZcol[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                  startZval[b][2*nrows + 2*v + 1] = SCIPvarGetUbLocal(vars[v]) - SCIPgetSolVal(scip, dualsol, vars[v]);
                  if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[b][2*nrows + 2*v + 1], relaxdata->warmstartprojminevdual) )
                     startZval[b][2*nrows + 2*v + 1] = relaxdata->warmstartprojminevdual;
                  else if ( relaxdata->warmstartiptype == 1 && SCIPisLT(scip, startZval[b][2*nrows + 2*v + 1], 1.0) )
                  {
                     /* since we want the value to be strictly positive, if the original entry is negative we just set it to warmstartipfactor */
                     if ( SCIPisLT(scip, startZval[b][2*nrows + 2*v + 1], 0.0) )
                        startZval[b][2*nrows + 2*v + 1] = relaxdata->warmstartipfactor;
                     else
                        startZval[b][2*nrows + 2*v + 1] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*nrows + 2*v + 1] + relaxdata->warmstartipfactor;
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     startZval[b][2*nrows + 2*v + 1] = (1 - relaxdata->warmstartipfactor) * startZval[b][2*nrows + 2*v] + relaxdata->warmstartipfactor * relaxdata->ipZval[b][2*nrows + 2*v + 1];

                     /* if this is non-positive, we shift it to a strictly positive value */
                     if ( SCIPisLT(scip, startZval[b][2*nrows + 2*v + 1], WARMSTART_MINVAL) )
                        startZval[b][2*nrows + 2*v + 1] = WARMSTART_MINVAL;
                  }

                  if ( relaxdata->warmstartpreoptsol && startZval[b][2*nrows + 2*v + 1] < WARMSTART_PREOPT_MIN_Z_LPVAL )
                     startZval[b][2*nrows + 2*v + 1] = WARMSTART_PREOPT_MIN_Z_LPVAL;

                  if ( relaxdata->warmstartprimaltype == 1 && relaxdata->warmstartiptype == 1 )
                  {
                     startXrow[b][2*nrows + 2*v] = 2*nrows + 2*v;
                     startXcol[b][2*nrows + 2*v] = 2*nrows + 2*v;
                     startXval[b][2*nrows + 2*v] = maxprimalentry;
                     startXrow[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                     startXcol[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                     startXval[b][2*nrows + 2*v + 1] = maxprimalentry;
                  }
                  else if ( relaxdata->warmstartprimaltype == 2 )
                  {
                     startXrow[b][2*nrows + 2*v] = startZrow[b][2*nrows + 2*v];
                     startXcol[b][2*nrows + 2*v] = startZcol[b][2*nrows + 2*v];
                     startXval[b][2*nrows + 2*v] = 1 / startZval[b][2*nrows + 2*v];
                     startXrow[b][2*nrows + 2*v + 1] = startZrow[b][2*nrows + 2*v + 1];
                     startXcol[b][2*nrows + 2*v + 1] = startZcol[b][2*nrows + 2*v + 1];
                     startXval[b][2*nrows + 2*v + 1] = 1 / startZval[b][2*nrows + 2*v + 1];
                  }
                  else if ( relaxdata->warmstartprimaltype != 3 && relaxdata->warmstartproject == 1 )
                  {
                     SCIPerrorMessage("Unknown value %d for warmstartprimaltype.\n", relaxdata->warmstartprimaltype);
                     SCIPABORT();
                  }
               }
            }

            /* Solve the primal rounding problem
             * \f{eqnarray*}{
             *    \max & & \sum_{k \in K} A_0^{(k)} \bullet (V^{(k)} \text{diag}(\lambda^{(k)}) (V^{(k)})^T) + \sum_{j \in J} c_j x_j - \sum_{i \in I_u} u_i v_i + \sum_{i \in I_\ell} \ell_i w_i \\
             *    \mbox{s.t.} & & \sum_{k \in K} A_i^{(k)} \bullet (V^{(k)} \text{diag}(\lambda^{(k)}) (V^{(k)})^T) + \sum_{j \in J} d_{ij} x_j - 1_{\{u_i < \infty\}} v_i + 1_{\{\ell_i > -\infty\}} w_i = b_i \quad \forall \ i \in I,\\
             *    & & \lambda^{(k)}_i \geq 0 \quad \forall \ k \in K, i \leq n \\
             *    & & x_j \geq 0 \quad \forall \ j \in J,\\
             *    & & v_i \geq 0 \quad \forall \ i \in I_u,\\
             *    & & w_i \geq 0 \quad \forall \ i \in I_\ell,
             * \f}
             * where \f$ V^{(k)} \text{diag}(\bar{\lambda}^{(k)}) (V^{(k)})^T \f$ is an eigenvector decomposition of the optimal primal solution
             * of the parent node, as well as the dual rounding problem
             * \f{eqnarray*}{
             *    \min & & b^T y \\
             *    \mbox{s.t.} & & \sum_{i \in I} d_{ij} y_i \geq c_j \quad \forall \ j \in J, \\
             *    & & \sum_{i \in I} A_i^{(k)} y_i - A_0^{(k)} = V^{(k)} \text{diag}(\lambda^{(k)}) (V^{(k)})^T \quad \forall \ k \in K, \\
             *    & & \ell_i \leq y_i \leq u_i \quad \forall \ i \in I, \\
             *    & & \lambda^{(k)}_i \geq 0 \quad \forall \ k \in K, i \leq n \\
             * \f}
             * where \f$ V^{(k)} \text{diag}(\bar{\lambda}^{(k)}) (V^{(k)})^T \f$ is now an eigenvector decomposition of the optimal solution
             * of the parent node for the dual problem. The matrix equation is reformulated as blocksize * (blocksize + 1) /2 linear constraints
             * over the lower triangular entries.
             */
            if ( relaxdata->warmstartproject == 4 )
            {
               SCIP_VAR** blockvars;
               SCIP_LPI* lpi;
               SCIP_ROW* row;
               struct timeval starttime;
               struct timeval currenttime;
               SCIP_Real** blockval;
               SCIP_Real** blockeigenvalues;
               SCIP_Real** blockeigenvectors;
               SCIP_Real** blockrowvals;
               SCIP_Real* obj;
               SCIP_Real* lb;
               SCIP_Real* ub;
               SCIP_Real* lhs;
               SCIP_Real* rhs;
               SCIP_Real* val;
               SCIP_Real* blockconstval;
               SCIP_Real* scaledeigenvectors;
               SCIP_Real* fullXmatrix;
               SCIP_Real* fullZmatrix;
               SCIP_Real rowlhs;
               SCIP_Real rowrhs;
               SCIP_Real varobj;
               SCIP_Real epsilon;
               SCIP_Real primalroundobj;
               SCIP_Real dualroundobj;
               int** blockcol;
               int** blockrow;
               int** blockrowcols;
               int* beg;
               int* ind;
               int* blocknvarnonz;
               int* blockconstcol;
               int* blockconstrow;
               int* blocksizes;
               int* nblockrownonz;
               int* rowinds;
               int pos;
               int indpos;
               int startpos;
               int blocknvars;
               int blocknnonz;
               int arraylength;
               int blockconstnnonz;
               int varind;
               int roundingvars;
               int matrixsize;
               int evind;
               int c;
               int matrixpos;
               int nroundingrows;
               int j;
               int nremovedentries;

               /* since the main purpose of the rounding problem approach is detecting infeasibility through the restricted primal problem,
                * it doesn't make sense to use this approach unless the whole primal solution is saved */
               if ( relaxdata->warmstartprimaltype != 3 )
               {
                  SCIPerrorMessage("Invalid parameter combination, use relax/warmstartproject = 4 only with relax/warmstartprimaltype = 3.\n");
                  return SCIP_PARAMETERWRONGVAL;
               }

               TIMEOFDAY_CALL( gettimeofday(&starttime, NULL) );/*lint !e438, !e550, !e641 */

               /* since we cannot compute the number of nonzeros of the solution of the rounding problem beforehand, we allocate the maximum possible (blocksize * (blocksize + 1) / 2 */
               for (b = 0; b < nblocks; b++)
               {
                  matrixsize = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);
                  matrixsize *= (matrixsize + 1) * 0.5;
                  startXnblocknonz[b] = matrixsize;

                  SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], matrixsize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], matrixsize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], matrixsize) );
               }

               /* allocate memory for LP and variable bound block of warmstart matrix (note that we need to allocate the maximum, since additional
                * entries may be generated through the convex combination */
               SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[nblocks], 2 * nvars + 2 * nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[nblocks], 2 * nvars + 2 * nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startXval[nblocks], 2 * nvars + 2 * nrows) );
               startXnblocknonz[nblocks] = 2 * nvars + 2 * nrows;

               SCIP_CALL( SCIPconsSavesdpsolGetPrimalMatrix(scip, conss[parentconsind], nblocks + 1, startXnblocknonz, startXrow, startXcol, startXval) );

               lpi = relaxdata->lpi;

               /* clear the old LP */
               SCIP_CALL( SCIPlpiClear(lpi) );

               /* we want to maximize for the primal rounding problem */
               SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MAXIMIZE) );

               /* if the varmapper has a different number of variables than SCIP, we might get problems with empty spots in the obj/lb/ub arrays */
               assert( SCIPsdpVarmapperGetNVars(relaxdata->varmapper) == nvars );

               /* initialize the rows */
               SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nvars) );
               for (v = 0; v < nvars; v++)
               {
                  varobj = SCIPvarGetObj(vars[v]);
                  lhs[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] = varobj;
                  rhs[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] = varobj;
               }

               /* this vector will be used to later map the solution to the correspoding blocks */
               SCIP_CALL( SCIPallocBufferArray(scip, &blocksizes, nblocks + 2) );
               SCIP_CALL( SCIPallocBufferArray(scip, &blockeigenvalues, nblocks) );
               SCIP_CALL( SCIPallocBufferArray(scip, &blockeigenvectors, nblocks) );

               /* initialize rows; columns will be added blockwise later (note that we reenter the whole problem each time
                * and do not allow warmstarts, we could store the LP- and bound-constraints between nodes (TODO, but would
                * need to check for added/removed cuts then), but warmstarting doesn't seem to make sense since the whole
                * SDP part is changed, which should form the main difficulty when solving (MI)SDPs */
               SCIP_CALL( SCIPlpiAddRows(lpi, nvars, lhs, rhs, NULL, 0, NULL, NULL, NULL) );

               SCIPfreeBufferArray(scip, &rhs);
               SCIPfreeBufferArray(scip, &lhs);

               /* add columns corresponding to bound constraints (at most we get two entries [lb & ub] per variable) */
               SCIP_CALL( SCIPallocBufferArray(scip, &obj, 2*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &lb, 2*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &ub, 2*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &beg, 2*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &ind, 2*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &val, 2*nvars) );

               /* iterate over all finite variable bounds and set the corresponding entries */
               pos = 0;
               for (v = 0; v < nvars; v++)
               {
                  if ( ! SCIPisInfinity(scip, -1 * SCIPvarGetLbLocal(vars[v])) )
                  {
                     obj[pos] = SCIPvarGetLbLocal(vars[v]);
                     lb[pos] = 0.0;
                     ub[pos] = SCIPlpiInfinity(lpi);
                     beg[pos] = pos;
                     ind[pos] = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v]);
                     val[pos] = 1.0;
                     pos++;
                  }
                  if ( ! SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[v])) )
                  {
                     obj[pos] = -1 * SCIPvarGetUbLocal(vars[v]);
                     lb[pos] = 0.0;
                     ub[pos] = SCIPlpiInfinity(lpi);
                     beg[pos] = pos;
                     ind[pos] = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v]);
                     val[pos] = -1.0;
                     pos++;
                  }
               }

               SCIP_CALL( SCIPlpiAddCols(lpi, pos, obj, lb, ub, NULL, pos, beg, ind, val) );
               blocksizes[0] = pos;
               roundingvars = pos;

               SCIPfreeBufferArray(scip, &val);
               SCIPfreeBufferArray(scip, &ind);
               SCIPfreeBufferArray(scip, &beg);
               SCIPfreeBufferArray(scip, &ub);
               SCIPfreeBufferArray(scip, &lb);
               SCIPfreeBufferArray(scip, &obj);

               /* add columns corresponding to linear constraints (at most we get two columns per ranged row) */
               SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

               SCIP_CALL( SCIPallocBufferArray(scip, &obj, 2*nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &lb, 2*nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &ub, 2*nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &beg, 2*nrows) );
               SCIP_CALL( SCIPallocBufferArray(scip, &ind, 2*nrows*nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &val, 2*nrows*nvars) );

               pos = 0;
               indpos = 0;
               /* iterate over all LP-ranged-rows and add corresponding entries if lhs and/or rhs is finite */
               for (r = 0; r < nrows; r++)
               {
                  row = rows[r];
                  assert( row != 0 );
                  rownnonz = SCIProwGetNNonz(row);

                  rowvals = SCIProwGetVals(row);
                  rowcols = SCIProwGetCols(row);
                  rowlhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
                  rowrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

                  if ( ! SCIPisInfinity(scip, -1 * rowlhs) )
                  {
                     obj[pos] = rowlhs;
                     lb[pos] = 0.0;
                     ub[pos] = SCIPlpiInfinity(lpi);
                     beg[pos] = indpos;
                     for (i = 0; i < rownnonz; i++)
                     {
                        ind[indpos] = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, SCIPcolGetVar(rowcols[i]));
                        val[indpos] = rowvals[i];
                        indpos++;
                     }
                     pos++;
                  }
                  /* since we assume >= constraints in the dual, we have to multiply all entries of <= constraints by -1 */
                  if ( ! SCIPisInfinity(scip, rowrhs) )
                  {
                     obj[pos] = -1 * rowrhs;
                     lb[pos] = 0.0;
                     ub[pos] = SCIPlpiInfinity(lpi);
                     beg[pos] = indpos;
                     for (i = 0; i < rownnonz; i++)
                     {
                        ind[indpos] = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, SCIPcolGetVar(rowcols[i]));
                        val[indpos] = -1 * rowvals[i];
                        indpos++;
                     }
                     pos++;
                  }
               }

               SCIP_CALL( SCIPlpiAddCols(lpi, pos, obj, lb, ub, NULL, indpos, beg, ind, val) );
               blocksizes[1] = pos;
               roundingvars += pos;

               SCIPfreeBufferArray(scip, &val);
               SCIPfreeBufferArray(scip, &ind);
               SCIPfreeBufferArray(scip, &beg);
               SCIPfreeBufferArray(scip, &ub);
               SCIPfreeBufferArray(scip, &lb);
               SCIPfreeBufferArray(scip, &obj);

               /* finally add columns corresponding to SDP-constraints */
               for (b = 0; b < nblocks; b++)
               {
                  /* get data for this SDP block */
                  SCIP_CALL( SCIPconsSdpGetNNonz(scip, sdpblocks[b], NULL, &blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blocknvarnonz, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockcol, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockrow, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockval, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstcol, blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstrow, blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstval, blockconstnnonz) );

                  arraylength = nvars;
                  SCIP_CALL( SCIPconsSdpGetData(scip, sdpblocks[b], &blocknvars, &blocknnonz, &blocksize, &arraylength, blocknvarnonz,
                                 blockcol, blockrow, blockval, blockvars, &blockconstnnonz, blockconstcol, blockconstrow, blockconstval) );
                  assert( arraylength == nvars ); /* arraylength should alwys be sufficient */

                  matrixsize = blocksize * blocksize;

                  SCIP_CALL( SCIPallocBufferArray(scip, &fullXmatrix, matrixsize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockeigenvalues[b], blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockeigenvectors[b], matrixsize) );

                  SCIP_CALL( expandSparseMatrix(startXnblocknonz[b], blocksize, startXrow[b], startXcol[b], startXval[b], fullXmatrix) );

                  SCIP_CALL( SCIPlapackComputeEigenvectorDecomposition(SCIPbuffer(scip), blocksize, fullXmatrix, blockeigenvalues[b], blockeigenvectors[b]) );

                  /* compute coefficients for rounding problems, we get blocksize many variables corresponding to the eigenvalues of X* */
                  SCIP_CALL( SCIPallocBufferArray(scip, &obj, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &lb, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &ub, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &beg, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &ind, blocksize*nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &val, blocksize*nvars) );

                  /* initialize arrays */
                  for (i = 0; i < blocksize; i++)
                  {
                     obj[i] = 0.0;
                     beg[i] = i * nvars;

                     for (v = 0; v < nvars; v++)
                     {
                        ind[i * nvars + v] = v;
                        val[i * nvars + v] = 0.0;
                     }

                     /* make all eigenvalues non-negative, so that matrix stays positive semidefinite */
                     lb[i] = 0.0;
                     ub[i] = SCIPlpiInfinity(lpi);
                  }

                  /* iterate over constant entries to compute objective coefficients */
                  for (i = 0; i < blockconstnnonz; i++)
                  {
                     /* for every constant matrix entry (k,l) and every eigenvector i, we get an entry A_kl * V_ki *V_li
                      * entry V_ki corresponds to entry k of the i-th eigenvector, which is given as the i-th row of the eigenvectors array
                      * note that we need to mulitply by two for non-diagonal entries to also account for entry (l,k) */
                     for (evind = 0; evind < blocksize; evind++)
                     {
                        if ( blockconstrow[i] == blockconstcol[i] )
                           obj[evind] += blockconstval[i] * blockeigenvectors[b][evind * blocksize + blockconstrow[i]] * blockeigenvectors[b][evind * blocksize + blockconstcol[i]];
                        else
                           obj[evind] += 2 * blockconstval[i] * blockeigenvectors[b][evind * blocksize + blockconstrow[i]] * blockeigenvectors[b][evind * blocksize + blockconstcol[i]];
                     }
                  }

                  SCIPfreeBufferArray(scip, &blockconstval);
                  SCIPfreeBufferArray(scip, &blockconstrow);
                  SCIPfreeBufferArray(scip, &blockconstcol);

                  /* compute constraint coefficients */
                  for (v = 0; v < blocknvars; v++)
                  {
                     varind = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, blockvars[v]);
                     for (i = 0; i < blocknvarnonz[v]; i++)
                     {
                        /* for every matrix entry (k,l) and every eigenvector i, we get an entry A_kl * V_kj *V_lj
                         * entry V_kj corresponds to entry k of the j-th eigenvector, which is given as the j-th row of the eigenvectors array
                         * note that we need to mulitply by two for non-diagonal entries to also account for entry (l,k) */
                        for (evind = 0; evind < blocksize; evind++)
                        {
                           if ( blockrow[v][i] == blockcol[v][i] )
                              val[evind * nvars + varind] += blockval[v][i] * blockeigenvectors[b][evind * blocksize + blockrow[v][i]] * blockeigenvectors[b][evind * blocksize + blockcol[v][i]];
                           else
                              val[evind * nvars + varind] += 2 * blockval[v][i] * blockeigenvectors[b][evind * blocksize + blockrow[v][i]] * blockeigenvectors[b][evind * blocksize + blockcol[v][i]];
                        }
                     }
                  }

                  SCIPfreeBufferArray(scip, &blockvars);
                  SCIPfreeBufferArray(scip, &blockval);
                  SCIPfreeBufferArray(scip, &blockrow);
                  SCIPfreeBufferArray(scip, &blockcol);
                  SCIPfreeBufferArray(scip, &blocknvarnonz);
                  SCIPfreeBufferArray(scip, &fullXmatrix);

                  /* remove zero entries */
                  nremovedentries = 0;
                  for (i = 0; i < blocksize; i++)
                  {
                     beg[i] = beg[i] - nremovedentries;
                     for (v = 0; v < nvars; v++)
                     {
                        if ( REALABS(val[i * nvars + v]) < SCIPepsilon(scip) )
                        {
                           nremovedentries++;
                        }
                        else
                        {
                           val[i * nvars + v - nremovedentries] = val[i * nvars + v];
                           ind[i * nvars + v - nremovedentries] = ind[i * nvars + v];
                        }
                     }
                  }

                  SCIP_CALL( SCIPlpiAddCols(lpi, blocksize, obj, lb, ub, NULL, blocksize*nvars - nremovedentries, beg, ind, val) );

                  blocksizes[2 + b] = blocksize;
                  roundingvars += blocksize;

                  SCIPfreeBufferArray(scip, &val);
                  SCIPfreeBufferArray(scip, &ind);
                  SCIPfreeBufferArray(scip, &beg);
                  SCIPfreeBufferArray(scip, &ub);
                  SCIPfreeBufferArray(scip, &lb);
                  SCIPfreeBufferArray(scip, &obj);
               }

               /* solve the problem (since we may encounter unboundedness, which sometimes causes trouble for the CPLEX interior-point solver,
                * we use the dual Simplex solver) */
               SCIP_CALL( SCIPlpiSolveDual(lpi) );

               /* get optimal objective value of the primal rounding problem (will be -infinity if infeasible) */
               SCIP_CALL( SCIPlpiGetObjval(lpi, &primalroundobj) );

               /* if the restricted primal problem is already dual infeasible, then the original primal has to be dual infeasible as
                * well, so the dual we actually want to solve is infeasible and we can cut the node off
                * the same is true by weak duality if the restricted primal already has a larger objective value than the current cutoff-bound */
               if ( SCIPlpiIsDualInfeasible(lpi) || SCIPisGE(scip, primalroundobj, SCIPgetCutoffbound(scip)) )
               {
                  if ( SCIPlpiIsDualInfeasible(lpi) )
                  {
                     SCIPdebugMsg(scip, "Infeasibility of node %lld detected through primal rounding problem during warmstarting\n",
                           SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

                     relaxdata->roundingprobinf++;
                  }
                  else if ( SCIPisGT(scip, primalroundobj, SCIPgetCutoffbound(scip)) )
                  {
                     SCIPdebugMsg(scip, "Suboptimality of node %lld detected through primal rounding problem during warmstarting:"
                           "lower bound = %f > %f = cutoffbound\n", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), primalroundobj, SCIPgetCutoffbound(scip));

                     relaxdata->roundingcutoff++;
                  }

                  /* free memory */
                  SCIPfreeBufferArrayNull(scip, &startXval[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXcol[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXrow[nblocks]);
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIPfreeBufferArrayNull(scip,&blockeigenvectors[b]);
                     SCIPfreeBufferArrayNull(scip,&blockeigenvalues[b]);
                     SCIPfreeBufferArrayNull(scip, &startXval[b]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                     SCIPfreeBufferArrayNull(scip, &startZval[b]);
                     SCIPfreeBufferArrayNull(scip, &startZcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startZrow[b]);
                  }
                  SCIPfreeBufferArray(scip, &blocksizes);
                  SCIPfreeBufferArray(scip, &blockeigenvectors);
                  SCIPfreeBufferArray(scip, &blockeigenvalues);
                  SCIPfreeBufferArrayNull(scip, &startXval);
                  SCIPfreeBufferArrayNull(scip, &startXcol);
                  SCIPfreeBufferArrayNull(scip, &startXrow);
                  SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
                  SCIPfreeBufferArrayNull(scip, &startZval);
                  SCIPfreeBufferArrayNull(scip, &startZcol);
                  SCIPfreeBufferArrayNull(scip, &startZrow);
                  SCIPfreeBufferArrayNull(scip, &startZnblocknonz);
                  SCIPfreeBufferArray(scip, &starty);

                  TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
                  relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;

                  relaxdata->feasible = FALSE;
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               else if ( relaxdata->warmstartroundonlyinf )
               {
                  /* in this case we only cared about checking infeasibility and coldstart now */
                  SCIPfreeBufferArrayNull(scip, &startXval[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXcol[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXrow[nblocks]);
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIPfreeBufferArrayNull(scip, &blockeigenvectors[b]);
                     SCIPfreeBufferArrayNull(scip, &blockeigenvalues[b]);
                     SCIPfreeBufferArrayNull(scip, &startXval[b]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                     SCIPfreeBufferArrayNull(scip, &startZval[b]);
                     SCIPfreeBufferArrayNull(scip, &startZcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startZrow[b]);
                  }
                  SCIPfreeBufferArray(scip, &blocksizes);
                  SCIPfreeBufferArray(scip, &blockeigenvectors);
                  SCIPfreeBufferArray(scip, &blockeigenvalues);
                  SCIPfreeBufferArrayNull(scip, &startXval);
                  SCIPfreeBufferArrayNull(scip, &startXcol);
                  SCIPfreeBufferArrayNull(scip, &startXrow);
                  SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
                  SCIPfreeBufferArrayNull(scip, &startZval);
                  SCIPfreeBufferArrayNull(scip, &startZcol);
                  SCIPfreeBufferArrayNull(scip, &startZrow);
                  SCIPfreeBufferArrayNull(scip, &startZnblocknonz);
                  SCIPfreeBufferArray(scip, &starty);

                  TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
                  relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;

                  SCIP_CALL(SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, startsetting, enforceslater, timelimit));
                  goto solved;
               }
               else if ( ! SCIPlpiIsOptimal(lpi) )
               {
                  SCIPdebugMsg(scip, "Solving without warmstart since solving of the primal rounding problem failed with status %d!\n", SCIPlpiGetInternalStatus(lpi));
                  relaxdata->primalroundfails++;

                  /* since warmstart computation failed, we solve without warmstart, free memory and skip the remaining warmstarting code */
                  SCIPfreeBufferArrayNull(scip, &startXval[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXcol[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXrow[nblocks]);
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIPfreeBufferArrayNull(scip, &blockeigenvectors[b]);
                     SCIPfreeBufferArrayNull(scip, &blockeigenvalues[b]);
                     SCIPfreeBufferArrayNull(scip, &startXval[b]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                     SCIPfreeBufferArrayNull(scip, &startZval[b]);
                     SCIPfreeBufferArrayNull(scip, &startZcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startZrow[b]);
                  }
                  SCIPfreeBufferArray(scip, &blocksizes);
                  SCIPfreeBufferArray(scip, &blockeigenvectors);
                  SCIPfreeBufferArray(scip, &blockeigenvalues);
                  SCIPfreeBufferArrayNull(scip, &startXval);
                  SCIPfreeBufferArrayNull(scip, &startXcol);
                  SCIPfreeBufferArrayNull(scip, &startXrow);
                  SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
                  SCIPfreeBufferArrayNull(scip, &startZval);
                  SCIPfreeBufferArrayNull(scip, &startZcol);
                  SCIPfreeBufferArrayNull(scip, &startZrow);
                  SCIPfreeBufferArrayNull(scip, &startZnblocknonz);
                  SCIPfreeBufferArray(scip, &starty);

                  TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
                  relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;

                  SCIP_CALL(SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, startsetting, enforceslater, timelimit));
                  goto solved;
               }
               else
               {
                  SCIP_Real* optev;
                  int evpos;

                  /* the problem was solved to optimality: we construct the primal matrix using the computed eigenvalues */
                  SCIP_CALL( SCIPallocBufferArray(scip, &optev, roundingvars) );

                  SCIP_CALL( SCIPlpiGetSol(lpi, NULL, optev, NULL, NULL, NULL) );

                  /* build varbound block */
                  pos = blocksizes[1]; /* to save some sorting later, the startX arrays should start with the LP block */
                  evpos = 0;
                  for (v = 0; v < nvars; v++)
                  {
                     startXrow[nblocks][pos] = 2 * nrows + 2 * v;
                     startXcol[nblocks][pos] = 2 * nrows + 2 * v;
                     if ( ! SCIPisInfinity(scip, -1 * SCIPvarGetLbLocal(vars[v])) )
                     {
                        startXval[nblocks][pos] = optev[evpos];
                        evpos++;
                     }
                     else
                        startXval[nblocks][pos] = SCIPinfinity(scip);
                     pos++;

                     startXrow[nblocks][pos] = 2 * nrows + 2 * v + 1;
                     startXcol[nblocks][pos] = 2 * nrows + 2 * v + 1;
                     if ( ! SCIPisInfinity(scip, SCIPvarGetUbLocal(vars[v])) )
                     {

                        startXval[nblocks][pos] = optev[evpos];
                        evpos++;
                     }
                     else
                        startXval[nblocks][pos] = SCIPinfinity(scip);
                     pos++;
                  }
                  assert( evpos == blocksizes[0] );

                  /* build LP block */
                  pos = 0;
                  evpos = blocksizes[0];
                  for (r = 0; r < nrows; r++)
                  {
                     row = rows[r];

                     rowlhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
                     rowrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

                     startXrow[nblocks][pos] = 2 * r;
                     startXcol[nblocks][pos] = 2 * r;
                     if ( ! SCIPisInfinity(scip, -1 * rowlhs) )
                     {
                        startXval[nblocks][pos] = optev[evpos];
                        evpos++;
                     }
                     else
                        startXval[nblocks][pos] = SCIPinfinity(scip);
                     pos++;

                     startXrow[nblocks][pos] = 2 * r + 1;
                     startXcol[nblocks][pos] = 2 * r + 1;
                     if ( ! SCIPisInfinity(scip, rowrhs) )
                     {
                        startXval[nblocks][pos] = optev[evpos];
                        evpos++;
                     }
                     else
                        startXval[nblocks][pos] = SCIPinfinity(scip);
                     pos++;
                  }
                  assert( evpos == blocksizes[0] + blocksizes[1] );

                  startXnblocknonz[nblocks] = blocksizes[0] + blocksizes[1];

                  /* build SDP blocks */
                  pos = blocksizes[0] + blocksizes[1];
                  for (b = 0; b < nblocks; b++)
                  {
                     blocksize = blocksizes[2 + b];
                     matrixsize = blocksize * blocksize;

                     /* duplicate memory of eigenvectors to compute diag(lambda_i^*) * U^T */
                     SCIP_CALL( SCIPduplicateBufferArray(scip, &scaledeigenvectors, blockeigenvectors[b], matrixsize) );

                     /* compute diag(lambda_i_+) * U^T */
                     SCIP_CALL( scaleTransposedMatrix(blocksize, scaledeigenvectors, &(optev[pos])) );

                     /* allocate memory for full X matrix */
                     SCIP_CALL( SCIPallocBufferArray(scip, &fullXmatrix, matrixsize) );

                     /* compute U * [diag(lambda_i_+) * U^T] (note that transposes are switched because LAPACK/Fortran uses column-first-format) */
                     SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksizes[2 + b], blocksizes[2 + b], blockeigenvectors[b], TRUE, blocksizes[2 + b], blocksizes[2 + b],
                           scaledeigenvectors, FALSE, fullXmatrix) );

                     /* extract sparse matrix */
                     startXnblocknonz[b] = 0;
                     epsilon = SCIPepsilon(scip);
                     for (r = 0; r < blocksize; r++)
                     {
                        for (c = r; c < blocksize; c++)
                        {
                           matrixpos = r * blocksize + c;
                           if ( REALABS(fullXmatrix[matrixpos]) > epsilon )
                           {
                              startXrow[b][startXnblocknonz[b]] = r;
                              startXcol[b][startXnblocknonz[b]] = c;
                              startXval[b][startXnblocknonz[b]] = fullXmatrix[matrixpos];
                              startXnblocknonz[b]++;
                           }
                        }
                     }
                     pos += blocksizes[2 + b];
                     SCIPfreeBufferArray(scip, &fullXmatrix);
                     SCIPfreeBufferArray(scip, &scaledeigenvectors);
                  }
                  SCIPfreeBufferArray(scip, &optev);
               }

               /* solve dual rounding problem */

               /* allocate memory for LP and varbound-block (we take a worst-case guess here) */
               matrixsize = 2 * nvars + 2 * nrows;
               SCIP_CALL( SCIPallocBufferArray(scip, &startZrow[nblocks], matrixsize) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZcol[nblocks], matrixsize) );
               SCIP_CALL( SCIPallocBufferArray(scip, &startZval[nblocks], matrixsize) );

               /* clear the old LP */
               SCIP_CALL( SCIPlpiClear(lpi) );

               /* we want to maximize for the primal rounding problem */
               SCIP_CALL( SCIPlpiChgObjsen(lpi, SCIP_OBJSEN_MINIMIZE) );

               /* compute number of columns */
               roundingvars = nvars;
               for (b = 0; b < nblocks; b++)
                  roundingvars += blocksizes[2 + b];

               SCIP_CALL( SCIPallocBufferArray(scip, &obj, roundingvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &lb, roundingvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &ub, roundingvars) );

               for (v = 0; v < nvars; v++)
               {
                  obj[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] = SCIPvarGetObj(vars[v]);
                  lb[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] = SCIPvarGetLbLocal(vars[v]);
                  ub[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] = SCIPvarGetUbLocal(vars[v]);
               }

               /* the eigenvalue-variables should all be non-negative, but do not appear in the objective */
               for (v = nvars; v < roundingvars; v++)
               {
                  obj[v] = 0.0;
                  lb[v] = 0.0;
                  ub[v] = SCIPlpiInfinity(lpi);
               }

               /* initialize columns; rows will be added blockwise later (note that we reenter the whole problem each time
                * and do not allow warmstarts, we could store the LP- and bound-constraints between nodes (TODO, but would
                * need to check for added/removed cuts then), but warmstarting doesn't seem to make sense since the whole
                * SDP part is changed, which should form the main difficulty when solving (MI)SDPs */
               SCIP_CALL( SCIPlpiAddCols(lpi, roundingvars, obj, lb, ub, NULL, 0, NULL, NULL, NULL) );

               SCIPfreeBufferArray(scip, &ub);
               SCIPfreeBufferArray(scip, &lb);
               SCIPfreeBufferArray(scip, &obj);

               /* add LP rows */
               SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

               for (r = 0; r < nrows; r++)
               {
                  row = rows[r];
                  assert( row != NULL );
                  rownnonz = SCIProwGetNNonz(row);

                  rowvals = SCIProwGetVals(row);
                  rowcols = SCIProwGetCols(row);
                  rowlhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
                  rowrhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

                  SCIP_CALL( SCIPallocBufferArray(scip, &rowinds, rownnonz) );

                  /* iterate over rowcols and get corresponding indices */
                  for (i = 0; i < rownnonz; i++)
                     rowinds[i] = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, SCIPcolGetVar(rowcols[i]));

                  pos = 0;

                  SCIP_CALL( SCIPlpiAddRows(lpi, 1, &rowlhs, &rowrhs, NULL, rownnonz, &pos, rowinds, rowvals) );

                  SCIPfreeBufferArray(scip, &rowinds);
               }

               /* for each SDP-block add constraints linking y-variables to eigenvalues of Z matrix */
               startpos = nvars;
               for (b = 0; b < nblocks; b++)
               {
                  /* get data for this SDP block */
                  SCIP_CALL( SCIPconsSdpGetNNonz(scip, sdpblocks[b], NULL, &blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blocknvarnonz, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockcol, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockrow, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockval, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nvars) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstcol, blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstrow, blockconstnnonz) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockconstval, blockconstnnonz) );
                  blocksize = blocksizes[2 + b];
                  matrixsize = blocksize * blocksize;

                  SCIP_CALL( SCIPallocBufferArray(scip, &fullZmatrix, matrixsize) );

                  SCIP_CALL( expandSparseMatrix(startZnblocknonz[b], blocksize, startZrow[b], startZcol[b], startZval[b], fullZmatrix) );

                  SCIP_CALL( SCIPlapackComputeEigenvectorDecomposition(SCIPbuffer(scip), blocksize, fullZmatrix, blockeigenvalues[b], blockeigenvectors[b]) );

                  arraylength = nvars;
                  SCIP_CALL( SCIPconsSdpGetData(scip, sdpblocks[b], &blocknvars, &blocknnonz, &blocksize, &arraylength, blocknvarnonz,
                                 blockcol, blockrow, blockval, blockvars, &blockconstnnonz, blockconstcol, blockconstrow, blockconstval) );

                  nroundingrows = blocksize * (blocksize + 1) * 0.5;

                  SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nroundingrows) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nroundingrows) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &nblockrownonz , nroundingrows) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockrowcols , nroundingrows) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &blockrowvals , nroundingrows) );

                  /* initialize lhs and rhs with 0 and allocate memory for nonzeros (we take the worst case of nvars y entries and blocksize lambda entries) */
                  for (i = 0; i < nroundingrows; i++)
                  {
                     lhs[i] = 0.0;
                     rhs[i] = 0.0;

                     nblockrownonz[i] = 0;
                     SCIP_CALL( SCIPallocBufferArray(scip, &blockrowcols[i], nvars + blocksize) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &blockrowvals[i], nvars + blocksize) );
                  }

                  /* iterate over constant matrix to compute right-hand sides of equality constraints */
                  for (i = 0; i < blockconstnnonz; i++)
                  {
                     assert( lhs[SCIPconsSdpCompLowerTriangPos(blockconstrow[i], blockconstcol[i])] == 0.0 ); /* there should only be a single entry per index-pair */
                     lhs[SCIPconsSdpCompLowerTriangPos(blockconstrow[i], blockconstcol[i])] = blockconstval[i];
                     rhs[SCIPconsSdpCompLowerTriangPos(blockconstrow[i], blockconstcol[i])] = blockconstval[i];
                  }

                  /* iterate over all nonzeros to add them to the roundingrows */
                  for (v = 0; v < blocknvars; v++)
                  {
                     varind = SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, blockvars[v]);

                     for (i = 0; i < blocknvarnonz[v]; i++)
                     {
                        pos = SCIPconsSdpCompLowerTriangPos(blockrow[v][i], blockcol[v][i]);
                        blockrowcols[pos][nblockrownonz[pos]] = varind;
                        blockrowvals[pos][nblockrownonz[pos]] = blockval[v][i];
                        nblockrownonz[pos]++;
                     }
                  }

                  /* add entries corresponding to eigenvalues */
                  for (evind = 0; evind < blocksize; evind++)
                  {
                     for (i = 0; i < blocksize; i++)
                     {
                        for (j = 0; j <= i; j++)
                        {
                           /* for index (i,j) and every eigenvector v, we get an entry -V_iv *V_jv (we get the -1 by transferring this to the left-hand side of the equation)
                            * entry V_iv corresponds to entry i of the v-th eigenvector, which is given as the v-th row of the eigenvectors array */
                           if ( SCIPisGT(scip, REALABS(-1 * blockeigenvectors[b][evind * blocksize + i] * blockeigenvectors[b][evind * blocksize + j]), 0.0) )
                           {
                              pos = SCIPconsSdpCompLowerTriangPos(i, j);
                              blockrowcols[pos][nblockrownonz[pos]] = startpos + evind;
                              blockrowvals[pos][nblockrownonz[pos]] = -1 * blockeigenvectors[b][evind * blocksize + i] * blockeigenvectors[b][evind * blocksize + j];
                              nblockrownonz[pos]++;
                           }
                        }
                     }
                  }
                  startpos += blocksize;

                  pos = 0;
                  /* add the rows one by one */
                  for (r = 0; r < nroundingrows; r++)
                  {
                     SCIP_CALL( SCIPlpiAddRows(lpi, 1, &lhs[r], &rhs[r], NULL, nblockrownonz[r], &pos, blockrowcols[r], blockrowvals[r]) );
                     SCIPfreeBufferArray(scip, &blockrowvals[r]);
                     SCIPfreeBufferArray(scip, &blockrowcols[r]);
                  }

                  SCIPfreeBufferArray(scip, &blockrowvals);
                  SCIPfreeBufferArray(scip, &blockrowcols);
                  SCIPfreeBufferArray(scip, &nblockrownonz);
                  SCIPfreeBufferArray(scip, &rhs);
                  SCIPfreeBufferArray(scip, &lhs);
                  SCIPfreeBufferArray(scip, &fullZmatrix);
                  SCIPfreeBufferArray(scip, &blockconstval);
                  SCIPfreeBufferArray(scip, &blockconstrow);
                  SCIPfreeBufferArray(scip, &blockconstcol);
                  SCIPfreeBufferArray(scip, &blockvars);
                  SCIPfreeBufferArray(scip, &blockval);
                  SCIPfreeBufferArray(scip, &blockrow);
                  SCIPfreeBufferArray(scip, &blockcol);
                  SCIPfreeBufferArray(scip, &blocknvarnonz);
               }

               /* solve the problem (for some reason dual simplex seems to work better here) */
               SCIP_CALL( SCIPlpiSolveDual(lpi) );

               if ( ! SCIPlpiIsOptimal(lpi) )
               {
                  SCIPdebugMsg(scip, "Solution of dual rounding problem failed with status %d, continuing without warmstart\n", SCIPlpiGetInternalStatus(lpi));
                  relaxdata->dualroundfails++;

                  /* free memory */
                  SCIPfreeBufferArrayNull(scip, &startZval[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startZcol[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startZrow[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXval[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXcol[nblocks]);
                  SCIPfreeBufferArrayNull(scip, &startXrow[nblocks]);
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIPfreeBufferArrayNull(scip, &blockeigenvectors[b]);
                     SCIPfreeBufferArrayNull(scip, &blockeigenvalues[b]);
                     SCIPfreeBufferArrayNull(scip, &startZval[b]);
                     SCIPfreeBufferArrayNull(scip, &startZcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startZrow[b]);
                     SCIPfreeBufferArrayNull(scip, &startXval[b]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                  }
                  SCIPfreeBufferArray(scip, &blocksizes);
                  SCIPfreeBufferArray(scip, &blockeigenvectors);
                  SCIPfreeBufferArray(scip, &blockeigenvalues);
                  SCIPfreeBufferArrayNull(scip, &startXval);
                  SCIPfreeBufferArrayNull(scip, &startXcol);
                  SCIPfreeBufferArrayNull(scip, &startXrow);
                  SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
                  SCIPfreeBufferArrayNull(scip, &startZval);
                  SCIPfreeBufferArrayNull(scip, &startZcol);
                  SCIPfreeBufferArrayNull(scip, &startZrow);
                  SCIPfreeBufferArrayNull(scip, &startZnblocknonz);
                  SCIPfreeBufferArray(scip, &starty);

                  TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
                  relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;

                  /* since warmstart computation failed, we solve without warmstart, free memory and skip the remaining warmstarting code */
                  SCIP_CALL(SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, startsetting, enforceslater, timelimit));

                  goto solved;
               }
               else
               {
                  SCIP_Real* optev;

                  relaxdata->roundstartsuccess++;

                  /* the problem was solved to optimality: we construct the dual vector and matrix using the computed eigenvalues */
                  SCIP_CALL( SCIPallocBufferArray(scip, &optev, nvars + roundingvars) );
                  SCIP_CALL( SCIPlpiGetSol(lpi, &dualroundobj, optev, NULL, NULL, NULL) );

                  /* if the objective values of the primal and dual rounding problem agree, the problem has been solved to optimality,
                   * since both of them are respective restrictions of the original primal and dual problem */
                  if ( SCIPisEQ(scip, primalroundobj, dualroundobj) )
                  {
                     SCIP_SOL* scipsol; /* TODO: eliminate this */
                     SCIP_CONS* savedcons;

                     SCIPdebugMsg(scip, "Node %lld solved to optimality through rounding problems with optimal objective %f\n",
                           SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), dualroundobj);

                     relaxdata->roundingoptimal++;

                     TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
                     relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;

                     /* create SCIP solution (first nvars entries of optev correspond to y variables) */
                     SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
                     SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, optev) );

                     *lowerbound = dualroundobj;
                     relaxdata->objval = dualroundobj;

                     /* copy solution */
                     SCIP_CALL( SCIPsetRelaxSolValsSol(scip, scipsol, TRUE) );

                     relaxdata->feasible = TRUE;
                     *result = SCIP_SUCCESS;

                     /* save solution for warmstarts */
                     if ( relaxdata->warmstart )
                     {
                        char consname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
                        int snprintfreturn; /* this is used to assert that the SCIP string concatenation works */
#endif

#ifndef NDEBUG
                        snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "saved_relax_sol_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
                        assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether name fit into string */
#else
                        (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "saved_relax_sol_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif
                        SCIP_CALL( createConsSavesdpsol(scip, &savedcons, consname, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), scipsol,
                                 maxprimalentry, nblocks + 1, startXnblocknonz, startXrow, startXcol, startXval) );

                        SCIP_CALL( SCIPaddCons(scip, savedcons) );
                        SCIP_CALL( SCIPreleaseCons(scip, &savedcons) );
                     }

                     SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

                     /* free memory */
                     SCIPfreeBufferArray(scip, &optev);
                     SCIPfreeBufferArrayNull(scip, &startZval[nblocks]);
                     SCIPfreeBufferArrayNull(scip, &startZcol[nblocks]);
                     SCIPfreeBufferArrayNull(scip, &startZrow[nblocks]);
                     SCIPfreeBufferArrayNull(scip, &startXval[nblocks]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[nblocks]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[nblocks]);
                     for (b = 0; b < nblocks; b++)
                     {
                        SCIPfreeBufferArrayNull(scip,&blockeigenvectors[b]);
                        SCIPfreeBufferArrayNull(scip,&blockeigenvalues[b]);
                        SCIPfreeBufferArrayNull(scip, &startXval[b]);
                        SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                        SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                        SCIPfreeBufferArrayNull(scip, &startZval[b]);
                        SCIPfreeBufferArrayNull(scip, &startZcol[b]);
                        SCIPfreeBufferArrayNull(scip, &startZrow[b]);
                     }
                     SCIPfreeBufferArray(scip, &blocksizes);
                     SCIPfreeBufferArray(scip, &blockeigenvectors);
                     SCIPfreeBufferArray(scip, &blockeigenvalues);
                     SCIPfreeBufferArrayNull(scip, &startXval);
                     SCIPfreeBufferArrayNull(scip, &startXcol);
                     SCIPfreeBufferArrayNull(scip, &startXrow);
                     SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
                     SCIPfreeBufferArrayNull(scip, &startZval);
                     SCIPfreeBufferArrayNull(scip, &startZcol);
                     SCIPfreeBufferArrayNull(scip, &startZrow);
                     SCIPfreeBufferArrayNull(scip, &startZnblocknonz);
                     SCIPfreeBufferArray(scip, &starty);

                     return SCIP_OKAY;
                  }

                  /* adjust dual vector */
                  for (v = 0; v < nvars; v++)
                  {
                     if ( relaxdata->warmstartiptype == 1 )
                     {
                        /* we take a convex combination with 0, so we just scale */
                        starty[v] = (1 - relaxdata->warmstartipfactor) * optev[v];
                     }
                     else if ( relaxdata->warmstartiptype == 2 )
                     {
                        /* take the convex combination with the saved analytic center */
                        starty[v] = (1 - relaxdata->warmstartipfactor) * optev[v] + relaxdata->warmstartipfactor * SCIPgetSolVal(scip, relaxdata->ipy, vars[v]);
                     }
                  }

                  /* build Z matrix */
                  pos = nvars;
                  for (b = 0; b < nblocks; b++)
                  {
                     blocksize = blocksizes[2 + b];
                     matrixsize = blocksize * blocksize;
                     /* duplicate memory of eigenvectors to compute diag(lambda_i^*) * U^T */
                     SCIP_CALL( SCIPduplicateBufferArray(scip, &scaledeigenvectors, blockeigenvectors[b], matrixsize) );

                     /* compute diag(lambda_i_+) * U^T */
                     SCIP_CALL( scaleTransposedMatrix(blocksize, scaledeigenvectors, &(optev[pos])) );

                     /* allocate memory for full Z matrix */
                     SCIP_CALL( SCIPallocBufferArray(scip, &fullZmatrix, matrixsize) );

                     /* compute U * [diag(lambda_i_+) * U^T] (note that transposes are switched because LAPACK uses column-first-format) */
                     SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksize, blocksize, blockeigenvectors[b], TRUE, blocksize, blocksize,
                           scaledeigenvectors, FALSE, fullZmatrix) );

                     /* extract sparse matrix */
                     startZnblocknonz[b] = 0;
                     epsilon = SCIPepsilon(scip);
                     for (r = 0; r < blocksize; r++)
                     {
                        for (c = r; c < blocksize; c++)
                        {
                           matrixpos = r * blocksize + c;
                           if ( REALABS(fullZmatrix[matrixpos]) > epsilon )
                           {
                              startZrow[b][startZnblocknonz[b]] = r;
                              startZcol[b][startZnblocknonz[b]] = c;
                              startZval[b][startZnblocknonz[b]] = fullZmatrix[matrixpos];
                              startZnblocknonz[b]++;
                           }
                        }
                     }
                     pos += blocksize;
                     SCIPfreeBufferArray(scip, &fullZmatrix);
                     SCIPfreeBufferArray(scip, &scaledeigenvectors);
                     SCIPfreeBufferArray(scip, &blockeigenvectors[b]);
                     SCIPfreeBufferArray(scip, &blockeigenvalues[b]);
                  }

                  SCIPfreeBufferArray(scip, &optev);
                  SCIPfreeBufferArray(scip, &blockeigenvectors);
                  SCIPfreeBufferArray(scip, &blockeigenvalues);
                  SCIPfreeBufferArray(scip, &blocksizes);

                  /* build LP and varbound block of Z matrix */

                  /* to get a positive definite matrix, all the entries need to be strictly positive */
                  startZnblocknonz[b] = 2 * nrows + 2 * nvars;

                  /* fill LP-block */
                  for (r = 0; r < nrows; r++)
                  {
                     /* compute row value for current solution */
                     rowval = 0.0;
                     rownnonz = SCIProwGetNNonz(rows[r]);
                     rowvals = SCIProwGetVals(rows[r]);
                     rowcols = SCIProwGetCols(rows[r]);

                     for (i = 0; i < rownnonz; i++)
                        rowval += starty[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, SCIPcolGetVar(rowcols[i]))] * rowvals[i];

                     startZrow[b][2*r] = 2*r;
                     startZcol[b][2*r] = 2*r;
                     startZval[b][2*r] = rowval - (SCIProwGetLhs(rows[r]) - SCIProwGetConstant(rows[r]));

                     startZrow[b][2*r + 1] = 2*r + 1;
                     startZcol[b][2*r + 1] = 2*r + 1;
                     startZval[b][2*r + 1] = SCIProwGetRhs(rows[r]) - SCIProwGetConstant(rows[r]) - rowval;
                  }

                  /* fill varbound block */
                  for (v = 0; v < nvars; v++)
                  {
                     startZrow[b][2*nrows + 2*v] = 2*nrows + 2*v;
                     startZcol[b][2*nrows + 2*v] = 2*nrows + 2*v;
                     startZval[b][2*nrows + 2*v] = starty[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])] - SCIPvarGetLbLocal(vars[v]);

                     startZrow[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                     startZcol[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
                     startZval[b][2*nrows + 2*v + 1] = SCIPvarGetUbLocal(vars[v]) - starty[SCIPsdpVarmapperGetSdpIndex(relaxdata->varmapper, vars[v])];
                  }
               }

               TIMEOFDAY_CALL( gettimeofday(&currenttime, NULL) );/*lint !e438, !e550, !e641 */
               relaxdata->roundingprobtime += (SCIP_Real) currenttime.tv_sec + (SCIP_Real) currenttime.tv_usec / 1e6 - (SCIP_Real) starttime.tv_sec - (SCIP_Real) starttime.tv_usec / 1e6;
            }

            /* if we saved the whole primal solution before, we can set it at once */
            if ( relaxdata->warmstartprimaltype == 3 )
            {
               /* if we wanted to use the rounding problem, we already created the solution and just want to do the convex combination */
               if ( relaxdata->warmstartproject != 4 )
               {
                  SCIP_CALL( SCIPconsSavesdpsolGetPrimalMatrixNonzeros(scip, conss[parentconsind], nblocks + 1, startXnblocknonz) );

                  /* allocate sufficient memory (memory for LP-block was already allocated); we allocate an extra blocksize for adding the diagonal matrix or analytic center */
                  if ( relaxdata->warmstartproject == 3 )
                  {
                     int matrixsize;

                     /* since we cannot compute the number of nonzeros of the projection beforehand, we allocate the maximum possible (blocksize * (blocksize + 1) / 2) */
                     for (b = 0; b < nblocks; b++)
                     {
                        matrixsize = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);
                        matrixsize *= (matrixsize + 1) * 0.5;

                        SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], matrixsize) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], matrixsize) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], matrixsize) );
                     }
                  }
                  else if ( relaxdata->warmstartiptype == 1 )
                  {
                     for (b = 0; b < nblocks; b++)
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b] + SCIPconsSdpGetBlocksize(scip, sdpblocks[b])) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b] + SCIPconsSdpGetBlocksize(scip, sdpblocks[b])) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b] + SCIPconsSdpGetBlocksize(scip, sdpblocks[b])) );
                     }
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     for (b = 0; b < nblocks; b++)
                     {
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b] + relaxdata->ipXnblocknonz[b]) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b] + relaxdata->ipXnblocknonz[b]) );
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b] + relaxdata->ipXnblocknonz[b]) );
                     }
                  }

                  SCIP_CALL( SCIPconsSavesdpsolGetPrimalMatrix(scip, conss[parentconsind], nblocks + 1, startXnblocknonz, startXrow, startXcol, startXval) );
               }
            }

            /* iterate over all blocks again to compute convex combination / projection */

            /* in case of warmstartprojpdsame first compute a single maximum value for all primal and dual blocks */
            if ( relaxdata->warmstartprojpdsame && SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) )
            {
               identitydiagonal = 1.0;

               for (b = 0; b < nblocks; b++)
               {
                  for (i = 0; i < startZnblocknonz[b]; i++)
                  {
                     if ( REALABS(startZval[b][i]) > identitydiagonal )
                        identitydiagonal = REALABS(startZval[b][i]);
                  }
                  for (i = 0; i < startXnblocknonz[b]; i++)
                  {
                     if ( REALABS(startXval[b][i]) > identitydiagonal )
                        identitydiagonal = REALABS(startXval[b][i]);
                  }
               }
               identitydiagonal *= relaxdata->warmstartipfactor; /* the diagonal entries of the scaled identity matrix */
            }

            for (b = 0; b < nblocks; b++)
            {
               /* compute projection onto psd cone (computed as U * diag(lambda_i_+) * U^T where U consists of the eigenvectors of the matrix) */
               if ( relaxdata->warmstartproject == 3 )
               {
                  SCIP_Real* fullXmatrix;
                  SCIP_Real* eigenvalues;
                  SCIP_Real* eigenvectors;
                  SCIP_Real* scaledeigenvectors;
                  SCIP_Real matrixsize;
                  SCIP_Real epsilon;
                  int c;
                  int matrixpos;

                  matrixsize = blocksize * blocksize;

                  SCIP_CALL( SCIPallocBufferArray(scip, &fullXmatrix, matrixsize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &eigenvectors, matrixsize) );

                  SCIP_CALL( expandSparseMatrix(startXnblocknonz[b], blocksize, startXrow[b], startXcol[b], startXval[b], fullXmatrix) );

                  SCIP_CALL( SCIPlapackComputeEigenvectorDecomposition(SCIPbuffer(scip), blocksize, fullXmatrix, eigenvalues, eigenvectors) );

                  /* duplicate memory of eigenvectors to compute diag(lambda_i_+) * U^T */
                  SCIP_CALL( SCIPduplicateBufferArray(scip, &scaledeigenvectors, eigenvectors, matrixsize) );

                  /* set all negative eigenvalues to zero (using the property that LAPACK returns them in ascending order) */
                  i = 0;
                  while (i < blocksize && SCIPisLT(scip, eigenvalues[i], relaxdata->warmstartprojminevprimal) )
                  {
                     eigenvalues[i] = relaxdata->warmstartprojminevprimal;
                     i++;
                  }

                  /* compute diag(lambda_i_+) * U^T */
                  SCIP_CALL( scaleTransposedMatrix(blocksize, scaledeigenvectors, eigenvalues) );

                  /* compute U * [diag(lambda_i_+) * U^T] (note that transposes are switched because LAPACK uses column-first-format) */
                  SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksize, blocksize, eigenvectors, TRUE, blocksize, blocksize, scaledeigenvectors,
                        FALSE, fullXmatrix) );

                  /* extract sparse matrix from projection */
                  startXnblocknonz[b] = 0;
                  epsilon = SCIPepsilon(scip);
                  for (r = 0; r < blocksize; r++)
                  {
                     for (c = r; c < blocksize; c++)
                     {
                        matrixpos = r * blocksize + c;
                        if ( REALABS(fullXmatrix[matrixpos]) > epsilon )
                        {
                           startXrow[b][startXnblocknonz[b]] = r;
                           startXcol[b][startXnblocknonz[b]] = c;
                           startXval[b][startXnblocknonz[b]] = fullXmatrix[matrixpos];
                           startXnblocknonz[b]++;
                        }
                     }
                  }

                  /* free memory */
                  SCIPfreeBufferArray(scip, &scaledeigenvectors);
                  SCIPfreeBufferArray(scip, &eigenvectors);
                  SCIPfreeBufferArray(scip, &eigenvalues);
                  SCIPfreeBufferArray(scip, &fullXmatrix);
               }

               /* use convex combination between X and scaled identity matrix to move solution to the interior */
               if ( SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) )
               {
                  if ( relaxdata->warmstartiptype == 1 )
                  {
                     /* compute maxprimalentry (should be at least one or warmstartprojminevprimal) */
                     if ( ! relaxdata->warmstartprojpdsame )
                     {
                        if ( relaxdata->warmstartproject == 3 )
                           maxprimalentry = relaxdata->warmstartprojminevprimal;
                        else
                           maxprimalentry = 1.0;
                        for (i = 0; i < startXnblocknonz[b]; i++)
                        {
                           if ( REALABS(startXval[b][i]) > maxprimalentry )
                              maxprimalentry = REALABS(startXval[b][i]);
                        }
                        identitydiagonal = relaxdata->warmstartipfactor * maxprimalentry; /* the diagonal entries of the scaled identity matrix */
                     }
                     blocksize = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);

                     SCIP_CALL( SCIPallocBufferArray(scip, &diagentryexists, blocksize) ); /* TODO: could allocate this once for Z and X with max-blocksize */
                     for (i = 0; i < blocksize; i++)
                        diagentryexists[i] = FALSE;

                     for (i = 0; i < startXnblocknonz[b]; i++)
                     {
                        if ( startXrow[b][i] == startXcol[b][i] )
                        {
                           startXval[b][i] = startXval[b][i] * (1 - relaxdata->warmstartipfactor) + identitydiagonal; /* add identity for diagonal entries */
                           assert( startXval[b][i] >= 0.0 ); /* since the original matrix should have been positive semidefinite, diagonal entries should be >= 0 */
                           diagentryexists[startXrow[b][i]] = TRUE;
                        }
                        else
                           startXval[b][i] *= (1 - relaxdata->warmstartipfactor); /* since this is an off-diagonal entry, we scale towards zero */
                     }

                     /* if a diagonal entry was missing (because we had a zero row before), we have to add it to the end */
                     for (i = 0; i < blocksize; i++)
                     {
                        if ( ! diagentryexists[i] )
                        {
                           startXrow[b][startXnblocknonz[b]] = i;
                           startXcol[b][startXnblocknonz[b]] = i;
                           startXval[b][startXnblocknonz[b]] = identitydiagonal;
                           startXnblocknonz[b]++;
                        }
                     }
                     SCIPfreeBufferArrayNull(scip, &diagentryexists);
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     /* iterate once over all entries to multiply them with (1 - warmstartipfactor) */
                     for (i = 0; i < startXnblocknonz[b]; i++)
                        startXval[b][i] *= 1 - relaxdata->warmstartipfactor;

                     /* merge the scaled interior point array into the warmstart array */
                     SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), relaxdata->ipXrow[b], relaxdata->ipXcol[b], relaxdata->ipXval[b], relaxdata->ipXnblocknonz[b], TRUE,
                           relaxdata->warmstartipfactor, startXrow[b], startXcol[b], startXval[b], &(startXnblocknonz[b]), startXnblocknonz[b] + relaxdata->ipXnblocknonz[b]) );
                  }
               }
            }

            for (b = 0; b < nblocks; b++)
            {
               /* use convex combination between Z and scaled identity matrix to move solution to the interior */
               if ( SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) )
               {
                  if ( relaxdata->warmstartiptype == 1 )
                  {
                     if ( ! relaxdata->warmstartprojpdsame )
                     {
                        /* compute maxdualentry (should be at least one) */
                        maxdualentry = 1.0;
                        for (i = 0; i < startZnblocknonz[b]; i++)
                        {
                           if ( REALABS(startZval[b][i]) > maxdualentry )
                              maxdualentry = REALABS(startZval[b][i]);
                        }
                        identitydiagonal = relaxdata->warmstartipfactor * maxdualentry; /* the diagonal entries of the scaled identity matrix */
                     }

                     SCIP_CALL( SCIPallocBufferArray(scip, &diagentryexists, blocksize) ); /* TODO: could allocate this once for Z and X with max-blocksize */
                     for (i = 0; i < blocksize; i++)
                        diagentryexists[i] = FALSE;

                     for (i = 0; i < startZnblocknonz[b]; i++)
                     {
                        if ( startZrow[b][i] == startZcol[b][i] )
                        {
                           startZval[b][i] = startZval[b][i] * (1 - relaxdata->warmstartipfactor) + identitydiagonal; /* add identity for diagonal entries */
                           assert( startZval[b][i] >= 0 ); /* since the original matrix should have been positive semidefinite, diagonal entries should be >= 0 */
                           diagentryexists[startZrow[b][i]] = TRUE;
                        }
                        else
                           startZval[b][i] *= (1 - relaxdata->warmstartipfactor); /* since this is an off-diagonal entry, we scale towards zero */
                     }

                     /* if a diagonal entry was missing (because we had a zero row before), we have to add it to the end */
                     for (i = 0; i < blocksize; i++)
                     {
                        if ( ! diagentryexists[i] )
                        {
                           startZrow[b][startZnblocknonz[b]] = i;
                           startZcol[b][startZnblocknonz[b]] = i;
                           startZval[b][startZnblocknonz[b]] = identitydiagonal;
                           startZnblocknonz[b]++;
                        }
                     }

                     SCIPfreeBufferArrayNull(scip, &diagentryexists);
                  }
                  else if ( relaxdata->warmstartiptype == 2 )
                  {
                     /* iterate once over all entries to multiply them with (1 - warmstartipfactor) */
                     for (i = 0; i < startZnblocknonz[b]; i++)
                        startZval[b][i] *= 1 - relaxdata->warmstartipfactor;

                     /* merge the scaled interior point array into the warmstart array */
                     SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), relaxdata->ipZrow[b], relaxdata->ipZcol[b], relaxdata->ipZval[b], relaxdata->ipZnblocknonz[b], TRUE,
                           relaxdata->warmstartipfactor, startZrow[b], startZcol[b], startZval[b], &(startZnblocknonz[b]), startZnblocknonz[b] + relaxdata->ipZnblocknonz[b]) );
                  }
               }
            }

            /* compute projection for LP block */
            if ( relaxdata->warmstartproject == 3 )
            {
               int nsavedentries;
               int lastentry;
               int j;

               /* sort indices by row/col; TODO: check if this is necessary */
               SCIPsortIntIntReal(startXrow[nblocks], startXcol[nblocks], startXval[nblocks], startXnblocknonz[nblocks]);

               /* iterate over all entries */
               nsavedentries = startXnblocknonz[nblocks];
               lastentry = 0;

               for (i = 0; i < nsavedentries; i++)
               {
                  assert( startXrow[nblocks][i] == startXcol[nblocks][i] ); /* this is the LP-block */
                  /* if some entries are missing, we add them at the end */
                  for (j = lastentry + 1; j < startXrow[nblocks][i]; j++)
                  {
                     assert( startXnblocknonz[nblocks] < 2 * nrows + 2 * nvars );
                     startXrow[nblocks][startXnblocknonz[nblocks]] = j;
                     startXcol[nblocks][startXnblocknonz[nblocks]] = j;
                     startXval[nblocks][startXnblocknonz[nblocks]] = relaxdata->warmstartprojminevprimal;
                     startXnblocknonz[nblocks]++;
                  }
                  if ( SCIPisLT(scip, startXval[b][i], 1.0) )
                     startXval[b][i] = relaxdata->warmstartprojminevprimal;

                  lastentry = startXrow[nblocks][i];
               }
               /* add missing entries at the end */
               for (j = lastentry + 1; j < 2 * nrows + 2 * nvars; j++)
               {
                  assert( startXnblocknonz[nblocks] < 2 * nrows + 2 * nvars );
                  startXrow[nblocks][startXnblocknonz[nblocks]] = j;
                  startXcol[nblocks][startXnblocknonz[nblocks]] = j;
                  startXval[nblocks][startXnblocknonz[nblocks]] = relaxdata->warmstartprojminevprimal;
                  startXnblocknonz[nblocks]++;
               }
            }

            /* take convex combination for LP block */
            if ( SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) )
            {
               if ( relaxdata->warmstartiptype == 1 )
               {
                  int nsavedentries;
                  int lastentry;
                  int j;

                  /* if warmstartprojpdsame is true we use the computed value, otherwise we use the maximum of this block, which is one */
                  if ( ! relaxdata->warmstartprojpdsame )
                     identitydiagonal = relaxdata->warmstartipfactor;

                  /* sort indices by row/col; TODO: check if this is necessary */
                  SCIPsortIntIntReal(startXrow[nblocks], startXcol[nblocks], startXval[nblocks], startXnblocknonz[nblocks]);

                  /* iterate over all entries */
                  nsavedentries = startXnblocknonz[nblocks];
                  lastentry = 0;

                  for (i = 0; i < nsavedentries; i++)
                  {
                     assert( startXrow[nblocks][i] == startXcol[nblocks][i] ); /* this is the LP-block */
                     /* if some entries are missing, we add them at the end */
                     for (j = lastentry + 1; j < startXrow[nblocks][i]; j++)
                     {
                        assert( startXnblocknonz[nblocks] < 2 * nrows + 2 * nvars );
                        startXrow[nblocks][startXnblocknonz[nblocks]] = j;
                        startXcol[nblocks][startXnblocknonz[nblocks]] = j;
                        /* if warmstartprojpdsame is true we use the computed value, otherwise we use the maximum of this block, which is one */
                        if ( relaxdata->warmstartprojpdsame )
                           startXval[nblocks][startXnblocknonz[nblocks]] = identitydiagonal;
                        else
                           startXval[nblocks][startXnblocknonz[nblocks]] = relaxdata->warmstartipfactor;
                        startXnblocknonz[nblocks]++;
                     }
                     /* we only take the convex combination if the value is less than one, since the maxblockentry is equal to the value
                      * otherwise, so taking the convex combination doesn't change anything in that case (unless warmstarprojpdsame)
                      */
                     if ( relaxdata->warmstartprojpdsame )
                        startXval[b][i] = (1 - relaxdata->warmstartipfactor) * startXval[b][i] + identitydiagonal;
                     else if ( SCIPisLT(scip, startXval[b][i], 1.0) )
                        startXval[b][i] = (1 - relaxdata->warmstartipfactor) * startXval[b][i] + relaxdata->warmstartipfactor;

                     lastentry = startXrow[nblocks][i];
                  }
                  /* add missing entries at the end */
                  for (j = lastentry + 1; j < 2 * nrows + 2 * nvars; j++)
                  {
                     assert( startXnblocknonz[nblocks] < 2 * nrows + 2 * nvars );
                     startXrow[nblocks][startXnblocknonz[nblocks]] = j;
                     startXcol[nblocks][startXnblocknonz[nblocks]] = j;
                     startXval[nblocks][startXnblocknonz[nblocks]] = identitydiagonal;
                     startXnblocknonz[nblocks]++;
                  }
               }
               else if ( relaxdata->warmstartiptype == 2 )
               {
                  /* iterate once over all entries to multiply them with (1 - warmstartipfactor) */
                  for (i = 0; i < startXnblocknonz[nblocks]; i++)
                     startXval[nblocks][i] *= 1 - relaxdata->warmstartipfactor;

                  /* merge the scaled interior point array into the warmstart array */
                  SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), relaxdata->ipXrow[nblocks],
                        relaxdata->ipXcol[nblocks], relaxdata->ipXval[nblocks], relaxdata->ipXnblocknonz[nblocks], TRUE,
                        relaxdata->warmstartipfactor, startXrow[nblocks], startXcol[nblocks], startXval[nblocks],
                        &(startXnblocknonz[nblocks]), startXnblocknonz[nblocks] + relaxdata->ipXnblocknonz[nblocks]) );
               }

               /* take the convex combination for the dual (in this case we do not need to check for missing entries since we added all of them ourselves) */
               if ( relaxdata->warmstartiptype == 1 )
               {
                  for (r = 0; r < nrows; r++)
                  {
                     /* for the project we just set all values smaller than minev to minev */
                     if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[b][2*r], relaxdata->warmstartprojminevdual) )
                        startZval[nblocks][2*r] = relaxdata->warmstartprojminevdual;
                     /* we only take the convex combination if the value is less than one, since the maxblockentry is equal to the value
                      * otherwise, so taking the convex combination doesn't change anything in that case (unless projpdsame)
                      */
                     else if ( relaxdata->warmstartiptype == 1 && (SCIPisLT(scip, startZval[nblocks][2*r], 1.0) || relaxdata->warmstartprojpdsame) )
                     {
                        /* since we want the value to be strictly positive, if the original entry is negative we just set it to identitydiagonal */
                        if ( SCIPisLT(scip, startZval[nblocks][2*r], 0.0) )
                           startZval[nblocks][2*r] = identitydiagonal;
                        else
                           startZval[nblocks][2*r] = (1 - relaxdata->warmstartipfactor) * startZval[nblocks][2*r] + identitydiagonal;
                     }

                     if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[nblocks][2*r + 1], relaxdata->warmstartprojminevdual) )
                        startZval[nblocks][2*r + 1] = relaxdata->warmstartprojminevdual;
                     else if ( relaxdata->warmstartiptype == 1 && (SCIPisLT(scip, startZval[nblocks][2*r + 1], 1.0) || relaxdata->warmstartprojpdsame) )
                     {
                        /* since we want the value to be strictly positive, if the original entry is negative we just set it to identitydiagonal */
                        if ( SCIPisLT(scip, startZval[nblocks][2*r + 1], 0.0) )
                           startZval[nblocks][2*r + 1] = identitydiagonal;
                        else
                           startZval[nblocks][2*r + 1] = (1 - relaxdata->warmstartipfactor) * startZval[nblocks][2*r + 1] + identitydiagonal;
                     }
                  }

                  for (v = 0; v < nvars; v++)
                  {
                     if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v], relaxdata->warmstartprojminevdual) )
                        startZval[nblocks][2*nrows + 2*v] = relaxdata->warmstartprojminevdual;
                     else if ( relaxdata->warmstartiptype == 1 && (SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v], 1.0) || relaxdata->warmstartprojpdsame) )
                     {
                        /* since we want the value to be strictly positive, if the original entry is negative we just set it to identitydiagonal */
                        if ( SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v], 0.0) )
                           startZval[nblocks][2*nrows + 2*v] = identitydiagonal;
                        else
                           startZval[nblocks][2*nrows + 2*v] = (1 - relaxdata->warmstartipfactor) * startZval[nblocks][2*nrows + 2*v] + identitydiagonal;
                     }

                     if ( relaxdata->warmstartiptype == 1 && relaxdata->warmstartproject == 3 && SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v + 1], relaxdata->warmstartprojminevdual) )
                        startZval[nblocks][2*nrows + 2*v + 1] = relaxdata->warmstartprojminevdual;
                     else if ( relaxdata->warmstartiptype == 1 && (SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v + 1], 1.0) || relaxdata->warmstartprojpdsame) )
                     {
                        /* since we want the value to be strictly positive, if the original entry is negative we just set it to identitydiagonal */
                        if ( SCIPisLT(scip, startZval[nblocks][2*nrows + 2*v + 1], 0.0) )
                           startZval[nblocks][2*nrows + 2*v + 1] = identitydiagonal;
                        else
                           startZval[nblocks][2*nrows + 2*v + 1] = (1 - relaxdata->warmstartipfactor) * startZval[nblocks][2*nrows + 2*v + 1] + identitydiagonal;
                     }
                  }
               }
               else if ( relaxdata->warmstartiptype == 2 )
               {
                  /* iterate once over all entries to multiply them with (1 - warmstartipfactor) */
                  for (i = 0; i < startZnblocknonz[nblocks]; i++)
                     startZval[nblocks][i] *= 1 - relaxdata->warmstartipfactor;

                  /* merge the scaled interior point array into the warmstart array */
                  SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), relaxdata->ipZrow[nblocks],
                        relaxdata->ipZcol[nblocks], relaxdata->ipZval[nblocks], relaxdata->ipZnblocknonz[nblocks], TRUE,
                        relaxdata->warmstartipfactor, startZrow[nblocks], startZcol[nblocks], startZval[nblocks],
                        &(startZnblocknonz[nblocks]), startZnblocknonz[nblocks] + relaxdata->ipZnblocknonz[nblocks]) );
               }
            }
         }

#ifdef SCIP_PRINT_WARMSTART
         SCIPdebugMessage("warmstart using the following point:\n");
         nblocks = SCIPconshdlrGetNConss(SCIPfindConshdlr(scip, "SDP"));
         for (i = 0; i < nvars; i++)
            SCIPdebugMessage("y[%d]=%f\n", i, starty[i]);

         if ( SCIPsdpiDoesWarmstartNeedPrimal() )
         {
            for (b = 0; b < nblocks + 1; b++)
            {
               SCIPdebugMessage("dual block %d\n", b);
               for (i = 0; i < startZnblocknonz[b]; i++)
               {
                  SCIPdebugMessage("Z(%d,%d)=%f\n", startZrow[b][i], startZcol[b][i], startZval[b][i]);
               }
            }
            for (b = 0; b < nblocks + 1; b++)
            {
               SCIPdebugMessage("primal block %d\n", b);
               for (i = 0; i < startXnblocknonz[b]; i++)
               {
                  SCIPdebugMessage("X(%d,%d)=%f\n", startXrow[b][i], startXcol[b][i], startXval[b][i]);
               }
            }
         }
#endif

         /* solve with given starting point */
         SCIP_CALL(SCIPsdpiSolve(sdpi, starty, startZnblocknonz, startZrow, startZcol, startZval, startXnblocknonz, startXrow,
                  startXcol, startXval, startsetting, enforceslater, timelimit));

         if ( SCIPsdpiDoesWarmstartNeedPrimal() )
         {
            SCIP_CONSHDLR* sdpconshdlr;

            sdpconshdlr = SCIPfindConshdlr(scip, "SDP");
            nblocks = SCIPconshdlrGetNConss(sdpconshdlr) + 1; /* +1 for LP block */

            /* free memory */
            for (b = 0; b < nblocks; b++)
            {
               SCIPfreeBufferArray(scip, &startXval[b]);
               SCIPfreeBufferArray(scip, &startXcol[b]);
               SCIPfreeBufferArray(scip, &startXrow[b]);
               SCIPfreeBufferArray(scip, &startZval[b]);
               SCIPfreeBufferArray(scip, &startZcol[b]);
               SCIPfreeBufferArray(scip, &startZrow[b]);
            }

            SCIPfreeBufferArray(scip, &startXval);
            SCIPfreeBufferArray(scip, &startXcol);
            SCIPfreeBufferArray(scip, &startXrow);
            SCIPfreeBufferArray(scip, &startXnblocknonz);
            SCIPfreeBufferArray(scip, &startZval);
            SCIPfreeBufferArray(scip, &startZcol);
            SCIPfreeBufferArray(scip, &startZrow);
            SCIPfreeBufferArray(scip, &startZnblocknonz);
         }
         SCIPfreeBufferArray(scip, &starty);
      }
   }

   solved:

   relaxdata->lastsdpnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* update calls, iterations and stability numbers (only if the SDP-solver was actually called) */
   relaxdata->sdpinterfacecalls++;
   naddedsdpcalls = 0;
   SCIP_CALL( SCIPsdpiGetSdpCalls(relaxdata->sdpi, &naddedsdpcalls) );
   usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   if ( naddedsdpcalls )
   {
      relaxdata->sdpcalls += naddedsdpcalls;
      naddediters = 0;
      SCIP_CALL( SCIPsdpiGetIterations(relaxdata->sdpi, &naddediters) );
      relaxdata->sdpiterations += naddediters;

      SCIP_CALL( SCIPsdpiSettingsUsed(relaxdata->sdpi, &usedsetting) );

      switch( usedsetting )/*lint --e{788}*/
      {
      case SCIP_SDPSOLVERSETTING_PENALTY:
         relaxdata->solvedpenalty++;
         break;
      case SCIP_SDPSOLVERSETTING_FAST:
         relaxdata->solvedfast++;
         break;
      case SCIP_SDPSOLVERSETTING_MEDIUM:
         relaxdata->solvedmedium++;
         break;
      case SCIP_SDPSOLVERSETTING_STABLE:
         relaxdata->solvedstable++;
         break;
      case SCIP_SDPSOLVERSETTING_UNSOLVED:
         relaxdata->unsolved++;
         break;
      default:
         break;
      }
      primalslater = SCIP_SDPSLATER_NOINFO;
      dualslater = SCIP_SDPSLATER_NOINFO;
      SCIP_CALL( SCIPsdpiSlater(relaxdata->sdpi, &primalslater, &dualslater) );
      switch( primalslater )/*lint --e{788}*/
      {
         case SCIP_SDPSLATER_NOINFO:
            relaxdata->npslatercheckfailed++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
            }
            break;
         case SCIP_SDPSLATER_NOT:
            relaxdata->npnoslater++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nnoslater++;
                  break;
            }
            break;
         case SCIP_SDPSLATER_HOLDS:
            relaxdata->npslaterholds++;
            switch( dualslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->ndnoslater++;
                  relaxdata->nnoslater++;
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->ndslaterholds++;
                  relaxdata->nslaterholds++;
                  break;
               case SCIP_SDPSLATER_INF:
                  relaxdata->nslaterinfeasible++;
                  break;
               default:
                  relaxdata->ndslatercheckfailed++;
                  relaxdata->nslatercheckfailed++;
                  break;
            }
            break;
            default:
               relaxdata->npslatercheckfailed++;
               relaxdata->ndslatercheckfailed++;
               relaxdata->nslatercheckfailed++;
               break;
      }
      slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
      SCIP_CALL( SCIPsdpiSlaterSettings(relaxdata->sdpi, &slatersetting) );
      switch( slatersetting )/*lint --e{788}*/
      {
         case SCIP_SDPSLATERSETTING_STABLEWSLATER:
            relaxdata->stablewslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLEWSLATER:
            relaxdata->unstablewslater++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYWSLATER:
            relaxdata->penaltywslater++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDWSLATER:
            relaxdata->boundedwslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER:
            relaxdata->unsolvedwslater++;
            break;
         case SCIP_SDPSLATERSETTING_STABLENOSLATER:
            relaxdata->stablenoslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLENOSLATER:
            relaxdata->unstablenoslater++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYNOSLATER:
            relaxdata->penaltynoslater++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER:
            relaxdata->boundednoslater++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER:
            relaxdata->unsolvednoslater++;
            break;
         case SCIP_SDPSLATERSETTING_STABLEINFEASIBLE:
            relaxdata->stableinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE:
            relaxdata->unstableinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE:
            relaxdata->penaltyinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE:
            relaxdata->boundedinfeasible++;
            break;
         case SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE:
            relaxdata->unsolvedinfeasible++;
            break;
         default:
            break;
      }
   }

   /* remember settings */
   if ( ! (strcmp(SCIPsdpiGetSolverName(), "DSDP") == 0) && ! (strstr(SCIPsdpiGetSolverName(), "MOSEK") != NULL) )
   {
      (void) SCIPsnprintf(saveconsname, SCIP_MAXSTRLEN, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
      SCIP_CALL( createConsSavedsdpsettings(scip, &savedsetting, saveconsname, usedsetting) );
      SCIP_CALL( SCIPaddCons(scip, savedsetting) );
      SCIP_CALL( SCIPreleaseCons(scip, &savedsetting) );
   }

   if ( ! SCIPsdpiWasSolved(sdpi) )
      relaxdata->feasible = FALSE;

   if ( SCIPinProbing(scip) )
      relaxdata->probingsolved = SCIPsdpiWasSolved(sdpi);
   else
      relaxdata->origsolved = SCIPsdpiSolvedOrig(sdpi);

   if ( SCIPsdpiIsAcceptable(sdpi) )
   {
#ifdef SCIP_MORE_DEBUG /* print the optimal solution */
      {
         int sollength;
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         sollength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &sollength) ); /* get both the objective and the solution from the SDP solver */

         assert( sollength == nvars ); /* If this isn't true any longer, the getSol-call was unsuccessfull, because the given array wasn't long enough,
                                        * but this can't happen, because the array has enough space for all SDP variables. */

         if ( SCIPsdpiFeasibilityKnown(sdpi) )
         {
            SCIPdebugMessage("optimal solution: objective = %f, dual feasible: %u, primal feasible: %u.\n",
                  objforscip, SCIPsdpiIsDualFeasible(sdpi), SCIPsdpiIsPrimalFeasible(sdpi));
         }
         else
         {
            SCIPdebugMessage("The solver could not determine feasibility ! ");
         }

         /* output solution */
         for (i = 0; i < nvars; ++i)
         {
            SCIPdebugMessage("<%s> = %f\n", SCIPvarGetName(vars[i]), solforscip[i]);
         }
         SCIPfreeBufferArray(scip, &solforscip);
      }
#endif

      if ( SCIPsdpiIsDualInfeasible(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to infeasibility.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsObjlimExc(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to objective limit.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsDualUnbounded(sdpi) )
      {
         SCIPdebugMessage("Node unbounded.");
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsPrimalFeasible(sdpi) && SCIPsdpiIsDualFeasible(sdpi) )
      {
         SCIP_SOL* scipsol; /* TODO: eliminate this from warmstart and save array instead */
         SCIP_SOL* preoptimalsol;
         SCIP_CONS* savedcons;
         int slength;
         SCIP_Bool preoptimalsolsuccess;

         /* get solution w.r.t. SCIP variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         slength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

         assert( slength == nvars ); /* If this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                      * but this can't happen, because the array has enough space for all sdp variables. */

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         *lowerbound = objforscip;
         relaxdata->objval = objforscip;

         /* copy solution */
         SCIP_CALL( SCIPsetRelaxSolVals(scip, nvars, vars, solforscip, TRUE) );

         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         preoptimalsolsuccess = FALSE;

         /* save solution for warmstarts (only if we did not use the penalty formulation, since this would change the problem structure) */
         if ( relaxdata->warmstart && SCIPsdpiSolvedOrig(relaxdata->sdpi) )
         {
            SCIP_Real maxprimalentry;
            int* startXnblocknonz;
            int** startXrow;
            int** startXcol;
            SCIP_Real** startXval;
            char consname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
            int snprintfreturn; /* this is used to assert that the SCIP string concatenation works */
#endif

            startXnblocknonz = NULL;
            startXrow = NULL;
            startXcol = NULL;
            startXval = NULL;

            /* use preoptimal solution if using DSDP and parameter is set accordingly */
            if ( relaxdata->warmstartpreoptsol )
            {
               if ( strcmp(SCIPsdpiGetSolverName(), "DSDP") == 0.0 || strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0.0 )
               {
                  SCIP_Real* preoptimalvec;
                  int nvarsgiven;

                  SCIP_CALL( SCIPallocBufferArray(scip, &preoptimalvec, nvars) );
                  nvarsgiven = nvars;

                  if ( SCIPsdpiDoesWarmstartNeedPrimal() )
                  {
                     maxprimalentry = 0.0;
                     if ( relaxdata->warmstartprimaltype == 3 )
                     {
                        nblocks = SCIPconshdlrGetNConss(SCIPfindConshdlr(scip, "SDP")) + 1; /* +1 for the LP part */
                        SCIP_CALL( SCIPallocBufferArray(scip, &startXnblocknonz, nblocks) );

                        /* get amount of memory to allocate for row/col/val from sdpi */
                        SCIP_CALL( SCIPsdpiGetPreoptimalPrimalNonzeros(relaxdata->sdpi, nblocks, startXnblocknonz) );

                        /* check if the primal matrix exists, otherwise skip creation of the savedsol contraint */
                        if ( startXnblocknonz[0] > -1 )
                        {
                           preoptimalsolsuccess = TRUE;

                           SCIP_CALL( SCIPallocBufferArray(scip, &startXrow, nblocks) );
                           SCIP_CALL( SCIPallocBufferArray(scip, &startXcol, nblocks) );
                           SCIP_CALL( SCIPallocBufferArray(scip, &startXval, nblocks) );

                           /* allocate memory for different blocks in row/col/val */
                           for (b = 0; b < nblocks; b++)
                           {
                              SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b]) );
                              SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b]) );
                              SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b]) );
                           }

                           SCIP_CALL( SCIPsdpiGetPreoptimalSol(relaxdata->sdpi, &preoptimalsolsuccess, preoptimalvec, &nvarsgiven,
                                 nblocks, startXnblocknonz, startXrow, startXcol, startXval) );
                        }
                        else
                           preoptimalsolsuccess = FALSE;
                     }
                     else
                     {
                        nblocks = 0;
                        maxprimalentry = SCIPsdpiGetMaxPrimalEntry(relaxdata->sdpi);
                     }
                  }
                  else
                  {
                     maxprimalentry = 0.0;
                     SCIP_CALL( SCIPsdpiGetPreoptimalSol(relaxdata->sdpi, &preoptimalsolsuccess, preoptimalvec, &nvarsgiven,
                           -1, NULL, NULL, NULL, NULL) );
                     nblocks = 0;
                  }

                  if ( preoptimalsolsuccess )
                  {
                     assert( nvarsgiven == nvars ); /* length of solution should be nvars */

                     /* create SCIP solution */
                     SCIP_CALL( SCIPcreateSol(scip, &preoptimalsol, NULL) );
                     SCIP_CALL( SCIPsetSolVals(scip, preoptimalsol, nvars, vars, preoptimalvec) );
                  }

                  SCIPfreeBufferArray(scip, &preoptimalvec);
               }
               else
               {
                  SCIPerrorMessage("Warmstarting with preoptimal solutions currently only supported for DSDP and SDPA \n");
                  return SCIP_LPERROR;
               }
            }
            else if ( SCIPsdpiDoesWarmstartNeedPrimal() )
            {
               maxprimalentry = 0.0;
               if ( relaxdata->warmstartprimaltype == 3 )
               {
                  nblocks = SCIPconshdlrGetNConss(SCIPfindConshdlr(scip, "SDP")) + 1; /* +1 for the LP part */
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXnblocknonz, nblocks) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXrow, nblocks) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXcol, nblocks) );
                  SCIP_CALL( SCIPallocBufferArray(scip, &startXval, nblocks) );

                  /* get amount of memory to allocate for row/col/val from sdpi */
                  SCIP_CALL( SCIPsdpiGetPrimalNonzeros(relaxdata->sdpi, nblocks, startXnblocknonz) );

                  /* allocate memory for different blocks in row/col/val */
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXrow[b], startXnblocknonz[b]) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXcol[b], startXnblocknonz[b]) );
                     SCIP_CALL( SCIPallocBufferArray(scip, &startXval[b], startXnblocknonz[b]) );
                  }
                  SCIP_CALL( SCIPsdpiGetPrimalMatrix(sdpi, nblocks, startXnblocknonz, startXrow, startXcol, startXval) );
               }
               else
               {
                  nblocks = 0;
                  maxprimalentry = SCIPsdpiGetMaxPrimalEntry(relaxdata->sdpi);
               }
            }
            else
            {
               maxprimalentry = 0.0;
               nblocks = 0;
            }
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "saved_relax_sol_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
            assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether name fit into string */
#else
   (void) SCIPsnprintf(consname, SCIP_MAXSTRLEN, "saved_relax_sol_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif
            if ( relaxdata->warmstartpreoptsol )
            {
               /* only create constraint if the preoptimal solution exists, otherwise we don't want to warmstart at all */
               if ( preoptimalsolsuccess )
               {
                  SCIP_CALL( createConsSavesdpsol(scip, &savedcons, consname, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), preoptimalsol,
                        maxprimalentry, nblocks, startXnblocknonz, startXrow, startXcol, startXval) );

                  SCIP_CALL( SCIPaddCons(scip, savedcons) );
                  SCIP_CALL( SCIPreleaseCons(scip, &savedcons) );
               }
            }
            else
            {
               SCIP_CALL( createConsSavesdpsol(scip, &savedcons, consname, SCIPnodeGetNumber(SCIPgetCurrentNode(scip)), scipsol,
                     maxprimalentry, nblocks, startXnblocknonz, startXrow, startXcol, startXval) );

               SCIP_CALL( SCIPaddCons(scip, savedcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &savedcons) );
            }

            if ( SCIPsdpiDoesWarmstartNeedPrimal() && relaxdata->warmstartprimaltype == 3 )
            {
               /* free memory for primal matrix */
               if ( startXnblocknonz[0] > 1 ) /* no memory was allocated if computation of preoptimal solution failed */
               {
                  for (b = 0; b < nblocks; b++)
                  {
                     SCIPfreeBufferArrayNull(scip, &startXval[b]);
                     SCIPfreeBufferArrayNull(scip, &startXcol[b]);
                     SCIPfreeBufferArrayNull(scip, &startXrow[b]);
                  }
                  SCIPfreeBufferArrayNull(scip, &startXval);
                  SCIPfreeBufferArrayNull(scip, &startXcol);
                  SCIPfreeBufferArrayNull(scip, &startXrow);
               }
               SCIPfreeBufferArrayNull(scip, &startXnblocknonz);
            }
         }

         SCIPfreeBufferArray(scip, &solforscip);
         SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
         if ( preoptimalsolsuccess )
         {
            SCIP_CALL( SCIPfreeSol(scip, &preoptimalsol) );
         }
      }
   }
   else
   {
      SCIP_Real objlb;

      if ( SCIPsdpiIsTimelimExc(relaxdata->sdpi) )
      {
         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      /* if we used the penalty approach, we might have calculated a good lower bound, even if we did not produce a feasible solution, otherwise we
       * keep the current bound, if the current bound is -infty, we abort */
      objlb = -SCIPinfinity(scip);
      SCIP_CALL( SCIPsdpiGetLowerObjbound(relaxdata->sdpi, &objlb) );
      if ( ! SCIPisInfinity(scip, objlb) )
      {
         *lowerbound = objlb;
         SCIPdebugMessage("The relaxation could not be solved, using best computed bound from penalty formulation.\n");
      }
      else if ( ! SCIPisInfinity(scip, -1 * SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip))) )
      {
         *lowerbound = SCIPnodeGetLowerbound(SCIPgetCurrentNode(scip));
         SCIPdebugMessage("The relaxation could not be solved, keeping old bound.\n");
      }
      else
      {
         *result = SCIP_SUSPENDED;
         SCIPerrorMessage("The relaxation of the root node could not be solved, so there is no hope to solve this instance.\n");
         return SCIP_ERROR;
      }

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** checks whether all variables are fixed */
static
SCIP_Bool allVarsFixed(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR** vars;
   int i;

   assert( scip != NULL );

   vars = SCIPgetVars(scip);

   /* try to find a variable that is not fixed */
   for (i = 0; i < SCIPgetNVars(scip); i++)
   {
      if ( SCIPisLT(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
         return FALSE;
   }

   /* if no variable with lower bound strictly lower than upper bound has been found, all variables are fixed */
   return TRUE;
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecSdp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_VAR** vars;
   SCIP_Real* ubs;
   SCIP_Bool cutoff;
   SCIP_SOL* scipsol; /* TODO: eliminate this */
   int nconss;
   int nvars;
   int i;
#ifdef SCIP_EVEN_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
#endif

   SCIPdebugMessage("Calling relaxExecSdp.\n");

   relaxdata = SCIPrelaxGetData(relax);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* don't run again if we already solved the current node (except during probing), and we solved the correct problem */
   if ( (relaxdata->lastsdpnode == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) && ( ! SCIPinProbing(scip) ) ) && relaxdata->origsolved && ! relaxdata->resolve )
   {
      SCIP_COL** cols;
      SCIP_Real objforscip;
      SCIP_Real* solforscip;
      int ncols;
      int slength;

      SCIPdebugMessage("Already solved SDP-relaxation for node %ld, returning with SCIP_SUCCESS so that no other relaxator is called.\n",
            SCIPrelaxGetData(relax)->lastsdpnode);

      if ( SCIPsdpiIsDualUnbounded(relaxdata->sdpi) )
      {
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }

      /* get solution w.r.t. SCIP variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
      slength = nvars;
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

      assert( slength == nvars ); /* If this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                   * but this can't happen, because the array has enough space for all sdp variables. */

      /* create SCIP solution */
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetRelaxSolVals(scip, nvars, vars, solforscip, TRUE) );

      *lowerbound = objforscip;

      /* copy solution */
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      for (i = 0; i < ncols; i++)
      {
         SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );
      }

      SCIP_CALL( SCIPmarkRelaxSolValid(scip, TRUE) );
      *result = SCIP_SUCCESS;

      SCIPfreeBufferArray(scip, &solforscip);
      SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* if we are solving a probing SDP, remember that we didn't solve the original problem */
   relaxdata->origsolved = FALSE;

   /* construct the lp and make sure, that everything is where it should be */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   if ( cutoff )
   {
      relaxdata->feasible = FALSE;
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* very important to call flushLP */
   SCIP_CALL( SCIPflushLP(scip) );

   /* get varmapper */
   nconss = SCIPgetNConss(scip);

#ifdef SCIP_EVEN_MORE_DEBUG
   for (i = 0; i < nvarsfordebug; i++)
   {
      SCIPdebugMessage("variable %s: status = %u, integral = %u, bounds = [%f, %f] \n", SCIPvarGetName(varsfordebug[i]), SCIPvarGetStatus(varsfordebug[i]),
         SCIPvarIsIntegral(varsfordebug[i]), SCIPvarGetLbLocal(varsfordebug[i]), SCIPvarGetUbLocal(varsfordebug[i]));
   }
#endif

   if ( nconss == 0 )
   {
      /* if there are no constraints, there is nothing to do */
      relaxdata->feasible = TRUE;
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( allVarsFixed(scip) )
   {
      SCIP_Bool feasible;

      /* if all variables, really all, are fixed, we give this fixed solution to SCIP */

      SCIP_CALL( SCIPallocBufferArray(scip, &ubs, nvars) );

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
      }

      SCIPdebugMessage("EVERYTHING IS FIXED, objective value = %f\n", *lowerbound);

      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      /* set the relaxation solution */
      for (i = 0; i < nvars; i++)
      {
         SCIP_CALL( SCIPsetRelaxSolVal(scip, vars[i], SCIPvarGetLbLocal(vars[i])) );
      }
      SCIP_CALL( SCIPmarkRelaxSolValid(scip, TRUE) );

      /* check if the solution really is feasible */
      SCIP_CALL( SCIPcheckSol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );

      relaxdata->feasible = feasible;

      SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

      SCIPfreeBufferArray(scip, &ubs);

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* update LP Data in Interface */
   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, TRUE, TRUE) );

   SCIP_CALL( calcRelax(scip, relaxdata, result, lowerbound));

   return SCIP_OKAY;
}


/** this method is called after presolving is finished, at this point the varmapper is prepared and the SDP Interface is initialized and gets
 *  the SDP information from the constraints */
static
SCIP_DECL_RELAXINITSOL(relaxInitSolSdp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RETCODE retcode;
   SCIP_VAR** vars;
   SCIP_Real gaptol;
   SCIP_Real feastol;
   SCIP_Real penaltyparam;
   SCIP_Real maxpenaltyparam;
   int npenaltyincr;
   SCIP_Bool sdpinfo;
   SCIP_Real givenpenaltyparam;
   SCIP_Real projminevprimal;
   SCIP_Real projminevdual;
   SCIP_Real preoptgap;
   int nthreads;
   int slatercheck;
   int nvars;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpinterfacecalls = 0;
   relaxdata->sdpiterations = 0;
   relaxdata->solvedfast = 0;
   relaxdata->solvedmedium = 0;
   relaxdata->solvedstable = 0;
   relaxdata->solvedpenalty = 0;
   relaxdata->stablewslater = 0;
   relaxdata->unstablewslater = 0;
   relaxdata->boundedwslater = 0;
   relaxdata->unsolvedwslater = 0;
   relaxdata->stablenoslater = 0;
   relaxdata->unsolvednoslater = 0;
   relaxdata->boundednoslater = 0;
   relaxdata->unsolvednoslater = 0;
   relaxdata->nslaterholds = 0;
   relaxdata->nnoslater = 0;
   relaxdata->nslatercheckfailed = 0;
   relaxdata->npslaterholds = 0;
   relaxdata->npnoslater = 0;
   relaxdata->npslatercheckfailed = 0;
   relaxdata->ndslaterholds = 0;
   relaxdata->ndnoslater = 0;
   relaxdata->ndslatercheckfailed = 0;
   relaxdata->nslaterinfeasible = 0;
   relaxdata->stableinfeasible = 0;
   relaxdata->unstableinfeasible = 0;
   relaxdata->penaltyinfeasible = 0;
   relaxdata->boundedinfeasible = 0;
   relaxdata->unsolvedinfeasible = 0;
   relaxdata->roundingprobinf = 0;
   relaxdata->primalroundfails = 0;
   relaxdata->dualroundfails = 0;
   relaxdata->roundstartsuccess = 0;
   relaxdata->roundingoptimal = 0;
   relaxdata->roundingcutoff = 0;
   relaxdata->roundingprobtime = 0.0;
   relaxdata->unsolved = 0;
   relaxdata->feasible = FALSE;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* all SCIPvars will be added to this list, and 3/4 seems like a good load factor (java uses this factor) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(relaxdata->varmapper), (int) ceil(1.33 * nvars)) );
   SCIP_CALL( SCIPsdpVarmapperAddVars(scip, relaxdata->varmapper, nvars, vars) );

   if ( nvars > 0 )
   {
      SCIP_CALL( putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, TRUE, FALSE) );
   }

   /* set the parameters of the SDP-Solver */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolvergaptol", &gaptol) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_GAPTOL, gaptol);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: gaptol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolverfeastol", &feastol) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_SDPSOLVERFEASTOL, feastol);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: sdpsolverfeastol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_EPSILON, SCIPepsilon(scip));
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: epsilon setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_FEASTOL, SCIPfeastol(scip));
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: feastol setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* set/compute the starting penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/penaltyparam", &penaltyparam) );
   if ( SCIPisGE(scip, penaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, penaltyparam);
      givenpenaltyparam = penaltyparam;
      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: penaltyparam setting not available -- SCIP parameter has no effect\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }
   }
   else
   {
      SCIP_Real maxcoeff;
      int v;

      /* compute the maximum coefficient in the objective */
      maxcoeff = 0.0;
      for (v = 0; v < nvars; v++)
      {
         if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(vars[v])), maxcoeff) )
            maxcoeff = REALABS(SCIPvarGetObj(vars[v]));
      }

      SCIP_CALL( SCIPsdpiComputePenaltyparam(relaxdata->sdpi, maxcoeff, &givenpenaltyparam) );
   }

   /* set/compute the maximum penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/maxpenaltyparam", &maxpenaltyparam) );
   if ( SCIPisGE(scip, maxpenaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_MAXPENALTYPARAM, maxpenaltyparam);

      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: maxpenaltyparam setting not available -- SCIP parameter has no effect.\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }

      /* check if the starting value is not bigger than the maximum one, otherwise update it */
      if ( SCIPisLT(scip, givenpenaltyparam, maxpenaltyparam) )
      {
         SCIPdebugMessage("Penalty parameter %f overwritten by maxpenaltyparam %f! \n", givenpenaltyparam, maxpenaltyparam);
         SCIP_CALL( SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, maxpenaltyparam) );
      }
   }
   else
   {
      SCIP_Real givenmaxpenaltyparam;

      SCIP_CALL( SCIPsdpiComputeMaxPenaltyparam(relaxdata->sdpi, givenpenaltyparam, &givenmaxpenaltyparam) );
   }

   /* set maximum number of penalty increasing rounds */
   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/npenaltyincr", &npenaltyincr) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_NPENALTYINCR, npenaltyincr);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: npenaltyincr setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }


   /* set/compute lambda star if SDPA is used as the SDP-Solver */
   if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0.0 )
   {
      SCIP_Real lambdastar;

      SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/lambdastar", &lambdastar) );
      if ( SCIPisGE(scip, lambdastar, 0.0) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, lambdastar);
      }
      else
      {
         SCIP_Real guess;
         SCIP_Real maxguess;
         SCIP_CONS** conss;
         int nconss;
         int c;

         /* iterate over all SDP-constraints to compute the biggest guess for lambdastar  */
         conss = SCIPgetConss(scip);
         nconss = SCIPgetNConss(scip);
         maxguess = 0.0;

         for (c = 0; c < nconss; c++)
         {
            /* only check the SDP constraints */
            if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
            {
               SCIP_CALL( SCIPconsSdpGuessInitialPoint(scip, conss[c], &guess) );
               if ( (! SCIPisInfinity(scip, maxguess) ) && SCIPisGT(scip, guess, maxguess) )
                  maxguess = guess;
            }
         }

         SCIP_CALL( SCIPsdpiComputeLambdastar(relaxdata->sdpi, maxguess) );
      }
      retcode = SCIPsdpiGetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, &relaxdata->lambdastar);
   }
   else
      relaxdata->lambdastar = 1.0;

   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: lambdastar setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* set/compute minimum eigenvalue for projecting warmstarting points */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/warmstartprminevpri", &projminevprimal) );
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/warmstartprminevdu", &projminevdual) );

   if ( SCIPisGE(scip, projminevprimal, 0.0) && SCIPisGE(scip, projminevdual, 0.0) )
   {
      relaxdata->warmstartprojminevprimal = projminevprimal;
      relaxdata->warmstartprojminevdual = projminevdual;
   }
   else if ( SCIPisGE(scip, projminevprimal, 0.0) && relaxdata->warmstartprojpdsame )
   {
      relaxdata->warmstartprojminevprimal = projminevprimal;
      relaxdata->warmstartprojminevdual = projminevprimal;
   }
   else if ( SCIPisGE(scip, projminevdual, 0.0) && relaxdata->warmstartprojpdsame )
   {
      relaxdata->warmstartprojminevprimal = projminevdual;
      relaxdata->warmstartprojminevdual = projminevdual;
   }
   else
   {
      SCIP_CONSHDLR* sdpconshdlr;
      SCIP_CONS** sdpblocks;
      int nsdpblocks;
      int b;
      int v;
      SCIP_Real maxsdprhs; /* note that we only take the maximum value of the SDP constraints, since these tend to be the most problematic */
      SCIP_Real maxobj;
      SCIP_Real maxsdpcoef; /* note that we only take the maximum value of the SDP constraints, since these tend to be the most problematic */
      SCIP_Real maxval;
      SCIP_Real sdpcoef;

      /* compute value as WARMSTART_PROJ_FACTOR * max{maxrhs, maxobj, maxsdpcoef} */
      sdpconshdlr = SCIPfindConshdlr(scip, "SDP");
      nsdpblocks = SCIPconshdlrGetNConss(sdpconshdlr);
      sdpblocks = SCIPconshdlrGetConss(sdpconshdlr);

      /* compute maxsdpcoef */
      maxsdpcoef = WARMSTART_PROJ_MINRHSOBJ;
      for (b = 0; b < nsdpblocks; b++)
      {
         sdpcoef = SCIPconsSdpGetMaxSdpCoef(scip, sdpblocks[b]);
         if ( SCIPisGT(scip, sdpcoef, maxsdpcoef) )
            maxsdpcoef = sdpcoef;
      }
      maxsdpcoef *= WARMSTART_PROJ_FACTOR_LHS; /* multiply by additional factor to account for summation of lhs entries */

      /* compute maxsdprhs */
      maxsdprhs = WARMSTART_PROJ_MINRHSOBJ;
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( SCIPisGT(scip, SCIPconsSdpGetMaxConstEntry(scip, sdpblocks[b]), maxsdprhs) )
            maxsdprhs = SCIPconsSdpGetMaxConstEntry(scip, sdpblocks[b]);
      }

      /* compute maxobj */
      maxobj = WARMSTART_PROJ_MINRHSOBJ;
      for (v = 0; v < nvars; v++)
      {
         if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(vars[v])), maxobj) )
            maxobj = REALABS(SCIPvarGetObj(vars[v]));
      }

      if ( relaxdata->warmstartprojpdsame )
      {
         maxval = SCIPisGT(scip, maxsdprhs, maxobj) ? maxsdprhs : maxobj;
         maxval = SCIPisGT(scip, maxsdpcoef, maxval) ? maxsdpcoef : maxval;

         relaxdata->warmstartprojminevprimal = WARMSTART_PROJ_FACTOR * maxval;
         relaxdata->warmstartprojminevdual = WARMSTART_PROJ_FACTOR * maxval;

         SCIPdebugMsg(scip, "Setting warmstartprojminev to %f\n", relaxdata->warmstartprojminevdual);
      }
      else
      {
         if ( ! SCIPisGE(scip, projminevprimal, 0.0) )
         {
            relaxdata->warmstartprojminevprimal = WARMSTART_PROJ_FACTOR_PRIMAL * (SCIPisGT(scip, maxobj, maxsdpcoef) ? maxobj : maxsdpcoef);

            SCIPdebugMsg(scip, "Setting warmstartprojminevprimal to %f\n", relaxdata->warmstartprojminevprimal);
         }

         if ( ! SCIPisGE(scip, projminevdual, 0.0) )
         {
            relaxdata->warmstartprojminevdual = WARMSTART_PROJ_FACTOR_PRIMAL * (SCIPisGT(scip, maxsdprhs, maxsdpcoef) ? maxsdprhs : maxsdpcoef);

            SCIPdebugMsg(scip, "Setting warmstartprojminevdual to %f\n", relaxdata->warmstartprojminevdual);
         }
      }
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/sdpinfo", &sdpinfo) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SDPINFO, (int) sdpinfo);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: sdpinfo setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/sdpsolverthreads", &nthreads) );
   /* only try to set nthreads if the value differs from the default to prevent unnecessary warning messages for unknown parameter */
   if ( nthreads != DEFAULT_SDPSOLVERTHREADS )
   {
      retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_NTHREADS, nthreads);
      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: nthreads setting not available -- SCIP parameter has no effect.\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }
   }

   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/slatercheck", &slatercheck) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SLATERCHECK, slatercheck);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: slatercheck setting not available -- SCIP parameter has no effect.\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* initialize objective limit in case it was set in an earlier optimize call */
   SCIP_CALL( SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPsdpiInfinity(relaxdata->sdpi)) );

   /* set warmstartpreoptimal gap if DSDP is used as the SDP-Solver and preoptimal solutions should be saved */
   if ( relaxdata->warmstartpreoptsol && (strcmp(SCIPsdpiGetSolverName(), "DSDP") == 0.0 || strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0.0) )
   {
      SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/warmstartpreoptgap", &preoptgap) );
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_WARMSTARTPOGAP, preoptgap);
      if ( retcode == SCIP_PARAMETERUNKNOWN )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "SDP Solver <%s>: warmstartpreoptgap setting not available -- SCIP parameter has no effect.\n",
            SCIPsdpiGetSolverName());
      }
      else
      {
         SCIP_CALL( retcode );
      }
   }

   return SCIP_OKAY;
}

/** copy method for SDP-relaxation handler (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopySdp)
{
   assert( scip != NULL );
   assert( relax != NULL );
   assert( strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0 );

   SCIP_CALL( SCIPincludeRelaxSdp(scip) );

   return SCIP_OKAY;
}

/** reset the relaxator's data */
static
SCIP_DECL_RELAXEXIT(relaxExitSdp)
{
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   SCIPdebugMessage("Exiting Relaxation Handler.\n");

   if ( relaxdata->displaystat && SCIPgetSubscipDepth(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "\nSDP iterations:\t\t\t\t%6d\n", relaxdata->sdpiterations);
      if ( relaxdata->sdpcalls )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Average SDP-iterations:\t\t\t%6.2f \n", (SCIP_Real) relaxdata->sdpiterations / (SCIP_Real) relaxdata->sdpcalls );
      }
      if ( relaxdata->sdpinterfacecalls )
      {
         if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedfast / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'medium settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedmedium / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedstable / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         else
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'default formulation' solved:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedfast / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage penalty formulation used:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->solvedpenalty / (SCIP_Real) relaxdata->sdpinterfacecalls);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage unsolved even with penalty:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->unsolved / (SCIP_Real) relaxdata->sdpinterfacecalls);
      }
      if ( relaxdata->slatercheck )
      {
         if ( relaxdata->sdpinterfacecalls )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater condition held:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npslaterholds / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater condition did not hold:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npnoslater / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage primal Slater check failed:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->npslatercheckfailed / (SCIP_Real) relaxdata->sdpinterfacecalls);

            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater condition held:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndslaterholds / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater condition did not hold:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndnoslater / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater check failed:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->ndslatercheckfailed / (SCIP_Real) relaxdata->sdpinterfacecalls);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage dual Slater check detected infeasibility:\t%6.2f \n", 100.0 * (SCIP_Real) relaxdata->nslaterinfeasible / (SCIP_Real) relaxdata->sdpinterfacecalls);
         }
         if ( relaxdata->nslaterholds )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stablewslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstablewslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltywslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundedwslater / (SCIP_Real) relaxdata->nslaterholds);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with primal and dual slater holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvedwslater / (SCIP_Real) relaxdata->nslaterholds);
         }
         if ( relaxdata->nnoslater )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stablenoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstablenoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltynoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundednoslater / (SCIP_Real) relaxdata->nnoslater);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with either primal or dual slater not holding:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvednoslater / (SCIP_Real) relaxdata->nnoslater);
         }
         if ( relaxdata->nslaterinfeasible )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'fastest settings' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->stableinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'stable settings' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unstableinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'penalty' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->penaltyinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'computed infeasible lower bound' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->boundedinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Percentage 'unsolved' with slater check showing infeasibility:\t%6.2f \n",
                  100.0 * (SCIP_Real) relaxdata->unsolvedinfeasible / (SCIP_Real) relaxdata->nslaterinfeasible);
         }
#ifdef SLATERSOLVED_ABSOLUTE
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with primal and dual slater holding:\t%d \n", relaxdata->nslaterholds);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'fastest settings' and primal and dual slater holding:\t%d \n", relaxdata->stablewslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'stable settings' and primal and dual slater holding:\t%d \n", relaxdata->unstablewslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'penalty' and primal and dual slater holding:\t%d \n", relaxdata->penaltywslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'computed infeasible lower bound' and primal and dual slater holding:\t%d \n", relaxdata->boundedwslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'unsolved' and primal and dual slater holding:\t%d \n", relaxdata->unsolvedwslater);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with either primal or dual slater not holding:\t%d \n", relaxdata->nnoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'fastest settings' and either primal or dual slater not holding:\t%d \n", relaxdata->stablenoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'stable settings' and either primal or dual slater not holding:\t%d \n", relaxdata->unstablenoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'penalty' and either primal or dual slater not holding:\t%d \n", relaxdata->penaltynoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'computed infeasible lower bound' and either primal or dual slater not holding:\t%d \n", relaxdata->boundednoslater);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes with 'unsolved' and either primal or dual slater not holding:\t%d \n", relaxdata->unsolvednoslater);

         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes:\t%d \n", relaxdata->nslaterinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'fastest settings':\t%d \n", relaxdata->stableinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'stable settings':\t%d \n", relaxdata->unstableinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'penalty':\t%d \n", relaxdata->penaltyinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'computed infeasible lower bound':\t%d \n", relaxdata->boundedinfeasible);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of infeasible nodes with 'unsolved':\t%d \n", relaxdata->unsolvedinfeasible);
#endif
      }
      if ( relaxdata->warmstart && relaxdata->warmstartproject == 4 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes detected infeasible through primal rounding problem:\t%d \n", relaxdata->roundingprobinf);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes that were successfully warmstarted using the rounding problems:\t%d \n", relaxdata->roundstartsuccess);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes where the primal rounding problem failed:\t%d \n", relaxdata->primalroundfails);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes where the dual rounding problem failed:\t%d \n", relaxdata->dualroundfails);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes where the optimal solution was determined by the rounding problem:\t%d \n", relaxdata->roundingoptimal);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Number of nodes cut off through bounding by the rounding problem:\t%d \n", relaxdata->roundingcutoff);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Time spent in rounding problems for warmstarting / detecting infeasibility:\t%f s \n", relaxdata->roundingprobtime);
      }
   }

   if ( relaxdata->varmapper != NULL )
   {
      SCIP_CALL( SCIPsdpVarmapperFree(scip, &(relaxdata->varmapper)) );
   }

   /* free warmstart data (the nblocks > 0 check is only needed in case the parameter was changed after initsol) */
   if ( relaxdata->warmstart && SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) && relaxdata->nblocks > 0 )
   {
      int b;

      for (b = 0; b < relaxdata->nblocks; b++)
      {
         if ( relaxdata->warmstartprimaltype != 2 && SCIPsdpiDoesWarmstartNeedPrimal() && relaxdata->ipXnblocknonz[b] > 0 )
         {
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipXval[b]), relaxdata->ipXnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipXcol[b]), relaxdata->ipXnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipXrow[b]), relaxdata->ipXnblocknonz[b]);
         }
         if ( relaxdata->ipZnblocknonz[b] > 0 )
         {
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipZval[b]), relaxdata->ipZnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipZcol[b]), relaxdata->ipZnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(relaxdata->ipZrow[b]), relaxdata->ipZnblocknonz[b]);
         }
      }
      if ( relaxdata->warmstartprimaltype != 2 && SCIPsdpiDoesWarmstartNeedPrimal() )
      {
         SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipXval, relaxdata->nblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipXcol, relaxdata->nblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipXrow, relaxdata->nblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipXnblocknonz, relaxdata->nblocks);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipZval, relaxdata->nblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipZcol, relaxdata->nblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipZrow, relaxdata->nblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &relaxdata->ipZnblocknonz, relaxdata->nblocks);
      SCIP_CALL( SCIPfreeSol(scip, &relaxdata->ipy) );
   }

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->feasible = FALSE;
   relaxdata->sdpiterations = 0;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpinterfacecalls = 0;
   relaxdata->lastsdpnode = 0;
   SCIP_CALL( SCIPsdpiClear(relaxdata->sdpi) );

   return SCIP_OKAY;
}

/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSdp)
{/*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if ( relaxdata->sdpi != NULL )
   {
      SCIP_CALL( SCIPsdpiFree(&(relaxdata->sdpi)) );
   }
   if ( relaxdata->lpi != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&(relaxdata->lpi)) );
   }

   SCIPfreeMemory(scip, &relaxdata);

   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}

/** creates the SDP-relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax;
   SCIP_SDPI* sdpi;
   SCIP_LPI* lpi;

   assert( scip != NULL );

   /* create SDP-relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   SCIP_CALL( SCIPsdpiCreate(&sdpi, SCIPgetMessagehdlr(scip), SCIPblkmem(scip), SCIPbuffer(scip)) );
   SCIP_CALL( SCIPlpiCreate(&lpi, SCIPgetMessagehdlr(scip), "SDProundingProb", SCIP_OBJSEN_MINIMIZE) );

   relaxdata->sdpi = sdpi;
   relaxdata->lpi = lpi;
   relaxdata->lastsdpnode = -1;
   relaxdata->nblocks = 0;
   relaxdata->ipXexists = FALSE;
   relaxdata->ipZexists = FALSE;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, relaxExecSdp, relaxdata) );
   assert( relax != NULL );

   /* include additional callbacks */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitSolSdp) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitSdp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeSdp) );
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopySdp) );

   /* add parameters for SDP-solver */
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolvergaptol",
         "the stopping criterion for the duality gap the sdpsolver should use",
         &(relaxdata->sdpsolvergaptol), TRUE, SCIPsdpiGetDefaultSdpiSolverGaptol(), 1e-20, 0.001, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolverfeastol",
         "the feasibility tolerance for the SDP solver",
         &(relaxdata->sdpsolverfeastol), TRUE, SCIPsdpiGetDefaultSdpiSolverFeastol(), 1e-17, 0.001, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/penaltyparam",
         "the starting value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->penaltyparam),
         TRUE, DEFAULT_PENALTYPARAM, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/maxpenaltyparam",
         "the maximum value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge; set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->maxpenaltyparam),
         TRUE, DEFAULT_MAXPENALTYPARAM, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/warmstartipfactor",
         "factor for interior point in convexcombination of IP and parent solution, if warmstarts are enabled", &(relaxdata->warmstartipfactor),
         TRUE, DEFAULT_WARMSTARTIPFACTOR, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/warmstartprimaltype",
         "how to warmstart the primal problem? 1: scaled identity/analytic center, 2: elementwise reciprocal, 3: saved primal sol", &(relaxdata->warmstartprimaltype), TRUE,
         DEFAULT_WARMSTARTPRIMALTYPE, 1, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/warmstartiptype",
         "which interior point to use for convex combination for warmstarts? 1: scaled identity, 2: analytic center", &(relaxdata->warmstartiptype), TRUE,
         DEFAULT_WARMSTARTIPTYPE, 1, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/warmstartproject",
         "how to update dual matrix for new bounds? 1: use old bounds, 2: use new bounds, 3: use new bounds and project on psd cone, 4: use new bounds and solve rounding problem", &(relaxdata->warmstartproject), TRUE,
         DEFAULT_WARMSTARTPROJECT, 1, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/warmstartprminevpri",
         "minimum eigenvalue to allow when projecting primal matrices onto the positive (semi-)definite cone for warmstarting; -1 to compute automatically", &(relaxdata->warmstartpmevprimalpar),
         TRUE, DEFAULT_WARMSTARTPROJMINEV, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/warmstartprminevdu",
         "minimum eigenvalue to allow when projecting dual matrices onto the positive (semi-)definite cone for warmstarting; -1 to compute automatically", &(relaxdata->warmstartpmevdualpar),
         TRUE, DEFAULT_WARMSTARTPROJMINEV, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/warmstartprojpdsame",
         "Should one shared minimum eigenvalue respectively maximum entry be computed for primal and dual problem instead of different ones for primal and dual and each block for projection or convex combination ?",
         &(relaxdata->warmstartprojpdsame), TRUE, DEFAULT_WARMSTARTPROJPDSAME, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/warmstartpreoptsol",
         "Should a preoptimal solution (with higher epsilon) instead of the optimal solution be used for warmstarts",
         &(relaxdata->warmstartpreoptsol), TRUE, DEFAULT_WARMSTART_PREOPTIMAL_SOL, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/warmstartpreoptgap",
         "If warmstartpreoptsol is TRUE, this is the gap where the preoptimal solution will be saved", &(relaxdata->warmstartpreoptgap),
         TRUE, DEFAULT_WARMSTARTPREOPTGAP, 0.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/warmstartroundonlyinf",
         "Only use solution of roundingproblem to detect infeasibility (only has an effect for warmstartproject = 4)",
         &(relaxdata->warmstartroundonlyinf), TRUE, DEFAULT_WARMSTARTROUNDONLYINF, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/npenaltyincr",
         "maximum number of times the penalty parameter will be increased if the penalty formulation failed", &(relaxdata->npenaltyincr), TRUE,
         SCIPsdpiGetDefaultSdpiSolverNpenaltyIncreases(), 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/lambdastar",
         "the parameter lambda star used by SDPA to set the initial point;"
         "set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->lambdastar),
         TRUE, DEFAULT_LAMBDASTAR, -1.0, 1e+20, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/slatercheck",
         "Should the Slater condition for the primal and dual problem be checked ahead of solving each SDP? 0: no, 1: yes but only for statistics, 2: yes and print warning for "
         "every problem not satisfying primal and dual Slater condition", &(relaxdata->slatercheck), TRUE, DEFAULT_SLATERCHECK, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/sdpinfo",
         "Should the SDP solver output information to the screen?",
         &(relaxdata->sdpinfo), TRUE, DEFAULT_SDPINFO, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/warmstart",
         "Should the SDP solver try to use warmstarts?",
         &(relaxdata->warmstart), TRUE, DEFAULT_WARMSTART, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/objlimit",
         "Should an objective limit be given to the SDP-Solver?",
         &(relaxdata->objlimit), TRUE, DEFAULT_OBJLIMIT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/resolve",
         "Should the relaxation be resolved after bound-tightenings were found during propagation (outside of probing)?",
         &(relaxdata->resolve), TRUE, DEFAULT_RESOLVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/tightenvb",
         "Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ?",
         &(relaxdata->tightenvb), TRUE, DEFAULT_TIGHTENVB, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/displaystatistics",
         "Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ?",
         &(relaxdata->displaystat), TRUE, DEFAULT_DISPLAYSTAT, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetfreq",
         "frequency for resetting parameters in SDP solver and trying again with fastest settings (-1: never, 0: only at depth settingsresetofs);"
         "currently only supported for SDPA",
         &(relaxdata->settingsresetfreq), TRUE, DEFAULT_SETTINGSRESETFREQ, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetofs",
         "frequency offset for resetting parameters in SDP solver and trying again with fastest settings; currently only supported for SDPA",
         &(relaxdata->settingsresetofs), TRUE, DEFAULT_SETTINGSRESETOFS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/sdpsolverthreads",
         "number of threads the SDP solver should use (-1 = number of cores); currently only supported for MOSEK",
         &(relaxdata->sdpsolverthreads), TRUE, DEFAULT_SDPSOLVERTHREADS, -1, INT_MAX, NULL, NULL) );


   /* add description of SDP-solver */
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPsdpiGetSolverName(), SCIPsdpiGetSolverDesc()) );

   return SCIP_OKAY;
}


/* external functions */

/** computes analytic centers of primal and dual feasible set and saves them in relaxdata
 * @note This function should be called at the end of the root node (or at least after the solving stage starts and before the first non-root node).
 */
SCIP_RETCODE SCIPrelaxSdpComputeAnalyticCenters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< SDP-relaxator to compute analytic centers for */
   )
{
   SCIP_CONSHDLR* sdpconshdlr;
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_SOLVING );

   SCIPdebugMsg(scip, "computing analytic centers for warmstarting\n");

   relaxdata = SCIPrelaxGetData(relax);

   /* this function should only be executed once */
   if ( relaxdata->ipXexists || relaxdata->ipZexists )
   {
      SCIPdebugMsg(scip, "aborting SCIPrelaxSdpComputeAnalyticCenters since analytic centers have already been computed\n");
      return SCIP_OKAY;
   }

   sdpconshdlr = SCIPfindConshdlr(scip, "SDP");

   /* if we want to warmstart using the analytic center, compute it now */
   if ( relaxdata->warmstart && (relaxdata->warmstartiptype == 2) && SCIPisGT(scip, relaxdata->warmstartipfactor, 0.0) && (SCIPconshdlrGetNConss(sdpconshdlr) + SCIPgetNLPRows(scip) > 0 ) )
   {
      int b;

      relaxdata->nblocks = SCIPgetNLPRows(scip) + SCIPgetNVars(scip) > 0 ? SCIPconshdlrGetNConss(sdpconshdlr) + 1 : SCIPconshdlrGetNConss(sdpconshdlr);

      if ( SCIPgetNVars(scip) > 0 )
      {
         SCIP_ROW** rows;
         SCIP_COL** rowcols;
         SCIP_CONS** sdpblocks;
         SCIP_Real* solforscip;
         SCIP_Real* rowvals;
         SCIP_Real timelimit;
         SCIP_Real rowval;
         int slength;
         int arraylength;
         int nrows;
         int rownnonz;
         int i;
         int r;
         int v;
         SCIP_SDPSOLVERSETTING usedsetting;
         SCIP_SDPSLATERSETTING slatersetting;
         SCIP_SDPSLATER primalslater;
         SCIP_SDPSLATER dualslater;
         int naddediters;
         int naddedsdpcalls;

         /* first solve SDP with primal objective (dual constant part) set to zero to compute analytic center of primal feasible set */
         if ( relaxdata->warmstartprimaltype != 2 && SCIPsdpiDoesWarmstartNeedPrimal() )
         {
            SCIP_CALL( putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, FALSE, TRUE) );
            SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, FALSE, TRUE) );

            /* set time limit */
            SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
            if ( ! SCIPisInfinity(scip, timelimit) )
            {
               timelimit -= SCIPgetSolvingTime(scip);
               if ( timelimit <= 0.0 )
                  return SCIP_OKAY;
            }

            /* TODO: might want to add an additional parameter to solve to disable penalty, since we cannot use that here anyways */
            SCIP_CALL(SCIPsdpiSolve(relaxdata->sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, timelimit));

            /* update calls, iterations and stability numbers (only if the SDP-solver was actually called) */
            relaxdata->sdpinterfacecalls++;
            naddedsdpcalls = 0;
            SCIP_CALL( SCIPsdpiGetSdpCalls(relaxdata->sdpi, &naddedsdpcalls) );
            usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
            if ( naddedsdpcalls )
            {
               relaxdata->sdpcalls += naddedsdpcalls;
               naddediters = 0;
               SCIP_CALL( SCIPsdpiGetIterations(relaxdata->sdpi, &naddediters) );
               relaxdata->sdpiterations += naddediters;

               SCIP_CALL( SCIPsdpiSettingsUsed(relaxdata->sdpi, &usedsetting) );

               switch( usedsetting )/*lint --e{788}*/
               {
               case SCIP_SDPSOLVERSETTING_PENALTY:
                  relaxdata->solvedpenalty++;
                  break;
               case SCIP_SDPSOLVERSETTING_FAST:
                  relaxdata->solvedfast++;
                  break;
               case SCIP_SDPSOLVERSETTING_MEDIUM:
                  relaxdata->solvedmedium++;
                  break;
               case SCIP_SDPSOLVERSETTING_STABLE:
                  relaxdata->solvedstable++;
                  break;
               case SCIP_SDPSOLVERSETTING_UNSOLVED:
                  relaxdata->unsolved++;
                  break;
               default:
                  break;
               }
               primalslater = SCIP_SDPSLATER_NOINFO;
               dualslater = SCIP_SDPSLATER_NOINFO;
               SCIP_CALL( SCIPsdpiSlater(relaxdata->sdpi, &primalslater, &dualslater) );
               switch( primalslater )/*lint --e{788}*/
               {
                  case SCIP_SDPSLATER_NOINFO:
                     relaxdata->npslatercheckfailed++;
                     switch( dualslater )/*lint --e{788}*/
                     {
                        case SCIP_SDPSLATER_NOINFO:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nslatercheckfailed++;
                           break;
                        case SCIP_SDPSLATER_NOT:
                           relaxdata->ndnoslater++;
                           relaxdata->nnoslater++;
                           break;
                        case SCIP_SDPSLATER_HOLDS:
                           relaxdata->ndslaterholds++;
                           relaxdata->nslatercheckfailed++;
                           break;
                        case SCIP_SDPSLATER_INF:
                           relaxdata->nslaterinfeasible++;
                           break;
                        default:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nslatercheckfailed++;
                           break;
                     }
                  break;
                  case SCIP_SDPSLATER_NOT:
                     relaxdata->npnoslater++;
                     switch( dualslater )/*lint --e{788}*/
                     {
                        case SCIP_SDPSLATER_NOINFO:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nnoslater++;
                           break;
                        case SCIP_SDPSLATER_NOT:
                           relaxdata->ndnoslater++;
                           relaxdata->nnoslater++;
                           break;
                        case SCIP_SDPSLATER_HOLDS:
                           relaxdata->ndslaterholds++;
                           relaxdata->nnoslater++;
                           break;
                        case SCIP_SDPSLATER_INF:
                           relaxdata->nslaterinfeasible++;
                           break;
                        default:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nnoslater++;
                           break;
                     }
                     break;
                  case SCIP_SDPSLATER_HOLDS:
                     relaxdata->npslaterholds++;
                     switch( dualslater )/*lint --e{788}*/
                     {
                        case SCIP_SDPSLATER_NOINFO:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nslatercheckfailed++;
                           break;
                        case SCIP_SDPSLATER_NOT:
                           relaxdata->ndnoslater++;
                           relaxdata->nnoslater++;
                           break;
                        case SCIP_SDPSLATER_HOLDS:
                           relaxdata->ndslaterholds++;
                           relaxdata->nslaterholds++;
                           break;
                        case SCIP_SDPSLATER_INF:
                           relaxdata->nslaterinfeasible++;
                           break;
                        default:
                           relaxdata->ndslatercheckfailed++;
                           relaxdata->nslatercheckfailed++;
                           break;
                     }
                     break;
                     default:
                        relaxdata->npslatercheckfailed++;
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nslatercheckfailed++;
                        break;
               }
               slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
               SCIP_CALL( SCIPsdpiSlaterSettings(relaxdata->sdpi, &slatersetting) );
               switch( slatersetting )/*lint --e{788}*/
               {
                  case SCIP_SDPSLATERSETTING_STABLEWSLATER:
                     relaxdata->stablewslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSTABLEWSLATER:
                     relaxdata->unstablewslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_PENALTYWSLATER:
                     relaxdata->penaltywslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_BOUNDEDWSLATER:
                     relaxdata->boundedwslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER:
                     relaxdata->unsolvedwslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_STABLENOSLATER:
                     relaxdata->stablenoslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSTABLENOSLATER:
                     relaxdata->unstablenoslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_PENALTYNOSLATER:
                     relaxdata->penaltynoslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER:
                     relaxdata->boundednoslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER:
                     relaxdata->unsolvednoslater++;
                     break;
                  case SCIP_SDPSLATERSETTING_STABLEINFEASIBLE:
                     relaxdata->stableinfeasible++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE:
                     relaxdata->unstableinfeasible++;
                     break;
                  case SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE:
                     relaxdata->penaltyinfeasible++;
                     break;
                  case SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE:
                     relaxdata->boundedinfeasible++;
                     break;
                  case SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE:
                     relaxdata->unsolvedinfeasible++;
                     break;
                  default:
                     break;
               }
            }

            if ( SCIPsdpiWasSolved(relaxdata->sdpi) && SCIPsdpiSolvedOrig(relaxdata->sdpi) && SCIPsdpiIsPrimalFeasible(relaxdata->sdpi) )
            {
               int npenaltybounds = 0;

               relaxdata->ipXexists = TRUE;

               /* allocate memory (for the different blocks the neccessary anount first needs to be computed) */
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXnblocknonz, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXrow, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXcol, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXval, relaxdata->nblocks) );

               SCIP_CALL(SCIPsdpiGetPrimalNonzeros(relaxdata->sdpi, relaxdata->nblocks, relaxdata->ipXnblocknonz));
               for (b = 0; b < relaxdata->nblocks; b++)
               {
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXrow[b], relaxdata->ipXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXcol[b], relaxdata->ipXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXval[b], relaxdata->ipXnblocknonz[b]) );
               }

               /* get primal solution */
               SCIP_CALL( SCIPsdpiGetPrimalMatrix(relaxdata->sdpi, relaxdata->nblocks, relaxdata->ipXnblocknonz,
                     relaxdata->ipXrow, relaxdata->ipXcol, relaxdata->ipXval) );

               /* count the number of primal entries corresponding to bounds of the penalty variable and remove them */
               for (i = 0; i < relaxdata->ipXnblocknonz[relaxdata->nblocks - 1]; i++)
               {
                  if ( relaxdata->ipXrow[relaxdata->nblocks - 1][i] == SCIPsdpVarmapperGetNVars(relaxdata->varmapper) )
                     npenaltybounds++;
               }

               if ( npenaltybounds > 0 )
               {
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipXrow[relaxdata->nblocks - 1],
                        relaxdata->ipXnblocknonz[relaxdata->nblocks - 1], relaxdata->ipXnblocknonz[relaxdata->nblocks - 1] - npenaltybounds) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipXcol[relaxdata->nblocks - 1],
                        relaxdata->ipXnblocknonz[relaxdata->nblocks - 1], relaxdata->ipXnblocknonz[relaxdata->nblocks - 1] - npenaltybounds) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipXval[relaxdata->nblocks - 1],
                        relaxdata->ipXnblocknonz[relaxdata->nblocks - 1], relaxdata->ipXnblocknonz[relaxdata->nblocks - 1] - npenaltybounds) );
                  relaxdata->ipXnblocknonz[relaxdata->nblocks - 1] = relaxdata->ipXnblocknonz[relaxdata->nblocks - 1] - npenaltybounds;
               }

               for (b = 0; b < relaxdata->nblocks; b++)
               {
                  /* TODO check if they are already sorted (sorting is needed since they will later be merged into warmstart arrays) */
                  SCIPsdpVarfixerSortRowCol(relaxdata->ipXrow[b], relaxdata->ipXcol[b], relaxdata->ipXval[b], relaxdata->ipXnblocknonz[b]);
               }

#ifdef SCIP_PRINT_WARMSTART
               SCIPdebugMessage("Computed primal analytic center:\n");
               for (b = 0; b < relaxdata->nblocks; b++)
               {
                  SCIPdebugMessage("primal matrix, block %d:\n", b);
                  for (i = 0; i < relaxdata->ipXnblocknonz[b]; i++)
                  {
                     SCIPdebugMessage("X_%d[%d,%d]: %f\n", b, relaxdata->ipXrow[b][i], relaxdata->ipXcol[b][i], relaxdata->ipXval[b][i]);
                  }
               }
#endif
            }
            else
            {
               relaxdata->ipXexists = TRUE;

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXnblocknonz, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXrow, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXcol, relaxdata->nblocks) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXval, relaxdata->nblocks) );

               sdpblocks = SCIPconshdlrGetConss(sdpconshdlr);

               for (b = 0; b < relaxdata->nblocks; b++)
               {
                  if ( b < relaxdata->nblocks - 1 )
                  {
                     /* SDP block */
                     relaxdata->ipXnblocknonz[b] = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);
                  }
                  else
                  {
                     /* LP block */
                     relaxdata->ipXnblocknonz[b] = SCIPgetNLPRows(scip) + 2 * SCIPgetNVars(scip);
                  }
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXrow[b], relaxdata->ipXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXcol[b], relaxdata->ipXnblocknonz[b]) );
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipXval[b], relaxdata->ipXnblocknonz[b]) );

                  for (i = 0; i < relaxdata->ipXnblocknonz[b]; i++)
                  {
                     relaxdata->ipXrow[b][i] = i;
                     relaxdata->ipXcol[b][i] = i;
                     relaxdata->ipXval[b][i] = relaxdata->lambdastar;
                  }
               }

               SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Failed to compute analytic center of primal feasible set, using scaled identity instead.\n");
            }
         }

         /* set dual objective coefficients to zero to compute analytic center of dual feasible set */
         SCIP_CALL( putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, TRUE, FALSE) );
         SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper, TRUE, FALSE) );

         /* set time limit */
         SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
         if ( ! SCIPisInfinity(scip, timelimit) )
         {
            timelimit -= SCIPgetSolvingTime(scip);
            if ( timelimit <= 0.0 )
               return SCIP_OKAY;
         }

         /* TODO: might want to add an additional parameter to solve to disable penalty, since we cannot use that here anyways */
         SCIP_CALL(SCIPsdpiSolve(relaxdata->sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, timelimit));

         /* update calls, iterations and stability numbers (only if the SDP-solver was actually called) */
         naddedsdpcalls = 0;
         SCIP_CALL( SCIPsdpiGetSdpCalls(relaxdata->sdpi, &naddedsdpcalls) );
         usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
         if ( naddedsdpcalls )
         {
            relaxdata->sdpinterfacecalls++;
            relaxdata->sdpcalls += naddedsdpcalls;
            naddediters = 0;
            SCIP_CALL( SCIPsdpiGetIterations(relaxdata->sdpi, &naddediters) );
            relaxdata->sdpiterations += naddediters;

            SCIP_CALL( SCIPsdpiSettingsUsed(relaxdata->sdpi, &usedsetting) );

            switch( usedsetting )/*lint --e{788}*/
            {
            case SCIP_SDPSOLVERSETTING_PENALTY:
               relaxdata->solvedpenalty++;
               break;
            case SCIP_SDPSOLVERSETTING_FAST:
               relaxdata->solvedfast++;
               break;
            case SCIP_SDPSOLVERSETTING_MEDIUM:
               relaxdata->solvedmedium++;
               break;
            case SCIP_SDPSOLVERSETTING_STABLE:
               relaxdata->solvedstable++;
               break;
            case SCIP_SDPSOLVERSETTING_UNSOLVED:
               relaxdata->unsolved++;
               break;
            default:
               break;
            }
            primalslater = SCIP_SDPSLATER_NOINFO;
            dualslater = SCIP_SDPSLATER_NOINFO;
            SCIP_CALL( SCIPsdpiSlater(relaxdata->sdpi, &primalslater, &dualslater) );
            switch( primalslater )/*lint --e{788}*/
            {
               case SCIP_SDPSLATER_NOINFO:
                  relaxdata->npslatercheckfailed++;
                  switch( dualslater )/*lint --e{788}*/
                  {
                     case SCIP_SDPSLATER_NOINFO:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nslatercheckfailed++;
                        break;
                     case SCIP_SDPSLATER_NOT:
                        relaxdata->ndnoslater++;
                        relaxdata->nnoslater++;
                        break;
                     case SCIP_SDPSLATER_HOLDS:
                        relaxdata->ndslaterholds++;
                        relaxdata->nslatercheckfailed++;
                        break;
                     case SCIP_SDPSLATER_INF:
                        relaxdata->nslaterinfeasible++;
                        break;
                     default:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nslatercheckfailed++;
                        break;
                  }
               break;
               case SCIP_SDPSLATER_NOT:
                  relaxdata->npnoslater++;
                  switch( dualslater )/*lint --e{788}*/
                  {
                     case SCIP_SDPSLATER_NOINFO:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nnoslater++;
                        break;
                     case SCIP_SDPSLATER_NOT:
                        relaxdata->ndnoslater++;
                        relaxdata->nnoslater++;
                        break;
                     case SCIP_SDPSLATER_HOLDS:
                        relaxdata->ndslaterholds++;
                        relaxdata->nnoslater++;
                        break;
                     case SCIP_SDPSLATER_INF:
                        relaxdata->nslaterinfeasible++;
                        break;
                     default:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nnoslater++;
                        break;
                  }
                  break;
               case SCIP_SDPSLATER_HOLDS:
                  relaxdata->npslaterholds++;
                  switch( dualslater )/*lint --e{788}*/
                  {
                     case SCIP_SDPSLATER_NOINFO:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nslatercheckfailed++;
                        break;
                     case SCIP_SDPSLATER_NOT:
                        relaxdata->ndnoslater++;
                        relaxdata->nnoslater++;
                        break;
                     case SCIP_SDPSLATER_HOLDS:
                        relaxdata->ndslaterholds++;
                        relaxdata->nslaterholds++;
                        break;
                     case SCIP_SDPSLATER_INF:
                        relaxdata->nslaterinfeasible++;
                        break;
                     default:
                        relaxdata->ndslatercheckfailed++;
                        relaxdata->nslatercheckfailed++;
                        break;
                  }
                  break;
                  default:
                     relaxdata->npslatercheckfailed++;
                     relaxdata->ndslatercheckfailed++;
                     relaxdata->nslatercheckfailed++;
                     break;
            }
            slatersetting = SCIP_SDPSLATERSETTING_NOINFO;
            SCIP_CALL( SCIPsdpiSlaterSettings(relaxdata->sdpi, &slatersetting) );
            switch( slatersetting )/*lint --e{788}*/
            {
               case SCIP_SDPSLATERSETTING_STABLEWSLATER:
                  relaxdata->stablewslater++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSTABLEWSLATER:
                  relaxdata->unstablewslater++;
                  break;
               case SCIP_SDPSLATERSETTING_PENALTYWSLATER:
                  relaxdata->penaltywslater++;
                  break;
               case SCIP_SDPSLATERSETTING_BOUNDEDWSLATER:
                  relaxdata->boundedwslater++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSOLVEDWSLATER:
                  relaxdata->unsolvedwslater++;
                  break;
               case SCIP_SDPSLATERSETTING_STABLENOSLATER:
                  relaxdata->stablenoslater++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSTABLENOSLATER:
                  relaxdata->unstablenoslater++;
                  break;
               case SCIP_SDPSLATERSETTING_PENALTYNOSLATER:
                  relaxdata->penaltynoslater++;
                  break;
               case SCIP_SDPSLATERSETTING_BOUNDEDNOSLATER:
                  relaxdata->boundednoslater++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSOLVEDNOSLATER:
                  relaxdata->unsolvednoslater++;
                  break;
               case SCIP_SDPSLATERSETTING_STABLEINFEASIBLE:
                  relaxdata->stableinfeasible++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSTABLEINFEASIBLE:
                  relaxdata->unstableinfeasible++;
                  break;
               case SCIP_SDPSLATERSETTING_PENALTYINFEASIBLE:
                  relaxdata->penaltyinfeasible++;
                  break;
               case SCIP_SDPSLATERSETTING_BOUNDEDINFEASIBLE:
                  relaxdata->boundedinfeasible++;
                  break;
               case SCIP_SDPSLATERSETTING_UNSOLVEDINFEASIBLE:
                  relaxdata->unsolvedinfeasible++;
                  break;
               default:
                  break;
            }
         }

         if ( SCIPsdpiWasSolved(relaxdata->sdpi) && SCIPsdpiSolvedOrig(relaxdata->sdpi) && SCIPsdpiIsDualFeasible(relaxdata->sdpi) )
         {
            int nvars;
            SCIP_VAR** vars;

            relaxdata->ipZexists = TRUE;

            nvars = SCIPgetNVars(scip);
            vars = SCIPgetVars(scip);

            /* allocate memory */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZnblocknonz, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZrow, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZcol, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZval, relaxdata->nblocks) );

            /* get solution w.r.t. SCIP variables */
            SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
            slength = nvars;

            SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, NULL, solforscip, &slength) ); /* get the solution from the SDP solver */

            assert( slength == nvars ); /* If this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                         * but this can't happen, because the array has enough space for all sdp variables. */

            /* create SCIP solution */
            SCIP_CALL( SCIPcreateSol(scip, &relaxdata->ipy, NULL) );
            SCIP_CALL( SCIPsetSolVals(scip, relaxdata->ipy, nvars, vars, solforscip) );
#ifdef SCIP_PRINT_WARMSTART
            SCIPdebugMessage("Computed dual analytic center:\n");
            for (i = 0; i < nvars; i++)
            {
               SCIPdebugMessage("y[%d] = %f\n", i, solforscip[i]);
            }
#endif

            SCIPfreeBufferArray(scip, &solforscip);

            /* compute SDP blocks of dual analytic center */
            sdpblocks = SCIPconshdlrGetConss(sdpconshdlr);
            for (b = 0; b < relaxdata->nblocks - 1; b++)
            {
               relaxdata->ipZnblocknonz[b] = SCIPconsSdpComputeUbSparseSdpMatrixLength(scip, sdpblocks[b]);
               arraylength = relaxdata->ipZnblocknonz[b];

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZrow[b], relaxdata->ipZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZcol[b], relaxdata->ipZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZval[b], relaxdata->ipZnblocknonz[b]) );

               /* compute Z matrix */
               SCIP_CALL( SCIPconsSdpComputeSparseSdpMatrix(scip, sdpblocks[b], relaxdata->ipy, &(relaxdata->ipZnblocknonz[b]), relaxdata->ipZrow[b], relaxdata->ipZcol[b], relaxdata->ipZval[b]) );

               assert( relaxdata->ipZnblocknonz[b] <= arraylength );

               if ( relaxdata->ipZnblocknonz[b] < arraylength )
               {
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipZrow[b], arraylength, relaxdata->ipZnblocknonz[b]) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipZcol[b], arraylength, relaxdata->ipZnblocknonz[b]) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &relaxdata->ipZval[b], arraylength, relaxdata->ipZnblocknonz[b]) );
               }

               /* TODO check if they are already sorted (sorting is needed since they will later be merged into warmstart arrays) */
               SCIPsdpVarfixerSortRowCol(relaxdata->ipZrow[b], relaxdata->ipZcol[b], relaxdata->ipZval[b], relaxdata->ipZnblocknonz[b]);
            }

            /* compute LP block */
            SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZrow[b], 2 * nrows + 2 * nvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZcol[b], 2 * nrows + 2 * nvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZval[b], 2 * nrows + 2 * nvars) );

            /* for the analytic center all the entries should be strictly positive */
            relaxdata->ipZnblocknonz[b] = 2 * nrows + 2 * nvars;

            for (r = 0; r < nrows; r++)
            {
               /* compute row value for current solution */
               rowval = 0.0;
               rownnonz = SCIProwGetNNonz(rows[r]);
               rowvals = SCIProwGetVals(rows[r]);
               rowcols = SCIProwGetCols(rows[r]);
               for (i = 0; i < rownnonz; i++)
                  rowval += SCIPgetSolVal(scip, relaxdata->ipy, SCIPcolGetVar(rowcols[i])) * rowvals[i];

               relaxdata->ipZrow[b][2*r] = 2*r;
               relaxdata->ipZcol[b][2*r] = 2*r;
               relaxdata->ipZval[b][2*r] = rowval - (SCIProwGetLhs(rows[r]) - SCIProwGetConstant(rows[r]));
               relaxdata->ipZrow[b][2*r + 1] = 2*r + 1;
               relaxdata->ipZcol[b][2*r + 1] = 2*r + 1;
               relaxdata->ipZval[b][2*r + 1] = SCIProwGetRhs(rows[r]) - SCIProwGetConstant(rows[r]) - rowval;
            }

            for (v = 0; v < nvars; v++)
            {
               relaxdata->ipZrow[b][2*nrows + 2*v] = 2*nrows + 2*v;
               relaxdata->ipZcol[b][2*nrows + 2*v] = 2*nrows + 2*v;
               relaxdata->ipZval[b][2*nrows + 2*v] = SCIPgetSolVal(scip, relaxdata->ipy, vars[v]) - SCIPvarGetLbLocal(vars[v]);
               relaxdata->ipZrow[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
               relaxdata->ipZcol[b][2*nrows + 2*v + 1] = 2*nrows + 2*v + 1;
               relaxdata->ipZval[b][2*nrows + 2*v + 1] = SCIPvarGetUbLocal(vars[v]) - SCIPgetSolVal(scip, relaxdata->ipy, vars[v]);
            }
#ifdef SCIP_PRINT_WARMSTART
            for (b = 0; b < relaxdata->nblocks - 1; b++)
            {
               SCIPdebugMessage("dual matrix, block %d:\n", b);
               for (i = 0; i < relaxdata->ipZnblocknonz[b]; i++)
               {
                  SCIPdebugMessage("Z_%d[%d,%d]: %f\n", b, relaxdata->ipZrow[b][i], relaxdata->ipZcol[b][i], relaxdata->ipZval[b][i]);
               }
            }
            SCIPdebugMessage("dual matrix, LP constraints:\n");
            for (r = 0; r < nrows; r++)
            {
               SCIPdebugMessage("Z_%d[%d,%d]: %f\n", relaxdata->nblocks, relaxdata->ipZrow[b][2*r], relaxdata->ipZcol[b][2*r], relaxdata->ipZval[b][2*r]);
               SCIPdebugMessage("Z_%d[%d,%d]: %f\n", relaxdata->nblocks, relaxdata->ipZrow[b][2*r+1], relaxdata->ipZcol[b][2*r+1], relaxdata->ipZval[b][2*r+1]);
            }
            for (v = 0; v < nvars; v++)
            {
               SCIPdebugMessage("Z_%d[%d,%d]: %f\n", relaxdata->nblocks,
                     relaxdata->ipZrow[b][2*nrows + 2*v], relaxdata->ipZcol[b][2*nrows + 2*v], relaxdata->ipZval[b][2*nrows + 2*v]);
               SCIPdebugMessage("Z_%d[%d,%d]: %f\n", relaxdata->nblocks,
                     relaxdata->ipZrow[b][2*nrows + 2*v + 1], relaxdata->ipZcol[b][2*nrows + 2*v + 1], relaxdata->ipZval[b][2*nrows + 2*v + 1]);
            }
#endif
         }
         else
         {
            /* use a scaled identity matrix (and y=0) if the computation of the dual analytic center failed */
            relaxdata->ipZexists = TRUE;

            /* y is set to the zero vector */
            SCIP_CALL( SCIPcreateSol(scip, &relaxdata->ipy, NULL) );

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZnblocknonz, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZrow, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZcol, relaxdata->nblocks) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZval, relaxdata->nblocks) );

            sdpblocks = SCIPconshdlrGetConss(sdpconshdlr);

            for (b = 0; b < relaxdata->nblocks; b++)
            {
               if ( b < relaxdata->nblocks - 1 )
               {
                  /* SDP block */
                  relaxdata->ipZnblocknonz[b] = SCIPconsSdpGetBlocksize(scip, sdpblocks[b]);
               }
               else
               {
                  /* LP block */
                  relaxdata->ipZnblocknonz[b] = SCIPgetNLPRows(scip) + 2 * SCIPgetNVars(scip);
               }
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZrow[b], relaxdata->ipZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZcol[b], relaxdata->ipZnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &relaxdata->ipZval[b], relaxdata->ipZnblocknonz[b]) );

               for (i = 0; i < relaxdata->ipXnblocknonz[b]; i++)
               {
                  relaxdata->ipZrow[b][i] = i;
                  relaxdata->ipZcol[b][i] = i;
                  relaxdata->ipZval[b][i] = relaxdata->lambdastar;
               }
            }

            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Failed to compute analytic center of dual feasible set, using scaled identity instead.\n");
         }
      }
   }
   return SCIP_OKAY;
}

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent
 *  primal variable.
 */
SCIP_RETCODE SCIPrelaxSdpGetPrimalBoundVars(
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get information for */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   SCIP_CALL( SCIPsdpiGetPrimalBoundVars(relaxdata->sdpi, lbvars, ubvars, arraylength) );

   return SCIP_OKAY;
}

/** returns optimal objective value of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get objective value for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            objval              /**< pointer to store the optimal objective value of the SDP-relaxation */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( objval != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;
   *objval = relaxdata->objval;

   return SCIP_OKAY;
}

/** returns values of all variables in the solution of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_RETCODE SCIPrelaxSdpGetRelaxSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get solution for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            solarray,           /**< pointer to store the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( solarray != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;

   if ( *sollength >= SCIPgetNVars(scip) )
   {
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, NULL, solarray, sollength) );
   }
   else
   {
      SCIPdebugMessage("Called SCIPrelaxSdpGetRelaxSol with an array that wasn't big enough, needed length %d, given %d!\n", SCIPgetNVars(scip), *sollength);
      *sollength = SCIPgetNVars(scip);
   }

   return SCIP_OKAY;
}

/** get the number of the SCIP-node which the current SDP solution belongs to */
long int SCIPrelaxSdpGetSdpNode(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->lastsdpnode;
}

/** Was the original problem solved for the last SDP-node (or a penalty or probing formulation) ? */
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->origsolved && SCIPsdpiSolvedOrig(relaxdata->sdpi);
}

/** Was the last probing SDP solved successfully ? */
SCIP_Bool SCIPrelaxSdpSolvedProbing(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->probingsolved;
}

/** returns whether the last solved problem was feasible */
SCIP_Bool SCIPrelaxSdpIsFeasible(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get feasibility for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->feasible );
}

/** returns total number of SDP-iterations */
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the iterations for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->sdpiterations;
}

/** returns number of SDPs solved by SDP-solver (including multiple calls for penalty formulation etc.) */
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->sdpcalls );
}

/** returns number of solved SDP-relaxations */
int SCIPrelaxSdpGetNSdpInterfaceCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->sdpinterfacecalls );
}

/** returns number of SDP-relaxations solved with fast settings */
int SCIPrelaxSdpGetNSdpFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedfast );
}

/** returns number of SDP-relaxations solved with medium settings */
int SCIPrelaxSdpGetNSdpMedium(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedmedium );
}

/** returns number of SDP-relaxations solved with stable settings */
int SCIPrelaxSdpGetNSdpStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedstable );
}

/** returns number of SDP-relaxations solved with penalty formulation */
int SCIPrelaxSdpGetNSdpPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedpenalty );
}

/** returns number of SDP-relaxations unsolved even when using a penalty formulation */
int SCIPrelaxSdpGetNSdpUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unsolved );
}

/** returns number of SDP-relaxations for which dual Slater condition held */
int SCIPrelaxSdpGetNdualSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->ndslaterholds );
}

/** returns number of SDP-relaxations for which dual Slater condition failed */
int SCIPrelaxSdpGetNdualSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->ndnoslater );
}

/** returns number of SDP-relaxations for which dual Slater condition showed infeasibility */
int SCIPrelaxSdpGetNdualSlaterInfeasible(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->nslaterinfeasible );
}

/** returns number of SDP-relaxations for which dual Slater condition could not be determined */
int SCIPrelaxSdpGetNdualSlaterUnknown(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->ndslatercheckfailed );
}

/** returns number of SDP-relaxations for which primal Slater condition held */
int SCIPrelaxSdpGetNprimalSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->npslaterholds );
}

/** returns number of SDP-relaxations for which primal Slater condition failed */
int SCIPrelaxSdpGetNprimalSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->npnoslater );
}

/** returns number of SDP-relaxations for which primal Slater condition could not be determined */
int SCIPrelaxSdpGetNprimalSlaterUnknown(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->npslatercheckfailed );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual */
int SCIPrelaxSdpGetNSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->nslaterholds );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with fastest settings */
int SCIPrelaxSdpGetNSlaterHoldsFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->stablewslater );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with stable settings */
int SCIPrelaxSdpGetNSlaterHoldsStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unstablewslater );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with penalty formulation */
int SCIPrelaxSdpGetNSlaterHoldsPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->penaltywslater );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, for which an infeasible lower bound could be computed */
int SCIPrelaxSdpGetNSlaterHoldsBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->boundedwslater );
}

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, unsolved even when using a penalty formulation */
int SCIPrelaxSdpGetNSlaterHoldsUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unsolvedwslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual */
int SCIPrelaxSdpGetNSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->nnoslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with fast settings */
int SCIPrelaxSdpGetNSlaterFailsFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->stablenoslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with stable settings */
int SCIPrelaxSdpGetNSlaterFailsStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unstablenoslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with penalty formulation */
int SCIPrelaxSdpGetNSlaterFailsPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->penaltynoslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, for which an infeasible lower bound could be computed */
int SCIPrelaxSdpGetNSlaterFailsBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->boundednoslater );
}

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, unsolved even when using a penalty formulation */
int SCIPrelaxSdpGetNSlaterFailsUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unsolvednoslater );
}

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with fast settings */
int SCIPrelaxSdpGetNSlaterInfeasibleFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->stableinfeasible );
}

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with stable settings */
int SCIPrelaxSdpGetNSlaterInfeasibleStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unstableinfeasible );
}

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with penalty formulation */
int SCIPrelaxSdpGetNSlaterInfeasiblePenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->penaltyinfeasible );
}

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, for which an infeasible lower bound could be computed */
int SCIPrelaxSdpGetNSlaterInfeasibleBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->boundedinfeasible );
}

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, unsolved even when using a penalty formulation */
int SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->unsolvedinfeasible );
}
