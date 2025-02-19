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

/**@file   heur_sdpinnerlp.c
 * @brief  set up inner approximation LP formulation and run heuristics
 * @author Marc Pfetsch
 *
 * The idea of this heuristic is to check whether there is a diagonally dominant solution via an
 * LP. This is described in@n
 * Optimization over structured subsets of positive semidefinite matrices via column generation@n
 * Amir Ali Ahmadi, Sanjeeb Dash, and Georgina Hall
 * Discrete Optimization 24 (2017) 129-151
 *
 * Assume that an SDP-constraint is given as \f$ \sum_{j=1}^n A_j y_j - A_0 \succeq 0 \f$. We then consider a set of
 * rank-1 SDP matrices \f$B_1, \dots, B_t\f$ and write
 * \f[
 *   \sum_{j=1}^n A_j y_j - A_0 = \sum_{i=1}^t B_i \alpha_i,\; \alpha \geq 0.
 * \f]
 * We use all \f$n^2\f$ matrices \f$B_i\f$ that arise as \f$u u^T\f$, where \f$u\f$ has at most 2 nonzero entries \f$\pm
 * 1\f$. It can be shown that this choice suffices to represent all diagonally dominant matrices. The given constraints
 * are actually linear in the \f$y\f$-variables. We set up a MIP and solve it. If it is feasible, it should provide a
 * feasible solution for the original problem. However, often the problem is infeasible or only provides a weak
 * solution.
 *
 * We currently do not use the column generation aspect of the paper mentioned above.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_sdpinnerlp.h"
#include "cons_sdp.h"
#include "scip/cons_linear.h"

#define HEUR_NAME             "sdpinnerlp"
#define HEUR_DESC             "inner approximation LP heuristic for SDPs"
#define HEUR_DISPCHAR         '!'
#define HEUR_PRIORITY         -1001000
#define HEUR_FREQ             -1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      TRUE  /* does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_STALLNODELIMIT          100L      /**< limit on number of nodes since last improving incumbent solutions */
#define DEFAULT_MAXSIZE                10000      /**< maximal size of the inner problem */


/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_Longint          stallnodelimit;     /**< limit on number of nodes since last improving incumbent solutions */
   int                   maxsize;            /**< maximal size of the inner problem */
};


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySdpInnerlp)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSdpInnerlp(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSdpInnerlp)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );

   /* free heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL);

   SCIPfreeBlockMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSdpInnerlp)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_HASHMAP* varmapfw;
   SCIP_CONSHDLR* conshdlrsdp;
   SCIP_CONS** conss;
   SCIP_VAR** subvars;
   SCIP_VAR** vars;
   SCIP_Bool success;
   SCIP_Real timelimit;
   SCIP_SOL** subsols;
   SCIP* subscip;
   int totalsize = 0;
   int nsubsols;
   int nconss;
   int nvars;
   int c;
   int i;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic if node was already detected to be infeasible */
   if ( nodeinfeasible )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* compute time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      timelimit = MAX(0, timelimit - SCIPgetSolvingTime(scip) );
   }
   if ( timelimit <= 0.01 )
      return SCIP_OKAY;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* estimate size of problem */
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   /* find SDP constraint handler */
   conshdlrsdp = SCIPfindConshdlr(scip, "SDP");
   if ( conshdlrsdp == NULL )
      return SCIP_OKAY;

   for (c = 0; c < nconss && totalsize < heurdata->maxsize; ++c)
   {
      int blocksize;

      /* skip non-SDP constraints */
      assert( conss[c] != NULL );
      if ( SCIPconsGetHdlr(conss[c]) != conshdlrsdp )
         continue;

      blocksize = SCIPconsSdpGetBlocksize(scip, conss[c]);
      totalsize += (blocksize * (blocksize - 1))/2;
   }

   if ( totalsize >= heurdata->maxsize )
   {
      SCIPdebugMsg(scip, "Skipping <%s>, because size would be too large.\n", SCIPheurGetName(heur));
      return SCIP_OKAY;
   }

   /* get original variable data */
   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );

   /* create subscip */
   SCIP_CALL( SCIPcreate(&subscip) );

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPhashmapCreate(&varmapfw, SCIPblkmem(subscip), nvars) );

   /* create a problem copy as sub SCIP */
   SCIP_CALL( SCIPcopyLargeNeighborhoodSearch(scip, subscip, varmapfw, "sdpinnerlp", NULL, NULL, 0, FALSE,
         FALSE, &success, NULL) );

   /* find SDP constraint handler */
   conshdlrsdp = SCIPfindConshdlr(subscip, "SDP");
   if ( conshdlrsdp == NULL )
   {
      SCIP_CALL( SCIPfree(&subscip) );
      return SCIP_OKAY;
   }

   /* copy subproblem variables into the same order as the source SCIP variables */
   SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
   for( i = 0; i < nvars; i++ )
      subvars[i] = (SCIP_VAR*) SCIPhashmapGetImage(varmapfw, vars[i]);

   /* free hash map */
   SCIPhashmapFree(&varmapfw);

   /* get orginal constraints: we have to copy them, because we are adding constraints which possibly leads to a reallocation */
   nconss = SCIPgetNConss(subscip);
   SCIP_CALL( SCIPduplicateBufferArray(scip, &conss, SCIPgetConss(subscip), nconss) );

   /* loop through all constraints */
   for (c = 0; c < nconss; ++c)
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_Real** matrices = NULL;
      SCIP_Real* constmatrix;
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_VAR** rayvars;
      SCIP_CONS* cons;
      SCIP_VAR** sdpvars;
      SCIP_Real val;
      int blocksize;
      int nsdpvars;
      int s;
      int t;

      assert( conss[c] != NULL );

      /* skip non-SDP constraints */
      if ( SCIPconsGetHdlr(conss[c]) != conshdlrsdp )
         continue;

      blocksize = SCIPconsSdpGetBlocksize(subscip, conss[c]);
      nsdpvars = SCIPconsSdpGetNVars(subscip, conss[c]);
      sdpvars = SCIPconsSdpGetVars(subscip, conss[c]);

      /* get matrices */
      SCIP_CALL( SCIPallocBufferArray(subscip, &constmatrix, blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(subscip, conss[c], constmatrix) );

      SCIP_CALL( SCIPallocBufferArray(subscip, &matrices, nsdpvars) );
      for (i = 0; i < nsdpvars; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(subscip, &matrices[i], blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullAj(subscip, conss[c], i, matrices[i]) );
      }

      SCIP_CALL( SCIPallocBufferArray(subscip, &consvars, nsdpvars + 2 * blocksize) );
      SCIP_CALL( SCIPallocBufferArray(subscip, &consvals, nsdpvars + 2 * blocksize) );
      SCIP_CALL( SCIPallocBufferArray(subscip, &rayvars, blocksize * blocksize) );

      /* Create ray variables: Variable rayvars[s * blocksize + t] corresponds to a rank-1 matrix. The submatrix indexed
       * by (s,t) is (1,1;1,1) or (1,-1;-1,1) depending on whether s < t or s > t. If s = t, we have a diagonal matrix
       * with a 1. */
      for (s = 0; s < blocksize; ++s)
      {
         for (t = 0; t <= s; ++t)
         {
            SCIP_Bool nonzero = FALSE;

            /* check whether entry (s,t) has a nonzero somewhere - otherwise we do not need variables or constraints */
            if ( SCIPisZero(subscip, constmatrix[s * blocksize + t]) )
            {
               for (i = 0; i < nsdpvars; ++i)
               {
                  val = matrices[i][s * blocksize + t];
                  if ( ! SCIPisZero(subscip, val) )
                  {
                     nonzero = TRUE;
                     break;
                  }
               }
            }
            else
               nonzero = TRUE;

            if ( nonzero || s == t )
            {
               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ray%d#%d", s, t);
               SCIP_CALL( SCIPcreateVarBasic(subscip, &rayvars[s * blocksize + t], name, 0.0, SCIPinfinity(subscip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
               SCIP_CALL( SCIPaddVar(subscip, rayvars[s * blocksize + t]) );

               if ( s != t )
               {
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "ray%d#%d", t, s);
                  SCIP_CALL( SCIPcreateVarBasic(subscip, &rayvars[t * blocksize + s], name, 0.0, SCIPinfinity(subscip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
                  SCIP_CALL( SCIPaddVar(subscip, rayvars[t * blocksize + s]) );
               }
            }
            else
            {
               rayvars[s * blocksize + t] = NULL;
               rayvars[t * blocksize + s] = NULL;
            }
         }
      }

      /* loop over all possible entries */
      for (s = 0; s < blocksize; ++s)
      {
         for (t = 0; t <= s; ++t)
         {
            int cnt = 0;

            /* skip 0-entries */
            if ( rayvars[s * blocksize + t] == NULL )
               continue;

            assert( rayvars[t * blocksize + s] != NULL );

            /* add entries for matrices */
            for (i = 0; i < nsdpvars; ++i)
            {
               val = matrices[i][s * blocksize + t];
               if ( ! SCIPisZero(subscip, val) )
               {
                  consvars[cnt] = sdpvars[i];
                  consvals[cnt++] = val;
               }
            }

            /* add ray variables: -1 because we have to bring the variables to the LHS */
            consvars[cnt] = rayvars[s * blocksize + t];
            consvals[cnt++] = -1.0;

            if ( s != t )
            {
               consvars[cnt] = rayvars[t * blocksize + s];
               consvals[cnt++] = +1.0;
            }
            else
            {
               /* add all off-diagonal variables */
               for (i = 0; i < blocksize; ++i)
               {
                  if ( i != s )
                  {
                     if ( rayvars[s * blocksize + i] != NULL )
                     {
                        assert( rayvars[i * blocksize + s] != NULL );

                        consvars[cnt] = rayvars[s * blocksize + i];
                        consvals[cnt++] = -1.0;

                        consvars[cnt] = rayvars[i * blocksize + s];
                        consvals[cnt++] = -1.0;
                     }
                  }
               }
            }

            /* add linear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin%d#%d", s, t);
            SCIP_CALL( SCIPcreateConsLinear(subscip, &cons, name, cnt, consvars, consvals, constmatrix[s * blocksize + t], constmatrix[s * blocksize + t],
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(subscip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIP_CALL( SCIPprintCons(subscip, cons, NULL) );
#endif
            SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         }
      }

      /* delete SDP constraint */
      SCIP_CALL( SCIPdelCons(subscip, conss[c]) );

      for (i = 0; i < blocksize * blocksize; ++i)
      {
         if ( rayvars[i] != NULL )
         {
            SCIP_CALL( SCIPreleaseVar(subscip, &rayvars[i]) );
         }
      }
      SCIPfreeBufferArray(subscip, &rayvars);
      SCIPfreeBufferArray(subscip, &consvals);
      SCIPfreeBufferArray(subscip, &consvars);

      for (i = nsdpvars-1; i >= 0; --i)
         SCIPfreeBufferArray(subscip, &matrices[i]);
      SCIPfreeBufferArray(subscip, &matrices);
      SCIPfreeBufferArray(subscip, &constmatrix);
   }

#if 0
   SCIP_CALL( SCIPwriteOrigProblem(subscip, "debug.lp", "lp", FALSE) );
#endif

   /* set individual time limit */
   if ( ! SCIPisInfinity(scip, timelimit) )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   }

   /* set stall node limit */
   SCIP_CALL( SCIPsetLongintParam(subscip, "limits/stallnodes", heurdata->stallnodelimit) );

#ifdef SCIP_MORE_DEBUG
   SCIPinfoMessage(scip, NULL, "\nSolving inner LP subproblem ...\n");
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 5) );
#else
   SCIPdebugMsg(scip, "Solving inner LP subproblem ...\n");
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );
#endif

   /* turn off recursive use */
   SCIP_CALL( SCIPsetIntParam(subscip, "heuristics/sdpinnerlp/freq", -1) );

   SCIP_CALL( SCIPsolve(subscip) );

   /* check, whether a solution was found */
   nsubsols = SCIPgetNSols(subscip);
   subsols = SCIPgetSols(subscip);
   success = FALSE;
   for (i = 0; i < nsubsols && ! success; ++i)
   {
      SCIP_SOL* newsol;

      SCIP_CALL( SCIPtranslateSubSol(scip, subscip, subsols[i], heur, subvars, &newsol) );

      SCIP_CALL( SCIPtrySolFree(scip, &newsol, FALSE, FALSE, TRUE, TRUE, TRUE, &success) );
      if ( success )
      {
         SCIPdebugMsg(scip, "Found solution.\n");
         *result = SCIP_FOUNDSOL;
      }
   }

   SCIP_CALL( SCIPfree(&subscip) );

   SCIPfreeBufferArray(scip, &subvars);
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the innerlp heuristic for SDPs and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSdpInnerlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create innerlp primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSdpInnerlp, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySdpInnerlp) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSdpInnerlp) );

   SCIP_CALL( SCIPaddLongintParam(scip, "heuristics/" HEUR_NAME "/stallnodelimit",
         "limit on number of nodes since last improving incumbent solutions",
         &heurdata->stallnodelimit, FALSE, DEFAULT_STALLNODELIMIT, -1LL, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "heuristics/" HEUR_NAME "/maxsize",
         "maximal size of the inner problem",
         &heurdata->maxsize, FALSE, DEFAULT_MAXSIZE, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
