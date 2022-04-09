/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_sdpfracround.c
 * @brief  fractional rounding heuristic for SDPs
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_sdpfracround.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define HEUR_NAME             "sdpfracround"
#define HEUR_DESC             "fractional rounding heuristic for SDPs"
#define HEUR_DISPCHAR         '^'
#define HEUR_PRIORITY         -1000500
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /* does the heuristic use a secondary SCIP instance? */


/*
 * Default parameter settings
 */

#define DEFAULT_RUNFORLP              FALSE  /**< Should fractional rounding be applied if we are solving LPs? */

/* locally defined heuristic data */
struct SCIP_HeurData
{
   SCIP_SOL*             sol;                /**< working solution */
   SCIP_Bool             runforlp;           /**< Should fractional rounding be applied if we are solving LPs? */
};


/*
 * Callback methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopySdpfracround)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurSdpFracround(scip) );

   return SCIP_OKAY;
}

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
static
SCIP_DECL_HEURFREE(heurFreeSdpfracround)
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

/** initialization method of primal heuristic (called after problem was transformed) */
static
SCIP_DECL_HEURINIT(heurInitSdpfracround)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* create working solution */
   SCIP_CALL( SCIPcreateSol(scip, &heurdata->sol, heur) );

   return SCIP_OKAY;
}

/** deinitialization method of primal heuristic (called before transformed problem is freed) */
static
SCIP_DECL_HEUREXIT(heurExitSdpfracround)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* free working solution */
   SCIP_CALL( SCIPfreeSol(scip, &heurdata->sol) );

   return SCIP_OKAY;
}

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecSdpfracround)
{  /*lint --e{715}*/
   SCIP_HEURDATA* heurdata;
   SCIP_CONSHDLR* conshdlrsdp;
   SCIP_RELAX* relaxsdp;
   SCIP_Real* sdpcandssol;
   SCIP_Real* sdpcandsfrac;
   SCIP_VAR** vars;
   SCIP_VAR** sdpcands;
   SCIP_SOL* relaxsol = NULL;
   SCIP_Bool usesdp = TRUE;
   SCIP_Bool cutoff = FALSE;
   int nsdpcands = 0;
   int ncontvars;
   int freq = -1;
   int nfixed = 0;
   int nrounded = 0;
   int nvars;
   int v;

   assert( heur != NULL );
   assert( strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0 );
   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DELAYED;

   /* do not call heuristic if node was already detected to be infeasible */
   if ( nodeinfeasible )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTRUN;

   /* get heuristic data */
   heurdata = SCIPheurGetData(heur);
   assert( heurdata != NULL );

   /* do not run if relaxation solution is not available and we do not want to run for LPs or no LP solution is available */
   if ( ! SCIPisRelaxSolValid(scip) )
   {
      /* exit if we do not want to run for LPs */
      if ( ! heurdata->runforlp )
         return SCIP_OKAY;

      /* exit if LP is not solved */
      if ( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
         return SCIP_OKAY;

      /* avoid solving for sub-SCIPs */
      if ( SCIPgetSubscipDepth(scip) > 0 )
         return SCIP_OKAY;

      usesdp = FALSE;
   }

   /* get relaxator - exit if not found */
   relaxsdp = SCIPfindRelax(scip, "SDP");
   if ( relaxsdp == NULL )
      return SCIP_OKAY;

   conshdlrsdp = SCIPfindConshdlr(scip, "SDP");
   if ( conshdlrsdp == NULL )
      return SCIP_OKAY;

   /* exit if there are no SDP constraints */
   if ( SCIPconshdlrGetNConss(conshdlrsdp) <= 0 )
      return SCIP_OKAY;

   /* get number of continuous variables */
   ncontvars = SCIPgetNContVars(scip) +  SCIPgetNImplVars(scip);

   /* decide which solution to use */
   if ( usesdp )
   {
      SCIP_CALL( SCIPcreateRelaxSol(scip, &relaxsol, heur) );
   }

   /* prepare probing mode */
   SCIP_CALL( SCIPstartProbing(scip) );

   /* get SDP/LP solution */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcands, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcandssol, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcandsfrac, nvars) );

   /* prepare solution to be changed */
   SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );

   /* collect fractional unfixed values */
   for (v = 0; v < nvars; ++v)
   {
      SCIP_Real val;
      SCIP_VAR* var;

      var = vars[v];
      if ( SCIPvarIsIntegral(var) )
      {
         val = SCIPgetSolVal(scip, relaxsol, var);
         if ( ! SCIPisFeasIntegral(scip, val) )
         {
            sdpcandssol[nsdpcands] = val;
            sdpcandsfrac[nsdpcands] = SCIPfeasFrac(scip, val);
            sdpcands[nsdpcands++] = var;
         }
         else
         {
            /* make sure value is really integral */
            val = SCIPfeasRound(scip, val);

            /* fixing variable to integral value */
            if ( SCIPisGT(scip, val, SCIPvarGetLbLocal(var)) )
            {
               SCIP_CALL( SCIPchgVarLbProbing(scip, var, val) );
               ++nfixed;
            }
            if ( SCIPisLT(scip, val, SCIPvarGetUbLocal(var)) )
            {
               SCIP_CALL( SCIPchgVarUbProbing(scip, var, val) );
               ++nfixed;
            }

            /* to avoid numerical noise, make sure variable integral */
            SCIP_CALL( SCIPsetSolVal(scip, heurdata->sol, var, val) );
         }
      }
   }

   /* possibly free relaxation (LP or SDP) solution */
   if ( relaxsol != NULL )
   {
      SCIP_CALL( SCIPfreeSol(scip, &relaxsol) );
   }

   /* do not proceed, if there are no fractional variables */
   if ( nsdpcands == 0 )
   {
      SCIP_CALL( SCIPendProbing(scip) );

      SCIPfreeBufferArray(scip, &sdpcandsfrac);
      SCIPfreeBufferArray(scip, &sdpcandssol);
      SCIPfreeBufferArray(scip, &sdpcands);

      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;

   SCIPdebugMsg(scip, "Node %"SCIP_LONGINT_FORMAT": executing SDP fractional rounding heuristic (depth %d, %d fractionals).\n",
      SCIPgetNNodes(scip), SCIPgetDepth(scip), nsdpcands);
   SCIPdebugMsg(scip, "Fixed %d bounds of variables with integral values.\n", nfixed);

   if ( ! usesdp )
   {
      /* temporarily change relaxator frequency, since otherwise relaxation will not be solved */
      freq = SCIPrelaxGetFreq(relaxsdp);
      SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );
   }

   /* sort according to increasing fractinal parts */
   SCIPsortRealRealPtr(sdpcandsfrac, sdpcandssol, (void*) sdpcands, nsdpcands);

   /* perform rounding loop */
   nfixed = 0;
   for (v = 0; v < nsdpcands && ! cutoff; ++v)
   {
      SCIP_Longint ndomreds;
      SCIP_VAR* var;
      SCIP_Real newval;
      SCIP_Real val;
      SCIP_Real ceilval;
      SCIP_Real floorval;
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Bool lbadjust;
      SCIP_Bool ubadjust;

      /* get next variable from permuted candidate array */
      var = sdpcands[v];
      val = sdpcandssol[v];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      assert( SCIPvarIsIntegral(var) );
      assert( ! SCIPisFeasIntegral(scip, val) );

      ceilval = SCIPfeasCeil(scip, val);
      floorval = SCIPfeasFloor(scip, val);

      /* Abort if rounded ceil and floor value lie outside the variable domain. Otherwise, check if bounds allow only
       * one rounding direction, anyway */
      if ( lb > ceilval + 0.5 || ub < floorval - 0.5 )
      {
         cutoff = TRUE;
         break;
      }
      else if ( SCIPisFeasEQ(scip, lb, ceilval) )
         newval = ceilval;
      else if ( SCIPisFeasEQ(scip, ub, floorval) )
         newval = floorval;
      else
      {
         newval = SCIPfeasRound(scip, val);
         ++nrounded;
      }

      lbadjust = SCIPisGT(scip, newval, lb);
      ubadjust = SCIPisLT(scip, newval, ub);

      if ( lbadjust || ubadjust )
      {
         SCIP_CALL( SCIPnewProbingNode(scip) );

         /* tighten the bounds to fix the variable for the probing node */
         if ( lbadjust )
         {
            SCIP_CALL( SCIPchgVarLbProbing(scip, var, newval) );
         }

         if ( ubadjust )
         {
            SCIP_CALL( SCIPchgVarUbProbing(scip, var, newval) );
         }
         ++nfixed;

         /* call propagation routines for the reduced problem */
         SCIP_CALL( SCIPpropagateProbing(scip, 1, &cutoff, &ndomreds) );
      }

      /* change solution value */
      if ( ! cutoff )
      {
         SCIP_CALL( SCIPsetSolVal(scip, heurdata->sol, var, newval) );
      }
   }

   /* check solution */
   if ( ! cutoff )
   {
      SCIP_Bool success;

      SCIPdebugMsg(scip, "Rounded %d variables.\n", nrounded);

      /* if there are no continuous variables, we can just try the solution */
      if ( ncontvars == 0 )
      {
#ifndef NDEBUG
         /* assert that solution is really integral */
         for (v = 0; v < nvars; ++v)
            assert( ! SCIPvarIsIntegral(vars[v]) || SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, heurdata->sol, vars[v])) );
#endif

         /* try to add solution to SCIP - do not need to check integrality here */
         SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

         if ( success )
         {
            SCIPdebugMsg(scip, "Found solution for full integral instance.\n");
            *result = SCIP_FOUNDSOL;
         }
         else
            SCIPdebugMsg(scip, "Solution not feasible for full integral instance.\n");
      }
      else if ( nfixed > 0 )
      {
         /* if there are continuous variables, we need to solve a final SDP */
         SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

         /* if solving was successfull */
         if ( SCIPrelaxSdpSolvedProbing(relaxsdp) )
         {
            if ( SCIPrelaxSdpIsFeasible(relaxsdp) )
            {
               /* check solution */
               SCIP_CALL( SCIPlinkRelaxSol(scip, heurdata->sol) );

               /* try to add solution to SCIP: check all constraints, including integrality */
               SCIP_CALL( SCIPtrySol(scip, heurdata->sol, FALSE, TRUE, TRUE, TRUE, TRUE, &success) );

               /* check, if solution was feasible and good enough */
               if ( success )
               {
                  SCIPdebugMsg(scip, "Solution was feasible and good enough.\n");
                  *result = SCIP_FOUNDSOL;
               }
               else
                  SCIPdebugMsg(scip, "Solution was not feasible.\n");
            }
            else
               SCIPdebugMsg(scip, "Problem was infeasible.\n");
         }
      }
      else
         SCIPdebugMsg(scip, "No fixings have been performed.\n");
   }
   else
      SCIPdebugMsg(scip, "Reached cutoff after %d roundings.\n", nrounded);

   /* free local problem */
   SCIP_CALL( SCIPendProbing(scip) );

   /* reset frequency of relaxator */
   if ( ! usesdp )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", freq) );
   }

   SCIPfreeBufferArray(scip, &sdpcandsfrac);
   SCIPfreeBufferArray(scip, &sdpcandssol);
   SCIPfreeBufferArray(scip, &sdpcands);

   SCIPdebugMsg(scip, "Finished fractional rounding heuristic.\n");

   return SCIP_OKAY;
}


/*
 * heuristic specific interface methods
 */

/** creates the fractional rounding heuristic for SDPs and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurSdpFracround(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create fractional rounding primal heuristic data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &heurdata) );

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecSdpfracround, heurdata) );

   assert( heur != NULL );

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopySdpfracround) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeSdpfracround) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitSdpfracround) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitSdpfracround) );

   /* fractional rounding heuristic parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "heuristics/sdpfracround/runforlp",
         "Should fractional rounding be applied if we are solving LPs?",
         &heurdata->runforlp, FALSE, DEFAULT_RUNFORLP, NULL, NULL) );

   return SCIP_OKAY;
}
