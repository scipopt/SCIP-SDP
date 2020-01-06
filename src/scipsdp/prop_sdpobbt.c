/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2019 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2019 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_sdpobbt.c
 * @brief  optimization-based bound tightening propagator for semidefinite programs
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/

#include <assert.h>
#include <string.h>

#include "prop_sdpobbt.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/* fundamental propagator properties */
#define PROP_NAME             "sdp-obbt"
#define PROP_DESC             "optimization-based bound tightening for SDPs"
#define PROP_PRIORITY         -1100000 /**< propagator priority */
#define PROP_FREQ             -1 /**< propagator frequency */
#define PROP_DELAY            FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING           SCIP_PROPTIMING_AFTERLPLOOP /**< propagation timing mask */

#define DEFAULT_PROPBIN       FALSE /**< should obbt be done for binary variables ? */
#define DEFAULT_PROPCONT      TRUE /**< should obbt be done for continuous variables ? */

/* TODO: maybe make this a parameter and/or have different values for integral and continuous variables */
#define TOLERANCE_FACTOR      2000 /**< only if the improvement is larger than TOLERANCE_FACTOR * gaptol it will be accepted */

/* TODO: fix memory error */


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             propbin;            /**< should obbt be done for binary variables ? */
   SCIP_Bool             propcont;           /**< should obbt be done for continuous variables ? */
   SCIP_Bool             delayed;            /**< did we delay the last call? */
   long long int         lastnode;           /**< the last node we ran for */
   SCIP_Real             lastcufoffbound;    /**< the last cutoffbound we ran for */
   SCIP_Real             sdpsolvergaptol;    /**< gap tolerance of the underlying SDP solver */
   SCIP_Bool             propenabled;        /**< whether propagator is enabled */
};


/*
 * Local methods
 */

/* the current bugfix branch (3.2.1) does not have SCIPsolveProbingRelax() -> do nothing */
#if ( (SCIP_VERSION > 321 || SCIP_SUBVERSION > 0) )
static
SCIP_RETCODE addObjCutoff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_ROW* row;
   SCIP_VAR** vars;
   char rowname[SCIP_MAXSTRLEN];
   int nvars;
   int v;

   assert( scip != NULL );
   assert( SCIPinProbing(scip) );

   SCIPdebugMsg(scip, "create objective cutoff and add it to the LP-constraints\n");

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* create objective cutoff row; set local flag to FALSE since primal cutoff is globally valid */
   (void) SCIPsnprintf(rowname, SCIP_MAXSTRLEN, "obbtsdp_objcutoff");
   SCIP_CALL( SCIPcreateEmptyRowUnspec(scip, &row, rowname, -SCIPinfinity(scip), SCIPgetCutoffbound(scip), FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

   for( v = 0; v < nvars; v++ )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, row, vars[v], SCIPvarGetObj(vars[v])) );
   }
   SCIP_CALL( SCIPflushRowExtensions(scip, row) );

   /* add row to the LP-constraints */
   SCIP_CALL( SCIPaddRowProbing(scip, row) );

   SCIP_CALL( SCIPreleaseRow(scip, &row) );

   return SCIP_OKAY;
}
#endif

/*
 * Callback methods of propagator
 */


/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopySdpObbt)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( prop != NULL );
   assert( strcmp(SCIPpropGetName(prop), PROP_NAME) == 0 );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropSdpObbt(scip) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSdpObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitSdpObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);

   propdata->lastnode = -1; /* we reset this to be able to run again if a new problem is read */

   return SCIP_OKAY;
}

/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolSdpObbt)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);

   if ( SCIPfindRelax(scip, "SDP") == NULL )
      propdata->propenabled = FALSE;
   else
   {
      propdata->propenabled = TRUE;
      SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolvergaptol", &(propdata->sdpsolvergaptol)) );
   }

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSdpObbt)
{  /*lint --e{715}*/
   /* the current bugfix branch (3.2.1) does not have SCIPsolveProbingRelax() -> do nothing */
#if ( (SCIP_VERSION > 321 || SCIP_SUBVERSION > 0) )
   int nvars;
   SCIP_VAR** vars;
   int v;
   SCIP_PROPDATA* propdata;
   SCIP_Real relaxval;
   SCIP_Bool cutoff;
   SCIP_Bool oldobjlimitparam;
   SCIP_Real probingval;
   SCIP_Bool success;
   SCIP_RELAX* relaxsdp;
   /* newbounds and newboundinds save the bound tightenings that should be inserted after probing ends, newboundinds saves the bounds the entries of newbounds
    * belong to, for this the variables are sorted 1 to nvars (and entry i means vars[i-1]), with negative entry for the lower bound and positive for the upper */
   SCIP_Real* newbounds;
   int* newboundinds;
   int nnewbounds;
   int i;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );

   propdata = SCIPpropGetData(prop);

   assert( propdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( ! propdata->propenabled )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Executing propExecSdpObbt! \n");

   /* do not run in: presolving, repropagation, probing mode, subscips, if no objective propagation is allowed */
#if ( SCIP_VERSION >= 700 || (SCIP_VERSION >= 602 && SCIP_SUBVERSION > 0) )
   if ( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) || !SCIPallowWeakDualReds(scip) || (SCIPgetSubscipDepth(scip) > 0) )
#else
   if ( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) || !SCIPallowObjProp(scip) || (SCIPgetSubscipDepth(scip) > 0) )
#endif
   {
      SCIPdebugMsg(scip, "Aborting propExecSdpObbt because we are in presolving, repropagation, probing mode, a subscip or no objective "
            "propagation is allowed!\n");
      return SCIP_OKAY;
   }

   /* delay if cutoffbound is infinite or no relaxation solution exists */
   if ( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) || (! SCIPisRelaxSolValid(scip)) )
   {
      /* if we already delayed in the last call, abort to prevent an infinite loop */
      if ( propdata->delayed )
      {
         SCIPdebugMsg(scip, "Aborting propExecSdpObbt since still cutoffbound is infinite or no relaxation solution exists\n");
         return SCIP_OKAY;
      }
      *result = SCIP_DELAYED;
      propdata->delayed = TRUE;
      SCIPdebugMsg(scip, "Delaying propExecSdpObbt since cutoffbound is infinite or no relaxation solution exists\n");
      return SCIP_OKAY;
   }

   /* delay if best solution was found by trivial heuristic (since in this case the cutoffbound will generally not be good enough) or objective propagation onky */
   if ( (SCIPgetBestSol(scip) == NULL) || ((SCIPsolGetHeur(SCIPgetBestSol(scip)) != NULL) && (strcmp(SCIPheurGetName(SCIPsolGetHeur(SCIPgetBestSol(scip))), "trivial") == 0)) )
   {
      /* if we already delayed in the last call, abort to prevent an infinite loop */
      if ( propdata->delayed )
      {
         SCIPdebugMsg(scip, "Aborting propExecSdpObbt since still best solution was found by trivial heuristic or simple objective propagation, which will not be good enough\n");
         return SCIP_OKAY;
      }
      *result = SCIP_DELAYED;
      propdata->delayed = TRUE;
      SCIPdebugMsg(scip, "Delaying propExecSdpObbt since best solution was found by trivial heuristic or simple objective propagation, which will not be good enough\n");
      return SCIP_OKAY;
   }

   if ( (SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode) && (SCIPisEQ(scip, SCIPgetCutoffbound(scip), propdata->lastcufoffbound)) )
   {
      SCIPdebugMsg(scip, "Not running again for node %lld with cutoffbound &f!\n", propdata->lastnode, propdata->lastcufoffbound);
      return SCIP_OKAY;
   }
   else
   {
      propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
      propdata->lastcufoffbound = SCIPgetCutoffbound(scip);
   }

   propdata->delayed = FALSE;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIPdebugMsg(scip, "start probing\n");

   SCIP_CALL( addObjCutoff(scip) );

   /* make sure that we don't use the objective cutoff for the changed objective */
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/objlimit", &oldobjlimitparam) );
   SCIP_CALL( SCIPsetBoolParam(scip, "relaxing/SDP/objlimit", FALSE) );

   /* allocate memory to save bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &newbounds, 2*nvars) );/*lint !e647*/
   SCIP_CALL( SCIPallocBufferArray(scip, &newboundinds, 2*nvars) );/*lint !e647*/

   *result = SCIP_DIDNOTFIND;

   /* set objective coefficients to zero */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[v], 0.0) );
   }

   nnewbounds = 0;

   for (v = 0; v < nvars; v++)
   {
      /* do not propagate binary or continous variables if the corresponding flag is set to false */
      if ( (( ! propdata->propbin ) && SCIPvarIsBinary(vars[v])) || (( ! propdata->propcont ) && ( ! SCIPvarIsIntegral(vars[v]))) )
      {
#ifdef SCIP_MORE_DEBUG
         if ( SCIPvarIsBinary(vars[v]) )
         {
            SCIPdebugMsg(scip, "Skipping binary variable %s\n", SCIPvarGetName(vars[v]));
         }
         else
         {
            SCIPdebugMsg(scip, "Skipping continuous variable %s\n", SCIPvarGetName(vars[v]));
         }
#endif
         continue;
      }

      /* get the value of this variable for the current relaxation */
      relaxval = SCIPgetRelaxSolVal(scip, vars[v]);

      /* only try obbt for the lower bound if it is not tight for the current relaxation's solution */
      if ( SCIPisFeasGT(scip, relaxval, SCIPvarGetLbLocal(vars[v])) )
      {
         /* set the objective to minimize y_v */
         SCIP_CALL( SCIPchgVarObjProbing(scip, vars[v], 1.0) );

         /* solve the probing problem */
         SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

         /* as cutoff doesn't work for relax sdp, we have to check ourselves, if we didn't manage to solve successfully, we abort (as this will
          * probably not get better for the other variables as we only change the objective */
         relaxsdp = SCIPfindRelax(scip, "SDP");

         if (! SCIPrelaxSdpSolvedProbing(relaxsdp))
         {
            SCIPdebugMsg(scip, "Aborting sdp-obbt, as we were unable to solve a probing sdp!\n");
            if ( *result != SCIP_REDUCEDDOM )
               *result = SCIP_DIDNOTRUN;
            break;
         }

         /* if the problem is infeasible, return with cutoff */
         if ( ! SCIPrelaxSdpIsFeasible(relaxsdp) )
         {
            SCIPdebugMsg(scip, "Probing sdp infeasible, so there can't be a better solution for this problem!\n");
            *result = SCIP_CUTOFF;
            break;
         }

         /* only check objective value if problem was bounded */
         if ( ! SCIPrelaxSdpIsUnbounded(relaxsdp) )
         {
            /* check if we managed to tighten the bound */
            success = FALSE; /* this will be ignored, we check solvedProbing instead */
            SCIP_CALL( SCIPrelaxSdpRelaxVal(relaxsdp, &success, &probingval) );

            /* only update if we improved the bound by at least gaptol, everything else might be inexactness of the solver */
            if ( SCIPisGT(scip, probingval - TOLERANCE_FACTOR * propdata->sdpsolvergaptol, SCIPvarGetLbLocal(vars[v])) )
            {
               /* update bound */
               SCIPdebugMsg(scip, "Obbt-Sdp tightened lower bound of variable %s from %f to %f !\n",
                     SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), probingval - propdata->sdpsolvergaptol);

               newbounds[nnewbounds] = probingval;
               newboundinds[nnewbounds] = -1 * (v+1);
               nnewbounds++;
               *result = SCIP_REDUCEDDOM;
            }
   #ifdef SCIP_MORE_DEBUG
            else
            {
               SCIPdebugMsg(scip, "Obbt-Sdp found lower bound of %f for variable %s, worse than old bound %f !\n",
                     probingval, SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]));
            }
   #endif
         }
#ifdef SCIP_MORE_DEBUG
         else
         {
            SCIPdebugMsg(scip, "Obbt-Sdp problem unbounded for variable %s!\n", SCIPvarGetName(vars[v]));
         }
#endif
      }
#ifdef SCIP_MORE_DEBUG
      else
      {
         SCIPdebugMsg(scip, "Skipping obbt for lower bound %f of variable %s, as current relaxation's solution is tight.\n",
               SCIPvarGetLbLocal(vars[v]), SCIPvarGetName(vars[v]));
      }
#endif

      /* only try obbt for the upper bound if it is not tight for the current relaxation's solution */
      if ( SCIPisFeasLT(scip, relaxval, SCIPvarGetUbLocal(vars[v])) )
      {
         /* set the objective to maximize y_v (minimize -y_v) */
         SCIP_CALL( SCIPchgVarObjProbing(scip, vars[v], -1.0) );

         /* solve the probing problem */
         SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

         /* as cutoff doesn't work for relax sdp, we have to check ourselves, if we didn't manage to solve successfully, we abort (as this will
          * probably not get better for the other variables as we only change the objective */
         relaxsdp = SCIPfindRelax(scip, "SDP");

         if (! SCIPrelaxSdpSolvedProbing(relaxsdp))
         {
            SCIPdebugMsg(scip, "Aborting sdp-obbt, as we were unable to solve a probing sdp!\n");
            if ( *result != SCIP_REDUCEDDOM )
               *result = SCIP_DIDNOTRUN;
            goto ENDPROBING;
         }

         /* if the problem is infeasible, return with cutoff */
         if ( ! SCIPrelaxSdpIsFeasible(relaxsdp) )
         {
            SCIPdebugMsg(scip, "Probing sdp infeasible, so there can't be a better solution for this problem!\n");
            *result = SCIP_CUTOFF;
            nnewbounds = 0;
            goto ENDPROBING;
         }

         /* only check objective value if problem was bounded */
         if ( ! SCIPrelaxSdpIsUnbounded(relaxsdp) )
         {
            /* check if we managed to tighten the bound */
            success = FALSE; /* this will be ignored, we check solvedProbing instead */
            SCIP_CALL( SCIPrelaxSdpRelaxVal(relaxsdp, &success, &probingval) );

            /* only update if we improved the bound by at least gaptol, everything else might be inexactness of the solver */
            if ( SCIPisLT(scip, -probingval + TOLERANCE_FACTOR * propdata->sdpsolvergaptol, SCIPvarGetUbLocal(vars[v])) )
            {
               SCIPdebugMsg(scip, "Obbt-Sdp tightened upper bound of variable %s from %f to %f !\n",
                     SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), -probingval + propdata->sdpsolvergaptol);

               newbounds[nnewbounds] = -probingval;
               newboundinds[nnewbounds] = v + 1;
               nnewbounds++;
            }
   #ifdef SCIP_MORE_DEBUG
            else
            {
               SCIPdebugMsg(scip, "Obbt-Sdp found upper bound of %f for variable %s, worse than old bound %f !\n",
                     -probingval, SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]));
            }
   #endif
         }
#ifdef SCIP_MORE_DEBUG
         else
         {
            SCIPdebugMsg(scip, "Obbt-Sdp problem unbounded for variable %s!\n", SCIPvarGetName(vars[v]));
         }
#endif
      }
#ifdef SCIP_MORE_DEBUG
      else
      {
         SCIPdebugMsg(scip, "Skipping obbt for upper bound %f of variable %s, as current relaxation's solution is tight.\n",
               SCIPvarGetUbLocal(vars[v]), SCIPvarGetName(vars[v]));
      }
#endif

      /* reset the objective coefficient to zero for the next variable */
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[v], 0.0) );
   }

   ENDPROBING:
   SCIP_CALL( SCIPendProbing(scip) );
   SCIPdebugMsg(scip, "end probing\n");
   SCIP_CALL( SCIPsetBoolParam(scip, "relaxing/SDP/objlimit", oldobjlimitparam) );

   for (i = 0; i < nnewbounds; i++)
   {
      if ( newboundinds[i] < 0)
      {
         SCIP_CALL( SCIPchgVarLb(scip, vars[-1 * newboundinds[i] - 1], newbounds[i]) ); /*lint !e679*/
         *result = SCIP_REDUCEDDOM;
      }
      else
      {
         /* Check if the (rounded) new upper bound is smaller than the updated lower bound, in that case return cutoff.
          * Note that this can only happen for integer variables, since by construction the lower bound computed by obbt
          * has to be smaller than the upper bound, but it can happen through rounding. */
         if ( SCIPvarIsBinary(vars[newboundinds[i] - 1]) && SCIPisLT(scip, SCIPfeasFloor(scip, newbounds[i]),
               SCIPvarGetLbLocal(vars[newboundinds[i] - 1]) ))
         {
            SCIPdebugMsg(scip, "Probing sdp founded conflicting bounds for integer variable %s -> cutoff!\n",
                  SCIPvarGetName(vars[newboundinds[i] - 1]));
            *result = SCIP_CUTOFF;
            break;
         }
         SCIP_CALL( SCIPchgVarUb(scip, vars[newboundinds[i] - 1], newbounds[i]) );
         *result = SCIP_REDUCEDDOM;
      }
   }

   SCIPfreeBufferArray(scip, &newboundinds);
   SCIPfreeBufferArray(scip, &newbounds);

   return SCIP_OKAY;

#else
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
#endif
}



/*
 * propagator specific interface methods
 */

/** creates the sdp-obbt propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSdpObbt(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create SdpObbt propagator data */
   propdata = NULL;
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->lastnode = -1;
   propdata->propenabled = TRUE;

   /* include propagator */
   /* use SCIPincludePropBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSdpObbt, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopySdpObbt) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSdpObbt) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSdpObbt) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolSdpObbt) );

   /* add SdpObbt propagator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propbin",
         "Should optimization-based bound tightening be performed for binary variables?",
         &propdata->propbin, TRUE, DEFAULT_PROPBIN, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propcont",
         "Should optimization-based bound tightening be performed for continuous variables?",
         &propdata->propcont, TRUE, DEFAULT_PROPCONT, NULL, NULL) );

   return SCIP_OKAY;
}
