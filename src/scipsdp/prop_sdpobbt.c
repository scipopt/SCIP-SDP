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

/**@file   prop_SdpObbt.c
 * @brief  SdpObbt propagator
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG*/

#include <assert.h>
#include <string.h>

#include "prop_sdpobbt.h"
#include "relax_sdp.h"

/* fundamental propagator properties */
#define PROP_NAME              "obbt-sdp"
#define PROP_DESC              "optimization-based bound tightening for sdps"
#define PROP_PRIORITY                 -1100000 /**< propagator priority */
#define PROP_FREQ                     -1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_AFTERLPLOOP/**< propagation timing mask */

#define DEFAULT_PROPBIN                FALSE /**< should obbt be done for binary variables ? */
#define DEFAULT_PROPCONT               TRUE /**< should obbt be done for continuous variables ? */


/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             propbin;            /**< should obbt be done for binary variables ? */
   SCIP_Bool             propcont;           /**< should obbt be done for continuous variables ? */
   long long int         lastnode;           /**< the last node we ran for */
};


/*
 * Local methods
 */

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

   SCIPdebugMessage("create objective cutoff and add it to the LP-constraints\n");

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

   return SCIP_OKAY;
}

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
{
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


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSdpObbt)
{  /*lint --e{715}*/
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

   SCIPdebugMessage("Executing propExecSdpObbt! \n");

   /* if there is no cutoff-bound, we don't want to run */
   if ( SCIPisInfinity(scip, SCIPgetCutoffbound(scip)) )
   {
      SCIPdebugMessage("Aborting propExecSdpObbt because of lack of cutoff-bound!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /* do not run in: presolving, repropagation, probing mode, if no objective propagation is allowed  */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING || SCIPinRepropagation(scip) || SCIPinProbing(scip) || !SCIPallowObjProp(scip) )
   {
      SCIPdebugMessage("Aborting propExecSdpObbt because we are in presolving, repropagation, probing mode or no objective propagation is allowed!\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   propdata = SCIPpropGetData(prop);

   assert( propdata != NULL );

   if ( SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == propdata->lastnode )
   {
      SCIPdebugMessage("Not running again for node %lld!\n", propdata->lastnode);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   else
   {
      propdata->lastnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   }

   /* start probing */
   SCIP_CALL( SCIPstartProbing(scip) );
   SCIPdebugMessage("start probing\n");

   SCIP_CALL( addObjCutoff(scip) );

   /* make sure that we don't use the objective cutoff for the changed objective */
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/objlimit", &oldobjlimitparam) );
   SCIP_CALL( SCIPsetBoolParam(scip, "relaxing/SDP/objlimit", FALSE) );

   /* allocate memory to save bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &newbounds, 2*nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newboundinds, 2*nvars) );

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
            SCIPdebugMessage("Skipping binary variable %s\n", SCIPvarGetName(vars[v]));
         }
         else
         {
            SCIPdebugMessage("Skipping continuous variable %s\n", SCIPvarGetName(vars[v]));
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
            SCIPdebugMessage("Aborting sdp-obbt, as we were unable to solve a probing sdp!\n");
            if ( *result != SCIP_REDUCEDDOM )
               *result = SCIP_DIDNOTRUN;
            goto TERMINATE;
         }

         /* if the problem is infeasible, return with cutoff */
         if ( ! SCIPrelaxSdpIsFeasible(relaxsdp) )
         {
            SCIPdebugMessage("Probing sdp infeasible, so there can't be a better solution for this problem!\n");
            *result = SCIP_CUTOFF;
            goto TERMINATE;
         }

         /* check if we managed to tighten the bound */
         success = FALSE; /* this will be ignored, we check solvedProbing instead */
         SCIP_CALL( SCIPrelaxSdpRelaxVal(relaxsdp, &success, &probingval) );

         if ( SCIPisGT(scip, probingval, SCIPvarGetLbLocal(vars[v])) )
         {
            /* update bound TODO: check if this works or needs to be done outside of probing */
            SCIPdebugMessage("Obbt-Sdp tightened lower bound of variable %s from %f to %f !\n",
                  SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]), probingval);

            newbounds[nnewbounds] = probingval;
            newboundinds[nnewbounds] = -1 * (v+1);
            nnewbounds++;
            *result = SCIP_REDUCEDDOM;
         }
#ifdef SCIP_MORE_DEBUG
         else
         {
            SCIPdebugMessage("Obbt-Sdp found lower bound of %f for variable %s, worse than old bound %f !\n",
                  probingval, SCIPvarGetName(vars[v]), SCIPvarGetLbLocal(vars[v]));
         }
#endif
      }
#ifdef SCIP_MORE_DEBUG
      else
      {
         SCIPdebugMessage("Skipping obbt for lower bound %f of variable %s, as current relaxation's solution is tight.\n",
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
            SCIPdebugMessage("Aborting sdp-obbt, as we were unable to solve a probing sdp!\n");
            if ( *result != SCIP_REDUCEDDOM )
               *result = SCIP_DIDNOTRUN;
            goto TERMINATE;
         }

         /* if the problem is infeasible, return with cutoff */
         if ( ! SCIPrelaxSdpIsFeasible(relaxsdp) )
         {
            SCIPdebugMessage("Probing sdp infeasible, so there can't be a better solution for this problem!\n");
            *result = SCIP_CUTOFF;
            goto TERMINATE;
         }

         /* check if we managed to tighten the bound */
         success = FALSE; /* this will be ignored, we check solvedProbing instead */
         SCIP_CALL( SCIPrelaxSdpRelaxVal(relaxsdp, &success, &probingval) );

         if ( SCIPisLT(scip, -probingval, SCIPvarGetUbLocal(vars[v])) )
         {
            /* update bound TODO: check if this works or needs to be done outside of probing */
            SCIPdebugMessage("Obbt-Sdp tightened upper bound of variable %s from %f to %f !\n",
                  SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]), -probingval);

            newbounds[nnewbounds] = -probingval;
            newboundinds[nnewbounds] = v + 1;
            nnewbounds++;
            *result = SCIP_REDUCEDDOM;
         }
#ifdef SCIP_MORE_DEBUG
         else
         {
            SCIPdebugMessage("Obbt-Sdp found upper bound of %f for variable %s, worse than old bound %f !\n",
                  -probingval, SCIPvarGetName(vars[v]), SCIPvarGetUbLocal(vars[v]));
         }
#endif
      }
#ifdef SCIP_MORE_DEBUG
      else
      {
         SCIPdebugMessage("Skipping obbt for upper bound %f of variable %s, as current relaxation's solution is tight.\n",
               SCIPvarGetUbLocal(vars[v]), SCIPvarGetName(vars[v]));
      }
#endif

      /* reset the objective coefficient to zero for the next variable */
      SCIP_CALL( SCIPchgVarObjProbing(scip, vars[v], 0.0) );
   }

   TERMINATE:
   SCIP_CALL( SCIPendProbing(scip) );
   SCIPdebugMessage("end probing\n");
   SCIP_CALL( SCIPsetBoolParam(scip, "relaxing/SDP/objlimit", oldobjlimitparam) );

   for (i = 0; i < nnewbounds; i++)
   {
      if ( newboundinds[i] < 0)
      {
         SCIP_CALL( SCIPchgVarLb(scip, vars[-1 * newboundinds[i] - 1], newbounds[i]) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUb(scip, vars[newboundinds[i] - 1], newbounds[i]) );
      }
   }

   SCIPfreeBufferArray(scip, &newbounds);
   SCIPfreeBufferArray(scip, &newboundinds);

   return SCIP_OKAY;
}



/*
 * propagator specific interface methods
 */

/** creates the SdpObbt propagator and includes it in SCIP */
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

   /* add SdpObbt propagator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propbin",
         " Should obbt be done for binary variables ? ",
         &propdata->propbin, TRUE, DEFAULT_PROPBIN, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/" PROP_NAME "/propcont",
         " Should obbt be done for continuous variables ? ",
         &propdata->propcont, TRUE, DEFAULT_PROPCONT, NULL, NULL) );

   return SCIP_OKAY;
}
