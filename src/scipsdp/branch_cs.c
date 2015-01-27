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

/**@file   branch_cs.c
 * @brief  special branching rule for compressed sensing problems in SCIPSDP
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <assert.h>

#include "branch_cs.h"
#include "relax_sdp.h"
#include "nodesel_prio.h"


#define BRANCHRULE_NAME            "cs"
#define BRANCHRULE_DESC            "special branching rule for compressed sensing problems in SCIPSDP"
#define BRANCHRULE_PRIORITY        -5000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0
#define USE_NODESELPRIO

/*
 * Data structures
 */

/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of branching rule
 */

/** copy method for branchrule plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_BRANCHCOPY(branchCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz branching rule not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define branchCopyXyz NULL
#endif

/** branching execution method for external candidates */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextCs)
{
   int i;
   int ncands;
   SCIP_VAR** cands = NULL;
   SCIP_Real* candssol; /* solution values of all candidates */
   SCIP_Real* candsscore; /* scores of all candidates */
   SCIP_VAR* bestobjvar = NULL; /* variable with currently best objective */
   SCIP_Real bestobjobj; /* objective of the current candidate with best objective */
   SCIP_Real bestobjval; /* value of the current candidate with best objective */
   SCIP_Real bestobjscore; /* score of the current candidate with best objective */
   int nvars;
   SCIP_VAR** vars;
   int nbinvars; /* the number of binary variables, as vars is sorted by variable type, vars[0] - vars[nbinvars - 1] will be the binary variables */
   int nconss;
   SCIP_CONS** conss;
   int nvarsincons;
   SCIP_VAR** varsincons;
   SCIP_Bool** coupledvars; /* is there a constraint coupling binary variable i and variable j ? */
   SCIP_Bool** singlecoupledvars; /* is variable j coupled with binary variable i AND with no other binary variable */
   int** binsincons; /* binsincons[i] gives a list of all binary variables (indexed as in the variable array) appearing in cons i */
   int* nbinsincons; /* nbinsincons[i] gives the length of binsincons[i] */
   int j;
   int c;
   int v;
   int binv;
   SCIP_Bool success;
   int binpos;
   int coupledbin;
   int cand;
   SCIP_Real currentobj;
   SCIP_Real oldobj;
   SCIP_Real* oldsol;
   int neededlength;
   SCIP_Real* newsol;
   int* candtobinmapper;
   SCIP_Real diagval;
#ifdef SCIP_DEBUG
   int diagvar;
#endif
   SCIP_RELAX* relaxsdp;
   SCIP_NODE* onechild;
   SCIP_NODE* zerochild;
#ifdef USE_NODESELPRIO
   SCIP_NODESEL* nodeselprio;
#endif

   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Executing External Branching method of SDP-CS-branching-rule!\n");

   /* get the external candidates, as we use the score only as a tiebreaker, we aren't interested in the number of variables of different types with maximal
    * score, so these return values are set to NULL */
   SCIPgetExternBranchCands(scip, &cands, &candssol, &candsscore, &ncands, NULL, NULL, NULL, NULL);

   assert( ncands > 0 ); /* branchExecext should only be called if the list of extern branching candidate is non-empty */

   /* look for other variables that are coupled with the binary variables */

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);
   assert( vars != NULL );
   nbinvars = SCIPgetNBinVars( scip );
   nconss = SCIPgetNConss(scip);
   conss = SCIPgetConss(scip);
   assert( conss != NULL );

   /* allocate memory to save the coupled variables and initialize the arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &nbinsincons, nconss) );
   for (i = 0; i < nconss; i++)
      nbinsincons[i] = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &binsincons, nconss) );
   for (i = 0; i < nconss; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(binsincons[i]), nbinvars) );
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &coupledvars, nbinvars) );
   for (i = 0; i < nbinvars; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(coupledvars[i]), nvars) );
      for (j = 0; j < nvars; j++)
         coupledvars[i][j] = FALSE;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &singlecoupledvars, nbinvars) );
   for (i = 0; i < nbinvars; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &(singlecoupledvars[i]), nvars) );
      for (j = 0; j < nvars; j++)
         singlecoupledvars[i][j] = FALSE;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );

   /* find all variables that are coupled to a binary variable */
   for (c = 0; c < nconss; c++)
   {
      /* first check which candidates appear in which constraints */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nvarsincons, &success) );
      if ( ! success )
      {
         SCIPdebugMessage("couldn't get variable information from constraint %s, so ignoring it for computing coupled variables\n", SCIPconsGetName(conss[c]));
         continue; /* if we can't get the variables of this constraint, we can't include variables coupled through this constraint */
      }
      assert( nvarsincons > 0 );
      SCIPgetConsVars(scip, conss[c], varsincons, nvarsincons, &success);
      assert( success ); /* we allocated enough memory */
      assert( varsincons != NULL );
      for (v = 0; v < nvarsincons; v++)
      {
         for (binv = 0; binv < nbinvars; binv++)
         {
            if ( varsincons[v] == vars[binv] )
            {
               binsincons[c][nbinsincons[c]] = binv;
               nbinsincons[c]++;
            }
         }
      }

      /* now save which variables are coupled to each binary variable by adding all those that appear in this constraint to all binary variables appearing
       * in this constraint */
      for (binpos = 0; binpos < nbinsincons[c]; binpos++)
      {
         for (v = 0; v < nvarsincons; v++)
         {
            /* the coupledvars-index corresponding to a variable is its variable index - nvars, because we work on the transformed variables which
             * have indices nvars to 2*nvars - 1, as their indices start after those of the original variables */
            coupledvars[binsincons[c][binpos]][SCIPvarGetIndex(varsincons[v]) - nvars] = TRUE;
         }
      }
   }

   /* finally find all variables, that are only coupled to multiple binary variables */
   for (v = 0; v < nvars; v++)
   {
      /* we use coupledcand to save if we have already found a binary variable this variable is coupled with (otherwise coupledbin = -1), if we find one, we
       * set coupledbin to that index, to easily set the corresponding entry to TRUE if we don't find another binary variable it is coupled with */
      coupledbin = -1;
      for (binv = 0; binv < nbinvars; binv++)
      {
         if ( coupledvars[binv][v] )
         {
            /* check if this is the first or the second found binary variable for this variable */
            if ( coupledbin == -1 )
            {
               /* this is the first binary variable this is coupled with, so it might be the only one and we save it to potentially later set singlecoupledvars
                * to true */
               coupledbin = binv;
            }
            else
            {
               /* we found a second binary variable, so we set coupledbin to -2 to not set the corresponding entry in singlecoupledvars to TRUE and continue
                * with the next variable */
               coupledbin = -2;
               continue;
            }
         }
      }
      if ( coupledbin > -1 )
      {
         /* as we found exactly one binary variable this variable is coupled with, we set the corresponding singlecoupledvars-entry to TRUE */
         singlecoupledvars[coupledbin][v] = TRUE;
      }
   }

   /* as we will now work on the candidates, we have to map them to the binary variables to find the right entries in the arrays */
   SCIP_CALL( SCIPallocBufferArray(scip, &candtobinmapper, ncands) );
   for (cand = 0; cand < ncands; cand++)
   {
      for (binv = 0; binv < nbinvars; binv++)
      {
         if ( vars[binv] == cands[cand] )
         {
            candtobinmapper[cand] = binv;
            break;
         }
      }
      assert( candtobinmapper[cand] < nbinvars ); /* all branching candidates should be binary variables */
   }

   /* get the current solution to compute a new solution from this */
   SCIP_CALL( SCIPallocBufferArray(scip, &oldsol, nvars) );
   neededlength = nvars;
   relaxsdp = SCIPfindRelax(scip, "SDP");
   SCIP_CALL( SCIPrelaxSdpGetRelaxSol(scip, relaxsdp, &success, oldsol, &neededlength) );
   assert( success );
   assert( neededlength == nvars );
   SCIP_CALL( SCIPrelaxSdpRelaxVal(relaxsdp, &success, &oldobj) );

#ifdef SCIP_DEBUG
   printf("current solution (omitting binary variables): ");
   for (v = nbinvars; v < nvars; v++)
   {
      printf("%s = %f, ", SCIPvarGetName(vars[v]), oldsol[v]);
   }
   printf(", old objective = %f\n", oldobj);

#endif

   /* allocate memory for computing the new solution */
   SCIP_CALL( SCIPallocBufferArray(scip, &newsol, nvars) );

   /* set maxobjobj to a value that is worse than that of every possible solution */
   if ( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
      bestobjobj = SCIPinfinity(scip);
   else
      bestobjobj = -SCIPinfinity(scip);

   /* iterate over all candidates and produce an approximation for the solution after fixing this variable to zero and therefore also an approximation for the
    * objective value after this fixing */
   for (cand = 0; cand < ncands; cand++)
   {
      currentobj = 0.0;
      /* find the corresponding diagonal entry, which will be the only variable that is only coupled to this binary variable */
      for (v = 0; v < nvars; v++)
      {
         if ( singlecoupledvars[cand][v] )
         {
            diagval = oldsol[v];
#ifdef SCIP_DEBUG
            diagvar = v;
#endif
            break;
         }
      }

      /* generate the new solution, we ommitt the binary variables, as they don't influence the objective */
      for (v = nbinvars; v < nvars; v++)
      {
         if ( coupledvars[candtobinmapper[cand]][v])
         {
            /* this variable belongs to this candidates row or column, so it will be fixed to zero */
            newsol[v] = 0.0;
         }
         else
         {
            newsol[v] = oldsol[v] * (1 / (1 - diagval) );
            currentobj += SCIPvarGetObj(vars[v]) * newsol[v];
         }
      }

#ifdef SCIP_DEBUG
      printf("candidate %s, corresponding diagonal entry %s with value %f, new solution (ommitting binary variables): ",
            SCIPvarGetName(cands[cand]), SCIPvarGetName(vars[diagvar]), diagval);
      for (v = nbinvars; v < nvars; v++)
      {
         printf("%s = %f, ", SCIPvarGetName(vars[v]), newsol[v]);
      }
      printf(", new objective = %f\n", currentobj);

#endif

      /* if the new objective for current candidate is better than before (smaller for a maximization problem, bigger for minimization, as we want to reduce
       * the gap), we update the currently best candidate */
      if ( (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE && currentobj < bestobjobj) ||
           (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MINIMIZE && currentobj > bestobjobj) )
      {
         bestobjvar = cands[cand];
         bestobjobj = currentobj;
         bestobjval = candssol[cand];
         bestobjscore = candsscore[cand];
      }
   }

   /* free Memory */
   SCIPfreeBufferArray(scip, &newsol);
   SCIPfreeBufferArray(scip, &oldsol);
   SCIPfreeBufferArray(scip, &candtobinmapper);
   SCIPfreeBufferArray(scip, &varsincons);
   for (i = 0; i < nbinvars; i++)
   {
      SCIPfreeBufferArray(scip, &(singlecoupledvars[i]));
   }
   SCIPfreeBufferArray(scip, &singlecoupledvars);
   for (i = 0; i < nbinvars; i++)
   {
      SCIPfreeBufferArray(scip, &(coupledvars[i]));
   }
   SCIPfreeBufferArray(scip, &coupledvars);
   for (i = 0; i < nconss; i++)
   {
      SCIPfreeBufferArray(scip, &(binsincons[i]));
   }
   SCIPfreeBufferArray(scip, &binsincons);
   SCIPfreeBufferArray(scip, &nbinsincons);

   /* branch */
   SCIPdebugMessage("branching on variable %s with value %f, predicted objective %f and score %f\n",
                     SCIPvarGetName(bestobjvar), bestobjval, bestobjobj, bestobjscore);
   assert( ! SCIPisFeasIntegral(scip, bestobjval) );
   SCIP_CALL( SCIPbranchVarVal(scip, bestobjvar, bestobjval, &zerochild, NULL, &onechild) );

   assert( (zerochild != NULL && onechild != NULL) );

   /* set the node selection priority of the child node, for the child node with the candidate set to zero we use the computed value, for the one with candidate
    * set to one we use the objective value of the current node (as it will usually not change much after fixing a variable to one), if we are minimizing, we
    * take the negative values, as in the case a lower objective value means higher priority */
#ifdef USE_NODESELPRIO
   nodeselprio = SCIPfindNodesel(scip, "prio");
   if (SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE)
   {
      SCIP_CALL( SCIPnodeselPrioInsertNodePrio(scip, nodeselprio, onechild, oldobj) );
      SCIP_CALL( SCIPnodeselPrioInsertNodePrio(scip, nodeselprio, zerochild, bestobjobj) );
   }
   else
   {
      SCIP_CALL( SCIPnodeselPrioInsertNodePrio(scip, nodeselprio, onechild, -oldobj) );
      SCIP_CALL( SCIPnodeselPrioInsertNodePrio(scip, nodeselprio, zerochild, -bestobjobj) );
   }
#else
   SCIP_CALL( SCIPchgChildPrio(scip, zerochild, -1.0) );
   SCIP_CALL( SCIPchgChildPrio(scip, onechild, 1.0) );
#endif


   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the SDP highest absolute objective rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleCs(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create empty branching rule data */
   branchruledata = NULL;

   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY,
         BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert(branchrule != NULL);

   /* set non fundamental callbacks via setter functions */
   /* SCIP_CALL( SCIPsetBranchruleCopy(scip, branchrule, branchCopyXyz) ); */
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextCs) );

   return SCIP_OKAY;
}
