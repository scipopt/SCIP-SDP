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

/**@file   branch_sdpobjective.c
 * @brief  highest absolute objective branching rule for SCIPSDP
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <assert.h>

#include "branch_sdpobjective.h"


#define BRANCHRULE_NAME            "sdpobjective"
#define BRANCHRULE_DESC            "branch on variable with highest absolute objective of the SDP"
#define BRANCHRULE_PRIORITY        7500000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0


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
SCIP_DECL_BRANCHEXECEXT(branchExecextSdpobjective)
{
   int i;
   int ncands;
   SCIP_VAR** cands = NULL;
   SCIP_Real* candssol; /* solution values of all candidates */
   SCIP_Real* candsscore; /* scores of all candidates */
   SCIP_VAR* maxobjvar = NULL; /* variable with currently highest absolute objective */
   SCIP_Real maxobjobj; /* objective of the current candidate with highest absolute objective */
   SCIP_Real maxobjval; /* value of the current candidate with highest absolute objective */
   SCIP_Real maxobjscore; /* score of the current candidate with highest absolute objective */

   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMessage("Executing External Branching method of SDP-objective!\n");

   /* get the external candidates, as we use the score only as a tiebreaker, we aren't interested in the number of variables of different types with maximal
    * score, so these return values are set to NULL */
   SCIPgetExternBranchCands(scip, &cands, &candssol, &candsscore, &ncands, NULL, NULL, NULL, NULL);

   assert( ncands > 0 ); /* branchExecext should only be called if the list of extern branching candidate is non-empty */

#ifdef SCIP_DEBUG
   printf("branching candidates for SDP-objective:\n");
   for (i = 0; i < ncands; i++)
      printf("%s, value = %f, objective = %f, score = %f\n", SCIPvarGetName(cands[i]), candssol[i], SCIPvarGetObj(cands[i]), candsscore[i]);
#endif

   maxobjobj = -1.0;

   /* iterate over all candidates and find the one with the highest absolute objective, use score as tiebreaker */
   for (i = 0; i < ncands; i++)
   {
      /* a candidate is better than the current one if:
       * - the absolute objective is (epsilon-)bigger than before or
       * - the absolute objective is (epsilon-)equal and the score is (epsilon-)bigger or
       * - the score is (epsilon-)equal and the absolute objective is (less than epsilon) bigger
       * - the absolute objective is (exactly) equal and the score is (less than epsilon) bigger
       */
      if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(cands[i])), maxobjobj) ||
          (SCIPisEQ(scip, REALABS(SCIPvarGetObj(cands[i])), maxobjobj) && SCIPisGT(scip, candsscore[i], maxobjscore)) ||
          (SCIPisEQ(scip, candsscore[i], maxobjscore) && SCIPvarGetObj(cands[i]) > maxobjobj) ||
          (REALABS(SCIPvarGetObj(cands[i])) == maxobjobj && candsscore[i] > maxobjscore) )
      {
         maxobjvar = cands[i];
         maxobjobj = REALABS(SCIPvarGetObj(cands[i]));
         maxobjval = candssol[i];
         maxobjscore = candsscore[i];
      }
   }

   assert( SCIPisGE(scip, maxobjobj, 0.0) );

   /* if all candidates have objective zero, we look for other variables that are coupled with the candidates and check their objective values */
   if ( SCIPisEQ(scip, maxobjobj, 0.0) )
   {
      int j;
      int c;
      int v;
      int cand;
      int candpos;
      int nvars;
      SCIP_VAR** vars;
      int nconss;
      SCIP_CONS** conss;
      int nvarsincons;
      SCIP_VAR** varsincons;
      SCIP_Bool** coupledvars; /* is there a constraint coupling candidate i and variable j ? */
      SCIP_Bool** singlecoupledvars; /* is variable j coupled with candidate i AND with no other candidate */
      int** candsincons; /* candsincons[i] gives a list of all candidates (indexed as in cands) appearing in cons i */
      int* ncandsincons; /* ncandsincons[i] gives the length of candsincons[i] */
      SCIP_Real currentobj;
      SCIP_Bool success;
      int coupledcand;

      SCIPdebugMessage("All branching candidates have objective 0.0, objective branching proceeds to check coupled variables, updated values for candidates: \n");

      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);
      assert( vars != NULL );
      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      assert( conss != NULL );

      /* allocate memory to save the coupled variables and initialize the arrays */
      SCIP_CALL( SCIPallocBufferArray(scip, &ncandsincons, nconss) );
      for (i = 0; i < nconss; i++)
         ncandsincons[i] = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &candsincons, nconss) );
      for (i = 0; i < nconss; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(candsincons[i]), ncands) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &coupledvars, ncands) );
      for (i = 0; i < ncands; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(coupledvars[i]), nvars) );
         for (j = 0; j < nvars; j++)
            coupledvars[i][j] = FALSE;
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &singlecoupledvars, ncands) );
      for (i = 0; i < ncands; i++)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(singlecoupledvars[i]), nvars) );
         for (j = 0; j < nvars; j++)
            singlecoupledvars[i][j] = FALSE;
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &varsincons, nvars) );

      /* find all variables that are coupled to a candidate */
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
            for (cand = 0; cand < ncands; cand++)
            {
               if ( varsincons[v] == cands[cand] )
               {
                  candsincons[c][ncandsincons[c]] = cand;
                  ncandsincons[c]++;
               }
            }
         }

         /* now save which variables are coupled to each candidate by adding all those that appear in this constraint to all candidates appearing in this constraint */
         for (candpos = 0; candpos < ncandsincons[c]; candpos++)
         {
            for (v = 0; v < nvarsincons; v++)
            {
               /* the coupledvars-index corresponding to a variable is its variable index - nvars, because we work on the transformed variables which
                * have indices nvars to 2*nvars - 1, as their indices start after those of the original variables */
               coupledvars[candsincons[c][candpos]][SCIPvarGetIndex(varsincons[v]) - nvars] = TRUE;
            }
         }
      }

      /* finally remove all variables, that are coupled to multiple candidates */
      for (v = 0; v < nvars; v++)
      {
         /* we use coupledcand to save if we have already found a candidate this variable is coupled with (otherwise coupledcand = -1), if we find one, we set
          * coupledcand to that index, to easily set the corresponding entry to TRUE if we don't find another candidate it is coupled with */
         coupledcand = -1;
         for (cand = 0; cand < ncands; cand++)
         {
            if ( coupledvars[cand][v] )
            {
               /* check if this is the first or the second found candidate for this variable */
               if ( coupledcand == -1 )
               {
                  /* this is the first candidate this is coupled with, so it might be the only one and we save it to potentially later set singlecoupledvars
                   * to true */
                  coupledcand = cand;
               }
               else
               {
                  /* we found a second candidate, so this variable won't be taken into account for the branching rule, so we set coupledcand to -2 to not set
                   * the corresponding entry in singlecoupledvars to TRUE and continue with the next variable */
                  coupledcand = -2;
                  break;
               }
            }
         }
         if ( coupledcand > -1 )
         {
            /* as we found exactly one candidate this variable is coupled with, we set the corresponding singlecoupledvars-entry to TRUE */
            singlecoupledvars[coupledcand][v] = TRUE;
         }
      }

      /* iterate over all variables and compute the total absolute objective of all coupled variables */
      for (cand = 0; cand < ncands; cand++)
      {
         currentobj = 0.0;
         for (v = 0; v < nvars; v++)
         {
            if (singlecoupledvars[cand][v])
               currentobj += REALABS(SCIPvarGetObj(vars[v]));
         }

         assert( SCIPisGE(scip, currentobj, 0.0) );

#ifdef SCIP_DEBUG
         printf("candidate %s, coupled with ", SCIPvarGetName(cands[cand]));
         for (v = 0; v < nvars; v++)
         {
            if (coupledvars[cand][v])
               printf("%s, ", SCIPvarGetName(vars[v]));
         }
         printf("out of those ");
         for (v = 0; v < nvars; v++)
         {
            if (singlecoupledvars[cand][v])
               printf("%s, ", SCIPvarGetName(vars[v]));
         }
         printf("are only coupled with this candidate, total objective = %f, score = %f\n", currentobj, candsscore[cand]);

#endif

         /* a candidate is better than the current one if:
          * - the total absolute objective is (epsilon-)bigger than before or
          * - the total absolute objective is (epsilon-)equal and the score is (epsilon-)bigger or
          * - the score is (epsilon-)equal and the total absolute objective is (less than epsilon) bigger
          * - the total absolute objective is (exactly) equal and the score is (less than epsilon) bigger
          */
         if ( SCIPisGT(scip, currentobj, maxobjobj) ||
             (SCIPisEQ(scip, currentobj, maxobjobj) && SCIPisGT(scip, candsscore[cand], maxobjscore)) ||
             (SCIPisEQ(scip, candsscore[cand], maxobjscore) && currentobj > maxobjobj) ||
             (currentobj == maxobjobj && candsscore[i] > maxobjscore) )
         {
            maxobjvar = cands[cand];
            maxobjobj = currentobj;
            maxobjval = candssol[cand];
            maxobjscore = candsscore[cand];
         }
      }

      /* free Memory */
      SCIPfreeBufferArray(scip, &varsincons);
      for (i = 0; i < ncands; i++)
      {
         SCIPfreeBufferArray(scip, &(singlecoupledvars[i]));
      }
      SCIPfreeBufferArray(scip, &singlecoupledvars);
      for (i = 0; i < ncands; i++)
      {
         SCIPfreeBufferArray(scip, &(coupledvars[i]));
      }
      SCIPfreeBufferArray(scip, &coupledvars);
      for (i = 0; i < nconss; i++)
      {
         SCIPfreeBufferArray(scip, &(candsincons[i]));
      }
      SCIPfreeBufferArray(scip, &candsincons);
      SCIPfreeBufferArray(scip, &ncandsincons);
   }

   /* branch */
   SCIPdebugMessage("branching on variable %s with value %f, absolute objective %f and score %f\n", SCIPvarGetName(maxobjvar), maxobjval, maxobjobj, maxobjscore);
   SCIP_CALL( SCIPbranchVarVal(scip, maxobjvar, maxobjval, NULL, NULL, NULL) );

   *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/*
 * branching rule specific interface methods
 */

/** creates the SDP highest absolute objective rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleSdpobjective(
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
   SCIP_CALL( SCIPsetBranchruleExecExt(scip, branchrule, branchExecextSdpobjective) );

   return SCIP_OKAY;
}