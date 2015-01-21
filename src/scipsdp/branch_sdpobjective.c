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
          (SCIPisEQ(scip, candsscore[i], maxobjscore) && SCIPfrac(scip, REALABS(SCIPvarGetObj(cands[i])) > maxobjobj)) ||
          (REALABS(SCIPvarGetObj(cands[i])) == maxobjobj && candsscore[i] > maxobjscore) )
      {
         maxobjvar = cands[i];
         maxobjobj = REALABS(SCIPvarGetObj(cands[i]));
         maxobjval = candssol[i];
         maxobjscore = candsscore[i];
      }
   }

   assert( maxobjobj >= 0 );

   /* if all candidates have objective zero, we look for other variables that are coupled with the candidates and check their objective values */

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
