/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_sdp.h
 * @ingroup RELAXATORS
 * @brief  SDP relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

#include "relax_sdp.h"

#include <cassert>                      // for assert
#include <cstdio>                       // for NULL, printf
#include "SdpInterface.h"               // for SdpInterface
#include "SdpSolverFactory.h"
#include "SdpProblem.h"                 // for SdpProblem
#include "SdpVarMapper.h"               // for SdpVarMapper



#define RELAX_NAME             "SDP"
#define RELAX_DESC             "SDP relaxator"
#define RELAX_PRIORITY         1
#define RELAX_FREQ             1

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
};



/** check if variable bounds are fulfilled */
static
SCIP_RETCODE check_bounds(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            var,                /**< variables */
   SCIP_SOL*             scipsol,            /**< solution of scip-type to check bounds for */
   bool*                 sol_is_feas         /**< pointer to store if solution is feasible */
   )
{
   *sol_is_feas = TRUE;
   for (int v = 0; v < nvars; ++v )
   {
      const SCIP_Real solval = SCIPgetSolVal(scip, scipsol, var[v]);

      if( solval != SCIP_UNKNOWN )
      {
         const SCIP_Real lb = SCIPvarGetLbGlobal(var[v]);
         const SCIP_Real ub = SCIPvarGetUbGlobal(var[v]);
         *sol_is_feas = *sol_is_feas && SCIPisFeasGE(scip, solval, lb) && SCIPisFeasLE(scip, solval, ub);
      }
   }
   return SCIP_OKAY;
}



/** call solve */
static
SCIP_RETCODE calc_relax(
   SdpInterface*         sdpsolver,          /**< sdpsolver class object */
   SdpProblem*           problemdata,        /**< data structure with problem-data of a specific node */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound,         /**< pointer to store lowerbound */
   SdpVarMapper*         varmapper           /**< varmapper class data */
   )
{
   SCIP_VAR** vars;
   vars = SCIPgetVars (scip);
   const int nvars = SCIPgetNVars(scip);
   char status;
   int solutiontype;
   double *sol_for_scip;

   SCIP_CALL( SCIPallocBufferArray(scip, &sol_for_scip, nvars));
   SCIP_CALL( sdpsolver->sdp_solve(varmapper, &status, &solutiontype, sol_for_scip) );

   double obj_val = 0;
   for (int i = 0; i < nvars; i++)
   {
      if (sol_for_scip[i] != 0)
      {
         obj_val += sol_for_scip[i] * SCIPvarGetObj(vars[i]);
      }

   }

   bool go_on_and_try = FALSE;
   assert(sol_for_scip != NULL);

   if ((status == 'c') && (solutiontype == 3) )
   {
      if (! varmapper->get_intsfixed())
      {  //we don't know anything about the problem status, so we solve again with a penalty formulation, that is able to find a feasible point if there is one
         SCIP_CALL( sdpsolver->again_with_penalty(result, lowerbound, varmapper, problemdata) );

         SCIPfreeBufferArray(scip, &sol_for_scip);
         return SCIP_OKAY;
      }
      else
      {
         //if all ints are fixed, we want to check feasbility and maybe use the solution
         go_on_and_try = TRUE;
      }

   }

   if ( (status == 's') && (solutiontype == 2) )
   {
      *result = SCIP_SUCCESS;
      *lowerbound = SCIPinfinity(scip);
      SCIPfreeBufferArray(scip, &sol_for_scip);
      return SCIP_OKAY;
   }

   if ( (status == 's') && (solutiontype == 1) )
   {
      *result = SCIP_CUTOFF;
      SCIPfreeBufferArray(scip, &sol_for_scip);
      return SCIP_OKAY;
   }

   if ((solutiontype == 0) || ((status == 's') && (solutiontype == 3)) || (go_on_and_try) )
   {
      bool sol_is_feas = TRUE;
      //create a scip_sol object and set the values
      SCIP_SOL* scipsol;
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_Bool stored;
      SCIP_CALL( SCIPsetSolVals (scip, scipsol, nvars, vars, sol_for_scip) );

      SCIP_CALL(check_bounds(scip, nvars, vars, scipsol, &sol_is_feas));

      if (sol_is_feas == 1)
      {
         *lowerbound = obj_val;


         SCIP_Bool  delayed;
         SCIP_Bool  cutoff_forsep;

         SCIP_CALL( SCIPseparateSol (scip, scipsol, FALSE, FALSE, &delayed, &cutoff_forsep) );
         SCIP_CALL( SCIPtrySol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, &stored) );

         SCIP_COL** cols;
         int ncols;
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

         for (int i = 0; i < ncols; i++)
         {
            SCIP_CALL(SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))));
         }

         SCIP_CALL(SCIPmarkRelaxSolValid(scip));
         SCIP_CALL(SCIPfreeSol(scip, &scipsol));
      }

      int ncuts = SCIPgetNCuts(scip);


      if ( (status == 's') && (solutiontype == 0) )
      {
         if (sol_is_feas == 0)
         {

            SCIP_CALL( sdpsolver->again_with_penalty(result, lowerbound, varmapper, problemdata) );

            SCIPfreeBufferArray(scip, &sol_for_scip);
            return SCIP_OKAY;
         }
      }

      if (SCIPgetNIntVars(scip) == 0 && SCIPgetNBinVars(scip) == 0) //there are no integer and binary vars, we are done
      {
         *result = SCIP_SUCCESS;
      }
      else
      {
         if (ncuts == 0)
         {
            *result = SCIP_SUCCESS;
         }
         else
         {
            *result = SCIP_SEPARATED;
         }
      }

      for (int i = 0; i < nvars; i++)
      {
         if (!SCIPisIntegral(scip, sol_for_scip[i]) && SCIPvarIsIntegral(vars[i]) && (SCIPvarGetLbLocal(vars[i]) != SCIPvarGetUbLocal(vars[i])))
            SCIP_CALL( SCIPaddExternBranchCand(scip, vars[i], 10000, sol_for_scip[i]) );
      }

      SCIPfreeBufferArray(scip, &sol_for_scip);

      return SCIP_OKAY;
   }
   if ( (status == 'c') && (solutiontype == 1 || solutiontype == 2))
   {
      //DSDP is not converged and the solution it produced was not feasible.


      if (varmapper->get_intsfixed())
      {
         *result = SCIP_CUTOFF;
         SCIPfreeBufferArray(scip, &sol_for_scip);
         return SCIP_OKAY;
      }
      else
      {
         *result = SCIP_DIDNOTRUN;
         //we are not converged, so we can't trust the infeasible or unbounded result, therefore we return didnotrun

         SCIPfreeBufferArray(scip, &sol_for_scip);

         return SCIP_OKAY;
      }
   }
   return SCIP_ERROR; //Cannot be reached
}

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecSDP)
{
   // construct the lp and make sure, that everything is where it should be
   SCIP_Bool cutoff;
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   if (cutoff)
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   // very important to call flusLP
   SCIP_CALL( SCIPflushLP(scip) );

   SdpVarMapper* varmapper;
   varmapper = new SdpVarMapper(scip);
   varmapper->init();

   SdpProblem* problemdata;
   problemdata = new SdpProblem(scip, varmapper);

   // it is possible to call this function for writing the problem of every node in sdpa-format to a file per node
   // SCIP_CALL(write_sdpafile(scip, problemdata, varmapper));

   int nlprows;
   nlprows = SCIPgetNCuts(scip) + SCIPgetNLPRows(scip);

   if ( (nlprows == 0) && (problemdata->get_nsdpcones() == 0) )
   {
      //if there are no constraints, there is nothing to do
      *result = SCIP_DIDNOTRUN;
      SCIP_CALL(varmapper->exit());
      delete varmapper;
      delete problemdata;
      return SCIP_OKAY;
   }

   SdpInterface* sdpsolver;
   char* value;
   SCIP_CALL(SCIPgetStringParam(scip, "sdpsolver", &value));
   sdpsolver = SdpSolverFactory::createSdpSolver(scip, "dsdp");

   SCIP_CALL( sdpsolver->put_data_in(problemdata, varmapper) );

   if (varmapper->get_allfixed())
   {
      // if all variables, really all, are fixed, I can't solve an sdp, because there is no interior point in this case, result is success and I'm separating the solution (the upper or lower bounds on a variable
      SCIPdebugMessage("EVERYTHING IS FIXED\n");
      SCIP_VAR** vars = SCIPgetVars(scip);
      const int nvars = SCIPgetNVars(scip);

      SCIP_Real* ubs;
      SCIP_CALL(SCIPallocBufferArray(scip, &ubs, nvars));

      *lowerbound = 0.0;
      for (int i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
      }
      int sense = SCIPgetObjsense(scip);
      if (sense == -1)
      {
         *lowerbound *= -1;
      }

      SCIP_SOL* scipsol;
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_Bool stored ;
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      SCIP_CALL( SCIPtrySolFree(scip, &scipsol, FALSE, TRUE, TRUE, TRUE, &stored) );

      if (stored == 1)
      {
         *result = SCIP_SUCCESS;
      }
      else
      {
         *result = SCIP_CUTOFF;
      }
      SCIPfreeBufferArray(scip, &ubs);
      delete sdpsolver;
      delete problemdata;
      SCIP_CALL(varmapper->exit());
      delete varmapper;

      return SCIP_OKAY;
   }

   SCIP_CALL( calc_relax(sdpsolver, problemdata, scip, result, lowerbound, varmapper));

   SCIP_CALL(varmapper->exit());
   delete varmapper;
   delete sdpsolver;
   delete problemdata;

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** creates the SDP relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSDP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax = NULL;

   /* create SDP relaxator data */

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecSDP, relaxdata) );
   assert( relax != NULL );

   /* add xyz relaxator parameters */

   return SCIP_OKAY;
}
