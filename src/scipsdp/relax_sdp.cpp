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

/**@file   relax_sdp.h
 * @ingroup RELAXATORS
 * @brief  SDP relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

//#define SCIP_DEBUG
//#define SCIP_MORE_DEBUG /* displays complete solution for each relaxation */
//#define SCIP_EVEN_MORE_DEBUG /* shows number of deleted empty cols/rows for every relaxation and variable status & bounds as well as all constraints in the beginning */

#include "relax_sdp.h"

#include <cassert>                      // for assert
#include <cstdio>                       // for NULL, printf
#include <cstring>                      // for NULL, strcmp

#include "SdpVarmapper.h"               // for SdpVarmapper
#include "sdpi/sdpi.h"                  // for SDP-Interface
#include "scipsdp/cons_sdp.h"           // for cons_check


#define RELAX_NAME                  "SDP"
#define RELAX_DESC                  "SDP relaxator"
#define RELAX_PRIORITY              1
#define RELAX_FREQ                  1

#define DEFAULT_SDPSOLVEREPSILON    1e-3     /**< the stopping criterion for the duality gap the sdpsolver should use */
#define DEFAULT_SDPSOLVERFEASTOL    1e-6     /**< the feasibility tolerance the SDP solver should use for the SDP constraints */


/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_SDPI*            sdpi;               /**< general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SdpVarmapper*         varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Real             ojbval;             /**< objective value of the last SDP relaxation */
   SCIP_Bool             origsolved;         /**< solved original problem to optimality (not only a penalty formulation */
   SCIP_Real             sdpsolverepsilon;   /**< the stopping criterion for the duality gap the sdpsolver should use */
   SCIP_Real             sdpsolverfeastol;   /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
   int                   sdpiterations;      /**< saves the total number of sdp-iterations */
   int                   sdpcalls;           /**< number of solved SDPs (used to compute average SDP iterations) */
};

/** inserts all the SDP data into the corresponding SDP Interface */
static
SCIP_RETCODE putSdpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper           /**< maps SCIP variables to their global SDP indices and vice versa */
   )
{
   int i;
   int j;
   int nvars;
   SCIP_VAR ** vars;
   SCIP_VAR ** blockvars;
   SCIP_CONS** conss;
   int nconss;
   int ind;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int nsdpblocks;
   int* sdpblocksizes;
   int sdpconstnnonz;
   int** constrow;
   int** constcol;
   SCIP_Real** constval;
   int sdpnnonz;
   int constnnonzcounter;
   int*** row;
   int*** col;
   SCIP_Real*** val;
   SCIP_CONSHDLR* conshdlr;
   int blocknnonz;
   int* nblockvars;
   int** nblockvarnonz;
   int* nconstblocknonz;
   int constlength;
   int** sdpvar;

   SCIP_Real param;
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDPRelax/sdpsolverepsilon", &param) );

   SCIPdebugMessage("Putting SDP Data in general interface! \n");

   assert( scip != NULL );
   assert( sdpi != NULL );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* prepare arrays of objective values and bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );

   for (i = 0; i < nvars; i++)
   {
      obj[i] = SCIPvarGetObj(vars[i]);
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
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

      const char* hdlrName = SCIPconshdlrGetName(conshdlr);

#ifdef SCIP_EVEN_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
#endif

      if ( strcmp(hdlrName, "SDP") == 0 )
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
      SCIP_CALL( SCIPallocBufferArray(scip, &(nblockvarnonz[i]), nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &col[i], nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &row[i], nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &val[i], nvars) );
   }

   /* get the SDP-data */
   ind = 0; /* index of the current sdp block in the complete sdp */
   SCIP_CALL( SCIPallocBufferArray(scip, &blockvars, nvars) );

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );

      const char* hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") == 0 )
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

         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar[ind]), nblockvars[ind]) );

         /* get global variable indices */
         for (j = 0; j < nblockvars[ind]; j++)
            sdpvar[ind][j] = SdpVarmapperGetSdpIndex(varmapper, blockvars[j]);

         ind++;
      }
   }

   /* free the memory that is no longer needed */
   SCIPfreeBufferArray(scip, &blockvars);

   SCIP_CALL(SCIPsdpiLoadSDP(sdpi, nvars,  obj,  lb,  ub, nsdpblocks,
                            sdpblocksizes, nblockvars, sdpconstnnonz, nconstblocknonz, constrow,
                            constcol, constval, sdpnnonz, nblockvarnonz, sdpvar,
                            row, col,  val, 0,
                            NULL, 0, NULL, NULL, NULL)); /* insert the SDP part, add an empty LP part */

   /* free the remaining memory */
   for (i = 0; i < nsdpblocks; i++)
   {
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

/** inserts all the LP data into the corresponding SDP Interface */
static
SCIP_RETCODE putLpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper           /**< data about fixed variables */
   )
{
   int i;
   int j;
   int nvars;
   int nconss;
   int scipnnonz;
   SCIP_Real* rhs;
   int nnonz;
   int* rowind;
   int* colind;
   SCIP_Real* val;
   SCIP_ROW** rows;
   int nrows;
   int rownnonz;
   SCIP_Real* rowvals;
   SCIP_COL** rowcols;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   int nrowssdpi;
   SCIP_VAR** vars;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* inds;

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   nvars = SCIPgetNVars(scip); /* or get these from the varmapper, depending on how it is done in the SDP-part */

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIPdebugMessage("inserting %d LPRows into the interface \n", nrows);

   /* compute the total number of LP nonzeroes in SCIP */
   scipnnonz = 0;
   for (i = 0; i < nrows; i++)
   {
      assert( rows[i] != NULL );
      scipnnonz += SCIProwGetNNonz(rows[i]);
   }

   /* allocate memory */
   /* the arrays need to be twice the size of those given be scipped because of the lack of left-hand-sides (which means rows could be duplicated, with
    * one constraint for the lhs and one for the rhs) */
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, 2 * nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowind, 2 * scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colind, 2 * scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, 2 * scipnnonz) );

   /* insert the nonzeroes */
   nnonz = 0; /* this is recomputed for the sdpi, because of the possible duplication of non-zeroes for lhs and rhs */
   nconss = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      rownnonz = SCIProwGetNNonz(rows[i]);

      rowvals = SCIProwGetVals(rows[i]);
      rowcols = SCIProwGetCols(rows[i]);
      sciplhs = SCIProwGetLhs(rows[i]) - SCIProwGetConstant(rows[i]);
      sciprhs = SCIProwGetRhs(rows[i]) - SCIProwGetConstant(rows[i]);

      /* add row >= lhs if lhs is finite */
      if ( sciplhs > -SCIPsdpiInfinity(sdpi) )
      {
         for (j = 0; j < rownnonz; j++)
         {
            colind[nnonz] = SdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = nconss;
            val[nnonz] = rowvals[j];
            nnonz++;
         }
         rhs[nconss] = sciplhs;
         nconss++;
      }

      /* add -row >= -rhs if rhs is finite */
      if (sciprhs < SCIPsdpiInfinity(sdpi))
      {
         for (j = 0; j < rownnonz; j++)
         {
            colind[nnonz] = SdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = nconss;
            val[nnonz] = -rowvals[j];
            nnonz++;
         }
         rhs[nconss] = -sciprhs;
         nconss++;
      }
   }

   /* delete the old LP-block from the sdpi */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &nrowssdpi) );
   if ( nrowssdpi > 0 )
   {
      SCIP_CALL( SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1) );
   }

   /* add the LP-block to the sdpi */
   SCIP_CALL( SCIPsdpiAddLPRows(sdpi, nconss, rhs, nnonz, (const int*)rowind, (const int*)colind, val) );

   /* free the remaining arrays */
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &colind);
   SCIPfreeBufferArray(scip, &rowind);
   SCIPfreeBufferArray(scip, &rhs);

   /* update the bounds */

   /* get the variables */
   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars(scip);
   assert( vars != NULL );

   /* prepare arrays of bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nvars) );

   /* get new bounds */
   for (i = 0; i < nvars; i++)
   {
      assert( vars[i] != NULL );
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
      inds[i] = i; /* we want to change all bounds, so all indices are included in inds */
   }

   /* inform interface */
   SCIP_CALL( SCIPsdpiChgBounds(sdpi, nvars, inds, lb, ub) );

   /* free the bounds-arrays */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}

/** checks the feasibility of the problem if the solver returned some ambiguous solution by calling it again with a
 *  formulation that only has the LP-part as constraints and tries to minimize the minimal eigenvalue of the SDP-constraint */
static
SCIP_RETCODE relaxIsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP-Interface structure */
   SCIP_RELAXDATA*       relaxdata,          /**< pointer to the data of the SDP relaxator */
   bool&                 success,            /**< could feasibility be determined */
   bool&                 feasible            /**< whether we obtained a feasible solution */
   )
{
   SCIP_Real obj;

   assert( sdpi != NULL );
   assert( scip != NULL );

   /* solve with penalty without objective */
   SCIP_CALL( SCIPsdpiSolvePenalty(sdpi, 1.0, FALSE, NULL, &(relaxdata->sdpiterations)) );

   SCIP_CALL( SCIPsdpiGetObjval(sdpi, &obj) );

   if ( SCIPsdpiIsAcceptable(sdpi) )
   {
      /* if solution is feasible and objective is <= 0, then there is a solution which is actually psd */
      if ( SCIPsdpiIsDualFeasible(sdpi) && SCIPisLE(scip, obj, 0.0) )
      {
         SCIPdebugMessage("Verified that a problem is feasible with a penalty-only-formulation.\n");
         success = true;
         feasible = true;
      }
      /* Now the objective is > 0, so it is feasible w.r.t. the LP-part, but there is no psd solution. */
      else if ( SCIPsdpiIsDualFeasible(sdpi) || SCIPsdpiIsDualInfeasible(sdpi) )
      {
         SCIPdebugMessage("Verified that a problem is infeasible with a penalty-only-formulation.\n");
         success = true;
         feasible = false;
      }
      else
      {
         SCIPdebugMessage("Even when using a penalty-only-formulation the SDP-Solver couldnot decide whether subproblem is feasible!");
         success = false;
         feasible = false;
      }
   }
   else
   {
      SCIPerrorMessage("Even when using a penalty-only-formulation the SDP-Solver didnot converge!");
      success= false;
      feasible = false;
   }

   return SCIP_OKAY;
}

/** calculate relaxation and process the relaxation results
 *
 *  May call itself recursively once to try again with a penalty formulation (or more often if the penalty formulation
 *  turns out to be unbounded).
 */
static
SCIP_RETCODE calc_relax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< data of the relaxator */
   SCIP_Bool             withpenalty,        /**< should a penalty formulation be used */
   SCIP_Real             penaltyparam,       /**< parameter for penalty formulation, if 0 the normal SDP is solved */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound          /**< pointer to store lowerbound */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;
   int v;
   SCIP_CONS** conss;
   int nconss;
   SCIP_SDPI* sdpi;
   SdpVarmapper* varmapper;

   SCIPdebugMessage("calc_relax called\n");

   assert( scip != NULL );
   assert( result != NULL );
   assert( lowerbound != NULL );
   assert( penaltyparam >= 0.0 );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars (scip);

   sdpi = relaxdata->sdpi;
   assert( sdpi != NULL );
   varmapper = relaxdata->varmapper;
   assert( varmapper != NULL );

   if ( withpenalty )
   {
      SCIP_CALL( SCIPsdpiSolvePenalty(sdpi, penaltyparam, TRUE, NULL, &(relaxdata->sdpiterations)) );
   }
   else
   {
      SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, &(relaxdata->sdpiterations)) );
      if ( SCIPsdpiIsAcceptable(sdpi) )
      {
         SCIP_CALL( SCIPsdpiGetObjval(relaxdata->sdpi, &(relaxdata->ojbval)) );
         relaxdata->origsolved = TRUE;
      }
      else
         relaxdata->origsolved = FALSE;
   }

#ifdef SCIP_MORE_DEBUG /* print the optimal solution */
   SCIP_Real objforscip;
   SCIP_Real* solforscip;
   SCIP_Bool allint;
   int sollength;

   SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
   sollength = nvars;
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &sollength) ); /* get both the objective and the solution from the SDP solver */

   assert( sollength == nvars ); /* if this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                   * but this can't happen, because the array has enough space for all sdp variables */

   SCIPdebugMessage("optimal solution: objective = %f, ", objforscip);
   if ( SCIPsdpiFeasibilityKnown(sdpi) )
   {
      SCIPdebugMessage("dual feasible: %d, ", SCIPsdpiIsDualFeasible(sdpi));
      SCIPdebugMessage("primal feasible: %d, ", SCIPsdpiIsPrimalFeasible(sdpi));
   }
   else
   {
      SCIPdebugMessage("The solver could not determine feasibility ! ");
   }
   for (i = 0; i < nvars; ++i)
   {
      SCIPinfoMessage(scip, NULL, "%s = %f, ", SCIPvarGetName(vars[i]), solforscip[i]);
   }
   SCIPdebugMessage("\n");

   SCIPfreeBufferArray(scip, &solforscip);
#endif

   if ( SCIPsdpiIsAcceptable(sdpi) && SCIPsdpiFeasibilityKnown(sdpi) )
   {
      if ( SCIPsdpiIsDualInfeasible(sdpi) )
      {
         *result = SCIP_CUTOFF;
         /* need to set lowerbound? */
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsDualUnbounded(sdpi) )
      {
         if ( (! withpenalty) || SCIPsdpiIsGEMaxPenParam(sdpi, penaltyparam) )
         {
            *result = SCIP_SUCCESS;
            *lowerbound = -SCIPinfinity(scip);
            return SCIP_OKAY;
         }
         else
         {
            /* the problem might be unbounded because of a too small penalty param */
            SCIPdebugMessage("calc_relax is called again with penaltyparameter %f because of unboundedness!\n", 10 * penaltyparam);

            /* recursive call - return result from there */
            SCIP_CALL( calc_relax(scip, relaxdata, TRUE, 10 * penaltyparam, result, lowerbound) );

            return SCIP_OKAY;
         }
      }
      else if ( SCIPsdpiIsPrimalFeasible(sdpi) && SCIPsdpiIsDualFeasible(sdpi) )
      {
#ifndef SCIP_MORE_DEBUG       /* with MORE_DEBUG these were created when accessing solution information to print it to the console */
         SCIP_Real objforscip;
         SCIP_Real* solforscip;
         SCIP_Bool allint;
#endif
         SCIP_SOL* scipsol;
         SCIP_Bool solisfeas = TRUE;
         SCIP_COL** cols;
         int ncols;
         int slength;

         /* get solution w.r.t. SCIP variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         slength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

         assert( slength == nvars ); /* if this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                      * but this can't happen, because the array has enough space for all sdp variables */

         /* check if the solution is integral */
         allint = TRUE;
         for (v = 0; v < nvars; v++)
         {
            if ( SCIPvarIsIntegral(SdpVarmapperGetSCIPvar(varmapper, v)) && ! (SCIPisIntegral(scip, solforscip[v])) )
            {
               allint = FALSE;
               break;
            }
         }

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         /* if called with penalty formulation check if the solution really is feasible */
         if ( withpenalty )
         {
            SCIP_RESULT conefeas;

            nconss = SCIPgetNConss(scip);
            conss = SCIPgetConss(scip);

            for (i = 0; i < nconss; ++i)
            {
               SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], scipsol, FALSE, TRUE, FALSE, &conefeas) );
               if ( conefeas == SCIP_INFEASIBLE )
               {
                  solisfeas = FALSE;
                  break;
               }
            }
         }

         /* this was initialized as true [and thus is always true if called without a penalty formulation], for a penalty formulation
          * the sdp-constraint was checked because this was relaxed during solving */
         if ( solisfeas )
         {
            SCIP_Bool stored;
            SCIP_Bool allfeas;

            *lowerbound = objforscip;

            if ( allint ) /* if the solution is integer, we might have found a new best solution for the MISDP */
            {
               SCIP_CALL( SCIPcheckSol(scip, scipsol, TRUE, FALSE, FALSE, FALSE, &allfeas) ); /* is this really needed ? */
               if ( allfeas )
               {
                  SCIP_CALL( SCIPtrySol(scip, scipsol, TRUE, FALSE, FALSE, FALSE, &stored) );
                  if (stored)
                     SCIPdebugMessage("feasible solution for MISDP found, cut node off, solution is new best solution \n");
                  else
                     SCIPdebugMessage("feasible solution for MISDP found, cut node off, solution is worse than earlier one \n");

                  SCIPfreeBufferArray(scip, &solforscip);
                  SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               SCIPdebugMessage("WARNING!!! Found a solution that is feasible for SDP and integrality, but infeasible for SCIP, this will probably not properly get enforced ! \n");
            }

            /* copy solution */
            SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
            for (i = 0; i < ncols; i++)
            {
               SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );
            }

            SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
            *result = SCIP_SUCCESS;

            /* if all int and binary vars are integral, nothing else needs to be done */
            if ( ! allint )
            {
               int oldncuts = SCIPgetNCuts(scip);
               /* ????????????? Should this be called from relaxator ?? */
               //SCIP_CALL( SCIPseparateSol(scip, scipsol, FALSE, FALSE, &delayed, &cutoff_forsep) );

               if ( SCIPgetNCuts(scip) > oldncuts )
                  *result = SCIP_SEPARATED;

               for (i = 0; i < nvars; ++i)
               {
                  SCIP_VAR* var = vars[i];
                  if ( SCIPvarIsIntegral(var) && ! SCIPisIntegral(scip, solforscip[i]) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                  {
                     SCIP_Real frac;

                     /* use current integer-infeasibility as score */
                     frac = SCIPfrac(scip, solforscip[i]);
                     // SCIP_Real inf = frac < 0.5 ? frac : 1 - frac;

                     SCIP_CALL( SCIPaddExternBranchCand(scip, var, frac, solforscip[i]) );
                  }
               }
            }
            SCIPfreeBufferArray(scip, &solforscip);
            SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
         }
         else
         {    /* solver returned feasible solution to relaxed problem with penalty formulation, but check for psd failed */
            SCIPfreeBufferArray(scip, &solforscip);
            if ( ! SCIPsdpiIsGEMaxPenParam(sdpi, penaltyparam) )
            {
               /* the penalty parameter was too small to create a feasible solution */
               SCIPdebugMessage("calc_relax is called again with penaltyparameter %f because the solution of the penalty problem was infeasible in the original problem!\n", 10.0 * penaltyparam);

               /* recursive call */
               SCIP_CALL( calc_relax(scip, relaxdata, TRUE, 10.0 * penaltyparam, result, lowerbound) );
            }
            else
            {
               /* A penalty-only-formulation showed that the problem is feasible, but we weren't able to produce a feasible solution. */
               /* We try to reuse the relaxation result of the parent node (if one exists): */
               SCIP_NODE* node = SCIPnodeGetParent(SCIPgetCurrentNode(scip));
               if ( node == 0 )
               {
                  *result = SCIP_SUSPENDED;
                  SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, as there is no parent node the relaxation is suspended. \n");
                  return SCIP_OKAY;
               }

               *lowerbound = SCIPnodeGetLowerbound(node);
               *result = SCIP_SUCCESS;
               SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, so the relaxation result from the parent node was copied. \n");
               SCIP_CALL( SCIPupdateLocalLowerbound(scip, *lowerbound) );
            }
         }
      }

      return SCIP_OKAY;
   }

   /* the solver either didnot converge or couldnot determine whether the problem is feasible */
   if ( withpenalty && (! SCIPsdpiIsGEMaxPenParam(sdpi, penaltyparam)) )
   {
      /* the penalty parameter was too small to make the SDP solver more stable */
      SCIPdebugMessage("calc_relax is called again with penaltyparameter %f because of non-convergence!\n", 10.0 * penaltyparam);
      SCIP_CALL( calc_relax(scip, relaxdata, TRUE, 10.0 * penaltyparam, result, lowerbound) );

      return SCIP_OKAY;
   }
   else if ( SCIPsdpiIsGEMaxPenParam(sdpi, penaltyparam) )
   {
      /* A penalty-only-formulation showed that the problem is feasible, but no feasible solution could be produced,
       * so we reuse the relaxation result of the parent node (if one exists) */
      SCIP_NODE* node = SCIPnodeGetParent(SCIPgetCurrentNode(scip));
      if ( node == 0 )
      {
         *result = SCIP_SUSPENDED;
         SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, as there is no parent node the relaxation is suspended. \n");
         return SCIP_OKAY;
      }

      *lowerbound = SCIPnodeGetLowerbound(node);
      *result = SCIP_SUCCESS;
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, *lowerbound) );
      SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, so the relaxation result from the parent node was copied. \n");
      return SCIP_OKAY;
   }

   /* because of earlier ifs and returns this is only done if penaltyparam == 0 and the solver didnont converge */

   /* check for feasibility via penalty-only-formulation */
   bool success;
   bool feasible;
   SCIP_CALL( relaxIsFeasible(scip, sdpi, relaxdata, success, feasible) );

   if ( success )
   {
      if ( feasible )
      {
         /* try again with penalty formulation */
         SCIP_CALL( calc_relax(scip, relaxdata, TRUE, 1.0, result, lowerbound) ); /* TODO: think about penalty parameter */
      }
      else
      {
         /* penalty-only-formulation showed, that the problem is infeasible */
         *result = SCIP_CUTOFF;
      }
   }
   else
   {
      /* even with penalty-only-formulation the solver didnot converge or couldn't determine feasibility,
      * so we reuse the relaxation result of the parent node (if one exists) */
      SCIP_NODE* node = SCIPnodeGetParent(SCIPgetCurrentNode(scip));
      if ( node == 0 )
      {
         *result = SCIP_SUSPENDED;
         SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, as there is no parent node the relaxation is suspended. \n");
         return SCIP_OKAY;
      }

      *lowerbound = SCIPnodeGetLowerbound(node);
      *result = SCIP_SUCCESS;
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, *lowerbound) );
      SCIPdebugMessage("The problem was shown to be feasible by a penalty formulation, but no solution was found, so the relaxation result from the parent node was copied. \n");
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
   int i;
   SCIP_VAR** vars;

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
SCIP_DECL_RELAXEXEC(relaxExecSDP)
{
   SCIP_RELAXDATA* relaxdata;
   int nconss;
   int i;
   int nvars;
   SCIP_VAR** vars;

   // construct the lp and make sure, that everything is where it should be
   SCIP_Bool cutoff;
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   if ( cutoff )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   // very important to call flushLP
   SCIP_CALL( SCIPflushLP(scip) );

   /* get varmapper */
   nconss = SCIPgetNConss(scip);

   // it is possible to call this function for writing the problem of every node in sdpa-format to a file per node
   // SCIP_CALL(write_sdpafile(scip, problemdata, varmapper));

#ifdef SCIP_EVEN_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
   for (i = 0; i < nvarsfordebug; i++)
   {
      SCIPdebugMessage("variable %s: status = %u, integral = %u, bounds = [%f, %f] \n", SCIPvarGetName(varsfordebug[i]), SCIPvarGetStatus(varsfordebug[i]),
         SCIPvarIsIntegral(varsfordebug[i]), SCIPvarGetLbLocal(varsfordebug[i]), SCIPvarGetUbLocal(varsfordebug[i]));
   }
#endif

   if ( nconss == 0 )
   {
      //if there are no constraints, there is nothing to do
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( allVarsFixed(scip) )
   {
      SCIP_Bool feasible;

      // if all variables, really all, are fixed, I can't solve an sdp, because there is no interior point in this case,
      // result is success and I'm separating the solution (the upper or lower bounds on a variable)
      SCIPdebugMessage("EVERYTHING IS FIXED\n");
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_Real* ubs;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ubs, nvars) );

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert( SCIPisEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
      }
      int sense = SCIPgetObjsense(scip);
      if ( sense == -1 )
         *lowerbound *= -1;

      SCIP_SOL* scipsol;
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_Bool stored ;
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      /* check if the solution really is feasible */
      SCIP_CALL( SCIPcheckSol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, &feasible) );

      if ( feasible )
         SCIP_CALL( SCIPtrySolFree(scip, &scipsol, FALSE, FALSE, FALSE, FALSE, &stored) );

      if (feasible && stored == 1)
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_CUTOFF;

      SCIPfreeBlockMemoryArray(scip, &ubs, nvars);

      return SCIP_OKAY;
   }

   relaxdata = SCIPrelaxGetData(relax);

   /* update LP Data in Interface */
   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );

   SCIP_CALL( calc_relax(scip, relaxdata, FALSE, 0.0, result, lowerbound));
   relaxdata->sdpcalls++;

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** this method is called after presolving is finished, at this point the varmapper is prepared and the SDP Interface is initialized and gets
 *  the SDP information from the constraints */
static
SCIP_DECL_RELAXINIT(relaxInitSolSDP)
{
   SCIP_RELAXDATA* relaxdata;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real epsilon;
   SCIP_Real feastol;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   relaxdata->ojbval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpiterations = 0;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SdpVarmapperCreate(scip, &(relaxdata->varmapper), ceil(1.33 * nvars)) ); /* all SCIPvars will be added to this list, and 3/4 seems
                                                                                        * like a good load factor (java uses this factor) */
   SCIP_CALL( SdpVarmapperAddVars(scip, relaxdata->varmapper, nvars, vars) );

   if ( SCIPgetNVars(scip) > 0 )
   {
      SCIP_CALL( putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );
   }

   /* set the parameters of the SDP-Solver */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDPRelax/sdpsolverepsilon", &epsilon) );
   SCIP_CALL( SCIPsdpiSetEpsilon(relaxdata->sdpi, epsilon) );

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDPRelax/sdpsolverfeastol", &feastol) );
   SCIP_CALL( SCIPsdpiSetFeastol(relaxdata->sdpi, feastol) );

   return SCIP_OKAY;
}

/*
* copy method for sdp relaxation handler (called when SCIP copies plugins)
static
SCIP_DECL_RELAXCOPY(relaxCopySDP)
{
   SCIP_RELAXDATA* sourcedata;
   SCIP_RELAXDATA* targetdata;

   SCIP_CALL ( SCIPallocBlockMemory(scip, &targetdata) );

   sourcedata = SCIPrelaxGetData(relax);

   //TODO: copying (needs a copy method for the sdpi_general), or is this only needed for local (i.e. in every node) heuristics ?

   return SCIP_OKAY;
}*/

/** free the relaxator's data */
static
SCIP_DECL_RELAXEXIT(relaxExitSDP)
{
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );

   SCIPdebugMessage("Exiting Relaxation Handler \n");

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   SCIP_CALL( SdpVarmapperFree(scip, &(relaxdata->varmapper)) );
   SCIP_CALL( SCIPsdpiFree(&(relaxdata->sdpi)) );
   SCIPfreeBlockMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}

/** creates the SDP relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSDP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax;
   SCIP_SDPI* sdpi;

   assert( scip != NULL );

   /* create SDP relaxator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &relaxdata) );
   SCIP_CALL( SCIPsdpiCreate(&sdpi, NULL, SCIPblkmem(scip)) );

   relaxdata->sdpi = sdpi;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecSDP, relaxdata) );
   assert( relax != NULL );

   /* include additional callbacks */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitSolSDP) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitSDP) );

   /* add parameters for SDP-solver */
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDPRelax/sdpsolverepsilon", "the stopping criterion for the duality gap the sdpsolver should use",
         &(relaxdata->sdpsolverepsilon), TRUE, DEFAULT_SDPSOLVEREPSILON, 1e-20, 0.001, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDPRelax/sdpsolverfeastol", "the feasibility tolerance the SDP solver should use for the SDP constraints",
         &(relaxdata->sdpsolverfeastol), TRUE, DEFAULT_SDPSOLVERFEASTOL, 1e-17, 0.001, NULL, NULL) );

   /* add description of SDP-solver */
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPsdpiGetSolverName(), SCIPsdpiGetSolverDesc()) );

   return SCIP_OKAY;
}


/* external functions */

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown. Note: if a variable is either fixed or unbounded in the dual
 *  problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_RETCODE SCIPrelaxGetPrimalBoundVars(
   SCIP_RELAX*           relax,              /**< SDP relaxator to information for */
   SCIP_Real*            lbvars,             /**< returns the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< returns the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars
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

/** returns optimal objective value of the current SDP relaxation, if the last SDP relaxation was successfully solved*/
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP relaxator to get objective value for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP relaxation solved successfully */
   SCIP_Real*            objval              /**< returns the optimal objective value of the SDP relaxation */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( objval != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;
   *objval = relaxdata->ojbval;

   return SCIP_OKAY;
}

/** returns total number of SDP iterations */
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*            relax               /**< SDP relaxator to get the iterations for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->sdpiterations;
}

/** returns number of solved SDP relaxations */
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*            relax               /**< SDP relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->sdpcalls );
}
