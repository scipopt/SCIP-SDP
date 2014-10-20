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


#define SCIP_DEBUG
#define SCIP_MORE_DEBUG /* shows number of deleted empty cols/rows and complete solution for every relaxation and variable status & bounds */

#include "relax_sdp.h"

#include <cassert>                      // for assert
#include <cstdio>                       // for NULL, printf
#include <cstring>                      // for NULL, strcmp

#include "SdpVarmapper.h"               // for SdpVarmapper
#include "sdpi/sdpi_general.h"          // for SDP-Interface
#include "scipsdp/cons_sdp.h"           // for cons_check

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
   SCIP_SDPI*            sdpi;               /* general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SdpVarmapper*         varmapper;          /* maps SCIP variables to their global SDP indices and vice versa */
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

   assert ( scip != NULL );
   assert ( sdpi != NULL );

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* prepare arrays of objective values and bounds */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &obj, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lb, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ub, nvars));

   for (i = 0; i < nvars; i++)
   {
      obj[i] = SCIPvarGetObj(vars[i]);
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
   }

   nconss = SCIPgetNConss(scip);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conss, nconss)) ;
   conss = SCIPgetConss(scip);

   /* count the number of sdpblocks and compute the number of nonzeros */
   nsdpblocks = 0;
   sdpnnonz = 0;
   sdpconstnnonz = 0;

   for (i=0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) ); // TODO: remove

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         nsdpblocks++;

         SCIP_CALL(SCIPconsSdpGetNNonz(scip, conss[i], &blocknnonz, &constnnonzcounter));
         sdpnnonz += blocknnonz;
         sdpconstnnonz += constnnonzcounter;
      }
   }

   /* create the sdp- and sdpconst-arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sdpblocksizes, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nblockvarnonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nconstblocknonz, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &col, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &row, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &val, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constcol, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constrow, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constval, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nblockvars, nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sdpvar, nsdpblocks) );

   for (i = 0; i < nsdpblocks; i++)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nblockvarnonz[i]), nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &col[i], nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &row[i], nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &val[i], nvars) );
   }

   /* get the SDP-data */
   ind = 0; /* index of the current sdp block in the complete sdp */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &blockvars, nvars));

   for (i = 0; i < nconss; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         assert( ind < nsdpblocks );

         /* allocate memory for the constant nonzeros */
         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[i], NULL, &constlength) );
         nconstblocknonz[ind] = constlength;
         if (constlength > 0)
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(constcol[ind]), constlength) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(constrow[ind]), constlength) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(constval[ind]), constlength) );
         }

         /* get the data */
         SCIPconsSdpGetData(scip, conss[i], &nblockvars[ind], &blocknnonz, &sdpblocksizes[ind], &nvars, nblockvarnonz[ind], col[ind],
            row[ind], val[ind], blockvars, &nconstblocknonz[ind], constcol[ind], constrow[ind], constval[ind]);

         /* nvars and nconstblocknonz[ind] would have been overwritten if the space in the given arrays hadn't been sufficient */
         assert ( nvars == SCIPgetNVars(scip) );
         assert ( nconstblocknonz[ind] <= constlength);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpvar[ind]), nblockvars[ind]));

         /* get global variable indices */
         for (j = 0; j < nblockvars[ind]; j++)
            sdpvar[ind][j] = SdpVarmapperGetSdpIndex(varmapper, blockvars[j]);

         ind++;
      }
   }

   /* free the memory that is no longer needed */
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);
   SCIPfreeBlockMemoryArray(scip, &blockvars, nvars);

   SCIP_CALL(SCIPsdpiLoadSDP(sdpi, nvars,  obj,  lb,  ub, nsdpblocks,
                            sdpblocksizes, nblockvars, sdpconstnnonz, nconstblocknonz, constrow,
                            constcol, constval, sdpnnonz, nblockvarnonz, sdpvar,
                            row, col,  val, 0,
                            NULL, 0, NULL, NULL, NULL)); /* insert the SDP part, add an empty LP part */

   /* free the remaining memory */


   for (i = 0; i < nsdpblocks; i++)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpvar[i]), nblockvars[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &val[i], nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &row[i], nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &col[i], nvars);
      SCIPfreeBlockMemoryArrayNull(scip, &(nblockvarnonz[i]), nvars);
      if (constlength > 0)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &(constval[i]), nconstblocknonz[ind]);
         SCIPfreeBlockMemoryArrayNull(scip, &(constrow[i]), nconstblocknonz[ind]);
         SCIPfreeBlockMemoryArrayNull(scip, &(constcol[i]), nconstblocknonz[ind]);
      }
   }

   SCIPfreeBlockMemoryArrayNull(scip, &sdpvar, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &nblockvars, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &constval, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &constrow, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &constcol, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &nconstblocknonz, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &nblockvarnonz, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &row, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &col, nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &val, sdpnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &sdpblocksizes, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &obj, nvars);
   SCIPfreeBlockMemoryArray(scip, &lb, nvars);
   SCIPfreeBlockMemoryArray(scip, &ub, nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &sdpblocksizes, nsdpblocks);

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
   int blocknnonz;
   SCIP_Real* rowvals;
   SCIP_COL** rowcols;
   SCIP_Real sciplhs;
   SCIP_Real sciprhs;
   int nrowssdpi;
   SCIP_VAR** vars;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int* inds;

   assert ( scip != NULL );
   assert ( sdpi != NULL );
   assert ( varmapper != NULL );

   nvars = SCIPgetNVars(scip); /* or get these from the varmapper, depending on how it is done in the SDP-part */

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rows, SCIPgetNLPRows(scip)) );
   SCIP_CALL_ABORT( SCIPgetLPRowsData(scip, &rows, &nrows) );

   SCIPdebugMessage("inserting %d LPRows into the interface \n", nrows);

   /* compute the total number of LP nonzeroes in SCIP */
   scipnnonz = 0;
   for (i = 0; i < nrows; i++)
      scipnnonz += SCIProwGetNNonz(rows[i]);

   /* allocate memory */
   /* the arrays need to be twice the size of those given be scipped because of the lack of left-hand-sides (which means rows could be duplicated, with
    * one constraint for the lhs and one for the rhs), these will be reallocated to the right size later */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &rhs, 2 * nrows));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &rowind, 2 * scipnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &colind, 2 * scipnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &val, 2 * scipnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &rowvals, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &rowcols, nvars));

   /* insert the nonzeroes */
   nnonz = 0; /* this is recomputed for the sdpi, because of the possible duplication of non-zeroes for lhs and rhs */
   nconss = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      blocknnonz = SCIProwGetNNonz(rows[i]);

      rowvals = SCIProwGetVals(rows[i]);
      rowcols = SCIProwGetCols(rows[i]);
      sciplhs = SCIProwGetLhs(rows[i])- SCIProwGetConstant(rows[i]);
      sciprhs = SCIProwGetRhs(rows[i])- SCIProwGetConstant(rows[i]);

      /* add row >= lhs if lhs is finite */
      if (sciplhs > -SCIPsdpiInfinity(sdpi))
      {
         for (j = 0; j < blocknnonz; j++)
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
         for (j = 0; j < blocknnonz; j++)
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

   /* these arrays are no longer needed */
   SCIPfreeBlockMemoryArray(scip, &rows, nrows);
   SCIPfreeBlockMemoryArray(scip, &rowvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &rowcols, nvars);

   /* reallocate some arrays depending on the number of lpnnonz and lhs/rhs */
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &rhs, 2 * nrows, nconss));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &rowind, 2 * scipnnonz, nnonz));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &colind, 2 * scipnnonz, nnonz));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &val, 2 * scipnnonz, nnonz));

   /* delete the old LP-block from the sdpi */
   SCIP_CALL(SCIPsdpiGetNLPRows(sdpi, &nrowssdpi));
   if (nrowssdpi > 0)
      SCIP_CALL(SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1));

   /* add the LP-block to the sdpi */
   SCIP_CALL(SCIPsdpiAddLPRows(sdpi, nconss, rhs, nnonz, (const int*)rowind, (const int*)colind, val));

   /* free the remaining arrays */
   SCIPfreeBlockMemoryArray(scip, &rhs, nconss);
   SCIPfreeBlockMemoryArray(scip, &rowind, nnonz);
   SCIPfreeBlockMemoryArray(scip, &colind, nnonz);
   SCIPfreeBlockMemoryArray(scip, &val, nnonz);

   /* update the bounds */

   /* get the variables */
   nvars = SCIPgetNVars(scip);
   assert ( nvars > 0 );
   vars = SCIPgetVars(scip);
   assert ( vars != NULL );

   /* prepare arrays of bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nvars) );

   /* get new bounds */
   for (i = 0; i < nvars; i++)
   {
      assert ( vars[i] != NULL );
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
      inds[i] = 1; /* we want to change all bounds */
   }

   /* inform interface */
   SCIP_CALL( SCIPsdpiChgBounds(sdpi, nvars, inds, lb, ub) );

   /* free the bounds-arrays */
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}

/** checks the feasibility of the problem if the solver returned some ambigous solution by calling it again with a
 *  formulation that only has the LP-part as constraints and tries to minimize the minimal eigenvalue of the SDP-constraint */
static
SCIP_RETCODE relaxIsFeasible(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP-Interface structure */
   bool&                 success,            /**< could feasibility be determined */
   bool&                 feasible            /**< whether we obtained a feasible solution */
   )
{
   SCIP_Real obj;

   assert ( sdpi != NULL );
   assert ( scip != NULL );

   /* solve with penalty without objective */
   SCIP_CALL( SCIPsdpiSolvePenalty(sdpi, 1.0, FALSE) );

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

/** calculate relaxation and process the relaxation results, may call itself recursively once to try again
 *  with a penalty formulation (or more often if the penalty formulation turns out to be unbounded) */
static
SCIP_RETCODE calc_relax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP-Interface structure */
   SdpVarmapper*         varmapper,          /**< varmapper class data */
   SCIP_Bool             withpenalty,        /**< should a penalty formulation be used */
   SCIP_Real             penaltyparam,       /**< parameter for penalty formulation, if 0 the normal SDP is solved */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound         /**< pointer to store lowerbound */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;
   int v;
   SCIP_CONS** conss;
   int nconss;

   assert ( sdpi != NULL );
   assert ( scip != NULL );
   assert ( result != NULL );
   assert ( lowerbound != NULL );
   assert ( varmapper != NULL );
   assert ( penaltyparam >= 0.0 );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars (scip);

   if ( withpenalty )
   {
      SCIP_CALL(SCIPsdpiSolvePenalty(sdpi, penaltyparam, TRUE));
   }
   else
   {
      SCIP_CALL(SCIPsdpiSolve(sdpi));
   }

#ifdef SCIP_MORE_DEBUG /* print the optimal solution */
   SCIP_Real objforscip;
   SCIP_Real* solforscip;
   SCIP_Bool allint;
   int sollength;

   SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
   sollength = nvars;
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &sollength) ); /* get both the objective and the solution from the SDP solver */

   assert ( sollength == nvars ); /* if this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                   * but this can't happen, because the array has enough space for all sdp variables */

   SCIPdebugMessage("optimal solution: objective = %f, ", objforscip);
   if( SCIPsdpiFeasibilityKnown(sdpi) )
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
      SCIPdebugMessage("y_%d = %f, ", i, solforscip[i]);
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
            SCIP_CALL( calc_relax(scip, sdpi, varmapper, TRUE, 10 * penaltyparam, result, lowerbound) );

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

         /* get solution w.r.t. SCIP variables */
         SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
         sollength = nvars;
         SCIP_CALL( SCIPsdpiGetSol(sdpi, &objforscip, solforscip, &sollength) ); /* get both the objective and the solution from the SDP solver */

         assert ( sollength == nvars ); /* if this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                         * but this can't happen, because the array has enough space for all sdp variables */

         /* check if the solution is integral */
         allint = TRUE;
         for (v = 0; v < nvars; v++)
         {
            if (SCIPvarIsIntegral(SdpVarmapperGetSCIPvar(varmapper, v)) && ! (SCIPisIntegral(scip, solforscip[v])))
            {
               allint = FALSE;
               break;
            }
         }

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         /* if called with penalty formulation check if the solution really is feasible */
         if( withpenalty )
         {
            SCIP_RESULT conefeas;

            nconss = SCIPgetNConss(scip);
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &conss, nconss));
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
            SCIPfreeBlockMemoryArray(scip, &conss, nconss);
         }

         /* this was initialized as true [and thus is always true if called without a penalty formulation], for a penalty formulation
          * the sdp-constraint was checked because this was relaxed during solving */
         if ( solisfeas )
         {
            SCIP_Bool delayed;
            SCIP_Bool cutoff_forsep;
            SCIP_Bool stored;

            *lowerbound = objforscip;

            /* test whether solution is feasible */
            SCIP_CALL( SCIPtrySol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, &stored) );

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
               SCIP_CALL( SCIPseparateSol(scip, scipsol, FALSE, FALSE, &delayed, &cutoff_forsep) );

               if ( SCIPgetNCuts(scip) > oldncuts )
                  *result = SCIP_SEPARATED;

               for (i = 0; i < nvars; ++i)
               {
                  SCIP_VAR* var = vars[i];
                  if ( SCIPvarIsIntegral(var) && ! SCIPisIntegral(scip, solforscip[i]) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                  {
                     /* @todo Check for better score or make it parameter/define */
                     SCIP_CALL( SCIPaddExternBranchCand(scip, var, 10000.0, solforscip[i]) );
                  }
               }
            }
            SCIPfreeBufferArray(scip, &solforscip);
            SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
         }
         else
         {    /* solver returned feasible solution to relaxed problem with penalty formulation, but check for psd failed */
            SCIPfreeBufferArray(scip, &solforscip);

            if (! SCIPsdpiIsGEMaxPenParam(sdpi, penaltyparam) )
            {
               /* the penalty parameter was too small to create a feasible solution */
               SCIPdebugMessage("calc_relax is called again with penaltyparameter %f because the solution of the penalty problem was infeasible in the original problem!\n", 2.0 * penaltyparam);

               /* recursive call */
               SCIP_CALL(calc_relax(scip, sdpi, varmapper, TRUE, 10.0 * penaltyparam, result, lowerbound));
            }
            else
            {   /* A penalty-only-formulation showed that the problem is feasible, but we weren't able to produce a feasible solution. */

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
      SCIPdebugMessage("calc_relax is called again with penaltyparameter %f because of non-convergence!\n", 10 * penaltyparam);
      SCIP_CALL(calc_relax(scip, sdpi, varmapper, TRUE, 10 * penaltyparam, result, lowerbound));

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
   SCIP_CALL( relaxIsFeasible(scip, sdpi, success, feasible) );

   if ( success )
   {
      if ( feasible )
      {
         /* try again with penalty formulation */
         SCIP_CALL(calc_relax(scip, sdpi, varmapper, TRUE, 1.0, result, lowerbound)); /* TODO: think about penalty parameter */
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

static
SCIP_Bool allVarsFixed(
    SCIP*                 scip                /**< SCIP data structure */
    )
{
   int i;
   SCIP_VAR** vars;

   assert ( scip != NULL );

   vars = SCIPgetVars(scip);

   /* try to find a variable that is not fixed */
   for (i = 0; i < SCIPgetNVars(scip); i++)
   {
      if (SCIPisLT(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])))
         return FALSE;
   }

   /* if no variable with lower bound strictly lower than upper bound has been found, all variables are fixed */
   return TRUE;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecSDP)
{
   SCIP_Cons** conss;
   SCIP_RELAXDATA* relaxdata;
   int nconss;
   int i;
   int nvars;
   SCIP_VAR** vars;
   int* indsforsdpi;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;

   // construct the lp and make sure, that everything is where it should be
   SCIP_Bool cutoff;

   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );
   if (cutoff)
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   // very important to call flushLP
   SCIP_CALL( SCIPflushLP(scip) );

   /* get varmapper */
   nconss = SCIPgetNConss(scip);
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &conss, nconss));
   conss = SCIPgetConss(scip);

   // it is possible to call this function for writing the problem of every node in sdpa-format to a file per node
   // SCIP_CALL(write_sdpafile(scip, problemdata, varmapper));

#ifdef SCIP_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
   for (i = 0; i < nvarsfordebug; i++)
   {
      SCIPdebugMessage("variable %d: status = %u, bounds = [%f, %f] \n", i, SCIPvarGetStatus(varsfordebug[i]), SCIPvarGetLbLocal(varsfordebug[i]), SCIPvarGetUbLocal(varsfordebug[i]));
   }
#endif

   if ( nconss == 0 )
   {
      //if there are no constraints, there is nothing to do
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if (allVarsFixed(scip))
   {
      // if all variables, really all, are fixed, I can't solve an sdp, because there is no interior point in this case, result is success and I'm separating the solution (the upper or lower bounds on a variable)
      SCIPdebugMessage("EVERYTHING IS FIXED\n");
      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_Real* ubs;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ubs, nvars));

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert ( SCIPisEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
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
      SCIPfreeBlockMemoryArray(scip, &ubs, nvars);

      return SCIP_OKAY;
   }

   relaxdata = SCIPrelaxGetData(relax);

   /* update objective and bounds and lp data */
   nvars= SCIPgetNVars(scip);
   assert ( nvars > 0 );
   vars= SCIPgetVars(scip);
   assert ( vars != NULL );

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &indsforsdpi, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &obj, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lb, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ub, nvars));

   for (i = 0; i < nvars; i++)
   {
      indsforsdpi[i] = i;
      obj[i] = SCIPvarGetObj(vars[i]);
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
   }

   SCIP_CALL(SCIPsdpiChgObj(relaxdata->sdpi, nvars, indsforsdpi, obj));
   SCIP_CALL(SCIPsdpiChgBounds(relaxdata->sdpi, nvars, indsforsdpi, lb, ub));

   SCIPfreeBlockMemoryArray(scip, &indsforsdpi, nvars);
   SCIPfreeBlockMemoryArray(scip, &obj, nvars);
   SCIPfreeBlockMemoryArray(scip, &lb, nvars);
   SCIPfreeBlockMemoryArray(scip, &ub, nvars);

   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );


   SCIP_CALL( calc_relax(scip, relaxdata->sdpi, relaxdata->varmapper, FALSE, 0.0, result, lowerbound));


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

   assert ( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   SdpVarmapperCreate(scip, &(relaxdata->varmapper), ceil(1.33 * nvars)); /* all SCIPvars will be added to this list, and 3/4 seems like a good
                                                                           * load factor (java uses this factor) */
   SdpVarmapperAddVars(scip, relaxdata->varmapper, nvars, vars);

   SCIP_CALL(putSdpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper));

   SCIPrelaxSetData(relax, relaxdata);

   return SCIP_OKAY;
}

/** copy method for sdp relaxation handler (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopySDP)
{
   SCIP_RELAXDATA* sourcedata;
   SCIP_RELAXDATA* targetdata;

   SCIP_CALL ( SCIPallocBlockMemory(scip, &targetdata) );

   sourcedata = SCIPrelaxGetData(relax);

   //TODO: copying (needs a copy method for the sdpi_general), or is this only needed for local (i.e. in every node) heuristics ?

   return SCIP_OKAY;
}


/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSDP)
{
SCIP_RELAXDATA* relaxdata;

relaxdata = SCIPrelaxGetData(relax);
assert(relaxdata != NULL);

SCIP_CALL( SdpVarmapperFree(scip, &(relaxdata->varmapper)) );
SCIPfreeMemory(scip, &relaxdata);
SCIPrelaxSetData(relax, NULL);

return SCIP_OKAY;
}

/** creates the SDP relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSDP(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;
   SCIP_SDPI* sdpi;

   assert ( scip != NULL );

   /* create SDP relaxator data */
   SCIP_CALL(SCIPallocBlockMemory(scip, &relaxdata));
   SCIP_CALL(SCIPsdpiCreate(&sdpi, NULL, SCIPblkmem(scip)));

   relaxdata->sdpi = sdpi;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecSDP, relaxdata) );
   assert( relax != NULL );

   /* include additional callbacks */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitSolSDP) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeSDP) );

   return SCIP_OKAY;
}
