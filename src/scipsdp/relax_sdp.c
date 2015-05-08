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

/**@file   relax_sdp.c
 * @ingroup RELAXATORS
 * @brief  SDP relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG   *//* displays complete solution for each relaxation */
/* #define SCIP_EVEN_MORE_DEBUG  *//* shows number of deleted empty cols/rows for every relaxation and variable status & bounds as well as all constraints in the beginning */

#include "relax_sdp.h"

#include "assert.h"                     /* for assert */
#include "string.h"                     /* for strcmp */

#include "SdpVarmapper.h"
#include "sdpi/sdpi.h"
#include "scipsdp/cons_sdp.h"


#define RELAX_NAME                  "SDP"
#define RELAX_DESC                  "SDP relaxator"
#define RELAX_PRIORITY              1
#define RELAX_FREQ                  1

#define DEFAULT_SDPSOLVEREPSILON    1e-5     /**< the stopping criterion for the duality gap the sdpsolver should use */
#define DEFAULT_SDPSOLVERFEASTOL    0.9e-4     /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
#define DEFAULT_THREADS             1        /**< number of threads used for SDP solving */
#define DEFAULT_OBJLIMIT            FALSE    /**< should an objective limit be given to the SDP-Solver ? */

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_SDPI*            sdpi;               /**< general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SdpVarmapper*         varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Real             objval;             /**< objective value of the last SDP relaxation */
   SCIP_Bool             origsolved;         /**< solved original problem to optimality (not only a penalty formulation */
   SCIP_Real             sdpsolverepsilon;   /**< the stopping criterion for the duality gap the sdpsolver should use */
   SCIP_Real             sdpsolverfeastol;   /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
   int                   sdpiterations;      /**< saves the total number of sdp-iterations */
   int                   threads;            /**< number of threads used for SDP solving */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             objlimit;           /**< Should an objective limit be given to the SDP solver? */
   int                   sdpcalls;           /**< number of solved SDPs (used to compute average SDP iterations) */
   long int              lastsdpnode;        /**< number of the SCIP node the current SDP-solution belongs to */
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
   const char* hdlrName;

   SCIP_Real param;
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolverepsilon", &param) );

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

      hdlrName = SCIPconshdlrGetName(conshdlr);

#ifdef SCIP_EVEN_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
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

      hdlrName = SCIPconshdlrGetName(conshdlr);

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
                            NULL, NULL, 0, NULL, NULL, NULL)); /* insert the SDP part, add an empty LP part */

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
   SCIP_Real* lhs;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &lhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rhs, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rowind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colind, scipnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &val, scipnnonz) );

   /* insert the nonzeroes */
   nnonz = 0; /* this is recomputed for the sdpi, because of the possible duplication of non-zeroes for lhs and rhs */
   nconss = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      SCIP_ROW* row;

      row = rows[i];
      assert( row != 0 );
      rownnonz = SCIProwGetNNonz(row);

      rowvals = SCIProwGetVals(row);
      rowcols = SCIProwGetCols(row);
      sciplhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      sciprhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

      for (j = 0; j < rownnonz; j++)
      {
         assert( SCIPcolGetVar(rowcols[j]) != 0 );
         colind[nnonz] = SdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
         rowind[nnonz] = nconss;
         val[nnonz] = rowvals[j];
         nnonz++;
      }
      lhs[nconss] = sciplhs;
      rhs[nconss] = sciprhs;
      nconss++;
   }

   /* delete the old LP-block from the sdpi */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &nrowssdpi) );
   if ( nrowssdpi > 0 )
   {
      SCIP_CALL( SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1) );
   }

   /* add the LP-block to the sdpi */
   SCIP_CALL( SCIPsdpiAddLPRows(sdpi, nconss, lhs, rhs, nnonz, (const int*)rowind, (const int*)colind, val) );

   /* free the remaining arrays */
   SCIPfreeBufferArray(scip, &val);
   SCIPfreeBufferArray(scip, &colind);
   SCIPfreeBufferArray(scip, &rowind);
   SCIPfreeBufferArray(scip, &rhs);
   SCIPfreeBufferArray(scip, &lhs);

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

/** calculate relaxation and process the relaxation results
 *
 *  May call itself recursively once to try again with a penalty formulation (or more often if the penalty formulation
 *  turns out to be unbounded).
 */
static
SCIP_RETCODE calc_relax(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAXDATA*       relaxdata,          /**< data of the relaxator */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound          /**< pointer to store lowerbound */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;
   int v;
   SCIP_SDPI* sdpi;
   SdpVarmapper* varmapper;

   SCIPdebugMessage("calc_relax called\n");

   assert( scip != NULL );
   assert( result != NULL );
   assert( lowerbound != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );
   vars = SCIPgetVars (scip);

   sdpi = relaxdata->sdpi;
   assert( sdpi != NULL );
   varmapper = relaxdata->varmapper;
   assert( varmapper != NULL );

   SCIP_CALL( SCIPsdpiSetIntpar(sdpi, SCIP_SDPPAR_SDPINFO, relaxdata->sdpinfo) );

   if ( relaxdata->objlimit )
   {
      /* set the objective limit */
      assert( SCIPgetUpperbound(scip) > -SCIPsdpiInfinity(sdpi) );
      SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPgetUpperbound(scip)) );
   }

   /* solve the problem */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, &(relaxdata->sdpiterations)) );
   relaxdata->lastsdpnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   if ( SCIPsdpiWasSolved(sdpi) && SCIPsdpiSolvedOrig(sdpi) )
      relaxdata->origsolved = TRUE;
   else if ( ! SCIPsdpiWasSolved(sdpi) )
   {
      /* We couldn't solve the problem, not even with a penalty formulation, so we reuse the relaxation result of the parent node (if one exists) */
      SCIP_NODE* node = SCIPnodeGetParent(SCIPgetCurrentNode(scip));
      if ( node == 0 )
      {
         /* TODO: if we could generate a feasible solution via penalty-only and a lower bound via penalty-with-objective, we could use those two together here */
         *result = SCIP_SUSPENDED;
         SCIPdebugMessage("The relaxation of the root node could not be solved, as there is no parent node the relaxation is suspended. \n");
         return SCIP_OKAY;
      }

      *lowerbound = SCIPnodeGetLowerbound(node);
      *result = SCIP_SUCCESS;
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, *lowerbound) );
      SCIPdebugMessage("The relaxation couldn't be solved, so the relaxation result from the parent node was copied. \n");
      return SCIP_OKAY;
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
      printf("%s = %f, ", SCIPvarGetName(vars[i]), solforscip[i]);
   }
   SCIPdebugMessage("\n");

   SCIPfreeBufferArray(scip, &solforscip);
#endif

   if ( SCIPsdpiIsAcceptable(sdpi) )
   {
      if ( SCIPsdpiIsDualInfeasible(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to infeasibility.");
         *result = SCIP_CUTOFF;
         /* TODO: need to set lowerbound? */
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsObjlimExc(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to objective limit.");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsDualUnbounded(sdpi) )
      {
         SCIPdebugMessage("Node unbounded.");
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsPrimalFeasible(sdpi) && SCIPsdpiIsDualFeasible(sdpi) )
      {
#ifndef SCIP_MORE_DEBUG       /* with MORE_DEBUG these were created when accessing solution information to print it to the console */
         SCIP_Real objforscip;
         SCIP_Real* solforscip;
         SCIP_Bool allint;
#endif
         SCIP_SOL* scipsol;
         SCIP_COL** cols;
         int ncols;
         int slength;
         SCIP_Bool stored;
         SCIP_Bool allfeas;

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
            if ( SCIPvarIsIntegral(SdpVarmapperGetSCIPvar(varmapper, v)) && ! (SCIPisFeasIntegral(scip, solforscip[v])) )
            {
               allint = FALSE;
               break;
            }
         }

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         *lowerbound = objforscip;

         if ( allint ) /* if the solution is integer, we might have found a new best solution for the MISDP */
         {
            SCIP_CALL( SCIPcheckSol(scip, scipsol, TRUE, FALSE, FALSE, FALSE, &allfeas) ); /* is this really needed ? */
            if ( allfeas )
            {
               SCIP_CALL( SCIPtrySol(scip, scipsol, TRUE, FALSE, FALSE, FALSE, &stored) );
               if (stored)
                  SCIPdebugMessage("feasible solution for MISDP found, cut node off, solution is stored \n");
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
            SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );

         SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
         *result = SCIP_SUCCESS;

         /* if all int and binary vars are integral, nothing else needs to be done */
         if ( ! allint )
         {
            int oldncuts = SCIPgetNCuts(scip);

            if ( SCIPgetNCuts(scip) > oldncuts )
               *result = SCIP_SEPARATED;

            for (i = 0; i < nvars; ++i)
            {
               SCIP_VAR* var = vars[i];
               if ( SCIPvarIsIntegral(var) && ! SCIPisFeasIntegral(scip, solforscip[i]) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
               {
                  /* we don't set a true score, we will just let the heuristic decide */
                  SCIP_CALL( SCIPaddExternBranchCand(scip, var, 10000, solforscip[i]) );
               }
            }
         }
         SCIPfreeBufferArray(scip, &solforscip);
         SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
      }

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
SCIP_DECL_RELAXEXEC(relaxExecSdp)
{
   SCIP_RELAXDATA* relaxdata;
   int nconss;
   int i;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real* ubs;
   SCIP_Bool cutoff;
   int sense;
   SCIP_SOL* scipsol;
   SCIP_Bool stored;
#ifdef SCIP_EVEN_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
#endif

   /* construct the lp and make sure, that everything is where it should be */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   if ( cutoff )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   /* very important to call flushLP */
   SCIP_CALL( SCIPflushLP(scip) );

   /* get varmapper */
   nconss = SCIPgetNConss(scip);

#ifdef SCIP_EVEN_MORE_DEBUG
   for (i = 0; i < nvarsfordebug; i++)
   {
      SCIPdebugMessage("variable %s: status = %u, integral = %u, bounds = [%f, %f] \n", SCIPvarGetName(varsfordebug[i]), SCIPvarGetStatus(varsfordebug[i]),
         SCIPvarIsIntegral(varsfordebug[i]), SCIPvarGetLbLocal(varsfordebug[i]), SCIPvarGetUbLocal(varsfordebug[i]));
   }
#endif

   if ( nconss == 0 )
   {
      /* if there are no constraints, there is nothing to do */
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( allVarsFixed(scip) )
   {
      SCIP_Bool feasible;

      /* if all variables, really all, are fixed, we give this fixed solution to SCIP */

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ubs, nvars) );

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert( SCIPisEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
      }
      sense = SCIPgetObjsense(scip);
      if ( sense == -1 )
         *lowerbound *= -1;

      SCIPdebugMessage("EVERYTHING IS FIXED, objective value = %f\n", *lowerbound);

      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      /* check if the solution really is feasible */
      SCIP_CALL( SCIPcheckSol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, &feasible) );

      if ( feasible )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &scipsol, FALSE, FALSE, FALSE, FALSE, &stored) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
      }

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

   SCIP_CALL( calc_relax(scip, relaxdata, result, lowerbound));
   relaxdata->sdpcalls++;

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** this method is called after presolving is finished, at this point the varmapper is prepared and the SDP Interface is initialized and gets
 *  the SDP information from the constraints */
static
SCIP_DECL_RELAXINIT(relaxInitSolSdp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RETCODE retcode;
   int nvars;
   SCIP_VAR** vars;
   SCIP_Real epsilon;
   SCIP_Real feastol;
   int threads;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );

   relaxdata->objval = 0.0;
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
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolverepsilon", &epsilon) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_EPSILON, epsilon);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: epsilon setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/sdpsolverfeastol", &feastol) );
   retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_FEASTOL, feastol);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: feastol setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/threads", &threads) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_THREADS, threads);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: threads setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/* copy method for sdp relaxation handler (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopySdp)
{
   assert( scip != NULL );
   assert( relax != NULL );
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   SCIPincludeRelaxSdp(scip);

   return SCIP_OKAY;
}

/** reset the relaxator's data */
static
SCIP_DECL_RELAXEXIT(relaxExitSdp)
{
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );

   SCIPdebugMessage("Exiting Relaxation Handler \n");

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   if ( relaxdata->varmapper != NULL )
   {
      SCIP_CALL( SdpVarmapperFree(scip, &(relaxdata->varmapper)) );
   }

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->sdpiterations = 0;
   relaxdata->sdpcalls = 0;
   relaxdata->lastsdpnode = 0;
   SCIP_CALL( SCIPsdpiClear(relaxdata->sdpi) );

   return SCIP_OKAY;
}

/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSdp)
{
   SCIP_RELAXDATA* relaxdata;

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   if ( relaxdata->sdpi != NULL )
   {
      SCIP_CALL( SCIPsdpiFree(&(relaxdata->sdpi)) );
   }

   SCIPfreeMemory(scip, &relaxdata);

   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** creates the SDP relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax;
   SCIP_SDPI* sdpi;

   assert( scip != NULL );

   /* create SDP relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );
   SCIP_CALL( SCIPsdpiCreate(&sdpi, NULL, SCIPblkmem(scip)) );

   relaxdata->sdpi = sdpi;
   relaxdata->lastsdpnode = -1;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecSdp, relaxdata) );
   assert( relax != NULL );

   /* include additional callbacks */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitSolSdp) );
   SCIP_CALL( SCIPsetRelaxExit(scip, relax, relaxExitSdp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeSdp) );
   SCIP_CALL( SCIPsetRelaxCopy(scip, relax, relaxCopySdp) );

   /* add parameters for SDP-solver */
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolverepsilon", "the stopping criterion for the duality gap the sdpsolver should use",
         &(relaxdata->sdpsolverepsilon), TRUE, DEFAULT_SDPSOLVEREPSILON, 1e-20, 0.001, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/sdpsolverfeastol", "the feasibility tolerance the SDP solver should use for the SDP constraints",
         &(relaxdata->sdpsolverfeastol), TRUE, DEFAULT_SDPSOLVERFEASTOL, 1e-17, 0.001, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/threads", "number of threads used for SDP solving",
         &(relaxdata->threads), TRUE, DEFAULT_THREADS, 1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/sdpinfo", "Should the SDP solver output information to the screen?",
         &(relaxdata->sdpinfo), TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/objlimit", "Should an objective limit be given to the SDP-Solver?",
         &(relaxdata->objlimit), TRUE, DEFAULT_OBJLIMIT, NULL, NULL) );

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
   *objval = relaxdata->objval;

   return SCIP_OKAY;
}

/** returns values of all variables in the solution of the current SDP relaxation, if the last SDP relaxation was successfully solved */
SCIP_RETCODE SCIPrelaxSdpGetRelaxSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax,              /**< SDP relaxator to get solution for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP relaxation solved successfully */
   SCIP_Real*            solarray,           /**< array to insert the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray, if this is less than nvars, it will be overwritten with the needed length and a
                                               *  debug message is thrown */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );
   assert( success != NULL );
   assert( solarray != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

   *success = relaxdata->origsolved;

   if ( *sollength >= SCIPgetNVars(scip) )
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, NULL, solarray, sollength) );
   else
   {
      SCIPdebugMessage("Called SCIPrelaxSdpGetRelaxSol with an array that wasn't big enough, needed length %d, given %d!\n", SCIPgetNVars(scip), *sollength);
      *sollength = SCIPgetNVars(scip);
   }

   return SCIP_OKAY;
}

/** get the number of the SCIP-node to which the current SDP solution belongs */
long int SCIPrelaxSdpGetSdpNode(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return SCIPrelaxGetData(relax)->lastsdpnode;
}

/** was the original problem solved for the last SDP-Node (or a penalty formulation) ? */
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return SCIPsdpiSolvedOrig(relaxdata->sdpi);
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
