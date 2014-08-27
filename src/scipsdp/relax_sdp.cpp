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

#include "SdpProblem.h"                 // for SdpProblem
#include "SdpVarmapper.h"               // for SdpVarMapper
#include "sdpi/sdpi_general.h"          // for SDP-Interface
#include "SdpCone.h"                    // for Iterators
#include "cons_sdp.h"                   // for cons_check

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
};

/** removes all indices j from an SDP block for which both row j and column j are completely empty/zero (for this all entries of sdpval
 *  really need to be nonzero, because these will not be checked) */
static
SCIP_RETCODE removeEmptySdpRowCols( //TODO: move this to sdpi and also check for empty lp rows (or add another function to do this)
   SCIP*                 scip,               /**< SCIP instance */
   const int             block,              /**< index of the SDP-block this is called for */
   const int             oldblocksize,       /**< old block size of the SDP-block */
   int*                  newblocksize,       /**< block size after removing redundant indices */
   const int             blockstartind,      /**< starting index of the block in the sdp-nonzero arrays */
   const int             nextblockstartind,  /**< starting index of the next block in the sdp-nonzero arrays */
   const int             constblockstartind, /**< starting index of the block in the constant sdp-nonzero arrays */
   const int             constnextblstartind,/**< starting index of the next block in the sdp-nonzero arrays */
   int*                  sdprowind,          /**< row-indices of the sdp-nonzeroes */
   int*                  sdpcolind,          /**< column-indices of the sdp-nonzeroes */
   int*                  sdpconstrowind,     /**< row-indices of the sdp-const-nonzeroes */
   int*                  sdpconstcolind      /**< column-indices of the sdp-const-nonzeroes */
   )
{
   int i;
   int j;
   int* deleted; /* if deleted[i] = -1, row and column i are both empty and will thus be deleted, otherwise this gives the number of
                  * indices before the current one that were deleted and therefore the number of positions this index needs to be shifted */
   int ndeleted; /* after the for-queue below this will give the total number of indices that will be deleted (this is also used to
                  * compute the positive entries of the deleted-array) */

   assert ( scip != NULL );
   assert ( block >= 0 );
   assert ( oldblocksize >= 0 );
   assert ( newblocksize != NULL );
   assert ( blockstartind >= 0 );
   assert ( nextblockstartind >= blockstartind );
   assert ( sdprowind != NULL );
   assert ( sdpcolind != NULL );

   /* compute the deleted-array */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &deleted, oldblocksize));
   ndeleted = 0;

   for (i = 0; i < oldblocksize; i++)
   {
      deleted[i] = -1; /* initialize this with -1, it will be overwritten if the index is found */
      /* go over the sdp-arrays */
      for (j = blockstartind; j < nextblockstartind; j++)
      {
         if (sdprowind[j] == i || sdpcolind[j] == i)
         {
            deleted[i] = ndeleted;  /* there exists an entry for this index, so it will not be deleted, but the number of positions
                                     * for shifting is saved */
            break;
         }
      }
      /* if nothing was found in the sdp-arrays go over the sdp-const arrays */
      if (deleted[i] == -1)
      {
         for (j = constblockstartind; j < constnextblstartind; j++)
         {
            if (sdpconstrowind[j] == i || sdpconstcolind[j] == i)
            {
               deleted[i] = ndeleted;  /* there exists an entry for this index, so it will not be deleted, but the number of positions
                                        * for shifting is saved */
               break;
            }
         }
      }
      /* the index wasn't found in the arrays, so it can be deleted */
      if (deleted[i] == -1)
      {
         ndeleted++;
#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("deleted the following index in block %d becase of empty row & col: %d\n", block, i);
#endif
      }
   }

   /* now shift all row- & column-indices according to the deleted-array, if at least one index was deleted */
   if (ndeleted > 0)
   {
      for (j = blockstartind; j < nextblockstartind; j++)
      {
         assert ( deleted[sdprowind[j]] > -1 );
         sdprowind[j] = sdprowind[j] - deleted[sdprowind[j]];

         assert ( deleted[sdpcolind[j]] > -1 );
         sdpcolind[j] = sdpcolind[j] - deleted[sdpcolind[j]];
      }
      for (j = constblockstartind; j < constnextblstartind; j++)
      {
         assert ( deleted[sdpconstrowind[j]] > -1 );
         sdpconstrowind[j] = sdpconstrowind[j] - deleted[sdpconstrowind[j]];

         assert ( deleted[sdpconstcolind[j]] > -1 );
         sdpconstcolind[j] = sdpconstcolind[j] - deleted[sdpconstcolind[j]];
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMessage("number of deleted indices in block %d because of empty rows & cols: %d\n", block, ndeleted);
#endif
   /* finally update the blocksize */
   *newblocksize = oldblocksize - ndeleted;

   SCIPfreeBlockMemoryArray(scip, &deleted, oldblocksize);

   return SCIP_OKAY;
}

/** inserts all the SDP data into the corresponding SDP Interface */
static
SCIP_RETCODE putSdpDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   int i;
   int j;
   int nvars;
   SCIP_VAR ** vars;
   SCIP_VAR ** blockvars;
   SCIP_CONS** conss;
   int ncons;
   int ind = 0;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int nsdpblocks;
   int* sdpblocksizes;
   int sdpconstnnonz;
   int* constrow;
   int* constcol;
   SCIP_Real* constval;
   int sdpnnonz;
   int constnnonzcounter;
   int** row;
   int** col;
   SCIP_Real** val;
   int** blockcol;
   int** blockrow;
   SCIP_Real** blockval;
   SCIP_CONSHDLR* conshdlr;
   int blocknnonz;
   int varind;
   int blocknvars;
   int* nblockvarnonz;
   int* nvarnonz;
   int* nconstblocknonz;

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

   ncons = SCIPgetNConss(scip);
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &conss, ncons));
   conss = SCIPgetConss(scip);

   /* count the number of sdpblocks and compute the number of nonzeros */
   nsdpblocks = 0;
   sdpnnonz = 0;
   sdpconstnnonz = 0;

   for (i=0; i < ncons; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         nsdpblocks++;

         SCIP_CALL(SCIPconsSdpGetNNonz(scip, conss[i], &blocknnonz, &constnnonzcounter));
         sdpnnonz += blocknnonz;
         sdpconstnnonz += constnnonzcounter;
      }
   }

   /* create the sdp- and sdpconst-arrays */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &nvarnonz, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &nblockvarnonz, nvars * nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &nconstblocknonz, nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &col, nvars * nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &row, nvars * nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &val, nvars * nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &blockcol, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &blockrow, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &blockval, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &constcol, nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &constrow, nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &constval, nsdpblocks));

   /* get the SDP-data */
   ind = 0; /* index of the current sdp block in the complete sdp */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &blockvars, nvars));

   for (i = 0; i < ncons; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         varind = 0;

         /* ?????????????????? */
         SCIPconsSdpGetData(scip, conss[i], &blocknvars, &blocknnonz, &sdpblocksizes[i], &nvars, nvarnonz, blockcol,
            blockrow, blockval, blockvars, &(nconstblocknonz[ind]), &(constcol[ind]), &(constrow[ind]), &(constval[ind]));

         /* nvars would have been overwritten if the space in the given arrays hadn't been sufficient */
         assert ( nvars == SCIPgetNVars(scip) );

         /* ??????????????? loop over variables in constraint ! */

         /* update the variable indices of the current block to the global variable indices */
         for (j = 0; j < nvars; j++)
         {
            if ( vars[j] == blockvars[varind])  /* this variable does exist in this block */
            {
               nblockvarnonz[ind * nvars + j] = nvarnonz[varind]; /* take the value given by the constraint */
               col[ind * nvars + j] = blockcol[varind];
               row[ind * nvars + j] = blockrow[varind];
               val[ind * nvars + j] = blockval[varind];
               varind++;
            }
            else                                /* this variable does not exist in this block */
               nblockvarnonz[ind * nvars + j] = 0; /* then there are no nonzeros for this variable in this block */
         }

         ind++;
      }
   }

   SCIPfreeBlockMemoryArray(scip, &conss, ncons);
   SCIPfreeBlockMemoryArray(scip, &blockvars, nvars);
   SCIPfreeBlockMemoryArray(scip, &nvarnonz, nvars);

   SCIP_CALL(SCIPsdpiLoadSDP(sdpi, nvars,  obj,  lb,  ub, nsdpblocks,
                            sdpblocksizes, sdpconstnnonz, nconstblocknonz , constrow,
                            constcol, constval, sdpnnonz, nblockvarnonz,
                            row, col,  val, 0,
                            NULL, 0, NULL, NULL, NULL)); /* insert the SDP part, add an empty LP part */

   SCIPfreeBlockMemoryArray(scip, &obj, nvars);
   SCIPfreeBlockMemoryArray(scip, &lb, nvars);
   SCIPfreeBlockMemoryArray(scip, &ub, nvars);
   SCIPfreeBlockMemoryArray(scip, &sdpblocksizes, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &constval, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &constrow, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &constcol, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &blockval, nvars);
   SCIPfreeBlockMemoryArray(scip, &blockrow, nvars);
   SCIPfreeBlockMemoryArray(scip, &blockcol, nvars);
   SCIPfreeBlockMemoryArray(scip, &nconstblocknonz, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &nblockvarnonz, nsdpblocks * nvars);
   SCIPfreeBlockMemoryArray(scip, &row, nsdpblocks * nvars);
   SCIPfreeBlockMemoryArray(scip, &col, nsdpblocks * nvars);
   SCIPfreeBlockMemoryArray(scip, &val, sdpnnonz);

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
   int ncons;
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
   ncons = 0; /* this will be increased for each finite lhs and rhs */

   for (i = 0; i < nrows; i++)
   {
      blocknnonz = SCIProwGetNNonz(rows[i]);

      rowvals = SCIProwGetVals(rows[i]);
      rowcols = SCIProwGetCols(rows[i]);
      sciplhs = SCIProwGetLhs(rows[i])- SCIProwGetConstant(rows[i]);
      sciprhs = SCIProwGetRhs(rows[i])- SCIProwGetConstant(rows[i]);

      /* add row >= lhs if lhs is finite */
      if (! (SCIPsdpiIsInfinity(sdpi, sciplhs))) /* TODO: this will also ignore lhs=+infty which renders the problem infeasible, can this happen/is this caught elsewhre ? */
      {
         for (j = 0; j < blocknnonz; j++)
         {
            colind[nnonz] = SdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = ncons;
            val[nnonz] = rowvals[j];
            nnonz++;
         }
         rhs[ncons] = sciplhs;
         ncons++;
      }

      /* add -row >= -rhs if rhs is finite */
      if (! (SCIPsdpiIsInfinity(sdpi, sciprhs))) /* TODO: this will also ignore rhs=-infty which renders the problem infeasible, can this happen/is this caught elsewhre ? */
      {
         for (j = 0; j < blocknnonz; j++)
         {
            colind[nnonz] = SdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = ncons;
            val[nnonz] = -rowvals[j];
            nnonz++;
         }
         rhs[ncons] = -sciprhs;
         ncons++;
      }
   }

   /* these arrays are no longer needed */
   SCIPfreeBlockMemoryArray(scip, &rows, nrows);
   SCIPfreeBlockMemoryArray(scip, &rowvals, nvars);
   SCIPfreeBlockMemoryArray(scip, &rowcols, nvars);

   /* reallocate some arrays depending on the number of lpnnonz and lhs/rhs */
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &rhs, 2 * nrows, ncons));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &rowind, 2 * scipnnonz, nnonz));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &colind, 2 * scipnnonz, nnonz));
   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &val, 2 * scipnnonz, nnonz));

   /* delete the old LP-block from the sdpi */
   SCIP_CALL(SCIPsdpiGetNLPRows(sdpi, &nrowssdpi));
   SCIP_CALL(SCIPsdpiDelLPRows(sdpi, 0, nrowssdpi - 1));

   /* add the LP-block to the sdpi */
   SCIP_CALL(SCIPsdpiAddLPRows(sdpi, ncons, rhs, nnonz, (const int*)rowind, (const int*)colind, val));

   /* free the remaining arrays */
   SCIPfreeBlockMemoryArray(scip, &rhs, ncons);
   SCIPfreeBlockMemoryArray(scip, &rowind, nnonz);
   SCIPfreeBlockMemoryArray(scip, &colind, nnonz);
   SCIPfreeBlockMemoryArray(scip, &val, nnonz);

   return SCIP_OKAY;
}

/** gets the solution from sdpi and adds the fixed variables to the solution of the SDP Solver and transforms it to the variable indices used in SCIP */
static
SCIP_RETCODE getTransformedSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarmapper*         varmapper,          /**< data about fixed variables */
   SCIP_Real*            solforscip,         /**< the solution indexed by SCIP variable indices */
   SCIP_Real*            objforscip,         /**< the objective value including the fixed variables */
   SCIP_Bool*            allint              /**< do all variables that should be integer have integer values */
   )
{
   int nsdpvars;
   SCIP_Real* sdpsol;
   SCIP_Real sdpobj;
   SCIP_VAR** vars;
   int nvars;
   int ind;
   int i;

   assert ( scip != NULL );
   assert ( sdpi != NULL );
   assert ( varmapper != NULL );
   assert ( solforscip != NULL );
   assert ( objforscip != NULL );
   assert ( allint != NULL );

   SCIP_CALL( SCIPsdpiGetNVars(sdpi, &nsdpvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sdpsol, nsdpvars) );

   SCIP_CALL( SCIPsdpiGetSol(sdpi, &sdpobj, sdpsol, nsdpvars) ); /* get both the objective and the solution from the SDP solver */

   *objforscip = sdpobj; /* initialize the objective with that of the SDP, later the objective parts of the fixed variables will be added to this */
   *allint = TRUE; /* initialize this as true until a variable is found that is non-integer but should be */

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars (scip);

   for (i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var = vars[i];
      assert( var != 0 );

      ind = SdpVarmapperGetSdpIndex(varmapper, var);
      if ( ind == -1 )
      {
         assert (SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

         solforscip[i] = SCIPvarGetUbLocal(var); /* this variable was fixed, so it is equal to its upper bound (which equals the lower bound) */
         SCIP_Real obj = SCIPvarGetObj(var);
         if ( ! SCIPisZero(scip, solforscip[i]) && ! SCIPisZero(scip, obj) )
            *objforscip += solforscip[i] * obj; /* add this fixed variables contribution to the objective value */
      }
      else
      {
         solforscip[i] = sdpsol[ind]; /* get the value of the corresponding variable in the SDP */
         if( *allint && SCIPvarIsIntegral(var) && !(SCIPisIntegral(scip, sdpsol[ind])) )
            *allint = FALSE;   /* a variable that should be integral but isnot was found, so allint is set to false */
         /* nothing needs to be done here with the objective value, because it was already initialized with the value from the SDP solver */
      }
   }

   SCIPfreeBlockMemoryArray(scip, &sdpsol, nsdpvars);

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
   SCIP_CONS** conss;
   int ncons;

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

   SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
   SCIP_CALL( getTransformedSol(scip, sdpi, varmapper, solforscip, &objforscip, &allint) );

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
         SCIP_CALL( getTransformedSol(scip, sdpi, varmapper, solforscip, &objforscip, &allint) );

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         /* if called with penalty formulation check if the solution really is feasible */
         if( withpenalty )
         {
            SCIP_RESULT conefeas;

            ncons = SCIPgetNConss(scip);
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &conss, ncons));
            conss = SCIPgetConss(scip);

            for (i = 0; i < ncons; ++i)
            {
               SCIP_CALL( consCheckSdp(scip, conss[i], scipsol, FALSE, TRUE, FALSE, &conefeas) );
               if ( conefeas == SCIP_INFEASIBLE )
               {
                  solisfeas = FALSE;
                  break;
               }
            }
            SCIPfreeBlockMemoryArray(scip, &conss, ncons);
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_RELAXDATA* relaxdata;
   int ncons;
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
   ncons = SCIPgetNConss(scip);
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &conss, ncons));
   conss = SCIPgetConss(scip);

   SdpVarmapper* varmapper;

   for (i=0; i < ncons; i++)
   {
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         varmapper = &((SCIPconshdlrGetData(conshdlr))->varmapper);
         break;
      }
   }

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

   int nlprows;
   nlprows = SCIPgetNCuts(scip) + SCIPgetNLPRows(scip);

   if ( ncons == 0 )
   {
      //if there are no constraints, there is nothing to do
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if (allVarsFixed(scip))
   {
      // if all variables, really all, are fixed, I can't solve an sdp, because there is no interior point in this case, result is success and I'm separating the solution (the upper or lower bounds on a variable)
      SCIPdebugMessage("EVERYTHING IS FIXED\n");
      SCIP_VAR** vars = SCIPgetVars(scip);
      const int nvars = SCIPgetNVars(scip);

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

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &indsforsdpi, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &obj, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lb, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ub, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vars, nvars));

   vars = SCIPgetVars(scip);

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
   SCIPfreeBlockMemoryArray(scip, &vars, nvars);

   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, varmapper) );


   SCIP_CALL( calc_relax(scip, relaxdata->sdpi, varmapper, FALSE, 0.0, result, lowerbound));


   return SCIP_OKAY;
}

/** this method is called after presolving is finished, at this point the SDP Interface is initialized and gets the SDP information from the constraints */
static
SCIP_DECL_RELAXINIT(relaxInitSolSDP)
{
   SCIP_RELAXDATA* relaxdata;

   assert ( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   SCIP_CALL(putSdpDataInInterface(scip, relaxdata->sdpi));

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
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;
   SCIP_SDPI sdpi;

   assert ( scip != NULL );

   /* create SDP relaxator data */
   SCIP_CALL(SCIPallocBlockMemory(scip, &relaxdata));
   SCIP_CALL(SCIPallocBlockMemory(scip, &sdpi));
   SCIP_CALL(SCIPsdpiCreate(sdpi, NULL, SCIPblkmem(scip)));

   relaxdata->sdpi = sdpi;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecSDP, relaxdata) );
   assert( relax != NULL );

   /* add xyz relaxator parameters */

   return SCIP_OKAY;
}

/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSDP)
{
SCIP_RELAXDATA* relaxdata;

relaxdata = SCIPrelaxGetData(relax);
assert(relaxdata != NULL);

SCIPfreeMemory(scip, &relaxdata);
SCIPrelaxSetData(relax, NULL);

return SCIP_OKAY;
}

