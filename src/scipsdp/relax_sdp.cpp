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
//#define SCIP_MORE_DEBUG /* shows number of deleted empty cols/rows and complete solution for every relaxation */

#include "relax_sdp.h"

#include <cassert>                      // for assert
#include <cstdio>                       // for NULL, printf

#include "SdpProblem.h"                 // for SdpProblem
#include "SdpVarMapper.h"               // for SdpVarMapper
#include "sdpi/sdpi.h"                  // for SDP-Interface
#include "SdpCone.h"                    // for Iterators
#include "objconshdlr_sdp.h"          // for cons_check

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


/*
* check if variable bounds are fulfilled
static
SCIP_RETCODE check_bounds(
   SCIP*                 scip,               *< SCIP data structure
   int                   nvars,              *< number of variables
   SCIP_VAR**            var,                *< variables
   SCIP_SOL*             scipsol,            *< solution of scip-type to check bounds for
   SCIP_Bool*            sol_is_feas         *< pointer to store if solution is feasible
   )
{
   assert ( scip != NULL );
   assert ( nvars >= 0 );
   assert ( var != NULL );
   assert ( scipsol != NULL );
   assert ( sol_is_feas != NULL );

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
}*/

/** removes all indices j from an SDP block for which both row j and column j are completely empty/zero (for this all entries of sdpval
 *  really need to be nonzero, because these will not be checked) */
static
SCIP_RETCODE removeEmptyRowCols(
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

/** inserts all the data into the corresponding SDP Interface */
static
SCIP_RETCODE putDataInInterface(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpProblem*           problemdata,        /**< data structure with problem-data of a specific node */
   SdpVarMapper*         varmapper           /**< data about fixed variables */
   )
{
   int i;
   int nvars;
   SCIP_VAR ** vars;
   SdpCone* sdpcone;
   int ind;
   int newblocksize;
   SCIP_VAR** fixedvars;
   int nfixedvars;
   double* fixedvalues;
   SdpCone::element el;
   SCIP_Real* sdpvar; /* this could as well be int, but SCIP only knows SCIPsortRealRealIntInt, but not IntRealIntInt or IntIntIntReal */
   int endindex;
   int constendindex;
   int nextindaftervar;
   int formatindsize;
   const int* formatind;
   const int* forconstraint;
   const SCIP_Real* forvals;
   SCIP_Real* obj;
   SCIP_Real* lb;
   SCIP_Real* ub;
   int nsdpblocks;
   int* sdpblocksizes;
   int sdpconstnnonz;
   int* sdpconstbegblock;
   int* sdpconstrowind;
   int* sdpconstcolind;
   SCIP_Real* sdpconstval;
   int sdpnnonz;
   int* sdpbegvarblock;
   int* sdprowind;
   int* sdpcolind;
   SCIP_Real* sdpval;
   int nlpcons;
   SCIP_Real* lprhs;
   int lpnnonz;
   int* lprowind;
   int* lpcolind;
   SCIP_Real* lpval;

   assert ( scip != NULL );
   assert ( sdpi != NULL );
   assert ( problemdata != NULL );
   assert ( varmapper != NULL );

   nsdpblocks = problemdata->get_nsdpcones();
   nvars = varmapper->get_sdp_nvars();

   vars = SCIPgetVars(scip);

   /* prepare arrays of objective values and bounds */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &obj, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lb, nvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ub, nvars));

   for (i = 0; i < nvars; i++)
   {
         obj[i] = SCIPvarGetObj(varmapper->get_scip_var(i));
         lb[i] = SCIPvarGetLbLocal(varmapper->get_scip_var(i));//-SCIPsdpiInfinity(sdpi); /* this should be changed to SCIPvarGetLbLocal after SdpProblem was changed to add bounds instead of converting them to LP cons (or just not add them at all) */
         ub[i] = SCIPvarGetUbLocal(varmapper->get_scip_var(i));//SCIPsdpiInfinity(sdpi); /* "---------------------------------------------------------------------------------------------------------" */
   }

   /* get SDPBlocksizes */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpblocksizes, nsdpblocks));

   for (i = 0; i < nsdpblocks; i++)
   {
      sdpblocksizes[i] = problemdata->get_sdpcone(i)->get_blocksize();
   }

   /* get all fixed variables */
   ind = 0;
   nfixedvars = varmapper->get_nfixed();
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &fixedvars, nfixedvars));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &fixedvalues, nfixedvars));

   for (i = 0; i < SCIPgetNVars(scip); i++)
   {
      if (varmapper->get_sdp_index(vars[i]) == -1 )
      {
         fixedvars[ind] = vars[i];
         fixedvalues[ind] = SCIPvarGetUbLocal(vars[i]);
         ind++;
         assert ( SCIPisEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])) );
      }
   }

   /* compute SDPconstbegblock and SDPconstnnonz and SDPbegvarblock (only the begblock part) and SDPnnonz
    * these need to be computed now to know how much space to allocate for the other arrays
    * for this the Iterator is needed, because these numbers depend on the number of nonzeroes of fixed variables */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpconstbegblock, nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpbegvarblock, nsdpblocks * nvars));

   /* the constant part */
   ind = 0;

   for (i = 0; i < nsdpblocks; i++)
   {
      sdpconstbegblock[i] = ind;

      sdpcone = problemdata->get_sdpcone(i);

      sdpcone->fix_vars(); //TODO: is this the right place to call it ?!?!?!?!?

      for (SdpCone::RhsIterator it = sdpcone->rhs_begin(fixedvars, nfixedvars, fixedvalues); it != sdpcone->rhs_end(); ++it)
      {
         ind++;
      }
   }

   sdpconstnnonz = ind;

   /* the non-constant part */
   ind = 0;

   for (i = 0; i < nsdpblocks; i++)
   {
      sdpbegvarblock[i * nvars] = ind;

      sdpcone = problemdata->get_sdpcone(i);

      for (SdpCone::LhsIterator it = sdpcone->lhs_begin(fixedvars, nfixedvars); it != sdpcone->lhs_end(); ++it)
      {
         ind++;
      }
   }

   sdpnnonz = ind;


   /* prepare sdpconst-arrays */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpconstrowind, sdpconstnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpconstcolind, sdpconstnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpconstval, sdpconstnnonz));
   ind = 0;

   for (i = 0; i < nsdpblocks; i++)
   {
      sdpcone = problemdata->get_sdpcone(i);

      for (SdpCone::RhsIterator it = sdpcone->rhs_begin(fixedvars, nfixedvars, fixedvalues); it != sdpcone->rhs_end(); ++it)
      {
         el = *it;
         sdpconstrowind[ind] = el.row;
         sdpconstcolind[ind] = el.col;
         sdpconstval[ind] = -1 * el.val; /* these are saved in SDPCone with the same sign as the other nonzeroes to be able to combine
                                          * them for fixed variables, but the interface assumes that -A_0 is added to the SDP-Blocks */
         ind++;
      }

      if (i < nsdpblocks - 1)
      {
         assert ( sdpconstbegblock[i+1] == ind );
      }
      else
      {
         assert ( sdpconstnnonz == ind );
      }

   }

   /* prepare sdp-arrays */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdprowind, sdpnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpcolind, sdpnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpval, sdpnnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpvar, sdpnnonz)); /* in this array the variables for the entries will be stored, later this will
                                                                   * be used to compute the sdpbegvarblock-array */

   ind = 0;

   for (i = 0; i < nsdpblocks; i++)
   {
      int j;

      sdpcone = problemdata->get_sdpcone(i);

      for (SdpCone::LhsIterator it = sdpcone->lhs_begin(fixedvars, nfixedvars); it != sdpcone->lhs_end(); ++it)
      {
         el = *it;
         sdprowind[ind] = el.row;
         sdpcolind[ind] = el.col;
         sdpval[ind] = el.val;
         sdpvar[ind] = varmapper->get_sdp_index(sdpcone->get_var(el.vidx));
         ind++;
      }

      if (i < nsdpblocks - 1)
      {
         assert ( sdpbegvarblock[(i+1) * nvars] == ind );
      }
      else
      {
         assert ( sdpnnonz == ind );
      }

      /* now sort all entries belonging to this block by their variables to compute the remaining entries of sdpbegvarblock */
      if (i < nsdpblocks - 1)
      {
         endindex = sdpbegvarblock[(i+1) * nvars];
         constendindex = sdpconstbegblock[i+1];
      }
      else
      {
         endindex = sdpnnonz;
         constendindex = sdpconstnnonz;
      }

      SCIPsortRealRealIntInt(sdpvar + sdpbegvarblock[i*nvars], sdpval + sdpbegvarblock[i*nvars], sdpcolind + sdpbegvarblock[i*nvars],
                                sdprowind + sdpbegvarblock[i*nvars], endindex - sdpbegvarblock[i * nvars]); /* + sdpbegvarblock[i*nvars] makes sure that sorting
                                                                                                             * starts at the beginning of this block */

      assert ( sdpvar[sdpbegvarblock[i*nvars]] >= 0 );

      nextindaftervar = sdpbegvarblock[i*nvars];
      for (j = 0; j < nvars; j++)
      {
         while (nextindaftervar < endindex && sdpvar[nextindaftervar] == j) /* get the first index that doesn't belong to this variable */
         {
            nextindaftervar++;
         }
         if (j < nvars - 1)
         {
            sdpbegvarblock[i * nvars + j + 1] = nextindaftervar;       /* set the remaining sdbegvarblock entries equal to the computed values */
         }
         else
         {
            assert (nextindaftervar == endindex );
         }
      }

      /* remove all empty rows and columns */
      SCIP_CALL(removeEmptyRowCols(scip, i, sdpblocksizes[i], &newblocksize, sdpbegvarblock[i * nvars], endindex,
            sdpconstbegblock[i], constendindex, sdprowind, sdpcolind, sdpconstrowind, sdpconstcolind));
   }

   SCIPfreeBlockMemoryArray(scip, &sdpvar, sdpnnonz);
   SCIPfreeBlockMemoryArray(scip, &fixedvars, nfixedvars);
   SCIPfreeBlockMemoryArray(scip, &fixedvalues, nfixedvars);

   /* prepare LP arrays */
   nlpcons = problemdata->get_size_lpblock();
   formatindsize = problemdata->get_for_matind_size();

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lprowind, formatindsize)); /* this is the worst case length if all right hand sides are zero,
                                                                                  * otherwise the last (number of rhs-nonzeroes) entries will stay empty */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lpcolind, formatindsize));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lpval, formatindsize));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lprhs, nlpcons));


   /* initialize rhs with zero TODO: change this if not only nonzeroes are saved in SDPProblem */
   for (i = 0; i < nlpcons; i++)
   {
      lprhs[i] = 0;
   }

   /*  partition the nonzeroes into right hand sides and "traditional" nonzeroes */
   formatind = problemdata->get_for_matind();
   forconstraint = problemdata->get_for_constraint();
   forvals = problemdata->get_for_vals();


   ind = 0;
   for (i = 0; i < formatindsize; i++)
   {
      if (forconstraint[i] == 0) /* this means it's a right hand side */
      {
         lprhs[formatind[i]] = -1 * forvals[i]; /* in SDPCone these are saved as lhs-nonzeroes with negative sign and variable 0, but
                                                 * the interface expects them as right hand sides, so the multiplication with -1 needs
                                                 * to be reversed */
      }
      else /* this means it's a "real" nonzero */
      {
         lprowind[ind] = formatind[i];
         lpcolind[ind] = forconstraint[i] - 1; /* these start at 1 in the SDP cone but should start at 0 for the interface */
         lpval[ind] = forvals[i];
         ind++;
      }
   }


   lpnnonz = ind;

   SCIP_CALL(SCIPsdpiLoadSDP(sdpi, nvars, (const SCIP_Real*) obj, (const SCIP_Real*) lb, (const SCIP_Real*) ub, nsdpblocks,
                            (const int*) sdpblocksizes, sdpconstnnonz, (const int*) sdpconstbegblock, (const int*) sdpconstrowind,
                            (const int*) sdpconstcolind, (const SCIP_Real*) sdpconstval, sdpnnonz, (const int*) sdpbegvarblock,
                            (const int*) sdprowind, (const int*) sdpcolind, (const SCIP_Real*) sdpval, nlpcons,
                            (const SCIP_Real*) lprhs, lpnnonz, (const int*) lprowind, (const int*) lpcolind,
                            (const SCIP_Real*) lpval));

   SCIPfreeBlockMemoryArray(scip, &obj, nvars);
   SCIPfreeBlockMemoryArray(scip, &lb, nvars);
   SCIPfreeBlockMemoryArray(scip, &ub, nvars);
   SCIPfreeBlockMemoryArray(scip, &sdpblocksizes, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &sdpconstbegblock, nsdpblocks);
   SCIPfreeBlockMemoryArray(scip, &sdpconstrowind, sdpconstnnonz);
   SCIPfreeBlockMemoryArray(scip, &sdpconstcolind, sdpconstnnonz);
   SCIPfreeBlockMemoryArray(scip, &sdpconstval, sdpconstnnonz);
   SCIPfreeBlockMemoryArray(scip, &sdpbegvarblock, nsdpblocks * nvars);
   SCIPfreeBlockMemoryArray(scip, &sdprowind, sdpnnonz);
   SCIPfreeBlockMemoryArray(scip, &sdpcolind, sdpnnonz);
   SCIPfreeBlockMemoryArray(scip, &sdpval, sdpnnonz);
   SCIPfreeBlockMemoryArray(scip, &lprhs, nlpcons);
   SCIPfreeBlockMemoryArray(scip, &lprowind, formatindsize);
   SCIPfreeBlockMemoryArray(scip, &lpcolind, formatindsize);
   SCIPfreeBlockMemoryArray(scip, &lpval, formatindsize);

   return SCIP_OKAY;
}

/** gets the solution from sdpi and adds the fixed variables to the solution of the SDP Solver and transforms it to the variable indices used in SCIP */
static
SCIP_RETCODE getTransformedSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SdpVarMapper*         varmapper,          /**< data about fixed variables */
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

      ind = varmapper->get_sdp_index(var);
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
         feasible = true;
      }
      /* Now the objective > 0, so it is feasible w.r.t. the LP-part, but there is no psd solution. */
      else if ( SCIPsdpiIsDualFeasible(sdpi) || SCIPsdpiIsDualInfeasible(sdpi) )
      {
         SCIPdebugMessage("Verified that a problem is infeasible with a penalty-only-formulation.\n");
         feasible = false;
      }
      else
      {
         SCIPerrorMessage("Cannot decide whether subproblem is feasible!"); /* can't be reached */
         SCIPABORT();
         return SCIP_ERROR;
      }
   }
   else
   {
      SCIPerrorMessage("Even when using a penalty-only-formulation the SDP-Solver didnot converge!");
      SCIPABORT();
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** calculate relaxation and process the relaxation results, may call itself recursively once to try again
 *  with a penalty formulation (or more often if the penalty formulation turns out to be unbounded) */
static
SCIP_RETCODE calc_relax(
   SCIP_SDPI*            sdpi,               /**< SDP-Interface structure */
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result,             /**< pointer to store result of relaxation process */
   SCIP_Real*            lowerbound,         /**< pointer to store lowerbound */
   SdpProblem*           problemdata,        /**< data structure with problem-data of a specific node */
   SdpVarMapper*         varmapper,          /**< varmapper class data */
   SCIP_Bool             withpenalty,        /**< should a penalty formulation be used */
   SCIP_Real             penaltyparam        /**< parameter for penalty formulation, if 0 the normal SDP is solved */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int i;

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
      SCIP_CALL(SCIPsdpiSolve(sdpi));
   }
   else
   {
      SCIP_CALL(SCIPsdpiSolvePenalty(sdpi, penaltyparam, TRUE));
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
      SCIPdebugMessage("DSDP returned primaldualunknown ! ");
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
            SCIP_CALL( calc_relax(sdpi, scip, result, lowerbound, problemdata, varmapper, TRUE, 10 * penaltyparam) );

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

            for (i = 0; i < problemdata->get_nsdpcones(); ++i)
            {
               SCIP_CALL( cons_check(scip, problemdata->get_sdpcone(i), scipsol, FALSE, TRUE, FALSE, &conefeas) );
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
               SCIP_CALL(calc_relax(sdpi, scip, result, lowerbound, problemdata, varmapper, TRUE, 10.0 * penaltyparam));
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
      SCIP_CALL(calc_relax(sdpi, scip, result, lowerbound, problemdata, varmapper, TRUE, 10 * penaltyparam));

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
   bool feasible;
   SCIP_CALL( relaxIsFeasible(scip, sdpi, feasible) );

   if ( feasible )
   {
      /* try again with penalty formulation */
      SCIP_CALL(calc_relax(sdpi, scip, result, lowerbound, problemdata, varmapper, TRUE, 1.0)); /* TODO: think about penalty parameter */
   }
   else
   {
      /* penalty-only-formulation showed, that the problem is infeasible */
      *result = SCIP_CUTOFF;
   }

   return SCIP_OKAY;
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

   // very important to call flushLP
   SCIP_CALL( SCIPflushLP(scip) );

   SdpVarMapper* varmapper;
   varmapper = new SdpVarMapper(scip);
   SCIP_CALL(varmapper->init());

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

   if (varmapper->get_allfixed())
   {
      // if all variables, really all, are fixed, I can't solve an sdp, because there is no interior point in this case, result is success and I'm separating the solution (the upper or lower bounds on a variable
      SCIPdebugMessage("EVERYTHING IS FIXED\n");
      SCIP_VAR** vars = SCIPgetVars(scip);
      const int nvars = SCIPgetNVars(scip);

      SCIP_Real* ubs;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &ubs, nvars));

      *lowerbound = 0.0;
      for (int i = 0; i < nvars; i++)
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
      delete problemdata;
      SCIP_CALL(varmapper->exit());
      delete varmapper;

      return SCIP_OKAY;
   }

   SCIP_SDPI* sdpi;
   BMS_BLKMEM* blkmem = SCIPblkmem(scip);
   SCIP_CALL(SCIPsdpiCreate(&sdpi, NULL, blkmem));

   SCIP_CALL( putDataInInterface(scip, sdpi,problemdata, varmapper) );

   SCIP_CALL( calc_relax(sdpi, scip, result, lowerbound, problemdata, varmapper, FALSE, 0.0));

   SCIP_CALL(varmapper->exit());

   delete varmapper;

   SCIP_CALL(SCIPsdpiFree(&sdpi));

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
