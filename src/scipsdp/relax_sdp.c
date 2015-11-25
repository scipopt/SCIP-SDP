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

/*#define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG   *//* displays complete solution for each relaxation */
/* #define SCIP_EVEN_MORE_DEBUG  *//* shows number of deleted empty cols/rows for every relaxation and variable status &
 * bounds as well as all constraints in the beginning */

#include "relax_sdp.h"

#include "assert.h"                     /*lint !e451*/
#include "string.h"                     /* for strcmp */

#include "SdpVarmapper.h"
#include "sdpi/sdpi.h"
#include "scipsdp/cons_sdp.h"
#include "scipsdp/cons_savedsdpsettings.h"


#define RELAX_NAME                  "SDP"
#define RELAX_DESC                  "SDP relaxator"
#define RELAX_PRIORITY              1
#define RELAX_FREQ                  1

#define DEFAULT_SDPSOLVEREPSILON    1e-5     /**< the stopping criterion for the duality gap the sdpsolver should use */
#define DEFAULT_SDPSOLVERFEASTOL    1e-4     /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
#define DEFAULT_PENALTYPARAM        -1       /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#define DEFAULT_LAMBDASTAR          -1       /**< the parameter lambda star used by SDPA to set the initial point */
#define DEFAULT_MAXPENALTYPARAM     -1       /**< the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
#if 0
#define DEFAULT_THREADS             1        /**< number of threads used for SDP solving */
#endif
#define DEFAULT_OBJLIMIT            FALSE    /**< should an objective limit be given to the SDP-Solver ? */
#define DEFAULT_RESOLVE             TRUE     /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
#define DEFAULT_TIGHTENVB           FALSE    /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
#define DEFAULT_SETTINGSRESETFREQ   -1       /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
#define DEFAULT_SETTINGSRESETOFS    0        /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */

#define MIN_PENALTYPARAM            1e5      /**< if the penalty parameter is to be computed, this is the minimum value it will take */
#define MAX_PENALTYPARAM            1e12     /**< if the penalty parameter is to be computed, this is the maximum value it will take */
#define PENALTYPARAM_FACTOR_DSDP    1e4      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define PENALTYPARAM_FACTOR_SDPA    1e1      /**< if the penalty parameter is to be computed, the maximal objective coefficient will be multiplied by this */
#define MAX_MAXPENALTYPARAM         1e15     /**< if the maximum penaltyparameter is to be computed, this is the maximum value it will take */
#define MAXPENALTYPARAM_FACTOR      1e6      /**< if the maximum penaltyparameter is to be computed, it will be set to penaltyparam * this */
#define MIN_LAMBDASTAR              1e0      /**< if lambda star is to be computed, this is the minimum value it will take */
#define MAX_LAMBDASTAR              1e8     /**< if lambda star is to be computed, this is the maximum value it will take */
#define LAMBDASTAR_FACTOR           1e0      /**< if lambda star is to be computed, the biggest guess of the SDP blocks is multiplied by this value */

#define PRINT_STATISTICS /* uncomment this to print additional statistics after the computation is finished */

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_SDPI*            sdpi;               /**< general SDP Interface that is given the data to presolve the SDP and give it so a solver specific interface */
   SdpVarmapper*         varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */
   SCIP_Real             objval;             /**< objective value of the last SDP relaxation */
   SCIP_Bool             origsolved;         /**< solved original problem to optimality (not only a penalty or probing formulation) */
   SCIP_Bool             probingsolved;      /**< was the last probing SDP solved successfully? */
   SCIP_Real             sdpsolverepsilon;   /**< the stopping criterion for the duality gap the sdpsolver should use */
   SCIP_Real             sdpsolverfeastol;   /**< the feasibility tolerance the SDP solver should use for the SDP constraints */
   SCIP_Real             penaltyparam;       /**< the starting penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   SCIP_Real             maxpenaltyparam;    /**< the maximum penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge */
   SCIP_Real             lambdastar;         /**< the parameter lambda star used by SDPA to set the initial point */
   int                   sdpiterations;      /**< saves the total number of sdp-iterations */
   int                   solvedfast;         /**< number of SDPs solved with fast settings */
   int                   solvedmedium;       /**< number of SDPs solved with medium settings */
   int                   solvedstable;       /**< number of SDPs solved with stable settings */
   int                   solvedpenalty;      /**< number of SDPs solved using penalty formulation */
#if 0
   int                   threads;            /**< number of threads used for SDP solving */
#endif
   SCIP_Bool             slatercheck;        /**< Should the Slater condition for the dual problem be check ahead of solving every SDP ? */
   SCIP_Bool             sdpinfo;            /**< Should the SDP solver output information to the screen? */
   SCIP_Bool             objlimit;           /**< Should an objective limit be given to the SDP solver? */
   SCIP_Bool             resolve;            /**< Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing) ? */
   SCIP_Bool             tightenvb;          /**< Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ? */
   int                   settingsresetfreq;  /**< frequency for resetting parameters in SDP solver and trying again with fastest settings */
   int                   settingsresetofs;   /**< frequency offset for resetting parameters in SDP solver and trying again with fastest settings */
   int                   sdpcalls;           /**< number of solved SDPs (used to compute average SDP iterations) */
   long int              lastsdpnode;        /**< number of the SCIP node the current SDP-solution belongs to */
   SCIP_Bool             feasible;           /**< was the last solved SDP feasible */
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
            sdpvar[ind][j] = SCIPsdpVarmapperGetSdpIndex(varmapper, blockvars[j]);

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

/** inserts all the LP data (including bounds and objective) into the corresponding SDP Interface */
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
   SCIP_Real* obj;
   int* objinds;
   SCIP_Bool tightenvb;

   assert( scip != NULL );
   assert( sdpi != NULL );
   assert( varmapper != NULL );

   nvars = SCIPgetNVars(scip);
   assert( nvars > 0 );

   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/tightenvb", &tightenvb) );

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
      SCIP_Bool tightened;
      SCIP_Bool swapped;
      SCIP_Real tightenedval;

      row = rows[i];
      assert( row != 0 );
      rownnonz = SCIProwGetNNonz(row);
      tightened = FALSE;

      rowvals = SCIProwGetVals(row);
      rowcols = SCIProwGetCols(row);
      sciplhs = SCIProwGetLhs(row) - SCIProwGetConstant(row);
      sciprhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);

      /* check whether we have a variable bound and can strenghten the big-M */
      if ( tightenvb && rownnonz == 2 && (SCIPisZero(scip, sciplhs) || SCIPisZero(scip, sciprhs) ) )
      {
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         SCIP_Real val1;
         SCIP_Real val2;

         val1 = rowvals[0];
         val2 = rowvals[1];

         assert( rowcols[0] != NULL );
         assert( rowcols[1] != NULL );
         var1 = SCIPcolGetVar(rowcols[0]);
         var2 = SCIPcolGetVar(rowcols[1]);
         assert( var1 != NULL );
         assert( var2 != NULL );

         /* check that variables are not locally fixed */
         if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
         {
            /* one coefficient must be 1 and the other negative */
            if ( (SCIPisEQ(scip, val1, 1.0) || SCIPisEQ(scip, val2, 1.0)) && ( SCIPisNegative(scip, val1) || SCIPisNegative(scip, val2) ) )
            {
               swapped = FALSE;
               /* We want x - a z <= 0 or x - a z >= 0, where var1 = x and var2 = z; possibly swap variables otherwise */
               if ( ! SCIPisEQ(scip, val1, 1.0) || ! SCIPisNegative(scip, val2) )
               {
                  SCIP_Real aval;

                  SCIPswapPointers((void**) &var1, (void**) &var2);

                  aval = val1;
                  val1 = val2;
                  val2 = aval;
                  swapped = TRUE;
               }

               /* var2 needs to be binary */
               if ( SCIPvarIsBinary(var2) )
               {
                  if ( SCIPisZero(scip, sciprhs) )
                  {
                     if ( SCIPisLT(scip, SCIPvarGetUbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetUbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }

                  if ( SCIPisZero(scip, sciplhs) )
                  {
                     if ( SCIPisGT(scip, SCIPvarGetLbLocal(var1), REALABS(val2)) )
                     {
                        SCIPdebugMessage("Big-M in %s changed from %f to %f\n", SCIProwGetName(row), REALABS(val2), SCIPvarGetLbLocal(var1));

                        tightened = TRUE;
                        tightenedval = -SCIPvarGetUbLocal(var1); /* negative sign because the coefficient needs to be negative */
                     }
                  }
               }
            }
         }
      }

      for (j = 0; j < rownnonz; j++)
      {
         /* if the Big-M was tightened, we use the new value (the position where this new value is used is dependant on wheter we needed to swap) */
         if ( tightened && ( (swapped && (j == 0)) || ((! swapped) && (j == 1)) ) ) /* use the tightened value */
         {
            if ( SCIPisFeasGT(scip, REALABS(tightenedval), 0.0) )
            {
               assert( SCIPcolGetVar(rowcols[j]) != 0 );
               colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
               rowind[nnonz] = nconss;
               val[nnonz] = tightenedval;
               nnonz++;
            }
         }
         else if ( SCIPisFeasGT(scip, REALABS(rowvals[j]), 0.0))
         {
            assert( SCIPcolGetVar(rowcols[j]) != 0 );
            colind[nnonz] = SCIPsdpVarmapperGetSdpIndex(varmapper, SCIPcolGetVar(rowcols[j]));
            rowind[nnonz] = nconss;
            val[nnonz] = rowvals[j];
            nnonz++;
         }
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
   vars = SCIPgetVars(scip);
   assert( vars != NULL );

   /* prepare arrays of bounds */
   SCIP_CALL( SCIPallocBufferArray(scip, &lb, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ub, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &obj, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objinds, nvars) );

   /* get new bounds and objective coefficients */
   for (i = 0; i < nvars; i++)
   {
      assert( vars[i] != NULL );
      lb[i] = SCIPvarGetLbLocal(vars[i]);
      ub[i] = SCIPvarGetUbLocal(vars[i]);
      inds[i] = i; /* we want to change all bounds, so all indices are included in inds */
      obj[i] = SCIPvarGetObj(vars[i]);
      objinds[i] = i;
   }

   /* inform interface */
   SCIP_CALL( SCIPsdpiChgBounds(sdpi, nvars, inds, lb, ub) );
   SCIP_CALL( SCIPsdpiChgObj(sdpi, nvars, objinds, obj) );

   /* free the bounds-arrays */
   SCIPfreeBufferArray(scip, &objinds);
   SCIPfreeBufferArray(scip, &obj);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &ub);
   SCIPfreeBufferArray(scip, &lb);

   return SCIP_OKAY;
}

/** calculate relaxation and process the relaxation results
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
   SCIP_Bool rootnode;
   int naddediters;
   SCIP_SDPSOLVERSETTING startsetting;
   SCIP_SDPSOLVERSETTING usedsetting;
   SCIP_CONS* savedsetting;
   char* saveconsname;
   SCIP_CONS** conss;
#ifndef NDEBUG
   int snprintfreturn; /* this is used to assert that the SCIP string concatenation works */
#endif
#ifdef SCIP_MORE_DEBUG
   SCIP_Real objforscip;
   SCIP_Real* solforscip;
   SCIP_Bool allint;
   int sollength;
#endif

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

   if ( relaxdata->objlimit )
   {
      /* set the objective limit */
      assert( SCIPgetUpperbound(scip) > -SCIPsdpiInfinity(sdpi) );
      SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_OBJLIMIT, SCIPgetUpperbound(scip)) );
   }

   /* if this is the root node and we cannot solve the problem, we want to check for the Slater condition independent of the SCIP parameter */
   rootnode = ! SCIPnodeGetParent(SCIPgetCurrentNode(scip));

   /* find settings to use for this relaxation */
   if ( rootnode || (SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == relaxdata->settingsresetofs) ||
         (relaxdata->settingsresetfreq > 0 && ((SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) - relaxdata->settingsresetofs) % relaxdata->settingsresetfreq == 0)) )
      startsetting = SCIP_SDPSOLVERSETTING_UNSOLVED; /* in the root node we have no information, at each multiple of resetfreq we reset */
   else
   {
      SCIP_CONSHDLR* conshdlr;
      int lastconsind;

      /* get constraint handler */
      conshdlr = SCIPfindConshdlr(scip, "Savedsdpsettings");
      if ( conshdlr == NULL )
      {
         SCIPerrorMessage("Savedsdpsettings constraint handler not found!\n");
         return SCIP_PLUGINNOTFOUND;
      }

      /* get constraints */
      conss = SCIPconshdlrGetConss(conshdlr);
      lastconsind = SCIPconshdlrGetNConss(conshdlr) - 1;


      assert ( conss != NULL );
      assert ( conss[lastconsind] != NULL ); /* we always use the last information we got (important e.g. in fracdiving) */

      /* start with the settings of the parentnode */
      startsetting = SCIPconsSavedsdpsettingsGetSettings(scip, conss[lastconsind]);

   }

   /* solve the problem */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, startsetting, rootnode) );
   relaxdata->lastsdpnode = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));

   /* update calls, iterations and stability numbers */
   relaxdata->sdpcalls++;
   naddediters = 0;
   SCIP_CALL( SCIPsdpiGetIterations(relaxdata->sdpi, &naddediters) );
   relaxdata->sdpiterations += naddediters;
   usedsetting = SCIP_SDPSOLVERSETTING_UNSOLVED;
   SCIP_CALL( SCIPsdpiSettingsUsed(relaxdata->sdpi, &usedsetting) );
   switch( usedsetting )/*lint --e{788}*/
   {
      case SCIP_SDPSOLVERSETTING_PENALTY:
         relaxdata->solvedpenalty++;
         break;
      case SCIP_SDPSOLVERSETTING_FAST:
         relaxdata->solvedfast++;
         break;
      case SCIP_SDPSOLVERSETTING_MEDIUM:
         relaxdata->solvedmedium++;
         break;
      case SCIP_SDPSOLVERSETTING_STABLE:
         relaxdata->solvedstable++;
         break;
      default:
         break;
   }

   /* remember settings */
   SCIP_CALL( SCIPallocBufferArray(scip, &saveconsname, 255));
#ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(saveconsname, 255, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
         assert( snprintfreturn < 256 ); /* this is the number of positions needed, we gave 255 */
#else
         SCIPsnprintf(saveconsname, 255, "savedsettings_node_%d", SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));
#endif
   SCIP_CALL( createConsSavedsdpsettings(scip, &savedsetting, saveconsname, usedsetting) );
   SCIP_CALL( SCIPaddCons(scip, savedsetting) );
   SCIP_CALL( SCIPreleaseCons(scip, &savedsetting) );
   SCIPfreeBufferArray(scip, &saveconsname);

   if ( SCIPsdpiWasSolved(sdpi) && SCIPsdpiSolvedOrig(sdpi) )
   {
      if ( SCIPinProbing(scip) )
         relaxdata->probingsolved = TRUE;
      else
         relaxdata->origsolved = TRUE;
   }
   else if ( ! SCIPsdpiWasSolved(sdpi) )
   {
      SCIP_Real objlb;
      SCIP_NODE* node;

      /* We couldn't solve the problem, not even with a penalty formulation, so we reuse the relaxation result of the parent node (if one exists) */
      node = SCIPnodeGetParent(SCIPgetCurrentNode(scip));

      relaxdata->origsolved = FALSE;
      if ( SCIPinProbing(scip) )
         relaxdata->probingsolved = FALSE;

      if ( node == NULL )
      {
         relaxdata->feasible = FALSE;
         *result = SCIP_SUSPENDED;
         SCIPerrorMessage("The relaxation of the root node could not be solved, so there is no hope to solve this instance. \n");
         return SCIP_ERROR;
      }

      relaxdata->feasible = FALSE;

      /* if we used the penalty approach, we might have calculated a good lower bound, even if we did not produce a feasible solution */
      objlb = -SCIPinfinity(scip);
      SCIP_CALL( SCIPsdpiGetLowerObjbound(relaxdata->sdpi, &objlb) );
      if ( ! SCIPisInfinity(scip, objlb) )
         *lowerbound = objlb;
      else
         *lowerbound = SCIPnodeGetLowerbound(node);

      *result = SCIP_SUCCESS;
      SCIP_CALL( SCIPupdateLocalLowerbound(scip, *lowerbound) );
      SCIPdebugMessage("The relaxation couldn't be solved, so the relaxation result from the parent node was copied. \n");
      return SCIP_OKAY;
   }

#ifdef SCIP_MORE_DEBUG /* print the optimal solution */

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
         SCIPdebugMessage("Node cut off due to infeasibility.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsObjlimExc(sdpi) )
      {
         SCIPdebugMessage("Node cut off due to objective limit.\n");
         relaxdata->feasible = FALSE;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else if ( SCIPsdpiIsDualUnbounded(sdpi) )
      {
         SCIPdebugMessage("Node unbounded.");
         relaxdata->feasible = TRUE;
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
            if ( SCIPvarIsIntegral(SCIPsdpVarmapperGetSCIPvar(varmapper, v)) && ! (SCIPisFeasIntegral(scip, solforscip[v])) )
            {
               allint = FALSE;
               break;
            }
         }

         /* create SCIP solution */
         SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
         SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

         *lowerbound = objforscip;
         relaxdata->objval = objforscip;

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

               /* set relax sol */
               SCIP_CALL( SCIPsetRelaxSolVals(scip, nvars, vars, solforscip) );
               SCIP_CALL( SCIPmarkRelaxSolValid(scip) );

               SCIPfreeBufferArray(scip, &solforscip);
               SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

               relaxdata->feasible = TRUE;
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            SCIPdebugMessage("Found a solution that is feasible for the SDP-solver and integrality, but infeasible for SCIP! \n");
         }

         /* copy solution */
         SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
         for (i = 0; i < ncols; i++)
            SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );

         SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;

         /* if all int and binary vars are integral, nothing else needs to be done */
         if ( ! allint )
         {
            for (i = 0; i < nvars; ++i)
            {
               SCIP_VAR* var = vars[i];
               if ( SCIPvarIsIntegral(var) && ! SCIPisFeasIntegral(scip, solforscip[i]) && ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
               {
                  /* we don't set a true score, we will just let the branching rule decide */
                  SCIP_CALL( SCIPaddExternBranchCand(scip, var, 10000.0, solforscip[i]) );
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
      if ( SCIPisFeasLT(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) )
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
   SCIP_SOL* scipsol;
   SCIP_Bool stored;
#ifdef SCIP_EVEN_MORE_DEBUG
   SCIP_VAR** varsfordebug = SCIPgetVars(scip);
   const int nvarsfordebug = SCIPgetNVars(scip);
#endif

   SCIPdebugMessage("Calling relaxExecSdp.\n");

   relaxdata = SCIPrelaxGetData(relax);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* don't run again if we already solved the current node (except during probing), and we solved the correct problem */
   if ( ( relaxdata->lastsdpnode == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) && ( ! SCIPinProbing(scip) ) )
         && relaxdata->origsolved && ( ! relaxdata->resolve) )
   {
      SCIP_COL** cols;
      int ncols;
      int slength;
      SCIP_Real objforscip;
      SCIP_Real* solforscip;

      SCIPdebugMessage("Already solved SDP-relaxation for node %ld, returning with SCIP_SUCCESS so that no other relaxator is called.\n",
            SCIPrelaxGetData(relax)->lastsdpnode);
      if ( SCIPsdpiIsDualUnbounded(relaxdata->sdpi) )
      {
         relaxdata->feasible = TRUE;
         *result = SCIP_SUCCESS;
         *lowerbound = -SCIPinfinity(scip);
         return SCIP_OKAY;
      }
      /* get solution w.r.t. SCIP variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &solforscip, nvars) );
      slength = nvars;
      SCIP_CALL( SCIPsdpiGetSol(relaxdata->sdpi, &objforscip, solforscip, &slength) ); /* get both the objective and the solution from the SDP solver */

      assert( slength == nvars ); /* if this isn't true any longer, the getSol-Call was unsuccessfull, because the given array wasn't long enough,
                                   * but this can't happen, because the array has enough space for all sdp variables */

      /* create SCIP solution */
      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, solforscip) );

      *lowerbound = objforscip;

      /* copy solution */
      SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
      for (i = 0; i < ncols; i++)
            SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPcolGetVar(cols[i]), SCIPgetSolVal(scip, scipsol, SCIPcolGetVar(cols[i]))) );

      SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
      *result = SCIP_SUCCESS;

      SCIPfreeBufferArray(scip, &solforscip);
      SCIP_CALL( SCIPfreeSol(scip, &scipsol) );

      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* if we are solving a probing SDP, remember that we didn't solve the original problem */
   relaxdata->origsolved = FALSE;

   /* construct the lp and make sure, that everything is where it should be */
   SCIP_CALL( SCIPconstructLP(scip, &cutoff) );

   if ( cutoff )
   {
      relaxdata->feasible = FALSE;
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
      relaxdata->feasible = TRUE;
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( allVarsFixed(scip) )
   {
      SCIP_Bool feasible;

      /* if all variables, really all, are fixed, we give this fixed solution to SCIP */

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &ubs, nvars) );

      *lowerbound = 0.0;
      for (i = 0; i < nvars; i++)
      {
         ubs[i] = SCIPvarGetUbLocal(vars[i]);
         *lowerbound += SCIPvarGetObj(vars[i]) * ubs[i];
         assert( SCIPisFeasEQ(scip, SCIPvarGetUbLocal(vars[i]), SCIPvarGetLbLocal(vars[i])));
      }
      if ( SCIPgetObjsense(scip) == -1 ) /*lint !e641*/
         *lowerbound *= -1;

      SCIPdebugMessage("EVERYTHING IS FIXED, objective value = %f\n", *lowerbound);

      SCIP_CALL( SCIPcreateSol(scip, &scipsol, NULL) );
      SCIP_CALL( SCIPsetSolVals(scip, scipsol, nvars, vars, ubs) );

      /* set the relaxation solution */
      for (i = 0; i < nvars; i++)
         SCIP_CALL( SCIPsetRelaxSolVal(scip, vars[i], SCIPvarGetLbLocal(vars[i])) );
      SCIP_CALL( SCIPmarkRelaxSolValid(scip) );

      /* check if the solution really is feasible */
      SCIP_CALL( SCIPcheckSol(scip, scipsol, FALSE, TRUE, TRUE, TRUE, &feasible) );

      stored = FALSE;
      if ( feasible )
      {
         SCIP_CALL( SCIPtrySolFree(scip, &scipsol, FALSE, FALSE, FALSE, FALSE, &stored) );
      }
      else
      {
         SCIP_CALL( SCIPfreeSol(scip, &scipsol) );
      }

      relaxdata->feasible = feasible;

      if (feasible && stored == 1)
      {
         *result = SCIP_CUTOFF;
         SCIPdebugMessage("New solution was stored, node is cut off !\n");
      }
      else
      {
         *result = SCIP_CUTOFF;
         SCIPdebugMessage("Fixed solution either infeasible or not good enough for storage, node cut off !\n");
      }

      SCIPfreeBlockMemoryArray(scip, &ubs, nvars);

      return SCIP_OKAY;
   }

   /* update LP Data in Interface */
   SCIP_CALL( putLpDataInInterface(scip, relaxdata->sdpi, relaxdata->varmapper) );

   SCIP_CALL( calc_relax(scip, relaxdata, result, lowerbound));

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
   SCIP_Real penaltyparam;
   SCIP_Real maxpenaltyparam;
#if 0
   int threads;
#endif
   SCIP_Bool sdpinfo;
   SCIP_Bool slatercheck;
   SCIP_Real givenpenaltyparam;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->sdpcalls = 0;
   relaxdata->sdpiterations = 0;
   relaxdata->solvedfast = 0;
   relaxdata->solvedmedium = 0;
   relaxdata->solvedstable = 0;
   relaxdata->solvedpenalty = 0;
   relaxdata->feasible = FALSE;

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   /* all SCIPvars will be added to this list, and 3/4 seems like a good load factor (java uses this factor) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(relaxdata->varmapper), (int) ceil(1.33 * nvars)) );
   SCIP_CALL( SCIPsdpVarmapperAddVars(scip, relaxdata->varmapper, nvars, vars) );

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

   /* set/compute the starting penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/penaltyparam", &penaltyparam) );
   if ( SCIPisGE(scip, penaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, penaltyparam);
      givenpenaltyparam = penaltyparam;
   }
   else
   {
      SCIP_Real maxcoeff;
      int v;
      SCIP_Real compval;

      /* we set the value to min{max{MIN_PENALTYPARAM, PENALTYPARAM_FACTOR * max_objective_coefficient}, MAX_PENALTYPARAM} */

      /* compute the maximum coefficient in the objective */
      maxcoeff = 0.0;
      for (v = 0; v < nvars; v++)
      {
         if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(vars[v])), maxcoeff) )
            maxcoeff = REALABS(SCIPvarGetObj(vars[v]));
      }

      /* compute the value we would like to set the penaltyparameter to */
      if ( strcmp(SCIPsdpiGetSolverName(), "DSDP") == 0 )
         compval = PENALTYPARAM_FACTOR_DSDP * maxcoeff;
      else if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0 )
         compval = PENALTYPARAM_FACTOR_SDPA * maxcoeff;
      else
      {
         SCIPdebugMessage("unknown SDP-Solver %s when setting penaltyparam !\n", SCIPsdpiGetSolverName());
         compval = SCIPinfinity(scip);
      }

      if ( SCIPisLT(scip, compval, MIN_PENALTYPARAM) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, MIN_PENALTYPARAM);
         SCIPdebugMessage("Setting penaltyparameter to %f.\n", MIN_PENALTYPARAM);
         givenpenaltyparam = MIN_PENALTYPARAM;
      }
      else if ( SCIPisGT(scip, compval, MAX_PENALTYPARAM) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, MAX_PENALTYPARAM);
         SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_PENALTYPARAM);
         givenpenaltyparam = MAX_PENALTYPARAM;
      }
      else
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, compval);
         SCIPdebugMessage("Setting penaltyparameter to %f.\n", compval);
         givenpenaltyparam = compval;
      }
   }
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: penaltyparam setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* set/compute the maximum penalty parameter */
   SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/maxpenaltyparam", &maxpenaltyparam) );
   if ( SCIPisGE(scip, maxpenaltyparam, 0.0) )
   {
      retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_MAXPENALTYPARAM, maxpenaltyparam);

      /* check if the starting value is not bigger than the maximum one, otherwise update it */
      if ( SCIPisLT(scip, givenpenaltyparam, maxpenaltyparam) )
      {
         SCIPdebugMessage("Penalty parameter %f overwritten by maxpenaltyparam %f! \n", givenpenaltyparam, maxpenaltyparam);
         SCIP_CALL( SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_PENALTYPARAM, maxpenaltyparam) );

      }
   }
   else
   {
      SCIP_Real compval;

      /* we set the value to min{MAX_MAXPENALTYPARAM, MAXPENALTYPARAM_FACTOR * penaltyparam} */

      assert( SCIPisLE(scip, givenpenaltyparam, MAX_MAXPENALTYPARAM) );

      /* compute the value we would like to set the penaltyparameter to */
      compval = givenpenaltyparam * MAXPENALTYPARAM_FACTOR;

      if ( SCIPisLT(scip, compval, MAX_MAXPENALTYPARAM) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_MAXPENALTYPARAM, compval);
         SCIPdebugMessage("Setting maximum penaltyparameter to %f.\n", compval);
      }
      else
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_MAXPENALTYPARAM, MAX_MAXPENALTYPARAM);
         SCIPdebugMessage("Setting penaltyparameter to %f.\n", MAX_MAXPENALTYPARAM);
      }
   }
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: maxpenaltyparam setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   /* set/compute lambda star if SDPA is used as the SDP-Solver */
   if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0.0 )
   {
      SCIP_Real lambdastar;

      SCIP_CALL( SCIPgetRealParam(scip, "relaxing/SDP/lambdastar", &lambdastar) );
      if ( SCIPisGE(scip, lambdastar, 0.0) )
      {
         retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, lambdastar);
      }
      else
      {
         SCIP_Real guess;
         SCIP_Real maxguess;
         SCIP_CONS** conss;
         int nconss;
         int c;
         SCIP_Real compval;

         /* we set the value to min{max{MIN_LAMBDASTAR, LAMBDASTAR_FACTOR * MAX_GUESS}, MAX_LAMBDASTAR}, where MAX_GUESS is the maximum of the guesses
          * of the SDP-Blocks */

         /* compute the maximum guess */
         conss = SCIPgetConss(scip);
         nconss = SCIPgetNConss(scip);
         maxguess = 0.0;

         for (c = 0; c < nconss; c++)
         {
            /* only check the SDP constraints */
            if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
            {
               SCIP_CALL( SCIPconsSdpGuessInitialPoint(scip, conss[c], &guess) );
               if ( SCIPisGT(scip, guess, maxguess) )
               {
                  maxguess = guess;
               }
            }
         }

         compval = LAMBDASTAR_FACTOR * maxguess;

         if ( SCIPisLT(scip, compval, MIN_LAMBDASTAR) )
         {
            retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, MIN_LAMBDASTAR);
            SCIPdebugMessage("Setting lambdastar to %f.\n", MIN_LAMBDASTAR);
         }
         else if ( SCIPisGT(scip, compval, MAX_LAMBDASTAR) )
         {
            retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, MAX_LAMBDASTAR);
            SCIPdebugMessage("Setting lambdastar to %f.\n", MAX_LAMBDASTAR);
         }
         else
         {
            retcode = SCIPsdpiSetRealpar(relaxdata->sdpi, SCIP_SDPPAR_LAMBDASTAR, compval);
            SCIPdebugMessage("Setting lambdastar to %f.\n", compval);
         }

      }
   }
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: lambdastar setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

#if 0
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
#endif

   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/sdpinfo", &sdpinfo) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SDPINFO, (int) sdpinfo);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: sdpinfo setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/SDP/slatercheck", &slatercheck) );
   retcode = SCIPsdpiSetIntpar(relaxdata->sdpi, SCIP_SDPPAR_SLATERCHECK, (int) slatercheck);
   if ( retcode == SCIP_PARAMETERUNKNOWN )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
         "SDP Solver <%s>: slatercheck setting not available -- SCIP parameter has no effect\n",
         SCIPsdpiGetSolverName());
   }
   else
   {
      SCIP_CALL( retcode );
   }

   return SCIP_OKAY;
}

/** copy method for sdp relaxation handler (called when SCIP copies plugins) */
static
SCIP_DECL_RELAXCOPY(relaxCopySdp)
{
   assert( scip != NULL );
   assert( relax != NULL );
   assert(strcmp(SCIPrelaxGetName(relax), RELAX_NAME) == 0);

   SCIP_CALL( SCIPincludeRelaxSdp(scip) );

   return SCIP_OKAY;
}

/** reset the relaxator's data */
static
SCIP_DECL_RELAXEXIT(relaxExitSdp)
{
   SCIP_RELAXDATA* relaxdata;

   assert( scip != NULL );
   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert( relaxdata != NULL );

#ifdef PRINT_STATISTICS
   SCIPinfoMessage(scip, NULL, "SDP-Iterations:%d \n", relaxdata->sdpiterations);
   SCIPinfoMessage(scip, NULL, "Average-SDP-Iterations:%f \n", (double) relaxdata->sdpiterations / (double) relaxdata->sdpcalls );
   SCIPinfoMessage(scip, NULL, "Fastest-Settings-Solved:%f \n", 100.0 * (double) relaxdata->solvedfast / (double) relaxdata->sdpcalls);
   SCIPinfoMessage(scip, NULL, "Medium-Settings-Solved:%f \n", 100.0 * (double) relaxdata->solvedmedium / (double) relaxdata->sdpcalls);
   SCIPinfoMessage(scip, NULL, "Stable-Settings-Solved:%f \n", 100.0 * (double) relaxdata->solvedstable / (double) relaxdata->sdpcalls);
   SCIPinfoMessage(scip, NULL, "Penalty-Percent:%f \n", 100.0 * (double) relaxdata->solvedpenalty / (double) relaxdata->sdpcalls);
#endif

   SCIPdebugMessage("Exiting Relaxation Handler \n");

   if ( relaxdata->varmapper != NULL )
   {
      SCIP_CALL( SCIPsdpVarmapperFree(scip, &(relaxdata->varmapper)) );
   }

   relaxdata->objval = 0.0;
   relaxdata->origsolved = FALSE;
   relaxdata->probingsolved = FALSE;
   relaxdata->feasible = FALSE;
   relaxdata->sdpiterations = 0;
   relaxdata->sdpcalls = 0;
   relaxdata->lastsdpnode = 0;
   SCIP_CALL( SCIPsdpiClear(relaxdata->sdpi) );

   return SCIP_OKAY;
}

/** free the relaxator's data */
static
SCIP_DECL_RELAXFREE(relaxFreeSdp)
{/*lint --e{715}*/
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
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/penaltyparam", "the starting value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge, set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->penaltyparam),
         TRUE, DEFAULT_PENALTYPARAM, -1.0, 1e+20, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/maxpenaltyparam", "the maximum value of the penalty parameter Gamma used for the penalty formulation if the "
         "SDP solver didn't converge, set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->maxpenaltyparam),
         TRUE, DEFAULT_MAXPENALTYPARAM, -1.0, 1e+20, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "relaxing/SDP/lambdastar", "the parameter lambda star used by SDPA to set the initial point , "
         "set this to a negative value to compute the parameter depending on the given problem", &(relaxdata->lambdastar),
            TRUE, DEFAULT_LAMBDASTAR, -1.0, 1e+20, NULL, NULL) );
#if 0
   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/threads", "number of threads used for SDP solving",
         &(relaxdata->threads), TRUE, DEFAULT_THREADS, 1, INT_MAX, NULL, NULL) );
#endif
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/slatercheck", "should the Slater condition for the dual problem be check ahead of solving each SDP?",
         &(relaxdata->slatercheck), TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/sdpinfo", "should the SDP solver output information to the screen?",
         &(relaxdata->sdpinfo), TRUE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/objlimit", "should an objective limit be given to the SDP-Solver?",
         &(relaxdata->objlimit), TRUE, DEFAULT_OBJLIMIT, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/resolve", "Are we allowed to solve the relaxation of a single node multiple times in a row"
         " (outside of probing) ?", &(relaxdata->resolve), TRUE, DEFAULT_RESOLVE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "relaxing/SDP/tightenvb", "Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ?",
         &(relaxdata->tightenvb), TRUE, DEFAULT_TIGHTENVB, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetfreq",
         "frequency for resetting parameters in SDP solver and trying again with fastest settings (-1: never, 0: only at depth settingsresetofs)",
         &(relaxdata->settingsresetfreq), TRUE, DEFAULT_SETTINGSRESETFREQ, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "relaxing/SDP/settingsresetofs",
         "frequency offset for resetting parameters in SDP solver and trying again with fastest settings",
         &(relaxdata->settingsresetofs), TRUE, DEFAULT_SETTINGSRESETOFS, 0, INT_MAX, NULL, NULL) );

   /* add description of SDP-solver */
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPsdpiGetSolverName(), SCIPsdpiGetSolverDesc()) );

   return SCIP_OKAY;
}


/* external functions */

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note if a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent
 *  primal variable.
 */
SCIP_RETCODE SCIPrelaxSdpGetPrimalBoundVars(
   SCIP_RELAX*           relax,              /**< SDP relaxator to information for */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
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
   SCIP_Real*            objval              /**< pointer to store the optimal objective value of the SDP relaxation */
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
   SCIP_Real*            solarray,           /**< pointer to store the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray */
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

/** was the original problem solved for the last SDP-Node (or a penalty or probing formulation) ? */
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->origsolved && SCIPsdpiSolvedOrig(relaxdata->sdpi);
}

/** was the last probing SDP solved successfully ? */
SCIP_Bool SCIPrelaxSdpSolvedProbing(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   )
{
   SCIP_RELAXDATA* relaxdata;

   assert( relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);

   assert( relaxdata != NULL );
   assert( relaxdata->sdpi != NULL );

   return relaxdata->probingsolved;
}

/** returns whether the last solved problem was feasible */
SCIP_Bool SCIPrelaxSdpIsFeasible(
   SCIP_RELAX*           relax               /**< SDP relaxator to get feasibility for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->feasible );
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

/** returns number of SDP relaxation solved with fast settings */
int SCIPrelaxSdpGetNSdpFast(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedfast );
}

/** returns number of SDP relaxation solved with medium settings */
int SCIPrelaxSdpGetNSdpMedium(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedmedium );
}

/** returns number of SDP relaxation solved with stable settings */
int SCIPrelaxSdpGetNSdpStable(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedstable );
}

/** returns number of SDP relaxation solved with penalty formulation */
int SCIPrelaxSdpGetNSdpPenalty(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the number of calls for */
   )
{
   assert( relax != NULL );
   assert( SCIPrelaxGetData(relax) != NULL );

   return ( SCIPrelaxGetData(relax)->solvedpenalty );
}
