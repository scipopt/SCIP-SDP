/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpsymmetry.c
 * @brief  routines for handling/detecting symmetries in SDPs
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scipsdp/sdpsymmetry.h"
#include "scipsdp/cons_sdp.h"

/** sorts real numbers
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortReal)
{
   SCIP_Real diffvals;
   SCIP_Real* vals;

   vals = (SCIP_Real*) dataptr;
   diffvals = vals[ind1] - vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
}

/** sorts integer numbers
 *
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
static
SCIP_DECL_SORTINDCOMP(SYMsortInt)
{
   SCIP_Real diffvals;
   int* vals;

   vals = (int*) dataptr;
   diffvals = vals[ind1] - vals[ind2];

   if ( diffvals < 0.0 )
      return -1;
   else if ( diffvals > 0.0 )
      return 1;

   return 0;
}

/** stores information about SDP constraints */
SCIP_RETCODE storeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata             /**< pointer to store SDP symmetry data */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   SCIP_Real** sdpval;
   int* sdpnvarnonz;
   int maxsdpnnonz;
   int maxsdpconstnnonz;
   int sdpconstnnonz;
   int sdpnnonz;
   int nconss;
   int nvars;
   int c;

   int** sdpcol;
   int** sdprow;
   SCIP_VAR** sdpvars;
   int* sdpconstcol;
   int* sdpconstrow;
   SCIP_Real* sdpconstval;

   assert( scip != NULL );
   assert( sdpdata != NULL );

   /* (partially) intialize SDP data */
   sdpdata->nsdpconss = 0;

   /* find SDP conshdlr */
   conshdlr = SCIPfindConshdlr(scip, "SDP");
   assert( conshdlr != NULL );

   nconss = SCIPconshdlrGetNConss(conshdlr);
   if ( nconss == 0 )
      return SCIP_OKAY;

   conss = SCIPconshdlrGetConss(conshdlr);
   assert( conss != NULL );

   nvars = SCIPgetNVars(scip);

   /* initialize SDP data */
   sdpdata->nsdpconss = nconss;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->blocksizes), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->nvars), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->vals), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->vars), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->valsbegins), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->cols), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->rows), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->colors), nconss) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->colors2), nconss) );

   /* temporary memory for copying constraints */
   maxsdpnnonz = 0;
   maxsdpconstnnonz = 0;
   for (c = 0; c < nconss; c++)
   {
      assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 ||
         strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0 );

      SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz) );

      if ( sdpnnonz > maxsdpnnonz )
         maxsdpnnonz = sdpnnonz;
      if ( sdpconstnnonz > maxsdpconstnnonz )
         maxsdpconstnnonz = sdpconstnnonz;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &sdpnvarnonz, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpcol, maxsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdprow, maxsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpval, maxsdpnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstcol, maxsdpconstnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstrow, maxsdpconstnnonz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstval, maxsdpconstnnonz) );

   /* fill data for each constraint */
   for (c = 0; c < nconss; ++c)
   {
      int sdpnvars;
      int sdpblocksize;
      int sdparraylength;
      int totalnnonz;
      int pos;
      int v;
      int i;

      sdparraylength = maxsdpnnonz;
      sdpconstnnonz = maxsdpconstnnonz;

      /* collect information from constraint */
      SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength,
            sdpnvarnonz, sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow, sdpconstval,
            NULL, NULL, NULL) );
      assert( sdparraylength <= maxsdpnnonz );
      assert( sdpconstnnonz <= maxsdpconstnnonz );

      sdpdata->blocksizes[c] = sdpblocksize;
      sdpdata->nvars[c] = sdpnvars;

      totalnnonz = sdpnnonz + sdpconstnnonz;

      /* allocate memory */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->valsbegins[c]), sdpnvars + 2) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->vals[c]), totalnnonz) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->vars[c]), sdpnvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->cols[c]), totalnnonz) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->rows[c]), totalnnonz) );

      /* copy data for variable matrices */
      pos = 0;
      for (v = 0; v < sdpnvars; ++v)
      {
         assert( SCIPvarGetProbindex(sdpvars[v]) >= 0 ); /* adapt graph construction if vars can be aggregated */
         sdpdata->valsbegins[c][v] = pos;
         sdpdata->vars[c][v] = sdpvars[v];
         for (i = 0; i < sdpnvarnonz[v]; ++i)
         {
            sdpdata->vals[c][pos] = sdpval[v][i];
            sdpdata->cols[c][pos] = sdpcol[v][i];
            sdpdata->rows[c][pos++] = sdprow[v][i];
         }
      }

      /* copy data for constant matrix */
      sdpdata->valsbegins[c][sdpnvars] = pos;
      for (i = 0; i < sdpconstnnonz; ++i)
      {
         sdpdata->vals[c][pos] = sdpconstval[i];
         sdpdata->cols[c][pos] = sdpconstcol[i];
         sdpdata->rows[c][pos++] = sdpconstrow[i];
      }
      sdpdata->valsbegins[c][sdpnvars + 1] = pos;

      /* allocate memory for colors */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->colors[c]), totalnnonz) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpdata->colors2[c]), totalnnonz) );
   }

   SCIPfreeBufferArray(scip, &sdpconstval);
   SCIPfreeBufferArray(scip, &sdpconstrow);
   SCIPfreeBufferArray(scip, &sdpconstcol);
   SCIPfreeBufferArray(scip, &sdpvars);
   SCIPfreeBufferArray(scip, &sdpval);
   SCIPfreeBufferArray(scip, &sdprow);
   SCIPfreeBufferArray(scip, &sdpcol);
   SCIPfreeBufferArray(scip, &sdpnvarnonz);

   return SCIP_OKAY;
}

/** frees information about SDP constraints */
SCIP_RETCODE freeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata             /**< pointer to store SDP symmetry data */
   )
{
   int c;

   assert( scip != NULL );
   assert( sdpdata != NULL );

   /* if there are no SDP constraints, there is nothing to be done */
   if ( sdpdata->nsdpconss == 0 )
      return SCIP_OKAY;

   for (c = 0; c < sdpdata->nsdpconss; ++c)
   {
      int sdpnvars;
      int totalnnonz;

      sdpnvars = sdpdata->nvars[c];
      totalnnonz = sdpdata->valsbegins[c][sdpdata->nvars[c] + 1] - sdpdata->valsbegins[c][0];

      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->valsbegins[c]), sdpnvars + 2);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->vals[c]), totalnnonz);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->vars[c]), sdpnvars);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->cols[c]), totalnnonz);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->rows[c]), totalnnonz);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->colors[c]), totalnnonz);
      SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->colors2[c]), totalnnonz);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->blocksizes), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->nvars), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->vals), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->vars), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->valsbegins), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->cols), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->rows), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->colors), sdpdata->nsdpconss);
   SCIPfreeBlockMemoryArrayNull(scip, &(sdpdata->colors2), sdpdata->nsdpconss);

   return SCIP_OKAY;
}

/** finds colors for symmetry detection graph */
SCIP_RETCODE findColorsSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata,            /**< pointer to store SDP symmetry data */
   int                   mincolorval         /**< value of smallest color */
   )
{
   int* consperm;
   int* blockperm;
   int oldblocksize = 0;
   int curblocksize;
   int* considx;
   int* valinconsidx;
   SCIP_Real* blockvals;
   int maxnblockvals;
   int nblockvals;
   int nconss;
   int curcolor;
   int cidx;
   int c;
   int v;
   int i;

   assert( scip != NULL );
   assert( sdpdata != NULL );
   assert( mincolorval >= 0 );

   nconss = sdpdata->nsdpconss;
   if ( nconss <= 0 )
      return SCIP_OKAY;

   assert( sdpdata->valsbegins != NULL );
   assert( sdpdata->blocksizes != NULL );
   assert( sdpdata->vals != NULL );
   assert( sdpdata->rows != NULL );
   assert( sdpdata->cols != NULL );
   assert( sdpdata->colors != NULL );
   assert( sdpdata->colors2 != NULL );

   /* sort SDP constraints based on their block size */
   SCIP_CALL( SCIPallocBufferArray(scip, &consperm, nconss) );
   SCIPsort(consperm, SYMsortInt, (void*) sdpdata->blocksizes, nconss);

   /* allocate memory to store all coefficients of SDP constraints of same block size, use block memory
    * since this can become large */
   maxnblockvals = 0;
   for (c = 0; c < nconss; ++c)
      maxnblockvals += sdpdata->valsbegins[c][sdpdata->nvars[c] + 1] - sdpdata->valsbegins[c][0];
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blockvals, maxnblockvals) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &considx, maxnblockvals) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &valinconsidx, maxnblockvals) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &blockperm, maxnblockvals) );

   /* iterate through SDP constraints with same block size, get their coefficients, and assign colors */
   curcolor = mincolorval;
   for (c = 0; c < nconss; ++c)
   {
      cidx = consperm[c];
      curblocksize = sdpdata->blocksizes[cidx];

      /* we have detected a new group of constraints */
      if ( curblocksize > oldblocksize )
      {
         oldblocksize = curblocksize;
         nblockvals = 0;
      }

      /* store coefficients (including constant block) and their relation to the constraints */
      for (v = 0; v <= sdpdata->nvars[cidx]; ++v)
      {
         for (i = sdpdata->valsbegins[cidx][v]; i < sdpdata->valsbegins[cidx][v + 1]; ++i)
         {
            blockvals[nblockvals] = sdpdata->vals[cidx][i];
            considx[nblockvals] = cidx;
            valinconsidx[nblockvals++] = i;
         }
      }

      /* store colors of SDP constraints in case the group of constraints ends here */
      if ( c == nconss - 1 || sdpdata->blocksizes[consperm[c + 1]] > curblocksize )
      {
         /* sort coefficients */
         SCIPsort(blockperm, SYMsortReals, (void*) blockvals, nblockvals);

         /* iterate over coefficients and store their colors */
         ++curcolor;
         sdpdata->colors[considx[blockperm[0]]][valinconsidx[blockperm[0]]] = curcolor;
         for (i = 1; i < nblockvals; ++i)
         {
            /* if we have found a new color */
            if ( SCIPisGT(scip, blockvals[blockperm[i]], blockvals[blockperm[i - 1]]) )
               ++curcolor;
            sdpdata->colors[considx[blockperm[i]]][valinconsidx[blockperm[i]]] = curcolor;
         }

         /* iterate over coefficients and store their second colors */
         ++curcolor;
         sdpdata->colors2[considx[blockperm[0]]][valinconsidx[blockperm[0]]] = curcolor;
         for (i = 1; i < nblockvals; ++i)
         {
            /* if we have found a new color */
            if ( SCIPisGT(scip, blockvals[blockperm[i]], blockvals[blockperm[i - 1]]) )
               ++curcolor;
            sdpdata->colors2[considx[blockperm[i]]][valinconsidx[blockperm[i]]] = curcolor;
         }
      }
   }

   SCIPfreeBlockMemoryArrayNull(scip, &blockperm, maxnblockvals);
   SCIPfreeBlockMemoryArrayNull(scip, &valinconsidx, maxnblockvals);
   SCIPfreeBlockMemoryArrayNull(scip, &considx, maxnblockvals);
   SCIPfreeBlockMemoryArrayNull(scip, &blockvals, maxnblockvals);
   SCIPfreeBufferArray(scip, &consperm);

   printf("Data visualization\n");
   for (c = 0; c < nconss; ++c)
   {
      int nblocks;

      nblocks = sdpdata->blocksizes[c];

      printf("constraint %d (block size %d)\n", c, nblocks);
      for (i = 0; i < sdpdata->nvars[c]; ++i)
      {
         int j;

         printf("\tVar %s:", SCIPvarGetName(sdpdata->vars[c][i]));
         for (j = sdpdata->valsbegins[c][i]; j < sdpdata->valsbegins[c][i+1]; ++j)
            printf(" [(%d,%d) %f -> (%d,%d)]", sdpdata->rows[c][j], sdpdata->cols[c][j], sdpdata->vals[c][j], sdpdata->colors[c][j], sdpdata->colors2[c][j]);
         printf("\n");
      }
      printf("\tConst:");
      for (i = sdpdata->valsbegins[c][sdpdata->nvars[c]]; i < sdpdata->valsbegins[c][sdpdata->nvars[c]+1]; ++i)
         printf(" [(%d,%d) %f -> (%d,%d)]", sdpdata->rows[c][i], sdpdata->cols[c][i], sdpdata->vals[c][i], sdpdata->colors[c][i], sdpdata->colors2[c][i]);
      printf("\n");
   }


   return SCIP_OKAY;
}
