/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpsolchecker.c
 * @brief  checks a given SDP solution for feasibility
 * @author Tristan Gally
 * @author Frederic Matter
 */

/* #define SCIP_DEBUG */
/* #define SCIP_MORE_DEBUG         /\* shows all values of variables and constraint and their violations *\/ */

#include "sdpi/sdpsolchecker.h"
#include "sdpi/lapack_interface.h"
#include "scip/pub_message.h"                /* for debug and error message */

/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY. */
#define BMS_CALL(x)   do                                                                                     \
                      {                                                                                      \
                         if( NULL == (x) )                                                                   \
                         {                                                                                   \
                            SCIPerrorMessage("No memory in function call.\n");                               \
                            return SCIP_NOMEMORY;                                                            \
                         }                                                                                   \
                      }                                                                                      \
                      while( FALSE )

#define INF 1e+16

/** Given a solution, an SDP instance and a feasibility tolerance, checks whether
 *  the smallest eigenvalue is >= -feastol for a given feasibility tolerance.
 */
SCIP_RETCODE SCIPsdpSolcheckerCheck(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constrain tmmatrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            solvector,          /**< values of all variables (including fixed ones) in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
   )
{/*lint --e{818}*/
   int i;
   int j;
   int b;
   int v;
   int ind;

   assert( bufmem != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nsdpblocks >= 0 );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( nsdpblocks == 0 || sdpnblockvars != NULL );
   assert( sdpconstnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpnblockvarnonz != NULL );
   assert( nsdpblocks == 0 || sdpvar != NULL );
   assert( nsdpblocks == 0 || sdprow != NULL );
   assert( nsdpblocks == 0 || sdpcol != NULL );
   assert( nsdpblocks == 0 || sdpval != NULL );
   assert( nsdpblocks == 0 || indchanges != NULL );
   assert( nsdpblocks == 0 || nremovedinds != NULL );
   assert( nsdpblocks == 0 || blockindchanges != NULL );
   assert( nlpcons >= 0 );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lprow != NULL );
   assert( nlpcons == 0 || lpcol != NULL );
   assert( nlpcons == 0 || lpval != NULL );
   assert( solvector != NULL );
   assert( feastol >= 0 );
   assert( infeasible != NULL );

   /* check variable bounds */
   for (i = 0; i < nvars; i++)
   {
      if ( solvector[i] < lb[i] - feastol || solvector[i] > ub[i] + feastol )
      {
         SCIPdebugMessage("solution found infeasible (feastol=%g) for dual variable bounds: x[%d] = %g <|= [%g, %g]\n",
            feastol, i, solvector[i], lb[i], ub[i]);
         *infeasible = TRUE;
         return SCIP_OKAY;
      }
   }

   /* check linear constraints (since DSDP sorts the lp-arrays by cols, we cannot expect them to be sorted) */
   if ( nlpcons > 0 )
   {
      SCIP_Real* lpconsvals;

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &lpconsvals, nlpcons) );

      /* initialize all rows with zero */
      for (i = 0; i < nlpcons; i++)
         lpconsvals[i] = 0;

      /* compute the values of all rows */
      for (i = 0; i < lpnnonz; i++)
      {
         if ( lb[lpcol[i]] < ub[lpcol[i]] - epsilon ) /* fixed variables are already included in lhs/rhs */
            lpconsvals[lprow[i]] += solvector[lpcol[i]] * lpval[i];
      }

      /* check all active constraints for feasibility */
      ind = 0; /* used to iterate over active constraints */
      for (i = 0; i < nlpcons; i++)
      {
         if ( lpconsvals[i] < lplhs[ind] - feastol || lpconsvals[i] > lprhs[ind] + feastol)
         {
            SCIPdebugMessage("solution found infeasible (feastol=%g) for dual lp constraint: LP-%d = %g <|= [%g,%g]\n",
               feastol, i, lpconsvals[i], lplhs[ind], lprhs[ind]);
            BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
            *infeasible = TRUE;
            return SCIP_OKAY;
         }

         ind++;
      }
      BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
   }

   /* check sdp constraints */
   if ( nsdpblocks > 0 )
   {
      SCIP_Real* fullsdpmatrix;
      SCIP_Real eigenvalue;
      int maxblocksize = 0;

      /* allocate memory */
      if ( nsdpblocks == 1 )
         maxblocksize = sdpblocksizes[0] - nremovedinds[0];
      else
      {
         /* calculate maximum size of any SDP block to not have to reallocate memory in between */
         for (b = 0; b < nsdpblocks; b++)
         {
            if ( (sdpblocksizes[b] - nremovedinds[b]) > maxblocksize )
               maxblocksize = sdpblocksizes[b] - nremovedinds[b];
         }
      }

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &fullsdpmatrix, maxblocksize * maxblocksize) );/*lint !e647*/

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            /* initialize lower triangular part of fullsdpmatrix with zero */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j <= i; j++)
               {
                  fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j] = 0.0; /*lint !e679*/
               }
            }

            /* iterate over all non-fixed variables and add the corresponding nonzeros */
            for (v = 0; v < sdpnblockvars[b]; v++)
            {
               if ( lb[sdpvar[b][v]] < ub[sdpvar[b][v]] - epsilon )
               {
                  for (i = 0; i < sdpnblockvarnonz[b][v]; i++)
                  {
                     fullsdpmatrix[((sdprow[b][v][i] - indchanges[b][sdprow[b][v][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                                   sdpcol[b][v][i] - indchanges[b][sdpcol[b][v][i]]] += solvector[sdpvar[b][v]] * sdpval[b][v][i];/*lint !e679*/
                  }
               }
            }

            /* add constant matrix */
            if ( sdpconstnblocknonz != NULL )
            {
               for (i = 0; i < sdpconstnblocknonz[b]; i++)
               {
                  fullsdpmatrix[((sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                     sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]]] -= sdpconstval[b][i];/*lint !e679*/
               }
            }

            /* extend to full symmetric matrix for LAPACK */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j < i; j++)
               {
                  fullsdpmatrix[j * (sdpblocksizes[b] - nremovedinds[b]) + i] = fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j];/*lint !e679*/
               }
            }

            /* compute smallest eigenvalue using LAPACK */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, FALSE, sdpblocksizes[b] - nremovedinds[b], fullsdpmatrix, 1, &eigenvalue, NULL) );

            if ( eigenvalue < - feastol )
            {
               SCIPdebugMessage("solution found infeasible (feastol=%g) for dual sdp constraint %d, smallest eigenvector %.10g\n",
                  feastol, b, eigenvalue);
               BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
               *infeasible = TRUE;
               return SCIP_OKAY;
            }
         }
      }

      BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
   }

   *infeasible = FALSE;
   return SCIP_OKAY;
}


/** Given a solution, an SDP instance and a feasibility tolerance, checks whether the smallest eigenvalue is >= -feastol
 *  for a given feasibility tolerance and returns maximal absolute violation and sum of absolute violations for bounds,
 *  linear constraints and SDP constraints for the corresponding dual problem.
 *
 * Note: Should not be called if solution is a certificate of primal infeasiblity!
 */
SCIP_RETCODE SCIPsdpSolcheckerCheckAndGetViolDual(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            solvector,          /**< values of all variables (including fixed ones) in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Real*            maxabsviolbnds,     /**< pointer to store maximal absolute violation of variable bounds  */
   SCIP_Real*            sumabsviolbnds,     /**< pointer to store sum of absolute violations of variable bounds */
   SCIP_Real*            maxabsviolcons,     /**< pointer to store maximal absolute violation of linear constraints */
   SCIP_Real*            sumabsviolcons,     /**< pointer to store sum of absolute violations of linear constraints */
   SCIP_Real*            maxabsviolsdp,      /**< pointer to store maximal absolute violation of SDP constraints */
   SCIP_Real*            sumabsviolsdp,      /**< pointer to store sum of absolute violations of SDP constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
   )
{/*lint --e{818}*/
   int i;
   int j;
   int b;
   int v;
   int ind;
   SCIP_Real viol;

   assert( bufmem != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nsdpblocks >= 0 );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( nsdpblocks == 0 || sdpnblockvars != NULL );
   assert( sdpconstnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpnblockvarnonz != NULL );
   assert( nsdpblocks == 0 || sdpvar != NULL );
   assert( nsdpblocks == 0 || sdprow != NULL );
   assert( nsdpblocks == 0 || sdpcol != NULL );
   assert( nsdpblocks == 0 || sdpval != NULL );
   assert( nsdpblocks == 0 || indchanges != NULL );
   assert( nsdpblocks == 0 || nremovedinds != NULL );
   assert( nsdpblocks == 0 || blockindchanges != NULL );
   assert( nlpcons >= 0 );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lprow != NULL );
   assert( nlpcons == 0 || lpcol != NULL );
   assert( nlpcons == 0 || lpval != NULL );
   assert( solvector != NULL );
   assert( feastol >= 0 );
   assert( maxabsviolbnds != NULL );
   assert( sumabsviolbnds != NULL );
   assert( maxabsviolcons != NULL );
   assert( sumabsviolcons != NULL );
   assert( maxabsviolsdp != NULL );
   assert( sumabsviolsdp != NULL );
   assert( infeasible != NULL );

   /* initialize violations with 0.0 */
   *maxabsviolbnds = 0.0;
   *sumabsviolbnds = 0.0;
   *maxabsviolcons = 0.0;
   *sumabsviolcons = 0.0;
   *maxabsviolsdp  = 0.0;
   *sumabsviolsdp  = 0.0;

   *infeasible = FALSE;

   /* check variable bounds */
   for (i = 0; i < nvars; i++)
   {
      viol = MAX3(lb[i] - solvector[i], solvector[i] - ub[i], 0.0);
      *sumabsviolbnds += viol;

#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("Dual variable %d: lb = %g, ub = %g, val = %g, viol = %g\n", i, lb[i], ub[i], solvector[i], viol);
#endif

      if ( viol > *maxabsviolbnds )
         *maxabsviolbnds = viol;

      if ( solvector[i] < lb[i] - feastol || solvector[i] > ub[i] + feastol )
      {
         SCIPdebugMessage("solution found infeasible (feastol=%f) for dual variable bounds: x[%d] = %g <|= [%g, %g]\n",
            feastol, i, solvector[i], lb[i], ub[i]);
         *infeasible = TRUE;
      }
   }

   /* check linear constraints (since DSDP sorts the lp-arrays by cols, we cannot expect them to be sorted) */
   if ( nlpcons > 0 )
   {
      SCIP_Real* lpconsvals;

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &lpconsvals, nlpcons) );

      /* initialize all rows with zero */
      for (i = 0; i < nlpcons; i++)
         lpconsvals[i] = 0;

      /* compute the values of all rows */
      for (i = 0; i < lpnnonz; i++)
      {
         if ( lb[lpcol[i]] < ub[lpcol[i]] - epsilon ) /* fixed variables are already included in lhs/rhs */
            lpconsvals[lprow[i]] += solvector[lpcol[i]] * lpval[i];
      }

      /* check all active constraints for feasibility */
      ind = 0; /* used to iterate over active constraints */
      for (i = 0; i < nlpcons; i++)
      {
         viol = MAX3(lplhs[ind] - lpconsvals[i], lpconsvals[i] - lprhs[ind], 0.0);
         *sumabsviolcons += viol;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Dual linear constraint %d: lhs = %g, rhs = %g, val = %.15g, viol = %.15g\n", i, lplhs[i], lprhs[i], lpconsvals[i], viol);
#endif

         if ( viol > *maxabsviolcons )
            *maxabsviolcons = viol;

         if ( lpconsvals[i] < lplhs[ind] - feastol || lpconsvals[i] > lprhs[ind] + feastol)
         {
            SCIPdebugMessage("solution found infeasible (feastol=%g) for dual lp constraint: LP-%d = %g <|= [%g,%g].\n",
               feastol, i, lpconsvals[i], lplhs[ind], lprhs[ind]);
            BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
            *infeasible = TRUE;
         }

         ind++;
      }
      BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
   }

   /* check sdp constraints */
   if ( nsdpblocks > 0 )
   {
      SCIP_Real* fullsdpmatrix;
      SCIP_Real eigenvalue;
      int maxblocksize = 0;

      /* allocate memory */
      if ( nsdpblocks == 1 )
         maxblocksize = sdpblocksizes[0] - nremovedinds[0];
      else
      {
         /* calculate maximum size of any SDP block to not have to reallocate memory in between */
         for (b = 0; b < nsdpblocks; b++)
         {
            if ( (sdpblocksizes[b] - nremovedinds[b]) > maxblocksize )
               maxblocksize = sdpblocksizes[b] - nremovedinds[b];
         }
      }

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &fullsdpmatrix, maxblocksize * maxblocksize) );/*lint !e647*/

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            /* initialize lower triangular part of fullsdpmatrix with zero */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j <= i; j++)
               {
                  fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j] = 0.0; /*lint !e679*/
               }
            }

            /* iterate over all non-fixed variables and add the corresponding nonzeros */
            for (v = 0; v < sdpnblockvars[b]; v++)
            {
               if ( lb[sdpvar[b][v]] < ub[sdpvar[b][v]] - epsilon )
               {
                  for (i = 0; i < sdpnblockvarnonz[b][v]; i++)
                  {
                     fullsdpmatrix[((sdprow[b][v][i] - indchanges[b][sdprow[b][v][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                                   sdpcol[b][v][i] - indchanges[b][sdpcol[b][v][i]]] += solvector[sdpvar[b][v]] * sdpval[b][v][i];/*lint !e679*/
                  }
               }
            }

            /* add constant matrix */
            if ( sdpconstnblocknonz != NULL )
            {
               for (i = 0; i < sdpconstnblocknonz[b]; i++)
               {
                  fullsdpmatrix[((sdpconstrow[b][i] - indchanges[b][sdpconstrow[b][i]]) * (sdpblocksizes[b] - nremovedinds[b])) +
                     sdpconstcol[b][i] - indchanges[b][sdpconstcol[b][i]]] -= sdpconstval[b][i];/*lint !e679*/
               }
            }

            /* extend to full symmetric matrix for LAPACK */
            for (i = 0; i < sdpblocksizes[b] - nremovedinds[b]; i++)
            {
               for (j = 0; j < i; j++)
               {
                  fullsdpmatrix[j * (sdpblocksizes[b] - nremovedinds[b]) + i] = fullsdpmatrix[i * (sdpblocksizes[b] - nremovedinds[b]) + j];/*lint !e679*/
               }
            }

            /* compute smallest eigenvalue using LAPACK */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, FALSE, sdpblocksizes[b] - nremovedinds[b], fullsdpmatrix, 1, &eigenvalue, NULL) );

            viol = MAX(-eigenvalue, 0.0);
            *sumabsviolsdp += viol;

#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("Dual SDP constraint %d: lambda_min = %.15g, viol = %.15g\n", b, eigenvalue, viol);
#endif

            if ( viol > *maxabsviolsdp )
               *maxabsviolsdp = viol;

            if ( eigenvalue < - feastol )
            {
               SCIPdebugMessage("solution found infeasible (feastol=%.10g) for dual SDP constraint %d, smallest eigenvector %.10g\n",
                     feastol, b, eigenvalue);
               BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
               *infeasible = TRUE;
            }
         }
      }

      BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
   }

   return SCIP_OKAY;
}


/** Given a solution, an SDP instance and a feasibility tolerance, returns maximal absolute violation and sum of
 *  absolute violations for bounds, linear constraints and SDP constraints for the corresponding primal problem. The
 *  solution should be given for the actual problem given to the SDP solver, i.e., after all fixed variables have been
 *  eliminated.
 *
 * Note: Should not be called if solution is a certificate of dual infeasiblity!
 */
SCIP_RETCODE SCIPsdpSolcheckerCheckAndGetViolPrimal(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables (in dual problem) */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int*                  inputtomosekmapper, /**< entry i gives the index of input variable i in MOSEK (starting from 0) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            solvector,          /**< values of all scalar variables in the solution that should be checked */
   SCIP_Real**           solmatrices,        /**< values of all matrix variables in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Real*            maxabsviolbnds,     /**< pointer to store maximal absolute violation of variable bounds  */
   SCIP_Real*            sumabsviolbnds,     /**< pointer to store sum of absolute violations of variable bounds */
   SCIP_Real*            maxabsviolcons,     /**< pointer to store maximal absolute violation of linear constraints */
   SCIP_Real*            sumabsviolcons,     /**< pointer to store sum of absolute violations of linear constraints */
   SCIP_Real*            maxabsviolsdp,      /**< pointer to store maximal absolute violation of SDP constraints */
   SCIP_Real*            sumabsviolsdp,      /**< pointer to store sum of absolute violations of SDP constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
   )
{/*lint --e{818}*/
   int i;
   int j;
   int b;
   int v;
   SCIP_Real viol;
   int nactivevars;
   int nvarbounds;
   int nfixedvars;
   SCIP_Real* objcoefs;
   int* varboundpos;
   int nlpvars;
   int nprimalvars;
   int nprimalmatrixvars;
   int maxblocksize;
   int* mosekblocksizes;
   int nprimallpcons;
   int blocksize;
   int mosekind;
   int mosekrow;
   int mosekcol;
   int blockvar;
   int k;
   int nnonz;
   int c;

   assert( bufmem != NULL );
   assert( lb != NULL );
   assert( ub != NULL );
   assert( nsdpblocks >= 0 );
   assert( nsdpblocks == 0 || sdpblocksizes != NULL );
   assert( nsdpblocks == 0 || sdpnblockvars != NULL );
   assert( sdpconstnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstnblocknonz != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( nsdpblocks == 0 || sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz >= 0 );
   assert( nsdpblocks == 0 || sdpnblockvarnonz != NULL );
   assert( nsdpblocks == 0 || sdpvar != NULL );
   assert( nsdpblocks == 0 || sdprow != NULL );
   assert( nsdpblocks == 0 || sdpcol != NULL );
   assert( nsdpblocks == 0 || sdpval != NULL );
   assert( nsdpblocks == 0 || indchanges != NULL );
   assert( nsdpblocks == 0 || nremovedinds != NULL );
   assert( nsdpblocks == 0 || blockindchanges != NULL );
   assert( 0 <= nremovedblocks && nremovedblocks <= nsdpblocks );
   assert( nlpcons >= 0 );
   assert( nlpcons == 0 || lplhs != NULL );
   assert( nlpcons == 0 || lprhs != NULL );
   assert( lpnnonz >= 0 );
   assert( nlpcons == 0 || lprow != NULL );
   assert( nlpcons == 0 || lpcol != NULL );
   assert( nlpcons == 0 || lpval != NULL );
   assert( solvector != NULL );
   assert( feastol >= 0 );
   assert( maxabsviolbnds != NULL );
   assert( sumabsviolbnds != NULL );
   assert( maxabsviolcons != NULL );
   assert( sumabsviolcons != NULL );
   assert( maxabsviolsdp != NULL );
   assert( sumabsviolsdp != NULL );
   assert( infeasible != NULL );

   /* initialize violations with 0.0 */
   *maxabsviolbnds = 0.0;
   *sumabsviolbnds = 0.0;
   *maxabsviolcons = 0.0;
   *sumabsviolcons = 0.0;
   *maxabsviolsdp  = 0.0;
   *sumabsviolsdp  = 0.0;

   *infeasible = FALSE;

   nactivevars = 0;
   nvarbounds = 0;
   nfixedvars = 0;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &objcoefs, nvars) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &varboundpos, 2 * nvars) );

   /* find fixed variables in dual problem */
   for (i = 0; i < nvars; i++)
   {
      if ( ub[i] - lb[i] <= epsilon )
         nfixedvars++;
      else
      {
         objcoefs[nactivevars] = obj[i];

         if ( lb[i] > -INF )
            varboundpos[nvarbounds++] = -(nactivevars + 1); /* negative sign means lower bound */

         if ( ub[i] < INF )
            varboundpos[nvarbounds++] = +(nactivevars + 1); /* positive sign means upper bound */

         nactivevars++;
      }
   }
   assert( nactivevars + nfixedvars == nvars );

   /* determine total number of sides in LP-constraints */
   nlpvars = 0;
   if ( nlpcons > 0 )
   {
      for (i = 0; i < nlpcons; i++)
      {
         if ( lplhs[i] > - INF )
            ++nlpvars;

         if ( lprhs[i] < INF )
            ++nlpvars;
      }
      assert( nlpvars <= 2 * nlpcons ); /* factor 2 comes from left- and right-hand-sides */
   }

   /* compute number of primal scalar variables */
   nprimalvars = nlpvars + nvarbounds;

   /* check variable bounds (they should all be nonnegative) */
   for (i = 0; i < nprimalvars; i++)
   {
      viol = MAX(-solvector[i], 0.0);
      *sumabsviolbnds += viol;

      if ( viol > *maxabsviolbnds )
         *maxabsviolbnds = viol;

#ifdef SCIP_MORE_DEBUG
      SCIPdebugMessage("Primal variable %d: val = %g, viol = %g\n", i, solvector[i], viol);
#endif

      if ( solvector[i] < - feastol )
      {
         SCIPdebugMessage("primal solution found infeasible (feastol=%f) for primal variable bounds: x[%d] = %g >= 0\n",
            feastol, i, solvector[i]);
         *infeasible = TRUE;
      }
   }

   /* compute sizes of primal matrix variables and save maximal size */
   nprimalmatrixvars = nsdpblocks - nremovedblocks;
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &mosekblocksizes, nprimalmatrixvars) );
   maxblocksize = 0;

   for (b = 0; b < nsdpblocks; b++)
   {
      if ( blockindchanges[b] > -1 )
      {
         assert( 0 <= blockindchanges[b] && blockindchanges[b] <= b && (b - blockindchanges[b]) <= (nsdpblocks - nremovedblocks) );
         mosekblocksizes[b - blockindchanges[b]] = sdpblocksizes[b] - nremovedinds[b];/*lint !e679*/

         if ( sdpblocksizes[b] - nremovedinds[b] > maxblocksize )
            maxblocksize = sdpblocksizes[b] - nremovedinds[b];
      }
   }


   /* compute number of primal linear constraints */
   nprimallpcons = nactivevars;

   /* check linear constraints */
   if ( nprimallpcons > 0 )
   {
      SCIP_Real* lpconsvals;

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &lpconsvals, nprimallpcons) );

      /* initialize all rows with zero */
      for (i = 0; i < nprimallpcons; i++)
         lpconsvals[i] = 0.0;

      /* first add the matrices A_i */
      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            for (blockvar = 0; blockvar < sdpnblockvars[b]; blockvar++)
            {
               v = inputtomosekmapper[sdpvar[b][blockvar]];

               /* check if the variable is active */
               if ( v > -1 )
               {
                  assert( v < nactivevars );

                  /* if there are removed indices, we have to adjust the column and row indices accordingly */
                  if ( nremovedinds[b] > 0 )
                  {
                     assert( sdpnblockvarnonz[b][blockvar] <= sdpnnonz );
                     for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
                     {
                        /* rows and cols with active nonzeros should not be removed */
                        assert( 0 <= indchanges[b][sdprow[b][blockvar][k]] && indchanges[b][sdprow[b][blockvar][k]] <= sdprow[b][blockvar][k] );
                        assert( 0 <= indchanges[b][sdpcol[b][blockvar][k]] && indchanges[b][sdpcol[b][blockvar][k]] <= sdpcol[b][blockvar][k] );

                        assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                        assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                        mosekrow = sdprow[b][blockvar][k] - indchanges[b][sdprow[b][blockvar][k]];
                        mosekcol = sdpcol[b][blockvar][k] - indchanges[b][sdpcol[b][blockvar][k]];

                        /* MOSEK stores matrices in column-first format?!? */
                        blocksize = mosekblocksizes[b - blockindchanges[b]];
                        mosekind = mosekcol * blocksize + mosekrow - mosekcol * (mosekcol + 1) / 2;
                        /* mosekind = mosekrow * (mosekrow + 1) / 2 + mosekcol; */

                        if ( mosekrow == mosekcol )
                           lpconsvals[v] += sdpval[b][blockvar][k] * solmatrices[b - blockindchanges[b]][mosekind];
                        else
                           lpconsvals[v] += 2 * sdpval[b][blockvar][k] * solmatrices[b - blockindchanges[b]][mosekind];
                     }
                     assert( k == sdpnblockvarnonz[b][blockvar] );
                  }
                  else
                  {
                     assert( sdpnblockvarnonz[b][blockvar] <= sdpnnonz );
                     for (k = 0; k < sdpnblockvarnonz[b][blockvar]; k++)
                     {
                        assert( 0 <= sdprow[b][blockvar][k] && sdprow[b][blockvar][k] < sdpblocksizes[b] );
                        assert( 0 <= sdpcol[b][blockvar][k] && sdpcol[b][blockvar][k] < sdpblocksizes[b] );

                        mosekrow = sdprow[b][blockvar][k];
                        mosekcol = sdpcol[b][blockvar][k];

                        /* MOSEK stores matrices in column-first format?!? */
                        blocksize = mosekblocksizes[b - blockindchanges[b]];
                        mosekind = mosekcol * blocksize + mosekrow - mosekcol * (mosekcol + 1) / 2;
                        /* mosekind = mosekrow * (mosekrow + 1) / 2 + mosekcol; */

                        if ( mosekrow == mosekcol )
                           lpconsvals[v] += sdpval[b][blockvar][k] * solmatrices[b - blockindchanges[b]][mosekind];
                        else
                           lpconsvals[v] += 2 * sdpval[b][blockvar][k] * solmatrices[b - blockindchanges[b]][mosekind];
                     }
                     assert( k == sdpnblockvarnonz[b][blockvar] );
                  }
               }
            }
         }
      }

      /* add the entries corresponding to the lp-constraints in the dual problem */
      if ( lpnnonz > 0 )
      {
         int currentrow;
         int varcnt = 0;

         currentrow = lprow[0];
         for (nnonz = 0; nnonz < lpnnonz; ++nnonz)
         {
            assert( nnonz == 0 || lprow[nnonz-1] <= lprow[nnonz] );  /* rows should be sorted */
            assert( lprow[nnonz] == currentrow );

            v = inputtomosekmapper[lpcol[nnonz]];
            if ( v >= 0 )
            {
               assert( v < nactivevars );

               /* treat left hand side */
               if ( lplhs[lprow[nnonz]] > - INF )
                  lpconsvals[v] += lpval[nnonz] * solvector[varcnt];

               /* treat right hand side */
               if ( lprhs[lprow[nnonz]] < INF )
                  lpconsvals[v] -= lpval[nnonz] * solvector[varcnt + 1];
            }

            /* we finished a row */
            if ( nnonz == lpnnonz - 1 || lprow[nnonz + 1] > currentrow )
            {
               /* treat left hand side */
               if ( lplhs[currentrow] > - INF )
                  varcnt++;

               /* treat right hand side */
               if ( lprhs[currentrow] < INF )
                  varcnt++;

               /* reset counters */
               if ( nnonz < lpnnonz )
                  currentrow = lprow[nnonz+1];
            }
         }
         assert( varcnt == nlpvars );
      }

      /* finally add the entries corresponding to varbounds in the dual problem */
      for (i = 0; i < nvarbounds; i++)
      {
         if ( varboundpos[i] < 0 )
         {
            /* lower bound */
            mosekrow = - varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
            assert( 0 <= mosekrow && mosekrow < nactivevars );
            lpconsvals[mosekrow] += solvector[nlpvars + i];
         }
         else
         {
            /* upper bound */
            assert( varboundpos[i] > 0 ); /* we should not have a zero as we wanted a clear differentiation between positive and negative */
            mosekrow = varboundpos[i] - 1; /* minus one because we added one to get strictly positive/negative values */
            assert( 0 <= mosekrow && mosekrow < nactivevars );
            lpconsvals[mosekrow] -= solvector[nlpvars + i];
         }
      }

      /* check all linear constraints for feasibility (they should equal the objective coefficient of the corresponding
       * dual variable) */
      for (c = 0; c < nprimallpcons; c++)
      {
         viol = REALABS(lpconsvals[c] - objcoefs[c]);
         *sumabsviolcons += viol;

         if ( viol > *maxabsviolcons )
            *maxabsviolcons = viol;

#ifdef SCIP_MORE_DEBUG
         SCIPdebugMessage("Primal linear constraint %d: lhs = %g, rhs = %g, viol = %g\n", c, lpconsvals[c], objcoefs[c], viol);
#endif

         if ( lpconsvals[c] < objcoefs[c] - feastol || lpconsvals[c] > objcoefs[c] + feastol )
         {
            SCIPdebugMessage("primal solution found infeasible (feastol=%g) for primal lp constraint: LP-%d = %g == [%g], viol = %g\n",
               feastol, c, lpconsvals[c], objcoefs[c], viol);
            *infeasible = TRUE;
         }

      }
      BMSfreeBufferMemoryArray(bufmem, &lpconsvals);
   }

   /* check sdp constraints */
   if ( nprimalmatrixvars > 0 )
   {
      SCIP_Real* fullsdpmatrix;
      SCIP_Real eigenvalue;

      assert( maxblocksize > 0 );

      BMS_CALL( BMSallocBufferMemoryArray(bufmem, &fullsdpmatrix, maxblocksize * maxblocksize) );/*lint !e647*/

      for (b = 0; b < nsdpblocks; b++)
      {
         if ( blockindchanges[b] > -1 )
         {
            blocksize = mosekblocksizes[b - blockindchanges[b]];

            /* iterate over all entries in the solution and create corresponding full matrix for LAPACK (in column-first format) */
            /* NOTE: MOSEK stores matrices in column-first format?!? */
            for (i = 0; i < blocksize; i++)
            {
               for (j = 0; j <= i; j++)
               {
                  mosekind = j * blocksize + i - j * (j + 1) / 2;
                  fullsdpmatrix[i * blocksize + j] = solmatrices[b - blockindchanges[b]][mosekind];
                  fullsdpmatrix[j * blocksize + i] = solmatrices[b - blockindchanges[b]][mosekind];
               }
            }

            /* compute smallest eigenvalue using LAPACK */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, FALSE, blocksize, fullsdpmatrix, 1, &eigenvalue, NULL) );

            viol = MAX(-eigenvalue, 0.0);
            *sumabsviolsdp += viol;

            if ( viol > *maxabsviolsdp )
               *maxabsviolsdp = viol;

#ifdef SCIP_MORE_DEBUG
            SCIPdebugMessage("Primal SDP variable %d: lambda_min = %g, viol = %g\n", b, eigenvalue, viol);
#endif

            if ( eigenvalue < - feastol )
            {
               SCIPdebugMessage("primal solution found infeasible (feastol=%.10g) for primal SDP variable %d, smallest eigenvector %.10g\n",
                     feastol, b - blockindchanges[b], eigenvalue);
               *infeasible = TRUE;
            }
         }
      }

      BMSfreeBufferMemoryArray(bufmem, &fullsdpmatrix);
   }

   BMSfreeBufferMemoryArray(bufmem, &mosekblocksizes);
   BMSfreeBufferMemoryArray(bufmem, &varboundpos);
   BMSfreeBufferMemoryArray(bufmem, &objcoefs);

   return SCIP_OKAY;
}
