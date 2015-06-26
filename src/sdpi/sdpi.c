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

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG*/

/**@file   sdpi.c
 * @brief  General interface methods for SDP-preprocessing (mainly fixing variables and removing empty rows/cols)
 * @author Tristan Gally
 */

#include <assert.h>

#include "sdpi/sdpisolver.h"
#include "sdpi/sdpi.h"
#include "scipsdp/SdpVarfixer.h"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */


/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMRY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** this will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem */
#define CHECK_IF_SOLVED(sdpi)  do                                                                             \
                      {                                                                                       \
                         if ( ! (sdpi->solved) )                                                              \
                         {                                                                                    \
                            SCIPerrorMessage("Tried to access solution information ahead of solving! \n");    \
                            return SCIP_LPERROR;                                                              \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** same as CHECK_IF_SOLVED, but this will be used in functions returning a boolean value */
#define CHECK_IF_SOLVED_BOOL(sdpi)  do                                                                        \
                      {                                                                                       \
                         if ( ! (sdpi->solved) )                                                              \
                         {                                                                                    \
                            SCIPerrorMessage("Tried to access solution information ahead of solving! \n");    \
                            return FALSE;                                                                     \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** duplicate an array that might be null (in that case null is returned, otherwise BMSduplicateMemory is called) */
#define DUPLICATE_ARRAY_NULL(blkmem, target, source, size) do                                                 \
                      {                                                                                       \
                         if (size > 0)                                                                        \
                            BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, target, source, size) );           \
                         else                                                                                 \
                            *target = NULL;                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** same as SCIP_CALL, but gives a SCIP_PARAMETERUNKNOWN error if it fails */
#define SCIP_CALL_PARAM(x)   do                                                                               \
                      {                                                                                       \
                         SCIP_RETCODE _restat_;                                                               \
                         if ( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                         {                                                                                    \
                            if ( _restat_ != SCIP_PARAMETERUNKNOWN )                                          \
                            {                                                                                 \
                               SCIPerrorMessage("Error <%d> in function call\n", _restat_);                   \
                               SCIPABORT();                                                                   \
                            }                                                                                 \
                            return _restat_;                                                                  \
                         }                                                                                    \
                      }                                                                                       \
                      while( FALSE )

/** data for SDPI */
struct SCIP_SDPi
{
   SCIP_SDPISOLVER*      sdpisolver;         /**< pointer to the interface for the SDP Solver */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            obj;                /**< objective function values of variables */
   SCIP_Real*            lb;                 /**< lower bounds of variables */
   SCIP_Real*            ub;                 /**< upper bounds of variables */
   int                   nsdpblocks;         /**< number of SDP-blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP-blocks */
   int*                  sdpnblockvars;      /**< number of variables in each SDP-block */

   /* constant SDP data: */
   int                   sdpconstnnonz;      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                               *  number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */

   /* non-constant SDP data: */
   int                   sdpnnonz;           /**< number of nonzero elements in the SDP-constraint matrices */
   int**                 sdpnblockvarnonz;   /**< sdpnblockvarnonz[i][j] gives the number of nonzeros for the j-th variable (not necessarly
                                               *  variable j) in the i-th block, this is also the length of row/col/val[i][j] */
   int**                 sdpvar;             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                               *  in the i-th block */
   int***                sdprow;             /**< pointer to the row-indices for each block and variable in this block, so row[i][j][k] gives
                                               *  the k-th nonzero of the j-th variable (not necessarly variable j) in the i-th block */
   int***                sdpcol;             /**< pointer to the column-indices for each block and variable in this block */
   SCIP_Real***          sdpval;             /**< pointer to the values of the nonzeros for each block and variable in this block */

   /* lp data: */
   int                   nlpcons;            /**< number of LP-constraints */
   SCIP_Real*            lplhs;              /**< left hand sides of LP rows */
   SCIP_Real*            lprhs;              /**< right hand sides of LP rows */
   int                   lpnnonz;            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow;              /**< row-index for each entry in lpval-array */
   int*                  lpcol;              /**< column-index for each entry in lpval-array */
   SCIP_Real*            lpval;              /**< values of LP-constraint matrix entries */

   /* other data */
   SCIP_Bool             slatercheck;        /**< should the Slater condition for the dual problem be checked ahead of each solving process */
   SCIP_Bool             solved;             /**< was the problem solved since the last change */
   SCIP_Bool             penalty;            /**< was the last solved problem a penalty formulation */
   SCIP_Bool             infeasible;         /**< was infeasibility detected in presolving */
   int                   sdpid;              /**< counter for the number of SDPs solved */
   SCIP_Real             epsilon;            /**< this is used for checking if primal and dual objective are equal */
   SCIP_Real             feastol;            /**< this is used to check if the SDP-Constraint is feasible */
};

/*
 * Local Functions
 */

/** For given row and column (i,j) checks if i >= j, so that i and j give a position in the lower
 *  triangular part, otherwise i and j will be switched. This function will be called whenever a position in a symmetric matrix
 *  is given, to prevent problems if position (i,j) is given but later (j,i) should be changed.
 */
static
void ensureLowerTriangular(
  int*                   i,                  /**< row index */
  int*                   j                   /**< column index */
  )
{
   if ( *i < *j )
   {
      int temp;
      temp = *i;
      *i = *j;
      *j = temp;
   }
}

#ifndef NDEBUG
/** tests if for a given variable the lower bound is in an epsilon neighborhood of the upper bound */
static
SCIP_Bool isFixed(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   int                   v                   /**< global index of the variable to check this for */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert ( sdpi != NULL );
   assert ( v < sdpi->nvars );
   assert ( sdpi->lb != NULL );
   assert ( sdpi->ub != NULL );

   lb = sdpi->lb[v];
   ub = sdpi->ub[v];

   assert( lb < ub + sdpi->feastol || sdpi->infeasible );

   return ( REALABS(ub-lb) <= sdpi->feastol );
}
#else
#define isFixed(sdpi, v) (REALABS(sdpi->ub[v] - sdpi->lb[v]) <= sdpi->epsilon)
#endif

/** Computes the constant Matrix after all variables with lb=ub have been fixed and their nonzeros were moved to the constant part. The five variables
 *  other than sdpi are used to return the matrix.
 *
 *  The size of sdpconstnblocknonz and the first pointers of sdpconst row/col/val should be equal to sdpi->nsdpblocks,
 *  the size of sdpconst row/col/val [i], which is given in sdpconstblocknnonz, needs to be sufficient, otherwise the
 *  needed length will be returned in sdpconstnblocknonz and a debug message will be thrown.
 */
static
SCIP_RETCODE compConstMatAfterFixings(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   int*                  sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval         /**< pointers to the values of the nonzeros for each block */
   )
{
   int i;
   int v;
   int block;
   int* nfixednonz;
   int** fixedrows;
   int** fixedcols;
   SCIP_Real** fixedvals;

   assert ( sdpi != NULL );
   assert ( sdpconstnnonz != NULL );
   assert ( sdpconstnblocknonz != NULL );
   assert ( sdpconstrow != NULL );
   assert ( sdpconstcol != NULL );
   assert ( sdpconstval != NULL );
#ifndef NDEBUG
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      assert ( sdpconstrow[block] != NULL );
      assert ( sdpconstcol[block] != NULL );
      assert ( sdpconstval[block] != NULL );
   }
#endif

   fixedrows = NULL;
   fixedcols = NULL;
   fixedvals = NULL;

   /* allocate memory for the nonzeros that need to be fixed, as this is only temporarly needed, we allocate as much as theoretically possible */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &nfixednonz, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedrows, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedcols, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &fixedvals, sdpi->nsdpblocks) );

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      /* compute the number of fixed nonzeros in this block */
      nfixednonz[block] = 0;
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if (isFixed(sdpi, sdpi->sdpvar[block][v]))
            nfixednonz[block] += sdpi->sdpnblockvarnonz[block][v];
      }

      fixedrows[block] = NULL;
      fixedcols[block] = NULL;
      fixedvals[block] = NULL;

      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedrows[block]), nfixednonz[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedcols[block]), nfixednonz[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(fixedvals[block]), nfixednonz[block]) );

      /* set nfixednonz to 0 to use it for indexing later (at the end of the next for-block it will again have the same value) */
      nfixednonz[block] = 0;
   }

   /* iterate over all variables, saving the nonzeros of the fixed ones */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if (isFixed(sdpi, sdpi->sdpvar[block][v]))
         {
            for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; i++)
            {
               fixedrows[block][nfixednonz[block]] = sdpi->sdprow[block][v][i];
               fixedcols[block][nfixednonz[block]] = sdpi->sdpcol[block][v][i];
               /* this is the final value to add, so we no longer have to remember, from which variable this nonzero comes,
                * the -1 comes from +y_iA_i but -A_0 */
               fixedvals[block][nfixednonz[block]] = - sdpi->sdpval[block][v][i] * sdpi->lb[sdpi->sdpvar[block][v]];
               nfixednonz[block]++;
            }
         }
      }
   }

   /* compute the constant matrix */
   *sdpconstnnonz = 0;
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      SCIP_CALL( SCIPsdpVarfixerMergeArraysIntoNew(sdpi->blkmem, sdpi->feastol, sdpi->sdpconstrow[block], sdpi->sdpconstcol[block], sdpi->sdpconstval[block],
                                               sdpi->sdpconstnblocknonz[block], fixedrows[block], fixedcols[block], fixedvals[block], nfixednonz[block],
                                               sdpconstrow[block], sdpconstcol[block], sdpconstval[block], &sdpconstnblocknonz[block]) );
      *sdpconstnnonz += sdpconstnblocknonz[block];
   }

   /* free all memory */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedvals[block]), nfixednonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedcols[block]), nfixednonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(fixedrows[block]), nfixednonz[block]);
   }
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedvals, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedcols, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &fixedrows, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &nfixednonz, sdpi->nsdpblocks);

   return SCIP_OKAY;
}

/** This takes the sdpi and the computed constant matrix after fixings as input and checks for empty rows and columns in each block, which should be
 *  removed to not harm the Slater condition. It also removes SDP blocks with no entries left, these are returned in blockindchanges and nremovedblocks.
 */
static
SCIP_RETCODE findEmptyRowColsSDP(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int**                 indchanges,         /**< pointer to store the changes needed to be done to the indices, if indchange[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this, i.e.
                                              *   the value to decrease this index by, this array should have memory allocated in the size
                                              *   sdpi->nsdpblocks times sdpi->sdpblocksizes[block] */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< pointer to store index change for each block, system is the same as for indchanges */
   int*                  nremovedblocks      /**< pointer to store the number of blocks to be removed from the SDP */
   )
{
   int block;
   int v;
   int i;
   int nfoundinds;

   assert( sdpi != NULL );
   assert( sdpconstnblocknonz != NULL );
   assert( sdpconstrow != NULL );
   assert( sdpconstcol != NULL );
   assert( sdpconstval != NULL );
   assert( indchanges != NULL );
   assert( nremovedinds != NULL );
   assert( blockindchanges != NULL );
   assert( nremovedblocks != NULL );

   /* initialize indchanges with -1 */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (i = 0; i < sdpi->sdpblocksizes[block]; i++)
         indchanges[block][i] = -1;
   }
   *nremovedblocks = 0;

   /* iterate over all active nonzeros, setting the values of indchange for their row and col to 1 (this is an intermediate value to save that the
    * index is still needed, it will later be set to the number of rows/cols deleted earlier) */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      /* the number of indices already found in this block, saved for prematurely stopping the loops */
      nfoundinds = 0;
      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         if ( ! (isFixed(sdpi, sdpi->sdpvar[block][v])) )
         {
            for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; i++)
            {
               assert ( REALABS(sdpi->sdpval[block][v][i]) > sdpi->epsilon); /* this should really be a nonzero */
               if ( indchanges[block][sdpi->sdprow[block][v][i]] == -1 )
               {
                  indchanges[block][sdpi->sdprow[block][v][i]] = 1;
                  nfoundinds++;
               }
               if ( indchanges[block][sdpi->sdpcol[block][v][i]] == -1 )
               {
                  indchanges[block][sdpi->sdpcol[block][v][i]] = 1;
                  nfoundinds++;
               }
               if ( nfoundinds == sdpi->sdpblocksizes[block] )
                  break;   /* we're done for this block */
            }
         }
         if (nfoundinds == sdpi->sdpblocksizes[block])
            break;   /* we're done for this block */
      }

      if ( nfoundinds < sdpi->sdpblocksizes[block] )
      {
         /* if some indices haven't been found yet, look in the constant part for them */
         for (i = 0; i < sdpconstnblocknonz[block]; i++)
         {
            assert ( REALABS(sdpconstval[block][i]) > sdpi->epsilon); /* this should really be a nonzero */
            if ( indchanges[block][sdpconstrow[block][i]] == -1 )
            {
               indchanges[block][sdpconstrow[block][i]] = 1;
               nfoundinds++;
            }
            if ( indchanges[block][sdpconstcol[block][i]] == -1 )
            {
               indchanges[block][sdpconstcol[block][i]] = 1;
               nfoundinds++;
            }
            if ( nfoundinds == sdpi->sdpblocksizes[block] )
               break;   /* we're done for this block */
         }
      }

      /* now iterate over all indices to compute the final values of indchanges, all 0 are set to -1, all 1 are changed to the number of -1 before it */
      nremovedinds[block] = 0;
      for (i = 0; i < sdpi->sdpblocksizes[block]; i++)
      {
         if ( indchanges[block][i] == -1 )
         {
            SCIPdebugMessage("empty row and col %d were removed from block %d of SDP %d\n", i, block, sdpi->sdpid);
            /* this index wasn't found (indchanges was initialized with 0), so it can be removed */
            nremovedinds[block]++;
         }
         else
         {
            /* this index has been found, so set the value to the number of removed inds before it */
            indchanges[block][i] = nremovedinds[block];
         }
      }

      /* check if the block became empty */
      if ( nremovedinds[block] == sdpi->sdpblocksizes[block] )
      {
         SCIPdebugMessage("empty block %d detected in SDP %d, this will be removed", block, sdpi->sdpid);
         blockindchanges[block] = -1;
         (*nremovedblocks)++;
      }
      else
         blockindchanges[block] = *nremovedblocks;
   }

   return SCIP_OKAY;
}

/** computes the number of active variables for each constraint, thereby detecting constraints that
 *  may be removed, and computes the LP-left- and right-hand-sides after including all locally fixed variables
 *  for all constraints with at least two remaining active variables
 */
static
SCIP_RETCODE computeLpLhsRhsAfterFixings(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   int*                  nactivelpcons,      /**< output: number of active LP-constraints */
   SCIP_Real*            lplhsafterfix,      /**< output: first nlpcons (output) entries give left-hand sides of
                                              *           active lp-constraints after fixing variables, these are
                                              *           in the same relative order as before (with non-active rows
                                              *           removed) */
   SCIP_Real*            lprhsafterfix,      /**< output: first nlpcons (output) entries give right-hand sides of
                                              *  		    active lp-constraints after fixing variables, these are
                                              *  		    in the same relative order as before (with non-active rows
                                              *  		    removed) */
   int*                  rownactivevars,     /**< output: number of active variable for every row */
   SCIP_Bool*            fixingsfound        /**< output: returns true if a variable was fixed during this function call */
   )
{
   int i;
   int c;
   int lastrow = -1;
   int nonzind = -1;
   int nonzcol = -1;
   SCIP_Real nonzval;

   assert( sdpi != NULL );
   assert( nactivelpcons != NULL );
   assert( sdpi->nlpcons == 0 || lplhsafterfix != NULL );
   assert( sdpi->nlpcons == 0 || lprhsafterfix != NULL );
   assert( sdpi->nlpcons == 0 || rownactivevars != NULL );
   assert( sdpi->nlpcons == 0 || fixingsfound != NULL );

   /* if there is no LP-part, there is nothing to do */
   if ( sdpi->nlpcons == 0 || sdpi->lpnnonz == 0 )
   {
      *nactivelpcons = 0;
      return SCIP_OKAY;
   }

   /* initialize rownactivevars */
   for (c = 0; c < sdpi->nlpcons; c++)
      rownactivevars[c] = 0;
   *nactivelpcons = 0;

   /* sort the lp arrays by row indices */
#if 0 /* TODO thoroughly test this */
   SCIPsortIntIntReal(sdpi->lprow, sdpi->lpcol, sdpi->lpval, sdpi->lpnnonz);
#endif

   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      assert( i == 0 || sdpi->lprow[i-1] <= sdpi->lprow[i] );

      /* we reached a new row */
      if ( sdpi->lprow[i] > lastrow )
      {
         /* if the last row had at least two active variables, we keep the lhs- and rhs-value */
         if ( lastrow >= 0 && rownactivevars[lastrow] > 1 )
            (*nactivelpcons)++;
         else if ( lastrow >= 0 && rownactivevars[lastrow] == 1 )
         {
            assert( 0 <= nonzind && nonzind < sdpi->lpnnonz );

            nonzcol = sdpi->lpcol[nonzind];
            assert( 0 <= nonzcol && nonzcol < sdpi->nvars );

            nonzval = sdpi->lpval[nonzind];
            assert( REALABS(nonzval) > sdpi->feastol );

            /* we have to check if this is an improvement of the current bound */
            if ( nonzval < 0.0 ) /* we have to compare with the upper bound for lhs and lower bound for rhs */
            {
               /* check for the left-hand-side */
               if ( (lplhsafterfix[*nactivelpcons] > - SCIPsdpiInfinity(sdpi)) &&
                    ((lplhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->feastol) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
                  sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }
                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
               /* check for the right-hand-side */
               if ( (lprhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
                    ((lprhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->feastol) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
                  sdpi->lb[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }
                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
            }
            else	/* we compare with the lower bound for lhs and upper bound for rhs */
            {
               /* check for the left-hand-side */
               if ( (lplhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
                    ((lplhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->feastol) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
                  sdpi->lb[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }
                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
               /* check for the right-hand-side */
               if ( (lprhsafterfix[*nactivelpcons] > - SCIPsdpiInfinity(sdpi)) &&
                    ((lprhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->feastol) )
               {
                  /* this bound is sharper than the original one */
                  SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
                     "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
                  sdpi->ub[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

                  /* check if this leads to a fixing of this variable */
                  if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
                  {
                     *fixingsfound = TRUE;
                     SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                        nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
                  }
                  /* check if this makes the problem infeasible */
                  if (sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol)
                  {
                     sdpi->infeasible = TRUE;
                     SCIPdebugMessage("We found an upper bound %f that is lower than the lower bound %f for variable %d, so the problem is infeasible !\n",
                           sdpi->ub[nonzcol], sdpi->lb[nonzcol], nonzcol);
                     return SCIP_OKAY;
                  }
               }
            }
         }
         else if ( lastrow >= 0 ) /* because of earlier ifs we have rownactivevars = 0 */
         {
            assert( lastrow == -1 || rownactivevars[lastrow] == 0 );
            /* we have a constraint lhs <= 0 <= rhs, so lhs should be non-positive and rhs non-negative, otherwise the problem is infeasible */
            if ( lplhsafterfix[*nactivelpcons] > sdpi->feastol || lprhsafterfix[*nactivelpcons] < -sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found a constraint which with given fixings reads %f <= 0 <= %f, so the current problem is infeasible !\n",
                     lplhsafterfix[*nactivelpcons], lprhsafterfix[*nactivelpcons] );
               return SCIP_OKAY;
            }
         }

         /* update lastrow for new row */
         lastrow = sdpi->lprow[i];

         /* start the next lhr & rhs with the original value */
         lplhsafterfix[*nactivelpcons] = sdpi->lplhs[lastrow];
         lprhsafterfix[*nactivelpcons] = sdpi->lprhs[lastrow];
      }

      /* if the variable is active, we increase rownactivevars */
      if ( ! isFixed(sdpi, sdpi->lpcol[i]) )
      {
         rownactivevars[lastrow]++;
         nonzind = i;
      }
      else
      {
         /* otherwise we add the value (coefficient * value of fixed variable) to the lhs and rhs, the minus comes from +A_i but -A_0 */
         lplhsafterfix[*nactivelpcons] -= sdpi->lpval[i] * sdpi->lb[sdpi->lpcol[i]];
         lprhsafterfix[*nactivelpcons] -= sdpi->lpval[i] * sdpi->lb[sdpi->lpcol[i]];
      }
   }
   /* for the last row of the lp we have to check if it is active, as in the above for-queue we only do so when the next row start */
   if ( rownactivevars[lastrow] > 1 )
      (*nactivelpcons)++;
   else if ( rownactivevars[lastrow] == 1 )
   {
      assert( 0 <= nonzind && nonzind < sdpi->lpnnonz );

      nonzcol = sdpi->lpcol[nonzind];
      assert( 0 <= nonzcol && nonzcol < sdpi->nvars );

      nonzval = sdpi->lpval[nonzind];
      assert( REALABS(nonzval) > sdpi->epsilon );

      /* we have to check if this is an improvement of the current bound */
      if ( nonzval < 0.0 ) /* we have to compare with the upper bound for lhs and lower bound for rhs */
      {
         /* check for the left-hand-side */
         if ( (lplhsafterfix[*nactivelpcons] > SCIPsdpiInfinity(sdpi)) &&
              ((lplhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->feastol) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
            sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }
            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
         /* check for the right-hand-side */
         if ( (lprhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
              ((lprhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] - sdpi->feastol) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
            sdpi->lb[nonzcol] = lprhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }
            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
      }
      else  /* we compare with the lower bound for lhs and upper bound for rhs */
      {
         /* check for the left-hand-side */
         if ( (lplhsafterfix[*nactivelpcons] < SCIPsdpiInfinity(sdpi)) &&
              ((lplhsafterfix[*nactivelpcons] / nonzval) > sdpi->lb[nonzcol] + sdpi->feastol) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, lower bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lplhsafterfix[*nactivelpcons] / nonzval, sdpi->lb[nonzcol]);
            sdpi->lb[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }
            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found a lower bound that is bigger than the upper bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
         /* check for the right-hand-side */
         if ( (lprhsafterfix[*nactivelpcons] > SCIPsdpiInfinity(sdpi)) &&
              ((lprhsafterfix[*nactivelpcons] / nonzval) < sdpi->ub[nonzcol] - sdpi->feastol) )
         {
            /* this bound is sharper than the original one */
            SCIPdebugMessage("empty LP-row %d has been removed from SDP %d, upper bound of variable %d has been sharpened to %f "
               "(originally %f)\n", lastrow, sdpi->sdpid, nonzcol, lprhsafterfix[*nactivelpcons] / nonzval, sdpi->ub[nonzcol]);
            sdpi->ub[nonzcol] = lplhsafterfix[*nactivelpcons] / nonzval;

            /* check if this leads to a fixing of this variable */
            if ( REALABS(sdpi->lb[nonzcol] - sdpi->ub[nonzcol]) < sdpi->feastol )
            {
               *fixingsfound = TRUE;
               SCIPdebugMessage("computeLpLhsRhsAfterFixings fixed variable %d to value %f in SDP %d\n",
                  nonzcol, sdpi->lb[nonzcol], sdpi->sdpid);
            }
            /* check if this makes the problem infeasible */
            if ( sdpi->ub[nonzcol] < sdpi->lb[nonzcol] - sdpi->feastol )
            {
               sdpi->infeasible = TRUE;
               SCIPdebugMessage("We found an upper bound that is lower than the lower bound, so the problem is infeasible !\n");
               return SCIP_OKAY;
            }
         }
      }
   }
   else
   {
      assert( lastrow == -1 || rownactivevars[lastrow] == 0 );
      /* we have a constraint lhs <= 0 <= rhs, so lhs should be non-positive and rhs non-negative, otherwise the problem is infeasible */
      if ( lplhsafterfix[*nactivelpcons] > sdpi->feastol || lprhsafterfix[*nactivelpcons] < -sdpi->feastol )
      {
         sdpi->infeasible = TRUE;
         SCIPdebugMessage("We found a constraint which with given fixings reads %f <= 0 <= %f, so the current problem is infeasible !\n",
               lplhsafterfix[*nactivelpcons], lprhsafterfix[*nactivelpcons] );
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver */
const char* SCIPsdpiGetSolverName(
   void
   )
{
   return SCIPsdpiSolverGetSolverName();
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   void
   )
{
   return SCIPsdpiSolverGetSolverDesc();
}

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return SCIPsdpiSolverGetSolverPointer(sdpi->sdpisolver);
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert ( sdpi != NULL );
   assert ( blkmem != NULL );

   SCIPdebugMessage("Calling SCIPsdpiCreate\n");

   BMS_CALL(BMSallocBlockMemory(blkmem, sdpi));

   SCIP_CALL( SCIPsdpiSolverCreate(&((*sdpi)->sdpisolver), messagehdlr, blkmem) );

   (*sdpi)->messagehdlr = messagehdlr;
   (*sdpi)->blkmem = blkmem;
   (*sdpi)->sdpid = 1;
   (*sdpi)->nvars = 0;
   (*sdpi)->nsdpblocks = 0;
   (*sdpi)->sdpconstnnonz = 0;
   (*sdpi)->sdpnnonz = 0;
   (*sdpi)->nlpcons = 0;
   (*sdpi)->lpnnonz = 0;
   (*sdpi)->slatercheck = FALSE;
   (*sdpi)->solved = FALSE;
   (*sdpi)->penalty = FALSE;
   (*sdpi)->infeasible = FALSE;

   (*sdpi)->obj = NULL;
   (*sdpi)->lb = NULL;
   (*sdpi)->ub = NULL;
   (*sdpi)->sdpblocksizes = NULL;
   (*sdpi)->sdpnblockvars = NULL;
   (*sdpi)->sdpconstnblocknonz = NULL;
   (*sdpi)->sdpconstrow = NULL;
   (*sdpi)->sdpconstcol = NULL;
   (*sdpi)->sdpconstval = NULL;
   (*sdpi)->sdpnblockvarnonz = NULL;
   (*sdpi)->sdpvar = NULL;
   (*sdpi)->sdprow = NULL;
   (*sdpi)->sdpcol = NULL;
   (*sdpi)->sdpval = NULL;
   (*sdpi)->lplhs = NULL;
   (*sdpi)->lprhs = NULL;
   (*sdpi)->lprow = NULL;
   (*sdpi)->lpcol = NULL;
   (*sdpi)->lpval = NULL;

   (*sdpi)->epsilon = 1e-5;
   (*sdpi)->feastol = 1e-4;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP interface structure */
   )
{
   int i;
   int j;

   SCIPdebugMessage("Calling SCIPsdpiFree \n");
   assert ( sdpi != NULL );
   assert ( *sdpi != NULL );

   /* free the LP part */
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lpval), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lpcol), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lprow), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lprhs), (*sdpi)->nlpcons);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->lplhs), (*sdpi)->nlpcons);

   /* free the individual nonzeros */
   for (i = 0; i < (*sdpi)->nsdpblocks; i++)
   {
      for (j = 0; j < (*sdpi)->sdpnblockvars[i]; j++)
      {
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
         BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol[i][j]), (*sdpi)->sdpnblockvarnonz[i][j]);
      }
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpvar[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonz[i]), (*sdpi)->sdpnblockvars[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstval[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstrow[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstcol[i]), (*sdpi)->sdpconstnblocknonz[i]);
   }

   /* free the rest */
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonz), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstnblocknonz), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpval), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpcol), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdprow), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpvar), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstval), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstcol), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpconstrow), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpnblockvars), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArrayNull((*sdpi)->blkmem, &((*sdpi)->sdpblocksizes), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->ub), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lb), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->obj), (*sdpi)->nvars);

   /* free the solver */
   SCIP_CALL( SCIPsdpiSolverFree(&((*sdpi)->sdpisolver)) );

   BMSfreeBlockMemory((*sdpi)->blkmem, sdpi);

   return SCIP_OKAY;
}

/** cloning method of the general SDP-Interface
 *  @note the solver specific interface is created anew and not copied */
SCIP_RETCODE SCIPsdpiClone(
   SCIP_SDPI*            oldsdpi,            /**< pointer to the SDP interface structure that should be cloned */
   SCIP_SDPI*            newsdpi             /**< pointer to an SDP interface structure to clone into */
   )
{
   BMS_BLKMEM* blkmem;
   int nvars;
   int nsdpblocks;
   int lpnnonz;
   int b;
   int v;

   assert( oldsdpi != NULL );

   SCIPdebugMessage("Cloning SDPI %d\n", oldsdpi->sdpid);

   /* general data */
   blkmem = oldsdpi->blkmem;
   nvars = oldsdpi->nvars;
   nsdpblocks = oldsdpi->nsdpblocks;
   lpnnonz = oldsdpi->lpnnonz;

   BMS_CALL( BMSallocBlockMemory(blkmem, &newsdpi) );

   SCIP_CALL( SCIPsdpiSolverCreate(&(newsdpi->sdpisolver), oldsdpi->messagehdlr, oldsdpi->blkmem) ); /* create new SDP-Solver Interface */

   newsdpi->messagehdlr = oldsdpi->messagehdlr;
   newsdpi->blkmem = blkmem;
   newsdpi->nvars = nvars;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->obj), oldsdpi->obj, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lb), oldsdpi->lb, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->ub), oldsdpi->ub, nvars) );

   newsdpi->nsdpblocks = nsdpblocks;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpblocksizes), oldsdpi->sdpblocksizes, nsdpblocks) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvars), oldsdpi->sdpnblockvars, nsdpblocks) );

   /* constant SDP data */
   newsdpi->sdpconstnnonz = oldsdpi->sdpconstnnonz;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstnblocknonz), oldsdpi->sdpconstnblocknonz, nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstrow), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstcol), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpconstval), nsdpblocks) );

   for (b = 0; b < nsdpblocks; b++)
   {
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstrow[b]), oldsdpi->sdpconstrow[b], oldsdpi->sdpconstnblocknonz[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstcol[b]), oldsdpi->sdpconstcol[b], oldsdpi->sdpconstnblocknonz[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpconstval[b]), oldsdpi->sdpconstval[b], oldsdpi->sdpconstnblocknonz[b]) );
   }

   /* SDP data */
   newsdpi->sdpnnonz = oldsdpi->sdpnnonz;

   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvarnonz), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpvar), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdprow), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpcol), nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpval), nsdpblocks) );

   for (b = 0; b < nsdpblocks; b++)
   {
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpnblockvarnonz[b]), oldsdpi->sdpnblockvarnonz[b], oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpvar[b]), oldsdpi->sdpvar[b], oldsdpi->sdpnblockvars[b]) );

      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdprow[b]), oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpcol[b]), oldsdpi->sdpnblockvars[b]) );
      BMS_CALL( BMSallocBlockMemoryArray(blkmem, &(newsdpi->sdpval[b]), oldsdpi->sdpnblockvars[b]) );

      for (v = 0; v < oldsdpi->sdpnblockvars[b]; v++)
      {
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdprow[b][v]), oldsdpi->sdprow[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpcol[b][v]), oldsdpi->sdpcol[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
         BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->sdpval[b][v]), oldsdpi->sdpval[b][v], oldsdpi->sdpnblockvarnonz[b][v]) );
      }
   }

   /* LP data */
   newsdpi->nlpcons = oldsdpi->nlpcons;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lplhs), oldsdpi->lplhs, oldsdpi->nlpcons) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lprhs), oldsdpi->lprhs, oldsdpi->nlpcons) );

   newsdpi->lpnnonz = lpnnonz;

   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lprow), oldsdpi->lprow, lpnnonz) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lpcol), oldsdpi->lpcol, lpnnonz) );
   BMS_CALL( BMSduplicateBlockMemoryArray(blkmem, &(newsdpi->lpval), oldsdpi->lpval, lpnnonz) );

   /* other data */
   newsdpi->solved = FALSE; /* as we don't copy the sdpisolver, this needs to be set to false */
   newsdpi->penalty = FALSE; /* all things about SDP-solutions are set to false as well, as we didn't solve the problem */
   newsdpi->infeasible = FALSE;
   newsdpi->sdpid = 1000000 + oldsdpi->sdpid; /* this is only used for debug output, setting it to this value should make it clear, that it is a new sdpi */
   newsdpi->epsilon = oldsdpi->epsilon;
   newsdpi->feastol = oldsdpi->feastol;

   return SCIP_OKAY;
}

/**@} */


/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data into SDP solver
 *
 *  @note as the SDP-constraint matrices are symmetric, only the upper triangular part of them must be specified
 *  @note there must be at least one variable, the SDP- and/or LP-part may be empty
 */
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables in each SDP block (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                               *  number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointer to row-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int**                 sdpconstcol,        /**< pointer to column-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   SCIP_Real**           sdpconstval,        /**< pointer to values of entries of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrices */
   int**                 sdpnblockvarnonz,   /**< sdpnblockvarnonz[i][j] gives the number of nonzeros for the j-th variable (not necessarly
                                               *  variable j) in the i-th block, this is also the length of row/col/val[i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                               *  in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable in this block, so row[i][j][k] gives
                                               *  the k-th nonzero of the j-th variable (not necessarly variable j) in the i-th block
                                               *  (may be NULL if sdptnnonz = 0)*/
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable in this block (may be NULL if sdptnnonz = 0)*/
   SCIP_Real***          sdpval,             /**< pointer to the values of the nonzeros for each block and variable in this
                                               *  block (may be NULL if sdptnnonz = 0)*/
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lplhs,              /**< left hand sides of LP rows (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   )
{
   int i;
   int v;
   int block;

   SCIPdebugMessage("Calling SCIPsdpiLoadSDP (%d)\n",sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( nvars > 0 );
   assert ( obj != NULL );
   assert ( lb != NULL );
   assert ( ub != NULL );

#ifdef SCIP_DEBUG
   if (sdpconstnnonz > 0 || sdpnnonz > 0 || nsdpblocks > 0)
   {
      assert ( sdpblocksizes != NULL );
      assert ( sdpnblockvars != NULL );
      assert ( nsdpblocks > 0 );
      assert ( sdpconstnblocknonz != NULL );
      assert ( sdpnblockvarnonz != NULL );

      if (sdpconstnnonz > 0)
      {
         assert ( sdpconstrow != NULL );
         assert ( sdpconstcol != NULL );
         assert ( sdpconstval != NULL );

         for (i = 0; i < nsdpblocks; i++)
         {
            if (sdpconstnblocknonz[i] > 0)
            {
               assert ( sdpconstrow[i] != NULL );
               assert ( sdpconstcol[i] != NULL );
               assert ( sdpconstval[i] != NULL );
            }
         }
      }

      if (sdpnnonz > 0)
      {
         assert ( sdprow != NULL );
         assert ( sdpcol != NULL );
         assert ( sdpval != NULL );

         for ( i = 0; i < nsdpblocks; i++ )
         {
            assert ( sdpcol[i] != NULL );
            assert ( sdprow[i] != NULL );
            assert ( sdpval[i] != NULL );

            for ( v = 0; v < sdpnblockvars[i]; v++)
            {
               if (sdpnblockvarnonz[i][v] > 0)
               {
                  assert ( sdpcol[i][v] != NULL );
                  assert ( sdprow[i][v] != NULL );
                  assert ( sdpval[i][v] != NULL );
               }
            }
         }
      }
   }
#endif

   assert ( nlpcons == 0 || lplhs != NULL );
   assert ( nlpcons == 0 || lprhs != NULL );
   assert ( lpnnonz == 0 || lprow != NULL );
   assert ( lpnnonz == 0 || lpcol != NULL );
   assert ( lpnnonz == 0 || lpval != NULL );

   /* memory allocation */

   /* first free the old arrays */
   for (block = sdpi->nsdpblocks - 1; block >= 0; block--)
   {
      for (v = sdpi->sdpnblockvars[block] - 1; v >= 0; v--)
      {
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpval[block][v]), sdpi->sdpnblockvarnonz[block][v]);
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdprow[block][v]), sdpi->sdpnblockvarnonz[block][v]);
         BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpcol[block][v]), sdpi->sdpnblockvarnonz[block][v]);
      }

      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpval[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdprow[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpcol[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpnblockvarnonz[block]), sdpi->sdpnblockvars[block]);
      BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpvar[block]), sdpi->sdpnblockvars[block]);
   }

   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->ub), sdpi->nvars);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lb), sdpi->nvars);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->obj), sdpi->nvars);

   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpnblockvars), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->sdpconstnblocknonz), sdpi->nsdpblocks);

   /* duplicate some arrays */
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), obj, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), lb, nvars) );
   BMS_CALL( BMSduplicateBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), ub, nvars) );
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpblocksizes, nsdpblocks);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpnblockvars), sdpnblockvars, nsdpblocks);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstnblocknonz), sdpconstnblocknonz, nsdpblocks);

   /* allocate memory for the sdp arrays & duplicate them */
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonz), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcol), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrow), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpvar), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow), sdpi->nsdpblocks, nsdpblocks) );
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->nsdpblocks, nsdpblocks) );

   for (block = 0; block < nsdpblocks; block++)
   {
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpnblockvarnonz[block]), sdpnblockvarnonz[block], sdpnblockvars[block]);

      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpconstcol[block], sdpconstnblocknonz[block]);
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpconstrow[block], sdpconstnblocknonz[block]);
      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpconstval[block], sdpconstnblocknonz[block]);

      /* make sure that we have a lower triangular matrix */
      for (i = 0; i < sdpi->sdpconstnblocknonz[block]; ++i)
         ensureLowerTriangular(&(sdpconstrow[block][i]), &(sdpconstcol[block][i]));

      DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpvar[block]), sdpvar[block], sdpnblockvars[block]);

      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol[block]), sdpnblockvars[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow[block]), sdpnblockvars[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval[block]), sdpnblockvars[block]) );

      for (v = 0; v < sdpi->sdpnblockvars[block]; v++)
      {
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpcol[block][v]), sdpcol[block][v], sdpnblockvarnonz[block][v]);
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdprow[block][v]), sdprow[block][v], sdpnblockvarnonz[block][v]);
         DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->sdpval[block][v]), sdpval[block][v], sdpnblockvarnonz[block][v]);

         /* make sure that we have a lower triangular matrix */
         for (i = 0; i < sdpi->sdpnblockvarnonz[block][v]; ++i)
            ensureLowerTriangular(&(sdprow[block][v][i]), &(sdpcol[block][v][i]));
      }
   }

   /* free old and duplicate new arrays for the LP part */
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);
   BMSfreeBlockMemoryArrayNull(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons);

   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lplhs), lplhs, nlpcons);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lprhs), lprhs, nlpcons);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lprow), lprow, lpnnonz);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lpcol), lpcol, lpnnonz);
   DUPLICATE_ARRAY_NULL(sdpi->blkmem, &(sdpi->lpval), lpval, lpnnonz);

   /* set the general information */
   sdpi->nvars = nvars;
   sdpi->nsdpblocks = nsdpblocks;

   sdpi->sdpconstnnonz = sdpconstnnonz;
   sdpi->sdpnnonz = sdpnnonz;

   /* LP part */
   sdpi->lpnnonz = lpnnonz;
   sdpi->nlpcons = nlpcons;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;

   return SCIP_OKAY;
}

/** adds rows to the LP-Block
 *
 *  @note arrays are not checked for duplicates, problems may appear if indices are added more than once
 */
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            row,                /**< row indices of constraint matrix entries, going from 0 to nrows - 1, these will be changed
                                               *  to nlpcons + i */
   const int*            col,                /**< column indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   int i;

   SCIPdebugMessage("Adding %d LP-Constraints to SDP %d.\n", nrows, sdpi->sdpid);

   assert ( sdpi != NULL );

   if ( nrows == 0 )
      return SCIP_OKAY; /* nothing to do in this case */

   assert ( lhs != NULL );
   assert ( rhs != NULL );
   assert ( nnonz >= 0 );
   assert ( row != NULL );
   assert ( col != NULL );
   assert ( val != NULL );

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons, sdpi->nlpcons + nrows) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons + nrows) ); /*lint !e776*/

   for (i = 0; i < nrows; i++)
   {
      sdpi->lplhs[sdpi->nlpcons + i] = lhs[i]; /*lint !e679*/
      sdpi->lprhs[sdpi->nlpcons + i] = rhs[i]; /*lint !e679*/
   }

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + nnonz) ); /*lint !e776*/

   for (i = 0; i < nnonz; i++)
   {
      assert ( 0 <= row[i] && row[i] < nrows );
      /* the new rows are added at the end, so the row indices are increased by the old number of LP-constraints */
      sdpi->lprow[sdpi->lpnnonz + i] = row[i] + sdpi->nlpcons; /*lint !e679*/

      assert ( 0 <= col[i] && col[i] < sdpi->nvars ); /* only existing vars should be added to the LP-constraints */
      sdpi->lpcol[sdpi->lpnnonz + i] = col[i]; /*lint !e679*/

      sdpi->lpval[sdpi->lpnnonz + i] = val[i]; /*lint !e679*/
   }

   sdpi->nlpcons = sdpi->nlpcons + nrows;
   sdpi->lpnnonz = sdpi->lpnnonz + nnonz;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;

   return SCIP_OKAY;
}

/** deletes all rows in the given range from the LP-Block */
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int i;
   int deletedrows;
   int firstrowind;
   int lastrowind;
   int deletednonz;

   SCIPdebugMessage("Deleting rows %d to %d from SDP %d.\n", firstrow, lastrow, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( firstrow >= 0 );
   assert ( firstrow <= lastrow );
   assert ( lastrow < sdpi->nlpcons );

   /* shorten the procedure if the whole LP-part is to be deleted */
   if (firstrow == 0 && lastrow == sdpi->nlpcons - 1)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons);

      sdpi->lplhs = NULL;
      sdpi->lprhs = NULL;
      sdpi->lpcol = NULL;
      sdpi->lprow = NULL;
      sdpi->lpval = NULL;

      sdpi->nlpcons = 0;
      sdpi->lpnnonz = 0;

      sdpi->solved = FALSE;
      sdpi->infeasible = FALSE;

      return SCIP_OKAY;
   }

   deletedrows = lastrow - firstrow + 1; /*lint !e834*/
   deletednonz = 0;

   /* first delete the left- and right-hand-sides */
   for (i = lastrow + 1; i < sdpi->nlpcons; i++) /* shift all rhs after the deleted rows */
   {
      sdpi->lplhs[i - deletedrows] = sdpi->lplhs[i]; /*lint !e679*/
      sdpi->lprhs[i - deletedrows] = sdpi->lprhs[i]; /*lint !e679*/
   }
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lplhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows) ); /*lint !e776*/

   /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
   SCIPsortIntIntReal(sdpi->lprow, sdpi->lpcol, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

   firstrowind = -1;
   /*iterate over the lprowind array to find the first index belonging to a row that should be deleted */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprow[i] >= firstrow && sdpi->lprow[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
      {
         firstrowind = i;
         lastrowind = i;
         i++;
         break;
      }
   }

   if (firstrowind > -1) /* if this is still 0 there are no nonzeroes for the given rows */
   {
      /* now find the last occurence of one of the rows (as these are sorted all in between also belong to deleted rows and will be removed) */
      while (i < sdpi->lpnnonz && sdpi->lprow[i] <= lastrow)
      {
         lastrowind++; /*lint !e644*/
         i++;
      }
      deletednonz = lastrowind - firstrowind + 1; /*lint !e834*/

      /* finally shift all LP-array-entries after the deleted rows */
      for (i = lastrowind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcol[i - deletednonz] = sdpi->lpcol[i]; /*lint !e679*/
         /* all rowindices after the deleted ones have to be lowered to still have ongoing indices from 0 to nlpcons-1 */
         sdpi->lprow[i - deletednonz] = sdpi->lprow[i] - deletedrows;  /*lint !e679*/
         sdpi->lpval[i - deletednonz] = sdpi->lpval[i]; /*lint !e679*/
      }
   }

   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcol), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprow), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz) ); /*lint !e776*/
   sdpi->nlpcons = sdpi->nlpcons - deletedrows;
   sdpi->lpnnonz = sdpi->lpnnonz - deletednonz;

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;

   return SCIP_OKAY;
}

/** deletes LP rows from SCIP_SDPI */
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of LP rows <br>
                                              *   input:  1 if row should be deleted, 0 otherwise <br>
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int i;
   int oldnlpcons;
   int deletedrows;

   SCIPdebugMessage("Calling SCIPsdpiDelLPRowset for SDP %d.\n", sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( dstat != NULL );

   oldnlpcons = sdpi->nlpcons;
   deletedrows = 0;

   for (i = 0; i < oldnlpcons; i++)
   {
      if (dstat[i] == 1)
      {
         /* delete this row, it is shifted by - deletedrows, because in this problem the earlier rows have already been deleted */
         SCIP_CALL( SCIPsdpiDelLPRows(sdpi, i - deletedrows, i - deletedrows) );
         dstat[i] = -1;
         deletedrows++;
      }
      else
         dstat[i] = i - deletedrows;
   }

   sdpi->solved = FALSE;
   sdpi->infeasible = FALSE;

   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );

   SCIPdebugMessage("Called SCIPsdpiClear in SDP %d.\n", sdpi->sdpid);

   /* we reset all counters */
   sdpi->sdpid = 1;
   SCIP_CALL( SCIPsdpiSolverResetCounter(sdpi->sdpisolver) );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of variables */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   int i;

   SCIPdebugMessage("Changing %d variable bounds in SDP %d\n", nvars, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( ind != NULL );
   assert ( lb != NULL );
   assert ( ub != NULL );

   for (i = 0; i < nvars; i++)
   {
      assert ( 0 <= ind[i] && ind[i] < sdpi->nvars );
      sdpi->lb[ind[i]] = lb[i];
      sdpi->ub[ind[i]] = ub[i];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes left and right hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPLhRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("Changing %d left and right hand sides of SDP %d\n", nrows, sdpi->sdpid);

   assert( sdpi != NULL );
   assert( 0 <= nrows && nrows <= sdpi->nlpcons );
   assert( ind != NULL );
   assert( lhs != NULL );
   assert( rhs != NULL );

   for (i = 0; i < nrows; i++)
   {
      assert ( ind[i] >= 0 );
      assert ( ind[i] < sdpi->nlpcons );
      sdpi->lplhs[ind[i]] = lhs[i];
      sdpi->lprhs[ind[i]] = rhs[i];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of LP-rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   )
{
   assert( sdpi != NULL );

   *nlprows = sdpi->nlpcons;

   return SCIP_OKAY;
}

/** gets the number of SDP-Blocks in the SDP */
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of blocks */
   )
{
   assert( sdpi != NULL );

   *nsdpblocks = sdpi->nsdpblocks;

   return SCIP_OKAY;
}

/** gets the number of variables in the SDP */
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   assert( sdpi != NULL );

   *nvars = sdpi->nvars;   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP constraint matrices */
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP constraint matrcies */
   )
{
   assert( sdpi != NULL );

   *nnonz = sdpi->sdpnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   )
{
   assert( sdpi != NULL );

   *nnonz = sdpi->sdpconstnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP Matrix */
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   )
{
   assert( sdpi != NULL );

   *nnonz = sdpi->lpnnonz;

   return SCIP_OKAY;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< pointer to store objective coefficients (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstvar >= 0 );
   assert( firstvar <= lastvar );
   assert( lastvar < sdpi->nvars);
   assert( vals != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++) /*lint !e834*/
      vals[i] = sdpi->obj[firstvar + i]; /*lint !e679*/

   return SCIP_OKAY;
}

/** gets current variable bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< pointer to store lower bound values, or NULL (memory of size lastvar - firstvar + 1 needs to be allocated) */
   SCIP_Real*            ubs                 /**< pointer to store upper bound values, or NULL (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstvar >= 0 );
   assert( firstvar <= lastvar );
   assert( lastvar < sdpi->nvars);
   assert( lbs != NULL );
   assert( ubs != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++) /*lint !e834*/
   {
      if (lbs != NULL)
         lbs[i] = sdpi->lb[firstvar + i]; /*lint !e679*/
      if (ubs != NULL)
         ubs[i] = sdpi->ub[firstvar + i]; /*lint !e679*/
   }
   return SCIP_OKAY;
}

/** gets current left hand sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetLhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss                /**< pointer to store left hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstrow >= 0 );
   assert( firstrow <= lastrow );
   assert( lastrow < sdpi->nlpcons);
   assert( lhss != NULL );

   for (i = 0; i < lastrow - firstrow + 1; i++) /*lint !e834*/
      lhss[firstrow + i] = sdpi->lplhs[i]; /*lint !e679*/

   return SCIP_OKAY;
}

/** gets current right hand sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< pointer to store right hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   )
{
   int i;

   assert( sdpi != NULL );
   assert( firstrow >= 0 );
   assert( firstrow <= lastrow );
   assert( lastrow < sdpi->nlpcons);
   assert( rhss != NULL );

   for (i = 0; i < lastrow - firstrow + 1; i++) /*lint !e834*/
      rhss[firstrow + i] = sdpi->lprhs[i]; /*lint !e679*/

   return SCIP_OKAY;
}


/**@} */



/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP, as start optionally a starting point for the solver may be given, if it is NULL, the solver will start from scratch */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            start,              /**< NULL or a starting point for the solver, this should have length nvars */
   int*                  totalsdpiterations, /**< the number of sdpiterations needed will be added to the int this points to */
   SCIP_Bool             enforceslatercheck  /**< always check for Slater condition in case the problem could not be solved and printf the solution
                                                  of this check */
   )
{
   int block;
   int sdpconstnnonz;
   int* sdpconstnblocknonz;
   int** sdpconstrow;
   int** sdpconstcol;
   SCIP_Real** sdpconstval;
   int** indchanges;
   int* nremovedinds;
   int newiterations;
   SCIP_Real* lplhsafterfix;
   SCIP_Real* lprhsafterfix;
   int* rowsnactivevars;
   SCIP_Bool fixingfound;
   int nactivelpcons;
   int* blockindchanges;
   int nremovedblocks;

   assert ( sdpi != NULL );
   assert ( totalsdpiterations != NULL );

   SCIPdebugMessage("Forwarding SDP %d to solver!\n", sdpi->sdpid);

   sdpconstnblocknonz = NULL;
   sdpconstrow = NULL;
   sdpconstcol = NULL;
   sdpconstval = NULL;
   indchanges = NULL;
   nremovedinds = NULL;
   nremovedblocks = 0;

   /* allocate memory for computing the constant matrix after fixings and finding empty rows and columns, this is as much as might possibly be
    * needed, this will be shrinked again before solving */
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstnblocknonz, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstrow, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstcol, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &sdpconstval, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &indchanges, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &nremovedinds, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &blockindchanges, sdpi->nsdpblocks) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &lplhsafterfix, sdpi->nlpcons) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &lprhsafterfix, sdpi->nlpcons) );
   BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &rowsnactivevars, sdpi->nlpcons) );

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      sdpconstrow[block] = NULL;
      sdpconstcol[block] = NULL;
      sdpconstval[block] = NULL;
      indchanges[block] = NULL;
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(indchanges[block]), sdpi->sdpblocksizes[block]) );
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
      BMS_CALL( BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz) ); /*lint !e776*/
   }

   /* compute the lplphss and lprhss, detect empty rows and check for additional variable fixings caused by boundchanges from
    * lp rows with a single active variable */
   do
   {
      fixingfound = FALSE;
      SCIP_CALL( computeLpLhsRhsAfterFixings(sdpi, &nactivelpcons, lplhsafterfix, lprhsafterfix, rowsnactivevars, &fixingfound) );
   }
   while (fixingfound);

   /* initialize sdpconstnblocknonz */
   for (block = 0; block < sdpi->nsdpblocks; block++)
      sdpconstnblocknonz[block] = sdpi->sdpnnonz + sdpi->sdpconstnnonz;

   SCIP_CALL( compConstMatAfterFixings(sdpi, &sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval) );

   /* shrink the constant arrays after the number of fixed nonzeros is known */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      assert ( sdpconstnblocknonz[block] <= sdpi->sdpnnonz + sdpi->sdpconstnnonz ); /* otherwise the memory wasn't sufficient,
                                                                                     * but we allocated more than enough */
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
      BMS_CALL( BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpi->sdpnnonz + sdpi->sdpconstnnonz, sdpconstnblocknonz[block]) ); /*lint !e776*/
   }

   SCIP_CALL (findEmptyRowColsSDP(sdpi, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval, indchanges, nremovedinds, blockindchanges, &nremovedblocks) );

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("SDP %d not given to solver, as infeasibility was detected during presolving!\n", sdpi->sdpid++);
      SCIP_CALL( SCIPsdpiSolverIncreaseCounter(sdpi->sdpisolver) );

      sdpi->solved = TRUE;
   }
   else
   {
      if ( sdpi->slatercheck )
      {
         SCIP_Real objval;
         SCIP_Bool origfeas;

         /* we solve the problem with a slack variable times identity added to the constraints and trying to minimize this slack variable r, if we are
          * still feasible for r < feastol, then we have an interior point with smallest eigenvalue > feastol, otherwise the Slater condition is harmed */
         SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, 1.0, FALSE, FALSE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
               sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
               sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
               sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
               rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, &origfeas) );

         /* if we didn't succeed, then probably the primal problem is troublesome */
         if ( (! SCIPsdpiSolverIsOptimal(sdpi->sdpisolver)) && (! SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver)) )
            printf("Unable to check Slater condition for dual problem of SDP %d, could mean that the Slater conidition for the primal problem"
                  " is not fullfilled.\n", sdpi->sdpid);
         else
         {
            SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objval) );

            if ( objval < - sdpi->feastol )
               SCIPdebugMessage("Slater condition for SDP %d is fullfilled for dual problem with smallest eigenvalue %f.\n", sdpi->sdpid, -1.0 * objval);
            else
               printf("Slater condition for SDP %d not fullfilled for dual problem as smallest eigenvalue was %f, expect numerical trouble.\n",
                  sdpi->sdpid, -1.0 * objval);
         }
      }

      /* try to solve the problem */
      SCIP_CALL( SCIPsdpiSolverLoadAndSolve(sdpi->sdpisolver, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
            sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
            sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
            sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
            sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
            rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start) );

      sdpi->penalty = FALSE;
      sdpi->solved = TRUE;

      /* if the solver didn't produce a satisfactory result, we have to try with penalty formulations */
      if ( ! SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) )
      {
         SCIP_Real penaltyparam;
         SCIP_Bool feasorig;

         penaltyparam = 1.0;
         feasorig = FALSE;

         /* if we didn't converge, first try to check feasibility with a penalty formulation */
         SCIPdebugMessage("Solver did not produce an acceptable result, trying SDP %d again with penalty formulation\n", sdpi->sdpid);

         SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, 1.0, FALSE, TRUE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
               sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
               sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
               sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
               sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
               rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, &feasorig) );

         /* if the solver converged and the solution is feasible for our original problem, the problem is feasible and we can continue to search for an
          * optimal solution, if the solver converged but the solution isn't feasible for our original problem, the problem is infeasible, if we still
          * didn't converge, we are out of luck */
         if ( SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && feasorig)
         {
            /* increase the penalty parameter until we are able to solve the problem and get a solution that is feasible for our original problem and
             * the problem isn't unbounded (this can be caused by a too small penalty parameter) or the penalty parameter gets too large */
            do
            {
               SCIPdebugMessage("Solver did not produce an acceptable result, trying SDP %d again with penaltyparameter %f\n", sdpi->sdpid, penaltyparam);

               SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, penaltyparam, TRUE, TRUE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
                     sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
                     sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
                     sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
                     sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
                     rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, &feasorig) );

               penaltyparam *= 10.0;
            }
            while ( (! SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) || ! feasorig || SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver))  &&
                     ! SCIPsdpiSolverIsGEMaxPenParam(sdpi->sdpisolver, penaltyparam) );

            /* check if we were able to solve the problem in the end */
            if ( SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && feasorig )
               sdpi->penalty = TRUE;
            else
            {
               /* TODO: Do we want to give the computed lower bound? Do we want to tell that we showed that the problem is feasible? */
               sdpi->solved = FALSE;
               sdpi->penalty = TRUE;
            }
         }
         else if ( SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver) && ! feasorig )
         {
            SCIPdebugMessage("Problem was found to be infeasible using a penalty formulation \n");
            sdpi->infeasible = TRUE;
            sdpi->penalty = TRUE;
         }
         else
         {
            SCIPdebugMessage("SDP-Solver could not solve the problem even after using a penalty formulation \n");
            sdpi->solved = FALSE;
            sdpi->penalty = TRUE;
         }
         /* if we still didn't succeed and enforceslatercheck was set, we finally test for the Slater condition to give a reason for failure */
         if ( sdpi->solved == FALSE && enforceslatercheck)
         {
            SCIP_Real objval;
            SCIP_Bool origfeas;

            /* we solve the problem with a slack variable times identity added to the constraints and trying to minimize this slack variable r, if we are
             * still feasible for r < feastol, then we have an interior point with smallest eigenvalue > feastol, otherwise the Slater condition is harmed */
            SCIP_CALL( SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, 1.0, FALSE, FALSE, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
                  sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpnblockvars, sdpconstnnonz,
                  sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
                  sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdpvar, sdpi->sdprow, sdpi->sdpcol,
                  sdpi->sdpval, indchanges, nremovedinds, blockindchanges, nremovedblocks, nactivelpcons, sdpi->nlpcons, lplhsafterfix, lprhsafterfix,
                  rowsnactivevars, sdpi->lpnnonz, sdpi->lprow, sdpi->lpcol, sdpi->lpval, start, &origfeas) );

            /* if we didn't succeed, then probably the primal problem is troublesome */
            if ( (! SCIPsdpiSolverIsOptimal(sdpi->sdpisolver)) && (! SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver)) )
               printf("SDP-solver could not solve root node relaxation, unable to check Slater condition for dual problem of SDP %d, could mean that the "
                     "Slater conidition for the primal problem is not fullfilled.\n", sdpi->sdpid);
            else
            {
               SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, &objval) );

               if ( objval < - sdpi->feastol )
                  printf("SDP-solver could not solve root node relaxation even though the Slater condition is fullfilled for the dual problem with smallest eigenvalue %f.\n",
                           -1.0 * objval);
               else
                  printf("SDP-solver could not solve root node relaxation, Slater condition is not fullfilled for the dual problem as smallest eigenvalue was %f.\n",
                        -1.0 * objval);
            }
         }
         else if ( sdpi->solved == FALSE )
            printf("Numerical trouble\n");
      }
   }

   /* empty the memory allocated here */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstval[block]), sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstcol[block]), sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpconstrow[block]), sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(indchanges[block]), sdpi->sdpblocksizes[block]);
   }
   BMSfreeBlockMemoryArray(sdpi->blkmem, &rowsnactivevars, sdpi->nlpcons);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &lprhsafterfix, sdpi->nlpcons);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &lplhsafterfix, sdpi->nlpcons);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &blockindchanges, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &nremovedinds, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &indchanges, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstval, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstcol, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstrow, sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sdpconstnblocknonz, sdpi->nsdpblocks);

   /* add the iterations needed to solve this SDP */
   if ( ! sdpi->infeasible )
   {
      SCIP_CALL( SCIPsdpiSolverGetIterations(sdpi->sdpisolver, &newiterations) );
      *totalsdpiterations += newiterations;
   }

   sdpi->sdpid++;

   return SCIP_OKAY;
}




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );

   return sdpi->solved;
}

/** returns whether the original problem was solved, if SCIPsdpiWasSolved = true and SCIPsdpiSolvedOrig = false, then a penalty formulation was solved */
SCIP_Bool SCIPsdpiSolvedOrig(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );

   return ( ! sdpi->penalty );
}

/** returns true if the solver could determine whether the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver doesn't know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be used */
SCIP_Bool SCIPsdpiFeasibilityKnown(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if (sdpi->infeasible)
      return TRUE;

   return SCIPsdpiSolverFeasibilityKnown(sdpi->sdpisolver);
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store the primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store the dual feasibility status */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      *dualfeasible = FALSE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetSolFeasibility(sdpi->sdpisolver, primalfeasible, dualfeasible) );

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal unboundedness not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalUnbounded(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be primal infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalInfeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be primal feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, primal feasibility not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsPrimalFeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, dual unboundedness not available\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsDualUnbounded(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing\n");
      return TRUE;
   }

   return SCIPsdpiSolverIsDualInfeasible(sdpi->sdpisolver);
}

/** returns TRUE iff SDP is proven to be dual feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsDualFeasible(sdpi->sdpisolver);
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiIsConverged(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, this counts as converged.\n");
      return TRUE;
   }

   return SCIPsdpiSolverIsConverged(sdpi->sdpisolver);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no objective limit available.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsObjlimExc(sdpi->sdpisolver);
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no iteration limit available.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsIterlimExc(sdpi->sdpisolver);
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no time limit available.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsTimelimExc(sdpi->sdpisolver);
}

/** returns the internal solution status of the solver, which has the following meaning:<br>
 * -1: solver wasn't started<br>
 *  0: converged<br>
 *  1: infeasible start<br>
 *  2: numerical problems<br>
 *  3: objective limit reached<br>
 *  4: iteration limit reached<br>
 *  5: time limit reached<br>
 *  6: user termination<br>
 *  7: other */
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );

   if ( ! sdpi->solved )
   {
      SCIPdebugMessage("Problem wasn't solved yet.\n");
      return -1;
   }
   else if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no internal status available.\n");
      return 0;
   }

   return SCIPsdpiSolverGetInternalStatus(sdpi->sdpisolver);
}

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, therefore there is no optimal solution.\n");
      return FALSE;
   }

   return SCIPsdpiSolverIsOptimal(sdpi->sdpisolver);
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 * that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiIsAcceptable(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED_BOOL(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, this is acceptable in a B&B context.\n");
      return TRUE;
   }

   return SCIPsdpiSolverIsAcceptable(sdpi->sdpisolver);
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   )
{
   assert( sdpi != NULL );
   assert( objval != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no objective value available.\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetObjval(sdpi->sdpisolver, objval) );

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return the needed length and
 *  a debug message */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dualsol vector, must be 0 if dualsol is NULL, if this is less than the number
                                               *  of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   )
{
   assert( sdpi != NULL );
   assert( dualsollength != NULL );
   assert( *dualsollength == 0 || dualsol != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no solution available.\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetSol(sdpi->sdpisolver, objval, dualsol, dualsollength) );

   return SCIP_OKAY;
}

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem, the last input should specify the length
 *  of the arrays, if this is less than the number of variables, the needed length will be returned and a debug message thrown
 *  @note if a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable */
SCIP_RETCODE SCIPsdpiGetPrimalBoundVars(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars<br>
                                                  output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   )
{
   assert( sdpi != NULL );
   assert( lbvars != NULL );
   assert( ubvars != NULL );
   assert( arraylength != NULL );
   assert( *arraylength >= 0 );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no primal variables available.\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetPrimalBoundVars(sdpi->sdpisolver, lbvars, ubvars, arraylength) );

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( sdpi->infeasible )
   {
      SCIPdebugMessage("Problem was found infeasible during preprocessing, no iterations needed.\n");
      *iterations = 0;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPsdpiSolverGetIterations(sdpi->sdpisolver, iterations) );

   return SCIP_OKAY;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP solver */
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   assert( sdpi != NULL  );

   return SCIPsdpiSolverInfinity(sdpi->sdpisolver);
}

/** checks if given value is treated as (plus or minus) infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   )
{
   assert( sdpi != NULL );

   return ((val <= -SCIPsdpiInfinity(sdpi)) || (val >= SCIPsdpiInfinity(sdpi)));
}

/** returns highest penalty parameter to be used */
SCIP_Real SCIPsdpiMaxPenParam(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   assert( sdpi != NULL );

   return SCIPsdpiSolverMaxPenParam(sdpi->sdpisolver);
}

/** checks if given value is greater or equal to the highest penalty parameter to be used */
SCIP_Bool SCIPsdpiIsGEMaxPenParam(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be compared to maximum penalty parameter */
   )
{
   assert( sdpi != NULL );

   return ((val <= -SCIPsdpiMaxPenParam(sdpi)) || (val >= SCIPsdpiMaxPenParam(sdpi)));
}

/** gets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );
   assert( dval != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_EPSILON:
      *dval = sdpi->epsilon;
      break;
   case SCIP_SDPPAR_FEASTOL:
      *dval = sdpi->feastol;
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

#ifndef NDEBUG
   {
      SCIP_Real val;
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, &val) );
      assert( REALABS(*dval - val) < sdpi->epsilon );
   }
#endif

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
   case SCIP_SDPPAR_EPSILON:
      sdpi->epsilon = dval;
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_FEASTOL:
      sdpi->feastol = dval;
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   case SCIP_SDPPAR_OBJLIMIT:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetRealpar(sdpi->sdpisolver, type, dval) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

#ifndef NDEBUG
   {
      SCIP_Real val;
      SCIP_CALL_PARAM( SCIPsdpiSolverGetRealpar(sdpi->sdpisolver, type, &val) );
      assert( REALABS(dval - val) < sdpi->epsilon );
   }
#endif

   return SCIP_OKAY;
}

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );
   assert( ival != NULL );

   switch( type )/*lint --e{788}*/
   {
#if 0
   case SCIP_SDPPAR_THREADS:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetIntpar(sdpi->sdpisolver, type, ival) );
      break;
#endif
   case SCIP_SDPPAR_SDPINFO:
      SCIP_CALL_PARAM( SCIPsdpiSolverGetIntpar(sdpi->sdpisolver, type, ival) );
      break;
   case SCIP_SDPPAR_SLATERCHECK:
      *ival = (int) sdpi->slatercheck;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

#ifndef NDEBUG
   {
      if ( type != SCIP_SDPPAR_SLATERCHECK )
      {
         int val;
         SCIP_CALL_PARAM( SCIPsdpiSolverGetIntpar(sdpi->sdpisolver, type, &val) );
         assert( *ival == val );
      }
   }
#endif

   return SCIP_OKAY;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   assert( sdpi != NULL );
   assert( sdpi->sdpisolver != NULL );

   switch( type )/*lint --e{788}*/
   {
#if 0
   case SCIP_SDPPAR_THREADS:
      SCIP_CALL_PARAM( SCIPsdpiSolverSetIntpar(sdpi->sdpisolver, type, ival) );
      break;
#endif
   case SCIP_SDPPAR_SDPINFO:
      assert( ival == 0 || ival == 1 ); /* this is a boolean parameter */
      SCIP_CALL_PARAM( SCIPsdpiSolverSetIntpar(sdpi->sdpisolver, type, ival) );
      break;
   case SCIP_SDPPAR_SLATERCHECK:
      assert( ival == 0 || ival == 1 ); /* this is a boolean parameter */
      sdpi->slatercheck = (SCIP_Bool) ival;
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }

#ifndef NDEBUG
   {
      if ( type != SCIP_SDPPAR_SLATERCHECK )
      {
         int val;
         SCIP_CALL_PARAM( SCIPsdpiSolverGetIntpar(sdpi->sdpisolver, type, &val) );
         assert( ival == val );
      }
   }
#endif

   return SCIP_OKAY;
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{/*lint --e{715}*/
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_LPERROR;
}

/**@} */
