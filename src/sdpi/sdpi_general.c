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

#define SCIP_DEBUG
#define SCIP_MORE_DEBUG

/**@file   sdpi_dsdp.c
 * @brief  interface for dsdp
 * @author Tristan Gally
 */

#include <assert.h>

#include "sdpi/sdpi.h"
#include "sdpi/sdpi_general.h"
#include "scipsdp/SdpVarfixer.h"

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */

struct SCIP_SDPi
{
   SCIP_SDPISOLVER*      sdpisolver;         /**< pointer to the interface for the SDP Solver */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   int                   nvars;              /**< number of variables */
   int                   nvarsorig;          /**< number of variables in original problem */
   SCIP_Real*            obj;                /**< objective function values of variables */
   SCIP_Real*            objorig;            /**< objective function values of variables in original problem */
   SCIP_Real*            lb;                 /**< lower bounds of variables */
   SCIP_Real*            lborig;             /**< lower bounds of variables in original problem */
   SCIP_Real*            ub;                 /**< upper bounds of variables */
   SCIP_Real*            uborig;             /**< upper bounds of variables in original problem */
   int                   nsdpblocks;         /**< number of SDP-blocks */
   int                   nsdpblocksorig;     /**< number of SDP-blocks in original problem */
   int*                  sdpblocksizes;      /**< sizes of the SDP-blocks */
   int*                  sdpblocksizesorig;  /**< sizes of the SDP-blocks in original problem */

   /* constant SDP data: */
   int                   sdpconstnnonz;      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int                   sdpconstnnonzorig;  /**< number of nonzero elements in the constant matrices of the SDP-Blocks in original problem */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                               *  number of entries  of sdpconst row/col/val [i] */
   int*                  sdpconstnblocknonzorig;/**< number of nonzeros for each variable in the constant part of the original problem,
                                                     *  also the i-th entry gives the number of entries  of sdpconst row/col/val orig [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstroworig;    /**< pointers to row-indices for each block in the original problem */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   int**                 sdpconstcolorig;    /**< pointers to column-indices for each block in the original problem */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */
   SCIP_Real**           sdpconstvalorig;    /**< pointers to the values of the nonzeros for each block in the original problem */

   /* non-constant SDP data: */
   int                   sdpnnonz;           /**< number of nonzero elements in the SDP-constraint matrices */
   int                   sdpnnonzorig;       /**< number of nonzero elements in the SDP-constraint matrices in original problem */
   int*                  sdpnblockvarnonz;   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i * nvars + j] */
   int*                  sdpnblockvarnonzorig;/**< entry i * nvars + j gives the number of nonzeros for block i and variable j in the original
                                                *  problem, this is exactly the number of entries of sdp row/col/val orig[i * nvars + j] */
   int**                 sdprow;             /**< pointer to the row-indices for each block and variable */
   int**                 sdproworig;         /**< pointer to the row-indices for each block and variable in the original problem */
   int**                 sdpcol;             /**< pointer to the column-indices for each block and variable */
   int**                 sdpcolorig;         /**< pointer to the column-indices for each block and variable in the original problem */
   SCIP_Real**           sdpval;             /**< pointer to the values of the nonzeros for each block and variable */
   SCIP_Real**           sdpvalorig;         /**< pointer to the values of the nonzeros for each block and variable in the original problem */

   /* lp data: */
   int                   nlpcons;            /**< number of LP-constraints */
   int                   nlpconsorig;        /**< number of LP-constraints in original problem */
   SCIP_Real*            lprhs;              /**< right hand sides of LP rows */
   SCIP_Real*            lprhsorig;          /**< right hand sides of LP rows in original problem */
   int                   lpnnonz;            /**< number of nonzero elements in the LP-constraint matrix */
   int                   lpnnonzorig;        /**< number of nonzero elements in the LP-constraint matrix in original problem */
   int*                  lprowind;           /**< row-index for each entry in lpval-array */
   int*                  lprowindorig;       /**< row-index for each entry in lpval-array in original problem */
   int*                  lpcolind;           /**< column-index for each entry in lpval-array */
   int*                  lpcolindorig;       /**< column-index for each entry in lpval-array in original problem */
   SCIP_Real*            lpval;               /**< values of LP-constraint matrix entries */
   SCIP_Real*            lpvalorig;           /**< values of LP-constraint matrix entries in original problem */

   /* other data */
   int                   sdpid;              /**< counter for the number of SDPs solved */
   int*                  sdptoinputmapper;   /**< entry i gives the original index of the i-th variable of the current sdp */
   int*                  inputtosdpmapper;   /**< entry i gives the current sdp index of the i-th input variable (or -1 if it is fixed) */
   SCIP_Real*            fixedvals;          /**< values the variables are fixed to, this is indexed by original variables */
   SCIP_Bool*            varstochgfixing;    /**< entry i is TRUE, if variable i is fixed and should be unfixed or is unfixed and should be fixed */
   SCIP_Bool             fixingstodo;        /**< TRUE if any entry of varstochgfixing is TRUE, FALSE if no fixings or unfixings have to be done */
};

static double epsilon    = 1e-6;             /**< this is used for checking if primal and dual objective are equal */

/*
 * Local Functions
 */

/**
 * For given row and column (i,j) checks if i >= j, so that i and j give a position in the lower
 * triangular part, otherwise i and j will be switched. This function will be called whenever a position in a symmetric matrix
 * is given, to prevent problems if position (i,j) is given but later (j,i) should be changed.
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
/**
 * Test if a lower bound lb really is smaller than an upper bound ub, meaning that lb <= ub - epsilon */
static
SCIP_Bool isFixed(
   SCIP_Real             lb,                 /**< lower bound */
   SCIP_Real             ub                  /**< upper bound */
   )
{
   assert( lb < ub + epsilon );

   return ( REALABS(ub-lb) <= epsilon);
}
#else
#define isFixed(lb,ub,epsilon) (REALABS(ub-lb) <= epsilon)
#endif


/** fix and unfix variables according to the varstochgfixing-array */
static
SCIP_RETCODE performFixing(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   )
{
   int i;
   int var;
   int block;
   int ind;
   int sdpind;
   int begvarblockind;
   int endindex;
   int nfixedvars;   /* this is a net value, if it is negative than more were unfixed than fixed */
   int nfixedsdpnonz;
   int nblockfixedsdpnonz; /* number of fixed or unfixed nonzeros */
   int nfixedlpnonz;
   int naddedblocknonz;
   int naddedlpnonz;
   int nblockvarnonz;
   int nleftshiftsvars; /* how many spots do the variable-arrays have to be moved to the left */
   int nleftshiftsdpnonz; /* "-----------------" sdpnonzero-arrays "------------------------" */
   int nleftshiftconstnonz; /* "----------------" sdpconstnonzero-arrays "------------------" */
   int nleftshiftsdparrays; /* "--------------" sdpbegvarblock-array "--------------------" */
   int nleftshiftlp; /* "-----------------------" lpnonz-array "----------------------------" */
   int nrightshiftsvarsdone; /* how many spots have the variable-arrays already moved to the right */
   int nrightshiftsdpnonzdone; /* "-------------------" sdpnonzero-arrays "----------------------" */
   int nrightshiftconstnonzdone; /* "-----------------" sdpconstnonzero-arrays "-----------------" */
   int nrigthshiftsdparraysdone; /* "---------------" sdpbegvarblock-array "-------------------" */
   int sdpind;
   int* colstoadd;   /* these arrays are used to change the const array after fixing/unfixing variables, each fixing/unfixing of a variable leads to
                      * a change in the value of an entry of the constmatrix, these are saved here to then check if this entry exists and should be
                      * changed, or if a new nonzero entry needs to be created, both fixings and unfixings are saved togehter, they only differ in the
                      * sign of the values */
   int* rowstoadd;
   SCIP_Real* valstoadd;
   int vararraylength;
   int sdpnonzarraylength;
   int sdpconstarraylength;
   int sdparraylength;
   int lastactivesdpvar;
   int nblockvarnonz;
   int newind;
   SCIP_Bool resetshift;
   SCIP_Bool fixedvar;
   SCIP_Bool unfixedvar;


   assert ( sdpi != NULL );

   /* check if there is anything to do */
   if (!(sdpi->fixingstodo))
      return SCIP_OKAY;

   nfixedvars = 0;
   nfixedsdpnonz = 0;
   nleftshiftvar = 0;
   nleftshiftsdpnonz = 0;
   nleftshiftconstnonz = 0;
   nleftshiftsdparrays = 0;
   nleftshiftlp = 0;
   nrightshiftvardone = 0;
   nrightshiftsdpnonzdone = 0;
   nrightshiftconstnonzdone = 0;
   nrigthshiftsdparraysdone = 0;
   vararraylength = sdpi->nvars;
   sdpnonzarraylength = sdpi->sdpnnonz;
   sdpconstarraylength = sdpi->sdpconstnnonz;
   lastactivesdpvar = -1;
   naddedlpnonz = 0;
   fixedvar = FALSE;
   unfixedvar = FALSE;

   /* add/remove obj/lb/ub etc. */
   for (var = 0; var < sdpi->nvarsorig; var++)
   {
      if (sdpi->varstochgfixing[var])
      {
         /* the status of this variable changes */

         if (sdpi->inputtosdpmapper[var] != -1)
         {
            /* this variable was active, so it will now be fixed */

            fixedvar = TRUE; /* remember that a variable was fixed */

            assert ( !(IsFixed(sdpi->lborig[var], sdpi->uborig[var]) ) );

            nfixedvars++;
            nleftshiftvar++;
            fixedvars[var] = sdpi->lb[var];
         }
         else
         {
            /* this variable was fixed, so it will now become active again */

            unfixedvar = TRUE; /* remember that a variable was unfixed */

            if (nleftshiftvar >= 1)
            {
               /* because of earlier fixings there is an empty spot in the arrays where this can be inserted */

               assert ( !(IsFixed(sdpi->lborig[var], sdpi->uborig[var])) );

               sdpi->obj[lastactivesdpvar + 1] = sdpi->objorig[var];
               sdpi->lb[lastactivesdpvar + 1] = sdpi->lborig[var];
               sdpi->ub[lastactivesdpvar + 1] = sdpi->uborig[var];

               /* one of the empty spots was taken */
               nleftshiftvar--;
               nfixedvars--;
               lastactivesdpvar++; /* another spot in the sdp arrays was taken */
            }
            else
            {
               assert ( !(IsFixed(sdpi->lborig[var], sdpi->uborig[var])) );

               /* increase the length of the arrays */
               BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), arraylength, arraylength + 1);
               BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), arraylength, arraylength + 1);
               BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), arraylength, arraylength + 1);
               arraylength++;

               /* shift all later variables in the arrays to the right */
               for (i = arraylength; i > lastactivesdpvar; i--)   /* going downwards to not overwrite later entries */
               {
                  sdpi->obj[i + 1] = sdpi->obj[i];
                  sdpi->lb[i + 1] = sdpi->lb[i];
                  sdpi->ub[i + 1] = sdpi->ub[i];
               }

               nrightshiftvardone++;
               nfixedvars--;

               sdpi->obj[lastactivesdpvar + 1] = sdpi->objorig[var];
               sdpi->lb[lastactivesdpvar + 1] = sdpi->lborig[var];
               sdpi->ub[lastactivesdpvar + 1] = sdpi->uborig[var];

               lastactivesdpvar++; /* another spot in the sdp arrays was taken */
            }
         }
      }
      else if (sdpi->inputtosdpmapper[var] != -1 && nleftshiftvars > 0)
      {
         /* this variable was active and stays active, so its entries have to be shifted according to fixings already done */

         lastactivesdpvar = sdpi->inputtosdpmapper[var] - nleftshiftvars  + nrightshiftvarsdone; /* because of the rightshifts after insertig unfixed
                                                                                                  * variables both the current and the target position
                                                                                                  * are moved that many positions to the right, the target
                                                                                                  * position is then redurced by the number of leftshifts
                                                                                                  * that are still to be done*/

         sdpi->obj[lastactivesdpvar] = sdpi->obj[lastactivesdpvar + nleftshiftvars];
         sdpi->lb[lastactivesdpvar] = sdpi->lb[lastactivesdpvar + nleftshiftvars];
         sdpi->ub[lastactivesdpvar] = sdpi->ub[lastactivesdpvar + nleftshiftvars];
      }
      else if (sdpi->inputtosdpmapper[var] != -1)
      {
         /* this variable was active and stays active */
         lastactivesdpvar = sdpi->inputtosdpmapper[var];
      }
      /* the else case would be is fixed and stays fixed, in this case nothing has to be done */
   }

   /* shrink the arrays according to the number of fixings (for the added entries they were already enlarge, so we only need to substract) */
   if (nleftshiftvars > 0)
   {
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), arraylength, arraylength - nleftshiftvars);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), arraylength, arraylength - nleftshiftvars);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), arraylength, arraylength - nleftshiftvars);
      arraylength -= leftshiftvars; // TODO: is this used again ?!? or is updating meaningless
   }

   /* prepare arrays to save the fixed/unfixed nonzeros to add them to the constant part */
   BMSallocBlockMemoryArray(sdpi->blkmem, &colstoadd, sdpi->sdpnnonzorig);
   BMSallocBlockMemoryArray(sdpi->blkmem, &rowstoadd, sdpi->sdpnnonzorig);
   BMSallocBlockMemoryArray(sdpi->blkmem, &valstoadd, sdpi->sdpnnonzorig);

   lastactivesdpvar = -1;

   /* SDP nonzeros */
   for (block = 0; block < sdpi->nsdpblocksorig; block++)
   {

      nblockfixedsdpnonz = 0;
      nblockaddednonz = 0;

      for (var = 0; var < sdpi->nvarsorig; var++)
      {
         if (sdpi->varstochgfixing[var])
         {
            /* the status of this variable changes */
            if (sdpi->inputtosdpmapper[var] != -1)
            {
               /* this variable was active, so it will now be fixed */

               sdpind = sdpi->inputtosdpmapper[var];

               /* copy the nonzeros to the toadd arrays */
               sdparrayind = block * sdpi->nvars + sdpind - nleftshiftsdparrays + nrigthshiftsdparraysdone;
               for (i = 0; i < sdpi->sdpnblockvarnonz[sdparrayind]; i++)
               {
                  colstoadd[nblockfixedsdpnonz] = sdpi->sdpcolind[i];
                  rowstoadd[nblockfixedsdpnonz] = sdpi->sdprowind[i];
                  valstoadd[nblockfixedsdpnonz] = -sdpi->sdpval[i] * sdpi->lborig[i]; /* minus because of +A_i but -A_0 */
               }

               nleftshiftsdparrays++;
            }
            else
            {
               /* this variable was fixed, so it will now become active again */
               sdparrayind = block * sdpi->nvarsorig + var;
               newind = block * sdpi->nvars + lastactivesdpvar - nleftshiftsdparrays + nrigthshiftsdparraysdone; /* the new index for the sdp arrays */
               lastactivesdpvar++; /* another spot in the sdp-arrays has to be taken */

               if (nleftshiftsdparrays == 0)
               {
                  /* add another entry to the sdp arrays */
                  BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonz), sdparraylength, sdparraylength + 1);
                  BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->col), sdparraylength, sdparraylength + 1);
                  BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->row), sdparraylength, sdparraylength + 1);
                  BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->val), sdparraylength, sdparraylength + 1);
                  sdparraylength++;

                  /* move all later entries one spot to the right */
                  for (i = sdparraylength - 1; i > newind; i--) /* going down to not overwrite something not yet used */
                  {
                     sdpi->sdpnblockvarnonz[i + 1] = sdpi->sdpnblockvarnonz[i];
                     sdpi->col[i + 1] = sdpi->col[i];
                     sdpi->row[i + 1] = sdpi->row[i];
                     sdpi->val[i + 1] = sdpi->val[i];
                  }
               }
               else
                  nleftshiftsdparrays--; /* one of the empty spots is used, so the rest is shifted one spot less */

               /* allocate memory for the nonzeros */
               BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->col[newind]), sdpi->sdpnblockvarnonz[sdparrayind]);
               BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->row[newind]), sdpi->sdpnblockvarnonz[sdparrayind]);
               BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->val[newind]), sdpi->sdpnblockvarnonz[sdparrayind]);

               /* copy the nonzero-entries to the new sdp arrays(just copying the pointers is not sufficient, because we want to be able to
                * change these entries without changing the original ones) and their negative values to the toadd arrays, to substract them
                * from the constant part */
               sdpi->sdpnblockvarnonz[newind] = sdpi->sdpnblockvarnonzorig[sdparrayind];
               for (i = 0; i < sdpi->sdpnblockvarnonzorig[sdparrayind]; i++)
               {
                  sdpi->sdpcol[newind][i] = sdpi->sdpcolorig[sdparrayind][i];
                  colstoadd[nblockfixedsdpnonz] = sdpi->sdpcolorig[sdparrayind][i];
                  sdpi->sdprow[newind][i] = sdpi->sdproworig[sdparrayind][i];
                  rowstoadd[nblockfixedsdpnonz] = sdpi->sdproworig[sdparrayind][i];
                  sdpi->sdpval[newind][i] = sdpi->sdpvalorig[sdparrayind][i];
                  valstoadd[nblockfixedsdpnonz] = +sdpi->sdpvalorig[sdparrayind][i] * sdpi->fixedvals[var]; /* plus = minus * minus because of +A_i but -A_0
                                                                                                             * and we want to decrease the corresponding
                                                                                                             * entry as it is moved to A_i */
                  nblockfixedsdponz++;
               }
            }
         }
         else if (sdpi->inputtosdpmapper[var] != -1 && nleftshiftsdpnonz > 0)
         {
            /* this variable was active and stays active, so the array-entries have to be shifted according to the fixings */

            sdpind = sdpi->inputtosdpmapper[var];

            sdpi->sdpnblockvarnonz[sdpind - nleftshiftsdparrays + nrightshiftsbegvarblockdone] = sdpi->sdpnblockvarnonz[sdpind + nrightshiftsbegvarblockdone];
            sdpi->sdpcol[sdpind - nleftshiftsdparrays + nrightshiftsbegvarblockdone] = sdpi->sdpcol[sdpind + nrightshiftsbegvarblockdone];
            sdpi->sdprow[sdpind - nleftshiftsdparrays + nrightshiftsbegvarblockdone] = sdpi->sdprow[sdpind + nrightshiftsbegvarblockdone];
            sdpi->sdpval[sdpind - nleftshiftsdparrays + nrightshiftsbegvarblockdone] = sdpi->sdpval[sdpind + nrightshiftsbegvarblockdone];

            lastactivesdpvar = sdpind - nleftshiftsdparrays + nrightshiftsbegvarblockdone;
         }
         else if (sdpi->inputtosdpmapper[var] != -1)
         {
            /* this variable was active and stays active, so the array-entries have to be shifted according to the fixings */
            lastactivesdpvar = sdpi->inputtosdpmapper[var];
         }
      }

      /* add the fixed/unfixed nonzeros and shift the constant nonzeros*/
      if (nblockfixedsdpnonz > 0)
      {
         SCIP_CALL( SdpVarfixerMergeArrays(sdpi->blkmem, rowstoadd, colstoadd, valstoadd, nblockfixedsdpnonz, FALSE, 1.0, sdpi->sdpconstrow[block],
                    sdpi->sdpconstcol[block], sdpi->sdpconstval[block], &(sdpi->sdpconstnblocknonz[block])) );
      }
   }

   /* free the arrays used to save the fixed/unfixed nonzeros to add them to the constant part */
   BMSfreeBlockMemoryArray(sdpi->blkmem, &colstoadd, sdpi->sdpnnonzorig);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &rowstoadd, sdpi->sdpnnonzorig);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &valstoadd, sdpi->sdpnnonzorig);

   /* shrink the sdp arrays if possible */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol), sdparraylength, sdparraylength - nleftshiftsdparrays);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow), sdparraylength, sdparraylength - nleftshiftsdparrays);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdparraylength, sdparraylength - nleftshiftsdparrays);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonz), sdparraylength, sdparraylength - nleftshiftsdparrays);

   /* compute the new number of total nonzeros */
   sdpi->sdpnnonz = 0;
   sdpi->sdpconstnnonz = 0;

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      sdpi->sdpconstnnonz += sdpi->sdpconstnblocknonz[block];
      for (var = 0; var < sdpi->nvars; var++)
         sdpi->sdpnnonz += sdpi->sdpnblockvarnonz[block * sdpi->nvars + var];
   }


   /* LP nonzeros */

   /*remove the fixed ones */
   if (fixedvar)
   {
      for (i = 0; i < sdpi->lpnnonz; i++)
      {
         if (sdpi->varstochgfixing[sdpi->lpcolind[i]])
         {
            /* remove this entry and decrease the rhs */
            nleftshiftlp++;
            sdpi->lprhs[sdpi->lprowind[i]] -= sdpi->lpval[i];
         }
         else if (nleftshiftlp > 0)
         {
            /* shift according to fixings */
            sdpi->lprowind[i - nleftshiftlp] = sdpi->lprowind[i];
            sdpi->lpcolind[i - nleftshiftlp] = sdpi->lpcolind[i];
            sdpi->lpval[i - nleftsthiftlp] = sdpi->lpval[i];
         }
      }
   }

   /* add the lp entries of unfixed variables */

   /* allocate memory for the insertions, this will be adjusted again later */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonzorig - nleftshiftlp);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonzorig - nleftshiftlp);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonzorig - nleftshiftlp);

   /* iterate over all original lp nonzeros to find those belonging to unfixed vars */
   if (unfixedvar)
   {
      for (i = 0; i < sdpi->lpnnonzorig; i++)
      {
         if (sdpi->varstochgfixing[sdpi->lpcolindori[i]] && sdpi->inputtosdpmapper[sdpi->lpcolindori[i]] == -1)
         {
            /* this lp nonzero has to be unfixed */
            sdpi->lprowind[sdpi->lpnnonz - nleftshiftlp + naddedlpnonz] = sdpi->lprowindorig[i];
            sdpi->lpcolind[sdpi->lpnnonz - nleftshiftlp + naddedlpnonz] = sdpi->lpcolindorig[i];
            sdpi->lpval[sdpi->lpnnonz - nleftshiftlp + naddedlpnonz] = sdpi->lpval[i];
            naddedlpnonz++;
         }
      }
   }

   /* adjust lp-arraylength */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonzorig - nleftshiftlp, sdpi->lpnnonz - nleftshiftlp + naddedlpnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonzorig - nleftshiftlp, sdpi->lpnnonz - nleftshiftlp + naddedlpnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonzorig - nleftshiftlp, sdpi->lpnnonz - nleftshiftlp + naddedlpnonz);

   sdpi->lpnnonz = sdpi->lpnnonz - nleftshiftlp + naddedlpnonz;

   /* update the varmappers and set varstochgfixing to false as the fixings are done*/
   nleftshiftvars = 0;
   nrightshiftvarsdone = 0;
   arraylength = sdpi->nvars;
   lastactivesdpvar = -1;

   for (var = 0; var < sdpi->nvarsorig; var++)
   {
      if (sdpi->varstochgfixing[var])
      {
         /* the status of this variable changes */
         if (sdpi->inputtosdpmapper[var] != -1)
         {
            /* this variable was active, so it will now be fixed */
            nleftshiftvar++;

            sdpi->inputtosdpmapper[var] = -1;
         }
         else
         {
            /* this variable was fixed, so it will now become active again */
            if (nleftshiftvars >= 1)
            {
               /* because of earlier fixing there are enough empty spots to insert this variable into sdptoinputmapper */
               sdpi->sdptoinputmapper[lastactivesdpvar + 1] = var;
               sdpi->inputtosdpmapper[var] = lastactivesdpvar + 1;

               /* one of the empty spots was taken */
               nleftshiftvar--;
               lastactivesdpvar++; /* another position in the sdp-arrays was filled */
            }
            else
            {
               BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdptoinputmapper), arraylength, arraylength + 1);
               arraylength++;

               /* move all entries behind this one spot to the right */
               for (i = arraylength; i > lastactivesdpvar; i--)   /* this has to go downwards to not overwrite folowind entries */
               {
                  sdpi->sdptoinputmapper[i + 1] = lastactivesdpvar[i];
               }

               sdpi->sdptoinputmapper[lastactivesdpvar + 1] = var;;
               sdpi->inputtosdpmapper[var] = lastactivesdpvar + 1;

               lastactivesdpvar++;
            }
         }
      }
      else if (sdpi->inputtosdpmapper[var] != -1 && nleftshiftvars > 0)
      {
         lastactivesdpvar = sdpi->inputtosdpmapper[var] - nleftshiftvars + nrightshiftvarsdone;

         /* this variable was active and stays active, so its entries have to be shifted according to fixings already done */
         sdpi->sdptoinputmapper[lastactivesdpvar] = sdpi->sdptoinputmapper[lastactivevar + nleftshiftvars];
         sdpi->inputtosdpmapper[var] = lastactivesdpvar;
      }
      else if (sdpi->inputtosdpmapper[var] != -1)
      {
         /* this variable was active and stays active */
         lastactivesdpvar = sdpi->inputtosdpmapper[var];
      }

      sdpi->varstochgfixing[var] = FALSE;
   }

   /* shrink the size of the sdptoinputmapper */
   if (nleftshiftvars > 0)
   {
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdptoinputmapper), arraylength, arraylength - nleftshiftvars);
      arraylength -= nleftshiftvars;
   }

   sdpi->fixingstodo = FALSE;
   sdpi->nvars -= nfixedvars;
   assert ( sdpi->nvars == arraylength );
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver, getting version doesn't seem to be supported by DSDP */
const char* SCIPsdpiGetSolverName(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   )
{
   return SCIPsdpiSolverGetSolverName;
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   )
{
   return SCIPsdpiSolverGetSolverDesc;
}

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   )
{
   return SCIPsdpiSolverGetSolverPointer;
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

   SCIPdebugMessage("Calling SCIPsdpiCreate (%d)\n",nextsdpid);

   BMSallocBlockMemory(blkmem, sdpi);

   (*sdpi)->messagehdlr = messagehdlr;
   (*sdpi)->blkmem = blkmem;
   (*sdpi)->sdpid = 1;
   (*sdpi)->nvars = 0;
   (*sdpi)->nvarsorig = 0;
   (*sdpi)->nsdpblocks = 0;
   (*sdpi)->nsdpblocksorig = 0;
   (*sdpi)->sdpconstnnonz = 0;
   (*sdpi)->sdpconstnnonzorig = 0;
   (*sdpi)->sdpnnonz = 0;
   (*sdpi)->sdpnnonzorig = 0;
   (*sdpi)->nlpcons = 0;
   (*sdpi)->nlpconsorig = 0;
   (*sdpi)->lpnnonz = 0;
   (*sdpi)->lpnnonzorig = 0;

   (*sdpi)->obj = NULL;
   (*sdpi)->objorig = NULL;
   (*sdpi)->lb = NULL;
   (*sdpi)->lborig = NULL;
   (*sdpi)->ub = NULL;
   (*sdpi)->uborig = NULL;
   (*sdpi)->sdpblocksizes = NULL;
   (*sdpi)->sdpblocksizesorig = NULL;
   (*sdpi)->sdpconstnblocknonz = NULL;
   (*sdpi)->sdpconstnblocknonzorig = NULL;
   (*sdpi)->sdpconstrow = NULL;
   (*sdpi)->sdpconstroworig = NULL;
   (*sdpi)->sdpconstcol = NULL;
   (*sdpi)->sdpconstcolorig = NULL;
   (*sdpi)->sdpconstval = NULL;
   (*sdpi)->sdpconstvalorig = NULL;
   (*sdpi)->sdpnblockvarnonz = NULL;
   (*sdpi)->sdpnblockvarnonzorig = NULL;
   (*sdpi)->sdprow = NULL;
   (*sdpi)->sdproworig = NULL;
   (*sdpi)->sdpcol = NULL;
   (*sdpi)->sdpcolorig = NULL;
   (*sdpi)->sdpval = NULL;
   (*sdpi)->sdpvalorig = NULL;
   (*sdpi)->lprhs = NULL;
   (*sdpi)->lprhsorig = NULL;
   (*sdpi)->lprowind = NULL;
   (*sdpi)->lprowindorig = NULL;
   (*sdpi)->lpcolind = NULL;
   (*sdpi)->lpcolindorig = NULL;
   (*sdpi)->lpval = NULL;
   (*sdpi)->lpvalorig = NULL;

   (*sdpi)->sdptoinputmapper = NULL;
   (*sdpi)->inputtosdpmapper = NULL;
   (*sdpi)->fixedvals = NULL;
   (*sdpi)->varstochgfixing = NULL;
   (*sdpi)->fixingstodo = FALSE;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP interface structure */
   )
{
   int i;
   int j;
   int pos;

   SCIPdebugMessage("Calling SCIPsdpiFree (%d)\n",(*sdpi)->sdpid);
   assert ( sdpi != NULL );
   assert ( *sdpi != NULL );

   if (((*sdpi)->dsdp) != NULL)
   {
   DSDP_CALL(DSDPDestroy((*sdpi)->dsdp));
   }

   /* free the individual nonzeros */
   for (i = 0; i < (*sdpi)->nsdpblocks; i++)
   {
      BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstcol[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstrow[i]), (*sdpi)->sdpconstnblocknonz[i]);
      BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstval[i]), (*sdpi)->sdpconstnblocknonz[i]);

      for (j = 0; j < (*sdpi)->nvars; j++)
      {
         pos = i * (*sdpi)->nvars + j;
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpcol[pos]), (*sdpi)->sdpnblockvarnonz[pos]);
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdprow[pos]), (*sdpi)->sdpnblockvarnonz[pos]);
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpval[pos]), (*sdpi)->sdpnblockvarnonz[pos]);
      }
   }

   /* and the same again for the original arrays */
   for (i = 0; i < (*sdpi)->nsdpblocksorig; i++)
      {
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstcolorig[i]), (*sdpi)->sdpconstnblocknonzorig[i]);
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstroworig[i]), (*sdpi)->sdpconstnblocknonzorig[i]);
         BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstvalorig[i]), (*sdpi)->sdpconstnblocknonzorig[i]);

         for (j = 0; j < (*sdpi)->nvarsorig; j++)
         {
            pos = i * (*sdpi)->nvarsorig + j;
            BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpcolorig[pos]), (*sdpi)->sdpnblockvarnonzorig[pos]);
            BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdproworig[pos]), (*sdpi)->sdpnblockvarnonzorig[pos]);
            BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpvalorig[pos]), (*sdpi)->sdpnblockvarnonzorig[pos]);
         }
      }

   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->obj), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->objorig), (*sdpi)->nvarsorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lb), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lborig), (*sdpi)->nvarsorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->ub), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->uborig), (*sdpi)->nvarsorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpblocksizes), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpblocksizesorig), (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstnblocknonz), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstnblocknonzorig), (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstrow), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstroworig), (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstcol), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstcolorig), (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstval), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstvalorig), (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonz), (*sdpi)->nvars * (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpnblockvarnonzorig), (*sdpi)->nvarsorig * (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdprow), (*sdpi)->nvars * (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdproworig), (*sdpi)->nvarsorig * (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpcol), (*sdpi)->nvars * (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpcolorig), (*sdpi)->nvarsorig * (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpval),(*sdpi)->nvars * (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpvalorig), (*sdpi)->nvarsorig * (*sdpi)->nsdpblocksorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprhs), (*sdpi)->nlpcons);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprhsorig), (*sdpi)->nlpconsorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprowind), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprowindorig), (*sdpi)->lpnnonzorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpcolind), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpcolindorig), (*sdpi)->lpnnonzorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpval), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpvalorig), (*sdpi)->lpnnonzorig);

   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdptoinputmapper), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->inputtosdpmapper), (*sdpi)->nvarsorig);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->varstochgfixing), (*sdpi)->nvarsorig);

   BMSfreeBlockMemory((*sdpi)->blkmem, sdpi);

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
 */
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                  *  number of entries  of sdpconst row/col/val [i] */
   int const* const*     sdpconstrow,        /**< pointer to row-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   const int**           sdpconstcol,        /**< pointer to column-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   const SCIP_Real**     sdpconstval,        /**< pointer to values of entries of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpnblockvarnonz,   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                                  *  the number of entries of sdp row/col/val [i * nvars + j] */
   const int**           sdprow,             /**< pointer to row indices for each block and variable (may be NULL if sdpnnonz = 0) */
   const int**           sdpcol,             /**< pointer to column indices for each block and variable (may be NULL if sdpnnonz = 0) */
   const SCIP_Real*      sdpval,             /**< pointer to values of nonzeros for each block and variable (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   )
{
   int i;
   int col;
   int row;
   int nfixedvars;
   int nfixednonz;
   int var;
   int block;
   int endindex;
   int naddednonz;
   int nblockaddednonz;
   int ind;
   int* coltoadd;
   int* rowtoadd;
   SCIP_Real* valtoadd;
   int nblockvarnonz;
   int nfixedlpnonz;
   int nblockfixednonz;

   SCIPdebugMessage("Calling SCIPsdpiLoadSDP (%d)\n",sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( obj != NULL );
   assert ( lb != NULL );
   assert ( ub != NULL );

#ifdef SCIP_DEBUG
   if (sdpconstnnonz > 0 || sdpnnonz > 0 || nsdpblocks > 0)
   {
      assert ( nsdpblocks > 0 );
      assert ( sdpconstbegblock != NULL );
      assert ( sdpbegvarblock != NULL );

      if (sdpconstnnonz > 0)
      {
         assert ( sdpconstrowind != NULL );
         assert ( sdpconstcolind != NULL );
         assert ( sdpconstval != NULL );
      }

      if (sdpnnonz > 0)
      {
         assert ( sdprowind != NULL );
         assert ( sdpcolind != NULL );
         assert ( sdpval != NULL );
      }
   }
#endif

   assert ( nlpcons == 0 || lprhs != NULL );
   assert ( lpnnonz == 0 || lprowind != NULL );
   assert ( lpnnonz == 0 || lpcolind != NULL );
   assert ( lpnnonz == 0 || lpval != NULL );

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->objorig), sdpi->nvarsorig, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lborig), sdpi->nvarsorig, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->uborig), sdpi->nvarsorig, nvars);

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdptoinputmapper), sdpi->nvars, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->inputtosdpmapper), sdpi->nvarsorig, nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->varstochgfixing), sdpi->nvarsorig, nvars);

   /* initialize varstochgfixing with FALSE */
   for (i = 0; i < nvars; i++)
   {
      sdpi->varstochgfixing[i] = FALSE;
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizesorig), sdpi->nsdpblocksorig, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstnblocknonz), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstnblocknonzorig), sdpi->nsdpblocksorig, nsdpblocks);

   nfixedvars = 0;

   /* set obj, lb, ub, check if they need to be fixed and set the varmapper arrays */
   for (i = 0; i < nvars; i++)
   {

      if ( !(IsFixed(lb[i], ub[i])))
      {
         (sdpi->obj)[i - nfixedvars] = obj[i];
         (sdpi->objorig)[i] = obj[i];
         (sdpi->lb)[i - nfixedvars] = lb[i];
         (sdpi->lborig)[i] = lb[i];
         (sdpi->ub)[i - nfixedvars] = ub[i];
         (sdpi->uborig)[i] = ub[i];

         sdpi->sdptoinputmapper[i - nfixedvars] = i;
         sdpi->inputtosdpmapper[i] = i - nfixedvars;
      }
      else
      {
         /* this variable should be fixed */
         sdpi->varstochgfixing[i] = TRUE;

         (sdpi->objorig)[i] = obj[i];
         (sdpi->lborig)[i] = lb[i];
         (sdpi->uborig)[i] = ub[i];

         nfixedvars++;

         for (block = 0; block < nsdpblocks; block++)
         {
            nfixednonz += sdpnblockvarnonz[block * nsdpblocks + i];
         }

         sdpi->inputtosdpmapper[i] = -1;
      }
   }

   /* shrink the initialized arrays according to the number of fixed variables */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), nvars, nvars - nfixedvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), nvars, nvars - nfixedvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), nvars, nvars - nfixedvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdptoinputmapper), nvars, nvars - nfixedvars);

   for (i = 0; i < nsdpblocks; i++)
   {
      (sdpi->sdpblocksizes)[i] = sdpblocksizes[i];
      (sdpi->sdpblocksizesorig)[i] = sdpblocksizes[i];
      (sdpi->sdpconstnblocknonz)[i] = sdpconstnblocknonz[i];
   }

   /* prepare the sdpnblockvarnonz arrays */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonz), sdpi->nvars * sdpi->nsdpblocks, (nvars - nfixedvars) * nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpnblockvarnonzorig), sdpi->nvarsorig * sdpi->nsdpblocksorig, nvars * nsdpblocks);

   for (block = 0; block < nsdpblocks; block++)
   {
      for (var = 0; var < nvars; var++)
      {
         sdpi->sdpnblockvarnonzorig[block * nvars + var] = sdpnblockvarnonz[block * nvars + var];
         if (! sdpi->varstochgfixing[var])
            sdpi->sdpnblockvarnonz[block * sdpi->nvars + sdpi->inputtosdpmapper[var]] = sdpnblockvarnonz[block * nvars + var];
      }
   }

   /* sdp nonzeroes */

   /* free all old nonzero arrays */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpi->sdpconstnblocknonz[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpi->sdpconstnblocknonz[block]);
      for (var = 0; var < sdpi->nvars; var++)
      {
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow[block][var]), sdpi->sdpnblockvarnonz[block * sdpi->nvars + var]);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol[block][var]), sdpi->sdpnblockvarnonz[block * sdpi->nvars + var]);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval[block][var]), sdpi->sdpnblockvarnonz[block * sdpi->nvars + var]);
      }
   }
   for (block = 0; block < sdpi->nsdpblocksorig; block++)
   {
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstroworig[block]), sdpi->sdpconstnblocknonzorig[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolorig[block]), sdpi->sdpconstnblocknonzorig[block]);
      BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstvalorig[block]), sdpi->sdpconstnblocknonzorig[block]);
      for (var = 0; var < sdpi->nvarsorig; var++)
      {
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdproworig[block][var]), sdpi->sdpnblockvarnonzorig[block * sdpi->nvarsorig + var]);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolorig[block][var]), sdpi->sdpnblockvarnonzorig[block * sdpi->nvarsorig + var]);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpvalorig[block][var]), sdpi->sdpnblockvarnonzorig[block * sdpi->nvarsorig + var]);
      }
   }

   /* allocate the memory for the constant pointer arrays */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrow), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstroworig), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcol), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolorig), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->nsdpblocks, nsdpblocks);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstvalorig), sdpi->nsdpblocks, nsdpblocks);

   /* allocate the memory for the non-constant pointer arrays */
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow), sdpi->nsdpblocks * sdpi->nvar, nsdpblocks * (nvars-nfixedvars));
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdproworig), sdpi->nsdpblocksorig * sdpi->nvarsorig, nsdpblocks * nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol), sdpi->nsdpblocks * sdpi->nvar, nsdpblocks * (nvars-nfixedvars));
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolorig), sdpi->nsdpblocksorig * sdpi->nvarsorig, nsdpblocks * nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->nsdpblocks * sdpi->nvar, nsdpblocks * (nvars-nfixedvars));
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpvalorig), sdpi->nsdpblocksorig * sdpi->nvarsorig, nsdpblocks * nvars);

   /* allocate the memory for the nonzero arrays */
   for (block = 0; block < nsdpblocks; block++)
   {
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrow[block]), sdpconstnblocknonz[block]);
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcol[block]), sdpconstnblocknonz[block]);
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval[block]), sdpconstnblocknonz[block]);
      for (var = 0; var < nvars; var++)
      {
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprow[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcol[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
      }
   }
   for (block = 0; block < nsdpblocks; block++)
   {
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstroworig[block]), sdpconstnblocknonz[block]);
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolorig[block]), sdpconstnblocknonz[block]);
      BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstvalorig[block]), sdpconstnblocknonz[block]);
      for (var = 0; var < nvars - nfixedvars; var++)
      {
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdproworig[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolorig[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
         BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpvalorig[block * nvars + var]), sdpnblockvarnonz[block * nvars + var]);
      }
   }


   sdpi->nvars = nvars - nfixedvars;
   sdpi->nvarsorig = nvars;
   sdpi->nsdpblocks = nsdpblocks;
   sdpi->nsdpblocksorig = nsdpblocks;

   /* insert the constant nonzeroes */
   for (block = 0; block < nsdpblocks; block++)
   {
      for (i = 0; i < sdpconstnblocknonz[block]; i++)
      {
         sdpi->sdpconstcol[block][i] = sdpconstcol[block][i];
         sdpi->sdpconstcolorig[block][i] = sdpconstcol[block][i];
         sdpi->sdpconstrow[block][i] = sdpconstrow[block][i];
         sdpi->sdpconstroworig[block][i] = sdpconstrow[block][i];
         sdpi->sdpconstval[block][i] = sdpconstval[block][i];
         sdpi->sdpconstvalorig[block][i] = sdpconstval[block][i];
      }
   }

   /* allocate memory for nonzeroes that need to be fixed */
   BMSallocBlockMemoryArray(sdpi->blkmem, &colstoadd, nfixednonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &rowstoadd, nfixednonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &valstoadd, nfixednonz);

   naddednonz = 0;

   /* insert the sdp-nonzeroes (constant and non-constant) and add the fixed ones to the const-arrays */
   for (block = 0; block < nsdpblocks; block++)
   {
      for (var = 0; var < nvars; var++)
      {
         if (sdpi->varstochgfixing[var])
         {
            assert ( IsFixed(lb[var], ub[var]));

            /* save the nonzeroes that need to be fixed in the toadd-arrays */
            for (i = 0; i < sdpnblockvarnonz[block * nvars + var]; i++)
            {
               row = sdprow[block * nvars + var][i];
               col = sdpcol[block * nvars + var][i];
               ensureLowerTriangular(&row, &col);

               sdpi->sdpcolorig[block * nvars + var][i] = col;
               colstoadd[nblockfixednonz] = col;
               sdpi->sdprowindorig[block * nvars + var] = row;
               rowstoadd[nblockfixednonz] = row;
               sdpi->sdpvalorig[block * nvars + var] = sdpval[i];
               valstoadd[nblockfixednonz] = -sdpval[i] * lborig[var]; /* this is the final contribution to the constant array, there
                                                                       * is no need to remember the variable that caused this entry
                                                                       * any longer, the minus comes from the fact that we have +A_i but -A_0 */
               nblockfixednonz++;
            }
         }
         else
         {
            assert ( !IsFixed(lb[var], ub[var]) );

            /* insert the sdp nonzeroes for this variable and block */
            for (i = 0; i < sdpnblockvarnonz[block * nvars + var]; i++)
            {
               row = sdprow[block * nvars + var][i];
               col = sdpcol[block * nvars + var][i];
               ensureLowerTriangular(&row, &col);

               sdpi->sdpcol[block * (nvars-nfixedvar) + var][i] = col;
               sdpi->sdpcolindorig[block * nvars + var][i] = col;
               sdpi->sdprowind[block * (nvars-nfixedvar) + var][i] = row;
               sdpi->sdprowindorig[block * nvars + var][i] = row;
               sdpi->sdpval[block * (nvars-nfixedvar) + var][i] = sdpval[block * nvars + var][i];
               sdpi->sdpvalorig[block * nvars + var][i] = sdpval[block * nvars + var][i];
            }
         }
      }

      /* after they have all been assembled, finally take care of the fixed nonzeroes */
      if (nblockfixednonz > 0 )
      {
         SdpVarfixerMergeArrays(sdpi->blkmem, rowstoadd, colstoadd, valstoadd, nblockfixednonz, FALSE, 1.0, sdpi->sdpconstrow[block], sdpi->sdpconstcol[block],
                                sdpi->sdpconstval[block], &sdpi->sdpconstnblocknonz[block]);
      }
   }

   /* free the memory for nonzeroes that need to be fixed */
   BMSfreeBlockMemoryArray(sdpi->blkmem, &colstoadd, nfixednonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &rowstoadd, nfixednonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &valstoadd, nfixednonz);

   /* compute sdpnnonz and sdpconstnnonz */
   sdpi->sdpnnonzorig = sdpnnonz;
   sdpi->sdpconstnnonzorig = sdpconstnnonz;

   sdpi->sdpnnonz = 0;
   sdpi->sdpconstnnonz = 0;

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      sdpi->sdpconstnnonz += sdpi->sdpconstnblocknonz[block];
      for (var = 0; var < sdpi->nvars; var++)
         sdpi->sdpnblockvarnonz[block * sdpi->nvars + var];
   }


   /* LP part */

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, nlpcons);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhsorig), sdpi->nlpconsorig, nlpcons);
   sdpi->nlpcons = nlpcons;
   sdpi->nlpconsorig = nlpcons;

   for (i = 0; i < nlpcons; i++)
   {
      (sdpi->lprhs)[i] = lprhs[i];
      (sdpi->lprhsorig)[i] = lprhs[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, lpnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowindorig), sdpi->lpnnonzorig, lpnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, lpnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolindorig), sdpi->lpnnonzorig, lpnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, lpnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpvalorig), sdpi->lpnnonzorig, lpnnonz);
   sdpi->lpnnonzorig = lpnnonz;

   nfixedlpnonz = 0;

   for (i = 0; i < lpnnonz; i++)
   {
      assert ( lpcolind[i] < nvars );
      if (sdpi->varstochgfixing[lpcolind[i]])
      {
         nfixedlpnonz++;

         (sdpi->lprowindorig)[i] = lprowind[i];
         (sdpi->lpcolindorig)[i] = lpcolind[i];
         (sdpi->lpvalorig)[i] = lpval[i];

         /* substract this value from this rows rhs */
         (sdpi->lprhs)[i] -= lpval[i] * lb[i];
      }
      else
      {
         (sdpi->lprowindorig)[i] = lprowind[i];
         (sdpi->lpcolindorig)[i] = lpcolind[i];
         (sdpi->lpvalorig)[i] = lpval[i];

         /* shift the entries by the number of fixedlpnonzeroes in front of them */
         (sdpi->lprowind)[i - nfixedlpnonz] = lprowind[i];
         (sdpi->lpcolind)[i - nfixedlpnonz] = lpcolind[i];
         (sdpi->lpval)[i - nfixedlpnonz] = lpval[i];
      }
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), lpnnonz, lpnnonz - nfixedlpnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), lpnnonz, lpnnonz - nfixedlpnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), lpnnonz, lpnnonz - nfixedlpnonz);
   sdpi->lpnnonz = lpnnonz - nfixedlpnonz;

   /* as all fixing have been performed, reset varstochgfixing to FALSE */
   for (i = 0; i < nvars; i++)
   {
      sdpi->varstochgfixing[i] = FALSE;
   }

   sdpi->fixingstodo = FALSE;

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** adds another SDP-Block to the problem
 *
 *  @note as \f A_i^j \f is symmetric, only the lower triangular part of it must be specified
 */
SCIP_RETCODE SCIPsdpiAddSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   blocksize,          /**< dimension of the matrices \f A_i^j \f */
   int                   constnnonz,         /**< sum of non-zeroes in the lower triagonal parts of \f A_0^j \f */
   const int*            constrowind,        /**< the row-indices of the non-zero entries of \f A_i^0 \f (may be NULL if constnnonz = 0) */
   const int*            constcolind,        /**< the column-indices of the non-zero entries of \f A_i^0 \f (may be NULL if constnnonz = 0) */
   const SCIP_Real*      constval,           /**< the values of \f A_i^0 \f as specified by constbegrow and constcolind (may be NULL if constnnonz = 0) */
   int                   nnonz,              /**< sum of non-zeroes in the lower triagonal parts of the \f A_i^j \f */
   const int*            begvar,             /**< start index of the matrix \f A_i^j \f for each i */
   const int*            rowind,             /**< the row-indices of the non-zero entries of \f A_i^j \f (may be NULL if nnonz = 0) */
   const int*            colind,             /**< the column-indices of the non-zero entries of \f A_i^j \f (may be NULL if nnonz = 0) */
   const SCIP_Real*      val                 /**< the values of of \f A_i^j \f as specified by begvar, rowind and colind (may be NULL if nnonz = 0) */
   )
{
   int i;
   int row;
   int col;

   SCIPdebugMessage("Adding a block to SDP %d\n",nextsdpid);

   assert ( sdpi != NULL );
   assert ( blocksize >= 0 );
   assert ( (constnnonz == 0 && nnonz == 0) || blocksize > 0 );
   assert ( constnnonz == 0 || constrowind != NULL );
   assert ( constnnonz == 0 || constcolind != NULL );
   assert ( constnnonz == 0 || constval != NULL );
   assert ( nnonz == 0 || rowind != NULL );
   assert ( nnonz == 0 || colind != NULL );
   assert ( nnonz == 0 || val != NULL );

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks + 1);
   (sdpi->sdpblocksizes)[sdpi->nsdpblocks] = blocksize; /* new SDP-Block will be added as the last block of the new SDP */

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks + 1);
   (sdpi->sdpconstbegblock)[sdpi->nsdpblocks] = sdpi->sdpconstnnonz; /* new SDP-Block starts after all the old ones in the arrays */

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);

   for (i = 0; i < constnnonz; i++)
   {
      assert ( constrowind[i] >= 0 );
      assert ( constrowind[i] < blocksize );
      assert ( constcolind[i] >= 0 );
      assert ( constcolind[i] < blocksize );

      row = constrowind[i];
      col = constcolind[i];
      ensureLowerTriangular(&row, &col);

      (sdpi->sdpconstrowind)[sdpi->sdpconstnnonz + i] = row;
      (sdpi->sdpconstcolind)[sdpi->sdpconstnnonz + i] = col;
      (sdpi->sdpconstval)[sdpi->sdpconstnnonz + i] = constval[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nsdpblocks * sdpi->nvars, (sdpi->nsdpblocks + 1) * sdpi->nvars);
   for (i = 0; i < sdpi->nvars; i++)
   {
      /* new SDP-Block starts after all the old ones in the arrays */
      (sdpi->sdpbegvarblock)[(sdpi->nsdpblocks * sdpi->nvars) + i] = sdpi->sdpnnonz + begvar[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);

   for (i = 0; i < nnonz; i++)
   {
      assert ( rowind[i] >= 0 );
      assert ( rowind[i] < blocksize );
      assert ( colind[i] >= 0 );
      assert ( colind[i] < blocksize );

      row = rowind[i];
      col = colind[i];
      ensureLowerTriangular(&row, &col);

      (sdpi->sdprowind)[sdpi->sdpnnonz + i] = row;
      (sdpi->sdpcolind)[sdpi->sdpnnonz + i] = col;
      (sdpi->sdpval)[sdpi->sdpnnonz + i] = val[i];
   }

   sdpi->nsdpblocks++;
   sdpi->sdpconstnnonz = sdpi->sdpconstnnonz + constnnonz;
   sdpi->sdpnnonz = sdpi->sdpnnonz + nnonz;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** deletes a SDP-Block */
SCIP_RETCODE SCIPsdpiDelSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block              /**< index of the SDP-Block (indexing starting at 1) that should be deleted */
   )
{
   int movingblock;
   int var;
   int i;
   int deletedconstnnonz;
   int newsdpconstnnonz;
   int deletednnonz;
   int newsdpnnonz;

   SCIPdebugMessage("Deleting block %d from SDP %d\n",block, nextsdpid);

   assert ( sdpi != NULL );
   assert ( block >= 0 );
   assert ( block < sdpi->nsdpblocks );

   if (block == sdpi->nsdpblocks - 1) /* the block can simply be deleted */
   {
      deletedconstnnonz = sdpi->sdpconstnnonz - sdpi->sdpconstbegblock[block]; /* sdpconstbegblock[block] gives the first index belonging to the deleted block,
                                                                               * all thereafter need to be deleted*/
      newsdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnnonz;
      deletednnonz = sdpi->sdpnnonz - sdpi->sdpbegvarblock[block * sdpi->nvars]; /* sdpbegvarblock[block*nvars] gives the first index of the deleted block */
      newsdpnnonz = sdpi->sdpnnonz - deletednnonz;

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, sdpi->nvars * (sdpi->nsdpblocks - 1));
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, newsdpnnonz);

      sdpi->nsdpblocks--;
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnnonz;
      sdpi->sdpnnonz = sdpi->sdpnnonz - deletednnonz;
   }
   else /* all blocks after the deleted block need to be shifted in the arrays */
   {
      /* compute the new numbers of nonzeroes, these need to be computed before the begvar-arrays are updated, but the old values are still needed for iterating */
      deletedconstnnonz =  sdpi->sdpconstbegblock[block + 1] - sdpi->sdpconstbegblock[block]; /* starting index of the next block minus starting index of the deleted block gives the number of
                                                                                              * nonzeroes of the deleted block, which is then substracted from the old value */
      newsdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnnonz;
      deletednnonz = sdpi->sdpbegvarblock[(block + 1) * sdpi->nvars] - sdpi->sdpbegvarblock[block * sdpi->nvars]; /* same as above, but because of the structure of the sdpbegvarblock-arrays
                                                                                                                  * block*nvars gives the first index of that block */
      newsdpnnonz = sdpi->sdpnnonz - deletednnonz;


      /* all later blocks need to be moved to the left in the arrays to fill the spot of the deleted block */
      for (movingblock = block + 1; movingblock < sdpi->nsdpblocks; movingblock++)
      {
         sdpi->sdpblocksizes[movingblock - 1] = sdpi->sdpblocksizes[movingblock];
         sdpi->sdpconstbegblock[movingblock - 1] = sdpi->sdpconstbegblock[movingblock] - deletedconstnnonz;
         for (var = 0; var < sdpi->nvars; var++)
         {
            sdpi->sdpbegvarblock[movingblock * sdpi->nvars + var - sdpi->nvars] =  sdpi->sdpbegvarblock[movingblock * sdpi->nvars + var] - deletednnonz; /* these are shifted nvars spaces
                                                                                                                                        * to the left, because there are nvars entries
                                                                                                                                        * in sdpbegvarblock belonging to the
                                                                                                                                        * deleted block */
         }
      }
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, sdpi->nvars * (sdpi->nsdpblocks - 1));

      /* shift all nonzeroes to the left by a number of spots equal to the number of nonzeroes in the deleted block */
      for (i = sdpi->sdpconstbegblock[block + 1]; i<sdpi->sdpconstnnonz; i++)
      {
         sdpi->sdpconstrowind[i - deletedconstnnonz] = sdpi->sdpconstrowind[i];
         sdpi->sdpconstcolind[i - deletedconstnnonz] = sdpi->sdpconstcolind[i];
         sdpi->sdpconstval[i - deletedconstnnonz] = sdpi->sdpconstval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, newsdpconstnnonz);

      for (i = sdpi->sdpbegvarblock[(block + 1) * sdpi->nvars]; i < sdpi->sdpnnonz; i++)
      {
         sdpi->sdprowind[i - deletednnonz] = sdpi->sdprowind[i];
         sdpi->sdpcolind[i - deletednnonz] = sdpi->sdpcolind[i];
         sdpi->sdpval[i - deletednnonz] = sdpi->sdpval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, newsdpnnonz);

      sdpi->nsdpblocks--;
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnnonz;
      sdpi->sdpnnonz = sdpi->sdpnnonz - deletednnonz;
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** adds additional variables to the SDP
 *
 *  @note arrays are not checked for duplicates, problems may appear if indices are added more than once
 */
SCIP_RETCODE SCIPsdpiAddVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to be added */
   const SCIP_Real*      obj,                /**< objective function values of new variables */
   const SCIP_Real*      lb,                 /**< lower bounds of new variables */
   const SCIP_Real*      ub,                 /**< upper bounds of new variables */
   int                   sdpnnonz,           /**< number of nonzero elements to be added to the SDP constraint matrices */
   const int*            sdpbegvarblock,     /**< start index of each new variable for each block in sdpval-array, or NULL if sdpnnonz == 0 */
   const int*            sdprowind,          /**< row indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const int*            sdpcolind,          /**< col indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const SCIP_Real*      sdpval,             /**< values of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   int                   lpnnonz,            /**< number of nonzero elements to be added to the LP constraint matrices */
   const int*            lprowind,           /**< row indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const int*            lpcolind,           /**< col indices of LP constraint matrix entries, indexed from 1 to the number of added vars, or NULL if lpnnonz == 0 */
   const SCIP_Real*      lpval               /**< values of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   )
{
   int i;
   int block;
   int toInsert;

   SCIPdebugMessage("Adding %d variables to SDP %d.\n", nvars, nextsdpid);

   assert ( sdpi != NULL );
   assert ( obj != NULL );
   assert ( lb != NULL );
   assert ( ub != NULL );
   assert ( sdpnnonz == 0 || sdpbegvarblock != NULL );
   assert ( sdpnnonz == 0 || sdprowind != NULL );
   assert ( sdpnnonz == 0 || sdpcolind != NULL );
   assert ( sdpnnonz == 0 || sdpval != NULL );
   assert ( lpnnonz == 0 || lprowind != NULL );
   assert ( lpnnonz == 0 || lpcolind != NULL );
   assert ( lpnnonz == 0 || lpval != NULL );

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars, sdpi->nvars + nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars, sdpi->nvars + nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars, sdpi->nvars + nvars);

   for (i = 0; i < nvars; i++)
   {
      sdpi->obj[sdpi->nvars + i] = obj[i];
      sdpi->lb[sdpi->nvars + i] = lb[i];
      sdpi->ub[sdpi->nvars + i] = ub[i];
   }

   if (sdpnnonz > 0)
   {
      int row;
      int col;
      int* begblockold; /* these are the old starting indices of the blocks (only for the first variable in that block), with extra entry equal to spdnnonz */

      BMSallocBlockMemoryArray(sdpi->blkmem, &begblockold, sdpi->nsdpblocks + 1);

      begblockold[sdpi->nsdpblocks] = sdpi->sdpnnonz; /* often begblockold[block+1] is needed, so this extra entry removes some additional if-clauses */

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, (sdpi->nvars + nvars) * sdpi->nsdpblocks);

      for (block = sdpi->nsdpblocks - 1; block > -1; block--)
      {
         begblockold[block] = sdpi->sdpbegvarblock[block * sdpi->nvars];

         for (i = sdpi->nvars - 1; i > -1; i--)
         {
            /* all the old entries have to be shifted to the right to get the needed space for inserting the new variables for all the blocks before it, also the
             * number of nonzeroes that are added for the new variables for the earlier blocks have to be added, for not overwriting needed entries this iteration
             *  has to go from right to left
             */
            sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + i] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] + sdpbegvarblock[block*nvars];
         }

         for (i = 0; i < nvars; i++)
         {
            /* insert the new values, add to them the start-index of the next block in the original problem (as all nonzeroes in earlier and this block of the
             * original problem as well as the new ones in the earlier blocks have to be added before the new nonzeroes in this block), for the last block the
             * number of nonzeroes is used, as sdpbegvarblock[nblocks*nvars] doesn't exist */
               sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars + i] = sdpbegvarblock[block * nvars + i] + begblockold[block + 1];
         }
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);

      /* now insert the nonzero-entries at the right positions of the arrays */
      for (block=sdpi->nsdpblocks - 1; block > -1; block--)
      {
         for (i = begblockold[block + 1]; i > begblockold[block] - 1; i--)
         {
            /* all the entries belonging to block j have to be shifted sdpbegvarblock[block*nvars] entries to the left, as this many entries belonging to
             * earlier blocks and new variables have to be inserted in front of them
             */
            sdpi->sdprowind[i + sdpbegvarblock[block * nvars]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i + sdpbegvarblock[block * nvars]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i + sdpbegvarblock[block * nvars]] = sdpi->sdpval[i];
         }

         /* compute the number of new nonzeroes to be inserted into this block */
         if (block == sdpi->nsdpblocks - 1)
         {
            toInsert = sdpnnonz - sdpbegvarblock[block * nvars];
         }
         else
         {
            toInsert = sdpbegvarblock[(block+1) * nvars] - sdpbegvarblock[block * nvars];
         }
         for (i = 0; i < toInsert; i++)
         {
            /* insert the new values for this block, they are inserted where the first new variable for the specific block starts */
            assert ( sdprowind[sdpbegvarblock[block * nvars]+i] < sdpi->sdpblocksizes[block] );
            assert ( sdpcolind[sdpbegvarblock[block * nvars]+i] < sdpi->sdpblocksizes[block] ); /* the row and column indices shouldn't exceed blocksizes */

            row = sdprowind[sdpbegvarblock[block * nvars]+i];
            col = sdpcolind[sdpbegvarblock[block * nvars]+i];

            ensureLowerTriangular(&row, &col);

            sdpi->sdprowind[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = row;
            sdpi->sdpcolind[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = col;
            sdpi->sdpval[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = sdpval[sdpbegvarblock[block * nvars]+i];
         }
      }

      BMSfreeBlockMemoryArray(sdpi->blkmem, &begblockold, sdpi->nsdpblocks);

      sdpi->sdpnnonz = sdpi->sdpnnonz + sdpnnonz;
   }

   if (lpnnonz > 0)
   {
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);

      for (i = 0; i < lpnnonz; i++)
      {
         assert ( lprowind[i] < sdpi->nlpcons ); /* only insert into existing LP constraints */

         sdpi->lprowind[sdpi->lpnnonz + i] = lprowind[i]; /* just add these at the end, they will be sorted before solving */
         sdpi->lpcolind[sdpi->lpnnonz + i] = lpcolind[i] + sdpi->nvars; /* the columns are added to the right of the old ones, so the column indices have to be shifted
                                                                         * by the number of old variables */
         sdpi->lpval[sdpi->lpnnonz + i] = lpval[i];
      }

      sdpi->lpnnonz = sdpi->lpnnonz + lpnnonz;
   }

   sdpi->nvars = sdpi->nvars + nvars;

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** deletes all variables in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to be deleted */
   int                   lastvar             /**< last variable to be deleted */
   )
{
   int block;
   int i;
   int* deletedsdpnonz; /* entry i will give the number of deleted sdp nonzeroes in block 1 till i+1 (with indexing starting at 1) */
   int deletedvars;
   int firstvarlpind;
   int lastvarlpind;
   int deletedlpnonz;
   int lastindexforshifting;

   SCIPdebugMessage("Deleting vars %d to %d from SDP %d.\n", firstvar, lastvar, nextsdpid);

   assert ( sdpi != NULL );
   assert ( firstvar >= 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar < sdpi->nvars );

   deletedvars = lastvar - firstvar + 1;

   BMSallocBlockMemoryArray(sdpi->blkmem, &deletedsdpnonz, sdpi->nsdpblocks);
   deletedsdpnonz[0] = sdpi->sdpbegvarblock[lastvar + 1] - sdpi->sdpbegvarblock[firstvar]; /* begvarblock[lastvar] gives the first index of the first
                                                                                            * non-deleted block */
   for (block = 1; block < sdpi->nsdpblocks; block++)
   {
      if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars - 1)
         {
         deletedsdpnonz[block] = deletedsdpnonz[block-1] + sdpi->sdpnnonz - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar];
         }
      else
      {
         deletedsdpnonz[block] = deletedsdpnonz[block-1] + sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar + 1]
                                                                             - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar];
      }
   }

   for (i=lastvar + 1; i < sdpi->nvars; i++)
   {
      sdpi->obj[i - deletedvars] = sdpi->obj[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars, sdpi->nvars - deletedvars);

   for (i=lastvar + 1; i < sdpi->nvars; i++)
   {
      sdpi->lb[i - deletedvars] = sdpi->lb[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars, sdpi->nvars - deletedvars);

   for (i=lastvar + 1; i < sdpi->nvars; i++)
   {
      sdpi->ub[i - deletedvars] = sdpi->ub[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars, sdpi->nvars - deletedvars);

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      if (block > 0)
      {
         /* first look at all nonzeroes in this given block before the deleted vars, for the first block there's
          * nothing to do, as no entries before those were deleted */
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars]; i < sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar]; i++)
         {
            sdpi->sdprowind[i - deletedsdpnonz[block - 1]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i - deletedsdpnonz[block - 1]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i - deletedsdpnonz[block - 1]] = sdpi->sdpval[i];
         }
      }
      /* then look at all nonzeroes in this given block after the deleted vars, here they are shifted by deletsdpnnonz[block] instead of block - 1 */
      if (lastvar < sdpi->nvars)  /* if the deleted var is the last one then there aren't any left after it in the same block, this extra if is needed,
                                   * because otherwise the start index of the for loop would be outside the bounds of sdpbegvarblock for the last block*/
      {
         if (block == sdpi->nsdpblocks - 1)
         {
            lastindexforshifting = sdpi->sdpnnonz;
         }
         else
         {
            lastindexforshifting = sdpi->sdpbegvarblock[(block+1) * sdpi->nvars];
         }
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar]; i < lastindexforshifting; i++)
         {
            sdpi->sdprowind[i - deletedsdpnonz[block]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i - deletedsdpnonz[block]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i - deletedsdpnonz[block]] = sdpi->sdpval[i];
         }
      }
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);

   /* sdpbegvarblock should be updated last to still be able to find the variables which should be moved or deleted */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (i = 0; i < sdpi->nvars; i++)
      {
         if (i < firstvar && block > 0) /* for block 0 nothing needs to be done prior to the first deleted var */
         {
            /* the entry will be moved to the corresponding position with the decreased number of variables and the number of deleted nonzeroes in earlier blocks
             * will be substracted */
            sdpi->sdpbegvarblock[block * (sdpi->nvars - deletedvars) + i] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] - deletedsdpnonz[block - 1];
         }
         else if (i > lastvar)
         {
            sdpi->sdpbegvarblock[block * (sdpi->nvars - deletedvars) + i - deletedvars] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] -
                  deletedsdpnonz[block]; /* this time it is moved even further left because in this block the variables were also deleted, and the entry also
                                           * also has to be decreased further because now also the deleted nonzeroes from this block must be substracted */
         }
      }
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, (sdpi->nvars - deletedvars) * sdpi->nsdpblocks);

   sdpi->sdpnnonz = sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1];

   BMSfreeBlockMemoryArray(sdpi->blkmem, &(deletedsdpnonz), sdpi->nsdpblocks);

   /* now the LP arrays are cleared of the deleted vars, for this they will have to be sorted first to have those belonging to the deleted vars together */
   SCIPsortIntIntReal(sdpi->lpcolind, sdpi->lprowind, sdpi->lpval, sdpi->lpnnonz); /* now the arrays should be sorted by nondecreasing column indices */

   /*iterate over the lpcolind array to find the first index belonging to a deleted var */
   firstvarlpind = -1; /* if this stays at -1 the deleted variables weren't part of the LP block */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lpcolind[i] >= firstvar && sdpi->lpcolind[i] <= lastvar)
      {
         firstvarlpind = i;
         lastvarlpind = i;
         i++;
         break;
      }
   }

   if (firstvarlpind > -1) /* if this is still -1 nothing has to be done for the LP part */
   {
   /* now find the last occurence of a deleted variable (as these are sorted all in between also belong to deleted vars and will be removed) */
      while (i < sdpi->lpnnonz && sdpi->lpcolind[i] <= lastvar)
      {
         lastvarlpind++;
         i++;
      }

      deletedlpnonz = lastvarlpind - firstvarlpind + 1;

      /* finally shift all LP-array-entries after the deleted variables */
      for (i = lastvarlpind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcolind[i - deletedlpnonz] = sdpi->lpcolind[i] - deletedvars; /* the column indices have to be decreased by the number of vars deleted
          * before that var */
         sdpi->lprowind[i - deletedlpnonz] = sdpi->lprowind[i];
         sdpi->lpval[i - deletedlpnonz] = sdpi->lpval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);

      sdpi->lpnnonz = sdpi->lpnnonz - deletedlpnonz;
   }
   sdpi->nvars = sdpi->nvars - deletedvars;

   sdpi->solved = FALSE;

   /* at this point there could be checked if any SDP-blocks or LP-rows have become empty (no variables left), but this isn't done,
    * because then the indices of all blocks/rows behind it would change, possibly creating problems if the user wanted to insert
    * variables into them (or even the deleted block/row) afterwards
    */

   return SCIP_OKAY;
}

/** deletes variables from SCIP_SDPI; the new position of a variable must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelVarset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of variables
                                              *   input:  1 if variable should be deleted, 0 if not
                                              *   output: new position of variable, -1 if variable was deleted */
   )
{
   int i;
   int deletedvars;
   int oldnvars;

   SCIPdebugMessage("Calling SCIPsdpiDelColset for sdpi %d.\n", sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( dstat != NULL );

   deletedvars = 0;
   oldnvars = sdpi->nvars;

   for (i=0; i < oldnvars; i++)
   {
      if (dstat[i] == 1)
      {
         SCIPsdpiDelVars(sdpi, i - deletedvars, i - deletedvars); /* delete this variable, as earliers deletions are already applied to the problem,
                                                                   * the index also has to be lowered by deletedvars */
         dstat[i] = -1;
         deletedvars++;
      }
      else
      {
         dstat[i] = i - deletedvars;
      }
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** adds rows to the LP-Block
 *
 *  @note arrays are not checked for duplicates, problems may appear if indices are added more than once
 */
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            rowind,             /**< row indices of constraint matrix entries, going from 0 to nrows - 1, these will be changed to nlpcons + i */
   const int*            colind,             /**< column indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   int i;

   SCIPdebugMessage("Adding %d LP-Constraints to SDP %d.\n", nrows, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( rhs != NULL );
   assert ( rowind != NULL );
   assert ( colind != NULL );
   assert ( val != NULL );

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons + nrows);
   for (i=0; i < nrows; i++)
   {
      sdpi->lprhs[sdpi->nlpcons + i] = rhs[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);

   for (i=0; i < nnonz; i++)
   {
      assert ( rowind[i] < nrows );
      sdpi->lprowind[sdpi->lpnnonz + i] = rowind[i] + sdpi->nlpcons; /* the new rows are added at the end, so the row indices are increased by the old
                                                                      * number of LP-constraints */

      assert ( colind[i] < sdpi->nvars ); /* only existing vars should be added to the LP-constraints */
      sdpi->lpcolind[sdpi->lpnnonz + i] = colind[i];

      sdpi->lpval[sdpi->lpnnonz + i] = val[i];
   }

   sdpi->nlpcons = sdpi->nlpcons + nrows;
   sdpi->lpnnonz = sdpi->lpnnonz + nnonz;

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP-Block */
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

   deletedrows = lastrow - firstrow + 1;
   deletednonz = 0;

   /* first delete the right-hand-sides */
   for (i = lastrow + 1; i < sdpi->nlpcons; i++) /* shift all rhs after the deleted rows */
   {
      sdpi->lprhs[i - deletedrows] = sdpi->lprhs[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows);

   /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
   SCIPsortIntIntReal(sdpi->lprowind, sdpi->lpcolind, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

   firstrowind = -1;
   /*iterate over the lprowind array to find the first index belonging to a row that should be deleted */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprowind[i] >= firstrow && sdpi->lprowind[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
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
      while (i < sdpi->lpnnonz && sdpi->lprowind[i] <= lastrow)
      {
         lastrowind++;
         i++;
      }
      deletednonz = lastrowind - firstrowind + 1;

      /* finally shift all LP-array-entries after the deleted rows */
      for (i = lastrowind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcolind[i - deletednonz] = sdpi->lpcolind[i];
         sdpi->lprowind[i - deletednonz] = sdpi->lprowind[i] - deletedrows; /* all rowindices after the deleted ones have to be lowered to still have ongoing
                                                                             * indices from 0 to nlpcons-1 */
         sdpi->lpval[i - deletednonz] = sdpi->lpval[i];
      }
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   sdpi->nlpcons = sdpi->nlpcons - deletedrows;
   sdpi->lpnnonz = sdpi->lpnnonz - deletednonz;

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** deletes LP rows from SCIP_SDPI; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of LP rows
                                              *   input:  1 if row should be deleted, 0 if not
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
         SCIPsdpiDelLPRows(sdpi, i - deletedrows, i - deletedrows); /* delete this row, it is shifted by - deletedrows, because in this
                                                                             * problem the earlier rows have already been deleted */
         dstat[i] = -1;
         deletedrows++;
      }
      else
         dstat[i] = i - deletedrows;
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Calling SCIPsdpiClear for SDPI (%d) \n", sdpi->sdpid);

   assert ( sdpi != NULL );

   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);

   sdpi->nvars = 0;
   sdpi->nsdpblocks = 0;
   sdpi->sdpconstnnonz = 0;
   sdpi->sdpnnonz = 0;
   sdpi->nlpcons = 0;
   sdpi->lpnnonz = 0;

   sdpi->solved = FALSE;

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
      assert ( ind[i] >= 0 );
      assert ( ind[i] < sdpi->nvars );
      sdpi->lb[ind[i]] = lb[i];
      sdpi->ub[ind[i]] = ub[i];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes right hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("Changing %d right hand sides of SDP %d\n", nrows, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( ind != NULL );
   assert ( rhs != NULL );

   for (i = 0; i < nrows; i++)
   {
      assert ( ind[i] >= 0 );
      assert ( ind[i] < sdpi->nlpcons );
      sdpi->lprhs[ind[i]] = rhs[i];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes a single coefficient in LP constraint matrix */
SCIP_RETCODE SCIPsdpiChgLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of LP-coefficient to change */
   int                   col,                /**< column number of LP-coefficient to change */
   SCIP_Real             newval              /**< new value of LP-coefficient */
   )
{
   int i;
   SCIP_Bool found;

   SCIPdebugMessage("Changed the LP Coefficient in row %d and colum %d of SDP %d.\n", row, col, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( row >= 0 );
   assert ( row < sdpi->nlpcons );
   assert ( col >= 0 );
   assert ( col < sdpi->nvars );

   /* check if that entry already exists */
   found = FALSE;
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprowind[i] == row && sdpi->lpcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->lpval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("An LP Coefficient in row %d and colum %d of SDP %d didn't exist so far or was zero, it is now added.\n", row, col, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + 1);

      sdpi->lprowind[sdpi->lpnnonz] = row;
      sdpi->lpcolind[sdpi->lpnnonz] = col;
      sdpi->lpval[sdpi->lpnnonz] = newval;
      sdpi->lpnnonz = sdpi->lpnnonz + 1;
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes objective values of variables in the SDP */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to change objective value for */
   int*                  ind,                /**< variable indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for variables */
   )
{
   int i;

   SCIPdebugMessage("Changing %d objective values of SDP %d.\n", ncols, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( ind != NULL );
   assert ( obj != NULL );

   for (i=0; i < nvars; i++)
   {
      assert ( ind[i] >= 0 );
      assert ( ind[i] < sdpi->nvars );
      sdpi->obj[ind[i]] = obj[ind[i]];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes a single coefficient in constant matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPConstCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   rowind,             /**< row index */
   int                   colind,             /**< column index*/
   const SCIP_Real      newval               /**< new value of given entry of constant matrix in given SDP-Block */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;
   SCIP_Bool found;

   SCIPdebugMessage("Changing a SDP coefficient in block %d of SDP %d.\n", block, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( block >= 0 );
   assert ( block < sdpi->nsdpblocks );
   assert ( rowind >= 0 );
   assert ( rowind < sdpi->sdpblocksizes[block] );
   assert ( colind >= 0 );
   assert ( colind < sdpi->sdpblocksizes[block] );

   row = rowind; /* pointers are needed to give these to the function checking if this is a position in the lower triangular part */
   col = colind;
   ensureLowerTriangular(&row, &col); /* make sure that this is a lower triangular position, otherwise it could happen that an upper
                                  * triangular position is added if the same entry is already filled in the lower triangular part */

   /* check if that entry already exists */
   found = FALSE;
   if (block == sdpi->nsdpblocks)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpconstbegblock[block + 1];
   }
   for (i = sdpi->sdpconstbegblock[block]; i < lastiterationindex; i++)
   {
      if (sdpi->sdpconstrowind[i] == row && sdpi->sdpconstcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->sdpconstval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("A constant SDP Coefficient in row %d and colum %d of block %d of SDP %d didn't exist so far or was zero, it is now added.\n",
            rowind, colind, block, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);

      /* move all sdpconstnonzeroes of later blocks one space in the arrays to be able to insert this one at the right position */
      for (i = sdpi->sdpconstnnonz - 1; i >= sdpi->sdpconstbegblock[block + 1]; i--)
      {
         sdpi->sdpconstrowind[i + 1] = sdpi->sdpconstrowind[i];
         sdpi->sdpconstcolind[i + 1] = sdpi->sdpconstcolind[i];
         sdpi->sdpconstval[i + 1] = sdpi->sdpconstval[i];
      }

      /* insert the new entries at the right position (namely what was originally the first position of the next block) */
      sdpi->sdpconstrowind[sdpi->sdpconstbegblock[block + 1]] = row;
      sdpi->sdpconstcolind[sdpi->sdpconstbegblock[block + 1]] = col;
      sdpi->sdpconstval[sdpi->sdpconstbegblock[block + 1]] = newval;

      /* update other information */
      for (i = block + 1; i < sdpi->nsdpblocks; i++)
         {
         sdpi->sdpconstbegblock[i]++; /* all later blocks start one spot later in the arrays */
         }
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz + 1;

   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** changes a single coefficient in a constraint matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   var,                /**< variable index */
   int                   rowind,             /**< row index between */
   int                   colind,             /**< column index between */
   const SCIP_Real      newval              /**< new value of given entry of the give constraint matrix in specified SDP-Block */
   )
{
   int i;
   int lastiterationindex;
   SCIP_Bool found;
   int row;
   int col;

   SCIPdebugMessage("Changing a Coefficient of Matrix A_%d^%d in SDP %d\n", var, block, sdpi->sdpid);

   assert ( sdpi != NULL );
   assert ( block >= 0 );
   assert ( block < sdpi->nsdpblocks );
   assert ( var >= 0 );
   assert ( var < sdpi->nvars );
   assert ( rowind >= 0 );
   assert ( rowind < sdpi->sdpblocksizes[block - 1] ); /* index shift because arrays start at 0 */
   assert ( colind >= 0 );
   assert ( colind < sdpi->sdpblocksizes[block - 1] );

   row = rowind;
   col = colind;
   ensureLowerTriangular(&row, &col); /* make sure that this is a lower triangular position, otherwise it could happen that an upper
                                    * triangular position is added if the same entry is already filled in the lower triangular part */

   /* check if that entry already exists */
   found = FALSE;
   if (block == sdpi->nsdpblocks - 1 && var == sdpi->nvars - 1)
   {
      lastiterationindex = sdpi->sdpnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpbegvarblock[block * sdpi->nvars + var + 1];
   }
   for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + var]; i < lastiterationindex; i++)
   {
      if (sdpi->sdprowind[i] == row && sdpi->sdpcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->sdpval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("A SDP Coefficient in row %d and colum %d of of Matrix A_%d^%d of SDP %d didn't exist so far or was zero, it is now added.\n",
            rowind, colind, var, block, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);

      /* move all sdpnonzeroes of later blocks and vars one space in the arrays to be able to insert this one at the right position */
      for (i = sdpi->sdpnnonz - 1; i >= sdpi->sdpbegvarblock[block * sdpi->nvars + var + 1]; i--)
      {
         sdpi->sdprowind[i + 1] = sdpi->sdprowind[i];
         sdpi->sdpcolind[i + 1] = sdpi->sdpcolind[i];
         sdpi->sdpval[i + 1] = sdpi->sdpval[i];
      }

      /* insert the new entries at the right position (namely what was originally the first position of the next block) */
      sdpi->sdprowind[sdpi->sdpbegvarblock[block * sdpi->nvars + var + 1]] = row;
      sdpi->sdpcolind[sdpi->sdpbegvarblock[block * sdpi->nvars + var + 1]] = col;
      sdpi->sdpval[sdpi->sdpbegvarblock[block * sdpi->nvars + var + 1]] = newval;

      /* update other information */
      for (i = block * sdpi->nvars + var + 1; i < sdpi->nsdpblocks * sdpi->nvars; i++)
         {
         sdpi->sdpbegvarblock[i]++; /* all later blocks start one spot later in the arrays */
         }
      sdpi->sdpnnonz = sdpi->sdpnnonz + 1;

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
   assert ( sdpi != NULL );
   *nlprows = sdpi->nlpcons;
   return SCIP_OKAY;
}

/** gets the number of SDP-Blocks in the SDP */
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of rows */
   )
{
   assert ( sdpi != NULL );
   *nsdpblocks = sdpi->nsdpblocks;
   return SCIP_OKAY;
}

/** gets the number of variables in the SDP */
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   assert ( sdpi != NULL );
   *nvars = sdpi->nvars;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP constraint matrices */
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP constraint matrcies */
   )
{
   assert ( sdpi != NULL );
   *nnonz = sdpi->sdpnnonz;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   )
{
   assert ( sdpi != NULL );
   *nnonz = sdpi->sdpconstnnonz;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP Matrix */
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   )
{
   assert ( sdpi != NULL );
   *nnonz = sdpi->lpnnonz;
   return SCIP_OKAY;
}

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either sdparraylength, sdpnnonz, sdpbegblock, sdprowind, sdpcolind and sdpval have to be 0 or NULL,
 *  or all of them have to be bigger than zero or non-NULL, the same is true for the lp-part.
 *  returns all information if the arrays were long enough or overwrites sdparraylength or lparraylength with the lentgh needed to store the information.
 */
SCIP_RETCODE SCIPsdpiGetVarInfos(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to extract information for */
   int                   lastvar,            /**< last variable to extract information for */
   SCIP_Real*            obj,                /**< buffer to store objective in, or NULL */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  sdparraylength,     /**< pointer to length of sdparrays, if this is less than sdpnnonz the needed length will be returned */
   int*                  sdpnnonz,           /**< pointer to store the number of nonzero elements for the given variables in the sdp-constraints returned, or NULL */
   int*                  sdpbegvarblock,     /**< buffer to store start index of each block in sdpval-array (length (lastvar-firstvar + 1) * nblocks), or NULL */
   int*                  sdprowind,          /**< buffer to store row indices of sdp constraint matrix entries, or NULL */
   int*                  sdpcolind,          /**< buffer to store column indices of sdp constraint matrix entries, or NULL */
   SCIP_Real*            sdpval,             /**< buffer to store values of sdp constraint matrix entries, or NULL */
   int*                  lparraylength,      /**< pointer to length of lparrays, if this is less than lpnnonz the needed length will be returned */
   int*                  lpnnonz,            /**< pointer to store the number of nonzero elements of the lp-constraints returned, or NULL */
   int*                  lprowind,           /**< buffer to store row indices of lp constraint matrix entries, or NULL */
   int*                  lpcolind,           /**< buffer to store column indices of lp constraint matrix entries (these will refer to the original variable indices), or NULL */
   SCIP_Real*            lpval               /**< buffer to store values of lp constraint matrix entries, or NULL */
   )
{
   int i;
   int numvars;
   int block;
   int lastiterationindex;
   int ind;
   int firstvarlpind;
   int lastvarlpind;

   assert ( sdpi != NULL );
   assert ( firstvar >= 0 );
   assert ( firstvar < sdpi->nvars );
   assert ( lastvar >= 0 );
   assert ( lastvar < sdpi->nvars );
   assert ( sdparraylength != NULL );
   assert ( lparraylength != NULL );

   numvars = lastvar - firstvar + 1;

   if (obj != NULL)
   {
      for (i = firstvar; i <= lastvar; i++)
      {
         obj[i] = sdpi->obj[i];
      }
   }

   if (lb != NULL)
   {
      assert ( ub != NULL );
      for (i = firstvar; i <= lastvar; i++)
      {
         lb[i] = sdpi->lb[i];
         ub[i] = sdpi->ub[i];
      }
   }

   if (*sdparraylength > 0)
   {
      assert ( sdpnnonz != NULL );
      assert ( sdpbegvarblock != NULL );
      assert ( sdprowind != NULL );
      assert ( sdpcolind != NULL );
      assert ( sdpval != NULL );

      /* count the number of nonzeroes for the given variables */
      for (block = 0; block < sdpi->nsdpblocks; block++)
      {
         if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars - 1)
         {
            *sdpnnonz = *sdpnnonz + sdpi->sdpnnonz - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar];
         }
         else
         {
            *sdpnnonz = *sdpnnonz + sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar + 1]
                                                         - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar];
         }
      }

      /* check if the arrays are long enough to store all information */
      if (*sdpnnonz > *sdparraylength)
      {
         SCIPdebugMessage("In SCIPsdpiGetVarInfos for vars %d to %d, the given sdp-arrays only had length %d, while length %d would have been needed.\n",
               firstvar, lastvar, *sdparraylength, *sdpnnonz);

         *sdparraylength = *sdpnnonz;
         return SCIP_OKAY;
      }

      ind = 0;
      for (block = 0; block < sdpi->nsdpblocks; block ++)
      {
         /* compute sdpbegvarblock */
         sdpbegvarblock[0] = 0;
         for (i = 0; i < numvars; i++)
         {
            if (block < sdpi->nsdpblocks - 1 || i < numvars + 1) /* for the last block-var-combination nothing has to be done, as this would only result
                                                                  * in sdpnnonz */
            {
               sdpbegvarblock[block * numvars + i + 1] = sdpbegvarblock[block * numvars + i] + sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar + i + 1]
                                                                                             - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar + i];
            }
         }

         /* copy the nonzeroes in the corresponding arrays */
         if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars - 1)
         {
            lastiterationindex = sdpi->sdpnnonz;
         }
         else
         {
            lastiterationindex = sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar + 1];
         }
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar]; i < lastiterationindex; i++)
         {
            sdprowind[ind] = sdpi->sdprowind[i];
            sdpcolind[ind] = sdpi->sdpcolind[i];
            sdpval[ind] = sdpi->sdpval[i];
            ind++;
         }
      }
   }

   if (*lparraylength > 0)
   {
      assert ( lpnnonz != 0 );
      assert ( lprowind != 0 );
      assert ( lpcolind != 0 );
      assert ( lpval != 0 );

      /* for copying the lp-nonzeroes first sort the arrays by columns, so that the entries corresponding to the asked for variables are in one block */
      SCIPsortIntIntReal(sdpi->lpcolind, sdpi->lprowind, sdpi->lpval, sdpi->lpnnonz);

      /*iterate over the lpcolind array to find the first index belonging to a one of the variables */
      firstvarlpind = -1; /* if this stays at -1 the variables weren't part of the LP block */
      for (i = 0; i < sdpi->lpnnonz; i++)
      {
         if (sdpi->lpcolind[i] >= firstvar && sdpi->lpcolind[i] <= lastvar)
         {
            firstvarlpind = i;
            lastvarlpind = i;
            i++;
            break;
         }
      }

      if (firstvarlpind > -1) /* if this is still -1 there are no entries */
      {
      /* now find the last occurence of one of the variable (as these are sorted all in between also belong to these vars) */
         while (i < sdpi->lpnnonz && sdpi->lpcolind[i] <= lastvar)
         {
            lastvarlpind++;
            i++;
         }

         *lpnnonz = lastvarlpind - firstvarlpind + 1;

         /* check if the provided arrays are sufficient for copying everything to them */
         if (*lparraylength < *lpnnonz)
         {
            SCIPdebugMessage("In SCIPsdpiGetVarInfos for vars %d to %d, the given lp-arrays only had length %d, while length %d would have been needed.\n",
                  firstvar, lastvar, *lparraylength, *lpnnonz);
            *lparraylength = *lpnnonz;
            return SCIP_OKAY;
         }

         /* now copy the lpnonzeroes */
         ind = 0;
         for (i = firstvarlpind; i <= lastvarlpind; i++)
         {
            lprowind[ind] = sdpi->lprowind[i];
            lpcolind[ind] = sdpi->lpcolind[i];
            lpval[ind] = sdpi->lpval[i];
            ind++;
         }
      }
   }

   return SCIP_OKAY;
}

/** gets LP rows from SDP problem object; the arrays have to be large enough to store all values.
 *  rhs can be null, otherwise it should contain an entry for each row that should be returned,
 *  either nnonz, rowind, colind, and val have to be NULL and *arraylength needs to be 0, or all of them have to be non-NULL.
 *  This will either return all wanted information or overwrite arraylength by the needed length to give this information.
 */
SCIP_RETCODE SCIPsdpiGetLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first LP row to get from SDP */
   int                   lastrow,            /**< last LP row to get from SDP */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  arraylength,        /**< pointer to length of the rowind and colind arrays, if this is less than nnonz it will be overwritten by nnonz */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int nrows;
   int i;
   int firstrowind;
   int lastrowind;
   int ind;

   assert ( sdpi != NULL );
   assert ( firstrow >= 0 );
   assert ( lastrow >= firstrow );
   assert ( lastrow < sdpi->nlpcons );
   assert ( arraylength != NULL );


   nrows = lastrow - firstrow + 1;

   if (rhs != NULL)
   {
      for (i = 0; i < nrows; i++)
      {
         rhs[firstrow + i] = sdpi->lprhs[i];
      }
   }

   if (*arraylength > 0)
   {
      assert ( nnonz != NULL );
      assert ( rowind != NULL );
      assert ( colind != NULL );
      assert ( val != NULL );

      /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
      SCIPsortIntIntReal(sdpi->lprowind, sdpi->lpcolind, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

      firstrowind = -1;
      /*iterate over the lprowind array to find the first index belonging to one of the rows */
      for (i = 0; i < sdpi->lpnnonz; i++)
      {
         if (sdpi->lprowind[i] >= firstrow && sdpi->lprowind[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
         {
            firstrowind = i;
            lastrowind = i;
            i++;
            break;
         }
      }

      if (firstrowind > -1) /* if this is still -1 there are no nonzeroes for the given rows */
      {
         /* now find the last occurence of one of the rows (as these are sorted all in between also belong to these rows) */
         while (i < sdpi->lpnnonz && sdpi->lprowind[i] <= lastrow)
         {
            lastrowind++;
            i++;
         }

         *nnonz = lastrowind - firstrowind + 1;

         /* check if given arrays are long enough */
         if (*nnonz > *arraylength)
         {
            SCIPdebugMessage("In SCIPsdpiGetLPRows for rows %d to %d, the given arrays only had length %d, while length %d would have been needed.\n",
                  firstrow, lastrow, *arraylength, *nnonz);
            *arraylength = *nnonz;
            return SCIP_OKAY;
         }

         /* copy the entries */
         ind = 0;
         for (i = firstrowind; i <= lastrowind; i++)
         {
            rowind[ind] = sdpi->lprowind[i];
            colind[ind] = sdpi->lpcolind[i];
            val[ind] = sdpi->lpval[i];
            ind++;
         }
      }
   }

   return SCIP_OKAY;
}

/** gets a number SDP blocks; the arrays have to be large enough to store all values.
 *  either constnnonz, constbegblock, constrow, constcolind, and constval have to be NULL, or all of them have to be non-NULL, same for the non-constant parts.
 *  This will either return the wanted information or overwrite constarraylength or arraylength by the length needed to store these informations.
 */
SCIP_RETCODE SCIPsdpiGetSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstblock,         /**< first SDP block to get from SDP */
   int                   lastblock,          /**< last SDP block to get from SDP */
   int*                  constarraylength,   /**< pointer to the length of the const arrays, if this is less than constnnonz this will be overwritten by it */
   int*                  constnnonz,         /**< pointer to store the number of nonzero elements returned for the constant part, or NULL */
   int*                  constbegblock,      /**< buffer to store the start indices of the different blocks in the constval-array, or NULL */
   int*                  constrowind,        /**< buffer to store row indices of entries of constant matrices, or NULL */
   int*                  constcolind,        /**< buffer to store column indices of entries of constant matrices, or NULL */
   SCIP_Real*            constval,           /**< buffer to store values of entries of constant matrices, or NULL */
   int*                  arraylength,        /**< pointer to the length of the sdp-arrays, if this is less than sdpnnonz this will be overwritten by it */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begvarblock,        /**< buffer to store the start indices of the different block/var-combinations in the val-array, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;

   assert ( sdpi != NULL );
   assert ( 0 <= firstblock );
   assert ( firstblock <= lastblock );
   assert ( lastblock < sdpi->nsdpblocks);

   if (*constarraylength > 0)
   {

      assert ( constnnonz != NULL );
      assert ( constbegblock != NULL );
      assert ( constrowind != NULL );
      assert ( constcolind != NULL );
      assert ( constval != NULL );

      if (lastblock == sdpi->nsdpblocks - 1)
      {
         *constnnonz = sdpi->sdpconstnnonz - sdpi->sdpconstbegblock[firstblock];
      }
      else
      {
         *constnnonz = sdpi->sdpconstbegblock[lastblock + 1] - sdpi->sdpconstbegblock[firstblock];
      }

      /* check if given arrays are sufficiently long */
      if (*constnnonz > *constarraylength)
      {
         SCIPdebugMessage("In SCIPsdpiGetSDPBlocks for blocks %d to %d, the given constant-arrays only had length %d, while length %d would have been needed.\n",
               firstblock, lastblock, *constarraylength, *constnnonz);
         *constarraylength = *constnnonz;
         return SCIP_OKAY;
      }

      /* compute constbegblock */
      for (i = 0; i <= lastblock - firstblock; i++)
      {
         constbegblock[i] = sdpi->sdpconstbegblock[firstblock + i] - sdpi->sdpconstbegblock[firstblock]; /* starting index of each block is the starting index in the
                                                                                                          * original problem minus that of the first block taken */
      }

      /* copy nonzeroes */
      for (i = 0; i < *constnnonz; i++)
      {
         constrowind[i] = sdpi->sdpconstrowind[sdpi->sdpconstbegblock[firstblock] + i];
         constcolind[i] = sdpi->sdpconstcolind[sdpi->sdpconstbegblock[firstblock] + i];
         constval[i] = sdpi->sdpconstval[sdpi->sdpconstbegblock[firstblock] + i];
      }
   }

   if (*arraylength > 0)
   {
      assert ( nnonz != NULL );
      assert ( begvarblock != NULL );
      assert ( rowind != NULL );
      assert ( colind != NULL );
      assert ( val != NULL );

      if (lastblock == sdpi->nsdpblocks - 1)
      {
         *nnonz = sdpi->sdpnnonz - sdpi->sdpbegvarblock[firstblock * sdpi->nvars];
      }
      else
      {
         *nnonz = sdpi->sdpbegvarblock[(lastblock + 1) * sdpi->nvars] - sdpi->sdpbegvarblock[firstblock * sdpi->nvars];
      }

      /* check if given arrays are long enough */
      if (*nnonz > *arraylength)
      {
         SCIPdebugMessage("In SCIPsdpiGetSDPBlocks for rows %d to %d, the given sdp-arrays only had length %d, while length %d would have been needed.\n",
               firstblock, lastblock, *arraylength, *nnonz);
         *arraylength = *nnonz;
         return SCIP_OKAY;
      }

      /* compute begvarblock */
      for (i = 0; i < (lastblock - firstblock +1) * sdpi->nvars; i++)
      {
         begvarblock[i] = sdpi->sdpbegvarblock[firstblock * sdpi->nvars + i] - sdpi->sdpbegvarblock[firstblock * sdpi->nvars];
      }

      /* copy nonzeroes */
      for (i = 0; i < *nnonz; i++)
      {
         rowind[i] = sdpi->sdprowind[sdpi->sdpbegvarblock[firstblock * sdpi->nvars] + i];
         colind[i] = sdpi->sdpcolind[sdpi->sdpbegvarblock[firstblock * sdpi->nvars] + i];
         val[i] = sdpi->sdpval[sdpi->sdpbegvarblock[firstblock * sdpi->nvars] + i];
      }
   }

   return SCIP_OKAY;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   int i;

   assert ( sdpi != NULL );
   assert ( firstvar >= 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar < sdpi->nvars);
   assert ( vals != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++)
   {
      vals[i] = sdpi->obj[firstvar + i];
   }
   return SCIP_OKAY;
}

/** gets current variable bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   assert ( sdpi != NULL );
   assert ( firstvar >= 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar < sdpi->nvars);
   assert ( lbs != NULL );
   assert ( ubs != NULL );

   for (i = 0; i < lastvar - firstvar + 1; i++)
   {
      if (lbs != NULL)
      {
         lbs[i] = sdpi->lb[firstvar + i];
      }
      if (ubs != NULL)
      {
         ubs[i] = sdpi->ub[firstvar + i];
      }
   }
   return SCIP_OKAY;
}

/** gets current right hand sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< array to store right hand side values */
   )
{
   int i;

   assert ( sdpi != NULL );
   assert ( firstrow >= 0 );
   assert ( firstrow <= lastrow );
   assert ( lastrow < sdpi->nlpcons);

   for (i = 0; i < lastrow - firstrow + 1; i++)
   {
      rhss[firstrow + i] = sdpi->lprhs[i];
   }

   return SCIP_OKAY;
}

/** gets a single coefficient of LP constraint matrix */
SCIP_RETCODE SCIPsdpiGetLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;

   assert ( sdpi != NULL );
   assert ( row >= 0 );
   assert ( row < sdpi->nlpcons );
   assert ( col >= 0 );
   assert ( col < sdpi->nvars);
   assert ( val != NULL );

   /* search for the entry */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lpcolind[i] == col && sdpi->lprowind[i] == row)
      {
         *val = sdpi->lpval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetLPCoef no entry for the given column = %d and row =%d was found.", col, row);
   SCIPABORT();
   return SCIP_ERROR;
}

/** gets a single coefficient of constant SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPConstCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   rowind,             /**< row number of coefficient */
   int                   colind,             /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;

   assert ( sdpi != NULL );
   assert ( block >= 0 );
   assert ( block < sdpi->nsdpblocks );
   assert ( rowind >= 0 );
   assert ( rowind < sdpi->sdpblocksizes[block] );
   assert ( colind >= 0 );
   assert ( colind < sdpi->sdpblocksizes[block] );

   row = rowind;
   col = colind;
   ensureLowerTriangular(&row, &col); /* Because the matrices are symmetric it doesn't matter if a position in the upper or lower triangular
                                    * was given, but only positions in the lower triangular path are saved in the corresponding arrays */

   /* search for the entry */
   if (block == sdpi->nsdpblocks - 1)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpconstbegblock[block + 1];
   }
   for (i = sdpi->sdpconstbegblock[block]; i < lastiterationindex; i++)
   {
      if (sdpi->sdpconstcolind[i] == col && sdpi->sdpconstrowind[i] == row)
      {
         *val = sdpi->sdpconstval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetSDPConstCoef no entry for the given column = %d and row = %d was found.", colind, rowind);
   SCIPABORT();
   return SCIP_ERROR;
}

/** gets a single coefficient of SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   var,                /**< variable index of coefficient, meaning the i in \f A_i^j \f */
   int                   rowind,             /**< row number of coefficient */
   int                   colind,             /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;

   assert ( sdpi != NULL );
   assert ( block >= 0 );
   assert ( block < sdpi->nsdpblocks );
   assert ( var >= 0 );
   assert ( var < sdpi->nvars );
   assert ( rowind >= 0 );
   assert ( rowind < sdpi->sdpblocksizes[block - 1] ); /* indexshift */
   assert ( colind >= 0 );
   assert ( colind < sdpi->sdpblocksizes[block - 1] ); /* indexshift again */

   row = rowind;
   col = colind;
   ensureLowerTriangular(&row, &col); /* Because the matrices are symmetric it doesn't matter if a position in the upper or lower triangular
                                    * was given, but only positions in the lower triangular path are saved in the corresponding arrays */

   /* search for the entry */
   if (block == sdpi->nsdpblocks - 1 && var == sdpi->nvars - 1)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpbegvarblock[block * sdpi->nvars + var];
   }
   for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + var - 1]; i < lastiterationindex; i++)
   {
      if (sdpi->sdpcolind[i] == col && sdpi->sdprowind[i] == row)
      {
         *val = sdpi->sdpval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetSDPConstCoef no entry for the given column=%d and row=%d was found.", colind, rowind);
   SCIPABORT();
   return SCIP_ERROR;
}


/**@} */



/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert ( sdpi != NULL );

   if (sdpi->fixingstodo)
      SCIP_CALL(performFixings(sdpi));

   /* remove empty rows/cols TODO*/

   return SCIPsdpiSolverLoadAndSolve(sdpi->sdpisolver, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub, sdpi->nsdpblocks, sdpi->sdpblocksizes,
                                            sdpi->sdpconstnnonz, sdpi->sdpconstnblocknonz, sdpi->sdpconstrowind, sdpi->sdpconstcolind,
                                            sdpi->sdpconstval, sdpi->sdpnnonz, sdpi->sdpnblockvarnonz, sdpi->sdprowind, sdpi->sdpcolind,
                                            sdpi->sdpval, sdpi->nlpcons, sdpi->lprhs, sdpi->lpnnonz, sdpi->lprowind, sdpi->lpcolind, sdpi->lpval);
}

/** solves the following penalty formulation of the SDP:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u}
 *   \f
 *   alternatively withObj can be to false to set \f b \f to false and only check for feasibility (if the optimal
 *   objective value is bigger than 0 the problem is infeasible, otherwise it's feasible) */
SCIP_RETCODE SCIPsdpiSolvePenalty(
      SCIP_SDPI*            sdpi,               /**< SDP interface structure */
      SCIP_Real             penaltyParam,       /**< the penalty parameter \f \Gamma \f above, needs to be >= 0 */
      SCIP_Bool             withObj             /**< if this is false, the objective is set to 0 */
   )
{
   assert ( sdpi != NULL );
   assert ( penaltyParam >= 0.0 );

   if (sdpi->fixingstodo)
      SCIP_CALL(performFixings(sdpi));

   /* remove empty rows/cols TODO*/

   if (SCIPsdpiSolverKnowsPenalty())
   {
      return SCIPsdpiSolverLoadAndSolveWithPenalty(sdpi->sdpisolver, penaltyParam, withObj, sdpi->nvars, sdpi->obj, sdpi->lb, sdpi->ub,
                                                    sdpi->nsdpblocks, sdpi->sdpblocksizes, sdpi->sdpconstnnonz, sdpi->sdpconstnblocknonz,
                                                    sdpi->sdpconstrowind, sdpi->sdpconstcolind, sdpi->sdpconstval, sdpi->sdpnnonz,
                                                    sdpi->sdpnblockvarnonz, sdpi->sdprowind, sdpi->sdpcolind, sdpi->sdpval, sdpi->nlpcons,
                                                    sdpi->lprhs, sdpi->lpnnonz, sdpi->lprowind, sdpi->lpcolind, sdpi->lpval);
   }
   else
   {
      SCIPerrorMessage("penalty formulation not yet supported for this solver");
      return SCIP_ERROR;
   }


}
/**@} */




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
   assert ( sdpi != NULL );
   return sdpi->solved;
}

/** returns true if the solver could determine whether or not the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver doesn't know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be used */
SCIP_Bool SCIPsdpiFeasibilityKnown(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   assert ( primalfeasible != NULL );
   assert ( dualfeasible != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));

   switch ( *pdfeasible)
   {
      case DSDP_PDFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      case DSDP_UNBOUNDED:
         *primalfeasible = FALSE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      case DSDP_INFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = FALSE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      default: /* should only include DSDP_PDUNKNOWN */
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
         SCIPABORT();
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiExistsPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}


/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiHasPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiExistsDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiHasDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiIsConverged(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_CONVERGED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_UPPERBOUND)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_MAX_IT)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return SCIP_ERROR;
}

/** returns the internal solution status of the solver */
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   switch ( *reason)
   {
      case DSDP_CONVERGED:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 0;
      }
      case DSDP_INFEASIBLE_START:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 1;
      }
      case DSDP_SMALL_STEPS:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 2;
      }
      case DSDP_INDEFINITE_SCHUR_MATRIX:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 2;
      }
      case DSDP_MAX_IT:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 4;
      }
      case DSDP_NUMERICAL_ERROR:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 2;
      }
      case DSDP_UPPERBOUND:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 3;
      }
      case DSDP_USER_TERMINATION:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 6;
      }
      default:
      {
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         return 7;
      }
   }
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return (SCIPsdpiIsConverged(sdpi) && SCIPsdpiIsPrimalFeasible(sdpi) && SCIPsdpiIsDualFeasible(sdpi));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 * that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiIsAcceptable(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   if (SCIPsdpiIsConverged(sdpi))
   {
      return TRUE;
   }
   else
   {
      double* pobj;
      double* dobj;
      double gap;

      printf("Numerical Trouble in DSDP!\n");

      /* if it didn't converge check the optimality gap */
      BMSallocBlockMemory(sdpi->blkmem, &pobj);
      BMSallocBlockMemory(sdpi->blkmem, &dobj);

      DSDP_CALL(DSDPGetPObjective(sdpi->dsdp, pobj));
      DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, dobj));

      gap = abs(*pobj - *dobj);

      if ((gap < epsilon) || ((gap / (0.5 * (abs(*pobj) + abs(*dobj)))) < epsilon)) /* this is the duality gap used in SDPA */
      {
         BMSfreeBlockMemory(sdpi->blkmem, &pobj);
         BMSfreeBlockMemory(sdpi->blkmem, &dobj);
         return TRUE;
      }
      else
      {
         BMSfreeBlockMemory(sdpi->blkmem, &pobj);
         BMSfreeBlockMemory(sdpi->blkmem, &dobj);
         return FALSE;
      }
   }
/* TODO: also check for primal feasibility, as this is also needed for optimality */
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiIgnoreInstability(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert ( sdpi != NULL );
   assert ( objval != NULL );
   CHECK_IF_SOLVED(sdpi);

   DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, objval));
   *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int                   dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL */
   )
{
   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   if ( objval != NULL )
   {
      DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, objval));
      *objval *= -1; /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   }

   if (dualsollength > 0)
   {
      assert(dualsol != NULL);
      DSDP_CALL(DSDPGetY(sdpi->dsdp, dualsol, dualsollength)); /*last entry needs to be the number of variables, will return an error otherwise */
   }
   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpi);

   DSDP_CALL(DSDPGetIts(sdpi->dsdp, iterations));
   return SCIP_OKAY;
}

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPsdpiGetRealSolQuality(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
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
   return 1E+20; /* default infinity from SCIP */
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   )
{
   return ((val <= -SCIPsdpiInfinity(sdpi)) || (val >= SCIPsdpiInfinity(sdpi)));
}

/** returns highest penalty parameter to be used */
SCIP_Real SCIPsdpiMaxPenParam(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   return 1E+10;  /* DSDP will start with penalty param 10^10 if called normally */
}

/** checks if given value is greater or equal to the highest penalty parameter to be used */
SCIP_Bool SCIPsdpiIsGEMaxPenParam(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be compared to maximum penalty parameter */
   )
{
   return ((val <= -SCIPsdpiMaxPenParam(sdpi)) || (val >= SCIPsdpiMaxPenParam(sdpi)));
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
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */
