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

/**@file   SdpVarfixer.h
 * @brief  adds the main functionality to fix/unfix/(multi-)aggregate variables by merging two three-tuple-arrays of row/col/val together
 * @author Tristan Gally
 */

#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__

#include "scip/type_misc.h"
#include "scip/def.h"
#include "SdpVarfixer.h"

/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMRY */
 #define BMS_CALL(x) do \
  { \
  if( NULL == (x) ) \
  { \
  SCIPerrorMessage("No memory in function call\n"); \
  return SCIP_NOMEMORY; \
  } \
  } \
  while( FALSE )

/**
 * sort the given row, col and val arrays first by non-decreasing row-indices, than for those by identical row-indices by non-increasing val-indices
 */
void SdpVarfixerSortRowCol(
   int*                  row,                /* row indices */
   int*                  col,                /* column indices */
   int*                  val,                /* values */
   int                   length              /* length of the given arrays */
   )
{
   int firstentry;
   int nextentry;

   /* first sort by row indices */
   SCIPsortIntIntReal(row, col, val, length);

   /* for those with identical row-indices now sort by non-decreasing col-index, first find all entries with the same row-index */
   nextentry = 0;
   while (nextentry < length)
   {
      firstentry = nextentry; /* the next row starts where the last one ended*/

      while (nextentry < length && targetrow[nextentry] == i) /* as long as the row still matches, increase nextentry */
      {
         nextentry++;
      }

      /* now sort all entries between firstentry and nextentry-1 by their col-indices */
      SCIPsortIntReal(col + firstentry, val + firstentry, nextentry - firstentry);
   }
}

/**
 * Merges two three-tuple-arrays together. The original arrays (which may have multiple entries for the same row and col) will be mulitplied with
 * scalar and then merged into the target arrays (which may not have multiple entries for the same row and col). If there is already an entry for
 * a row/col combination, these two entries will be combined (their values added together), if they cancel each other out the nonzero entry will
 *  be removed. If you think of the matrices described by the two arrays, this is a matrix addition (but only working on the nonzeros for efficiency).
 */
SCIP_RETCODE SdpVarfixerMergeArrays(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  originrow,          /** original row-index-array that is going to be merged */
   int*                  origincol,          /** original column-index-array that is going to be merged */
   SCIP_Real*            originval,          /** original nonzero-values-array that is going to be merged */
   int                   originlength,       /** length of the original arrays */
   SCIP_Bool             originsorted,       /** are the origin arrays already sorted by non-decreasing row and in case of ties col */
   SCIP_Real             scalar,             /** scalar that the original nonzero-values will be multiplied with before merging */
   int*                  targetrow,          /** row-index-array the original array will be merged into */
   int*                  targetcol,          /** column-index-array the original array will be merged into */
   SCIP_real*            targetval,          /** nonzero-values-array the original array will be merged into */
   int*                  targetlength        /** length of the target arrays the original arrays will be merged into, this will be updated to the
                                               * new length after the mergings */
   )
{
   int ind;
   int i;
   int nleftshifts; /* if some nonzeros of the target arrays get deleted, this saves the number of spots the following entries have to be moved
                     * to the left */
   int naddednonz; /* this gives the number of nonzeros that were added to the end of the arrays (this does NOT include those that were added in
                       the middle of the arrays be decreasing the number of leftshifts) */
   int insertionpos;

   assert ( blkmem != NULL );
   assert ( originrow != NULL );
   assert ( origincol != NULL );
   assert ( originval != NULL );
   assert ( origlength >= 0 );
   assert ( targetrow != NULL );
   assert ( targetcol != NULL );
   assert ( targetval != NULL );
   assert ( targetlength != NULL );
   assert ( *targetlength >= 0 );

   /* sort the target and origin arrays first by row and then by col to make searching for entries easier */
   sortRowCol(targetrow, targetcol, targetval, targetlength);

   if (! (originsorted))
      sortRowCol(originrow, origincol, originval, originlength);

   /* allocate memory for the maximum possible size of the target arrays, they will be decreased again afterwards after the number
    * of added nonzeros is known */
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetrow, targetlength, targetlength + originlength));
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetcol, targetlength, targetlength + originlength));
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetval, targetlength, targetlength + originlength));

   ind = 0; /* this will be used to traverse the nonzeros of the target arrays */
   naddednonz = 0;

   /* iterate over all nonzeroes */
   for (i = 0; i < originlength; i++)
   {
      /* search the target arrays for an entry at this position, as both the origin and the target arrays are sorted, we go on until
       * we find an entry that is not < as the current entry in the origin arrays according to this sorting, if this has equal row/col,
       * we have found the entry we have to edit, if it is >, then we know, that there is no identical entry, and we can just add a new
       * entry for this row and col */
      while (ind < targetlength && (targetrow[ind] < originrow[i] || (targetrow[ind] == originrow[i] && targetcol[ind] < origincol[i])))
      {
         /* shift the target nonzeros to the left if needed */
         if (nleftshifts > 0)
         {
            targetrow[ind - nleftshifts] = targetrow[ind];
            targetcol[ind - nleftshifts] = targetcol[ind];
            targetval[ind - nleftshifts] = targetval[ind];
         }
         ind++;
      }

      if (ind < targetlength && (targetrow[ind] == originrow[i] && targetcol[ind] == origincol[i]))
      {
         /* add to the old entry */

         /* shift the entry to the left if needed and change the value */
         if (nleftshifts > 0)
         {
            targetrow[ind - nleftshifts] = targetrow[ind];
            targetcol[ind - nleftshifts] = targetcol[ind];
         }
         targetval[ind - nleftshifts] += scalar * originval[i];

         /* there could be multiple entries to add with identical row and col, so look for further ones in the next entries until there
          * are no more */
         while (i + 1 < originlength && originrow[i + 1] == targetrow[ind - nleftshifts] && origincol[i + 1] == targetcol[ind - nleftshifts])
         {
            targetval[ind - nleftshifts] += scalar * originval[i + 1];
            i++;
         }

         if (REALABS(sdpi->sdpconstval[ind - nleftshiftconstnonz]) < epsilon)
         {
            /* the nonzero became zero */
            nleftshiftconstnonz++;
         }
         ind++; /* as we already added all origin-entries belonging to this row/col and also shifted the entry, we can continue with the next one */
      }
      else  /* create a new entry */
      {
         if (nleftshifts > 0)
         {
            /* we can add the nonzero at one of the empty spots */
            insertionpos = ind - nleftshifts;
            nleftshifts--; /* as one empty spot was filled, all remaining nonzeros should be moved one less position to the left */
         }
         else
         {
            /* add it to the end */
            insertionpos = targetlength + naddednonz;
            naddednonz++;
         }

         /* add the nonzero to the computed position */
         targetrow[insertionpos] = originrow[i];
         targetcol[insertionpos] = origincol[i];
         targetval[insertionpos] = originval[i];

         /* there could be multiple entries to add with identical row and col, so look for further ones in the next entries until there are no more */
         while (i + 1 < originlength && originrow[i + 1] == targetrow[insertionpos] && origincol[i + 1] == targetcol[insertionpos])
         {
            targetval[insertionpos] += originval[i + 1];
            i++;
         }

         /* if there were indeed multiple entries, check if they did cancel each other out, in that case remove the entry */
         if (REALABS(targetval[insertionpos]) < epsilon)
         {
            /* depending on where this actually zero nonzero was added, either add another leftshift to overwrite it or decrease the number of addednonz */
            if (insertionpos <= ind)
               nleftshifts++;
            else
               naddednonz--;
         }
      }
   }
   /* shift the remaining constnonzeros */
   if (nleftshifts > 0)
   {
      while (ind < originlength)
      {
         targetrow[ind - nleftshifts] = targetrow[ind];
         targetcol[ind - nleftshifts] = targetrow[ind];
         targetval[ind - nleftshifts] = targetrow[ind];
         ind++;
      }
   }

   /* shrink the targetarrays to the size that is really needed */
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetrow + originlength, targetlength, targetlength + naddednonz - nleftshifts));
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetcol + originlength, targetlength, targetlength + naddednonz - nleftshifts));
   BMS_CALL(BMSreallocBlockMemoryArray(blkmem, &targetval + originlength, targetlength, targetlength + naddednonz - nleftshifts));

   *targetlength = *targetlength + nadednonz - nleftshifts;

   return SCIP_OKAY;
}
