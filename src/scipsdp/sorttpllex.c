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

/**@file   sorttpllex.c
 * @ingroup OTHER_CFILES
 * @brief  template functions for lexicographic sorting
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/def.h"

#define SORTTPL_SHELLSORTMAX    25 /* maximal size for shell sort */
#define SORTTPL_MINSIZENINTHER 729 /* minimum input size to use ninther (median of nine) for pivot selection */

#ifndef SORTTPL_NAMEEXT
#error You need to define SORTTPL_NAMEEXT.
#endif
#ifndef SORTTPL_KEYTYPE
#error You need to define SORTTPL_KEYTYPE.
#endif

#ifdef SORTTPL_EXPANDNAME
#undef SORTTPL_EXPANDNAME
#endif
#ifdef SORTTPL_NAME
#undef SORTTPL_NAME
#endif

/* enabling and disabling additional lines in the code */
#ifdef SORTTPL_FIELD1TYPE
#define SORTTPL_HASFIELD1(x)    x
#define SORTTPL_HASFIELD1PAR(x) x,
#else
#define SORTTPL_HASFIELD1(x)    /**/
#define SORTTPL_HASFIELD1PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD2TYPE
#define SORTTPL_HASFIELD2(x)    x
#define SORTTPL_HASFIELD2PAR(x) x,
#else
#define SORTTPL_HASFIELD2(x)    /**/
#define SORTTPL_HASFIELD2PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD3TYPE
#define SORTTPL_HASFIELD3(x)    x
#define SORTTPL_HASFIELD3PAR(x) x,
#else
#define SORTTPL_HASFIELD3(x)    /**/
#define SORTTPL_HASFIELD3PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD4TYPE
#define SORTTPL_HASFIELD4(x)    x
#define SORTTPL_HASFIELD4PAR(x) x,
#else
#define SORTTPL_HASFIELD4(x)    /**/
#define SORTTPL_HASFIELD4PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD5TYPE
#define SORTTPL_HASFIELD5(x)    x
#define SORTTPL_HASFIELD5PAR(x) x,
#else
#define SORTTPL_HASFIELD5(x)    /**/
#define SORTTPL_HASFIELD5PAR(x) /**/
#endif
#ifdef SORTTPL_FIELD6TYPE
#define SORTTPL_HASFIELD6(x)    x
#define SORTTPL_HASFIELD6PAR(x) x,
#else
#define SORTTPL_HASFIELD6(x)    /**/
#define SORTTPL_HASFIELD6PAR(x) /**/
#endif

/* the two-step macro definition is needed, such that macro arguments
 * get expanded by prescan of the C preprocessor (see "info cpp",
 * chapter 3.10.6: Argument Prescan)
 */
#define SORTTPL_EXPANDNAME(method, methodname) \
   method ## methodname
#define SORTTPL_NAME(method, methodname) \
  SORTTPL_EXPANDNAME(method, methodname)

/** default comparer for integers */
static
int SORTTPL_CMP(
   SORTTPL_KEYTYPE       x1,                 /**< primary key */
   SORTTPL_KEYTYPE       x2,                 /**< secondary key */
   SORTTPL_KEYTYPE       y1,                 /**< primary key */
   SORTTPL_KEYTYPE       y2                  /**< secondary key */
   )
{
   if ( x1 < y1 )
      return -1;

   if ( x1 > y1 )
      return 1;
   assert( x1 == y1 );

   if ( x2 < y2 )
      return -1;

   if ( x2 > y2 )
      return 1;

   return 0;
}


#ifdef SORTTPL_BACKWARDS
#define SORTTPL_ISLEXBETTER(x1,x2,y1,y2) (SORTTPL_CMP(x1,x2,y1,y2) > 0)
#define SORTTPL_ISLEXWORSE(x1,x2,y1,y2)  (SORTTPL_CMP(x1,x2,y1,y2) < 0)
#else
#define SORTTPL_ISLEXBETTER(x1,x2,y1,y2) (SORTTPL_CMP(x1,x2,y1,y2) < 0)
#define SORTTPL_ISLEXWORSE(x1,x2,y1,y2)  (SORTTPL_CMP(x1,x2,y1,y2) > 0)
#endif

/* swapping two variables */
#define SORTTPL_SWAP(T,x,y) \
   {                \
      T temp = x;   \
      x = y;        \
      y = temp;     \
   }


/** shell-lex-sort an array of data elements; use it only for arrays smaller than 25 entries */
static
void SORTTPL_NAME(sorttpl_lexShellSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   int                   start,              /**< starting index */
   int                   end                 /**< ending index */
   )
{
   static const int incs[3] = {1, 5, 19}; /* sequence of increments */
   int k;

   assert( key1 != NULL );
   assert( key2 != NULL );
   assert( start <= end );

   for (k = 2; k >= 0; --k)
   {
      int h = incs[k];
      int first = h + start;
      int i;

      for (i = first; i <= end; ++i)
      {
         int j;
         SORTTPL_KEYTYPE tempkey1 = key1[i];
         SORTTPL_KEYTYPE tempkey2 = key2[i];

         SORTTPL_HASFIELD1( SORTTPL_FIELD1TYPE tempfield1 = field1[i]; )
         SORTTPL_HASFIELD2( SORTTPL_FIELD2TYPE tempfield2 = field2[i]; )
         SORTTPL_HASFIELD3( SORTTPL_FIELD3TYPE tempfield3 = field3[i]; )
         SORTTPL_HASFIELD4( SORTTPL_FIELD4TYPE tempfield4 = field4[i]; )
         SORTTPL_HASFIELD5( SORTTPL_FIELD5TYPE tempfield5 = field5[i]; )
         SORTTPL_HASFIELD6( SORTTPL_FIELD6TYPE tempfield6 = field6[i]; )

         j = i;
         while (j >= first && SORTTPL_ISLEXBETTER(tempkey1, tempkey2, key1[j-h], key2[j-h]))
         {
            key1[j] = key1[j-h];
            key2[j] = key2[j-h];

            SORTTPL_HASFIELD1( field1[j] = field1[j-h]; )
            SORTTPL_HASFIELD2( field2[j] = field2[j-h]; )
            SORTTPL_HASFIELD3( field3[j] = field3[j-h]; )
            SORTTPL_HASFIELD4( field4[j] = field4[j-h]; )
            SORTTPL_HASFIELD5( field5[j] = field5[j-h]; )
            SORTTPL_HASFIELD6( field6[j] = field6[j-h]; )
            j -= h;
         }

         key1[j] = tempkey1;
         key2[j] = tempkey2;

         SORTTPL_HASFIELD1( field1[j] = tempfield1; )
         SORTTPL_HASFIELD2( field2[j] = tempfield2; )
         SORTTPL_HASFIELD3( field3[j] = tempfield3; )
         SORTTPL_HASFIELD4( field4[j] = tempfield4; )
         SORTTPL_HASFIELD5( field5[j] = tempfield5; )
         SORTTPL_HASFIELD6( field6[j] = tempfield6; )
      }
   }
}

/** returns the index a, b, or c of the lexicographic median element among key[a], key[b], and key[c] */
static
int SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   int                   a,                  /**< first index of the key array to consider */
   int                   b,                  /**< second index of the key array to consider */
   int                   c                   /**< third index of the array to consider */
   )
{
   assert( a >= 0 );
   assert( b >= 0 );
   assert( c >= 0 );
   assert( a != b );
   assert( b != c );
   assert( c != a );

   /* let the elements in the unsorted order be a, b, c at positions start, mid, and end */
   if ( SORTTPL_ISLEXBETTER(key1[a], key2[a], key1[b], key2[b]) ) /* a <= b */
   {
      if ( SORTTPL_ISLEXBETTER(key1[b], key2[b], key1[c], key2[c]) ) /* b <= c */
         /* resulting permutation: a b c */
         return b;
      else /* b > c */
      {
         if ( SORTTPL_ISLEXBETTER(key1[a], key2[a], key1[c], key2[c]) ) /* a <= c */
            /* resulting permutation: a c b */
            return c;
         else
            /* resulting permutation: c a b */
            return a;
      }
   }
   else /* a > b */
   {
      if ( SORTTPL_ISLEXBETTER(key1[b], key2[b], key1[c], key2[c]) )
      {
         if ( SORTTPL_ISLEXBETTER(key1[a], key2[a], key1[c], key2[c]) )
            /* resulting permutation: b a c */
            return a;
         else
            /* resulting permutation: b c a */
            return c;
      }
      else
         /* resulting permutation: c b a */
         return b;
   }
}

/** guess a lexicographic median for the key array [start, ..., end] by using the median of the first, last, and middle element */
static
int SORTTPL_NAME(sorttpl_lexSelectPivotIndex, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   int                   start,              /**< first index of the key array to consider */
   int                   end                 /**< last index of the key array to consider */
   )
{
   int pivotindex;

   /* use the middle index on small arrays */
   if ( end - start + 1 <= SORTTPL_SHELLSORTMAX )
      pivotindex = (start + end) / 2;
   else if ( end - start + 1 < SORTTPL_MINSIZENINTHER )
   {
      /* select the median of the first, last, and middle element as pivot element */
      int mid = (start + end) / 2;
      pivotindex = SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)(key1, key2, start, mid, end);
   }
   else
   {
      /* use the median of medians of nine evenly distributed elements of the key array */
      int gap = (end - start + 1) / 9;
      int median1;
      int median2;
      int median3;

      /* this should always hold */
      assert( start + 8 * gap <= end );

      /* collect 3 medians evenly distributed over the array */
      median1 = SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)(key1, key2, start, start + gap, start + 2 * gap);

      median2 = SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)(key1, key2, start + 3 * gap, start + 4 * gap, start + 5 * gap);
      median3 = SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)(key1, key2, start + 6 * gap, start + 7 * gap, start + 8 * gap);

      /* compute and return the median of the medians */
      pivotindex = SORTTPL_NAME(sorttpl_lexMedianThree, SORTTPL_NAMEEXT)(key1, key2, median1, median2, median3);
   }

   return pivotindex;
}

/** lexicographic quick-sort an array of pointers; pivot is the medial element */
static
void SORTTPL_NAME(sorttpl_lexQSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   int                   start,              /**< starting index */
   int                   end,                /**< ending index */
   SCIP_Bool             type                /**< TRUE, if quick-sort should start with with key[lo] < pivot <= key[hi], key[lo] <= pivot < key[hi] otherwise */
   )
{
   assert( start <= end );

   /* use quick-sort for long lists */
   while ( end - start >= SORTTPL_SHELLSORTMAX )
   {
      SORTTPL_KEYTYPE pivotkey1;
      SORTTPL_KEYTYPE pivotkey2;
      int lo;
      int hi;
      int mid;

      /* select pivot element */
      mid = SORTTPL_NAME(sorttpl_lexSelectPivotIndex, SORTTPL_NAMEEXT)(key1, key2, start, end);
      pivotkey1 = key1[mid];
      pivotkey2 = key2[mid];

      /* partition the array into elements < pivot [start,hi] and elements >= pivot [lo,end] */
      lo = start;
      hi = end;
      for ( ;; )
      {
         if ( type )
         {
            while ( lo < end && SORTTPL_ISLEXBETTER(key1[lo], key2[lo], pivotkey1, pivotkey2) )
               lo++;
            while ( hi > start && ! SORTTPL_ISLEXBETTER(key1[hi], key2[hi], pivotkey1, pivotkey2) )
               hi--;
         }
         else
         {
            while ( lo < end && ! SORTTPL_ISLEXWORSE(key1[lo], key2[lo], pivotkey1, pivotkey2) )
               lo++;
            while ( hi > start && SORTTPL_ISLEXWORSE(key1[hi], key2[hi], pivotkey1, pivotkey2) )
               hi--;
         }

         if ( lo >= hi )
            break;

         SORTTPL_SWAP(SORTTPL_KEYTYPE, key1[lo], key1[hi]);
         SORTTPL_SWAP(SORTTPL_KEYTYPE, key2[lo], key2[hi]);
         SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[hi]); )
         SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[hi]); )
         SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[hi]); )
         SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[lo], field4[hi]); )
         SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[lo], field5[hi]); )
         SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[lo], field6[hi]); )

         lo++;
         hi--;
      }
      assert( (hi == lo-1) || (type && hi == start) || (!type && lo == end) );

      /* skip entries which are equal to the pivot element (three partitions, <, =, > than pivot)*/
      if ( type )
      {
         while ( lo < end && ! SORTTPL_ISLEXBETTER(pivotkey1, pivotkey2, key1[lo], key2[lo]) )
            lo++;

         /* make sure that we have at least one element in the smaller partition */
         if ( lo == start )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert( ! SORTTPL_ISLEXBETTER(key1[mid], key2[mid], pivotkey1, pivotkey2) ); /* the pivot element did not change its position */
            assert( ! SORTTPL_ISLEXBETTER(pivotkey1, pivotkey2, key1[mid], key2[mid]) );
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key1[lo], key1[mid]);
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key2[lo], key2[mid]);
            SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[lo], field1[mid]); )
            SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[lo], field2[mid]); )
            SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[lo], field3[mid]); )
            SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[lo], field4[mid]); )
            SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[lo], field5[mid]); )
            SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[lo], field6[mid]); )
            lo++;
         }
      }
      else
      {
         while ( hi > start && ! SORTTPL_ISLEXWORSE(pivotkey1, pivotkey2, key1[hi], key2[hi]) )
            hi--;

         /* make sure that we have at least one element in the smaller partition */
         if ( hi == end )
         {
            /* everything is greater or equal than the pivot element: move pivot to the left (degenerate case) */
            assert( ! SORTTPL_ISLEXBETTER(key1[mid], key2[mid], pivotkey1, pivotkey2) ); /* the pivot element did not change its position */
            assert( ! SORTTPL_ISLEXBETTER(pivotkey1, pivotkey2, key1[mid], key2[mid]) );
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key1[hi], key1[mid]);
            SORTTPL_SWAP(SORTTPL_KEYTYPE, key2[hi], key2[mid]);
            SORTTPL_HASFIELD1( SORTTPL_SWAP(SORTTPL_FIELD1TYPE, field1[hi], field1[mid]); )
            SORTTPL_HASFIELD2( SORTTPL_SWAP(SORTTPL_FIELD2TYPE, field2[hi], field2[mid]); )
            SORTTPL_HASFIELD3( SORTTPL_SWAP(SORTTPL_FIELD3TYPE, field3[hi], field3[mid]); )
            SORTTPL_HASFIELD4( SORTTPL_SWAP(SORTTPL_FIELD4TYPE, field4[hi], field4[mid]); )
            SORTTPL_HASFIELD5( SORTTPL_SWAP(SORTTPL_FIELD5TYPE, field5[hi], field5[mid]); )
            SORTTPL_HASFIELD6( SORTTPL_SWAP(SORTTPL_FIELD6TYPE, field6[hi], field6[mid]); )
            hi--;
         }
      }

      /* sort the smaller partition by a recursive call, sort the larger part without recursion */
      if ( hi - start <= end - lo )
      {
         /* sort [start,hi] with a recursive call */
         if ( start < hi )
         {
            SORTTPL_NAME(sorttpl_lexQSort, SORTTPL_NAMEEXT)
               (key1, key2,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASFIELD4PAR(field4)
                SORTTPL_HASFIELD5PAR(field5)
                SORTTPL_HASFIELD6PAR(field6)
                  start, hi, !type);
         }

         /* now focus on the larger part [lo,end] */
         start = lo;
      }
      else
      {
         if( lo < end )
         {
            /* sort [lo,end] with a recursive call */
            SORTTPL_NAME(sorttpl_lexQSort, SORTTPL_NAMEEXT)
               (key1, key2,
                SORTTPL_HASFIELD1PAR(field1)
                SORTTPL_HASFIELD2PAR(field2)
                SORTTPL_HASFIELD3PAR(field3)
                SORTTPL_HASFIELD4PAR(field4)
                SORTTPL_HASFIELD5PAR(field5)
                SORTTPL_HASFIELD6PAR(field6)
                  lo, end, !type);
         }

         /* now focus on the larger part [start,hi] */
         end = hi;
      }
      type = !type;
   }

   /* use shell sort on the remaining small list */
   if ( end - start >= 1 )
   {
      SORTTPL_NAME(sorttpl_lexShellSort, SORTTPL_NAMEEXT)
         (key1, key2,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            start, end);
   }
}

#ifndef NDEBUG
/** verifies that an array is indeed sorted */
static
void SORTTPL_NAME(sorttpl_lexCheckSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   int                   len                 /**< length of the array */
   )
{
   int i;

   for ( i = 0; i < len-1; i++ )
   {
      assert( ! SORTTPL_ISLEXBETTER(key1[i+1], key2[i+1], key1[i], key2[i]) );
   }
}
#endif

/** SCIPsort...(): lexicographically sorts array 'key1/key2' and performs the same permutations on the additional 'field' arrays */
static
void SORTTPL_NAME(SCIPlexSort, SORTTPL_NAMEEXT)
(
   SORTTPL_KEYTYPE*      key1,               /**< pointer to data array that defines the primary order */
   SORTTPL_KEYTYPE*      key2,               /**< pointer to data array that defines the secondary order */
   SORTTPL_HASFIELD1PAR(  SORTTPL_FIELD1TYPE*    field1 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD2PAR(  SORTTPL_FIELD2TYPE*    field2 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD3PAR(  SORTTPL_FIELD3TYPE*    field3 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD4PAR(  SORTTPL_FIELD4TYPE*    field4 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD5PAR(  SORTTPL_FIELD5TYPE*    field5 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   SORTTPL_HASFIELD6PAR(  SORTTPL_FIELD6TYPE*    field6 )      /**< additional field that should be sorted in the same way */ /*lint !e665*/
   int                   len                 /**< length of arrays */
   )
{
   /* ignore the trivial cases */
   if ( len <= 1 )
      return;

   /* use shell sort on the remaining small list */
   if ( len <= SORTTPL_SHELLSORTMAX )
   {
      SORTTPL_NAME(sorttpl_lexShellSort, SORTTPL_NAMEEXT)
         (key1, key2,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            0, len-1);
   }
   else
   {
      SORTTPL_NAME(sorttpl_lexQSort, SORTTPL_NAMEEXT)
         (key1, key2,
            SORTTPL_HASFIELD1PAR(field1)
            SORTTPL_HASFIELD2PAR(field2)
            SORTTPL_HASFIELD3PAR(field3)
            SORTTPL_HASFIELD4PAR(field4)
            SORTTPL_HASFIELD5PAR(field5)
            SORTTPL_HASFIELD6PAR(field6)
            0, len-1, TRUE);
   }
#ifndef NDEBUG
   SORTTPL_NAME(sorttpl_lexCheckSort, SORTTPL_NAMEEXT)(key1, key2, len);
#endif
}

/* undefine template parameters and local defines */
#undef SORTTPL_NAMEEXT
#undef SORTTPL_KEYTYPE
#undef SORTTPL_FIELD1TYPE
#undef SORTTPL_FIELD2TYPE
#undef SORTTPL_FIELD3TYPE
#undef SORTTPL_FIELD4TYPE
#undef SORTTPL_FIELD5TYPE
#undef SORTTPL_FIELD6TYPE
#undef SORTTPL_HASFIELD1
#undef SORTTPL_HASFIELD2
#undef SORTTPL_HASFIELD3
#undef SORTTPL_HASFIELD4
#undef SORTTPL_HASFIELD5
#undef SORTTPL_HASFIELD6
#undef SORTTPL_HASFIELD1PAR
#undef SORTTPL_HASFIELD2PAR
#undef SORTTPL_HASFIELD3PAR
#undef SORTTPL_HASFIELD4PAR
#undef SORTTPL_HASFIELD5PAR
#undef SORTTPL_HASFIELD6PAR
#undef SORTTPL_ISLEXBETTER
#undef SORTTPL_ISLEXWORSE
#undef SORTTPL_CMP
#undef SORTTPL_SWAP
#undef SORTTPL_SHELLSORTMAX
#undef SORTTPL_MINSIZENINTHER
