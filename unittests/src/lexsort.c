/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2021 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2021 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lexsort.c
 * @brief  unit testing lexicographic sorting
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "include/scip_test.h"
#include "scip/scip.h"

/* defines for lexicographic sorting macros */
#define SORTTPL_NAMEEXT     IntIntInt
#define SORTTPL_KEYTYPE     int
#define SORTTPL_FIELD1TYPE  int
#include "scipsdp/sorttpllex.c"

/* global variables */
static SCIP_RANDNUMGEN* randgen;
static SCIP* scip;
static unsigned int randomseed = 42;

#define LEN 1000

/* test suite */
static
void setup(void)
{
   SCIPcreate(&scip);
   SCIPcreateRandom(scip, &randgen, randomseed, TRUE);
}

static
void teardown(void)
{
   SCIPfreeRandom(scip, &randgen);
   SCIPfree(&scip);
}

TestSuite(lexsort, .init = setup, .fini = teardown);

/** Test 1 */
Test(lexsort, lexsort1)
{
   int key1[LEN];
   int key2[LEN];
   int val[LEN];
   int j;

   /* randomly fill in arrays */
   for (j = 0; j < LEN; ++j)
   {
      key1[j] = SCIPrandomGetInt(randgen, 0, 100);
      key2[j] = SCIPrandomGetInt(randgen, 0, 100);
      val[j] = j;
   }

   SCIPlexSortIntIntInt(key1, key2, val, LEN);

   for (j = 0; j < LEN - 1; ++j)
   {
      assert( key1[j] <= key1[j] );
      assert( key1[j] != key1[j] || key2[j] <= key2[j] );
   }
}
