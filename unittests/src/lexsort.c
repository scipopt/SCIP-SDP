/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2023 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2023 Zuse Institute Berlin                             */
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
