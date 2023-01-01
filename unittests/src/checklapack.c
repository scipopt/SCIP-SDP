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

/**@file   checklapack.c
 * @brief  unit test for checking Lapack routines
 * @author Frederic Matter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include "include/scip_test.h"
#include "scipsdp/scipsdpdefplugins.h"
#include "sdpi/lapack_interface.h"

/* global SCIP data structure */
SCIP* scipsdp;

#define EPS  1e-6


/** setup of test suite */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scipsdp) );

   /* include default SCIP-SDP plugins */
   SCIP_CALL( SCIPSDPincludeDefaultPlugins(scipsdp) );
}

/** deinitialization method of test */
static
void teardown(void)
{
   SCIP_CALL( SCIPfree(&scipsdp) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(checklapack, .init = setup, .fini = teardown);



/** TESTS **/

/** Test 1
 *
 * Check matrix matrix multiplication:
 *
 * [1.0 3.0] * [5.0 7.0]^T = [26.0 30.0]
 * [2.0 4.0]   [6.0 8.0]     [38.0 44.0]
 *
 */
Test(checklapack, test1)
{
   /* data with fixed values: */
   int blocksize = 2;
   SCIP_Real* matrixA;
   SCIP_Real* matrixB;
   SCIP_Real* matrixC;

   SCIP_CALL( SCIPallocBufferArray(scipsdp, &matrixA, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scipsdp, &matrixB, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scipsdp, &matrixC, blocksize * blocksize) );

   /* Note: Lapack uses column-first format! */
   matrixA[0] = 1.0;
   matrixA[1] = 2.0;
   matrixA[2] = 3.0;
   matrixA[3] = 4.0;

   matrixB[0] = 5.0;
   matrixB[1] = 6.0;
   matrixB[2] = 7.0;
   matrixB[3] = 8.0;

   matrixC[0] = 1.0;
   matrixC[1] = 1.0;
   matrixC[2] = 1.0;
   matrixC[3] = 1.0;

   SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksize, blocksize, matrixA, FALSE, blocksize, blocksize, matrixB,
         TRUE, matrixC) );

   cr_assert_float_eq(matrixC[0], 26.0, EPS, "C[0]: %g != %g\n", matrixC[0], 26.0);
   cr_assert_float_eq(matrixC[1], 38.0, EPS, "C[0]: %g != %g\n", matrixC[1], 38.0);
   cr_assert_float_eq(matrixC[2], 30.0, EPS, "C[0]: %g != %g\n", matrixC[2], 30.0);
   cr_assert_float_eq(matrixC[3], 44.0, EPS, "C[0]: %g != %g\n", matrixC[3], 44.0);

   SCIPfreeBufferArray(scipsdp, &matrixA);
   SCIPfreeBufferArray(scipsdp, &matrixB);
   SCIPfreeBufferArray(scipsdp, &matrixC);
}
