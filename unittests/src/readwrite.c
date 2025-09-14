/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2025 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2025 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   readwrite.cpp
 * @brief  unit test for checking reading and writing of MISDPs in CBF and SDPA format
 * @author Tim Schmidt
 * @author Frederic Matter
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scipsdp/scipsdpdefplugins.h"
#include "include/scip_test.h"

/* global SCIP data structure */
SCIP* scipsdp;

#define EPS       1e-6
#define RANK1EPS  1e-5


//! macro to check for a SCIP error and possibly exit
#define SCIP_CALL_STOP(x)  do                                  \
                       {                                       \
                          SCIP_RETCODE _restat_;               \
                          if( (_restat_ = (x)) != SCIP_OKAY )  \
                          {                                    \
                             SCIPprintError(_restat_);         \
                             abort();                          \
                           }                                   \
                       }                                       \
                       while( FALSE )


/** run test for CBF and DAT-S */
static
SCIP_RETCODE runTests(
   const char*           path,               /**< path to testfile */
   const char*           basename,           /**< basename of testfile */
   const char*           extension,          /**< extension of testfile */
   int                   idx,                /**< index of test */
   SCIP_Real             eps,                /**< epsilon for testing */
   SCIP_Real             objsense,           /**< objective sense for SDPA (can only write minimization problems) */
   SCIP_Bool             testpresol          /**< whether writing/reading the presolved problem should be tested */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;
   SCIP_Real obj4;

   /* read problem and solve it */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s/%s.%s", path, basename, extension);
   SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );

   /* write problem in CBF format */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d.cbf", idx);
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, name, "cbf", FALSE) );

   /* write problem in SDPA formata */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d.dat-s", idx);
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, name, "dat-s", FALSE) );

   /* write presolved problem in CBF format */
   if ( testpresol )
   {
      SCIP_CALL( SCIPpresolve(scipsdp) );
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d-pre.cbf", idx);
      SCIP_CALL( SCIPwriteTransProblem(scipsdp, name, "cbf", FALSE) );
   }

   /* now solve */
   SCIP_CALL( SCIPsolve(scipsdp) );
   obj1 = SCIPgetDualbound(scipsdp);

   /* read CBF problem again */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d.cbf", idx);
   SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );
   SCIP_CALL( SCIPsolve(scipsdp) );
   obj2 = SCIPgetDualbound(scipsdp);

   /* read SDPA problem again */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d.dat-s", idx);
   SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );
   SCIP_CALL( SCIPsolve(scipsdp) );
   obj3 = SCIPgetDualbound(scipsdp);

   /* read CBF presolved problem again */
   if ( testpresol )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d-pre.cbf", idx);
      SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );
      SCIP_CALL( SCIPsolve(scipsdp) );
      obj4 = SCIPgetDualbound(scipsdp);
   }

   cr_assert_float_eq(obj1, obj2, eps, "Optimal values differ: %g (SDPA original) != %g (CBF written)\n", obj1, obj2);
   cr_assert_float_eq(obj1, objsense * obj3, eps, "Optimal values differ: %g (SDPA original) != %g (SDPA written)\n", obj1, obj3);
   if ( testpresol )
   {
      cr_assert_float_eq(obj1, obj4, eps, "Optimal values differ: %g (SDPA original) != %g (presolved CBF written)\n", obj1, obj4);
   }

   return SCIP_OKAY;
}


/** run test for CIP */
static
SCIP_RETCODE runTestsCIP(
   const char*           path,               /**< path to testfile */
   const char*           basename,           /**< basename of testfile */
   const char*           extension,          /**< extension of testfile */
   int                   idx                 /**< index of test */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem and solve it */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s/%s.%s", path, basename, extension);
   SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );
   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CIP format */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "test%d.cip", idx);
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, name, "cip", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, name, NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );
   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (SDPA original) != %g (CIP written)\n", obj1, obj2);

   return SCIP_OKAY;
}


/** setup of test suite */
static
void setup(void)
{
   SCIP_CALL_STOP( SCIPcreate(&scipsdp) );

   /* include default SCIP-SDP plugins */
   SCIP_CALL_STOP( SCIPSDPincludeDefaultPlugins(scipsdp) );
}

/** deinitialization method of test */
static
void teardown(void)
{
   /* deinitialization */
   SCIP_CALL_STOP( SCIPfree(&scipsdp) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(readwrite, .init = setup, .fini = teardown);


/** TESTS **/

/** Test 1 */
Test(readwrite, test1)
{
   SCIP_CALL_STOP( runTests("../instances", "example_small", "dat-s", 1, EPS, 1, TRUE) );
}

/** Test 2 */
Test(readwrite, test2)
{
   SCIP_CALL_STOP( runTests("../instances", "example_small_cbf", "cbf", 2, EPS, 1, TRUE) );
}

/** Test 3 */
Test(readwrite, test3)
{
   SCIP_CALL_STOP( runTests("../instances", "example_inf", "dat-s", 3, EPS, 1, TRUE) );
}

/** Test 4 */
Test(readwrite, test4)
{
   SCIP_CALL_STOP( runTests("../instances", "example_TT", "dat-s.gz", 4, EPS, 1, TRUE) );
}

/** Test 5 */
Test(readwrite, test5)
{
   SCIP_CALL_STOP( runTests("../instances", "example_CLS", "dat-s.gz", 5, EPS, 1, TRUE) );
}

/** Test 6 */
Test(readwrite, test6)
{
   SCIP_CALL_STOP( runTests("../instances", "example_MkP", "dat-s.gz", 6, EPS, 1, TRUE) );
}

/** Test 7 */
Test(readwrite, test7)
{
   SCIP_CALL_STOP( runTests("../instances", "example_cbf_primal", "cbf", 7, EPS, 1, TRUE) );
}

/** Test 8 */
Test(readwrite, test8)
{
   SCIP_CALL_STOP( runTests("../instances", "example_cbf_mix", "cbf", 8, EPS, 1, TRUE) );
}

/** Test 9 */
Test(readwrite, test9)
{
   SCIP_CALL_STOP( runTests("../instances", "example_cbf_dual", "cbf", 9, EPS, 1, TRUE) );
}

/** Test 10 */
Test(readwrite, test10)
{
   SCIP_CALL_STOP( runTests("../instances", "example_multaggr", "cbf", 10, EPS, 1, TRUE) );
}

/** Test 11 */
Test(readwrite, test11)
{
   SCIP_CALL_STOP( runTests("../instances", "example_rank1_primal", "cbf", 11, RANK1EPS, -1, FALSE) );
}

/** Test 12 */
Test(readwrite, test12)
{
   SCIP_CALL_STOP( runTests("../instances", "example_rank1_dual", "cbf", 12, RANK1EPS, -1, FALSE) );
}

/** Test 13 */
Test(readwrite, test13)
{
   SCIP_CALL_STOP( runTests("../instances", "example_diagzeroimpl", "cbf", 13, EPS, 1, TRUE) );
}

/** Test 14 */
Test(readwrite, test14)
{
   SCIP_CALL_STOP( runTests("instances", "nolincons_dual", "cbf", 14, EPS, 1, TRUE) );
}

/** Test 15 */
Test(readwrite, test15)
{
   SCIP_CALL_STOP( runTests("instances", "nolincons_primal", "cbf", 15, EPS, 1, TRUE) );
}

/** Test 16 */
Test(readwrite, test16)
{
   SCIP_CALL_STOP( runTests("instances", "nopsdcons", "cbf", 16, EPS, 1, TRUE) );
}

/** Test 17 */
Test(readwrite, test17)
{
   SCIP_CALL_STOP( runTestsCIP("../instances", "example_small_ind", "dat-s", 17) );
}

/** Test 18 */
Test(readwrite, test18)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "../lib/scip/check/instances/MIP/stein27.fzn", NULL) );

   /* write problem in CBF format */
   SCIP_CALL_STOP( SCIPwriteOrigProblem(scipsdp, "test18.cbf", "cbf", FALSE) );

   /* write presolved problem */
   SCIP_CALL_STOP( SCIPpresolve(scipsdp) );
   SCIP_CALL_STOP( SCIPwriteTransProblem(scipsdp, "test18-pre.cbf", "cbf", FALSE) );

   /* now solve */
   SCIP_CALL_STOP( SCIPsolve(scipsdp) );
   obj1 = SCIPgetDualbound(scipsdp);

   /* read CBF problem again */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "test18.cbf", NULL) );
   SCIP_CALL_STOP( SCIPsolve(scipsdp) );
   obj2 = SCIPgetDualbound(scipsdp);

   /* read CBF presolved problem again */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "test18-pre.cbf", NULL) );
   SCIP_CALL_STOP( SCIPsolve(scipsdp) );
   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (SDPA original) != %g (CBF written)\n", obj1, obj2);
   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (SDPA original) != %g (presolved CBF written)\n", obj1, obj3);
}

/** Test 19 */
Test(readwrite, sign)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem in CBF format with L+ */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "instances/example_small_L-.cbf", NULL) );

   /* write problem in SDPA format */
   SCIP_CALL_STOP( SCIPwriteOrigProblem(scipsdp, "example_small_L-.dat-s", "dat-s", FALSE) );

   /* read problem in CBF format with L- */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "instances/example_small_L+.cbf", NULL) );

   /* write problem in SDPA format */
   SCIP_CALL_STOP( SCIPwriteOrigProblem(scipsdp, "example_small_L+.dat-s", "dat-s", FALSE) );

   /* read problem with L- in SDPA format and solve it */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "example_small_L-.dat-s", NULL) );

   SCIP_CALL_STOP( SCIPsolve(scipsdp) );
   obj1 = SCIPgetDualbound(scipsdp);

   /* read problem with L+ in SDPA format and solve it */
   SCIP_CALL_STOP( SCIPreadProb(scipsdp, "example_small_L+.dat-s", NULL) );

   SCIP_CALL_STOP( SCIPsolve(scipsdp) );
   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (SDPA from CBF with L-) != %g (SDPA from CBF with L+)\n", obj1, obj2);
}
