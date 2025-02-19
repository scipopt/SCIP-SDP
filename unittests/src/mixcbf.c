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

/**@file   mixcbf.cpp
 * @brief  unit test for checking reading MISDPs in primal, dual and mixed form
 * @author Marc Pfetsch
 * @author Frederic Matter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scipsdp/scipsdpdefplugins.h"
#include "include/scip_test.h"

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
   /* deinitialization */
   SCIP_CALL( SCIPfree(&scipsdp) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(readmix, .init = setup, .fini = teardown);

/** TESTS **/

/** Test 1 */
Test(readmix, readCBFmixreadCBFdual)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem in mixed form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_cbf_mix.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* read the same problem in dual form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_cbf_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (CBF Mixed) != %g (CBF dual)\n", obj1, obj2);
}
