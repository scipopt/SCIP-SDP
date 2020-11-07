/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2019 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2019 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   readwrite.cpp
 * @brief  unit test for checking reading and writing of MISDPs
 * @author Marc Pfetsch
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

TestSuite(readwrite, .init = setup, .fini = teardown);


/** TESTS **/

/** Test 1 */
Test(readwrite, readSDPAwriteCBF)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_small.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test1.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test1.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test1.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test1.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (SDPA) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (SDPA) != %g (CBF)\n", obj1, obj2);
}

/** Test 2 */
Test(readwrite, short2)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_small_cbf.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test2.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test2.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test2.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test2.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 3 */
Test(readwrite, short3)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_inf.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test3.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test3.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test3.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test3.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}

/** Test 4 */
Test(readwrite, short4)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_TT.dat-s.gz", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test4.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test4.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test4.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test4.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 5 */
Test(readwrite, short5)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_CLS.dat-s.gz", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test5.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test5.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test5.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test5.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 6 */
Test(readwrite, short6)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_MkP.dat-s.gz", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test6.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test6.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test6.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test6.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 7 */
Test(readwrite, short7)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_cbf_primal.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test7.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test7.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test7.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test7.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 8 */
Test(readwrite, short8)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_cbf_mix.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test8.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test8.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test8.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test8.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 9 */
Test(readwrite, short9)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_cbf_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test9.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test9.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test9.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test9.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 10 */
Test(readwrite, short10)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_multaggr.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test10.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test10.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test10.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test10.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 11 */
Test(readwrite, short11)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_primal.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test11.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test11.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test11.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test11.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, -obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 12 */
Test(readwrite, short12)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test12.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test12.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test12.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test12.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, -obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 13 */
Test(readwrite, short13)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_diagzeroimpl.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test13.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test13.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test13.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test13.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}


/** Test 14 */
Test(readwrite, nolinconsdual)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "nolincons_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test14.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test14.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test14.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test14.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}

/** Test 15 */
Test(readwrite, nolinconsprimal)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "nolincons_primal.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test15.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test15.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test15.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test15.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}



/** Test 16 */
Test(readwrite, nopsdcons)
{
   SCIP_Real obj1;
   SCIP_Real obj2;
   SCIP_Real obj3;

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "nopsdcons.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test16.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test16.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);
   
   /* write problem in SDPA format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "test16.dat-s", "dat-s", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "test16.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj3 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj3, EPS, "Optimal values differ: %g (Read Orig) != %g (SDPA)\n", obj1, obj3);
   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (Read Orig) != %g (Read Written)\n", obj1, obj2);
}
