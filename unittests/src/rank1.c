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

/**@file   rank1.c
 * @brief  unit test for checking reading of rank1-information in CBF-files
 * @author Marc Pfetsch
 * @author Frederic Matter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scipsdp/scipsdpdefplugins.h"
#include "include/scip_test.h"

/* global SCIP data structure */
SCIP* scipsdp;

#define EPS  1e-5

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

TestSuite(rank1, .init = setup, .fini = teardown);


/** TESTS **/

/** Test 1 */
Test(rank1, readPrimalDualRank1)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem in primal form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_primal.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* read problem in dual form and solve it again */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (CBF primal) != %g (CBF dual)\n", obj1, obj2);
}

/** Test 2 */
Test(rank1, readWritePrimalRank)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem in primal form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_primal.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "example_rank1.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "example_rank1.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (CBF primal original) != %g (CBF primal after writing)\n", obj1, obj2);
}

/** Test 3 */
Test(rank1, readWriteDualRank)
{
   SCIP_Real obj1;
   SCIP_Real obj2;

   /* read problem in dual form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_rank1_dual.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "example_rank1.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "example_rank1.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (CBF dual original) != %g (CBF dual after writing)\n", obj1, obj2);
}
