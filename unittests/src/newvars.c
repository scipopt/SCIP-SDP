/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2020 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2020 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   newvars.c
 * @brief  unit test for checking handling of new variables
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
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

TestSuite(newvars, .init = setup, .fini = teardown);


/** TESTS **/

/** Test 1 */
Test(newvars, addvars)
{
   SCIP_CONS* newcons;
   SCIP_VAR* newvar1;
   SCIP_VAR* newvar2;
   SCIP_VAR* newvar3;
   SCIP_Real trueobj = -8.0;
   SCIP_Real obj;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3] = {1.0, 1.0, 1.0};

   /* read problem in primal form and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_small.dat-s", NULL) );

   /* presolve */
   SCIP_CALL( SCIPpresolve(scipsdp) );

   /* add a new variable */
   SCIP_CALL( SCIPcreateVarBasic(scipsdp, &newvar1, "newvar1", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scipsdp, newvar1) );

   SCIP_CALL( SCIPcreateVarBasic(scipsdp, &newvar2, "newvar2", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scipsdp, newvar2) );

   SCIP_CALL( SCIPcreateVarBasic(scipsdp, &newvar3, "newvar3", 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scipsdp, newvar3) );

   vars[0] = newvar1;
   vars[1] = newvar2;
   vars[2] = newvar3;

   SCIP_CALL( SCIPcreateConsBasicLinear(scipsdp, &newcons, "newcons", 3, vars, vals, 1.0, 1.0) );

   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar3) );
   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar2) );
   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar1) );
   SCIP_CALL( SCIPreleaseCons(scipsdp, &newcons) );

   /* presolve again */
   SCIP_CALL( SCIPpresolve(scipsdp) );

   /* solve */
   SCIP_CALL( SCIPsolve(scipsdp) );

   obj = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj, trueobj, EPS, "Optimal values differ: %g != %g (true value)\n", obj, trueobj);
}
