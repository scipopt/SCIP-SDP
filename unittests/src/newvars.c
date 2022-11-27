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

/**@file   newvars.c
 * @brief  unit test for checking handling of new variables
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/cons_linear.h>
#include "scipsdp/scipsdpdefplugins.h"
#include "include/scip_test.h"
#include "scipsdp/cons_sdp.h"

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
   SCIP_CALL( SCIPaddCons(scipsdp, newcons) );

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

/** Test 2 */
Test(newvars, addSDPcons)
{
   SCIP_CONS* newcons;
   SCIP_CONS* newSDPcons;
   SCIP_VAR* newvar1;
   SCIP_VAR* newvar2;
   SCIP_VAR* newvar3;
   SCIP_Real* sdpvals[3];
   SCIP_Real trueobj = -8.0;
   SCIP_Real obj;
   SCIP_Real constval = -1.0;
   SCIP_VAR* vars[3];
   SCIP_Real vals[3] = {1.0, 1.0, 1.0};
   int nvarnonz[3] = {1, 1, 1};
   int constcol = 0;
   int constrow = 0;
   int* cols[3];
   int* rows[3];
   int i;

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
   SCIP_CALL( SCIPaddCons(scipsdp, newcons) );

   for (i = 0; i < 3; i++)
   {
      SCIP_CALL( SCIPallocBufferArray(scipsdp, &cols[i], 1) );
      SCIP_CALL( SCIPallocBufferArray(scipsdp, &rows[i], 1) );
      SCIP_CALL( SCIPallocBufferArray(scipsdp, &sdpvals[i], 1) );
   }

   cols[0][0] = 0;
   cols[1][0] = 0;
   cols[2][0] = 1;

   rows[0][0] = 0;
   rows[1][0] = 1;
   rows[2][0] = 1;

   sdpvals[0][0] = 1.0;
   sdpvals[1][0] = 1.0;
   sdpvals[2][0] = 1.0;


   SCIP_CALL( SCIPcreateConsSdp(scipsdp, &newSDPcons, "newSDPcons", 3, 3, 2, nvarnonz,
         cols, rows, sdpvals, vars, 1, &constcol, &constrow, &constval, FALSE) );
   SCIP_CALL( SCIPaddCons(scipsdp, newSDPcons) );

   for (i = 0; i < 3; i++)
   {
      SCIPfreeBufferArray(scipsdp, &sdpvals[i]);
      SCIPfreeBufferArray(scipsdp, &rows[i]);
      SCIPfreeBufferArray(scipsdp, &cols[i]);
   }

   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar3) );
   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar2) );
   SCIP_CALL( SCIPreleaseVar(scipsdp, &newvar1) );
   SCIP_CALL( SCIPreleaseCons(scipsdp, &newcons) );
   SCIP_CALL( SCIPreleaseCons(scipsdp, &newSDPcons) );

   /* presolve again */
   SCIP_CALL( SCIPpresolve(scipsdp) );

   /* solve */
   SCIP_CALL( SCIPsolve(scipsdp) );

   obj = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj, trueobj, EPS, "Optimal values differ: %g != %g (true value)\n", obj, trueobj);
}


