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

/**@file   readwrite.c
 * @brief  unit test for checking reading and writing of MISDPs
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "objscip/objscipdefplugins.h"

#include "include/scip_test.h"
#include "scipsdp/cons_sdp.h"
#include "scipsdp/cons_savesdpsol.h"
#include "scipsdp/cons_savedsdpsettings.h"
#include "scipsdp/relax_sdp.h"
#include "scipsdp/objreader_sdpa.h"
#include "scipsdp/reader_cbf.h"
#include "scipsdp/prop_sdpredcost.h"
#include "scipsdp/disp_sdpiterations.h"
#include "scipsdp/disp_sdpavgiterations.h"
#include "scipsdp/disp_sdpfastsettings.h"
#include "scipsdp/disp_sdppenalty.h"
#include "scipsdp/disp_sdpunsolved.h"
#include "scipsdp/branch_sdpmostfrac.h"
#include "scipsdp/branch_sdpmostinf.h"
#include "scipsdp/branch_sdpobjective.h"
#include "scipsdp/branch_sdpinfobjective.h"
#include "scipsdp/heur_sdpfracdiving.h"
#include "scipsdp/heur_sdprand.h"
#include "scipsdp/prop_sdpobbt.h"
#include "scipsdp/prop_companalcent.h"
#include "scipsdp/table_relaxsdp.h"
#include "scipsdp/table_sdpsolversuccess.h"
#include "scipsdp/table_slater.h"

/* global SCIP data structure */
SCIP* scipsdp;

#define EPS  1e-6

using namespace scip;


/** setup of test suite */
static
void setup(void)
{
   SCIP_CALL( SCIPcreate(&scipsdp) );

   /* include new plugins */
   SCIP_CALL( SCIPincludeObjReader(scipsdp, new ObjReaderSDPA(scipsdp), TRUE) );
   SCIP_CALL( SCIPincludeReaderCbf(scipsdp) );
   SCIP_CALL( SCIPincludeConshdlrSdp(scipsdp) );
   SCIP_CALL( SCIPincludeConshdlrSavesdpsol(scipsdp) );
   SCIP_CALL( SCIPincludeConshdlrSavedsdpsettings(scipsdp) );
   SCIP_CALL( SCIPincludeRelaxSdp(scipsdp) );
   SCIP_CALL( SCIPincludePropSdpredcost(scipsdp) );
   SCIP_CALL( SCIPincludeBranchruleSdpmostfrac(scipsdp) );
   SCIP_CALL( SCIPincludeBranchruleSdpmostinf(scipsdp) );
   SCIP_CALL( SCIPincludeBranchruleSdpobjective(scipsdp) );
   SCIP_CALL( SCIPincludeBranchruleSdpinfobjective(scipsdp) );
   SCIP_CALL( SCIPincludeHeurSdpFracdiving(scipsdp) );
   SCIP_CALL( SCIPincludeHeurSdpRand(scipsdp) );
   SCIP_CALL( SCIPincludePropSdpObbt(scipsdp) );
   SCIP_CALL( SCIPincludePropCompAnalCent(scipsdp) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scipsdp) );

   /* Choose between LP and SDP relaxations */
   SCIP_CALL( SCIPsetIntParam(scipsdp, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(scipsdp, "relaxing/SDP/freq", 1) );

   /* change epsilons for numerical stability */
   SCIP_CALL( SCIPsetRealParam(scipsdp, "numerics/epsilon", 1e-9) );
   SCIP_CALL( SCIPsetRealParam(scipsdp, "numerics/sumepsilon", 1e-6) );
   SCIP_CALL( SCIPsetRealParam(scipsdp, "numerics/feastol", 1e-6) );

   /* parameters for separation */
   SCIP_CALL( SCIPsetBoolParam(scipsdp, "lp/cleanuprows", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scipsdp, "lp/cleanuprowsroot", FALSE) );

   SCIP_CALL( SCIPsetIntParam(scipsdp, "nodeselection/hybridestim/stdpriority", 1000000) );
   SCIP_CALL( SCIPsetIntParam(scipsdp, "nodeselection/hybridestim/maxplungedepth", 0) );
   SCIP_CALL( SCIPsetRealParam(scipsdp, "nodeselection/hybridestim/estimweight", 0.0) );
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

   /* read problem and solve it */
   SCIP_CALL( SCIPreadProb(scipsdp, "../instances/example_small.dat-s", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj1 = SCIPgetDualbound(scipsdp);

   /* write problem in CBF format */
   SCIP_CALL( SCIPwriteOrigProblem(scipsdp, "example_small.cbf", "cbf", FALSE) );

   /* read problem again */
   SCIP_CALL( SCIPreadProb(scipsdp, "example_small.cbf", NULL) );

   SCIP_CALL( SCIPsolve(scipsdp) );

   obj2 = SCIPgetDualbound(scipsdp);

   cr_assert_float_eq(obj1, obj2, EPS, "Optimal values differ: %g (SDPA) != %g (CBF)\n", obj1, obj2);
}
