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

/**@file   readerrors.c
 * @brief  unit test for checking error while reading MISDPs in SDPA format
 * @author Tim Schmidt
 * @author Frederic Matter
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <signal.h>
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
Test(readerrors, blocks_col, setup, teardown)
{
   SCIP_RETCODE retcode;

   /* read problem, this should fail with a SCIP_READERROR. */
   retcode = SCIPreadProb(scipsdp, "./instances/blocks_col.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_invalidline, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_invalidline.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocksizes_0, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocksizes_0.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocksizes_invalidsymb, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocksizes_invalidsymb.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocksizes_LPblocks, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocksizes_LPblocks.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocksizes_toofew, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocksizes_toofew.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_LPnononz, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_LPnononz.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_row, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_row.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_sdpblock, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_sdpblock.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_SDPnononz, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_SDPnononz.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, blocks_var, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/blocks_var.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, int_invalid, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/int_invalid.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, int_noast, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/int_noast.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, int_var, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/int_var.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, LPblock_LPcons, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/LPblock_LPcons.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, LPblock_nondiag, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/LPblock_nondiag.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, LPblock_var, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/LPblock_var.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, nblocks_invalidsymb, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/nblocks_invalidsymb.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, nblocks_neg, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/nblocks_neg.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, nvars_invalidsymb, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/nvars_invalidsymb.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, nvars_neg, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/nvars_neg.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, objcoeff_invalidsymb, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/objcoeff_invalidsymb.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, objcoeff_toofew, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/objcoeff_toofew.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, rnk1_before_int, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/rnk1_before_int.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, rnk1_block, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/rnk1_block.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, rnk1_forLP, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/rnk1_forLP.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, rnk1_invalid, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/rnk1_invalid.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}

Test(readerrors, rnk1_noast, setup, teardown)
{
   SCIP_RETCODE retcode;

   retcode = SCIPreadProb(scipsdp, "./instances/rnk1_noast.dat-s", NULL);
   cr_expect(retcode == SCIP_READERROR);
}
