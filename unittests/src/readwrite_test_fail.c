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

/**@file   readwrite_test_fail.c
 * @brief  unit test for checking reading and writing of MISDPs in SDPA format
 * @author Marc Pfetsch
 *
 *  @todo: In order for this test to pass correctly, the calls to SCIPABORT() in reader_sdpa.c (and also reader_cbf.c)
 *  need to be removed. Otherwise SCIPABORT() breaks the program before the retcode SCIP_READERROR is returned.
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
Test(readwrite, readSDPAfail)
{
   /* SCIP_Real obj1; */
   /* SCIP_Real obj2; */
   SCIP_RETCODE retcode = 0;

   /* read problem, this should fail with a SCIP_READERROR. */
   retcode = SCIPreadProb(scipsdp, "example_small_fail.dat-s", NULL);

   cr_expect(retcode == SCIP_READERROR);
}
