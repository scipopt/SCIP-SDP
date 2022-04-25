/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2021 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2021 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define SCIP_DEBUG */
/* #define SCIP_MORE_DEBUG */        /* enable full solver output and expected as well as actual feas status */

/**@file   check1dsdp.c
 * @brief  unit test for checking 1 variable SDP solving
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <sdpi/sdpi.h>
#include <scip/message.h>
#include <scip/message_default.h>

#include "include/scip_test.h"
#include "sdpi/solveonevarsdp.h"

#define EPS 1e-6

/* global variables */
static SCIP_SDPI* sdpi = NULL;
static BMS_BLKMEM* blockmem = NULL;
static BMS_BUFMEM* buffermem = NULL;
static SCIP_MESSAGEHDLR* messagehdlr = NULL;


/* macro for parameters */
#define SCIP_CALL_PARAM(x) /*lint -e527 */ do                                                   \
{                                                                                               \
   SCIP_RETCODE _restat_;                                                                       \
   if ( (_restat_ = (x)) != SCIP_OKAY && (_restat_ != SCIP_PARAMETERUNKNOWN) )                  \
   {                                                                                            \
      SCIPerrorMessage("[%s:%d] Error <%d> in function call\n", __FILE__, __LINE__, _restat_);  \
      abort();                                                                                  \
   }                                                                                            \
}                                                                                               \
while ( FALSE )

/** setup of test suite */
static
void setup(void)
{
   blockmem = BMScreateBlockMemory(1, 10);
   buffermem = BMScreateBufferMemory(SCIP_DEFAULT_MEM_ARRAYGROWFAC, SCIP_DEFAULT_MEM_ARRAYGROWINIT, FALSE);
   cr_assert( blockmem != NULL );
   cr_assert( buffermem != NULL );

   /* create message handler */
   SCIP_CALL( SCIPcreateMessagehdlrDefault(&messagehdlr, TRUE, NULL, FALSE) );

   /* create SDPI */
   SCIP_CALL( SCIPsdpiCreate(&sdpi, messagehdlr, blockmem, buffermem) );

   /* set feasibility tolerance */
   SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_SDPSOLVERFEASTOL, EPS) );
   SCIP_CALL( SCIPsdpiSetRealpar(sdpi, SCIP_SDPPAR_GAPTOL, EPS) );

#ifdef SCIP_DEBUG
   /* turn on output */
   SCIP_CALL( SCIPsdpiSetIntpar(sdpi, SCIP_SDPPAR_SDPINFO, 1) );
#endif
}

/** deinitialization method of test */
static
void teardown(void)
{
   /* release message handler */
   SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );

   SCIP_CALL( SCIPsdpiFree(&sdpi) );

   BMSdestroyBufferMemory(&buffermem);
   BMSdestroyBlockMemory(&blockmem);

   cr_assert_eq(BMSgetMemoryUsed(), 0, "There is a memory leak!");
}

TestSuite(checksdpi, .init = setup, .fini = teardown);



/** TESTS **/

/** Test 1
 *
 * min x1
 *     1 - x1 <= 0    <=>  x1 - 1 >= 0
 *     0 <= x1 <= 1
 *
 * with optimal solution (1).
 */
Test(checksdpi, test1)
{
   /* data with fixed values: */
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.0;
   SCIP_Real rhs = 1.0;
   SCIP_Real lhs = 1;
   int row = 0;
   int col = 0;
   SCIP_Real val = 1.0;
   int ncols = 1;
   int nrows = 1;
   int nnonz = 1;

   int sdpconstrow[1] = {0};
   int sdpconstcol[1] = {0};
   SCIP_Real sdpconstval[1] = {1.0};
   int sdprow[1] = {0};
   int sdpcol[1] = {0};
   SCIP_Real sdpval[1] = {1.0};

   /* checking result */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   SCIP_Real objval;
   SCIP_Real dualsol;

   rhs = SCIPsdpiInfinity(sdpi);

   /* load LP data, but leave SDP block empty */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, ncols, &obj, &lb, &ub, 0, NULL, NULL, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL,
         nrows, &lhs, &rhs, nnonz, &row, &col, &val) );

   /* solve problem: no Slater-check, no time limit */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, 1e20, NULL, NULL) );

   /* check feasibility status */
   SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, &primalfeasible, &dualfeasible) );
   cr_assert( primalfeasible );
   cr_assert( dualfeasible );

   /* get solution */
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objval, &dualsol, &ncols) );
   cr_assert( SCIPsdpiIsDualFeasible(sdpi) );

   cr_assert_float_eq(dualsol, 1.0, EPS, "Violation of dual solution: %g != %g\n", dualsol, 1.0);

   /* directly solve problem */
   SCIP_CALL( SCIPsolveOneVarSDP(buffermem, obj, lb, ub, 2, 1, sdpconstrow, sdpconstcol, sdpconstval, 1, sdprow, sdpcol, sdpval, 1e20, EPS, &objval, &dualsol) );

   cr_assert_float_eq(dualsol, 1.0, EPS, "Violation of 1d solution: %g != %g\n", dualsol, 1.0);
}

/** Test 2
 *
 * min x1
 *     x1 >= 1          <=>  x1 -1 >= 0     <=>   (1  0) x1 - (1  0) >= 0
 *     x1 <= 1              -x1 +1 >= 0           (0 -1)      (0 -1)
 *     0 <= x1 <= 1.5
 *
 * with optimal solution (1).
 */
Test(checksdpi, test2)
{
   /* data with fixed values: */
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 1.5;
   SCIP_Real rhs[2];
   SCIP_Real lhs[2];
   int row[2] = {0, 1};
   int col[2] = {0, 0};
   SCIP_Real val[2] = {1.0, 1.0};
   int ncols = 1;
   int nrows = 2;
   int nnonz = 2;

   int sdpconstrow[2] = {0, 1};
   int sdpconstcol[2] = {0, 1};
   SCIP_Real sdpconstval[2] = {1.0, -1.0};
   int sdprow[2] = {0, 1};
   int sdpcol[2] = {0, 1};
   SCIP_Real sdpval[2] = {1.0, -1.0};

   /* checking result */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   SCIP_Real objval;
   SCIP_Real dualsol;

   rhs[0] = SCIPsdpiInfinity(sdpi);
   rhs[1] = 1.0;
   lhs[0] = 1.0;
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   /* load LP data, but leave SDP block empty */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, ncols, &obj, &lb, &ub, 0, NULL, NULL, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL,
         nrows, lhs, rhs, nnonz, row, col, val) );

   /* solve problem: no Slater-check, no time limit */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, 1e20, NULL, NULL) );

   /* check feasibility status */
   SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, &primalfeasible, &dualfeasible) );
   cr_assert( primalfeasible );
   cr_assert( dualfeasible );

   /* get solution */
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objval, &dualsol, &ncols) );
   cr_assert( SCIPsdpiIsDualFeasible(sdpi) );

   cr_assert_float_eq(dualsol, 1.0, EPS, "Violation of dual solution: %g != %g\n", dualsol, 1.0);

   /* directly solve problem */
   SCIP_CALL( SCIPsolveOneVarSDP(buffermem, obj, lb, ub, 2, 2, sdpconstrow, sdpconstcol, sdpconstval, 2, sdprow, sdpcol, sdpval, 1e20, EPS, &objval, &dualsol) );

   cr_assert_float_eq(dualsol, 1.0, EPS, "Violation of 1d solution: %g != %g\n", dualsol, 1.0);
}

/** Test 3
 *
 * min x1
 *     x1 >= 1             <=>   x1 - 1    >= 0   <=>  (1  0) x1 - (1  0   ) >= 0
 *     x1 <= 0.99               -x1 + 0.99 >= 0        (0 -1)      (0 -0.99)
 *     0 <= x1 <= 2
 *
 * which is infeasible.
 */
Test(checksdpi, test3)
{
   /* data with fixed values: */
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 2.0;
   SCIP_Real rhs[2];
   SCIP_Real lhs[2];
   int row[2] = {0, 1};
   int col[2] = {0, 0};
   SCIP_Real val[2] = {1.0, 1.0};
   int ncols = 1;
   int nrows = 2;
   int nnonz = 2;

   int sdpconstrow[2] = {0, 1};
   int sdpconstcol[2] = {0, 1};
   SCIP_Real sdpconstval[2] = {1.0, -0.99};
   int sdprow[2] = {0, 1};
   int sdpcol[2] = {0, 1};
   SCIP_Real sdpval[2] = {1.0, -1.0};

   /* checking result */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   SCIP_Real objval;
   SCIP_Real dualsol;

   rhs[0] = SCIPsdpiInfinity(sdpi);
   rhs[1] = 0.99;
   lhs[0] = 1.0;
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   /* load LP data, but leave SDP block empty */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, ncols, &obj, &lb, &ub, 0, NULL, NULL, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL,
         nrows, lhs, rhs, nnonz, row, col, val) );

   /* solve problem: no Slater-check, no time limit */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, 1e20, NULL, NULL) );

   /* check feasibility status */
   SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, &primalfeasible, &dualfeasible) );
   cr_assert( ! dualfeasible );

   /* get solution */
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objval, &dualsol, &ncols) );

   cr_assert( SCIPsdpiIsDualInfeasible(sdpi) );

   /* directly solve problem */
   SCIP_CALL( SCIPsolveOneVarSDP(buffermem, obj, lb, ub, 2, 2, sdpconstrow, sdpconstcol, sdpconstval, 2, sdprow, sdpcol, sdpval, 1e20, EPS, &objval, &dualsol) );

   cr_assert( objval >= 1e20 );
}

/** Test 4
 *
 *  Same as Test 3, but the matrices are multiplied by the orthonormal matrix
 *  V = [  0.97325  -0.22975 ]
 *      [ -0.22975  -0.97325 ]
 */
Test(checksdpi, test4)
{
   /* data with fixed values: */
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 2.0;

   int sdpconstrow[3] = {0, 1, 1};
   int sdpconstcol[3] = {0, 0, 1};
   SCIP_Real sdpconstval[3] = {0.89496, -0.44498, -0.88496};
   int sdprow[3] = {0, 1, 1};
   int sdpcol[3] = {0, 0, 1};
   SCIP_Real sdpval[3] = {0.89443, -0.44721, -0.89443};

   /* checking result */
   SCIP_Real objval;
   SCIP_Real dualsol;

   /* directly solve problem */
   SCIP_CALL( SCIPsolveOneVarSDP(buffermem, obj, lb, ub, 2, 3, sdpconstrow, sdpconstcol, sdpconstval, 3, sdprow, sdpcol, sdpval, 1e20, EPS, &objval, &dualsol) );

   cr_assert(objval >= 1e20);
}

/** Test 5
 *
 *  A = [1,-2;-2,5] with eigenvalue (0.17157, 5.82843)
 *  B = [-2,1;1,3] with eigenvalues (-2.1926, 3.1926)
 */
Test(checksdpi, test5)
{
   /* data with fixed values: */
   SCIP_Real obj = 1.0;
   SCIP_Real lb = 0.0;
   SCIP_Real ub = 2.0;

   int sdpconstrow[3] = {0, 1, 1};
   int sdpconstcol[3] = {0, 0, 1};
   SCIP_Real sdpconstval[3] = {-2.0, 1.0, 3.0};
   int sdprow[3] = {0, 1, 1};
   int sdpcol[3] = {0, 0, 1};
   SCIP_Real sdpval[3] = {1.0, -2.0, 5.0};

   int sdpblocksizes = 2;
   int sdpnblockvars = 1;
   int sdpconstnblocknonz = 3;
   int* sdpconstrows;
   int* sdpconstcols;
   SCIP_Real* sdpconstvals;
   int sdpnblockvarnonz = 3;
   int* sdpnblockvarnonzs;
   int sdpvar = 0;
   int* sdpvars;
   int* sdprows;
   int* sdpcols;
   SCIP_Real* sdpvals;
   int** sdprowss;
   int** sdpcolss;
   SCIP_Real** sdpvalss;

   /* checking result */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   SCIP_Real objval;
   SCIP_Real dualsol;
   int ncols = 1;

   sdpconstrows = (int*) sdpconstrow;
   sdpconstcols = (int*) sdpconstcol;
   sdpconstvals = (SCIP_Real*) sdpconstval;
   sdpnblockvarnonzs = &sdpnblockvarnonz;

   sdpvars = &sdpvar;
   sdprows = &sdprow[0];
   sdpcols = &sdpcol[0];
   sdpvals = &sdpval[0];

   sdprowss = &sdprows;
   sdpcolss = &sdpcols;
   sdpvalss = &sdpvals;

   /* load SDP data */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, 1, &obj, &lb, &ub, 1, &sdpblocksizes, &sdpnblockvars, 3, &sdpconstnblocknonz, &sdpconstrows, &sdpconstcols, &sdpconstvals,
         3, &sdpnblockvarnonzs, &sdpvars, &sdprowss, &sdpcolss, &sdpvalss, 0, NULL, NULL, 0, NULL, NULL, NULL) );

   /* solve problem: no Slater-check, no time limit */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, 1e20, NULL, NULL) );

   /* check feasibility status */
   SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, &primalfeasible, &dualfeasible) );
   cr_assert( primalfeasible );
   cr_assert( dualfeasible );

   /* get solution */
   SCIP_CALL( SCIPsdpiGetSol(sdpi, &objval, &dualsol, &ncols) );
   cr_assert_float_eq(dualsol, 1.541381, EPS, "Violation of dual solution: %g != %g\n", dualsol,  1.541381);

   /* directly solve problem */
   SCIP_CALL( SCIPsolveOneVarSDP(buffermem, obj, lb, ub, 2, 3, sdpconstrow, sdpconstcol, sdpconstval, 3, sdprow, sdpcol, sdpval, 1e20, EPS, &objval, &dualsol) );
   cr_assert_float_eq(dualsol, 1.541381, EPS, "Violation of dual solution: %g != %g\n", dualsol, 1.541381);
}
