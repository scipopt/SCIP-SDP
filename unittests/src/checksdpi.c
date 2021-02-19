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

/* #define SCIP_MORE_DEBUG         /\* enable full solver output and expected as well as actual feas status *\/ */

/**@file   checksdpi.c
 * @brief  unit test for checking SDPI
 * @author Marc Pfetsch
 * @author Frederic Matter
 *
 * We perform tests with solving several examples. They are inspired by the tests for the LPI in SCIP.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <sdpi/sdpi.h>
#include <scip/message.h>
#include <scip/message_default.h>

#include "include/scip_test.h"

#define EPS 1e-6

/** expected feasibility status for primal or dual problem */
enum SCIPfeasStatus
{
   SCIPfeas      = 0,    /**< the problem is feasible */
   SCIPunbounded = 1,    /**< the problem is unbounded (and feasible) */
   SCIPray       = 2,    /**< there exists a ray, but it is not necessarily feasible */
   SCIPinfeas    = 3     /**< the problem is infeasible */
};
typedef enum SCIPfeasStatus SCIPFEASSTATUS;

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

/* local functions */

/** solve problem */
static
SCIP_RETCODE solveTest(
   int                   ncols,              /**< number of columns */
   int                   nrows,              /**< number of rows */
   SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
   SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
   SCIP_Real*            exp_dualsol         /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
   )
{
   /* solution data */
   SCIP_Real objval;
   SCIP_Real* dualsol;

   /* auxiliary data */
   SCIP_Bool primalfeasible;
   SCIP_Bool dualfeasible;
   int ntmprows;
   int ntmpcols;
   int j;

#ifdef SCIP_MORE_DEBUG
   int sdpinfo;

   SCIPmessageFPrintInfo(messagehdlr, NULL, "Expected primal feas status: %d\n", exp_primalfeas);
   SCIPmessageFPrintInfo(messagehdlr, NULL, "Expected dual feas status: %d\n", exp_dualfeas);

   /* enable full output */
   SCIP_CALL( SCIPsdpiSetIntpar(sdpi, 5, 1) );
   SCIP_CALL( SCIPsdpiGetIntpar(sdpi, 5, &sdpinfo) );
#endif

   /* check size */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &ntmprows) );
   SCIP_CALL( SCIPsdpiGetNVars(sdpi, &ntmpcols) );
   cr_assert( nrows == ntmprows );
   cr_assert( ncols == ntmpcols );

   /* solve problem: no Slater-check, no time limit */
   SCIP_CALL( SCIPsdpiSolve(sdpi, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, SCIP_SDPSOLVERSETTING_UNSOLVED, FALSE, 1e20) );

   /* check status */
   cr_assert( SCIPsdpiWasSolved(sdpi) );
   cr_assert( ! SCIPsdpiIsObjlimExc(sdpi) );
   cr_assert( ! SCIPsdpiIsIterlimExc(sdpi) );
   cr_assert( ! SCIPsdpiIsTimelimExc(sdpi) );

   /* check feasibility status */
   SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, &primalfeasible, &dualfeasible) );

#ifdef SCIP_MORE_DEBUG
   SCIPmessageFPrintInfo(messagehdlr, NULL, "Primal feasible?: %d\n", primalfeasible);
   SCIPmessageFPrintInfo(messagehdlr, NULL, "Dual feasible?: %d\n", dualfeasible);
#endif

   /* if we are feasible, we should be optimal */
   if ( exp_primalfeas == SCIPfeas && exp_dualfeas == SCIPfeas )
   {
      cr_assert( SCIPsdpiIsOptimal(sdpi) );

      cr_assert( SCIPsdpiIsDualFeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsDualInfeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsDualUnbounded(sdpi) );

      cr_assert( SCIPsdpiIsPrimalFeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsPrimalInfeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsPrimalUnbounded(sdpi) );
   }

   /* check more primal statuses */
   switch ( exp_primalfeas )
   {
   case SCIPfeas:
      cr_assert( primalfeasible );
      cr_assert( SCIPsdpiIsPrimalFeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsPrimalInfeasible(sdpi) );

      cr_assert( ! SCIPsdpiIsDualUnbounded(sdpi) );
      break;

   case SCIPunbounded:
      /* cr_assert( SCIPsdpiIsPrimalFeasible(sdpi) ); */ /* we do not know whether we can prove primal feasibility */
      cr_assert( ! SCIPsdpiIsPrimalInfeasible(sdpi) );
      /* cr_assert( SCIPsdpiIsPrimalUnbounded(sdpi) ); */ /* we do not know whether we can prove unboundedness in the primal */

      cr_assert( ! SCIPsdpiIsDualFeasible(sdpi) );
      cr_assert( SCIPsdpiIsDualInfeasible(sdpi) );
      break;

   case SCIPray:
      cr_assert( ! SCIPsdpiIsDualFeasible(sdpi) );
      cr_assert( SCIPsdpiIsDualInfeasible(sdpi) );
      break;

   case SCIPinfeas:
      cr_assert( ! primalfeasible );
      cr_assert( ! SCIPsdpiIsPrimalFeasible(sdpi) );
      /* cr_assert( SCIPsdpiIsPrimalInfeasible(sdpi) ); */ /* we do not know whether we can determine primal infeasibility in the primal */
      cr_assert( ! SCIPsdpiIsPrimalUnbounded(sdpi) );
      break;

   default:
      abort();
   }

   /* check more dual statuses */
   switch ( exp_dualfeas )
   {
   case SCIPfeas:
      cr_assert( dualfeasible );
      cr_assert( SCIPsdpiIsDualFeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsDualInfeasible(sdpi) );

      cr_assert( ! SCIPsdpiIsPrimalUnbounded(sdpi) );
      break;

   case SCIPunbounded:
      /* Since Mosek solves the primal problem instead of the dual problem and cannot detect unboundedness, we cannot
       * determine the status of the dual problem purely based on the solution status returned by Mosek. We only know
       * that the primal problem should be infeasible, since otherwise dual unboundedness is not possible due to weak
       * duality.
       */

      /* cr_assert( SCIPsdpiIsDualFeasible(sdpi) ); */
      /* cr_assert( ! SCIPsdpiIsDualInfeasible(sdpi) ); */
      /* cr_assert( SCIPsdpiIsDualUnbounded(sdpi) ); */

      cr_assert( ! SCIPsdpiIsPrimalFeasible(sdpi) );
      cr_assert( SCIPsdpiIsPrimalInfeasible(sdpi) );
      break;

   case SCIPray:
      cr_assert( ! SCIPsdpiIsPrimalFeasible(sdpi) );
      cr_assert( SCIPsdpiIsPrimalInfeasible(sdpi) );
      break;

   case SCIPinfeas:
      cr_assert( ! dualfeasible );
      cr_assert( ! SCIPsdpiIsDualFeasible(sdpi) );
      cr_assert( SCIPsdpiIsDualInfeasible(sdpi) );
      cr_assert( ! SCIPsdpiIsDualUnbounded(sdpi) );
      break;

   default:
      abort();
   }

   /* allocate storage for solution */
   BMSallocMemoryArray(&dualsol, nrows);

   if ( exp_dualfeas == SCIPfeas )
   {
      /* get solution */
      SCIP_CALL( SCIPsdpiGetSol(sdpi, &objval, dualsol, &ncols) );

      for (j = 0; j < ncols; ++j)
      {
         cr_assert_float_eq(dualsol[j], exp_dualsol[j], EPS, "Violation of dual solution %d: %g != %g\n", j, dualsol[j], exp_dualsol[j]);
      }
   }

   BMSfreeMemoryArray(&dualsol);

   return SCIP_OKAY;
}

/** perform basic test for the given LP problem */
static
SCIP_RETCODE performLPTest(
   int                   ncols,              /**< number of columns */
   SCIP_Real*            obj,                /**< objective function values of columns */
   SCIP_Real*            lb,                 /**< lower bounds of columns */
   SCIP_Real*            ub,                 /**< upper bounds of columns */
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs,                /**< left hand sides of rows */
   SCIP_Real*            rhs,                /**< right hand sides of rows */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  row,                /**< row-indices of constraint-matrix entries */
   int*                  col,                /**< column-indices of constraint-matrix entries */
   SCIP_Real*            val,                /**< values of constraint-matrix entries */
   SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
   SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
   SCIP_Real*            exp_dualsol         /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
   )
{
   /* load LP data, but leave SDP block empty */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, ncols, obj, lb, ub, 0, NULL, NULL, 0, NULL, NULL, NULL, NULL, 0, NULL, NULL, NULL, NULL, NULL,
         nrows, lhs, rhs, nnonz, row, col, val) );

   cr_assert( ! SCIPsdpiWasSolved(sdpi) );

   /* solve problem */
   SCIP_CALL( solveTest(ncols, nrows, exp_primalfeas, exp_dualfeas, exp_dualsol) );

   return SCIP_OKAY;
}

/** perform basic test for the given SDP problem */
static
SCIP_RETCODE performSDPTest(
   int                   ncols,              /**< number of columns */
   SCIP_Real*            obj,                /**< objective function values of columns */
   SCIP_Real*            lb,                 /**< lower bounds of columns */
   SCIP_Real*            ub,                 /**< upper bounds of columns */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables in each SDP-block (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int**                 sdpconstrow,        /**< pointer to row-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int**                 sdpconstcol,        /**< pointer to column-indices of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   SCIP_Real**           sdpconstval,        /**< pointer to values of entries of constant matrix for each block (may be NULL if sdpconstnnonz = 0) */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrices */
   int**                 sdpnblockvarnonz,   /**< sdpnblockvarnonz[i][j] gives the number of nonzeros for the j-th variable (not necessarly
                                              *   variable j) in the i-th block, this is also the length of row/col/val[i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the global index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable in this block, so row[i][j][k] gives
                                              *   the k-th nonzero of the j-th variable (not necessarly variable j) in the i-th block
                                              *   (may be NULL if sdpnnonz = 0)*/
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable in this block (may be NULL if sdpnnonz = 0)*/
   SCIP_Real***          sdpval,             /**< pointer to the values of the nonzeros for each block and variable in this
                                              *   block (may be NULL if sdpnnonz = 0)*/
   int                   nrows,              /**< number of rows */
   SCIP_Real*            lhs,                /**< left hand sides of rows */
   SCIP_Real*            rhs,                /**< right hand sides of rows */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   int*                  row,                /**< row-indices of constraint-matrix entries */
   int*                  col,                /**< column-indices of constraint-matrix entries */
   SCIP_Real*            val,                /**< values of constraint-matrix entries */
   SCIPFEASSTATUS        exp_primalfeas,     /**< expected primal feasibility status */
   SCIPFEASSTATUS        exp_dualfeas,       /**< expected primal feasibility status */
   SCIP_Real*            exp_dualsol         /**< expected dual optimal solution or dual ray if dual is unbounded or NULL */
   )
{
   /* load LP data, but leave SDP block empty */
   SCIP_CALL( SCIPsdpiLoadSDP(sdpi, ncols, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz, sdpconstnblocknonz, sdpconstrow, sdpconstcol, sdpconstval,
         sdpnnonz, sdpnblockvarnonz, sdpvar, sdprow, sdpcol, sdpval, nrows, lhs, rhs, nnonz, row, col, val) );

   cr_assert( ! SCIPsdpiWasSolved(sdpi) );

   /* solve problem */
   SCIP_CALL( solveTest(ncols, nrows, exp_primalfeas, exp_dualfeas, exp_dualsol) );

   return SCIP_OKAY;
}

/** check whether data in LP solver aggrees with original data */
static
SCIP_RETCODE checkData(
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   int                   nnonz               /**< number of nonzero elements in the constraint matrix */
   )
{
   SCIP_Real* lplb;
   SCIP_Real* lpub;
   SCIP_Real* lpobj;
   SCIP_Real* lplhs;
   SCIP_Real* lprhs;
   int lpncols;
   int lpnrows;
   int lpnonz;
   int i;
   int j;

   /* check number of rows and columns */
   SCIP_CALL( SCIPsdpiGetNLPRows(sdpi, &lpnrows) );
   SCIP_CALL( SCIPsdpiGetNVars(sdpi, &lpncols) );
   cr_assert( lpnrows == nrows );
   cr_assert( lpncols == ncols );

   /* get number of nonzeros in matrix */
   SCIP_CALL( SCIPsdpiGetLPNNonz(sdpi, &lpnonz) );
   cr_assert( lpnonz == nnonz );

   /* allocate storage for data */
   BMSallocMemoryArray(&lplb, ncols);
   BMSallocMemoryArray(&lpub, ncols);
   BMSallocMemoryArray(&lpobj, ncols);

   /* get LP data */
   SCIP_CALL( SCIPsdpiGetBounds(sdpi, 0, ncols-1, lplb, lpub) );
   SCIP_CALL( SCIPsdpiGetObj(sdpi, 0, ncols-1, lpobj) );

   /* compare data */
   for (j = 0; j < ncols; ++j)
   {
      cr_assert_float_eq(lplb[j], lb[j], EPS, "Violation of lower bound %d: %g != %g\n", j, lplb[j], lb[j]);
      cr_assert_float_eq(lpub[j], ub[j], EPS, "Violation of upper bound %d: %g != %g\n", j, lpub[j], ub[j]);

      cr_assert_float_eq(lpobj[j], obj[j], EPS, "Violation of objective coefficient %d: %g != %g\n", j, lpobj[j], obj[j]);
   }

   BMSfreeMemoryArray(&lpobj);
   BMSfreeMemoryArray(&lpub);
   BMSfreeMemoryArray(&lplb);

   /* compare lhs/rhs */
   BMSallocMemoryArray(&lplhs, nrows);
   BMSallocMemoryArray(&lprhs, nrows);

   SCIP_CALL( SCIPsdpiGetLhSides(sdpi, 0, nrows - 1, lplhs) );
   SCIP_CALL( SCIPsdpiGetRhSides(sdpi, 0, nrows - 1, lprhs) );

   for (i = 0; i < nrows; ++i)
   {
      cr_assert_float_eq(lplhs[i], lhs[i], EPS, "Violation of lhs %d: %g != %g\n", i, lplhs[i], lhs[i]);
      cr_assert_float_eq(lplhs[i], lhs[i], EPS, "Violation of rhs %d: %g != %g\n", i, lprhs[i], rhs[i]);
   }

   BMSfreeMemoryArray(&lprhs);
   BMSfreeMemoryArray(&lplhs);

   return SCIP_OKAY;
}


/** TESTS **/

/** Test 1
 *
 * min -3 x1 -  x2
 *      2 x1 +   x2 <= 10
 *        x1 + 3 x2 <= 15
 *        x1,    x2 >= 0
 *
 * with optimal solution (5, 0).
 */
Test(checksdpi, test1)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-3, -1};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real rhs[2] = {10, 15};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real ub[2];
   SCIP_Real lhs[2];

   /* expected solutions */
   SCIP_Real exp_dualsol[2] = {5, 0};

   /* fill variable data */
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPfeas, SCIPfeas, exp_dualsol) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 2
 *
 * min -3 x1 -   x2
 *      2 x1 +   x2 <= 10
 *        x1 + 3 x2 <= 15
 *        x1, x2 free
 *
 * which is unbounded.
 */
Test(checksdpi, test2)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-3, -1};
   SCIP_Real rhs[2] = {10, 15};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real lb[2];
   SCIP_Real ub[2];
   SCIP_Real lhs[2];

   /* fill variable data */
   lb[0] = -SCIPsdpiInfinity(sdpi);
   lb[1] = -SCIPsdpiInfinity(sdpi);
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPinfeas, SCIPunbounded, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 3
 *
 * min 10 y1 + 15 y2
 *      2 y1 +   y2 == 3
 *        y1 + 3 y2 == 1
 *        y1,    y2 >= 0
 *
 * which is dual unbounded (this is the dual of the problem in Test 2).
 *
 * Then use rhs (1, 1) with primal optimal solution (0.4,0.2), dual optimal solution (3, 4), activity (0, 0), and redcost (0, 0).
 */
Test(checksdpi, test3)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {10, 15};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real lhs[2] = {3, 1};
   SCIP_Real rhs[2] = {3, 1};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real ub[2];

   /* fill variable data */
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPunbounded,  SCIPinfeas, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 4
 *
 * min -x1 - x2
 *      x1 - x2 <= 0
 *    - x1 + x2 <= -1
 *      x1,  x2 free
 *
 * which primal and dual infeasible.
 */
Test(checksdpi, test4)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-1, -1};
   SCIP_Real rhs[2] = {0, -1};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {1, -1, -1, 1};

   /* data to be filled */
   SCIP_Real lhs[2];
   SCIP_Real lb[2];
   SCIP_Real ub[2];

   /* fill variable data */
   lb[0] = -SCIPsdpiInfinity(sdpi);
   lb[1] = -SCIPsdpiInfinity(sdpi);
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPinfeas,  SCIPinfeas, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 5
 *
 * min -3 x1 -  x2
 *      2 x1 +   x2 <= 10
 *        x1 + 3 x2 <= 15
 *      0 <= x1 <= 0
 *      0 <= x2 <= 0
 *
 * with fixed variables, which is feasible.
 */
Test(checksdpi, test5)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-3, -1};
   SCIP_Real lb[2] = {0, 0};
   SCIP_Real ub[2] = {0, 0};
   SCIP_Real rhs[2] = {10, 15};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real lhs[2];

   /* expected solutions */
   SCIP_Real exp_dualsol[2] = {0, 0};

   /* fill data */
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPfeas, SCIPfeas, exp_dualsol) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 6
 *
 * min -3 x1 -  x2
 *      2 x1 +   x2 <= 10
 *        x1 + 3 x2 <= 15
 *      4 <= x1 <= 4
 *      3 <= x2 <= 3
 *
 * with fixed variables, which is infeasible.
 */
Test(checksdpi, test6)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-3, -1};
   SCIP_Real lb[2] = {4, 3};
   SCIP_Real ub[2] = {4, 3};
   SCIP_Real rhs[2] = {10, 15};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real lhs[2];

   /* fill data */
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPunbounded, SCIPinfeas, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 7
 *
 * min -3 x1 -  x2
 *      2 x1 +   x2 <= 10
 *        x1 + 3 x2 <= 15
 *      1 <= x1 <= 0
 *      0 <= x2 <= 10
 *
 * with conflicting bounds (infeasible)
 */
Test(checksdpi, test7)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-3, -1};
   SCIP_Real lb[2] = {1, 0};
   SCIP_Real ub[2] = {0, 10};
   SCIP_Real rhs[2] = {10, 15};
   int row[4] = {0, 0, 1, 1};
   int col[4] = {0, 1, 0, 1};
   SCIP_Real val[4] = {2, 1, 1, 3};

   /* data to be filled */
   SCIP_Real lhs[2];

   /* fill data */
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   lhs[1] = -SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performLPTest(2, obj, lb, ub, 2, lhs, rhs, 4, row, col, val, SCIPunbounded, SCIPinfeas, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 2, lhs, rhs, 4) );
}

/** Test 8
 *
 *  The following example is motived by an example of Fridberg:
 *  inf   x3
 *        x1 + x2 <= 0        [y1]
 *        x3      >= -1       [y2]
 *        [x1, x2, x3]
 *        [x2, x1,  0]  >= 0  [X]
 *        [x3,  0, x1]
 *
 *  An optimal solution is x = (t, -t, 0) with value 0: Indeed the SDP constraint implies:
 *  - x1 >= 0
 *  - x1^3 - x1 x3^2 - x1 x2^2 >= 0    <=>  x1 (x1^2 - x3^2 - x2^2) >= 0
 *  - x1^2 - x2^2 >= 0
 *  - x1^2 - x3^2 >= 0
 *  Write x2 = -x1 - s for t >= 0. Then x1^2 >= (-x1 - s)^2, which is only possible for s = 0.
 *  Thus, x2 = -x1. Moreover, x1^2 >= x2^2 + x3^2 = x1^2 + x3^2, which implies x3 = 0.
 *  Note, however, that there is a whole line that is optimal.
 *
 *  The dual is (see the brackets above for the dual variables):
 *  sup  -y2
 *       -y1 + X_11 + X_22 + X_33 == 0
 *       -y1 + X_21 == 0
 *        y2 + X_31 == 1
 *       X psd, y1, y2 >= 0.
 *  This problem is feasible, since y2 = 1, y2 = 0, X = 0 is feasible.
 *
 *  Note: We expect this test to crash!
 */
Test(checksdpi, test8)
{
   /* data with fixed values: */
   SCIP_Real obj[3] = {0, 0, 1};
   int row[3] = {0, 0, 1};
   int col[3] = {0, 1, 2};
   SCIP_Real val[3] = {1, 1, 1};
   int nsdpblocks = 1;
   int sdpblocksizes[1] = {3};
   int sdpnblockvars[1] = {3};
   int sdpconstnblocknonz[1] = {0};
   int sdpnnonz = 5;
   int* sdpnblockvarnonz;
   int* sdpvar;
   int** sdprow;
   int** sdpcol;
   SCIP_Real** sdpval;

   int sdpnblockvarnonzs[3] = {3, 1, 1};
   int sdpvars[3] = {0, 1, 2};
   int sdprowss[5] = {0, 1, 2, 1, 2};
   int sdpcolss[5] = {0, 1, 2, 0, 0};
   SCIP_Real sdpvalss[5] = {1.0, 1.0, 1.0, 1.0, 1.0};

   /* data to be filled */
   int* sdprows[3];
   int* sdpcols[3];
   SCIP_Real* sdpvals[3];
   SCIP_Real lhs[2];
   SCIP_Real rhs[2];
   SCIP_Real lb[3];
   SCIP_Real ub[3];

   /* expected solutions */
   SCIP_Real exp_dualsol[3] = {0, 0, 0};

   sdpnblockvarnonz = &sdpnblockvarnonzs[0];
   sdpvar = &sdpvars[0];
   sdprow = &sdprows[0];
   sdpcol = &sdpcols[0];
   sdpval = &sdpvals[0];
   sdprows[0] = &sdprowss[0];
   sdprows[1] = &sdprowss[3];
   sdprows[2] = &sdprowss[4];
   sdpcols[0] = &sdpcolss[0];
   sdpcols[1] = &sdpcolss[3];
   sdpcols[2] = &sdpcolss[4];
   sdpvals[0] = &sdpvalss[0];
   sdpvals[1] = &sdpvalss[3];
   sdpvals[2] = &sdpvalss[4];

   /* fill data */
   lhs[0] = -SCIPsdpiInfinity(sdpi);
   rhs[0] = 0.0;
   lhs[1] = -1.0;
   rhs[1] = SCIPsdpiInfinity(sdpi);
   lb[0] = -SCIPsdpiInfinity(sdpi);
   lb[1] = -SCIPsdpiInfinity(sdpi);
   lb[2] = -SCIPsdpiInfinity(sdpi);
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);
   ub[2] = SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performSDPTest(3, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, 0, sdpconstnblocknonz,
         NULL, NULL, NULL, sdpnnonz, &sdpnblockvarnonz, &sdpvar, &sdprow, &sdpcol, &sdpval,
         2, lhs, rhs, 3, row, col, val, SCIPfeas, SCIPfeas, exp_dualsol) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(3, obj, lb, ub, 2, lhs, rhs, 5) );
}

/** Test 9
 *
 *  inf   -x1
 *        x1 >= -1 [y1]
 *        x1 <= 1  [y2]
 *        x2 >= -1 [y3]
 *        x2 <= 1  [y4]
 *        [x1,      1]
 *        [1, 0.75*x2]  >= 0  [X]
 *
 *  This problem is infeasible, since the SDP constraint implies
 *  x1 * x2 >= 4/3, but x1 * x2 <= 1 due to the variable bounds.
 *
 *  The dual is (see the brackets above for the dual variables):
 *  sup  -2*X21 - y2 - y4 - y1 - y3
 *       y1 - y2 + X_11 == -1
 *       y3 - y4 + 0.75*X_22 == 0
 *       X psd, y1, y2, y3, y4 >= 0.
 *
 *  This problem is feasible, since X = 0, y1 = 0, y2 = 1, y3 = 0, y4 = 0
 *  is feasible. Moreover, choosing
 *  X = [t,     -t] >= 0
 *      [-t, t + 1]
 *  y1 = y3 = 0, y2 = t + 1, y4 = 0.75*(t + 1)
 *  is also feasible and yields a solution value
 *  0.25*t - 1.75
 *  which shows that the problem is in fact unbounded.
 */
Test(checksdpi, test9)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-1, 0};
   int row[4] = {0, 1, 2, 3};
   int col[4] = {0, 0, 1, 1};
   SCIP_Real val[4] = {1, -1, 1, -1};
   int nsdpblocks = 1;
   int sdpblocksizes[1] = {2};
   int sdpnblockvars[1] = {2};
   int sdpconstnblocknonz[1] = {1};
   int sdpconstnnonz = 1;
   int sdpnnonz = 2;
   int* sdpnblockvarnonz;
   int* sdpvar;
   int** sdprow;
   int** sdpcol;
   SCIP_Real** sdpval;
   int* sdpconstrow;
   int* sdpconstcol;
   SCIP_Real* sdpconstval;

   int sdpnblockvarnonzs[2] = {1, 1};
   int sdpvars[2] = {0, 1};
   int sdprowss[2] = {0, 1};
   int sdpcolss[2] = {0, 1};
   int sdpconstrowss[1] = {1};
   int sdpconstcolss[1] = {0};
   SCIP_Real sdpconstvalss[1] = {-1.0};
   SCIP_Real sdpvalss[2] = {1.0, 0.75};

   /* data to be filled */
   int* sdprows[2];
   int* sdpcols[2];
   SCIP_Real* sdpvals[2];
   SCIP_Real lhs[4];
   SCIP_Real rhs[4];
   SCIP_Real lb[2];
   SCIP_Real ub[2];

   sdpnblockvarnonz = &sdpnblockvarnonzs[0];
   sdpvar = &sdpvars[0];
   sdprow = &sdprows[0];
   sdpcol = &sdpcols[0];
   sdpval = &sdpvals[0];
   sdprows[0] = &sdprowss[0];
   sdprows[1] = &sdprowss[1];
   sdpcols[0] = &sdpcolss[0];
   sdpcols[1] = &sdpcolss[1];
   sdpvals[0] = &sdpvalss[0];
   sdpvals[1] = &sdpvalss[1];

   sdpconstrow = &sdpconstrowss[0];
   sdpconstcol = &sdpconstcolss[0];
   sdpconstval = &sdpconstvalss[0];

   /* fill data */
   lhs[0] = -1.0;
   lhs[1] = -1.0;
   lhs[2] = -1.0;
   lhs[3] = -1.0;
   rhs[0] = SCIPsdpiInfinity(sdpi);
   rhs[1] = SCIPsdpiInfinity(sdpi);
   rhs[2] = SCIPsdpiInfinity(sdpi);
   rhs[3] = SCIPsdpiInfinity(sdpi);
   lb[0] = -SCIPsdpiInfinity(sdpi);
   lb[1] = -SCIPsdpiInfinity(sdpi);
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performSDPTest(2, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz, sdpconstnblocknonz,
         &sdpconstrow, &sdpconstcol, &sdpconstval, sdpnnonz, &sdpnblockvarnonz, &sdpvar, &sdprow, &sdpcol, &sdpval,
         4, lhs, rhs, 4, row, col, val, SCIPunbounded, SCIPinfeas, NULL) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 4, lhs, rhs, 4) );
}


/** Test 10
 *
 *  inf   -x1 - x2
 *        x1 >= -1 [y1]
 *        x1 <= 1  [y2]
 *        x2 >= -1 [y3]
 *        x2 <= 1  [y4]
 *        [x1,  0]
 *        [0,  x2]  >= 0  [X]
 *
 *  This problem is feasible with optimal solution (x1,x2) = (1,1) and
 *  optimal value -2.
 *
 *  The dual is (see the brackets above for the dual variables):
 *  sup  -y2 - y4 - y1 - y3
 *       y1 - y2 + X_11 == -1
 *       y3 - y4 + X_22 == -1
 *       X psd, y1, y2, y3, y4 >= 0.
 *
 *  This problem is also feasible with optimal solution X = 0, y1 = y3 = 0, y2 = y4 = 1,
 *  and optimal value -2.
 */
Test(checksdpi, test10)
{
   /* data with fixed values: */
   SCIP_Real obj[2] = {-1.0, -1.0};
   int row[4] = {0, 1, 2, 3};
   int col[4] = {0, 0, 1, 1};
   SCIP_Real val[4] = {1, -1, 1, -1};
   int nsdpblocks = 1;
   int sdpblocksizes[1] = {2};
   int sdpnblockvars[1] = {2};
   int sdpconstnblocknonz[1] = {0};
   int sdpconstnnonz = 0;
   int sdpnnonz = 2;
   int* sdpnblockvarnonz;
   int* sdpvar;
   int** sdprow;
   int** sdpcol;
   SCIP_Real** sdpval;

   int sdpnblockvarnonzs[2] = {1, 1};
   int sdpvars[2] = {0, 1};
   int sdprowss[2] = {0, 1};
   int sdpcolss[2] = {0, 1};
   SCIP_Real sdpvalss[2] = {1.0, 1.0};

   /* data to be filled */
   int* sdprows[2];
   int* sdpcols[2];
   SCIP_Real* sdpvals[2];
   SCIP_Real lhs[4];
   SCIP_Real rhs[4];
   SCIP_Real lb[2];
   SCIP_Real ub[2];

   /* expected solutions */
   SCIP_Real exp_dualsol[2] = {1.0, 1.0};

   sdpnblockvarnonz = &sdpnblockvarnonzs[0];
   sdpvar = &sdpvars[0];
   sdprow = &sdprows[0];
   sdpcol = &sdpcols[0];
   sdpval = &sdpvals[0];
   sdprows[0] = &sdprowss[0];
   sdprows[1] = &sdprowss[1];
   sdpcols[0] = &sdpcolss[0];
   sdpcols[1] = &sdpcolss[1];
   sdpvals[0] = &sdpvalss[0];
   sdpvals[1] = &sdpvalss[1];

   /* fill data */
   lhs[0] = -1.0;
   lhs[1] = -1.0;
   lhs[2] = -1.0;
   lhs[3] = -1.0;
   rhs[0] = SCIPsdpiInfinity(sdpi);
   rhs[1] = SCIPsdpiInfinity(sdpi);
   rhs[2] = SCIPsdpiInfinity(sdpi);
   rhs[3] = SCIPsdpiInfinity(sdpi);
   lb[0] = -SCIPsdpiInfinity(sdpi);
   lb[1] = -SCIPsdpiInfinity(sdpi);
   ub[0] = SCIPsdpiInfinity(sdpi);
   ub[1] = SCIPsdpiInfinity(sdpi);

   SCIP_CALL( performSDPTest(2, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpnblockvars, sdpconstnnonz, sdpconstnblocknonz,
         NULL, NULL, NULL, sdpnnonz, &sdpnblockvarnonz, &sdpvar, &sdprow, &sdpcol, &sdpval,
         4, lhs, rhs, 4, row, col, val, SCIPfeas, SCIPfeas, exp_dualsol) );

   /* check that data stored in sdpi is still the same */
   SCIP_CALL( checkData(2, obj, lb, ub, 4, lhs, rhs, 4) );
}
