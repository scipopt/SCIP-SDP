/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpi.h
 * @brief  General interface methods for SDP-preprocessing (mainly fixing variables and removing empty rows/cols)
 * @author Tristan Gally
 *
 * This file specifies a generic SDP-solver interface used by SCIP to create, modify, and solve semidefinite programs of
 * the (dual) form
 * \f{align*}{
 *    \min\quad & b^T y \\
 *    \mbox{s.t.} & \sum_{j \in J} A_j^{(k)} y_j - A_0^{(k)} \succeq 0 & \forall \ k \in K, \\
 *     & \sum_{j \in J} d_{ij}\, y_j \geq c_i & \forall \ i \in I, \\
 *     & \ell_j \leq y_j \leq u_j & \forall \ j \in J,
 * \f}
 * for symmetric matrices \f$ A_i^{(k)} \in S_{n_k} \f$ and a matrix \f$ D \in \mathbb{R}^{I \times J} \f$ and query
 * information about the solution.
 * The code refers to this problem as the @em dual.
 *
 * We consider a problem (primal or dual) to be unbounded if there exists a ray and it is feasible.
 *
 * All indexing (rows, columns, blocks and variables) starts at 0.
 *
 * Although it includes a few SCIP header files, e.g., because it uses SCIP's return codes, it can be used independently of
 * any SCIP instance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPI_H__
#define __SCIP_SDPI_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "sdpi/type_sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and potentially version of SDP-solver */
SCIP_EXPORT
const char* SCIPsdpiGetSolverName(
   void
   );

/** gets description of SDP-solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPsdpiGetSolverDesc(
   void
   );

/** gets pointer for SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
SCIP_EXPORT
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
SCIP_EXPORT
int SCIPsdpiGetDefaultSdpiSolverNpenaltyIncreases(
   void
   );

/** Should primal solution values be saved for warmstarting purposes? */
SCIP_EXPORT
SCIP_Bool SCIPsdpiDoesWarmstartNeedPrimal(
   void
   );

/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an sdpi object */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   );

/** deletes an sdpi object */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP-interface structure */
   );

/** cloning method of the general SDP-Interface
 *
 *  @note The solver specific interface is created anew and not copied.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiClone(
   SCIP_SDPI*            oldsdpi,            /**< pointer to the SDP-interface structure that should be cloned */
   SCIP_SDPI*            newsdpi             /**< pointer to an SDP-interface structure to clone into */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data into SDP-solver
 *
 *  @note As the SDP-constraint-matrices are symmetric, only the lower triangular part of them must be specified.
 * @note It is assumed that the matrices are in lower triangular form.
 *  @note There must be at least one variable, the SDP- and/or LP-part may be empty.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective function values of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
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
                                              *   (may be NULL if sdpnnonz = 0) */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable in this block (may be NULL if sdptnnonz = 0) */
   SCIP_Real***          sdpval,             /**< pointer to the values of the nonzeros for each block and variable in this block (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of LP rows (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   );

/** adds rows to the LP-Block
 *
 *  @note Arrays are not checked for duplicates, problems may appear if indices are added more than once.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left-hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right-hand sides of new rows */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            row,                /**< row-indices of constraint-matrix entries, going from 0 to nrows - 1, these will be changed
                                              *   to nlpcons + i */
   const int*            col,                /**< column-indices of constraint-matrix entries */
   const SCIP_Real*      val                 /**< values of constraint-matrix entries */
   );

/** deletes all rows in the given range from the LP-Block */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes LP-rows from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  dstat               /**< deletion status of LP rows <br>
                                              *   input:  1 if row should be deleted, 0 otherwise <br>
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole SDP */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** changes objective coefficients of variables */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change objective coefficients for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      obj                 /**< new objective coefficients */
   );

/** changes lower and upper bounds of variables */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   );

/** changes left- and right-hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPLhRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      lhs,                /**< new values for left-hand sides */
   const SCIP_Real*      rhs                 /**< new values for right-hand sides */
   );


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** returns the currently installed sdpi message handler, or NULL if messages are currently suppressed */
SCIP_MESSAGEHDLR* SCIPsdpiGetMessagehdlr(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets the number of LP-rows in the SDP */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   );

/** gets the number of SDP-Blocks in the SDP */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of blocks */
   );

/** gets the number of variables in the SDP */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   );

/** gets the number of nonzero elements in the SDP-constraint-matrices */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP-constraint-matrices */
   );

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   );

/** gets the number of nonzero elements in the LP-Matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   );

/** gets SDP data from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetSDPdata(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int**                 sdpblocksizes,      /**< sizes of the SDP-blocks */
   int**                 sdpnblockvars,      /**< number of variables in each SDP-block */
   int***                sdpnblockvarnonz,   /**< sdpnblockvarnonz[i][j] = nonzeros of j-th variable in i-th block (length of row/col/val[i][j]) */
   int****               sdprow,             /**< sdprow[b][v][j] = row of j-th nonzero of variable v in block b */
   int****               sdpcol,             /**< sdprow[b][v][j] = column of j-th nonzero of variable v in block b */
   SCIP_Real****         sdpval,             /**< sdpval[i][j][k] = value of j-th nonzero of variable v in block b */
   int**                 sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] */
   int***                sdpconstrow,        /**< pointers to row-indices for each block */
   int***                sdpconstcol,        /**< pointers to column-indices for each block */
   SCIP_Real***          sdpconstval         /**< pointers to the values of the nonzeros for each block */
   );

/** gets objective coefficients from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< pointer to store objective coefficients (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );

/** gets current variable lower and/or upper bounds from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< pointer to store lower bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   SCIP_Real*            ubs                 /**< pointer to store upper bound values (memory of size lastvar - firstvar + 1 needs to be allocated), or NULL */
   );

/** gets current left-hand sides from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetLhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss                /**< pointer to store left-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );

/** gets current right-hand sides from SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< pointer to store right-hand side values (memory of size lastvar - firstvar + 1 needs to be allocated) */
   );


/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP, as start optionally a starting point for the solver may be given, if it is NULL, the solver will start from scratch
 *
 *  @note starting point needs to be given with original indices (before any local presolving), last block should be the LP block with indices
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independent of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            starty,             /**< NULL or dual vector y as starting point for the solver, this should have length nvars */
   int*                  startZnblocknonz,   /**< dual matrix Z = sum Ai yi as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startZrow,          /**< dual matrix Z = sum Ai yi as starting point for the solver: row indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int**                 startZcol,          /**< dual matrix Z = sum Ai yi as starting point for the solver: column indices for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   SCIP_Real**           startZval,          /**< dual matrix Z = sum Ai yi as starting point for the solver: values for each block;
                                              *   may be NULL if startZnblocknonz = NULL */
   int*                  startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                              *   also length of corresponding row/col/val-arrays; or NULL */
   int**                 startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   int**                 startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_Real**           startXval,          /**< primal matrix X as starting point for the solver: values for each block;
                                              *   may be NULL if startXnblocknonz = NULL */
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP or MOSEK, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Bool             enforceslatercheck, /**< always check for Slater condition in case the problem could not be solved and printf the solution
                                              *   of this check */
   SCIP_Real             timelimit           /**< after this many seconds solving will be aborted (currently only implemented for DSDP and MOSEK) */
   );

/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was successfully called after the last modification of the SDP */
SCIP_EXPORT
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns whether the original problem was solved, if SCIPsdpiWasSolved = true and SCIPsdpiSolvedOrig = false, then a penalty formulation was solved */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolvedOrig(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns true if the solver could determine whether the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver does not know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be used
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiFeasibilityKnown(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets information about proven primal and dual feasibility of the current SDP-solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Bool*            primalfeasible,     /**< pointer to store the proven primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< pointer to store the proven dual feasibility status */
   );

/** returns TRUE iff SDP is proven to be primal unbounded */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be primal infeasible */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be primal feasible */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual unbounded */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual infeasible */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP is proven to be dual feasible */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the solver converged */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsConverged(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the objective limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff the time limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns the internal solution status of the solver, which has the following meaning:<br>
 *  -1: solver was not started<br>
 *  0: converged<br>
 *  1: infeasible start<br>
 *  2: numerical problems<br>
 *  3: objective limit reached<br>
 *  4: iteration limit reached<br>
 *  5: time limit reached<br>
 *  6: user termination<br>
 *  7: other
 */
SCIP_EXPORT
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 * that is still acceptable inside a Branch & Bound framework
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsAcceptable(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** gets objective value of solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   );

/** gets the best lower bound on the objective (this is equal to objval, if the problem was solved successfully, but can also give a bound
 *  if we did not get a feasible solution using the penalty approach)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetLowerObjbound(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objlb               /**< pointer to store the lower bound on the objective value */
   );

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return the needed length and
 *  a debug message
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dualsol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a debug-message will be thrown and this is set to the needed value */
   );

/** return number of nonzeros for each block of the primal solution matrix X for the preoptimal solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPreoptimalPrimalNonzeros(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block
                                              *   or first entry -1 if no primal solution is available */
   );

/** gets preoptimal dual solution vector and primal matrix for warmstarting purposes
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity
 *  @note If dualsollength isn't equal to the number of variables this will return the needed length and a debug message is thrown.
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPreoptimalSol(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Bool*            success,            /**< could a preoptimal solution be returned ? */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength,      /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) or -1 if no primal matrix should be returned */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block (or NULL if nblocks = -1)
                                              *   output: number of nonzeros in each block (or NULL if nblocks = -1) */
   int**                 startXrow,          /**< pointer to store row indices of X or first entry -1 if no primal solution is available */
   int**                 startXcol,          /**< pointer to store column indices of X (or NULL if nblocks = -1) */
   SCIP_Real**           startXval           /**< pointer to store values of X (or NULL if nblocks = -1) */
   );

/** gets the primal solution corresponding to the lower and upper variable-bounds in the primal problem
 *
 *  The arrays should have size nvars.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPrimalBoundVars(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real*            lbvals,             /**< array to store the values of the variables corresponding to lower bounds in the primal problem */
   SCIP_Real*            ubvals,             /**< array to store the values of the variables corresponding to upper bounds in the primal problem */
   SCIP_Bool*            success             /**< pointer to store whether values could be retrieved */
   );

/** gets the primal variables corresponding to the LP sidex
 *
 *  @note If an LP row was removed, we return a value of 0.0. This can happen if the row is redundant, e.g., all
 *  involved variables are fixed, or it contains variable a single variable only.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPrimalLPSides(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real*            lhsvals,            /**< array to store the values of the variables corresponding to LP lhs */
   SCIP_Real*            rhsvals,            /**< array to store the values of the variables corresponding to LP rhs */
   SCIP_Bool*            success             /**< pointer to store whether values could be retrieved */
   );

/** return number of nonzeros for each block of the primal solution matrix X */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPrimalNonzeros(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   );

/** returns the primal matrix X
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independent of some lhs/rhs being infinity
 *
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPrimalMatrix(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block
                                              *   output: number of nonzeros in each block */
   int**                 startXrow,          /**< pointer to store row indices of X */
   int**                 startXcol,          /**< pointer to store column indices of X */
   SCIP_Real**           startXval           /**< pointer to store values of X */
   );

/** returns the primal solution matrix (without LP rows) */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetPrimalSolutionMatrix(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP-interface structure */
   SCIP_Real**           primalmatrices,     /**< pointer to store values of the primal matrix */
   SCIP_Bool*            success             /**< pointer to store whether the call was successfull */
   );

/** return the maximum absolute value of the optimal primal matrix */
SCIP_EXPORT
SCIP_Real SCIPsdpiGetMaxPrimalEntry(
   SCIP_SDPI*            sdpi                /**< pointer to an SDP-interface structure */
   );

/** gets the time for the last SDP optimization call of solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetTime(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real*            opttime             /**< pointer to store the time for optimization of the solver */
   );

/** gets the number of SDP-iterations of the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetSdpCalls(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   );

/** returns which settings the SDP-solver used in the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSettingsUsed(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
   );

/** returns which settings the SDP-solver used in the last solve call and whether primal and dual Slater condition were fullfilled */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSlaterSettings(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATERSETTING* slatersetting      /**< the combination of Slater conditions and successfull settings */
   );

/** returns whether primal and dual Slater condition held for last solved SDP */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSlater(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPSLATER*       primalslater,       /**< pointer to save whether primal Slater condition held */
   SCIP_SDPSLATER*       dualslater          /**< pointer to save whether dual Slater condition held */
   );

/** returns some statistcs */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetStatistics(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int*                  ninfeasible,        /**< pointer to store the total number of times infeasibility was detected in presolving */
   int*                  nallfixed,          /**< pointer to store the total number of times all variables were fixed */
   int*                  nonevarsdp          /**< pointer to store the total number of times a one variable SDP was solved */
   );

/**@} */




/*
 * SDPi State Methods
 */



/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP-solver */
SCIP_EXPORT
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*            sdpi                /**< SDP-interface structure */
   );

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_EXPORT
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   );

/** gets floating point parameter of SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< pointer to store the parameter value */
   );

/** sets floating point parameter of SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets integer parameter of SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );

/** sets integer parameter of SDP-interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** compute and set lambdastar (only used for SDPA) */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiComputeLambdastar(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   );

/** compute and set the penalty parameter */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiComputePenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   );

/** compute and set the maximum penalty parameter */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiComputeMaxPenaltyparam(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   );

/** sets the type of the clock */
SCIP_EXPORT
void SCIPsdpiClockSetType(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   int                   clocktype           /**< type of clock (1 = CPU, 2 = Wall) */
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   );

/** writes SDP to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP-interface structure */
   const char*           fname               /**< file name */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
