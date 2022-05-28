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

/**@file   sdpisolver.h
 * @brief  interface methods for specific SDP-solvers
 * @author Marc Pfetsch
 * @author Tristan Gally
 *
 * This file specifies a generic SDP-solver interface used by SCIP to create, modify, and solve semidefinite programs of
 * the (dual) form
 *
 *   \f{eqnarray*}{
 *   	  \min & & b^T y \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u
 *   \f}
 * for symmetric matrices \f$ A_j^i \in S_{k_i} \f$, a matrix \f$ D \in \mathbb{R}^{k_0 \times n} \f$ and query information about the solution.
 *
 * All indexing (rows, columns, blocks and variables) starts at 0.
 *
 * Although it includes a few SCIP header files, e.g., because it uses SCIP's return codes, it can be used independently of
 * any SCIP instance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPISOLVER_H__
#define __SCIP_SDPISOLVER_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "sdpi/type_sdpi.h"
#include "sdpi/sdpiclock.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_SDPiSolver SCIP_SDPISOLVER;                 /**< solver dependent SDP interface */

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version (if available) of SDP-solver */
SCIP_EXPORT
const char* SCIPsdpiSolverGetSolverName(
   void
   );

/** gets description of SDP-solver (developer, webpage, ...) */
SCIP_EXPORT
const char* SCIPsdpiSolverGetSolverDesc(
   void
   );

/** gets pointer to SDP-solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP-solver object.
 */
SCIP_EXPORT
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** gets default number of increases of penalty parameter for SDP-solver in SCIP-SDP */
SCIP_EXPORT
int SCIPsdpiSolverGetDefaultSdpiSolverNpenaltyIncreases(
   void
   );

/** Should primal solution values be saved for warmstarting purposes? */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverDoesWarmstartNeedPrimal(
   void
   );

/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP solver interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< SDP-solver interface */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem              /**< buffer memory */
   );

/** deletes an SDP solver interface */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< SDP-solver interface */
   );

/** increases the SDP-Counter */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverIncreaseCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** reset the SDP-Counter to zero */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverResetCounter(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** loads and solves an SDP
 *
 *  For the non-constant SDP- and the LP-part, the original arrays before fixings should be given, for the constant
 *  SDP-part the arrays AFTER fixings should be given. In addition, an array needs to be given, that for every block and
 *  every row/col index within that block either has value -1, meaning that this index should be deleted, or a
 *  non-negative integer stating the number of indices before it that are to be deleated, meaning that this index will
 *  be decreased by that number, in addition to that the total number of deleted indices for each block should be given.
 *  Optionally an array start may be given with a starting point for the solver (if this is NULL then the solver should
 *  start from scratch).
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 *  @note starting point needs to be given with original indices (before any local presolving), last block should be the LP block with indices
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independant of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
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
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP) */
   SDPI_CLOCK*           usedsdpitime        /**< clock to measure how much time has been used for the current solve */
);

/** loads and solves an SDP using a penalty formulation
 *
 *  The penalty formulation of the SDP is:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy + r \cdot \mathbb{I} \geq d \\
 *      & & l \leq y \leq u \\
 *      & & r \geq 0.\f}
 *  Alternatively withobj can be set to false to set b to 0 and only check for feasibility (if the optimal objective value is
 *  bigger than 0 the problem is infeasible, otherwise it's feasible), and rbound can be set to false to remove the non-negativity condition on r.
 *  For the non-constant SDP- and the LP-part the original arrays before fixings should be given, for the constant SDP-part the arrays AFTER fixings
 *  should be given. In addition, an array needs to be given, that for every block and every row/col index within that block either has value
 *  -1, meaning that this index should be deleted, or a non-negative integer stating the number of indices before it that are to be deleated,
 *  meaning that this index will be decreased by that number. Moreover, the total number of deleted indices for each block should be given.
 *  An optional starting point for the solver may be given; if it is NULL, the solver will start from scratch.
 *
 *  @warning Depending on the solver, the given lp arrays might get sorted in their original position.
 *  @note starting point needs to be given with original indices (before any local presolving), last block should be the LP block with indices
 *  lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ... independant of some lhs/rhs being infinity (the starting point
 *  will later be adjusted accordingly)
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the Gamma above, needs to be >= 0 */
   SCIP_Bool             withobj,            /**< if this is false the objective is set to 0 */
   SCIP_Bool             rbound,             /**< should r be non-negative ? */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            obj,                /**< objective coefficients of variables */
   SCIP_Real*            lb,                 /**< lower bounds of variables */
   SCIP_Real*            ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int*                  sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   int*                  sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int**                 sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS */
   int**                 sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real**           sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int**                 sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int**                 sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int***                sdprow,             /**< pointer to the row-indices for each block and variable */
   int***                sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real***          sdpval,             /**< values of SDP-constraintmmatrix entries (may be NULL if sdpnnonz = 0) */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   SCIP_Real*            lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   SCIP_Real*            lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   int*                  lprow,              /**< row-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   int*                  lpcol,              /**< column-index for each entry in lpval-array, might get sorted (may be NULL if lpnnonz = 0) */
   SCIP_Real*            lpval,              /**< values of LP-constraint-matrix entries, might get sorted (may be NULL if lpnnonz = 0) */
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
   SCIP_SDPSOLVERSETTING startsettings,      /**< settings used to start with in SDPA, currently not used for DSDP, set this to
                                              *   SCIP_SDPSOLVERSETTING_UNSOLVED to ignore it and start from scratch */
   SCIP_Real             timelimit,          /**< after this many seconds solving will be aborted (currently only implemented for DSDP) */
   SDPI_CLOCK*           usedsdpitime,       /**< clock to measure how much time has been used for the current solve */
   SCIP_Bool*            feasorig,           /**< pointer to store if the solution to the penalty-formulation is feasible for the original problem
                                              *   (may be NULL if penaltyparam = 0) */
   SCIP_Bool*            penaltybound        /**< pointer to store if the primal solution reached the bound Tr(X) <= penaltyparam in the primal problem,
                                              *   this is also an indication of the penalty parameter being to small (may be NULL if not needed) */
);



/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns true if the solver could determine whether the problem is feasible
 *
 *  So it returns true if the solver knows that the problem is feasible/infeasible/unbounded, it returns false if the
 *  solver does not know anything about the feasibility status and thus the functions IsPrimalFeasible etc. should not be
 *  used.
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff SDP is proven to be primal unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP is proven to be primal infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP is proven to be primal feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP is proven to be dual unbounded,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP is proven to be dual infeasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP is proven to be dual feasible,
 *  returns FALSE with a debug-message if the solver could not determine feasibility
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff the solver converged
 *
 * @todo Define this more precisely and adjust SDP-solvers (currently Mosek returns true for infeasibility and SDPA false).
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff the objective limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff the iteration limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff the time limit was reached */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns the internal solution status of the solver, which has the following meaning:<br>
 * -1: solver was not started<br>
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
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** returns TRUE iff SDP was solved to optimality or some other status was reached
 *  that is still acceptable inside a Branch & Bound framework
 */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** tries to reset the internal status of the SDP-solver in order to ignore an instability of the last solving call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real*            objval              /**< pointer to store the objective value */
   );

/** gets dual solution vector for feasible SDPs
 *
 *  If dualsollength isn't equal to the number of variables this will return the needed length and a debug message is thrown.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            objval,             /**< pointer to store the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   );

/** return number of nonzeros for each block of the primal solution matrix X for the preoptimal solution */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPreoptimalPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block
                                              *   or first entry -1 if no primal solution is available */
   );

/** gets preoptimal dual solution vector and primal matrix for warmstarting purposes
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independant of some lhs/rhs being infinity
 *  @note If dualsollength isn't equal to the number of variables this will return the needed length and a debug message is thrown.
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPreoptimalSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Bool*            success,            /**< could a preoptimal solution be returned ? */
   SCIP_Real*            dualsol,            /**< pointer to store the dual solution vector, may be NULL if not needed */
   int*                  dualsollength,      /**< length of the dual sol vector, must be 0 if dualsol is NULL, if this is less than the number
                                              *   of variables in the SDP, a DebugMessage will be thrown and this is set to the needed value */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) or -1 if no primal matrix should be returned */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block (or NULL if nblocks = -1)
                                              *   output: number of nonzeros in each block or first entry -1 if no primal solution is available */
   int**                 startXrow,          /**< pointer to store row indices of X (or NULL if nblocks = -1) */
   int**                 startXcol,          /**< pointer to store column indices of X (or NULL if nblocks = -1) */
   SCIP_Real**           startXval           /**< pointer to store values of X (or NULL if nblocks = -1) */
   );

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note If a variable is either fixed or unbounded in the dual problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPrimalBoundVars(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars <br>
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   );

/** gets the primal solution corresponding to the LP row sides */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPrimalLPSides(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int                   nlpcons,            /**< number of LP rows */
   SCIP_Real*            lplhs,              /**< lhs of LP rows */
   SCIP_Real*            lprhs,              /**< rhs of LP rows */
   SCIP_Real*            lhsvals,            /**< array to store the values of the variables corresponding to LP lhs */
   SCIP_Real*            rhsvals             /**< array to store the values of the variables corresponding to LP rhs */
   );

/** return number of nonzeros for each block of the primal solution matrix X (including lp block) */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPrimalNonzeros(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz    /**< pointer to store number of nonzeros for row/col/val-arrays in each block */
   );

/** returns the primal matrix X
 *
 *  @note last block will be the LP block (if one exists) with indices lhs(row0), rhs(row0), lhs(row1), ..., lb(var1), ub(var1), lb(var2), ...
 *  independant of some lhs/rhs being infinity
 *  @note If the allocated memory for row/col/val is insufficient, a debug message will be thrown and the neccessary amount is returned in startXnblocknonz
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPrimalMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nblocks,            /**< length of startXnblocknonz (should be nsdpblocks + 1) */
   int*                  startXnblocknonz,   /**< input: allocated memory for row/col/val-arrays in each block
                                              *   output: number of nonzeros in each block */
   int**                 startXrow,          /**< pointer to store row indices of X */
   int**                 startXcol,          /**< pointer to store column indices of X */
   SCIP_Real**           startXval           /**< pointer to store values of X */
   );

/** returns the primal solution matrix (without LP rows) */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetPrimalSolutionMatrix(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   int                   nsdpblocks,         /**< number of blocks */
   int*                  sdpblocksizes,      /**< sizes of the blocks */
   int**                 indchanges,         /**< changes needed to be done to the indices, if indchanges[block][nonz]=-1, then
                                              *   the index can be removed, otherwise it gives the number of indices removed before this */
   int*                  nremovedinds,       /**< pointer to store the number of rows/cols to be fixed for each block */
   int*                  blockindchanges,    /**< pointer to store index change for each block, system is the same as for indchanges */
   SCIP_Real**           primalmatrices      /**< pointer to store values of the primal matrix */
   );

/** return the maximum absolute value of the optimal primal matrix */
SCIP_EXPORT
SCIP_Real SCIPsdpiSolverGetMaxPrimalEntry(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to an SDP-solver interface */
   );

/** gets the time for the last SDP optimization call of solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetTime(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real*            opttime             /**< pointer to store the time for optimization of the solver */
   );

/** gets the number of SDP iterations of the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/** gets the number of calls to the SDP-solver for the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetSdpCalls(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   int*                  calls               /**< pointer to store the number of calls to the SDP-solver for the last solve call */
   );

/** gets the settings used by the SDP solver for the last solve call */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverSettingsUsed(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPSOLVERSETTING* usedsetting        /**< the setting used by the SDP-solver */
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
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver          /**< SDP-solver interface */
   );

/** checks if given value is treated as (plus or minus) infinity in the SDP-solver */
SCIP_EXPORT
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_Real             val                 /**< value to be checked for infinity */
   );

/** gets floating point parameter of SDP-Solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of SDP-Solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverSetRealpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/** gets integer parameter of SDP-Solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< pointer to store the parameter value */
   );

/** sets integer parameter of SDP-Solver */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** compute and set lambdastar (only used for SDPA) */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverComputeLambdastar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxguess            /**< maximum guess for lambda star of all SDP-constraints */
   );

/** compute and set the penalty parameter */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverComputePenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             maxcoeff,           /**< maximum objective coefficient */
   SCIP_Real*            penaltyparam        /**< the computed penalty parameter */
   );

/** compute and set the maximum penalty parameter */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverComputeMaxPenaltyparam(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP-solver interface */
   SCIP_Real             penaltyparam,       /**< the initial penalty parameter */
   SCIP_Real*            maxpenaltyparam     /**< the computed maximum penalty parameter */
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   const char*           fname               /**< file name */
   );

/** writes SDP to a file */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< SDP-solver interface */
   const char*           fname               /**< file name */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
