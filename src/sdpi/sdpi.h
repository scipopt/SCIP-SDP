/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpi.h
 * @brief  interface methods for specific SDP solvers
 * @author Marc Pfetsch
 * @author Leif Naundorf
 * @author Tristan Gally
 *
 * This file specifies a generic SDP solver interface used by SCIP to create, modify, and solve semidefinite programs of
 * the (dual) form
 *
 *   \f{eqnarray*}{
 *   	  \min & & b^T y \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u
 *   \f}
 *
 * for symmetric matrices \f A_j^i \in S_{k_i} \f and a matrix \f D \in \mathbb{R}^{k_0 \times n} \f and query information about the solution.
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

#include "type_sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_SDPiSolver SCIP_SDPISOLVER;                 /**< solver dependent SDP interface */

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of SDP solver */
EXTERN
const char* SCIPsdpiSolverGetSolverName(
   void
   );

/** gets description of SDP solver (developer, webpage, ...) */
EXTERN
const char* SCIPsdpiSolverGetSolverDesc(
   void
   );

/** Does the solver have a way to solve a penalty formulation on its own or must one be provided */
EXTERN
SCIP_Bool SCIPsdpiSolverKnowsPenalty(
   void
   );

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
EXTERN
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,          /**< SDP interface solver structure */
   SCIP_MESSAGEHDLR*     messagehdlr,         /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem              /**< block memory */
   );

/** deletes an SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver,          /**< SDP interface solver structure */
   );

/**@} */

/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** loads an SDP and immediatly solves it */
EXTERN
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                *  number of entries  of sdpconst row/col/val [i] */
   const int**          sdpconstrow,        /**< pointers to row-indices for each block */
   const int**          sdpconstcol,        /**< pointers to column-indices for each block */
   const SCIP_Real**     sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpnblockvarnonz,   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i * nvars + j] */
   const int**          sdprow,             /**< pointer to the row-indices for each block and variable */
   const int**          sdpcol,             /**< pointer to the column-indices for each block and variable */
   const SCIP_Real**     sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   );

/** loads a SDP solves the following penalty formulation of the SDP:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u}
 *   \f
 *   alternatively withObj can be to false to set \f b \f to false and only check for feasibility (if the optimal
 *   objective value is bigger than 0 the problem is infeasible, otherwise it's feasible)
 *
 *   attention: this only works for some solvers, check with SCIPsdpiKnowsPenalty first, otherwise this returns an error (in that case you should form
 *   the penalty formulation yourself and pass it via LoadAndSolve */
EXTERN
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Real             gamma,              /**< the penalty parameter above, needs to be >= 0 */
   SCIP_Bool             withObj,            /**< if this is false, the objective is set to 0 */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                *  number of entries  of sdpconst row/col/val [i] */
   const int**          sdpconstrow,        /**< pointers to row-indices for each block */
   const int**          sdpconstcol,        /**< pointers to column-indices for each block */
   const SCIP_Real**     sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpnblockvarnonz,   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i * nvars + j] */
   const int**          sdprow,             /**< pointer to the row-indices for each block and variable */
   const int**          sdpcol,             /**< pointer to the column-indices for each block and variable */
   const SCIP_Real**     sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
);



/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
EXTERN
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns true if the solver could determine whether or not the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver doesn't know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be used */
EXTERN
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** gets information about primal and dual feasibility of the current SDP solution */
EXTERN
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverExistsPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverHasPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be primal unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be primal infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be primal feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverExistsDualRay(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverHasDualRay(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be dual unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be dual infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP is proven to be dual feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
EXTERN
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff the solver converged */
EXTERN
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff the objective limit was reached */
EXTERN
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff the iteration limit was reached */
EXTERN
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff the time limit was reached */
EXTERN
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns the internal solution status of the solver, which has the following meaning:
 *  0: converged
 *  1: infeasible start
 *  2: numerical problems
 *  3: objective limit reached
 *  4: iteration limit reached
 *  5: time limit reached
 *  6: user termination
 *  7: other */
EXTERN
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP was solved to optimality, meaning the solver converged and returned primal and dual feasible solutions */
EXTERN
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 * that is still acceptable inside a Branch & Bound framework */
EXTERN
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
EXTERN
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
EXTERN
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Real*            objval              /**< stores the objective value */
   );

/** gets dual solution vector for feasible SDPs, if dualsollength isn't equal to the number of variables this will return an error */
EXTERN
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int                   dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL */
   );

/** gets the number of SDP iterations of the last solve call */
EXTERN
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
EXTERN
SCIP_RETCODE SCIPsdpiSolverGetRealSolQuality(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
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

/** returns value treated as infinity in the SDP solver */
EXTERN
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** checks if given value is treated as infinity in the SDP solver */
EXTERN
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   );

/** returns highest penalty parameter to be used */
EXTERN
SCIP_Real SCIPsdpiSolverMaxPenParam(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   );

/** checks if given value is greater or equal to the highest penalty parameter to be used */
EXTERN
SCIP_Bool SCIPsdpiSolverIsGEMaxPenParam(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   SCIP_Real            val                 /**< value to be compared to maximum penalty parameter */
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
EXTERN
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   const char*           fname               /**< file name */
   );

/** writes SDP to a file */
EXTERN
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   const char*           fname               /**< file name */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
