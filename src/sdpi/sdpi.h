/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of SDP solver */
EXTERN
const char* SCIPsdpiGetSolverName(
   void
   );

/** gets description of SDP solver (developer, webpage, ...) */
EXTERN
const char* SCIPsdpiGetSolverDesc(
   void
   );

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
EXTERN
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   );

/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP interface structure */
   int                   nvars,              /**< number of variables (needed in DSDP to create solver) */
   int                   nblocks,            /**< number of SDP-blocks (needed in DSDP to create SDPCone) */
   SCIP_MESSAGEHDLR*     messagehdlr         /**< message handler to use for printing messages, or NULL */
   );

/** deletes an SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP interface structure */
   );

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data into SDP solver
 *
 *  @note as the SDP-constraint matrices are symmetric, only the upper triangular part of them must be specified
 */
EXTERN
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstbegblock,   /**< start index of each block in sdpconstval-array */
   const int*            sdpconstrowind,     /**< row-index for each entry in sdpconstval-array */
   const int*            sdpconstcolind,     /**< column-index for each entry in sdpconstval-array */
   const SCIP_Real*      sdpconstval,        /**< values of entries of constant matrices in SDP-Block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   const int*            sdpbegvarblock,     /**< entry j*nvars + i is the start index of matrix \f A_i^j \f in sdpval, particularly entry i*nvars gives the starting point of block j */
   const int*            sdprowind,          /**< row-index for each entry in sdpval-array */
   const int*            sdpcolind,          /**< column-index for each entry in sdpval-array */
   const SCIP_Real*      sdpval,             /**< values of SDP-constraint matrix entries */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (going from 1 to nlpcons) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (going from 1 to nvars) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries */
   );

/** adds another SDP-Block to the problem
 *
 *  @note as \f A_i^j \f is symmetric, only the lower triangular part of it must be specified
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrices \f A_i^j \f */
   int                   nnonz,              /**< sum of non-zeroes in the lower triagonal parts of the \f A_i^j \f */
   const int*            constbegrow,        /**< start index of the rows of the matrix \f A_i^0 \f */
   const int*            constcolind,        /**< the column-indices of the non-zero entries of \f A_i^0 \f */
   const SCIP_Real*      constval,           /**< the values of \f A_i^0 \f as specified by constbegrow and constcolind */
   const int*            begvar,             /**< start index of the matrix \f A_i^j \f for each j */
   const int*            rowind,             /**< the row-indices of the non-zero entries of \f A_i^j \f */
   const int*            colind,             /**< the column-indices of the non-zero entries of \f A_i^j \f */
   const SCIP_Real*      val                 /**< the values of of \f A_i^j \f as specified by begvar, rowind and colind */
   );

/** deletes a SDP-Block */
EXTERN
SCIP_RETCODE SCIPsdpiDelSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block              /**< index of the SDP-Block that should be deleted */
   );

/** adds additional variables to the SDP
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to be added */
   const SCIP_Real*      obj,                /**< objective function values of new variables */
   const SCIP_Real*      lb,                 /**< lower bounds of new variables */
   const SCIP_Real*      ub,                 /**< upper bounds of new variables */
   char**                varnames,           /**< variable names, or NULL */
   int                   sdpnnonz,           /**< number of nonzero elements to be added to the SDP constraint matrices */
   const int*            sdpbegvar,          /**< start index of each variable in sdpval-array, or NULL if sdpnnonz == 0 */
   const int*            sdprowind,          /**< row indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const int*            sdpcolind,          /**< col indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const SCIP_Real*      sdppval,            /**< values of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   int                   lpnnonz,            /**< number of nonzero elements to be added to the LP constraint matrices */
   const int*            lpbegvar,           /**< start index of each variable in lpval-array, or NULL if lpnnonz == 0 */
   const int*            lprowind,           /**< row indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const int*            lpcolind,           /**< col indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const SCIP_Real*      lpval               /**< values of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   );

/** deletes all variables in the given range from SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to be deleted */
   int                   lastvar             /**< last variable to be deleted */
   );

/** deletes variables from SCIP_SDPI; the new position of a variable must not be greater that its old position */
EXTERN
SCIP_RETCODE SCIPsdpiDelColset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of variables
                                              *   input:  1 if variable should be deleted, 0 if not
                                              *   output: new position of variable, -1 if variable was deleted */
   );

/** adds rows to the LP-Block
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all rows in the given range from LP-Block */
EXTERN
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes LP rows from SCIP_SDPI; the new position of a row must not be greater that its old position */
EXTERN
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of LP rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole SDP */
EXTERN
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** changes lower and upper bounds of variables */
EXTERN
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   );

/** changes right hand sides of LP rows */
EXTERN
SCIP_RETCODE SCIPsdpiChgLPRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient in LP constraint matrix */
EXTERN
SCIP_RETCODE SCIPsdpiChgLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of LP-coefficient to change */
   int                   col,                /**< column number of LP-coefficient to change */
   SCIP_Real             newval              /**< new value of LP-coefficient */
   );

/** changes objective values of variables in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of variables to change objective value for */
   int*                  ind,                /**< variable indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for variables */
   );

/** changes a single coefficient in constant matrix of given SDP-Block */
EXTERN
SCIP_RETCODE SCIPsdpiChgSDPConstCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   rowind,             /**< row index */
   int                   colind,             /**< column index */
   const SCIP_Real*      val                 /**< new value of given entry of constant matrix in given SDP-Block */
   );

/** changes a single coefficient in a constraint matrix of given SDP-Block */
EXTERN
SCIP_RETCODE SCIPsdpiChgSDPCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   var,                /**< variable index */
   int                   rowind,             /**< row index */
   int                   colind,             /**< column index */
   const SCIP_Real*      val                 /**< new value of given entry of the give constraint matrix in specified SDP-Block */
   );




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/* @todo: # rows -> SDP constraints
 *         return dimension
 */
/** gets the number of LP-rows in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   );

/** gets the number of SDP-Blocks in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of rows */
   );

/** gets the number of variables in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   );

/** gets the number of nonzero elements in the SDP constraint matrices */
EXTERN
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP constraint matrcies */
   );

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
EXTERN
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   );

/** gets the number of nonzero elements in the LP Matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   );

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either sdpnnonz, sdpbegblock, sdprowind, sdpcolind and sdpval have to be NULL, or all of them have to be non-NULL, the same is true for the lp-part.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetVarInfos(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to extract information for */
   int                   lastvar,            /**< last variable to extract information for */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  sdpnnonz,           /**< pointer to store the number of nonzero elements of the sdp-constraints returned, or NULL */
   int*                  sdpbegblock,        /**< buffer to store start index of each block in sdpval-array, or NULL */
   int*                  sdprowind,          /**< buffer to store row indices of sdp constraint matrix entries, or NULL */
   int*                  sdpcolind,          /**< buffer to store column indices of sdp constraint matrix entries, or NULL */
   SCIP_Real*            sdpval,             /**< buffer to store values of sdp constraint matrix entries, or NULL */
   int*                  lpnnonz,            /**< pointer to store the number of nonzero elements of the lp-constraints returned, or NULL */
   int*                  lpbegrow,           /**< buffer to store start index of each row in lpval-array, or NULL */
   int*                  lpcolind,           /**< buffer to store column indices of lp constraint matrix entries, or NULL */
   SCIP_Real*            lpval               /**< buffer to store values of lp constraint matrix entries, or NULL */
   );

/** gets LP rows from SDP problem object; the arrays have to be large enough to store all values.
 *  rhs can be null,
 *  either nnonz, begrow, colind, and val have to be NULL, or all of them have to be non-NULL.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first LP row to get from SDP */
   int                   lastrow,            /**< last LP row to get from SDP */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begrow,             /**< buffer to store start index of each row in colind- and val-array, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets a number SDP blocks; the arrays have to be large enough to store all values.
 *  either constnnonz, constbegblock, constrow, constcolind, and constval have to be NULL, or all of them have to be non-NULL, same for the non-constant parts.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstblock,         /**< first SDP block to get from SDP */
   int                   lastblock,          /**< last SDP block to get from SDP */
   int*                  constnnonz,         /**< pointer to store the number of nonzero elements returned for the constant part, or NULL */
   int*                  constbegblock,      /**< buffer to store the start indices of the different blocks in the constval-array, or NULL */
   int*                  constrowind,        /**< buffer to store row indices of entries of constant matrices, or NULL */
   int*                  constcolind,        /**< buffer to store column indices of entries of constant matrices, or NULL */
   SCIP_Real*            constval,           /**< buffer to store values of entries of constant matrices, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begblock,           /**< buffer to store the start indices of the different blocks in the val-array, or NULL */
   int*                  varind,             /**< buffer to store variable indices of constraint matrix entries, meaning the i in \f A_i^j \f, in val-array, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets objective coefficients from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   );

/** gets current variable bounds from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   );

/** gets current right hand sides from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   );

/** gets a single coefficient of LP constraint matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   );

/** gets a single coefficient of constant SDP constraint matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetSDPConstCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   );

/** gets a single coefficient of SDP constraint matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetSDPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   var,                /**< variable index of coefficient, meaning the i in \f A_i^j \f, in val-array, or NULL */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   );


/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
EXTERN
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** gets information about primal and dual feasibility of the current SDP solution */
EXTERN
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   );

/** returns TRUE iff SDP is proven to be primal unbounded */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to be primal infeasible */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to be primal feasible */
EXTERN
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to be dual unbounded */
EXTERN
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to be dual infeasible */
EXTERN
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to be dual feasible */
EXTERN
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP was solved to optimality */
EXTERN
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff the objective limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff the iteration limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff the time limit was reached */
EXTERN
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns the internal solution status of the solver */
EXTERN
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
EXTERN
SCIP_RETCODE SCIPsdpiIgnoreInstability(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   );

/** gets objective value of solution */
EXTERN
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   );

/** gets dual solution vector for feasible SDPs */
EXTERN
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int                   dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL */
   );

/** gets the number of SDP iterations of the last solve call */
EXTERN
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   );

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetRealSolQuality(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   );

/**@} */




/*
 * SDPi State Methods
 */

/**@name SDPi State Methods */
/**@{ */

/** stores SDPi state (like solve status since last data manipulation) into sdpistate object */
EXTERN
SCIP_RETCODE SCIPsdpiGetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like solve status since last data manipulation) */
   );

/** loads SDPi state into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetState()
 */
EXTERN
SCIP_RETCODE SCIPsdpiSetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE*       sdpistate           /**< SDPi state information */
   );

/** clears current SDPi state (like basis information) of the solver */
EXTERN
SCIP_RETCODE SCIPsdpiClearState(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** frees SDPi state information */
EXTERN
SCIP_RETCODE SCIPsdpiFreeState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information */
   );

/** reads SDPi state from a file */
EXTERN
SCIP_RETCODE SCIPsdpiReadState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/** writes SDPi state to a file */
EXTERN
SCIP_RETCODE SCIPsdpiWriteState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */


/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   );

/** sets integer parameter of SDP */
EXTERN
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   );

/** gets floating point parameter of SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   );

/** sets floating point parameter of SDP */
EXTERN
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   );

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP solver */
EXTERN
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   );

/** checks if given value is treated as infinity in the SDP solver */
EXTERN
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   );

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
EXTERN
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/** writes SDP to a file */
EXTERN
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
