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
 *
 * This file specifies a generic SDP solver interface used by SCIP to create, modify, and solve semidefinite programs of
 * the form
 *
 *   min/max   obj * x + \sum_j <C_j,X_j>
 *      lhs <=   A * x   + \sum_j <A_ij,X_j> <= rhs
 *      lb  <=       x  <= ub
 *      X_j psd
 *
 * and query information about the solution.  See also http://docs.mosek.com/7.0/capi/Semidefinite_optimization.html.
 *
 * Although it includes a few SCIP header files, e.g., because it uses SCIP's return codes, it can be used independently of
 * any SCIP instance.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPI_H__
#define __SCIP_SDPI_H__


#include "scip/def.h"
#include "scip/blockmemshell/memory.h"
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
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
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

/** copies SDP data with column matrix into SDP solver */
EXTERN
SCIP_RETCODE SCIPsdpiLoadColSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   );

/** adds a number of positive-semidefinite variables X_j to the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiAddSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of semidefinite variables to be added */
   const int*            dims                /**< dimensions of the semidefinite variables to be added */
   );

/** deletes a number of positive-semidefinite variables X_j from the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of SDP variables
                                              *   input:  1 if SDP variable should be deleted, 0 if not
                                              *   output: new indices of SDP variable, -1 if variable was deleted */
   );

/** adds a term of the form <C,X_j> to the objective funtion */
EXTERN
SCIP_RETCODE SCIPsdpiAddSDPTermObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix C */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of C */
   const int*            subi,               /**< the row-indices of the non-zero entries of C */
   const int*            subj,               /**< the column-indices of the non-zero entries of C */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   ind                 /**< index of the symmetric variable X */
   );

/** adds a term of the form <A,X_j> to a row of the SDP
 *
 *  @note as A is symmetric, only the lower triangular part of it must be specified
 *  @note If there is already a matrix specified at the given position (ind,row), the old matrix will be overwritten
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix A */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of A */
   const int*            subi,               /**< the row-indices of the non-zero entries of A */
   const int*            subj,               /**< the column-indices of the non-zero entries of A */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   ind,                /**< index of the symmetric variable X */
   int                   row                 /**< row of the SDP where this term should be added */
   );

/** deletes a term of the form <A,X_j> from a row of the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ind,                /**< index of the symmetric variable X whose term should be deleted */
   int                   row                 /**< row of the SDP where this term should be deleted */
   );

/** adds columns to the SDP
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of new columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all columns in the given range from SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   );

/** deletes columns from SCIP_SDPI; the new position of a column must not be greater that its old position */
EXTERN
SCIP_RETCODE SCIPsdpiDelColset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   );

/** adds rows to the SDP
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
EXTERN
SCIP_RETCODE SCIPsdpiAddRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   );

/** deletes all rows in the given range from SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   );

/** deletes rows from SCIP_SDPI; the new position of a row must not be greater that its old position */
EXTERN
SCIP_RETCODE SCIPsdpiDelRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   );

/** clears the whole SDP */
EXTERN
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** changes lower and upper bounds of columns */
EXTERN
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   );

/** changes left and right hand sides of rows */
EXTERN
SCIP_RETCODE SCIPsdpiChgSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   );

/** changes a single coefficient */
EXTERN
SCIP_RETCODE SCIPsdpiChgCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   );

/** changes the objective sense */
EXTERN
SCIP_RETCODE SCIPsdpiChgObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   );

/** changes objective values of columns in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   );

/**@todo add change methods for SDP variables and matrices */
/** Leif: MOSEK currently has no function to directly change the SDP variables.
 *        The only way is by using SCIPsdpiDelSDPVars() and then SCIPsdpiAddSDPVars().
 *        We can implement this by an extra function, do we really need this?
 *
 *        To change a matrix term of the form <A,X_j> in a row is supported by using
 *        the function SCIPsdpiAddSDPTerm(). If we define a position which is already specified,
 *        the old matrix at this position will be overwritten by the new matrix.
 *        Again, do we need to implement a special SCIPsdpiChgSDPTerm() function for this?
 */

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/* @todo: # rows -> SDP constraints
 *         return dimension
 */
/** gets the number of rows in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   );

/** gets the number of columns in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   );

/** gets the objective sense of the SDP */
SCIP_RETCODE SCIPsdpiGetObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   );

/** gets the number of nonzero elements in the SDP constraint matrix */
EXTERN
SCIP_RETCODE SCIPsdpiGetNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   );

/** gets the number of nonzero semidefinite terms of the form <A,X_j> in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPTerms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpterms           /**< pointer to store the number of SDP terms */
   );

/** gets the number of semidefinite variables */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpvars            /**< pointer to store the number of SDP variables */
   );

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get from SDP */
   int                   lastcol,            /**< last column to get from SDP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets rows from SDP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
EXTERN
SCIP_RETCODE SCIPsdpiGetRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get from SDP */
   int                   lastrow,            /**< last row to get from SDP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   );

/** gets column names */
EXTERN
SCIP_RETCODE SCIPsdpiGetColNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get name from SDP */
   int                   lastcol,            /**< last column to get name from SDP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   );

/** gets row names */
EXTERN
SCIP_RETCODE SCIPsdpiGetRowNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get name from SDP */
   int                   lastrow,            /**< last row to get name from SDP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   );

/** gets objective coefficients from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   );

/** gets current bounds from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   );

/** gets current row sides from SDP problem object */
EXTERN
SCIP_RETCODE SCIPsdpiGetSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   );

/** gets a single coefficient */
EXTERN
SCIP_RETCODE SCIPsdpiGetCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   );

/** gets a symmetric matrix of the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetSymmetricMatrix(
   SCIP_SDPI*           sdpi,                /**< SDP interface structure */
   int                  index,               /**< index of the semidefinite matrix */
   int                  length,              /**< length of the output arrays subi, subj and valij */
   int*                 subi,                /**< Row indices of nonzero entries */
   int*                 subj,                /**< Column indices of nonzero entries */
   SCIP_Real*           valij                /**< values of the nonzero entries defined by the arrays subi and subj */
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

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
EXTERN
SCIP_Bool SCIPsdpiExistsPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
EXTERN
SCIP_Bool SCIPsdpiHasPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
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

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
EXTERN
SCIP_Bool SCIPsdpiExistsDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   );

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
EXTERN
SCIP_Bool SCIPsdpiHasDualRay(
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

/** returns TRUE iff current SDP basis is stable */
EXTERN
SCIP_Bool SCIPsdpiIsStable(
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

/** gets primal and dual solution vectors for feasible SDPs */
EXTERN
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   );

/** gets the primal solution for a semidefinite variable */
EXTERN
SCIP_RETCODE SCIPsdpiGetSolSDPVar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ind,                /**< index of semidefinite variable */
   SCIP_Real*            barxj               /**< primal semidefinite solution vector, must have size dim*(dim+1)/2 */
   );

/** gets primal ray for unbounded SDPs */
EXTERN
SCIP_RETCODE SCIPsdpiGetPrimalRay(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   );

/** gets dual Farkas proof for infeasibility */
EXTERN
SCIP_RETCODE SCIPsdpiGetDualfarkas(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
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

/** stores SDPi state (like basis information) into sdpistate object */
EXTERN
SCIP_RETCODE SCIPsdpiGetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like basis information) */
   );

/** loads SDPi state (like basis information) into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetState()
 */
EXTERN
SCIP_RETCODE SCIPsdpiSetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE*       sdpistate           /**< SDPi state information (like basis information) */
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
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like basis information) */
   );

/** checks, whether the given SDPi state contains simplex basis information */
EXTERN
SCIP_Bool SCIPsdpiHasStateBasis(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPISTATE*       sdpistate           /**< SDPi state information (like basis information) */
   );

/** reads SDPi state (like basis information from a file */
EXTERN
SCIP_RETCODE SCIPsdpiReadState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/** writes SDPi state (like basis information) to a file */
EXTERN
SCIP_RETCODE SCIPsdpiWriteState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   );

/**@} */


/*
 * SDPi Pricing Norms Methods
 */

/**@name SDPi Pricing Norms Methods */
/**@{ */

/** stores SDPi pricing norms into sdpinorms object */
extern
SCIP_RETCODE SCIPsdpiGetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
   );

/** loads SDPi pricing norms into solver; note that the SDP might have been extended with additional
 *  columns and rows since the norms were stored with SCIPsdpiGetNorms()
 */
extern
SCIP_RETCODE SCIPsdpiSetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS*       sdpinorms           /**< SDPi pricing norms information */
   );

/** frees SDPi pricing norms information */
extern
SCIP_RETCODE SCIPsdpiFreeNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
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
