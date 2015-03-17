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

/**@file   sdpi_none.c
 * @ingroup SDPIS
 * @brief  dummy interface for the case no SDP solver is needed
 * @author Stefan Heinz
 * @author Leif Naundorf
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "sdpi/sdpi.h"
#include "scip/pub_message.h"

#define SDPINAME          "NONE"              /**< name of the SDPI interface */
#define SDPIINFINITY       1e20               /**< infinity value */

/** SDP interface
 *
 *  Store several statistic values about the SDP. These values are only needed in order to provide a rudimentary
 *  communication, e.g., there are asserts that check the number of rows and columns.
 */
struct SCIP_SDPi
{
   int                   nrows;              /**< number of rows */
   int                   ncols;              /**< number of columns */
   int                   nsdpvars;           /**< number of semidefinite variables */
};


/*
 * Local Methods
 */

/** error handling method */
static
void errorMessageAbort(
   void
   )
{
   SCIPerrorMessage("No SDP solver available (SDPS=none).\n");
   SCIPerrorMessage("Ensure <sdp/solvefreq = -1>; note that continuous variables might require an SDP-solver.\n");
   SCIPABORT();
}

/** error handling method */
static
void errorMessage(
   void
   )
{
   SCIPerrorMessage("No SDP solver available (SDPS=none).\n");
   SCIPerrorMessage("Ensure <sdp/solvefreq = -1>; note that continuous variables might require an SDP-solver.\n");
}

/*
 * SDP Interface Methods
 */


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of SDP solver */
const char* SCIPsdpiGetSolverName(
   void
   )
{
   return SDPINAME;
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   void
   )
{
   return "dummy SDP solver interface which solely purpose is to resolve references at linking";
}

/** gets pointer for SDP solver - use only with great care */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*             sdpi                 /**< pointer to an SDP interface structure */
   )
{  /*lint --e{715}*/
   return (void*) NULL;
}
/**@} */




/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{  /*lint --e{715}*/
   assert(sdpi != NULL);
   SCIPdebugMessage("SCIPsdpiCreate()\n");
   SCIPdebugMessage("Note that there is no SDP solver linked to the binary\n");

   /* create empty SDPI */
   SCIP_ALLOC( BMSallocMemory(sdpi) );
   (*sdpi)->nrows = 0;
   (*sdpi)->ncols = 0;
   (*sdpi)->nsdpvars = 0;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**            sdpi                 /**< pointer to an SDP interface structure */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   SCIPdebugMessage("SCIPsdpiFree()\n");

   BMSfreeMemory(sdpi);

   return SCIP_OKAY;
}

/**@} */




/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data with column matrix into SDP solver */
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
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );

   sdpi->nrows = nrows;
   sdpi->ncols = ncols;
   assert( sdpi->nrows >= 0 );
   assert( sdpi->ncols >= 0 );

   return SCIP_OKAY;
}

/** adds a term of the form <C,X_j> to the objective funtion */
SCIP_RETCODE SCIPsdpiAddSDPTermObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix C */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of C */
   const int*            subi,               /**< the row-indices of the non-zero entries of C */
   const int*            subj,               /**< the column-indices of the non-zero entries of C */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   ind                 /**< index of the symmetric variable X */
   )
{
   return SCIP_OKAY;
}

/** adds a term of the form <A,X_j> to a row of the SDP
 *
 *  @note as A is symmetric, only the lower triangular part of it must be specified
 */
SCIP_RETCODE SCIPsdpiAddSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix A */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of A */
   const int*            subi,               /**< the row-indices of the non-zero entries of A */
   const int*            subj,               /**< the column-indices of the non-zero entries of A */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   ind,                /**< index of the symmetric variable X */
   int                   row                 /**< row of the SDP where this term should be added */
   )
{
   return SCIP_OKAY;
}

/** deletes a term of the form <A,X_j> from a row of the SDP */
SCIP_RETCODE SCIPsdpiDelSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ind,                /**< index of the symmetric variable X whose term should be deleted */
   int                   row                 /**< row of the SDP where this term should be deleted */
   )
{
   return SCIP_OKAY;
}

/** adds a number of positive-semidefinite variables X_j to the SDP */
SCIP_RETCODE SCIPsdpiAddSDPVars(
   SCIP_SDPI*           sdpi,                /**< SDP interface structure */
   int                  nvars,               /**< number of semidefinite variables to be added */
   const int*           dims                 /**< dimensions of the semidefinite variables to be added */
   )
{
   assert( sdpi != NULL );
   assert( sdpi-> nsdpvars >= 0 );

   sdpi->nsdpvars += nvars;

   return SCIP_OKAY;
}

/** deletes a number of positive-semidefinite variables X_j from the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiDelSDPVars(
   SCIP_SDPI*           sdpi,                /**< SDP interface structure */
   int*                 dstat                /**< deletion status of SDP variables
                                              *   input:  1 if SDP variable should be deleted, 0 if not
                                              *   output: new indices of SDP variable, -1 if variable was deleted */
   )
{
   int count = 0;
   int i;

   assert( sdpi != NULL );
   assert( dstat != NULL );

   for( i = 0; i < sdpi->nsdpvars; i++ )
   {
      if( dstat[i] == 1 )
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = count;
         count++;
      }
   }
   sdpi->nsdpvars = count;
   assert( sdpi-> nsdpvars >= 0 );

   return SCIP_OKAY;
}

/** adds columns to the SDP */
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
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( sdpi->ncols >= 0 );

   sdpi->ncols += ncols;

   return SCIP_OKAY;
}

/** deletes all columns in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( sdpi->ncols >= 0 );

   sdpi->ncols -= lastcol - firstcol + 1;
   assert( sdpi->ncols >= 0 );

   return SCIP_OKAY;
}

/** deletes columns from SCIP_SDP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelColset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{  /*lint --e{715}*/
   int cnt = 0;
   int j;

   assert( sdpi != NULL );
   assert( dstat != NULL );
   assert( sdpi->ncols >= 0 );

   for (j = 0; j < sdpi->ncols; ++j)
   {
      if ( dstat[j] )
      {
         dstat[j] = -1;
      }
      else
      {
         cnt++;
         dstat[j] = cnt;
      }
   }
   assert( sdpi->ncols >= cnt );
   sdpi->ncols = cnt;
   assert( sdpi->ncols >= 0 );

   return SCIP_OKAY;
}

/** adds rows to the SDP */
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
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( sdpi->nrows >= 0 );

   sdpi->nrows += nrows;

   return SCIP_OKAY;
}

/** deletes all rows in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( sdpi->nrows >= 0 );

   sdpi->nrows -= lastrow - firstrow + 1;
   assert( sdpi->nrows >= 0 );

   return SCIP_OKAY;
}

/** deletes rows from SCIP_SDP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{  /*lint --e{715}*/
   int cnt = 0;
   int i;

   assert( sdpi != NULL );
   assert( dstat != NULL );
   assert( sdpi->nrows >= 0 );

   for (i = 0; i < sdpi->nrows; ++i)
   {
      if ( dstat[i] )
      {
         ++cnt;
         dstat[i] = -1;
      }
      else
         dstat[i] = cnt;
   }
   sdpi->nrows -= cnt;
   assert( sdpi->nrows >= 0 );

   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( sdpi->nrows >= 0 );
   assert( sdpi->ncols >= 0 );

   sdpi->nrows = 0;
   sdpi->ncols = 0;
   sdpi->nsdpvars = 0;

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPsdpiChgSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPsdpiChgCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPsdpiChgObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** changes objective values of columns in the SDP */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPsdpiScaleRow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPsdpiScaleCol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/**@} */




/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( nrows != NULL );
   assert( sdpi->nrows >= 0 );

   *nrows = sdpi->nrows;

   return SCIP_OKAY;
}

/** gets the number of columns in the SDP */
SCIP_RETCODE SCIPsdpiGetNCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL );
   assert( ncols != NULL );
   assert( sdpi->ncols >= 0 );

   *ncols = sdpi->ncols;

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{  /*lint --e{715}*/
   assert(nnonz != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets the number of nonzero semidefinite terms of the form <A,X_j> in the SDP */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPTerms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpterms           /**< pointer to store the number of SDP terms */
   )
{
   return SCIP_PLUGINNOTFOUND;
}

/** gets the number of semidefinite variables */
EXTERN
SCIP_RETCODE SCIPsdpiGetNSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpvars            /**< pointer to store the number of SDP variables */
   )
{
   assert( sdpi != NULL);
   assert( sdpi->nsdpvars >= 0 );

   *nsdpvars = sdpi->nsdpvars;

   return SCIP_OKAY;
}

/** gets columns from SDP problem object; the arrays have to be large enough to store all values
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets rows from SDP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
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
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets column names */
SCIP_RETCODE SCIPsdpiGetColNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get name from SDP */
   int                   lastcol,            /**< last column to get name from SDP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets row names */
SCIP_RETCODE SCIPsdpiGetRowNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get name from SDP */
   int                   lastrow,            /**< last row to get name from SDP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets the objective sense of the SDP */
SCIP_RETCODE SCIPsdpiGetObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets current bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets current row sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{  /*lint --e{715}*/
   assert(firstrow <= lastrow);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPsdpiGetCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets a symmetric matrix of the SDP */
SCIP_RETCODE SCIPsdpiGetSymmetricMatrix(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   index,              /**< index of the semidefinite matrix */
   int                   length,             /**< length of the output arrays subi, subj and valij */
   int*                  subi,               /**< Row indices of nonzero entries */
   int*                  subj,               /**< Column indices of nonzero entries */
   SCIP_Real*            valij               /**< values of the nonzero entries defined by the arrays subi and subj */
   )
{
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** calls primal simplex to solve the SDP */
SCIP_RETCODE SCIPsdpiSolvePrimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** calls dual simplex to solve the SDP */
SCIP_RETCODE SCIPsdpiSolveDual(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** calls barrier or interior point algorithm to solve the SDP with crossover to simplex basis */
SCIP_RETCODE SCIPsdpiSolveBarrier(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPsdpiStartStrongbranch(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL);
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPsdpiEndStrongbranch(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   assert( sdpi != NULL);
   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPsdpiStrongbranchFrac(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{  /*lint --e{715}*/
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPsdpiStrongbranchesFrac(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{  /*lint --e{715}*/
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPsdpiStrongbranchInt(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{  /*lint --e{715}*/
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPsdpiStrongbranchesInt(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{  /*lint --e{715}*/
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{  /*lint --e{715}*/
   assert(primalfeasible != NULL);
   assert(dualfeasible != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPsdpiExistsPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPsdpiHasPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal feasible */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPsdpiExistsDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPsdpiHasDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual feasible */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff current SDP basis is stable */
SCIP_Bool SCIPsdpiIsStable(
   SCIP_SDPI*            sdpi               /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** returns the internal solution status of the solver */
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiIgnoreInstability(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{  /*lint --e{715}*/
   assert(success != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{  /*lint --e{715}*/
   assert(objval != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets the primal solution for a semidefinite variable */
SCIP_RETCODE SCIPsdpiGetSolSDPVar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ind,                /**< index of semidefinite variable */
   SCIP_Real*            barxj               /**< primal semidefinite solution vector, must have size dim*(dim+1)/2 */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets primal ray for unbounded SDPs */
SCIP_RETCODE SCIPsdpiGetPrimalRay(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPsdpiGetDualfarkas(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{  /*lint --e{715}*/
   assert(dualfarkas != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{  /*lint --e{715}*/
   assert(iterations != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPsdpiGetRealSolQuality(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   assert(sdpi != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/**@} */




/*
 * SDP Basis Methods
 */

/**@name SDP Basis Methods */
/**@{ */

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPsdpiGetBase(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPsdpiSetBase(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{  /*lint --e{715}*/
   assert(cstat != NULL);
   assert(rstat != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
extern
SCIP_RETCODE SCIPsdpiGetBasisInd(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPsdpiGetBInvRow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   r,                  /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPsdpiGetBInvCol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the SDP;
                                              *   you have to call SCIPsdpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the SDP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPsdpiGetBInvARow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   r,                  /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPsdpiGetBInvRow(), or NULL */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPsdpiGetBInvACol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * SDP State Methods
 */

/**@name SDP State Methods */
/**@{ */

/** stores SDPi state (like basis information) into sdpistate object */
SCIP_RETCODE SCIPsdpiGetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(blkmem != NULL);
   assert(sdpistate != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** loads SDPi state (like basis information) into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetState()
 */
SCIP_RETCODE SCIPsdpiSetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE*       sdpistate           /**< SDPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(blkmem != NULL);
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** clears current SDPi state (like basis information) of the solver */
SCIP_RETCODE SCIPsdpiClearState(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   assert(sdpi != NULL);
   return SCIP_OKAY;
}

/** frees SDPi state information */
SCIP_RETCODE SCIPsdpiFreeState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like basis information) */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** checks, whether the given SDP state contains simplex basis information */
SCIP_Bool SCIPsdpiHasStateBasis(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPISTATE*       sdpistate           /**< SDP state information (like basis information) */
   )
{  /*lint --e{715}*/
   errorMessageAbort();
   return FALSE;
}

/** reads SDP state (like basis information from a file */
SCIP_RETCODE SCIPsdpiReadState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** writes SDP state (like basis information) to a file */
SCIP_RETCODE SCIPsdpiWriteState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * SDP Pricing Norms Methods
 */

/**@name SDP Pricing Norms Methods */
/**@{ */

/** stores SDPi pricing norms information
 *  @todo should we store norm information?
 */
SCIP_RETCODE SCIPsdpiGetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** loads SDPi pricing norms into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetNorms()
 */
SCIP_RETCODE SCIPsdpiSetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS*       sdpinorms           /**< SDPi pricing norms information */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPsdpiFreeNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */




/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   assert(ival != NULL);
   return SCIP_PARAMETERUNKNOWN;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/** gets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{  /*lint --e{715}*/
   assert(dval != NULL);
   return SCIP_PARAMETERUNKNOWN;
}

/** sets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{  /*lint --e{715}*/
   return SCIP_PARAMETERUNKNOWN;
}

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverGetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< parameter value */
   )
{
   return SCIP_PARAMETERUNKNOWN;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSolverSetIntpar(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   return SCIP_PARAMETERUNKNOWN;
}

/**@} */

/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP solver */
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   return SDPIINFINITY;
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   if( val >= SDPIINFINITY )
      return TRUE;
   return FALSE;
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{  /*lint --e{715}*/
   errorMessage();
   return SCIP_PLUGINNOTFOUND;
}

/**@} */
