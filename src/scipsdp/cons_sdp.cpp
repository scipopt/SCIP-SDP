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
//#define SCIP_DEBUG
/**@file   cons_sdp.cpp
 * @brief  constraint handler for sdp-constraints
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 */

#include "cons_sdp.h"

#include <cassert>                      // for assert
//#include <cmath>                        // for floor //TODO: lint says it's not needed
#include <cstring>                      // for NULL, strcmp

#include "config.h"                     // for F77_FUNC
#include "stdlib.h"                     /* for fabs */

#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/SdpVarfixer.h"

#include "scip/cons_linear.h"           // for SCIPcreateConsLinear
#include "scip/scip.h"                  // for SCIPallocBufferArray, etc

#define CONSHDLR_NAME          "SDP"
#define CONSHDLR_DESC          "SDP constraints of the form \\sum_{j} A_j y_j - A_0 psd"
#define CONSHDLR_SEPAPRIORITY  +1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/** constraint data for sdp constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in this SDP constraint */
   int                   nnonz;              /**< number of nonzeroes in this SDP constraint */
   int                   blocksize;          /**< size of this SDP-block */
   int*                  nvarnonz;           /**< length of the arrays pointed to by col/row/val, number of nonzeros for each variable */
   int**                 col;                /**< pointers to the column indices of the nonzeroes for each variable */
   int**                 row;                /**< pointers to the row indices of the nonzeroes for each variable */
   SCIP_Real**           val;                /**< pointers to the values of the nonzeroes for each variable*/
   //SCIP_Var**            vars;               /**< SCIP_Variables present in this SDP constraint, ordered by their begvar-indices */
   SdpVarmapper          varmapper;          /**< maps SCIP_Variables to their position in begvar (if they are present in this constraint) and vice versa */
   int                   constnnonz;         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol;           /**< column indices of the constant nonzeroes */
   int*                  constrow;           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval;           /**< values of the constant nonzeroes */
};

struct SCIP_ConshdlrData
{
   SdpVarmapper          varmapper;          /**< maps SCIP variables to their global SDP indices and vice versa */
};


extern "C" {
/** BLAS Fortran subroutine DGEMV */
void F77_FUNC(dgemv, DGEMV)(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
}

/** call matrix-vector multipication
 *
 *  @note all memory must be allocated outside
 */
static
SCIP_RETCODE Blas_DGEMV(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   SCIP_Real             alpha,              /**< scaling parameter */
   double*               matrix,             /**< the matrix we want to multiply */
   double*               vector,             /**< vector we want to multiply with the matrix */
   double                beta,               /**< scaling parameter */
   double*               result              /**< vector where the result is put in */
   )
{
   /* store everything in local variables????????? */
   char TRANS = 'N';
   int M = nrows;
   int N = ncols;
   double ALPHA = alpha;
   double* A = matrix;
   int LDA = nrows;
   double* X = vector;
   int INCX = 1;
   double BETA = beta;
   double* Y = result;
   int INCY = 1;

   F77_FUNC(dgemv, DGEMV)(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);

   return SCIP_OKAY;
}

extern "C" {
/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)(
   char* JOBZ, char* RANGE, char* UPLO,
   int* N, double* A, int* LDA,
   double* VL, double* VU,
   int* IL, int* IU,
   double* ABSTOL, int* M, double* W, double* Z,
   int* LDZ, int* ISUPPZ, double* WORK,
   int* LWORK, int* IWORK, int* LIWORK,
   int* INFO );
}

/** computes the i-th eigenvalue, where 1 is the smallest and n the largest */
static
SCIP_RETCODE computeIthEigenvalue(
   SCIP*                 scip,               /**< SCIP data structure*/
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   double*               A,                  /**< matrix for which eigenvalues should be computed */
   int                   i,                  /**< index of eigenvalue to be computed */
   double*               eigenvalue,         /**< pointer to store eigenvalue */
   double*               eigenvector         /**< pointer to array to store eigenvector */
   )
{
   /* store everything in local variables????????? */
   int     N = n;
   int     INFO;
   char    JOBZ = geteigenvectors ? 'V' : 'N';
   char    RANGE = 'I';
   char    UPLO = 'L';
   int     LDA  = N;
   double* WORK;
   int     LWORK;
   int*    IWORK;
   int     LIWORK;
   //    int*    ISUPPZ;
   double* WTMP;
   double  ABSTOL = 0.0;
   int     IL = i;
   int     IU = i;
   int     M = 1;
   int     LDZ = LDA;
   double  WSIZE;
   int     WISIZE;

   // standard LAPACK workspace query
   LWORK = -1;
   LIWORK = -1;

   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      &IL, &IU,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO );

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   // allocate workspace
   LWORK = (int) WSIZE + 1;
   LIWORK = WISIZE;

   SCIP_CALL( SCIPallocBufferArray(scip, &WORK, LWORK) );
   SCIP_CALL( SCIPallocBufferArray(scip, &IWORK, LIWORK) );

   SCIP_CALL( SCIPallocBufferArray(scip, &WTMP, N) );

   // call the function
   double VL = -1e20;
   double VU = 1e20;
   int ISUPPZ[2];

   /* @todo Can NULL be passed via eigenvector? */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      &IL, &IU,
      &ABSTOL, &M, WTMP, eigenvector,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO );

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   // handle output
   *eigenvalue = WTMP[0];

   SCIPfreeBufferArray(scip, &WORK);
   SCIPfreeBufferArray(scip, &IWORK);
   SCIPfreeBufferArray(scip, &WTMP);

   return SCIP_OKAY;
}

/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static
int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}

/**
 * takes an 0.5*n*(n+1) array of a symmetric matrix and expands it to a n*n array of the full matrix to input into LAPACK
 */
static
SCIP_RETCODE expandSymMatrix(
   int                   size,               /**< size of the matrix, named n below */
   SCIP_Real*            symMat,             /**< symmetric matrix indexed via compLowerTriangPos that should be expanded */
   SCIP_Real*            fullMat             /**< n*n matrix, that is the symmetric expansion of symMat */
   )
{
   int i;
   int j;
   int ind;

   assert ( size >= 0 );
   assert ( symMat != NULL );
   assert ( fullMat != NULL );

   /* traverse the lower triangular part in the order of the indices and copy the values to both lower and upper triangular part */
   for (i = 0; i < size; i++)
   {
      for (j = 0; j <= i; j++)
      {
         assert ( ind == compLowerTriangPos(i,j) );
         fullMat[i*blocksize + j] = symMat[ind];
         fullMat[j*blocksize + i] = symMat[ind];
         ind++;
      }
   }

   return SCIP_OKAY;
}

/**
 * For a given vector \f y \f computes the SDP-Matrix \f \sum_{j=1}^n A_j y_j - A_0 \f for this SDP block.
 * Length of the matrix array needs to be (length of y) * (length of y + 1) /2, this will be indexed by compLowerTriangPos
 */
static
SCIP_RETCODE computeSdpMatrix(
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             y,                  /**< solution to separate */
   SCIP_Real*            matrix              /**< the SDP-Matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int ind;
   int nvars;
   int endindex;

   assert ( cons != NULL );
   assert ( y != NULL );
   assert ( matrix != NULL );

   consdata = SCIPconsGetData(cons);
   nvars = constdata->nvars;

   /* initialize the matrix with 0 */
   for (i = 0; i < 0.5 * nvars * (nvars + 1); i++)
      matrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      for (ind = 0; i < consdata->nvarnonz[i]; i++)
         matrix[compLowerTriangPos(consdata->row[ind], consdata->col[ind])] += y[i] * consdata->val[i][ind];
   }

   /* substract the constant part */
   for (ind = 0; ind < constnnonz; ind++)
      matrix[compLowerTriangPos(consdata->constrow[ind], consdata->constcol[ind])] -= consdata->constval[ind];

   return SCIP_OKAY;
}

/**
 * For a given variable-index j and a Vector v computes \f v^T A_j v \f.
 */
static
SCIP_RETCODE multiplyConstraintMatrix(
   SCIP_CONS*            cons,               /**< the SDP constraint that includes the Matrix \f A_j \f */
   int                   j,                  /**< variable-index of the matrix to multiply with */
   SCIP_Real*            v,                  /**< vector to multiply with */
   SCIP_Real*            vAv                 /**< the resulting scalar \f v^T A_j v \f */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert ( cons != NULL );
   assert ( j >= 0 );
   assert ( vav != NULL );

   consdata = SCIPconsGetData(cons);

   assert ( j < consdata->nvars );

   /* initialize the product with 0 */
   vAv = 0.0;

   for (i = 0; i < consdata->nvarnonz[i]; i++)
   {
      if (consdata->col[j][i] == consdata->row[j][i])
         vAv += v[consdata->col[j][i]] * consdata->val[i] * v[consdata->row[j][i]];
      else
      {
         vAv += 2.0 * v[consdata->col[j][i]] * consdata->val[i] * v[consdata->row[j][i]]; /* *2 because the matrix is symmetric and there is one identical
                                                                                           * contribution each from lower and upper triangular part */
      }
   }

   return SCIP_OKAY;
}

/**
 * Get the maximum absolute value of an entry of the constant matrix.
 */
static
SCIP_Real* getMaxConstEntry(
   SCIP_CONS*            cons                /**< the SDP constraint that includes the Matrix \f A_j \f */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   SCIP_Real* max;

   consdata = SCIPconsGetData(cons);

   /* initialize max with the absolute value of the first entry of the constant matrix */
   max = fabs(consdata->constval[0]);

   /* iterate over the remaining arrays, updating max if a higher absolute value is found */
   for (i = 1; i < consdata->constnnonz; i++)
   {
      if (fabs(consdata->constval[i]) > max)
         max = fabs(consdata->constval[i]);
   }

   return max;
}


/** separate current solution using a cut, with the eigenvectors and -values of the solution matrix
 *
 * This function computes the eigenvectors of the matrix, takes the first one and multiplies the matrix with it
 * \f x^T*A_i*x = coeff[i], x^T*A_0*x =lhs \f
 */
static
SCIP_RETCODE cutUsingEigenvector(
   SCIP*                 scip,               /**<SCIP data structure*/
   SCIP_CONS*            cons,               /**<the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**<solution to separate*/
   SCIP_Real*            coeff,              /**<coefficients of the computed cut*/
   SCIP_Real*            lhs                 /**<lhs of the computed cut*/
   )
{
   SCIP_CONSDATA* consdata;
   consdata = SCIPconsGetData(cons);

   *lhs = 0.0;
   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real* fullconstmatrix = NULL;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real* output_vector = NULL;
   int blocksize = consdata->blocksize;


   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1) ));
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ));
   SCIP_CALL( SCIPallocBufferArray(scip, &fullconstmatrix, blocksize * blocksize));
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvector, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &output_vector, blocksize) );

   //compute the matrix \f \sum_j A_j y_j \f
   SCIP_CALL( computeSdpMatrix(cons, sol, matrix));

   //expand it because LAPACK wants the full matrix instead of the lower triangular part
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix));

   SCIP_CALL( computeIthEigenvalue(scip, TRUE, blocksize, matrix, 1, eigenvalues, eigenvector) );

   //get full constant matrix
   SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, fullconstmatrix));

   /* multiply eigenvector with constant matrix to get lhs (after multiplying again with eigenvector from the left */
   SCIP_CALL(Blas_DGEMV(blocksize, blocksize, 1.0, const_matrix, eigenvector, 0.0, output_vector));

   for (int j = 0; j < blocksize; ++j)
   {
      *lhs += eigenvector[j] * output_vector[j];
   }


   /* compute \f v^T A_j v \f for eigenvector v and each Matrix \f A_j \f to get the coefficients of the LP cut */
   for (int j = 0; j < consdata->nvars; ++j)
   {
      multiplyConstraintMatrix(cons, j, eigenvector, &coeff[j]);
   }

   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &eigenvalues);
   SCIPfreeBufferArray(scip, &fullconstmatrix);
   SCIPfreeBufferArray(scip, &eigenvector);
   SCIPfreeBufferArray(scip, &output_vector);

   return SCIP_OKAY;
}

/** checks feasibility for a single SDP-Cone */
SCIP_RETCODE checkSdpCon(
   SCIP*                scip,               /**< SCIP data structure */
   SCIP_CONS*           cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*            sol,                /**< the solution to check feasibility for */
   SCIP_Bool            checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool            checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool            printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*         result              /**< pointer to store the result of the feasibility checking call */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   double check_value;

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * blocksize) );

   SCIP_CALL( computeSdpMatrix(cons, sol, matrix) );
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, fullmatrix, 1, eigenvalues, NULL) );

   //we are going to use one of the dimacs error norms for checking feasiblity.
   //we use the second one: err=max{0, -lambda_min(x)/(1+maximumentry of rhs}

   check_value = (-eigenvalues[0]) / (1 + getMaxConstEntry(cons));

   *result = SCIPisLE(scip, check_value, 0.0) ? SCIP_FEASIBLE : SCIP_INFEASIBLE;

#ifdef SCIP_DEBUG
   if( *result == SCIP_INFEASIBLE)
   {
      SCIPdebugMessage("In cons_check a matrix was found not to be sdp because of eigenvalue = %f, dimacs error norm = %f \n",
            eigenvalues[0], check_value);
   }
#endif

   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &eigenvalues);

   // @todo checkintegrality ?!?, checklprows ?!?, printreason ?!?

   return SCIP_OKAY;
}

/** separates the current solution */
static
SCIP_RETCODE separateSol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         conss,              /**< array of constraints to process */
   SCIP_SOL*          sol,                /**< primal solution that should be separated */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */)
{
   SCIP_CONSDATA* consdata;
   int nvars;
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff = NULL;

   consdata = SCIPconsGetData(cons);

   nvars = consdata->nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars ) );

   SCIP_CALL(cutUsingEigenvector(scip, sdpcone, sol, coeff, &lhs));

   SCIP_COL** cols;
   SCIP_Real* rcoeff = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rcoeff, nvars ) );

   int len = 0;

   for (int j = 0; j < nvars; ++j)
   {
      if (SCIPisZero(scip, coeff[j]))
      {
         continue;
      }

      cols[len] = SCIPvarGetCol(SCIPvarGetProbvar(sdpcone->get_var(j)));
      rcoeff[len] = coeff[j];
      ++len;
   }

   SCIP_ROW* row;
   SCIP_CALL( SCIPcreateRowCons(scip, &row, conshdlr, "sepa_eig_sdp", len, cols, rcoeff, lhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );

   if (SCIPisCutEfficacious(scip, sol, row) )
   {
      SCIP_Bool infeasible;
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );
      if ( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;
      SCIP_CALL( SCIPresetConsAge(scip, conss) );
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeBufferArray(scip, &cols );
   SCIPfreeBufferArray(scip, &rcoeff );
   SCIPfreeBufferArray(scip, &coeff );

   return SCIP_OKAY;

}

/** approximates the sdpcone using the fact that every diagonal entry must be non-negative, so it adds the LP-cut
 *  \f \sum_{j = 1}^m (A_j)_{kk} y_j - (A_0)_{kk} \geq 0 \quad \forall k \leq n \f */
static
SCIP_RETCODE diagGEzero(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraitns*/
   int               nconss,     /**<number of constraints*/
   int*              naddconss   /**<pointer to store how many constraints were added*/
   )
{
   int blocksize;
   SCIP_CONSDATA* consdata;
   int nvars;
   int i;
   int j;
   int k;
   int block;
   SCIP_Real rhs = SCIPinfinity(scip);

   for (block = 0; block < nconss; ++block)
   {
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      consdata = SCIPconsGetData(conss[block]);

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      rhs = SCIPinfinity(scip);

      SCIP_Real* matrix;
      SCIP_CALL(SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1)));
      SCIP_CALL(SCIPgetLowerTriangConstMatrix(scip, conss[block], matrix));

      SCIP_Real* cons_array;
      SCIP_Real* lhs_array;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cons_array, blocksize * nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &lhs_array, blocksize) );

      /* the lhs is the (k,k)-entry of the constant matrix */
      for (k = 0; k < blocksize; ++k)
      {
         lhs_array[k] = matrix[compLowerTriangPos(k,k)];
      }

      /* get the (k,k)-entry of every matrix A_j */
      for (j = 0; j < nvars; ++j)
      {
         for (k = 0; k < blocksize; ++k)
         {
            /* initialize these with 0 */
            cons_array[k * nvars + j] = 0.0;
         }

         /* go through the nonzeroes of A_j and look for diagonal entries */
         for (i = 0; i < consdata->nvarnonz[j]; i++)
         {
            if (consdata->col[j][i] == consdata->row[j][i])
               cons_array[consdata->col[j][i] * nvars + j] = consdata->val[j][i];
         }
      }

      /* get the corresponding SCIP variables */
      SCIP_VAR** vars;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vars, nvars));

      for (int j = 0; j < nvars; ++j)
      {
         vars[j] = SdpVarmapperGetSCIPvar(&(consdata->varmapper), j);
      }

      /* add the LP-cuts to SCIP */
      for (int k = 0; k < blocksize; ++k)
      {
         SCIP_CONS* cons;

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "cl", consdata->nvars, vars, cons_array + k * consdata->nvars, lhs_array[k], rhs,
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         ++(*naddconss);
      }

      SCIPfreeBlockMemoryArray(scip, &vars, nvars );
      SCIPfreeBufferArray(scip, &matrix);
      SCIPfreeBlockMemoryArray(scip, &cons_array, blocksize * nvars);
      SCIPfreeBlockMemoryArray(scip, &lhs_array, blocksize);
   }
   return SCIP_OKAY;
}

/** presolve-routine that adds some constraints for approximation of the sdpcone, if there is an entry (i,j) in the constant
 *  matrix, than some variable k with \f (A_k)_{ii} \f needs to be >0 because sdp-matrices are diagonal-dominant (and the same for j)
 *
 * not clear if this is really true for all SDPs, probably only works if A_i and A_0 are all semidefinite (or at least
 * have positive diagonal entries) and all variables appearing in the SDP constraint are integer, then sum_{A_i_kk >0}
 * 1*y_i >= 1 is feasible, because it means that (sum A_i)_kk > 0 because all diagonal entries are positive (they can't
 * cancel each other) and at least one variable needs to be >=1 becaue this is equal to >0 for integers
 */
static
SCIP_RETCODE diagDominant(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraints*/
   int               nconss,     /**<number of constraints*/
   int*              naddconss   /**<pointer to store how many constraints were added*/
   )
{
   SCIP_Bool* nonzerorows; /* entry i will be 1 if there is an entry \f (A_0)_ij \f for some \f j \neq i \f */
   int blocksize;
   int i;
   int j;
   int nvars;
   int var;
   int endindex;
   SCIP_CONS* cons;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );
   assert ( naddconss != NULL );

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);
      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      SCIP_CALL(SCIPallocBufferArray(scip, &nonzerorows, blocksize));

      /* initialize nonzerorows with FALSE */
      for (j = 0; j < blocksize; j++)
      {
         nonzerorows[j] = FALSE;
      }

      /* iterate over all nonzeroes of the constant matrix and set all corresponding rows/cols to true */
      for (j = 0; j < consdata->constnnonz; j++)
      {
         if (consdata->constcol[j] != consdata->constrow[j])
         {
            assert ( ! ( SCIPisEQ(consdata->constval[j], 0.0) ) );
            nonzerorows[consdata->constcol[j]] = TRUE;
            nonzerorows[consdata->constrow[j]] = TRUE;
         }
      }

      /* diagvars[i] is an array with all variables with a diagol entry (i,i) in the corresponding matrix, if nonzerorows[i] is true or NULL otherwise
       * the outer array goes over all rows to ease the access, but only for those that are really needed memory will be allocated
       */
      int** diagvars;
      SCIP_CALL(SCIPallocBufferArray(scip, &diagvars, blocksize));
      int* ndiagvars;
      SCIP_CALL(SCIPallocBufferArray(scip, &ndiagvars, blocksize));
      for (j = 0; j < blocksize; ++j)
      {
         ndiagvars[j] = 0;
         if (nonzerorows[j])
            SCIP_CALL(SCIPallocBufferArray(scip, &diagvars[j], nvars));
      }

      /* find all variables with corresponding diagonal entries */
      for (var = 0; var < nvars; var++)
      {
         endindex = (var == nvars - 1) ? consdata->nnonz : consdata->begvar[var + 1];
         for (j = 0; j < nvarnonz[var]; j++)
         {
            if (consdata->col[var][j] == consdata->row[var][j])
            {
               diagvars[consdata->col[var][j]][ndiagvars[consdata->col[var][j]]] = var;
               ndiagvars[consdata->col[var][j]]++;
            }
         }
      }

      SCIP_VAR** vars;
      SCIP_Real* vals;

      for (j = 0; j < blocksize; ++j)
      {
         if (nonzerorows[j])
         {
            SCIP_CALL(SCIPallocBufferArray(scip, &vals, ndiagvars[j]));
            SCIP_CALL(SCIPallocBufferArray(scip, &vars, ndiagvars[j]));

            /* get the corresponding SCIP variables and set all coefficients to 1 */
            for (var = 0; var < ndiagvars[j]; ++var)
            {
               vars[var] = SdpVarmapperGetSCIPvar(&(consdata->varmapper), diagvars[j][var]);
               vals[var] = 1.0;
            }

            /* add the linear constraint sum_j 1.0 * diagvars[j] >= 1.0 */
            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "sum_diag_geq_1", ndiagvars[j], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
            (*naddconss)++;

            SCIPfreeBufferArray(scip, &vars);
            SCIPfreeBufferArray(scip, &vals);
            SCIPfreeBufferArray(scip, &diagvars[j]);
         }
      }

      SCIPfreeBufferArray(scip, &diagvars);
      SCIPfreeBufferArray(scip, &ndiagvars);
      SCIPfreeBufferArray(scip, &nonzerorows);
   }
   return SCIP_OKAY;
}

/** detects if there are blocks with size one and transfers it to a lp-row */
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraints to check*/
   int               nconss,     /**<number of constraints to check*/
   int*              naddconss,  /**<pointer to store how many constraints were added*/
   int*              ndelconss   /**<pointer to store how many constraints were deleted*/
   )
{
   SCIP_CONSHDLR* hdlr;
   int nnonz;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   int nvars;
   int i;
   int j;
   SCIP_Real rhs;
   int constnnonz;
   int count;
   int var;
   int endindex;

   for (i = 0; i < nconss; ++i)
   {
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);

      /* if there is a 1x1 SDP-Block */
      if (consdata->blocksize == 1)
      {
         nvars = consdata->nvars;
         nnonz = consdata->nnonz;
         constnnonz = consdata->constnnonz;
         SCIP_CALL(SCIPallocBufferArray(scip, &vars, nvars));
         SCIP_CALL(SCIPallocBufferArray(scip, &coeffs, nnonz));


         // get all lhs-entries
         count = 0;

         for (var = 0; var < nvars; var++)
         {
            for (j = 0; j < consdata->nvarnonz[var]; j++)
            {
               assert ( consdata->col[var][j] == 0 && consdata->row[var][j] == 0 ); /* if the block is size one, all entries should have row and col equal to 0 */
               coeffs[count] = consdata->val[var][j];
               vars[count] = SdpVarmapperGetSCIPvar(&(consdata->varmapper), var);
               count++;
            }
         }

         //get rhs
         assert ( consdata->constnnonz <= 1 ); /* the 1x1 constant matrix may only have one entry */

         rhs = (consdata->constnnonz == 1) ? consdata->constval[0] : 0.0; /* if this one entry is not 0, than this is the rhs, otherwise it's 0 */

         //add new linear cons
         SCIP_CONS* cons;

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "1x1", consdata->nvars, vars, coeffs, rhs, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL(SCIPaddCons(scip, cons));
         SCIP_CALL(SCIPreleaseCons(scip, &cons));

         (*naddconss)++;

         //delete old 1x1 sdpcone
         SCIP_CALL(SCIPdelCons(scip, conss[i]));
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &coeffs);
      }
   }
   return SCIP_OKAY;
}

/** presolve routine that looks through the data and eliminates fixed variables */
static
SCIP_RETCODE fixVars(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_CONS**       conss,   /**< array with constraints to check */
   int               nconss,  /**< number of constraints to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int block;
   int var;
   int i;
   int ind;
   int endindex;
   int nfixedvars;
   int nfixednonz;
   int newindex;
   int nvarnonz;
   int* savedcol;
   int* savedrow;
   int* savedval;
   int naddednonz;
   int nleftshifts;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );

   /* allocate memory to save nonzeros that need to be fixed */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedcol, consdata->nnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedrow, consdata->nnonz));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedval, consdata->nnonz));

   for (block = 0; i < nconss; ++i)
   {
      conshdlr = SCIPconsGetHdlr(conss[block]);
      assert( conshdlr != NULL );;

      assert ( strcmp(SCIPconshdlrGetName(conshdlr), "SDP") == 0);

      consdata = SCIPconsGetData(conss[block]);

      /* initialize this with zero for each block */
      nfixedvars = 0;
      nfixednonz = 0;

      for (var = 0; var < consdata->nvars; var++)
      {
         /* check if the variable is fixed in SCIP */
         if (SCIPvarGetStatus(SCIPvarGetProbvar(SdpVarmapper(&(consdata->varmapper), var))) == SCIP_VARSTATUS_FIXED)
         {
            assert ( SCIPisEQ(scip, SCIPvarGetLbGlobal(SCIPvarGetProbvar(SdpVarmapper(&(consdata->varmapper), var))),
                                    SCIPvarGetUbGlobal(SCIPvarGetProbvar(SdpVarmapper(&(consdata->varmapper), var)))) );

            for (i = 0; i < consdata->nvarnonz[i]; i++)
            {
               savedcol[nfixednonz] = consdata->col[i];
               savedrow[nfixednonz] = consdata->row[i];
               savedval[nfixednonz] = -consdata->val[i] * SCIPvarGetLbGlobal(SCIPvarGetProbvar(SdpVarmapper(&(consdata->varmapper), var)));
               /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */;
            }
            nfixedvars++;
            nfixednonz++;
         }
         /* if it isn't fixed shift the row/col/val/nvarnonz entries */
         else if (nfixedvars > 0)
         {
            consdata->col[var - nfixedvars] = consdata->col[var];
            consdata->row[var - nfixedvars] = consdata->row[var];
            consdata->val[var - nfixedvars] = consdata->val[var];
            consdata->nvarnonz[var - nfixedvars] = consdata->nvarnonz[var];
         }
      }

      /* free memory for the sdp-arrays equal to the number of fixed variables */
      SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &(consdata->col), consdata->nvars, consdata->nvars - nfixedvars));
      SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &(consdata->row), consdata->nvars, consdata->nvars - nfixedvars));
      SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &(consdata->val), consdata->nvars, consdata->nvars - nfixedvars));
      SCIP_CALL(SCIPreallocBlockMemoryArray(scip, &(consdata->nvarnonz), consdata->nvars, consdata->nvars - nfixedvars));

      /* insert the fixed variables into the constant arrays */
      SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, nfixednonz, FALSE,  1.0, consdata->constrow, consdata->constcol,
                                       consdata->constval, &(consdata->constnnonz)) );


      /* free the saved arrays */
      SCIPfreeBlockMemoryArray(scip, &savedcol, consdata->nnonz);
      SCIPfreeBlockMemoryArray(scip, &savedrow, consdata->nnonz);
      SCIPfreeBlockMemoryArray(scip, &savedval, consdata->nnonz);

      consdata->nvars -= nfixedvars;

      /* remove all fixed variables from the varmapper */
      for (var = 0; var < consdata->nvars; var++)
      {
         /* check if the variable is fixed in SCIP */
         if (SCIPvarGetStatus(SCIPvarGetProbvar(SdpVarmapper(&(consdata->varmapper), var))) == SCIP_VARSTATUS_FIXED)
            SCIP_CALL( SdpVarmapperRemoveSdpIndex(scip, consdata->varmapper, var) );
      }
   }

   return SCIP_OKAY;
}

/** presolve routine that looks through the data and handles aggregated and multiaggregated variables
 */
static
SCIP_RETCODE multiaggrVars(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_CONS**       conss,   /**< array with constraints to check */
   int               nconss,  /**< number of constraints to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   SdpVarmapper* varmapper;
   SdpVarmapper* globalvarmapper;
   int block;
   int var;
   int i;
   int j;
   int endindex;
   int naddednonz;
   int newind;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   int ind;
   SCIP_VAR** aggrvars;
   SCIP_Real* scalars;
   int naggrvars;
   SCIP_Real constant;
   int requiredsize;
   int nvars;
   int globalnvars;
   int aggrind;
   int aggrconsind;
   int aggrglobalind;
   int* coltoadd;
   int* rowtoadd;
   SCIP_Real* valtoadd;
   int naggrnonz;
   int nremovedvars; /* this is a net value, if it is negative, more were added than removed */
   int nleftshifts;
   int posleftshiftsdone;  /* all variables after this position still need to be move nfestshifts positions to the left */
   int naggrnonz;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );

   for (block = 0; block < nconss; ++block)
   {
      conshdlr = SCIPconsGetHdlr(conss[block]);
      assert( conshdlr != NULL );

      assert( strcmp(SCIPconshdlrGetName(conshdlr), "SDP") == 0 );

      consdata = SCIPconsGetData(conss[block]);
      varmapper = &(consdata->varmapper); /* gets indices of variables in current block */
      globalvarmapper = &(SCIPconshdlrGetData(conshdlr)->varmapper); /* gets indices of variables in the whole sdp */
      globalnvars = SCIPgetNVars(scip);
      nvars = consdata->nvars;
      nremovedvars = 0;
      nleftshifts = 0;

      for (var = 0; var < nvars; var++)
      {
         if (SCIPvarGetStatus(SdpVarmapper(varmapper, var)) == SCIP_VARSTATUS_AGGREGATED ||
             SCIPvarGetStatus(SdpVarmapper(varmapper, var)) == SCIP_VARSTATUS_MULTAGGR)
         {
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &aggrvars, globalnvars));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &scalars, globalnvars));

            aggrvars[0] = consdata->vars[var];
            naggrvars = 1;

            /* get the variables this var was aggregated to */
            SCIP_CALL(SCIPgetProbvarLinearSum(scip, aggrvars, scalars, &naggrvars, globalnvars, &constant, &requiredsize, TRUE));
            assert( requiredsize <= globalnvars );

            /* save the nonzeroes of the (multi)aggregated var */
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedcol, consdata->nvarnonz[var]));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedrow, consdata->nvarnonz[var]));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedval, consdata->nvarnonz[var]));

            naggrnonz = 0;

            for (i = 0; i < consdata->nvarnonz[var]; i++)
            {
               savedcol[naggrnonz] = consdata->col[i];
               savedrow[naggrnonz] = consdata->row[i];
               savedval[naggrnonz] = consdata->val[i];
               naggrnonz++;
            }

            assert ( naggrnonz == consdata->nvarnonz[var] );

            /* sort them by nondecreasing row and then col to make the search for already existing entries easier (this is done here, because it
             * only needs to be done once and not for each variable this is multiaggregated to) */
            SdpVarfixerSortRowsCols(savedrow, savedcol, savedval, naggrnonz);

            nremovedvars++;
            nleftshifts++;
            /* remove the variable from the varmapper */
            SCIP_CALL( SdpVarmapperRemoveSdpIndex(scip, consdata->varmapper, var) );

            /* iterate over all variables this was aggregated to and insert the corresponding nonzeroes */
            for (aggrind = 0; aggrind < naggrvars; aggrind++)
            {

               /* check if the variable allready exists in this block */
               if (SdpVarmapperExistsSCIPvar(varmapper, aggrvars[aggrind]))
               {
                  /* get the variable index in the current constraint */
                  aggrconsind = SdpVarmapperGetSdpIndex(varmapper, aggrvars[aggrind]);

                  /* if the index is bigger than posleftshiftsdone we have to add nleftshifts, because the array-entries weren't yet shifted (but the
                   * varmapper was already updated) */
                  if (aggrconsind > posleftshiftsdone)
                     aggrconsind += nleftshifts;

                  SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, naggrnonz, TRUE, scalars[aggrind],
                             consdata->row[aggrconsind], consdata->col[aggrconsind], consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind])) );
               }
               else
               {
                  /* the variable has to be added to this constraint */

                  /* find the right position to insert this variable, to not have to rearrange later when combining with other blocks */
                  aggrglobalind = SdpVarmapperGetSdpIndex(globalvarmapper, aggrvars[aggrind]);
                  aggrconsind = 0;

                  /* find the right index, such that all variables before it have a lower global SDP-index and all behind have a higher one */
                  while (aggrconsind < consdata->nvars - nremovedvars &&
                          SdpVarmapperGetSdpIndex(globalvarmapper, SdpVarmapperGetSCIPvar(varmapper, aggrconsind)) < aggrglobalind)
                     aggrconsind++;

                  assert (SdpVarmapperGetSdpIndex(globalvarmapper, SdpVarmapperGetSCIPvar(varmapper, aggrconsind)) > aggrglobalind); /* = would be the if-part*/

                  SdpVarmapperInsertVar(scip, varmapper, aggrvars[aggrind], aggrconsind);

                  /* check if there is an empty spot in the sdp-arrays to insert the new variable */
                  if (nleftshifts > 0)
                  {
                     if (aggrconsind <= posleftshiftsdone + 1)
                     {
                        /* if the position where the new variable should be inserted is located ahead of the empty spots, all in between have to be
                         * shifted one position to the right, also if it should be inserted immediatly after last variable ahead of the empty spots,
                         * we may immediatly insert it at the first position of the empty spots */
                        for (i = posleftshiftsdone; i >= aggrconsind; i--)
                        {
                           consdata->row[i + 1] = consdata->row[i];
                           consdata->col[i + 1] = consdata->col[i];
                           consdata->val[i + 1] = consdata->val[i];
                           consdata->nvarnonz[i + 1] = consdata->nvarnonz[i];
                        }

                        /* allocate memory for the nonzeros */
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->row[aggrconsind], naggrnonz) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->col[aggrconsind], naggrnonz) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->val[aggrconsind], naggrnonz) );

                        /* add the new entries and nonzeros to the arrays */
                        consdata->nvarnonz[aggrconsind] = naggrnonz;
                        for (i = 0; i < naggrnonz; i++)
                        {
                           consdata->row[aggrconsind][i] = savedrow[i];
                           consdata->col[aggrconsind][i] = savedcol[i];
                           consdata->val[aggrconsind][i] = savedval[i];
                        }
                     }
                     else
                     {
                        /* it should be inserted after at least one variable that hasn't been shifted yet, so we have to move all of these at least
                         * one position, if we have to move them already, we might as well move them to their right position in the hope that we
                         * don't have to move them again (if we do have to move them again, we have to touch them again either way) */
                        for (i = nposleftshiftsdone; i < aggrconsind; i++)
                        {
                           consdata->row[i] = consdata->row[i + nleftshifts];
                           consdata->col[i] = consdata->col[i + nleftshifts];
                           consdata->val[i] = consdata->val[i + nleftshifts];
                           consdata->nvarnonz[i] = consdata->nvarnonz[i + nleftshifts];
                        }

                        /* remember, that we have moved these */
                        nposleftshiftsdone = aggrconsind;

                        /* allocate memory for the nonzeros */
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->row[aggrconsind], naggrnonz) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->col[aggrconsind], naggrnonz) );
                        SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->val[aggrconsind], naggrnonz) );

                        /* add the new entries and nonzeros to the arrays */
                        consdata->nvarnonz[aggrconsind] = naggrnonz;
                        for (i = 0; i < naggrnonz; i++)
                        {
                           consdata->row[aggrconsind][i] = savedrow[i];
                           consdata->col[aggrconsind][i] = savedcol[i];
                           consdata->val[aggrconsind][i] = savedval[i];
                        }
                     }
                     nleftshifts--;
                  }
                  else
                  {
                     /* we have to enlarge the arrays and move all in between one position to the right */
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->nvarnonz), consdata->nvars, consdata->nvars + 1));
                     //TODO: could go for length of globalvarmapper to only have to do this once (then we need to save this e.g. as a bool),
                     //but that also increases the likelyhood of actually having to move the data
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col), consdata->nvars, consdata->nvars + 1));
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row), consdata->nvars, consdata->nvars + 1));
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val), consdata->nvars, consdata->nvars + 1));

                     /* allocate memory for the nonzeros */
                     SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->row[aggrconsind], naggrnonz) );
                     SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->col[aggrconsind], naggrnonz) );
                     SCIP_CALL( SCIPallocBlockMemoryArray(scip, consdata->val[aggrconsind], naggrnonz) );

                     /* add the new entries and nonzeros to the arrays */
                     consdata->nvarnonz[aggrconsind] = naggrnonz;
                     for (i = 0; i < naggrnonz; i++)
                     {
                        consdata->row[aggrconsind][i] = savedrow[i];
                        consdata->col[aggrconsind][i] = savedcol[i];
                        consdata->val[aggrconsind][i] = savedval[i];
                     }
                  }
                  nremovedvars--;
               }
            }

            SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, naggrnonz, TRUE, constant,
                                              consdata->constrow[aggrconsind], consdata->constcol[aggrconsind],
                                              consdata->constval[aggrconsind], &(consdata->constnnonz)) );

            /* free all arrays that are no longer needed */
            SCIPfreeBlockMemoryArray(scip, &savedcol, nvarnonz);
            SCIPfreeBlockMemoryArray(scip, &savedrow, nvarnonz);
            SCIPfreeBlockMemoryArray(scip, &savedval, nvarnonz);
            SCIPfreeBlockMemoryArray(scip, &aggrvars, globalnvars);
            SCIPfreeBlockMemoryArray(scip, &scalars, globalnvars);
         }
         else if (nleftshifts > 0 && var > posleftshiftsdone)
         {
            /* move the entries in all sdp-arrays */
            consdata->col[var - nfixedvars] = consdata->col[var];
            consdata->row[var - nfixedvars] = consdata->row[var];
            consdata->val[var - nfixedvars] = consdata->val[var];
            consdata->nvarnonz[var - nfixedvars] = consdata->nvarnonz[var];
         }
      }

      consdata->nvars = consdata->nvars - nfixedvars;

      /* recompute sdpnnonz */
      consdata->sdpnnonz = 0;
      for (var = 0; var < consdata->nvars; var++)
      {
         consdata->sdpnnonz += consdata->nvarnonz[var];
      }
   }

   return SCIP_OKAY;
}


/** after the problem is transformed swap all variables in this constraint for the transformed ones */
static
SCIP_DECL_CONSINIT(consInitSdp)
{
   /* get transformed variables, if we are in the transformed problem */
   SCIP_CONSDATA* consdata;
   int i;
   int k;

   assert ( SCIPisTransformed(scip) );

   for (i = 0; i < nconss; ++i)
   {
      /* ?????? turn this code into an assert ??????? */
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      consdata = SCIPconsGetData(conss[i]);

      SdpVarmapperTransform(scip, &(consdata->varmapper));
      SdpVarmapperTransform(scip, &(SCIPconshdlrGetData(conshdlr)->varmapper));
   }
   return SCIP_OKAY;
}

/** locks a variable up if the corresponding constraint matrix is not positive semidefinite, locks it down if it is not negative semidefinite */
static
SCIP_DECL_CONSLOCK(consLockSdp)
{
   SCIP_Real* Aj;
   SCIP_CONSDATA* consdata;
   SdpVarmapper* varmapper;
   int blocksize;
   int var;
   int nvars;
   SCIP_Real eigenvalue;

   consdata = SCIPconsGetData(conss[block]);
   varmapper = &(consdata->varmapper);
   blocksize = consdata->blocksize;
   nvars = consdata->nvars;

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &Aj, blocksize * blocksize));

   for (var = 0; var < nvars; v++)
   {
      SCIP_CALL(SCIPconsSdpGetFullAj(scip, cons, var, Aj));

      /* compute the smallest eigenvalue */
      SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
      {
         /* as the lowest eigenvalue is negative, the matrix is not positive semidefinite, so adding more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, SdpVarmapperGetSCIPvar(varmapper, var), nlocksneg, nlockspos) );
      }

      /* compute the biggest eigenvalue */
      SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
      if ( SCIPisPositive(scip, eigenvalue) )
      {
         /* as the biggest eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, SdpVarmapperGetSCIPvar(varmapper, var), nlockspos, nlocksneg) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &constmatrix, blocksize * blocksize);

   return SCIP_OKAY;
}


/** after presolving the Constraint Handler Data with the varmapper is initialized */
static
SCIP_DECL_CONSEXITPRE(consExitpreSdp)
{
   int nvars;
   SCIP_Var** vars;
   SCIP_ConshdlrData conshdlrdata;

   assert ( scip != NULL );
   assert ( conshdlr != NULL );

   nvars = SCIPgetNVars(scip);
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vars, nvars));
   vars = SCIPgetVars(scip);

   SCIP_CALL(SCIPallocBlockMemory(scip, &conshdlrdata));

   SdpVarmapperCreate(scip, conshdlrdata->varmapper);
   SdpVarmapperAddVars(scip, conshdlrdata->varmapper, nvars, vars);

   SCIPconshdlrSetData(conshdlr, conshdlrdata);

   SCIP_CALL(fixVars(scip, conss, nconss));
   SCIP_CALL(multiaggrVars(scip, conss, nconss));

   SCIPfreeBlockMemoryArray(scip, &vars, nvars);

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSdp)
{
   assert( result != 0 );

   if ( nrounds == 0 )
   {
      SCIP_CALL( diagGEzero(scip, conss, nconss, naddconss) );
   }

   SCIP_CALL( move_1x1_blocks_to_lp(scip, conss, nconss, naddconss, ndelconss) );

   SCIP_CALL( fixVars(scip, conss, nconss) );

   if ( nrounds == 0 )
   {
      SCIP_CALL( diagDominant(scip, conss, nconss, naddconss) );  /*TODO: could be activated for some problem classes
      but doesn't seem to work in the general case */
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** creates transformed constraint */
static
SCIP_DECL_CONSTRANS(consTransSdp)
{
   SCIP_CONSDATA* sourcedata = NULL;
   SCIP_CONSDATA* targetdata = NULL;
   SdpVarmapper targetmapper;
   int i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocMemory(scip, &targetdata) );
   // TODO: copy complete arrays instead of pointers
   targetdata->nvars = sourcedata->nvars;
   targetdata->nnonz = sourcedata->nnonz;
   targetdata->blocksize = sourcedata->blocksize;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->begvar), sourcedata->nvars));
   for (i = 0; i < sourcedata->nvars; i++)
      targetdata->begvar[i] = sourcedata->begvar[i];

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->col), sourcedata->nnonz));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->row), sourcedata->nnonz));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->val), sourcedata->nnonz));
   for (i = 0; i < sourcedata->nnonz; i++)
   {
      targetdata->col[i] = sourcedata->col[i];
      targetdata->row[i] = sourcedata->row[i];
      targetdata->val[i] = sourcedata->val[i];
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &newmapper));
   SCIP_CALL(SdpVarmapperClone(scip, &(sourcedata->varmapper), &newmapper));

   targetdata->constnnonz = sourcedata->constnnonz;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->constcol), sourcedata->constnnonz));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->constrow), sourcedata->constnnonz));
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->constval), sourcedata->constnnonz));
   for (i = 0; i < sourcedata->nnonz; i++)
   {
      targetdata->constcol[i] = sourcedata->constcol[i];
      targetdata->constrow[i] = sourcedata->constrow[i];
      targetdata->constval[i] = sourcedata->constval[i];
   }
   assert( targetdata != NULL );


   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/**checks feasiblity of constraint, e.g. the positive semidefiniteness*/
static
SCIP_DECL_CONSCHECK(consCheckSdp)
{

   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CALL( checkSdpCon(scip, conss[i], sol, checkintegrality, checklprows, printreason, result) );
      if (*result == SCIP_INFEASIBLE)
      {
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/** enforce pseudo solution method, returns didnotrun, if objinfeasible, computes cut otherwise*/
static
SCIP_DECL_CONSENFOPS(consEnfopsSdp)
{
   SCIP_CONSDATA* consdata;
   int nvars;
   SCIP_Real* coeff;
   SCIP_Real lhs;

   if (objinfeasible)
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CALL( checkSdpCon(scip, conss[i], NULL, 0, 0, 0, result) );

      lhs = 0.0;
      consdata = SCIPconsGetData(conss[i]);
      nvars = consdata->nvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );
      SCIP_CALL(cutUsingEigenvector(scip, conss[i], NULL, coeff, &lhs));

      if (*result != SCIP_INFEASIBLE)
      {
         lhs = floor(lhs);
      }

      SCIPfreeBufferArray(scip, &coeff);

   }

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions
 *
 *  enforce lp solution method, must be implemented, but there is no lp in the sdp-case, so returns SCIP_ERROR.
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSdp)
{
   SCIP_CONSDATA* consdata;
   bool all_feasible = TRUE;
   bool separated = FALSE;
   for (int i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      SCIP_CALL( checkSdpCon(scip, conss[i], NULL, 0, 0, 0, result) );
      if (*result == SCIP_FEASIBLE)
      {
         continue;
      }
      all_feasible = FALSE;

      int nvars = consdata->nvars;


      SCIP_Real lhs = 0.0;
      SCIP_Real* coeff = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );

      SCIP_CALL(cutUsingEigenvector(scip, conss[i], NULL, coeff, &lhs));

      if (*result == !SCIP_INFEASIBLE)
      {
         lhs = floor(lhs);
      }
      SCIP_ROW* row;
      SCIP_Real rhs = SCIPinfinity(scip); //local modifiable, removable

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, "eigenvectorcut_enfolp", lhs, rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for (int j = 0; j < nvars; ++j)
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, SdpVarmapperGetSCIPvar(consdata->varmapper, j), coeff[j]) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      SCIP_Bool infeasible;
      SCIP_CALL(SCIPaddCut(scip, NULL, row, FALSE, &infeasible));

      if ( infeasible )
         *result = SCIP_CUTOFF;
      else
      {
         SCIP_CALL(SCIPaddPoolCut(scip, row));

         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );
         *result = SCIP_SEPARATED;
         separated = TRUE;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      SCIPfreeBufferArray(scip, &coeff);
   }
   if (all_feasible)
   {
      return SCIP_OKAY;
   }
   if (separated)
   {
      *result = SCIP_SEPARATED;
   }


   SCIP_VAR** vars;
   vars = SCIPgetVars(scip);
   int count = 0;
   for (int i = 0; i < SCIPgetNVars(scip); ++i)
   {
      if( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) && SCIPvarIsIntegral(vars[i]))
      {

         SCIP_CALL( SCIPaddExternBranchCand(scip, vars[i], 10000, SCIP_INVALID) );
         count++;
      }
   }

   return SCIP_OKAY;
}


/** separates a solution using constraint specific ideas, gives cuts to scip */
static
SCIP_DECL_CONSSEPASOL(consSepasolSdp)
{
   *result = SCIP_DIDNOTFIND;
   for (int i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL(separateSol(scip, conshdlr, conss[i], sol, result));
   }

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSdp)
{
   *result = SCIP_DIDNOTFIND;
   for (int i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL(separateSol(scip, conshdlr, conss[i], NULL, result));
   }

   return SCIP_OKAY;

}


/** delete method of sdp constrainthandler */
static
SCIP_DECL_CONSDELETE(consDeleteSdp)
{
   int i;

   assert(consdata != NULL);

   for (i = 0; i < consdata->nvars; i++)
   {
      SCIPfreeBlockMemoryArray(scip, &consdata->col[i], consdata->nvarnonz[i]);
      SCIPfreeBlockMemoryArray(scip, &consdata->row[i], consdata->nvarnonz[i]);
      SCIPfreeBlockMemoryArray(scip, &consdata->val[i], consdata->nvarnonz[i]);
   }
   SCIPfreeBlockMemoryArray(scip, &consdata->col, consdata->nvars);
   SCIPfreeBlockMemoryArray(scip, &consdata->row, consdata->nvars);
   SCIPfreeBlockMemoryArray(scip, &consdata->val, consdata->nvars);
   SCIPfreeBlockMemoryArray(scip, &consdata->nvarnonz, consdata->nvars);
   SCIPfreeBlockMemoryArray(scip, &consdata->constcol, consdata->nnonz);
   SCIPfreeBlockMemoryArray(scip, &consdata->constrow, consdata->nnonz);
   SCIPfreeBlockMemoryArray(scip, &consdata->constval, consdata->nnonz);
   SdpVarmapperFree(scip, consdata->varmapper);
   SCIPfreeMemory(scip, consdata);
   return SCIP_OKAY;
}

/** creates the handler for sdp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   assert( scip != 0 );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSdp, consEnfopsSdp, consCheckSdp, consLockSdp, 0) );

   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSdp) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitSdp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );

   return SCIP_OKAY;
}


/** get the data belonging to a single SDP-constraint
 *  in arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient for the data that
 *  needs to be inserted, a debug message will be thrown and this variable will be set to the needed length */
EXTERN
SCIP_RETCODE SCIPconsSdpGetData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nvars,              /**< number of variables in this SDP constraint */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  blocksize,          /**< size of this SDP-block */
   int*                  arraylength,        /**< length of the given nvarnonz, col, row and val arrays, if this is too short this will return the needed length*/
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val are pointing to */
   int**                 col,                /**< pointers to column indices of the nonzeroes for each variable */
   int**                 row,                /**< pointers to row indices of the nonzeroes for each variable */
   SCIP_Real**           val,                /**< pointers to values of the nonzeroes for each variable */
   SdpVarmapper*         varmapper,          /**< varmapper mapping the variables to positions in the sdp-arrays and vice versa */
   int*                  constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint */
   int**                 constcol,           /**< pointer to column indices of the constant nonzeroes */
   int**                 constrow,           /**< pointer to row indices of the constant nonzeroes */
   SCIP_Real**           constval,           /**< pointer to values of the constant nonzeroes */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   char* name

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( nvars != NULL );
   assert ( nnonz != NULL );
   assert ( blocksize != NULL );
   assert ( nvarnonzlength != NULL );
   assert ( nvarnonz != NULL );
   assert ( arraylength != NULL );
   assert ( col != NULL );
   assert ( row != NULL );
   assert ( val != NULL );
   assert ( varmapper != NULL );
   assert ( constnnonz != NULL );
   assert ( constarraylength != NULL );
   assert ( constcol != NULL );
   assert ( constrow != NULL );
   assert ( constval != NULL );

   consdata = SCIPconsGetData(cons);
   name = SCIPconsGetName(cons);

   *nvars = consdata->nvars;
   *nnonz = consdata->nnonz;
   *constnnonz = consdata->constnnonz;
   *blocksize = consdata->blocksize;

   varmapper = &(consdata->varmapper);

   /* check that the sdp-arrays are long enough to store the information */
   if (arraylength < consdata->nvars)
   {
      SCIPdebugMessage("nvarnonz, col, row and val arrays were not long enough to store the information for cons %s, they need to be at least size %d! \n", name, consdata->nnonz);
      *arraylength = consdata->nvars;
   }
   else
   {
      for (i=0; i < nnonz; i++)
      {
         nvarnonz[i] = consdata->nvarnonz[i];
         col[i] = &(consdata->col[i]);
         row[i] = &(consdata->row[i]);
         val[i] = &(consdata->val[i]);
      }
   }

   /* set the constant pointers (if a constant part exists) */
   if (consdata->constnnonz > 0)
   {
      *constcol = &(consdata->constcol[0]);
      *constrow = &(consdata->constrow[0]);
      *consdtval = &(consdata->constval[0]);
   }

   return SCIP_OKAY;
}

/** gets the number of nonzeroes and constant nonzeroes for this SDP constraint */
EXTERN
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( nnonz != NULL );
   assert ( constnnonz != NULL );

   consdata = SCIPconsGetData(cons);
   *nnonz = consdata->nnonz;
   *constnnonz = constdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the full constraint Matrix \f A_j \f for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int                   j,                  /**< the variable j to get the corresponding matrix \f A_j \f for */
   SCIP_Real*            Aj,                 /**< pointer to store the full matrix \f A_j \f */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int blocksize;
   int endindex;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( j >= 0 );
   assert ( mat != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   assert ( j < blockize );

   for (int i = 0; i < consdata->nvarnonz[j]; i++)
   {
      Aj[consdata->col[j][i] * blocksize + consdata->row[j][i]] = consdata->val[j][i];
      Aj[consdata->row[j][i] * blocksize + consdata->col[j][i]] = consdata->val[j][i];
   }

   return SCIP_OKAY;
}

/** gives an n*n-long array with the full constant matrix */
SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   SCIP_Real*            mat                 /**< pointer to store the full constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int blocksize;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( mat != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   for (int i = 0; i < consdata->constnnonz; i++)
   {
      mat[consdata->constcol[i] * blocksize + consdata->constrow[i]] = consdata->constval[i];
      mat[consdata->constrow[i] * blocksize + consdata->constcol[i]] = consdata->constval[i];
   }

   return SCIP_OKAY;
}

/** gives an 0.5*n*(n+1)-long array with the lower triangular part of the constant matrix indexed by compLowerTriangPos */
SCIP_RETCODE SCIPconsSdpGetLowerTriangConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   SCIP_Real*            mat                 /**< pointer to store the lower triangular part of the constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( mat != NULL );

   consdata = SCIPconsGetData(cons);

   for (int i = 0; i < consdata->constnnonz; i++)
   {
      mat[compLowerTriangPos(consdata->constcol[i], consdata->constrow[i])] = consdata->constval[i];
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSFREE(consFreeSDP)
{
SCIP_CONSHDLRDATA* conshdlrdata;
conshdlrdata = SCIPconshdlrGetData(conshdlr);
assert(conshdlrdata != NULL);
SdpVarmapperFree(scip, &(conshdlrdata->varmapper));
SCIPfreeBlockMemory(scip, &conshdlrdata);
SCIPconshdlrSetData(conshdlr, NULL);

return SCIP_OKAY;
}

/**creates an sdp-constraint*/
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in this SDP constraint */
   int                   nnonz,              /**< number of nonzeroes in this SDP constraint */
   int                   blocksize,          /**< size of this SDP-block */
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val point to */
   int**                 col,                /**< pointer to column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to values of the nonzeroes for each variable */
   SCIP_Var**            vars,               /**< SCIP_Variables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeroes */
   int*                  constrow,           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< values of the constant nonzeroes */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   int j;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( name != NULL );
   assert ( nvars >= 0 );
   assert ( nnonz >= 0 );
   assert ( blocksize >= 0 );
   assert ( constnnonz >= 0 );
   assert ( nvars == 0 || vars != NULL );
   assert ( nnonz == 0 || (bevar != NULL && col != NULL && row != NULL && val != NULL ));
   assert ( constnnonz == 0 || (constcol != NULL && constrow != NULL && constval != NULL ));

   conshdlr = SCIPfindConshdlr(scip, "SDP");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, nnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, nnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, nnonz) );
   for (i = 0; i < nvars; i++)
   {
      assert ( nvarnonz[i] >= 0 );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i]));
   }

   consdata->nvars = nvars;
   consdata->nnonz = nnonz;
   consdata->constnnonz = constnnonz;
   consdata->blocksize = blocksize;

   for (i = 0; i < nvars; i++)
   {
      consdata->nvarnonz[i] = nvarnonz[i];

      for (j = 0; j < nvarnonz[i]; j++)
      {
         assert ( col[i][j] >= 0 );
         assert ( col[i][j] < blocksize );
         assert ( row[i][j] >= 0 );
         assert ( row[i][j] < blocksize );
         consdata->col[i][j] = col[i][j];
         consdata->row[i][j] = row[i][j];
         consdata->val[i][j] = val[i][j];
      }
   }

   for (i = 0; i < constnnonz; i++)
   {
      consdata->constcol[i] = constcol[i];
      consdata->constrow[i] = constrow[i];
      consdata->constval[i] = constval[i];
   }

   SdpVarmapperCreate(scip, consdata->varmapper);
   SdpVarmapperAddVars(scip, consdata->varmapper, nvars, vars);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );


   return SCIP_OKAY;
}
