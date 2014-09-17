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
   SCIP_Var**            vars;               /**< SCIP_Variables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz;         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol;           /**< column indices of the constant nonzeroes */
   int*                  constrow;           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval;           /**< values of the constant nonzeroes */
};

static int neigvalcuts = 0; /* this is used to give the eigenvalue-cuts distinguishable names */
static int neigveccuts = 0; /* this is used to give the eigenvector-cuts distinguishable names */
static int ndiaggezerocuts = 0; /* this is used to give the diagGEzero-cuts distinguishable names */
static int ndiagdomcuts = 0; /* this is used to give the diagDominant-cuts distinguishable names */
static int n1x1blocks = 0; /* this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */


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

   SCIPfreeBufferArray(scip, &WTMP);
   SCIPfreeBufferArray(scip, &IWORK);
   SCIPfreeBufferArray(scip, &WORK);

   return SCIP_OKAY;
}

#ifndef NDEBUG
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
#else
#define compLowerTriangPos(i, j) (i*(i+1)/2 + j)
#endif

/** takes an 0.5*n*(n+1) array of a symmetric matrix and expands it to a n*n array of the full matrix to input into LAPACK */
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

   ind = 0;

   /* traverse the lower triangular part in the order of the indices and copy the values to both lower and upper triangular part */
   for (i = 0; i < size; i++)
   {
      for (j = 0; j <= i; j++)
      {
         assert ( ind == compLowerTriangPos(i,j) );
         fullMat[i*size + j] = symMat[ind];
         fullMat[j*size + i] = symMat[ind];
         ind++;
      }
   }

   return SCIP_OKAY;
}

/** For a given vector \f y \f computes the SDP-Matrix \f \sum_{j=1}^m A_j y_j - A_0 \f for this SDP block.
 *
 *  The length of the matrix array needs to be (length of y) * (length of y + 1) /2, this will be indexed by compLowerTriangPos.
 */
static
SCIP_RETCODE computeSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             y,                  /**< solution to separate */
   SCIP_Real*            matrix              /**< the SDP-Matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int ind;
   int nvars;
   SCIP_Real yval;

   assert ( cons != NULL );
   assert ( y != NULL );
   assert ( matrix != NULL );

   consdata = SCIPconsGetData(cons);
   nvars = consdata->nvars;

   /* initialize the matrix with 0 */
   for (i = 0; i < 0.5 * nvars * (nvars + 1); i++)
      matrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      yval = SCIPgetSolVal(scip, y, consdata->vars[i]);
      for (ind = 0; i < consdata->nvarnonz[i]; i++)
         matrix[compLowerTriangPos(consdata->row[i][ind], consdata->col[i][ind])] += yval * consdata->val[i][ind];
   }

   /* substract the constant part */
   for (ind = 0; ind < consdata->constnnonz; ind++)
      matrix[compLowerTriangPos(consdata->constrow[ind], consdata->constcol[ind])] -= consdata->constval[ind];

   return SCIP_OKAY;
}

/** For a given variable-index j and a Vector v computes \f$ v^T A_j v \f$. */
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
   assert ( vAv != NULL );

   consdata = SCIPconsGetData(cons);

   assert ( j < consdata->nvars );

   /* initialize the product with 0 */
   *vAv = 0.0;

   for (i = 0; i < consdata->nvarnonz[i]; i++)
   {
      if (consdata->col[j][i] == consdata->row[j][i])
         *vAv += v[consdata->col[j][i]] * consdata->val[j][i] * v[consdata->row[j][i]];
      else
      {
         /* Multiply by 2, because the matrix is symmetric and there is one identical contribution each from lower and upper triangular part. */
         *vAv += 2.0 * v[consdata->col[j][i]] * consdata->val[j][i] * v[consdata->row[j][i]];
      }
   }

   return SCIP_OKAY;
}

/** Get the maximum absolute value of an entry of the constant matrix. */
static
SCIP_Real getMaxConstEntry(
   SCIP_CONS*            cons                /**< the SDP constraint that includes the Matrix \f A_j \f */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real max;
   int i;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* initialize max with the absolute value of the first entry of the constant matrix */
   max = REALABS(consdata->constval[0]);

   /* iterate over the remaining arrays, updating max if a higher absolute value is found */
   for (i = 1; i < consdata->constnnonz; i++)
   {
      if ( REALABS(consdata->constval[i]) > max)
         max = REALABS(consdata->constval[i]);
   }

   return max;
}


/** separate current solution using a cut, with the eigenvectors and -values of the solution matrix
 *
 *  This function computes the eigenvectors of the matrix, takes the first one and multiplies the matrix with it
 *  \f$ x^T*A_i*x = coeff[i], x^T*A_0*x =lhs \f$.
 */
static
SCIP_RETCODE cutUsingEigenvector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real*            coeff,              /**< coefficients of the computed cut */
   SCIP_Real*            lhs                 /**< lhs of the computed cut */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real* fullconstmatrix = NULL;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real* output_vector = NULL;
   int blocksize;
   int j;

   assert( cons != NULL );
   assert( lhs != NULL );

   *lhs = 0.0;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1) ));
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ));
   SCIP_CALL( SCIPallocBufferArray(scip, &fullconstmatrix, blocksize * blocksize));
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvector, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &output_vector, blocksize) );

   /* compute the matrix \f$ \sum_j A_j y_j \f$ */
   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix));

   /* expand it because LAPACK wants the full matrix instead of the lower triangular part */
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( computeIthEigenvalue(scip, TRUE, blocksize, matrix, 1, eigenvalues, eigenvector) );

   /* get full constant matrix */
   SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, fullconstmatrix) );

   /* multiply eigenvector with constant matrix to get lhs (after multiplying again with eigenvector from the left) */
   SCIP_CALL( Blas_DGEMV(blocksize, blocksize, 1.0, fullconstmatrix, eigenvector, 0.0, output_vector) );

   for (j = 0; j < blocksize; ++j)
      *lhs += eigenvector[j] * output_vector[j];

   /* compute \f$ v^T A_j v \f$ for eigenvector v and each matrix \f$ A_j \f$ to get the coefficients of the LP cut */
   for (j = 0; j < consdata->nvars; ++j)
      multiplyConstraintMatrix(cons, j, eigenvector, &coeff[j]);

   SCIPfreeBufferArray(scip, &output_vector);
   SCIPfreeBufferArray(scip, &eigenvector);
   SCIPfreeBufferArray(scip, &fullconstmatrix);
   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &eigenvalues);

   return SCIP_OKAY;
}

/** checks feasibility for a single SDP-Cone */
SCIP_RETCODE SCIPconsSdpCheckSdpCons(
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
   SCIP_Real* fullmatrix = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize) );

   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix) );
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, fullmatrix, 1, eigenvalues, NULL) );

   //we are going to use one of the dimacs error norms for checking feasiblity.
   //we use the second one: err=max{0, -lambda_min(x)/(1+maximumentry of rhs}

   check_value = (-eigenvalues[0]) / (1.0 + getMaxConstEntry(cons));

   if ( SCIPisLE(scip, check_value, 0.0) )
      *result = SCIP_FEASIBLE;
   else
   {
      *result = SCIP_INFEASIBLE;
      if ( printreason )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "non psd matrix (eigenvalue %f, dimacs error norm = %f).\n", eigenvalues[0], check_value);
      }
   }

   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &eigenvalues);

   return SCIP_OKAY;
}

/** separates the current solution */
static
SCIP_RETCODE separateSol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         cons,               /**< constraint to process */
   SCIP_SOL*          sol,                /**< primal solution that should be separated */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff = NULL;
   int nvars;
   char* cutname;
   int snprintfreturn;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars ) );

   SCIP_CALL( cutUsingEigenvector(scip, cons, sol, coeff, &lhs) );

   SCIP_VAR** vars = NULL;
   SCIP_Real* vals = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars ) );

   int len = 0;

   for (int j = 0; j < nvars; ++j)
   {
      if ( SCIPisZero(scip, coeff[j]) )
         continue;

      vars[len] = SCIPvarGetProbvar(consdata->vars[j]);
      vals[len] = coeff[j];
      ++len;
   }

   SCIP_ROW* row;
   SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255) );
   snprintfreturn = SCIPsnprintf(cutname, 255, "sepa_eig_sdp_%d", ++neigvalcuts);
   assert ( snprintfreturn < 256 ); /* it returns the number of positions that would have been needed, if that is more than 255, it failed */
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, cutname , lhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, row, len, vars, vals) );

   if ( SCIPisCutEfficacious(scip, sol, row) )
   {
      SCIP_Bool infeasible;
      SCIP_CALL( SCIPaddCut(scip, sol, row, FALSE, &infeasible) );
      if ( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;
      SCIP_CALL( SCIPresetConsAge(scip, cons) );
   }

   SCIP_CALL( SCIPreleaseRow(scip, &row) );
   SCIPfreeBufferArray(scip, &cutname);

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coeff);

   return SCIP_OKAY;
}

/** approximates the sdpcone using the fact that every diagonal entry must be non-negative, so it adds the LP-cut
 *  \f$ \sum_{j = 1}^m (A_j)_{kk} y_j - (A_0)_{kk} \geq 0 \quad \forall k \leq n \f$
 */
static
SCIP_RETCODE diagGEzero(
   SCIP*             scip,       /**< SCIP data structure */
   SCIP_CONS**       conss,      /**< array of constraitns */
   int               nconss,     /**< number of constraints */
   int*              naddconss   /**< pointer to store how many constraints were added */
   )
{
   int blocksize;

   SCIP_CONSDATA* consdata;
   int nvars;
   int i;
   int j;
   int k;
   int c;
   SCIP_Real rhs = SCIPinfinity(scip);
   char* cutname;
   int snprintfreturn;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSHDLR* conshdlr;
      conshdlr = SCIPconsGetHdlr(conss[c]);
      assert( conshdlr != NULL );
      const char* conshdlrName;
      conshdlrName = SCIPconshdlrGetName(conshdlr);

      assert ( strcmp(conshdlrName, "SDP") == 0);

      consdata = SCIPconsGetData(conss[c]);

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      rhs = SCIPinfinity(scip);

      SCIP_Real* matrix;
      SCIP_CALL(SCIPallocBufferArray(scip, &matrix, 0.5 * blocksize * (blocksize+1)));
      SCIP_CALL(SCIPconsSdpGetLowerTriangConstMatrix(scip, conss[c], matrix));

      SCIP_Real* cons_array;
      SCIP_Real* lhs_array;
      SCIP_CALL( SCIPallocBufferArray(scip, &cons_array, blocksize * nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &lhs_array, blocksize) );

      /* the lhs is the (k,k)-entry of the constant matrix */
      for (k = 0; k < blocksize; ++k)
         lhs_array[k] = matrix[compLowerTriangPos(k,k)];

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
            if ( consdata->col[j][i] == consdata->row[j][i] )
               cons_array[consdata->col[j][i] * nvars + j] = consdata->val[j][i];
         }
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255));

      /* add the LP-cuts to SCIP */
      for (k = 0; k < blocksize; ++k)
      {
         SCIP_CONS* cons;
         snprintfreturn = SCIPsnprintf(cutname, 255, "diag_ge_zero_%d", ++ndiaggezerocuts);
         assert ( snprintfreturn < 256 ); /* this is the number of positions needed, we gave 255 */

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, consdata->vars, cons_array + k * consdata->nvars, lhs_array[k], rhs,
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         ++(*naddconss);
      }

      SCIPfreeBufferArray(scip, &cutname);
      SCIPfreeBufferArray(scip, &lhs_array);
      SCIPfreeBufferArray(scip, &cons_array);
      SCIPfreeBufferArray(scip, &matrix);
   }

   return SCIP_OKAY;
}

/** presolve-routine that adds some constraints for approximation of the sdpcone, if there is an entry (i,j) in the constant
 *  matrix, than some variable k with \f (A_k)_{ii} \f needs to be >0 because sdp-matrices are diagonal-dominant (and the same for j)
 *
 * not clear if this is really true for all SDPs, probably only works if A_i and A_0 are all semidefinite (or for A_0 at least
 * have positive diagonal entries) and all variables appearing in the SDP constraint are integer, then sum_{A_i_kk >0}
 * 1*y_i >= 1 is feasible, because it means that (sum A_i)_kk > 0 because all diagonal entries are positive (they can't
 * cancel each other) and at least one variable needs to be >=1 becaue this is equal to >0 for integers
 */
static
SCIP_RETCODE diagDominant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   SCIP_Bool* nonzerorows; /* entry i will be 1 if there is an entry \f (A_0)_ij \f for some \f j \neq i \f */
   int blocksize;
   int i;
   int j;
   int nvars;
   int var;
   SCIP_CONS* cons;
   char* cutname;
   int snprintfreturn;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );
   assert ( naddconss != NULL );

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CONSHDLR* conshdlr;
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* conshdlrName;
      conshdlrName = SCIPconshdlrGetName(conshdlr);

      assert ( strcmp(conshdlrName, "SDP") == 0);

      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &nonzerorows, blocksize) );

      /* initialize nonzerorows with FALSE */
      for (j = 0; j < blocksize; j++)
         nonzerorows[j] = FALSE;

      /* iterate over all nonzeroes of the constant matrix and set all corresponding rows/cols to true */
      for (j = 0; j < consdata->constnnonz; j++)
      {
         if ( consdata->constcol[j] != consdata->constrow[j] )
         {
            assert ( ! SCIPisEQ(scip, consdata->constval[j], 0.0) );
            nonzerorows[consdata->constcol[j]] = TRUE;
            nonzerorows[consdata->constrow[j]] = TRUE;
         }
      }

      /* diagvars[i] is an array with all variables with a diagonal entry (i,i) in the corresponding matrix, if nonzerorows[i] is true or NULL otherwise
       * the outer array goes over all rows to ease the access, but only for those that are really needed memory will be allocated */
      int** diagvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &diagvars, blocksize) );
      int* ndiagvars;
      SCIP_CALL( SCIPallocBufferArray(scip, &ndiagvars, blocksize) );
      for (j = 0; j < blocksize; ++j)
      {
         ndiagvars[j] = 0;
         if (nonzerorows[j])
            SCIP_CALL( SCIPallocBufferArray(scip, &(diagvars[j]), nvars) );
      }

      /* find all variables with corresponding diagonal entries for a row with nonzero non-diagonal constant entry */
      for (var = 0; var < nvars; var++)
      {
         for (j = 0; j < consdata->nvarnonz[var]; j++)
         {
            if ( consdata->col[var][j] == consdata->row[var][j] && nonzerorows[consdata->row[var][j]])
            {
               assert ( ! SCIPisEQ(scip, consdata->val[var][j], 0.0) );
               diagvars[consdata->col[var][j]][ndiagvars[consdata->col[var][j]]] = var;
               ndiagvars[consdata->col[var][j]]++;
            }
         }
      }

      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255));

      for (j = 0; j < blocksize; ++j)
      {
         if ( nonzerorows[j] )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &vals, ndiagvars[j]) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, ndiagvars[j]) );

            /* get the corresponding SCIP variables and set all coefficients to 1 */
            for (var = 0; var < ndiagvars[j]; ++var)
            {
               vars[var] = consdata->vars[diagvars[j][var]];
               vals[var] = 1.0;
            }

            snprintfreturn = SCIPsnprintf(cutname, 255, "diag_dom_%d", ++ndiagdomcuts);
            assert ( snprintfreturn < 256 ); /* the return is the number of spots needed, we gave 255 */

            /* add the linear constraint sum_j 1.0 * diagvars[j] >= 1.0 */
            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, cutname , ndiagvars[j], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));
            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));
            (*naddconss)++;

            SCIPfreeBufferArray(scip, &vars);
            SCIPfreeBufferArray(scip, &vals);
            SCIPfreeBufferArray(scip, &diagvars[j]);
         }
      }

      SCIPfreeBufferArray(scip, &cutname);
      SCIPfreeBufferArray(scip, &ndiagvars);
      SCIPfreeBufferArray(scip, &diagvars);
      SCIPfreeBufferArray(scip, &nonzerorows);
   }

   return SCIP_OKAY;
}

/** detects if there are blocks with size one and transfers them to lp-rows */
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*                 scip,               /**<SCIP data structure*/
   SCIP_CONS**           conss,              /**<array of constraints to check*/
   int                   nconss,             /**<number of constraints to check*/
   int*                  naddconss,          /**<pointer to store how many constraints were added*/
   int*                  ndelconss           /**<pointer to store how many constraints were deleted*/
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
   int count;
   int var;
   char* cutname;
   int snprintfreturn;

   for (i = 0; i < nconss; ++i)
   {
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      SCIP_CONSDATA* consdata = SCIPconsGetData(conss[i]);

      /* if there is a 1x1 SDP-Block */
      if ( consdata->blocksize == 1 )
      {
         nvars = consdata->nvars;
         nnonz = consdata->nnonz;
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
               vars[count] = consdata->vars[var];
               count++;
            }
         }

         //get rhs
         assert ( consdata->constnnonz <= 1 ); /* the 1x1 constant matrix may only have one entry */

         rhs = (consdata->constnnonz == 1) ? consdata->constval[0] : 0.0; /* if this one entry is not 0, than this is the rhs, otherwise it's 0 */

         //add new linear cons
         SCIP_CONS* cons;
         SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255) );
         snprintfreturn = SCIPsnprintf(cutname, 255, "1x1block_%d", ++n1x1blocks);
         assert ( snprintfreturn < 256 ); /* the return is the number of spots needed, we gave 255 */

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, vars, coeffs, rhs, SCIPinfinity(scip),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL(SCIPaddCons(scip, cons));
         SCIP_CALL(SCIPreleaseCons(scip, &cons));

         (*naddconss)++;

         //delete old 1x1 sdpcone
         SCIP_CALL(SCIPdelCons(scip, conss[i]));
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &cutname);
         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &coeffs);
      }
   }
   return SCIP_OKAY;
}

/** presolve routine that looks through the data and eliminates fixed variables */
static
SCIP_RETCODE fixVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array with constraints to check */
   int                   nconss              /**< number of constraints to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int i;
   int nfixednonz;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   int c;
   int v;
   int oldnvars;
   int arraylength;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );

   for (c = 0; c < nconss; ++c)
   {
      conshdlr = SCIPconsGetHdlr(conss[c]);
      assert( conshdlr != NULL );

      assert ( strcmp(SCIPconshdlrGetName(conshdlr), "SDP") == 0);

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* allocate memory to save nonzeros that need to be fixed */
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedcol, consdata->nnonz));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedrow, consdata->nnonz));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedval, consdata->nnonz));

      /* initialize this with zero for each block */
      nfixednonz = 0;

      oldnvars = consdata->nvars;

      for (v = 0; v < consdata->nvars; v++)
      {
         SCIP_VAR* var;

         var = SCIPvarGetProbvar(consdata->vars[v]);

         /* check if the variable is fixed in SCIP */
         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
         {
            assert ( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) );

            /* the nonzeros are saved to later be inserted into the constant part (this is only done after all nonzeros of fixed variables have been
             * assembled, because we need to sort the constant nonzeros and loop over them, which we only want to do once and not once for each fixed
             * variable) */
            for (i = 0; i < consdata->nvarnonz[v]; i++)
            {
               savedcol[nfixednonz] = consdata->col[v][i];
               savedrow[nfixednonz] = consdata->row[v][i];
               savedval[nfixednonz] = -consdata->val[v][i] * SCIPvarGetLbGlobal(var);
               /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */
               nfixednonz++;
               consdata->nnonz--;
            }
            /* as the variables don't need to be sorted, we just put the last variable into the empty spot and decrease sizes by one (at the end) */
            consdata->col[v] = consdata->col[consdata->nvars - 1];
            consdata->row[v] = consdata->row[consdata->nvars - 1];
            consdata->val[v] = consdata->val[consdata->nvars - 1];
            consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
            consdata->vars[v] = consdata->vars[consdata->nvars - 1];
            consdata->nvars--;
            v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
         }
      }

      /* free memory for the sdp-arrays equal to the number of fixed variables */
      if (consdata->nvars < oldnvars)
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col), oldnvars, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row), oldnvars, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val), oldnvars, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->nvarnonz), oldnvars, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->vars), oldnvars, consdata->nvars) );
      }

      /* allocate the maximally needed memory */
      arraylength = consdata->constnnonz + nfixednonz;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), consdata->constnnonz, arraylength) );

      /* insert the fixed variables into the constant arrays */
      SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, nfixednonz, FALSE, 1.0, consdata->constrow, consdata->constcol,
            consdata->constval, &(consdata->constnnonz), arraylength) );

      /* shrink the arrays if nonzeros could be combined */
      assert ( consdata->constnnonz <= arraylength );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), arraylength, consdata->constnnonz) );

      /* free the saved arrays */
      SCIPfreeBlockMemoryArray(scip, &savedval, consdata->nnonz);
      SCIPfreeBlockMemoryArray(scip, &savedrow, consdata->nnonz);
      SCIPfreeBlockMemoryArray(scip, &savedcol, consdata->nnonz);
   }

   return SCIP_OKAY;
}

/** presolve routine that looks through the data and handles aggregated and multiaggregated variables */
static
SCIP_RETCODE multiaggrVars(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_CONS**       conss,   /**< array with constraints to check */
   int               nconss  /**< number of constraints to check */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;
   int block;
   int var;
   int i;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   SCIP_VAR** aggrvars;
   SCIP_Real* scalars;
   int naggrvars;
   SCIP_Real constant;
   int requiredsize;
   int globalnvars;
   int aggrind;
   int aggrconsind;
   int naggrnonz;
   int vararraylength;
   int aggrtargetlength;

   assert ( scip != NULL );
   assert ( conss != NULL );
   assert ( nconss >= 0 );

   for (block = 0; block < nconss; ++block)
   {
      conshdlr = SCIPconsGetHdlr(conss[block]);
      assert( conshdlr != NULL );

      assert( strcmp(SCIPconshdlrGetName(conshdlr), "SDP") == 0 );

      consdata = SCIPconsGetData(conss[block]);
      globalnvars = SCIPgetNVars(scip);
      vararraylength = consdata->nvars;

      for (var = 0; var < consdata->nvars; var++)
      {
         if (SCIPvarGetStatus(consdata->vars[var]) == SCIP_VARSTATUS_AGGREGATED ||
             SCIPvarGetStatus(consdata->vars[var]) == SCIP_VARSTATUS_MULTAGGR)
         {
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &aggrvars, globalnvars));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &scalars, globalnvars));

            aggrvars[0] = consdata->vars[var];
            naggrvars = 1;

            /* get the variables this var was aggregated to */
            SCIP_CALL(SCIPgetProbvarLinearSum(scip, aggrvars, scalars, &naggrvars, globalnvars, &constant, &requiredsize, TRUE));
            assert( requiredsize <= globalnvars ); /* requiredsize is the number of empty spots in aggrvars needed, globalnvars is the number
                                                    * of spots we provided */

            /* save the nonzeroes of the (multi)aggregated var */
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedcol, consdata->nvarnonz[var]));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedrow, consdata->nvarnonz[var]));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &savedval, consdata->nvarnonz[var]));

            naggrnonz = 0;

            for (i = 0; i < consdata->nvarnonz[var]; i++)
            {
               savedcol[naggrnonz] = consdata->col[var][i];
               savedrow[naggrnonz] = consdata->row[var][i];
               savedval[naggrnonz] = consdata->val[var][i];
               naggrnonz++;
            }

            assert ( naggrnonz == consdata->nvarnonz[var] );

            /* sort them by nondecreasing row and then col to make the search for already existing entries easier (this is done here, because it
             * only needs to be done once and not for each variable this is multiaggregated to) */
            SdpVarfixerSortRowCol(savedrow, savedcol, savedval, naggrnonz);

            /* fill the empty spot of the (mutli-)aggregated variable with the last variable of this constraint (as they don't have to be sorted) */
            consdata->col[var] = consdata->col[consdata->nvars - 1];
            consdata->row[var] = consdata->row[consdata->nvars - 1];
            consdata->val[var] = consdata->val[consdata->nvars - 1];
            consdata->nvarnonz[var] = consdata->nvarnonz[consdata->nvars - 1];
            consdata->vars[var] = consdata->vars[consdata->nvars - 1];
            consdata->nvars--;
            var--; /* we need to check again if the variable we just shifted to this position also needs to be (multi-)aggregated */

            /* iterate over all variables this was aggregated to and insert the corresponding nonzeroes */
            for (aggrind = 0; aggrind < naggrvars; aggrind++)
            {
               /* check if the variable already exists in this block */
               aggrconsind = -1;
               for (i = 0; i < consdata->nvars; i++)
               {
                  if (consdata->vars[i] == aggrvars[aggrind])
                  {
                     aggrconsind = i;
                     break;
                  }
               }

               if (aggrconsind > -1)
               {
                  /* if the varialbe to aggregate to is already part of this sdp-constraint, just add the nonzeros of the old variable to it */

                  /* resize the arrays to the maximally needed length */
                  aggrtargetlength = consdata->nvarnonz[aggrconsind] + naggrnonz;
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );

                  SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, naggrnonz, TRUE, scalars[aggrind],
                             consdata->row[aggrconsind], consdata->col[aggrconsind], consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]),
                             aggrtargetlength) );

                  /* shrink them again if nonzeros could be combined */
                  assert ( consdata->nvarnonz[aggrconsind] <= aggrtargetlength );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
               }
               else
               {
                  /* the variable has to be added to this constraint */

                  /* check if we have to enlarge the arrays */
                  if (consdata->nvars == vararraylength)
                  {
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, vararraylength, globalnvars) );
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, vararraylength, globalnvars) );
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, vararraylength, globalnvars) );
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, vararraylength, globalnvars) );
                     SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, vararraylength, globalnvars) );
                     vararraylength = globalnvars; /* we don't want to enlarge this by one for every variable added, so we immediately set it to the
                                                    * maximum possible size */
                  }

                  /* we insert this variable at the last position, as the ordering doesn't matter */
                  consdata->vars[consdata->nvars] = aggrvars[aggrind];
                  consdata->nvarnonz[consdata->nvars] = naggrnonz; /* as there were no nonzeros thus far, the number of nonzeros equals the number
                                                                    * of nonzeros of the aggregated variable */

                  /* as there were no nonzeros thus far, we can just duplicate the saved arrays to get the nonzeros for the new variable */
                  SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), savedcol, naggrnonz) );
                  SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), savedcol, naggrnonz) );
                  SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), savedcol, naggrnonz) );

                  consdata->nvars++;
               }
            }

            /* merge the aggregated nonzeros into the constant arrays */

            /* reallocate the constant arrays to the maximally needed size */
            aggrtargetlength = consdata->constnnonz + naggrnonz;
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), consdata->constnnonz, aggrtargetlength) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), consdata->constnnonz, aggrtargetlength) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), consdata->constnnonz, aggrtargetlength) );

            SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, naggrnonz, TRUE, constant, consdata->constrow,
                                               consdata->constcol, consdata->constval, &(consdata->constnnonz), aggrtargetlength) );

            /* shrink the arrays again if nonzeros could be combined */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), aggrtargetlength, consdata->constnnonz) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), aggrtargetlength, consdata->constnnonz) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), aggrtargetlength, consdata->constnnonz) );

            /* free all arrays that are no longer needed */
            SCIPfreeBlockMemoryArray(scip, &savedval, consdata->nvarnonz[var]);
            SCIPfreeBlockMemoryArray(scip, &savedrow, consdata->nvarnonz[var]);
            SCIPfreeBlockMemoryArray(scip, &savedcol, consdata->nvarnonz[var]);
            SCIPfreeBlockMemoryArray(scip, &scalars, globalnvars);
            SCIPfreeBlockMemoryArray(scip, &aggrvars, globalnvars);
         }
      }

      /* shrink the variable arrays if they were enlarged too much */
      assert ( consdata->nvars <= vararraylength );
      if (consdata->nvars < vararraylength)
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, vararraylength, consdata->nvars) );
      }


      /* recompute sdpnnonz */
      consdata->nnonz = 0;
      for (var = 0; var < consdata->nvars; var++)
         consdata->nnonz += consdata->nvarnonz[var];
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
   int v;
   SCIP_VAR* var;

   assert ( SCIPisTransformed(scip) );

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      assert ( strcmp(hdlrName, "SDP") == 0);

      consdata = SCIPconsGetData(conss[i]);

      for (v = 0; v < consdata->nvars; v++)
      {
         SCIP_CALL( SCIPgetTransformedVar(scip, consdata->vars[v], &var) );
         consdata->vars[v] = var;
      }
   }
   return SCIP_OKAY;
}

/** locks a variable up if the corresponding constraint matrix is not positive semidefinite, locks it down if it is not negative semidefinite */
static
SCIP_DECL_CONSLOCK(consLockSdp)
{
   SCIP_Real* Aj;
   SCIP_CONSDATA* consdata;
   int blocksize;
   int var;
   int nvars;
   SCIP_Real eigenvalue;

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;
   nvars = consdata->nvars;

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &Aj, blocksize * blocksize));

   for (var = 0; var < nvars; var++)
   {
      SCIP_CALL(SCIPconsSdpGetFullAj(scip, cons, var, Aj));

      /* compute the smallest eigenvalue */
      SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
      {
         /* as the lowest eigenvalue is negative, the matrix is not positive semidefinite, so adding more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlocksneg, nlockspos) );
      }

      /* compute the biggest eigenvalue */
      SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
      if ( SCIPisPositive(scip, eigenvalue) )
      {
         /* as the biggest eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlockspos, nlocksneg) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &Aj, blocksize * blocksize);

   return SCIP_OKAY;
}


/** after presolving variables are fixed and multiaggregated */
static
SCIP_DECL_CONSEXITPRE(consExitpreSdp)
{
   assert ( scip != NULL );
   assert ( conss != NULL );

   SCIP_CALL(fixVars(scip, conss, nconss));
   SCIP_CALL(multiaggrVars(scip, conss, nconss));

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSdp)
{
   assert( result != 0 );

   if ( nrounds == 0 )
      SCIP_CALL( diagGEzero(scip, conss, nconss, naddconss) );

   SCIP_CALL( move_1x1_blocks_to_lp(scip, conss, nconss, naddconss, ndelconss) );

   SCIP_CALL( fixVars(scip, conss, nconss) );

   if ( nrounds == 0 )
   {
      SCIP_CALL( diagDominant(scip, conss, nconss, naddconss) ); /*TODO: could be activated for some problem classes
      but doesn't seem to work in the general case */
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** creates transformed constraint */
static
SCIP_DECL_CONSTRANS(consTransSdp)
{
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   int i;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   /* copy some general data */
   targetdata->nvars = sourcedata->nvars;
   targetdata->nnonz = sourcedata->nnonz;
   targetdata->blocksize = sourcedata->blocksize;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->nvarnonz), sourcedata->nvarnonz, sourcedata->nvars) );

   /* copy the non-constant nonzeros */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->col), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->row), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->val), sourcedata->nvars) );

   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->col[i]), sourcedata->col[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->row[i]), sourcedata->row[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->val[i]), sourcedata->val[i], sourcedata->nvarnonz[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->vars), sourcedata->nvars));

   /* copy & transform the vars array */
   for (i = 0; i < sourcedata->nvars; i++)
      targetdata->vars[i] = SCIPvarGetTransVar(sourcedata->vars[i]);

   /* copy the constant nonzeros */
   targetdata->constnnonz = sourcedata->constnnonz;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constcol), sourcedata->constcol, sourcedata->constnnonz));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constrow), sourcedata->constrow, sourcedata->constnnonz));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constval), sourcedata->constval, sourcedata->constnnonz));

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
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], sol, checkintegrality, checklprows, printreason, result) );
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
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], NULL, 0, 0, 0, result) );

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
   char* cutname;
   int snprintfreturn;

   for (int i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], NULL, 0, 0, 0, result) );
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
      SCIP_CALL (SCIPallocBufferArray(scip, &cutname, 255) );
      snprintfreturn = SCIPsnprintf(cutname, 255, "eigenvectorcut_enfolp_%d", ++neigveccuts);
      assert ( snprintfreturn < 256 ); /* this is the number of spots needed, we gave 255 */
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, cutname , lhs, rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for (int j = 0; j < nvars; ++j)
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], coeff[j]) );
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
      SCIPfreeBufferArray(scip, &cutname);
      SCIPfreeBufferArray(scip, &coeff);
   }
   if (all_feasible)
      return SCIP_OKAY;
   if (separated)
      *result = SCIP_SEPARATED;


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
      SCIP_CALL(separateSol(scip, conshdlr, conss[i], sol, result));

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSdp)
{
   *result = SCIP_DIDNOTFIND;
   for (int i = 0; i < nusefulconss; ++i)
      SCIP_CALL(separateSol(scip, conshdlr, conss[i], NULL, result));

   return SCIP_OKAY;
}


/** delete method of sdp constrainthandler */
static
SCIP_DECL_CONSDELETE(consDeleteSdp)
{
   int i;

   assert(consdata != NULL);

   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->col[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->row[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->val[i], (*consdata)->nvarnonz[i]);
   }
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->col, (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->row, (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->val, (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->nvarnonz, (*consdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->constcol, (*consdata)->nnonz);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->constrow, (*consdata)->nnonz);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->constval, (*consdata)->nnonz);
   SCIPfreeMemory(scip, consdata);

   return SCIP_OKAY;
}

/* print a SDP constraint */
static
SCIP_DECL_CONSPRINT(consPrintSdp)
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* fullmatrix;
   int v;
   int i;
   int j;

   assert ( cons != NULL );

   consdata = SCIPconsGetData(cons);

   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, consdata->blocksize * consdata->blocksize) );

   /* print the non-constant matrices, for this they first have to be assembled in fullmatrix */
   for (v = 0; v < consdata->nvars; v++)
   {
      /* assemble the matrix */

      /* first initialize it with zero */
      for (i = 0; i < consdata->blocksize; i++)
      {
         for (j = 0; j < consdata->blocksize; j++)
            fullmatrix[i * consdata->blocksize + j] = 0.0;
      }

      /* then add the nonzeros */
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         fullmatrix[consdata->row[v][i] * consdata->blocksize + consdata->col[v][i]] = consdata->val[v][i]; /* lower triangular entry */
         fullmatrix[consdata->col[v][i] * consdata->blocksize + consdata->row[v][i]] = consdata->val[v][i]; /* upper triangular entry */
      }
      /* print it */
      SCIPinfoMessage(scip, file, "+\n");
      for (i = 0; i < consdata->blocksize; i++)
      {
         SCIPinfoMessage(scip, file, "( ");
         for (j = 0; j < consdata->blocksize; j++)
            SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]);
         SCIPinfoMessage(scip, file, ")\n");
      }
      SCIPinfoMessage(scip, file, "* %s\n", SCIPvarGetName(consdata->vars[v]));
   }

   /* print the constant-matrix */

   /* assemble the matrix */

   /* first initialize it with zero */
   for (i = 0; i < consdata->blocksize; i++)
   {
      for (j = 0; j < consdata->blocksize; j++)
         fullmatrix[i * consdata->blocksize + j] = 0.0;
   }

   /* then add the nonzeros */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      fullmatrix[consdata->constrow[i] * consdata->blocksize + consdata->constcol[i]] = consdata->constval[i]; /* lower triangular entry */
      fullmatrix[consdata->constcol[i] * consdata->blocksize + consdata->constrow[i]] = consdata->constval[i]; /* upper triangular entry */
   }
   /* print it */
   SCIPinfoMessage(scip, file, "-\n");
   for (i = 0; i < consdata->blocksize; i++)
   {
      SCIPinfoMessage(scip, file, "( ");
      for (j = 0; j < consdata->blocksize; j++)
         SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]);
      SCIPinfoMessage(scip, file, ")\n");
   }
   SCIPinfoMessage(scip, file, ">=0\n");

   SCIPfreeBufferArray(scip, &fullmatrix);
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
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSdp) );

   return SCIP_OKAY;
}

/** get the data belonging to a single SDP-constraint
 *  in arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown */
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
   SCIP_Var**            vars,               /**< the SCIP variables present in this constraint, indexing equals indices in col/row/val */
   int*                  constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint, also length of the const arrays */
   int*                  constcol,           /**< pointer to column indices of the constant nonzeroes */
   int*                  constrow,           /**< pointer to row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< pointer to values of the constant nonzeroes */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   const char* name;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( nvars != NULL );
   assert ( nnonz != NULL );
   assert ( blocksize != NULL );
   assert ( arraylength != NULL );
   assert ( nvarnonz != NULL );
   assert ( col != NULL );
   assert ( row != NULL );
   assert ( val != NULL );
   assert ( vars != NULL );
   assert ( constnnonz != NULL );
   assert ( constcol != NULL );
   assert ( constrow != NULL );
   assert ( constval != NULL );

   consdata = SCIPconsGetData(cons);
   name = SCIPconsGetName(cons);

   *nvars = consdata->nvars;
   *nnonz = consdata->nnonz;
   *blocksize = consdata->blocksize;

   vars = consdata->vars;

   /* check that the sdp-arrays are long enough to store the information */
   if (*arraylength < consdata->nvars)
   {
      SCIPdebugMessage("nvarnonz, col, row and val arrays were not long enough to store the information for cons %s, they need to be at least"
                       "size %d, given was only length %d! \n", name, consdata->nvars, *arraylength);
      *arraylength = consdata->nvars;
   }
   else
   {
      for (i = 0; i < consdata->nvars; i++)
      {
         nvarnonz[i] = consdata->nvarnonz[i];
         col[i] = consdata->col[i];
         row[i] = consdata->row[i];
         val[i] = consdata->val[i];
      }
   }

   /* set the constant pointers (if a constant part exists) */
   if (consdata->constnnonz > 0)
   {
      if (consdata->constnnonz > *constnnonz)
      {
         SCIPdebugMessage("The constant nonzeros arrays were not long enough to store the information for cons %s, they need to be at least"
                                "size %d, given was only length %d! \n", name, consdata->constnnonz, *constnnonz);
      }
      else
      {
         for (i = 0; i < consdata->constnnonz; i++)
         {
            constcol[i] = consdata->constcol[i];
            constrow[i] = consdata->constrow[i];
            constval[i] = consdata->constval[i];
         }
      }
      *constnnonz = consdata->constnnonz;
   }

   return SCIP_OKAY;
}

/** gets the number of nonzeroes and constant nonzeroes for this SDP constraint, either nnonz oder constnnonz may be NULL */
EXTERN
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  constnnonz          /**< number of nonzeroes in the constant part of this SDP constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert ( scip != NULL );
   assert ( cons != NULL );

   consdata = SCIPconsGetData(cons);

   if (nnonz != NULL)
      *nnonz = consdata->nnonz;

   if (constnnonz != NULL)
      *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the full constraint Matrix \f A_j \f for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int                   j,                  /**< the variable j to get the corresponding matrix \f A_j \f for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f A_j \f */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( j >= 0 );
   assert ( Aj != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   assert ( j < consdata->nvars );

   for (i = 0; i < consdata->nvarnonz[j]; i++)
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

   for (i = 0; i < consdata->constnnonz; i++)
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

   for (i = 0; i < consdata->constnnonz; i++)
   {
      mat[compLowerTriangPos(consdata->constcol[i], consdata->constrow[i])] = consdata->constval[i];
   }

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
   assert ( nnonz == 0 || (nvarnonz != NULL && col != NULL && row != NULL && val != NULL ));
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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars));

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

   for (i = 0; i < nvars; i++)
      consdata->vars[i] = vars[i];


   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );


   return SCIP_OKAY;
}
} //TODO: these two brackets are obviously too much, if a remove them the compiler says } missing, now it gives a warning for syntax error for both
