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

/**@file   objconshdlr_sdp.cpp
 * @brief  constraint handler for sdp-constraints
 * @author Sonja Mars, Lars Schewe, Tristan Gally
 */

//#define SCIP_DEBUG

#include "objconshdlr_sdp.h"

#include <cassert>                      // for assert
//#include <cmath>                        // for floor //TODO: lint says it's not needed
#include <cstring>                      // for NULL, strcmp

#include "SdpCone.h"                    // for SdpCone
#include "config.h"                     // for F77_FUNC

#include "scip/cons_linear.h"           // for SCIPcreateConsLinear
#include "scip/scip.h"                  // for SCIPallocBufferArray, etc

/** struct with consdata*/
struct SCIP_ConsData
{
   SdpCone *sdpcone;
};


extern "C" {
/** BLAS Fortran subroutine DGEMV */
void F77_FUNC(dgemv, DGEMV)(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);
}

/** call matrix-multipication,
 *@note all memory must be allocated outside
 */
static
SCIP_RETCODE Blas_DGEMV(
   int num_rows,                /**<numer of rows in matrix*/
   int num_cols,                /**<numver of cols in matrix*/
   SCIP_Real alpha,             /**<*/
   double* matrix,              /**<the matrix we want to multiply*/
   double* vector_to_multiwith, /**<vector we want to multiply with the matrix*/
   double beta,                 /**<*/
   double* output_vector        /**<vector where the result is put in*/)
{
   char TRANS = 'N';
   int M = num_rows;
   int N = num_cols;
   double ALPHA = alpha;
   double* A = matrix;
   int LDA = num_rows;
   double* X = vector_to_multiwith;
   int INCX = 1;
   double BETA = beta;
   double* Y = output_vector;
   int INCY = 1;

   F77_FUNC(dgemv, DGEMV)(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);

   return SCIP_OKAY;

}

extern "C" {
/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)( char* JOBZ, char* RANGE, char* UPLO,
   int* N, double* A, int* LDA,
   double* VL, double* VU,
   int* IL, int* IU,
   double* ABSTOL, int* M, double* W, double* Z,
   int* LDZ, int* ISUPPZ, double* WORK,
   int* LWORK, int* IWORK, int* LIWORK,
   int* INFO );
}

/**computes the i-th eigenvalue, where 1 is the smallest and n the largest*/
static
SCIP_RETCODE computeIthEigenvalue(
   SCIP* scip,               /**<SCIP data structure*/
   bool computeeigenvectors, /**<should also the eigenvectors be computed*/
   int n,                    /**<blocksize of matrix*/
   double* A,                /**<matrix from which eigenvalues should be computed*/
   int i,                    /**<this i-th eigenvalues is computed*/
   double* eigenvalue,       /**<pointer to store eigenvalue*/
   double* eigenvector       /**<pointer to store eigenvector*/
   )
{
   int     N = n;
   int     INFO;
   char    JOBZ = computeeigenvectors ? 'V' : 'N';
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

   /// standard LAPACK workspace query
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

   if( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   /// allocate workspace
   LWORK = (int) WSIZE + 1;
   LIWORK = WISIZE;

   SCIP_CALL( SCIPallocBufferArray(scip, &WORK, LWORK) );
   SCIP_CALL( SCIPallocBufferArray(scip, &IWORK, LIWORK) );

   SCIP_CALL( SCIPallocBufferArray(scip, &WTMP, N) );


   /// call the function

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

   if( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d\n", INFO);
      return SCIP_ERROR;
   }

   /// handle output

   *eigenvalue = WTMP[0];

   SCIPfreeBufferArray(scip, &WORK);
   SCIPfreeBufferArray(scip, &IWORK);
   SCIPfreeBufferArray(scip, &WTMP);

   return SCIP_OKAY;
}

/** separate current solution using a cut, with the eigenvectors and -values of the solution matrix
 * @note: this function computes the eigenvectors of the matrix, takes the first one and multiplies the matrix with it
 * x^T*A_i*x = coeff[i],  x^T*A_0*x =lhs

 */
static
SCIP_RETCODE cut_using_eigenvector(
   SCIP_SOL* sol,                   /**<solution to separate*/
   SCIP* scip,                      /**<SCIP data structure*/
   SdpCone*  sdpcone,               /**<sdpcone with data*/
   SCIP_Real* coeff,                /**<coefficients of the computed cut*/
   SCIP_Real *lhs                   /**<lhs of the computed cut*/
   )
{
   *lhs = 0.0;
   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   SCIP_Real* const_matrix = NULL;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real* output_vector = NULL;
   int blocksize = sdpcone->get_blocksize();
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &const_matrix, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvector, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &output_vector, blocksize) );

   //compute the matrix \sum_i A_i x_i
   SCIP_CALL( sdpcone->assemble_matrix_from_solution(matrix, sol));

   SCIP_CALL( computeIthEigenvalue(scip, TRUE, blocksize, matrix, 1, eigenvalues, eigenvector) );

   SCIP_CALL(sdpcone->get_constant_matrix(const_matrix));
   SCIP_CALL(Blas_DGEMV(blocksize, blocksize, 1.0, const_matrix, eigenvector, 0.0, output_vector));

   for (int j = 0; j < blocksize; ++j)
   {
      *lhs += eigenvector[j] * output_vector[j];
   }


   for (int j = 0; j < sdpcone->get_nvars(); ++j)
   {
      SCIP_CALL(sdpcone->form_inner_prod_with_constraint_matrix(j, eigenvector, &coeff[j]));
   }

   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &eigenvalues);
   SCIPfreeBufferArray(scip, &const_matrix);
   SCIPfreeBufferArray(scip, &eigenvector);
   SCIPfreeBufferArray(scip, &output_vector);

   return SCIP_OKAY;
}

SCIP_RETCODE cons_check(
   SCIP*              scip,               /**< SCIP data structure */
   SdpCone*           sdpcone,            /**< sdpcone to check positive semidefiniteness for */
   SCIP_SOL*          sol,                /**< the solution to check feasibility for */
   SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
   )
{
   SCIP_Real* eigenvalues = NULL;
   SCIP_Real* matrix = NULL;
   int blocksize = sdpcone->get_blocksize();
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * blocksize) );

   SCIP_CALL(sdpcone->assemble_matrix_from_solution(matrix, sol));

   SCIP_CALL( computeIthEigenvalue(scip, FALSE, blocksize, matrix, 1, eigenvalues, NULL) );

   //we are going to use one of the dimacs error norms for checking feasiblity.
   //we use the second one: err=max{0, -lambda_min(x)/(1+maximumentry of rhs}

   double check_value;
   double max_rhs = sdpcone->get_max_rhs();
   check_value = (-eigenvalues[0]) / (1 + max_rhs);

   *result = SCIPisLE(scip, check_value, 0.0) ? SCIP_FEASIBLE : SCIP_INFEASIBLE;

#ifdef SCIP_DEBUG
   if( *result == SCIP_INFEASIBLE)
   {
      SCIPdebugMessage("In cons_check a matrix was found not to be sdp because of eigenvalue = %f, dimacs error norm = %f \n",
            eigenvalues[0], check_value);
   }
#endif

   SCIPfreeBufferArray(scip, &matrix);
   SCIPfreeBufferArray(scip, &eigenvalues);

   // TODO: checkintegrality ?!?, checklprows ?!?, printreason ?!?

   return SCIP_OKAY;
}

/**separates the current solution*/
static
SCIP_RETCODE separate_curr_sol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         conss,              /**< array of constraints to process */
   SCIP_SOL*          sol,                /**< primal solution that should be separated */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */)
{
   SdpCone*  sdpcone;
   SCIP_CALL( getSdpCone(scip, conss, &sdpcone ) );
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff = NULL;

   int nvars = sdpcone->get_nvars();
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars ) );


   SCIP_CALL(cut_using_eigenvector(sol, scip, sdpcone, coeff, &lhs));

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

/**approximates the sdpcone using some trivial facts, i.e. every diagonal entry must be non-negative*/
static
SCIP_RETCODE trivial_approx(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraitns*/
   int               nconss,     /**<number of constraints*/
   int*              naddconss,  /**<pointer to store how many constraints were added*/
   SCIP_RESULT*      result      /**<pointer to store the result*/
   )
{
   SCIP_Real rhs = SCIPinfinity(scip);

   for (int i = 0; i < nconss; ++i)
   {
      SdpCone* sdpcone;
      //this approximates the sdpcone using the constraints that all diagonal entries must be >=0
      int blocksize;
      SCIP_CALL( getSdpCone(scip, conss[i], &sdpcone ) );

      blocksize = sdpcone->get_blocksize();
      rhs = SCIPinfinity(scip);

      double* matrix;
      SCIP_CALL(SCIPallocBufferArray(scip, &matrix, blocksize * blocksize));
      SCIP_CALL(sdpcone->get_constant_matrix(matrix));

      double* cons_array;
      double* lhs_array;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &cons_array, blocksize * sdpcone->get_nvars()));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &lhs_array, blocksize));

      for (int k = 0; k < blocksize; ++k)
      {
         lhs_array[k] = matrix[k * blocksize + k];
      }

      for (int j = 0; j < sdpcone->get_nvars(); ++j)
      {
         SCIP_CALL( sdpcone->get_constraint_matrix(matrix, j) );
         for (int k = 0; k < blocksize; ++k)
         {
            cons_array[k * sdpcone->get_nvars() + j] = matrix[k * blocksize + k];
         }
      }

      SCIP_VAR** vars;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vars, sdpcone->get_nvars()));

      for (int j = 0; j < sdpcone->get_nvars(); ++j)
      {
         vars[j] = sdpcone->get_var(j);
      }

      for (int k = 0; k < blocksize; ++k)
      {
         SCIP_CONS* cons;

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "cl", sdpcone->get_nvars(), vars, cons_array + k * sdpcone->get_nvars(), lhs_array[k], rhs, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL(SCIPaddCons(scip, cons));
         SCIP_CALL(SCIPreleaseCons(scip, &cons));
         ++(*naddconss);
      }

      SCIPfreeBlockMemoryArray(scip, &vars, sdpcone->get_nvars() );
      SCIPfreeBufferArray(scip, &matrix);
      SCIPfreeBlockMemoryArray(scip, &cons_array, blocksize * sdpcone->get_nvars());
      SCIPfreeBlockMemoryArray(scip, &lhs_array, blocksize);

      *result = SCIP_SUCCESS;
   }
   return SCIP_OKAY;
}

/* not clear if this is really true for all SDPs, probably only works if A_i and A_0 are all semidefinite (or at least have positive diagonal entries) and all variables appearing in the SDP constraint are integer, then sum_{A_i_kk >0} 1*y_i >= 1 is feasible, because it means that (sum A_i)_kk > 0 because all diagonal entries are positive (they can't cancel each other) and at least one variable needs to be >=1 becaue this is equal to >0 for integers */
/**presolve-routine that adds some constraints for approximation of the sdpcone, if there is an entry on the right hand side there must be a corresponding diagonal entry*/
static
SCIP_RETCODE trivial_ineq_from_rhs(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraints*/
   int               nconss,     /**<number of constraints*/
   int*              naddconss,  /**<pointer to store how many constraints were added*/
   SCIP_RESULT*      result      /**<pointer to store, if procedure was sucessfull*/
   )
{
   for (int i = 0; i < nconss; ++i)
   {
      /* ?????? turn this code into an assert ??????? */
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      if ( strcmp(hdlrName, "SDP") != 0)
         continue;

      SCIP_CONSDATA* data = SCIPconsGetData(conss[i]);
      int* const_row;
      SCIP_CALL(SCIPallocBufferArray(scip, &const_row, data->sdpcone->get_const_nnz()));
      int count = 0;
      for ( SdpCone::RhsIterator it = data->sdpcone->rhs_begin(NULL, 0, NULL); it != data->sdpcone->rhs_end(); ++it)
      {
         SdpCone::element el = *it;
         if (el.col != el.row && (count == 0 || el.row != const_row[count - 1]))
         {
            const_row[count] = el.row;
            count++;
         }
      }


      int** diag_var;
      SCIP_CALL(SCIPallocBufferArray(scip, &diag_var, count));
      int* len_diag_var;
      SCIP_CALL(SCIPallocBufferArray(scip, &len_diag_var, count));
      for (int j = 0; j < count; ++j)
      {
         len_diag_var[j] = 0;
         SCIP_CALL(SCIPallocBufferArray(scip, &diag_var[j], data->sdpcone->get_nvars()));
      }

      for ( SdpCone::LhsIterator it = data->sdpcone->lhs_begin(NULL, 0); it != data->sdpcone->lhs_end(); ++it)
      {
            SdpCone::element el = *it;

         for (int j = 0; j < count; ++j)
         {
            if (const_row[j] == el.row && const_row[j] == el.col && (len_diag_var[j] == 0 || el.vidx != diag_var[j][len_diag_var[j] - 1]))
            {
               len_diag_var[j]++;
               diag_var[j][len_diag_var[j] - 1] = el.vidx;

            }
         }
      }

      SCIP_VAR** vars;
      SCIP_Real* vals;

      for (int j = 0; j < count; ++j)
      {
         if(len_diag_var[j] > 0)
         {
            SCIP_CALL(SCIPallocBufferArray(scip, &vals, len_diag_var[j]));
            SCIP_CALL(SCIPallocBufferArray(scip, &vars, len_diag_var[j]));

            for (int k = 0; k < len_diag_var[j]; ++k)
            {
               vars[k] = data->sdpcone->get_var(diag_var[j][k]);
               vals[k] = 1.0;
            }


            //add new linear cons
            SCIP_CONS* cons;

            SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "sum_diag_geq_1", len_diag_var[j], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

            SCIP_CALL(SCIPaddCons(scip, cons));
            SCIP_CALL(SCIPreleaseCons(scip, &cons));


            //adds linear constraint, for summe_j diag_var[j] >=1
            (*naddconss)++;
            SCIPfreeBufferArray(scip, &vars);
            SCIPfreeBufferArray(scip, &vals);
         }
      }
      for (int j = 0; j < count; ++j)
      {
        SCIPfreeBufferArray(scip, &diag_var[j]);
      }

      SCIPfreeBufferArray(scip, &diag_var);
      SCIPfreeBufferArray(scip, &len_diag_var);
      SCIPfreeBufferArray(scip, &const_row);
   }
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**detects if there are blocks with size one and transfers it to a lp-row*/
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*             scip,       /**<SCIP data structure*/
   SCIP_CONS**       conss,      /**<array of constraints to check*/
   int               nconss,     /**<number of constraints to check*/
   int*              naddconss,  /**<pointer to store how many constraints were added*/
   int*              ndelconss,  /**<pointer to store how many constraints were deleted*/
   SCIP_RESULT*      result      /**<pointer to store the result*/
   )
{
   for (int i = 0; i < nconss; ++i)
   {
      /* ?????? turn this code into an assert ??????? */
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      if ( strcmp(hdlrName, "SDP") != 0)
         continue;

      SCIP_CONSDATA* data = SCIPconsGetData(conss[i]);
      if (data->sdpcone->get_blocksize() == 1)
      {
         //Check if there is a 1x1 block (row and col do not matter anymore)
         int nvars = data->sdpcone->get_nvars();
         SCIP_VAR** vars;
         double* coeffs;
         int nnz = data->sdpcone->get_nnz();
         double* rhs;
         int const_nnz = data->sdpcone->get_const_nnz();
         SCIP_CALL(SCIPallocBufferArray(scip, &vars, nvars));
         SCIP_CALL(SCIPallocBufferArray(scip, &coeffs, nnz));
         SCIP_CALL(SCIPallocBufferArray(scip, &rhs, const_nnz));
         for (int j = 0; j < nvars; ++j)
         {
            vars[j] = data->sdpcone->get_var(j);
         }
         //lhs-iterator
         int count = 0;

         for ( SdpCone::LhsIterator it = data->sdpcone->lhs_begin(NULL, 0); it != data->sdpcone->lhs_end(); ++it)
         {
            SdpCone::element el = *it;
            coeffs[count] = el.val;
            count++;
         }

         //rhs-iterator
         count = 0;
         for ( SdpCone::RhsIterator it = data->sdpcone->rhs_begin(NULL, 0, NULL); it != data->sdpcone->rhs_end(); ++it)
         {
            SdpCone::element el = *it;
            rhs[count] = -el.val;
            count++;
         }
         assert(count <= 1);

         //add new linear cons
         SCIP_CONS* cons;

         SCIP_CALL(SCIPcreateConsLinear(scip, &cons, "1x1", data->sdpcone->get_nvars(), vars, coeffs, *rhs, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE));

         SCIP_CALL(SCIPaddCons(scip, cons));
         SCIP_CALL(SCIPreleaseCons(scip, &cons));

         (*naddconss)++;

         //delete old 1x1 sdpcone
         SCIP_CALL(SCIPdelCons(scip, conss[i]));
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &vars);
         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &rhs);

      }
   }

   return SCIP_OKAY;
}

/** presolve routine that looks through the data and eliminates fixed or deleted or aggregated or negated variables
 *
 *  SDP-solver often have problems with fixed variables.
 */
static
SCIP_RETCODE fix_vars(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_CONS**       conss,   /**< array with constraints to check */
   int               nconss,  /**< number of constraints to check */
   SCIP_RESULT*      result   /**< pointer to store the result */
   )
{
   for (int i = 0; i < nconss; ++i)
   {
      /* ?????? turn this code into an assert ??????? */
      SCIP_CONSHDLR* conshdlr;
      conshdlr = SCIPconsGetHdlr(conss[i]);
      assert( conshdlr != NULL );
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(conshdlr);

      if ( strcmp(hdlrName, "SDP") != 0)
         continue;

      SCIP_CONSDATA* data = SCIPconsGetData(conss[i]);
      SCIP_CALL( data->sdpcone->fix_vars() );
   }

   return SCIP_OKAY;
}



/** initialization method of constraint handler (called after problem has been transformed) */
SCIP_RETCODE ObjConshdlrSdp::scip_init(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
   int                nconss              /**< number of constraints in transformed problem */
   )
{
   /* get transformed variables, if we are in the transformed problem */
   SCIP_CONSDATA* consdata = NULL;

   if( SCIPisTransformed(scip) )
   {
      for (int i = 0; i < nconss; ++i)
      {
         /* ?????? turn this code into an assert ??????? */
         SCIP_CONSHDLR* hdlr;
         hdlr = SCIPconsGetHdlr(conss[i]);
         assert(hdlr != NULL);
         const char* hdlrName;
         hdlrName = SCIPconshdlrGetName(hdlr);

         if ( strcmp(hdlrName, "SDP") != 0)
            continue;

         consdata = SCIPconsGetData(conss[i]);

         SCIP_CALL( consdata->sdpcone->transform_vars());
      }
   }

   return SCIP_OKAY;
}


/** locks variables in sdpcone up and down if necessary */
SCIP_RETCODE ObjConshdlrSdp::scip_lock(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the *   constraint handler does not need constraints */
   int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
   int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
   )
{
   SdpCone* sdpcone;

   SCIP_CALL( getSdpCone(scip, cons, &sdpcone) );

   SCIP_Real* matrix = NULL;
   int blocksize = sdpcone->get_blocksize();
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * blocksize) );

   int nvars = sdpcone->get_nvars();

   for (int i = 0; i < nvars; ++i)
   {
      int shrunk_blocksize;
      bool has_lock = false;

      SCIP_CALL( sdpcone->get_shrunk_constraint_matrix(matrix, &shrunk_blocksize, i) );

      /* ?????????? */
      if ( shrunk_blocksize == 0 )
         continue;

      /* determine whether matrix contains only positive/negative entries */
      bool only_neg = true;
      bool only_pos = true;
      for(int j = 0; j < shrunk_blocksize; ++j)
      {
         if ( SCIPisNegative(scip, matrix[j * shrunk_blocksize + j]) )
         {
            only_pos = false;
            if ( ! only_neg )
               break;
         }

         if ( SCIPisPositive(scip, matrix[j * shrunk_blocksize + j]) )
         {
            only_neg = false;
            if ( ! only_pos )
               break;
         }
      }

      /* up lock, if only negative entries */
      if ( only_neg )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, sdpcone->get_var(i), nlocksneg, nlockspos) );
         has_lock = true;
      }

      /* down lock, if only positive entries */
      if ( only_pos )
      {
         SCIP_CALL( SCIPaddVarLocks(scip, sdpcone->get_var(i), nlockspos, nlocksneg) );
         has_lock = true;
      }

      /* treate mixed cases */
      if ( ! only_pos && ! only_neg)
      {
         double eigenvalue;
         SCIP_CALL( computeIthEigenvalue(scip, false, shrunk_blocksize, matrix, 1, &eigenvalue, NULL) );
         if ( SCIPisNegative(scip, eigenvalue) )
         {
            /* if the smallest eigenvalue is negative, we lock the variable upwards, because increasing the variable
             * will make everything more negative, decreasing the variable will make everything more feasible */
            SCIP_CALL( SCIPaddVarLocks(scip, sdpcone->get_var(i), nlocksneg, nlockspos) );
            has_lock = true;
         }

         SCIP_CALL( computeIthEigenvalue(scip, false, shrunk_blocksize, matrix, shrunk_blocksize, &eigenvalue, NULL) );
         if ( SCIPisPositive(scip, eigenvalue) )
         {
            /* if the biggest eigenvalue is positive, we lock the variable downwards, because, increasing is ok, but
             * decreasing can probably make the problem infeasible */
            SCIP_CALL( SCIPaddVarLocks(scip, sdpcone->get_var(i), nlockspos, nlocksneg) );
            has_lock = true;
         }
      }

      assert( has_lock ); // this can fail, if 0 is a eigenvector, which is not possible
   }

   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
SCIP_RETCODE ObjConshdlrSdp::scip_initpre(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
   int                nconss             /**< number of constraints in transformed problem */
   )
{

   for (int j = 0; j < nconss; ++j)
   {
      /* ?????? turn this code into an assert ??????? */
      SCIP_CONSHDLR* hdlr;
      hdlr = SCIPconsGetHdlr(conss[j]);
      assert(hdlr != NULL);
      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      if ( strcmp(hdlrName, "SDP") != 0)
         continue;

      SdpCone* sdpcone;
      SCIP_CALL(getSdpCone(scip, conss[j], &sdpcone ));
      //do not multiaggregate variables
      for (int i = 0; i < sdpcone->get_nvars(); ++i)
      {
         SCIP_CALL(SCIPmarkDoNotMultaggrVar(scip, sdpcone->get_var(i)) );
      }
   }
   return SCIP_OKAY;
}


SCIP_RETCODE ObjConshdlrSdp::scip_presol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< no. of constraints to process */
   int                nrounds,            /**< no. of presolving rounds already done */
   int                nnewfixedvars,      /**< no. of variables fixed since last call to presolving method */
   int                nnewaggrvars,       /**< no. of variables aggregated since last call to presolving method */
   int                nnewchgvartypes,    /**< no. of variable type changes since last call to presolving method */
   int                nnewchgbds,         /**< no. of variable bounds tightend since last call to presolving method */
   int                nnewholes,          /**< no. of domain holes added since last call to presolving method */
   int                nnewdelconss,       /**< no. of deleted constraints since last call to presolving method */
   int                nnewaddconss,       /**< no. of added constraints since last call to presolving method */
   int                nnewupgdconss,      /**< no. of upgraded constraints since last call to presolving method */
   int                nnewchgcoefs,       /**< no. of changed coefficients since last call to presolving method */
   int                nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolving method */
   int*               nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
   int*               naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
   int*               nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
   int*               nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
   int*               naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
   int*               ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
   int*               naddconss,          /**< pointer to count total number of added constraints of all presolvers */
   int*               nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
   int*               nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
   int*               nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
   SCIP_RESULT*       result              /**< pointer to store the result of the presolving call */
   )
{
   assert( result != 0 );

   *result = SCIP_DIDNOTRUN;

   if ( nrounds == 0 )
   {
      SCIP_CALL( trivial_approx(scip, conss, nconss, naddconss, result) );
   }

   SCIP_CALL( move_1x1_blocks_to_lp(scip, conss, nconss, naddconss, ndelconss, result) );

   SCIP_CALL( fix_vars(scip, conss, nconss, result) );

   if ( nrounds == 0 )
   {
      SCIP_CALL( trivial_ineq_from_rhs(scip, conss, nconss, naddconss, result) );  /*TODO: could be activated for some problem classes
      but doesn't seem to work in the general case */
   }

   return SCIP_OKAY;
}


/** creates transformed constraint*/
SCIP_RETCODE ObjConshdlrSdp::scip_trans(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         sourcecons,         /**< source constraint to transform */
   SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
   )
{
   SCIP_CONSDATA* sourcedata = NULL;
   SCIP_CONSDATA* targetdata = NULL;

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocMemory(scip, &targetdata) );//122 Ist eins zu eins aus conshdlrsubtour
   targetdata->sdpcone = new SdpCone(*sourcedata->sdpcone);
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
SCIP_RETCODE ObjConshdlrSdp::scip_check(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   SCIP_SOL*          sol,                /**< the solution to check feasibility for */
   SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
   )
{
   SdpCone* sdpcone;

   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CALL(getSdpCone(scip, conss[i], &sdpcone));
      SCIP_CALL(cons_check(scip, sdpcone, sol, checkintegrality, checklprows, printreason, result));
      if (*result == SCIP_INFEASIBLE)
      {
         return SCIP_OKAY;
      }
   }

   return SCIP_OKAY;
}


/** enforce pseudo solution mehtod, returns didnotrun, if objinfeasible, computes cut otherwise*/
SCIP_RETCODE ObjConshdlrSdp::scip_enfops(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
   )
{
   SdpCone* sdpcone;
   if (objinfeasible)
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }
   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CALL(getSdpCone(scip, conss[i], &sdpcone));
      SCIP_CALL( cons_check(scip, sdpcone, NULL, 0, 0, 0, result) );

      SCIP_Real lhs = 0.0;
      SCIP_Real* coeff = NULL;
      int nvars = sdpcone->get_nvars();
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );
      SCIP_CALL(cut_using_eigenvector(NULL, scip, sdpcone, coeff, &lhs));

      if (*result != SCIP_INFEASIBLE)
      {
         lhs = floor(lhs);
      }

      // TODO: this if does nothing


      SCIPfreeBufferArray(scip, &coeff);

   }

   return SCIP_OKAY;
}


/**enforce lp solution method, must be implemented, but there is no lp in the sdp-case, so returns SCIP_ERROR*/
SCIP_RETCODE ObjConshdlrSdp::scip_enfolp(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
   SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
   )
{
   SdpCone* sdpcone;
   bool all_feasible = TRUE;
   bool separated = FALSE;
   for (int i = 0; i < nconss; ++i)
   {
      SCIP_CALL( getSdpCone(scip, conss[i], &sdpcone ) );
      SCIP_CALL( cons_check(scip, sdpcone, NULL, 0, 0, 0, result) );
      if (*result == SCIP_FEASIBLE)
      {
         continue;
      }
      all_feasible = FALSE;

      int nvars = sdpcone->get_nvars();


      SCIP_Real lhs = 0.0;
      SCIP_Real* coeff = NULL;
      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );

      SCIP_CALL(cut_using_eigenvector(NULL, scip, sdpcone, coeff, &lhs));

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
         SCIP_CALL( SCIPaddVarToRow(scip, row, sdpcone->get_var(j), coeff[j]) );
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


/**separates a solution using constraint specific ideas, gives cuts to scip*/
SCIP_RETCODE ObjConshdlrSdp::scip_sepasol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_SOL*          sol,                /**< primal solution that should be separated */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
   )
{
   *result = SCIP_DIDNOTFIND;
   for (int i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL(separate_curr_sol(scip, conshdlr, conss[i], sol, result));
   }

   return SCIP_OKAY;
}



/** separation method of constraint handler for LP solution */
SCIP_RETCODE ObjConshdlrSdp::scip_sepalp(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS**        conss,              /**< array of constraints to process */
   int                nconss,             /**< number of constraints to process */
   int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
   SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
   )
{
   *result = SCIP_DIDNOTFIND;
   for (int i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL(separate_curr_sol(scip, conshdlr, conss[i], NULL, result));
   }

   return SCIP_OKAY;

}


/**delete method of sdp constrainthandler */
SCIP_RETCODE ObjConshdlrSdp::scip_delete(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*         cons,               /**< the constraint belonging to the constraint data */
   SCIP_CONSDATA**    consdata            /**< pointer to the constraint data to free */
   )
{
   assert(consdata != NULL);

   delete (*consdata)->sdpcone;
   SCIPfreeMemory(scip, consdata);
   return SCIP_OKAY;
}


SCIP_RETCODE getSdpCone(
   SCIP* scip,
   SCIP_CONS* conss,
   SdpCone** sdpcone)
{
   SCIP_CONSDATA* consdata = NULL;
   consdata = SCIPconsGetData(conss);
   assert( consdata != NULL );

   *sdpcone = consdata->sdpcone;

   assert(sdpcone != NULL);

   return SCIP_OKAY;
}



/**creates an sdp-constraint*/
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   const SdpCone&              sdpcone             /**< the sdpcone*/
   )
{

   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSDATA* consdata = NULL;

   conshdlr = SCIPfindConshdlr(scip, "SDP");
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocMemory(scip, &consdata) );

   consdata->sdpcone = new SdpCone(sdpcone);


   assert(consdata->sdpcone != NULL);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );


   return SCIP_OKAY;
}
