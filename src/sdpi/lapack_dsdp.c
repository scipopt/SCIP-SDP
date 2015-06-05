/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG*/

/**@file   lapack_dsdp.c
 * @brief  interface methods for eigenvector computation and matrix multiplication using standard LAPACK and BLAS
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 */

#include <assert.h>

#include "lapack.h"
#include "config.h"                     /* for F77_FUNC */

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"


/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMRY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** transforms a double (that should be integer, but might be off by some numerical error) to an integer by adding an epsilon and rounding down */
#define DOUBLETOINT(x) ((int) (x + 0.5))

/*
 * BLAS/LAPACK Calls
 */

/**@name BLAS/LAPACK Calls */
/**@{ */

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


/** BLAS Fortran subroutine DGEMV */
void F77_FUNC(dgemv, DGEMV)(char* TRANS, int* M, int* N, double* ALPHA, double* A, int* LDA, double* X, int* INCX, double* BETA, double* Y, int* INCY);


/**@} */


/*
 * Functions
 */

/**@name Functions */
/**@{ */

/** computes the i-th eigenvalue, where 1 is the smallest and n the largest, matrix has to be given with all n^2 entries */
SCIP_RETCODE SCIPlapackComputeIthEigenvalue(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to array to store eigenvector */
   )
{
   int N;
   int INFO;
   char JOBZ;
   char RANGE;
   char UPLO;
   int LDA;
   double* WORK;
   int LWORK;
   int* IWORK;
   int LIWORK;
   double* WTMP;
   double ABSTOL;
   int IL;
   int IU;
   int M;
   int LDZ;
   double WSIZE;
   int WISIZE;
   double VL;
   double VU;
   int* ISUPPZ;


   assert( blkmem != NULL );
   assert( n > 0 );
   assert( A != NULL );
   assert( 0 < i && i <= n );
   assert( eigenvalue != NULL );
   assert( ( ! geteigenvectors) || eigenvector != NULL );

   N = n;
   JOBZ = geteigenvectors ? 'V' : 'N';
   RANGE = 'I';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;
   IL = i;
   IU = i;
   M = 1;
   LDZ = n;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1;
   LIWORK = -1;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
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

   /* allocate workspace */
   LWORK = DOUBLETOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &WORK, LWORK) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &IWORK, LIWORK) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &WTMP, N) );
   BMS_CALL( BMSallocBlockMemoryArray(blkmem, &ISUPPZ, 2) );

   /* call the function */
   VL = -1e20;
   VU = 1e20;

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

   /* handle output */
   *eigenvalue = WTMP[0];

   /* free memory */
   BMSfreeBlockMemoryArray(blkmem, &ISUPPZ, 2);
   BMSfreeBlockMemoryArray(blkmem, &WTMP, N);
   BMSfreeBlockMemoryArray(blkmem, &IWORK, LIWORK);
   BMSfreeBlockMemoryArray(blkmem, &WORK, LWORK);

   return SCIP_OKAY;
}

/** performs matrix-vector-multiplication using BLAS */
SCIP_RETCODE SCIPlapackMatrixVectorMult(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   SCIP_Real*            matrix,             /**< the matrix we want to multiply */
   SCIP_Real*            vector,             /**< vector we want to multiply with the matrix */
   SCIP_Real*            result              /**< pointer to store the resulting vector */
   )
{
   char TRANS;
   int M;
   int N;
   double ALPHA;
   double* A;
   int LDA;
   double* X;
   int INCX;
   double BETA;
   double* Y;
   int INCY;

   TRANS = 'N';
   M = nrows;
   N = ncols;
   ALPHA = 1.0;
   A = matrix;
   LDA = nrows;
   X = vector;
   INCX = 1;
   BETA = 0.0;
   Y = result;
   INCY = 1;

   F77_FUNC(dgemv, DGEMV)(&TRANS, &M, &N, &ALPHA, A, &LDA, X, &INCX, &BETA, Y, &INCY);

   return SCIP_OKAY;
}

/**@} */

