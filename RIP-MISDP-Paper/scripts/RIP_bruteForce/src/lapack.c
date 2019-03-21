/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2016 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2016 Zuse Institute Berlin                             */
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
#include <stdio.h>
#include <stdlib.h>

#include "lapack.h"
#include "config.h"                     /* for F77_FUNC */
#include "stddef.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/** transforms a SCIP_Real (that should be integer, but might be off by some numerical error) to an integer by adding an epsilon and rounding down */
#define RealTOINT(x) ((int) (x + 0.5))

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
void F77_FUNC(dgemm, DGEMM)(char* TRANSA, char* TRANSB, int* M, int* N, int* K, double* ALPHA,
      double* A, int* LDA, double* B, int* LDB, double* BETA, double* C, int* LDC);


/**@} */


/*
 * Functions
 */

/**@name Functions */
/**@{ */

/** computes smallest and largest eigenvalue of a given matrix using LAPACK */
void SCIPlapackComputeSmallestLargestEigenvalue(
   int                   n,                  /**< size of matrix */
   double*               A,                  /**< matrix for which eigenvalues should be computed */
   double*               smallestev,         /**< pointer to store smallest eigenvalue */
   double*               largestev           /**< pointer to store largest eigenvalue */
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

   assert( n > 0 );
   assert( A != NULL );
   assert( smallestev != NULL );
   assert( largestev != NULL );

   N = n;
   JOBZ = 'N';
   RANGE = 'I';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;
   IL = 1;
   IU = n;
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
      printf("There was an error when calling DSYEVR. INFO = %d\n", INFO);
   }

   /* allocate workspace */
   LWORK = RealTOINT(WSIZE);
   LIWORK = WISIZE;

   WORK = (double*) malloc(LWORK * sizeof(double));
   IWORK = (int*) malloc(LIWORK * sizeof(int));
   WTMP = (double*) malloc(N * sizeof(double));
   ISUPPZ = (int*) malloc(2 * sizeof(int));

   /* call the function */
   VL = -1e20;
   VU = 1e20;

   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      &IL, &IU,
      &ABSTOL, &M, WTMP, NULL,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO );

   if ( INFO != 0 )
   {
      printf("There was an error when calling DSYEVR. INFO = %d\n", INFO);
   }

   /* handle output */
   *smallestev = WTMP[0];
   *largestev = WTMP[n - 1];

   /* free memory */
   free(ISUPPZ);
   free(WTMP);
   free(IWORK);
   free(WORK);
}

/** for given A computes A^T A using BLAS */
void SCIPlapackMatrixTransposedMatrixMult(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   double*               matrix,             /**< the matrix we want to multiply */
   double*               result              /**< pointer to store the resulting matrix */
   )
{
   char TRANSA;
   char TRANSB;
   int M;
   int N;
   int K;
   double ALPHA;
   double* A;
   int LDA;
   double* B;
   int LDB;
   double BETA;
   double* C;
   int LDC;

   TRANSA = 'T';
   TRANSB = 'N';
   M = ncols;
   N = ncols;
   K = nrows;
   ALPHA = 1.0;
   A = matrix;
   LDA = K;
   B = matrix;
   LDB = K;
   BETA = 1.0;
   C = result;
   LDC = M;

   F77_FUNC(dgemm, DGEMM)(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
}

/**@} */

