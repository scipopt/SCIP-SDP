/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2019 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2019 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lapack_sdpa.c
 * @brief  interface methods for eigenvector computation and matrix multiplication using openblas
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 * @author Marc Pfetsch
 *
 * This file is used to call the LAPACK routine DSYEVR (double-symmetric-eigenvector computation) and the
 * BLAS routine DGEMV (double-general-matrix-vector multiplication).
 *
 * If a version of openblas/lapack with 64 bit integers is used in connection with SDPA, then one should define
 * LAPACK_LONGLONGINT. Note that one cannot detect the precision with which openblas/lapack is compiled, so this
 * define is essential for correct functioning (resulting in a segmentation fault otherwise).
 */

#include <assert.h>

#include "lapack_interface.h"
#include "config.h"                          /* for F77_FUNC */

#include "scip/def.h"
#include "scip/pub_message.h"                /* for debug and error message */
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/


/* if we use 64 bit integers then use long long int, otherwise int */
#ifdef LAPACK_LONGLONGINT
typedef long long int LAPACKINTTYPE;
#else
typedef int LAPACKINTTYPE;
#endif


/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )

/** transforms a SCIP_Real (that should be integer, but might be off by some numerical error) to an integer by adding an epsilon and rounding down */
#define SCIP_RealTOINT(x) ((LAPACKINTTYPE) (x + 0.5))

/*
 * BLAS/LAPACK Calls
 */

/**@name BLAS/LAPACK Calls */
/**@{ */

/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)(
   char* JOBZ, char* RANGE, char* UPLO,
   LAPACKINTTYPE* N, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* VL, SCIP_Real* VU,
   LAPACKINTTYPE* IL, LAPACKINTTYPE* IU,
   SCIP_Real* ABSTOL, LAPACKINTTYPE* M, SCIP_Real* W, SCIP_Real* Z,
   LAPACKINTTYPE* LDZ, LAPACKINTTYPE* ISUPPZ, SCIP_Real* WORK,
   LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK, LAPACKINTTYPE* LIWORK,
   int* INFO);


/** BLAS Fortran subroutine DGEMV */
void F77_FUNC(dgemv, DGEMV)(char* TRANS, LAPACKINTTYPE* M,
   LAPACKINTTYPE* N, SCIP_Real* ALPHA, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* X, LAPACKINTTYPE* INCX, SCIP_Real* BETA, SCIP_Real* Y, LAPACKINTTYPE* INCY);


/** BLAS Fortran subroutine DGEMM */
void F77_FUNC(dgemm, DGEMM)(char* TRANSA, char* TRANSB, LAPACKINTTYPE* M, LAPACKINTTYPE* N, LAPACKINTTYPE* K, SCIP_Real* ALPHA,
      SCIP_Real* A, LAPACKINTTYPE* LDA, SCIP_Real* B, LAPACKINTTYPE* LDB, SCIP_Real* BETA, SCIP_Real* C, LAPACKINTTYPE* LDC );


/**@} */


/*
 * Functions
 */

/**@name Functions */
/**@{ */

/** computes the i-th eigenvalue of a symmetric matrix using LAPACK, where 1 is the smallest and n the largest, matrix has to be given with all \f$n^2\f$ entries */
SCIP_RETCODE SCIPlapackComputeIthEigenvalue(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to array to store eigenvector */
   )
{
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   char JOBZ;
   char RANGE;
   char UPLO;
   LAPACKINTTYPE LDA;
   SCIP_Real* WORK;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE LIWORK;
   SCIP_Real* WTMP;
   SCIP_Real ABSTOL;
   LAPACKINTTYPE IL;
   LAPACKINTTYPE IU;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   SCIP_Real WSIZE;
   LAPACKINTTYPE WISIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   LAPACKINTTYPE* ISUPPZ;

   assert( bufmem != NULL );
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
#ifdef LAPACK_LONGLONGINT
   LWORK = -1LL;
   LIWORK = -1LL;
#else
   LWORK = -1;
   LIWORK = -1;
#endif

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      &IL, &IU,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO);

   /* for some reason this code seems to be called with INFO=0 within UG */
   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %lld.\n", INFO);
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WTMP, (int) N) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, 2) ); /*lint !e506*/

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
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", INFO);
      return SCIP_ERROR;
   }

   /* handle output */
   *eigenvalue = WTMP[0];

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &WTMP);/*lint !e737*/
   BMSfreeBufferMemoryArray(bufmem, &IWORK);/*lint !e737*/
   BMSfreeBufferMemoryArray(bufmem, &WORK);/*lint !e737*/

   return SCIP_OKAY;
}


/** computes the eigenvector decomposition of a symmetric matrix using LAPACK */
SCIP_RETCODE SCIPlapackComputeEigenvectorDecomposition(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which the decomposition should be computed */
   SCIP_Real*            eigenvalues,        /**< pointer to store eigenvalues (should be length n) */
   SCIP_Real*            eigenvectors        /**< pointer to store eigenvectors (should be length n*n), eigenvectors are given as rows  */
   )
{
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   char JOBZ;
   char RANGE;
   char UPLO;
   LAPACKINTTYPE LDA;
   SCIP_Real* WORK;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE LIWORK;
   SCIP_Real ABSTOL;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   SCIP_Real WSIZE;
   LAPACKINTTYPE WISIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   LAPACKINTTYPE* ISUPPZ;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( A != NULL );
   assert( eigenvalues != NULL );
   assert( eigenvectors != NULL );

   N = n;
   JOBZ = 'V';
   RANGE = 'A';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;
   M = n;
   LDZ = n;

   /* standard LAPACK workspace query, to get the amount of needed memory */
#ifdef LAPACK_LONGLONGINT
   LWORK = -1LL;
   LIWORK = -1LL;
#else
   LWORK = -1;
   LIWORK = -1;
#endif

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      NULL, NULL,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO );

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", INFO);
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, (int) 2 * N) );

   /* call the function */
   VL = -1e20;
   VU = 1e20;

   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      NULL, NULL,
      &ABSTOL, &M, eigenvalues, eigenvectors,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO );

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", INFO);
      return SCIP_ERROR;
   }

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);/*lint !e737*/
   BMSfreeBufferMemoryArray(bufmem, &WORK);/*lint !e737*/

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
   LAPACKINTTYPE M;
   LAPACKINTTYPE N;
   SCIP_Real ALPHA;
   SCIP_Real* A;
   LAPACKINTTYPE LDA;
   SCIP_Real* X;
   LAPACKINTTYPE INCX;
   SCIP_Real BETA;
   SCIP_Real* Y;
   LAPACKINTTYPE INCY;

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


/** performs matrix-matrix-multiplication A*B using BLAS */
SCIP_RETCODE SCIPlapackMatrixMatrixMult(
   int                   nrowsA,             /**< number of rows in matrix A */
   int                   ncolsA,             /**< number of cols in matrix A */
   SCIP_Real*            matrixA,            /**< matrix A given as nrowsA * ncolsA array */
   SCIP_Bool             transposeA,         /**< Should matrix A be transposed before multiplication? */
   int                   nrowsB,             /**< number of rows in matrix B */
   int                   ncolsB,             /**< number of cols in matrix B */
   SCIP_Real*            matrixB,            /**< matrix B given as ncolsA * ncolsB array */
   SCIP_Bool             transposeB,         /**< Should matrix B be transposed before multiplication? */
   SCIP_Real*            result              /**< pointer to nrowsA * nrowsB array to store the resulting matrix */
   )
{
   char TRANSA;
   char TRANSB;
   LAPACKINTTYPE M;
   LAPACKINTTYPE N;
   LAPACKINTTYPE K;
   SCIP_Real ALPHA;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE LDB;
   SCIP_Real BETA;
   LAPACKINTTYPE LDC;

   assert( (transposeA && transposeB && (nrowsA == ncolsB)) || (transposeA && !transposeB && (nrowsA == nrowsB))
      || (!transposeA && transposeB && (ncolsA == ncolsB)) || (!transposeA && !transposeB && (ncolsA == nrowsB)) );

   TRANSA = transposeA ? 'T' : 'N';
   TRANSB = transposeB ? 'T' : 'N';
   M = transposeA ? ncolsA : nrowsA;
   N = transposeB ? nrowsB : ncolsB;
   K = transposeA ? nrowsA : ncolsA;
   ALPHA = 1.0;
   LDA = transposeA ? K : M;
   LDB = transposeB ? N : K;
   BETA = 0.0;
   LDC = M;

   F77_FUNC(dgemm, DGEMM)(&TRANSA, &TRANSB, &M, &N, &K, &ALPHA, matrixA, &LDA, matrixB, &LDB, &BETA, result, &LDC);

   return SCIP_OKAY;
}

/**@} */
