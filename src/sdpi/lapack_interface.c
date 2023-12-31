/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2023 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2023 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   lapack_interface.c
 * @brief  interface methods for eigenvector computation and matrix multiplication using openblas
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 * @author Marc Pfetsch
 *
 * This file is used to call the LAPACK routine DSYEVR (double-symmetric-eigenvector computation) and the
 * BLAS routine DGEMV (double-general-matrix-vector multiplication).
 *
 * LAPACK can be built with 32- or 64-bit integers, which is not visible to the outside. This interface tries to work
 * around this issue. Since the Fortran routines are called by reference, they only get a pointer. We always use 64-bit
 * integers on input, but reduce the output to 32-bit integers. We assume that all sizes can be represented in 32-bit
 * integers.
 */

/* #define PRINTMATRICES     /\* Should all matrices appearing in best rank-1 approximation heuristic be printed? *\/ */

#include <assert.h>

#include "lapack_interface.h"
#include "config.h"                          /* for F77_FUNC */

#include "scip/def.h"
#include "scip/pub_message.h"                /* for debug and error message */
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/* use int type from Openblas if available */
#ifdef OPENBLAS
#include <cblas.h>
typedef blasint LAPACKINTTYPE;
#else
/* otherwise we use int as a default */
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

/** transforms a SCIP_Real (that should be integer, but might be off by some numerical error) to an integer by adding 0.5 and rounding down */
#define SCIP_RealTOINT(x) ((LAPACKINTTYPE) (x + 0.5))

/*
 * BLAS/LAPACK Calls
 */

/**@name BLAS/LAPACK Calls */
/**@{ */

/** LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevr, DSYEVR)(char* JOBZ, char* RANGE, char* UPLO,
   LAPACKINTTYPE* N, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* VL, SCIP_Real* VU,
   LAPACKINTTYPE* IL, LAPACKINTTYPE* IU,
   SCIP_Real* ABSTOL, LAPACKINTTYPE* M, SCIP_Real* W, SCIP_Real* Z,
   LAPACKINTTYPE* LDZ, LAPACKINTTYPE* ISUPPZ, SCIP_Real* WORK,
   LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK, LAPACKINTTYPE* LIWORK,
   LAPACKINTTYPE* INFO);

/** Alternative LAPACK Fortran subroutine DSYEVR */
void F77_FUNC(dsyevx, DSYEVX)(char* JOBZ, char* RANGE, char* UPLO,
   LAPACKINTTYPE* N, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* VL, SCIP_Real* VU,
   LAPACKINTTYPE* IL, LAPACKINTTYPE* IU,
   SCIP_Real* ABSTOL, LAPACKINTTYPE* M, SCIP_Real* W, SCIP_Real* Z,
   LAPACKINTTYPE* LDZ, SCIP_Real* WORK, LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK,
   LAPACKINTTYPE* IFAIL, LAPACKINTTYPE* INFO);

/** BLAS Fortran subroutine DGEMV */
void F77_FUNC(dgemv, DGEMV)(char* TRANS, LAPACKINTTYPE* M,
   LAPACKINTTYPE* N, SCIP_Real* ALPHA, SCIP_Real* A, LAPACKINTTYPE* LDA,
   SCIP_Real* X, LAPACKINTTYPE* INCX, SCIP_Real* BETA, SCIP_Real* Y, LAPACKINTTYPE* INCY);

/** BLAS Fortran subroutine DGEMM */
void F77_FUNC(dgemm, DGEMM)(char* TRANSA, char* TRANSB, LAPACKINTTYPE* M, LAPACKINTTYPE* N, LAPACKINTTYPE* K, SCIP_Real* ALPHA,
      SCIP_Real* A, LAPACKINTTYPE* LDA, SCIP_Real* B, LAPACKINTTYPE* LDB, SCIP_Real* BETA, SCIP_Real* C, LAPACKINTTYPE* LDC);

/* LAPACK Fortran subroutine DGELSD */
void F77_FUNC(dgelsd, DGELSD)(LAPACKINTTYPE* M, LAPACKINTTYPE* N, LAPACKINTTYPE* NRHS,
      SCIP_Real* A, LAPACKINTTYPE* LDA, SCIP_Real* b, LAPACKINTTYPE* LDB, SCIP_Real* S, SCIP_Real* RCOND, LAPACKINTTYPE* RANK,
      SCIP_Real* WORK, LAPACKINTTYPE* LWORK, LAPACKINTTYPE* IWORK, LAPACKINTTYPE* INFO);

/**@} */


/** converts a number stored in a long long int to an int, depending on big- or little endian machines
 *
 *  We assume that the number actually fits into an int. Thus, if more bits are used, we assume that the number is
 *  negative.
 */
static
int convertToInt(
   long long int         num                 /**< number to be converted */
   )
{
   long long int work;
   int checkval = 1;

   assert(sizeof(work) > sizeof(checkval)); /*lint !e506*/

   /* if we have a little-endian machine (e.g, x86), the sought value is in the bottom part */
   if ( *(int8_t*)&checkval != 0 ) /*lint !e774*/
   {
      /* if the top part is nonzero, we assume that the number is negative */
      if ( *((int8_t*)&num + 4) != 0 ) /*lint !e2662*/
      {
         work = -num;
         return -(*((int*)&work));
      }
      return *((int*)&num);
   }

   /* otherwise we have a big-endian machine (e.g., PowerPC); the sought value is in the top part */
   assert( *(int8_t*)&checkval == 0 );

   /* if the bottom part is nonzero, we assume that the number is negative */
   if ( *(int8_t*)&num != 0 ) /*lint !e774*/
   {
      work = -num;
      return -(*((int*)&work + 4)); /*lint !e2662*/
   }
   return *((int*)&num + 4);
}


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
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed - will be destroyed! */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to array to store eigenvector */
   )
{
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE* ISUPPZ;
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE WISIZE;
   LAPACKINTTYPE IL;
   LAPACKINTTYPE IU;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE LIWORK;
   SCIP_Real* WORK;
   SCIP_Real* WTMP;
   SCIP_Real ABSTOL;
   SCIP_Real WSIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   char JOBZ;
   char RANGE;
   char UPLO;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( 0 < i && i <= n );
   assert( eigenvalue != NULL );
   assert( ! geteigenvectors || eigenvector != NULL );

   N = n;
   JOBZ = geteigenvectors ? 'V' : 'N';
   RANGE = 'I';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0; /* we use abstol = 0, since some lapack return an error otherwise */
   VL = -1e20;
   VU = 1e20;
   IL = i;
   IU = i;
   M = 1;
   LDZ = n;
   INFO = 0LL;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;
   LIWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      &IL, &IU,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
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
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      &IL, &IU,
      &ABSTOL, &M, WTMP, eigenvector,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO);

   /* handle output */
   if ( convertToInt(INFO) == 0 )
     *eigenvalue = WTMP[0];

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &WTMP);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** computes i-th eigenvalue of a symmetric matrix using alternative algorithm in LAPACK, where 1 is the smallest and n the largest, matrix has to be given with all \f$n^2\f$ entries */
SCIP_RETCODE SCIPlapackComputeIthEigenvalueAlternative(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed - will be destroyed! */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to array to store eigenvector */
   )
{
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE* IFAIL;
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE IL;
   LAPACKINTTYPE IU;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   LAPACKINTTYPE LWORK;
   SCIP_Real* WORK;
   SCIP_Real* WTMP;
   SCIP_Real WSIZE;
   SCIP_Real ABSTOL;
   SCIP_Real VL;
   SCIP_Real VU;
   char JOBZ;
   char RANGE;
   char UPLO;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( 0 < i && i <= n );
   assert( eigenvalue != NULL );
   assert( ! geteigenvectors || eigenvector != NULL );

   N = n;
   JOBZ = geteigenvectors ? 'V' : 'N';
   RANGE = 'I';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;  /* we use abstol = 0, since some lapack return an error otherwise */
   VL = -1e20;
   VU = 1e20;
   IL = i;
   IU = i;
   M = 1;
   LDZ = n;
   INFO = 0LL;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE */
   F77_FUNC(dsyevx, DSYEVX)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      &IL, &IU,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, &WSIZE, &LWORK, NULL, NULL,
      &INFO);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) 5 * N) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WTMP, (int) N) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IFAIL, (int) N) ); /*lint !e506*/

   /* call the function */
   F77_FUNC(dsyevx, DSYEVX)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      &IL, &IU,
      &ABSTOL, &M, WTMP, eigenvector,
      &LDZ, WORK, &LWORK, IWORK, IFAIL,
      &INFO);

   /* handle output */
   if ( convertToInt(INFO) == 0 )
     *eigenvalue = WTMP[0];

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &IFAIL);
   BMSfreeBufferMemoryArray(bufmem, &WTMP);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVX. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** computes eigenvectors corresponding to negative eigenvalues of a symmetric matrix using LAPACK, matrix has to be given with all \f$n^2\f$ entries */
SCIP_RETCODE SCIPlapackComputeEigenvectorsNegative(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvectors should be computed - will be destroyed! */
   SCIP_Real             tol,                /**< tolerance; the eigenvalues will be in the interval (-1e20, -tol] */
   int*                  neigenvalues,       /**< pointer to store the number of negative eigenvalues */
   SCIP_Real*            eigenvalues,        /**< array for eigenvalues (should be length n) */
   SCIP_Real*            eigenvectors        /**< array for eigenvectors (should be length n*n), eigenvectors are given as rows  */
   )
{
   LAPACKINTTYPE* ISUPPZ;
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE LIWORK;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   LAPACKINTTYPE WISIZE;
   SCIP_Real* WORK;
   SCIP_Real ABSTOL;
   SCIP_Real WSIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   char JOBZ;
   char RANGE;
   char UPLO;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( tol >= 0 );
   assert( neigenvalues != NULL );
   assert( eigenvalues != NULL );
   assert( eigenvectors != NULL );

   N = n;
   JOBZ = 'V';
   RANGE = 'V';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;  /* we use abstol = 0, since some lapack return an error otherwise */
   LDZ = n;
   M = -1;
   INFO = 0LL;

   /* interval of allowed values */
   VL = -1e30;
   VU = -tol;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;
   LIWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      &VL, &VU,
      NULL, NULL,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO);

   /* for some reason this code seems to be called with INFO=0 within UG */
   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, (int) 2 * N) );

   /* call the function */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      NULL, NULL,
      &ABSTOL, &M, eigenvalues, eigenvectors,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO);

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   *neigenvalues = convertToInt(M);

   return SCIP_OKAY;
}


/** computes the eigenvector decomposition of a symmetric matrix using LAPACK */
SCIP_RETCODE SCIPlapackComputeEigenvectorDecomposition(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which the decomposition should be computed - will be destroyed! */
   SCIP_Real*            eigenvalues,        /**< pointer to store eigenvalues (should be length n) */
   SCIP_Real*            eigenvectors        /**< pointer to store eigenvectors (should be length n*n), eigenvectors are given as rows  */
   )
{
   LAPACKINTTYPE* ISUPPZ;
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE N;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE LIWORK;
   LAPACKINTTYPE M;
   LAPACKINTTYPE LDZ;
   LAPACKINTTYPE WISIZE;
   SCIP_Real* WORK;
   SCIP_Real ABSTOL;
   SCIP_Real WSIZE;
   SCIP_Real VL;
   SCIP_Real VU;
   char JOBZ;
   char RANGE;
   char UPLO;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( eigenvalues != NULL );
   assert( eigenvectors != NULL );

   N = n;
   JOBZ = 'V';
   RANGE = 'A';
   UPLO = 'L';
   LDA  = n;
   ABSTOL = 0.0;  /* we use abstol = 0, since some lapack return an error otherwise */
   M = n;
   LDZ = n;
   VL = -1e20;
   VU = 1e20;
   INFO = 0LL;

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;
   LIWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entries of [the 1x1 arrays]) WSIZE and WISIZE */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, NULL, &LDA,
      NULL, NULL,
      NULL, NULL,
      &ABSTOL, &M, NULL, NULL,
      &LDZ, NULL, &WSIZE,
      &LWORK, &WISIZE, &LIWORK,
      &INFO);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &ISUPPZ, (int) 2 * N) );

   /* call the function */
   F77_FUNC(dsyevr, DSYEVR)( &JOBZ, &RANGE, &UPLO,
      &N, A, &LDA,
      &VL, &VU,
      NULL, NULL,
      &ABSTOL, &M, eigenvalues, eigenvectors,
      &LDZ, ISUPPZ, WORK,
      &LWORK, IWORK, &LIWORK,
      &INFO);

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &ISUPPZ);
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSYEVR. INFO = %d.\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

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
   LAPACKINTTYPE M;
   LAPACKINTTYPE N;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE INCX;
   LAPACKINTTYPE INCY;
   SCIP_Real* A;
   SCIP_Real* X;
   SCIP_Real* Y;
   SCIP_Real BETA;
   SCIP_Real ALPHA;
   char TRANS;

   assert( nrows > 0 );
   assert( nrows < INT_MAX );
   assert( ncols > 0 );
   assert( ncols < INT_MAX );
   assert( matrix != NULL );
   assert( vector != NULL );
   assert( result != NULL );

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
   LAPACKINTTYPE M;
   LAPACKINTTYPE N;
   LAPACKINTTYPE K;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE LDB;
   LAPACKINTTYPE LDC;
   SCIP_Real ALPHA;
   SCIP_Real BETA;
   char TRANSA;
   char TRANSB;

   assert( nrowsA > 0 );
   assert( nrowsA < INT_MAX );
   assert( ncolsA > 0 );
   assert( ncolsA < INT_MAX );
   assert( nrowsB > 0 );
   assert( nrowsB < INT_MAX );
   assert( ncolsB > 0 );
   assert( ncolsB < INT_MAX );
   assert( matrixA != NULL );
   assert( matrixB != NULL );
   assert( result != NULL );

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


/** computes the minimum-norm solution to a real linear least squares problem: minimize 2-norm(| b - A*x |) using LAPACK
 *  (uses singular value decomposition of A). A is an M-by-N matrix which may be rank-deficient.
 */
SCIP_RETCODE SCIPlapackLinearSolve(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   m,                  /**< number of rows of A */
   int                   n,                  /**< number of columns of A */
   SCIP_Real*            A,                  /**< coefficient matrix of the linear system */
   SCIP_Real*            b,                  /**< right-hand side of the linear system (should be length max(m,n)) */
   SCIP_Real*            x                   /**< pointer to store values for x (should be length n) */
   )
{
   LAPACKINTTYPE* IWORK;
   LAPACKINTTYPE M;
   LAPACKINTTYPE N;
   LAPACKINTTYPE NRHS;
   LAPACKINTTYPE LDA;
   LAPACKINTTYPE LDB;
   LAPACKINTTYPE RANK;
   LAPACKINTTYPE LWORK;
   LAPACKINTTYPE LIWORK;
   LAPACKINTTYPE INFO;
   LAPACKINTTYPE WISIZE;
   SCIP_Real* WORK;
   SCIP_Real* S;
   SCIP_Real RCOND;
   SCIP_Real WSIZE;
   int i;

   assert( bufmem != NULL );
   assert( m > 0 );
   assert( m < INT_MAX );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( b != NULL );
   assert( x != NULL );

   M = m;
   N = n;
   NRHS = 1;
   LDA = m;
   LDB = MAX(M,N);
   RCOND = 0.0;
   INFO = 0LL;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &S, MIN(m,n)) );

   /* standard LAPACK workspace query, to get the amount of needed memory */
   LWORK = -1LL;

   /* this computes the internally needed memory and returns this as (the first entry of [the 1x1 array]) WSIZE */
   F77_FUNC(dgelsd, DGELSD)( &M, &N, &NRHS,
      NULL, &LDA, NULL, &LDB, NULL,
      &RCOND, &RANK, &WSIZE, &LWORK,
      &WISIZE, &INFO);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DGELSD. INFO = %d\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* allocate workspace */
   LWORK = SCIP_RealTOINT(WSIZE);
   LIWORK = WISIZE;

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORK, (int) LWORK) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &IWORK, (int) LIWORK) );

   /* call the function */
   F77_FUNC(dgelsd, DGELSD)( &M, &N, &NRHS,
      A, &LDA, b, &LDB, S, &RCOND, &RANK,
      WORK, &LWORK, IWORK, &INFO);

#ifdef PRINTMATRICES
   {
      SCIP_Real residual = 0.0;

      printf("LWORK = %d\n", LWORK);
      printf("LIWORK = %d\n", LIWORK);
      printf("A has size (%d,%d), is of rank %d\n", M, N, RANK);
      printf("Minimum l2-norm solution of linear equation system:\n");

      for (i = 0; i < n; ++i)
         printf("(%d, %f)   ", i, b[i]);

      for (i = n; i < m; ++i)
         residual += b[i] * b[i];
      printf("\n");

      printf("Residual sum-of-squares for the solution is %f\n", residual);
   }
#endif

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &IWORK);
   BMSfreeBufferMemoryArray(bufmem, &WORK);
   BMSfreeBufferMemoryArray(bufmem, &S);

   if ( convertToInt(INFO) != 0 )
   {
      SCIPerrorMessage("There was an error when calling DGELSD. INFO = %d\n", convertToInt(INFO));
      return SCIP_ERROR;
   }

   /* LAPACK overwrites the right-hand side with the result */
   for (i = 0; i < n; ++i)
      x[i] = b[i];

   return SCIP_OKAY;
}

/**@} */
