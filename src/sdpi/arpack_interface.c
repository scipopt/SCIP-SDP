/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   arpack_interface.c
 * @brief  interface methods for eigenvector computation using arpack
 * @author Marc Pfetsch
 */

#include <assert.h>

#include "arpack_interface.h"
#include "config.h"                          /* for F77_FUNC */

#include "scip/def.h"
#include "scip/pub_message.h"                /* for debug and error message */
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/


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

/*
 * ARPACK Calls
 */

/**@name ARPACK Calls */
/**@{ */

/** ARPACK Fortran subroutine dsaupd */
void F77_FUNC(dsaupd, DSAUPD)(int* IDO, char* BMAT, int* N, char* WHICH, int* NEV, SCIP_Real* TOL,
   SCIP_Real* RESID, int* NCV, SCIP_Real* V, int* LDV, int* IPARAM, int* IPNTR,
   SCIP_Real* WORKD, SCIP_Real* WORKL, int* LWORKL, int* INFO);

/** ARPACK Fortran subroutine dseupd */
void F77_FUNC(dseupd, DSEUPD)(int* RVEC, char* HOWMNY, int* SELECT, SCIP_Real* D, SCIP_Real* Z,
   int* LDZ, SCIP_Real* SIGMA, char* BMAT, int* N, char* WHICH, int* NEV, SCIP_Real* TOL,
   SCIP_Real* RESID, int* NCV, SCIP_Real* V, int* LDV, int* IPARAM, int* IPNTR,
   SCIP_Real* WORKD, SCIP_Real* WORKL, int* LWORKL, int* INFO);

/**@} */


/*
 * Functions
 */

/**@name Functions */
/**@{ */

#define MAXITER 2147483647

/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK */
SCIP_EXPORT
SCIP_RETCODE SCIParpackComputeSmallestEigenvector(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed in column-major form */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< array for eigenvector */
   )
{  /*lint --e{715}*/
#ifdef ARPACK
   int IPARAM[11];
   int IPNTR[11];
   int SELECT[2];
   int IDO;
   int N;
   int NEV;
   int NCV;
   int LDV;
   int LWORKL;
   int INFO;
   int RVEC;
   SCIP_Real* RESID;
   SCIP_Real* WORKD;
   SCIP_Real* WORKL;
   SCIP_Real* V;
   SCIP_Real TOL;
   SCIP_Real SIGMA;
   char BMAT;
   char WHICH[2];
   char HOWMNY;
   int i;
   int j;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( A != NULL );
   assert( eigenvalue != NULL );
   assert( eigenvector != NULL );

   IDO = 0;
   BMAT = 'I';
   N = n;
   WHICH[0] = 'S'; /* compute smallest algebraic eigenvalue */
   WHICH[1] = 'A';
   NEV = 1;
   TOL = 0.0;
   NCV = MIN(n,4);
   LDV = n;
   IPARAM[0] = 1;       /* exact shifts */
   IPARAM[2] = MAXITER; /* maximal number of iterations */
   IPARAM[6] = 1;       /* Mode 1: A*x = lambda*x, A symmetric, => OP = A  and  B = I. */
   LWORKL = NCV * (NCV + 8);  /* must be at least NCV**2 + 8*NCV */
   INFO = 0;            /* use random starting vector */

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &RESID, n) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &V, n * NCV) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORKD, 3 * n) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORKL, LWORKL) );

   /* enter loop with "reverse communication interface" */
   do
   {
      F77_FUNC(dsaupd, DSAUPD)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
         WORKD, WORKL, &LWORKL, &INFO);

      /* if maximal number of iterations have been reached - still try to extract eigenvalue/eigenvector */
      if ( INFO == 1 )
         break;

      if ( INFO != 0 )
      {
         SCIPerrorMessage("There was an error when calling DSAUPD. INFO = %d.\n", INFO);
         return SCIP_ERROR;
      }

      if ( IDO == 1 || IDO == -1 )
      {
         SCIP_Real* x;
         SCIP_Real* y;

         x = &WORKD[IPNTR[0] - 1];   /* input vector */
         y = &WORKD[IPNTR[1] - 1];   /* output vector */

         /* perform matrix vector multiplication */
         for (i = 0; i < n; ++i)
         {
            SCIP_Real sum = 0.0;
            for (j = 0; j < n; ++j)
               sum += A[i * n + j] * x[j];
            y[i] = sum;
         }
      }
   }
   while ( IDO == -1 || IDO == 1 );

   /* treat possible errors */
   if ( IPARAM[4] <= 0 )
   {
      SCIPerrorMessage("No Ritz vectors have been computed with ARPACK.\n");
      return SCIP_ERROR;
   }
   else if ( IPARAM[2] > MAXITER )
   {
      /* return error if we could not find a Ritz vector and ran into the iteration limit */
      SCIPerrorMessage("Reached iteration limit %d in ARPACK.\n", IPARAM[2]);
      return SCIP_ERROR;
   }

   /* post process */
   RVEC = 1;   /* compute eigenvectors */
   HOWMNY = 'A';
   SIGMA = 0.0;
   INFO = 0;

   F77_FUNC(dseupd, DSEUPD)(&RVEC, &HOWMNY, SELECT, eigenvalue, V, &LDV, &SIGMA,
      &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
      WORKD, WORKL, &LWORKL, &INFO);

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSEUPD. INFO = %d.\n", INFO);
      return SCIP_ERROR;
   }

   /* copy handle output */
   for (i = 0; i < n; ++i)
      eigenvector[i] = V[i];

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &WORKL);
   BMSfreeBufferMemoryArray(bufmem, &WORKD);
   BMSfreeBufferMemoryArray(bufmem, &V);
   BMSfreeBufferMemoryArray(bufmem, &RESID);
#endif

   return SCIP_OKAY;
}


/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK; specialized sparse version for mu A - B */
SCIP_EXPORT
SCIP_RETCODE SCIParpackComputeSmallestEigenvectorOneVar(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real             mu,                 /**< scaling factor for A matrix */
   int                   annonz,             /**< number of nonzero elements in A */
   int*                  arow,               /**< array of row-indices of A */
   int*                  acol,               /**< array of column-indices of A */
   SCIP_Real*            aval,               /**< array of nonzero values of entries of A */
   int                   bnnonz,             /**< number of nonzero elements in B */
   int*                  brow,               /**< array of row-indices of nonzero matrix entries in B */
   int*                  bcol,               /**< array of column-indices of nonzero matrix entries in B*/
   SCIP_Real*            bval,               /**< array of nonzero values in B */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< array for eigenvector */
   )
{  /*lint --e{715}*/
#ifdef ARPACK
   int IPARAM[11];
   int IPNTR[11];
   int SELECT[2];
   int IDO;
   int N;
   int NEV;
   int NCV;
   int LDV;
   int LWORKL;
   int INFO;
   int RVEC;
   SCIP_Real* RESID;
   SCIP_Real* WORKD;
   SCIP_Real* WORKL;
   SCIP_Real* V;
   SCIP_Real TOL;
   SCIP_Real SIGMA;
   char BMAT;
   char WHICH[2];
   char HOWMNY;
   int i;

   assert( bufmem != NULL );
   assert( n > 0 );
   assert( n < INT_MAX );
   assert( arow != NULL );
   assert( acol != NULL );
   assert( aval != NULL );
   assert( brow != NULL );
   assert( bcol != NULL );
   assert( bval != NULL );
   assert( eigenvalue != NULL );
   assert( eigenvector != NULL );

   IDO = 0;
   BMAT = 'I';
   N = n;
   WHICH[0] = 'S'; /* compute smallest algebraic eigenvalue */
   WHICH[1] = 'A';
   NEV = 1;
   TOL = 0.0;
   NCV = MIN(n,4);
   LDV = n;
   IPARAM[0] = 1;       /* exact shifts */
   IPARAM[2] = MAXITER; /* maximal number of iterations */
   IPARAM[6] = 1;       /* Mode 1: A*x = lambda*x, A symmetric, => OP = A  and  B = I. */
   LWORKL = NCV * (NCV + 8);  /* must be at least NCV**2 + 8*NCV */
   INFO = 0;            /* use random starting vector */

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &RESID, n) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &V, n * NCV) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORKD, 3 * n) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &WORKL, LWORKL) );

   /* enter loop with "reverse communication interface" */
   do
   {
      F77_FUNC(dsaupd, DSAUPD)(&IDO, &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
         WORKD, WORKL, &LWORKL, &INFO);

      /* if maximal number of iterations have been reached - still try to extract eigenvalue/eigenvector */
      if ( INFO == 1 )
         break;

      if ( INFO != 0 )
      {
         SCIPerrorMessage("There was an error when calling DSAUPD. INFO = %d.\n", INFO);
         return SCIP_ERROR;
      }

      if ( IDO == 1 || IDO == -1 )
      {
         SCIP_Real* x;
         SCIP_Real* y;
         int r;
         int c;

         x = &WORKD[IPNTR[0] - 1];   /* input vector */
         y = &WORKD[IPNTR[1] - 1];   /* output vector */

         /* perform matrix vector multiplication */
         for (i = 0; i < n; ++i)
            y[i] = 0.0;

         for (i = 0; i < annonz; ++i)
         {
            r = arow[i];
            c = acol[i];
            assert( 0 <= r && r < n );
            assert( 0 <= c && c < n );
            y[r] += mu * aval[i] * x[c];
            if ( r != c )
               y[c] += mu * aval[i] * x[r];
         }

         for (i = 0; i < bnnonz; ++i)
         {
            r = brow[i];
            c = bcol[i];
            assert( 0 <= r && r < n );
            assert( 0 <= c && c < n );
            y[r] -= bval[i] * x[c];
            if ( r != c )
               y[c] -= bval[i] * x[r];
         }
      }
   }
   while ( IDO == -1 || IDO == 1 );

   /* treat possible errors */
   if ( IPARAM[4] <= 0 )
   {
      SCIPerrorMessage("No Ritz vectors have been computed with ARPACK.\n");
      return SCIP_ERROR;
   }
   else if ( IPARAM[2] > MAXITER )
   {
      /* return error if we could not find a Ritz vector and ran into the iteration limit */
      SCIPerrorMessage("Reached iteration limit %d in ARPACK.\n", IPARAM[2]);
      return SCIP_ERROR;
   }

   /* post process */
   RVEC = 1;   /* compute eigenvectors */
   HOWMNY = 'A';
   SIGMA = 0.0;
   INFO = 0;

   F77_FUNC(dseupd, DSEUPD)(&RVEC, &HOWMNY, SELECT, eigenvalue, V, &LDV, &SIGMA,
      &BMAT, &N, WHICH, &NEV, &TOL, RESID, &NCV, V, &LDV, IPARAM, IPNTR,
      WORKD, WORKL, &LWORKL, &INFO);

   if ( INFO != 0 )
   {
      SCIPerrorMessage("There was an error when calling DSEUPD. INFO = %d.\n", INFO);
      return SCIP_ERROR;
   }

   /* copy handle output */
   for (i = 0; i < n; ++i)
      eigenvector[i] = V[i];

   /* free memory */
   BMSfreeBufferMemoryArray(bufmem, &WORKL);
   BMSfreeBufferMemoryArray(bufmem, &WORKD);
   BMSfreeBufferMemoryArray(bufmem, &V);
   BMSfreeBufferMemoryArray(bufmem, &RESID);
#endif

   return SCIP_OKAY;
}

/**@} */
