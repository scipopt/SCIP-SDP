#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>
#include "lapack.h"

long long int counter;
int timelimitprinted = 0;


/** Computes smallest and largest eigenvalue for any submatrix of given size */
static
void computeExtremalEigenvaluesForSubmatrices(
   int                   kused,              /* number of indices already chosen */
   int                   kopen,              /* number of indices still to be chosen */
   int*                  chosen,             /* list of already chosen indices (length = kused) */
   int                   smallest,           /* smallest possible index */
   int                   largest,            /* largest possible index */
   double*               matrix,             /* matrix to compute eigenvalues for */
   double*               allocedmem,         /* allocated memory of size (kused+kopen)^2 to save submatrix (will be overwritten) */
   double*               minev,              /* smallest eigenvalue of any submatrix of size at most kused+kopen */
   double*               maxev,              /* largest eigenvalue of any submatrix of size at most kused+kopen */
   double                timelimit           /* timeofday in seconds when the timelimit is reached (may be negative for no timelimit) */
   )
{
   int i;
   int j;
   struct timeval currenttime;
   double currentseconds;


   assert( kused >= 0 );
   assert( kopen >= 0 );
   assert( kused == 0 || chosen != NULL );
   assert( smallest >= 0 );
   assert( largest >= smallest || kopen == 0 );
   assert( matrix != NULL );
   assert( allocedmem != NULL || kused + kopen == 0 );

   *minev = INFINITY;
   *maxev = -INFINITY;

   /* check the current time, if we are behind the timelimit: return */
   if ( timelimit >= 0 )
   {
      gettimeofday(&currenttime, NULL);
      currentseconds = (double) currenttime.tv_sec + (double) currenttime.tv_usec / 1e6;

      if ( currentseconds >= timelimit )
      {
         if ( ! timelimitprinted )
         {
            timelimitprinted = 1;
            printf("time limit reached ! \n");
         }
         return;
      }
   }

   if ( kopen == 0 )
   {
      /* all indices are chosen, now the smallest/largest eigenvalue is computed */

      /* assemble the submatrix */
      for (i = 0; i < kused; i++)
      {
         for (j = 0; j < kused; j++)
            allocedmem[i * kused + j] = matrix[chosen[i] * (largest + 1) + chosen[j]];
      }

      /* compute the eigenvalues */
      SCIPlapackComputeSmallestLargestEigenvalue(kused, allocedmem, minev, maxev);
      counter++;

      return;
   }
   else
   {
      double smallestev;
      double largestev;

      for (i = smallest; i <= largest - kopen + 1; i++)
      {
         int* newchosen;

         /* we need to allocate new memory for every recursive call */
         newchosen = (int*) malloc((kused + 1) * sizeof(int));
         for (j = 0; j < kused; j++)
            newchosen[j] = chosen[j];
         newchosen[kused] = i;

         /* recursive call */
         computeExtremalEigenvaluesForSubmatrices(kused + 1, kopen - 1, newchosen, i + 1, largest, matrix, allocedmem, &smallestev, &largestev, timelimit);

         /* compare values, update optimum */
         if ( smallestev < *minev )
            *minev = smallestev;
         if ( largestev > *maxev )
            *maxev = largestev;
         free(newchosen);
      }
      return;
   }
}

/** main function
 * first argument: path to file
 * second (optional) argument: timelimit in seconds */
int main (
   int argc,
   char** argv
   )
{
   double* matrix;
   int matrixlength;
   int m;
   int n;
   int k;
   int i;
   int j;
   FILE* file;
   double* ATA; /* the matrix A^T A, which we will compute eigenvalues for */
   double* submatrix; /* submatrix of ATA for indices in given sets, needs to be rewritten each times cause dsyevr overwrites its entries */
   double minev;
   double maxev;
   double delta;
   double gamma;
   struct timeval starttime;
   struct timeval finishtime;
   double startseconds;
   double finishseconds;
   double timelimit;

   counter = 0;

   if ( argc <= 1 )
   {
      printf("please call this with the matrix file as input");
      return -1;
   }

   file = fopen(argv[1], "r");

   fscanf(file, "%*[^\n]\n"); /* skip the first line */

   fscanf(file, "%*c %*c %d %*c %*c %d %*c %*c %d %*[^\n] \n", &m, &n, &k);

   printf("m = %d, n = %d, k = %d \n", m, n, k);

   matrixlength = m * n;

   matrix = (double*) malloc(matrixlength * sizeof(double));

   for (i = 0; i < m; i++)
   {
      for (j = 0; j < n-1; j++)
         fscanf(file, "%lf \n", &matrix[j*m+i]);
      fscanf(file, "%lf \n", &matrix[j*m+i]); /* also read/remove the newline */
   }

   gettimeofday(&starttime, NULL);
   startseconds = (double) starttime.tv_sec + (double) starttime.tv_usec / 1e6;

   ATA = (double*) malloc(n*n*sizeof(double));

   SCIPlapackMatrixTransposedMatrixMult(m, n, matrix, ATA);

   free(matrix);

   /* get timelimit */
   if ( argc >= 3 )
      timelimit = startseconds + strtod(argv[2], NULL);
   else
      timelimit = -1;

   /* compute eigenvalues */
   submatrix = (double*) malloc(k * k * sizeof(double));

   computeExtremalEigenvaluesForSubmatrices(0, k, NULL, 0, n-1, ATA, submatrix, &minev, &maxev, timelimit);

   free(submatrix);
   free(ATA);

   fclose(file);

   delta = (1 - minev) > (maxev - 1) ? 1 - minev : maxev - 1;
   gamma = maxev / minev;

   gettimeofday(&finishtime, NULL);
   finishseconds = (double) finishtime.tv_sec + (double) finishtime.tv_usec / 1e6;

   printf("minev = %f, maxev = %f, delta = %f, gamma = %f\n", minev, maxev, delta, gamma);

   printf("eigenvalue computations = %lld, seconds = %f\n", counter, finishseconds - startseconds);

   return 0;
}
