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

/**@file   cons_sdp.c
 * @brief  constraint handler for sdp-constraints
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 */

/* #define SCIP_DEBUG*/
/* #define SCIP_MORE_DEBUG *//* shows all cuts added */

#include "cons_sdp.h"

#include <assert.h>                     /*lint !e451*/
#include <string.h>                     /* for NULL, strcmp */

#include "sdpi/lapack.h"

#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/SdpVarfixer.h"

#include "scip/cons_linear.h"           /* for SCIPcreateConsLinear */
#include "scip/scip.h"                  /* for SCIPallocBufferArray, etc */

#define CONSHDLR_NAME          "SDP"
#define CONSHDLR_DESC          "SDP constraints of the form \\sum_{j} A_j y_j - A_0 psd"
#define CONSHDLR_SEPAPRIORITY  +1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */


/** constraint data for sdp constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in this SDP constraint */
   int                   nnonz;              /**< number of nonzeroes in this SDP constraint */
   int                   blocksize;          /**< size of this SDP-block */
   int*                  nvarnonz;           /**< length of the arrays pointed to by col/row/val, number of nonzeros for each variable */
   int**                 col;                /**< pointers to the column indices of the nonzeroes for each variable */
   int**                 row;                /**< pointers to the row indices of the nonzeroes for each variable */
   SCIP_Real**           val;                /**< pointers to the values of the nonzeroes for each variable */
   SCIP_VAR**            vars;               /**< SCIP_VARiables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz;         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol;           /**< column indices of the constant nonzeroes */
   int*                  constrow;           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval;           /**< values of the constant nonzeroes */
   SCIP_Real             maxrhsentry;        /**< maximum entry of constant matrix (needed for DIMACS error norm) */
};

/** SDP constraint handler data */
struct SCIP_ConshdlrData
{
   int                   neigveccuts;        /**< this is used to give the eigenvector-cuts distinguishable names */
   int                   ndiaggezerocuts;    /**< this is used to give the diagGEzero-cuts distinguishable names */
   int                   ndiagdomcuts;       /**< this is used to give the diagDominant-cuts distinguishable names */
   int                   n1x1blocks;         /**< this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */
};

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
   int ind = 0;

   assert( size >= 0 );
   assert( symMat != NULL );
   assert( fullMat != NULL );

   /* traverse the lower triangular part in the order of the indices and copy the values to both lower and upper triangular part */
   for (i = 0; i < size; i++)
   {
      for (j = 0; j <= i; j++)
      {
         assert( ind == compLowerTriangPos(i,j) );
         fullMat[i*size + j] = symMat[ind]; /*lint !e679*/
         fullMat[j*size + i] = symMat[ind]; /*lint !e679*/
         ind++;
      }
   }

   return SCIP_OKAY;
}

/** For a given vector \f$ y \f$ computes the SDP-Matrix \f$ \sum_{j=1}^m A_j y_j - A_0 \f$ for this SDP block.
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
   int blocksize;
   SCIP_Real yval;

   assert( cons != NULL );
   assert( matrix != NULL );

   consdata = SCIPconsGetData(cons);
   nvars = consdata->nvars;
   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1))/2; i++)
      matrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      yval = SCIPgetSolVal(scip, y, consdata->vars[i]);
      for (ind = 0; ind < consdata->nvarnonz[i]; ind++)
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
   SCIP_CONS*            cons,               /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   int                   j,                  /**< variable-index of the matrix to multiply with */
   SCIP_Real*            v,                  /**< vector to multiply with */
   SCIP_Real*            vAv                 /**< the resulting scalar \f$ v^T A_j v \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert( cons != NULL );
   assert( j >= 0 );
   assert( vAv != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   assert( j < consdata->nvars );

   /* initialize the product with 0 */
   *vAv = 0.0;

   for (i = 0; i < consdata->nvarnonz[j]; i++)
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

/** Set the maximum absolute value of an entry of the constant matrix.
 *  This must be done before presolving, because otherwise this is influenced by variable fixings (which might lead
 *  to solutions being feasible in presolving no longer being feasible afterwards) */
static
SCIP_RETCODE setMaxRhsEntry(
   SCIP_CONS*            cons                /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real max;
   int i;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* initialize max with zero (this is used if there is no constant-matrix) */
   max = 0.0;

   /* iterate over the entries of the constant matrix, updating max if a higher absolute value is found */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      if ( REALABS(consdata->constval[i]) > max )
         max = REALABS(consdata->constval[i]);
   }

   consdata->maxrhsentry = max;

   return SCIP_OKAY;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize+1))/2 ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullconstmatrix, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvector, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &output_vector, blocksize) );

   /* compute the matrix \f$ \sum_j A_j y_j - A_0 \f$ */
   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix) );

   /* expand it because LAPACK wants the full matrix instead of the lower triangular part */
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPblkmem(scip), TRUE, blocksize, fullmatrix, 1, eigenvalues, eigenvector) );

   /* get full constant matrix */
   SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, fullconstmatrix) );

   /* multiply eigenvector with constant matrix to get lhs (after multiplying again with eigenvector from the left) */
   SCIP_CALL( SCIPlapackMatrixVectorMult(blocksize, blocksize, fullconstmatrix, eigenvector, output_vector) );

   for (j = 0; j < blocksize; ++j)
      *lhs += eigenvector[j] * output_vector[j];

   /* compute \f$ v^T A_j v \f$ for eigenvector v and each matrix \f$ A_j \f$ to get the coefficients of the LP cut */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_CALL( multiplyConstraintMatrix(cons, j, eigenvector, &coeff[j]) );
   }

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
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**< the solution to check feasibility for */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the feasibility checking call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int blocksize;
   double check_value;
   SCIP_Real eigenvalue;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( result != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), "SDP") == 0);

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize+1)) / 2) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize) );
   SCIP_CALL( computeSdpMatrix(scip, cons, sol, matrix) );
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPblkmem(scip), FALSE, blocksize, fullmatrix, 1, &eigenvalue, NULL) );

   /* This enables checking the second DIMACS Error Norm: err=max{0, -lambda_min(x)/(1+maximumentry of rhs} */
#ifdef DIMACS
   check_value = (-eigenvalue) / (1.0 + consdata->maxrhsentry);
#else
   check_value = -eigenvalue;
#endif

   if ( SCIPisFeasLE(scip, check_value, 0.0) )
      *result = SCIP_FEASIBLE;
   else
   {
      *result = SCIP_INFEASIBLE;
#ifdef SCIP_DEBUG
      if ( printreason )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
         SCIPinfoMessage(scip, NULL, "non psd matrix (eigenvalue %f, dimacs error norm = %f).\n", eigenvalue, check_value);
      }
#endif
   }

   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** separates the current solution */
static
SCIP_RETCODE separateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real lhs = 0.0;
   SCIP_Real* coeff = NULL;
   int nvars;
   char* cutname;
#ifndef NDEBUG
   int snprintfreturn; /* this is used to assert that the SCIP string concatenation works */
#endif
   int j;
   SCIP_COL** cols;
   SCIP_Real* vals;
   int len;
   SCIP_ROW* row;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;
   SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars ) );

   SCIP_CALL( cutUsingEigenvector(scip, cons, sol, coeff, &lhs) );

   SCIP_CALL( SCIPallocBufferArray(scip, &cols, nvars ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars ) );

   len = 0;
   for (j = 0; j < nvars; ++j)
   {
      if ( SCIPisZero(scip, coeff[j]) )
         continue;

      cols[len] = SCIPvarGetCol(SCIPvarGetProbvar(consdata->vars[j]));
      vals[len] = coeff[j];
      ++len;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255) );
#ifndef NDEBUG
   snprintfreturn = SCIPsnprintf(cutname, 255, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
   assert( snprintfreturn < 256 ); /* it returns the number of positions that would have been needed, if that is more than 255, it failed */
#else
   SCIPsnprintf(cutname, 255, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
#endif
   SCIP_CALL( SCIPcreateRowCons(scip, &row, conshdlr, cutname , len, cols, vals, lhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );

   if ( SCIPisCutEfficacious(scip, sol, row) )
   {
#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "Added cut %s: ", cutname);
      SCIPinfoMessage(scip, NULL, "%f <= ", lhs);
      for (j = 0; j < len; j++)
         SCIPinfoMessage(scip, NULL, "+ (%f)*%s", vals[j], SCIPvarGetName(SCIPcolGetVar(cols[j])));
      SCIPinfoMessage(scip, NULL, "\n");
#endif

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
   SCIPfreeBufferArray(scip, &cols);
   SCIPfreeBufferArray(scip, &coeff);

   return SCIP_OKAY;
}

/** approximates the sdpcone using the fact that every diagonal entry must be non-negative, so it adds the LP-cut
 *  \f$ \sum_{j = 1}^m (A_j)_{kk} y_j - (A_0)_{kk} \geq 0 \quad \forall k \leq n \f$
 */
static
SCIP_RETCODE diagGEzero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   int blocksize;
   SCIP_CONSDATA* consdata;
   int nvars;
   int i;
   int j;
   int k;
   int c;
   SCIP_Real rhs;
   char* cutname;
#ifndef NDEBUG
   int snprintfreturn; /* used to check if sdnprintf worked */
#endif
   SCIP_Real* matrix;
   SCIP_Real* cons_array;
   SCIP_Real* lhs_array;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSHDLR* conshdlr;
#ifndef NDEBUG
      const char* conshdlrName;
#endif

      conshdlr = SCIPconsGetHdlr(conss[c]);
      assert( conshdlr != NULL );
#ifndef NDEBUG
      conshdlrName = SCIPconshdlrGetName(conshdlr);
      assert( strcmp(conshdlrName, "SDP") == 0);
#endif

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;
      rhs = SCIPinfinity(scip);

      SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize + 1)) / 2) );
      SCIP_CALL( SCIPconsSdpGetLowerTriangConstMatrix(scip, conss[c], matrix) );

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
            cons_array[k * nvars + j] = 0.0; /*lint !e679*/
         }

         /* go through the nonzeroes of A_j and look for diagonal entries */
         for (i = 0; i < consdata->nvarnonz[j]; i++)
         {
            if ( consdata->col[j][i] == consdata->row[j][i] )
               cons_array[consdata->col[j][i] * nvars + j] = consdata->val[j][i]; /*lint !e679*/
         }
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255));

      /* add the LP-cuts to SCIP */
      for (k = 0; k < blocksize; ++k)
      {
         SCIP_CONS* cons;
         SCIP_CONSHDLRDATA* conshdlrdata;

         conshdlrdata = SCIPconshdlrGetData(conshdlr);
#ifndef NDEBUG
         snprintfreturn = SCIPsnprintf(cutname, 255, "diag_ge_zero_%d", ++(conshdlrdata->ndiaggezerocuts));
         assert( snprintfreturn < 256 ); /* this is the number of positions needed, we gave 255 */
#else
         SCIPsnprintf(cutname, 255, "diag_ge_zero_%d", ++(conshdlrdata->ndiaggezerocuts));
#endif

#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
         SCIPinfoMessage(scip, NULL, "%f <= ", lhs_array[k]);
         for (i = 0; i < consdata->nvars; i++)
         {
            if ( ! SCIPisZero(scip, cons_array[k * consdata->nvars + i]) )
               SCIPinfoMessage(scip, NULL, "+ (%f)*%s", cons_array[k * consdata->nvars + i], SCIPvarGetName(consdata->vars[i]));
         }
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, consdata->vars, cons_array + k * consdata->nvars, lhs_array[k], rhs,
               TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/

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

#if 0
/** presolve-routine that adds some constraints for approximation of the sdpcone, if there is an entry (i,j) in the constant
 *  matrix, than some variable k with \f$ (A_k)_{ii} \f$ needs to be >0 because sdp-matrices are diagonal-dominant (and the same for j)
 *
 *  Not clear if this is really true for all SDPs, probably only works if A_i and A_0 are all semidefinite (or for A_0 at least
 *  have positive diagonal entries) and all variables appearing in the SDP constraint are integer, then sum_{A_i_kk >0}
 *  1*y_i >= 1 is feasible, because it means that (sum A_i)_kk > 0 because all diagonal entries are positive (they can't
 *  cancel each other) and at least one variable needs to be >=1 because this is equal to >0 for integers.
 */
static
SCIP_RETCODE diagDominant(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   SCIP_Bool* nonzerorows;  /* entry i will be 1 if there is an entry (A_0)_ij \f$ for some \f$ j \neq i */
   int blocksize;
   int i;
   int j;
   int nvars;
   int var;
   SCIP_CONS* cons;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* cutname;
#ifndef NDEBUG
   int snprintfreturn;
#endif
   int** diagvars;
   int* ndiagvars;

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nconss >= 0 );
   assert( naddconss != NULL );

   for (i = 0; i < nconss; ++i)
   {
      assert( conss[i] != NULL );
      assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "SDP") == 0 );

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
            assert( ! SCIPisEQ(scip, consdata->constval[j], 0.0) );
            nonzerorows[consdata->constcol[j]] = TRUE;
            nonzerorows[consdata->constrow[j]] = TRUE;
         }
      }

      /* diagvars[i] is an array with all variables with a diagonal entry (i,i) in the corresponding matrix, if nonzerorows[i] is true or NULL otherwise
       * the outer array goes over all rows to ease the access, but only for those that are really needed memory will be allocated */
      SCIP_CALL( SCIPallocBufferArray(scip, &diagvars, blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ndiagvars, blocksize) );
      for (j = 0; j < blocksize; ++j)
      {
         ndiagvars[j] = 0;
         if ( nonzerorows[j] )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &(diagvars[j]), nvars) );
         }
      }

      /* find all variables with corresponding diagonal entries for a row with nonzero non-diagonal constant entry */
      for (var = 0; var < nvars; var++)
      {
         for (j = 0; j < consdata->nvarnonz[var]; j++)
         {
            if ( consdata->col[var][j] == consdata->row[var][j] && nonzerorows[consdata->row[var][j]])
            {
               assert( ! SCIPisEQ(scip, consdata->val[var][j], 0.0) );
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

            conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(conss[i]));
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(cutname, 255, "diag_dom_%d", ++(conshdlrdata->ndiagdomcuts));
            assert( snprintfreturn < 256 );  /* the return is the number of spots needed, we gave 255 */
#else
            SCIPsnprintf(cutname, 255, "diag_dom_%d", ++(conshdlrdata->ndiagdomcuts));
#endif

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
            SCIPinfoMessage(scip, NULL, "1 <= ");
            for (i = 0; i < ndiagvars[j]; i++)
               SCIPinfoMessage(scip, NULL, "+ (%f)*%s", vals[i], SCIPvarGetName(vars[i]));
            SCIPinfoMessage(scip, NULL, "\n");
#endif

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
#endif

/** detects if there are blocks with size one and transfers them to lp-rows */
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddconss,          /**< pointer to store how many constraints were added */
   int*                  ndelconss,          /**< pointer to store how many constraints were deleted */
   SCIP_RESULT*          result              /**< pointer to store if this routine was successfull or if it detected infeasibility */
   )
{
   SCIP_CONSHDLR* hdlr;
   SCIP_CONS* cons;
   SCIP_CONSHDLRDATA* conshdlrdata;
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
   SCIP_CONSDATA* consdata;
#ifndef NDEBUG
   int snprintfreturn; /* used to assert the return code of snprintf */
   const char* hdlrName;
#endif

   *result = SCIP_SUCCESS;

   for (i = 0; i < nconss; ++i)
   {
      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);

#ifndef NDEBUG
      hdlrName = SCIPconshdlrGetName(hdlr);
      assert( strcmp(hdlrName, "SDP") == 0);
#endif

      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* if there is a 1x1 SDP-Block */
      if ( consdata->blocksize == 1 )
      {
         nvars = consdata->nvars;
         nnonz = consdata->nnonz;
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nnonz) );

         /* get all lhs-entries */
         count = 0;

         for (var = 0; var < nvars; var++)
         {
            assert( consdata->nvarnonz[var] <= 1 ); /* in a 1x1 block there may be at most one entry per variable */

            for (j = 0; j < consdata->nvarnonz[var]; j++)
            {
               assert( consdata->col[var][j] == 0 && consdata->row[var][j] == 0 ); /* if the block is size one, all entries should have row and col equal to 0 */
               coeffs[count] = consdata->val[var][j];
               vars[count] = consdata->vars[var];
               count++;
            }
         }

         /* get rhs */
         assert( consdata->constnnonz <= 1 ); /* the 1x1 constant matrix may only have one entry */

         rhs = (consdata->constnnonz == 1) ? consdata->constval[0] : 0.0; /* if this one entry is not 0, than this is the rhs, otherwise it's 0 */

         /* if there is more than one left-hand-side-entry, add a linear constraint, otherwise update the variable bound */
         if ( count > 1 )
         {
            /* add new linear cons */
            conshdlrdata = SCIPconshdlrGetData(hdlr);
            SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255) );
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(cutname, 255, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
            assert( snprintfreturn < 256 ); /* the return is the number of spots needed, we gave 255 */
#else
            SCIPsnprintf(cutname, 255, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
#endif

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added lp-constraint %s: ", cutname);
            for (i = 0; i < consdata->nvars; i++)
               SCIPinfoMessage(scip, NULL, "+ (%f)*%s", coeffs[i], SCIPvarGetName(vars[i]));
            SCIPinfoMessage(scip, NULL, " <= %f \n", rhs);
#endif

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, consdata->nvars, vars, coeffs, rhs, SCIPinfinity(scip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            SCIPfreeBufferArray(scip, &cutname);

            (*naddconss)++;
         }
         else
         {
            /* we compare this new variable with the current (local) bounds, we don't do an epsilon check here, because 10^(-21)*x >= 10^(-19) for
             * epsilon = 10^(-20) would be infeasible then, even though it only says x >= 100 */
            if ( coeffs[0] > 0.0 )
            {
               /* this gives a lower bound, if it is bigger than the current one, we need to update it */
               if ( SCIPisFeasGT(scip, rhs / coeffs[0], SCIPvarGetLbLocal(vars[0])) )
               {
                  SCIPdebugMessage("Changing lower bound of variable %s from %f to %f because of 1x1-SDP-constraint %s!\n",
                        SCIPvarGetName(vars[0]), SCIPvarGetLbLocal(vars[0]), rhs / coeffs[0], SCIPconsGetName(conss[i]));
                  SCIP_CALL( SCIPchgVarLb(scip, vars[0], rhs / coeffs[0]) );
               }
               else
               {
                  SCIPdebugMessage("Deleting 1x1-SDP-constraint %s, new lower bound %f for variable %s no improvement over old bound %f!\n",
                        SCIPconsGetName(conss[i]), rhs / coeffs[0], SCIPvarGetName(vars[0]), SCIPvarGetLbLocal(vars[0]));
               }
            }
            else if ( coeffs[0] < 0.0 )
            {
               /* this gives an upper bound, if it is lower than the current one, we need to update it */
               if (SCIPisFeasLT(scip, -rhs / coeffs[0], SCIPvarGetUbLocal(vars[0])))
               {
                  SCIPdebugMessage("Changing upper bound of variable %s from %f to %f because of 1x1-SDP-constraint %s!\n",
                                          SCIPvarGetName(vars[0]), SCIPvarGetUbLocal(vars[0]), -rhs / coeffs[0], SCIPconsGetName(conss[i]));
                  SCIP_CALL( SCIPchgVarUb(scip, vars[0], -rhs / coeffs[0]) );
               }
               else
               {
                  SCIPdebugMessage("Deleting 1x1-SDP-constraint %s, new upper bound %f for variable %s no improvement over old bound %f!\n",
                        SCIPconsGetName(conss[i]), -rhs / coeffs[0], SCIPvarGetName(vars[0]), SCIPvarGetUbLocal(vars[0]));
               }
            }
            else
            {
               SCIPdebugMessage("Detected 1x1 SDP-block without any nonzero coefficients \n");
               if ( SCIPisFeasGT(scip, rhs, 0.0) )
               {
                  SCIPdebugMessage("Detected infeasibility in 1x1 SDP-block without any nonzero coefficients but with strictly positive rhs\n");
                  *result = SCIP_CUTOFF;
                  /* delete old 1x1 sdpcone */
                  SCIP_CALL( SCIPdelCons(scip, conss[i]) );
                  (*ndelconss)++;

                  SCIPfreeBufferArray(scip, &coeffs);
                  SCIPfreeBufferArray(scip, &vars);

                  return SCIP_OKAY; /* the node is infeasible, we don't care for the other constraints */
               }
            }
         }

         /* delete old 1x1 sdpcone */
         SCIP_CALL(SCIPdelCons(scip, conss[i]));
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &vars);
      }
   }
   return SCIP_OKAY;
}

/** local function to perform (parts of) multiaggregation of a single variable within fixAndAggrVars */
static
SCIP_RETCODE multiaggrVar(
   SCIP*                scip,                /* SCIP pointer */
   SCIP_CONS*           cons,                /* constraint to multiaggregate for */
   int*                 v,                   /* position of the variable that gets (multi-)aggregated */
   SCIP_VAR**           aggrvars,            /* variables this has to be (multi-)aggregated to */
   SCIP_Real*           scalars,             /* scalar parts to multiply with for each variable this is aggregated to */
   int                  naggrvars,           /* number of variables this is (multi-)aggregated to */
   SCIP_Real            constant,            /* the constant part for the (multi-)aggregation */
   int*                 savedcol,            /* array of columns for nonzeros that need to be added to the constant part */
   int*                 savedrow,            /* array of rows for nonzeros that need to be added to the constant part */
   SCIP_Real*           savedval,            /* array of values for nonzeros that need to be added to the constant part */
   int*                 nfixednonz,          /* length of the arrays of saved nonzeros for the constant part */
   int*                 vararraylength       /* length of the variable array */
   )
{
   int i;
   SCIP_CONSDATA* consdata;
   int startind;
   int aggrind;
   int aggrtargetlength;
   int globalnvars;
   int aggrconsind;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( scalars != NULL );
   assert( naggrvars > 0 );
   assert( savedcol != NULL );
   assert( savedrow != NULL );
   assert( savedval != NULL );
   assert( nfixednonz != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* save the current nfixednonz-index, all entries starting from here will need to be added to the variables this is aggregated to */
   startind = *nfixednonz;

   if ( SCIPisEQ(scip, constant, 0.0) )
   {
      /* if there is no constant part, we just save the nonzeros to add them to the variables they are aggregated to, we do this to be able to remove
       * the nonzero-arrays for this variable to be able to fill it with a newly inserted variable, as copying all variables, if we created an empty
       * gap in the variable arrays, will be more time consuming then copying all variables (as we will usually have more variables than nonzeros per
       * variable */
      for (i = 0; i < consdata->nvarnonz[*v]; i++)
      {
         savedcol[*nfixednonz] = consdata->col[*v][i];
         savedrow[*nfixednonz] = consdata->row[*v][i];
         savedval[*nfixednonz] = consdata->val[*v][i];
         (*nfixednonz)++;
      }
   }
   else
   {
      for (i = 0; i < consdata->nvarnonz[*v]; i++)
      {
         savedcol[*nfixednonz] = consdata->col[*v][i];
         savedrow[*nfixednonz] = consdata->row[*v][i];
         savedval[*nfixednonz] = consdata->val[*v][i] * constant;
         (*nfixednonz)++;
      }
   }
   assert( *nfixednonz - startind == consdata->nvarnonz[*v] );

   /* sort them by nondecreasing row and then col to make the search for already existing entries easier (this is done here, because it
    * only needs to be done once and not for each variable this is multiaggregated to), we add startind to the pointers to only start where we started
    * inserting, the number of elements added to the saved arrays for this variable is nfixednonz - startind */
   SdpVarfixerSortRowCol(savedrow + startind, savedcol + startind, savedval + startind, *nfixednonz - startind);

   /* fill the empty spot of the (multi-)aggregated variable with the last variable of this constraint (as they don't have to be sorted) */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[*v]) );
   consdata->col[*v] = consdata->col[consdata->nvars - 1];
   consdata->row[*v] = consdata->row[consdata->nvars - 1];
   consdata->val[*v] = consdata->val[consdata->nvars - 1];
   consdata->nvarnonz[*v] = consdata->nvarnonz[consdata->nvars - 1];
   consdata->vars[*v] = consdata->vars[consdata->nvars - 1];
   (consdata->nvars)--;
   (*v)--; /* we need to check again if the variable we just shifted to this position also needs to be (multi-)aggregated */

   /* iterate over all variables this was aggregated to and insert the corresponding nonzeroes */
   for (aggrind = 0; aggrind < naggrvars; aggrind++)
   {
      /* check if the variable already exists in this block */
      aggrconsind = -1;
      for (i = 0; i < consdata->nvars; i++)
      {
         if ( consdata->vars[i] == aggrvars[aggrind] )
         {
            aggrconsind = i;
            break;
         }
      }

      if ( aggrconsind > -1 )
      {
         /* if the varialbe to aggregate to is already part of this sdp-constraint, just add the nonzeros of the old variable to it */

         /* resize the arrays to the maximally needed length */
         aggrtargetlength = consdata->nvarnonz[aggrconsind] + *nfixednonz - startind;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );

         if ( SCIPisEQ(scip, constant, 0.0) )
         {
            /* in this case we saved the original values in savedval, we add startind to the pointers to only add those from
             * the current variable, the number of entries is the current position minus the position whre we started */
            SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow + startind, savedcol + startind, savedval + startind,
                        *nfixednonz - startind, TRUE, scalars[aggrind], consdata->row[aggrconsind], consdata->col[aggrconsind],
                        consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]), aggrtargetlength) );
         }
         else
         {
            /* in this case we saved the original values * constant, so we now have to divide by constant, we add startind to the pointers
             * to only add those from the current variable, the number of entries is the current position minus the position whre we started */
            SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow + startind, savedcol + startind, savedval + startind,
                        *nfixednonz - startind, TRUE, scalars[aggrind] / constant, consdata->row[aggrconsind], consdata->col[aggrconsind],
                        consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]), aggrtargetlength) );
         }

         /* shrink them again if nonzeros could be combined */
         assert( consdata->nvarnonz[aggrconsind] <= aggrtargetlength );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
      }
      else
      {
         /* the variable has to be added to this constraint */

         SCIPdebugMessage("adding variable %s to SDP constraint %s because of (multi-)aggregation\n", SCIPvarGetName(aggrvars[aggrind]), SCIPconsGetName(cons));

         /* check if we have to enlarge the arrays */
         if ( consdata->nvars == *vararraylength )
         {
            globalnvars = SCIPgetNVars(scip);

            /* we don't want to enlarge this by one for every variable added, so we immediately set it to the maximum possible size */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, *vararraylength, globalnvars) );
            *vararraylength = globalnvars;
         }

         /* we insert this variable at the last position, as the ordering doesn't matter */
         SCIP_CALL( SCIPcaptureVar(scip, aggrvars[aggrind]) );
         consdata->vars[consdata->nvars] = aggrvars[aggrind];

         /* as there were no nonzeros thus far, we can just duplicate the saved arrays to get the nonzeros for the new variable */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), savedcol + startind, *nfixednonz - startind) ); /*lint !e776*/
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), savedrow + startind, *nfixednonz - startind) ); /*lint !e776*/

         /* if scalars[aggrind] = constant, we would multiply with 1, if constant = 0, we didn't divide by constant, so in these cases, we can just
          * memcopy the array of nonzero-values */
         /* TODO: only checking scalar and constant for feas-equality might lead to big differences, if the nonzeros they are multiplied with are big,
          * e.g. scalar = 0, constant = 10^(-6), nonzero = 10^(10) leads to new nonzero of 10^4 instead of 0 */
         if ( SCIPisEQ(scip, scalars[aggrind], constant) || SCIPisEQ(scip, constant, 0.0) )
         {
            SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), savedval + startind, *nfixednonz - startind) ); /*lint !e776*/
            consdata->nvarnonz[consdata->nvars] = *nfixednonz - startind; /* as there were no nonzeros thus far, the number of nonzeros equals
                                                                           * the number of nonzeros of the aggregated variable */
         }
         else  /* we have to multiply all entries by scalar before inserting them */
         {
            SCIP_Real epsilon;

            SCIP_CALL( SCIPgetRealParam(scip, "numerics/epsilon", &epsilon) );

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), *nfixednonz - startind) );

            consdata->nvarnonz[consdata->nvars] = 0;

            for (i = 0; i < *nfixednonz - startind; i++)
            {
               /* if both scalar and savedval are small this might become too small */
               if ( (scalars[i] / constant) * savedval[startind + i] >= epsilon )  /*lint !e679*/
               {
                  consdata->val[consdata->nvars][consdata->nvarnonz[consdata->nvars]] = (scalars[i] / constant) * savedval[startind + i]; /*lint !e679*/
                  consdata->nvarnonz[consdata->nvars]++;
               }
            }
         }

         if ( consdata->nvarnonz[consdata->nvars] > 0 ) /* if scalar and all savedvals were to small */
            consdata->nvars++;
      }
   }

   /* if the constant part is zero, we may delete the nonzero-entries from the saved arrays (by resetting the nfixednonz entry to where
    * it started, so that these entries will be overwritten */
   if ( SCIPisEQ(scip, constant, 0.0) )
      *nfixednonz = startind;

   return SCIP_OKAY;
}


/** presolve routine that looks through the data and handles fixed, (multi-)aggregated and negated variables */
static
SCIP_RETCODE fixAndAggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array with constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_Bool             aggregate           /**< do we want to (mutli-)aggregate variables ? */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int nfixednonz;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   int c;
   int v;
   int arraylength;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_Real scalar;
   SCIP_Real* scalars;
   int naggrvars;
   SCIP_Real constant;
   int requiredsize;
   int globalnvars;
   int vararraylength;
   SCIP_Bool negated;


   /* loop over all variables once, add all fixed to savedrow/col/val, for all multiaggregated variables, if constant-scalar =!= 0, add
    * constant-scalar * entry to savedrow/col/val and call mergeArrays for all aggrvars for savedrow[startindex of this var] and scalar/constant-scalar,
    * if constant-scalar = 0, add 1*entry to savedrow/col/val, call mergeArrays for all aggrvars for savedrow[startindex of this var] and scalar and later
    * reduze the saved size of savedrow/col/val by the number of nonzeros of the mutliagrregated variable to not add them to the constant part later */

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nconss >= 0 );

   SCIPdebugMessage("Calling fixAndAggrVars with aggregate = %u\n", aggregate);

   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );
      assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0);

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* allocate memory to save nonzeros that need to be fixed */
      SCIP_CALL( SCIPallocBufferArray(scip, &savedcol, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedrow, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedval, consdata->nnonz) );

      /* initialize this with zero for each block */
      nfixednonz = 0;

      vararraylength = consdata->nvars;
      globalnvars = SCIPgetNVars(scip);

      for (v = 0; v < consdata->nvars; v++)/*lint --e{850}*/
      {
         negated = FALSE;
         /* if the variable is negated, get the negation var */
         if ( SCIPvarIsBinary(consdata->vars[v]) && SCIPvarIsNegated(consdata->vars[v]) )
         {
            negated = TRUE;
            var = SCIPvarGetNegationVar(consdata->vars[v]);
         }
         else
            var = consdata->vars[v];

         /* check if the variable is fixed in SCIP */
         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)))
         {
            assert( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) );

            SCIPdebugMessage("Globally fixing Variable %s to value %f !\n", SCIPvarGetName(var), SCIPvarGetLbGlobal(var));

            if ( ((! negated) && (! SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0))) || (negated && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0)) )
            {
               /* the nonzeros are saved to later be inserted into the constant part (this is only done after all nonzeros of fixed variables have been
                * assembled, because we need to sort the constant nonzeros and loop over them, which we only want to do once and not once for each fixed
                * variable) */
               for (i = 0; i < consdata->nvarnonz[v]; i++)
               {
                  savedcol[nfixednonz] = consdata->col[v][i];
                  savedrow[nfixednonz] = consdata->row[v][i];
                  /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */
                  if ( ! negated )
                     savedval[nfixednonz] = consdata->val[v][i] * SCIPvarGetLbGlobal(var);
                  else
                     savedval[nfixednonz] = consdata->val[v][i]; /* if it is the negation of a variable fixed to zero, this variable is fixed to one */

                  nfixednonz++;
                  consdata->nnonz--;
               }
            }
            else
            {
               /* if the variable is fixed to zero, the nonzeros will just vanish, so we only reduce the number of nonzeros */
               consdata->nnonz -= consdata->nvarnonz[v];
            }
            /* as the variables don't need to be sorted, we just put the last variable into the empty spot and decrease sizes by one (at the end) */
            SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[v])) );
            consdata->col[v] = consdata->col[consdata->nvars - 1];
            consdata->row[v] = consdata->row[consdata->nvars - 1];
            consdata->val[v] = consdata->val[consdata->nvars - 1];
            consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
            consdata->vars[v] = consdata->vars[consdata->nvars - 1];
            consdata->nvars--;
            v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
         }
         else if ( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED ||
                   SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR)
                  && aggregate )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, globalnvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &scalars, globalnvars) );

            /* this is how they should be initialized before calling SCIPgetProbvarLinearSum */
            if (! negated)
            {
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 0.0;
               scalars[0] = 1.0;
            }
            else
            {
               /* if this variable is the negation of var, than we look for a representation of 1-var */
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 1.0;
               scalars[0] = -1.0;
            }

            /* get the variables this var was aggregated to */
            SCIP_CALL( SCIPgetProbvarLinearSum(scip, aggrvars, scalars, &naggrvars, globalnvars, &constant, &requiredsize, TRUE) );
            assert( requiredsize <= globalnvars ); /* requiredsize is the number of empty spots in aggrvars needed, globalnvars is the number
                                                    * of spots we provided */

            /* Debugmessages for the (multi-)aggregation */
#ifdef SCIP_DEBUG
            if ( SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_AGGREGATED )
               SCIPdebugMessage("aggregating variable %s to ", SCIPvarGetName(var));
            else
               SCIPdebugMessage("multiaggregating variable %s to ", SCIPvarGetName(var));
            for (i = 0; i < naggrvars; i++)
               printf("+ (%f2) * %s ", scalars[i], SCIPvarGetName(aggrvars[i]));
            printf("+ (%f2) \n", constant);
#endif

            /* add the nonzeros to the saved-arrays for the constant part, remove the nonzeros for the old variables and add them to the variables this variable
             * was (multi-)aggregated to */
            SCIP_CALL( multiaggrVar(scip, conss[c], &v, aggrvars, scalars, naggrvars, constant, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );

            SCIPfreeBufferArray(scip, &aggrvars);
            SCIPfreeBufferArray(scip, &scalars);
         }
         else if ( negated && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN)
                  && aggregate)
         {
             /* if var1 is the negation of var2, then this is equivalent to it being aggregated to -var2 + 1 = 1 - var2 */

            SCIPdebugMessage("Changing variable %s to negation of variable %s !\n", SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(var));

            scalar = -1.0;

            SCIP_CALL( multiaggrVar(scip, conss[c], &v, &var, &scalar, 1, 1.0, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );
         }
      }

      /* shrink the variable arrays if they were enlarged too much (or more vars were removed than added) */
      assert( consdata->nvars <= vararraylength );
      if ( consdata->nvars < vararraylength )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, vararraylength, consdata->nvars) );
      }

      /* allocate the maximally needed memory for inserting the fixed variables into the constant part */
      arraylength = consdata->constnnonz + nfixednonz;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), consdata->constnnonz, arraylength) );

      /* insert the fixed variables into the constant arrays, as we have +A_i but -A_0 we mutliply them by -1 */
      SCIP_CALL( SdpVarfixerMergeArrays(SCIPblkmem(scip), savedrow, savedcol, savedval, nfixednonz, FALSE, -1.0, consdata->constrow, consdata->constcol,
                consdata->constval, &(consdata->constnnonz), arraylength) );

      assert( consdata->constnnonz <= arraylength ); /* the allocated memory should always be sufficient */

      /* shrink the arrays if nonzeros could be combined */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), arraylength, consdata->constnnonz) );

      /* free the saved arrays */
      SCIPfreeBufferArray(scip, &savedval);
      SCIPfreeBufferArray(scip, &savedrow);
      SCIPfreeBufferArray(scip, &savedcol);

      /* recompute sdpnnonz */
      consdata->nnonz = 0;
      for (v = 0; v < consdata->nvars; v++)
         consdata->nnonz += consdata->nvarnonz[v];
   }

   return SCIP_OKAY;
}

/** informs constraint handler that the presolving process is being started */
static
SCIP_DECL_CONSINITPRE(consInitpreSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   conshdlrdata->neigveccuts = 0; /* this is used to give the eigenvector-cuts distinguishable names */
   conshdlrdata->ndiaggezerocuts = 0; /* this is used to give the diagGEzero-cuts distinguishable names */
   conshdlrdata->ndiagdomcuts = 0; /* this is used to give the diagDominant-cuts distinguishable names */
   conshdlrdata->n1x1blocks = 0; /* this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */

   return SCIP_OKAY;
}

/** locks a variable up if the corresponding constraint matrix is not positive semidefinite, locks it down if it is not negative semidefinite */
static
SCIP_DECL_CONSLOCK(consLockSdp)
{/*lint --e{715}*/
   SCIP_Real* Aj;
   SCIP_CONSDATA* consdata;
   int blocksize;
   int var;
   int nvars;
   SCIP_Real eigenvalue;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;
   nvars = consdata->nvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) );

   for (var = 0; var < nvars; var++)
   {
      SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, var, Aj) );

      /* compute the smallest eigenvalue */
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPblkmem(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
      {
         /* as the lowest eigenvalue is negative, the matrix is not positive semidefinite, so adding more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlocksneg, nlockspos) );
      }

      /* if the smallest eigenvalue is already positive, we don't need to compute the biggest one */
      if ( SCIPisPositive(scip, eigenvalue) )
      {
         /* as an eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
          * semidefiniteness of the SDP-matrix */
         SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlockspos, nlocksneg) );
      }
      else
      {
         /* compute the biggest eigenvalue */
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPblkmem(scip), FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
         if ( SCIPisPositive(scip, eigenvalue) )
         {
            /* as the biggest eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
             * semidefiniteness of the SDP-matrix */
            SCIP_CALL( SCIPaddVarLocks(scip, consdata->vars[var], nlockspos, nlocksneg) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &Aj);

   return SCIP_OKAY;
}

/** after presolving variables are fixed and multiaggregated */
static
SCIP_DECL_CONSEXITPRE(consExitpreSdp)
{/*lint --e{715}*/
   assert( scip != NULL );
   assert( conss != NULL );

   SCIP_CALL( fixAndAggrVars(scip, conss, nconss, TRUE) );

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSdp)
{/*lint --e{715}*/
   assert( result != 0 );

   if ( nrounds == 0 )
   {
      SCIP_CALL( diagGEzero(scip, conss, nconss, naddconss) );
   }

   SCIP_CALL( move_1x1_blocks_to_lp(scip, conss, nconss, naddconss, ndelconss, result) );

   SCIP_CALL( fixAndAggrVars(scip, conss, nconss, FALSE) ); /* the FALSE means we only do fixings and not (multi-)aggregations or negations */

#if 0
   if ( nrounds == 0 )
   {
      SCIP_CALL( diagDominant(scip, conss, nconss, naddconss) ); /* TODO: could be activated for some problem classes but doesn't work in the general case */
   }
#endif

   return SCIP_OKAY;
}

/** creates transformed constraint */
static
SCIP_DECL_CONSTRANS(consTransSdp)
{/*lint --e{715}*/
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
   {
      targetdata->vars[i] = SCIPvarGetTransVar(sourcedata->vars[i]);
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->vars[i]) );
   }

   /* copy the constant nonzeros */
   targetdata->constnnonz = sourcedata->constnnonz;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constcol), sourcedata->constcol, sourcedata->constnnonz));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constrow), sourcedata->constrow, sourcedata->constnnonz));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constval), sourcedata->constval, sourcedata->constnnonz));

   /* copy the maxrhsentry */
   targetdata->maxrhsentry = sourcedata->maxrhsentry;

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** checks feasiblity of constraint, e.g. the positive semidefiniteness */
static
SCIP_DECL_CONSCHECK(consCheckSdp)
{/*lint --e{715}*/
   int i;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], sol, checkintegrality, checklprows, printreason, result) );
      if ( *result == SCIP_INFEASIBLE )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** enforce pseudo solution method
 *
 *  Returns didnotrun, if objinfeasible, computes cut otherwise.
 */
static
SCIP_DECL_CONSENFOPS(consEnfopsSdp)
{/*lint --e{715}*/
   int i;

   assert( scip != NULL );
   assert( result != NULL );
   assert( conss != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( objinfeasible )
   {
      SCIPdebugMessage("-> pseudo solution is objective infeasible, return.\n");
      return SCIP_OKAY;
   }

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], NULL, 0, 0, 0, result) );

      if (*result == SCIP_INFEASIBLE)
      {
         /* if it is infeasible for one SDP constraint, it is infeasible for the whole problem */
         SCIPdebugMessage("-> pseudo solution infeasible for SDP-constraint %s, return.\n", SCIPconsGetName(conss[i]));
         return SCIP_OKAY;
      }
   }

   SCIPdebugMessage("-> pseudo solution feasible for all SDP-constraints.\n");

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions
 *
 *  Enforce lp solution method, if some block is not psd an eigenvector cut is added.
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool all_feasible = TRUE;
   SCIP_Bool separated = FALSE;
   char* cutname;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif
   int i;
   int j;
   int nvars;
   SCIP_ROW* row;
   SCIP_Bool infeasible;
   SCIP_Real lhs;
   SCIP_Real* coeff;
   SCIP_Real rhs;
#if 0 /* TODO: see below */
   SCIP_VAR** vars;
   int count;
#endif

   assert( result != NULL );
   *result = SCIP_FEASIBLE;

   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conss[i], NULL, 0, 0, 0, result) );
      if ( *result == SCIP_FEASIBLE )
         continue;

      all_feasible = FALSE;

      nvars = consdata->nvars;
      lhs = 0.0;
      coeff = NULL;

      SCIP_CALL( SCIPallocBufferArray(scip, &coeff, nvars) );
      SCIP_CALL( cutUsingEigenvector(scip, conss[i], NULL, coeff, &lhs) );

      rhs = SCIPinfinity(scip);
      conshdlrdata = SCIPconshdlrGetData(conshdlr);

      SCIP_CALL( SCIPallocBufferArray(scip, &cutname, 255) );
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(cutname, 255, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
      assert( snprintfreturn < 256 ); /* this is the number of spots needed, we gave 255 */
#else
      SCIPsnprintf(cutname, 255, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
#endif

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conshdlr, cutname , lhs, rhs, FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for (j = 0; j < nvars; ++j)
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, consdata->vars[j], coeff[j]) );
      }

      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "Added cut %s: ", cutname);
      SCIPinfoMessage(scip, NULL, "%f <= ", lhs);
      for (j = 0; j < nvars; j++)
         SCIPinfoMessage(scip, NULL, "+ (%f)*%s", coeff[j], SCIPvarGetName(consdata->vars[j]));
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      SCIP_CALL( SCIPaddCut(scip, NULL, row, FALSE, &infeasible) );

      if ( infeasible )
         *result = SCIP_CUTOFF;
      else
      {
         SCIP_CALL( SCIPaddPoolCut(scip, row) );

         SCIP_CALL( SCIPresetConsAge(scip, conss[i]) );
         *result = SCIP_SEPARATED;
         separated = TRUE;
      }
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
      SCIPfreeBufferArray(scip, &cutname);
      SCIPfreeBufferArray(scip, &coeff);
   }
   if ( all_feasible )
      return SCIP_OKAY;

   if ( separated )
      *result = SCIP_SEPARATED;
#if 0 /* TODO: should this be done here or is it the task of the conshdlr_integer ? if we do it, we should use Feastol and also check for integrality of solution and the counter should obviously be removed */
   vars = SCIPgetVars(scip);
   count = 0;
   for (i = 0; i < SCIPgetNVars(scip); ++i)
   {
      if ( !SCIPisRelEQ(scip, SCIPvarGetLbLocal(vars[i]), SCIPvarGetUbLocal(vars[i])) && SCIPvarIsIntegral(vars[i]))
      {
         SCIP_CALL( SCIPaddExternBranchCand(scip, vars[i], 10000.0, SCIP_INVALID) );
         count++;
      }
   }
#endif

   return SCIP_OKAY;
}

/** separates a solution using constraint specific ideas, gives cuts to SCIP */
static
SCIP_DECL_CONSSEPASOL(consSepasolSdp)
{/*lint --e{715}*/
   int i;

   *result = SCIP_DIDNOTFIND;
   for (i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL( separateSol(scip, conshdlr, conss[i], sol, result) );
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSdp)
{/*lint --e{715}*/
   int i;

   *result = SCIP_DIDNOTFIND;
   for (i = 0; i < nusefulconss; ++i)
   {
      SCIP_CALL( separateSol(scip, conshdlr, conss[i], NULL, result) );
   }

   return SCIP_OKAY;
}

/** delete method of SDP constrainthandler */
static
SCIP_DECL_CONSDELETE(consDeleteSdp)
{/*lint --e{715}*/
   int i;

   assert( cons != NULL );
   assert( consdata != NULL );

   SCIPdebugMessage("deleting %s \n", SCIPconsGetName(cons));

   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col[i], (*consdata)->nvarnonz[i]);
   }

   /* release all variables */
   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[i])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constval, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constrow, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constcol, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->nvarnonz, (*consdata)->nvars);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** free method of sdp constrainthandler */
static
SCIP_DECL_CONSFREE(consFreeSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** copy an SDP constraint handler */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySdp)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConshdlrSdp(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** copy an SDP constraint*/
static
SCIP_DECL_CONSCOPY(consCopySdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_Bool success;
   SCIP_VAR** targetvars;
   SCIP_VAR* var;
   int i;

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );

   SCIPdebugMessage("Copying SDP constraint %s\n", SCIPconsGetName(sourcecons));

   *valid = TRUE;

   /* as we can only map active variables, we have to make sure, that the constraint contains no fixed or (multi-)aggregated vars, after
    * exitpresolve (stage 6) this should always be the case, earlier than that we need to call fixAndAggrVars */
   if ( SCIPgetStage(sourcescip)  <= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CALL( fixAndAggrVars(scip, &sourcecons, 1, TRUE) );
   }


   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, sourcedata->nvars) );

   /* map all variables in the constraint */
   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->vars[i], &var, varmap, consmap, global, &success) );
      if ( success )
      {
         targetvars[i] = var;
         SCIP_CALL( SCIPcaptureVar(scip, targetvars[i]) );
      }
      else
         *valid = FALSE;
   }

   /* create the new constraint */
   SCIP_CALL( SCIPcreateConsSdp( scip, cons, name, sourcedata->nvars, sourcedata->nnonz, sourcedata->blocksize, sourcedata->nvarnonz,
                                 sourcedata->col, sourcedata->row, sourcedata->val, targetvars, sourcedata->constnnonz,
                                 sourcedata->constcol, sourcedata->constrow, sourcedata->constval) );

   SCIPfreeBufferArray(scip, &targetvars);

   return SCIP_OKAY;
}

/** print an SDP constraint */
static
SCIP_DECL_CONSPRINT(consPrintSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* fullmatrix;
   int v;
   int i;
   int j;

   assert( cons != NULL );

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
            fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
      }

      /* then add the nonzeros */
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         fullmatrix[consdata->row[v][i] * consdata->blocksize + consdata->col[v][i]] = consdata->val[v][i]; /* lower triangular entry */ /*lint !e679*/
         fullmatrix[consdata->col[v][i] * consdata->blocksize + consdata->row[v][i]] = consdata->val[v][i]; /* upper triangular entry */ /*lint !e679*/
      }

      /* print it */
      SCIPinfoMessage(scip, file, "+\n");
      for (i = 0; i < consdata->blocksize; i++)
      {
         SCIPinfoMessage(scip, file, "( ");
         for (j = 0; j < consdata->blocksize; j++)
            SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
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
         fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
   }

   /* then add the nonzeros */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      fullmatrix[consdata->constrow[i] * consdata->blocksize + consdata->constcol[i]] = consdata->constval[i]; /* lower triangular entry */ /*lint !e679*/
      fullmatrix[consdata->constcol[i] * consdata->blocksize + consdata->constrow[i]] = consdata->constval[i]; /* upper triangular entry */ /*lint !e679*/
   }

   /* print it */
   SCIPinfoMessage(scip, file, "-\n");
   for (i = 0; i < consdata->blocksize; i++)
   {
      SCIPinfoMessage(scip, file, "( ");
      for (j = 0; j < consdata->blocksize; j++)
         SCIPinfoMessage(scip, file, "%f ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
      SCIPinfoMessage(scip, file, ")\n");
   }
   SCIPinfoMessage(scip, file, ">=0\n");

   SCIPfreeBufferArray(scip, &fullmatrix);

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables (if possible) */
static
SCIP_DECL_CONSGETVARS(consGetVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int i;
   int nvars;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars != NULL );
   assert( success != NULL );
   assert( varssize >= 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   if ( nvars > varssize )
   {
      SCIPdebugMessage("consGetVarsIndicator called for array of size %d, needed size %d.\n", varssize, nvars);
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (i = 0; i < nvars; i++)
   {
      vars[i] = consdata->vars[i];
   }

   *success = TRUE;
   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables (if possible) */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( success != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nvars = consdata->nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

/** creates the handler for sdp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != 0 );

   /* allocate memory for the conshdlrdata */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSdp, consEnfopsSdp, consCheckSdp, consLockSdp, conshdlrdata) );

   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSdp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSdp) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySdp, consCopySdp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr,consInitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_DELAYPRESOL) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSdp) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSdp) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSdp) );

   return SCIP_OKAY;
}

/** get the data belonging to a single SDP-constraint
 *
 *  In arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length.
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown.
 */
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
   SCIP_VAR**            vars,               /**< the SCIP variables present in this constraint, indexing equals indices in col/row/val */
   int*                  constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint, also length of the const arrays */
   int*                  constcol,           /**< pointer to column indices of the constant nonzeroes */
   int*                  constrow,           /**< pointer to row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< pointer to values of the constant nonzeroes */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   const char* name;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( nnonz != NULL );
   assert( blocksize != NULL );
   assert( arraylength != NULL );
   assert( nvarnonz != NULL );
   assert( col != NULL );
   assert( row != NULL );
   assert( val != NULL );
   assert( vars != NULL );
   assert( constnnonz != NULL );

   consdata = SCIPconsGetData(cons);
   name = SCIPconsGetName(cons);

   assert( consdata->constnnonz == 0 || ( constcol != NULL && constrow != NULL && constval != NULL ) );

   *nvars = consdata->nvars;
   *nnonz = consdata->nnonz;
   *blocksize = consdata->blocksize;

   for (i = 0; i < consdata->nvars; i++)
      vars[i] = consdata->vars[i];

   /* check that the sdp-arrays are long enough to store the information */
   if ( *arraylength < consdata->nvars )
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
         /* set the pointers for each variable */
         col[i] = consdata->col[i];
         row[i] = consdata->row[i];
         val[i] = consdata->val[i];
      }
   }

   /* set the constant pointers (if a constant part exists) */
   if ( consdata->constnnonz > 0 )
   {
      if ( consdata->constnnonz > *constnnonz )
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
   }

   *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the number of nonzeroes and constant nonzeroes for this SDP constraint
 *
 *  Either nnonz or constnnonz may be NULL.
 */
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  constnnonz          /**< number of nonzeroes in the constant part of this SDP constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( nnonz != NULL )
      *nnonz = consdata->nnonz;

   if ( constnnonz != NULL )
      *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the full constraint Matrix \f$ A_j \f$ for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int                   j,                  /**< the variable j to get the corresponding matrix \f$ A_j \f$ for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( j >= 0 );
   assert( Aj != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;

   assert( j < consdata->nvars );

   for (i = 0; i < blocksize * blocksize; i++)
      Aj[i] = 0;

   for (i = 0; i < consdata->nvarnonz[j]; i++)
   {
      Aj[consdata->col[j][i] * blocksize + consdata->row[j][i]] = consdata->val[j][i]; /*lint !e679*/
      Aj[consdata->row[j][i] * blocksize + consdata->col[j][i]] = consdata->val[j][i]; /*lint !e679*/
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
   int blocksize;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   for (i = 0; i < blocksize; i++)
   {
      for (j = 0; j < blocksize; j++)
         mat[i * blocksize + j] = 0.0; /*lint !e679*/
   }

   for (i = 0; i < consdata->constnnonz; i++)
   {
      mat[consdata->constcol[i] * blocksize + consdata->constrow[i]] = consdata->constval[i]; /*lint !e679*/
      mat[consdata->constrow[i] * blocksize + consdata->constcol[i]] = consdata->constval[i]; /*lint !e679*/
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
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1)) / 2; i++)
      mat[i] = 0.0;

   for (i = 0; i < consdata->constnnonz; i++)
      mat[compLowerTriangPos(consdata->constrow[i], consdata->constcol[i])] = consdata->constval[i];

   return SCIP_OKAY;
}

/** creates an SDP-constraint */
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
   SCIP_VAR**            vars,               /**< SCIP_VARiables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeroes */
   int*                  constrow,           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< values of the constant nonzeroes */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( name != NULL );
   assert( nvars >= 0 );
   assert( nnonz >= 0 );
   assert( blocksize >= 0 );
   assert( constnnonz >= 0 );
   assert( nvars == 0 || vars != NULL );
   assert( nnonz == 0 || (nvarnonz != NULL && col != NULL && row != NULL && val != NULL ));
   assert( constnnonz == 0 || (constcol != NULL && constrow != NULL && constval != NULL ));

   conshdlr = SCIPfindConshdlr(scip, "SDP");
   if ( conshdlr == NULL )
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
      assert( nvarnonz[i] >= 0 );

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
         assert( col[i][j] >= 0 );
         assert( col[i][j] < blocksize );
         assert( row[i][j] >= 0 );
         assert( row[i][j] < blocksize );

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
   {
      consdata->vars[i] = vars[i];
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

   return SCIP_OKAY;
}
