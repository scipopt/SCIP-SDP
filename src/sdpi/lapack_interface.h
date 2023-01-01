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

/**@file   lapack_interface.h
 * @brief  interface methods for eigenvector computation and matrix multiplication using openblas
 * @author Tristan Gally
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LAPACK_INTERFACE_H__
#define __SCIP_LAPACK_INTERFACE_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** computes the i-th eigenvalue of a symmetric matrix using LAPACK, where 1 is the smallest and n the largest, matrix has to be given with all \f$n^2\f$ entries */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackComputeIthEigenvalue(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed - will be destroyed! */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to store eigenvector */
   );

/** computes i-th eigenvalue of a symmetric matrix using alternative algorithm in LAPACK, where 1 is the smallest and n the largest, matrix has to be given with all \f$n^2\f$ entries */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackComputeIthEigenvalueAlternative(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Bool             geteigenvectors,    /**< Should also the eigenvectors be computed? */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed - will be destroyed! */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to array to store eigenvector */
   );

/** computes eigenvectors corresponding to negative eigenvalues of a symmetric matrix using LAPACK, matrix has to be given with all \f$n^2\f$ entries */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackComputeEigenvectorsNegative(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvectors should be computed - will be destroyed! */
   SCIP_Real             tol,                /**< tolerance; the eigenvalues will be in the interval (-1e20, -tol] */
   int*                  neigenvalues,       /**< pointer to store the number of negative eigenvalues */
   SCIP_Real*            eigenvalues,        /**< array for eigenvalues (should be length n) */
   SCIP_Real*            eigenvectors        /**< array for eigenvectors (should be length n*n), eigenvectors are given as rows  */
   );

/** computes the eigenvector decomposition of a symmetric matrix using LAPACK */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackComputeEigenvectorDecomposition(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which the decomposition should be computed - will be destroyed! */
   SCIP_Real*            eigenvalues,        /**< pointer to store eigenvalues (should be length n) */
   SCIP_Real*            eigenvectors        /**< pointer to store eigenvectors (should be length n*n), eigenvectors are given as rows  */
   );

/** performs matrix-vector-multiplication using BLAS */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackMatrixVectorMult(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   SCIP_Real*            matrix,             /**< the matrix we want to multiply */
   SCIP_Real*            vector,             /**< vector we want to multiply with the matrix */
   SCIP_Real*            result              /**< pointer to store the resulting vector */
   );

/** performs matrix-matrix-multiplication A*B using BLAS */
SCIP_EXPORT
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
   );

/** computes the minimum-norm solution to a real linear least squares problem: minimize 2-norm(| b - A*x |) using LAPACK
 (uses singular value decomposition of A). A is an M-by-N matrix which may be rank-deficient. */
SCIP_RETCODE SCIPlapackLinearSolve(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   m,                  /**< number of rows of A */
   int                   n,                  /**< number of columns of A */
   SCIP_Real*            A,                  /**< coefficient matrix of the linear system */
   SCIP_Real*            b,                  /**< right-hand side of the linear system (should be length n) */
   SCIP_Real*            x                   /**< pointer to store values for x (should be length n) */
   );

#ifdef __cplusplus
}
#endif

#endif
