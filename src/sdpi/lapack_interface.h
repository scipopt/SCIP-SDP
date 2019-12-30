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
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed */
   int                   i,                  /**< index of eigenvalue to be computed */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< pointer to store eigenvector */
   );

/** computes the eigenvector decomposition of a symmetric matrix using LAPACK */
SCIP_EXPORT
SCIP_RETCODE SCIPlapackComputeEigenvectorDecomposition(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which the decomposition should be computed */
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

#ifdef __cplusplus
}
#endif

#endif
