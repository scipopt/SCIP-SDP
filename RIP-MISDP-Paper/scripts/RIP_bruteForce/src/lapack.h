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

/**@file   lapack.h
 * @brief  interface methods for eigenvector computation and matrix multiplication using different versions of LAPACK and BLAS
 * @author Tristan Gally
 *
 * This file is used to call the LAPACK routine DSYEVR (double-symmetric-eigenvector computation) and the
 * BLAS routine DGEMV (double-general-matrix-vector multiplication). It is needed because different SDP-solvers
 * need different BLAS/LAPACK-versions with different data types (for example long long int for
 * Openblas/SDPA vs. int for ATLAS/DSDP).
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_LAPACK_H__
#define __SCIP_LAPACK_H__

#ifdef __cplusplus
extern "C" {
#endif

/** computes smallest and largest eigenvalue of a given matrix using LAPACK */
void SCIPlapackComputeSmallestLargestEigenvalue(
   int                   n,                  /**< size of matrix */
   double*               A,                  /**< matrix for which eigenvalues should be computed */
   double*               smallestev,         /**< pointer to store smallest eigenvalue */
   double*               largestev           /**< pointer to store largest eigenvalue */
   );

/** for given A computes A^T A using BLAS */
void SCIPlapackMatrixTransposedMatrixMult(
   int                   nrows,              /**< number of rows in matrix */
   int                   ncols,              /**< number of cols in matrix */
   double*               matrix,             /**< the matrix we want to multiply */
   double*               result              /**< pointer to store the resulting matrix */
   );

#ifdef __cplusplus
}
#endif

#endif
