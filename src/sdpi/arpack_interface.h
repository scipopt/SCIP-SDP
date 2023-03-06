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

/**@file   arpack_interface.h
 * @brief  interface methods for eigenvector computation using arpack
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_ARPACK_INTERFACE_H__
#define __SCIP_ARPACK_INTERFACE_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK */
SCIP_EXPORT
SCIP_RETCODE SCIParpackComputeSmallestEigenvector(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed in column-major form */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< array for eigenvector */
   );

/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK; specialized sparse version for mu A + B */
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
   );

#ifdef __cplusplus
}
#endif

#endif
