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

/**@file   compute_symmetry.h
 * @brief  interface for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_COMPUTE_SYMMETRY_H_
#define __SCIP_COMPUTE_SYMMETRY_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "symmetry/struct_sdpsymmetry.h"

/** return whether symmetry can be computed */
SCIP_EXPORT
SCIP_Bool SDPSYMcanComputeSymmetry(void);

/** return name of external program used to compute generators */
SCIP_EXPORT
const char* SDPSYMsymmetryGetName(void);

/** return description of external program used to compute generators */
SCIP_EXPORT
const char* SDPSYMsymmetryGetDesc(void);

/** compute generators of symmetry group */
SCIP_EXPORT
SCIP_RETCODE SDPSYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SDPSYM_MATRIXDATA*    matrixdata,         /**< data for MIP matrix */
   SDPSYM_EXPRDATA*      exprdata,           /**< data for nonlinear constraints */
   SDPSYM_SDPDATA*       sdpdata,            /**< data for SDP constraints */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of size of group */
   SCIP_HASHSET*         fixedvars,          /**< hash set storing variable that need to be fixed */
   int                   nfixedvars          /**< number of fixed variables */
   );

#ifdef __cplusplus
}
#endif

#endif
