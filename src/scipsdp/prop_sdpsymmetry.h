/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2024 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2024 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_sdpsymmetry.h
 * @ingroup PROPAGATORS
 * @brief  propagator for symmetry handling
 * @author Marc Pfetsch
 * @author Thomas Rehn
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_SDPSYMMETRY_H_
#define __SCIP_PROP_SDPSYMMETRY_H_

#include <scip/scip.h>

#include <symmetry/type_sdpsymmetry.h>

#ifdef __cplusplus
extern "C" {
#endif

/** include symmetry propagator */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropSdpSymmetry(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return currently available symmetry group information */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSdpSymmetry(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  npermvars,          /**< pointer to store number of variables for permutations */
   SCIP_VAR***           permvars,           /**< pointer to store variables on which permutations act */
   SCIP_HASHMAP**        permvarmap,         /**< pointer to store hash map of permvars (or NULL) */
   int*                  nperms,             /**< pointer to store number of permutations */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix (or NULL)*/
   int***                permstrans,         /**< pointer to store permutation generators as (npermvars x nperms) matrix (or NULL)*/
   SCIP_Real*            log10groupsize,     /**< pointer to store log10 of group size (or NULL) */
   SCIP_Bool*            binvaraffected,     /**< pointer to store whether binary variables are affected */
   int**                 components,         /**< pointer to store components of symmetry group (or NULL) */
   int**                 componentbegins,    /**< pointer to store begin positions of components in components array (or NULL) */
   int**                 vartocomponent,     /**< pointer to store assignment from variable to its component (or NULL) */
   int*                  ncomponents         /**< pointer to store number of components (or NULL) */
   );

/** return whether orbital fixing is enabled */
SCIP_EXPORT
SCIP_Bool SCIPSDPisOrbitalfixingEnabled(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** return number of the symmetry group's generators */
SCIP_EXPORT
int SCIPgetSdpSymmetryNGenerators(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
