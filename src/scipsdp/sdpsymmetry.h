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

/**@file   sdpsymmetry.h
 * @brief  routines for handling/detecting symmetries in SDPs
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPSYMMETRY_H__
#define __SCIP_SDPSYMMETRY_H__

#include "scip/scip.h"

#include "symmetry/struct_sdpsymmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

/** stores information about SDP constraints */
SCIP_EXPORT
SCIP_RETCODE storeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SDPSYM_SDPDATA*       sdpdata             /**< pointer to store SDP symmetry data */
   );

/** frees information about SDP constraints */
SCIP_EXPORT
SCIP_RETCODE freeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SDPSYM_SDPDATA*       sdpdata             /**< pointer to store SDP symmetry data */
   );

/** finds colors for symmetry detection graph */
SCIP_RETCODE findColorsSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SDPSYM_SDPDATA*       sdpdata,            /**< pointer to store SDP symmetry data */
   int                   mincolorval,        /**< value of smallest color */
   SCIP_HASHSET*         fixedvars           /**< hash set storing variable that need to be fixed */
   );

#ifdef __cplusplus
}
#endif

#endif
