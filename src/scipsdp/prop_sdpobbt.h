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

/**@file   prop_sdpobbt.h
 * @ingroup PROPAGATORS
 * @brief  optimization-based bound tightening propagator for semidefinite programs
 * @author Tristan Gally
 *
 * In Optimization-Based Bound Tightening (OBBT), we solve auxiliary SDPs of the form
 * \f[
 *      \min / \max \, \{ y_i \mid y \in SDP' \},
 * \f]
 *
 * where \f$SDP'\f$ is the current SDP relaxation restricted by the primal cutoff constraint \f$c^T y <= z\f$, \f$z\f$ the
 * current cutoff bound. Trivially, the optimal objective value of this LP provides a valid lower/upper bound on
 * variable \f$y_i\f$.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PROP_SDPOBBT_H__
#define __SCIP_PROP_SDPOBBT_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the sdp-obbt propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropSdpObbt(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
