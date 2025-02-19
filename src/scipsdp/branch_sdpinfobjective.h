/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2025 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2025 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_sdpinfobjective.h
 * @ingroup BRANCHINGRULES
 * @brief  combined infeasibility and absolute objective branching rule for SCIP-SDP
 * @author Tristan Gally
 *
 * Branch on variable with highest product of fractionality/integral-infeasibility and absolute objective value in the SDP.
 *
 * Will do nothing for continuous variables, since these are what the external callbacks of the SCIP branching rules are for.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_SDPINFOBJECTIVE_H__
#define __SCIP_BRANCH_SDPINFOBJECTIVE_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the SDP combined infeasibility and absolute objective branching rule and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBranchruleSdpinfobjective(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
