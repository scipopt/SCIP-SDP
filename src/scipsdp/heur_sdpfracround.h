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

/**@file heur_sdpfracround.h
 * @ingroup PRIMALHEURISTICS
 * @brief  fractional rounding heuristic for SDPs
 * @author Marc Pfetsch
 *
 * Fractional rounding heuristic for SDPs: Takes the solution of the SDP-relaxation and iteratively rounds all integer
 * variables. After each rounding the problem is propagated. If the SDP includes continuous variables, the remaining SDP
 * after the fixings is solved again.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_SDPFRACROUND_H__
#define __SCIP_HEUR_SDPFRACROUND_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the fractional rounding heuristic for SDPs and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeHeurSdpFracround(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
