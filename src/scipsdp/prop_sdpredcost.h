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

/**@file   prop_sdpredcost.h
 * @ingroup PROPAGATORS
 * @brief  reduced cost / dual fixing for SDPs
 * @author Tristan Gally
 *
 *  Propagates bounds
 *
 *  \f$ y_j \leq \ell_j + \frac{v_{CO} - \bar{v}}{\bar{X}_{n+m+j,n+m+j}} \f$,
 *
 *  \f$ y_j \geq u_j - \frac{v_{CO} - \bar{v}}{\bar{X}_{n+j,n+j}} \f$
 *
 *  where \f$\bar{v}\f$ is the value of the current SDP-relaxation, \f$v_{CO}\f$ is the cutoffbound and \f$\bar{X}_{n+m+j,n+m+j}\f$ the value of the
 *  corresponding primal solution.
 */

#ifndef __SCIP_PROP_SDPREDCOST_H_
#define __SCIP_PROP_SDPREDCOST_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the Sdpredcost propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropSdpredcost(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
