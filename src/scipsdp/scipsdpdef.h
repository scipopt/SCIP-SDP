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

/**@file   scipsdpdef.h
 * @brief  main definitions for SCIP-SDP
 * @author Marc Pfetsch
 */

#ifndef __SCIPSDP_DEF_H__
#define __SCIPSDP_DEF_H__

#include "scip/scip.h"


#ifdef __cplusplus
extern "C" {
#endif

#define SCIPSDP_VERSION               401 /**< SCIP-SDP version number (multiplied by 100 to get integer number) */
#define SCIPSDP_SUBVERSION              0 /**< SCIP-SDP sub version number */
#define SCIPSDP_APIVERSION              0 /**< SCIP-SDP API version number */
#define SCIPSDP_COPYRIGHT   "Copyright (C) 2011-2022 TU Darmstadt"

#define SCIPSDPmajorVersion SCIPSDP_VERSION / 100
#define SCIPSDPminorVersion (SCIPSDP_VERSION/10) % 10
#define SCIPSDPtechVersion  SCIPSDP_VERSION % 10

#ifdef __cplusplus
}
#endif

#endif
