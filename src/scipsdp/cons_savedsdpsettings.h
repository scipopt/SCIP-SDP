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

/**@file   cons_savedsdpsettings.h
 * @brief  constraint handler for saving SDP settings
 * @author Tristan Gally
 *
 * A constraint that is always feasible which can be used to save and recover settings used
 * to solve the SDP-relaxation at the current node.
 */

#ifndef __SCIP_CONS_SAVEDSDPSETTINGS_H_
#define __SCIP_CONS_SAVEDSDPSETTINGS_H_

#include "scip/scip.h"
#include "sdpi/type_sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** include Savedsdpsettings constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrSavedsdpsettings(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create a savedsdpsettings constraint, i.e. save the current settings for the SDP-relaxation of this node */
SCIP_EXPORT
SCIP_RETCODE createConsSavedsdpsettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_SDPSOLVERSETTING settings            /**< settings to save */
   );

/** get the settings used to solve the SDP relaxation in this node */
SCIP_EXPORT
SCIP_SDPSOLVERSETTING SCIPconsSavedsdpsettingsGetSettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get starting point for */
   );

#ifdef __cplusplus
}
#endif

#endif
