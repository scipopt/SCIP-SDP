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

/**@file   sdpiclock.h
 * @brief  methods for clocks and timing
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SDPICLOCK_H__
#define __SDPICLOCK_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "sdpi/type_sdpiclock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates a clock and initializes it */
SCIP_RETCODE SDPIclockCreate(
   SDPI_CLOCK**          clck                /**< pointer to clock timer */
   );

/** frees a clock */
void SDPIclockFree(
   SDPI_CLOCK**          clck                /**< pointer to clock timer */
   );

/** sets the type of the clock */
void SDPIclockSetType(
   SDPI_CLOCK*           clck,               /**< clock timer */
   SDPI_CLOCKTYPE        clocktype           /**< type of clock */
   );

/** starts measurement of time in the given clock, update the clock's type if it is bound to the default type */
void SDPIclockStart(
   SDPI_CLOCK*           clck                /**< clock timer */
   );

/** stops measurement of time in the given clock */
void SDPIclockStop(
   SDPI_CLOCK*           clck                /**< clock timer */
   );

/** gets the used time of this clock in seconds */
SCIP_Real SDPIclockGetTime(
   SDPI_CLOCK*           clck                /**< clock timer */
   );

#ifdef __cplusplus
}
#endif

#endif
