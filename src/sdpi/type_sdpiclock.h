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

/**@file   type_sdpiclock.h
 * @brief  type definitions for clocks and timing issues
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_SDPICLOCK_H__
#define __TYPE_SDPICLOCK_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

enum SDPI_ClockType
{
   SDPI_CLOCKTYPE_CPU     = 1,          /**< use CPU clock */
   SDPI_CLOCKTYPE_WALL    = 2           /**< use wall clock */
};
typedef enum SDPI_ClockType SDPI_CLOCKTYPE;       /**< clock type to use */

typedef struct SDPI_Clock SDPI_CLOCK;             /**< clock timer */
typedef struct SDPI_CPUClock SDPI_CPUCLOCK;       /**< CPU clock counter */
typedef struct SDPI_WallClock SDPI_WALLCLOCK;     /**< wall clock counter */

#ifdef __cplusplus
}
#endif

#endif
