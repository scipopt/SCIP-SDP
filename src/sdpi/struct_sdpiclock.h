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

/**@file   struct_sdpiclock.h
 * @brief  datastructures for clocks and timing issues
 * @author Tobias Achterberg
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SDPICLOCK_H__
#define __STRUCT_SDPICLOCK_H__


#if defined(_WIN32) || defined(_WIN64)
#include <time.h>
#else
#include <sys/times.h>
#endif

#include "scip/def.h"
#include "scip/type_clock.h"

#ifdef __cplusplus
extern "C" {
#endif

/** CPU clock counter */
struct SDPI_CPUClock
{
   clock_t               user;               /**< clock ticks for user CPU time */
};

/** wall clock counter */
struct SDPI_WallClock
{
   long                  sec;                /**< seconds counter */
   long                  usec;               /**< microseconds counter */
};

/** clock timer */
struct SDPI_Clock
{
   union
   {
      SDPI_CPUCLOCK      cpuclock;           /**< CPU clock counter */
      SDPI_WALLCLOCK     wallclock;          /**< wall clock counter */
   } data;
   int                   nruns;              /**< number of SCIPclockStart() calls without SCIPclockStop() calls */
   SDPI_CLOCKTYPE        clocktype;          /**< current type of clock used */
};

#ifdef __cplusplus
}
#endif

#endif
