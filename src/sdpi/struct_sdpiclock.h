/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2021 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2021 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
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
