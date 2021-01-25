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
