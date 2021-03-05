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
 *
 * This is a modified version of clock.c in SCIP in order to make it independent from SCIP_SET. Note that here the time
 * is restarted every time SDPIclockStart() is called.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>
#endif
#include <time.h>

#include "scip/pub_message.h"
#include "blockmemshell/memory.h"

#include "sdpi/sdpiclock.h"
#include "sdpi/struct_sdpiclock.h"


/** converts CPU clock ticks into seconds */
static
SCIP_Real cputime2sec(
   clock_t               cputime             /**< clock ticks for CPU time */
   )
{
   clock_t clocks_per_second;

#if defined(_WIN32) || defined(_WIN64)
   clocks_per_second = 100;
#else
#ifndef CLK_TCK
   clocks_per_second = sysconf(_SC_CLK_TCK);
#else
   clocks_per_second = CLK_TCK;
#endif
#endif

   return (SCIP_Real)cputime / (SCIP_Real)clocks_per_second;
}

/*lint -esym(*,timeval)*/
/*lint -esym(*,gettimeofday)*/

/** converts wall clock time into seconds */
static
SCIP_Real walltime2sec(
   long                  sec,                /**< seconds counter */
   long                  usec                /**< microseconds counter */
   )
{
   return (SCIP_Real)sec + 0.000001 * (SCIP_Real)usec;
}

/** creates a clock and initializes it */
SCIP_RETCODE SDPIclockCreate(
   SDPI_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   assert(clck != NULL);

   SCIP_ALLOC( BMSallocMemory(clck) );

   (*clck)->clocktype = SDPI_CLOCKTYPE_WALL;
   (*clck)->nruns = 0;

   return SCIP_OKAY;
}

/** frees a clock */
void SDPIclockFree(
   SDPI_CLOCK**          clck                /**< pointer to clock timer */
   )
{
   assert(clck != NULL);

   BMSfreeMemory(clck);
}

/** sets the type of the clock */
void SDPIclockSetType(
   SDPI_CLOCK*           clck,               /**< clock timer */
   SDPI_CLOCKTYPE        clocktype           /**< type of clock */
   )
{
   assert( clck != NULL );

   SCIPdebugMessage("setting type of clock %p to %d.\n", (void*)clck, clocktype);

   clck->clocktype = clocktype;
}

/** starts measurement of time in the given clock */
void SDPIclockStart(
   SDPI_CLOCK*           clck                /**< clock timer */
   )
{
#if defined(_WIN32) || defined(_WIN64)
   FILETIME creationtime;
   FILETIME exittime;
   FILETIME kerneltime;
   FILETIME usertime;
#else
   struct timeval tp; /*lint !e86*/
   struct tms now;
#endif

   assert( clck != NULL );
   assert( clck->nruns == 0 );

   SCIPdebugMessage("starting clock %p (type %d).\n", (void*)clck, clck->clocktype);

   switch ( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
      GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
      clck->data.cpuclock.user = - usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L;
#else
      (void)times(&now);
      clck->data.cpuclock.user = - now.tms_utime;
#endif
      break;

   case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
      clck->data.wallclock.sec = - time(NULL);
#else
      gettimeofday(&tp, NULL);
      if( tp.tv_usec > 0 ) /*lint !e115 !e40*/
      {
         clck->data.wallclock.sec = - (tp.tv_sec + 1); /*lint !e115 !e40*/
         clck->data.wallclock.usec = (1000000 - tp.tv_usec); /*lint !e115 !e40*/
      }
      else
      {
         clck->data.wallclock.sec = - tp.tv_sec; /*lint !e115 !e40*/
         clck->data.wallclock.usec = - tp.tv_usec; /*lint !e115 !e40*/
      }
#endif
      break;

   default:
      SCIPerrorMessage("invalid clock type\n");
      SCIPABORT();
   }

   ++clck->nruns;
}

/** stops measurement of time in the given clock */
void SDPIclockStop(
   SDPI_CLOCK*           clck                /**< clock timer */
   )
{
#if defined(_WIN32) || defined(_WIN64)
   FILETIME creationtime;
   FILETIME exittime;
   FILETIME kerneltime;
   FILETIME usertime;
#else
   struct timeval tp; /*lint !e86*/
   struct tms now;
#endif

   assert( clck != NULL );
   assert( clck->nruns == 1 );

   --clck->nruns;

   SCIPdebugMessage("stopping clock %p (type %d)\n", (void*)clck, clck->clocktype);

   switch ( clck->clocktype )
   {
   case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
      GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
      clck->data.cpuclock.user += usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L;
#else
      (void)times(&now);
      clck->data.cpuclock.user += now.tms_utime;
#endif
      break;

   case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
      clck->data.wallclock.sec += time(NULL);
#else
      gettimeofday(&tp, NULL);
      if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
      {
         clck->data.wallclock.sec += (tp.tv_sec + 1); /*lint !e115 !e40*/
         clck->data.wallclock.usec -= (1000000 - tp.tv_usec); /*lint !e115 !e40*/
      }
      else
      {
         clck->data.wallclock.sec += tp.tv_sec; /*lint !e115 !e40*/
         clck->data.wallclock.usec += tp.tv_usec; /*lint !e115 !e40*/
      }
#endif
      break;

   default:
      SCIPerrorMessage("invalid clock type\n");
      SCIPABORT();
   }
}

/** gets the used time of this clock in seconds */
SCIP_Real SDPIclockGetTime(
   SDPI_CLOCK*           clck                /**< clock timer */
   )
{
   SCIP_Real result = 0.0;

   assert( clck != NULL );

   SCIPdebugMessage("getting time of clock %p (type %d, nruns=%d)\n", (void*)clck, clck->clocktype, clck->nruns);

   if ( clck->nruns == 0 )
   {
      /* the clock is not running: convert the clocks timer into seconds */
      switch ( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
         result = cputime2sec(clck->data.cpuclock.user);
         break;
      case SCIP_CLOCKTYPE_WALL:
         result = walltime2sec(clck->data.wallclock.sec, clck->data.wallclock.usec);
         break;
      default:
         SCIPerrorMessage("invalid clock type\n");
         SCIPABORT();
      }
   }
   else
   {
#if defined(_WIN32) || defined(_WIN64)
      FILETIME creationtime;
      FILETIME exittime;
      FILETIME kerneltime;
      FILETIME usertime;
#else
      struct timeval tp; /*lint !e86*/
      struct tms now;
#endif

      /* the clock is currently running: we have to add the current time to the clocks timer */
      switch ( clck->clocktype )
      {
      case SCIP_CLOCKTYPE_CPU:
#if defined(_WIN32) || defined(_WIN64)
         GetProcessTimes(GetCurrentProcess(), &creationtime, &exittime, &kerneltime, &usertime);
         result = cputime2sec(clck->data.cpuclock.user + usertime.dwHighDateTime * 42950 + usertime.dwLowDateTime / 100000L);
#else
         (void)times(&now);
         result = cputime2sec(clck->data.cpuclock.user + now.tms_utime);
#endif
         break;
      case SCIP_CLOCKTYPE_WALL:
#if defined(_WIN32) || defined(_WIN64)
         result = walltime2sec(clck->data.wallclock.sec + time(NULL), 0);
#else
         gettimeofday(&tp, NULL);
         if( tp.tv_usec + clck->data.wallclock.usec > 1000000 ) /*lint !e115 !e40*/
            result = walltime2sec(clck->data.wallclock.sec + tp.tv_sec + 1, (clck->data.wallclock.usec - 1000000) + tp.tv_usec); /*lint !e115 !e40*/
         else
            result = walltime2sec(clck->data.wallclock.sec + tp.tv_sec, clck->data.wallclock.usec + tp.tv_usec); /*lint !e115 !e40*/
#endif
         break;
      default:
         SCIPerrorMessage("invalid clock type\n");
         SCIPABORT();
      }
   }

   return result;
}
