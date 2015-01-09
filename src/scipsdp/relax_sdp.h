/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_sdp.h
 * @ingroup RELAXATORS
 * @brief  SDP relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAXSDP_H__
#define __SCIP_RELAXSDP_H__

#include "scip/scip.h"
#include "sdpi/sdpi_general.h"          // for SDP-Interface

#ifdef __cplusplus
extern "C" {
#endif

/** creates the SDP relaxator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeRelaxSDP(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns pointer to SDP Interface structure */
EXTERN
SCIP_SDPI* SCIPrelaxSdpGetSdpi(
   SCIP_RELAX*           relax               /**< SDP relaxator to get sdpi for */
   );

/** returns optimal objective value of the current SDP relaxation, if the last SDP relaxation was successfully solved*/
EXTERN
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP relaxator to get objective value for */
   SCIP_Bool*            success,            /**< was the last SDP relaxation solved successfully? */
   SCIP_Real*            objval              /**< returns the optimal objective value of the SDP relaxation */
   );

/** returns total number of SDP iterations */
EXTERN
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*            relax               /**< SDP relaxator to get the iterations for */
   );

/** returns number of solved SDP relaxations */
EXTERN
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*            relax               /**< SDP relaxator to get the average iterations for */
   );

#ifdef __cplusplus
}
#endif

#endif
