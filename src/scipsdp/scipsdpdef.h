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

#define SCIPSDP_VERSION               400 /**< SCIP-SDP version number (multiplied by 100 to get integer number) */
#define SCIPSDP_SUBVERSION              0 /**< SCIP-SDP sub version number */
#define SCIPSDP_APIVERSION              0 /**< SCIP-SDP API version number */
#define SCIPSDP_COPYRIGHT   "Copyright (C) 2011-2021 TU Darmstadt"

#define SCIPSDPmajorVersion SCIPSDP_VERSION / 100
#define SCIPSDPminorVersion (SCIPSDP_VERSION/10) % 10
#define SCIPSDPtechVersion  SCIPSDP_VERSION % 10

#ifdef __cplusplus
}
#endif

#endif
