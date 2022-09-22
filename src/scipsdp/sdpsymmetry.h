/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sdpsymmetry.h
 * @brief  routines for handling/detecting symmetries in SDPs
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SDPSYMMETRY_H__
#define __SCIP_SDPSYMMETRY_H__

#include "scip/scip.h"

#include "scipsdp/struct_sdpsymmetry.h"

#ifdef __cplusplus
extern "C" {
#endif

/** stores information about SDP constraints */
SCIP_EXPORT
SCIP_RETCODE storeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata             /**< pointer to store SDP symmetry data */
   );

/** frees information about SDP constraints */
SCIP_EXPORT
SCIP_RETCODE freeSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata             /**< pointer to store SDP symmetry data */
   );

/** finds colors for symmetry detection graph */
SCIP_RETCODE findColorsSDPSymmetryData(
   SCIP*                 scip,               /**< SCIP data structure */
   SYM_SDPDATA*          sdpdata,            /**< pointer to store SDP symmetry data */
   int                   mincolorval         /**< value of smallest color */
   );

#ifdef __cplusplus
}
#endif

#endif
