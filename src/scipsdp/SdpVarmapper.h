/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   SdpVarmapper.h
 * @brief  class that maps SCIP variables to SDP indices (the SCIP variables are given SDP indices in the order in which they were inserted)
 * @author Tristan Gally
 */

#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__

#include "scip/scip.h"
#include "scip/type_misc.h" /* for SCIP Hashmap */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct Sdpvarmapper SdpVarmapper;

/** creates the SDP Varmapper */
EXTERN
SCIP_RETCODE SdpVarmapperCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper**        varmapper,          /**< Pointer to the Varmapper that should be created */
   int                   size                /**< initial size of the sciptosdp-hashmap */
   );

/** frees the SDP Varmapper */
EXTERN
SCIP_RETCODE SdpVarmapperFree(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper**        varmapper          /**< Pointer to the Varmapper that should be freed */
   );

/** adds the given variables (if not already existent) to the end of the Varmapper */
EXTERN
SCIP_RETCODE SdpVarmapperAddVars(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper,         /**< Varmapper to add variables to */
   int                   nvars,             /**< number of variables to add to the varmapper */
   SCIP_VAR**            vars               /**< SCIP variables to add to the varmapper */
   );

/** adds the given variable (if not already existant) to the Varmapper at the given position */
EXTERN
SCIP_RETCODE SdpVarmapperInsertVar(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper,         /**< Varmapper to add variables to */
   SCIP_VAR*            var,                /**< SCIP variable to add to the varmapper */
   int                  pos                 /**< position where the variable should be added */
   );

/** gets the number of variables */
EXTERN
int SdpVarmapperGetNVars(
   SdpVarmapper*         varmapper          /**< Varmapper to get number of variables for */
   );

/** Is the given SCIP variable included in the varmapper? */
EXTERN
SCIP_Bool SdpVarmapperExistsSCIPvar(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   SCIP_VAR*             var                /**< SCIP variables to get sdp index for */
   );

/** gets the sdp index for the given SCIP variable */
EXTERN
int SdpVarmapperGetSdpIndex(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   SCIP_VAR*             var                /**< SCIP variables to get sdp index for */
   );

/** gets the corresponding SCIP variable for the given sdp variable index */
EXTERN
SCIP_VAR* SdpVarmapperGetSCIPvar(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   int                   ind                /**< index of the sdp variable */
   );

/** removes the variable for the given Sdp index from the varmapper, decreasing the indices of all later variables by 1 */
EXTERN
SCIP_RETCODE SdpVarmapperRemoveSdpIndex(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   int                   ind                /**< index of the sdp variable */
   );

/** swaps all SCIP variables for their transformed counterparts */
EXTERN
SCIP_RETCODE SdpVarmapperTransform(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper          /**< pointer to the Varmapper that should be transformed */
   );

/** clones the varmapper in the second argument to the varmapper in the third argument */
EXTERN
SCIP_RETCODE SdpVarmapperClone(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         oldmapper,         /**< pointer to the Varmapper that should be cloned */
   SdpVarmapper*         newmapper          /**< pointer to the Varmapper that should become a clone of the other one */
   );

#ifdef __cplusplus
}
#endif

#endif
