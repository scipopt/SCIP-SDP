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

/**@file   cons_sdp.h
 * @brief  constraint handler for SDPs
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 */

#ifndef __SCIP_CONSHDLR_SDP_H__
#define __SCIP_CONSHDLR_SDP_H__

#include "scip/scip.h"

/*
#ifdef __cplusplus
extern "C" {
#endif
*/

class SdpCone;  // lines 20-20

/** creates the handler for sdp constraints and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** method for creating sdp constraints */
EXTERN
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   const SdpCone&        sdpcone             /**< the sdpcone */
   /* @todo pass C-style information instead of cone */
   );

/** gets (a pointer to) the sdpcone of a sdp-constraint */
EXTERN
SCIP_RETCODE getSdpCone(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            conss,              /**< constraint to get data for */
   SdpCone**             sdpcone             /**< output: pointer to the cone */
   );

/** checks feasibility for a single SDP-Cone */
EXTERN
SCIP_RETCODE consCheckSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpCone*              sdpcone,            /**< sdpcone to check positive semidefiniteness for */
   SCIP_SOL*             sol,                /**< the solution to check feasibility for */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the feasibility checking call */
   );

/*
#ifdef __cplusplus
}
#endif
*/

#endif
