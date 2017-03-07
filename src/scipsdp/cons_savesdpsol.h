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

#ifndef __SCIP_CONS_SAVEDSDPSOL_H_
#define __SCIP_CONS_SAVEDSDPSOL_H_

#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** include Savedsdpsol constraint handler */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSavesdpsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create a Savedsdpsol-Cons, i.e. save the current optimal solution for the SDP-relaxation of this node */
EXTERN
SCIP_RETCODE createConsSavesdpsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables and therefore length of sol */
   SCIP_Real*            sol                 /**< optimal solution for SDP-relaxation of this node */
   );

/** for the given cons of type Savedsdpsol returns the previous dual solution vector y, length should start with the length of the array, this
 *  needs to be atleast the number of variables in scip and will be overwritten by this value, if it wasn't sufficient a debugMessage will be thrown
 */
EXTERN
SCIP_RETCODE getDualVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get starting point for */
   SCIP_Real*            sol,                /**< output: previous dual solution vector y */
   int*                  length              /**< input: length of sol-array, output: number of entries in sol-array */
   );

#ifdef __cplusplus
}
#endif

#endif
