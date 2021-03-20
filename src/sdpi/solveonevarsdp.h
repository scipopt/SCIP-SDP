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

/**@file   solveonvarsdp.h
 * @brief  Solve SDP with one variable
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SOLVEONEVARSDP_H__
#define __SCIP_SOLVEONEVARSDP_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** solves SDP with one variable and one SDP block */
SCIP_EXPORT
SCIP_RETCODE SCIPsolveOneVarSDP(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real             obj,                /**< objective coefficient of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   int                   blocksize,          /**< size of the SDP-block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrix of the SDP-block */
   int*                  sdpconstrow,        /**< pointer to row-indices of constant matrix */
   int*                  sdpconstcol,        /**< pointer to column-indices of constant matrix */
   SCIP_Real*            sdpconstval,        /**< pointer to values of entries of constant matrix */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< pointer to the row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< pointer to the column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< pointer to the nonzero values */
   SCIP_Real             infinity,           /**< infinity value */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            objval,             /**< pointer to store optimal objective value */
   SCIP_Real*            optval              /**< pointer to store optimal value of variable */
   );

#ifdef __cplusplus
}
#endif

#endif
