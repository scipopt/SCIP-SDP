/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2024 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2024 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solveonevarsdp.h
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
   int*                  sdpconstrow,        /**< array of row-indices of constant matrix */
   int*                  sdpconstcol,        /**< array of column-indices of constant matrix */
   SCIP_Real*            sdpconstval,        /**< array of nonzero values of entries of constant matrix */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< array of row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< array of column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< array of nonzero values */
   SCIP_Real             infinity,           /**< infinity value */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            certificatevector,  /**< array to store a certificate (eigen)vector (or NULL if not required) */
   SCIP_Real*            certificatevalue,   /**< array to store a certificate value (or NULL if not required) */
   SCIP_Real*            objval,             /**< pointer to store optimal objective value */
   SCIP_Real*            optval              /**< pointer to store optimal value of variable */
   );

/** solves SDP with one variable and one SDP block - variant for dense constant matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPsolveOneVarSDPDense(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real             obj,                /**< objective coefficient of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   int                   blocksize,          /**< size of the SDP-block */
   SCIP_Real*            fullconstmatrix,    /**< dense full constant matrix */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< array of row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< array of column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< array of nonzero values */
   SCIP_Real             infinity,           /**< infinity value */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            objval,             /**< pointer to store optimal objective value */
   SCIP_Real*            optval              /**< pointer to store optimal value of variable */
   );

#ifdef __cplusplus
}
#endif

#endif
