/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_savesdpsol.h
 * @brief  constraint handler for saving SDP solutions in nodes
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

/** include Savesdpsol constraint handler */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrSavesdpsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** create a Savesdpsol-Cons, i.e. save the current optimal solution for the SDP-relaxation of this node */
SCIP_EXPORT
SCIP_RETCODE createConsSavesdpsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Longint          node,               /**< index of the node the solution belongs to */
   SCIP_SOL*             sol,                /**< optimal solution for SDP-relaxation of this node */
   SCIP_Real             maxprimalentry,     /**< maximum absolute value of primal matrix */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz,   /**< primal matrix X as starting point for the solver: number of nonzeros for each block,
                                               *  also length of corresponding row/col/val-arrays; or NULL if nblocks = 0 */
   int**                 startXrow,          /**< primal matrix X as starting point for the solver: row indices for each block;
                                               *  or NULL if nblocks = 0 */
   int**                 startXcol,          /**< primal matrix X as starting point for the solver: column indices for each block;
                                               *  or NULL if nblocks = 0 */
   SCIP_Real**           startXval          /**< primal matrix X as starting point for the solver: values for each block;
                                               *  or NULL if nblocks = 0 */
   );

/** for the given cons of type Savesdpsol returns the node the information belongs to */
SCIP_EXPORT
SCIP_Longint SCIPconsSavesdpsolGetNodeIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get starting point for */
   );

/** for the given cons of type Savesdpsol returns the previous dual solution vector y */
SCIP_EXPORT
SCIP_SOL* SCIPconsSavesdpsolGetDualVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get starting point for */
   );

/** for the given cons of type Savesdpsol returns the maximum entry of primal solution X */
SCIP_EXPORT
SCIP_Real SCIPconsSavesdpsolGetMaxPrimalEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint to get maximum primal entry for */
   );

/** for the given cons of type Savesdpsol returns the number of nonzeros for each block of previous primal solution X */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSavesdpsolGetPrimalMatrixNonzeros(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get maximum primal entry for */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz    /**< input: allocated memory for startXrow/col/val
                                               *  output: length of startXrow/col/val */
   );

/** for the given cons of type Savesdpsol returns the previous primal solution X */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSavesdpsolGetPrimalMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get maximum primal entry for */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz,   /**< input: allocated memory for startXrow/col/val
                                               *  output: length of startXrow/col/val */
   int**                 startXrow,          /**< pointer to store pointer to row indices of X */
   int**                 startXcol,          /**< pointer to store pointer to column indices of X */
   SCIP_Real**           startXval           /**< pointer to store pointer to values of X */
   );

#ifdef __cplusplus
}
#endif

#endif
