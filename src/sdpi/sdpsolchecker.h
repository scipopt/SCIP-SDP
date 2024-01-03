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

/**@file   sdpsolchecker.h
 * @brief  checks a given SDP solution for feasibility
 * @author Tristan Gally
 *
 * Given a solution, an SDP instance and a feasibility tolerance, checks whether
 * the smallest eigenvalue is >= -feastol for a given feasibility tolerance.
 */

#ifndef __SCIP_SDPSOLCHECKER_H__
#define __SCIP_SDPSOLCHECKER_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Given a solution, an SDP instance and a feasibility tolerance, checks whether
 *  the smallest eigenvalue is >= -feastol for a given feasibility tolerance.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpSolcheckerCheck(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   const int*            sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int* const*           sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int* const*           sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real* const*     sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int* const*           sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int* const*           sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int** const*          sdprow,             /**< pointer to the row-indices for each block and variable */
   int** const*          sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real** const*    sdpval,             /**< values of SDP-constraint mmatrix entries (may be NULL if sdpnnonz = 0) */
   int* const*           indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   const int*            nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   const int*            blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   const int*            lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   const SCIP_Real*      lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   const SCIP_Real*      lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   const int*            lpbeg,              /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            lpind,              /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      lpval,              /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      solvector,          /**< values of all variables (including fixed ones) in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
);

/** Given a solution, an SDP instance and a feasibility tolerance, checks whether the smallest eigenvalue is >= -feastol
 *  for a given feasibility tolerance and returns maximal absolute violation and sum of absolute violations for bounds,
 *  linear constraints and SDP constraints for the corresponding dual problem.
 *
 * Note: Should not be called if solution is a certificate of primal infeasiblity!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpSolcheckerCheckAndGetViolDual(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   const int*            sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int* const*           sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int* const*           sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real* const*     sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int* const*           sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                              *   the number of entries of sdp row/col/val [i][j] */
   int* const*           sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                              *   in the i-th block */
   int** const*          sdprow,             /**< pointer to the row-indices for each block and variable */
   int** const*          sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real** const*    sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int* const*           indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                              *   be removed, otherwise it gives the number of indices removed before this */
   const int*            nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   const int*            blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   const int*            lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   const SCIP_Real*      lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   const SCIP_Real*      lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   const int*            lpbeg,              /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            lpind,              /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      lpval,              /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      solvector,          /**< values of all variables (including fixed ones) in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Real*            maxabsviolbnds,     /**< pointer to store maximal absolute violation of variable bounds  */
   SCIP_Real*            sumabsviolbnds,     /**< pointer to store sum of absolute violations of variable bounds */
   SCIP_Real*            maxabsviolcons,     /**< pointer to store maximal absolute violation of linear constraints */
   SCIP_Real*            sumabsviolcons,     /**< pointer to store sum of absolute violations of linear constraints */
   SCIP_Real*            maxabsviolsdp,      /**< pointer to store maximal absolute violation of SDP constraints */
   SCIP_Real*            sumabsviolsdp,      /**< pointer to store sum of absolute violations of SDP constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
);

/** Given a solution, an SDP instance and a feasibility tolerance, returns maximal absolute violation and sum of
 *  absolute violations for bounds, linear constraints and SDP constraints for the corresponding primal problem. The
 *  solution should be given for the actual problem given to the SDP solver, i.e., after all fixed variables have been
 *  eliminated.
 *
 * Note: Should not be called if solution is a certificate of dual infeasiblity!
 */
SCIP_EXPORT
SCIP_RETCODE SCIPsdpSolcheckerCheckAndGetViolPrimal(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective coefficients of variables (in dual problem) */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   const int*            inputtomosekmapper, /**< entry i gives the index of input variable i in MOSEK (starting from 0) or
                                              *   -j (j=1, 2, ..., nvars - nactivevars) if the variable is fixed, the value and objective value of
                                              *   this fixed variable can be found in entry j-1 of fixedval/obj */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   const int*            sdpnblockvars,      /**< number of variables that exist in each block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-blocks AFTER FIXINGS */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries  of sdpconst row/col/val [i] AFTER FIXINGS */
   int* const*           sdpconstrow,        /**< pointers to row-indices for each block AFTER FIXINGS*/
   int* const*           sdpconstcol,        /**< pointers to column-indices for each block AFTER FIXINGS */
   SCIP_Real* const*     sdpconstval,        /**< pointers to the values of the nonzeros for each block AFTER FIXINGS */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int* const*           sdpnblockvarnonz,   /**< entry [i][j] gives the number of nonzeros for block i and variable j, this is exactly
                                        *   the number of entries of sdp row/col/val [i][j] */
   int* const*           sdpvar,             /**< sdpvar[i][j] gives the sdp-index of the j-th variable (according to the sorting for row/col/val)
                                        *   in the i-th block */
   int** const*          sdprow,             /**< pointer to the row-indices for each block and variable */
   int** const*          sdpcol,             /**< pointer to the column-indices for each block and variable */
   SCIP_Real** const*    sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int* const*           indchanges,         /**< changes needed to be done to the indices, if indchanges[block][ind]=-1, then the index can
                                        *   be removed, otherwise it gives the number of indices removed before this */
   const int*            nremovedinds,       /**< the number of rows/cols to be fixed for each block */
   const int*            blockindchanges,    /**< block indizes will be modified by these, see indchanges */
   int                   nremovedblocks,     /**< number of empty blocks that should be removed */
   int                   nlpcons,            /**< number of active (at least two nonzeros) LP-constraints */
   const int*            lpindchanges,       /**< array for the number of LP-constraints removed before the current one (-1 if removed itself) */
   const SCIP_Real*      lplhs,              /**< left-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   const SCIP_Real*      lprhs,              /**< right-hand sides of active LP-rows after fixings (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint-matrix */
   const int*            lpbeg,              /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            lpind,              /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      lpval,              /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      solvector,          /**< values of all scalar variables in the solution that should be checked */
   SCIP_Real* const*     solmatrices,        /**< values of all matrix variables in the solution that should be checked */
   SCIP_Real             feastol,            /**< feasibility tolerance to check feasibility for */
   SCIP_Real             epsilon,            /**< tolerance used to check for fixed variables */
   SCIP_Real*            maxabsviolbnds,     /**< pointer to store maximal absolute violation of variable bounds  */
   SCIP_Real*            sumabsviolbnds,     /**< pointer to store sum of absolute violations of variable bounds */
   SCIP_Real*            maxabsviolcons,     /**< pointer to store maximal absolute violation of linear constraints */
   SCIP_Real*            sumabsviolcons,     /**< pointer to store sum of absolute violations of linear constraints */
   SCIP_Real*            maxabsviolsdp,      /**< pointer to store maximal absolute violation of SDP constraints */
   SCIP_Real*            sumabsviolsdp,      /**< pointer to store sum of absolute violations of SDP constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether solution is feasible */
);

#ifdef __cplusplus
}
#endif

#endif
