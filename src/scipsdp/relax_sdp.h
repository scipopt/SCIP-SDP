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
#include "sdpi/sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the SDP relaxator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown. Note: if a variable is either fixed or unbounded in the dual
 *  problem, a zero will be returned for the non-existent primal variable.
 */
EXTERN
SCIP_RETCODE SCIPrelaxSdpGetPrimalBoundVars(
   SCIP_RELAX*           relax,              /**< SDP relaxator to information for */
   SCIP_Real*            lbvars,             /**< returns the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< returns the variables corresponding to upper bounds in the dual problems */
   int*                  arraylength         /**< input: length of lbvars and ubvars
                                              *   output: number of elements inserted into lbvars/ubvars (or needed length if it wasn't sufficient) */
   );

/** returns optimal objective value of the current SDP relaxation, if the last SDP relaxation was successfully solved */
EXTERN
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP relaxator to get objective value for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP relaxation solved successfully */
   SCIP_Real*            objval              /**< returns the optimal objective value of the SDP relaxation */
   );

/** returns values of all variables in the solution of the current SDP relaxation, if the last SDP relaxation was successfully solved */
EXTERN
SCIP_RETCODE SCIPrelaxSdpGetRelaxSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax,              /**< SDP relaxator to get solution for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP relaxation solved successfully */
   SCIP_Real*            solarray,           /**< array to insert the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray, if this is less than nvars, it will be overwritten with the needed length and a
                                               *  debug message is thrown */
   );

/** get the number of the SCIP-node to which the current SDP solution belongs */
EXTERN
long int SCIPrelaxSdpGetSdpNode(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   );

/** was the original problem solved for the last SDP-Node (or a penalty formulation) ? */
EXTERN
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP relaxator to get solution for */
   );

/** returns total number of SDP iterations */
EXTERN
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the iterations for */
   );

/** returns number of solved SDP relaxations */
EXTERN
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*           relax               /**< SDP relaxator to get the number of calls for */
   );

#ifdef __cplusplus
}
#endif

#endif
