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

/**@file   relax_sdp.h
 * @ingroup RELAXATORS
 * @brief  SDP-relaxator
 * @author Sonja Mars
 * @author Tristan Gally
 *
 * Relaxator to solve semidefinite programs of the form
 * \f{eqnarray*}{
 *    \min & & b^T y \\
 *    \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i \succeq 0 \quad \forall i \leq m \\
 *    & & Dy \geq d \\
 *    & & \ell \leq y \leq u
 * \f}
 * for symmetric matrices \f$ A_j^i \in S_{k_i} \f$ and a matrix \f$ D \in \mathbb{R}^{k_0 \times n}. \f$
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RELAXSDP_H__
#define __SCIP_RELAXSDP_H__

#include "scip/scip.h"
#include "sdpi/sdpi.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the SDP-relaxator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** computes analytic centers of primal and dual feasible set and saves them in relaxdata
 * @note This function should be called at the end of the root node (or at least after the solving stage starts and before the first non-root node).
 */
SCIP_EXPORT
SCIP_RETCODE SCIPrelaxSdpComputeAnalyticCenters(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RELAX*           relax               /**< SDP-relaxator to compute analytic centers for */
   );

/** gets the primal variables corresponding to the lower and upper variable-bounds in the dual problem
 *
 *  The last input should specify the length of the arrays. If this is less than the number of variables, the needed
 *  length will be returned and a debug message thrown.
 *
 *  @note If a variable is either fixed or unbounded in the dual
 *  problem, a zero will be returned for the non-existent primal variable.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPrelaxSdpGetPrimalBoundVars(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get information for */
   SCIP_VAR**            vars,               /**< variables to get bounds for */
   int                   nvars,              /**< number of variables */
   SCIP_Real*            lbvars,             /**< pointer to store the values of the variables corresponding to lower bounds in the dual problems */
   SCIP_Real*            ubvars,             /**< pointer to store the values of the variables corresponding to upper bounds in the dual problems */
   SCIP_Bool*            success             /**< pointer to store success (may fail if problem is infeasible or all variables are fixed) */
   );

/** returns optimal objective value of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_EXPORT
SCIP_RETCODE SCIPrelaxSdpRelaxVal(
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get objective value for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            objval              /**< pointer to store the optimal objective value of the SDP-relaxation */
   );

/** returns values of all variables in the solution of the current SDP-relaxation if the last SDP-relaxation was successfully solved */
SCIP_EXPORT
SCIP_RETCODE SCIPrelaxSdpGetRelaxSol(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax,              /**< SDP-relaxator to get solution for */
   SCIP_Bool*            success,            /**< pointer to store whether the last SDP-relaxation was solved successfully */
   SCIP_Real*            solarray,           /**< pointer to store the solution, this has to be at least length nvars */
   int*                  sollength           /**< length of the solarray, if this is less than nvars, it will be overwritten with the needed length and a
                                               *  debug message is thrown */
   );

/** get the number of the SCIP-node which the current SDP solution belongs to */
SCIP_EXPORT
SCIP_Longint SCIPrelaxSdpGetSdpNode(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   );

/** Was the original problem solved for the last SDP-node (or a penalty or probing formulation) ? */
SCIP_EXPORT
SCIP_Bool SCIPrelaxSdpSolvedOrig(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   );

/** Was the last probing SDP solved successfully ? */
SCIP_EXPORT
SCIP_Bool SCIPrelaxSdpSolvedProbing(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get solution for */
   );

/** returns whether the last solved problem was feasible */
SCIP_EXPORT
SCIP_Bool SCIPrelaxSdpIsFeasible(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get feasibility for */
   );

/** returns whether the last solved problem was unbounded */
SCIP_EXPORT
SCIP_Bool SCIPrelaxSdpIsUnbounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to check for unboundedness */
   );

/** returns time in optimization of solver */
SCIP_EXPORT
SCIP_Real SCIPrelaxSdpGetOptTime(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the iterations for */
   );

/** returns total number of SDP-iterations */
SCIP_EXPORT
int SCIPrelaxSdpGetNIterations(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the iterations for */
   );

/** returns number of SDPs solved by SDP-solver (including multiple calls for penalty formulation etc.) */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of solved SDP-relaxations */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpInterfaceCalls(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations solved with fastest settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations solved with medium settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpMedium(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations solved with stable settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations solved with penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations unsolved even when using a penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSdpUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get the number of calls for */
   );

/** returns number of SDP-relaxations for which dual Slater condition held */
SCIP_EXPORT
int SCIPrelaxSdpGetNdualSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which dual Slater condition failed */
SCIP_EXPORT
int SCIPrelaxSdpGetNdualSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which dual Slater condition showed infeasibility */
SCIP_EXPORT
int SCIPrelaxSdpGetNdualSlaterInfeasible(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which dual Slater condition could not be determined */
SCIP_EXPORT
int SCIPrelaxSdpGetNdualSlaterUnknown(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which primal Slater condition held */
SCIP_EXPORT
int SCIPrelaxSdpGetNprimalSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which primal Slater condition failed */
SCIP_EXPORT
int SCIPrelaxSdpGetNprimalSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations for which primal Slater condition could not be determined */
SCIP_EXPORT
int SCIPrelaxSdpGetNprimalSlaterUnknown(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHolds(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with fastest settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHoldsFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with stable settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHoldsStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, solved with penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHoldsPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, for which an infeasible lower bound could be computed */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHoldsBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition holding for primal and dual, unsolved even when using a penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterHoldsUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFails(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with fast settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFailsFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with stable settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFailsStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, solved with penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFailsPenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, for which an infeasible lower bound could be computed */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFailsBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slater condition failing for primal or dual, unsolved even when using a penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterFailsUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with fast settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterInfeasibleFast(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with stable settings */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterInfeasibleStable(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, solved with penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterInfeasiblePenalty(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, for which an infeasible lower bound could be computed */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterInfeasibleBounded(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns number of SDP-relaxations with Slatercheck showing infeasibility, unsolved even when using a penalty formulation */
SCIP_EXPORT
int SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(
   SCIP_RELAX*           relax               /**< SDP-relaxator to get number for */
   );

/** returns solving time in SDP solver */
SCIP_EXPORT
SCIP_Real SCIPrelaxSdpGetSolvingTime(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_RELAX*           relax               /**< SDP-relaxator to get timer for */
   );

#ifdef __cplusplus
}
#endif

#endif
