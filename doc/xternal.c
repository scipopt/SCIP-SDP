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

/**@file   xternal.c
 * @brief  main document page
 * @author Tristan Gally
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 *
 * @version 4.1.0
 * @author Marc Pfetsch; Sonja Mars, Lars Schewe, Tristan Gally, Frederic Matter
 * @date 2011-2022
 *
 * SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite programs (MISDPs) of the form
 *
 *   \f{equation*}{
 * 	\begin{aligned}
 *      \inf \quad & b^T y && \\
 *      \mbox{s.t.} \quad & \sum_{i = 1}^m A_i y_i - A_0 \succeq 0,&& \\
 *	& y_i \in \mathbb{Z} && \forall \ i \in \mathcal{I}.
 *	\end{aligned}
 *   \f}
 *
 * It combines the branch-and-bound framework of SCIP with interior-point SDP-solvers to solve MISDPs using either a
 * nonlinear branch-and-bound approach or an outer-approximation-based cutting-plane approach. In addition to providing
 * a constraint handler for SDP-constraints and a relaxator to solve continuous SDP-relaxations using interior-point
 * solvers, SCIPSDP adds several heuristics and propagators to SCIP. The MISDPs can be read in using either an extended
 * SDPA-format or the CBF-format. To use the nonlinear branch-and-bound approach one of the following SDP-solvers needs
 * to be installed:
 *
 * - DSDP
 * - SDPA
 * - MOSEK
 *
 * The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main
 * purposes of the code is handling cases where the Slater condition does not hold using a penalty approach. However, in
 * some cases the SDP-solvers may still fail because of numerical difficulties or even return wrong results, which cannot
 * be compensated. For this purpose there is the possibility to check the Slater condition for the primal and dual problem
 * before the solution of each SDP by setting a SCIP parameter, for details see the parameters tab.
 */
