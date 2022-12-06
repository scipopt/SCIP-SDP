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
 * SCIP-SDP allows to solve MISDPs using a nonlinear branch-and-bound approach or a linear programming cutting-plane
 * approach. In the first case (the default), the semidefinite programming (SDP) relaxations are solve using
 * interior-point SDP-solvers. In the second case, cutting planes based on eigenvector are generated.  SCIP-SDP is based
 * on the branch-and-cut framework <A HREF="https://scipopt.org">SCIP</A>.  In addition to providing a constraint
 * handler for SDP-constraints and a relaxator to solve continuous SDP-relaxations using interior-point solvers,
 * SCIP-SDP adds several heuristics and propagators to SCIP.
 *
 * The MISDPs can be read in using either an extended SDPA-format or the CBF-format. There is also an interface for
 * Matlab/Octave on <A HREF="https://github.com/scipopt/MatlabSCIPInterface">GitHub</A>.
 *
 * To use the nonlinear branch-and-bound approach one of the following SDP-solvers needs to be installed:
 *
 * - DSDP
 * - SDPA
 * - MOSEK
 *
 * Mixed-integer semidefinite programs are sometimes numerically challenging to solve. One reason is that the Slater
 * condition may not hold for the SDP-relaxations of some of the nodes. SCIP-SDP implements several methods that try to
 * recover from a failure to accurately solve the relaxation.
 */
