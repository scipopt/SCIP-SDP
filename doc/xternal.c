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

/**@file   xternal.c
 * @brief  main document page
 * @author Tristan Gally
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 *
 * @version 4.2.0
 * @author Marc Pfetsch; Sonja Mars, Lars Schewe, Tristan Gally, Frederic Matter, Christopher Hojny
 * @date 2011-2024
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
