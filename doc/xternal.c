/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2016 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2016 Zuse Institute Berlin                             */
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
 * @version 2.1.0
 * @author Tristan Gally, Marc Pfetsch; Sonja Mars, Lars Schewe
 * @date 2011-2016
 *
 * SCIP-SDP is a plugin for SCIP to solve mixed integer semidefinite programs (MISDPs). It combines the branch-and-bound
 * framework of SCIP with interior-point SDP-solvers. It provides the data handling, some presolving and propagation as
 * well as a reader for a modified sparse SDPA-format with additional lines for integrality constraints (see 
 * data_format.txt). It is possible to solve the resulting SDP-relaxations using a linear approximation procedure, but for
 * full functionality one of the following SDP-solvers needs to be installed:
 *
 * - DSDP
 * - SDPA
 *
 * The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main
 * purposes of the code is ensuring that the slater condition is not harmed by fixing variables in the branch-and-bound
 * process. However in some cases the combination of variable fixings and specific linear or semidefinite constraints might
 * still lead to relaxations for which the Slater condition no longer holds. In this case the SDP-solvers may be unable to
 * solve the relaxations or even return wrong results, which cannot be compensated. For this purpose there is the 
 * possibility to check the Slater condition for the primal and dual problem before the solution of each SDP by setting a 
 * SCIP parameter, for details see the parameters tab.
 */

/** @page PARAMETERS Additional Parameters
 * The following important parameters (with these default values) were added:
 *
 * <table>
 * <tr><td>relaxing/SDP/freq = 1</td> <td>set this to -1 and lp/solvefreq to 1 to solve LP relaxations with eigenvector cuts</td></tr>
 * <tr><td>constraints/SDP/threads = 1</td> <td>number of threads used for openblas; only available with OpenBLAS/OpenMP and compile options SDPS=sdpa and MKL=true (default for SDPS=sdpa)</td></tr>
 * <tr><td>relaxing/SDP/displaystatistics = TRUE</td> <td>Should statistics about SDP iterations and solver settings/success be printed after quitting SCIP-SDP ?</td></tr>
 * <tr><td>relaxing/SDP/slatercheck = 0</td> <td>Should the Slater condition for the primal and dual problem be checked ahead of solving each SDP? [0: no, 1: yes and output statistics, 2: yes and print warning for every problem not satisfying primal and dual Slater condition]</td></tr>
 * <tr><td>relaxing/SDP/sdpinfo = FALSE</td> <td>Should output of the SDP-Solver be printed to the console?</td></tr>
 * <tr><td>branching/sdpinfobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpinfobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>branching/sdpobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>display/sdpfastsettings/active = 0</td> <td>Should the number of SDP-relaxations solved with the fastest setting (SDPA) or the default formulation (DSDP) be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>display/sdppenalty/active = 0</td> <td>Should the number of SDP-relaxations solved using a penalty formulation be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>display/sdpunsolved/active = 1</td> <td>Should the number of SDP-relaxations that could not be solved be displayed in the console? [0: off, 1: auto, 2:on]</td></tr>
 * <tr><td>heuristics/sdpfracdiving/freq = -1</td> <td>set this to >= 0 to enable a fractional diving heuristic for SDPs</td></tr>
 * <tr><td>heuristics/sdprand/freq = 10</td> <td>set this to -1 to disable the randomized rounding heuristic</td></tr>
 * <tr><td>heuristics/sdprand/generalints = FALSE</td> <td>Should randomized rounding also be applied if there are general integer variables and not only binary variables ?</td></tr>
 * <tr><td>heuristics/sdprand/nrounds = 5</td> <td>number of rounding rounds</td></tr>
 * <tr><td>propagating/obbt-sdp/freq = -1</td> <td>set this to 0 or more to enable SDP-OBBT</td></tr>
 * <tr><td>propagating/obbt/tightcontboundsprobing = FALSE</td> <td>should continuous bounds be tightened during the probing mode?</td></tr>
 * <tr><td>propagating/obbt/tightintboundsprobing = TRUE</td> <td>should integral bounds be tightened during the probing mode?</td></tr>
 * <tr><td>propagating/sdpredcost/freq = 1</td> <td>set this to -1 to disable reduced cost fixing for SDPs</td></tr>
 * <tr><td>propagating/sdpredcost/forbins = TRUE</td> <td>should sdp reduced cost fixing be executed for binary variables?</td></tr>
 * <tr><td>propagating/sdpredcost/forintcons = TRUE</td> <td>should sdp reduced cost fixing be executed for integer and continuous variables?</td></tr>
 * <tr><td>relaxing/SDP/lambdastar = -1</td> <td>the parameter lambda star used by SDPA to set the initial point , set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>relaxing/SDP/maxpenaltyparam = -1</td> <td>the maximum value of the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge, set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>relaxing/SDP/objlimit = FALSE</td> <td>Should an objective limit be given to the SDP-Solver?</td></tr>
 * <tr><td>relaxing/SDP/penaltyparam = -1</td> <td>the starting value of the penalty parameter Gamma used for the penalty formulation if the SDP solver didn't converge, set this to a negative value to compute the parameter depending on the given problem</td></tr>
 * <tr><td>relaxing/SDP/resolve = TRUE</td> <td>Are we allowed to solve the relaxation of a single node multiple times in a row (outside of probing)?</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverepsilon = 1e-04</td> <td>sets the bound for the duality gap in the SDP-Solver</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverfeastol = 1e-06</td> <td>feasibility tolerance for the SDP-Solver (should be less or equal to numerics/feastol)</td></tr>
 * <tr><td>relaxing/SDP/settingsresetfreq = -1</td> <td>frequency for resetting parameters in SDP solver and trying again with fastest settings [-1: never, 0: only at depth settingsresetofs, n: all nodes with depth a multiple of n]</td></tr>
 * <tr><td>relaxing/SDP/settingsresetofs = 0</td> <td>frequency offset for resetting parameters in SDP solver and trying again with fastest settings</td></tr>
 * <tr><td>relaxing/SDP/tightenvb = TRUE</td> <td>Should Big-Ms in varbound-like constraints be tightened before giving them to the SDP-solver ?</td></tr>
 * </table>
 */
