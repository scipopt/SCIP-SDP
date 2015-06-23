/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Tristan Gally
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 *
 * @version 2.0.0
 * @author Tristan Gally, Marc Pfetsch; Sonja Mars, Lars Schewe
 * @date 2011-2015
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
 * Please note that the interface to SDPA is still in beta state. It works well for some instances and is faster than DSDP
 * for those, but currently fails for others because of numerical problems, as some parameters need further tuning.
 *
 * The solution process of interior-point methods for SDPs is highly dependent on the Slater condition. One of the main
 * purposes of the code is ensuring that the slater condition is not harmed by fixing variables in the branch-and-bound
 * process. However in some cases the combination of variable fixings and specific linear or semidefinite constraints might
 * still lead to relaxations for which the Slater condition no longer holds. In this case the SDP-solvers may be unable to
 * solve the relaxations or even return wrong results, which cannot be compensated. For this purpose there is the 
 * possibility to check the Slater condition for the dual problem (which still does not guarantee it for the primal)
 * before the solution of each SDP by setting a SCIP parameter, for details see the parameters tab.
 */

/** @page PARAMETERS Additional Parameters
 * The following important parameters (with these default values) were added:
 *
 * <table>
 * <tr><td>relaxing/SDP/freq = 1</td> <td>set this to -1 and lp/solvefreq to 1 to solve LP relaxations with eigenvector cuts</td></tr>
 * <tr><td>propagating/sdpredcost/freq = 1</td> <td>set this to -1 to disable reduced cost fixing for SDPs</td></tr>
 * <tr><td>propagating/sdpredcost/advanced/forbins = TRUE</td> <td>should sdp reduced cost fixing be executed for binary variables?</tr>
 * <tr><td>propagating/sdpredcost/advanced/forintcons = TRUE</td> <td>should sdp reduced cost fixing be executed for integer and continuous variables?</tr>
 * <tr><td>relaxing/SDP/sdpsolverepsilon = 0.0001</td> <td>sets the bound for the duality gap in the SDP-Solver</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverfeastol = 0.00001</td> <td>feasibility tolerance for the SDP-Solver (should be less or equal to numerics/feastol)</td></tr>
 * <tr><td>relaxing/SDP/sdpinfo = FALSE</td> <td>should output of the SDP-Solver be printed to the console?</td></tr>
 * <tr><td>relaxing/SDP/objlimit = FALSE</td> <td>should an objective limit be given to the SDP-Solver?</td></tr>
 * <tr><td>relaxing/SDP/slatercheck = FALSE</td> <td>should the Slater condition for the dual problem be checked ahead of solving each SDP?</td></tr>
 * <tr><td>branching/sdpobjective/coupledvars = FALSE</td> <td>if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpobjective/singlecoupledvars = FALSE</td> <td>if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>branching/sdpinfobjective/coupledvars = FALSE</td> <td>if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpinfobjective/singlecoupledvars = FALSE</td> <td>if all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * </table>
 */
