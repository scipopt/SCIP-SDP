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
 * @version 2.0
 * @author Sonja Mars, Lars Schewe, Tristan Gally, Marc Pfetsch
 * @date 2011-2015
 *
 * The SCIP-SDP-Package is a plug-in for the Software SCIP. It is able to solve Mixed-Integer SDPs using
 * interior-point SDP-solvers (currently DSDP and SDPA) and SCIP as a branch and bound framework. It provides a lot
 * of data handling, some presolve routines and a propagator for reduced cost/dual fixing for semidefinite problems
 * as well as a linear approximation procedure for MISDPs.
 *
 * The solution process of interior-point methods for SDPs is highly dependent on the slater condition. One of the 
 * main purposes of the code is ensuring that the slater condition is not harmed by fixing variables alone. Still the 
 * combination of variable fixings and specific linear or semidefinite constraints might lead to relaxations in which 
 * the slater condition no longer holds. In this case the SDP-solvers may be unable to solve the relaxations or even 
 * return wrong results, which cannot be compensated in all cases. For this purpose there is the possibility to check 
 * the slater condition before the solution of each SDP, for details see the parameters tab.
 */

/** @page PARAMETERS Additional Parameters
 * The following important parameters (with these default values) were added:
 *
 * <table>
 * <tr><td>relaxing/SDP/freq = 1</td> <td>set this to -1 and lp/solvefreq to 1 to solve LP relaxations and add eigenvector cuts</td></tr>
 * <tr><td>propagating/sdpredcost/freq = 1</td> <td>set this to -1 to disable reduced cost fixing for SDPs</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverepsilon = 0.0001</td> <td>sets the bound for the duality gap in the SDP-Solver and also the difference that is allowed when fixing variables before inserting into the SDP-solver</td></tr>
 * <tr><td>relaxing/SDP/sdpsolverfeastol = 0.00001</td> <td>a matrix is considered positive semidefinite if the smallest eigenvalue is bigger than -sdpsolverfeastol</td></tr>
<!-- * <tr><td>relaxing/SDP/threads = 1</td> <td>number of threads to use for the SDP-Solver, currently only supported for SDPA</td></tr> -->
 * <tr><td>relaxing/SDP/sdpinfo = FALSE</td> <td>Should output of the SDP-Solver be printed to the console?</td></tr>
 * <tr><td>relaxing/SDP/objlimit = FALSE</td> <td>Should an objective limit be given to the SDP-Solver?</td></tr>
 * <tr><td>relaxing/SDP/slatercheck = FALSE</td> <td>Should the slater condition for the dual problem be check ahead of solving each SDP?</td></tr>
 * <tr><td>branching/sdpobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * <tr><td>branching/sdpinfobjective/coupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints?</td></tr>
 * <tr><td>branching/sdpinfobjective/singlecoupledvars = FALSE</td> <td>If all branching candidates have objective zero, should we use the sum of the absolute objectives of all continuous variables coupled with the candidate through constraints in which no other candidate appears?</td></tr>
 * </table>
 */
