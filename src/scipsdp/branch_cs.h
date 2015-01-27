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

/**@file   branch_cs.h
 * @ingroup BRANCHINGRULES
 * @brief  special branching rule for compressed sensing problems in SCIPSDP
 * @author Tristan Gally
 *
 * This branching rule is used for a MISDP to compute exact left- and right-hand-side-restricted-isometry-constants. The idea is, that if we branch on a binary
 * variable, the objective of the problem with this binary variable set to one will most likely change only minimally, while we can compute a good approximation
 * on the solution (and therefore the objective) of the problem where it is set to zero. In this case, because of the Trace-Constraint, if we set a binary
 * variable to zero, thereby eliminating it, the whole column and row will be set to zero (all variables coupled with it), but because this decreases the trace
 * by the value of the corresponding diagonal entry (the only continuous variable, that is only coupled with this binary variable), we have to multiply all
 * other entries of the matrix by 1-(1/y), where y is the corresponding diagonal entry.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_CS_H__
#define __SCIP_BRANCH_CS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the compressed sensing branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleCs(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
