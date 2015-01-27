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

/**@file   nodesel_prio.h
 * @ingroup NODESELECTORS
 * @brief  nodeselector that chooses candidate with highest priority (dual bound is used as tiebreaker)
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_NODESEL_PRIO_H__
#define __SCIP_NODESEL_PRIO_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the priority node selector and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeNodeselPrio(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** inform the nodeselector about a new node being added by the branching rule and gives it the corresponding priority value */
EXTERN
SCIP_RETCODE SCIPnodeselPrioInsertNodePrio(
   SCIP*                  scip,              /**< scip-pointer */
   SCIP_NODESEL*          nodesel,           /**< pointer to the nodeselector */
   SCIP_NODE*             node,              /**< node for which priority is given */
   SCIP_Real              priority           /**< nodeselection-priority of the given node */
   );

#ifdef __cplusplus
}
#endif

#endif
