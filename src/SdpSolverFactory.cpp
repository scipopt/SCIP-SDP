/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2012 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
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

/**@file   SdpSolverFactory.cpp
 * @brief  turns solvers on and off
 * @author Lars Schewe, Sonja Mars
 */

#include "SdpSolverFactory.h"

#ifdef USE_DSDP
#include "DsdpInterface.h"
#endif

SdpInterface* SdpSolverFactory::createSdpSolver(SCIP* scip, const std::string& name)
{
#ifdef USE_DSDP
   if (name == "dsdp")
   {
      return (new DsdpInterface(scip));
   } 
#endif

   assert(false); // TODO(LS) error handling
}

SdpInterface* SdpSolverFactory::createSdpSolver(SCIP* scip)
{
#ifdef USE_DSDP
   return (new DsdpInterface(scip));
#endif

   assert(false); // TODO(LS) error handling
}
