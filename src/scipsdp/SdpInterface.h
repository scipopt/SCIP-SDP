/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2014 Discrete Optimization, TU Darmstadt               */
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

/**@file   Sdpinterface.cpp
 * @brief  interface to sdp-solver-functions
 * @author Lars Schewe, Sonja Mars
 */

#ifndef SDPINTERFACE_H
#define SDPINTERFACE_H

#include "scip/def.h"                   // for SCIP_Real
#include "scip/type_result.h"           // for SCIP_RESULT
#include "scip/type_retcode.h"          // for SCIP_RETCODE

class SdpProblem;
class SdpVarMapper;

class SdpInterface
{
 public:
   virtual ~SdpInterface() {}

   /**method that puts all the data into the solver object*/
   virtual SCIP_RETCODE put_data_in(
      SdpProblem* problemdata,        /**<problemdata class object*/ 
      SdpVarMapper* varmapper         /**<varmapper class object*/ ) = 0;
   
   /**calls dsdps solving routine and analyses the output*/
   virtual SCIP_RETCODE sdp_solve(
      SdpVarMapper* varmapper,   /**<varmapper class object*/
      char* status,               /**<pointer to store convergenst status of dsdp*/
      int* solutiontype,          /**<pointer to store solutiontype of dsdp*/
      double* sol_for_scip        /**<pointer to store solution in scip format*/) = 0;

   /**also calls dsdps solving routine but after adding a penalty term*/
   virtual SCIP_RETCODE again_with_penalty(
      SCIP_RESULT* result,             /**<pointer to store result*/
      SCIP_Real* lowerbound,           /**<pointer to store lowerbound for scip*/
      SdpVarMapper* varmapper,        /**<varmapper class object*/
      SdpProblem* problemdata
      ) = 0;
};

#endif
