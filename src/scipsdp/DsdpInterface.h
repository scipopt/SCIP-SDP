/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   DsdpInterface.h
 * @brief  interface to dsdp solver
 * @author Sonja Mars
 */

#ifndef _DsdpInterface_h
#define _DsdpInterface_h

#include "SdpInterface.h"

#include "dsdpbasictypes.h"             // for DSDP, DSDP_C
//#include "scip/def.h"                   // for SCIP_Real
//#include "scip/type_result.h"           // for SCIP_RESULT
//#include "scip/type_retcode.h"          // for SCIP_RETCODE
//#include "scip/type_scip.h"             // for SCIP
#include "scip/scip.h"

class SdpProblem;
class SdpVarMapper;

/**class for communicating with dsdp*/
class DsdpInterface : public SdpInterface
{
 public:
   /**constructor*/
   DsdpInterface(SCIP* scip);

   /**destructor*/
   virtual ~DsdpInterface();

   /**method that puts all the data into the solver object*/
   virtual SCIP_RETCODE put_data_in(
      SdpProblem* problemdata,
      SdpVarMapper* varmapper         /**<varmapper class object*/ );

   /**calls dsdps solving routine and analyses the output*/
   virtual SCIP_RETCODE sdp_solve(
      SdpVarMapper* varmapper,   /**<varmapper class object*/
      char* status,               /**<pointer to store convergenst status of dsdp*/
      int* solutiontype,          /**<pointer to store solutiontype of dsdp*/
      double* sol_for_scip        /**<pointer to store solution in scip format*/);

   /**also calls dsdps solving routine but after adding a penalty term*/
   virtual SCIP_RETCODE again_with_penalty(
      SCIP_RESULT* result,             /**<pointer to store result*/
      SCIP_Real* lowerbound,           /**<pointer to store lowerbound for scip*/
      SdpVarMapper* varmapper,        /**<varmapper class object*/
      SdpProblem* problemdata
      );


 private:
   SCIP*       scip_;      /**<SCIP data structure*/
   //for DSDP
   DSDP        dsdp_;      /**<DSDP-object*/
   int*        matind_;    /**<array needed for dsdp, the place (calculated using cols and rows) where the nnz are written is saved here*/
   double*     nnz_;       /**<nonzeros in a specific sdpcone*/
   int*        ittt_;      /**<array needed for dsdp*/
};

#endif
