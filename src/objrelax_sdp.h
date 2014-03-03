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

/**@file   objrelax_sdp.h
 * @brief  sdp-relaxator
 * @author Sonja Mars
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJRELAXSDP_H__
#define __SCIP_OBJRELAXSDP_H__

#include "objscip/objrelax.h"           // for ObjRelax
#include "scip/def.h"                   // for SCIP_Real
#include "scip/type_relax.h"            // for SCIP_RELAX
#include "scip/type_result.h"           // for SCIP_RESULT
#include "scip/type_retcode.h"          // for SCIP_RETCODE
#include "scip/type_scip.h"             // for SCIP

namespace scip
{
   /**class of the objective sdp-relaxator*/
   class ObjRelaxSdp : public ObjRelax
   {

   public:
      /*lint --e{1540}*/

      /** default constructor */
      ObjRelaxSdp(SCIP * scip) ;

      /** destructor */
      virtual ~ObjRelaxSdp() {}



      /** execution method of relaxator
       *
       *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
       *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
       *
       *  possible return values for *result (if more than one applies, the first in the list should be used):
       *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
       *  - SCIP_CONSADDED  : an additional constraint was generated, and the relaxator should not be called again on the
       *                      same relaxation
       *  - SCIP_REDUCEDDOM : a variable's domain was reduced, and the relaxator should not be called again on the same
       *                      relaxation
       *  - SCIP_SEPARATED  : a cutting plane was generated, and the relaxator should not be called again on the same
       *                      relaxation
       *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
       *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
       *                      planes); however, it is able to continue the solving in order to improve the dual bound
       *  - SCIP_DIDNOTRUN  : the relaxator was skipped
       */
      virtual SCIP_RETCODE scip_exec(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_RELAX*        relax,              /**< the relaxator itself */
         SCIP_Real*         lowerbound,         /**< pointer to store a lowerbound for the current node */
         SCIP_RESULT*       result              /**< pointer to store the result of the relaxation call */
         );

  
   };

} /* namespace scip */



#endif
