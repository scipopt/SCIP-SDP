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

/**@file   SdpVarMapper.cpp
 * @brief  class that maps scip variables to sdp indices
 * @author Sonja Mars
 */

#include "SdpVarMapper.h"

#include <cassert>                     // for assert

//#include "scip/def.h"                   // for FALSE, SCIP_CALL, SCIP_Real, etc
//#include "scip/pub_message.h"           // for SCIPdebugMessage
//#include "scip/type_prop.h"
//#include "scip/type_lp.h"
//#include "scip/pub_var.h"               // for SCIPvarGetName, etc
#include "scip/scip.h"                  // for SCIPinfinity, SCIPisEQ, etc

/**constructor, initializes private members*/
SdpVarMapper::SdpVarMapper(SCIP* scip): sdp_nvars_(0), allfixed_(FALSE), intsfixed_(FALSE), nfixed_(0), scip_(scip)
{
}


/**init-function of the SsdpVarMapper class,
 *@note must be called by the user for initialise the values of the map and the vectors, maps variable names to the index the variable has for calling sdp, -1 if variable is fixed
 */
SCIP_RETCODE SdpVarMapper::init()
{
   SCIP_VAR** vars;
   int nvars;

   vars = SCIPgetVars(scip_);
   nvars = SCIPgetNVars(scip_);
   scip_var_in_sdp_order_.reserve(nvars);
   int num_fixed_bins_ints = 0;

   //take all the variables and compare their bounds
   for (int i = 0; i < nvars; i++)
   {
      SCIP_Real lb;
      SCIP_Real ub;
      lb = SCIPvarGetLbLocal(vars[i]);
      ub = SCIPvarGetUbLocal(vars[i]);

      if(!SCIPisEQ(scip_, lb, ub))
      {
         //if the variable is not fixed, this variable will be a variable for sdp, so put it in the vector and count _sdp_nvars +1 (zero is also a valid index
         scip_var_in_sdp_order_.push_back(vars[i]);
         name_index_map_[SCIPvarGetName(vars[i])] = sdp_nvars_;
         sdp_nvars_++;


         SCIP_CALL(SCIPcaptureVar(scip_, vars[i]));

      }
      else
      {
         //if the variable is fixed, it won't be given to sdp, so we save -1 als index
         name_index_map_[SCIPvarGetName(vars[i])] = -1;

      }
      if(SCIPisGT(scip_, lb, -SCIPinfinity(scip_)) && SCIPisLT(scip_, ub , SCIPinfinity(scip_)) && SCIPisEQ(scip_, ub, lb) )
      {
         if (SCIPvarIsIntegral(vars[i]))
         {
            num_fixed_bins_ints++;//we count how many of the integer and binary variables are fixed
         }
         nfixed_++;
         //we count how many variables are fixed
      }
   }

   int num_ints_and_bins = SCIPgetNIntVars(scip_) + SCIPgetNBinVars(scip_);

   if (num_fixed_bins_ints == num_ints_and_bins && num_ints_and_bins != 0)
   {
      SCIPdebugMessage("INTS FIXED\n");
      intsfixed_ = TRUE;
   }

   if (nfixed_ == nvars)
   {
      allfixed_ = TRUE;
   }

   assert(sdp_nvars_ == scip_var_in_sdp_order_.size());

   return SCIP_OKAY;
}

/**exit-method, needs to release the variables*/
SCIP_RETCODE SdpVarMapper::exit()
{
   //take all the variables and compare their bounds
   for (unsigned int i = 0; i < scip_var_in_sdp_order_.size(); i++)
   {
      SCIP_CALL(SCIPreleaseVar(scip_, &scip_var_in_sdp_order_[i]));
   }
   return SCIP_OKAY;
}

/**return the maximum index of a variable in sdp*/
int SdpVarMapper::get_sdp_nvars()
{
   return sdp_nvars_;
}

/**return the maximum index of a variable in sdp*/
bool SdpVarMapper::get_allfixed()
{
   return allfixed_;
}
/**return the maximum index of a variable in sdp*/
bool SdpVarMapper::get_intsfixed()
{
   return intsfixed_;
}

/**return the number of fixed variables*/
int SdpVarMapper::get_nfixed()
{
   return nfixed_;
}

/**converts a scip var into a sdp-index*/
int SdpVarMapper::get_sdp_index(SCIP_VAR* var /**<SCIP variable to be converted*/)
{
   return name_index_map_[SCIPvarGetName(var)];
}

/**takes a sdp-index and gets the corresponding scip_var*/
SCIP_VAR* SdpVarMapper::get_scip_var(int idx/**<idx of scip variable*/)
{
   return scip_var_in_sdp_order_[idx];
}


