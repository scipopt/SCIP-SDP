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

/**@file   SdpVarMapper.h
 * @brief  class that maps scip variables to sdp indices
 * @author Sonja Mars
 */

#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__

#include <vector>
#include <map>
#include <string>
#include "scip/type_scip.h"
#include "scip/type_retcode.h"
#include "scip/type_var.h"

/**class for mapping scip variables and sdp indices*/
class SdpVarMapper
{
 public:
   /**constructor*/
   SdpVarMapper(
   SCIP* scip);

   /**init method that creates all the vectors and number and captures variables*/
   SCIP_RETCODE init();

   /**frees memory and releases variables*/
   SCIP_RETCODE exit();

   /**return the indes of sdp of var*/
   int get_sdp_index(
   SCIP_VAR* var/**<scip variable for which index is calculated for*/);

   /**returns the scip variable for which index idx stands for*/
   SCIP_VAR* get_scip_var(
      int idx/**<index of the varialbe which should be returned*/);

   /**return the number of non fixed variable, given to sdp*/
   int get_sdp_nvars();

   /**return wheather all variables are fixed*/
   bool get_allfixed();

   /**return wheather all ints are fixed*/
   bool get_intsfixed();

   /**return number of fixed variables*/
   int get_nfixed();


 private:

   unsigned int sdp_nvars_;                        /**<number of variables in sdp*/
   bool allfixed_;                                 /**<bool for indicating if all variables are fixed*/
   bool intsfixed_;                                /**<bool for indicating if all integer and binary variables are fixed*/
   int nfixed_;                                    /**<number of currently fixed variables*/

   SCIP* scip_;                                    /**<SCIP data structure*/
   std::vector<SCIP_VAR*> scip_var_in_sdp_order_;  /**<vector with non-fixed scip varialbes in the order they are given to sdp*/
   std::map<std::string, int> name_index_map_;     /**<maps a variable name to an index, all variables are in here, if a variable is fixed, it gets indes -1*/

};



#endif
