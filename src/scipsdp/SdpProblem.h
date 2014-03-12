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

/**@file   SdpProblem.h
 * @brief  Class, where the sdp-data and lp-data is stored
 * @author Sonja Mars
 */

#ifndef SDPPROBLEM_H
#define SDPPROBLEM_H

#include <vector>                       // for vector
//#include "scip/def.h"                   // for SCIP_Real
//#include "scip/type_lp.h"               // for SCIP_ROW
//#include "scip/type_retcode.h"          // for SCIP_RETCODE
//#include "scip/type_scip.h"             // for SCIP
#include "scip/scip.h"

class SdpCone;  // lines 10-10
class SdpVarMapper;  // lines 11-11

class SdpProblem
{
 public:
   /**Constructor*/
   SdpProblem(SCIP* scip, SdpVarMapper* varmapper);

   /**Destructor*/
   ~SdpProblem();

   /** gets the ith sdpcone*/
   SdpCone* get_sdpcone(int i) const;

   /** gets the number of sdpcones*/
   int get_nsdpcones() const ;

   /**gets the vector for_vals*/
   const double* get_for_vals() const;

   /**gets the vector for_constraints*/
   const int* get_for_constraint() const;

   /** gets the vector for_matind*/
   const int* get_for_matind() const;

   /**gets the size of the vector for_vals*/
   int get_for_vals_size() const;

   /**gets the size of the vector for_constraints*/
   int get_for_constraint_size() const;

   /**gets the size of the vector for_matind*/
   int get_for_matind_size() const;

   /**gets the size of the lp-block*/
   int get_size_lpblock() const;

   /**gets the number of nonzeros in the lp-block*/
   int get_lp_nnz() const;

   /**Method for adding linear constraints to the structure we need for dsdp
   *the arrays for_* will later be added to dsdp as lp cone*/
   SCIP_RETCODE addconstraint(
   SdpVarMapper* varmapper,          /**<varmapper class data*/
   int *position,                    /**<diagonal-position where constraint should be added*/
   double rhs,                       /**<rhs of constraint*/
   int* lininds,                     /**<indices of variables in constraint*/
   int nlininds,                     /**<number of variables in constraint*/
   SCIP_Real* vals,                  /**<coefficients of variables in constraint*/
   SCIP_Real* ubs                    /**upper bound of variable<*/);

   /** adds bound for variables */
   SCIP_RETCODE addbound(
      int*               position,           /**< position where bound constraint is added */
      double             rhs,                /**< rhs of variable */
      int                lininds,            /**< variable index */
      double             linvals,            /**< coefficient of variable */
      SdpVarMapper*      varmapper           /**< varmapper class object */
      );

   /**adding lp data out of scip to the vector ding->for_matind, ding->for_block, ding->for_constraint, ding->for_vals*/
   SCIP_RETCODE get_rows_data(
   SdpVarMapper* varmapper,          /**<varmapper class data*/
   SCIP_ROW** rows,                  /**<rows to add to lp-block*/
   int nrows,                        /**<number of rows*/
   int* position                     /**<pointer to store position of diagonal entries of lp block*/);

 private:
   SCIP* scip_;                       /**scip-data structure*/
   SdpCone**  sdpcones_;             /**<sdpcones*/
   int        nsdpcones_;            /**<number of sdpcones*/
   std::vector <double> for_vals_;   /**<vector to store nonzeros*/
   std::vector <int> for_constraint_;/**<vector to store number of variables of nonzero*/
   std::vector <int> for_matind_;    /**<vector to store position in array where nonzero is saved*/
   int size_lpblock_;

};

#endif
