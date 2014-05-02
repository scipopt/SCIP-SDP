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

/**@file   SdpCone.h
 * @brief  Class, where the sdp-data is stored
 * @author Sonja Mars, Lars Schewe
 */

#ifndef SDPCONE_H
#define SDPCONE_H

#include <iterator>

//#include "scip/type_var.h"
//#include "scip/type_retcode.h"
//#include "scip/type_scip.h"
//#include "scip/type_sol.h"
#include "scip/scip.h"

class SdpCone
{
public:

   // we assume that all memory is allocated via SCIP and we take control over all arrays
   SdpCone(SCIP* scip,
      int blocksize,
      SCIP_VAR** vars,
      int* col,
      int* row,
      double* vals,
      int nnz,
      int* const_col,
      int* const_row,
      double* const_vals,
      int const_nnz);

   SdpCone(const SdpCone& sdpcone);

   void swap(SdpCone& other);

   SdpCone& operator=(SdpCone other);

   ~SdpCone();

   /**returns the blocksize of this sdpcone*/
   int get_blocksize() const;

   /**returns the maximal value of the coefficients of the rhs*/
   double get_max_rhs() const;

   /**writes the constraints matrix vor variable vidx in matrix*/
   SCIP_RETCODE get_constraint_matrix(double* matrix, int vidx) const;

   /**
    * computes a version of the constraint matrix #vidx where all zero rows/cols are removed
    *
    * @param[output] matrix preallocated array containing the entries of the constraint matrix in col-major order in entries [0, (new_size^2 -1)]
    * @param[output] new_size contains size of shrunk constraint matrix on output
    * @param[input] vidx constraint matrix to return
    */
   SCIP_RETCODE get_shrunk_constraint_matrix(double* matrix, int* new_size, int vidx) const;

   /**writes a product of constraint matrix with vector vector into out*/
   SCIP_RETCODE form_inner_prod_with_constraint_matrix(int vidx, double* vector, double* out) const;

   /**calculates sum_i A_i x_i for a given solution x_i*/
   SCIP_RETCODE assemble_matrix_from_solution(double* matrix, SCIP_SOL* sol) const;

   /**writes A_0 into matrix*/
   SCIP_RETCODE get_constant_matrix(double* matrix) const;

   /**updated sdpcone-data after presolve and handles negated, fixed, aggregated and deleted vars*/
   SCIP_RETCODE fix_vars();

   /**transforms variables for SCIP to have an array of tranformed variables in the sdpcone*/
   SCIP_RETCODE transform_vars();

   /**returns the number of variables that belong to this sdpcone*/
   int get_nvars() const;

   /**returns the Scip_var for the idx-th index*/
   SCIP_VAR* get_var(int idx) const;

   /**returns the number of nonzeros in this sdpcone*/
   int get_nnz() const;

   /**returns the number of nonzeros in the constant matrix*/
   int get_const_nnz() const;

   /**returns the number of fixed variables in the sdpcone*/
   int get_sdpcone_nfixed();

   /**takes the scip_vars and writes their postion into the array n_fixed_vars*/
   SCIP_RETCODE transform_vars_to_pos(SCIP* scip, SCIP_VAR** fixed_vars, int* n_fixed_vars);

   /**
    * representation of cone elements for outside view
    */
   struct element
   {
      int row;  /// row index
      int col;  /// col index
      int vidx; /// variable index
      int eidx; /// element index in packed format
      double val; /// value
   };

   /**
    * Iterator class for elements of the left hand side
    */
   class LhsIterator : public std::iterator<std::input_iterator_tag,
                                            const SdpCone::element>
   {
   public:

      /**
       * Constructor for valid lhs iterator
       */
      LhsIterator(SdpCone* c, int n_fixed_vars);

      /**
       * constructs an iterator which is already at end
       */
      LhsIterator();

      /**
       * standard iterator operations
       */
      LhsIterator& operator++();
      LhsIterator operator++(int);
      const SdpCone::element& operator*() const;
      const SdpCone::element* operator->() const;
      bool operator==(const LhsIterator& other);
      bool operator!=(const LhsIterator& other);

   private:
      SdpCone* sdpcone_; /// back pointer to owning SdpCone
      int vpos_; /// position in the variable array
      int epos_; /// position in the col/row array
      int n_fixed_vars_;
      SdpCone::element curr_element_; /// current position
      bool end_; /// is at end?
   };

   /**iterator for going through the lhs and handling fixed variables*/
   LhsIterator lhs_begin(SCIP_VAR** fixed_vars, int n_fixed_vars);
   LhsIterator lhs_end();

   /**
    * Iterator class for elements of the right hand side
    */
   class RhsIterator : public std::iterator<std::input_iterator_tag,
                                            const SdpCone::element>
   {
   public:
      /**
       * Constructor for valid rhs iterator (including fixings)
       */
      RhsIterator(SdpCone* c, int n_fixed_vars, double* fixed_values);

      /**
       * constructs an iterator which is already at end
       */
      RhsIterator();

      ~RhsIterator();

      /**
       * standard iterator operations
       */
      RhsIterator& operator++();
      RhsIterator operator++(int);
      const SdpCone::element& operator*() const;
      const SdpCone::element* operator->() const;
      bool operator==(const RhsIterator& other);
      bool operator!=(const RhsIterator& other);

   private:
      SdpCone* sdpcone_; /// back pointer to owning SdpCone
      int pos_;
      int* epos_; /// position in the col/row array
      SdpCone::element curr_element_;  /// current element

      int n_fixed_vars_; /// number of fixings
      double* fixed_values_; /// value of fixings
      bool end_; /// is at end?
   };


   /**iterator for going through the rhs and handling fixed variables*/
   RhsIterator rhs_begin(SCIP_VAR** fixed_vars, int n_fixed_vars, double* fixed_values);
   RhsIterator rhs_end();


private:
   /**
    * @internal helper function to convert input into CSR
    *
    * @param vars variable array
    * @param col column array
    * @param row row array
    * @param vals values array
    */
   void compress_representation(
      SCIP_VAR** vars,
      int* col,
      int* row,
      double* vals);

   void sort_const(int* const_col, int* const_row, double* const_vals);


   SCIP* scip_;
   bool captured_;      /**<internal helper variable for memory management*/

   int blocksize_;      /**<number of rows or number of cols*/

   int nvars_;          /**<number of unique vars in constraint*/
   SCIP_VAR** uvars_;   /**<array of unique vars in constraint*/
   int* vbeg_;          /**<array indicating where matrix i begins (with sentinel)*/

   int * col_;          /**<columns, where the nonzero entries are located*/
   int * row_;          /**<rows, where the nonzero entries are located*/

   double* vals_;       /**<nonzero entries in the matrix*/
   int nnz_;            /**<number of nonzero entries in the block*/

   int * const_col_;    /**<columns, where the nonzero entries of the rhs are located*/
   int * const_row_;    /**<rows, where the nonzero entries of the rhs are located*/
   double* const_vals_; /**<nonzero entries of the constant matrix*/
   int const_nnz_;      /**<numer of nonzero entries in the constant matrix*/

   int* pos_fixed_vars_;
   int sdpcone_nfixed_ ;
};

#endif
