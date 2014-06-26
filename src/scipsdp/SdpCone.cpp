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

/**@file   SdpCone.cpp
 * @brief  Class, where the sdp-data is stored
 * @author Sonja Mars, Lars Schewe, Tristan Gally
 */

#define SCIP_DEBUG

#include "SdpCone.h"

#include <cassert>                      // for assert
#include <algorithm>                    // for swap, max, sort, unique

#include "scip/scip.h"                  // for SCIPallocMemoryArray, etc


/** compare variables with less than*/
static
bool varLT (
   SCIP_VAR *var1,  /**<SCIP_Var one */
   SCIP_VAR *var2   /**<SCIP_VAR two*/)
{
   return (SCIPvarCompare(var1, var2) == -1);
}

/** compare if two variables are equal*/
static
bool varEQ (
   SCIP_VAR *var1,  /**<SCIP_Var one */
   SCIP_VAR *var2   /**<SCIP_Var two */)
{
   return (SCIPvarCompare(var1, var2) == 0);
}


SdpCone::SdpCone(SCIP* scip,
   int blocksize,
   SCIP_VAR** vars,
   int* col,
   int* row,
   double* vals,
   int nnz,
   int* const_col,
   int* const_row,
   double* const_vals,
   int const_nnz) :
   scip_(scip),
   captured_(TRUE),
   blocksize_(blocksize),
   nvars_(0),
   uvars_(NULL),
   vbeg_(NULL),
   col_(NULL),
   row_(NULL),
   vals_(NULL),
   nnz_(nnz),
   const_col_(NULL),
   const_row_(NULL),
   const_vals_(NULL),
   const_nnz_(const_nnz),
   pos_fixed_vars_(NULL),
   sdpcone_nfixed_(0)
{
   SCIP_VAR** variables;

   SCIPduplicateBufferArray(scip_, &variables, vars, nnz_);

   std::sort(variables, variables + nnz_, varLT);

   SCIP_VAR** new_end = std::unique(variables, variables + nnz_, varEQ);

   nvars_ = new_end - variables;

   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip_, &uvars_, nvars_));
   for (int i = 0; i < nvars_; ++i)
   {
      uvars_[i] = variables[i];
   }

   for (int i = 0; i < nvars_; ++i)
   {
      SCIP_CALL_ABORT(SCIPcaptureVar(scip_, uvars_[i]));
   }

   SCIPfreeBufferArray(scip_, &variables);

   /// copy constant matrix into the cone

   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip, &const_col_, const_nnz_));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip, &const_row_,  const_nnz_));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip, &const_vals_, const_nnz_));


   sort_const(const_col, const_row, const_vals);

   compress_representation(vars, col, row, vals);

   /// free all input arrays

   SCIPfreeBlockMemoryArray(scip_, &vars, nnz_);
   SCIPfreeBlockMemoryArray(scip_, &row, nnz_);
   SCIPfreeBlockMemoryArray(scip_, &col, nnz_);
   SCIPfreeBlockMemoryArray(scip_, &vals, nnz_);
   SCIPfreeBlockMemoryArray(scip, &const_col, const_nnz_);
   SCIPfreeBlockMemoryArray(scip, &const_row, const_nnz_);
   SCIPfreeBlockMemoryArray(scip, &const_vals, const_nnz_);

}

SdpCone::~SdpCone()
{
   if (captured_)
   {
      for (int i = 0; i < nvars_; ++i)
      {
         SCIP_CALL_ABORT(SCIPreleaseVar(scip_, &uvars_[i]));
      }
   }
   if (pos_fixed_vars_)
   {
      SCIPfreeMemoryArray(scip_, &pos_fixed_vars_);
   }

   SCIPfreeMemoryArray(scip_, &uvars_);
   SCIPfreeMemoryArray(scip_, &vbeg_);
   SCIPfreeMemoryArray(scip_, &row_);
   SCIPfreeMemoryArray(scip_, &col_);
   SCIPfreeMemoryArray(scip_, &vals_);
   SCIPfreeMemoryArray(scip_, &const_vals_);
   SCIPfreeMemoryArray(scip_, &const_row_);
   SCIPfreeMemoryArray(scip_, &const_col_);
}

SdpCone::SdpCone(const SdpCone& sdpcone) :
   scip_(sdpcone.scip_),
   captured_(FALSE),
   blocksize_(sdpcone.blocksize_),
   nvars_(sdpcone.nvars_),
   uvars_(NULL),
   vbeg_(NULL),
   col_(NULL),
   row_(NULL),
   vals_(NULL),
   nnz_(sdpcone.nnz_),
   const_col_(NULL),
   const_row_(NULL),
   const_vals_(NULL),
   const_nnz_(sdpcone.const_nnz_),
   pos_fixed_vars_(NULL)
{
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &uvars_, sdpcone.uvars_, nvars_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &vbeg_, sdpcone.vbeg_, nvars_ + 1));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &row_, sdpcone.row_, nnz_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &col_, sdpcone.col_, nnz_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &vals_, sdpcone.vals_, nnz_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &const_vals_, sdpcone.const_vals_, const_nnz_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &const_row_, sdpcone.const_row_, const_nnz_));
   SCIP_CALL_ABORT(SCIPduplicateMemoryArray(scip_, &const_col_, sdpcone.const_col_, const_nnz_));

   SCIP_STAGE stage;
   stage = SCIPgetStage(scip_);
   if (stage != SCIP_STAGE_FREETRANS)
   {
      captured_ = TRUE;
      for (int i = 0; i < nvars_; ++i)
      {
         SCIP_CALL_ABORT(SCIPcaptureVar(scip_, uvars_[i]));
      }
   }
}

void SdpCone::swap(SdpCone& other)
{
   std::swap(scip_, other.scip_);
   std::swap(captured_, other.captured_);
   std::swap(blocksize_, other.blocksize_);
   std::swap(nvars_, other.nvars_);
   std::swap(uvars_, other.uvars_);
   std::swap(vbeg_, other.vbeg_);
   std::swap(col_, other.col_);
   std::swap(row_, other.row_);
   std::swap(vals_, other.vals_);
   std::swap(nnz_, other.nnz_);
   std::swap(const_col_, other.const_col_);
   std::swap(const_row_, other.const_row_);
   std::swap(const_vals_, other.const_vals_);
   std::swap(const_nnz_, other.const_nnz_);
}

SdpCone& SdpCone::operator=(SdpCone other)
{
   swap(other);

   return *this;
}

int SdpCone::get_nvars() const
{
   return nvars_;
}

SCIP_VAR* SdpCone::get_var(int idx) const
{
   assert(0 <= idx);
   assert(idx < nvars_);

   return uvars_[idx];
}

int SdpCone::get_blocksize() const
{
   return blocksize_;
}

double SdpCone::get_max_rhs() const
{
   double max_rhs;

   max_rhs = 0.0;

   for (int j = 0; j < const_nnz_; ++j)
   {
      max_rhs = std::max(max_rhs, fabs(const_vals_[j]));
   }
   return max_rhs;
}


SCIP_RETCODE SdpCone::form_inner_prod_with_constraint_matrix(int vidx, double* vector, double* out) const
{
   *out = 0.0;

   for (int i = vbeg_[vidx]; i < vbeg_[vidx + 1]; ++i)
   {

      *out += vals_[i] * vector[row_[i] - 1] * vector[col_[i] - 1];

      if (( col_[i] - 1) != (row_[i] - 1))
      {
         *out += vals_[i] * vector[col_[i] - 1] * vector[row_[i] - 1];
      }
   }

   return SCIP_OKAY;
}

SCIP_RETCODE SdpCone::get_constraint_matrix(double* matrix,
   int vidx) const
{
   //matrix must already have the right dimensions
   assert(matrix != NULL);
   for (int j = 0; j < (blocksize_ * blocksize_); ++j)
   {
      matrix[j] = 0.0;
   }

   for (int i = vbeg_[vidx]; i < vbeg_[vidx + 1]; ++i)
   {
      matrix[ (col_[i] - 1) * blocksize_ + (row_[i] - 1) ] += vals_[i] ;

      if (( col_[i] - 1) != (row_[i] - 1))
      {
         matrix[ (row_[i] - 1) * blocksize_ + (col_[i] - 1) ] += vals_[i];
      }
   }
   return SCIP_OKAY;
}

SCIP_RETCODE SdpCone::assemble_matrix_from_solution(double* matrix,
   SCIP_SOL* sol) const

{
   assert(matrix != NULL);
   for (int j = 0; j < (blocksize_ * blocksize_); ++j)
   {
      matrix[j] = 0.0;
   }

   for (int k = 0; k < nvars_; ++k)
   {
      double solval = SCIPgetSolVal(scip_, sol, uvars_[k]);

      for (int i = vbeg_[k]; i < vbeg_[k + 1]; ++i)
      {
         matrix[ (col_[i] - 1) * blocksize_ + (row_[i] - 1) ] += vals_[i] *  solval;

         if (( col_[i] - 1) != (row_[i] - 1))
         {
            matrix[ (row_[i] - 1) * blocksize_ + (col_[i] - 1) ] += vals_[i] *  solval;
         }
      }
   }

   ///@attention we subtract the constant part from the result matrix
   for (int i = 0; i < const_nnz_; ++i)
   {
      matrix[ (const_col_[i] - 1) * blocksize_ + (const_row_[i] - 1) ] -= const_vals_[i];

      if (( const_col_[i] - 1) != (const_row_[i] - 1))
      {
         matrix[ (const_row_[i] - 1) * blocksize_ + (const_col_[i] - 1) ] -= const_vals_[i];
      }
   }

   return SCIP_OKAY;
}


SCIP_RETCODE SdpCone::get_constant_matrix(double* matrix) const
{
   //matrix must already have the right dimensions
   assert(matrix != NULL);
   for (int j = 0; j < (blocksize_ * blocksize_); ++j)
   {
      matrix[j] = 0.0;
   }

   for (int k = 0; k < (const_nnz_); ++k)
   {
      matrix[ (const_col_[k] - 1) * blocksize_ + (const_row_[k] - 1) ] += const_vals_[k] ;

      if (( const_col_[k] - 1) != (const_row_[k] - 1))
      {
         matrix[ (const_row_[k] - 1) * blocksize_ + (const_col_[k] - 1) ] += const_vals_[k];
      }
   }
   return SCIP_OKAY;
}

SCIP_RETCODE SdpCone::fix_vars()
{
   int how_many_deleted = 0;
   int deleted_nz = 0;
   for (int j = 0; j < nvars_; ++j)
   {
      if ( (SCIPvarGetStatus(SCIPvarGetProbvar(uvars_[j])) == SCIP_VARSTATUS_FIXED) || (SCIPvarGetStatus(SCIPvarGetProbvar(uvars_[j])) == SCIP_VARSTATUS_AGGREGATED) || SCIPisEQ(scip_, SCIPvarGetLbLocal(SCIPvarGetProbvar(uvars_[j])),SCIPvarGetUbLocal(SCIPvarGetProbvar(uvars_[j]))) )
      {
         how_many_deleted++; //number of deleted and aggregated vars
         deleted_nz = deleted_nz + vbeg_[j + 1] - vbeg_[j];
#ifdef SCIP_DEBUG
         SCIPdebugMessage("variable %s has been fixed to value %f with varstatus %u \n", SCIPvarGetName(SCIPvarGetProbvar(uvars_[j])), SCIPvarGetLbLocal(SCIPvarGetProbvar(uvars_[j])), SCIPvarGetStatus(SCIPvarGetProbvar(uvars_[j])));
#endif
      }
   }

   int new_nnz = nnz_;

   int* new_const_col;
   int* new_const_row;
   double* new_const_vals;
   int new_const_nnz;

   new_const_nnz = const_nnz_;

   SCIP_CALL( SCIPallocBufferArray(scip_, &new_const_col, const_nnz_ + deleted_nz));
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_const_row, const_nnz_ + deleted_nz));
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_const_vals, const_nnz_ + deleted_nz));

   for (int i = 0; i < const_nnz_; ++i)
   {
      new_const_col[i] = const_col_[i];
      new_const_row[i] = const_row_[i];
      new_const_vals[i] = const_vals_[i];
   }

   int *aggr_row;
   int *aggr_col;
   double *aggr_vals;
   SCIP_VAR **aggr_vars;
   int* aggr_vbeg;


   SCIP_CALL( SCIPallocBufferArray(scip_, &aggr_row, new_nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &aggr_col, new_nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &aggr_vals, new_nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &aggr_vbeg, nvars_) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &aggr_vars, nvars_) );

   aggr_vbeg[0] = 0;

   int count_aggr = 0;
   int count_nnz_aggr = 0;
   for (int i = 0; i < nvars_; ++i)
   {

      SCIP_VAR* temp_prob_var = SCIPvarGetProbvar (uvars_[i]);

      if (SCIPvarGetProbindex(temp_prob_var) == -1 || SCIPisEQ(scip_, SCIPvarGetLbLocal(temp_prob_var), SCIPvarGetUbLocal(temp_prob_var)))
      {  // TODO: latter part seems to cause problems
         SCIP_VARSTATUS status;
         status = SCIPvarGetStatus(uvars_[i]);

         if (status == SCIP_VARSTATUS_AGGREGATED)
         {
            count_aggr++;
         }
         int save_position = -1;

         assert( SCIPisEQ(scip_, SCIPvarGetLbLocal(temp_prob_var), SCIPvarGetUbLocal(temp_prob_var)) );

         SCIP_Real constant = 0;
         SCIP_Real scalar = 1;
         SCIP_VAR* var;
         var = uvars_[i];

         for (int k = vbeg_[i]; k < vbeg_[i + 1]; ++k)
         {
            double val = 0;
            if (status == SCIP_VARSTATUS_NEGATED)
            {
               val = -1;
            }
            if (status == SCIP_VARSTATUS_FIXED )
            {
               val = -vals_[k] * SCIPvarGetLbLocal(temp_prob_var);
            }

            if (SCIPisEQ(scip_, SCIPvarGetLbLocal(temp_prob_var), SCIPvarGetUbLocal(temp_prob_var)) && status == SCIP_VARSTATUS_COLUMN)
            {  //no need to tell SCIP that something needs to be fixed, because at this point SCIP would only change bounds accordingly, which it already has
               val = -vals_[k] * SCIPvarGetLbLocal(temp_prob_var);
            }

            if (status == SCIP_VARSTATUS_AGGREGATED)
            {
               SCIP_CALL(SCIPgetProbvarSum(scip_, &var, &scalar, &constant));
               val = constant;
            }

            for (int j = 0; j < new_const_nnz; ++j)
            {

               if (new_const_row[j] == row_[k] && new_const_col[j] == col_[k])
               {
                  //there is already an entry at this position on the rhs
                  save_position = j;
               }
            }
            if (save_position == -1 && !SCIPisEQ(scip_, val, 0))
            {
               new_const_col[new_const_nnz] = col_[k];
               new_const_row[new_const_nnz] = row_[k];
               new_const_vals[new_const_nnz] = val;
               new_const_nnz++;
            } else if (save_position > -1)
            {
               new_const_vals[save_position] += val;
            }

            if (status == SCIP_VARSTATUS_NEGATED)
            {
               vals_[k] = -vals_[k];
            }
            if (status == SCIP_VARSTATUS_AGGREGATED)
            {
               aggr_row[count_nnz_aggr] = row_[k];
               aggr_col[count_nnz_aggr] = col_[k];
               aggr_vals[count_nnz_aggr] = vals_[k] / scalar;

               count_nnz_aggr++;
            }
         }
         aggr_vars[i] = var;
         aggr_vbeg[i + 1] = count_nnz_aggr;
      }
   }

   SCIPfreeMemoryArray(scip_, &const_col_);
   SCIPfreeMemoryArray(scip_, &const_row_);
   SCIPfreeMemoryArray(scip_, &const_vals_);
   const_nnz_ = new_const_nnz;
   SCIP_CALL( SCIPallocMemoryArray(scip_, &const_col_, const_nnz_));
   SCIP_CALL( SCIPallocMemoryArray(scip_, &const_row_, const_nnz_));
   SCIP_CALL( SCIPallocMemoryArray(scip_, &const_vals_, const_nnz_));

   for (int i = 0; i < const_nnz_; ++i)
   {
      const_col_[i] = new_const_col[i];
      const_row_[i] = new_const_row[i];
      const_vals_[i] = new_const_vals[i];
   }

   SCIPfreeBufferArray(scip_, &new_const_col);
   SCIPfreeBufferArray(scip_, &new_const_row);
   SCIPfreeBufferArray(scip_, &new_const_vals);

   int *new_row;
   int *new_col;
   double *new_vals;
   SCIP_VAR **new_vars;
   int* new_vbeg;


   SCIP_CALL( SCIPallocBufferArray(scip_, &new_row, new_nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_col, new_nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_vals, new_nnz) );

   int save_position = -1;
   int count_all = 0;
   int count_vars = 0;

   int no_more_there = 0;
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_vbeg, nvars_ + 1 - how_many_deleted) );
   SCIP_CALL( SCIPallocBufferArray(scip_, &new_vars, nvars_ - how_many_deleted) );


   for (int k = 0; k < nvars_; ++k)
   {
      if ((SCIPvarGetStatus(SCIPvarGetProbvar(uvars_[k])) != SCIP_VARSTATUS_FIXED) && (SCIPvarGetStatus(SCIPvarGetProbvar(uvars_[k])) != SCIP_VARSTATUS_AGGREGATED) && ! SCIPisEQ(scip_, SCIPvarGetLbLocal(SCIPvarGetProbvar(uvars_[k])),SCIPvarGetUbLocal(SCIPvarGetProbvar(uvars_[k]))))
      {
         new_vars[count_vars] = uvars_[k];
         new_vbeg[count_vars] = vbeg_[k] - no_more_there;
         count_vars++;
         for (int l = vbeg_[k]; l < vbeg_[k + 1]; ++l)
         {
            new_col[count_all] = col_[l];
            new_row[count_all] = row_[l];
            new_vals[count_all] = vals_[l];
            count_all++;
         }

         // look for other variables that are aggregated to this one (only needs to be done if it isn't deleted)
         for (int i = 0; i < count_aggr; ++i)
         {
            if (aggr_vars[i] == uvars_[k])
            {
               save_position = -1;
               for (int j = aggr_vbeg[i]; j < aggr_vbeg[i + 1]; ++j)
               {
                  for (int l = new_vbeg[count_vars]; l < count_all; ++l)
                  {
                     if (new_row[l] == aggr_row[j] && new_col[l] == aggr_col[j])
                     {
                        //there is already an entry at this position of the matrix for this variable
                        save_position = l;
                     }
                  }
                  if (save_position == -1)
                  {   //no entry at this postion
                     new_col[count_all] = aggr_col[j];
                     new_row[count_all] = aggr_row[j];
                     new_vals[count_all] = aggr_vals[j];
                     count_all++;
                  }
                  else
                  {
                     new_vals[save_position] += aggr_vals[j];
                  }
                  new_vbeg[count_vars] = count_all;
               }
            }
         }
      }
      else
      {
         no_more_there += vbeg_[k + 1] - vbeg_[k];
         SCIP_CALL(SCIPreleaseVar(scip_, &uvars_[k] ));
      }
   }

   new_vbeg[count_vars] = vbeg_[nvars_] - no_more_there;
   SCIPfreeBufferArray(scip_, &aggr_row);
   SCIPfreeBufferArray(scip_, &aggr_col);
   SCIPfreeBufferArray(scip_, &aggr_vals);
   SCIPfreeBufferArray(scip_, &aggr_vbeg);
   SCIPfreeBufferArray(scip_, &aggr_vars);

   SCIPfreeMemoryArray(scip_, &row_);
   SCIPfreeMemoryArray(scip_, &col_);
   SCIPfreeMemoryArray(scip_, &vals_);
   SCIPfreeMemoryArray(scip_, &uvars_);
   SCIPfreeMemoryArray(scip_, &vbeg_);
   nnz_ = new_nnz;
   nvars_ = nvars_ - how_many_deleted;
   SCIP_CALL( SCIPallocMemoryArray(scip_, &row_, new_nnz));
   SCIP_CALL( SCIPallocMemoryArray(scip_, &col_, new_nnz));
   SCIP_CALL( SCIPallocMemoryArray(scip_, &vals_, new_nnz));

   for (int i = 0; i < new_nnz; ++i)
   {
      row_[i] = new_row[i];
      vals_[i] = new_vals[i];
      col_[i] = new_col[i];
   }

   SCIP_CALL( SCIPallocMemoryArray(scip_, &uvars_, nvars_));
   SCIP_CALL( SCIPallocMemoryArray(scip_, &vbeg_, nvars_ + 1));
   for (int i = 0; i < nvars_; ++i)
   {
      vbeg_[i] = new_vbeg[i];
      uvars_[i] = new_vars[i];
   }
   vbeg_[nvars_] = new_vbeg[nvars_];

   SCIPfreeBufferArray(scip_, &new_vars);
   SCIPfreeBufferArray(scip_, &new_row);
   SCIPfreeBufferArray(scip_, &new_col);
   SCIPfreeBufferArray(scip_, &new_vals);
   SCIPfreeBufferArray(scip_, &new_vbeg);

   return SCIP_OKAY;
}

SCIP_RETCODE SdpCone::transform_vars()
{
   for (int k = 0; k < nvars_; ++k)
   {
      SCIP_VAR* tvar;
      SCIP_CALL(SCIPgetTransformedVar(scip_, uvars_[k], &tvar));
      SCIP_CALL(SCIPcaptureVar(scip_, tvar));
      SCIP_CALL(SCIPreleaseVar(scip_, &uvars_[k]));
      uvars_[k] = tvar;
   }

   return SCIP_OKAY;
}

void SdpCone::sort_const(int* const_col, int* const_row, double* const_vals)
{
   int* rctr;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &rctr, blocksize_));
   int* cctr;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &cctr, blocksize_));
   int* rbeg;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &rbeg, blocksize_ + 1));
   int* cbeg;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &cbeg, blocksize_ + 1));
   for (int k = 0; k < blocksize_; ++k)
   {
      cctr[k] = 0;
      rctr[k] = 0;
      rbeg[k] = 0;
      cbeg[k] = 0;
   }

   rbeg[blocksize_] = 0;
   cbeg[blocksize_] = 0;


   //count the number of entries for every row
   for (int i = 0; i < const_nnz_; ++i)
   {
      int crow = const_row[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (crow == k + 1)
         {
            rctr[k]++;
         }
      }
   }

   int ctr = 0;
   for (int k = 0; k < blocksize_; ++k)
   {
      rbeg[k] = ctr;
      ctr += rctr[k];
   }

   assert (ctr == const_nnz_);
   rbeg[blocksize_] = ctr;

   //count the number of entries for every col
   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
   }

   for (int i = 0; i < const_nnz_; ++i)
   {
      int ccol = const_col[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (ccol == k + 1)
         {
            cctr[k]++;
         }
      }
   }

   ctr = 0;
   for (int k = 0; k < blocksize_; ++k)
   {
      cbeg[k] = ctr;
      ctr += cctr[k];
   }
   assert (ctr == const_nnz_);
   cbeg[blocksize_] = ctr;

   // reset rctr
   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
      cctr[k] = 0;
   }
   int* tmp_row;
   int* tmp_col;
   double* tmp_vals;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_row, const_nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_col, const_nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_vals, const_nnz_));


   //sorting by column
   for (int i = 0; i < const_nnz_; ++i)
   {
      int ccol = const_col[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (ccol != k + 1)
            continue;

         int idx = cbeg[k] + cctr[k];
         assert(idx < cbeg[k + 1]);

         tmp_row[idx] = const_row[i];
         tmp_col[idx] = const_col[i];
         tmp_vals[idx] = const_vals[i];

         cctr[k]++;
         break;
      }
   }

   //sorting by row
   for (int i = 0; i < const_nnz_; ++i)
   {
      int crow = tmp_row[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (crow != k + 1)
            continue;

         int idx = rbeg[k] + rctr[k];
         assert(idx < rbeg[k + 1]);

         const_row_[idx] = tmp_row[i];
         const_col_[idx] = tmp_col[i];
         const_vals_[idx] = tmp_vals[i];

         rctr[k]++;
         break;
      }
   }

   SCIPfreeBufferArray(scip_, &tmp_col);
   SCIPfreeBufferArray(scip_, &tmp_row);
   SCIPfreeBufferArray(scip_, &tmp_vals);
   SCIPfreeBufferArray(scip_, &rbeg);
   SCIPfreeBufferArray(scip_, &cbeg);
   SCIPfreeBufferArray(scip_, &rctr);
   SCIPfreeBufferArray(scip_, &cctr);

}

void SdpCone::compress_representation(SCIP_VAR** vars,
   int* col,
   int* row,
   double* vals)
{
   // compute the vbeg array
   // make one pass through data and count occurences of vars
   // in tmp array vctr

   int* vctr;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &vctr, nvars_));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip_, &vbeg_, nvars_ + 1));
   int* rctr;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &rctr, blocksize_));
   int* rbeg;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &rbeg, blocksize_ + 1));
   int* cbeg;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &cbeg, blocksize_ + 1));

   for (int k = 0; k < nvars_ + 1; ++k)
   {
      vbeg_[k] = 0;
   }

   for (int k = 0; k < blocksize_ + 1; ++k)
   {
      rbeg[k] = 0;
      cbeg[k] = 0;
   }

   for (int k = 0; k < nvars_; ++k)
   {
      vctr[k] = 0;
   }

   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
   }

   //count the number of entries for every var
   for (int i = 0; i < nnz_; ++i)
   {
      SCIP_VAR* cvar = vars[i];
      for (int k = 0; k < nvars_; ++k)
      {
         if (cvar != uvars_[k])
         {
            continue;
         }

         vctr[k]++;
         break;
      }
   }

   int ctr = 0;
   for (int k = 0; k < nvars_; ++k)
   {
      vbeg_[k] = ctr;
      ctr += vctr[k];
   }
   assert (ctr == nnz_);
   vbeg_[nvars_] = ctr;


   //count the number of entries for every row
   for (int i = 0; i < nnz_; ++i)
   {
      int crow = row[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (crow == k + 1)
         {
            rctr[k]++;
         }
      }
   }

   ctr = 0;
   for (int k = 0; k < blocksize_; ++k)
   {
      rbeg[k] = ctr;
      ctr += rctr[k];
   }

   assert (ctr == nnz_);
   rbeg[blocksize_] = ctr;


   //count the number of entries for every col
   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
   }

   for (int i = 0; i < nnz_; ++i)
   {
      int ccol = col[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (ccol == k + 1)
         {
            rctr[k]++;
         }
      }
   }

   ctr = 0;
   for (int k = 0; k < blocksize_; ++k)
   {
      cbeg[k] = ctr;
      ctr += rctr[k];
   }
   assert (ctr == nnz_);
   cbeg[blocksize_] = ctr;

   // reset vctr
   for (int k = 0; k < nvars_; ++k)
   {
      vctr[k] = 0;
   }

   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
   }

   // now fill vars_, row_, col_ and vals_ according to vbeg
   // using vctr as local pointers

   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip_, &row_, nnz_));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip_, &col_, nnz_));
   SCIP_CALL_ABORT(SCIPallocMemoryArray(scip_, &vals_, nnz_));

   int* tmp_row;
   int* tmp_col;
   double* tmp_vals;
   SCIP_VAR** tmp_vars;
   SCIP_VAR** d_tmp_vars;
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_row, nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_col, nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_vals, nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &tmp_vars, nnz_));
   SCIP_CALL_ABORT(SCIPallocBufferArray(scip_, &d_tmp_vars, nnz_));

   //we need to sort out data
   //sorting by column
   for (int i = 0; i < nnz_; ++i)
   {
      int ccol = col[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (ccol != k + 1)
            continue;

         int idx = cbeg[k] + rctr[k];
         assert(idx < cbeg[k + 1]);

         row_[idx] = row[i];
         col_[idx] = col[i];
         vals_[idx] = vals[i];
         d_tmp_vars[idx] = vars[i];

         rctr[k]++;

         break;
      }
   }

   for (int k = 0; k < blocksize_; ++k)
   {
      rctr[k] = 0;
   }

   //sorting by row
   for (int i = 0; i < nnz_; ++i)
   {
      int crow = row_[i];
      for (int k = 0; k < blocksize_; ++k)
      {
         if (crow != k + 1)
            continue;

         int idx = rbeg[k] + rctr[k];
         assert(idx < rbeg[k + 1]);

         tmp_row[idx] = row_[i];
         tmp_col[idx] = col_[i];
         tmp_vals[idx] = vals_[i];
         tmp_vars[idx] = d_tmp_vars[i];

         rctr[k]++;
         break;
      }
   }


   //sorting by vars
   for (int i = 0; i < nnz_; ++i)
   {
      SCIP_VAR* cvar = tmp_vars[i];
      for (int k = 0; k < nvars_; ++k)
      {
         if (cvar != uvars_[k])
            continue;

         int idx = vbeg_[k] + vctr[k];
         assert(idx < vbeg_[k + 1]);

         row_[idx] = tmp_row[i];
         col_[idx] = tmp_col[i];
         vals_[idx] = tmp_vals[i];

         vctr[k]++;

         break;
      }
   }


   int k = 0;
   for (int i = 0; i < nnz_; ++i)
   {
      if (i>=vbeg_[k+1]) {
         k++;
      }


   }

   SCIPfreeBufferArray(scip_, &rbeg);
   SCIPfreeBufferArray(scip_, &cbeg);
   SCIPfreeBufferArray(scip_, &tmp_col);
   SCIPfreeBufferArray(scip_, &tmp_row);
   SCIPfreeBufferArray(scip_, &tmp_vals);
   SCIPfreeBufferArray(scip_, &tmp_vars);
   SCIPfreeBufferArray(scip_, &d_tmp_vars);
   SCIPfreeBufferArray(scip_, &rctr);
   SCIPfreeBufferArray(scip_, &vctr);
}



SCIP_RETCODE SdpCone::get_shrunk_constraint_matrix(double* matrix, int* actual_size, int vidx) const
{
   assert( matrix != 0 );
   assert( actual_size != 0 );

   // alloc array to store row/col sparsity pattern
   int* pattern;

   SCIP_CALL(SCIPallocBufferArray(scip_, &pattern, blocksize_ + 1));

   for (int j = 0; j < blocksize_ + 1; ++j)
   {
      pattern[j] = 0;
   }

   // find rows/cols that have a nonzero entry
   for (int i = vbeg_[vidx]; i < vbeg_[vidx + 1]; ++i)
   {
      pattern[row_[i]] = 1;
      pattern[col_[i]] = 1;
   }

   // update pattern array
   // it now contains the index of the row-col in the new matrix
   int new_size = 0;
   for (int j = 0; j < blocksize_ + 1; ++j)
   {
      int tmp = pattern[j];
      pattern[j] = new_size;
      new_size += tmp;
   }

   // new_size is now the actual_size
   *actual_size = new_size;

   // zero matrix

   for (int i = 0; i < blocksize_ * blocksize_; ++i)
   {
      matrix[i] = 0.0;
   }

   // fill matrix
   for (int i = vbeg_[vidx]; i < vbeg_[vidx + 1]; ++i)
   {
      int new_row = pattern[row_[i]]; // NB: one based indexing in row/col
      int new_col = pattern[col_[i]];

      matrix[new_col * new_size + new_row] += vals_[i] ;

      if (new_row != new_col)
      {
         matrix[new_row * new_size + new_col] += vals_[i];
      }
   }

   // free memory
   SCIPfreeBufferArray(scip_, &pattern);

   return SCIP_OKAY;
}

int SdpCone::get_nnz() const
{
   return nnz_;
}

int SdpCone::get_const_nnz() const
{
   return const_nnz_;
}


SCIP_RETCODE SdpCone::transform_vars_to_pos(SCIP* scip, SCIP_VAR** fixed_vars, int *n_fixed_vars)
{
   if (*n_fixed_vars > 0)
   {
      bool found = FALSE;
      if (pos_fixed_vars_)
      {
         SCIPfreeMemoryArray(scip_, &pos_fixed_vars_);
      }
      SCIP_CALL(SCIPallocMemoryArray(scip, &pos_fixed_vars_, *n_fixed_vars));
      int tmp_nfixed_vars = *n_fixed_vars;
      int count = 0;


      for (int i = 0; i < *n_fixed_vars; ++i)
      {
         found = FALSE;
         for (int j = 0; j < nvars_; ++j)
         {
            if (varEQ(uvars_[j], fixed_vars[i]))
            {
               pos_fixed_vars_[count] = j;
               found = TRUE;
               count++;
               break;
            }
         }
         if (!found)
         {
            tmp_nfixed_vars--;
         }
      }
      *n_fixed_vars = tmp_nfixed_vars;
   }


   return SCIP_OKAY;
}

int SdpCone::get_sdpcone_nfixed()
{

   sdpcone_nfixed_ = 0;
   for (int i = 0;  i < nvars_; ++i)
   {

      if (SCIPisEQ(scip_, SCIPvarGetLbLocal(uvars_[i]),SCIPvarGetUbLocal(uvars_[i])))
      {

         sdpcone_nfixed_++;
      }
   }
   return sdpcone_nfixed_;
}

SdpCone::LhsIterator SdpCone::lhs_begin(SCIP_VAR** fixed_vars, int n_fixed_vars)
{
   SCIP_CALL_ABORT(transform_vars_to_pos(scip_, fixed_vars, &n_fixed_vars));
   return LhsIterator(this, n_fixed_vars);
}

SdpCone::LhsIterator SdpCone::lhs_end()
{
   return SdpCone::LhsIterator();
}

SdpCone::RhsIterator SdpCone::rhs_begin(SCIP_VAR** fixed_vars, int n_fixed_vars, double* fixed_values)
{
   SCIP_CALL_ABORT(transform_vars_to_pos(scip_, fixed_vars, &n_fixed_vars));
   return SdpCone::RhsIterator(this, n_fixed_vars, fixed_values);

}

SdpCone::RhsIterator SdpCone::rhs_end()
{
   return SdpCone::RhsIterator();
}


SdpCone::RhsIterator::RhsIterator(SdpCone* c, int n_fixed_vars, double* fixed_values) : sdpcone_(c), pos_(0), epos_(NULL), curr_element_(), n_fixed_vars_(n_fixed_vars), fixed_values_(fixed_values), end_(FALSE)
{
   SCIP_CALL_ABORT(SCIPallocMemoryArray(sdpcone_->scip_, &epos_, n_fixed_vars_));
   //initialize epos
   for (int i = 0; i < n_fixed_vars_; ++i)
   {
      epos_[i] = sdpcone_->vbeg_[sdpcone_->pos_fixed_vars_[i]];
   }

   if (sdpcone_->const_nnz_ == 0 && n_fixed_vars_ == 0)
   {
      end_ = TRUE;
      return;
   }
   else if (sdpcone_->const_nnz_ == 0 && n_fixed_vars_ != 0)
   {
      pos_ = -1;
   }

   int crow = -1;
   int ccol = -1;
   double val = 0.0;

   bool stop = FALSE;
   //look for rows
   for (int k = 0; k < sdpcone_->blocksize_; k++)
   {//look for cols
      for (int j = k; j < sdpcone_->blocksize_; j++)
      {
         crow = j; //row and cols are swaped
         ccol = k;
         val = 0.0;

         if (pos_ != -1 && (crow == sdpcone_->const_col_[pos_] - 1) && (ccol == sdpcone_->const_row_[pos_] - 1)) // because it was checked that there is a nonzero, the check pos_ < const_nonzeroes isn't needed
         {
            val = -sdpcone_->const_vals_[pos_];
            pos_++;
            stop = TRUE;
         }
         for (int i = 0; i < n_fixed_vars_; ++i)
         {
            if (epos_[i] != -1 && crow == (sdpcone_->col_[epos_[i]] - 1) && ccol == (sdpcone_->row_[epos_[i]] - 1))
            {
               val += sdpcone_->vals_[epos_[i]] * fixed_values_[i];
               (epos_[i])++;
               stop = TRUE;
            }
         }
         if (stop)
         {
            break;
         }
      }
      if (stop)
      {
         break;
      }
   }
   curr_element_.row = crow;
   curr_element_.col = ccol;
   curr_element_.val = val;
   curr_element_.vidx = -1;
   curr_element_.eidx = (crow * (crow + 1) / 2 + ccol);
}


SdpCone::RhsIterator::RhsIterator() : sdpcone_(NULL), pos_(-1), epos_(NULL), curr_element_(), n_fixed_vars_(0), fixed_values_(NULL), end_(TRUE) {}


SdpCone::RhsIterator::~RhsIterator()
{
   if (epos_)
   {
      SCIPfreeMemoryArray(sdpcone_->scip_, &epos_);
   }
   sdpcone_ = NULL;
   epos_ = NULL;
}


SdpCone::RhsIterator& SdpCone::RhsIterator::operator++()
{
   bool done = TRUE;
   for (int i = 0; i < n_fixed_vars_; ++i)
   {
      if (epos_[i] == -1)
         continue;

      if ((epos_[i] != -1) && (epos_[i] >= sdpcone_->vbeg_[sdpcone_->pos_fixed_vars_[i] + 1]))
      {
         epos_[i] = -1;
      }
      else
      {
         done = FALSE;
      }
   }

   if (pos_ != -1 && pos_ >= sdpcone_->const_nnz_)
   {
      pos_ = -1;
   }

   // set to end iterator if we are done
   if ((done || n_fixed_vars_ == 0) && pos_ == -1)
   {
      end_ = TRUE;
      return *this;
   }

   int crow;
   int ccol;
   double val = 0.0;

   bool stop = FALSE;
   //look for cols
   for (ccol = 0; ccol < sdpcone_->blocksize_; ccol++)
   {//look for rows
      for (crow = ccol; crow < sdpcone_->blocksize_; crow++)
      {
         val = 0.0;

         if (pos_ != -1 && (crow == sdpcone_->const_col_[pos_] - 1) && (ccol == sdpcone_->const_row_[pos_] - 1))
         {
            val = -sdpcone_->const_vals_[pos_];
            pos_++;
            stop = TRUE;
         }
         for (int i = 0; i < n_fixed_vars_; ++i)
         {
            if (epos_[i] != -1 && crow == (sdpcone_->col_[epos_[i]] - 1) && ccol == (sdpcone_->row_[epos_[i]] - 1))
            {
               val += sdpcone_->vals_[epos_[i]] * fixed_values_[i];
               (epos_[i])++;
               stop = TRUE;
            }
         }
         if (stop)
         {
            break;
         }
      }
      if (stop)
      {
         break;
      }
   }

   curr_element_.row = crow;
   curr_element_.col = ccol;
   curr_element_.val = val;
   curr_element_.vidx = -1;
   curr_element_.eidx = (crow * (crow + 1) / 2 + ccol);

   return *this;
}

SdpCone::RhsIterator SdpCone::RhsIterator::operator++(int) {
   SdpCone::RhsIterator tmp(*this);
   ++(*this);
   return tmp;
}

const SdpCone::element& SdpCone::RhsIterator::operator*() const
{
   return curr_element_;
}

const SdpCone::element* SdpCone::RhsIterator::operator->() const
{
   return &curr_element_;
}

bool SdpCone::RhsIterator::operator==(const SdpCone::RhsIterator& other)
{
   bool both_at_end = end_ && other.end_;
   bool both_not_at_end = !end_ && !other.end_;
   bool same_position = (epos_ == other.epos_) && (pos_ == other.pos_);

   bool result = (both_at_end || (both_not_at_end && same_position));
   return result;
}

bool SdpCone::RhsIterator::operator!=(const SdpCone::RhsIterator& other)
{
   return !(*this == other);
}



SdpCone::LhsIterator::LhsIterator(SdpCone* c, int n_fixed_vars) : sdpcone_(c), vpos_(0), epos_(0), n_fixed_vars_(n_fixed_vars), curr_element_(), end_(FALSE)
{
   if (sdpcone_->nnz_ == 0 || n_fixed_vars_ == sdpcone_->get_nvars())
   {
      end_ = TRUE;
   }
   else
   {
      int i = -1;
      bool find_first_unfixed = TRUE;
      if (n_fixed_vars > 0)
      {
         while (find_first_unfixed)
         {
            ++i;

            if (i >= n_fixed_vars_ || (i != -1  && sdpcone_->pos_fixed_vars_[i] != i)) // this works for the first unfixed, because they are entered in pos_fixed with ascending index, so if the first vars are fixed then var 0 is in pos_fixed[0], var[1] in pos_fixed[1], ...
            {
               epos_ = sdpcone_->vbeg_[i];
               vpos_ = i;
               find_first_unfixed = FALSE;
            }
         }
      }
      else
      {
         epos_ = 0;
      }

      int row = sdpcone_->col_[epos_] - 1;
      int col = sdpcone_->row_[epos_] - 1;
      double val = sdpcone_->vals_[epos_];

      curr_element_.row = row;
      curr_element_.col = col;
      curr_element_.val = val;
      curr_element_.vidx = vpos_;
      curr_element_.eidx = (row * (row + 1) / 2 + col);
   }
}


SdpCone::LhsIterator::LhsIterator() : sdpcone_(NULL), vpos_(-1), epos_(-1),  n_fixed_vars_(0),curr_element_(), end_(TRUE){}

SdpCone::LhsIterator& SdpCone::LhsIterator::operator++()
{
   // increment element
   ++epos_;

   // check whether we need to increment the variable
   if (epos_ >= sdpcone_->vbeg_[vpos_ + 1])
   {
      ++vpos_;

      bool var_is_fixed = FALSE;
      for (int i = 0; i < n_fixed_vars_; ++i)
      {
         if (sdpcone_->pos_fixed_vars_[i] == vpos_)
         {
            var_is_fixed = TRUE;
            break;
         }
      }

      if (var_is_fixed)
      {
         bool find_next_unfixed = TRUE;
         while (find_next_unfixed && vpos_ <= sdpcone_->nvars_)
         {
            ++vpos_;
            epos_ = sdpcone_->vbeg_[vpos_];
            for (int i = 0; i < n_fixed_vars_; ++i)
            {
               find_next_unfixed = FALSE;

               if (sdpcone_->pos_fixed_vars_[i] == vpos_)
               {
                  find_next_unfixed = TRUE;
                  break;
               }
            }
         }
      }
   }

   // set to end iterator if we are done
   if (vpos_ == sdpcone_->nvars_)
   {
      end_ = TRUE;
   }
   else
   {
      int row = -1;
      int col = -1;
      double val = 0.0;
      if (epos_ < sdpcone_->get_nnz())
      {
         row = sdpcone_->col_[epos_] - 1;
         col = sdpcone_->row_[epos_] - 1;
         val = sdpcone_->vals_[epos_];
      }

      curr_element_.row = row;
      curr_element_.col = col;
      curr_element_.val = val;
      curr_element_.vidx = vpos_;
      curr_element_.eidx = (row * (row + 1) / 2 + col);
   }

   return *this;
}

SdpCone::LhsIterator SdpCone::LhsIterator::operator++(int) {
   SdpCone::LhsIterator tmp(*this);
   ++(*this);
   return tmp;
}

const SdpCone::element& SdpCone::LhsIterator::operator*() const
{
   return curr_element_;
}

const SdpCone::element* SdpCone::LhsIterator::operator->() const
{
   return &curr_element_;
}

bool SdpCone::LhsIterator::operator==(const SdpCone::LhsIterator& other)
{
   bool both_at_end = end_ && other.end_;
   bool both_not_at_end = !end_ && !other.end_;
   bool same_position = (epos_ == other.epos_) && (vpos_ == other.vpos_);
   bool result = (both_at_end || (both_not_at_end && same_position));
   return result;
}

bool SdpCone::LhsIterator::operator!=(const SdpCone::LhsIterator& other)
{
   return !(*this == other);
}
