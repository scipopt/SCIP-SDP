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

/**@file   SdpProblem.cpp
 * @brief  Class, where the sdp-data and lp-data is stored
 * @author Sonja Mars
 */

#include "SdpProblem.h"

#include <cassert>                      // for assert
#include <cstring>                      // for NULL, strcmp
#include <vector>                       // for vector

#include "SdpVarMapper.h"               // for SdpVarMapper
#include "objconshdlr_sdp.h"            // for getSdpCone
#include "scip/pub_cons.h"              // for SCIPconsGetHdlr, etc
#include "scip/pub_lp.h"                // for SCIPcolGetVar, etc
#include "scip/pub_message.h"           // for SCIPdebugMessage
#include "scip/pub_var.h"               // for SCIPvarGetUbLocal, etc
#include "scip/scip.h"                  // for SCIPisEQ, SCIPinfinity, etc
#include "scip/type_cons.h"             // for SCIP_CONS, SCIP_CONSHDLR
#include "scip/type_var.h"              // for SCIP_VAR

/**Method for adding linear constraints to the structure we need for dsdp
 *the arrays for_* will later be added to dsdp as lp cone*/
SCIP_RETCODE SdpProblem::addconstraint(
   SdpVarMapper* varmapper,         /**<varmapper class data*/
   int *position,                    /**<diagonal-position where constraint should be added*/
   double rhs,                       /**<rhs of constraint*/
   int* lininds,                     /**<indices of variables in constraint*/
   int nlininds,                     /**<number of variables in constraint*/
   SCIP_Real* vals,                  /**<coefficients of variables in constraint*/
   SCIP_Real* ubs                    /**upper bound of variable<*/)
{
   bool increase_position = FALSE;
   int nmat = 0;
   int remember_size = for_vals_.size() - 1;
   bool something_over = TRUE;
   int l = 0;

   //rhs
   if (varmapper->get_nfixed() == 0)
   {
      if (!SCIPisEQ(scip_, rhs, 0.0))
      {
         for_vals_.push_back(rhs);
         for_matind_.push_back(*position);
         for_constraint_.push_back(nmat);
         remember_size = for_vals_.size() - 1;
         something_over = TRUE;
      }
   }
   else  //there are fixed variables
   {
      while (l < nlininds)
      {
         if (lininds[l] != -1)
         {
            //if there exists at least one variables
            if (!SCIPisEQ(scip_, rhs, 0.0))  //add everything if the rhs is not zero
            {
               for_vals_.push_back(rhs);
               for_matind_.push_back(*position);
               for_constraint_.push_back(nmat);
               remember_size = for_vals_.size() - 1;
               something_over = TRUE;
            }
            break;//stop if we found one variable for every equality
         }
         else
         {
            l++;
            something_over = FALSE;
         }
      }
   }

   //vars and indices...
   for (int k = 0; k < nlininds; k++)
   {
      if (lininds[k] == -1)
      {
         //lininds[k] is gefixed
         if (SCIPisEQ(scip_, rhs, 0.0))
         {
            //the rhs of the inequality was zero, so there is no entry and we can just add it
            if (!SCIPisEQ(scip_, vals[k]*ubs[k], 0.0))  //we only add, if the new value is not zero
            {
               for_matind_.push_back(*position);
               for_constraint_.push_back(0);
               for_vals_.push_back(vals[k]*ubs[k]);
            }
         }
         else //the rhs was not zero to it is possible that an enty exists
         {
            if (SCIPisEQ(scip_, ubs[k], 0.0))
            {
               //do nothing
            }
            else  //the variable is fixed, so we have to do something
            {
               if (something_over == TRUE) //maybe there is nothing left from our inequality, only to something, if anything is left
               {
                  //the rhs has another sign so do +=
                  for_vals_[remember_size] += vals[k] * ubs[k];
               }
            }
         }
      }
      else  //add everything if nothing was fixed
      {
         increase_position = TRUE;
         nmat = lininds[k] + 1;
         for_matind_.push_back(*position);
         for_constraint_.push_back(nmat);
         for_vals_.push_back(vals[k]);
      }
   }

   if (something_over &&  remember_size >= 0 && SCIPisEQ(scip_, for_vals_[remember_size], 0.0) )
   {
      for_vals_[remember_size] = 0.0;
      for_matind_[remember_size] = -1;
      for_constraint_[remember_size] = -1;
   }
   if (increase_position)
   {
      (*position)++;
   }
   return SCIP_OKAY;
}

/**adds bound for variables*/
SCIP_RETCODE SdpProblem::addbound(
   int *position,                     /**<position where bound constraint is added*/
   double rhs,                        /**<rhs of variable*/
   int lininds,                       /**<variable index*/
   double linvals,                    /**<coefficient of variable*/
   SdpVarMapper* varmapper            /**<varmapper class object*/)
{
   int nmat = 0;
   if (!SCIPisEQ(scip_, rhs, 0.0))
   {
      for_vals_.push_back(rhs);
      for_matind_.push_back(*position);
      for_constraint_.push_back(nmat);
   }
   nmat = lininds + 1;
   for_matind_.push_back(*position);
   for_constraint_.push_back(nmat);
   for_vals_.push_back(linvals);
   (*position)++;

   return SCIP_OKAY;
}

/**adding lp data out of scip to the vector ding->for_matind, ding->for_block, ding->for_constraint, ding->for_vals*/
SCIP_RETCODE SdpProblem::get_rows_data(
   SdpVarMapper* varmapper,          /**<varmapper class data*/
   SCIP_ROW** rows,                  /**<rows to add to lp-block*/
   int nrows,                        /**<number of rows*/
   int* position                     /**<pointer to store position of diagonal entries of lp block*/
  )
{
   for (int i = 0 ; i < nrows; i++)
   {
      const SCIP_Real lhs = SCIProwGetLhs(rows[i])- SCIProwGetConstant(rows[i]);
      const SCIP_Real rhs = SCIProwGetRhs(rows[i])- SCIProwGetConstant(rows[i]);

      const int nlininds = SCIProwGetNNonz(rows[i]);

      const SCIP_Real* vals;
      SCIP_COL* const* row_cols;

      vals = SCIProwGetVals(rows[i]);
      row_cols = SCIProwGetCols(rows[i]);

      int* lininds;
      SCIP_CALL( SCIPallocBufferArray(scip_, &lininds, nlininds) );
      SCIP_Real* ubs;
      SCIP_CALL( SCIPallocBufferArray(scip_, &ubs, nlininds) );

      for (int j = 0; j < nlininds; j++)
      {
         SCIP_VAR* var = SCIPcolGetVar(row_cols[j]);
         lininds[j] = varmapper->get_sdp_index(var);
         ubs[j] = SCIPvarGetUbLocal(var);
      }

      SCIP_Real* vals_for_change;
      SCIP_CALL( SCIPallocBufferArray(scip_, &vals_for_change, nlininds) );

      if (lhs > -SCIPinfinity(scip_))
      {
         for (int j = 0; j < nlininds; j++)
         {
            vals_for_change[j] = vals[j];
         }
         SCIP_CALL( addconstraint(varmapper, position, -lhs, lininds, nlininds, vals_for_change, ubs) );
      }

      if (rhs < SCIPinfinity(scip_))
      {
         for (int j = 0; j < nlininds; j++)
         {
            vals_for_change[j] = -vals[j];
         }
         SCIP_CALL( addconstraint(varmapper, position, rhs, lininds, nlininds, vals_for_change, ubs) );
      }
      SCIPfreeBufferArray(scip_, &vals_for_change);
      SCIPfreeBufferArray(scip_, &lininds);
      SCIPfreeBufferArray(scip_, &ubs);
   }
   return SCIP_OKAY;
}


SdpProblem::SdpProblem(SCIP* scip, SdpVarMapper* varmapper) :
   scip_(scip),
   sdpcones_(NULL),
   nsdpcones_(0),
   for_vals_(),
   for_constraint_(),
   for_matind_(),
   size_lpblock_(0)
{
   //get constraints from scip
   const int nconss   = SCIPgetNConss(scip);
   SCIP_CONS** conss;
   SCIP_CONSHDLR* hdlr;

   conss = SCIPgetConss(scip);
   //find out how many sdpcones exists
   for (int i = 0; i < nconss; i++)
   {
      assert(conss[i] != NULL);

      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);

      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      if ( strcmp(hdlrName, "SDP") == 0)
      {
         nsdpcones_ = SCIPconshdlrGetNConss(hdlr);
         break;
      }
   }

   int position = 0;
   //get all the lp data
   //get cols
   SCIP_COL** cols;
   int        ncols;
   SCIP_CALL_ABORT( SCIPgetLPColsData(scip, &cols, &ncols) );

   //get rows
   SCIP_ROW** rows;
   int        nrows;
   SCIP_CALL_ABORT( SCIPgetLPRowsData(scip, &rows, &nrows) );


   SCIPdebugMessage("nrows: %d  ncols: %d  \n", nrows, ncols);

   //write lp data in vectors of sdpproblem
   //A_0 is u and -l in the linear block and A_i
   //LP
   SCIP_CALL_ABORT( get_rows_data(varmapper, rows, nrows, &position));


   for (int i = 0; i < varmapper->get_sdp_nvars(); i++)
   {

      const SCIP_Real ub = SCIPvarGetUbLocal(varmapper->get_scip_var(i));
      const SCIP_Real lb = SCIPvarGetLbLocal(varmapper->get_scip_var(i));

      assert ( ub > lb);

      if (SCIPisGT(scip, lb, -SCIPinfinity(scip) ) )
      {
         double rhs = -lb;
         double linvals = 1.0;
         addbound(&position, rhs, i, linvals, varmapper);
      }
      if (SCIPisLT(scip, ub, SCIPinfinity(scip) ) )
      {
         double rhs = ub;
         double linvals = - 1.0;
         addbound(&position, rhs, i, linvals, varmapper);
      }
   }
   size_lpblock_ = position;

   SCIP_CALL_ABORT(SCIPallocBlockMemoryArray(scip, &sdpcones_, nsdpcones_));
   int count = 0;
   SdpCone* tmp_sdpcone;
   //get all the sdp constraints and add them to the sdpdata
   for (int i = 0; i < nconss; i++)
   {
      tmp_sdpcone = NULL;
      assert(conss[i] != NULL);

      hdlr = SCIPconsGetHdlr(conss[i]);
      assert(hdlr != NULL);

      const char* hdlrName;
      hdlrName = SCIPconshdlrGetName(hdlr);

      if ( strcmp(hdlrName, "SDP") != 0)
         continue;

      SCIP_CALL_ABORT(getSdpCone(scip, conss[i] , &sdpcones_[count]));
      count++;
   }
}

SdpProblem::~SdpProblem()
{
}

SdpCone* SdpProblem::SdpProblem::get_sdpcone(int i) const
{
   return sdpcones_[i];
}

int SdpProblem::get_nsdpcones() const
{
   return nsdpcones_;
}

const double* SdpProblem::get_for_vals() const
{
   return &for_vals_[0];
}

int SdpProblem::get_for_vals_size() const
{
   return for_vals_.size();
}

const int* SdpProblem::get_for_constraint() const
{
   return &for_constraint_[0];
}

int SdpProblem::get_for_constraint_size() const
{
   return for_constraint_.size();
}

const int* SdpProblem::get_for_matind() const
{
   return &for_matind_[0];
}

int SdpProblem::get_for_matind_size() const
{
   return for_matind_.size();
}

int SdpProblem::get_size_lpblock() const
{
   return size_lpblock_;
}

int SdpProblem::get_lp_nnz() const
{
   return for_matind_.size();
}



