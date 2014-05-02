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

/**@file   DsdpInterface.cpp
 * @brief  inteface to dsdp-solver
 * @author Sonja Mars
 */
#include "DsdpInterface.h"
#include <cstdio>                       // for printf, NULL

#include "dsdp5.h"                      // for DSDPUsePenalty, etc
#include "dsdpmem.h"                    // for DSDPCALLOC2, DSDPFREE

#include "SdpCone.h"                    // for SdpCone, SdpCone::element, etc
#include "SdpProblem.h"                 // for SdpProblem
#include "SdpVarMapper.h"               // for SdpVarMapper
#include "blockmemshell/memory.h"       // for BMSallocMemoryArray
//#include "scip/type_prop.h"             // for SCIP_PROP, etc
//#include "scip/pub_message.h"           // for SCIPerrorMessage, etc
//#include "scip/pub_tree.h"              // for SCIPnodeGetLowerbound, etc
//#include "scip/pub_var.h"               // for SCIPvarGetUbLocal, etc
#include "scip/scip.h"                  // for SCIPgetNVars, etc
//#include "scip/type_tree.h"             // for SCIP_NODE
//#include "scip/type_var.h"              // for SCIP_VAR

/**Method to convert a info, which will be returned by DSDP into an SCIP return-code*/
#define INFO_TO_SCIPCALL(x)   do                                \
   {                                                            \
      SCIP_RETCODE _restat_;                                    \
      if( ((x)) != 0 )                                          \
      {                                                         \
         _restat_ = SCIP_ERROR;                                 \
         SCIPerrorMessage("Error <%d> in function call\n", x);  \
         return _restat_;                                       \
      }                                                         \
   }                                                            \
   while( FALSE )

/**The following methods are taken completely or partly from DSDP5.8 file readspda.c written by     Steven J. Benson and Yinyu Ye
 *(C) COPYRIGHT 2004 UNIVERSITY OF CHICAGO*/

/**fuction from dsdp: gets marker of arrays*/
static int GetMarkers(int block, int constraint, int blockn[], int constraintn[],
   int*m3)
{
   int i = 0;
   while (blockn[i] == block && constraintn[i] == constraint)
   {
      i++;
   }
   *m3 = i;
   return 0;
}

/**function from dsdp: counts nonzeros*/
static int CountNonzeroMatrices(int block, int blockn[], int constraintn[], int*m3)
{
   int i = 0, cvar = -1, nnzmats = 0;
   while (blockn[i] == block)
   {
      if (constraintn[i] > cvar)
      {
         cvar = constraintn[i];
         nnzmats++;
      }
      i++;
   }
   *m3 = nnzmats;
   return 0;
}

/**function from dsdp: checks if there is a constant matrix*/
static int CheckForConstantMat(double v[], int nnz, int n)
{
   int i;
   double vv;
   if (n <= 1)
   {
      return 0;
   }
   if (nnz != (n * n + n) / 2)
   {
      return 0;
   }
   for (vv = v[0], i = 1; i < nnz; i++)
   {
      if (v[i] != vv)
      {
         return 0;
      }
   }
   return 1;
}

/**function from dsdp: partitions a set*/
static int partition(int list1[], int list2[], int list3[], double list5[], int lstart, int lend)
{
   int k = lend;
   int pivot1 = list1[k], pivot2 = list2[k], pivot3 = list3[k];
   double pivot5 = list5[k];
   int bottom = lstart - 1, top = lend;
   int done = 0;
   int ordered = 1;
   while (!done)
   {

      while (!done)
      {
         bottom = bottom + 1;

         if (bottom == top)
         {
            done = 1;
            break;
         }
         if ( list1[bottom] > pivot1 ||
            (list1[bottom] == pivot1 && list2[bottom] > pivot2) ||
            (list1[bottom] == pivot1 && list2[bottom] == pivot2 &&
               list3[bottom] > pivot3) )
         {
            list1[top] = list1[bottom];
            list2[top] = list2[bottom];
            list3[top] = list3[bottom];
            list5[top] = list5[bottom];
            ordered = 0;
            break;
         }
      }
      while (!done)
      {
         top = top - 1;

         if (top == bottom)
         {
            done = 1;
            break;
         }
         if ( list1[top] < pivot1 ||
            (list1[top] == pivot1 && list2[top] < pivot2) ||
            (list1[top] == pivot1 && list2[top] == pivot2 && list3[top] < pivot3))
         {
            list1[bottom] = list1[top];
            list2[bottom] = list2[top];
            list3[bottom] = list3[top];
            list5[bottom] = list5[top];
            ordered = 0;
            break;
         }
      }
   }
   list1[top] = pivot1;
   list2[top] = pivot2;
   list3[top] = pivot3;
   list5[top] = pivot5;

   ordered = 0;

   if (bottom == lend)
   {
      ordered = 1;
      for (k = lstart; k < lend - 1; k++)
      {
         if ( (list1[k] > list1[k + 1]) ||
            (list1[k] == list1[k + 1] && list2[k] > list2[k + 1]) ||
            (list1[k] == list1[k + 1] && list2[k] == list2[k + 1] && list3[k] > list3[k + 1]) )
         {
            ordered = 0;
            break;
         }
      }
   }
   if (ordered && lend - lstart > 2)
   {
      top = (lend + lstart) / 2;
   }
   return top;
}

/**function from dsdp: sorts the arrays*/
static int qusort(int list1[], int list2[], int list3[], double list5[], int lstart, int lend)
{
   int split;
   if (lstart < lend)
   {
      split = partition(list1, list2, list3, list5, lstart, lend);
      qusort(list1, list2, list3, list5, lstart, split - 1);
      qusort(list1, list2, list3, list5, split + 1, lend);
   }
   else
   {
      return 0;
   }
   return 0;
}


/**sets the objectieve with the right variable order*/
static SCIP_RETCODE set_objective(
   DSDP dsdp,               /**<DSDP object*/
   SCIP_VAR** vars,         /**<scip variables*/
   int nvars,               /**<number of variables*/
   SdpVarMapper* varmapper  /**<varmapper class object*/)
{
   ////Set Dual Objective vector
   for (int i = 0; i < nvars; ++i)
   {
      //if variable is not fixed
      int index = varmapper->get_sdp_index(vars[i]);
      if (index != -1)
      {
         const double obj = SCIPvarGetObj(vars[i]);
         int info = 0;
         info = DSDPSetDualObjective(dsdp, index + 1, obj);
         INFO_TO_SCIPCALL(info);
      }
   }
   return SCIP_OKAY;
}

/**converts a dsdp solution with only non-fixed variables to a scip solution, changes order of variables and puts fixed variables in it*/
static SCIP_RETCODE convert_sol(
   SCIP* scip,             /**<SCIP data structure*/
   DSDP dsdp,              /**<DSDP object*/
   double* sol_for_scip,   /**<pointer to store solution in scip format*/
   SdpVarMapper* varmapper /**<varmapper class object*/)
{
   SCIP_VAR** vars;
   vars = SCIPgetVars (scip);
   const int nvars = SCIPgetNVars(scip);
   int info = 0;


   double* sol;
   DSDPCALLOC2(&sol, double, varmapper->get_sdp_nvars(), &info);
   INFO_TO_SCIPCALL(info);
   DSDPGetY(dsdp, sol, varmapper->get_sdp_nvars());

   for (int i = 0; i < nvars; i++)
   {
      int index = varmapper->get_sdp_index(vars[i]);
      if (index != -1)
      {
         sol_for_scip[i] = -sol[index];
      }
      else
      {
         sol_for_scip[i] = SCIPvarGetUbLocal(vars[i]);
      }
   }


   DSDPFREE(&sol, &info);
   return SCIP_OKAY;
}


static
SCIP_RETCODE transform_data(
   SCIP* scip,              /**<SCIP data structure*/
   SdpProblem* problemdata, /**<class with problemdata of a node*/
   int** block,             /**<pointer to store blocknumbers of nonzero entries*/
   int** matind,            /**<pointer to store indices for position in the matrix of nonzero entries*/
   int** constraint,        /**<pointer to store constraint number, that is the number of the variable this nonzero entry belongs to*/
   double** nnz,            /**<pointer to store nonzero entry*/
   int** blocksizes,        /**<pointer to store size of blocks*/
   char** conetypes,        /**<pointer to store type of cones (sdp or lp)*/
   SdpVarMapper* varmapper  /**<varmapper class object*/
)
{
   int nsdpcones = problemdata->get_nsdpcones();

   int info;
   SdpCone* sdpcone;
   int len_of_arrays = 0;
   for (int i = 0; i < nsdpcones; i++)
   {
      len_of_arrays += problemdata->get_sdpcone(i)->get_nnz() + problemdata->get_sdpcone(i)->get_const_nnz();
   }

   DSDPCALLOC2(blocksizes, int, (problemdata->get_nsdpcones() + 1), &info);
   INFO_TO_SCIPCALL(info);
   DSDPCALLOC2(conetypes, char, (problemdata->get_nsdpcones() + 1), &info);
   INFO_TO_SCIPCALL(info);

   int* matind_tmp;
   int* block_tmp;
   int* constraint_tmp;
   double* vals_tmp;
   int* row_tmp;
   int* col_tmp;

   SCIP_CALL(SCIPallocBufferArray(scip, &matind_tmp, len_of_arrays));
   SCIP_CALL(SCIPallocBufferArray(scip, &vals_tmp, len_of_arrays));
   SCIP_CALL(SCIPallocBufferArray(scip, &block_tmp, len_of_arrays));
   SCIP_CALL(SCIPallocBufferArray(scip, &constraint_tmp, len_of_arrays));
   SCIP_CALL(SCIPallocBufferArray(scip, &row_tmp, len_of_arrays));
   SCIP_CALL(SCIPallocBufferArray(scip, &col_tmp, len_of_arrays));

   (*conetypes)[nsdpcones] = 'L';
   (*blocksizes)[nsdpcones] = problemdata->get_size_lpblock();


   SCIP_VAR** vars = SCIPgetVars(scip);
   SCIP_VAR** fixed_vars;
   int n_fixed_vars = varmapper->get_nfixed();
   SCIP_CALL(SCIPallocBufferArray(scip, &fixed_vars, n_fixed_vars));
   double* fixed_values;
   SCIP_CALL(SCIPallocBufferArray(scip, &fixed_values, n_fixed_vars));

   int count = 0;

   for (int i = 0; i < SCIPgetNVars(scip); ++i)
   {
      if (varmapper->get_sdp_index(vars[i]) == -1 )
      {
         fixed_vars[count] = vars[i];
         fixed_values[count] = SCIPvarGetUbLocal(vars[i]);
         count++;
         assert (SCIPvarGetUbLocal(vars[i]) == SCIPvarGetLbLocal(vars[i]));

      }
   }

   count = 0;
      int save_ctr = 0;

   for (int i = 0; i < nsdpcones; ++i)
   {
      (*blocksizes)[i] = problemdata->get_sdpcone(i)->get_blocksize();
      (*conetypes)[i] = 'S';

      sdpcone = problemdata->get_sdpcone(i);

      //A_0 for SDP
      for ( SdpCone::RhsIterator it = sdpcone->rhs_begin(fixed_vars, n_fixed_vars, fixed_values); it != sdpcone->rhs_end(); ++it)
      {
         SdpCone::element el = *it;
         if (el.val != 0)
         {
            matind_tmp[count] = el.eidx;
            row_tmp[count] =  el.row;
            col_tmp[count] = el.col;
            //packed format for only saving upper triangle matrix as rows (or as cols)
            //this is zero based
            block_tmp[count] = i + 1;
            constraint_tmp[count] = 0;
            vals_tmp[count] = el.val;//wenn hier kein - steht, geht truss1 schief

            count++;
         }
      }

      //all A_i, i=1,...n for SDP

      for ( SdpCone::LhsIterator it = sdpcone->lhs_begin(fixed_vars, n_fixed_vars); it != sdpcone->lhs_end(); ++it)
      {
         SdpCone::element el = *it;
         SCIP_VAR* var = sdpcone->get_var(el.vidx);

         matind_tmp[count] = el.eidx;
         row_tmp[count] =  el.row;
         col_tmp[count] = el.col;
         //packed format for only saving upper triangle matrix as rows (or as cols)
         //this is zero based
         block_tmp[count] = i + 1;
         constraint_tmp[count] = varmapper->get_sdp_index(var) + 1;
         vals_tmp[count] = el.val;

         count++;
      }

      //delete all indices j for which both row j and column j are completely zero
      int* found;
      SCIP_CALL(SCIPallocBufferArray(scip, &found, sdpcone->get_blocksize()));
      for (int k = 0; k < sdpcone->get_blocksize(); ++k)
      {
         found[k] = 0;
      }

      for (int j = 0; j < sdpcone->get_blocksize(); ++j)
      {
         for (int k = save_ctr ; k < count; ++k)
         {
            if (row_tmp[k] == j || col_tmp[k] == j)
            {
               found[j] = 1;
               break;
            }
         }
      }


      int num_not_deleted = 0;
      for (int j = 0; j < sdpcone->get_blocksize(); ++j)
      {
         num_not_deleted += found[j];
      }


      (*blocksizes)[i] = num_not_deleted;

      int sum_del = 0;

      int row_and_col_to_del = sdpcone->get_blocksize() + 5;

      if (num_not_deleted != sdpcone->get_blocksize())
      {
         for (int j = 0; j < sdpcone->get_blocksize(); ++j)
         {
            if (found[j] == 0)
            {
               row_and_col_to_del = j - sum_del;
               sum_del++;
               for (int k = save_ctr; k < count; k++)
               {
                  if (row_tmp[k] >= row_and_col_to_del)
                  {
                     row_tmp[k] = row_tmp[k] - 1;
                  }
                  if (col_tmp[k] >= row_and_col_to_del)
                  {
                     col_tmp[k] = col_tmp[k] - 1;
                  }
               }
            }
         }
      }
      save_ctr = count;
      SCIPfreeBufferArray(scip, &found);
   }


   SCIPfreeBufferArray(scip, &fixed_vars);
   SCIPfreeBufferArray(scip, &fixed_values);



   len_of_arrays = count + problemdata->get_lp_nnz();
   DSDPCALLOC2(matind, int, len_of_arrays + 1, &info);
   INFO_TO_SCIPCALL(info);
   DSDPCALLOC2(nnz, double, len_of_arrays + 1, &info);
   INFO_TO_SCIPCALL(info);
   DSDPCALLOC2(block, int, len_of_arrays + 1, &info);
   INFO_TO_SCIPCALL(info);
   DSDPCALLOC2(constraint, int, len_of_arrays + 1, &info);
   INFO_TO_SCIPCALL(info);

   int sdp_count = 0;
   //insert sdp-block
   for (int i = 0; i < count; ++i)
   {
      if (vals_tmp[i] != 0.0) //In vals_tmp I wrote explicit 0.0, so its ok to compare with it
      {
         (*matind)[sdp_count] = (row_tmp[i] * (row_tmp[i] + 1) / 2 + col_tmp[i]);//matind_tmp[i];
         assert ( (*matind)[sdp_count] >= 0);
         (*block)[sdp_count] = block_tmp[i];
         (*constraint)[sdp_count] = constraint_tmp[i];
         (*nnz)[sdp_count] = vals_tmp[i];
         sdp_count++;
      }
   }
   const int* for_mat = problemdata->get_for_matind();
   const int* cons = problemdata->get_for_constraint();
   const double* values = problemdata->get_for_vals();


   count = sdp_count;
   //   insert lp-block
   for (int i = 0; i < problemdata->get_for_matind_size(); ++i)
   {
      if (values[i] != 0.0) //In values I wrote explicit 0.0, so its ok to compare with it
      {
         (*matind)[count] = for_mat[i];

         (*block)[count] = nsdpcones + 1;
         (*constraint)[count] = cons[i];
         (*nnz)[count] = values[i];
         count++;
      }
   }


   len_of_arrays = count;

   (*block)[len_of_arrays] = nsdpcones + 2;
   (*constraint)[len_of_arrays] = (varmapper->get_sdp_nvars() + 2);
   (*matind)[len_of_arrays] = 10000000;
   (*nnz)[len_of_arrays] = 0.0;
   qusort(*block, *constraint, *matind, *nnz, 0, len_of_arrays - 1);


   SCIPfreeBufferArray(scip, &vals_tmp);
   SCIPfreeBufferArray(scip, &matind_tmp);
   SCIPfreeBufferArray(scip, &constraint_tmp);
   SCIPfreeBufferArray(scip, &block_tmp);
   SCIPfreeBufferArray(scip, &row_tmp);
   SCIPfreeBufferArray(scip, &col_tmp);
   return SCIP_OKAY;
}

/**constructor*/
DsdpInterface::DsdpInterface(SCIP* scip): scip_(scip), dsdp_(NULL), matind_(NULL), nnz_(NULL), ittt_(NULL)
{
}


/**method that puts all the data into the solver object, some parts of this function are copied from dsdp*/
SCIP_RETCODE DsdpInterface::put_data_in(
   SdpProblem* problemdata,        /**<class with problemdata of a node*/
   SdpVarMapper* varmapper         /**<varmapper class object*/ )
{
   int     info        = 0;
   int* blocksizes;
   char* conetypes;
   int* constraint;
   int* block;
   int nblocks = problemdata->get_nsdpcones();
   int n;


   SCIP_CALL(transform_data(scip_, problemdata, &block, &matind_ , &constraint, &nnz_, &blocksizes, &conetypes,varmapper));

   const int nvars = SCIPgetNVars(scip_);
   SCIP_VAR ** vars = SCIPgetVars(scip_);

   int nvars_wo_fixed = varmapper->get_sdp_nvars();


   if (dsdp_)
   {
      //delete and free in DSDP everything existing so far
      dsdp_ = NULL;
   }
   //create dsdp-object
   info = DSDPCreate(nvars_wo_fixed, &(dsdp_));
   INFO_TO_SCIPCALL(info);



   SDPCone dsdp_sdpcone = 0;
   LPCone lpcone = 0;

   int spot, ijnnz, nzmats, np, sdpnmax, stat1;
   int sspot;

   SCIPdebugMessage("  Blocksize(linear): %d\n", blocksizes[nblocks]);
   if (nblocks > 0)
   {
      info = DSDPCreateSDPCone (dsdp_, nblocks, &dsdp_sdpcone);
      INFO_TO_SCIPCALL(info);
   }


   SCIP_CALL (set_objective(dsdp_, vars, nvars, varmapper));


   double* y0;
   DSDPCALLOC2(&y0, double, nvars_wo_fixed, &info);
   INFO_TO_SCIPCALL(info);

   for (int i = 0; i < nvars_wo_fixed; i++)
   {
      y0[i] = 0.0;
   }
   for (int i = 0; i < nvars_wo_fixed; i++)
   {
      info = DSDPSetY0(dsdp_, i + 1, y0[i]);
   }

   spot = 0;
   ijnnz = 0;
   np = 0;
   sdpnmax = 1;
   stat1 = 1;
   // Insert the SDP data

   for (int j = 0; j < nblocks + 1; j++)
   {
      if (conetypes[j] == 'S')
      {

         n = blocksizes[j];
         info = CountNonzeroMatrices(j + 1, block + spot, constraint + spot, &nzmats);
         INFO_TO_SCIPCALL(info);
         info = SDPConeSetBlockSize(dsdp_sdpcone, j, n);
         INFO_TO_SCIPCALL(info);
         info = SDPConeSetSparsity(dsdp_sdpcone, j, nzmats);
         INFO_TO_SCIPCALL(info);
         info = SDPConeSetStorageFormat(dsdp_sdpcone, j, 'P'); //P for packed storage format
         INFO_TO_SCIPCALL(info);
         np += n;
         if (sdpnmax < n)
         {
            sdpnmax = n;
         }
         if (stat1 < nzmats)
         {
            stat1 = nzmats;
         }

         for (int i = 0; i <= nvars_wo_fixed; i++)
         {
            info = GetMarkers(j + 1, i, block + spot, constraint + spot, &ijnnz);
            INFO_TO_SCIPCALL(info);
            if ( ijnnz == 0 )
            {
            }
            else if (CheckForConstantMat((nnz_ + spot), ijnnz, n))
            {
               info = SDPConeSetConstantMat(dsdp_sdpcone, j, i, n, nnz_[spot + 1]);
               INFO_TO_SCIPCALL(info);

               info = SDPConeSetXArray(dsdp_sdpcone, j, n, nnz_ + spot, (n * (n + 1) / 2));
               INFO_TO_SCIPCALL(info);

            }
            else if (ijnnz == (n * (n + 1) / 2) )  // check for dense matrix
            {
               info = SDPConeSetADenseVecMat(dsdp_sdpcone, j, i, n, 1.0, nnz_ + spot, ijnnz);
               INFO_TO_SCIPCALL(info);
            }
            else         // sparse matrix
            {
               info = SDPConeSetASparseVecMat(dsdp_sdpcone, j, i, n, 1.0, 0, matind_ + spot, nnz_ + spot, ijnnz);
            }
            spot += ijnnz;

         }

      }

      else if (conetypes[j] == 'L')
      {
         info = DSDPCreateLPCone(dsdp_, &lpcone);
         INFO_TO_SCIPCALL(info);
         n = blocksizes[j];
         np += n;
         sspot = spot;
         DSDPCALLOC2(&ittt_, int, (nvars_wo_fixed + 2), &info);
         INFO_TO_SCIPCALL(info);
         for (int i = 0; i <= nvars_wo_fixed; i++)
         {
            ittt_[i] = 0;
         }
         for (int i = 0; i <= nvars_wo_fixed; i++)
         {
            info = GetMarkers(j + 1, i, block + spot, constraint + spot, &ijnnz);
            INFO_TO_SCIPCALL(info);
            ittt_[i + 1] = ijnnz;
            spot += ijnnz;
         }
         for (int i = 1; i <= nvars_wo_fixed; i++)
         {
            ittt_[i + 1] += ittt_[i];
         }

         info = LPConeSetData(lpcone, n, ittt_, matind_ + sspot, nnz_ + sspot);
         INFO_TO_SCIPCALL(info);

      }
   }


   int its = (nvars_wo_fixed - 2) / sdpnmax;
   if (np < 100 && its == 0)
   {
      its = 1;
   }
   if (its >= 1)
   {
      its++;
   }
   its = its * its;
   if (nvars_wo_fixed < 2000 && its > 10)
   {
      its = 10;
   }
   if (its > 12)
   {
      its = 12;
   }

   info = DSDPReuseMatrix(dsdp_, its);
   INFO_TO_SCIPCALL(info);

   //WARNING!!! DON'T FREE ittt_, nnz_, matind_ in here, dsdp needs
   //them and uses them somewhere beyond this point. Freeing causes segmentation faults.
   double   dnorm[3];
   info = DSDPGetDataNorms(dsdp_, dnorm);
   INFO_TO_SCIPCALL(info);
   if (dnorm[0] == 0)
   {
      info = DSDPSetR0(dsdp_, np);
      INFO_TO_SCIPCALL(info);
      info = DSDPSetGapTolerance(dsdp_, 1e-3);
      INFO_TO_SCIPCALL(info);
      info = DSDPSetYBounds(dsdp_, -1.0, 1.0);
      INFO_TO_SCIPCALL(info);
   }

   DSDPFREE(&constraint, &info);
   INFO_TO_SCIPCALL(info);

   DSDPFREE(&block, &info);
   INFO_TO_SCIPCALL(info);

   DSDPFREE(&blocksizes, &info);
   INFO_TO_SCIPCALL(info);

   DSDPFREE(&y0, &info);
   INFO_TO_SCIPCALL(info);

   DSDPFREE(&conetypes, &info);
   INFO_TO_SCIPCALL(info);

   return SCIP_OKAY;
}



/**calls dsdps solving routine and analyses the output, some parts of this function are copied from dsdp*/
SCIP_RETCODE DsdpInterface::sdp_solve(
   SdpVarMapper* varmapper,    /**<varmapper class object*/
   char* status,               /**<pointer to store convergenst status of dsdp*/
   int* solutiontype,          /**<pointer to store solutiontype of dsdp*/
   double* sol_for_scip        /**<pointer to store solution in scip format*/)
{
   int      info = 0;
   double   ddobj, ppobj, dobj, pobj;
   DSDPSolutionType pdfeasible;
   pdfeasible = DSDP_INFEASIBLE;
   DSDPTerminationReason reason;

   info = DSDPSetup(dsdp_);

   if (info)
   {
      printf("\nProblem Setting problem.  Likely insufficient memory\n");
      return SCIP_NOMEMORY;
   }

   info = DSDPSolve(dsdp_);

   if (info)
   {
      printf("\nNumerical errors encountered in DSDPSolve(). \n");
      return SCIP_ERROR;
   }

   info = DSDPStopReason(dsdp_, &reason);
   INFO_TO_SCIPCALL(info);

   if (reason != DSDP_INFEASIBLE_START)
   {
      info = DSDPComputeX(dsdp_);//computes the solution and sets the feasibility parameter
      INFO_TO_SCIPCALL(info);
   }

   info = DSDPGetSolutionType(dsdp_, &pdfeasible);
   INFO_TO_SCIPCALL(info);
   info = DSDPGetDDObjective(dsdp_, &ddobj);
   INFO_TO_SCIPCALL(info);
   info = DSDPGetPPObjective(dsdp_, &ppobj);
   INFO_TO_SCIPCALL(info);
   info = DSDPGetDObjective(dsdp_, &dobj);
   INFO_TO_SCIPCALL(info);
   info = DSDPGetPObjective(dsdp_, &pobj);
   INFO_TO_SCIPCALL(info);

   *status = 'c';

   if (reason == DSDP_CONVERGED)
   {
      SCIPdebugMessage("DSDP Converged. \n");
      *status = 's';
      if ( !SCIPisFeasEQ(scip_, ppobj, ddobj) && pdfeasible == DSDP_PDFEASIBLE)
      {
         SCIPdebugMessage("BE CAREFUL!: pobj: %g, dobj: %g \n", ppobj, ddobj);
         //This should never happen
         *status = 'c';
         return SCIP_ERROR;
      }
   }
   else if ( reason == DSDP_UPPERBOUND )
   {
      SCIPdebugMessage("DSDP Terminated Because Dual Objective Exceeded its Bound\n");
      *status = 'c';
   }
   else if ( reason == DSDP_SMALL_STEPS )
   {
      SCIPdebugMessage("DSDP Terminated Due to Small Steps\n");
      *status = 'c';
   }
   else if ( reason == DSDP_MAX_IT)
   {
      SCIPdebugMessage("DSDP Terminated Due Maximum Number of Iterations\n");
      *status = 'c';
   }
   else if ( reason == DSDP_INFEASIBLE_START)
   {
      SCIPdebugMessage("DSDP Terminated Due to Infeasible Starting Point\n");
      *status = 'c';
   }
   else if ( reason == DSDP_INDEFINITE_SCHUR_MATRIX)
   {
      SCIPdebugMessage("DSDP Terminated Due to Indefinite Schur Complement\n");
      *status = 'c';
   }
   else if ( reason == DSDP_NUMERICAL_ERROR)
   {
      SCIPdebugMessage("Another numerical error occurred. Check solution \n");
      *status = 'c'; //solution is checkt automatically in sciptrysol
   }
   else if ( reason == DSDP_USER_TERMINATION)
   {
      SCIPdebugMessage("We stopped the process, dsdp did not do anything\n");
      *status = 'c';
   }
   else
   {
      SCIPdebugMessage("DSDP Finished: this must mean, that the reason flag is continue_iterating\n");
      *status = 'c';
      return SCIP_ERROR;
   }
   if (pdfeasible == DSDP_PDFEASIBLE)
   {
      *solutiontype = 0;
      SCIPdebugMessage("DSDP tells us the solution is feasible\n");
   }
   if (pdfeasible == DSDP_UNBOUNDED )
   {
      *solutiontype = 2;
      SCIPdebugMessage("DSDP Dual Unbounded, Primal Infeasible\n");
   }
   else if ( pdfeasible == DSDP_INFEASIBLE )
   {
      *solutiontype = 1;
      SCIPdebugMessage("DSDP Primal Unbounded, Dual Infeasible\n");
   }
   else if ( pdfeasible == DSDP_PDUNKNOWN )
   {
      *solutiontype = 3;
      SCIPdebugMessage("UNKNOWN: wer kann schon sagen, was hier passiert. Dsdp dazu:  Hmm.  Not clear whether either solution is feasible.\n");
   }


   SCIPdebugMessage("P Objective  : %16.8e     %16.8e\n", ppobj, pobj);
   SCIPdebugMessage("DSDP Solution: %16.8e     %16.8e\n\n", ddobj, dobj);

   double   derr[6];
   info = DSDPGetFinalErrors(dsdp_, derr);
   convert_sol(scip_, dsdp_, sol_for_scip, varmapper);
   return SCIP_OKAY;
}

/**free the memory we allocated for dsdp*/
DsdpInterface::~DsdpInterface()
{
   int info;
   if (dsdp_)
   {
      DSDPDestroy(dsdp_);
   }

   DSDPFREE(&(matind_), &info);
   DSDPFREE(&(nnz_), &info);
   DSDPFREE(&(ittt_), &info);

}

/**also calls dsdps solving routine but after adding a penalty term
 * @note if a first solve was not successful, we start solving again, but with changed
 *objective. Now we use a penalty term. So this is just a feasibility check.
 */
SCIP_RETCODE DsdpInterface::again_with_penalty(
   SCIP_RESULT* result,             /**<pointer to store result*/
   SCIP_Real* lowerbound,           /**<pointer to store lowerbound for scip*/
   SdpVarMapper* varmapper,         /**<varmapper class object*/
   SdpProblem* problemdata          /**<class with problemdata of a node*/
   )
{
   char     status;
   double   bound;
   int solutiontype;

   DSDPUsePenalty(dsdp_, 1);

   SCIP_CALL(put_data_in(problemdata, varmapper));

   if (varmapper->get_intsfixed() != 0)
   {
      //isallfixed = 1 if all the variables are fixed and we don't call DSPD_solve
      *result = SCIP_DIDNOTRUN;
      return SCIP_ERROR;
   }

   //Use a penalty term, set the old objective to zero and solve again
   for (int i = 0; i < varmapper->get_sdp_nvars() ; i++)
   {
      DSDPSetDualObjective(dsdp_, i + 1, 0.0);
   }
   int nvars = SCIPgetNVars(scip_);
   double *sol_for_scip;
   SCIP_ALLOC( BMSallocMemoryArray(&sol_for_scip, nvars));

   SCIP_CALL( sdp_solve(varmapper, &status, &solutiontype, sol_for_scip) );
   SCIPfreeMemoryArray(scip, &sol_for_scip);
   DSDPGetDObjective(dsdp_, &bound);

   if (SCIPisGT(scip_, bound, 0.0))
   {
      printf("if DSDP with penalty converges and the solution is feasible. Then the objective should be negative or zero, but it is > 0.0.");
      return SCIP_ERROR;
   }    //if DSDP with penalty converges and the solution is feasible. Then the objective
   //should be negative or zero. If this is the case, then there exists a feasible
   //solution. We return DIDNOTRUN, because we don't know what else to return.
   //Returning cut off is wrong.

   //if DSDP with a penalty term did not converge or can't produce a feasible
   //solution, then there is no feasible solution. In this case we cut off the node.
   if ( status == 'c' || solutiontype != 0 )
   {
      SCIPdebugMessage("\n DSDP with penalty told us the solution is not feasible ->cutoff\n");
      DSDPUsePenalty(dsdp_, 0);
      *result = SCIP_CUTOFF;
      return  SCIP_OKAY;
   }

   if (status == 's' && solutiontype == 0)
   {
      SCIP_NODE* node = SCIPnodeGetParent(SCIPgetCurrentNode(scip_));
      if (!node)
      {
         *result = SCIP_SUSPENDED;
         DSDPUsePenalty(dsdp_, 0);
         return SCIP_OKAY;
      }
      else
      {
         *lowerbound = SCIPnodeGetLowerbound(node);
         *result = SCIP_SUCCESS;
         SCIP_CALL( SCIPupdateLocalLowerbound(scip_, *lowerbound) );
         DSDPUsePenalty(dsdp_, 0);
         return SCIP_OKAY;
      }
   }

   DSDPUsePenalty(dsdp_, 0);
   return SCIP_OKAY;
}

