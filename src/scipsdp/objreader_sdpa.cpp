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

/**@file   objreader_sdpa.cpp
 * @brief  Reader for SDPA-files
 * @author Jakob Schelbert, Sonja Mars, Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "objreader_sdpa.h"

#include <cassert>                     // for assert
#include <cctype>                      // for isspace
#include <cstdio>                      // for printf
#include <cstdlib>                     // for abs
#include <istream>                      // for istream, etc
#include <string>                       // for getline, string
#include "BlockMemoryAllocator.h"       // for BlockMemoryAllocator
#include "ScipStreamBuffer.h"           // for ScipStreamBuffer
#include "SdpCone.h"                    // for SdpCone
#include "objconshdlr_sdp.h"            // for SCIPcreateConsSdp
#include "scip/cons_linear.h"           // for SCIPaddCoefLinear, etc
#include "scip/def.h"                   // for SCIP_CALL, FALSE, TRUE
#include "scip/pub_fileio.h"            // for SCIPfopen, SCIP_FILE
#include "scip/pub_message.h"           // for SCIPdebugMessage
#include "scip/pub_misc.h"              // for SCIPsnprintf
#include "scip/scip.h"                  // for SCIPinfinity, etc
#include "scip/type_cons.h"             // for SCIP_CONS
#include "scip/type_var.h"              // for SCIP_VAR, etc

#include "SdpProblem.h"
#include "SdpVarMapper.h"
#include <fstream>
#include <sstream>

namespace
{

   inline void drop_space(std::istream& line)
   {
      while(std::isspace(line.peek()))
      {
         line.ignore(1);
      }
      return;
   }

   inline void drop_rest_line (std::istream& s)
   {
      std::string tmp;
      std::getline(s, tmp);
      return;
   }

}

namespace scip
{


   /** problem reading method of reader
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
    *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
    *
    *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
    */
   SCIP_RETCODE ObjReaderSDPA::scip_read(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_READER*       reader,             /**< the file reader itself */
      const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
      SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
      )
   {
      *result = SCIP_DIDNOTRUN;

      int numvars;                        //Number of variables
      int numblocks;                      //Number of all blocks (SDP + LP)
      int numsdpblocks;                   //Number of SDP-blocks
      int numlpblocks;                    //Number of LP-blocks
      int alllpblocksize;                 //Size of all LP-blocks added

      std::vector<int, BlockMemoryAllocator<int> > blockpattern =
      std::vector<int, BlockMemoryAllocator<int> >(BlockMemoryAllocator<int>(scip));      //Vector with the sizes of all blocks
      std::vector<double> object;         //Objectivevector
      std::vector<SDPBlock> blockstruct;	//Blockstructure
      LPBlock LPData;                     //LP Data
      std::vector<bool> blockislp;         //Is the block an LP block?
      std::vector<int> intvars;           //Indices of integer variables
      std::vector <int> lp_block_num;
      std::vector <int> lp_block_size;
      int new_row_index;
      bool lp_block_already_done;

      SCIP_FILE* scip_file = SCIPfopen(filename, "r");
      if (!scip_file)
         return SCIP_READERROR;

      // setup buffer
      ScipStreamBuffer scip_buffer(scip, scip_file, true);

      // setup our stream from the new buffer
      std::istream file(&scip_buffer);

      if( !file )
         return SCIP_READERROR;
      file.clear();

      // drop comments
      char fst_col('"');
      fst_col = file.peek();
      while (fst_col == '"' || fst_col == '*')
      {
         drop_rest_line(file);
         fst_col = file.peek();
      }

      //  read numvar
      drop_space(file);
      file >> numvars;
      drop_rest_line(file);

      // read numblocks
      drop_space(file);
      file >> numblocks;

      drop_rest_line(file);

      numlpblocks = 0;
      numsdpblocks = 0;
      alllpblocksize = 0;

      // read block pattern
      blockpattern = std::vector<int, BlockMemoryAllocator<int> >(numblocks, 0, BlockMemoryAllocator<int>(scip));
      blockislp = std::vector<bool>(numblocks, FALSE);
      lp_block_num = std::vector<int>(numblocks, 0);
      lp_block_size = std::vector<int>(numblocks, 0);

      for (int j = 0; j < numblocks; ++j)
      {
         drop_space(file);
         file >> blockpattern[j];
         if (blockpattern[j] > 0)
         {
            numsdpblocks++;
            blockstruct.push_back(SDPBlock(blockpattern[j]));

         }
         else if (blockpattern[j] < 0)
         {
            //LP block has a negative coefficient!
            numlpblocks++;
            alllpblocksize += abs(blockpattern[j]);
            blockislp[j] = TRUE;
            blockstruct.push_back(SDPBlock(0));
            lp_block_num[j] = numlpblocks;
            lp_block_size[numlpblocks] = abs(blockpattern[j]);

         }
         else
         {
            printf("Blocklength 0 seems a bit odd, don't you think!\n");
         }
      }

      assert(numblocks == numsdpblocks + numlpblocks);

      // read objective
      object = std::vector<double>(numvars, 0.0);
      drop_rest_line(file);
      drop_space(file);
      for (int i = 0; i < numvars; ++i)
      {
         file >> object[i];
         drop_space(file);
      }

      SCIPdebugMessage("Number of variables: %d \n", numvars);
      SCIPdebugMessage("Number of blocks: %d \n", numblocks);
      SCIPdebugMessage("Number of SDP- and LP-cones: %d, %d \n", numsdpblocks, numlpblocks);


      // construct blocks

      //construct LP block
      LPData.rows = std::vector<LProw>(alllpblocksize);
      LPData.numrows = alllpblocksize;

      std::vector<int> for_indices;
      std::string commentline;

      // read data
      while(!file.eof())
      {
      	if(file.peek() == '*') // comment
      	{
      		std::getline(file, commentline);
      		if (commentline.find("*INT") == 0) // if current line starts with *INT then go to Integer definitions
      		{
      			drop_space(file); // drop \newline
      			break;
      		}
      		else // usual comment line
      		{
      			drop_rest_line(file);
      			drop_space(file);
      		}
      	}
      	else
      	{
      		int var_index, block_index; // block id
      		int row_index, col_index; // position in matrix
      		double val;
      		drop_space(file);

      		file >> var_index;
      		drop_space(file);
      		file >> block_index;
      		drop_space(file);
      		file >> row_index;
      		drop_space(file);
      		file >> col_index;
      		drop_space(file);
      		file >> val;

      		if (SCIPisEQ(scip, val, 0.0))
      		{
      			drop_rest_line(file);
      			drop_space(file);
      			continue;
      		}

      		//sdp-block
      		if (!blockislp[block_index - 1])
      		{
      			int row_leq_col = FALSE;
      			if (row_index <= col_index)
      			{
      				row_leq_col = TRUE;
      			}
      			else
      			{
      				if (row_leq_col)
      				{
      					return SCIP_ERROR;
      				}
      				int save_row = row_index;
      				row_index = col_index;
      				col_index = save_row;
      			}

      			if (var_index == 0)
      			{
      				blockstruct[block_index - 1].constcolumns.push_back(col_index);
      				blockstruct[block_index - 1].constrows.push_back(row_index);
      				blockstruct[block_index - 1].constvalues.push_back(val);
      				blockstruct[block_index - 1].constnum_nonzeros++;
      			}
      			else
      			{
      				blockstruct[block_index - 1].columns.push_back(col_index);
      				blockstruct[block_index - 1].rows.push_back(row_index);
      				blockstruct[block_index - 1].values.push_back(val);
      				blockstruct[block_index - 1].variables.push_back(var_index);
      				blockstruct[block_index - 1].num_nonzeros++;
      			}
      		}
      		//lp-block
      		else if (blockislp[block_index - 1])
      		{
      			assert(row_index == col_index);
      			new_row_index = (row_index - 1) + (lp_block_num[block_index - 1] - 1) * lp_block_size[lp_block_num[block_index - 1] - 1];
      			LPData.rows[new_row_index].data.push_back(std::make_pair(var_index, val));
      			SCIPdebugMessage("block_index: %d, row: %d, var: %d, val: %g\n", block_index, new_row_index, var_index,val );

      		}

      		drop_rest_line(file);
      		drop_space(file);
      	}
      }

      //read integer variables
      intvars = std::vector<int>(numvars, 0);

      while(file.peek() == '*')
      {
         int index;
         file.ignore(1);
         file >> index;
         //in the SDPA-file the variable numbers start at 1!
         intvars[index - 1] = 1;
         SCIPdebugMessage("Variable %d is integer.\n", index - 1);
         drop_rest_line(file);
         drop_space(file);
      }


      /************************/
      /* create empty problem */
      /************************/

      SCIP_CALL( SCIPcreateProb(scip, filename, 0, 0, 0, 0, 0, 0, 0) );

      /*****************/
      /* add variables */
      /*****************/

      std::vector<SCIP_VAR*> VariablesX ( numvars );

      for ( int i = 0; i < numvars; ++i)
      {
         SCIP_VAR* var;
         char      var_name[255];
         SCIPsnprintf(var_name, 255, "X_%d", i);


         if (intvars[i] == 1)
         {
            SCIP_CALL( SCIPcreateVar(scip, &var, var_name, -SCIPinfinity(scip), SCIPinfinity(scip), object[i], SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));
         }
         else
         {
            SCIP_CALL( SCIPcreateVar(scip, &var, var_name,  -SCIPinfinity(scip), SCIPinfinity(scip), object[i],  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
         }

         SCIP_CALL( SCIPaddVar(scip, var) );
         VariablesX[i] = var;

         /* release variable for the reader. */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }


      /*********************************/
      /* create SDP and LP constraints */
      /*********************************/

      lp_block_already_done = FALSE;
      for (int bindex = 0; bindex < numblocks; ++bindex)
      {
         if (!blockislp[bindex])
         {

            SCIPdebugMessage("Begin construction of SDP constraint for block %d.\n", bindex);

            int blocksize = blockpattern[bindex];
            int nnz = blockstruct[bindex].num_nonzeros;
            int const_nnz = blockstruct[bindex].constnum_nonzeros;

            SCIP_VAR ** vars;
            int * col;
            int * row;
            double* vals;

            int * const_col;
            int * const_row;
            double* const_vals;

            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vars, blockstruct[bindex].num_nonzeros));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &col, blockstruct[bindex].num_nonzeros));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &row, blockstruct[bindex].num_nonzeros));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &vals, blockstruct[bindex].num_nonzeros));

            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &const_col, blockstruct[bindex].constnum_nonzeros));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &const_row, blockstruct[bindex].constnum_nonzeros));
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &const_vals, blockstruct[bindex].constnum_nonzeros));
            int idx = 0;
            for (int k = 0; k < const_nnz; ++k)
            {
               if (blockstruct[bindex].constvalues[k] > 10e-6 || blockstruct[bindex].constvalues[k] < -10e-6)
               {
                  const_col[idx] = blockstruct[bindex].constcolumns[k];
                  const_row[idx] = blockstruct[bindex].constrows[k];
                  const_vals[idx] = blockstruct[bindex].constvalues[k];
                  idx++;
               }
            }
            const_nnz = idx;

            idx = 0;
            for (int k = 0; k < nnz; ++k)
            {
               if (blockstruct[bindex].values[k] > 10e-6 || blockstruct[bindex].values[k] < -10e-6)
               {
                  vars[idx] = VariablesX[blockstruct[bindex].variables[k] - 1];
                  col[idx] = blockstruct[bindex].columns[k];
                  row[idx] = blockstruct[bindex].rows[k];
                  vals[idx] = blockstruct[bindex].values[k];
                  idx++;
               }
            }
            nnz = idx;

            SdpCone sdpcone(scip, blocksize, vars, col, row, vals, nnz, const_col, const_row, const_vals, const_nnz);

            SCIP_CONS* sdpcon;
            char       sdpcon_name[255];
            SCIPsnprintf(sdpcon_name, 255, "SDP-Constraint-%d", bindex);
            SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcon, sdpcon_name, sdpcone) );
            SCIP_CALL( SCIPaddCons(scip, sdpcon) );

            SCIP_CALL( SCIPreleaseCons(scip, &sdpcon) );

            SCIPdebugMessage("Construction of SDP constraint for block %d completed.\n", bindex);

         }
         else
         {
            //construct lp-block only once
            if (!lp_block_already_done)
            {
               lp_block_already_done = TRUE;
               SCIPdebugMessage("Begin construction of LP (block %d).\n", bindex);

               for (int row_i = 0; row_i < LPData.numrows; ++row_i)
               {
                  SCIP_CONS* LPcon;
                  char       LPcon_name[255];
                  SCIPsnprintf(LPcon_name, 255, "LP-Con-%d", row_i);

                  //Get right hand side of the constraint
                  double LPlhs = 0.0;

                  for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)
                  {
                     if (LPData.rows[row_i].data[var_i].first == 0)
                     {
                        LPlhs = LPData.rows[row_i].data[var_i].second;
                     }
                  }

                  //Create constraint
                  SCIP_CALL( SCIPcreateConsLinear(scip, &LPcon, LPcon_name, 0, 0, 0, LPlhs, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

                  SCIP_CALL( SCIPaddCons(scip, LPcon) );

                  //Insert variables into constraint:
                  for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)
                  {
                     if (LPData.rows[row_i].data[var_i].first != 0)
                     {
                        SCIP_CALL( SCIPaddCoefLinear(scip, LPcon, VariablesX[LPData.rows[row_i].data[var_i].first - 1], LPData.rows[row_i].data[var_i].second) );
                     }
                  }
                  SCIP_CALL( SCIPreleaseCons(scip, &LPcon) );
               }
            }
         }
      }

      *result = SCIP_SUCCESS;

      return SCIP_OKAY;
   }

   /** method for writing a problem in a sdpa-file, for example the problem in a branching node with all node specific variable fixing*/
   SCIP_RETCODE write_sdpafile(
      SCIP*             scip,                /**< SCIP data structure */
      SdpProblem*       problemdata,         /**< node specific problem data */
      SdpVarMapper*     varmapper            /**<varmapper class object*/
   )
   {
      int nvars = SCIPgetNVars(scip);
      SCIP_VAR** vars = SCIPgetVars(scip);
      int count_nvars_notfixed = 0;
      for (int i = 0; i < nvars; ++i)
      {
         //if variable is not fixed
         if (varmapper->get_sdp_index(vars[i]) != -1)
         {
            count_nvars_notfixed++;
         }
      }

      SCIP_Longint num_nodes = SCIPgetNNodes (scip);
      int depth = SCIPgetDepth(scip);
      std::string s;
      std::stringstream out;
      out << num_nodes;
      s = out.str();


      std::string filename ("Files/test_" );
      filename.append(s);
      filename.append("_depth_");
      std::stringstream ding;
      ding << depth;
      s = ding.str();
      filename.append(s);
      filename.append(".dat-s");


      std::ofstream fs(filename.c_str());
      //write nvars and nblocks
      fs << count_nvars_notfixed << std::endl;
      //we write the lp data in one lp-block, therefore we have one more block, than sdpcones
      fs << problemdata->get_nsdpcones() + 1  << std::endl;


      SCIP_VAR** fixed_vars;
      int n_fixed_vars = varmapper->get_nfixed();
      SCIP_CALL(SCIPallocBufferArray(scip, &fixed_vars, n_fixed_vars));
      double* fixed_values;
      SCIP_CALL(SCIPallocBufferArray(scip, &fixed_values, n_fixed_vars));
      int* blocksizes;
      SCIP_CALL(SCIPallocBufferArray(scip, &blocksizes, problemdata->get_nsdpcones()));

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

      std::vector <int> row;
      std::vector <int> col;
      std::vector <double> val;
      std::vector <int> block;
      std::vector <int> var;

      int save_ctr = 0;

      //SDP-blocks
      for (int i = 0 ; i < problemdata->get_nsdpcones(); ++i)
      {
         SdpCone* sdpcone = problemdata->get_sdpcone(i);
         //A_0
         if (sdpcone->get_const_nnz() > 0 || varmapper->get_nfixed() > 0)
         {
            for ( SdpCone::RhsIterator it = sdpcone->rhs_begin(fixed_vars, n_fixed_vars, fixed_values); it != sdpcone->rhs_end(); ++it)
            {
               SdpCone::element el = *it;
               if (el.val != 0)
               {
                  var.push_back(0);
                  block.push_back(i + 1);
                  col.push_back(el.col + 1);
                  row.push_back(el.row + 1);
                  val.push_back(-el.val);
               }
            }
         }

         //A_i i!=0
         for ( SdpCone::LhsIterator it = sdpcone->lhs_begin(fixed_vars, n_fixed_vars); it != sdpcone->lhs_end(); ++it)
         {
            SdpCone::element el = *it;
            if (el.val != 0)
            {
               var.push_back(varmapper->get_sdp_index(sdpcone->get_var(el.vidx)) + 1);
               block.push_back(i + 1);
               col.push_back(el.col + 1);
               row.push_back(el.row + 1);
               val.push_back(el.val);
            }
         }

         //we need a array that tells us, if a there is an entry for a specific row
         int* found;
         SCIP_CALL(SCIPallocBufferArray(scip, &found, sdpcone->get_blocksize()));
         for (int k = 0; k < sdpcone->get_blocksize(); ++k)
         {
            found[k] = 0;
         }

         for (int j = 0; j < sdpcone->get_blocksize(); ++j)
         {
            for (unsigned int k = save_ctr ; k < row.size(); ++k)
            {
               if (row[k] == j || col[k] == j)
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


         blocksizes[i] = num_not_deleted + 1;

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
                  for (unsigned int k = save_ctr; k < row.size(); k++)
                  {
                     if (row[k] >= row_and_col_to_del)
                     {
                        row[k] = row[k] - 1;
                     }
                     if (col[k] >= row_and_col_to_del)
                     {
                        col[k] = col[k] - 1;
                     }
                  }
               }
            }
         }

         save_ctr = col.size();
         SCIPfreeBufferArray(scip, &found);
      }

      //write blocksizes
      for (int i = 0; i < problemdata->get_nsdpcones(); ++i)
      {
         fs << blocksizes[i] << " ";
      }

      fs << -problemdata->get_size_lpblock() << std::endl;

      //write objective
      for (int i = 0; i < nvars; ++i)
      {
         //if variable is not fixed
         if (varmapper->get_sdp_index(vars[i]) != -1)
         {
            fs << SCIPvarGetObj(vars[i]) << " ";
         }
      }
      fs << std::endl;

      //now write sdp-data
      for (unsigned int i = 0; i < row.size(); ++i)
      {
         fs << var[i] << " " << block[i] << " " << col[i]+1 << " " << row[i]+1 << " " << val[i] << std::endl;
      }

      SCIPfreeBufferArray(scip, &blocksizes);
      SCIPfreeBufferArray(scip, &fixed_vars);
      SCIPfreeBufferArray(scip, &fixed_values);

      const int* cons = problemdata->get_for_constraint();
      const double* vals = problemdata->get_for_vals();
      const int* matind = problemdata->get_for_matind();

      //write lp-block
      int lp_block_number = problemdata->get_nsdpcones() + 1;

      for (int i = 0;  i < problemdata->get_for_vals_size(); ++i)
      {
         if (vals[i] != 0)
         {
            if (cons[i] == 0)
            {
               fs << cons[i] << " " << lp_block_number << " " << matind[i] + 1 << " " << matind[i] + 1 << " " << -vals[i]<< std::endl;
            }
            else
            {
               fs << cons[i] << " " << lp_block_number << " " << matind[i] + 1 << " " << matind[i] + 1 << " " << vals[i]<< std::endl;
            }
         }
      }
      fs << std::endl;

      fs.close();
      return SCIP_OKAY;
   }
}//end of namespace scip

