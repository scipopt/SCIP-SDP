/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objreader_sdpa.cpp
 * @brief  Reader for SDPA-files
 * @author Jakob Schelbert
 * @author Sonja Mars
 * @author Tristan Gally
 */

//#define SCIP_DEBUG
//#define SCIP_MORE_DEBUG

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "objreader_sdpa.h"

#include <cassert>                     // for assert
#include <cctype>                      // for isspace
#include <cstdio>                      // for printf
#include <cstdlib>                     // for abs                      <- lint says this needs an include guard
#include <istream>                      // for istream, etc
#include <string>                       // for getline, string

#include "BlockMemoryAllocator.h"       // for BlockMemoryAllocator
#include "ScipStreamBuffer.h"           // for ScipStreamBuffer
#include "scipsdp/cons_sdp.h"           // for SCIPcreateConsSdp

#include "scip/cons_linear.h"           // for SCIPaddCoefLinear, etc
#include "scip/scip.h"                  // for SCIPinfinity, etc

namespace
{

/* drop spaces and all brackets that are allowed within the blocks in the sdpa format */
   inline void drop_space(std::istream& line)
   {
      while(std::isspace(line.peek()) || line.peek() == '(' || line.peek() == '{' || line.peek() == ')' || line.peek() == '}' || line.peek() == ',')
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

   /** function for removing comments in between the variable & block definitions */
   static
   SCIP_RETCODE dropComments(
      std::istream*       file                /* the file instance that is read */
      )
   {
      char fst_col('"');
      fst_col = (*file).peek();
      while (fst_col == '"' || fst_col == '*')
      {
         drop_rest_line(*file);
         fst_col = (*file).peek();
      }
      return SCIP_OKAY;
   }


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
      int* nvarnonz;                      // nblockvarnonz[i] gives the number of nonzeros for variable i
      SCIP_Real feastol;                  // used to check for nonzeros

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

      // initialize feastol
      SCIP_CALL( SCIPgetRealParam(scip, "numerics/feastol", &feastol) );

      if( !file )
         return SCIP_READERROR;
      file.clear();

      SCIP_CALL(dropComments(&file));

      //  read numvar
      drop_space(file);
      file >> numvars;
      drop_rest_line(file);

      SCIP_CALL(dropComments(&file));

      // read numblocks
      drop_space(file);
      file >> numblocks;

      drop_rest_line(file);

      numlpblocks = 0;
      numsdpblocks = 0;
      alllpblocksize = 0;

      SCIP_CALL(dropComments(&file));

      // read block pattern
      blockpattern = std::vector<int, BlockMemoryAllocator<int> >(numblocks, 0, BlockMemoryAllocator<int>(scip));
      blockislp = std::vector<bool>(numblocks, false);
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
            blockislp[j] = true;
            blockstruct.push_back(SDPBlock(0));
            lp_block_num[j] = numlpblocks;
            lp_block_size[numlpblocks] = abs(blockpattern[j]);

         }
         else
            printf("Blocklength 0 seems a bit odd, don't you think!\n");
      }

      assert(numblocks == numsdpblocks + numlpblocks);

      drop_rest_line(file);
      drop_space(file);

      SCIP_CALL(dropComments(&file));

      // read objective
      object = std::vector<double>(numvars, 0.0);
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
      			if (row_index < col_index)
      			{
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
               SCIPdebugMessage("SDP entry: block_index: %d, row: %d, col: %d, var: %d, val: %g\n", block_index, row_index, col_index, var_index,val );
      		}
      		//lp-block
      		else if (blockislp[block_index - 1])
      		{
      			assert(row_index == col_index);
      			new_row_index = (row_index - 1) + (lp_block_num[block_index - 1] - 1) * lp_block_size[lp_block_num[block_index - 1] - 1];
      			LPData.rows[new_row_index].data.push_back(std::make_pair(var_index, val));
      			SCIPdebugMessage("LP entry: block_index: %d, row: %d, var: %d, val: %g\n", block_index, new_row_index, var_index,val );
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
            SCIP_CALL( SCIPcreateVar(scip, &var, var_name,  -SCIPinfinity(scip), SCIPinfinity(scip), object[i],
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));
         }

         SCIP_CALL( SCIPaddVar(scip, var) );
         VariablesX[i] = var;

         /* release variable for the reader. */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }


      /*********************************/
      /* create SDP and LP constraints */
      /*********************************/

      lp_block_already_done = false;
      for (int bindex = 0; bindex < numblocks; ++bindex)
      {
         if (!blockislp[bindex])
         {
            int nvars;
            int nnonz;
            int blocksize;
            int* varind; /* this is used to sort the nonzeroes by variable-indices and check which variables are actually included in this block */
            int* col;
            int* row;
            SCIP_Real* val;
            int** colpointer;
            int** rowpointer;
            SCIP_Real** valpointer;
            SCIP_Var** vars;
            int constnnonz;
            int* constcol;
            int* constrow;
            SCIP_Real* constval;
            int k;
            int ind;
            int nextindaftervar;
            int firstindforvar;
            SCIP_Bool varused;

            SCIPdebugMessage("Begin construction of SDP constraint for block %d.\n", bindex);

            blocksize = blockpattern[bindex];
            nnonz = blockstruct[bindex].num_nonzeros;
            ind = 0;

            /* allocate memory */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varind, blockstruct[bindex].num_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &col, blockstruct[bindex].num_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &row, blockstruct[bindex].num_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &val, blockstruct[bindex].num_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colpointer, numvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowpointer, numvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &valpointer, numvars) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, numvars) );

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constcol, blockstruct[bindex].constnum_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constrow, blockstruct[bindex].constnum_nonzeros) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constval, blockstruct[bindex].constnum_nonzeros) );

            /* allocate memory for nblockvarnonz & initialize it with zero */
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &nvarnonz, numvars));
            for (int i = 0; i < numvars; i++)
               nvarnonz[i] = 0;

            /* prepare the constant arrays */
            for (k = 0; k < blockstruct[bindex].constnum_nonzeros; ++k)
            {
               if (blockstruct[bindex].constvalues[k] > feastol || blockstruct[bindex].constvalues[k] < -feastol)
               {
                  constcol[ind] = blockstruct[bindex].constcolumns[k] - 1; /* the sdpa format counts from 1 to blocksize, we want to start from 0 */
                  constrow[ind] = blockstruct[bindex].constrows[k] - 1;
                  constval[ind] = blockstruct[bindex].constvalues[k];
                  ind++;
               }
            }
            constnnonz = ind;

            /* prepare the non-constant arrays */
            ind = 0;
            for (k = 0; k < nnonz; ++k)
            {
               if (blockstruct[bindex].values[k] > feastol || blockstruct[bindex].values[k] < -feastol)
               {
                  varind[ind] = blockstruct[bindex].variables[k] - 1;
                  col[ind] = blockstruct[bindex].columns[k] - 1;
                  row[ind] = blockstruct[bindex].rows[k] - 1;
                  val[ind] = blockstruct[bindex].values[k];
                  ind++;
               }
            }
            nnonz = ind;

            SCIPsortIntIntIntReal(varind, col, row, val, nnonz); /* sort the nonzeroes by non-decreasing variable indices */

            /* create the pointer arrays and insert used variables into vars-array */
            nextindaftervar = 0;
            ind = 0; /* sdp index of the current variable */
            for (k = 0; k < numvars; k++)
            {
               varused = FALSE;
               firstindforvar = nextindaftervar; /* this variable starts where the last one ended */
               while (nextindaftervar < nnonz && varind[nextindaftervar] == k) /* get the first index that doesn't belong to this variable */
               {
                  nextindaftervar++;
                  varused = TRUE;
                  nvarnonz[ind]++;
               }
               if (varused)
               {
                  vars[ind] = VariablesX[k]; /* if the variable is used, add it to the vars array */
                  colpointer[ind] = &col[firstindforvar]; /* save a pointer to the first nonzero belonging to this variable */
                  rowpointer[ind] = &row[firstindforvar];
                  valpointer[ind] = &val[firstindforvar];
                  ind++;
               }
            }

            assert (nextindaftervar == nnonz);

            /* this was only needed to compute the vars arrays */
            SCIPfreeBlockMemoryArray(scip, &varind, blockstruct[bindex].num_nonzeros);

            nvars = ind;

            SCIP_CONS* sdpcon;
            char       sdpcon_name[255];
            SCIPsnprintf(sdpcon_name, 255, "SDP-Constraint-%d", bindex);
            SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcon, sdpcon_name, nvars, nnonz, blocksize, nvarnonz, colpointer,
                  rowpointer, valpointer, vars, constnnonz, constcol, constrow, constval) );

#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, sdpcon, NULL) );
#endif

            SCIP_CALL( SCIPaddCons(scip, sdpcon) );

            SCIP_CALL( SCIPreleaseCons(scip, &sdpcon) );

            /* free the used arrays */
            SCIPfreeBlockMemoryArray(scip, &nvarnonz, numvars);
            SCIPfreeBlockMemoryArray(scip, &constval, blockstruct[bindex].constnum_nonzeros);
            SCIPfreeBlockMemoryArray(scip, &constrow, blockstruct[bindex].constnum_nonzeros);
            SCIPfreeBlockMemoryArray(scip, &constcol, blockstruct[bindex].constnum_nonzeros);
            SCIPfreeBlockMemoryArray(scip, &vars, numvars);
            SCIPfreeBlockMemoryArray(scip, &valpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &rowpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &colpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &val, blockstruct[bindex].num_nonzeros);
            SCIPfreeBlockMemoryArray(scip, &row, blockstruct[bindex].num_nonzeros);
            SCIPfreeBlockMemoryArray(scip, &col, blockstruct[bindex].num_nonzeros);

            SCIPdebugMessage("Construction of SDP constraint for block %d completed.\n", bindex);
         }
         else
         {
            //construct lp-block only once
            if (!lp_block_already_done)
            {
               lp_block_already_done = true;
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
                  SCIP_CALL( SCIPcreateConsLinear(scip, &LPcon, LPcon_name, 0, 0, 0, LPlhs, SCIPinfinity(scip), TRUE, TRUE,
                        TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

                  SCIP_CALL( SCIPaddCons(scip, LPcon) );

                  //Insert variables into constraint:
                  for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)
                  {
                     if (LPData.rows[row_i].data[var_i].first != 0)
                     {
                        SCIP_CALL( SCIPaddCoefLinear(scip, LPcon, VariablesX[LPData.rows[row_i].data[var_i].first - 1], LPData.rows[row_i].data[var_i].second) );
                     }
                  }
#ifdef SCIP_MORE_DEBUG
                  SCIP_CALL( SCIPprintCons(scip, LPcon, NULL) );
#endif
                  SCIP_CALL( SCIPreleaseCons(scip, &LPcon) );
               }
            }
         }
      }

      *result = SCIP_SUCCESS;

      return SCIP_OKAY;
   }

}//end of namespace scip
