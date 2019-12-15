/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2019 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2019 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objreader_sdpaind.cpp
 * @brief  Reader for SDPA-Files with indicator constraints (-var in linear constraint => indicator constraint with var as indicator variable)
 * @author Jakob Schelbert
 * @author Sonja Mars
 * @author Tristan Gally
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "objreader_sdpaind.h"

#include <cassert>                      // for assert
#include <cctype>                       // for isspace
#include <cstdio>                       // for printf
#include <cstdlib>                      // for abs                      /*lint !e10*//*lint !e129*/
#include <istream>                      // for istream, etc
#include <string>                       // for getline, string

#include "BlockMemoryAllocator.h"       // for BlockMemoryAllocator
#include "ScipStreamBuffer.h"           // for ScipStreamBuffer
#include "scipsdp/cons_sdp.h"           // for SCIPcreateConsSdp

#include "scip/cons_linear.h"           // for SCIPaddCoefLinear, etc
#include "scip/cons_indicator.h"        // for SCIPcreateConsIndicatorLinCons
#include "scip/scip.h"                  // for SCIPinfinity, etc

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

/** drop spaces and all brackets that are allowed within the blocks in the sdpa format */
inline void drop_space(std::istream& line)
{
   while ( std::isspace(line.peek()) || line.peek() == '(' || line.peek() == '{' || line.peek() == ')' || line.peek() == '}' || line.peek() == ',' )
   {
      (void) line.ignore(1);/*lint !e747*/
   }
}

/** drops the rest of the line */
inline void drop_rest_line(std::istream& s)
{
   std::string tmp;
   (void) std::getline(s, tmp);
}

/** function for removing comments in between the variable & block definitions */
static
void dropComments(
   std::istream*         file                /**< the file instance that is read */
   )
{
   char fst_col;
   fst_col = (char) (*file).peek();
   while (fst_col == '"' || fst_col == '*')
   {
      drop_rest_line(*file);
      fst_col = (char) (*file).peek();
   }
}

/** drop spaces and all brackets that are allowed within the blocks in the sdpa format, throws an error if it reaches a newline */
static
SCIP_RETCODE dropSpaceNewlineError(
   std::istream&         line                /**< the file instance that is read */
   )
{
   if ( line.peek() == '\n' )
   {
      SCIPerrorMessage("Input File invalid, SDP/LP-block rows need to consist of five entries each, see data_format.txt.\n");
      return SCIP_ERROR;
   }

   while ( std::isspace(line.peek()) || line.peek() == '(' || line.peek() == '{' || line.peek() == ')' || line.peek() == '}' || line.peek() == ',' )
   {
      (void) line.ignore(1);/*lint !e747*/
      if ( (char) line.peek() == '\n' )
      {
         SCIPerrorMessage("Input File invalid, SDP/LP-block rows need to consist of five entries each, see data_format.txt\n");
         return SCIP_ERROR;
      }
   }
   return SCIP_OKAY;
}

/** checks that only spaces, newlines or comments follow in the current line */
static
SCIP_RETCODE checkForLineEnd(
   std::istream&         line                /**< the file instance that is read */
   )
{
   if ( line.peek() == '\n' || line.peek() == '"' || line.peek() == '*' || line.peek() == EOF )
      return SCIP_OKAY;

   while ( std::isspace(line.peek()) || line.peek() == '(' || line.peek() == '{' || line.peek() == ')' || line.peek() == '}' || line.peek() == ',' )
   {
      (void) line.ignore(1);/*lint !e747*/

      if ( line.peek() == '\n' || line.peek() == '"' || line.peek() == '*' || line.peek() == EOF )
         return SCIP_OKAY;
      else if ( std::isspace(line.peek()) || line.peek() == '(' || line.peek() == '{' || line.peek() == ')' || line.peek() == '}' || line.peek() == ',' )
         continue;
      else
      {
         SCIPerrorMessage("Input File invalid, SDP/LP-block rows need to consist of five entries each, see data_format.txt\n");
         return SCIP_ERROR;
      }
   }
   return SCIP_OKAY;
}

/** function to test whether the next character in the input string is a digit (or a minus), if it isn't SCIP aborts with a corresponding error */
static
SCIP_RETCODE testDigit(
   std::istream*         file                /** the file instance that is read */
   )
{
   if ( (! isdigit((*file).peek())) && (! ((*file).peek() == '-')) )
   {
      SCIPerrorMessage("Input File invalid, got character '%c', but only numerals allowed in SDP/LP-block rows, see data_format.txt\n", (*file).peek());
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** function to check whether the given index is within the given bounds, if not an error message for the given string will be thrown */
static
SCIP_RETCODE checkIndex(
   const char*           indexname,          /**< name of the index that will be used in the error message */
   int                   value,              /**< value to check against the upper bound */
   int                   lb,                 /**< lower bound to check against */
   int                   ub                  /**< upper bound to check against */
   )
{
   if ( value < lb )
   {
      SCIPerrorMessage("In an SDP/LP-block-line %s index %d was smaller than %d.\n", indexname, value, lb);
      return SCIP_ERROR;
   }
   if ( value > ub )
   {
      SCIPerrorMessage("In an SDP/LP-block-line %s index %d was larger than given number of %ss %d.\n", indexname, value, indexname, ub);
      return SCIP_ERROR;
   }
   return SCIP_OKAY;
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
   SCIP_DECL_READERREAD(ObjReaderSDPAind::scip_read)
   {/*lint --e{715}*/
      int numvars;                        // Number of variables
      int numblocks;                      // Number of all blocks (SDP + LP)
      int numsdpblocks;                   // Number of SDP-blocks
      int numlpblocks;                    // Number of LP-blocks
      int alllpblocksize;                 // Size of all LP-blocks added
      int* nvarnonz;                      // nblockvarnonz[i] gives the number of nonzeros for variable i

      assert( scip != NULL );
      assert( filename != NULL );
      assert( result != NULL );

      *result = SCIP_DIDNOTRUN;

      std::vector<int, BlockMemoryAllocator<int> > blockpattern =
         std::vector<int, BlockMemoryAllocator<int> >(BlockMemoryAllocator<int>(scip));      // Vector with the sizes of all blocks
      std::vector<SCIP_Real> object;         // Objectivevector
      std::vector<SDPBlock> blockstruct;     // Blockstructure
      LPBlock LPData;                        // LP Data
      std::vector<bool> blockislp;           // Is the block an LP block?
      std::vector<int> intvars;              // Indices of integer variables
      std::vector<int> lp_block_num;
      std::vector<int> lp_block_size;
      int new_row_index;
      bool lp_block_already_done;

      SCIP_FILE* scip_file = SCIPfopen(filename, "r");
      if ( ! scip_file )
         return SCIP_READERROR;

      // setup buffer
      ScipStreamBuffer scip_buffer(scip, scip_file, true);

      // setup our stream from the new buffer
      std::istream file(&scip_buffer);

      if ( ! file )
         return SCIP_READERROR;
      file.clear();

      dropComments(&file);

      // read numvar
      drop_space(file);
      file >> numvars;
      if ( numvars < 0 )
      {
         SCIPerrorMessage("Number of variables is negative!\n");
         return SCIP_READERROR;
      }
      drop_rest_line(file);

      dropComments(&file);

      // read numblocks
      drop_space(file);
      file >> numblocks;
      if ( numblocks < 0 )
      {
         SCIPerrorMessage("Number of blocks is negative!\n");
         return SCIP_READERROR;
      }
      drop_rest_line(file);

      numlpblocks = 0;
      numsdpblocks = 0;
      alllpblocksize = 0;

      dropComments(&file);

      // read block pattern
      blockpattern = std::vector<int, BlockMemoryAllocator<int> >(numblocks, 0, BlockMemoryAllocator<int>(scip));
      blockislp = std::vector<bool>(numblocks, false);/*lint !e747*//*lint !e732*/
      lp_block_num = std::vector<int>(numblocks, 0);
      lp_block_size = std::vector<int>(numblocks, 0);

      for (int j = 0; j < numblocks; ++j)
      {
         drop_space(file);
         file >> blockpattern[j];/*lint !e747*//*lint !e732*/
         if ( blockpattern[j] > 0 )/*lint !e747*//*lint !e732*/
         {
            numsdpblocks++;
            blockstruct.push_back(SDPBlock(blockpattern[j]));/*lint !e747*//*lint !e732*/
         }
         else if ( blockpattern[j] < 0 )/*lint !e747*//*lint !e732*/
         {
            // LP block has a negative coefficient!
            numlpblocks++;
            alllpblocksize += abs(blockpattern[j]);/*lint !e747*//*lint !e732*/
            blockislp[j] = true; /*lint !e747*//*lint !e732*//*lint !e1793*/
            blockstruct.push_back(SDPBlock(0));
            lp_block_num[j] = numlpblocks;/*lint !e747*//*lint !e732*/
            lp_block_size[numlpblocks - 1] = abs(blockpattern[j]);/*lint !e747*//*lint !e732*/
         }
         else
            SCIPwarningMessage(scip, "A blocklength 0 seems a bit odd, don't you think!\n");
      }
      assert(numblocks == numsdpblocks + numlpblocks);

      drop_rest_line(file);
      drop_space(file);
      dropComments(&file);

      // read objective
      object = std::vector<SCIP_Real>(numvars, 0.0);/*lint !e747*//*lint !e732*/
      for (int i = 0; i < numvars; ++i)
      {
         file >> object[i];/*lint !e747*//*lint !e732*/
         drop_space(file);
      }

      SCIPdebugMsg(scip, "Number of variables: %d \n", numvars);
      SCIPdebugMsg(scip, "Number of blocks: %d \n", numblocks);
      SCIPdebugMsg(scip, "Number of SDP-cones: %d\n", numsdpblocks);
      SCIPdebugMsg(scip, "Number of LP-cones: %d\n", numlpblocks);

      // construct blocks

      // construct LP block
      LPData.rows = std::vector<LProw>(alllpblocksize);/*lint !e747*//*lint !e732*/
      LPData.numrows = alllpblocksize;
      SCIPdebugMsg(scip, "Number of LP constraints: %d\n", alllpblocksize);

      std::string commentline;

      // read data
      while ( ! file.eof() )
      {
         if ( file.peek() == '*' ) // comment
         {
            (void) std::getline(file, commentline);
            if ( commentline.find("*INT") == 0 ) // if current line starts with *INT then go to Integer definitions
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
            SCIP_Real val;

            drop_space(file);

            SCIP_CALL( testDigit(&file) );
            file >> var_index;
            SCIP_CALL( checkIndex("variable", var_index, 0, numvars) );
            SCIP_CALL( dropSpaceNewlineError(file) );

            SCIP_CALL( testDigit(&file) );
            file >> block_index;
            SCIP_CALL( checkIndex("block", block_index, 1, numblocks) );
            SCIP_CALL( dropSpaceNewlineError(file) );

            SCIP_CALL( testDigit(&file) );
            file >> row_index;
            SCIP_CALL( checkIndex("row", row_index, 1, (blockislp[block_index - 1] ? LPData.numrows : blockstruct[block_index - 1].blocksize)) );/*lint !e732*//*lint !e747*/
            SCIP_CALL( dropSpaceNewlineError(file) );

            SCIP_CALL( testDigit(&file) );
            file >> col_index;
            SCIP_CALL( checkIndex("column", col_index, 1, (blockislp[block_index - 1] ? LPData.numrows : blockstruct[block_index - 1].blocksize)) );/*lint !e732*//*lint !e747*/
            SCIP_CALL( dropSpaceNewlineError(file) );

            SCIP_CALL( testDigit(&file) );
            file >> val;
            SCIP_CALL( checkForLineEnd(file) );

            if ( SCIPisEQ(scip, val, 0.0) )
            {
               drop_rest_line(file);
               drop_space(file);
               continue;
            }

            // sdp-block
            if ( ! blockislp[block_index - 1] )/*lint !e732*//*lint !e747*/
            {
               if ( row_index < col_index )
               {
                  int save_row = row_index;
                  row_index = col_index;
                  col_index = save_row;
               }

               if ( var_index == 0 )
               {
                  blockstruct[block_index - 1].constcolumns.push_back(col_index);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].constrows.push_back(row_index);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].constvalues.push_back(val);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].constnum_nonzeros++;/*lint !e732*//*lint !e747*/
               }
               else
               {
                  blockstruct[block_index - 1].columns.push_back(col_index);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].rows.push_back(row_index);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].values.push_back(val);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].variables.push_back(var_index);/*lint !e732*//*lint !e747*/
                  blockstruct[block_index - 1].num_nonzeros++;/*lint !e732*//*lint !e747*/
               }
               SCIPdebugMsg(scip, "SDP entry: block_index: %d, row: %d, col: %d, var: %d, val: %g\n", block_index, row_index, col_index, var_index,val );/*lint !e525*/
            }
            // lp-block
            else if ( blockislp[block_index - 1] )/*lint !e732*//*lint !e747*/
            {
               assert( row_index == col_index );
               if ( lp_block_num[block_index - 1] == 1 )   /*lint !e732*//*lint !e747*/
                  new_row_index = row_index - 1;
               else // we combine all lp blocks to a single one, so we add the total number of rows of earlier blocks to the row index
               {
                  int rowoffset = 0;

                  for (int b = 0; b < lp_block_num[block_index - 1] - 1; b++)  /*lint !e732*//*lint !e747*/
                     rowoffset += lp_block_size[b];  /*lint !e732*//*lint !e747*/

                  new_row_index = rowoffset + row_index - 1;
               }
               LPData.rows[new_row_index].data.push_back(std::make_pair(var_index, val));/*lint !e732*//*lint !e747*/
               SCIPdebugMsg(scip, "LP entry: row: %d, var: %d, val: %g\n", new_row_index, var_index,val );
            }

            drop_rest_line(file);
            drop_space(file);
         }
      }

      // read integer variables
      intvars = std::vector<int>(numvars, 0);

      std::string str;
      (void) std::getline(file, str);
      while ( str[0] == '*' && isdigit(str[1]) )
      {
         // get index of integer variable
         int idx;
         (void) str.erase(0, 1);  /*lint !e747*/
         idx = atoi(str.c_str());
         SCIP_CALL( checkIndex("variable", idx, 1, numvars) );

         // in the SDPA-file the variable numbers start at 1!
         intvars[idx - 1] = 1;  /*lint !e732*//*lint !e747*/
         SCIPdebugMsg(scip, "Variable %d is integer.\n", idx - 1);

         (void) std::getline(file, str);
      }

      /************************/
      /* create empty problem */
      /************************/

      SCIP_CALL( SCIPcreateProb(scip, filename, 0, 0, 0, 0, 0, 0, 0) );

      /*****************/
      /* add variables */
      /*****************/

      std::vector<SCIP_VAR*> VariablesX ( numvars );/*lint !e732*//*lint !e747*/

      for (int i = 0; i < numvars; ++i)
      {
         SCIP_VAR* var;
         char      var_name[SCIP_MAXSTRLEN];
#ifndef NDEBUG
         int snprintfreturn = SCIPsnprintf(var_name, SCIP_MAXSTRLEN, "X_%d", i);
         assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
         (void) SCIPsnprintf(var_name, SCIP_MAXSTRLEN, "X_%d", i);
#endif

         if ( intvars[i] == 1 )/*lint !e732*//*lint !e747*/
         {
            SCIP_CALL( SCIPcreateVar(scip, &var, var_name, -SCIPinfinity(scip), SCIPinfinity(scip), object[i], SCIP_VARTYPE_INTEGER, TRUE, FALSE, 0, 0, 0, 0, 0));/*lint !e732*//*lint !e747*/
         }
         else
         {
            SCIP_CALL( SCIPcreateVar(scip, &var, var_name,  -SCIPinfinity(scip), SCIPinfinity(scip), object[i],
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));/*lint !e732*//*lint !e747*/
         }

         SCIP_CALL( SCIPaddVar(scip, var) );
         VariablesX[i] = var;/*lint !e732*//*lint !e747*/

         /* release variable for the reader. */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );
      }


      /*********************************/
      /* create SDP and LP constraints */
      /*********************************/

      lp_block_already_done = false;
      for (int bindex = 0; bindex < numblocks; ++bindex)
      {
         if ( ! blockislp[bindex] )/*lint !e732*//*lint !e747*/
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

            SCIPdebugMsg(scip, "Begin construction of SDP constraint for block %d.\n", bindex);

            blocksize = blockpattern[bindex];/*lint !e732*//*lint !e747*/
            nnonz = blockstruct[bindex].num_nonzeros;/*lint !e732*//*lint !e747*/
            ind = 0;

            /* allocate memory */
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varind, blockstruct[bindex].num_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &col, blockstruct[bindex].num_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &row, blockstruct[bindex].num_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &val, blockstruct[bindex].num_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &colpointer, numvars) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowpointer, numvars) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &valpointer, numvars) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vars, numvars) );/*lint !e732*//*lint !e530*/

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constcol, blockstruct[bindex].constnum_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constrow, blockstruct[bindex].constnum_nonzeros) );/*lint !e732*//*lint !e530*/
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &constval, blockstruct[bindex].constnum_nonzeros) );/*lint !e732*//*lint !e530*/

            /* allocate memory for nblockvarnonz & initialize it with zero */
            SCIP_CALL(SCIPallocBlockMemoryArray(scip, &nvarnonz, numvars));/*lint !e530*/
            for (int i = 0; i < numvars; i++)
               nvarnonz[i] = 0;

            /* prepare the constant arrays */
            for (k = 0; k < blockstruct[bindex].constnum_nonzeros; ++k)/*lint !e732*//*lint !e747*/
            {
               if ( ! SCIPisZero(scip, blockstruct[bindex].constvalues[k]) )/*lint !e732*//*lint !e747*/
               {
                  constcol[ind] = blockstruct[bindex].constcolumns[k] - 1;/*lint !e732*//*lint !e747*/ /* the sdpa format counts from 1 to blocksize, we want to start from 0 */
                  constrow[ind] = blockstruct[bindex].constrows[k] - 1;/*lint !e732*//*lint !e747*/
                  constval[ind] = blockstruct[bindex].constvalues[k];/*lint !e732*//*lint !e747*/
                  ind++;
               }
            }
            constnnonz = ind;

            /* prepare the non-constant arrays */
            ind = 0;
            for (k = 0; k < nnonz; ++k)
            {
               if ( ! SCIPisZero(scip, blockstruct[bindex].values[k]) )/*lint !e732*//*lint !e747*/
               {
                  varind[ind] = blockstruct[bindex].variables[k] - 1;/*lint !e732*//*lint !e747*/
                  col[ind] = blockstruct[bindex].columns[k] - 1;/*lint !e732*//*lint !e747*/
                  row[ind] = blockstruct[bindex].rows[k] - 1;/*lint !e732*//*lint !e747*/
                  val[ind] = blockstruct[bindex].values[k];/*lint !e732*//*lint !e747*/
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
               while ( nextindaftervar < nnonz && varind[nextindaftervar] == k ) /* get the first index that doesn't belong to this variable */
               {
                  nextindaftervar++;
                  varused = TRUE;
                  nvarnonz[ind]++;
               }
               if (varused)
               {
                  vars[ind] = VariablesX[k];/*lint !e732*//*lint !e747*/ /* if the variable is used, add it to the vars array */
                  colpointer[ind] = &col[firstindforvar]; /* save a pointer to the first nonzero belonging to this variable */
                  rowpointer[ind] = &row[firstindforvar];
                  valpointer[ind] = &val[firstindforvar];
                  ind++;
               }
            }

            assert (nextindaftervar == nnonz);

            /* this was only needed to compute the vars arrays */
            SCIPfreeBlockMemoryArray(scip, &varind, blockstruct[bindex].num_nonzeros);/*lint !e747*/

            nvars = ind;

            SCIP_CONS* sdpcon;
            char       sdpcon_name[SCIP_MAXSTRLEN];
#ifndef NDEBUG
            int snprintfreturn = SCIPsnprintf(sdpcon_name, SCIP_MAXSTRLEN, "SDP-Constraint-%d", bindex);
            assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
            (void) SCIPsnprintf(sdpcon_name, SCIP_MAXSTRLEN, "SDP-Constraint-%d", bindex);
#endif
            SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcon, sdpcon_name, nvars, nnonz, blocksize, nvarnonz, colpointer,
                  rowpointer, valpointer, vars, constnnonz, constcol, constrow, constval, FALSE) );

#ifdef SCIP_MORE_DEBUG
            SCIP_CALL( SCIPprintCons(scip, sdpcon, NULL) );
#endif

            SCIP_CALL( SCIPaddCons(scip, sdpcon) );

            SCIP_CALL( SCIPreleaseCons(scip, &sdpcon) );

            /* free the used arrays */
            SCIPfreeBlockMemoryArray(scip, &nvarnonz, numvars);
            SCIPfreeBlockMemoryArray(scip, &constval, blockstruct[bindex].constnum_nonzeros);/*lint !e747*/
            SCIPfreeBlockMemoryArray(scip, &constrow, blockstruct[bindex].constnum_nonzeros);/*lint !e747*/
            SCIPfreeBlockMemoryArray(scip, &constcol, blockstruct[bindex].constnum_nonzeros);/*lint !e747*/
            SCIPfreeBlockMemoryArray(scip, &vars, numvars);
            SCIPfreeBlockMemoryArray(scip, &valpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &rowpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &colpointer, numvars);
            SCIPfreeBlockMemoryArray(scip, &val, blockstruct[bindex].num_nonzeros);/*lint !e747*/
            SCIPfreeBlockMemoryArray(scip, &row, blockstruct[bindex].num_nonzeros);/*lint !e747*/
            SCIPfreeBlockMemoryArray(scip, &col, blockstruct[bindex].num_nonzeros);/*lint !e747*/

            SCIPdebugMsg(scip, "Construction of SDP constraint for block %d completed.\n", bindex);
         }
         else
         {
            // construct lp-block only once
            if ( ! lp_block_already_done )
            {
               int nindconss = 0;

               lp_block_already_done = true;
               SCIPdebugMsg(scip, "Begin construction of LP (block %d).\n", bindex);

               for (int row_i = 0; row_i < LPData.numrows; ++row_i)
               {
                  SCIP_CONS* LPcon;
                  char       LPcon_name[SCIP_MAXSTRLEN];
                  int        indicator = 0;
#ifndef NDEBUG
                  int snprintfreturn = SCIPsnprintf(LPcon_name, SCIP_MAXSTRLEN, "LP-Con-%d", row_i);
                  assert( snprintfreturn < SCIP_MAXSTRLEN );
#else
                  (void) SCIPsnprintf(LPcon_name, SCIP_MAXSTRLEN, "LP-Con-%d", row_i);
#endif

                  // Get right hand side of the constraint
                  SCIP_Real LPlhs = 0.0;
                  for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)/*lint !e732*//*lint !e747*/
                  {
                     if (LPData.rows[row_i].data[var_i].first == 0)/*lint !e732*//*lint !e747*/
                     {
                        LPlhs = LPData.rows[row_i].data[var_i].second;/*lint !e732*//*lint !e747*/
                     }
                  }

                  /* check for indicator constraint */
                  for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)/*lint !e732*//*lint !e747*/
                  {
                     if (LPData.rows[row_i].data[var_i].first < 0)/*lint !e732*//*lint !e747*/
                        indicator++;
                  }

                  if ( indicator == 0 )
                  {
                     /* regular linear constraint */

                     //Create constraint
                     SCIP_CALL( SCIPcreateConsLinear(scip, &LPcon, LPcon_name, 0, 0, 0, LPlhs, SCIPinfinity(scip), TRUE, TRUE,
                           TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

                     SCIP_CALL( SCIPaddCons(scip, LPcon) );

                     //Insert variables into constraint:
                     for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)/*lint !e732*//*lint !e747*/
                     {
                        if (LPData.rows[row_i].data[var_i].first != 0)/*lint !e732*//*lint !e747*/
                        {
                           SCIP_CALL( SCIPaddCoefLinear(scip, LPcon, VariablesX[LPData.rows[row_i].data[var_i].first - 1], LPData.rows[row_i].data[var_i].second) );/*lint !e732*//*lint !e747*/
                        }
                     }
   #ifdef SCIP_MORE_DEBUG
                     SCIP_CALL( SCIPprintCons(scip, LPcon, NULL) );
   #endif
                     SCIP_CALL( SCIPreleaseCons(scip, &LPcon) );
                  }
                  else
                  {
                     SCIP_VAR* indvar = 0;
                     SCIP_VAR* slackvar = 0;
                     SCIP_CONS* indcons;
                     char      cons_name[SCIP_MAXSTRLEN];

#ifndef NDEBUG
                     snprintfreturn = SCIPsnprintf(cons_name, SCIP_MAXSTRLEN, "cons_indicator_%d", nindconss);
                     assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
                     (void) SCIPsnprintf(cons_name, SCIP_MAXSTRLEN, "cons_indicator_%d", nindconss);
#endif

                     /* indicator constraint */

                     /* each indicator constraint should only have a single indicator variable */
                     assert( indicator <= 1 );

                     //Create constraint
                     SCIP_CALL( SCIPcreateConsLinear(scip, &LPcon, LPcon_name, 0, 0, 0, LPlhs, SCIPinfinity(scip), TRUE, TRUE,
                           TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

                     //Insert variables into constraint:
                     for (unsigned int var_i = 0; var_i < LPData.rows[row_i].data.size(); ++var_i)/*lint !e732*//*lint !e747*/
                     {
                        if ( LPData.rows[row_i].data[var_i].first > 0 )/*lint !e732*//*lint !e747*/
                        {
                           SCIP_CALL( SCIPaddCoefLinear(scip, LPcon, VariablesX[LPData.rows[row_i].data[var_i].first - 1], LPData.rows[row_i].data[var_i].second) );/*lint !e732*//*lint !e747*/
                        }
                        else if ( LPData.rows[row_i].data[var_i].first < 0 )/*lint !e732*//*lint !e747*/
                        {
                           char      var_name[SCIP_MAXSTRLEN];
                           SCIP_Bool infeasible;
#ifndef NDEBUG
                           snprintfreturn = SCIPsnprintf(var_name, SCIP_MAXSTRLEN, "s_%d", nindconss);
                           assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
                           (void) SCIPsnprintf(var_name, SCIP_MAXSTRLEN, "s_%d", nindconss);
#endif
                           nindconss++;

                           /* Store indvar, note that we need to make it binary now, since otherwise the indicator constraint handler
                            * will throw an assert later and we cannot verify that it is binary now, so we print a warning. */
                           indvar = VariablesX[-1 * LPData.rows[row_i].data[var_i].first - 1];
#ifdef SCIP_WARN_INDICATOR
                           SCIPwarningMessage(scip, "Changing type of variable %s to binary, since it appears as an indicator variable "
                                 "in constraint %s\n", SCIPvarGetName(indvar), LPcon_name);
#endif
                           assert( SCIPvarIsIntegral( VariablesX[-1 * LPData.rows[row_i].data[var_i].first - 1] ) );
                           SCIP_CALL( SCIPchgVarLbGlobal(scip, indvar, 0.0) );
                           SCIP_CALL( SCIPchgVarUbGlobal(scip, indvar, 1.0) );
                           SCIP_CALL( SCIPchgVarType(scip, indvar, SCIP_VARTYPE_BINARY, &infeasible) );

                           /* create slack variable and add it to the constraint*/
                           SCIP_CALL( SCIPcreateVar(scip, &slackvar, var_name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0));/*lint !e732*//*lint !e747*/

                           SCIP_CALL( SCIPaddVar(scip, slackvar) );

                           /* slack variable has coefficient +1 since we have a >= 0 constraint */
                           SCIP_CALL( SCIPaddCoefLinear(scip, LPcon, slackvar, +1.0) );/*lint !e732*//*lint !e747*/
                        }
                     }

                     /* add constraint */
                     assert( indvar != 0 );
                     assert( slackvar != 0 );
                     SCIP_CALL( SCIPcreateConsIndicatorLinCons( scip, &indcons, cons_name, indvar, LPcon, slackvar,
                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );
                     SCIP_CALL( SCIPaddCons(scip, LPcon) );
                     SCIP_CALL( SCIPaddCons(scip, indcons) );

#ifdef SCIP_MORE_DEBUG
                     SCIP_CALL( SCIPprintCons(scip, indcons, NULL) );
#endif

                     /* release both conss and the slackvar */
                     SCIP_CALL( SCIPreleaseCons(scip, &indcons) );
                     SCIP_CALL( SCIPreleaseVar(scip, &slackvar) );
                     SCIP_CALL( SCIPreleaseCons(scip, &LPcon) );
                  }
               }
            }
         }
      }

      *result = SCIP_SUCCESS;

      return SCIP_OKAY;
   }

}//end of namespace scip
