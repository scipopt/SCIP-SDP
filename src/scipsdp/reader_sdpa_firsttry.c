/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2020 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2020 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*									     */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* #define SCIP_MORE_DEBUG */
/* #define SCIP_DEBUG */

/**@file   reader_sdpa_firsttry.c
 * @brief  file reader for mixed-integer semidefinite programs in SDPA format
 * @author Tim Schmidt
 * @author Frederic Matter
 *
 * @todo Allow to write varbounds other than -infinity/infinity as linear constraints.
 * @todo Allow to write a transformed problem.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>                      /* for strcmp */

#include "scipsdp/reader_sdpa_firsttry.h"
#include "scipsdp/cons_sdp.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h" /* for SCIPcreateConsIndicatorLinCons */

#undef SCIPABORT
#define SCIPABORT() {}

#define READER_NAME             "sdpareader"
#define READER_DESC             "file reader and writer for MISDPs in sdpa format"
#define READER_EXTENSION        "dat-s"

#define SDPA_CHECK_NONNEG       TRUE      /**< when writing: check linear constraints and move nonnegativity(-positivity)
                                           *  constraints to definition of variables (which are now defined in non-negative
                                           *  orthant) */
                                          /*  TODO: currently doesn't work for ranged rows (which are not created by sdpa
                                           *  reader) */

/* Use SDPA_NAME_FORMAT instead of %s when parsing lines, to avoid buffer overflow. */
#define MACRO_STR_EXPAND(tok) #tok
#define MACRO_STR(tok) MACRO_STR_EXPAND(tok)
#define SDPA_NAME_FORMAT "%" MACRO_STR(SDPA_MAX_NAME) "s"
#define SDPA_MAX_LINE  512      /* Last 3 chars reserved for '\r\n\0' */
#define SDPA_MAX_NAME  512	/* Frage: sind die vielfachen von 2 wirklich nötig? (vorher 512) */


char SDPA_LINE_BUFFER[SDPA_MAX_LINE];
char SDPA_NAME_BUFFER[SDPA_MAX_NAME];
double OBJ_VALUE_BUFFER[SDPA_MAX_LINE];
int BLOCK_SIZE_BUFFER[SDPA_MAX_LINE];
int SDPA_IS_COMMENT = 0;

struct SDPA_Data{
   SCIP_Bool*            sdpblockrank1;      /**< rank-1 information for each SDP block (TRUE = should be rank 1) */
   int                   nsdpblocksrank1;    /**< number of SDP constraints/blocks that should be rank 1 */
   int                   nvars;              /**< number of variables and length of createdvars-array */
   SCIP_VAR**            createdvars;        /**< array of variables created by the SDPA reader */
   int                   nlinconss;          /**< number of constraints and length of createdconss-array */
   SCIP_CONS**           createdconss;       /**< array of constraints created by the SDPA reader */
   int                   nsdpblocks;         /**< number of SDP constraints/blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP blocks */
   int*                  sdpnblocknonz;      /**< number of nonzeros for each SDP block */
   int*                  sdpnblockvars;      /**< number of variables for each SDP block */
   int**                 nvarnonz;           /**< number of nonzeros for each block and each variable */
   SCIP_VAR***           sdpblockvars;       /**< SCIP variables appearing in each block */
   int**                 sdprow;             /**< array of all row indices for each SDP block */
   int**                 sdpcol;             /**< array of all column indices for each SDP block */
   SCIP_Real**           sdpval;             /**< array of all values of SDP nonzeros for each SDP block */
   int***                rowpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int***                colpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   SCIP_Real***          valpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries of sdpconst row/col/val [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */
   int                   nsdpaconstblock;    /**< number of constraint blocks specified by sdpa*/
   int*                  memorysizessdp;     /**< size of memory allocated for each sdp constraint */
   int*                  memorysizescon;     /**< size of memory allocated for each linear constraint */
   int                   locationConBlock;   /**< the index of the linear constraint description in the constraint block */

};

typedef struct SDPA_Data SDPA_DATA;

/*
 * Local methods
 */

/** find next line in the integer description */
static
SCIP_RETCODE fIntgets(
   SCIP_FILE*            pFile,              /**< file to read from */
   SCIP_Longint*         linecount           /**< current linecount */
   )
{
   assert( pFile != NULL );
   assert( linecount != NULL );

   /* Find first non-commentary line */
   while ( SCIPfgets(SDPA_LINE_BUFFER, (int) sizeof(SDPA_LINE_BUFFER), pFile) != NULL )
   {
      if(SDPA_LINE_BUFFER[0] != '\n')
      {
         ++(*linecount);
         return SCIP_OKAY;
      }
   }

   return SCIP_READERROR;
}


/** finds first non-commentary line in given file */
static
SCIP_RETCODE SDPAfgets(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pFile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   int                   readLongLine        /**< zero in case of long line */
   )
{
   assert( pFile != NULL );
   assert( linecount != NULL );

   /* set linebreak */
   SDPA_LINE_BUFFER[SDPA_MAX_LINE-2] = '\0';

   /* find first non-commentary line */
   while ( SCIPfgets(SDPA_LINE_BUFFER, (int) sizeof(SDPA_LINE_BUFFER), pFile) != NULL )
   {

      /*check for integer or rank1 section */
      if ( strncmp(SDPA_LINE_BUFFER, "*INTEGER", 8) == 0 ||strncmp(SDPA_LINE_BUFFER, "*RANK1", 5) == 0)
      {
         ++(*linecount);
         return SCIP_OKAY;
      }

      /* skip line if it is a comment */
      if ( SDPA_LINE_BUFFER[0] != '*' && SDPA_LINE_BUFFER[0] != '"'&& SDPA_LINE_BUFFER[0] != '\n' && SDPA_IS_COMMENT == 0 )
      {
         /* set SDPA_IS_COMMENT if the line ends with a comment */
         if ( ! (strrchr(SDPA_LINE_BUFFER, '*') == NULL) || ! (strrchr(SDPA_LINE_BUFFER, '"') == NULL)
            || ! (strrchr(SDPA_LINE_BUFFER, '=') == NULL) )
         {
            SDPA_IS_COMMENT = 1;
            ++(*linecount);
	  }

         /* set SDPA_IS_COMMENT to zero if the line ends */
         if ( ! (strrchr(SDPA_LINE_BUFFER, '\n') == NULL) )
         {
            if ( SDPA_IS_COMMENT == 0 )
               ++(*linecount);

            SDPA_IS_COMMENT = 0;
         }

         /* if line is too long for our buffer correct the buffer and correct position in file */
         if ( SDPA_LINE_BUFFER[SDPA_MAX_LINE - 2] != '\0' && readLongLine == 0 )
         {
            char* last;

            /* buffer is full; erase last token since it might be incomplete */
            last = strrchr(SDPA_LINE_BUFFER, ' ');

            if ( last == NULL )
            {
               SCIPwarningMessage(scip, "we read %d characters from the file; this might indicate a corrupted input file!",
                  SDPA_MAX_LINE - 2);
               SDPA_LINE_BUFFER[SDPA_MAX_LINE - 2] = '\0';
               SCIPdebugMsg(scip, "the buffer might be corrupted\n");
            }
            else
            {
               SCIPfseek(pFile, -(long) strlen(last), SEEK_CUR);
               SCIPinfoMessage(scip, NULL, "correct buffer, reread the last %ld characters\n", (long) strlen(last));
               *last = '\0';
            }
            return SCIP_OKAY;
         }
         else
            return SCIP_OKAY;
      }

      /* set SDPA_IS_COMMENT if the line ends with a comment */
      if ( (! (strrchr(SDPA_LINE_BUFFER, '*') == NULL) || ! (strrchr(SDPA_LINE_BUFFER, '"') == NULL)
            || ! (strrchr(SDPA_LINE_BUFFER, '=') == NULL)) && SDPA_IS_COMMENT == 0 )
      {
         SDPA_IS_COMMENT = 1;
         ++(*linecount);
      }

      /* set SDPA_IS_COMMENT to zero if the line ends */
      if ( ! (strrchr(SDPA_LINE_BUFFER, '\n') == NULL) )
         SDPA_IS_COMMENT = 0;
      else
         SDPA_LINE_BUFFER[SDPA_MAX_LINE - 2] = '\0';
   }
   return SCIP_READERROR;
}



/** method for reading a line of double values with fixed length */
static
int readLineDouble(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   int                   nvals,              /**< number of values to read */
   SCIP_Real*            values              /**< values that have been read */
   )
{
   int i = 0;
   char* token;
   SCIP_Real val;
   char* nonconstendptr;
   char* rest;

   do
   {
      SCIP_CALL( SDPAfgets(scip, pfile, linecount, 0) );
      rest = SDPA_LINE_BUFFER;
      while(( token = SCIPstrtok(rest, " ", &rest) ))
      {
         SCIPstrToRealValue(token, &val, &nonconstendptr);
         *(values + i) = val;
         i = i + 1;
         if ( i >= nvals )
            break;
      }
   }
   while ( i < nvals );

   assert( i == nvals );

   return i;
}


/** method for reading a line of integer values with fixed length*/
static
int readLineInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   int                   nvals,              /**< number of values to read */
   int*                  values              /**< values that have been read */
   )
{
   int i = 0;
   char* token;
   int val;
   char* nonconstendptr;
   char* rest;

   do
   {
      SCIP_CALL( SDPAfgets(scip, pfile, linecount, 0) );
      rest = SDPA_LINE_BUFFER;
      while(( token = SCIPstrtok(rest, " ", &rest) ))
      {
         SCIPstrToIntValue(token, &val, &nonconstendptr);
         *(values + i) = val;
         i = i + 1;
         if ( i >= nvals )
            break;
      }
   }
   while ( i < nvals );

   assert( i == nvals );

   return i;
}



/** reads the number of variables from given SDPA-file */
static
SCIP_RETCODE SDPAreadNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   int cnt = 0;
   int v;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( SDPAfgets(scip,pfile, linecount, 1) );

   if ( sscanf(SDPA_LINE_BUFFER, "%i", &(data->nvars)) != 1 )
   {
      SCIPerrorMessage("Could not read number of scalar variables in line %" SCIP_LONGINT_FORMAT ".\n", *linecount);
      SCIPABORT();
      return SCIP_READERROR;
   }

   if ( data->nvars < 0 )
   {
      SCIPerrorMessage("Number of scalar variables %d in line %" SCIP_LONGINT_FORMAT " should be non-negative!\n",
         data->nvars, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert( data->nvars >= 0 );

   /* loop through different variable types */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdvars), data->nvars) );

   /* create corresponding variables */
   for (v = 0; v < data->nvars; v++)
   {
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
      assert( snprintfreturn < SCIP_MAXSTRLEN );
#else
      (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
#endif

      SCIP_CALL( SCIPcreateVar(scip, &var, varname, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL,
            NULL, NULL) );

      SCIP_CALL( SCIPaddVar(scip, var) );
      data->createdvars[cnt++] = var; /*lint !e732*//*lint !e747*/

      /* release variable for the reader */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}


/** reads the number of constraint blocks from given SDPA-file */
static
SCIP_RETCODE SDPAreadNBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{
   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( SDPAfgets(scip, pfile, linecount, 1) );

   if ( sscanf(SDPA_LINE_BUFFER, "%i", &(data->nsdpaconstblock)) != 1 )
   {
      SCIPerrorMessage("Could not read number of SDP blocks in line %" SCIP_LONGINT_FORMAT ".\n", *linecount);
      SCIPABORT();
      return SCIP_READERROR;
   }

   if ( data->nsdpaconstblock < 0 )
   {
      SCIPerrorMessage("Number of SDP blocks %d in line %" SCIP_LONGINT_FORMAT " should be non-negative!\n",
         data->nsdpaconstblock, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert( data->nsdpaconstblock >= 0 );
   return SCIP_OKAY;
}


/** reads SDP-constraint sizes and number of linear constraints from given SDPA-file */
static
SCIP_RETCODE SDPAreadBlockSize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   int b;
   int cnt = 0;
   int nsdpblocks = 0;
   int* blockVals;
   int nblocks;
   int* blockValsPsd;

#ifndef NDEBUG
   int snprintfreturn;
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &blockVals, data->nsdpaconstblock) );
   SCIP_CALL( SCIPallocBufferArray(scip, &blockValsPsd, data->nsdpaconstblock) );

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   nblocks = readLineInt(scip, pfile, linecount, data->nsdpaconstblock, blockVals);

   if ( data->nsdpaconstblock != nblocks )
   {
      SCIPerrorMessage("Number of specified blocksizes %d in line %" SCIP_LONGINT_FORMAT
         " does not match number of blocks %d.\n", *linecount, nblocks, data->nsdpaconstblock);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   for (int i = 0; i < nblocks; i++)
   {
   /* if the entry is less than zero it describes the LP blocks */
      if ( *(blockVals + i) < 0 )
      {
         if(data->locationConBlock == -1)
         {
            data->locationConBlock = i;
         }
         else
         {
            SCIPerrorMessage("Only one LP block can be defined in line %" SCIP_LONGINT_FORMAT
               " but at least two blocksizes are negative.\n", *linecount);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }
         data->nlinconss = *(blockVals + i) * - 1;
      }
      else
      {
         if ( *(blockVals + i) == 0 )
         {
         SCIPerrorMessage("Encountered a block size of 0 in line %" SCIP_LONGINT_FORMAT " which is not valid.\n",
            *linecount);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }
         *(blockValsPsd + nsdpblocks) = *(blockVals + i);
         nsdpblocks ++;
      }
   }

   assert( data->locationConBlock < 0 || data->nlinconss > 0 );
   assert( data->locationConBlock >= 0 || data->nlinconss == 0 );

   if ( data->nlinconss < 0 )
   {
      SCIPerrorMessage("Number of linear constraints %d in line %" SCIP_LONGINT_FORMAT " should be non-negative!\n",
         data->nlinconss, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert( data->nlinconss >= 0 );
   assert( nsdpblocks >= 0 );

   data->nsdpblocks = nsdpblocks;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblocksizes), data->nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockrank1), data->nsdpblocks) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdconss), data->nlinconss) );


   for (b = 0; b < nsdpblocks; b++)
   {
      assert( blockValsPsd[b] > 0 );    
      data->sdpblocksizes[b] = *(blockValsPsd + b);

      /* initialize rank-1 information to FALSE, will eventually be changed in SDPAreadRank1 */
      data->sdpblockrank1[b] = FALSE;
   }


   /* create corresponding constraints */
   for (int c = 0; c < data->nlinconss; c++)
   {
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LP_%d", cnt);
      assert( snprintfreturn < SCIP_MAXSTRLEN );
#else
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "linear_%d", cnt);
#endif
      /* linear constraints are specified as 0 <= cons <= SCIPinfinity(scip) */
      SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, 0.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE,
            FALSE, FALSE, FALSE, FALSE, FALSE) );

      SCIP_CALL( SCIPaddCons(scip, cons) );
      data->createdconss[cnt++] = cons;

      /* release constraint for the reader. */
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   assert( cnt == data->nlinconss );

   SCIPfreeBufferArray(scip, &blockValsPsd);
   SCIPfreeBufferArray(scip, &blockVals);

   return SCIP_OKAY;
}

/** reads objective values for scalar variables from given SDPA-file */
static
SCIP_RETCODE SDPAreadObjVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int nzerocoef = 0;
   int v;
   int nValsRead = 0;
   SCIP_Real* objVals;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &objVals, data->nvars ) );

   if ( data->createdvars == NULL)
   {
      SCIPerrorMessage("Number of variables needs to be specified before objective values!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 );

   nValsRead = readLineDouble(scip, pfile, linecount, data->nvars, objVals);

   assert(data->nvars == nValsRead);

   for (v = 0; v < data->nvars; v++)
   {
      if ( SCIPisZero(scip, *(objVals + v)) )
         ++nzerocoef;
      else
         SCIP_CALL( SCIPchgVarObj(scip, data->createdvars[v], *(objVals + v)) );
   }

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Found %d objective coefficients with absolute value less than epsilon = %g.\n",
         nzerocoef, SCIPepsilon(scip));
   }

   SCIPfreeBufferArray(scip, &objVals);
   return SCIP_OKAY;
}


/** reads the SDP-constraint blocks and the linear constraint block */
static
SCIP_RETCODE SDPAreadBlocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data,               /**< data pointer to save the results in */
   const char*           filename   	      /**< name of the file that is currently read*/
   )
{
   SCIP_Real val;
   int** sdpvar;
   int b; /** current block */
   int v; /** current variable */
   int c; /** location of the linear constraint block */
   int row;
   int col;
   int firstindforvar;
   int nextindaftervar;
   int nzerocoef = 0;

   int emptySdpBlocks = 0;
   int emptyConBlocks = 0;

   int* currentEntriesSdp;
   int* currentEntriesCon;
   int* currentEntriesLinCon;
   int** sdprow_local;             /**< array of all row indices for each SDP block */
   int** sdpcol_local;             /**< array of all column indices for each SDP block */
   SCIP_Real** sdpval_local;       /**< array of all values of SDP nonzeros for each SDP block */
                                   /**   number of entries of sdpconst row/col/val [i] */
   int** sdpconstrow_local;        /**< pointers to row-indices for each block */
   int** sdpconstcol_local;        /**< pointers to column-indices for each block */
   SCIP_Real **sdpconstval_local;
   SCIP_VAR* indvar = 0;
   SCIP_Bool infeasible;
   int nindcons = 0;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );
   
   assert( data->nvars >= 0 );

   if ( data->nvars < 0 || data-> createdvars == NULL )
   {
      SCIPerrorMessage("Number of variables needs to be specified before entries of the blocks!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   if ( data->nlinconss > 0 )
      SCIP_CALL( SCIPallocBufferArray(scip, &currentEntriesLinCon, data->nlinconss) );

   if(data->nsdpblocks > 0)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->memorysizessdp), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->memorysizescon), data->nsdpblocks) ); 

      SCIP_CALL( SCIPallocBufferArray(scip, &currentEntriesSdp, data->nsdpblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &currentEntriesCon, data->nsdpblocks) );

      /* set initial memory size*/
      for (b = 0; b < data->nsdpblocks; b++) 
      {
         data->memorysizessdp[b] = 8; 
         data->memorysizescon[b] = 8;
         currentEntriesCon[b] = 0;
         currentEntriesSdp[b] = 0;
      }

      for (c = 0; c < data->nlinconss ; c++)
      {
         currentEntriesLinCon[c] = 0;
      }


      if ( data->nsdpblocks < 0 )
      {
         SCIPerrorMessage("Number of blocks needs to be specified before entries of the blocks!\n");
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      if ( data->sdpblocksizes == NULL ) 
      {
         SCIPerrorMessage("Sizes of the SDP blocks need to be specified before entries of the blocks!\n");
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
      assert( data->nlinconss >= 0 );

      /* initialize sdpnblocknonz with 0 */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblocknonz), data->nsdpblocks) ); 
      for (b = 0; b < data->nsdpblocks; b++)
         data->sdpnblocknonz[b] = 0;

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdprow_local), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpcol_local), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpval_local), data->nsdpblocks) );

      for (b = 0; b < data->nsdpblocks; b++) 
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpvar[b]), data->memorysizessdp[b]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdprow_local[b]), data->memorysizessdp[b]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpcol_local[b]), data->memorysizessdp[b]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpval_local[b]), data->memorysizessdp[b]) );
      }

      /* initialize sdpconstnblocknonz with 0 */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstnblocknonz), data->nsdpblocks) );
      for (b = 0; b < data->nsdpblocks; b++)
         data->sdpconstnblocknonz[b] = 0;

      /* allocate memory (constnnonz for each block, since we do not yet know the distribution) */
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstrow_local), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstcol_local), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstval_local), data->nsdpblocks) );


      for (b = 0; b < data->nsdpblocks; b++) 
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstrow_local[b]), data->memorysizescon[b]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstcol_local[b]), data->memorysizescon[b]) );
         SCIP_CALL( SCIPallocBufferArray(scip, &(sdpconstval_local[b]), data->memorysizescon[b]) );
      }  
   }

   while ( SDPAfgets(scip, pfile, linecount,1) == SCIP_OKAY )
   {
      if ( strncmp(SDPA_LINE_BUFFER, "*INTEGER", 8) == 0 || strncmp(SDPA_LINE_BUFFER, "*RANK1", 5) == 0 )
         break;

      if ( sscanf(SDPA_LINE_BUFFER, "%i %i %i %i %lf", &v, &b, &row, &col, &val) != 5 )
      {
         SCIPerrorMessage("Could not read block entry in line %" SCIP_LONGINT_FORMAT ".\n", *linecount);
         SCIPABORT();
         return SCIP_READERROR;
      }

      /* switch from sdpa counting (starting from 1) to scip counting (starting from 0) */
      v -= 1;
      b -= 1;
      row -= 1;
      col -= 1;

      /* check if this entry belongs to the LP block (FALSE) or to an SDP block (TRUE)*/
      if ( b != data->locationConBlock )
      {
      	 /* check if the LP block was already read and adjust the counter */
         if ( b > data->locationConBlock && data->locationConBlock >= 0 )
            b = b - 1;

         if ( v < - 1 || v >= data->nvars )
         {
            SCIPerrorMessage("Given SDP-coefficient in line %" SCIP_LONGINT_FORMAT
               " for variable %d which does not exist!\n", *linecount, v+1);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         /* check if this entry belongs to the constant part of the SDP block (v = -1) or not (v >= 0) */
         if ( v >= 0 )
         {
            if ( b < 0 || b >= data->nsdpblocks )
            {
               SCIPerrorMessage("Given SDP-coefficient in line %" SCIP_LONGINT_FORMAT
                  " for SDP-constraint %d which does not exist!\n", *linecount, b + 1);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( row < 0 || row >= data->sdpblocksizes[b] )
            {
               SCIPerrorMessage("Row index %d of given SDP coefficient in line %" SCIP_LONGINT_FORMAT
                  " is negative or larger than blocksize %d!\n", row +1, *linecount, data->sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( col < 0 || col >= data->sdpblocksizes[b] )
            {
               SCIPerrorMessage("Column index %d of given SDP coefficient in line %" SCIP_LONGINT_FORMAT
                  " is negative or larger than blocksize %d!\n", col + 1, *linecount, data->sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( SCIPisZero(scip, val) )
               ++nzerocoef;
            else
            {
               currentEntriesSdp[b]++;

               /* if the current memory is not sufficient reallocate*/
               if ( currentEntriesSdp[b] >= data->memorysizessdp[b] )
               {
                  int newsize = SCIPcalcMemGrowSize(scip, data->memorysizessdp[b] + 1);
                  assert( newsize > data->memorysizessdp[b] );
                  assert( newsize > currentEntriesSdp[b] );

                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpvar[b]), newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdprow_local[b]), newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpcol_local[b]), newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpval_local[b]), newsize) );
                  data->memorysizessdp[b] = newsize;
               }

               sdpvar[b][data->sdpnblocknonz[b]] = v;

               /* make sure matrix is in lower triangular form */
               if ( col > row )
               {
                  sdprow_local[b][data->sdpnblocknonz[b]] = col;
                  sdpcol_local[b][data->sdpnblocknonz[b]] = row;
               }
               else
               {
                  sdprow_local[b][data->sdpnblocknonz[b]] = row;
                  sdpcol_local[b][data->sdpnblocknonz[b]] = col;
               }

               sdpval_local[b][data->sdpnblocknonz[b]] = val;
               data->sdpnblocknonz[b]++;
            }
         }
         else /* constant part of SDP block*/
         {
            if ( b < 0 || b >= data->nsdpblocks )
            {
               SCIPerrorMessage("Given constant entry in line %" SCIP_LONGINT_FORMAT
                  " for SDP-constraint %d which does not exist!\n", *linecount, b + 1);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( row < 0 || row >= data->sdpblocksizes[b] )
            {
               SCIPerrorMessage("Row index %d of given constant SDP-entry in line %" SCIP_LONGINT_FORMAT
                  " is negative or larger than blocksize %d!\n", row + 1, *linecount, data->sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( col < 0 || col >= data->sdpblocksizes[b] )
            {
               SCIPerrorMessage("Column index %d of given constant SDP-entry in line %" SCIP_LONGINT_FORMAT
                  " is negative or larger than blocksize %d!\n", col + 1, *linecount, data->sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( SCIPisZero(scip, val) )
               ++nzerocoef;
            else
            {
               currentEntriesCon[b]++;

               /* if the current memory is not sufficient reallocate*/
               if ( currentEntriesCon[b] >= data->memorysizescon[b] )
               {
                  int newsize = SCIPcalcMemGrowSize(scip, data->memorysizescon[b] + 1);
                  assert( newsize > data->memorysizescon[b] );
                  assert( newsize > currentEntriesCon[b] );

                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpconstrow_local[b]), newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpconstcol_local[b]), newsize) );
                  SCIP_CALL( SCIPreallocBufferArray(scip, &(sdpconstval_local[b]), newsize) );

                  data->memorysizescon[b] = newsize;
               }
            }

            /* make sure matrix is in lower triangular form */
            if ( col > row )
            {
               sdpconstrow_local[b][data->sdpconstnblocknonz[b]] = col;
               sdpconstcol_local[b][data->sdpconstnblocknonz[b]] = row;
            }
            else
            {
               sdpconstrow_local[b][data->sdpconstnblocknonz[b]] = row;
               sdpconstcol_local[b][data->sdpconstnblocknonz[b]] = col;
            }
            sdpconstval_local[b][data->sdpconstnblocknonz[b]] = val;
            data->sdpconstnblocknonz[b]++;
         }
      }
      else /* LP block */
      {
         /* check if this entry belongs to the constant part of the LP block (v = -1) or not (v >= 0 || v < -1) the latter for indicator variables  */
         if ( v >= 0 )
         {
            /* linear constraints are specified on the diagonal of the LP block */
            if ( row != col )
            {
               SCIPerrorMessage("Given linear coefficient in line %" SCIP_LONGINT_FORMAT
                  " is not located on the diagonal!\n", *linecount);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            assert( row == col );
            c = row;

            if ( c < 0 || c >= data->nlinconss )
            {
               SCIPerrorMessage("Given linear coefficient in line %" SCIP_LONGINT_FORMAT
                  " for constraint %d which does not exist!\n", *linecount, c + 1);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( v < 0 || v >= data->nvars )
            {
               SCIPerrorMessage("Given linear coefficient in line %" SCIP_LONGINT_FORMAT
                  " for variable %d which does not exist!\n", *linecount, v + 1);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if ( SCIPisZero(scip, val) )
            {
               ++nzerocoef;
            }
            else
            {
               SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdvars[v],val) );/*lint !e732*//*lint !e747*/
               currentEntriesLinCon[c]++;
            }         
         }
         else /* constant part or indicator constraint*/
         {
            if( v < -1 )  /* indicator constraint*/
            {
               SCIP_CONS* indcons;
               SCIP_VAR* slackvar = 0;
               char indconsname[SCIP_MAXSTRLEN];
               char slackvarname[SCIP_MAXSTRLEN];
               char linearconsname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
               int snprintfreturn;
#endif

               v = -v - 2;  /* adjust variable index to be positive */

#ifndef NDEBUG
               snprintfreturn = SCIPsnprintf(indconsname, SCIP_MAXSTRLEN, "cons_indicator_%d", nindcons);
               assert( snprintfreturn < SCIP_MAXSTRLEN);
               snprintfreturn = SCIPsnprintf(slackvarname, SCIP_MAXSTRLEN, "indslack_%s", indconsname);
               assert( snprintfreturn < SCIP_MAXSTRLEN);
               snprintfreturn = SCIPsnprintf(linearconsname, SCIP_MAXSTRLEN, "indlin_%s", indconsname);
               assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
               (void) SCIPsnprintf(indconsname, SCIP_MAXSTRLEN, "cons_indicator_%d", nindcons);
               (void) SCIPsnprintf(slackvarname, SCIP_MAXSTRLEN, "indslack_%s", indconsname);
               (void) SCIPsnprintf(linearconsname, SCIP_MAXSTRLEN, "indlin_%s", indconsname);
#endif

               /* create slack variable and add it to the constraint */
               SCIP_CALL( SCIPcreateVar(scip, &slackvar, slackvarname, 0.0, SCIPinfinity(scip), 0.0, 
                  SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, 0, 0, 0, 0, 0)); 
               SCIP_CALL( SCIPaddVar(scip, slackvar) ); 
               
               SCIP_CALL( SCIPaddCoefLinear(scip,data->createdconss[c] , slackvar, +1.0) );/*lint !e732*//*lint !e747*/

               /* change name of the corresponding linear constraint */
               SCIP_CALL( SCIPchgConsName(scip, data->createdconss[c], linearconsname) );
               
               indvar= data->createdvars[v];	
               SCIP_CALL( SCIPchgVarLbGlobal(scip, indvar, 0.0) );
               SCIP_CALL( SCIPchgVarUbGlobal(scip, indvar, 1.0) );
               SCIP_CALL( SCIPchgVarType(scip, indvar, SCIP_VARTYPE_BINARY, &infeasible) );
               SCIP_CALL( SCIPcreateConsIndicatorLinCons( scip, &indcons, indconsname, indvar,data->createdconss[c], slackvar,
                           TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) ); 
               SCIP_CALL( SCIPaddCons(scip, indcons) );
               
               /* release both conss and the slackvar */
               SCIP_CALL( SCIPreleaseCons(scip, &indcons) );
               SCIP_CALL( SCIPreleaseVar(scip, &slackvar) );
               
               nindcons++;
            }
            else /* constant part */
            {
               c = row;

               if ( c < 0 || c >= data->nlinconss )
               {
                  SCIPerrorMessage("Given constant part in line %" SCIP_LONGINT_FORMAT
                     " for scalar constraint %d which does not exist!\n", *linecount, c + 1);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
                  ++nzerocoef;
               else
               {
                  assert( ! SCIPisInfinity(scip, - SCIPgetLhsLinear(scip, data->createdconss[c])) );
                  assert( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, data->createdconss[c])) );

                  /* All linear constraints are greater or equal constraints */
                  SCIP_CALL( SCIPchgLhsLinear(scip, data->createdconss[c], val) );
               }
            }
         }
      }
   }

   for (b = 0; b < data->nsdpblocks; b++)
   {
      if ( currentEntriesSdp[b] == 0 )
      {
         SCIPerrorMessage("SDP block number %d does not contain any nonzero entries!\n", b + 1);
         emptySdpBlocks++;
      }
   }

   if (emptySdpBlocks > 0)
   {
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   for (c = 0; c < data->nlinconss; c++)
   {
      if(currentEntriesLinCon[c] == 0) 
      {
         SCIPerrorMessage("Linear constraint number %d does not contain nonzero entries!\n", c + 1);
         emptyConBlocks++;
      }
   }
   
   if (emptyConBlocks > 0)
   {
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   if(data->nsdpblocks > 0)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval), data->nsdpblocks) );

      for (b = 0; b < data->nsdpblocks; b++) 
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow[b]), data->memorysizessdp[b]) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol[b]), data->memorysizessdp[b]) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval[b]), data->memorysizessdp[b]) );
      }

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstrow), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstcol), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstval), data->nsdpblocks) );

      for (b = 0; b < data->nsdpblocks; b++) 
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstrow[b]), data->memorysizescon[b]) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstcol[b]), data->memorysizescon[b]) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstval[b]), data->memorysizescon[b]) );
      }

      /* write the read blocks from buffer memory to the data object */
      for (b = 0; b < data->nsdpblocks; b++) 
      {
         int nonz = 0;

         for (nonz = 0; nonz < data->sdpnblocknonz[b]; nonz++)
         {
            if ( sdpcol_local[b][nonz] > sdprow_local[b][nonz] )
            {
               data->sdprow[b][nonz] = sdpcol_local[b][nonz];
               data->sdpcol[b][nonz] = sdprow_local[b][nonz];
            }
            else
            {
               data->sdprow[b][nonz] = sdprow_local[b][nonz];
               data->sdpcol[b][nonz] = sdpcol_local[b][nonz];
            }
            data->sdpval[b][nonz] = sdpval_local[b][nonz];
         }

         for (nonz = 0; nonz < data->sdpconstnblocknonz[b]; nonz++)
         {
            if (  sdpconstcol_local[b][nonz] > sdpconstrow_local[b][nonz] )
            {
               data->sdpconstrow[b][nonz] = sdpconstcol_local[b][nonz];
               data->sdpconstcol[b][nonz] = sdpconstrow_local[b][nonz];
            }
            else
            {
               data->sdpconstrow[b][nonz] = sdpconstrow_local[b][nonz];
               data->sdpconstcol[b][nonz] = sdpconstcol_local[b][nonz];
            }
            data->sdpconstval[b][nonz] = sdpconstval_local[b][nonz];
         }
      }

      /* construct pointer arrays */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblockvars), data->nsdpblocks) ); 
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer), data->nsdpblocks) );

      /* sdp blocks as specified in sdpa file */
      for (b = 0; b < data->nsdpblocks; b++)
      {
         /* sort the nonzeroes by non-decreasing variable indices */
         SCIPsortIntIntIntReal(sdpvar[b], data->sdprow[b], data->sdpcol[b], data->sdpval[b], data->sdpnblocknonz[b]);

         /* create the pointer arrays and insert used variables into vars-array */
         nextindaftervar = 0;
         data->sdpnblockvars[b] = 0;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars[b]), data->nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz[b]), data->nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer[b]), data->nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer[b]), data->nvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer[b]), data->nvars) );

         for (v = 0; v < data->nvars; v++)
         {
            SCIP_Bool varused;
            varused = FALSE;

            firstindforvar = nextindaftervar; /* this variable starts where the last one ended */
            data->nvarnonz[b][data->sdpnblockvars[b]] = 0;

            /* get the first index that doesn't belong to this variable */
            while ( nextindaftervar < data->sdpnblocknonz[b] && sdpvar[b][nextindaftervar] == v )
            {
               nextindaftervar++;
               varused = TRUE;
               data->nvarnonz[b][data->sdpnblockvars[b]]++;
            }

            if ( varused )
            {
               /* if the variable is used, add it to the vars array */
               data->sdpblockvars[b][data->sdpnblockvars[b]] = data->createdvars[v];
               /* save a pointer to the first nonzero belonging to this variable */
               data->rowpointer[b][data->sdpnblockvars[b]] = &(data->sdprow[b][firstindforvar]);
               data->colpointer[b][data->sdpnblockvars[b]] = &(data->sdpcol[b][firstindforvar]);
               data->valpointer[b][data->sdpnblockvars[b]] = &(data->sdpval[b][firstindforvar]);
               data->sdpnblockvars[b]++;
            }
         }
         assert( nextindaftervar == data->sdpnblocknonz[b] );
      }

      if ( nzerocoef > 0 )
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "Found %d block coefficients with absolute value less than epsilon = %g.\n", nzerocoef, SCIPepsilon(scip));

      /* free buffer memory */
      for (b = 0; b < data->nsdpblocks; b++)
      {
         SCIPfreeBufferArray(scip, &(sdpconstval_local[b]));
         SCIPfreeBufferArray(scip, &(sdpconstcol_local[b]));
         SCIPfreeBufferArray(scip, &(sdpconstrow_local[b]));
         SCIPfreeBufferArray(scip, &(sdpval_local[b]));
         SCIPfreeBufferArray(scip, &(sdpcol_local[b]));
         SCIPfreeBufferArray(scip, &(sdprow_local[b]));
         SCIPfreeBufferArray(scip, &(sdpvar[b]));
      }
      SCIPfreeBufferArray(scip, &(sdpconstval_local));
      SCIPfreeBufferArray(scip, &(sdpconstcol_local));
      SCIPfreeBufferArray(scip, &(sdpconstrow_local));
      SCIPfreeBufferArray(scip, &(sdpval_local));
      SCIPfreeBufferArray(scip, &(sdpcol_local));
      SCIPfreeBufferArray(scip, &(sdprow_local));
      SCIPfreeBufferArray(scip, &sdpvar);
      SCIPfreeBufferArray(scip, &currentEntriesCon);
      SCIPfreeBufferArray(scip, &currentEntriesSdp);
   }

   if ( data->nlinconss > 0 )
      SCIPfreeBufferArray(scip, &currentEntriesLinCon);

   return SCIP_OKAY;
}




/** reads integrality conditions from given SDPA-file */
static
SCIP_RETCODE SDPAreadInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int v;
   SCIP_Bool infeasible;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdvars == NULL )
   {
      SCIPerrorMessage("Number of variables needs to be specified before integer section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert( data->nvars >= 0 );

   /* read to end of file */
   while ( fIntgets(pfile, linecount) == SCIP_OKAY )
   {
      if (strncmp(SDPA_LINE_BUFFER, "*RANK1", 5) == 0 )
         break;

      if ( sscanf(SDPA_LINE_BUFFER + 1, "%i", &v) != 1 ) /* move the index by one to ignore the first character of the line */
      {
         SCIPerrorMessage("Could not read variable index in line %" SCIP_LONGINT_FORMAT ".\n", *linecount);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      if ( v < 1 || v > data->nvars )
      {
         SCIPerrorMessage("Given integrality in line %" SCIP_LONGINT_FORMAT " for variable %d which does not exist!\n",
            *linecount, v);
         SCIPABORT();
         return SCIP_READERROR;
      }

      v -= 1;
      
      if( SCIPvarGetType(data->createdvars[v]) != SCIP_VARTYPE_BINARY )
      { 
         SCIP_CALL( SCIPchgVarType(scip, data->createdvars[v], SCIP_VARTYPE_INTEGER, &infeasible) );
      }
      if ( infeasible )
      {
         SCIPerrorMessage("Infeasibility detected because of integrality of variable %s!\n",
            SCIPvarGetName(data->createdvars[v]));
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}


/** reads rank1 conditions from given SDPA-file */
static
SCIP_RETCODE SDPAreadRank1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int v;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->sdpblocksizes == NULL )
   {
      SCIPerrorMessage("SDP blocks need to be specified before rank-1 section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert( data->nvars >= 0 );

   /* read to end of file */
   while ( fIntgets(pfile, linecount) == SCIP_OKAY )
   {
      char* ps = SDPA_LINE_BUFFER + 1;
      if ( sscanf(ps, "%i", &v) != 1 )
      {
         SCIPerrorMessage("Could not read SDP block index in line %" SCIP_LONGINT_FORMAT ".\n", *linecount);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      if ( v < 1 || v > data-> nsdpblocks)
      {
         SCIPerrorMessage("Given rank1 in line %" SCIP_LONGINT_FORMAT " for SDP block %d which does not exist!\n",
            *linecount, v);
         SCIPABORT();
         return SCIP_READERROR;
      }

      v -= 1;

      data->sdpblockrank1[v] = TRUE;
      ++data->nsdpblocksrank1;
   }

   return SCIP_OKAY;
}



/** frees all data allocated for the SDPA-data-struct */
static
SCIP_RETCODE SDPAfreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   SDPA_DATA*            data                /**< data pointer to save the results in */
   )
{
   int b = 0;

   assert( scip != NULL );
   assert( data != NULL );
   if(data->nsdpblocks > 0)
   {
      assert( data->nvars > 0 );

      for (b = 0; b < data->nsdpblocks; b++)
      {
         assert( data->memorysizescon[b] > 0);

         SCIPfreeBlockMemoryArrayNull(scip, &(data->valpointer[b]), data->nvars);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->colpointer[b]), data->nvars);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->rowpointer[b]), data->nvars);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->nvarnonz[b]), data->nvars);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpblockvars[b]), data->nvars);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstval[b]), data->memorysizescon[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstcol[b]), data->memorysizescon[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstrow[b]), data->memorysizescon[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpval[b]), data->memorysizessdp[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpcol[b]), data->memorysizessdp[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdprow[b]), data->memorysizessdp[b]);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &data->valpointer, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->colpointer, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->rowpointer, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->nvarnonz, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblockvars, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblockvars, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstval, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstcol, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstrow, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpval, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpcol, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdprow, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstnblocknonz, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblocknonz, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->memorysizescon, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->memorysizessdp, data->nsdpblocks);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &data->createdconss, data->nlinconss);
   SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblockrank1, data->nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblocksizes, data->nsdpblocks);
   SCIPfreeBlockMemoryArrayNull(scip, &data->createdvars, data->nvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopySdpa)
{  /*lint --e{715,818}*/
   assert( scip != NULL );

   SCIP_CALL( SCIPincludeReaderSdpa(scip) );

   return SCIP_OKAY;
}


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSdpa)
{  /*lint --e{715,818}*/
   SCIP_FILE* scipfile;
   SCIP_Longint linecount = 0;
   SDPA_DATA* data;
   int b;

   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Reading file %s ...\n", filename);

   scipfile = SCIPfopen(filename, "r");

   if ( ! scipfile )
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBuffer(scip, &data) );
   data->nsdpblocks = -1;
   data->nsdpblocksrank1 = 0;
   data->nlinconss = 0;
   data->nvars = -1;

   data->sdpblockrank1 = NULL;
   data->sdpblocksizes = NULL;
   data->sdpnblocknonz = NULL;
   data->sdpnblockvars = NULL;
   data->sdpblockvars = NULL;
   data->sdprow = NULL;
   data->sdpcol = NULL;
   data->sdpval = NULL;
   data->sdpconstnblocknonz = NULL;
   data->sdpconstrow = NULL;
   data->sdpconstcol = NULL;
   data->sdpconstval = NULL;
   data->locationConBlock = -1;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* read the file */
   SCIPdebugMsg(scip, "Reading number of variables\n");
   SCIP_CALL( SDPAreadNVars(scip, scipfile, &linecount, data) );

   SCIPdebugMsg(scip, "Reading number of blocks\n");
   SCIP_CALL( SDPAreadNBlocks(scip, scipfile, &linecount, data) );

   SCIPdebugMsg(scip, "Reading blocksizes\n");
   SCIP_CALL( SDPAreadBlockSize(scip, scipfile, &linecount, data) );

   SCIPdebugMsg(scip, "Reading objective values\n");
   SCIP_CALL( SDPAreadObjVals(scip, scipfile, &linecount, data) );

   SCIPdebugMsg(scip, "Reading block entries\n");
   SCIP_CALL( SDPAreadBlocks(scip, scipfile, &linecount, data, filename) );

   if ( strncmp(SDPA_LINE_BUFFER, "*INTEGER", 8) == 0 )
   {
      SCIPdebugMsg(scip, "Reading integer section\n");
      SCIP_CALL( SDPAreadInt(scip, scipfile, &linecount, data) );
   }

   SCIPdebugMsg(scip, "Reading rank1 section\n");
   SCIP_CALL( SDPAreadRank1(scip, scipfile, &linecount, data) );


   /* close the file (and make sure SCIPfclose returns 0) */
   if ( SCIPfclose(scipfile) )
      return SCIP_READERROR;

#ifdef SCIP_MORE_DEBUG
   for (b = 0; b < SCIPgetNConss(scip); b++)
   {
      SCIP_CALL( SCIPprintCons(scip, SCIPgetConss(scip)[b], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
   }
#endif

   /* create SDP-constraints */
   for (b = 0; b < data->nsdpblocks; b++)
   {
      SCIP_CONS *sdpcons;
      char sdpconname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
      int snprintfreturn;
#endif
      assert( data->sdpblocksizes[b] > 0 );
      assert( (data->sdpnblockvars[b] > 0 && data->sdpnblocknonz[b] > 0) || (data->sdpconstnblocknonz[b] > 0) );
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
      assert( snprintfreturn < SCIP_MAXSTRLEN );
#else
      (void) SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
#endif

      /* special treatment of case without constant PSD blocks */
      if ( data->sdpconstnblocknonz == NULL )
      {
         if ( ! data->sdpblockrank1[b] )
         {
            SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
                  data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b],
                  data->valpointer[b], data->sdpblockvars[b], 0, NULL, NULL, NULL, TRUE) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsSdpRank1(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
                  data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b],
                  data->valpointer[b], data->sdpblockvars[b], 0, NULL, NULL, NULL, TRUE) );
         }
      }
      else
      {
         if ( ! data->sdpblockrank1[b] )
         {
            SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data->sdpnblockvars[b],data->sdpnblocknonz[b],
                  data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b],data->rowpointer[b], data->valpointer[b],
                  data->sdpblockvars[b], data->sdpconstnblocknonz[b],data->sdpconstcol[b], data->sdpconstrow[b],
                  data->sdpconstval[b], TRUE) );
         }
         else
         {
            SCIP_CALL( SCIPcreateConsSdpRank1(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
                  data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b], data->valpointer[b],
                  data->sdpblockvars[b], data->sdpconstnblocknonz[b], data->sdpconstcol[b], data->sdpconstrow[b],
                  data->sdpconstval[b], TRUE) );
         }
      }
#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, sdpcons, NULL) );
#endif
      SCIP_CALL( SCIPaddCons(scip, sdpcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &sdpcons) );
   }

   SCIP_CALL( SDPAfreeData(scip, data) );
   SCIPfreeBufferNull(scip, &data);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteSdpa)
{  /*lint --e{715,818}*/
   SCIP_VAR** linvars;
   SCIP_Real* linvals;
   SCIP_VAR** sdpvars;
   SCIP_Real** sdpval;
   SCIP_Real* sdpconstval;
   int** sdpcol;
   int** sdprow;
   int* sdpconstcol;
   int* sdpconstrow;
   int* sdpnvarnonz;
   int* varsenses;
   int* consssenses;
   int nsdpconss = 0;
   int sdpnvars;
   int sdpnnonz;
   int totalsdpnnonz = 0;
   int sdpblocksize;
   int sdparraylength;
   int totalsdpconstnnonz = 0;
   int sdpconstnnonz;
   int consind = 0;
   int linconsind = 0;
   int consmax = 0;
   int c;
   int i;
   int v;
   SCIP_Real val;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int blocks;
   int const_sign = 1;
   int nChangedConstraints = 0;
   int nvarbndslinconss = 0;
   int nlinconss = 0;
   int nrank1sdpblocks = 0;
   int objcoeff = 1;

   assert( scip != NULL );
   assert( result != NULL );
   assert( nvars > 0 );
   assert( nconss >= 0 );

   SCIPdebugMsg(scip, "Writing problem in SDPA format to file.\n");
   *result = SCIP_DIDNOTRUN;

   if ( transformed )
   {
      SCIPerrorMessage("SDPA reader currently only supports writing original problems!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

#ifndef NDEBUG
   for (v = 0; v < nvars; v++)
   {
      assert( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL );
   }
#endif
   /* NVARS */

   /* write number of variables */
   SCIPinfoMessage(scip, file, "%d\n", nvars);

   /* collect different variable senses */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsenses, nvars) );

   for (v = 0; v < nvars; v++)
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbOriginal(vars[v]);
      ub = SCIPvarGetUbOriginal(vars[v]);

      varsenses[v] = 0;
      if ( SCIPisZero(scip, lb) )
      {
         varsenses[v] = 1;
         nvarbndslinconss++;
      }
      else
      {
         if ( ! SCIPisInfinity(scip, -lb) )
         {
            SCIPerrorMessage("Can only handle variables with lower bound 0 or minus infinity.\n");
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }
      }

      if ( SCIPisZero(scip, ub) )
      {
         varsenses[v] = -1;
         nvarbndslinconss++;
      }
      else
      {
         if ( ! SCIPisInfinity(scip, ub) )
         {
            SCIPerrorMessage("Can only handle variables with upper bound 0 or infinity.\n");
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }
      }
   }

   /* collect different constraint senses */
   SCIP_CALL( SCIPallocBufferArray(scip, &consssenses, nconss) );

   /* check if all constraints are either linear or SDP*/
   for (c = 0; c < nconss; c++)
   {
      if ( (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0)
         && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0) &&
          (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0))
      {
         SCIPerrorMessage("SDPA reader currently only supports linear, SDP and SDPrank1 constraints!\n");
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      /* count number of rank1 sdp blocks */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0 )
         ++nrank1sdpblocks;

      /* only check linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") == 0 )
      {

         lhs = SCIPgetLhsLinear(scip, conss[c]);
         rhs = SCIPgetRhsLinear(scip, conss[c]);

         if ( SCIPisEQ(scip, lhs, rhs) )
         {
            assert( ! SCIPisInfinity(scip, -lhs) && ! SCIPisInfinity(scip, rhs) );
            nlinconss += 2;
            consssenses[c] = 0;
         }
         else
         {
            if ( ! SCIPisInfinity(scip, -lhs) && ! SCIPisInfinity(scip, rhs) )
            {
               SCIPerrorMessage("Cannot handle ranged rows.\n");
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }
            else
            {
               assert( SCIPisInfinity(scip, -lhs) || SCIPisInfinity(scip, rhs) );

               if ( ! SCIPisInfinity(scip, -lhs) )
                  consssenses[c] = 1;
               else if ( ! SCIPisInfinity(scip, rhs) )
                  consssenses[c] = -1;

               nlinconss++;
            }
         }
      }
      else
      {
         assert( (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0) ||
            (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0) );

         /* count number of SDP constraints (conshdlrGetNConss doesn't seem to work before transformation) */
         ++nsdpconss;
         ++consmax;

         /* count SDP nonzeros */
         SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz) );
         totalsdpnnonz += sdpnnonz;
         totalsdpconstnnonz += sdpconstnnonz;
      }
   }

   blocks = nsdpconss;
   
   if(blocks > 0 && totalsdpnnonz == 0)
   {
      SCIPerrorMessage("There are %d SDP blocks but no nonzero coefficients. \n", blocks);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   
   if ( nvarbndslinconss + nlinconss > 0 )
      blocks++;

   /* write number of blocks */
   SCIPinfoMessage(scip, file, "%d\n", blocks);

   /* write sizes of the SDP blocks and number of linear constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 && strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1")   )
         continue;

      SCIPinfoMessage(scip, file, "%d ", SCIPconsSdpGetBlocksize(scip, conss[c]));
   }

   if ( nvarbndslinconss + nlinconss > 0 )
      SCIPinfoMessage(scip, file, "-%d \n", nvarbndslinconss + nlinconss);
   else
      SCIPinfoMessage(scip, file, "\n");

   /* write the objective values */
   /* If objsense = maximize, multiply objective values with -1 */
   if ( objsense == SCIP_OBJSENSE_MAXIMIZE)
   {
      objcoeff = -1;
      SCIPinfoMessage(scip, NULL, "WARNING: Transforming original maximization problem to a minimization problem by multiplying all objective coefficients by -1. \n");
   }
      
  
   for (v = 0; v < nvars; v++)
   {
      SCIP_Real obj;

      obj = SCIPvarGetObj(vars[v]);

      if ( ! SCIPisZero(scip, obj) )
         SCIPinfoMessage(scip, file, "%.15g ", obj * objcoeff);
      else
         SCIPinfoMessage(scip, file, "%.15g ", 0.0);
   }
   SCIPinfoMessage(scip, file, "\n");

   
   if( nsdpconss > 0 )
   {
      /* allocate memory for SDPdata */
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpnvarnonz, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpcol, totalsdpnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdprow, totalsdpnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpval, totalsdpnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstcol, totalsdpconstnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstrow, totalsdpconstnnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstval, totalsdpconstnnonz) );

      sdparraylength = totalsdpnnonz;
      sdpconstnnonz = totalsdpconstnnonz;
   }
   
   /* write variable bounds as linear constraints */	
   for (c = 0; c < nvars; c++)
   {
      assert(varsenses[c] == 0 || varsenses[c] == -1 || varsenses[c] == 1 );

      if(varsenses[c] == 0)
         continue;

      if(varsenses[c] == -1)
      {
         ++linconsind;
         SCIPinfoMessage(scip, file, "%d %d %d %d -1.0\n", c + 1, consmax + 1, linconsind, linconsind);
      }
      else
      {
         ++linconsind;
         SCIPinfoMessage(scip, file, "%d %d %d %d 1.0\n", c + 1, consmax + 1, linconsind, linconsind);
      }
   }
  
   /* write SDP constraint blocks */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 ||
         strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0 )
      {

         /* write SDP nonzeros */
         if ( totalsdpnnonz > 0 )
         { 	
            /* coefficient matrices */

            /* initialization for SDPconsSDPGetData-call */
            sdparraylength = totalsdpnnonz;
            sdpconstnnonz = totalsdpconstnnonz;

            SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength,
                  sdpnvarnonz, sdprow, sdpcol, sdpval, sdpvars, &sdpconstnnonz,  sdpconstrow, sdpconstcol,
                  sdpconstval, NULL, NULL, NULL) );

            assert( sdpconstnnonz <= totalsdpconstnnonz );
            assert( sdparraylength <= totalsdpnnonz );

            for (v = 0; v < sdpnvars; v++)
            {
               for (i = 0; i < sdpnvarnonz[v]; i++)
               {
                  int ind;
                  ind = SCIPvarGetProbindex(sdpvars[v]);
                  assert( 0 <= ind && ind < nvars );
                  SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n",  ind + 1, consind + 1, sdprow[v][i]+ 1 , sdpcol[v][i] + 1,
                     sdpval[v][i]);
               } 
            }

            /* constant matrix */

            /* initialization for SDPconsSDPGetData-call */
            sdparraylength = totalsdpnnonz;

            assert( sdpconstnnonz <= totalsdpconstnnonz );
            assert( sdparraylength <= totalsdpnnonz );

            for (i = 0; i < sdpconstnnonz; i++)
            {
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n",0, consind + 1, sdpconstrow[i] + 1, sdpconstcol[i] + 1,
                  sdpconstval[i]);
            }
            consind++;           
         }
      }
      else
      {  
         linconsind++;
      
         lhs = SCIPgetLhsLinear(scip, conss[c]);
         rhs = SCIPgetRhsLinear(scip, conss[c]);
         const_sign = 1;

         /* in case of unconstrained left side and constrained right side swap the inequality by multipling with -1 */
         if ( ! SCIPisInfinity(scip, rhs) && SCIPisInfinity(scip, -lhs) )
         {
      	     const_sign = -1;
      	     nChangedConstraints++;
         }

         linvars = SCIPgetVarsLinear(scip, conss[c]);
         linvals = SCIPgetValsLinear(scip, conss[c]);

         if(!SCIPisEQ(scip, lhs, rhs))
         {

            for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
            {
               i = SCIPvarGetProbindex(linvars[v]);
               assert( 0 <= i && i < nvars );
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", i + 1, consmax + 1, linconsind, linconsind, linvals[v] * const_sign);
            }

            /* write the constant part of the LP block */

            if ( const_sign < 0 )
               val = SCIPgetRhsLinear(scip, conss[c]);
            else
               val = SCIPgetLhsLinear(scip, conss[c]);

            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", 0, consmax + 1, linconsind, linconsind, val * const_sign);
         }
         else  /* write linear constraint block */
         {
   
            for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
            {  
               i = SCIPvarGetProbindex(linvars[v]);
               assert( 0 <= i && i < nvars );
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", i + 1, consmax + 1, linconsind,linconsind, linvals[v] * const_sign);
            }

            /* write the constant part of the LP block */

            if ( const_sign < 0 )
               val = SCIPgetRhsLinear(scip, conss[c]);
            else
               val = SCIPgetLhsLinear(scip, conss[c]);

            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", 0, consmax + 1, linconsind, linconsind, val * const_sign);
         
         
            ++linconsind;  
         
            for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
            {
               i = SCIPvarGetProbindex(linvars[v]);
               assert( 0 <= i && i < nvars );
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", i + 1, consmax + 1, linconsind,linconsind, linvals[v] * const_sign*(-1));
            }

            /* write the constant part of the LP block */

            if ( const_sign < 0 )
               val = SCIPgetRhsLinear(scip, conss[c]);
            else
               val = SCIPgetLhsLinear(scip, conss[c]);

            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", 0, consmax + 1, linconsind, linconsind, val * const_sign*(-1));       
         }   
      }
   }
   
   
   if ( nChangedConstraints > 0 )
   	SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Changed the sign of %d constraints. \n", nChangedConstraints);

   /* write integrality constraints */
   if ( nbinvars + nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "*INTEGER\n");
      for (v = 0; v < nbinvars + nintvars; v++)
      {
         assert( SCIPvarIsIntegral(vars[v]) );
         SCIPinfoMessage(scip, file, "*%d\n", v + 1);
      }
   }

   /* write rank-1 SDP constraints (if existing) */
   if ( nrank1sdpblocks > 0 )
   {
      consind = 0;
      SCIPinfoMessage(scip, file, "*RANK1\n");
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") == 0 )
            continue;

	 if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0)
	    consind++;

	 if(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0)
	 {
            assert( SCIPconsSdpShouldBeRankOne(conss[c]) );
            SCIPinfoMessage(scip, file, "*%d\n", consind);
         }
      }
   }

   if( nsdpconss > 0 )
   {
      SCIPfreeBufferArray(scip, &sdpconstval);
      SCIPfreeBufferArray(scip, &sdpconstrow);
      SCIPfreeBufferArray(scip, &sdpconstcol);
      SCIPfreeBufferArray(scip, &sdpvars);
      SCIPfreeBufferArray(scip, &sdpval);
      SCIPfreeBufferArray(scip, &sdprow);
      SCIPfreeBufferArray(scip, &sdpcol);
      SCIPfreeBufferArray(scip, &sdpnvarnonz);
   }   
   SCIPfreeBufferArray(scip, &consssenses);
   SCIPfreeBufferArray(scip, &varsenses);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the SDPA file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSdpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata = NULL;
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert( reader != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopySdpa) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadSdpa) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteSdpa) );

   return SCIP_OKAY;
}
