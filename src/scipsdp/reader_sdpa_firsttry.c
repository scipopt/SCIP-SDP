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
/*                                                 				*/
/*	Test									*/
/*										*/
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

//#define SCIP_MORE_DEBUG 
//#define SCIP_DEBUG

/**@file   reader_cbf.c
 * @brief  file reader for mixed-integer semidefinite programs in SDPA format
 * @author Tristan Gally
 * @author Henrik A. Friberg
 * @author Marc Pfetsch
 * @author Frederic Matter
 *
 * As an extension of the original SDPA format it is possible to specify the constraint that a psd variable and/or a
 * SDP-constraint has rank 1. This is done by using the two new keywords PSDCONRANK1 and PSDVARRANK1 with the following
 * structure:
 *
 * PSDCONRANK1
 * r
 * m_1
 * ...
 * m_r,
 * where \f$r\f$ is the number of SDP-constraints with a rank-1 constraint, and \f$m_1, \dots, m_r\f$ are the indices of the
 * SDP-constraints with a rank-1 constraints.
 *
 * PSDVARRANK1
 * r
 * v_1
 * ...
 * v_r,
 * where \f$r\f$ is the number of psd variables with a rank-1 constraint, and \f$v_1, \dots, v_r\f$ are the indices of the psd
 * variables with a rank-1 constraints.
 *
 * @todo Allow to read SOC constraints in SDPA format.
 * @todo Allow to write varbounds as linear constraints.
 * @todo Allow to write a transformed problem.
 * @todo Allow to write a problem in primal form.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>                      /* for strcmp */

#include "scipsdp/reader_sdpa_firsttry.h"
#include "scipsdp/cons_sdp.h"
#include "scip/cons_linear.h"


#define READER_NAME             "sdpareader"
#define READER_DESC             "file reader and writer for MISDPs in sdpa format"
#define READER_EXTENSION        "dat-s"

#define SDPA_VERSION_NR         2         /**< version number for SDPA format */
#define SDPA_CHECK_NONNEG       TRUE      /**< when writing: check linear constraints and move nonnegativity(-positivity)
                                           *  constraints to definition of variables (which are now defined in non-negative
                                           *  orthant) */
/*  TODO: currently doesn't work for ranged rows (which are not created by sdpa
                                           *  reader) */

/* Use SDPA_NAME_FORMAT instead of %s when parsing lines, to avoid buffer overflow. */
#define MACRO_STR_EXPAND(tok) #tok
#define MACRO_STR(tok) MACRO_STR_EXPAND(tok)
#define SDPA_NAME_FORMAT "%" MACRO_STR(SDPA_MAX_NAME) "s"
#define SDPA_MAX_LINE  512       /* Last 3 chars reserved for '\r\n\0' */
#define SDPA_MAX_NAME  512

char SDPA_LINE_BUFFER[SDPA_MAX_LINE];
char SDPA_LINE_BUFFER2[SDPA_MAX_LINE];
char SDPA_NAME_BUFFER[SDPA_MAX_NAME];
double OBJ_VALUE_BUFFER[SDPA_MAX_LINE];
int BLOCK_SIZE_BUFFER[SDPA_MAX_LINE];

struct SDPA_Data{
   SCIP_Bool *sdpblockrank1;      /**< rank-1 information for each SDP block (TRUE = should be rank 1) */
   int nsdpblocksrank1;    /**< number of SDP constraints/blocks that should be rank 1 */

   int nvars;              /**< number of variables and length of createdvars-array */
   SCIP_VAR **createdvars;        /**< array of variables created by the SDPA reader */
   int nconss;             /**< number of constraints and length of createdconss-array */
   SCIP_CONS **createdconss;       /**< array of constraints created by the SDPA reader */

   int nsdpblocks;         /**< number of SDP constraints/blocks */
   int *sdpblocksizes;      /**< sizes of the SDP blocks */
   int *sdpnblocknonz;      /**< number of nonzeros for each SDP block */
   int *sdpnblockvars;      /**< number of variables for each SDP block */
   int **nvarnonz;           /**< number of nonzeros for each block and each variable */
   SCIP_VAR ***sdpblockvars;       /**< SCIP variables appearing in each block */
   int **sdprow;             /**< array of all row indices for each SDP block */
   int **sdpcol;             /**< array of all column indices for each SDP block */
   SCIP_Real **sdpval;             /**< array of all values of SDP nonzeros for each SDP block */
   int nnonz;              /**< number of nonzeros in blocks */
   int ***rowpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int ***colpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   SCIP_Real ***valpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int *sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries of sdpconst row/col/val [i] */
   int **sdpconstrow;        /**< pointers to row-indices for each block */
   int **sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real **sdpconstval;        /**< pointers to the values of the nonzeros for each block */
   int constnnonz;         /**< number of nonzeros in const block */

   int nsdpaconstblock;    /**< number of constraint blocks specified by sdpa*/

   int *memorysizessdp;
   int *memorysizescon;
   int locationConBlock;

};

typedef struct SDPA_Data SDPA_DATA;

static
int readLineDouble(
        char *LINE_BUFFER,              /** Line as text */
        double *values,
        int number_of_vars           /** Number of vars to read from line  */
){

   const char s[2] = " ";
   char *token;
   char *end;

   /* get the first token */
   token = strtok(LINE_BUFFER, s);
   int i = 0;
   /* walk through other tokens */
   while( token != NULL)
   {
      *(values + i) = strtod(token, &end);
      token = strtok(NULL, s);
      i = i + 1;
      if( i >= number_of_vars )
      {
         break;
      }

   }
   return i;

}

static
int readLineInt(
        char *LINE_BUFFER,              /** Line as text */
        int *values,
        int number_of_vars           /** Number of vars to read from line  */
){

   const char s[2] = " ";
   char *token;
   char *end;

   /* get the first token */
   token = strtok(LINE_BUFFER, s);
   int i = 0;
   /* walk through other tokens */
   while( token != NULL)
   {     
      *(values + i) = strtol(token, &end, 10);//atoi(token);
      token = strtok(NULL, s);
      i = i + 1;
      if( i >= number_of_vars )
      {
         break;
      }
   }
   return i;
}

static
SCIP_RETCODE fIntgets(
        SCIP_FILE *pFile,              /**< file to read from */
        SCIP_Longint *linecount           /**< current linecount */
){
   assert(pFile != NULL);
   assert(linecount != NULL);

   /* Find first non-commentary line */
   while( SCIPfgets(SDPA_LINE_BUFFER, (int) sizeof(SDPA_LINE_BUFFER), pFile) != NULL) 
   {
      ++ (*linecount);
      return SCIP_OKAY;
   }

   return SCIP_READERROR;
}

/*
 * Local methods
 */

/** finds first non-commentary line in given file */
static
SCIP_RETCODE SDPAfgets(
        SCIP_FILE *pFile,              /**< file to read from */
        SCIP_Longint *linecount           /**< current linecount */
){
   assert(pFile != NULL);
   assert(linecount != NULL);

   /* Find first non-commentary line */
   while( SCIPfgets(SDPA_LINE_BUFFER, (int) sizeof(SDPA_LINE_BUFFER), pFile) != NULL) 
   {
      ++ (*linecount);
      if( strncmp(SDPA_LINE_BUFFER, "*INTEGER", 8) == 0 )
      {
         return SCIP_OKAY;
      }
      if( SDPA_LINE_BUFFER[0] != '*' )
      {
         return SCIP_OKAY;
      }
   }

   return SCIP_READERROR;
}


/** finds first non-commentary line in given file */
static
SCIP_RETCODE SDPAfgets2(
        SCIP_FILE *pFile,              /**< file to read from */
        SCIP_Longint *linecount           /**< current linecount */
){
   assert(pFile != NULL);
   assert(linecount != NULL);

   /* Find first non-commentary line */
   while( SCIPfgets(SDPA_LINE_BUFFER2, (int) sizeof(SDPA_LINE_BUFFER2), pFile) !=
          NULL) 
   {
      ++ (*linecount);
      if( strncmp(SDPA_LINE_BUFFER2, "*INTEGER", 8) == 0 )
      {
         return SCIP_OKAY;
      }
      if( SDPA_LINE_BUFFER2[0] != '*' )
      {
         return SCIP_OKAY;
      }
   }
   return SCIP_READERROR;
}


static
SCIP_RETCODE SDPAreadVar(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *pfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data                /**< data pointer to save the results in */
){
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR *var;
   int cnt = 0;

   int v;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert(scip != NULL);
   assert(pfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   SCIP_CALL(SDPAfgets(pfile, linecount));

   if( sscanf(SDPA_LINE_BUFFER, "%i", &(data -> nvars)) != 1 )
   {
      SCIPerrorMessage("Could not read number of scalar variables and conic domains in line %"
      SCIP_LONGINT_FORMAT
      ".\n", *linecount);
      SCIPABORT();
      return SCIP_READERROR;
   }

   if( data -> nvars < 0 ){
      SCIPerrorMessage("Number of scalar variables %d in line %"
      SCIP_LONGINT_FORMAT
      " should be non-negative!\n", data -> nvars, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }


   assert(data -> nvars >= 0);

   /* loop through different variable types */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> createdvars), data -> nvars));

   SCIP_Real lb = - SCIPinfinity(scip);
   SCIP_Real ub = SCIPinfinity(scip);

   /* create corresponding variables */
   for( v = 0; v < data -> nvars; v ++ ){
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
      assert(snprintfreturn < SCIP_MAXSTRLEN);
#else
      (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
#endif

      SCIP_CALL(SCIPcreateVar(scip, &var, varname, lb, ub, 0.0,
                              SCIP_VARTYPE_CONTINUOUS,     
                              TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));

      SCIP_CALL(SCIPaddVar(scip, var));
      data -> createdvars[cnt ++] = var;/*lint !e732*//*lint !e747*/

      /* release variable for the reader */
      SCIP_CALL(SCIPreleaseVar(scip, &var));
   }


   return SCIP_OKAY;
}


/** reads the number and type of constraints from given SDPA-file */
static
SCIP_RETCODE SDPAreadCon(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *pfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data                /**< data pointer to save the results in */
){


#ifndef NDEBUG

#endif

   assert(scip != NULL);
   assert(pfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   SCIP_CALL(SDPAfgets(pfile, linecount));

   if( sscanf(SDPA_LINE_BUFFER, "%i", &(data -> nsdpaconstblock)) != 1 )
   {
      SCIPerrorMessage("Could not read number of scalar constraints and conic domains in line %"
      SCIP_LONGINT_FORMAT
      ".\n", *linecount);
      SCIPABORT();
      return SCIP_READERROR;
   }

   if( data -> nsdpaconstblock < 0 )
   {
      SCIPerrorMessage("Number of scalar constraints %d in line %"
      SCIP_LONGINT_FORMAT
      " should be non-negative!\n", data -> nsdpaconstblock, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   assert(data -> nsdpaconstblock >= 0);
   return SCIP_OKAY;
}


/** reads SDP-constraint sizes from given SDPA-file */
static
SCIP_RETCODE SDPAreadPsdCon(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *pfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data                /**< data pointer to save the results in */
){
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS *cons;
   int b;
   int cnt = 0;
   int ncbfsdpblocks = 0;
   int nlpblocks =0;

   int *blockVals;
   int nblocks;
   int *blockValsLp;
   int *blockValsPsd;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   SCIP_CALL(SCIPallocBufferArray(scip, &blockVals, SDPA_MAX_LINE/2));
   SCIP_CALL(SCIPallocBufferArray(scip, &blockValsLp, SDPA_MAX_LINE/2));
   SCIP_CALL(SCIPallocBufferArray(scip, &blockValsPsd, SDPA_MAX_LINE/2));


   assert(scip != NULL);
   assert(pfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   SCIP_CALL(SDPAfgets(pfile, linecount));

   nblocks = readLineInt(SDPA_LINE_BUFFER, blockVals,
                         data -> nsdpaconstblock); 

   if( data -> nsdpaconstblock != nblocks ){
      SCIPerrorMessage("Number of nblocks returned by readLineInt in line %"
      SCIP_LONGINT_FORMAT
      " should equal data->nsdpaconstblock. Expected: %i got: %i\n", *linecount,data->nsdpaconstblock, data -> nconss);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/ //FRAGE: was heißt dieser Komentar?
   }



   for( int i = 0; i < nblocks; i ++ )
   {
      if( *(blockVals + i) < 0 )
      {
      	
	nlpblocks++;
         data -> locationConBlock = i;
         data -> nconss = *(blockVals + i) * - 1;
         
      }
      else
      {
         *(blockValsPsd + ncbfsdpblocks) = *(blockVals + i);
         ncbfsdpblocks ++;
      }
   }
   if( data -> nconss < 0 )
   {
      SCIPerrorMessage("Number of scalar constraints %d in line %"
      SCIP_LONGINT_FORMAT
      " should be non-negative!\n", data -> nconss, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

	if(nlpblocks>1)
      	{
      	 SCIPerrorMessage("Only one LP can be defined in line %"
      SCIP_LONGINT_FORMAT
      " but %i detected.", *linecount,nlpblocks);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
 
      	}
      	
	
   assert(data -> nconss >= 0);


   if( ncbfsdpblocks < 0 )
   {
      SCIPerrorMessage("Number of SDP-blocks %d in line %"
      SCIP_LONGINT_FORMAT
      " should be non-negative!\n", ncbfsdpblocks, *linecount);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   data -> nsdpblocks = ncbfsdpblocks;


   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpblocksizes), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpblockrank1), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> createdconss), data -> nconss));


   for( b = 0; b < ncbfsdpblocks; b ++ )
   {
      data -> sdpblocksizes[b] = *(blockValsPsd + b);
      if( data -> sdpblocksizes[b] <= 0 )
      {
         SCIPerrorMessage("Size %d of SDP-block %d in line %"
         SCIP_LONGINT_FORMAT
         " should be positive! \n", data -> sdpblocksizes[b], b, *linecount );
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
      /* initialize rank-1 information to FALSE, will eventually be changed in PSDCONRANK1 */
      data -> sdpblockrank1[b] = FALSE;
   }

   SCIP_Real lhs;
   SCIP_Real rhs;

   rhs = SCIPinfinity(
           scip);      
   lhs = 0.0;
   /* create corresponding constraints */
   int c;
   for( c = 0; c < data -> nconss; c ++ )
   {
#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LP_%d", cnt);
      assert(snprintfreturn < SCIP_MAXSTRLEN);
#else
      (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "linear_%d", cnt);
#endif
      SCIP_CALL(SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhs, rhs,
                                     TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      SCIP_CALL(SCIPaddCons(scip, cons));
      data -> createdconss[cnt ++] = cons;
      /* release constraint for the reader. */
      SCIP_CALL(SCIPreleaseCons(scip, &cons));
   }

   if( cnt != data -> nconss )
   {
      SCIPerrorMessage(
              "Total number of constraints for different cone types not equal to total number of constraints!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }


   SCIPfreeBufferArray(scip, &blockValsPsd);
   SCIPfreeBufferArray(scip, &blockValsLp);
   SCIPfreeBufferArray(scip, &blockVals);

   return SCIP_OKAY;
}

/** reads objective values for scalar variables from given SDPA-file */
static
SCIP_RETCODE SDPAreadObjAcoord(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *pfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data                /**< data pointer to save the results in */
){  /*lint --e{818}*/

   int nzerocoef = 0;
   int i;
   int v;
   double *objVals;

   assert(scip != NULL);
   assert(pfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   SCIP_CALL(SCIPallocBufferArray(scip, &objVals, SDPA_MAX_LINE/2));

   if( data -> createdvars == NULL)
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'OBJACOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert(data -> nvars >= 0);

   SCIP_CALL(SDPAfgets(pfile, linecount));

   readLineDouble(SDPA_LINE_BUFFER, objVals, data -> nvars);

   for( v = 0; v <
               data -> nvars; v ++ )              
   {
      if( v < 0 || v >= data -> nvars )
      {
         SCIPerrorMessage("Given objective coefficient in line %"
         SCIP_LONGINT_FORMAT
         " for scalar variable %d which does not exist!\n", *linecount, v);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      if( SCIPisZero(scip, *(objVals + v)))
      {
         ++ nzerocoef;
      }
      else
      {
         SCIP_CALL(SCIPchgVarObj(scip, data -> createdvars[v], *(objVals + v)));
      }
   }

   if( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "OBJACOORD: Found %d coefficients with absolute value less than epsilon = %g.\n", nzerocoef,
                      SCIPepsilon(scip));
   }

   SCIPfreeBufferArray(scip, &objVals);
   return SCIP_OKAY;
}


/** reads nonzero coefficients of SDP-constraints from given SDPA-file */
static
SCIP_RETCODE SDPAreadHcoord(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *scipfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data,                /**< data pointer to save the results in */


        char *filename
){

   SCIP_Real val;
   int **sdpvar;


   int b;
   int v;
   int c; //Constraint für A/Bcoord
   int row;
   int col;
   int firstindforvar;
   int nextindaftervar;
   int nzerocoef = 0;
   int ncbfsdpblocks;
   int nauxnonz = 0;            /* TODO: finde number of nonzeros in each auxiliary sdp block for reformulating matrix variables using
                                 * scalar variables */

   assert(scip != NULL);
   assert(scipfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> memorysizessdp), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> memorysizescon), data -> nsdpblocks));

   SCIP_Longint linecount_copy = 0;
   SCIP_FILE *scip_file_copy = SCIPfopen(filename, "r");

   for( b = 0; b < data -> nsdpblocks; b ++ ) 
   {
      data -> memorysizessdp[b] = 0;
      data -> memorysizescon[b] = 0;
   }

   int s = 0;
   while( SDPAfgets2(scip_file_copy, &linecount_copy) == SCIP_OKAY )
   {
   
      s = s + 1;
      if( strncmp(SDPA_LINE_BUFFER2, "*INTEGER", 8) == 0 )
      {
         break;
      }

      int b2;
      int v2;
      int row2;
      int col2;
      double val2;

      if( s > 4 ){
         if( sscanf(SDPA_LINE_BUFFER2, "%i %i %i %i %lf", &v2, &b2, &row2, &col2, &val2) != 5 )
         { 
            SCIPerrorMessage("Could not read entry of HCOORD in line %"
            SCIP_LONGINT_FORMAT
            ". Wrong number of arguments. \n", linecount_copy);
            SCIPABORT();
            return SCIP_READERROR;
         }

         if( b2 - 1 != data -> locationConBlock ){
            if( b2 - 1 > data -> locationConBlock && data -> locationConBlock >= 0 )
            {
               b2 = b2 - 1;
            }

            if( v2 == 0 )
            {
               data -> memorysizescon[b2 - 1] += 1;
            }
            else
            {
               data -> memorysizessdp[b2 - 1] += 1;
            }
         }


      }
   }


   if( data -> nsdpblocks < 0 )
   {
      SCIPerrorMessage("Need to have 'PSDCON' section before 'HCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert(data -> nvars >= 0);

   if( data -> nvars < 0 )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'HCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   if( data -> createdvars ==
       NULL)  
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'ACOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }


   /* get number of sdp blocks specified by PSDCON (without auxiliary sdp blocks for reformulating matrix variables
     * using scalar variables), save number of nonzeros needed for the auxiliary sdp blocks in nauxnonz */


   ncbfsdpblocks = data -> nsdpblocks;
   nauxnonz = 0;


   if( data -> createdconss == NULL)
   {
      SCIPerrorMessage("Need to have 'CON' section before 'BCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert(data -> nconss >= 0);


   /* initialize sdpnblocknonz with 0 */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpnblocknonz), data -> nsdpblocks));
   for( b = 0; b < data -> nsdpblocks; b ++ )
   {
      data -> sdpnblocknonz[b] = 0;
   }

   /* allocate memory */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &sdpvar, data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdprow), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpcol), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpval), data -> nsdpblocks));

   for( b = 0; b < data -> nsdpblocks; b ++ ) 
   {
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(sdpvar[b]),
                                          data -> memorysizessdp[b]));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdprow[b]), data -> memorysizessdp[b]));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpcol[b]), data -> memorysizessdp[b]));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpval[b]), data -> memorysizessdp[b]));
   }

   /* initialize sdpconstnblocknonz with 0 */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstnblocknonz), data -> nsdpblocks));
   for( b = 0; b < data -> nsdpblocks; b ++ )
   { 
      data -> sdpconstnblocknonz[b] = 0;
   }


   /* allocate memory (constnnonz for each block, since we do not yet know the distribution) */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstrow), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstcol), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstval), data -> nsdpblocks));

   for( b = 0; b < data -> nsdpblocks; b ++ )
   {
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstrow[b]), data -> memorysizescon[b]));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstcol[b]), data -> memorysizescon[b]));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpconstval[b]), data -> memorysizescon[b]));
   }


   while( SDPAfgets(scipfile, linecount) == SCIP_OKAY )
   {
    
      if( strncmp(SDPA_LINE_BUFFER, "*INTEGER", 8) == 0 )
      {
         break;
      }

      if( sscanf(SDPA_LINE_BUFFER, "%i %i %i %i %lf", &v, &b, &row, &col, &val) != 5 )
      { 
         SCIPerrorMessage("Could not read entry of HCOORD in line %"
         SCIP_LONGINT_FORMAT
         ".\n", *linecount);
         SCIPABORT();
         return SCIP_READERROR;
      }
      
      
      
      
      v = v - 1;
      b = b - 1;
      row = row - 1;
      col = col - 1;

      if( b != data -> locationConBlock )
      { 
         if( b > data -> locationConBlock && data -> locationConBlock >= 0 )
         {
            b = b - 1;
         }


	if( v < -1 || v >= data -> nvars )
            {
               SCIPerrorMessage("Given SDP-coefficient in line %"
               SCIP_LONGINT_FORMAT
               " for variable %d which does not exist!\n", *linecount, v);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/}
            


         if( v >= 0 )
         {  //TIM: checke ob Hcoord oder Dcoord

            if( b < 0 || b >= ncbfsdpblocks )
            {
               SCIPerrorMessage("Given SDP-coefficient in line %" 
               SCIP_LONGINT_FORMAT
               " for SDP-constraint %d which does not exist!\n", *linecount, b);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            

            if( row < 0 || row >= data -> sdpblocksizes[b] )
            {
               SCIPerrorMessage("Row index %d of given SDP coefficient in line %"
               SCIP_LONGINT_FORMAT
               " is negative or larger than blocksize %d!\n",
                       row, *linecount, data -> sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( col < 0 || col >= data -> sdpblocksizes[b] )
            {
               SCIPerrorMessage("Column index %d of given SDP coefficient in line %"
               SCIP_LONGINT_FORMAT
               " is negative or larger than blocksize %d!\n",
                       col, *linecount, data -> sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( SCIPisZero(scip, val))
            {
               ++ nzerocoef;
            }else{
               sdpvar[b][data -> sdpnblocknonz[b]] = v;

               /* make sure matrix is in lower triangular form */
               if( col > row )
               {
                  data -> sdprow[b][data -> sdpnblocknonz[b]] = col;
                  data -> sdpcol[b][data -> sdpnblocknonz[b]] = row;
               }
               else
               {
                  data -> sdprow[b][data -> sdpnblocknonz[b]] = row;
                  data -> sdpcol[b][data -> sdpnblocknonz[b]] = col;
               }

               data -> sdpval[b][data -> sdpnblocknonz[b]] = val;
               data -> sdpnblocknonz[b] ++;
            }

         }
         else
         {

            if( b < 0 || b >= data -> nsdpblocks )
            {
               SCIPerrorMessage("Given constant entry in line %" 
               SCIP_LONGINT_FORMAT
               " for SDP-constraint %d which does not exist!\n", *linecount, b);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( row < 0 || row >= data -> sdpblocksizes[b] )
            {
               SCIPerrorMessage("Row index %d of given constant SDP-entry in line %"
               SCIP_LONGINT_FORMAT
               " is negative or larger than blocksize %d!\n",
                       row, *linecount, data -> sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( col < 0 || col >= data -> sdpblocksizes[b] )
            {
               SCIPerrorMessage("Column index %d of given constant SDP-entry in line %"
               SCIP_LONGINT_FORMAT
               " is negative or larger than blocksize %d!\n",
                       col, *linecount, data -> sdpblocksizes[b]);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( SCIPisZero(scip, val))
            {
               ++ nzerocoef;
            }
            else
            {
               /* make sure matrix is in lower triangular form */

               if( col > row )
               {
                  data -> sdpconstrow[b][data -> sdpconstnblocknonz[b]] = col;
                  data -> sdpconstcol[b][data -> sdpconstnblocknonz[b]] = row;
               }
               else
               {
                  data -> sdpconstrow[b][data -> sdpconstnblocknonz[b]] = row;
                  data -> sdpconstcol[b][data -> sdpconstnblocknonz[b]] = col;
               }
               data -> sdpconstval[b][data -> sdpconstnblocknonz[b]] = val;
               data -> sdpconstnblocknonz[b] ++;
            }


         }
      }
      else
      {
         if( v >= 0 )
         { //TIM: ob ACoord oder BCoord
            c = row; //TIM: macht einfach row zur constraint
            if( c < 0 || c >= data -> nconss )
            {
               SCIPerrorMessage("Given linear coefficient in line %"
               SCIP_LONGINT_FORMAT
               " for constraint %d which does not exist!\n", *linecount, c);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( v < 0 || v >= data -> nvars )
            {
               SCIPerrorMessage("Given linear coefficient in line %"
               SCIP_LONGINT_FORMAT
               " for variable %d which does not exist!\n", *linecount, v);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( SCIPisZero(scip, val))
            {
               ++ nzerocoef;
            }
            else
            {
               SCIP_CALL(SCIPaddCoefLinear(scip, data -> createdconss[c], data -> createdvars[v],
                                           val));/*lint !e732*//*lint !e747*/
            }
         }
         else
         {
            c = row;

            if( c < 0 || c >= data -> nconss )
            {
               SCIPerrorMessage("Given constant part in line %"
               SCIP_LONGINT_FORMAT
               " for scalar constraint %d which does not exist!\n", *linecount, c);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }

            if( SCIPisZero(scip, val))
            {
               ++ nzerocoef;
            }
            else
            {
               /* check type */
               if( ! SCIPisInfinity(scip, - SCIPgetLhsLinear(scip, data -> createdconss[c])))
               {
                  /* greater or equal constraint -> left-hand side*/  
                  SCIP_CALL(SCIPchgLhsLinear(scip, data -> createdconss[c], val));
               }

               if( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, data -> createdconss[c])))
               {
                  /* less or equal constraint -> right-hand side*/
                  SCIP_CALL(SCIPchgRhsLinear(scip, data -> createdconss[c], val));
               }
            }
         }
      }
   }

   /* construct pointer arrays */
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpnblockvars), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpblockvars), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> nvarnonz), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> rowpointer), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> colpointer), data -> nsdpblocks));
   SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> valpointer), data -> nsdpblocks));

   /* sdp blocks as specified in cbf file in HCOORD */
   for( b = 0; b < ncbfsdpblocks; b ++ ){
      /* sort the nonzeroes by non-decreasing variable indices */
      SCIPsortIntIntIntReal(sdpvar[b], data -> sdprow[b], data -> sdpcol[b], data -> sdpval[b],
                            data -> sdpnblocknonz[b]);

      /* create the pointer arrays and insert used variables into vars-array */
      nextindaftervar = 0;
      data -> sdpnblockvars[b] = 0;
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> sdpblockvars[b]), data -> nvars));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> nvarnonz[b]), data -> nvars));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> rowpointer[b]), data -> nvars));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> colpointer[b]), data -> nvars));
      SCIP_CALL(SCIPallocBlockMemoryArray(scip, &(data -> valpointer[b]), data -> nvars));

      for( v = 0; v < data -> nvars; v ++ )
      {
         SCIP_Bool varused = FALSE;

         firstindforvar = nextindaftervar; /* this variable starts where the last one ended */
         data -> nvarnonz[b][data -> sdpnblockvars[b]] = 0;

         while( nextindaftervar < data -> sdpnblocknonz[b] &&
                sdpvar[b][nextindaftervar] == v ) /* get the first index that doesn't belong to this variable */
         {
            nextindaftervar ++;
            varused = TRUE;
            data -> nvarnonz[b][data -> sdpnblockvars[b]] ++;
         }

         if( varused )
         {
            data -> sdpblockvars[b][data -> sdpnblockvars[b]] = data -> createdvars[v];/*lint !e732*//*lint !e747*/ /* if the variable is used, add it to the vars array */
            data -> rowpointer[b][data -> sdpnblockvars[b]] = &(data -> sdprow[b][firstindforvar]); /* save a pointer to the first nonzero belonging to this variable */
            data -> colpointer[b][data -> sdpnblockvars[b]] = &(data -> sdpcol[b][firstindforvar]);
            data -> valpointer[b][data -> sdpnblockvars[b]] = &(data -> sdpval[b][firstindforvar]);
            data -> sdpnblockvars[b] ++;
         }
      }

      assert(nextindaftervar == data -> sdpnblocknonz[b]);
   }


   /* free SDP-var array which is no longer needed */
   for( b = 0; b < data -> nsdpblocks; b ++ )
      SCIPfreeBlockMemoryArray(scip, &(sdpvar[b]), data -> memorysizessdp[b]);

   SCIPfreeBlockMemoryArray(scip, &sdpvar, data -> nsdpblocks);

   if( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                      "HCOORD: Found %d coefficients with absolute value less than epsilon = %g.\n", nzerocoef,
                      SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}


/** reads integrality conditions from given SDPA-file */
static
SCIP_RETCODE SDPAreadInt(
        SCIP *scip,               /**< SCIP data structure */
        SCIP_FILE *scipfile,              /**< file to read from */
        SCIP_Longint *linecount,          /**< current linecount */
        SDPA_DATA *data                /**< data pointer to save the results in */
){  /*lint --e{818}*/

   int v;
   SCIP_Bool infeasible;

   assert(scip != NULL);
   assert(scipfile != NULL);
   assert(linecount != NULL);
   assert(data != NULL);

   if( data -> createdvars == NULL)
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'INT' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert(data -> nvars >= 0);


   while( fIntgets(scipfile, linecount) == SCIP_OKAY )
   {
      char *ps = SDPA_LINE_BUFFER + 1;
      if( sscanf(ps, "%i", &v) != 1 )
      {
         SCIPerrorMessage("Could not read variable index in line %"
         SCIP_LONGINT_FORMAT
         ".\n", *linecount);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

	if( v < 1 || v > data -> nvars )
            {
               SCIPerrorMessage("Given integrality in line %"
               SCIP_LONGINT_FORMAT
               " for variable %d which does not exist!\n", *linecount, v);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/}	



      v = v - 1;

      SCIP_CALL(SCIPchgVarType(scip, data -> createdvars[v], SCIP_VARTYPE_INTEGER, &infeasible));


      if( infeasible )
      {
         SCIPerrorMessage("Infeasibility detected because of integrality of variable %s!\n",
                          SCIPvarGetName(data -> createdvars[v]));
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }

   return SCIP_OKAY;
}


/** frees all data allocated for the SDPA-data-struct */
static
SCIP_RETCODE SDPAfreeData(
        SCIP *scip,               /**< SCIP data structure */
        SDPA_DATA *data                /**< data pointer to save the results in */
){
   SCIP_Bool allocated = FALSE;
   int b = 0;


   int ncbfsdpblocks;

   assert(scip != NULL);
   assert(data != NULL);

   /* we only allocated memory for the const blocks if there were any nonzeros */
   /* TODO: could also think about saving this in the struct instead, which would cost one bool and save some (unimportant) time here */
   while( data -> sdpconstnblocknonz != NULL && allocated == FALSE && b < data -> nsdpblocks )
   {
      if( data -> sdpconstnblocknonz[b] > 0 )
         allocated = TRUE;
      b ++;
   }

   if( allocated )
   {
      for( b = 0; b < data -> nsdpblocks; b ++ )
      {
         SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpconstval[b]), data -> memorysizescon[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpconstcol[b]), data -> memorysizescon[b]);
         SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpconstrow[b]), data -> memorysizescon[b]);
      }

      SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpconstval, data -> nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpconstcol, data -> nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpconstrow, data -> nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpconstnblocknonz, data -> nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data -> memorysizescon, data -> nsdpblocks);

   }

   /* we only allocated memory for the sdpblocks if there were any nonzeros */
   b = 0;
   while( allocated == FALSE && b < data -> nsdpblocks )
   {
      if( data -> sdpnblocknonz[b] > 0 )
         allocated = TRUE;
      b ++;
   }

   if( allocated )
   {
      /* get number of sdp blocks specified by PSDCON (without auxiliary sdp blocks for reformulating matrix variables
       * using scalar variables), save number of nonzeros needed for the auxiliary sdp blocks in nauxnonz */

      ncbfsdpblocks = data -> nsdpblocks;

     
         /* some SDP constraints specified in the SDPA file! */
         assert(ncbfsdpblocks > 0);

         for( b = 0; b < ncbfsdpblocks; b ++ )
         {
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> valpointer[b]), data -> nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> colpointer[b]), data -> nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> rowpointer[b]), data -> nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpval[b]), data -> memorysizessdp[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpcol[b]), data -> memorysizessdp[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdprow[b]), data -> memorysizessdp[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> sdpblockvars[b]), data -> nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data -> nvarnonz[b]), data -> nvars);
         }

         SCIPfreeBlockMemoryArrayNull(scip, &data -> memorysizessdp, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> valpointer, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> colpointer, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> rowpointer, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpval, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpcol, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdprow, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpblockvars, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> nvarnonz, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpnblockvars, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpnblocknonz, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpblockrank1, data -> nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data -> sdpblocksizes, data -> nsdpblocks);
      
   }

   if( data -> nconss > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &data -> createdconss, data -> nconss);
   }

   if( data -> nvars > 0 )
      SCIPfreeBlockMemoryArrayNull(scip, &data -> createdvars, data -> nvars);

   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCbf)
        {  /*lint --e{715,818}*/
                assert(scip != NULL);

        SCIP_CALL( SCIPincludeReaderCbf(scip));

        return SCIP_OKAY;
        }


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCbf)
        {  /*lint --e{715,818}*/
                SCIP_FILE * scipfile;
        SCIP_Longint linecount = 0;
        SCIP_Bool versionread = FALSE;
        SCIP_Bool objread = FALSE;
        SDPA_DATA* data;
        int b;

        assert(result != NULL);

        *result = SCIP_DIDNOTRUN;

        SCIPdebugMsg(scip, "Reading file %s ...\n", filename);

        scipfile = SCIPfopen(filename, "r");

        if( ! scipfile )
        return SCIP_READERROR;

        SCIP_CALL( SCIPallocBuffer(scip, &data));
        data->nsdpblocks = -1;
        data->nsdpblocksrank1 = 0;
        data->nconss = 0;
        data->nvars = -1;
        data->constnnonz = 0;
        data->nnonz = 0;

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
        data->locationConBlock=-1;

        /* create empty problem */
        SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL));

        SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));

        SCIPdebugMsg(scip, "Reading VAR\n");             
        SCIP_CALL( SDPAreadVar(scip, scipfile, &linecount, data));

        SCIPdebugMsg(scip, "Reading CON\n");               
        SCIP_CALL( SDPAreadCon(scip, scipfile, &linecount, data));

        SCIPdebugMsg(scip, "Reading PSDCON\n");        
        SCIP_CALL( SDPAreadPsdCon(scip, scipfile, &linecount, data));

        SCIPdebugMsg(scip, "Reading OBJACOORD\n");           
        SCIP_CALL( SDPAreadObjAcoord(scip, scipfile, &linecount, data));

        SCIPdebugMsg(scip, "Reading HCOORD\n");          
        SCIP_CALL( SDPAreadHcoord(scip, scipfile, &linecount, data, filename));

        SCIPdebugMsg(scip, "Reading INT\n");                  
        SCIP_CALL( SDPAreadInt(scip, scipfile, &linecount, data));

/* close the file (and make sure SCIPfclose returns 0) */
        if( SCIPfclose(scipfile))
        return SCIP_READERROR;

#ifdef SCIP_MORE_DEBUG
   for (b = 0; b < SCIPgetNConss(scip); b++)
   {
      SCIP_CALL( SCIPprintCons(scip, SCIPgetConss(scip)[b], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
   }
#endif

/* create SDP-constraints */
        for(b = 0; b < data->nsdpblocks; b++)
        {
           SCIP_CONS *sdpcons;
           char sdpconname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
           int snprintfreturn;
#endif
           assert(data -> sdpblocksizes[b] > 0);
           assert((data -> sdpnblockvars[b] > 0 && data -> sdpnblocknonz[b] > 0) ||
                  (data -> sdpconstnblocknonz[b] > 0));
#ifndef NDEBUG
           snprintfreturn = SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
           assert(snprintfreturn < SCIP_MAXSTRLEN);
#else
           (void) SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
#endif

/* special treatment of case without constant PSD blocks */
           if( data -> sdpconstnblocknonz == NULL)
           {
              if( ! data -> sdpblockrank1[b] )
              {
                 SCIP_CALL(SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data -> sdpnblockvars[b],
                                             data -> sdpnblocknonz[b],
                                             data -> sdpblocksizes[b], data -> nvarnonz[b], data -> colpointer[b],
                                             data -> rowpointer[b], data -> valpointer[b],
                                             data -> sdpblockvars[b], 0, NULL, NULL, NULL));
              }
              else
              {
                 SCIP_CALL(SCIPcreateConsSdpRank1(scip, &sdpcons, sdpconname, data -> sdpnblockvars[b],
                                                  data -> sdpnblocknonz[b],
                                                  data -> sdpblocksizes[b], data -> nvarnonz[b], data -> colpointer[b],
                                                  data -> rowpointer[b], data -> valpointer[b],
                                                  data -> sdpblockvars[b], 0, NULL, NULL, NULL));
              }
           }
           else
           {
              if( ! data -> sdpblockrank1[b] )
              {
                  SCIP_CALL(SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data -> sdpnblockvars[b],
                                             data -> sdpnblocknonz[b],
                                             data -> sdpblocksizes[b], data -> nvarnonz[b], data -> colpointer[b],
                                             data -> rowpointer[b], data -> valpointer[b],
                                             data -> sdpblockvars[b], data -> sdpconstnblocknonz[b],
                                             data -> sdpconstcol[b], data -> sdpconstrow[b],
                                             data -> sdpconstval[b]));
              }
              else
              {
                 SCIP_CALL(SCIPcreateConsSdpRank1(scip, &sdpcons, sdpconname, data -> sdpnblockvars[b],
                                                  data -> sdpnblocknonz[b],
                                                  data -> sdpblocksizes[b], data -> nvarnonz[b], data -> colpointer[b],
                                                  data -> rowpointer[b], data -> valpointer[b],
                                                  data -> sdpblockvars[b], data -> sdpconstnblocknonz[b],
                                                  data -> sdpconstcol[b], data -> sdpconstrow[b],
                                                  data -> sdpconstval[b]));
              }
           }

#ifdef SCIP_MORE_DEBUG
           SCIP_CALL( SCIPprintCons(scip, sdpcons, NULL) );
#endif

           SCIP_CALL(SCIPaddCons(scip, sdpcons));

           SCIP_CALL(SCIPreleaseCons(scip, &sdpcons));
        }

        SCIP_CALL( SDPAfreeData(scip, data));
        SCIPfreeBufferNull(scip, &data);

        *result = SCIP_SUCCESS;

        return SCIP_OKAY;
        }


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCbf)
        {  /*lint --e{715,818}*/
        SCIP_VAR * *linvars;
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
        int nsdpconss;
        int sdpnvars;
        int sdpnnonz;
        int totalsdpnnonz;
        int sdpblocksize;
        int sdparraylength;
        int totalsdpconstnnonz;
        int sdpconstnnonz;
        int nobjnonz;
        int nnonz;
        int nbnonz;
        int nsenses;
        int nconsssenses;
        int lastsense;
        int consind;
        int c;
        int i;
        int v;
        int nrank1sdpblocks;

        assert(scip != NULL);
        assert(result != NULL);

        SCIPdebugMsg(scip, "Writing problem in SDPA format to file.\n");
        *result = SCIP_DIDNOTRUN;

        if( transformed )
        {
           SCIPerrorMessage("SDPA reader currently only supports writing original problems!\n");
           SCIPABORT();
           return SCIP_READERROR; /*lint !e527*/
        }

        for(c = 0; c < nconss; c++)
        {
           if((strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0)
              && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0)
              && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0)){
              SCIPerrorMessage("SDPA reader currently only supports linear and SDP constraints!\n");
              SCIPABORT();
              return SCIP_READERROR; /*lint !e527*/
           }
        }

#ifndef NDEBUG
        for(v = 0; v < nvars; v++)
        {
           assert(SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL);
        }
#endif

        /* write version number */
        SCIPinfoMessage(scip, file, "VER\n%d\n\n", SDPA_VERSION_NR);

        /* write objective sense */
        SCIPinfoMessage(scip, file, "OBJSENSE\n%s\n\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "MIN" : "MAX");

        /* collect different variable senses */
        SCIP_CALL( SCIPallocBufferArray(scip, &varsenses, nvars));
        for(v = 0; v < nvars; v++)
        {
           SCIP_Real lb;
           SCIP_Real ub;

           lb = SCIPvarGetLbOriginal(vars[v]);
           ub = SCIPvarGetUbOriginal(vars[v]);

           varsenses[v] = 0;
           if( SCIPisZero(scip, lb))
              varsenses[v] = 1;
           else{
              if( ! SCIPisInfinity(scip, - lb)){
                 SCIPerrorMessage("Can only handle variables with lower bound 0 or minus infinity.\n");
                 SCIPABORT();
                 return SCIP_READERROR; /*lint !e527*/
              }
           }

           if( SCIPisZero(scip, ub))
              varsenses[v] = - 1;
           else{
              if( ! SCIPisInfinity(scip, ub)){
                 SCIPerrorMessage("Can only handle variables with upper bound 0 or infinity.\n");
                 SCIPABORT();
                 return SCIP_READERROR; /*lint !e527*/
              }
           }
        }

        /* now determine senses of constraints - possibly disable constraints */
        SCIP_CALL( SCIPallocBufferArray(scip, &consssenses, nconss));
        for(c = 0; c < nconss; c++)
        {
           SCIP_Real lhs;
           SCIP_Real rhs;

           consssenses[c] = 2;

           /* only count linear constraints */
           if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
              continue;

           lhs = SCIPgetLhsLinear(scip, conss[c]);
           rhs = SCIPgetRhsLinear(scip, conss[c]);

#ifdef SDPA_CHECK_NONNEG
           /* check if there are constraints that would determine senses of variables */
           if( SCIPgetNVarsLinear(scip, conss[c]) == 1 &&
               ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))){
              /* the nonzero should be a true nonzero */
              assert(! SCIPisZero(scip, SCIPgetValsLinear(scip, conss[c])[0]));
              assert(SCIPgetVarsLinear(scip, conss[c]) != NULL);

              v = SCIPvarGetProbindex(SCIPgetVarsLinear(scip, conss[c])[0]);
              assert(0 <= v && v < nvars);

              if( SCIPgetValsLinear(scip, conss[c])[0] > 0.0 ){
                 if( SCIPisZero(scip, lhs)){
                    varsenses[v] = 1;
                    continue;
                 }else if( SCIPisZero(scip, rhs)){
                    varsenses[v] = - 1;
                    continue;
                 }
              }else{
                 if( SCIPisZero(scip, lhs)){
                    varsenses[v] = - 1;
                    continue;
                 }else if( SCIPisZero(scip, rhs)){
                    varsenses[v] = 1;
                    continue;
                 }
              }
           }
#endif

           if( SCIPisEQ(scip, lhs, rhs)){
              assert(! SCIPisInfinity(scip, - lhs));
              assert(! SCIPisInfinity(scip, rhs));
              consssenses[c] = 0;
           }else{
              if( ! SCIPisInfinity(scip, - lhs) && ! SCIPisInfinity(scip, rhs)){
                 SCIPerrorMessage("Cannot handle ranged rows.\n");
                 SCIPABORT();
                 return SCIP_READERROR; /*lint !e527*/
              }

              if( ! SCIPisInfinity(scip, - lhs))
                 consssenses[c] = 1;
              else if( ! SCIPisInfinity(scip, rhs))
                 consssenses[c] = - 1;
           }
        }

        /* compute different varsenses */
        nsenses = 0;
        lastsense = 2;
        for(v = 0; v < nvars; v++)
        {
           if( varsenses[v] != lastsense ){
              ++ nsenses;
              lastsense = varsenses[v];
           }
        }

        /* write variable senses */
        SCIPinfoMessage(scip, file, "VAR\n%d %d\n", nvars, nsenses);
        lastsense = varsenses[0];
        nsenses = 1;
        for(v = 1; v < nvars; v++)
        {
           if( varsenses[v] != lastsense ){
              if( lastsense == 0 )
                 SCIPinfoMessage(scip, file, "F %d\n", nsenses);
              else if( lastsense == - 1 )
                 SCIPinfoMessage(scip, file, "L- %d\n", nsenses);
              else{
                 assert(lastsense == 1);
                 SCIPinfoMessage(scip, file, "L+ %d\n", nsenses);
              }
              nsenses = 0;
              lastsense = varsenses[v];
           }
           ++ nsenses;
        }
        if( lastsense == 0 )
        SCIPinfoMessage(scip, file, "F %d\n\n", nsenses);
        else if( lastsense == -1 )
        SCIPinfoMessage(scip, file, "L- %d\n\n", nsenses);
        else
        {
           assert(lastsense == 1);
           SCIPinfoMessage(scip, file, "L+ %d\n\n", nsenses);
        }

        /* write integrality constraints */
        if( nbinvars + nintvars > 0 )
        {
           SCIPinfoMessage(scip, file, "INT\n%d\n", nbinvars + nintvars);

           for( v = 0; v < nbinvars + nintvars; v ++ ){
              assert(SCIPvarIsIntegral(vars[v]));
              SCIPinfoMessage(scip, file, "%d\n", v);
           }
           SCIPinfoMessage(scip, file, "\n");
        }

        /* compute different consssenses */
        nsenses = 0;
        lastsense = 3;
        i = 0;
        for(c = 0; c < nconss; c++)
        {
           if( consssenses[c] == 2 )
              continue;

           ++ i;
           if( consssenses[c] != lastsense ){
              ++ nsenses;
              lastsense = consssenses[c];
           }
        }

        /* write constraint senses */
        SCIPinfoMessage(scip, file, "CON\n%d %d\n", i, nsenses);
        nconsssenses = nsenses;
        c = 0;
        while(c < nconss && consssenses[c] == 2)
        ++c;
        if( c < nconss )
        {
           lastsense = consssenses[c];
           nsenses = 1;
           ++ c;
           for( ; c < nconss; ++ c ){
              if( consssenses[c] == 2 )
                 continue;

              if( consssenses[c] != lastsense ){
                 if( lastsense == 0 )
                    SCIPinfoMessage(scip, file, "L= %d\n", nsenses);
                 else if( lastsense == - 1 )
                    SCIPinfoMessage(scip, file, "L- %d\n", nsenses);
                 else{
                    assert(lastsense == 1);
                    SCIPinfoMessage(scip, file, "L+ %d\n", nsenses);
                 }
                 nsenses = 0;
                 lastsense = consssenses[c];
              }
              ++ nsenses;
           }
        }
        if( lastsense == 0 )
        SCIPinfoMessage(scip, file, "L= %d\n\n", nsenses);
        else if( lastsense == -1 )
        SCIPinfoMessage(scip, file, "L- %d\n\n", nsenses);
        else
        {
           assert(lastsense == 1);
           SCIPinfoMessage(scip, file, "L+ %d\n\n", nsenses);
        }

        /* count number of SDP constraints (conshdlrGetNConss doesn't seem to work before transformation) */
        nsdpconss = 0;
        for(c = 0; c < nconss; c++)
        {
           if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 ||
               strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0 )
              ++ nsdpconss;
        }

        /* write SDP constraints */
        SCIPinfoMessage(scip, file, "PSDCON\n%d\n", nsdpconss);

        for(c = 0; c < nconss; c++)
        {
           if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 &&
               strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0 )
              continue;

           SCIPinfoMessage(scip, file, "%d\n", SCIPconsSdpGetBlocksize(scip, conss[c]));
        }

        SCIPinfoMessage(scip, file, "\n");

        /* count number of rank-1 SDP constraints */
        nrank1sdpblocks = 0;
        for(c = 0; c < nconss; c++)
        {
           if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") == 0 )
              ++ nrank1sdpblocks;
        }

        /* write rank-1 SDP constraints (if existing) */
        if( nrank1sdpblocks > 0 )
        {
           SCIPinfoMessage(scip, file, "PSDCONRANK1\n%d\n", nrank1sdpblocks);
           consind = 0;

           for( c = 0; c < nconss; c ++ ){
              if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0 )
                 continue;

              assert(SCIPconsSdpShouldBeRankOne(conss[c]));
              SCIPinfoMessage(scip, file, "%d\n", consind);
              consind ++;
           }

           SCIPinfoMessage(scip, file, "\n");
        }

        /* count number of nonzero objective coefficients */
        nobjnonz = 0;
        for(v = 0; v < nvars; v++)
        {
           if( ! SCIPisZero(scip, SCIPvarGetObj(vars[v])))
              ++ nobjnonz;
        }

        /* write objective */
        SCIPinfoMessage(scip, file, "OBJACOORD\n%d\n", nobjnonz);

        for(v = 0; v < nvars; v++)
        {
           SCIP_Real obj;

           obj = SCIPvarGetObj(vars[v]);
           if( ! SCIPisZero(scip, obj)){
              SCIPinfoMessage(scip, file, "%d %.15g\n", v, obj);
           }
        }
        SCIPinfoMessage(scip, file, "\n");

        /* write coefficients of linear constraints */
        if( nconsssenses > 0 )
        {
           /* count number of nonzero coefficients in linear constraints */
           nnonz = 0;
           nbnonz = 0;
           for( c = 0; c < nconss; c ++ ){
              if( consssenses[c] == - 1 ){
                 assert(! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])));
                 nnonz += SCIPgetNVarsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])))
                    ++ nbnonz;
              }else if( consssenses[c] == 1 ){
                 assert(! SCIPisInfinity(scip, - SCIPgetLhsLinear(scip, conss[c])));
                 nnonz += SCIPgetNVarsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])))
                    ++ nbnonz;
              }else if( consssenses[c] == 0 ){
                 assert(SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])));
                 nnonz += SCIPgetNVarsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])))
                    ++ nbnonz;
              }
           }

           /* write linear nonzero coefficients */
           SCIPinfoMessage(scip, file, "ACOORD\n%d\n", nnonz);
           consind = 0;
           for( c = 0; c < nconss; c ++ ){
              if( consssenses[c] == - 1 || consssenses[c] == 0 || consssenses[c] == 1 ){
                 assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") == 0);

                 linvars = SCIPgetVarsLinear(scip, conss[c]);
                 linvals = SCIPgetValsLinear(scip, conss[c]);

                 for( v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v ++ ){
                    i = SCIPvarGetProbindex(linvars[v]);
                    assert(0 <= i && i < nvars);
                    SCIPinfoMessage(scip, file, "%d %d %.15g\n", consind, i, linvals[v]);
                 }
                 ++ consind;
              }
           }
           SCIPinfoMessage(scip, file, "\n");

           /* write constant part of linear constraints */
           SCIPinfoMessage(scip, file, "BCOORD\n%d\n", nbnonz);
           consind = 0;
           for( c = 0; c < nconss; c ++ ){
              SCIP_Real val;
              if( consssenses[c] == - 1 ){
                 val = SCIPgetRhsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, val))
                    SCIPinfoMessage(scip, file, "%d %.15g\n", consind, - val);
                 consind ++;
              }else if( consssenses[c] == 1 ){
                 val = SCIPgetLhsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, val))
                    SCIPinfoMessage(scip, file, "%d %.15g\n", consind, - val);
                 consind ++;
              }else if( consssenses[c] == 0 ){
                 val = SCIPgetLhsLinear(scip, conss[c]);
                 if( ! SCIPisZero(scip, val))
                    SCIPinfoMessage(scip, file, "%d %.15g\n", consind, - val);
                 consind ++;
              }
           }
           SCIPinfoMessage(scip, file, "\n");
        }

        /* count SDP nonzeros */
        totalsdpnnonz = 0;
        totalsdpconstnnonz = 0;
        for(c = 0; c < nconss; c++)
        {
           if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 &&
               strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0 )
              continue;

           SCIP_CALL(SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz));
           totalsdpnnonz += sdpnnonz;
           totalsdpconstnnonz += sdpconstnnonz;
        }

        /* allocate memory for SDPdata */
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpnvarnonz, nvars));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpcol, totalsdpnnonz));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdprow, totalsdpnnonz));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpval, totalsdpnnonz));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpvars, nvars));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstcol, totalsdpconstnnonz));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstrow, totalsdpconstnnonz));
        SCIP_CALL( SCIPallocBufferArray(scip, &sdpconstval, totalsdpconstnnonz));

        sdparraylength = totalsdpnnonz;
        sdpconstnnonz = totalsdpconstnnonz;

        /* write SDP nonzeros */
        if( totalsdpnnonz > 0 )
        {
           SCIPinfoMessage(scip, file, "HCOORD\n%d\n", totalsdpnnonz);
           consind = 0;
           for( c = 0; c < nconss; c ++ ){
              if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 &&
                  strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0 )
                 continue;

              /* initialization for SDPconsSDPGetData-call */
              sdparraylength = totalsdpnnonz;
              sdpconstnnonz = totalsdpconstnnonz;

              SCIP_CALL(SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength,
                                           sdpnvarnonz,
                                           sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow,
                                           sdpconstval, NULL, NULL, NULL));

              assert(sdpconstnnonz <= totalsdpconstnnonz);
              assert(sdparraylength <= totalsdpnnonz);

              for( v = 0; v < sdpnvars; v ++ ){
                 for( i = 0; i < sdpnvarnonz[v]; i ++ ){
                    int ind;
                    ind = SCIPvarGetProbindex(sdpvars[v]);
                    assert(0 <= ind && ind < nvars);
                    SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", consind, ind, sdprow[v][i], sdpcol[v][i],
                                    sdpval[v][i]);
                 }
              }
              consind ++;
           }
           SCIPinfoMessage(scip, file, "\n");
        }

        /* write nonzeros of constant part of SDP constraint */
        if( totalsdpnnonz > 0 )
        {
           SCIPinfoMessage(scip, file, "DCOORD\n%d\n", totalsdpconstnnonz);
           consind = 0;
           for( c = 0; c < nconss; c ++ ){
              if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 &&
                  strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDPrank1") != 0 )
                 continue;

              /* initialization for SDPconsSDPGetData-call */
              sdparraylength = totalsdpnnonz;
              sdpconstnnonz = totalsdpconstnnonz;

              SCIP_CALL(SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength,
                                           sdpnvarnonz,
                                           sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow,
                                           sdpconstval, NULL, NULL, NULL));

              assert(sdpconstnnonz <= totalsdpconstnnonz);
              assert(sdparraylength <= totalsdpnnonz);

              for( i = 0; i < sdpconstnnonz; i ++ ){
                 SCIPinfoMessage(scip, file, "%d %d %d %.15g\n", consind, sdpconstrow[i], sdpconstcol[i],
                                 - sdpconstval[i]);
              }
              consind ++;
           }
        }

        SCIPfreeBufferArray(scip, &sdpconstval);
        SCIPfreeBufferArray(scip, &sdpconstrow);
        SCIPfreeBufferArray(scip, &sdpconstcol);
        SCIPfreeBufferArray(scip, &sdpvars);
        SCIPfreeBufferArray(scip, &sdpval);
        SCIPfreeBufferArray(scip, &sdprow);
        SCIPfreeBufferArray(scip, &sdpcol);
        SCIPfreeBufferArray(scip, &sdpnvarnonz);
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
        SCIP *scip                /**< SCIP data structure */
){
   SCIP_READERDATA *readerdata = NULL;
   SCIP_READER *reader;

   /* include reader */
   SCIP_CALL(SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata));

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL(SCIPsetReaderCopy(scip, reader, readerCopyCbf));
   SCIP_CALL(SCIPsetReaderRead(scip, reader, readerReadCbf));
   SCIP_CALL(SCIPsetReaderWrite(scip, reader, readerWriteCbf));

   return SCIP_OKAY;
}
