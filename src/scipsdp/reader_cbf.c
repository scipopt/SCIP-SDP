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

/* #define SCIP_MORE_DEBUG */

/**@file   reader_cbf.c
 * @brief  file reader for mixed-integer semidefinite programs in CBF format
 * @author Tristan Gally
 * @author Henrik A. Friberg
 * @author Marc Pfetsch
 * @author Frederic Matter
 *
 * @todo Allow to read SOC constraints in CBF format.
 * @todo Allow to write varbounds as linear constraints.
 * @todo Allow to write a transformed problem.
 * @todo Allow to write a problem in primal form.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>                      /* for strcmp */

#include "scipsdp/reader_cbf.h"
#include "scipsdp/cons_sdp.h"
#include "scip/cons_linear.h"


#define READER_NAME             "cbfreader"
#define READER_DESC             "file reader and writer for MISDPs in cbf format"
#define READER_EXTENSION        "cbf"

#define CBF_VERSION_NR         2         /**< version number for CBF format */
#define CBF_CHECK_NONNEG       TRUE      /**< when writing: check linear constraints and move nonnegativity(-positivity)
                                           *  constraints to definition of variables (which are now defined in non-negative
                                           *  orthant) */
                                          /*  TODO: currently doesn't work for ranged rows (which are not created by sdpa
                                           *  reader) */

/* Use CBF_NAME_FORMAT instead of %s when parsing lines, to avoid buffer overflow. */
#define MACRO_STR_EXPAND(tok) #tok
#define MACRO_STR(tok) MACRO_STR_EXPAND(tok)
#define CBF_NAME_FORMAT "%" MACRO_STR(CBF_MAX_NAME) "s"
#define CBF_MAX_LINE  512       /* Last 3 chars reserved for '\r\n\0' */
#define CBF_MAX_NAME  512

char CBF_LINE_BUFFER[CBF_MAX_LINE];
char CBF_NAME_BUFFER[CBF_MAX_NAME];

struct CBF_Data
{
   int                   npsdvars;           /**< number of psd variables and length of createdpsdvars-array  */
   int*                  psdvarsizes;        /**< sizes of the psd variables */
   SCIP_VAR****          createdpsdvars;     /**< array of psd variables created by the CBF reader */
   SCIP_Bool             noorigsdpcons;      /**< Are there SDP constraints specified in the CBF file?  */

   int                   nvars;              /**< number of variables and length of createdvars-array */
   SCIP_VAR**            createdvars;        /**< array of variables created by the CBF reader */
   int                   nconss;             /**< number of constraints and length of createdconss-array */
   SCIP_CONS**           createdconss;       /**< array of constraints created by the CBF reader */

   int                   nsdpblocks;         /**< number of SDP constraints/blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP blocks */
   int*                  sdpnblocknonz;      /**< number of nonzeros for each SDP block */
   int*                  sdpnblockvars;      /**< number of variables for each SDP block */
   int**                 nvarnonz;           /**< number of nonzeros for each block and each variable */
   SCIP_VAR***           sdpblockvars;       /**< SCIP variables appearing in each block */
   int**                 sdprow;             /**< array of all row indices for each SDP block */
   int**                 sdpcol;             /**< array of all column indices for each SDP block */
   SCIP_Real**           sdpval;             /**< array of all values of SDP nonzeros for each SDP block */
   int                   nnonz;              /**< number of nonzeros in blocks */
   int***                rowpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int***                colpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   SCIP_Real***          valpointer;         /**< array of pointers to first entries in row-array for each block and variable */
   int*                  sdpconstnblocknonz; /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                              *   number of entries of sdpconst row/col/val [i] */
   int**                 sdpconstrow;        /**< pointers to row-indices for each block */
   int**                 sdpconstcol;        /**< pointers to column-indices for each block */
   SCIP_Real**           sdpconstval;        /**< pointers to the values of the nonzeros for each block */
   int                   constnnonz;         /**< number of nonzeros in const block */
};

typedef struct CBF_Data CBF_DATA;


/*
 * Local methods
 */

/** finds first non-commentary line in given file */
static
SCIP_RETCODE CBFfgets(
   SCIP_FILE*            pFile,              /**< file to read from */
   SCIP_Longint*         linecount           /**< current linecount */
   )
{
   assert( pFile != NULL );
   assert( linecount != NULL );

   /* Find first non-commentary line */
   while ( SCIPfgets(CBF_LINE_BUFFER, (int) sizeof(CBF_LINE_BUFFER), pFile) != NULL )
   {
      ++(*linecount);

      if ( CBF_LINE_BUFFER[0] != '#' )
         return SCIP_OKAY;
   }

   return SCIP_READERROR;
}

/** reads objective sense from given CBF-file */
static
SCIP_RETCODE CBFreadObjsense(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount           /**< current linecount */
   )
{
   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER) == 1 )
   {
      if ( strcmp(CBF_NAME_BUFFER, "MIN") == 0 )
      {
         SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
      }
      else if ( strcmp(CBF_NAME_BUFFER, "MAX") == 0 )
      {
         SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
      }
      else
      {
         SCIPerrorMessage("OBJSENSE should be either MIN or MAX.\n");
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
   {
      return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** reads the number and type of scalar variables from given CBF-file */
static
SCIP_RETCODE CBFreadVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   int nvartypevars;
   int nvartypes;
   int cnt = 0;
   int t;
   int v;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i %i", &(data->nvars), &nvartypes) != 2 )
      return SCIP_READERROR;

   if ( data->nvars < 0 )
   {
      SCIPerrorMessage("Number of scalar variables %d should be non-negative!\n", data->nvars);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   if ( nvartypes < 0 )
   {
      SCIPerrorMessage("Number of conic variable domains %d should be non-negative!\n", nvartypes);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 && nvartypes >= 0 );

   /* loop through different variable types */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdvars), data->nvars) );
   for (t = 0; t < nvartypes; t++)
   {
      SCIP_CALL( CBFfgets(pfile, linecount) );

      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %i", CBF_NAME_BUFFER, &nvartypevars) == 2 )
      {
         SCIP_Real lb;
         SCIP_Real ub;

         if ( nvartypevars <= 0 )
         {
            SCIPerrorMessage("Number of scalar variables %d should be positive!\n", nvartypevars);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         lb = -SCIPinfinity(scip);
         ub = SCIPinfinity(scip);

         if ( strcmp(CBF_NAME_BUFFER, "L+") == 0 )
         {
            lb = 0.0;
         }
         else if ( strcmp(CBF_NAME_BUFFER, "L-") == 0 )
         {
            ub = 0.0;
         }
         else if ( strcmp(CBF_NAME_BUFFER, "F") != 0 )
         {
            SCIPerrorMessage("CBF-Reader of SCIP-SDP currently only supports non-negative, non-positive and free variables!\n");
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         /* create corresponding variables */
         for (v = 0; v < nvartypevars; v++)
         {
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
            assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
            (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "x_%d", cnt);
#endif

            SCIP_CALL( SCIPcreateVar(scip, &var, varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS,
                  TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));/*lint !e732*//*lint !e747*/

            SCIP_CALL( SCIPaddVar(scip, var) );
            data->createdvars[cnt++] = var;/*lint !e732*//*lint !e747*/

            /* release variable for the reader */
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
      }
      else
         return SCIP_READERROR;
   }

   if ( cnt != data->nvars )
   {
      SCIPerrorMessage("Total number of scalar variables for different cone types not equal to total number of scalar variables!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** reads the number and type of psd variables from given CBF-file */
static
SCIP_RETCODE CBFreadPsdvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   char varname[SCIP_MAXSTRLEN];
   SCIP_VAR* var;
   int i;
   int j;
   int t;
   int sizepsdvar;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   /* PSDVAR need to be in front of PSDCON! */
   if ( data->nsdpblocks > -1 )
   {
      SCIPerrorMessage("Need to have 'PSDVAR' section before 'PSDCON' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &(data->npsdvars)) != 1 )
      return SCIP_READERROR;

   if ( data->npsdvars < 0 )
   {
      SCIPerrorMessage("Number of psd variables %d should be non-negative!\n", data->npsdvars);
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->npsdvars >= 0 );

   /* loop through different psd variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdpsdvars), data->npsdvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->psdvarsizes), data->npsdvars) );

   for (t = 0; t < data->npsdvars; t++)
   {
      SCIP_CALL( CBFfgets(pfile, linecount) );

      if ( sscanf(CBF_LINE_BUFFER, "%i", &sizepsdvar) == 1 )
      {
         SCIP_Real lb;
         SCIP_Real ub;
         int nscalarvars;
         int cnt = 0;

         if ( sizepsdvar <= 0 )
         {
            SCIPerrorMessage("Size %d of psd variable should be positive!\n", sizepsdvar);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         /* for each psd variable of size n_i create 1/2*n_i*(n_i+1) scalar variables */
         nscalarvars = sizepsdvar * (sizepsdvar + 1) * 0.5;
         data->psdvarsizes[t] = sizepsdvar;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdpsdvars[t]), sizepsdvar) );

         lb = -SCIPinfinity(scip);
         ub = SCIPinfinity(scip);

         for (i = 0; i < sizepsdvar; ++i)
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdpsdvars[t][i]), i+1) );
            for (j = 0; j <= i; ++j)
            {
#ifndef NDEBUG
               snprintfreturn = SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d%d%d", t, i, j);
               assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
               (void)SCIPsnprintf(varname, SCIP_MAXSTRLEN, "y_%d%d%d", t, i, j);
#endif

               SCIP_CALL( SCIPcreateVar(scip, &var, varname, lb, ub, 0.0, SCIP_VARTYPE_CONTINUOUS,
                     TRUE, FALSE, NULL, NULL, NULL, NULL, NULL));/*lint !e732*//*lint !e747*/

               SCIP_CALL( SCIPaddVar(scip, var) );
               data->createdpsdvars[t][i][j] = var;/*lint !e732*//*lint !e747*/
               ++cnt;

               /* release variable for the reader */
               SCIP_CALL( SCIPreleaseVar(scip, &var) );
            }
         }
         assert( cnt == nscalarvars );
      }
      else
         return SCIP_READERROR;
   }

   return SCIP_OKAY;
}

/** reads the number and type of constraints from given CBF-file */
static
SCIP_RETCODE CBFreadCon(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   char consname[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   int nconstypes;
   int nconstypeconss;
   int t;
   int c;
   int cnt = 0;
#ifndef NDEBUG
   int snprintfreturn;
#endif

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i %i", &(data->nconss), &nconstypes) != 2 )
      return SCIP_READERROR;
   else
   {
      if ( data->nconss < 0 )
      {
         SCIPerrorMessage("Number of scalar constraints %d should be non-negative!\n", data->nconss);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      if ( nconstypes < 0 )
      {
         SCIPerrorMessage("Number of conic constraint domains %d should be non-negative!\n", nconstypes);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
      assert( data->nconss >= 0 && nconstypes >= 0 );
   }

   /* loop through different constraint types */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->createdconss), data->nconss) );
   for (t = 0; t < nconstypes; t++)
   {
      SCIP_CALL( CBFfgets(pfile, linecount) );

      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT" %i", CBF_NAME_BUFFER, &nconstypeconss) == 2 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

         if ( nconstypeconss <= 0 )
         {
            SCIPerrorMessage("Number of constraints %d should be positive!\n", nconstypeconss);
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         lhs = -SCIPinfinity(scip);
         rhs = SCIPinfinity(scip);

         if ( strcmp(CBF_NAME_BUFFER, "L+") == 0 )
         {
            lhs = 0.0;
         }
         else if ( strcmp(CBF_NAME_BUFFER, "L-") == 0 )
         {
            rhs = 0.0;
         }
         else if ( strcmp(CBF_NAME_BUFFER, "L=") == 0 )
         {
            lhs = 0.0;
            rhs = 0.0;
         }
         else
         {
            SCIPerrorMessage("CBF-Reader of SCIP-SDP currently only supports linear greater or equal, less or equal and"
               "equality constraints!\n");
            SCIPABORT();
            return SCIP_READERROR; /*lint !e527*/
         }

         /* create corresponding constraints */
         for (c = 0; c < nconstypeconss; c++)
         {
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(consname, SCIP_MAXSTRLEN, "LP_%d", cnt);
            assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
            (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "linear_%d", cnt);
#endif

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname, 0, NULL, NULL, lhs, rhs,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));

            SCIP_CALL( SCIPaddCons(scip, cons) );
            data->createdconss[cnt++] = cons;

            /* release constraint for the reader. */
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
         }
      }
      else
         return SCIP_READERROR;
   }

   if ( cnt != data->nconss )
   {
      SCIPerrorMessage("Total number of constraints for different cone types not equal to total number of constraints!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   return SCIP_OKAY;
}

/** reads integrality conditions from given CBF-file */
static
SCIP_RETCODE CBFreadInt(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   int nintvars;
   int i;
   int v;
   SCIP_Bool infeasible;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdvars == NULL )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'INT' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nintvars) == 1 )
   {
      if ( nintvars >= 0 )
      {
         for (i = 0; i < nintvars; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i", &v) == 1 )
            {
               SCIP_CALL( SCIPchgVarType(scip, data->createdvars[v], SCIP_VARTYPE_INTEGER, &infeasible)   );

               if ( infeasible )
               {
                  SCIPerrorMessage("Infeasibility detected because of integrality of variable %s!\n", SCIPvarGetName(data->createdvars[v]));
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }
            }
            else
            {
               SCIPerrorMessage("Number of integrality constraints %d should be non-negative!\n", nintvars);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }
         }
      }
      else
      {
         SCIPerrorMessage("Number of integrality constraints %d should be non-negative!\n", nintvars);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads SDP-constraint sizes from given CBF-file */
static
SCIP_RETCODE CBFreadPsdcon(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   int b;
   int ncbfsdpblocks;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &ncbfsdpblocks) == 1 )
   {
      if ( ncbfsdpblocks >= 0 )
      {
         /* increase number of sdp blocks by number of psd variables that need to be transformed into a psd constraint
          * with scalar variables */
         if ( data->npsdvars > 0 )
            data->nsdpblocks = ncbfsdpblocks + data->npsdvars;
         else
            data->nsdpblocks = ncbfsdpblocks;

         data->noorigsdpcons = FALSE;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblocksizes), data->nsdpblocks) );

         for (b = 0; b < ncbfsdpblocks; b++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i", &(data->sdpblocksizes[b])) == 1 )
            {
               if ( data->sdpblocksizes[b] <= 0 )
               {
                  SCIPerrorMessage("Size %d of SDP-block %d should be positive!\n", data->sdpblocksizes[b], b);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }
            }
            else
               return SCIP_READERROR;
         }

         for (b = 0; b < data->npsdvars; b++)
            data->sdpblocksizes[ncbfsdpblocks + b] = data->psdvarsizes[b];
      }
      else
      {
         SCIPerrorMessage("Number of SDP-blocks %d should be non-negative!\n", data->nsdpblocks);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   return SCIP_OKAY;
}

/** reads objective values for matrix variables from given CBF-file */
static
SCIP_RETCODE CBFreadObjfcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   SCIP_Real val;
   int nobjcoefs;
   int nzerocoef = 0;
   int i;
   int v;
   int row;
   int col;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdpsdvars == NULL )
   {
      SCIPerrorMessage("Need to have 'PSDVAR' section before 'OBJFCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->npsdvars >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nobjcoefs) == 1 )
   {
      if ( nobjcoefs >= 0 )
      {
         for (i = 0; i < nobjcoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );

            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %lf", &v, &row, &col, &val) == 4 )
            {
               if ( v < 0 || v >= data->npsdvars )
               {
                  SCIPerrorMessage("Given objective coefficient for matrix variable %d which does not exist!\n", v);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( row < 0 || row >= data->psdvarsizes[v] )
               {
                  SCIPerrorMessage("Row index %d of given coefficient for matrix variable %d in objective function is negative or larger than varsize %d!\n",
                     row, v, data->psdvarsizes[v]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( col < 0 || col >= data->psdvarsizes[v] )
               {
                  SCIPerrorMessage("Column index %d of given coefficient for matrix variable %d in objective function is negative or larger than varsize %d!\n",
                     col, v, data->psdvarsizes[v]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  /* make sure matrix is in lower triangular form */
                  if ( row < col )
                  {
                     SCIP_CALL( SCIPchgVarObj(scip, data->createdpsdvars[v][col][row], 2*val) );
                  }
                  else if ( row == col )
                  {
                     SCIP_CALL( SCIPchgVarObj(scip, data->createdpsdvars[v][col][row], val) );
                  }
                  else
                  {
                     SCIP_CALL( SCIPchgVarObj(scip, data->createdpsdvars[v][row][col], 2*val) );
                  }
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of objective coefficients for matrix variables %d should be non-negative!\n", nobjcoefs);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "OBJFCOORD: Found %d coefficients with absolute value less than epsilon = %f.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads objective values for scalar variables from given CBF-file */
static
SCIP_RETCODE CBFreadObjacoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   SCIP_Real val;
   int nobjcoefs;
   int nzerocoef = 0;
   int i;
   int v;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdvars == NULL )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'OBJACOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nobjcoefs) == 1 )
   {
      if ( nobjcoefs >= 0 )
      {
         for (i = 0; i < nobjcoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %lf", &v, &val) == 2 )
            {
               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given objective coefficient for scalar variable %d which does not exist!\n", v);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  SCIP_CALL( SCIPchgVarObj(scip, data->createdvars[v], val) );
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of objective coefficients for scalar variables %d should be non-negative!\n", nobjcoefs);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "OBJACOORD: Found %d coefficients with absolute value less than epsilon = %g.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads matrix variable coefficients from given CBF-file */
static
SCIP_RETCODE CBFreadFcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   SCIP_Real val;
   int nzerocoef = 0;
   int ncoefs;
   int c;
   int i;
   int v;
   int row;
   int col;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdpsdvars == NULL )
   {
      SCIPerrorMessage("Need to have 'PSDVAR' section before 'FCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->npsdvars >= 0 );

   if ( data->createdconss == NULL )
   {
      SCIPerrorMessage("Need to have 'CON' section before 'FCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nconss >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &ncoefs) == 1 )
   {
      if ( ncoefs >= 0 )
      {
         for (i = 0; i < ncoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );

            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %i %lf", &c, &v, &row, &col, &val) == 5 )
            {
               if ( c < 0 || c >= data->nconss )
               {
                  SCIPerrorMessage("Given matrix variable coefficient for constraint %d which does not exist!\n", c);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( v < 0 || v >= data->npsdvars )
               {
                  SCIPerrorMessage("Given coefficient for matrix variable %d which does not exist!\n", v);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( row < 0 || row >= data->psdvarsizes[v] )
               {
                  SCIPerrorMessage("Row index %d of given coefficient for matrix variable %d in scalar constraint %d is negative or larger than varsize %d!\n",
                     row, v, c, data->psdvarsizes[v]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( col < 0 || col >= data->psdvarsizes[v] )
               {
                  SCIPerrorMessage("Column index %d of given coefficient for matrix variable %d in scalar constraint %d is negative or larger than varsize %d!\n",
                     col, v, c, data->psdvarsizes[v]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  /* make sure matrix is in lower triangular form */
                  if ( row < col )
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdpsdvars[v][col][row], 2*val) );/*lint !e732*//*lint !e747*/
                  }
                  else if ( row == col )
                     SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdpsdvars[v][row][col], val) );/*lint !e732*//*lint !e747*/
                  else
                  {
                     SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdpsdvars[v][row][col], 2*val) );/*lint !e732*//*lint !e747*/
                  }
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of matrix variable coefficients %d should be non-negative!\n", ncoefs);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "FCOORD: Found %d coefficients with absolute value less than epsilon = %f.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads linear coefficients from given CBF-file */
static
SCIP_RETCODE CBFreadAcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   SCIP_Real val;
   int nzerocoef = 0;
   int ncoefs;
   int c;
   int i;
   int v;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdvars == NULL )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'ACOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 );

   if ( data->createdconss == NULL )
   {
      SCIPerrorMessage("Need to have 'CON' section before 'ACOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nconss >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &ncoefs) == 1 )
   {
      if ( ncoefs >= 0 )
      {
         for (i = 0; i < ncoefs; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %i %lf", &c, &v, &val) == 3 )
            {
               if ( c < 0 || c >= data->nconss )
               {
                  SCIPerrorMessage("Given linear coefficient for constraint %d which does not exist!\n", c);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }
               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given linear coefficient for variable %d which does not exist!\n", v);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }
               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, data->createdconss[c], data->createdvars[v], val) );/*lint !e732*//*lint !e747*/
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of linear coefficients %d should be non-negative!\n", ncoefs);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "ACOORD: Found %d coefficients with absolute value less than epsilon = %g.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads left- and right-hand sides from given CBF-file */
static
SCIP_RETCODE CBFreadBcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{  /*lint --e{818}*/
   SCIP_Real val;
   int nzerocoef = 0;
   int nsides;
   int c;
   int i;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->createdconss == NULL )
   {
      SCIPerrorMessage("Need to have 'CON' section before 'BCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nconss >= 0 );

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nsides) == 1 )
   {
      if ( nsides >= 0 )
      {
         for (i = 0; i < nsides; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %lf", &c, &val) == 2 )
            {
               if ( c < 0 || c >= data->nconss )
               {
                  SCIPerrorMessage("Given constant part for constraint %d which does not exist!\n", c);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  /* check type */
                  if ( ! SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, data->createdconss[c])) )
                  {
                     /* greater or equal constraint -> left-hand side (minus since we have Ax + b >= 0) */
                     SCIP_CALL( SCIPchgLhsLinear(scip, data->createdconss[c], -val) );
                  }

                  if ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, data->createdconss[c]) ) )
                  {
                     /* less or equal constraint -> right-hand side (minus since we have Ax + b <= 0) */
                     SCIP_CALL( SCIPchgRhsLinear(scip, data->createdconss[c], -val) );
                  }
               }
            }
            else
               return SCIP_READERROR;
         }
      }
      else
      {
         SCIPerrorMessage("Number of left- and right-hand sides %d should be non-negative!\n", nsides);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "BCOORD: Found %d coefficients with absolute value less than epsilon = %f.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads nonzero coefficients of SDP-constraints from given CBF-file */
static
SCIP_RETCODE CBFreadHcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   SCIP_Real val;
   int** sdpvar;
   int nnonz;
   int i;
   int b;
   int v;
   int row;
   int col;
   int firstindforvar;
   int nextindaftervar;
   int nzerocoef = 0;
   int ncbfsdpblocks;
   int nauxnonz;                /* number of nonzeros in each auxiliary sdp block for reformulating matrix variables using
                                 * scalar variables */

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->nsdpblocks < 0 )
   {
      SCIPerrorMessage("Need to have 'PSDCON' section before 'HCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }
   assert( data->nvars >= 0 );

   if ( data->nvars < 0 )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'HCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   /* get number of sdp blocks specified by PSDCON (without auxiliary sdp blocks for reformulating matrix variables
    * using scalar variables), save number of nonzeros needed for the auxiliary sdp blocks in nauxnonz */
   if ( data->npsdvars > 0 )
   {
      ncbfsdpblocks = data->nsdpblocks - data->npsdvars;
      for (i = 0; i < data->npsdvars; i++)
         nauxnonz += data->psdvarsizes[i] * (data->psdvarsizes[i] + 1) * 0.5;
   }
   else
   {
      ncbfsdpblocks = data->nsdpblocks;
      nauxnonz = 0;
   }

   /* initialize sdpnblocknonz with 0 */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblocknonz), data->nsdpblocks) );
   for (b = 0; b < data->nsdpblocks; b++)
      data->sdpnblocknonz[b] = 0;

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &nnonz) == 1 )
   {
      data->nnonz = nnonz + nauxnonz;
      if ( nnonz >= 0 )
      {
         /* allocate memory (nnonz for each block, since we do not yet know the distribution) */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sdpvar, data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval), data->nsdpblocks) );

         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sdpvar[b]), data->nnonz) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow[b]), data->nnonz) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol[b]), data->nnonz) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval[b]), data->nnonz) );
         }

         for (i = 0; i < nnonz; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );

            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %i %lf", &b, &v, &row, &col, &val) == 5 )
            {
               if ( b < 0 || b >= ncbfsdpblocks )
               {
                  SCIPerrorMessage("Given SDP-coefficient for SDP-constraint %d which does not exist!\n", b);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( v < 0 || v >= data->nvars )
               {
                  SCIPerrorMessage("Given SDP-coefficient for variable %d which does not exist!\n", v);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( row < 0 || row >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Row index %d of given SDP coefficient is negative or larger than blocksize %d!\n",
                     row, data->sdpblocksizes[b]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( col < 0 || col >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Column index %d of given SDP coefficient is negative or larger than blocksize %d!\n",
                     col, data->sdpblocksizes[b]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  sdpvar[b][data->sdpnblocknonz[b]] = v;

                  /* make sure matrix is in lower triangular form */
                  if ( col > row )
                  {
                     data->sdprow[b][data->sdpnblocknonz[b]] = col;
                     data->sdpcol[b][data->sdpnblocknonz[b]] = row;
                  }
                  else
                  {
                     data->sdprow[b][data->sdpnblocknonz[b]] = row;
                     data->sdpcol[b][data->sdpnblocknonz[b]] = col;
                  }

                  data->sdpval[b][data->sdpnblocknonz[b]] = val;
                  data->sdpnblocknonz[b]++;
               }
            }
            else
               return SCIP_READERROR;
         }

         /* construct entries for auxiliary sdp blocks (reformulation of psdvars) */
         if ( data->npsdvars > 0 )
         {
            val = 1.0;
            for (v = 0; v < data->npsdvars; v++)
            {
               b = ncbfsdpblocks + v;
               for (row = 0; row < data->psdvarsizes[v]; row++)
               {
                  for (col = 0; col <= row; col++)
                  {
                     sdpvar[b][data->sdpnblocknonz[b]] = data->nvars + 1;
                     data->sdprow[b][data->sdpnblocknonz[b]] = row;
                     data->sdpcol[b][data->sdpnblocknonz[b]] = col;
                     data->sdpval[b][data->sdpnblocknonz[b]] = val;
                     data->sdpnblocknonz[b]++;
                  }
               }
               assert( data->sdpnblocknonz[b] == data->psdvarsizes[v] * (data->psdvarsizes[v] + 1) * 0.5 );
               assert( b == data->nsdpblocks - 1 );
            }
         }

         /* construct pointer arrays */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblockvars), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer), data->nsdpblocks) );

         /* sdp blocks as specified in cbf file in HCOORD */
         for (b = 0; b < ncbfsdpblocks; b++)
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
               SCIP_Bool varused = FALSE;

               firstindforvar = nextindaftervar; /* this variable starts where the last one ended */
               data->nvarnonz[b][data->sdpnblockvars[b]] = 0;

               while (nextindaftervar < data->sdpnblocknonz[b] && sdpvar[b][nextindaftervar] == v) /* get the first index that doesn't belong to this variable */
               {
                  nextindaftervar++;
                  varused = TRUE;
                  data->nvarnonz[b][data->sdpnblockvars[b]]++;
               }

               if ( varused )
               {
                  data->sdpblockvars[b][data->sdpnblockvars[b]] = data->createdvars[v];/*lint !e732*//*lint !e747*/ /* if the variable is used, add it to the vars array */
                  data->rowpointer[b][data->sdpnblockvars[b]] = &(data->sdprow[b][firstindforvar]); /* save a pointer to the first nonzero belonging to this variable */
                  data->colpointer[b][data->sdpnblockvars[b]] = &(data->sdpcol[b][firstindforvar]);
                  data->valpointer[b][data->sdpnblockvars[b]] = &(data->sdpval[b][firstindforvar]);
                  data->sdpnblockvars[b]++;
               }
            }

            assert( nextindaftervar == data->sdpnblocknonz[b] );
         }

         /* auxiliary sdp blocks (reformulation of psd variables) */
         if ( data->npsdvars > 0 )
         {
            assert( data->nsdpblocks == ncbfsdpblocks + data->npsdvars );
            for (v = 0; v < data->npsdvars; v++)
            {
               int varidx = 0;

               /* create the pointer arrays and insert used variables into vars-array */
               b = ncbfsdpblocks + v;

               data->sdpnblockvars[b] = data->sdpnblocknonz[b];

               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars[b]),  data->sdpnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz[b]),  data->sdpnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer[b]),  data->sdpnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer[b]),  data->sdpnblocknonz[b]) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer[b]),  data->sdpnblocknonz[b]) );

               for (row = 0; row < data->psdvarsizes[v]; row++)
               {
                  for (col = 0; col <= row; col++)
                  {
                     data->nvarnonz[b][varidx] = 1;
                     data->sdpblockvars[b][varidx] = data->createdpsdvars[v][row][col];
                     data->rowpointer[b][varidx] = &(data->sdprow[b][varidx]);
                     data->colpointer[b][varidx] = &(data->sdpcol[b][varidx]);
                     data->valpointer[b][varidx] = &(data->sdpval[b][varidx]);
                     varidx++;
                  }
               }
               assert( varidx == data->sdpnblocknonz[b] );
            }
         }

         /* free SDP-var array which is no longer needed */
         for (b = 0; b < data->nsdpblocks; b++)
            SCIPfreeBlockMemoryArray(scip, &(sdpvar[b]), data->nnonz);

         SCIPfreeBlockMemoryArray(scip, &sdpvar, data->nsdpblocks);
      }
      else
      {
         SCIPerrorMessage("Number of nonzero coefficients of SDP-constraints %d should be non-negative!\n", nnonz);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "HCOORD: Found %d coefficients with absolute value less than epsilon = %f.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** reads constant entries of SDP-constraints from given CBF-file */
static
SCIP_RETCODE CBFreadDcoord(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_FILE*            pfile,              /**< file to read from */
   SCIP_Longint*         linecount,          /**< current linecount */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   SCIP_Real val;
   int nzerocoef = 0;
   int constnnonz;
   int b;
   int i;
   int row;
   int col;

   assert( scip != NULL );
   assert( pfile != NULL );
   assert( linecount != NULL );
   assert( data != NULL );

   if ( data->nsdpblocks < 0 )
   {
      SCIPerrorMessage("Need to have 'PSDCON' section before 'DCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   if ( data->nvars < 0 )
   {
      SCIPerrorMessage("Need to have 'VAR' section before 'DCOORD' section!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   SCIP_CALL( CBFfgets(pfile, linecount) );

   if ( sscanf(CBF_LINE_BUFFER, "%i", &constnnonz) == 1 )
   {
      if ( constnnonz < 0 )
      {
         SCIPerrorMessage("Number of constant entries of SDP-constraints %d should be non-negative!\n", constnnonz);
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }

      data->constnnonz = constnnonz;
      if ( constnnonz > 0 )
      {
         /* initialize sdpconstnblocknonz with 0 */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstnblocknonz), data->nsdpblocks) );
         for (b = 0; b < data->nsdpblocks; b++)
            data->sdpconstnblocknonz[b] = 0;

         /* allocate memory (constnnonz for each block, since we do not yet know the distribution) */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstrow), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstcol), data->nsdpblocks) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstval), data->nsdpblocks) );

         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstrow[b]), constnnonz) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstcol[b]), constnnonz) );
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpconstval[b]), constnnonz) );
         }

         for (i = 0; i < constnnonz; i++)
         {
            SCIP_CALL( CBFfgets(pfile, linecount) );
            if ( sscanf(CBF_LINE_BUFFER, "%i %i %i %lf", &b, &row, &col, &val) == 4 )
            {
               if ( b < 0 || b >= data->nsdpblocks )
               {
                  SCIPerrorMessage("Given constant entry for SDP-constraint %d which does not exist!\n", b);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( row < 0 || row >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Row index %d of given constant SDP-entry is negative or larger than blocksize %d!\n",
                     row, data->sdpblocksizes[b]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( col < 0 || col >= data->sdpblocksizes[b] )
               {
                  SCIPerrorMessage("Column index %d of given constant SDP-entry is negative or larger than blocksize %d!\n",
                     col, data->sdpblocksizes[b]);
                  SCIPABORT();
                  return SCIP_READERROR; /*lint !e527*/
               }

               if ( SCIPisZero(scip, val) )
               {
                  ++nzerocoef;
               }
               else
               {
                  /* make sure matrix is in lower triangular form */
                  if ( col > row )
                  {
                     data->sdpconstrow[b][data->sdpconstnblocknonz[b]] = col;
                     data->sdpconstcol[b][data->sdpconstnblocknonz[b]] = row;
                  }
                  else
                  {
                     data->sdpconstrow[b][data->sdpconstnblocknonz[b]] = row;
                     data->sdpconstcol[b][data->sdpconstnblocknonz[b]] = col;
                  }
                  data->sdpconstval[b][data->sdpconstnblocknonz[b]] = -val;
                  data->sdpconstnblocknonz[b]++;
               }
            }
            else
               return SCIP_READERROR;
         }
      }
   }
   else
      return SCIP_READERROR;

   if ( nzerocoef > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
         "DCOORD: Found %d coefficients with absolute value less than epsilon = %g.\n", nzerocoef, SCIPepsilon(scip));
   }

   return SCIP_OKAY;
}

/** frees all data allocated for the CBF-data-struct */
static
SCIP_RETCODE CBFfreeData(
   SCIP*                 scip,               /**< SCIP data structure */
   CBF_DATA*             data                /**< data pointer to save the results in */
   )
{
   SCIP_Bool allocated = FALSE;
   int b = 0;
   int i;
   int t;
   int ncbfsdpblocks;

   assert( scip != NULL );
   assert( data != NULL );

   /* we only allocated memory for the const blocks if there were any nonzeros */
   /* TODO: could also think about saving this in the struct instead, which would cost one bool and save some (unimportant) time here */
   while ( data->sdpconstnblocknonz != NULL && allocated == FALSE && b < data->nsdpblocks)
   {
      if (data->sdpconstnblocknonz[b] > 0)
         allocated = TRUE;
      b++;
   }

   if ( allocated )
   {
      for (b = 0; b < data->nsdpblocks; b++)
      {
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstval[b]), data->constnnonz);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstcol[b]), data->constnnonz);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpconstrow[b]), data->constnnonz);
      }
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstval, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstcol, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstrow, data->nsdpblocks);
      SCIPfreeBlockMemoryArrayNull(scip, &data->sdpconstnblocknonz, data->nsdpblocks);
   }

   /* we only allocated memory for the sdpblocks if there were any nonzeros */
   b = 0;
   while (allocated == FALSE && b < data->nsdpblocks)
   {
      if (data->sdpnblocknonz[b] > 0)
         allocated = TRUE;
      b++;
   }

   if ( allocated )
   {
      /* get number of sdp blocks specified by PSDCON (without auxiliary sdp blocks for reformulating matrix variables
       * using scalar variables), save number of nonzeros needed for the auxiliary sdp blocks in nauxnonz */
      if ( data->npsdvars > 0 )
         ncbfsdpblocks = data->nsdpblocks - data->npsdvars;
      else
         ncbfsdpblocks = data->nsdpblocks;

      if ( data->noorigsdpcons )
      {
         /* no SDP constraints specified in the CBF file! */
         assert( ncbfsdpblocks == 0 );

         for (b = 0; b < data->nsdpblocks; b++)
         {
            SCIPfreeBlockMemoryArrayNull(scip, &(data->valpointer[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->colpointer[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->rowpointer[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpval[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpcol[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdprow[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpblockvars[b]), data->sdpnblocknonz[b]);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->nvarnonz[b]), data->sdpnblocknonz[b]);
         }

         SCIPfreeBlockMemoryArrayNull(scip, &data->valpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->colpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->rowpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpval, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpcol, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdprow, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblockvars, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->nvarnonz, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblockvars, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblocknonz, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblocksizes, data->nsdpblocks);
      }
      else
      {
         /* some SDP constraints specified in the CBF file! */
         assert( ncbfsdpblocks > 0 );

         for (b = 0; b < ncbfsdpblocks; b++)
         {
            SCIPfreeBlockMemoryArrayNull(scip, &(data->valpointer[b]), data->nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->colpointer[b]), data->nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->rowpointer[b]), data->nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpval[b]), data->nnonz);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpcol[b]), data->nnonz);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdprow[b]), data->nnonz);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpblockvars[b]), data->nvars);
            SCIPfreeBlockMemoryArrayNull(scip, &(data->nvarnonz[b]), data->nvars);
         }

         if ( data->npsdvars > 0 )
         {
            for (b = ncbfsdpblocks; b < data->nsdpblocks; b++)
            {
               SCIPfreeBlockMemoryArrayNull(scip, &(data->valpointer[b]), data->sdpnblocknonz[b]);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->colpointer[b]), data->sdpnblocknonz[b]);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->rowpointer[b]), data->sdpnblocknonz[b]);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpval[b]), data->nnonz);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpcol[b]), data->nnonz);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->sdprow[b]), data->nnonz);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->sdpblockvars[b]), data->sdpnblocknonz[b]);
               SCIPfreeBlockMemoryArrayNull(scip, &(data->nvarnonz[b]), data->sdpnblocknonz[b]);
            }
         }

         SCIPfreeBlockMemoryArrayNull(scip, &data->valpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->colpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->rowpointer, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpval, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpcol, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdprow, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblockvars, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->nvarnonz, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblockvars, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpnblocknonz, data->nsdpblocks);
         SCIPfreeBlockMemoryArrayNull(scip, &data->sdpblocksizes, data->nsdpblocks);
      }
   }

   if (data->nconss > 0)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &data->createdconss, data->nconss);
   }

   if (data->nvars > 0)
      SCIPfreeBlockMemoryArrayNull(scip, &data->createdvars, data->nvars);

   if ( data-> npsdvars > 0 )
   {
      for (t = 0; t < data->npsdvars; t++)
      {
         for (i = 0; i < data->psdvarsizes[t]; ++i)
            SCIPfreeBlockMemoryArrayNull(scip, &(data->createdpsdvars[t][i]), i+1);
         SCIPfreeBlockMemoryArrayNull(scip, &(data->createdpsdvars[t]), data->psdvarsizes[t]);
      }

      SCIPfreeBlockMemoryArrayNull(scip, &(data->psdvarsizes), data->npsdvars);
      SCIPfreeBlockMemoryArrayNull(scip, &(data->createdpsdvars), data->npsdvars);
   }

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

   SCIP_CALL( SCIPincludeReaderCbf(scip) );

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadCbf)
{  /*lint --e{715,818}*/
   SCIP_FILE* scipfile;
   SCIP_Longint linecount = 0;
   SCIP_Bool versionread = FALSE;
   SCIP_Bool objread = FALSE;
   CBF_DATA* data;
   int b;

   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   SCIPdebugMsg(scip, "Reading file %s ...\n", filename);

   scipfile = SCIPfopen(filename, "r");

   if ( ! scipfile )
      return SCIP_READERROR;

   SCIP_CALL( SCIPallocBuffer(scip, &data) );
   data->nsdpblocks = -1;
   data->nconss = -1;
   data->nvars = -1;
   data->npsdvars = -1;
   data->constnnonz = 0;
   data->nnonz = 0;
   data->noorigsdpcons= FALSE;

   data->psdvarsizes = NULL;
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

   /* create empty problem */
   SCIP_CALL( SCIPcreateProb(scip, filename, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   while( CBFfgets(scipfile, &linecount) == SCIP_OKAY )
   {
      /* Parse keyword on non-empty lines */
      if ( sscanf(CBF_LINE_BUFFER, CBF_NAME_FORMAT, CBF_NAME_BUFFER) == 1 )
      {
         /* first line should be version number */
         if ( ! versionread )
         {
            if ( strcmp(CBF_NAME_BUFFER, "VER") == 0 )
            {
               int ver;

               SCIP_CALL( CBFfgets(scipfile, &linecount) );

               if ( sscanf(CBF_LINE_BUFFER, "%i", &ver) == 1 )
               {
                  SCIPdebugMsg(scip, "file version %d.\n", ver);
                  if ( ver < 1 )
                  {
                     SCIPerrorMessage("Strange version number %d; need at least version 1.\n", ver);
                     SCIPABORT();
                     return SCIP_READERROR; /*lint !e527*/
                  }
                  else if ( ver > CBF_VERSION_NR )
                  {
                     SCIPerrorMessage("Version %d too new; only supported up to version %d.\n", CBF_VERSION_NR);
                     SCIPABORT();
                     return SCIP_READERROR; /*lint !e527*/
                  }
                  else
                     versionread = TRUE;
               }
               else
                  return SCIP_READERROR;
            }
            else
            {
               SCIPerrorMessage("First keyword should be VER.\n");
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }
         }
         else
         {
            if ( strcmp(CBF_NAME_BUFFER, "OBJSENSE") == 0 )
            {
               SCIPdebugMsg(scip, "Reading OBJSENSE\n");
               SCIP_CALL( CBFreadObjsense(scip, scipfile, &linecount) );
               objread = TRUE;
            }
            else if ( strcmp(CBF_NAME_BUFFER, "VAR") == 0 )
            {
               SCIPdebugMsg(scip, "Reading VAR\n");
               SCIP_CALL( CBFreadVar(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "CON") == 0 )
            {
               SCIPdebugMsg(scip, "Reading CON\n");
               SCIP_CALL( CBFreadCon(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "INT") == 0 )
            {
               SCIPdebugMsg(scip, "Reading INT\n");
               SCIP_CALL( CBFreadInt(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "PSDCON") == 0 )
            {
               SCIPdebugMsg(scip, "Reading PSDCON\n");
               SCIP_CALL( CBFreadPsdcon(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "PSDVAR") == 0 )
            {
               SCIPdebugMsg(scip, "Reading PSDVAR\n");
               SCIP_CALL( CBFreadPsdvar(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "OBJFCOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading OBJFCOORD\n");
               SCIP_CALL( CBFreadObjfcoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "OBJACOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading OBJACOORD\n");
               SCIP_CALL( CBFreadObjacoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "OBJBCOORD") == 0 )
            {
               SCIPerrorMessage("constant part in objective value not supported by SCIP!\n");
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }
            else if ( strcmp(CBF_NAME_BUFFER, "FCOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading FCOORD\n");
               SCIP_CALL( CBFreadFcoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "ACOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading ACOORD\n");
               SCIP_CALL( CBFreadAcoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "BCOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading BCOORD\n");
               SCIP_CALL( CBFreadBcoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "HCOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading HCOORD\n");
               SCIP_CALL( CBFreadHcoord(scip, scipfile, &linecount, data) );
            }
            else if ( strcmp(CBF_NAME_BUFFER, "DCOORD") == 0 )
            {
               SCIPdebugMsg(scip, "Reading DCOORD\n");
               SCIP_CALL( CBFreadDcoord(scip, scipfile, &linecount, data) );
            }
            else
            {
               SCIPerrorMessage("Keyword %s not recognized!\n", CBF_NAME_BUFFER);
               SCIPABORT();
               return SCIP_READERROR; /*lint !e527*/
            }
         }
      }
   }

   /* Psd vars are created using scalar vars and a corresponding psd constraint that gets created in the function
    * CBFreadPsdcon. If there no psd cons in the original problem, then the psd constraint for the reformulation of the
    * original psd vars needs to be created at this point! */
   if ( data->npsdvars > data->nsdpblocks )
   {
      int nauxvars;
      int nnonz;
      int varidx;
      int row;
      int col;
      int val;

      assert( data->nsdpblocks == -1 );
      assert( data->createdpsdvars != NULL );

      data->noorigsdpcons = TRUE;

      data->nsdpblocks = data->npsdvars;
      data->nnonz = 0;

      /* allocate memory for auxiliary sdp blocks */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblocksizes), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblocknonz), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpnblockvars), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer), data->nsdpblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer), data->nsdpblocks) );

      for (b = 0; b < data->nsdpblocks; b++)
      {
         nauxvars = data->psdvarsizes[b] * (data->psdvarsizes[b] + 1) * 0.5;
         nnonz += nauxvars;
         varidx = 0;

         data->sdpblocksizes[b] = data->psdvarsizes[b];
         data->sdpnblocknonz[b] = nauxvars;
         data->sdpnblockvars[b] = nauxvars;

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->nvarnonz[b]), nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpblockvars[b]), nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdprow[b]), nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpcol[b]), nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->sdpval[b]), nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->rowpointer[b]),  nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->colpointer[b]),  nauxvars) );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(data->valpointer[b]),  nauxvars) );

         val = 1;
         for (row = 0; row < data->psdvarsizes[b]; row++)
         {
            for (col = 0; col <= row; col++)
            {
               data->sdprow[b][varidx] = row;
               data->sdpcol[b][varidx] = col;
               data->sdpval[b][varidx] = val;

               data->nvarnonz[b][varidx] = 1;
               data->sdpblockvars[b][varidx] = data->createdpsdvars[b][row][col];
               data->rowpointer[b][varidx] = &(data->sdprow[b][varidx]);
               data->colpointer[b][varidx] = &(data->sdpcol[b][varidx]);
               data->valpointer[b][varidx] = &(data->sdpval[b][varidx]);
               varidx++;
            }
         }
         assert( varidx == nauxvars );
      }
   }

   if ( ! objread )
   {
      SCIPerrorMessage("Keyword OBJSENSE is missing!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

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
      SCIP_CONS* sdpcons;
      char sdpconname[SCIP_MAXSTRLEN];
#ifndef NDEBUG
      int snprintfreturn;
#endif

      assert( data->sdpblocksizes[b] > 0 );
      assert( (data->sdpnblockvars[b] > 0 && data->sdpnblocknonz[b] > 0) || (data->sdpconstnblocknonz[b] > 0) );

#ifndef NDEBUG
      snprintfreturn = SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
      assert( snprintfreturn < SCIP_MAXSTRLEN);
#else
      (void) SCIPsnprintf(sdpconname, SCIP_MAXSTRLEN, "SDP_%d", b);
#endif

      /* special treatment of case without constant PSD blocks */
      if ( data->sdpconstnblocknonz == NULL )
      {
         SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
               data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b], data->valpointer[b],
               data->sdpblockvars[b], 0, NULL, NULL, NULL) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsSdp(scip, &sdpcons, sdpconname, data->sdpnblockvars[b], data->sdpnblocknonz[b],
               data->sdpblocksizes[b], data->nvarnonz[b], data->colpointer[b], data->rowpointer[b], data->valpointer[b],
               data->sdpblockvars[b], data->sdpconstnblocknonz[b], data->sdpconstcol[b], data->sdpconstrow[b],
               data->sdpconstval[b]) );
      }

#ifdef SCIP_MORE_DEBUG
      SCIP_CALL( SCIPprintCons(scip, sdpcons, NULL) );
#endif

      SCIP_CALL( SCIPaddCons(scip, sdpcons) );

      SCIP_CALL( SCIPreleaseCons(scip, &sdpcons) );
   }

   SCIP_CALL( CBFfreeData(scip, data) );
   SCIPfreeBufferNull(scip, &data);

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCbf)
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

   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Writing problem in CBF format to file.\n");
   *result = SCIP_DIDNOTRUN;

   if ( transformed )
   {
      SCIPerrorMessage("CBF reader currently only supports writing original problems!\n");
      SCIPABORT();
      return SCIP_READERROR; /*lint !e527*/
   }

   for (c = 0; c < nconss; c++)
   {
      if ( (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0)
         && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 ) )
      {
         SCIPerrorMessage("CBF reader currently only supports linear and SDP constraints!\n");
         SCIPABORT();
         return SCIP_READERROR; /*lint !e527*/
      }
   }

#ifndef NDEBUG
   for (v = 0; v < nvars; v++)
   {
      assert( SCIPvarGetStatus(vars[v]) == SCIP_VARSTATUS_ORIGINAL );
   }
#endif

   /* write version number */
   SCIPinfoMessage(scip, file, "VER\n%d\n\n", CBF_VERSION_NR);

   /* write objective sense */
   SCIPinfoMessage(scip, file, "OBJSENSE\n%s\n\n", objsense == SCIP_OBJSENSE_MINIMIZE ? "MIN" : "MAX");

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
         varsenses[v] = 1;
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
         varsenses[v] = -1;
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

   /* now determine senses of constraints - possibly disable constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &consssenses, nconss) );
   for (c = 0; c < nconss; c++)
   {
      SCIP_Real lhs;
      SCIP_Real rhs;

      consssenses[c] = 2;

      /* only count linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      lhs = SCIPgetLhsLinear(scip, conss[c]);
      rhs = SCIPgetRhsLinear(scip, conss[c]);

#ifdef CBF_CHECK_NONNEG
      /* check if there are constraints that would determine senses of variables */
      if ( SCIPgetNVarsLinear(scip, conss[c]) == 1 && ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
      {
         /* the nonzero should be a true nonzero */
         assert( ! SCIPisZero(scip, SCIPgetValsLinear(scip, conss[c])[0]) );
         assert( SCIPgetVarsLinear(scip, conss[c]) != NULL );

         v = SCIPvarGetProbindex(SCIPgetVarsLinear(scip, conss[c])[0]);
         assert( 0 <= v && v < nvars );

         if ( SCIPgetValsLinear(scip, conss[c])[0] > 0.0 )
         {
            if ( SCIPisZero(scip, lhs) )
            {
               varsenses[v] = 1;
               continue;
            }
            else if ( SCIPisZero(scip, rhs) )
            {
               varsenses[v] = -1;
               continue;
            }
         }
         else
         {
            if ( SCIPisZero(scip, lhs) )
            {
               varsenses[v] = -1;
               continue;
            }
            else if ( SCIPisZero(scip, rhs) )
            {
               varsenses[v] = 1;
               continue;
            }
         }
      }
#endif

      if ( SCIPisEQ(scip, lhs, rhs) )
      {
         assert( ! SCIPisInfinity(scip, -lhs) );
         assert( ! SCIPisInfinity(scip, rhs) );
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

         if ( ! SCIPisInfinity(scip, -lhs) )
            consssenses[c] = 1;
         else if ( ! SCIPisInfinity(scip, rhs) )
            consssenses[c] = -1;
      }
   }

   /* compute different varsenses */
   nsenses = 0;
   lastsense = 2;
   for (v = 0; v < nvars; v++)
   {
      if ( varsenses[v] != lastsense )
      {
         ++nsenses;
         lastsense = varsenses[v];
      }
   }

   /* write variable senses */
   SCIPinfoMessage(scip, file, "VAR\n%d %d\n", nvars, nsenses);
   lastsense = varsenses[0];
   nsenses = 1;
   for (v = 1; v < nvars; v++)
   {
      if ( varsenses[v] != lastsense )
      {
         if ( lastsense == 0 )
            SCIPinfoMessage(scip, file, "F %d\n", nsenses);
         else if ( lastsense == -1 )
            SCIPinfoMessage(scip, file, "L- %d\n", nsenses);
         else
         {
            assert( lastsense == 1 );
            SCIPinfoMessage(scip, file, "L+ %d\n", nsenses);
         }
         nsenses = 0;
         lastsense = varsenses[v];
      }
      ++nsenses;
   }
   if ( lastsense == 0 )
      SCIPinfoMessage(scip, file, "F %d\n\n", nsenses);
   else if ( lastsense == -1 )
      SCIPinfoMessage(scip, file, "L- %d\n\n", nsenses);
   else
   {
      assert( lastsense == 1 );
      SCIPinfoMessage(scip, file, "L+ %d\n\n", nsenses);
   }

   /* write integrality constraints */
   if ( nbinvars + nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "INT\n%d\n", nbinvars + nintvars);

      for (v = 0; v < nbinvars + nintvars; v++)
      {
         assert( SCIPvarIsIntegral(vars[v]) );
         SCIPinfoMessage(scip, file, "%d\n", v);
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /* compute different consssenses */
   nsenses = 0;
   lastsense = 3;
   i = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( consssenses[c] == 2 )
         continue;

      ++i;
      if ( consssenses[c] != lastsense )
      {
         ++nsenses;
         lastsense = consssenses[c];
      }
   }

   /* write constraint senses */
   SCIPinfoMessage(scip, file, "CON\n%d %d\n", i, nsenses);
   nconsssenses = nsenses;
   c = 0;
   while (c < nconss && consssenses[c] == 2)
      ++c;
   if ( c < nconss )
   {
      lastsense = consssenses[c];
      nsenses = 1;
      ++c;
      for (; c < nconss; ++c)
      {
         if ( consssenses[c] == 2 )
            continue;

         if ( consssenses[c] != lastsense )
         {
            if ( lastsense == 0 )
               SCIPinfoMessage(scip, file, "L= %d\n", nsenses);
            else if ( lastsense == -1 )
               SCIPinfoMessage(scip, file, "L- %d\n", nsenses);
            else
            {
               assert( lastsense == 1 );
               SCIPinfoMessage(scip, file, "L+ %d\n", nsenses);
            }
            nsenses = 0;
            lastsense = consssenses[c];
         }
         ++nsenses;
      }
   }
   if ( lastsense == 0 )
      SCIPinfoMessage(scip, file, "L= %d\n\n", nsenses);
   else if ( lastsense == -1 )
      SCIPinfoMessage(scip, file, "L- %d\n\n", nsenses);
   else
   {
      assert( lastsense == 1 );
      SCIPinfoMessage(scip, file, "L+ %d\n\n", nsenses);
   }

   /* count number of SDP constraints (conshdlrGetNConss doesn't seem to work before transformation) */
   nsdpconss = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
         ++nsdpconss;
   }

   /* write SDP constraints */
   SCIPinfoMessage(scip, file, "PSDCON\n%d\n", nsdpconss);

   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 )
         continue;

      SCIPinfoMessage(scip, file, "%d\n", SCIPconsSdpGetBlocksize(scip, conss[c]));
   }

   SCIPinfoMessage(scip, file, "\n");

   /* count number of nonzero objective coefficients */
   nobjnonz = 0;
   for (v = 0; v < nvars; v++)
   {
      if ( ! SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
         ++nobjnonz;
   }

   /* write objective */
   SCIPinfoMessage(scip, file, "OBJACOORD\n%d\n", nobjnonz);

   for (v = 0; v < nvars; v++)
   {
      SCIP_Real obj;

      obj = SCIPvarGetObj(vars[v]);
      if ( ! SCIPisZero(scip, obj) )
      {
         SCIPinfoMessage(scip, file, "%d %.15g\n", v, obj);
      }
   }
   SCIPinfoMessage(scip, file, "\n");

   /* write coefficients of linear constraints */
   if ( nconsssenses > 0 )
   {
      /* count number of nonzero coefficients in linear constraints */
      nnonz = 0;
      nbnonz = 0;
      for (c = 0; c < nconss; c++)
      {
         if ( consssenses[c] == -1 )
         {
            assert( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) );
            nnonz += SCIPgetNVarsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
               ++nbnonz;
         }
         else if ( consssenses[c] == 1 )
         {
            assert( ! SCIPisInfinity(scip, -SCIPgetLhsLinear(scip, conss[c])) );
            nnonz += SCIPgetNVarsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
               ++nbnonz;
         }
         else if ( consssenses[c] == 0 )
         {
            assert( SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) );
            nnonz += SCIPgetNVarsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
               ++nbnonz;
         }
      }

      /* write linear nonzero coefficients */
      SCIPinfoMessage(scip, file, "ACOORD\n%d\n", nnonz);
      consind = 0;
      for (c = 0; c < nconss; c++)
      {
         if ( consssenses[c] == -1 || consssenses[c] == 0 || consssenses[c] == 1 )
         {
            assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") == 0 );

            linvars = SCIPgetVarsLinear(scip, conss[c]);
            linvals = SCIPgetValsLinear(scip, conss[c]);

            for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
            {
               i = SCIPvarGetProbindex(linvars[v]);
               assert( 0 <= i && i < nvars );
               SCIPinfoMessage(scip, file, "%d %d %.15g\n", consind, i, linvals[v]);
            }
            ++consind;
         }
      }
      SCIPinfoMessage(scip, file, "\n");

      /* write constant part of linear constraints */
      SCIPinfoMessage(scip, file, "BCOORD\n%d\n", nbnonz);
      consind = 0;
      for (c = 0; c < nconss; c++)
      {
         SCIP_Real val;
         if ( consssenses[c] == -1 )
         {
            val = SCIPgetRhsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %.15g\n", consind, -val);
            consind++;
         }
         else if ( consssenses[c] == 1 )
         {
            val = SCIPgetLhsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %.15g\n", consind, -val);
            consind++;
         }
         else if ( consssenses[c] == 0 )
         {
            val = SCIPgetLhsLinear(scip, conss[c]);
            if ( ! SCIPisZero(scip, val) )
               SCIPinfoMessage(scip, file, "%d %.15g\n", consind, -val);
            consind++;
         }
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /* count SDP nonzeros */
   totalsdpnnonz = 0;
   totalsdpconstnnonz = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])),"SDP") != 0 )
         continue;

      SCIP_CALL( SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz) );
      totalsdpnnonz += sdpnnonz;
      totalsdpconstnnonz += sdpconstnnonz;
   }

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

   /* write SDP nonzeros */
   if ( totalsdpnnonz > 0 )
   {
      SCIPinfoMessage(scip, file, "HCOORD\n%d\n", totalsdpnnonz);
      consind = 0;
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])),"SDP") != 0 )
            continue;

         /* initialization for SDPconsSDPGetData-call */
         sdparraylength = totalsdpnnonz;
         sdpconstnnonz = totalsdpconstnnonz;

         SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength, sdpnvarnonz,
               sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow, sdpconstval) );

         assert( sdpconstnnonz <= totalsdpconstnnonz );
         assert( sdparraylength <= totalsdpnnonz);

         for (v = 0; v < sdpnvars; v++)
         {
            for (i = 0; i < sdpnvarnonz[v]; i++)
            {
               int ind;
               ind = SCIPvarGetProbindex(sdpvars[v]);
               assert( 0 <= ind && ind < nvars );
               SCIPinfoMessage(scip, file, "%d %d %d %d %.15g\n", consind, ind, sdprow[v][i], sdpcol[v][i], sdpval[v][i]);
            }
         }
         consind++;
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /* write nonzeros of constant part of SDP constraint */
   if ( totalsdpnnonz > 0 )
   {
      SCIPinfoMessage(scip, file, "DCOORD\n%d\n", totalsdpconstnnonz);
      consind = 0;
      for (c = 0; c < nconss; c++)
      {
         if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 )
            continue;

         /* initialization for SDPconsSDPGetData-call */
         sdparraylength = totalsdpnnonz;
         sdpconstnnonz = totalsdpconstnnonz;

         SCIP_CALL( SCIPconsSdpGetData(scip, conss[c], &sdpnvars, &sdpnnonz, &sdpblocksize, &sdparraylength, sdpnvarnonz,
               sdpcol, sdprow, sdpval, sdpvars, &sdpconstnnonz, sdpconstcol, sdpconstrow, sdpconstval) );

         assert( sdpconstnnonz <= totalsdpconstnnonz );
         assert( sdparraylength <= totalsdpnnonz);

         for (i = 0; i < sdpconstnnonz; i++)
         {
            SCIPinfoMessage(scip, file, "%d %d %d %.15g\n", consind, sdpconstrow[i], sdpconstcol[i], -sdpconstval[i]);
         }
         consind++;
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

/** includes the CBF file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderCbf(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata = NULL;
   SCIP_READER* reader;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert( reader != NULL );

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCbf) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCbf) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCbf) );

   return SCIP_OKAY;
}
