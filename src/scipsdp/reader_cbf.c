/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2016 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2016 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_cbf.c
 * @brief  file reader for mixed-integer semidefinite programs in CBF format
 * @author Tristan Gally
 * @author Henrik A. Friberg TODO
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>                      /* for strcmp */

#include "scipsdp/reader_cbf.h"
#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/cons_sdp.h"
#include "scip/cons_linear.h"


#define READER_NAME             "cbfreader"
#define READER_DESC             "file reader and writer for MISDPs in cbf format"
#define READER_EXTENSION        "cbf"

#define CBF_VERSION_NR         1         /* version numberfor CBF format */
#define CBF_CHECK_NONNEG       TRUE      /* when writing: check linear constraints and move nonnegativity(-positivity)
                                          * constraints to definition of variables (which are now defined in non-negative
                                          * orthant)
                                          * TODO: currently doesn't work for ranged rows (which are not created by sdpa
                                          * reader) */


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyCbf)
{
   assert(scip != NULL);

   SCIP_CALL( SCIPincludeReaderCbf(scip) );

   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_READERFREE(readerFreeCbf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cbf reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerFreeCbf NULL
#endif


/** problem reading method of reader */
#if 0
static
SCIP_DECL_READERREAD(readerReadCbf)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of cbf reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerReadCbf NULL
#endif


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCbf)
{  /*lint --e{715}*/
   SdpVarmapper* varmapper;
   SCIP_VAR** linvars;
   SCIP_Real* linvals;
   SCIP_VAR** sdpvars;
   int nsdpconss;
   int sdpnvars;
   int sdpnnonz;
   int totalsdpnnonz;
   int sdpblocksize;
   int sdparraylength;
   int* sdpnvarnonz;
   int** sdpcol;
   int** sdprow;
   SCIP_Real** sdpval;
   int totalsdpconstnnonz;
   int sdpconstnnonz;
   int* sdpconstcol;
   int* sdpconstrow;
   SCIP_Real* sdpconstval;
   int c;
   int i;
   int v;
   int consind;
   int neqconss;
   int nleqconss;
   int ngeqconss;
   int nobjnonz;
   int nnonz;
   int nbnonz;
#ifdef CBF_CHECK_NONNEG
   int nconsdisabled;
   SCIP_Bool* consdisabled;
   int nposorth;
   SCIP_Bool* posorth;
   int nnegorth;
   SCIP_Bool* negorth;
#endif

   assert( scip != NULL );

   SCIPdebugMsg(scip, "Writing problem in CBF format to %s\n", file);
   *result = SCIP_DIDNOTRUN;

   if ( transformed )
   {
      SCIPerrorMessage("CBF reader currently only supports writing original problems!\n");
      SCIPABORT(); /*lint --e{527}*/
   }

   for (c = 0; c < nconss; c++)
   {
      if ( (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0)
            && (strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") != 0 ) )
      {
         SCIPerrorMessage("CBF reader currently only supports linear and SDP constraints!\n");
         SCIPABORT(); /*lint --e{527}*/
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

#ifdef CBF_CHECK_NONNEG
   /* check for constraints with only a single non-zero and lhs/rhs zero, in that case disable the constraint
    * and force variable to be in positive or negative orthant
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &consdisabled, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &posorth, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &negorth, nvars) );

   for (v = 0; v < nvars; v++)
   {
      posorth[v] = FALSE;
      negorth[v] = FALSE;
   }
   nposorth = 0;
   nnegorth = 0;
   nconsdisabled = 0;

   for (c = 0; c < nconss; c++)
   {
      consdisabled[c] = 0;

      /* only check linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      /* only check constraints with exactly one nonzero */
      if ( SCIPgetNVarsLinear(scip, conss[c]) != 1 )
         continue;

      /* the nonzero should be a true nonzero */
      assert( ! SCIPisZero(scip, SCIPgetValsLinear(scip, conss[c])[0]) );

      /* the variable should not be fixed */
      assert( SCIPisLT(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) );

      if ( SCIPgetValsLinear(scip, conss[c])[0] > 0.0 )
      {
         if ( SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            nconsdisabled++;
            posorth[v] = TRUE;
            nposorth++;

            if ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) )
            {
               SCIPerrorMessage("Detection of non-negativity constraints currently not supported for ranged rows!\n");
               SCIPABORT(); /*lint --e{527}*/
            }
         }
         else if ( SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            nconsdisabled++;
            negorth[v] = TRUE;
            nnegorth++;
         }
      }
      else
      {
         if ( SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            nconsdisabled++;
            negorth[v] = TRUE;
            nnegorth++;
         }
         else if ( SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])) )
         {
            consdisabled[c] = 1;
            nconsdisabled++;
            posorth[v] = TRUE;
            nposorth++;
         }
      }
   }

   assert( (0 <= nposorth) && (0 <= nnegorth) && (nposorth + nnegorth <= nvars) );

   /* prepare varmapper that maps SCIP variables to indices for CBF format (3/4 is the load factor java uses) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(varmapper), (int) ceil(1.33 * nvars)) );
   /* in this case we first need to add all variables in the positive orthant, then all in the negative orthant
    * and finally the free variables, since they will be input and counted in CBF this way
    */
   for (v = 0; v < nvars; v++)
   {
      if ( posorth[v] )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }
   for (v = 0; v < nvars; v++)
   {
      if ( negorth[v] )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }
   for (v = 0; v < nvars; v++)
   {
      if ( ( ! posorth[v] ) && ( ! negorth[v] ) )
      {
         SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, 1, &(vars[v])) );
      }
   }

   /* write variables */
   if ( (nposorth == 0) && (nnegorth == 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 1\nF %d\n\n", nvars, nvars);
   else if ( (nposorth > 0) && (nnegorth == 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 2\nL+ %d\nF %d\n\n", nvars, nposorth, nvars - nposorth);
   else if ( (nposorth == 0) && (nnegorth > 0) )
      SCIPinfoMessage(scip, file, "VAR\n%d 2\nL- %d\nF %d\n\n", nvars, nposorth, nvars - nposorth);
   else
      SCIPinfoMessage(scip, file, "VAR\n%d 3\nL+ %d\nL- %d\nF %d\n\n", nvars, nposorth, nnegorth, nvars - nposorth - nnegorth);

#else
   /* prepare varmapper that maps SCIP variables to indices for CBF format (3/4 is the load factor java uses) */
   SCIP_CALL( SCIPsdpVarmapperCreate(scip, &(varmapper), (int) ceil(1.33 * nvars)) );
   SCIP_CALL( SCIPsdpVarmapperAddVars(scip, varmapper, nvars, vars) );

   /* write variables */
   SCIPinfoMessage(scip, file, "VAR\n%d 1\nF %d\n\n", nvars, nvars);
#endif

   /* write integrality constraints */
   if ( nbinvars + nintvars > 0 )
   {
      SCIPinfoMessage(scip, file, "INT\n%d\n", nbinvars + nintvars);

      for (v = 0; v < nbinvars + nintvars; v++)
      {
         assert( SCIPvarIsIntegral(vars[v]) );
         SCIPinfoMessage(scip, file, "%d\n", SCIPsdpVarmapperGetSdpIndex(varmapper, vars[v]));
      }

      SCIPinfoMessage(scip, file, "\n");
   }

   /* count number of equality/geq/leq constraints */
   neqconss = 0;
   ngeqconss = 0;
   nleqconss = 0;
   for (c = 0; c < nconss; c++)
   {
      /* only count linear constraints */
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      if ( ( ! SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))) &&
            ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) ) &&
            SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         neqconss++;
      else
      {
         if ( ! SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c])) )
            ngeqconss++;
         if ( ! SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c])) )
            nleqconss++;
      }
   }

   /* write constraints */
   if ( neqconss + ngeqconss + nleqconss > 0 )
   {
      if ( ((neqconss > 0) && (ngeqconss == 0) && (nleqconss == 0))
            || ((neqconss == 0) && (ngeqconss > 0) && (nleqconss == 0))
            || ((neqconss == 0) && (ngeqconss == 0) && (nleqconss > 0)) )
         SCIPinfoMessage(scip, file, "CON\n%d 1\n", neqconss + ngeqconss + nleqconss);
      else if ( ((neqconss > 0) && (ngeqconss > 0) && (nleqconss == 0))
            || ((neqconss > 0) && (ngeqconss == 0) && (nleqconss > 0))
            || ((neqconss == 0) && (ngeqconss > 0) && (nleqconss > 0)) )
         SCIPinfoMessage(scip, file, "CON\n%d 2\n", neqconss + ngeqconss + nleqconss);
      else
      {
         assert( (neqconss > 0) && (ngeqconss > 0) && (nleqconss > 0) );
         SCIPinfoMessage(scip, file, "CON\n%d 3\n", neqconss + ngeqconss + nleqconss);
      }

      if ( ngeqconss > 0 )
         SCIPinfoMessage(scip, file, "L+ %d", ngeqconss);
      if ( nleqconss > 0 )
         SCIPinfoMessage(scip, file, "L- %d", nleqconss);
      if ( neqconss > 0 )
         SCIPinfoMessage(scip, file, "L= %d", neqconss);

      SCIPinfoMessage(scip, file, "\n\n");
   }

   /* count number of SDP constraints (conshdlrGetNConss doesn't seem to work before transformation) */
   nsdpconss = 0;
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "SDP") == 0 )
         nsdpconss++;
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
         nobjnonz++;
   }

   /* write objective */
   SCIPinfoMessage(scip, file, "OBJACOORD\n%d\n", nobjnonz);

   for (v = 0; v < nvars; v++)
   {
      if ( ! SCIPisZero(scip, SCIPvarGetObj(vars[v])) )
      {
         SCIPinfoMessage(scip, file, "%d %f\n", SCIPsdpVarmapperGetSdpIndex(varmapper, vars[v]), SCIPvarGetObj(vars[v]));
      }
   }

   SCIPinfoMessage(scip, file, "\n");

   /* count number of nonzero coefficients in linear constraints */
   nnonz = 0;
   nbnonz = 0;

   /* first iterate over all greater or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      nnonz += SCIPgetNVarsLinear(scip, conss[c]);
      nbnonz += ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])));
   }
   /* iterate over all less or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      nnonz += SCIPgetNVarsLinear(scip, conss[c]);
      nbnonz += ( ! SCIPisZero(scip, SCIPgetRhsLinear(scip, conss[c])));
   }
   /* finally iterate over all equality constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      nnonz += SCIPgetNVarsLinear(scip, conss[c]);
      nbnonz += ( ! SCIPisZero(scip, SCIPgetLhsLinear(scip, conss[c])));
   }

   /* write linear nonzero coefficients */
   SCIPinfoMessage(scip, file, "ACOORD\n%d\n", nnonz);
   consind = 0;

   /* first iterate over all greater or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      linvars = SCIPgetVarsLinear(scip, conss[c]);
      linvals = SCIPgetValsLinear(scip, conss[c]);

      for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
      {
         SCIPinfoMessage(scip, file, "%d %d %f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
      }
      consind++;
   }
   /* iterate over all less or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      linvars = SCIPgetVarsLinear(scip, conss[c]);
      linvals = SCIPgetValsLinear(scip, conss[c]);

      for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
      {
         SCIPinfoMessage(scip, file, "%d %d %f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
      }
      consind++;
   }
   /* finally iterate over all equality constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      linvars = SCIPgetVarsLinear(scip, conss[c]);
      linvals = SCIPgetValsLinear(scip, conss[c]);

      for (v = 0; v < SCIPgetNVarsLinear(scip, conss[c]); v++)
      {
         SCIPinfoMessage(scip, file, "%d %d %f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, linvars[v]), linvals[v]);
      }
      consind++;
   }

   SCIPinfoMessage(scip, file, "\n");

   /* write constant part of linear constraints */
   SCIPinfoMessage(scip, file, "BCOORD\n%d\n", nbnonz);
   consind = 0;

   /* first iterate over all greater or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      SCIPinfoMessage(scip, file, "%d %f\n", consind, SCIPgetLhsLinear(scip, conss[c]));
      consind++;
   }
   /* iterate over all less or equal constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c])) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      SCIPinfoMessage(scip, file, "%d %f\n", consind, SCIPgetRhsLinear(scip, conss[c]));
      consind++;
   }
   /* finally iterate over all equality constraints */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), "linear") != 0 )
         continue;

      if ( SCIPisInfinity(scip, -1 * SCIPgetLhsLinear(scip, conss[c]))
            || SCIPisInfinity(scip, SCIPgetRhsLinear(scip, conss[c]))
            || ( ! SCIPisEQ(scip, SCIPgetLhsLinear(scip, conss[c]), SCIPgetRhsLinear(scip, conss[c]))) )
         continue;
#ifdef CBF_CHECK_NONNEG
      if ( consdisabled[c] )
         continue;
#endif

      SCIPinfoMessage(scip, file, "%d %f\n", consind, SCIPgetLhsLinear(scip, conss[c]));
      consind++;
   }
   SCIPinfoMessage(scip, file, "\n");

   /*count SDP nonzeros */
   for (c = 0; c < nconss; c++)
   {
      if ( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])),"SDP") != 0 )
         continue;

      SCIPconsSdpGetNNonz(scip, conss[c], &sdpnnonz, &sdpconstnnonz);

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

      for (v = 0; v < nvars; v++)
      {
         for (i = 0; i < sdpnvarnonz[v]; i++)
         {
            SCIPinfoMessage(scip, file, "%d %d %d %d %f\n", consind, SCIPsdpVarmapperGetSdpIndex(varmapper, sdpvars[v]),
                  sdprow[v][i], sdpcol[v][i], sdpval[v][i]);
         }
      }
      consind++;
   }
   SCIPinfoMessage(scip, file, "\n");

   /* write nonzeros of constant part of SDP constraint */
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
         SCIPinfoMessage(scip, file, "%d %d %d %f\n", consind, sdpconstrow[i], sdpconstcol[i], sdpconstval[i]);
      }
      consind++;
   }

   SCIPfreeBufferArray(scip, &sdpconstval);
   SCIPfreeBufferArray(scip, &sdpconstrow);
   SCIPfreeBufferArray(scip, &sdpconstcol);
   SCIPfreeBufferArray(scip, &sdpvars);
   SCIPfreeBufferArray(scip, &sdpval);
   SCIPfreeBufferArray(scip, &sdprow);
   SCIPfreeBufferArray(scip, &sdpcol);
   SCIPfreeBufferArray(scip, &sdpnvarnonz);
   SCIP_CALL( SCIPsdpVarmapperFree(scip, &varmapper) );
#ifdef CBF_CHECK_NONNEG
   SCIPfreeBufferArray(scip, &negorth);
   SCIPfreeBufferArray(scip, &posorth);
   SCIPfreeBufferArray(scip, &consdisabled);
#endif

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
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyCbf) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadCbf) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteCbf) );

   return SCIP_OKAY;
}
