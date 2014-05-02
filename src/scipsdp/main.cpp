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

/**@file   main.cpp
 * @brief  driver-file for solving MISDPs
 * @author Sonja Mars
 */
// standard library includes
#include <string>

#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

#include "objconshdlr_sdp.h"
#include "relax_sdp.h"
#include "objreader_sdpa.h"

using namespace scip;

/**run scip and set some parameters*/
static
SCIP_RETCODE runSCIP(
   int argc,
   char** argv)
{

   printf("Starting solver.\n");
   SCIP* scip = NULL;


   SCIP_CALL( SCIPcreate(&scip) );
   SCIPprintVersion(scip, NULL);

   //include new plugins
   SCIP_CALL( SCIPincludeObjReader(scip, new ObjReaderSDPA(scip), TRUE) );

   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ObjConshdlrSdp(scip), TRUE) );

   SCIP_CALL( SCIPincludeRelaxSDP(scip) );

   const char* name = "sdpsolver";
   const char * 	desc = "which sdpsolver should be called";

   SCIP_PARAMDATA* 	paramdata = NULL;

   SCIP_CALL( SCIPaddStringParam	(scip, name, desc, NULL, FALSE, "dsdp" , NULL, paramdata)	);

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );

   //Choose between LP and SDP relaxations
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip,"relaxing/SDP/freq", 1) );

   //Do some stuff to be numerically stable
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-6) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-4) );

   SCIP_CALL( SCIPsetStringParam(scip, "sdpsolver", "dsdp") );

   //parameters for separation
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 10) );

   //# maximum age a cut can reach before it is deleted from the global cut pool, or -1 to keep all cuts
   //# [type: int, range: [-1,2147483647], default: 100]
   SCIP_CALL( SCIPsetIntParam(scip, "separating/cutagelimit", 10));


   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 20));

   //Parameters for node selection
/*   int relaxfreq;
   SCIP_CALL(SCIPgetIntParam(scip, "relaxing/SDPRelax/freq", &relaxfreq));
   if (relaxfreq==1) {                                                           also doesn't work, see line 81*/
      SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/stdpriority", 1000000));
      SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/maxplungedepth", 0));
      SCIP_CALL( SCIPsetRealParam(scip, "nodeselection/hybridestim/estimweight", 0.0));

      SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority",-2000000));
      SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/priority",-2000000));
  // }

   //turn off int-obj
   SCIP_CALL( SCIPsetIntParam(scip, "separating/intobj/freq", -1));

   //read parameter file
   if (argc > 2 )
   {
      SCIP_CALL( SCIPreadParams(scip,argv[2]));
   }




   printf("\n read problem\n");
   printf("============\n");

   if (argc > 1)
   {
      SCIP_CALL( SCIPreadProb(scip, argv[1], NULL) );
   }
   else
   {
      printf("Wrong call of main, not enough arguments. Call using: \
      ./main instancefile.dat-s    or    ./main instancefile.dat-s parameterfile.set");
   }

   /* solve problem */
   printf("\nsolve problem\n");
   printf("=============\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPprintStatistics(scip, NULL));
   printf("\nprimal solution:\n" );
   printf("================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );

   FILE* solfile;
   std::string solfilename = "Solution.sol";
   if(argc > 1)
   {
      solfilename = std::string(argv[1]) + "_misdp.sol";
   }


   solfile = fopen(solfilename.c_str(), "w");
   if( solfile == NULL )
   {
      SCIPerrorMessage("error creating file <%s>\n", solfilename.c_str());
   }
   else
   {
      SCIP_CALL( SCIPprintBestSol(scip, solfile, FALSE) );
      fclose(solfile);
   }

   /********************
    * Deinitialization *
    ********************/
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/**main function */
int main (
   int argc,
   char** argv)
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
