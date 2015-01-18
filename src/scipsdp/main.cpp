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
 * @brief  main file for solving MISDPs
 * @author Sonja Mars
 * @author Tristan Gally
 */

#include "objscip/objscipdefplugins.h"

#include "cons_sdp.h"
#include "relax_sdp.h"
#include "objreader_sdpa.h"
#include "prop_sdpredcost.h"
#include "disp_sdpiterations.h"
#include "disp_sdpavgiterations.h"

using namespace scip;

/** run scip and set some parameters */
static
SCIP_RETCODE runSCIP(
   int                   argc,               /**< number of command line arguments */
   char**                argv                /**< pointer to command line arguments */
   )
{
   SCIP* scip = NULL;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include new plugins */
   SCIP_CALL( SCIPincludeObjReader(scip, new ObjReaderSDPA(scip), TRUE) );
   SCIP_CALL( SCIPincludeConshdlrSdp(scip) );
   SCIP_CALL( SCIPincludeRelaxSDP(scip) );
   SCIP_CALL( SCIPincludePropSdpredcost(scip) );

   /* add parameter for SDP solver */
   SCIP_CALL( SCIPaddStringParam(scip, "sdpsolver", "SDP solver", 0, FALSE, "dsdp" , 0, 0) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* disable subscips - for the meantime */
   SCIP_CALL( SCIPsetSubscipsOff(scip, TRUE) );


   /* change certain paramters: */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );

   /* Choose between LP and SDP relaxations */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );
   SCIP_CALL( SCIPincludeDispSdpiterations(scip) );
   SCIP_CALL( SCIPincludeDispSdpavgiterations(scip) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/lpiterations/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/lpavgiterations/active", 0) );

   /* change epsilons for numerical stability */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-6) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-4) );

   SCIP_CALL( SCIPsetStringParam(scip, "sdpsolver", "dsdp") );

   /* parameters for separation */
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE) );
   SCIP_CALL( SCIPsetIntParam(scip, "lp/rowagelimit", 10) );

   /* maximum age a cut can reach before it is deleted from the global cut pool, or -1 to keep all cuts */
   SCIP_CALL( SCIPsetIntParam(scip, "separating/cutagelimit", 10) );

   SCIP_CALL( SCIPsetIntParam(scip, "separating/maxrounds", 20) );

   /* Parameters for node selection */

   /* Because in the SDP-world there are no warmstarts as for LPs, the main advantage for DFS (that the change in the
    * problem is minimal and therefore the Simplex can continue with the current Basis) is lost and best first search, which
    * provably needs the least number of nodes (see the Dissertation of Tobias Achterberg, the node selection rule with
    * the least number of nodes, allways has to be a best first search), is the optimal choice
    */
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/stdpriority", 1000000) );
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/maxplungedepth", 0) );
   SCIP_CALL( SCIPsetRealParam(scip, "nodeselection/hybridestim/estimweight", 0.0) );

   /* This theoretically leaves mostinf as the one with the highest priority, but this doesn't work properly for SDP (as do the SCIP implementations of mostinf
    * or leastinf), and just chooses the external candidate with the highest priority, and these priorities are chosen by fractionality, so what is done in
    * the end is a most infeasible branching
    */
   //SCIP_CALL( SCIPsetIntParam(scip, "branching/pscost/priority",-2000000));
   SCIP_CALL( SCIPsetIntParam(scip, "branching/relpscost/priority",-2000000) );
   SCIP_CALL( SCIPsetIntParam(scip, "branching/mostinf/priority",2000000) );

   /* turn off int-obj ????? */
   SCIP_CALL( SCIPsetIntParam(scip, "separating/intobj/freq", -1) );

   /* run interactive shell */
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, "scip.set") );

   /* deinitialization */
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main function */
int main (
   int                   argc,               /**< number of command line arguments */
   char**                argv                /**< pointer to command line arguments */
   )
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
