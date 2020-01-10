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

/**@file   scipsdpdefplugins.cpp
 * @brief  default SCIP-SCP plugins
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIPSDPVERSION              "3.1.2"


#include "scipsdp/scipsdpdefplugins.h"
#include "scip/scipdefplugins.h"

#include "objscip/objscipdefplugins.h"

#include "cons_sdp.h"
#include "scipsdp/cons_savesdpsol.h"
#include "cons_savedsdpsettings.h"
#include "relax_sdp.h"
#include "objreader_sdpa.h"
#include "objreader_sdpaind.h"
#include "reader_cbf.h"
#include "prop_sdpredcost.h"
#include "disp_sdpiterations.h"
#include "disp_sdpavgiterations.h"
#include "disp_sdpfastsettings.h"
#include "disp_sdppenalty.h"
#include "disp_sdpunsolved.h"
#include "branch_sdpmostfrac.h"
#include "branch_sdpmostinf.h"
#include "branch_sdpobjective.h"
#include "branch_sdpinfobjective.h"
#include "heur_sdpfracdiving.h"
#include "heur_sdprand.h"
#include "prop_sdpobbt.h"
#include "prop_companalcent.h"
#include "scipsdpgithash.c"
#include "table_relaxsdp.h"
#include "table_sdpsolversuccess.h"
#include "table_slater.h"

/* hack to allow to change the name of the dialog without needing to copy everything */
#include "scip/struct_dialog.h"

using namespace scip;

/** includes default SCIP-SDP plugins */
SCIP_RETCODE SCIPSDPincludeDefaultPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   char scipsdpname[SCIP_MAXSTRLEN];
   char scipsdpdesc[SCIP_MAXSTRLEN];
   SCIP_DIALOG* dialog;

   /* add description */
   (void) SCIPsnprintf(scipsdpname, SCIP_MAXSTRLEN, "SCIP-SDP %s", SCIPSDPVERSION);
   (void) SCIPsnprintf(scipsdpdesc, SCIP_MAXSTRLEN, "Mixed Integer Semidefinite Programming Plugin for SCIP "
         "[GitHash: %s] (www.opt.tu-darmstadt.de/scipsdp/)", SCIPSDP_GITHASH);
   SCIP_CALL( SCIPincludeExternalCodeInformation(scip, scipsdpname, scipsdpdesc) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include new plugins */
   SCIP_CALL( SCIPincludeObjReader(scip, new ObjReaderSDPAind(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjReader(scip, new ObjReaderSDPA(scip), TRUE) );
   SCIP_CALL( SCIPincludeReaderCbf(scip) );
   SCIP_CALL( SCIPincludeConshdlrSdp(scip) );
   SCIP_CALL( SCIPincludeConshdlrSdpRank1(scip) );
   SCIP_CALL( SCIPincludeConshdlrSavesdpsol(scip) );
   SCIP_CALL( SCIPincludeConshdlrSavedsdpsettings(scip) );
   SCIP_CALL( SCIPincludeRelaxSdp(scip) );
   SCIP_CALL( SCIPincludePropSdpredcost(scip) );
   SCIP_CALL( SCIPincludeBranchruleSdpmostfrac(scip) );
   SCIP_CALL( SCIPincludeBranchruleSdpmostinf(scip) );
   SCIP_CALL( SCIPincludeBranchruleSdpobjective(scip) );
   SCIP_CALL( SCIPincludeBranchruleSdpinfobjective(scip) );
   SCIP_CALL( SCIPincludeHeurSdpFracdiving(scip) );
   SCIP_CALL( SCIPincludeHeurSdpRand(scip) );
   SCIP_CALL( SCIPincludePropSdpObbt(scip) );
   SCIP_CALL( SCIPincludePropCompAnalCent(scip) );

   /* change name of dialog */
   dialog = SCIPgetRootDialog(scip);
   BMSfreeMemoryArrayNull(&dialog->name);
   SCIP_ALLOC( BMSallocMemoryArray(&dialog->name, 9) );
   (void) SCIPstrncpy(dialog->name, "SCIP-SDP", 9);

   /* include displays */
   SCIP_CALL( SCIPincludeDispSdpiterations(scip) );
   SCIP_CALL( SCIPincludeDispSdpavgiterations(scip) );
   SCIP_CALL( SCIPincludeDispSdpfastsettings(scip) );
   SCIP_CALL( SCIPincludeDispSdppenalty(scip) );
   SCIP_CALL( SCIPincludeDispSdpunsolved(scip) );

   /* include tables */
   SCIP_CALL( SCIPincludeTableRelaxSdp(scip) );
   SCIP_CALL( SCIPincludeTableSdpSolverSuccess(scip) );
   SCIP_CALL( SCIPincludeTableSlater(scip) );

   /* set clocktype to walltime to not add multiple threads together */
   SCIP_CALL( SCIPsetIntParam(scip, "timing/clocktype", 2) );

   /* Choose between LP and SDP relaxations */
   SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/lpiterations/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/lpavgiterations/active", 0) );

   /* display numerical problems in SDPs instead of current columns and strong branching */
   SCIP_CALL( SCIPsetIntParam(scip, "display/nfrac/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/curcols/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/strongbranchs/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/sdpfastsettings/active", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "display/sdppenalty/active", 0) );

   /* display SDP statistics instead of default relaxator statistics */
   SCIP_CALL( SCIPsetBoolParam(scip, "table/relaxator/active", FALSE) );

   /* change epsilons for numerical stability */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/epsilon", 1e-9) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/sumepsilon", 1e-6) );
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-6) );

   /* parameters for separation */
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE) );
   SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE) );

   /* Because in the SDP-world there are no warmstarts as for LPs, the main advantage for DFS (that the change in the
    * problem is minimal and therefore the Simplex can continue with the current Basis) is lost and best first search, which
    * provably needs the least number of nodes (see the Dissertation of Tobias Achterberg, the node selection rule with
    * the least number of nodes, allways has to be a best first search), is the optimal choice
    */
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/stdpriority", 1000000) );
   SCIP_CALL( SCIPsetIntParam(scip, "nodeselection/hybridestim/maxplungedepth", 0) );
   SCIP_CALL( SCIPsetRealParam(scip, "nodeselection/hybridestim/estimweight", 0.0) );

   return SCIP_OKAY;
}