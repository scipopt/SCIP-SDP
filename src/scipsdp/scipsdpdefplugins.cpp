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
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scipsdpdefplugins.cpp
 * @brief  default SCIP-SCP plugins
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIPSDPVERSION              "3.2.0"


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

/* hack to change default parameter values */
#include <scip/paramset.h>

using namespace scip;

/** reset some default parameter values */
SCIP_RETCODE SCIPSDPsetDefaultParams(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PARAM* param;

   /* turn off LP solving - note that the SDP relaxator is on by default */
   param = SCIPgetParam(scip, "lp/solvefreq");
   SCIPparamSetDefaultInt(param, -1);

   param = SCIPgetParam(scip, "lp/cleanuprows");
   SCIPparamSetDefaultBool(param, FALSE);

   param = SCIPgetParam(scip, "lp/cleanuprowsroot");
   SCIPparamSetDefaultBool(param, FALSE);

   /* change default display */
   param = SCIPgetParam(scip, "display/lpiterations/active");
   SCIPparamSetDefaultInt(param, 0);

   param = SCIPgetParam(scip, "display/lpavgiterations/active");
   SCIPparamSetDefaultInt(param, 0);

   param = SCIPgetParam(scip, "display/nfrac/active");
   SCIPparamSetDefaultInt(param, 0);

   param = SCIPgetParam(scip, "display/curcols/active");
   SCIPparamSetDefaultInt(param, 0);

   param = SCIPgetParam(scip, "display/strongbranchs/active");
   SCIPparamSetDefaultInt(param, 0);

   /* Because in the SDP-world there are no warmstarts as for LPs, the main advantage for DFS (that the change in the
    * problem is minimal and therefore the Simplex can continue with the current Basis) is lost and best first search, which
    * provably needs the least number of nodes (see the Dissertation of Tobias Achterberg, the node selection rule with
    * the least number of nodes, allways has to be a best first search), is the optimal choice
    */
   param = SCIPgetParam(scip, "nodeselection/hybridestim/stdpriority");
   SCIPparamSetDefaultInt(param, 1000000);

   param = SCIPgetParam(scip, "nodeselection/hybridestim/maxplungedepth");
   SCIPparamSetDefaultInt(param, 0);

   /* now set parameters to their default value */
   SCIP_CALL( SCIPresetParams(scip) );

   /* The function SCIPparamSetDefaultReal() does not yet exist. We therefore just set the parameter */
   SCIP_CALL( SCIPsetRealParam(scip, "nodeselection/hybridestim/estimweight", 0.0) );

   return SCIP_OKAY;
}

/** callback function to adapt further parameters once changed */
SCIP_DECL_PARAMCHGD(SCIPparamChgdSolvesdps)
{
   int value;

   value = SCIPparamGetInt(param);
   if ( value == 1 )
   {
      /* turn on SDP solving, turn off LP solving */
      SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );
      SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
      SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", FALSE) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Turning on SDP solving, turning off LP solving, cleanuprows(root) = FALSE.\n");
   }
   else
   {
      /* turn off SDP solving, turn on LP solving */
      SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", -1) );
      SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", 1) );
      SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprows", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "lp/cleanuprowsroot", TRUE) );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Turning on LP solving, turning off SDP solving, cleanuprows(root) = TRUE.\n");
   }

   return SCIP_OKAY;
}

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

   /* add parameter to decide whether SDPs or LPs are solved */
   SCIP_CALL( SCIPaddIntParam(scip, "misc/solvesdps", "solve SDPs (1) or LPs (0)", NULL, FALSE, 1, 0, 1, SCIPparamChgdSolvesdps, NULL) );

   /* change default parameter settings */
   SCIP_CALL( SCIPSDPsetDefaultParams(scip) );

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

   return SCIP_OKAY;
}
