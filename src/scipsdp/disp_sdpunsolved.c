/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_sdpunsolved.c
 * @brief  Column to display the percentage of SDP-relaxations that could not be solved even using a penalty formulation
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "disp_sdpunsolved.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define DISP_NAME             "sdpunsolved"
#define DISP_DESC             "percentage of unsolved SDPs"
#define DISP_HEADER           "SDP  uns"
#define DISP_WIDTH            8              /**< the width of the display column */
#define DISP_PRIORITY         40000          /**< the priority of the display column */
#define DISP_POSITION         1475           /**< the relative position of the display column */
#define DISP_STRIPLINE        TRUE           /**< default for displaying column separated with a line from its right neighbor */




/*
 * Data structures
 */

/** display column data */
struct SCIP_DispData
{
   SCIP_RELAX*           relaxSDP;           /**< pointer to the SDP relaxator whose fastsettings should be displayed */
};



/*
 * Callback methods of display column
 */

/** copy method for dialog plugins (called when SCIP copies plugins) */
static
SCIP_DECL_DISPCOPY(dispCopySdpunsolved)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(disp != NULL);

   /* call inclusion method of dialog */
   SCIP_CALL( SCIPincludeDispSdpunsolved(scip) );

   return SCIP_OKAY;
}

/** destructor of display column to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DISPFREE(dispFreeSdpunsolved)
{
   SCIP_DISPDATA* dispdata;

   assert( scip != NULL );
   assert( disp != NULL );
   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );

   SCIPfreeMemory(scip, &dispdata);
   SCIPdispSetData(disp, NULL);

   return SCIP_OKAY;
}

/** solving process initialization method of display column (called when branch and bound process is about to begin) */
static
SCIP_DECL_DISPINITSOL(dispInitsolSdpunsolved)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   assert( disp != NULL );
   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );

   dispdata->relaxSDP = SCIPfindRelax(scip, "SDP");

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputSdpunsolved)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   assert( scip != NULL );
   assert( disp != NULL );

   dispdata = SCIPdispGetData(disp);

   assert( dispdata != NULL );

   if ( dispdata->relaxSDP != NULL )
   {
      if ( SCIPrelaxSdpGetNSdpCalls(dispdata->relaxSDP) == 0 )
      {
         SCIPinfoMessage(scip, file, "   --   ");
      }
      else
      {
         SCIP_Real fastpercent;
         fastpercent = (SCIP_Real) SCIPrelaxSdpGetNSdpUnsolved(dispdata->relaxSDP) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(dispdata->relaxSDP);
         SCIPinfoMessage(scip, file, "%7.2f%%", 100.0 * fastpercent);
      }
   }

   return SCIP_OKAY;
}


/*
 * display column specific interface methods
 */

/** creates the SDP-fastsettings display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispSdpunsolved(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata = NULL;

   assert( scip != NULL );

   /* create display column data */
   SCIP_CALL( SCIPallocMemory(scip, &dispdata) );

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO,
         dispCopySdpunsolved,
         dispFreeSdpunsolved, NULL, NULL,
         dispInitsolSdpunsolved, NULL, dispOutputSdpunsolved,
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   return SCIP_OKAY;
}
