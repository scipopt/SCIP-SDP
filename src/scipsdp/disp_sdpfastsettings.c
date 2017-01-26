/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2017 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2017 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_sdpfastsettings.c
 * @brief  Column to display the percentage of SDP-relaxations that were solved using fast settings
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "disp_sdpfastsettings.h"
#include "relax_sdp.h"

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define DISP_NAME             "sdpfastsettings"
#define DISP_DESC             "percentage of fastsettings for SDP solver"
#define DISP_HEADER           "SDP fast"
#define DISP_WIDTH            8              /**< the width of the display column */
#define DISP_PRIORITY         1000           /**< the priority of the display column */
#define DISP_POSITION         1425           /**< the relative position of the display column */
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
SCIP_DECL_DISPCOPY(dispCopySdpfastsettings)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(disp != NULL);

   /* call inclusion method of dialog */
   SCIP_CALL( SCIPincludeDispSdpfastsettings(scip) );

   return SCIP_OKAY;
}

/** destructor of display column to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DISPFREE(dispFreeSdpfastsettings)
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
SCIP_DECL_DISPINITSOL(dispInitsolSdpfastsettings)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   assert( disp != NULL );
   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );

   dispdata->relaxSDP = SCIPfindRelax(scip, "SDP");
   assert( dispdata->relaxSDP != NULL );

   return SCIP_OKAY;
}

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputSdpfastsettings)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;


   assert( scip != NULL );
   assert( disp != NULL );

   dispdata = SCIPdispGetData(disp);

   assert( dispdata != NULL );
   assert( dispdata->relaxSDP != NULL );

   if ( SCIPrelaxSdpGetNSdpCalls(dispdata->relaxSDP) == 0 )
   {
      SCIPinfoMessage(scip, file, "   --   ");
   }
   else
   {
      SCIP_Real fastpercent;
      fastpercent = (SCIP_Real) SCIPrelaxSdpGetNSdpFast(dispdata->relaxSDP) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(dispdata->relaxSDP);
      SCIPinfoMessage(scip, file, "%7.2f%%", 100.0 * fastpercent);
   }

   return SCIP_OKAY;
}


/*
 * display column specific interface methods
 */

/** creates the SDP-fastsettings display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispSdpfastsettings(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata = NULL;

   assert( scip != NULL );

   /* create display column data */
   SCIP_CALL( SCIPallocMemory(scip, &dispdata) );

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO,
         dispCopySdpfastsettings,
         dispFreeSdpfastsettings, NULL, NULL,
         dispInitsolSdpfastsettings, NULL, dispOutputSdpfastsettings,
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   return SCIP_OKAY;
}
