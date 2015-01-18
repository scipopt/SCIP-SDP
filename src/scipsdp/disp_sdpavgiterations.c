/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp_sdpavgiterations.c
 * @brief  Column to display the average number of SDP iterations
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "disp_sdpavgiterations.h"
#include "relax_sdp.h"


#define DISP_NAME             "sdpavgiterations"
#define DISP_DESC             "average number of SDP iterations"
#define DISP_HEADER           "SDP it/n"
#define DISP_WIDTH            8              /**< the width of the display column */
#define DISP_PRIORITY         25001          /**< the priority of the display column */
#define DISP_POSITION         1400           /**< the relative position of the display column */
#define DISP_STRIPLINE        TRUE           /**< default for displaying column separated with a line from its right neighbor */




/*
 * Data structures
 */

/** display column data */
struct SCIP_DispData
{
   SCIP_RELAX*           relaxSDP;           /**< pointer to the SDP relaxator whose iterations should be displayed */
};


/*
 * Callback methods of display column
 */

/** copy method for dialog plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_DISPCOPY(dispCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispCopyXyz NULL
#endif

/** destructor of display column to free user data (called when SCIP is exiting) */
static
SCIP_DECL_DISPFREE(dispFreeSdpavgiterations)
{
   SCIP_DISPDATA* dispdata;

   assert( disp != NULL );

   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );

   SCIPfreeMemory(scip, &dispdata);
   SCIPdispSetData(disp, NULL);

   return SCIP_OKAY;
}

/** initialization method of display column (called after problem was transformed) */
#if 0
static
SCIP_DECL_DISPINIT(dispInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispInitXyz NULL
#endif

/** deinitialization method of display column (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_DISPEXIT(dispExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitXyz NULL
#endif

/** solving process initialization method of display column (called when branch and bound process is about to begin) */
static
SCIP_DECL_DISPINITSOL(dispInitsolSdpavgiterations)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   assert ( disp != NULL );

   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );

   dispdata->relaxSDP = SCIPfindRelax(scip, "SDP");
   assert( dispdata->relaxSDP != NULL );

   return SCIP_OKAY;
}

/** solving process deinitialization method of display column (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_DISPEXITSOL(dispExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz display column not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define dispExitsolXyz NULL
#endif

/** output method of display column to output file stream 'file' */
static
SCIP_DECL_DISPOUTPUT(dispOutputSdpavgiterations)
{  /*lint --e{715}*/
   SCIP_DISPDATA* dispdata;

   assert( scip != NULL );
   assert( disp != NULL );

   dispdata = SCIPdispGetData(disp);
   assert( dispdata != NULL );
   assert( dispdata->relaxSDP != NULL );

   if ( SCIPrelaxSdpGetNSdpCalls(dispdata->relaxSDP) == 0 )
      SCIPinfoMessage(scip, file, "      - ");
   else
      SCIPinfoMessage(scip, file, "%7.1f ", (double) SCIPrelaxSdpGetNIterations(dispdata->relaxSDP) / (double) SCIPrelaxSdpGetNSdpCalls(dispdata->relaxSDP) );

   return SCIP_OKAY;
}



/*
 * display column specific interface methods
 */

/** creates the SDP-average-iterations display column and includes it in SCIP */
SCIP_RETCODE SCIPincludeDispSdpavgiterations(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DISPDATA* dispdata = NULL;

   assert ( scip != NULL );

   /* create display column data */
   SCIP_CALL( SCIPallocMemory(scip, &dispdata) );

   /* include display column */
   SCIP_CALL( SCIPincludeDisp(scip, DISP_NAME, DISP_DESC, DISP_HEADER, SCIP_DISPSTATUS_AUTO,
         dispCopyXyz,
         dispFreeSdpavgiterations, dispInitXyz, dispExitXyz,
         dispInitsolSdpavgiterations, dispExitsolXyz, dispOutputSdpavgiterations,
         dispdata, DISP_WIDTH, DISP_PRIORITY, DISP_POSITION, DISP_STRIPLINE) );

   return SCIP_OKAY;
}
