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

/**@file   table_sdpsolversuccess.c
 * @brief  SDP solver success statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "string.h"                     /* for strcmp */

#include "table_sdpsolversuccess.h"
#include "relax_sdp.h"
#include "sdpi/sdpi.h"


#define TABLE_NAME              "sdpsolversuccess"
#define TABLE_DESC              "SDP solver success statistics table"
#define TABLE_POSITION          17100              /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_SOLVING /**< output of the statistics table is only printed from this stage onwards */


/*
 * Data structures
 */

/** statistics table data */
struct SCIP_TableData
{
   SCIP_RELAX*           relaxSDP;           /**< pointer to the SDP relaxator whose iterations should be displayed */
   SCIP_Bool             absolute;           /**< Should statistics be printed in absolute numbers (true) or percentages (false)? */
};


/*
 * Callback methods of statistics table
 */

/** copy method for statistics table plugins (called when SCIP copies plugins) */
static
SCIP_DECL_TABLECOPY(tableCopySdpSolverSuccess)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( table != NULL );

   SCIP_CALL( SCIPincludeTableSdpSolverSuccess(scip) );

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeSdpSolverSuccess)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );
   assert( table != NULL );
   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   SCIPfreeMemory(scip, &tabledata);
   SCIPtableSetData(table, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of statistics table (called when branch and bound process is about to begin) */
static
SCIP_DECL_TABLEINITSOL(tableInitsolSdpSolverSuccess)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;

   assert( table != NULL );
   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   tabledata->relaxSDP = SCIPfindRelax(scip, "SDP");
   assert( tabledata->relaxSDP != NULL );

   return SCIP_OKAY;
}


/** output method of statistics table to output file stream 'file' */
static
SCIP_DECL_TABLEOUTPUT(tableOutputSdpSolverSuccess)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;
   SCIP_RELAX* relaxsdp;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   relaxsdp = tabledata->relaxSDP;
   assert( relaxsdp != NULL );

   if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0 )
   {
      SCIPinfoMessage(scip, file, "    SDP-Solvers    :       Time    Opttime Fast     Medium     Stable    Penalty    Unsolved\n");
      if ( tabledata->absolute )
      {
         SCIPinfoMessage(scip, file, "     %-14.14s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
            "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
            SCIPsdpiGetSolverName(), SCIPrelaxSdpGetSolvingTime(scip, relaxsdp), SCIPrelaxSdpGetOptTime(relaxsdp),
            SCIPrelaxSdpGetNSdpFast(relaxsdp), SCIPrelaxSdpGetNSdpMedium(relaxsdp), SCIPrelaxSdpGetNSdpStable(relaxsdp), SCIPrelaxSdpGetNSdpPenalty(relaxsdp),
            SCIPrelaxSdpGetNSdpUnsolved(relaxsdp));
      }
      else
      {
         SCIPinfoMessage(scip, file, "     %-14.14s: %10.2f %10.2f %8.2f %% %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
            SCIPsdpiGetSolverName(), SCIPrelaxSdpGetSolvingTime(scip, relaxsdp), SCIPrelaxSdpGetOptTime(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpMedium(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpStable(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp));
      }
   }
   else
   {
      SCIPinfoMessage(scip, file, "    SDP-Solvers    :       Time    Opttime Default    Penalty   Unsolved\n");
      if ( tabledata->absolute )
      {
         SCIPinfoMessage(scip, file, "     %-14.14s: %10.2f %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
            "%10" SCIP_LONGINT_FORMAT "\n",
            SCIPsdpiGetSolverName(), SCIPrelaxSdpGetSolvingTime(scip, relaxsdp), SCIPrelaxSdpGetOptTime(relaxsdp),
            SCIPrelaxSdpGetNSdpFast(relaxsdp), SCIPrelaxSdpGetNSdpPenalty(relaxsdp), SCIPrelaxSdpGetNSdpUnsolved(relaxsdp));
      }
      else
      {
         if ( SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %10.2f %10.2f %8.2f %% %8.2f %% %8.2f %%\n",
               SCIPsdpiGetSolverName(), SCIPrelaxSdpGetSolvingTime(scip, relaxsdp), SCIPrelaxSdpGetOptTime(relaxsdp),
               100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
               100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
               100.0 * (SCIP_Real) SCIPrelaxSdpGetNSdpUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp));
         }
         else
         {
            SCIPinfoMessage(scip, file, "     %-14.14s:   %10s %10s %8s   %8s   %8s\n", SCIPsdpiGetSolverName(), "-", "-", "-", "-", "-");
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * statistics table specific interface methods
 */

/** creates the SDP solver success statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTableSdpSolverSuccess(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );

   /* create statistics table data */
   SCIP_CALL( SCIPallocMemory(scip, &tabledata) );

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         tableCopySdpSolverSuccess, tableFreeSdpSolverSuccess, NULL, NULL,
         tableInitsolSdpSolverSuccess, NULL, tableOutputSdpSolverSuccess,
         tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* add "absolute" parameter */
   SCIP_CALL( SCIPaddBoolParam( scip, "table/sdpsolversuccess/absolute", "Should statistics be printed in absolute numbers (true) or percentages (false)?",
         &(tabledata->absolute), FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
