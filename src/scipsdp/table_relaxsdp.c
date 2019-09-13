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

/**@file   table_relaxsdp.c
 * @brief  advanced SDP relaxator statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "table_relaxsdp.h"
#include "relax_sdp.h"


#define TABLE_NAME              "relaxsdp"
#define TABLE_DESC              "advanced SDP relaxator statistics table"
#define TABLE_POSITION          16000              /**< the position of the statistics table */
#define TABLE_EARLIEST_STAGE    SCIP_STAGE_SOLVING /**< output of the statistics table is only printed from this stage onwards */


/*
 * Data structures
 */

/** statistics table data */
struct SCIP_TableData
{
   SCIP_RELAX*           relaxSDP;           /**< pointer to the SDP relaxator whose iterations should be displayed */
};


/*
 * Callback methods of statistics table
 */

/** copy method for statistics table plugins (called when SCIP copies plugins) */
static
SCIP_DECL_TABLECOPY(tableCopyRelaxSdp)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( table != NULL );

   SCIP_CALL( SCIPincludeTableRelaxSdp(scip) );

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeRelaxSdp)
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
SCIP_DECL_TABLEINITSOL(tableInitsolRelaxSdp)
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
SCIP_DECL_TABLEOUTPUT(tableOutputRelaxSdp)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;
   SCIP_RELAX* relaxsdp;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   relaxsdp = tabledata->relaxSDP;
   assert( relaxsdp != NULL );

   SCIPinfoMessage(scip, file, "Relaxators         :       Time      Calls Iterations  Iter/call\n");

   if ( SCIPrelaxSdpGetNSdpCalls(relaxsdp) > 0 )
   {
      SCIPinfoMessage(scip, file, "  %-17.17s: %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10.2f \n",
         "SDP", SCIPrelaxGetTime(relaxsdp), SCIPrelaxGetNCalls(relaxsdp), SCIPrelaxSdpGetNIterations(relaxsdp),
         (SCIP_Real) SCIPrelaxSdpGetNIterations(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpCalls(relaxsdp));
   }
   else
   {
      SCIPinfoMessage(scip, file, "  %-17.17s: %10.2f %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10s \n",
         "SDP", SCIPrelaxGetTime(relaxsdp), SCIPrelaxGetNCalls(relaxsdp), SCIPrelaxSdpGetNIterations(relaxsdp), "-");
   }

   return SCIP_OKAY;
}


/*
 * statistics table specific interface methods
 */

/** creates the advanced SDP relaxator statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTableRelaxSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );

   /* create statistics table data */
   SCIP_CALL( SCIPallocMemory(scip, &tabledata) );

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         tableCopyRelaxSdp, tableFreeRelaxSdp, NULL, NULL,
         tableInitsolRelaxSdp, NULL, tableOutputRelaxSdp,
         tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   return SCIP_OKAY;
}
