/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
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

/**@file   table_slater.c
 * @brief  Slater statistics table
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "string.h"                     /* for strcmp */

#include "table_slater.h"
#include "relax_sdp.h"
#include "sdpi/sdpi.h"


#define TABLE_NAME              "slater"
#define TABLE_DESC              "Slater statistics table (needs relaxing/SDP/slatercheck > 0)"
#define TABLE_POSITION          16200              /**< the position of the statistics table */
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
SCIP_DECL_TABLECOPY(tableCopySlater)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( table != NULL );

   SCIP_CALL( SCIPincludeTableSlater(scip) );

   return SCIP_OKAY;
}


/** destructor of statistics table to free user data (called when SCIP is exiting) */
static
SCIP_DECL_TABLEFREE(tableFreeSlater)
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
SCIP_DECL_TABLEINITSOL(tableInitsolSlater)
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
SCIP_DECL_TABLEOUTPUT(tableOutputSlater)
{  /*lint --e{715}*/
   SCIP_TABLEDATA* tabledata;
   SCIP_RELAX* relaxsdp;
   int relaxslatercheck;

   assert( scip != NULL );
   assert( table != NULL );

   tabledata = SCIPtableGetData(table);
   assert( tabledata != NULL );

   relaxsdp = tabledata->relaxSDP;
   assert( relaxsdp != NULL );

   /* check if Slater statistics were stored in relaxator */
   SCIP_CALL( SCIPgetIntParam(scip, "relaxing/SDP/slatercheck", &relaxslatercheck) );

   if ( relaxslatercheck == 0 )
   {
      SCIPinfoMessage(scip, file, "    Slater:   no information available when relaxing/SDP/slatercheck = 0\n");
      return SCIP_OKAY;
   }
   else if ( relaxslatercheck < 0 || relaxslatercheck > 2 )
   {
      SCIPerrorMessage("Unknown parameter value %d for parameter relaxing/SDP/slatercheck in table/slater\n", relaxslatercheck);
      return SCIP_PARAMETERWRONGVAL;
   }

   /* Slater statistics */
   SCIPinfoMessage(scip, file, "    Slater         :      Holds      Fails Infeasible    Unknown\n");
   if ( tabledata->absolute )
   {
      SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
            "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
            "Dual Slater",
            SCIPrelaxSdpGetNdualSlaterHolds(relaxsdp), SCIPrelaxSdpGetNdualSlaterFails(relaxsdp),
            SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp), SCIPrelaxSdpGetNdualSlaterUnknown(relaxsdp));

      SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
            "         - %10" SCIP_LONGINT_FORMAT "\n",
            "Primal Slater",
            SCIPrelaxSdpGetNprimalSlaterHolds(relaxsdp), SCIPrelaxSdpGetNprimalSlaterFails(relaxsdp), SCIPrelaxSdpGetNprimalSlaterUnknown(relaxsdp));
   }
   else
   {
      SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
            "Dual Slater",
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNdualSlaterHolds(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNdualSlaterFails(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNdualSlaterUnknown(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp));

      SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %%        -   %8.2f %%\n",
            "Primal Slater",
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNprimalSlaterHolds(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNprimalSlaterFails(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp),
            100.0 * (SCIP_Real) SCIPrelaxSdpGetNprimalSlaterUnknown(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSdpInterfaceCalls(relaxsdp));
   }


   /* Slater solved statistics */

   if ( strcmp(SCIPsdpiGetSolverName(), "SDPA") == 0 )
   {
      SCIPinfoMessage(scip, file, "    Slater Solves  :       Fast     Stable    Penalty    Bounded    Unsolved\n");
      if ( tabledata->absolute )
      {
         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Slater holds",
               SCIPrelaxSdpGetNSlaterHoldsFast(relaxsdp), SCIPrelaxSdpGetNSlaterHoldsStable(relaxsdp), SCIPrelaxSdpGetNSlaterHoldsPenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterHoldsBounded(relaxsdp), SCIPrelaxSdpGetNSlaterHoldsUnsolved(relaxsdp));

         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Slater fails",
               SCIPrelaxSdpGetNSlaterFailsFast(relaxsdp), SCIPrelaxSdpGetNSlaterFailsStable(relaxsdp), SCIPrelaxSdpGetNSlaterFailsPenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterFailsBounded(relaxsdp), SCIPrelaxSdpGetNSlaterFailsUnsolved(relaxsdp));

         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Infeasible",
               SCIPrelaxSdpGetNSlaterInfeasibleFast(relaxsdp), SCIPrelaxSdpGetNSlaterInfeasibleStable(relaxsdp), SCIPrelaxSdpGetNSlaterInfeasiblePenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterInfeasibleBounded(relaxsdp), SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(relaxsdp));
      }
      else
      {
         if ( SCIPrelaxSdpGetNSlaterHolds(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Slater holds",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsStable(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp));
         }

         if ( SCIPrelaxSdpGetNSlaterFails(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Slater fails",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsStable(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp));
         }

         if ( SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Infeasible",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleStable(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasiblePenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp));
         }
      }
   }
   else
   {
      SCIPinfoMessage(scip, file, "    Slater         :       Fast    Penalty    Bounded    Unsolved\n");
      if ( tabledata->absolute )
      {
         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Slater holds",
               SCIPrelaxSdpGetNSlaterHoldsFast(relaxsdp), SCIPrelaxSdpGetNSlaterHoldsPenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterHoldsBounded(relaxsdp), SCIPrelaxSdpGetNSlaterHoldsUnsolved(relaxsdp));

         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Slater fails",
               SCIPrelaxSdpGetNSlaterFailsFast(relaxsdp), SCIPrelaxSdpGetNSlaterFailsPenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterFailsBounded(relaxsdp), SCIPrelaxSdpGetNSlaterFailsUnsolved(relaxsdp));

         SCIPinfoMessage(scip, file, "     %-14.14s: %10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT " "
               "%10" SCIP_LONGINT_FORMAT " %10" SCIP_LONGINT_FORMAT "\n",
               "Infeasible",
               SCIPrelaxSdpGetNSlaterInfeasibleFast(relaxsdp), SCIPrelaxSdpGetNSlaterInfeasiblePenalty(relaxsdp),
               SCIPrelaxSdpGetNSlaterInfeasibleBounded(relaxsdp), SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(relaxsdp));
      }
      else
      {
         if ( SCIPrelaxSdpGetNSlaterHolds(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Slater holds",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterHoldsUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterHolds(relaxsdp));
         }

         if ( SCIPrelaxSdpGetNSlaterFails(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Slater fails",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsPenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterFailsUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNSlaterFails(relaxsdp));
         }

         if ( SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp) > 0 )
         {
            SCIPinfoMessage(scip, file, "     %-14.14s: %8.2f %% %8.2f %% %8.2f %% %8.2f %%\n",
                  "Infeasible",
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleFast(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasiblePenalty(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleBounded(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp),
                  100.0 * (SCIP_Real) SCIPrelaxSdpGetNSlaterInfeasibleUnsolved(relaxsdp) / (SCIP_Real) SCIPrelaxSdpGetNdualSlaterInfeasible(relaxsdp));
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * statistics table specific interface methods
 */

/** creates the advanced SDP relaxator statistics table and includes it in SCIP */
SCIP_RETCODE SCIPincludeTableSlater(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_TABLEDATA* tabledata;

   assert( scip != NULL );

   /* create statistics table data */
   SCIP_CALL( SCIPallocMemory(scip, &tabledata) );

   /* include statistics table */
   SCIP_CALL( SCIPincludeTable(scip, TABLE_NAME, TABLE_DESC, TRUE,
         tableCopySlater, tableFreeSlater, NULL, NULL,
         tableInitsolSlater, NULL, tableOutputSlater,
         tabledata, TABLE_POSITION, TABLE_EARLIEST_STAGE) );

   /* add "absolute" parameter */
   SCIP_CALL( SCIPaddBoolParam( scip, "table/slater/absolute", "Should statistics be printed in absolute numbers (true) or percentages (false)?",
         &(tabledata->absolute), FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
