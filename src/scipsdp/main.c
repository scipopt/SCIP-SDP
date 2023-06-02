/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2023 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2023 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.c
 * @brief  main file for solving MISDPs
 * @author Sonja Mars
 * @author Tristan Gally
 */

#include "scipsdp/scipsdpdefplugins.h"
#include "scip/scipshell.h"

/** run scip and set some parameters */
static
SCIP_RETCODE runSCIP(
   int                   argc,               /**< number of command line arguments */
   char**                argv                /**< pointer to command line arguments */
   )
{
   SCIP* scip = NULL;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include plugins */
   SCIP_CALL( SCIPSDPincludeDefaultPlugins(scip) );

   /* change certain paramters: */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 5) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

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
