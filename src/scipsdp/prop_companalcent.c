/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2024 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2024 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_companalcent.c
 * @brief  compute analytic center propagator
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_companalcent.h"
#include "relax_sdp.h"

/* fundamental propagator properties */
#define PROP_NAME              "companalcent"
#define PROP_DESC              "computes analytic center and forwards it to SDP relaxation handler"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     0 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_AFTERLPLOOP/**< propagation timing mask */

/*
 * Callback methods of propagator
 */

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecCompAnalCent)
{  /*lint --e{715}*/
   SCIP_RELAX* relax;

   relax = SCIPfindRelax(scip, "SDP"); /* get SDP relaxation handler */
   assert( relax != NULL );

   SCIP_CALL( SCIPrelaxSdpComputeAnalyticCenters(scip, relax) );

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the compute analytic center propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropCompAnalCent(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata = NULL;
   SCIP_PROP* prop = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecCompAnalCent, propdata) );

   assert( prop != NULL );

   return SCIP_OKAY;
}
