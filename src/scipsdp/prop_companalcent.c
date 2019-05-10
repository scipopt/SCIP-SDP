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
 * Data structures
 */


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


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
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create xyz propagator data */
   propdata = NULL;
   prop = NULL;

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecCompAnalCent, propdata) );

   assert(prop != NULL);

   return SCIP_OKAY;
}
