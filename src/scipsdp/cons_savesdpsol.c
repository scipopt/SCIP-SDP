/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014      Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2014 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//#define SCIP_DEBUG
/**@file   cons_sdp.cpp
 * @brief  constraint handler for sdp-constraints
 * @author Tristan Gally
 */

//#define SCIP_DEBUG
//#define SCIP_MORE_DEBUG /* shows all cuts added */

#include "cons_savesdpsol.h"
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include <string.h>
#include <assert.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "Savesdpsol"
#define CONSHDLR_DESC          "saving the SDP solution at each node of the tree constraint handler"
#define CONSHDLR_SEPAPRIORITY         0 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_DELAYPRESOL      FALSE /**< should presolving method be delayed, if other presolvers found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/** constraint data to store optimal solution */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables and therefore length of sol */
   SCIP_Real*            sol;                /**< optimal solution for SDP-relaxation of this node */
};

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSavesdpsol)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( *consdata != NULL );

   SCIPdebugMessage("Deleting store node data constraint: <%s>.\n", SCIPconsGetName(cons));

   SCIPfreeBlockMemoryArrayNull(scip, &((*consdata)->sol), (*consdata)->nvars);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpSavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsSavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for primal solutions */
static
SCIP_DECL_CONSCHECK(consCheckSavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* do nothing */
   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockSavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* do nothing */
   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySavesdpsol)
{
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( valid != NULL );

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeConshdlrSavesdpsol(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopySavesdpsol)
{  /*lint --e{715}*/

   /* do not do anything: no Savesdpsol constraint should be present in the copy */
   return SCIP_OKAY;
}


/** include store bool constraint handler */
extern
SCIP_RETCODE SCIPincludeConshdlrSavesdpsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;

   /* include constraint handler */
   conshdlr = NULL;
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSavesdpsol, consEnfopsSavesdpsol, consCheckSavesdpsol, consLockSavesdpsol,
         NULL) );
   assert( conshdlr != NULL );

   /* set additional callbacks */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSavesdpsol) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySavesdpsol, consCopySavesdpsol) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSavesdpsol) );

   return SCIP_OKAY;
}



/*
 * External functions
 */

/** create a Savesdpsol-Cons, i.e. save the current optimal solution for the SDP-relaxation of this node */
SCIP_RETCODE createConsSavesdpsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables and therefore length of sol */
   SCIP_Real*            sol                 /**< optimal solution for SDP-relaxation of this node */
   )
{
   SCIP_CONSDATA* consdata = NULL;
   SCIP_CONSHDLR* conshdlr;
   int i;

   assert ( scip != NULL );
   assert ( name != NULL );
   assert ( nvars >= 0 );
   assert ( sol != NULL );

   /* find the node data constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "Savesdpsol");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Savesdpsol constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->sol), nvars) );
   consdata->nvars = nvars;
   for (i = 0; i < nvars; i++)
      consdata->sol[i] = sol[i];

   SCIPdebugMessage("Creating Savesdpsol constraint <%s>.\n", name);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, TRUE, FALSE, TRUE));

   return SCIP_OKAY;
}

/** for the given cons of type Savesdpsol returns the previous solution, length should start with the length of the array, this needs to be atleast
 *  the number of variables in scip and will be overwritten by this value, if it wasn't sufficient a debugMessage will be thrown */
SCIP_RETCODE getStartingPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to get starting point for */
   SCIP_Real*            sol,                /**< output: previous solution */
   int*                  length              /**< input: length of sol-array, output: number of entries in sol-array */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert ( scip != NULL );
   assert ( cons != NULL );
   assert ( sol != NULL );
   assert ( length != NULL );

   consdata = SCIPconsGetData(cons);

   assert ( consdata != NULL );

   if (*length < consdata->nvars)
   {
      SCIPdebugMessage("Not enough space for getStartingPoint in cons_savsdpsol, given %d, needed %d !\n", *length, consdata->nvars);
      *length = consdata->nvars;
      return SCIP_OKAY;
   }

   *length = consdata->nvars;
   for (i = 0; i < consdata->nvars; i++)
      sol[i] = consdata->sol[i];

   return SCIP_OKAY;
}
