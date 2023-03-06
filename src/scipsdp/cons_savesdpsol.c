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

/**@file   cons_savesdpsol.c
 * @brief  constraint handler for saving SDP solutions in nodes
 * @author Tristan Gally
 * @author Marc Pfetsch
 */

/*#define SCIP_DEBUG*/
/*#define SCIP_MORE_DEBUG *//* shows all cuts added */

#include "cons_savesdpsol.h"
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include <string.h>
#include <assert.h>

/* constraint handler properties */
#define CONSHDLR_NAME          "Savesdpsol"
#define CONSHDLR_DESC          "saving the SDP solution at each node of the tree constraint handler"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/** constraint data to store solution */
struct SCIP_ConsData
{
   SCIP_Longint          node;               /**< index of the node the solution belongs to */
   SCIP_SOL*             sol;                /**< optimal solution for SDP-relaxation of this node; TODO: change to array*/
   SCIP_Real             maxprimalentry;     /**< maximal absolute value of primal matrix */
   int                   nlpcons;            /**< number of LP constraints of solution */
   int                   nblocks;            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz;   /**< starting point primal matrix X: number of nonzeros for each block (or NULL if nblocks == 0) */
   int**                 startXrow;          /**< starting point primal matrix X: row indices for each block (or NULL if nblocks = 0) */
   int**                 startXcol;          /**< starting point primal matrix X: column indices for each block (or NULL if nblocks = 0) */
   SCIP_Real**           startXval;          /**< starting point primal matrix X: values for each block (or NULL if nblocks = 0) */
};

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSavesdpsol)
{  /*lint --e{715}*/
   int b;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( *consdata != NULL );

   SCIPdebugMsg(scip, "Deleting store node data constraint: <%s>.\n", SCIPconsGetName(cons));

   assert( (*consdata)->startXval != NULL );
   assert( (*consdata)->startXcol != NULL );
   assert( (*consdata)->startXrow != NULL );
   assert( (*consdata)->startXnblocknonz != NULL );
   assert( (*consdata)->sol != NULL );

   for (b = 0; b < (*consdata)->nblocks; b++)
   {
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXval[b], (*consdata)->startXnblocknonz[b]);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXcol[b], (*consdata)->startXnblocknonz[b]);
      SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXrow[b], (*consdata)->startXnblocknonz[b]);
   }
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXval, (*consdata)->nblocks);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXcol, (*consdata)->nblocks);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXrow, (*consdata)->nblocks);
   SCIPfreeBlockMemoryArray(scip, &(*consdata)->startXnblocknonz, (*consdata)->nblocks);

   SCIP_CALL( SCIPfreeSol(scip, &(*consdata)->sol) );
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSavesdpsol)
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   /* do nothing */
   return SCIP_OKAY;
}


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySavesdpsol)
{  /*lint --e{715}*/
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
   /* do not copy: no Savesdpsol constraint should be present in the copy */
   return SCIP_OKAY;
}


/** include Savesdpsol constraint handler */
SCIP_RETCODE SCIPincludeConshdlrSavesdpsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;

   /* include constraint handler */
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

/** create a Savesdpsol constraint, i.e., save solution for the SDP-relaxation */
SCIP_RETCODE createConsSavesdpsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_Longint          node,               /**< index of the node the solution belongs to */
   int                   nlpcons,            /**< number of LP constraints of solution */
   SCIP_SOL*             sol,                /**< optimal solution for SDP-relaxation of this node */
   SCIP_Real             maxprimalentry,     /**< maximal absolute value of primal matrix */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz,   /**< starting point primal matrix X: number of nonzeros for each block (or NULL if nblocks == 0) */
   int**                 startXrow,          /**< starting point primal matrix X: row indices for each block (or NULL if nblocks = 0) */
   int**                 startXcol,          /**< starting point primal matrix X: column indices for each block (or NULL if nblocks = 0) */
   SCIP_Real**           startXval           /**< starting point primal matrix X: values for each block (or NULL if nblocks = 0) */
   )
{
   SCIP_CONSDATA* consdata = NULL;
   SCIP_CONSHDLR* conshdlr;
   int b;

   assert( scip != NULL );
   assert( name != NULL );
   assert( sol != NULL );
   assert( nblocks >= 0 );
   assert( nblocks == 0 || startXnblocknonz != NULL );
   assert( nblocks == 0 || startXrow != NULL );
   assert( nblocks == 0 || startXcol != NULL );
   assert( nblocks == 0 || startXval != NULL );

   SCIPdebugMsg(scip, "Creating Savesdpsol constraint <%s>.\n", name);

   /* find the node data constraint handler */
   conshdlr = SCIPfindConshdlr(scip, "Savesdpsol");
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Savesdpsol constraint handler not found.\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );

   consdata->node = node;
   consdata->nlpcons = nlpcons;

   SCIP_CALL( SCIPcreateSolCopy(scip, &consdata->sol, sol) );
   SCIP_CALL( SCIPunlinkSol(scip, consdata->sol) );
   consdata->maxprimalentry = maxprimalentry;

   if ( nblocks > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->startXnblocknonz, startXnblocknonz, nblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->startXrow, nblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->startXcol, nblocks) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->startXval, nblocks) );

      for (b = 0; b < nblocks; b++)
      {
         assert( startXnblocknonz != NULL );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->startXrow[b], startXrow[b], startXnblocknonz[b]) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->startXcol[b], startXcol[b], startXnblocknonz[b]) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &consdata->startXval[b], startXval[b], startXnblocknonz[b]) );
      }

      consdata->nblocks = nblocks;
   }
   else
   {
      consdata->startXnblocknonz = NULL;
      consdata->startXrow = NULL;
      consdata->startXcol = NULL;
      consdata->startXval = NULL;
      consdata->nblocks = 0;
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, FALSE,
         TRUE, FALSE, TRUE, FALSE, TRUE));

   return SCIP_OKAY;
}

/** for the given Savesdpsol constraint returns the node the information belongs to */
SCIP_Longint SCIPconsSavesdpsolGetNodeIndex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< Savesdpsol constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->node;
}

/** for the given Savesdpsol constraint returns the previous dual solution y */
SCIP_SOL* SCIPconsSavesdpsolGetDualSolution(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< Savesdpsol constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->sol != NULL );

   return consdata->sol;
}

/** for the given Savesdpsol constraint returns the maximal entry of primal solution X */
SCIP_Real SCIPconsSavesdpsolGetMaxPrimalEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< Savesdpsol constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->maxprimalentry;
}

/** for the given Savesdpsol constraint returns the number of LP constraints */
int SCIPconsSavesdpsolGetNLPcons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< Savesdpsol constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->nlpcons;
}

/** for the given cons of type Savesdpsol returns the number of nonzeros for each block of previous primal solution X */
SCIP_RETCODE SCIPconsSavesdpsolGetPrimalMatrixNonzeros(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< Savesdpsol constraint */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz    /**< input: allocated memory for startXrow/col/val; output: length of startXrow/col/val */
   )
{
   SCIP_CONSDATA* consdata;
   int b;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( nblocks != consdata->nblocks )
   {
      SCIPerrorMessage("SCIPconsSavesdpsolGetPrimalMatrix expected nblocks = %d but got %d\n", consdata->nblocks, nblocks);
      return SCIP_ERROR;
   }

   for (b = 0; b < nblocks; b++)
      startXnblocknonz[b] = consdata->startXnblocknonz[b];

   return SCIP_OKAY;
}

/** for the given cons of type Savesdpsol returns the previous primal solution X */
SCIP_RETCODE SCIPconsSavesdpsolGetPrimalMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< Savesdpsol constraint */
   int                   nblocks,            /**< number of blocks INCLUDING lp-block */
   int*                  startXnblocknonz,   /**< input: allocated memory for startXrow/col/val; output: length of startXrow/col/val */
   int**                 startXrow,          /**< pointer to store pointer to row indices of X */
   int**                 startXcol,          /**< pointer to store pointer to column indices of X */
   SCIP_Real**           startXval           /**< pointer to store pointer to values of X */
   )
{
   SCIP_CONSDATA* consdata;
   int b;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( nblocks != consdata->nblocks )
   {
      SCIPerrorMessage("SCIPconsSavesdpsolGetPrimalMatrix expected nblocks = %d, but got %d.\n", consdata->nblocks, nblocks);
      return SCIP_ERROR;
   }

   for (b = 0; b < nblocks; b++)
   {
      assert( startXnblocknonz[b] >= consdata->startXnblocknonz[b] );
      startXnblocknonz[b] = consdata->startXnblocknonz[b];
      for (i = 0; i < consdata->startXnblocknonz[b]; i++)
      {
         startXrow[b][i] = consdata->startXrow[b][i];
         startXcol[b][i] = consdata->startXcol[b][i];
         startXval[b][i] = consdata->startXval[b][i];
      }
   }

   return SCIP_OKAY;
}
