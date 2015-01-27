/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programms based on SCIP.                                     */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2015 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2015 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   nodesel_prio.c
 * @brief  nodeselector that chooses candidate with highest priority (dual bound is used as tiebreaker)
 * @author Tristan Gally
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* #define SCIP_DEBUG */

#include <assert.h>

#include "nodesel_prio.h"
#include "scip/type_misc.h" /* for SCIP Hashmap */


#define NODESEL_NAME            "prio"
#define NODESEL_DESC            "nodeselector that chooses candidate with highest priority (dual bound is used as tiebreaker)"
#define NODESEL_STDPRIORITY     10000000
#define NODESEL_MEMSAVEPRIORITY 10000000
#define HASHMAP_MAXLOADFACTOR   0.75 /**< if this factor is exceeded, the hashmap for the node priorities is enlarged */
#define HASHMAP_INITIALSIZEFACTOR 100 /**< the initial size of the hashmap is this times the number of binary and integer variables */
#define HASHMAP_ENLARGEFACTOR   2 /**< if the MAXLOADFACTOR of the hasmpa is exceeded, it is enlarged by this factor */


/*
 * Data structures
 */
struct SCIP_NodeselData
{
   SCIP_HASHMAP*         nodeprios;          /**< hashmap mapping B&B-Nodes to their priority, this returns a position in savedvals where the priority is saved */
   SCIP_Real*            savedvals;          /**< priority values saved in the hashmap */
   int                   nnodes;             /**< number of nodes in the hashmap */
   int                   size;               /**< size of the hashmap */
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** If the hashmap gets too crowded, create a new one of larger size to replace the old one and copy all entries to it */
static
SCIP_RETCODE CheckHashmapSize(
   SCIP*                 scip,               /**< scip-pointer */
   SCIP_NODESELDATA*     nodeseldata         /**< pointer to the data of the nodeselector */
   )
{
   assert( scip != NULL );
   assert( nodeseldata != NULL );

   if ( nodeseldata->nnodes >= HASHMAP_MAXLOADFACTOR * nodeseldata->size)
   {
      SCIP_HASHMAP* newhashmap;
      SCIP_HASHMAPLIST* currentlist;
      int i;

      SCIPdebugMessage("Resizing priority-nodeselector hashmap from size %d to size %d. \n", nodeseldata->size, (int) HASHMAP_ENLARGEFACTOR * nodeseldata->size);

      SCIPhashmapCreate(&newhashmap, SCIPblkmem(scip), (int) HASHMAP_ENLARGEFACTOR * nodeseldata->size);

      /* copy all entries from the old hashmap to the new one */
      for (i = 0; i < SCIPhashmapGetNLists(nodeseldata->nodeprios); i++)
      {
         currentlist = SCIPhashmapGetList(nodeseldata->nodeprios, i);

         while (currentlist != NULL)
         {
            SCIPhashmapInsert(newhashmap, SCIPhashmapListGetOrigin(currentlist), SCIPhashmapListGetImage(currentlist));

            currentlist = SCIPhashmapListGetNext(currentlist);
         }
      }

      SCIPhashmapFree(&(nodeseldata->nodeprios));
      nodeseldata->nodeprios = newhashmap;

      /* also enlarge the savedvals-array */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(nodeseldata->savedvals), nodeseldata->size, (int) HASHMAP_ENLARGEFACTOR * nodeseldata->size) );

      nodeseldata->size = (int) HASHMAP_ENLARGEFACTOR * nodeseldata->size;
   }

   return SCIP_OKAY;
}

/** inform the nodeselector about a new node being added by the branching rule and gives it the corresponding priority value */
SCIP_RETCODE SCIPnodeselPrioInsertNodePrio(
   SCIP*                  scip,              /**< scip pointer */
   SCIP_NODESEL*          nodesel,           /**< pointer to the nodeselector */
   SCIP_NODE*             node,              /**< node for which priority is given */
   SCIP_Real              priority           /**< nodeselection-priority of the given node */
   )
{
   SCIP_NODESELDATA* nodeseldata;

   assert( scip != NULL );
   assert( nodesel != NULL );
   assert( node != NULL );

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert( nodeseldata != NULL );

   /* check if the hashmap needs to be enlarge */
   SCIP_CALL( CheckHashmapSize(scip, nodeseldata) );

   nodeseldata->savedvals[nodeseldata->nnodes] = priority;
   SCIP_CALL( SCIPhashmapInsert(nodeseldata->nodeprios, (void*) node, (void*) (size_t) nodeseldata->nnodes) );
   nodeseldata->nnodes++;

   return SCIP_OKAY;
}


/*
 * Callback methods of node selector
 */

/* TODO: Implement all necessary node selector methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for node selector plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_NODESELCOPY(nodeselCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselCopyXyz NULL
#endif

/** destructor of node selector to free user data (called when SCIP is exiting) */
static
SCIP_DECL_NODESELFREE(nodeselFreePrio)
{
  SCIP_NODESELDATA* nodeseldata;

  nodeseldata = SCIPnodeselGetData(nodesel);
  assert(nodeseldata != NULL);

  if ( nodeseldata->size > 0 )
  {
     SCIPhashmapFree(&(nodeseldata->nodeprios));
     SCIPfreeBlockMemoryArray(scip, &(nodeseldata->savedvals), nodeseldata->size);
  }
  SCIPfreeMemory(scip, &nodeseldata);

  SCIPnodeselSetData(nodesel, NULL);

  return SCIP_OKAY;
}


/** initialization method of node selector (called after problem was transformed) */
static
SCIP_DECL_NODESELINIT(nodeselInitPrio)
{
   SCIP_NODESELDATA* nodeseldata;

   assert( nodesel != NULL );
   nodeseldata = SCIPnodeselGetData(nodesel);

   SCIPhashmapCreate(&(nodeseldata->nodeprios), SCIPblkmem(scip), HASHMAP_INITIALSIZEFACTOR * ( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) ) );
   nodeseldata->nnodes = 0;
   nodeseldata->size = HASHMAP_INITIALSIZEFACTOR * ( SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nodeseldata->savedvals), nodeseldata->size) );

   return SCIP_OKAY;
}


/** deinitialization method of node selector (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_NODESELEXIT(nodeselExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitXyz NULL
#endif


/** solving process initialization method of node selector (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_NODESELINITSOL(nodeselInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselInitsolXyz NULL
#endif


/** solving process deinitialization method of node selector (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_NODESELEXITSOL(nodeselExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz node selector not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nodeselExitsolXyz NULL
#endif


/** node selection method of node selector */
static
SCIP_DECL_NODESELSELECT(nodeselSelectPrio)
{
#ifdef SCIP_DEBUG
   SCIP_NODESELDATA* nodeseldata;
#endif
   assert( scip != NULL );
   assert( selnode != NULL );

   *selnode = SCIPgetBestNode(scip);

#ifdef SCIP_DEBUG
   nodeseldata = SCIPnodeselGetData(nodesel);
   assert( nodeseldata != NULL );
   SCIPdebugMessage("priority nodeselector chose node with priority %f.\n",
         nodeseldata->savedvals[(int) (size_t) SCIPhashmapGetImage(nodeseldata->nodeprios, (void*) *selnode)]);
#endif
   return SCIP_OKAY;
}

/** node comparison method of node selector */
static
SCIP_DECL_NODESELCOMP(nodeselCompPrio)
{
   SCIP_NODESELDATA* nodeseldata;

   assert( nodesel != NULL );
   assert( node1 != NULL );
   assert( node2 != NULL );

   nodeseldata = SCIPnodeselGetData(nodesel);
   assert( nodeseldata != NULL );

   /* if node 1 has a higher priority, we return a negative value, if node 2's priority is higher, it will be positive, in case of equal priorities, we take
    * the one with better dual bound (which is called lower bound in SCIP) */

   if ( SCIPisGT(scip, nodeseldata->savedvals[(int) (size_t) SCIPhashmapGetImage(nodeseldata->nodeprios, (void*) node1)],
                       nodeseldata->savedvals[(int) (size_t) SCIPhashmapGetImage(nodeseldata->nodeprios, (void*) node2)]) )
      return -1;
   else if ( SCIPisLT(scip, nodeseldata->savedvals[(int) (size_t) SCIPhashmapGetImage(nodeseldata->nodeprios, (void*) node1)],
                            nodeseldata->savedvals[(int) (size_t) SCIPhashmapGetImage(nodeseldata->nodeprios, (void*) node2)]) )
      return +1;
   else if (SCIPisLT(scip, SCIPnodeGetLowerbound(node1), SCIPnodeGetLowerbound(node2)))
      return -1;
   else if (SCIPisGT(scip, SCIPnodeGetLowerbound(node1), SCIPnodeGetLowerbound(node2)))
      return +1;
   else
      return 0;
   //TODO: check if Lowerbound should be smaller as well if maximizing instead of minimizing, but bfs-nodeselector doesn't check for objsense as well
}


/*
 * node selector specific interface methods
 */

/** creates the priority node selector and includes it in SCIP */
SCIP_RETCODE SCIPincludeNodeselPrio(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NODESELDATA* nodeseldata;
   SCIP_NODESEL* nodesel;

   /* create prio node selector data */
   nodeseldata = NULL;
   /* create node selector specific data here */
   SCIP_CALL( SCIPallocMemory(scip, &nodeseldata) );

   nodesel = NULL;

   /* use SCIPincludeNodeselBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeNodeselBasic(scip, &nodesel, NODESEL_NAME, NODESEL_DESC, NODESEL_STDPRIORITY, NODESEL_MEMSAVEPRIORITY,
          nodeselSelectPrio, nodeselCompPrio, nodeseldata) );

   assert(nodesel != NULL);

   /* set non fundamental callbacks via setter functions */
   /*SCIP_CALL( SCIPsetNodeselCopy(scip, nodesel, nodeselCopyXyz) );*/
   SCIP_CALL( SCIPsetNodeselFree(scip, nodesel, nodeselFreePrio) );
   SCIP_CALL( SCIPsetNodeselInit(scip, nodesel, nodeselInitPrio) );
/*   SCIP_CALL( SCIPsetNodeselExit(scip, nodesel, nodeselExitXyz) );
   SCIP_CALL( SCIPsetNodeselInitsol(scip, nodesel, nodeselInitsolXyz) );
   SCIP_CALL( SCIPsetNodeselExitsol(scip, nodesel, nodeselExitsolXyz) );*/

   return SCIP_OKAY;
}
