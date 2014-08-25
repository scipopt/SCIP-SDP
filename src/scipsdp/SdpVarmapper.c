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

/**@file   SdpVarmapper.h
 * @brief  class that maps scip variables to sdp indices (the SCIP variables are given sdp indices in the order in which they were inserted)
 * @author Tristan Gally
 */

#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__

#include "scip/scip.h"
#include "scip/type_misc.h"
#include "SdpVarmapper.h"

struct Sdpvarmapper
{
   SCIP_VAR**            sdptoscip;          /**< array of SCIP variables indexed by their SDP indices */
   SCIP_HASHMAP*         sciptosdp;          /**< hashmap that maps SCIP variables to their SDP indices */
   int                   nvars;              /**< number of variables saved in this varmapper */
};

/** creates the SDP Varmapper */
SCIP_RETCODE SdpVarmapperCreate(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper          /**< Pointer to the Varmapper that should be created */
      )
{
   assert ( scip != NULL );
   assert ( varmapper != NULL );

   SCIP_CALL(SCIPallocBlockMemory(scip, &varmapper));
   SCIP_CALL(SCIPhashmapCreate(&varmapper->sciptosdp, SCIPblkmem(scip), 0));
   varmapper->nvars = 0;

   return SCIP_OKAY;
}

/** frees the SDP Varmapper */
SCIP_RETCODE SdpVarmapperFree(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper          /**< Pointer to the Varmapper that should be freed */
      )
{
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );

   /* release all vars */
   for (i = 0; i < varmapper->nvars; i++)
   {
      SCIP_CALL(SCIPreleaseVar(scip, &vars[i]));
   }

   SCIP_CALL(SCIPhashmapFree(&varmapper->sciptosdp));
   SCIPfreeBlockMemoryArray(scip, &varmapper->sdptoscip, varmapper->nvars);
   SCIPfreeBlockMemory(scip, &varmapper);

   return SCIP_OKAY;
}

/** adds the given variables (if not already existant) to the end of the Varmapper */
SCIP_RETCODE SdpVarmapperAddVars(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper,         /**< Varmapper to add variables to */
   int                   nvars,             /**< number of variables to add to the varmapper */
   SCIP_VAR**            vars               /**< SCIP variables to add to the varmapper */
   )
{
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( nvars >= 0 );
   assert ( vars != NULL );

   SCIP_CALL(SCIPreallocBlockMemoryArray(scip, varmapper->sciptosdp, varmapper->nvars, varmapper->nvars + nvars));

   for (i = 0; i < nvars; i++)
   {
      if (!(SCIPhashmapExists(varmapper->sciptosdp, vars[i]))) /* make sure, that there are no duplicates in the lists */
      {
         sdptoscip[varmapper->nvars] = vars[i];
         SCIP_CALL(SCIPhashmapInsert(varmapper->sciptosdp, vars[i], (void*) varmapper->nvars));
         varmapper->nvars++;
         SCIP_CALL(SCIPcaptureVar(scip, vars[i]));
      }
      else
         SCIPdebugMessage("variable %s was not added to the varmapper as it was allready part of it \n", SCIPvarGetName(vars[i]));
   }

   return SCIP_OKAY;
}

/** adds the given variable (if not already existant) to the Varmapper at the given position */
SCIP_RETCODE SdpVarmapperInsertVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SdpVarmapper*         varmapper,          /**< Varmapper to add variables to */
   SCIP_VAR*             var,                /**< SCIP variable to add to the varmapper */
   int                   pos                 /**< position where the variable should be added */
   )
{
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( var != NULL );
   assert ( pos >= 0 );
   assert ( pos <= varmapper->nvars );

   if (!(SCIPhashmapExists(varmapper->sciptosdp, var))) /* make sure, that there are no duplicates in the lists */
   {
      if ( pos == varmapper->nvars )   /* add it to the end */
      {
         SCIP_CALL(SdpVarmapperAddVars(scip, varmapper, 1, &var));
      }
      else
      {
         SCIP_CALL(SCIPreallocBlockMemoryArray(scip, varmapper->sciptosdp, varmapper->nvars, varmapper->nvars + 1));

         /* move all variables after pos one spot to the right to make room for the new one */
         for (i = nvars - 1; i >= pos; i--)
         {
            varmapper->sdptoscip[i + 1] = varmapper->sdptoscip[i];
            SCIP_CALL(SCIPhashmapSetImage(varmapper->sciptosdp, varmapper->sdptoscip[i + 1], (void*) i + 1));
         }

         varmapper->sdptoscip[pos] = var;
         SCIP_CALL(SCIPhashmapInsert(varmapper->sciptosdp, var, (void*) pos));
      }
   }
   else
      SCIPdebugMessage("variable %s was not added to the varmapper as it was allready part of it \n", SCIPvarGetName(var));

   return SCIP_OKAY;
}

/** gets the number of variables */
int SdpVarmapperGetNVars(
   SdpVarmapper*         varmapper          /**< Varmapper to get number of variables for */
   )
{
   assert ( varmapper != NULL );

   return varmapper->nvars;
}

/** is this SCIP variable included in the varmapper? */
SCIP_Bool SdpVarmapperExistsSCIPvar(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   SCIP_VAR*             var                /**< SCIP variables to get sdp index for */
   )
{
   assert ( varmapper != NULL );
   assert ( var != NULL );

   return SCIPhashmapExists(varmapper->sciptosdp, var);
}

/** gets the sdp index for the given SCIP variable */
int SdpVarmapperGetSdpIndex(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   SCIP_VAR*             var                /**< SCIP variables to get sdp index for */
   )
{
   assert ( varmapper != NULL );
   assert ( var != NULL );

   return (int) SCIPhashmapGetImage(varmapper->sciptosdp, var);
}

/** gets the corresponding SCIP variable for the given sdp variable index */
SCIP_VAR* SdpVarmapperGetSCIPvar(
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   int                   ind,               /**< index of the sdp variable */
   )
{
   assert ( varmapper != NULL );
   assert ( ind >= 0 );
   assert ( ind < varmapper->nvars );

   return varmapper->sdptoscip[i];
}

/** removes the variable for the given Sdp index from the varmapper, decreasing the indices of all later variables by 1 */
EXTERN
SCIP_RETCODE SdpVarmapperRemoveSdpIndex(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper,         /**< Varmapper to get variable index for */
   int                   ind                /**< index of the sdp variable */
   )
{
   SCIP_VAR* var;
   int i;

   assert ( scip != NULL );
   assert ( varmapper != NULL );
   assert ( ind >= 0 );
   assert ( ind < varmapper->nvars );

   var = varmapper->sdptoscip[ind];

   assert ( SCIPhashmapExists(varmapper->sciptosdp, var) );

   SCIP_CALL( SCIPhashmapRemove(varmapper->sciptosdp, var) );

   /* shift all entries of the sdptoscip-array behind ind one to the left and update their sciptosdp-entries */
   for (i = ind + 1; i < varmapper->nvars; i++)
   {
      varmapper->sdptoscip[i - 1] = varmapper->sdptoscip[i];
      SCIP_CALL( SCIPhashmapSetImage (varmapper->sciptosdp, varmapper->sdptoscip[i - 1], i - 1) );
   }

   /* reallocate memory */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, varmapper->sdptoscip, varmapper->nvars, varmapper->nvars - 1) );

   varmapper->nvars--;

   return SCIP_OKAY;
}

/** swaps all SCIP variables for their transformed counterparts */
SCIP_RETCODE SdpVarmapperTransform(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         varmapper          /**< Pointer to the Varmapper that should be transformed */
   )
{
   int k;
   SCIP_VAR* var;

   assert ( scip != NULL );
   assert ( varmapper != NULL );

   for (k = 0; k < varmapper->nvars; ++k)
   {
      SCIP_CALL(SCIPgetTransformedVar(scip, sdptoscip[k], &var));
      SCIP_CALL(SCIPcaptureVar(scip, var));

      SCIP_CALL(SCIPhashmapRemove(varmapper->sciptosdp, sdptoscip[k]));
      SCIP_CALL(SCIPhashmapInsert(varmapper->sciptosdp, var, (void*) k));

      SCIP_CALL(SCIPreleaseVar(scip, &sdptoscip[k]));

      sdptoscip[k] = var;
   }

   return SCIP_OKAY;
}

/** clones the varmapper in the second argument to the varmapper in the third argument */
EXTERN
SCIP_RETCODE SdpVarmapperClone(
   SCIP*                 scip,              /**< SCIP data structure */
   SdpVarmapper*         oldmapper,         /**< Pointer to the Varmapper that should be cloned */
   SdpVarmapper*         newmapper          /**< Pointer to the Varmapper that should become a clone of the other one */
   )
{
   int nvars;
   int i;

   nvars = oldmapper->nvars;

   newmapper->nvars = nvars;

   /* allocate memory */
   SCIPallocBlockMemoryArray(scip, &newmapper->sdptoscip, nvars);
   SCIPallocBlockMemoryArray(scip, &newmapper->sciptosdp, nvars);

   /* copy entries */
   for (i = 0; i < nvars; i++)
   {
      newmapper->sdptoscip[i] = oldmapper->sdptoscip[i];
      SCIP_CALL(SCIPhashmapInsert(newmapper->sciptosdp, oldmapper->sdptoscip[i], (void*) i));
   }
}
