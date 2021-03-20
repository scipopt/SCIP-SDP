/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2021 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2021 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solveonvarsdp.h
 * @brief  Solve SDP with one variable
 * @author Marc Pfetsch
 */

#include "sdpi/solveonevarsdp.h"
#include "sdpi/lapack_interface.h"

/** Checks if a BMSallocMemory-call was successfull, otherwise returns SCIP_NOMEMORY */
#define BMS_CALL(x)   do                                                                                      \
                      {                                                                                       \
                          if( NULL == (x) )                                                                   \
                          {                                                                                   \
                             SCIPerrorMessage("No memory in function call\n");                                \
                             return SCIP_NOMEMORY;                                                            \
                          }                                                                                   \
                      }                                                                                       \
                      while( FALSE )


/** determine whether linear combination with value alpha is feasible */
static
SCIP_RETCODE SCIPoneVarFeasible(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   blocksize,          /**< size of the SDP-block */
   SCIP_Real*            tmpmatrix,          /**< temporary matrix */
   SCIP_Real*            fullconstmatrix,    /**< constant matrix */
   SCIP_Real*            fullmatrix,         /**< constant matrix */
   SCIP_Real             alpha,              /**< variable value to test */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< corresponding eigenvector */
   )
{
   int i;

   assert( fullconstmatrix != NULL );
   assert( fullmatrix != NULL );
   assert( eigenvalue != NULL );

   /* compute linear combination */
   for (i = 0; i < blocksize * blocksize; ++i)
      tmpmatrix[i] = alpha * fullmatrix[i] - fullconstmatrix[i];

   if ( eigenvector != NULL )
   {
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, TRUE, blocksize, tmpmatrix, 1, eigenvalue, eigenvector) );
   }
   else
   {
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(bufmem, FALSE, blocksize, tmpmatrix, 1, eigenvalue, NULL) );
   }

   return SCIP_OKAY;
}


/** solves SDP with one variable and one SDP block */
SCIP_RETCODE SCIPsolveOneVarSDP(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real             obj,                /**< objective coefficient of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   int                   blocksize,          /**< size of the SDP-block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrix of the SDP-block */
   int*                  sdpconstrow,        /**< pointer to row-indices of constant matrix */
   int*                  sdpconstcol,        /**< pointer to column-indices of constant matrix */
   SCIP_Real*            sdpconstval,        /**< pointer to values of entries of constant matrix */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< pointer to the row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< pointer to the column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< pointer to the nonzero values */
   SCIP_Real             infinity,           /**< infinity value */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            objval,             /**< pointer to store optimal objective value */
   SCIP_Real*            optval              /**< pointer to store optimal value of variable */
   )
{
   SCIP_Real* fullconstmatrix;
   SCIP_Real* fullmatrix;
   SCIP_Real* tmpmatrix;
   SCIP_Real* eigenvector;
   SCIP_Real eigenvalue;
   SCIP_Real supergradient;
   int r;
   int c;
   int i;

   assert( sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz == 0 || sdprow != NULL );
   assert( sdpnnonz == 0 || sdpcol != NULL );
   assert( sdpnnonz == 0 || sdpval != NULL );
   assert( objval != NULL );
   assert( optval != NULL );

   *objval = SCIP_INVALID;
   *optval = SCIP_INVALID;

   /* can currently only treat the case with finite lower and upper bounds */
   if ( lb <= -infinity )
      return SCIP_OKAY;
   if ( ub >= infinity )
      return SCIP_OKAY;

   /* can currently only treat nonnegative objectie functions */
   if ( obj < 0.0 )
      return SCIP_OKAY;

   SCIPdebugMessage("Solve SDP with one variable (obj = %g, lb = %g, ub = %g).\n", obj, lb, ub);

   /* fill in full matrices */
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullconstmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &tmpmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &eigenvector, blocksize) );

   for (i = 0; i < sdpconstnnonz; ++i)
   {
      r = sdpconstrow[i];
      c = sdpconstcol[i];
      fullconstmatrix[r * blocksize + c] = sdpconstval[i];
      fullconstmatrix[c * blocksize + r] = sdpconstval[i];
   }

   for (i = 0; i < sdpnnonz; ++i)
   {
      r = sdprow[i];
      c = sdpcol[i];
      fullmatrix[r * blocksize + c] = sdpval[i];
      fullmatrix[c * blocksize + r] = sdpval[i];
   }

   /* check finite upper bound */
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, ub, feastol, &eigenvalue, eigenvector) );

   /* if combination is not psd */
   if ( eigenvalue < -feastol )
   {
      /* compute supergradient value */
      supergradient = 0.0;
      for (i = 0; i < sdpnnonz; ++i)
      {
         r = sdprow[i];
         c = sdpcol[i];

         supergradient += sdpval[i] * eigenvector[r] * eigenvector[c];
         if ( r != c )
            supergradient += sdpval[i] * eigenvector[c] * eigenvector[r];
      }

      /* if supergradient is positive combination, then the problem is infeasible, because the minimal eigenvalue is increasing and we are not psd */
      if ( supergradient > feastol )
      {
         SCIPdebugMessage("Problem is infeasible (eigenvalue: %g, supergradient = %g).\n", eigenvalue, supergradient);
         *objval = infinity;
         *optval = ub;

         goto TERMINATE;
      }
   }

   /* otherwise check lower bound */
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, lb, feastol, &eigenvalue, NULL) );

   /* if combination is psd, then the lower bound is optimal */
   if ( eigenvalue >= -feastol )
   {
      SCIPdebugMessage("Lower bound is optimal.\n");
      *objval = obj * lb;
      *optval = lb;

      goto TERMINATE;
   }

 TERMINATE:
   BMSfreeBufferMemoryArray(bufmem, &eigenvector);
   BMSfreeBufferMemoryArray(bufmem, &tmpmatrix);
   BMSfreeBufferMemoryArray(bufmem, &fullmatrix);
   BMSfreeBufferMemoryArray(bufmem, &fullconstmatrix);

   return SCIP_OKAY;
}
