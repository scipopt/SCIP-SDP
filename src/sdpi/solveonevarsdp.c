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


/** compute supergradient */
static
void computeSupergradient(
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< array of row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< array of column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< array of nonzero values */
   SCIP_Real*            eigenvector,        /**< array of eigenvector */
   SCIP_Real*            supergradient       /**< pointer to store the value of the supergradient */
   )
{
   int r;
   int c;
   int i;

   assert( supergradient != NULL );

   *supergradient = 0.0;
   for (i = 0; i < sdpnnonz; ++i)
   {
      r = sdprow[i];
      c = sdpcol[i];

      *supergradient += sdpval[i] * eigenvector[r] * eigenvector[c];
      if ( r != c )
         *supergradient += sdpval[i] * eigenvector[c] * eigenvector[r];
   }
}

/** solves SDP with one variable and one SDP block */
SCIP_RETCODE SCIPsolveOneVarSDP(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real             obj,                /**< objective coefficient of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   int                   blocksize,          /**< size of the SDP-block */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrix of the SDP-block */
   int*                  sdpconstrow,        /**< array of row-indices of constant matrix */
   int*                  sdpconstcol,        /**< array of column-indices of constant matrix */
   SCIP_Real*            sdpconstval,        /**< array of nonzero values of entries of constant matrix */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< array of row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< array of column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< array of nonzero values */
   SCIP_Real             infinity,           /**< infinity value */
   SCIP_Real             feastol,            /**< feasibility tolerance */
   SCIP_Real*            objval,             /**< pointer to store optimal objective value */
   SCIP_Real*            optval              /**< pointer to store optimal value of variable */
   )
{
   SCIP_Real* fullconstmatrix;
   SCIP_Real* fullmatrix;
   SCIP_Real* tmpmatrix;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real* eigenvectorub;
   SCIP_Real* eigenvectorlb;
   SCIP_Real eigenvalue;
   SCIP_Real eigenvaluelb;
   SCIP_Real eigenvalueub;
   SCIP_Real supergradient = SCIP_INVALID;
   SCIP_Real mu;
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

   /* can currently only treat nonnegative objective functions */
   if ( obj < 0.0 )
      return SCIP_OKAY;

   SCIPdebugMessage("Solve SDP with one variable (obj = %g, lb = %g, ub = %g).\n", obj, lb, ub);

   /* fill in full matrices */
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullconstmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &tmpmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &eigenvectorlb, blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &eigenvectorub, blocksize) );

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
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, ub, feastol, &eigenvalueub, eigenvectorub) );
   SCIPdebugMessage("ub = %g, eigenvalue: %g\n", ub, eigenvalueub);

   /* if matrix is not psd */
   if ( eigenvalueub < -feastol )
   {
      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvectorub, &supergradient);

      /* if supergradient is positive, then the problem is infeasible, because the minimal eigenvalue is increasing and we are not psd */
      if ( supergradient > feastol )
      {
         SCIPdebugMessage("Problem is infeasible (eigenvalue: %g, supergradient = %g).\n", eigenvalueub, supergradient);
         *objval = infinity;
         *optval = ub;

         goto TERMINATE;
      }
   }

   /* otherwise check lower bound */
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, lb, feastol, &eigenvaluelb, eigenvectorlb) );
   SCIPdebugMessage("lb = %g, eigenvalue: %g\n", lb, eigenvaluelb);

   /* if matrix is psd, then the lower bound is optimal */
   if ( eigenvaluelb >= -feastol )
   {
      SCIPdebugMessage("Lower bound is optimal.\n");
      *objval = obj * lb;
      *optval = lb;

      goto TERMINATE;
   }

   /* choose starting point depending on which eigenvalue is nearer to 0 */
   if ( REALABS(eigenvaluelb) < REALABS(eigenvalueub) )
   {
      mu = lb;
      eigenvector = eigenvectorlb;
      eigenvalue = eigenvaluelb;

      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);
   }
   else
   {
      mu = ub;
      eigenvector = eigenvectorub;
      eigenvalue = eigenvalueub;

      if ( eigenvalueub < -feastol )
         assert( supergradient != SCIP_INVALID ); /* supergradient has been computed already */
      else
      {
         /* compute supergradient value */
         computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);
      }
   }
   assert( eigenvector != NULL );
   assert( supergradient != SCIP_INVALID );

   /* we now in the case in which the lower bound is not psd, but the upper bound is */
   while ( ub - lb > feastol )
   {
      /* is supergradient is large enough */
      if ( REALABS(supergradient) > feastol )
      {
         /* compute estimate based on where the supergradient would reach 0 */
         mu = mu - eigenvalue / supergradient;

         /* if estimate is smaller than the lower bound, use bisection */
         if ( mu <= lb )
            mu = (lb + ub) / 2.0;
#if 0
         /* if the estimate is too close to lb or ub, use bisection */
         else if ( mu <= lb + 0.05 * (ub - lb) ||  mu >= ub - 0.05 * (ub - lb) )
            mu = (lb + ub) / 2.0;
#endif
      }
      else
         mu = (lb + ub) / 2.0;

      /* compute eigenvalue and eigenvector */
      SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, mu, feastol, &eigenvalue, eigenvector) );

      /* update supergradient */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      SCIPdebugMessage("mu = %.15g in [%.15g, %.15g] (delta: %g), eigenvalue: %g, supergradient: %g\n", mu, lb, ub, ub - lb, eigenvalue, supergradient);

      /* check early termination */
      if ( REALABS(eigenvalue) <= feastol / 10.0 && REALABS(supergradient) > feastol )
         break;

      if ( eigenvalue >= 0.0 )
         ub = mu;
      else
         lb = mu;
   }
   SCIPdebugMessage("Solution is %.15g in [%.15g, %.15g] (delta: %g), eigenvalue: %g, supergradient: %g\n", mu, lb, ub, ub - lb, eigenvalue, supergradient);

   *objval = obj * mu;
   *optval = mu;

 TERMINATE:
   BMSfreeBufferMemoryArray(bufmem, &eigenvectorub);
   BMSfreeBufferMemoryArray(bufmem, &eigenvectorlb);
   BMSfreeBufferMemoryArray(bufmem, &tmpmatrix);
   BMSfreeBufferMemoryArray(bufmem, &fullmatrix);
   BMSfreeBufferMemoryArray(bufmem, &fullconstmatrix);

   return SCIP_OKAY;
}
