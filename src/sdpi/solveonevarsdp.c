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

/**@file   solveonevarsdp.c
 * @brief  Solve SDP with one variable
 * @author Marc Pfetsch
 */

#include "scip/pub_misc.h"
#include "sdpi/solveonevarsdp.h"
#include "sdpi/lapack_interface.h"
#include "sdpi/arpack_interface.h"

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

#ifdef ARPACK

/** determine whether linear combination with value alpha is feasible */
static
SCIP_RETCODE SCIPoneVarFeasibleSparse(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   blocksize,          /**< size of the SDP-block */
   SCIP_Real             alpha,              /**< variable value to test */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint-matrix */
   int*                  sdprow,             /**< array of row-indices of nonzero matrix entries */
   int*                  sdpcol,             /**< array of column-indices of nonzero matrix entries */
   SCIP_Real*            sdpval,             /**< array of nonzero values */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrix of the SDP-block */
   int*                  sdpconstrow,        /**< array of row-indices of constant matrix */
   int*                  sdpconstcol,        /**< array of column-indices of constant matrix */
   SCIP_Real*            sdpconstval,        /**< array of nonzero values of entries of constant matrix */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< corresponding eigenvector */
   )
{
   assert( sdpconstnnonz == 0 || sdpconstrow != NULL );
   assert( sdpconstnnonz == 0 || sdpconstcol != NULL );
   assert( sdpconstnnonz == 0 || sdpconstval != NULL );
   assert( sdpnnonz == 0 || sdprow != NULL );
   assert( sdpnnonz == 0 || sdpcol != NULL );
   assert( sdpnnonz == 0 || sdpval != NULL );
   assert( eigenvalue != NULL );

   SCIP_CALL( SCIParpackComputeSmallestEigenvectorOneVar(bufmem, blocksize, alpha, sdpnnonz, sdprow, sdpcol, sdpval, sdpconstnnonz, sdpconstrow, sdpconstcol, sdpconstval, eigenvalue, eigenvector) );

   return SCIP_OKAY;
}

#else

/** determine whether linear combination with value alpha is feasible */
static
SCIP_RETCODE SCIPoneVarFeasible(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   blocksize,          /**< size of the SDP-block */
   SCIP_Real*            tmpmatrix,          /**< temporary matrix */
   SCIP_Real*            fullconstmatrix,    /**< constant matrix */
   SCIP_Real*            fullmatrix,         /**< constant matrix */
   SCIP_Real             alpha,              /**< variable value to test */
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

#endif

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
   SCIP_Real* fullconstmatrix = NULL;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real* tmpmatrix = NULL;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real eigenvalue;
   SCIP_Real supergradient = SCIP_INVALID;
   SCIP_Real mu;
#ifndef ARPACK
   int r;
   int c;
   int i;
#endif

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

   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &eigenvector, blocksize) );

   /* fill in full matrices */
#ifndef ARPACK
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullconstmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &tmpmatrix, blocksize * blocksize) );

   for (i = 0; i < sdpconstnnonz; ++i)
   {
      r = sdpconstrow[i];
      c = sdpconstcol[i];
      assert( fullconstmatrix[r * blocksize + c] == 0.0 );
      assert( fullconstmatrix[c * blocksize + r] == 0.0 );
      fullconstmatrix[r * blocksize + c] = sdpconstval[i];
      fullconstmatrix[c * blocksize + r] = sdpconstval[i];
   }

   for (i = 0; i < sdpnnonz; ++i)
   {
      r = sdprow[i];
      c = sdpcol[i];
      assert( fullmatrix[r * blocksize + c] == 0.0 );
      assert( fullmatrix[c * blocksize + r] == 0.0 );
      fullmatrix[r * blocksize + c] = sdpval[i];
      fullmatrix[c * blocksize + r] = sdpval[i];
   }
#endif

   /* check upper bound */
#ifdef ARPACK
   SCIP_CALL( SCIPoneVarFeasibleSparse(bufmem, blocksize, ub, sdpnnonz, sdprow, sdpcol, sdpval, sdpconstnnonz, sdpconstrow, sdpconstcol, sdpconstval, &eigenvalue, eigenvector) );
#else
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, ub, &eigenvalue, eigenvector) );
#endif

   SCIPdebugMessage("ub = %g, minimal eigenvalue: %g\n", ub, eigenvalue);

   /* if matrix is not psd */
   if ( eigenvalue < -feastol )
   {
      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);
      SCIPdebugMessage(" -> supergradient: %g\n", supergradient);

      /* if supergradient is positive, then the problem is infeasible, because the minimal eigenvalue is increasing and we are not psd */
      if ( supergradient > 0.0 )
      {
         SCIPdebugMessage("Problem is infeasible (minimal eigenvalue: %g, supergradient = %g).\n", eigenvalue, supergradient);
         *objval = infinity;
         *optval = ub;

         goto TERMINATE;
      }
   }

   /* otherwise check lower bound */
#ifdef ARPACK
   SCIP_CALL( SCIPoneVarFeasibleSparse(bufmem, blocksize, lb, sdpnnonz, sdprow, sdpcol, sdpval, sdpconstnnonz, sdpconstrow, sdpconstcol, sdpconstval, &eigenvalue, eigenvector) );
#else
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, lb, &eigenvalue, eigenvector) );
#endif

   /* if matrix is psd, then the lower bound is optimal */
   if ( eigenvalue >= -feastol )
   {
      SCIPdebugMessage("Lower bound is optimal.\n");
      *objval = obj * lb;
      *optval = lb;

      goto TERMINATE;
   }
   else
   {
      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      /* if supergradient not positive, then the problem is infeasible because the eigenvalue is decreasing and we are not psd */
      if ( supergradient <= 0.0 )
      {
         SCIPdebugMessage("Problem is infeasible (minimal eigenvalue: %g, supergradient: %g).\n", eigenvalue, supergradient);
         *objval = infinity;
         *optval = lb;

         goto TERMINATE;
      }
   }
   assert( supergradient != SCIP_INVALID );
   SCIPdebugMessage("lb = %g, minimal eigenvalue: %g, supergradient: %g\n", lb, eigenvalue, supergradient);

   /* Invariant of the loop: the lower bound is not psd and the supergradient is positive */
   assert( eigenvector != NULL );
   assert( eigenvalue < -feastol );
   assert( supergradient > 0.0 );
   mu = lb;
   while ( eigenvalue < -feastol && supergradient > 0.0 )
   {
      assert( eigenvalue < -feastol );
      assert( supergradient > 0.0 );

      /* Compute estimate based on where the supergradient would reach -feastol/2.0 based on the supergradient inequality
      *  f(mu) \leq f(muold) + (mu - muold) g. We use feastol/2.0 to avoid little rounding errors. */
      mu = mu - (feastol / 2.0 + eigenvalue) / supergradient;

      /* Stop if we are larger than ub. In this case the problem is infeasible, since eigenvalue < -feastol by the while
       * condition. Infeasibility is detected below. */
      if ( mu > ub )
         break;

      /* compute eigenvalue and eigenvector */
#ifdef ARPACK
      SCIP_CALL( SCIPoneVarFeasibleSparse(bufmem, blocksize, mu, sdpnnonz, sdprow, sdpcol, sdpval, sdpconstnnonz, sdpconstrow, sdpconstcol, sdpconstval, &eigenvalue, eigenvector) );
#else
      SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, mu, &eigenvalue, eigenvector) );
#endif

      /* update supergradient */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      SCIPdebugMessage("mu = %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
   }

   /* check whether we are infeasible */
   if ( eigenvalue < -feastol )
   {
      SCIPdebugMessage("Problem infeasible; detected at %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
      *objval = infinity;
      *optval = mu;
   }
   else
   {
      SCIPdebugMessage("Solution is %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
      *objval = obj * mu;
      *optval = mu;
   }

 TERMINATE:
   BMSfreeBufferMemoryArrayNull(bufmem, &tmpmatrix);
   BMSfreeBufferMemoryArrayNull(bufmem, &fullmatrix);
   BMSfreeBufferMemoryArrayNull(bufmem, &fullconstmatrix);
   BMSfreeBufferMemoryArray(bufmem, &eigenvector);

   return SCIP_OKAY;
}


/** solves SDP with one variable and one SDP block - variant for dense constant matrix */
SCIP_RETCODE SCIPsolveOneVarSDPDense(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_Real             obj,                /**< objective coefficient of variable */
   SCIP_Real             lb,                 /**< lower bound of variable */
   SCIP_Real             ub,                 /**< upper bound of variable */
   int                   blocksize,          /**< size of the SDP-block */
   SCIP_Real*            fullconstmatrix,    /**< dense full constant matrix */
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
   SCIP_Real* fullmatrix;
   SCIP_Real* tmpmatrix;
   SCIP_Real* eigenvector = NULL;
   SCIP_Real eigenvalue;
   SCIP_Real supergradient = SCIP_INVALID;
   SCIP_Real mu;
   int r;
   int c;
   int i;

   assert( fullconstmatrix != NULL );
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
   BMS_CALL( BMSallocClearBufferMemoryArray(bufmem, &fullmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &tmpmatrix, blocksize * blocksize) );
   BMS_CALL( BMSallocBufferMemoryArray(bufmem, &eigenvector, blocksize) );

   for (i = 0; i < sdpnnonz; ++i)
   {
      r = sdprow[i];
      c = sdpcol[i];
      assert( fullmatrix[r * blocksize + c] == 0.0 );
      assert( fullmatrix[c * blocksize + r] == 0.0 );
      fullmatrix[r * blocksize + c] = sdpval[i];
      fullmatrix[c * blocksize + r] = sdpval[i];
   }

   /* check upper bound */
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, ub, &eigenvalue, eigenvector) );
   SCIPdebugMessage("ub = %g, minimal eigenvalue: %g\n", ub, eigenvalue);

   /* if matrix is not psd */
   if ( eigenvalue < -feastol )
   {
      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      /* if supergradient is positive, then the problem is infeasible, because the minimal eigenvalue is increasing and we are not psd */
      if ( supergradient > 0.0 )
      {
         SCIPdebugMessage("Problem is infeasible (minimal eigenvalue: %g, supergradient = %g).\n", eigenvalue, supergradient);
         *objval = infinity;
         *optval = ub;

         goto TERMINATE;
      }
   }

   /* otherwise check lower bound */
   SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, lb, &eigenvalue, eigenvector) );

   /* if matrix is psd, then the lower bound is optimal */
   if ( eigenvalue >= -feastol )
   {
      SCIPdebugMessage("Lower bound is optimal.\n");
      *objval = obj * lb;
      *optval = lb;

      goto TERMINATE;
   }
   else
   {
      /* compute supergradient value */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      /* if supergradient not positive, then the problem is infeasible because the eigenvalue is decreasing and we are not psd */
      if ( supergradient <= 0.0 )
      {
         SCIPdebugMessage("Problem is infeasible (minimal eigenvalue: %g, supergradient: %g).\n", eigenvalue, supergradient);
         *objval = infinity;
         *optval = lb;

         goto TERMINATE;
      }
   }
   assert( supergradient != SCIP_INVALID );
   SCIPdebugMessage("lb = %g, minimal eigenvalue: %g, supergradient: %g\n", lb, eigenvalue, supergradient);

   /* Invariant of the loop: the lower bound is not psd and the supergradient is positive */
   assert( eigenvector != NULL );
   assert( eigenvalue < -feastol );
   assert( supergradient > 0.0 );
   mu = lb;
   while ( eigenvalue < -feastol && supergradient > 0.0 )
   {
      assert( eigenvalue < -feastol );
      assert( supergradient > 0.0 );

      /* Compute estimate based on where the supergradient would reach -feastol/2.0 based on the supergradient inequality
      *  f(mu) \leq f(muold) + (mu - muold) g. We use feastol/2.0 to avoid little rounding errors. */
      mu = mu - (feastol / 2.0 + eigenvalue) / supergradient;

      /* Stop if we are larger than ub. In this case the problem is infeasible, since eigenvalue < -feastol by the while
       * condition. Infeasibility is detected below. */
      if ( mu > ub )
         break;

      /* compute eigenvalue and eigenvector */
      SCIP_CALL( SCIPoneVarFeasible(bufmem, blocksize, tmpmatrix, fullconstmatrix, fullmatrix, mu, &eigenvalue, eigenvector) );

      /* update supergradient */
      computeSupergradient(sdpnnonz, sdprow, sdpcol, sdpval, eigenvector, &supergradient);

      SCIPdebugMessage("mu = %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
   }

   /* check whether we are infeasible */
   if ( eigenvalue < -feastol )
   {
      SCIPdebugMessage("Problem infeasible; detected at %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
      *objval = infinity;
      *optval = mu;
   }
   else
   {
      SCIPdebugMessage("Solution is %.15g in [%.15g, %.15g], minimal eigenvalue: %g, supergradient: %g.\n", mu, lb, ub, eigenvalue, supergradient);
      *objval = obj * mu;
      *optval = mu;
   }

 TERMINATE:
   BMSfreeBufferMemoryArray(bufmem, &eigenvector);
   BMSfreeBufferMemoryArray(bufmem, &tmpmatrix);
   BMSfreeBufferMemoryArray(bufmem, &fullmatrix);

   return SCIP_OKAY;
}
