/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compute_symmetry_none.cpp
 * @brief  interface for no symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_sdpsymmetry.h"

/** return whether symmetry can be computed */
SCIP_Bool SDPSYMcanComputeSymmetry(void)
{
   return FALSE;
}

/** return name of external program used to compute generators */
const char* SDPSYMsymmetryGetName(void)
{
   return "none";
}

/** return description of external program used to compute generators */
const char* SDPSYMsymmetryGetDesc(void)
{
   return "";
}

/** compute generators of symmetry group */ /*lint -e{715}*/
SCIP_RETCODE SDPSYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SDPSYM_MATRIXDATA*    matrixdata,         /**< data for MIP matrix */
   SDPSYM_EXPRDATA*      exprdata,           /**< data for nonlinear constraints */
   SDPSYM_SDPDATA*       sdpdata,            /**< data for SDP constraints */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store size of group */
   SCIP_Bool             fixsdpvars          /**< whether variables in SDP constraints shall be fixed */
   )
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( matrixdata != NULL );
   assert( exprdata != NULL );
   assert( sdpdata != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;

   return SCIP_OKAY;
}
