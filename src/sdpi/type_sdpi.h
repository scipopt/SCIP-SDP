/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_sdpi.h
 * @brief  type definitions for specific SDP solvers interface
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SDPI_H__
#define __SCIP_TYPE_SDPI_H__

#ifdef __cplusplus
extern "C" {
#endif

/* for now, we reuse the enums SCIP_OBJSEN, SCIP_PRICING, and SCIP_BASESTAT from the LPI */
#include "lpi/type_lpi.h"

/** SDP solver parameters */
enum SCIP_SDPParam
{
   SCIP_SDPPAR_FROMSCRATCH    =  0,      /**< solver should start from scratch at next call? */
   SCIP_SDPPAR_FASTMIP        =  1,      /**< fast mip setting of SDP solver */
   SCIP_SDPPAR_SCALING        =  2,      /**< should SDP solver use scaling? */
   SCIP_SDPPAR_PRESOLVING     =  3,      /**< should SDP solver use presolving? */
   SCIP_SDPPAR_PRICING        =  4,      /**< pricing strategy */
   SCIP_SDPPAR_SDPINFO        =  5,      /**< should SDP solver output information to the screen? */
   SCIP_SDPPAR_FEASTOL        =  6,      /**< feasibility tolerance for primal variables and slacks */
   SCIP_SDPPAR_DUALFEASTOL    =  7,      /**< feasibility tolerance for dual variables and reduced costs */
   SCIP_SDPPAR_BARRIERCONVTOL =  8,      /**< convergence tolerance used in barrier algorithm */
   SCIP_SDPPAR_LOBJLIM        =  9,      /**< lower objective limit */
   SCIP_SDPPAR_UOBJLIM        = 10,      /**< upper objective limit */
   SCIP_SDPPAR_SDPITLIM       = 11,      /**< SDP iteration limit */
   SCIP_SDPPAR_SDPTILIM       = 12,      /**< SDP time limit */
   SCIP_SDPPAR_MARKOWITZ      = 13,      /**< Markowitz tolerance */
   SCIP_SDPPAR_ROWREPSWITCH   = 14,      /**< simplex algorithm shall use row representation of the basis
                                          *   if number of rows divided by number of columns exceeds this value */
   SCIP_SDPPAR_THREADS        = 15,      /**< number of threads used to solve the SDP */
   SCIP_SDPPAR_CONDITIONLIMIT = 16       /**< maximum condition number of SDP basis counted as stable */
};
typedef enum SCIP_SDPParam SCIP_SDPPARAM;

/** SDP solution quality quantities */
enum SCIP_SDPSolQuality
{
   SCIP_SDPSOLQUALITY_ESTIMCONDITION,    /**< estimated condition number of (scaled) basis matrix (SCIP_Real) */
   SCIP_SDPSOLQUALITY_EXACTCONDITION     /**< exact condition number of (scaled) basis matrix (SCIP_Real) */
};
typedef enum SCIP_SDPSolQuality SCIP_SDPSOLQUALITY;

typedef struct SCIP_SDPi SCIP_SDPI;                 /**< solver dependent SDP interface */
typedef struct SCIP_SDPiState SCIP_SDPISTATE;       /**< complete SDP state (i.e. basis information) */
typedef struct SCIP_SDPiNorms SCIP_SDPINORMS;       /**< SDP pricing norms information */

#ifdef __cplusplus
}
#endif

#endif
