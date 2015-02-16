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
 * @author Tristan Gally
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
   SCIP_SDPPAR_EPSILON       = 0,      /* convergence tolerance */
   SCIP_SDPPAR_FEASTOL       = 1,      /* feasibility tolerance */
   SCIP_SDPPAR_OBJLIMIT      = 2       /* objective limit, if the SDP solver computes a lower bound for the minimzation problem that is bigger than this, it may stop */
};
typedef enum SCIP_SDPParam SCIP_SDPPARAM;

typedef struct SCIP_SDPi SCIP_SDPI;                 /**< solver independent SDP interface */

#ifdef __cplusplus
}
#endif

#endif
