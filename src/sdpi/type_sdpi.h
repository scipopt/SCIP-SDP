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

/**@file   type_sdpi.h
 * @brief  type definitions for specific SDP solver interfaces
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
   SCIP_SDPPAR_EPSILON       = 0,      /**< convergence tolerance */
   SCIP_SDPPAR_FEASTOL       = 1,      /**< feasibility tolerance */
   SCIP_SDPPAR_OBJLIMIT      = 2,      /**< objective limit, if the SDP solver computes a lower bound for the minimzation
                                         *  problem that is bigger than this, it may stop */
#if 0
   SCIP_SDPPAR_THREADS       = 3,      /**< numer of threads */
#endif
   SCIP_SDPPAR_SDPINFO       = 4,      /**< should the SDP solver output information to the screen? */
   SCIP_SDPPAR_SLATERCHECK   = 5       /**< should the slater condition for the dual problem be checked before solving each SDP ? */
};
typedef enum SCIP_SDPParam SCIP_SDPPARAM;

typedef struct SCIP_SDPi SCIP_SDPI;                 /**< solver independent SDP interface */

#ifdef __cplusplus
}
#endif

#endif
