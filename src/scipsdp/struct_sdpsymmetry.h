/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2022 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2022 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_sdpsymmetry.h
 * @brief  structs for symmetry detection in SDP
 * @author Christopher Hojny
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SDPSYMMETRY_H__
#define __SCIP_TYPE_SDPSYMMETRY_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data for symmetry group computation on SDP onstraints */
struct SYM_Sdpdata
{
   int                   nsdpconss;          /**< number of SDP constraints */
   int*                  blocksizes;         /**< for each SDP cons, the size of the SDP-block */
   int*                  nvars;              /**< for each SDP cons, the number of variables in this cons */
   SCIP_Real**           vals;               /**< for each SDP cons, the values of nonzeros for each variable and
                                              *   constant block */
   SCIP_VAR***           vars;               /**< for each SDP cons, the variables present in the cons */
   int**                 valsbegins;         /**< for each SDP cons, the begin positions of new matrix in vals */
   int**                 cols;               /**< for each SDP cons, the column index of value in vals */
   int**                 rows;               /**< for each SDP cons, the row index of value in vals */
   int**                 colors;             /**< for each SDP cons, the color of the coefficients */
   int**                 colors2;            /**< for each SDP cons, the second color of the coefficients */
};
typedef struct SYM_Sdpdata SYM_SDPDATA;

#ifdef __cplusplus
}
#endif

#endif
