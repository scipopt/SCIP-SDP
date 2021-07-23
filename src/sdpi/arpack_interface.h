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

/**@file   arpack_interface.h
 * @brief  interface methods for eigenvector computation using arpack
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_ARPACK_INTERFACE_H__
#define __SCIP_ARPACK_INTERFACE_H__

#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"

#ifdef __cplusplus
extern "C" {
#endif

/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK */
SCIP_EXPORT
SCIP_RETCODE SCIParpackComputeSmallestEigenvector(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real*            A,                  /**< matrix for which eigenvalues should be computed - will be destroyed! */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< array for eigenvector */
   );

/** computes an eigenvector for the smallest eigenvalue of a symmetric matrix using ARPACK; specialized sparse version for mu A + B */
SCIP_EXPORT
SCIP_RETCODE SCIParpackComputeSmallestEigenvectorOneVar(
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   int                   n,                  /**< size of matrix */
   SCIP_Real             mu,                 /**< scaling factor for A matrix */
   int                   annonz,             /**< number of nonzero elements in A */
   int*                  arow,               /**< array of row-indices of A */
   int*                  acol,               /**< array of column-indices of A */
   SCIP_Real*            aval,               /**< array of nonzero values of entries of A */
   int                   bnnonz,             /**< number of nonzero elements in B */
   int*                  brow,               /**< array of row-indices of nonzero matrix entries in B */
   int*                  bcol,               /**< array of column-indices of nonzero matrix entries in B*/
   SCIP_Real*            bval,               /**< array of nonzero values in B */
   SCIP_Real*            eigenvalue,         /**< pointer to store eigenvalue */
   SCIP_Real*            eigenvector         /**< array for eigenvector */
   );

#ifdef __cplusplus
}
#endif

#endif
