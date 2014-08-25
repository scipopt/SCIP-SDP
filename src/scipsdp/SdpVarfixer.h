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

/**@file   SdpVarfixer.h
 * @brief  adds the main functionality to fix/unfix/(multi-)aggregate variables by merging two three-tuple-arrays of row/col/val together
 * @author Tristan Gally
 */

#ifndef __SDPVARMAPPER_H__
#define __SDPVARMAPPER_H__

#include "scip/type_misc.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * sort the given row, col and val arrays first by non-decreasing row-indices, than for those by identical row-indices by non-increasing val-indices
 */
EXTERN
void SdpVarfixerSortRowCol(
   int*                  row,                /* row indices */
   int*                  col,                /* column indices */
   int*                  val,                /* values */
   int                   length              /* length of the given arrays */
   )

/**
 * Merges two three-tuple-arrays together. The original arrays will be mulitplied with scalar and then merged into the target arrays. If there
 * is already an entry for a row/col combination, these two entries will be combined (their values added together), if they cancel each other out
 * the nonzero entry will be removed.
 */
EXTERN
SCIP_RETCODE SdpVarfixerMergeArrays(
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int*                  originrow,          /** original row-index-array that is going to be merged */
   int*                  origincol,          /** original column-index-array that is going to be merged */
   SCIP_Real*            originval,          /** original nonzero-values-array that is going to be merged */
   int                   origlength,         /** length of the original arrays */
   SCIP_Bool             originsorted,       /** are the origin arrays already sorted by non-decreasing row and in case of ties col */
   SCIP_Real             scalar,             /** scalar that the original nonzero-values will be multiplied with before merging */
   int*                  targetrow,          /** row-index-array the original array will be merged into */
   int*                  targetcol,          /** column-index-array the original array will be merged into */
   SCIP_real*            targetval,          /** nonzero-values-array the original array will be merged into */
   int*                  targetlength        /** length of the target arrays the original arrays will be merged into, this will be updated to the
                                               * new length after the mergings */
   );
