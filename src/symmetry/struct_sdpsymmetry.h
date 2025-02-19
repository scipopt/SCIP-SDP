/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2025 Discrete Optimization, TU Darmstadt               */
/*                                                                           */
/*                                                                           */
/* Licensed under the Apache License, Version 2.0 (the "License");           */
/* you may not use this file except in compliance with the License.          */
/* You may obtain a copy of the License at                                   */
/*                                                                           */
/*     http://www.apache.org/licenses/LICENSE-2.0                            */
/*                                                                           */
/* Unless required by applicable law or agreed to in writing, software       */
/* distributed under the License is distributed on an "AS IS" BASIS,         */
/* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  */
/* See the License for the specific language governing permissions and       */
/* limitations under the License.                                            */
/*                                                                           */
/*                                                                           */
/* Based on SCIP - Solving Constraint Integer Programs                       */
/* Copyright (C) 2002-2025 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_symmetry.h
 * @brief  structs for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_SDPSYMMETRY_H_
#define __SCIP_STRUCT_SDPSYMMETRY_H_

#include "scip/scip.h"
#include "symmetry/type_sdpsymmetry.h"
#include "scip/type_expr.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of variables that are considered to be equivalent */
struct SDPSYM_Vartype
{
   SCIP_Real             obj;                /**< objective of variable */
   SCIP_Real             lb;                 /**< lower bound of variable */
   SCIP_Real             ub;                 /**< upper bound of variable */
   SCIP_VARTYPE          type;               /**< type of variable */
   int                   nconss;             /**< number of conss a variable is contained in */
   int                   color;              /**< store color */
};

/** data of operators that are considered to be equivalent */
struct SDPSYM_Optype
{
   SCIP_EXPR*            expr;               /**< the underlying expression */
   int                   level;              /**< level of operator in its expression tree */
   int                   color;              /**< store color */
};

/** data of constants that are considered to be equivalent */
struct SDPSYM_Consttype
{
   SCIP_Real             value;              /**< value of constant */
   int                   color;              /**< store color */
};

/** data of coefficients that are considered to be equivalent */
struct SDPSYM_Rhstype
{
   SCIP_Real             lhs;                /**< value of left-hand-side */
   SCIP_Real             rhs;                /**< value of right-hand-side */
   int                   color;              /**< store color */
};

/** data for symmetry group computation on linear constraints */
struct SDPSYM_Matrixdata
{
   SCIP_Real*            matcoef;            /**< nonzero coefficients appearing in the matrix */
   SCIP_Real*            rhscoef;            /**< rhs coefficients */
   SDPSYM_RHSSENSE*      rhssense;           /**< sense of rhs */
   int*                  matrhsidx;          /**< indices of rhs corresponding to matrix entries */
   int*                  matvaridx;          /**< indices of variables for matrix entries */
   int*                  matidx;             /**< indices in mat(rhs/var)idx array corresponding to matrix coefficients */
   int*                  rhsidx;             /**< indices in rhstype array corresponding to rhs coefficients */
   int*                  permvarcolors;      /**< array for storing the colors of the individual variables */
   int*                  matcoefcolors;      /**< array for storing the colors of all matrix coefficients */
   int*                  rhscoefcolors;      /**< array for storing the colors of all rhs coefficients */
   SCIP_VAR**            permvars;           /**< variables on which permutations act */
   int                   npermvars;          /**< number of variables for permutations */
   int                   nmatcoef;           /**< number of coefficients in matrix */
   int                   nrhscoef;           /**< number of coefficients in rhs */
   int                   nmaxmatcoef;        /**< maximal number of matrix coefficients (will be increase on demand) */
   int                   nuniquevars;        /**< number of unique variable types */
   int                   nuniquerhs;         /**< number of unique rhs types */
   int                   nuniquemat;         /**< number of unique matrix coefficients */
};

/** data for symmetry group computation on nonlinear constraints */
struct SDPSYM_Exprdata
{
   int                   nuniqueconstants;   /**< number of unique constants */
   int                   nuniqueoperators;   /**< number of unique operators */
   int                   nuniquecoefs;       /**< number of unique coefficients */
};

/** data for symmetry group computation on SDP onstraints */
struct SDPSYM_Sdpdata
{
   int                   nsdpconss;          /**< number of SDP constraints */
   int*                  blocksizes;         /**< for each SDP cons, the size of the SDP-block */
   int*                  nvars;              /**< for each SDP cons, the number of variables in this cons */
   SCIP_Real**           vals;               /**< for each SDP cons, the values of nonzeros for each variable block */
   SCIP_VAR***           vars;               /**< for each SDP cons, the variables present in the cons */
   int**                 valsbegins;         /**< for each SDP cons, the begin positions of new matrix in vals */
   int**                 cols;               /**< for each SDP cons, the column index of value in vals */
   int**                 rows;               /**< for each SDP cons, the row index of value in vals */
   SCIP_Real**           constvals;          /**< for each SDP cons, the values of nonzeros for each constant block */
   int*                  nconstvals;         /**< for each SDP cons, the number of nonzeros in the constant block */
   int**                 constcols;          /**< for each SDP cons, the column index of value in constvals */
   int**                 constrows;          /**< for each SDP cons, the row index of value in constvals */
   int**                 colors;             /**< for each SDP cons, the color of the coefficients */
   int**                 colors2;            /**< for each SDP cons, the second color of the coefficients */
   int**                 constcolors;        /**< for each SDP cons, the color of constant coefficients */
   int**                 constcolors2;       /**< for each SDP cons, the second color of constant coefficients */
   int                   lastcolorused;      /**< last used color */
   SCIP_Real*            minvals;            /**< for each variable block its minimum value in vals */
   SCIP_Real*            maxvals;            /**< for each variable blick its maximum value in vals */
   int                   nsdpvarblocks;      /**< number of variable blocks stored in minvals/maxvals */
   int                   maxsdpvarblocks;    /**< maximum number of variable blocks (length of minvals/maxvals) */
};


#ifdef __cplusplus
}
#endif

#endif
