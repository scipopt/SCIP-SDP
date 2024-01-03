/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt,              */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2024 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2024 Zuse Institute Berlin                             */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_sdp.h
 * @brief  Constraint handler for SDP-constraints
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 *
 * Constraint handler for semidefinite constraints of the form \f$ \sum_{j=1}^n A_j y_j - A_0 \succeq 0 \f$,
 * where the matrices \f$A_j\f$ and \f$A_0\f$ need to be symmetric. Only the nonzero entries of the matrices
 * are stored.
 *
 */

#ifndef __SCIP_CONSHDLR_SDP_H__
#define __SCIP_CONSHDLR_SDP_H__

#include "scip/scip.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for SDP constraints and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates the handler for rank 1 SDP constraints and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeConshdlrSdpRank1(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates an SDP-constraint
 *
 *  The matrices should be lower triangular.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in this SDP constraint */
   int                   nnonz,              /**< number of nonzeros in this SDP constraint */
   int                   blocksize,          /**< size of this SDP-block */
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val point to */
   int**                 col,                /**< pointer to column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< SCIP_VARiables present in this SDP constraint that correspond to the indices in col/row/val */
   int                   constnnonz,         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeros */
   int*                  constrow,           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< values of the constant nonzeros */
   SCIP_Bool             removeduplicates    /**< Should duplicate matrix entries be removed (then order of col/row/val might change)? */
   );

/** creates a rank 1 SDP-constraint
 *
 *  The matrices should be lower triangular.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateConsSdpRank1(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in this SDP constraint */
   int                   nnonz,              /**< number of nonzeros in this SDP constraint */
   int                   blocksize,          /**< size of this SDP-block */
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val point to */
   int**                 col,                /**< pointer to column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< SCIP_VARiables present in this SDP constraint that correspond to the indices in col/row/val */
   int                   constnnonz,         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeros */
   int*                  constrow,           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< values of the constant nonzeros */
   SCIP_Bool             removeduplicates    /**< Should duplicate matrix entries be removed (then order of col/row/val might change)? */
   );

/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
SCIP_EXPORT
int SCIPconsSdpCompLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   );

/** get the data belonging to a single SDP-constraint
 *
 *  In arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length.
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nvars,              /**< pointer to store the number of variables in this SDP constraint */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  blocksize,          /**< pointer to store the size of this SDP-block */
   int*                  arraylength,        /**< length of the given nvarnonz, col, row and val arrays, if this is too short this will return the needed length*/
   int*                  nvarnonz,           /**< pointer to store the number of nonzeros for each variable, also length of the arrays col/row/val are
                                               *  pointing to */
   int**                 col,                /**< pointer to store the column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to store the row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to store the values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< pointer to store the SCIP variables present in this constraint that correspond to the indices in col/row/val */
   int*                  constnnonz,         /**< pointer to store the number of nonzeros in the constant part of this SDP constraint, also length of
                                               *  the const arrays */
   int*                  constcol,           /**< pointer to store the column indices of the constant nonzeros */
   int*                  constrow,           /**< pointer to store the row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< pointer to store the values of the constant nonzeros */
   SCIP_Bool*            rankone,            /**< pointer to store if matrix should be rank one (or NULL, if information not necessary) */
   int**                 maxevsubmat,        /**< pointer to store two row indices of 2x2 subdeterminant with maximal eigenvalue [-1,-1 if not yet computed] (or NULL, if information not necessary) */
   SCIP_Bool*            addedquadcons       /**< pointer to store if the quadratic 2x2-minor constraints already added (in the rank1-case) (or NULL, if information not necessary) */
   );

/** gets the number of nonzeros and constant nonzeros for this SDP constraint
 *
 *  Either nnonz or constnnonz may be NULL if only the other one is needed.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get number of nonzeros for */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  constnnonz          /**< pointer to store the number of nonzeros in the constant part of this SDP constraint */
   );

/** gets the number of variables of the SDP constraint */
SCIP_EXPORT
int SCIPconsSdpGetNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get number of variables for */
   );

/** gets the variables of the SDP constraint */
SCIP_EXPORT
SCIP_VAR** SCIPconsSdpGetVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get variables for */
   );

/** gets the blocksize of the SDP constraint */
SCIP_EXPORT
int SCIPconsSdpGetBlocksize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get blocksize for */
   );

/** gets the full constraint Matrix \f$ A_j \f$ for a given variable j */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   int                   j,                  /**< the variable j to get the corresponding matrix \f A_j \f for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f A_j \f */
   );

/** gives an n*n-long array with the full constant matrix */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   SCIP_Real*            mat                 /**< pointer to store the full constant matrix */
   );

/** gives a 0.5*n*(n+1)-long array with the lower triangular part of the constant matrix indexed by SCIPconsSdpCompLowerTriangPos */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetLowerTriangConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   SCIP_Real*            mat                 /**< pointer to store the lower triangular part of the constant matrix */
   );

/** Compute a heuristic guess for a good starting solution \f$ \lambda ^* \cdot I \f$.
 *
 *  The solution is computed as
 *  \f[
 *  \lambda^* = \max \Bigg\{S \cdot \max_{i \in [m]} \{|u_i|, |l_i|\} \cdot \max_{i \in [m]} \|A_i\|_\infty + \|C\|_\infty,
 *  \frac{\max_{i \in [m]} b_i}{S \cdot \min_{i \in [m]} \min_{j, \ell \in [n]} (A_i)_{j\ell} } \Bigg\},
 *  \f]
 *  where \f$ S = \frac{ | \text{nonzero-entries of all } A_i | }{0.5 \cdot \text{ blocksize } (\text{ blocksize } + 1)} \f$
 *  measures the sparsity of the matrices.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGuessInitialPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint to guess an initial point for */
   SCIP_Real*            lambdastar          /**< pointer to store the guess for the initial point */
   );

/** Gets maximum absolute entry of constant matrix \f$ A_0 \f$ */
SCIP_EXPORT
SCIP_Real SCIPconsSdpGetMaxConstEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< the constraint to get the maximum constant matrix entry for */
   );

/** Gets maximum absolute entry of all matrices \f$ A_i \f$ */
SCIP_EXPORT
SCIP_Real SCIPconsSdpGetMaxSdpCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< the constraint to get the maximum constant matrix entry for */
   );

/** Computes an upper bound on the number of nonzeros of the (dual) SDP matrix \f$ Z = \sum_{j=1}^n A_j y_j - A_0 \f$,
 *  this should be used to allocate enough memory before calling SCIPconsSdpComputeSparseSdpMatrix.
 *
 *  Upper bound is computed as \f$ \min \{ \sum_{v \leq m} \text{nvarnonz}(v) + \text{constnnonz}, n \cdot (n+1) / 2 \} \f$.
 */
SCIP_EXPORT
int SCIPconsSdpComputeUbSparseSdpMatrixLength(
   SCIP_CONS*            cons                /**< the constraint for which the Matrix should be assembled */
   );

/** Computes (dual) SDP matrix \f$ Z = \sum_{j=1}^n A_j y_j - A_0 \f$ and returns it in sparse format
 *  @note row, col and val should have memory allocated equal to SCIPconsSdpComputeUbSparseSdpMatrixLength(),
 *        if the memory is not sufficient, length will be set to -1 and an error will be thrown
 */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpComputeSparseSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**< the solution to assemble the matrix for */
   int*                  length,             /**< input: allocated memory for row/col/val arrays
                                               *  output: number of nonzeros of the matrix / length of row/col/val arrays */
   int*                  row,                /**< pointer to store row indices of SDP-matrix */
   int*                  col,                /**< pointer to store column indices of SDP-matrix */
   SCIP_Real*            val                 /**< pointer to store values of SDP-matrix */
   );

/** returns wheter matrix should be rank one */
SCIP_EXPORT
SCIP_Bool SCIPconsSdpShouldBeRankOne(
   SCIP_CONS*            cons                /**< the constraint for which the existence of a rank one constraint should be checked */
   );

/** returns two row indices of 2x2 subdeterminant with maximal eigenvalue [or -1,-1 if not available] */
SCIP_EXPORT
SCIP_RETCODE SCIPconsSdpGetMaxEVSubmat(
   SCIP_CONS*            cons,               /**< the constraint for which the existence of a rank one constraint should be checked */
   int**                 maxevsubmat         /**< pointer to store the two row indices of 2x2 subdeterminant with
                                                maximal eigenvalue [or -1,-1 if not available] */
   );

/** returns whether the quadratic 2x2-minor constraints are already added (in the rank1-case) */
SCIP_EXPORT
SCIP_Bool SCIPconsSdpAddedQuadCons(
   SCIP_CONS*            cons                /**< the constraint for which it should be checked whether the quadratic 2x2-minor constraints are already added (in the rank1-case) */
   );

#ifdef __cplusplus
}
#endif

#endif
