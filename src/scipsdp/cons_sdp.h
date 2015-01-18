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

/**@file   cons_sdp.h
 * @brief  constraint handler for SDPs
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 */

#ifndef __SCIP_CONSHDLR_SDP_H__
#define __SCIP_CONSHDLR_SDP_H__

#include "scip/scip.h"


#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for SDP constraints and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sort given arrays by nondecreasing rows and then cols */
EXTERN
SCIP_RETCODE SCIPconsSDPsortRowsCols(
   int*                  row,                /** array of row indices, will be returned in nondecreasing order */
   int*                  col,                /** array of col indices, for equal rows this will be the tiebreaker */
   SCIP_Real*            val,                /** this will be returned in the order of row and col */
   int                   maxrow,             /** maximum value in the row array (values should go from 0 to maxrow */
   int                   length              /** length of the arrays that should be sorted */
   );

/** creates an SDP-constraint */
SCIP_RETCODE SCIPcreateConsSdp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   nvars,              /**< number of variables in this SDP constraint */
   int                   nnonz,              /**< number of nonzeroes in this SDP constraint */
   int                   blocksize,          /**< size of this SDP-block */
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val point to */
   int**                 col,                /**< pointer to column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to values of the nonzeroes for each variable */
   SCIP_Var**            vars,               /**< SCIP_Variables present in this SDP constraint, ordered by their begvar-indices */
   int                   constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeroes */
   int*                  constrow,           /**< row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< values of the constant nonzeroes */
   );

/** get the data belonging to a single SDP-constraint
 *
 *  In arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length.
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown.
 */
EXTERN
SCIP_RETCODE SCIPconsSdpGetData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nvars,              /**< number of variables in this SDP constraint */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  blocksize,          /**< size of this SDP-block */
   int*                  arraylength,        /**< length of the given nvarnonz, col, row and val arrays, if this is too short this will return the needed length*/
   int*                  nvarnonz,           /**< number of nonzeros for each variable, also length of the arrays col/row/val are pointing to */
   int**                 col,                /**< pointers to column indices of the nonzeroes for each variable */
   int**                 row,                /**< pointers to row indices of the nonzeroes for each variable */
   SCIP_Real**           val,                /**< pointers to values of the nonzeroes for each variable */
   SCIP_Var**            vars,               /**< the SCIP variables present in this constraint, indexing equals indices in col/row/val */
   int*                  constnnonz,         /**< number of nonzeroes in the constant part of this SDP constraint, also length of the const arrays */
   int*                  constcol,           /**< pointer to column indices of the constant nonzeroes */
   int*                  constrow,           /**< pointer to row indices of the constant nonzeroes */
   SCIP_Real*            constval            /**< pointer to values of the constant nonzeroes */
   );

/** gets the number of nonzeroes and constant nonzeroes for this SDP constraint
 *
 *  Either nnonz or constnnonz may be NULL.
 */
EXTERN
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nnonz,              /**< number of nonzeroes in this SDP constraint */
   int*                  constnnonz          /**< number of nonzeroes in the constant part of this SDP constraint */
   );

/** gets the full constraint Matrix \f A_j \f for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int                   j,                  /**< the variable j to get the corresponding matrix \f A_j \f for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f A_j \f */
   );

/** gives an n*n-long array with the full constant matrix */
EXTERN
SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   SCIP_Real*            mat                 /**< pointer to store the full constant matrix */
   );

/** gives an 0.5*n*(n+1)-long array with the lower triangular part of the constant matrix indexed by compLowerTriangPos */
SCIP_RETCODE SCIPconsSdpGetLowerTriangConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   SCIP_Real*            mat                 /**< pointer to store the lower triangular part of the constant matrix */
   );

/** checks feasibility for a single SDP-Cone */
SCIP_RETCODE SCIPconsSdpCheckSdpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**< the solution to check feasibility for */
   SCIP_Bool             checkintegrality,   /**< has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< have current LP rows to be checked? */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the feasibility checking call */
   );

#ifdef __cplusplus
}
#endif

#endif
