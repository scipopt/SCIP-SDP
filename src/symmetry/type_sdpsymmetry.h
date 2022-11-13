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

/**@file   type_symmetry.h
 * @brief  type definitions for symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SDPSYMMETRY_H_
#define __SCIP_TYPE_SDPSYMMETRY_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** symmetry type specification */
#define SDPSYM_SPEC_INTEGER             UINT32_C(0x00000001)  /**< need symmetries for integer variables only */
#define SDPSYM_SPEC_BINARY              UINT32_C(0x00000002)  /**< need symmetries for binary variables only */
#define SDPSYM_SPEC_REAL                UINT32_C(0x00000004)  /**< need symmetries also for continuous variables */

typedef uint32_t SDPSYM_SPEC;           /**< types of variables handled by symmetry */

/** symmetry timings */
#define SDPSYM_COMPUTETIMING_BEFOREPRESOL    0  /**< compute symmetries before presolving */
#define SDPSYM_COMPUTETIMING_DURINGPRESOL    1  /**< compute symmetries during presolving */
#define SDPSYM_COMPUTETIMING_AFTERPRESOL     2  /**< compute symmetries after presolving */

/** define sense of rhs */
enum SDPSYM_Rhssense
{
   SDPSYM_SENSE_UNKOWN     = 0,              /**< unknown sense */
   SDPSYM_SENSE_INEQUALITY = 1,              /**< linear inequality */
   SDPSYM_SENSE_EQUATION   = 2,              /**< linear equation */
   SDPSYM_SENSE_XOR        = 3,              /**< XOR constraint */
   SDPSYM_SENSE_AND        = 4,              /**< AND constraint */
   SDPSYM_SENSE_OR         = 5,              /**< OR constrant */
   SDPSYM_SENSE_BOUNDIS_TYPE_1 = 6,          /**< bounddisjunction type 1 */
   SDPSYM_SENSE_BOUNDIS_TYPE_2 = 7           /**< bounddisjunction type 2 */
};
typedef enum SDPSYM_Rhssense SDPSYM_RHSSENSE;

/* type of symmetry handling codes */
#define SDPSYM_HANDLETYPE_NONE             UINT32_C(0x00000000)  /**< no symmetry handling */
#define SDPSYM_HANDLETYPE_SYMBREAK         UINT32_C(0x00000001)  /**< symmetry breaking inequalities */
#define SDPSYM_HANDLETYPE_ORBITALFIXING    UINT32_C(0x00000002)  /**< orbital fixing */
#define SDPSYM_HANDLETYPE_SST              UINT32_C(0x00000004)  /**< Schreier Sims cuts */
#define SDPSYM_HANDLETYPE_SYMCONS (SDPSYM_HANDLETYPE_SYMBREAK | SDPSYM_HANDLETYPE_SST)

typedef uint32_t SDPSYM_HANDLETYPE;          /**< type of symmetry handling */

typedef struct SDPSYM_Vartype SDPSYM_VARTYPE;      /**< data of variables that are considered to be equivalent */
typedef struct SDPSYM_Optype SDPSYM_OPTYPE;        /**< data of operators that are considered to be equivalent */
typedef struct SDPSYM_Consttype SDPSYM_CONSTTYPE;  /**< data of constants that are considered to be equivalent */
typedef struct SDPSYM_Rhstype SDPSYM_RHSTYPE;      /**< data of constraint sides that are considered to be equivalent */
typedef struct SDPSYM_Matrixdata SDPSYM_MATRIXDATA;/**< data for symmetry group computation on linear constraints */
typedef struct SDPSYM_Exprdata SDPSYM_EXPRDATA;    /**< data for symmetry group computation on nonlinear constraints */
typedef struct SDPSYM_Sdpdata SDPSYM_SDPDATA;      /**< data for symmetry group computation on SDP constraints */


#ifdef __cplusplus
}
#endif

#endif
