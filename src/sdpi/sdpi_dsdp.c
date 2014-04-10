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

#define SCIP_DEBUG

/**@file   sdpi_dsdp.c
 * @brief  interface for dsdp
 * @author Tristan Gally
 */

#include <assert.h>
#include "sdpi/sdpi.h"
#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */
#include "dsdpmem.h"                         /* for DSDPCALLOC2, DSDPFREE */
#include "scip/scip.h"

/* calls a DSDP-Function and transforms the return-code to a SCIP_ERROR if needed */
#define DSDP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if( (_dsdperrorcode_ = (x)) != 0 )                                                  \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/* same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced */
#define DSDP_CALLM(x)   do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if( (_dsdperrorcode_ = (x)) != 0 )                                                  \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_NOMEMORY;                                                            \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

#define IS_POSINF(x) ((x) >= SCIP_DEFAULT_INFINITY)
#define IS_NEGINF(x) ((x) <= -SCIP_DEFAULT_INFINITY)

struct SCIP_SDPi
{
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object to add sdp-constraints to */
   LPCone                lpcone;             /**< lpcone-object to add lp-constraints to */
   BCone                 bcone;              /**< bcone-object to add variable bounds to */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   int                   sdpid;              /**< identifier for debug-messages */
   int                   nvars;              /**< number of variables */
   SCIP_Real*            obj;                /**< objective function values of variables */
   SCIP_Real*            lb;                 /**< lower bounds of variables */
   SCIP_Real*            ub;                 /**< upper bounds of variables */
   int                   nsdpblocks;         /**< number of SDP-blocks */
   int*                  sdpblocksizes;      /**< sizes of the SDP-blocks */
   int                   sdpconstnnonz;      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   int*                  sdpconstbegblock;   /**< start index of each block in sdpconstval-array */
   int*                  sdpconstrowind;     /**< row-index for each entry in sdpconstval-array */
   int*                  sdpconstcolind;     /**< column-index for each entry in sdpconstval-array */
   SCIP_Real*            sdpconstval;        /**< values of entries of constant matrices in SDP-Block */
   int                   sdpnnonz;           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpbegvarblock;     /**< entry j*nvars + i is the start index of matrix \f A_i^j \f in sdpval,
                                              *   particularly entry i*nvars gives the starting point of block j */
   int*                  sdprowind;          /**< row-index for each entry in sdpval-array */
   int*                  sdpcolind;          /**< column-index for each entry in sdpval-array */
   SCIP_Real*            sdpval;             /**< values of SDP-constraint matrix entries */
   int                   nlpcons;            /**< number of LP-constraints */
   SCIP_Real*            lprhs;              /**< right hand sides of LP rows */
   int                   lpnnonz;            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprowind;           /**< row-index for each entry in lpval-array (going from 1 to nlpcons) */
   int*                  lpcolind;           /**< column-index for each entry in lpval-array (going from 1 to nvars) */
   SCIP_Real*            lpval;               /**< values of LP-constraint matrix entries */

};

static int nextsdpid     =  1;               /**< used to give ids to the generated sdps for debugging messages */


/*
 * Local Functions
 */

/** for given row and column (i,j) (indices going from 1 to n) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1
 */
static int compLowerTriangPos(
   int         i,                            /**< row index */
   int         j                             /**< column index */
   )
{
   /* the formula is for a position in the lower triangular part, if a position in the upper triangular part is given,
    * switch row and column (all matrices must be symmetric) */
   if ( i < j )
   {
      int temp;
      temp = i;
      i = j;
      j = temp;
   }
   assert(i >= 1);
   assert(j >= 1);

   return i*(i-1)/2 + j - 1;
}


/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver, getting version doesn't seem to be supported by DSDP */
const char* SCIPsdpiGetSolverName(
   void
   )
{
   return "DSDP";
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   void
   )
{
   return "Dual-Scaling Interior Point Solver for Semidefinite Programming developed by Steve Benson, Yinyu Ye, and Xiong Zhang (http://www.mcs.anl.gov/hs/software/DSDP/)";
}

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                 /**< pointer to an SDP interface structure */
   )
{
   return (void*) sdpi->dsdp;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,         /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   DSDP newdsdp;
   SDPCone newsdpcone;
   LPCone newlpcone;
   BCone newbcone;

   /* these will be properly initialized only immediatly prior to solving because DSDP and the SDPCone need information about the number
    * of variables and sdpblocks during creation */
   newdsdp=0;
   newsdpcone=0;
   newlpcone=0;
   newbcone=0;

   SCIPdebugMessage("Calling SCIPsdpiCreate (%d)\n",nextsdpid);

   assert ( sdpi != NULL );

   BMSallocBlockMemory(blkmem, sdpi);

   (*sdpi)->messagehdlr = messagehdlr;
   (*sdpi)->blkmem = blkmem;
   (*sdpi)->dsdp = newdsdp;
   (*sdpi)->sdpcone = newsdpcone;
   (*sdpi)->lpcone = newlpcone;
   (*sdpi)->bcone = newbcone;
   (*sdpi)->sdpid = nextsdpid++;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP interface structure */
   )
{
   SCIPdebugMessage("Calling SCIPsdpiFree (%d)\n",(*sdpi)->sdpid);
   assert(sdpi != NULL);
   assert(*sdpi != NULL);

   DSDP_CALL(DSDPDestroy((*sdpi)->dsdp));

   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->obj), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lb), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->ub), (*sdpi)->nvars);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpblocksizes), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstbegblock), (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstrowind), (*sdpi)->sdpconstnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstcolind), (*sdpi)->sdpconstnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpconstval), (*sdpi)->sdpconstnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpbegvarblock), (*sdpi)->nvars * (*sdpi)->nsdpblocks);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdprowind), (*sdpi)->sdpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpcolind), (*sdpi)->sdpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->sdpval), (*sdpi)->sdpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprhs), (*sdpi)->nlpcons);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lprowind), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpcolind), (*sdpi)->lpnnonz);
   BMSfreeBlockMemoryArray((*sdpi)->blkmem, &((*sdpi)->lpval), (*sdpi)->lpnnonz);

   BMSfreeBlockMemory((*sdpi)->blkmem, sdpi);

   return SCIP_OKAY;
}

/**@} */


/*
 * Modification Methods
 */

/**@name Modification Methods */
/**@{ */

/** copies SDP data into SDP solver
 *
 *  @note as the SDP-constraint matrices are symmetric, only the upper triangular part of them must be specified
 */
SCIP_RETCODE SCIPsdpiLoadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   int*                  sdpblocksizes,      /**< sizes of the SDP-blocks */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstbegblock,   /**< start index of each block in sdpconstval-array */
   const int*            sdpconstrowind,     /**< row-index for each entry in sdpconstval-array */
   const int*            sdpconstcolind,     /**< column-index for each entry in sdpconstval-array */
   const SCIP_Real*      sdpconstval,        /**< values of entries of constant matrices in SDP-Block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   const int*            sdpbegvarblock,     /**< entry j*nvars + i is the start index of matrix \f A_i^j \f in sdpval,
                                              *   particularly entry j*nvars gives the starting point of block j (with numbering starting at 0) */
   const int*            sdprowind,          /**< row-index for each entry in sdpval-array (going from 1 to the blocksize of that block)*/
   const int*            sdpcolind,          /**< column-index for each entry in sdpval-array (going from 1 to the blocksize of that block)*/
   const SCIP_Real*      sdpval,             /**< values of SDP-constraint matrix entries */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (going from 1 to nlpcons) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (going from 1 to nvars) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries */
   )
{
   int i;

   /* copy all inputs into the corresponding sdpi-parameters to later put them into DSDP prior to solving when the final
    * number of blocks and variables are known
    */

   sdpi->nvars = nvars;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), nvars);
   for (i = 0; i < nvars; i++)
   {
      (sdpi->obj)[i] = obj[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), nvars);
   for (i = 0; i < nvars; i++)
   {
      (sdpi->lb)[i] = lb[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), nvars);
   for (i = 0; i < nvars; i++)
   {
      (sdpi->ub)[i] = ub[i];
   }
   sdpi->nsdpblocks = nsdpblocks;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), nsdpblocks);
   for (i = 0; i < nsdpblocks; i++)
   {
      (sdpi->sdpblocksizes)[i] = sdpblocksizes[i];
   }
   sdpi->sdpconstnnonz = sdpconstnnonz;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), nsdpblocks);
   for (i = 0; i < nsdpblocks; i++)
   {
      (sdpi->sdpconstbegblock)[i] = sdpconstbegblock[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpconstnnonz);
   for (i = 0; i < sdpconstnnonz; i++)
   {
      (sdpi->sdpconstrowind)[i] = sdpconstrowind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpconstnnonz);
   for (i = 0; i < sdpconstnnonz; i++)
   {
      (sdpi->sdpconstcolind)[i] = sdpconstcolind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpconstnnonz);
   for (i = 0; i < sdpconstnnonz; i++)
   {
      (sdpi->sdpconstval)[i] = sdpconstval[i];
   }
   sdpi->sdpnnonz = sdpnnonz;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), nvars * nsdpblocks);
   for (i = 0; i < nvars * nsdpblocks; i++)
   {
      (sdpi->sdpbegvarblock)[i] = sdpbegvarblock[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpnnonz);
   for (i = 0; i < sdpnnonz; i++)
   {
      (sdpi->sdprowind)[i] = sdprowind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpnnonz);
   for (i = 0; i < sdpnnonz; i++)
   {
      (sdpi->sdpcolind)[i] = sdpcolind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpnnonz);
   for (i = 0; i < sdpnnonz; i++)
   {
      (sdpi->sdpval)[i] = sdpval[i];
   }
   sdpi->nlpcons = nlpcons;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), nlpcons);
   for (i = 0; i < nlpcons; i++)
   {
      (sdpi->lprhs)[i] = lprhs[i];
   }
   sdpi->lpnnonz = lpnnonz;
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), lpnnonz);
   for (i = 0; i < lpnnonz; i++)
   {
      (sdpi->lprowind)[i] = lprowind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), lpnnonz);
   for (i = 0; i < lpnnonz; i++)
   {
      (sdpi->lpcolind)[i] = lpcolind[i];
   }
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), lpnnonz);
   for (i = 0; i < lpnnonz; i++)
   {
      (sdpi->lpval)[i] = lpval[i];
   }

   return SCIP_OKAY;
}

/** adds another SDP-Block to the problem
 *
 *  @note as \f A_i^j \f is symmetric, only the lower triangular part of it must be specified
 */
SCIP_RETCODE SCIPsdpiAddSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrices \f A_i^j \f */
   int                   nnonz,              /**< sum of non-zeroes in the lower triagonal parts of the \f A_i^j \f */
   const int*            constbegrow,        /**< start index of the rows of the matrix \f A_i^0 \f */
   const int*            constcolind,        /**< the column-indices of the non-zero entries of \f A_i^0 \f */
   const SCIP_Real*      constval,           /**< the values of \f A_i^0 \f as specified by constbegrow and constcolind */
   const int*            begvar,             /**< start index of the matrix \f A_i^j \f for each j */
   const int*            rowind,             /**< the row-indices of the non-zero entries of \f A_i^j \f */
   const int*            colind,             /**< the column-indices of the non-zero entries of \f A_i^j \f */
   const SCIP_Real*      val                 /**< the values of of \f A_i^j \f as specified by begvar, rowind and colind */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** deletes a SDP-Block */
SCIP_RETCODE SCIPsdpiDelSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block              /**< index of the SDP-Block that should be deleted */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** adds additional variables to the SDP
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
SCIP_RETCODE SCIPsdpiAddVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to be added */
   const SCIP_Real*      obj,                /**< objective function values of new variables */
   const SCIP_Real*      lb,                 /**< lower bounds of new variables */
   const SCIP_Real*      ub,                 /**< upper bounds of new variables */
   char**                varnames,           /**< variable names, or NULL */
   int                   sdpnnonz,           /**< number of nonzero elements to be added to the SDP constraint matrices */
   const int*            sdpbegvar,          /**< start index of each variable in sdpval-array, or NULL if sdpnnonz == 0 */
   const int*            sdprowind,          /**< row indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const int*            sdpcolind,          /**< col indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const SCIP_Real*      sdppval,            /**< values of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   int                   lpnnonz,            /**< number of nonzero elements to be added to the LP constraint matrices */
   const int*            lpbegvar,           /**< start index of each variable in lpval-array, or NULL if lpnnonz == 0 */
   const int*            lprowind,           /**< row indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const int*            lpcolind,           /**< col indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const SCIP_Real*      lpval               /**< values of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** deletes all variables in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to be deleted */
   int                   lastvar             /**< last variable to be deleted */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** deletes variables from SCIP_SDPI; the new position of a variable must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelColset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of variables
                                              *   input:  1 if variable should be deleted, 0 if not
                                              *   output: new position of variable, -1 if variable was deleted */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** adds rows to the LP-Block
 *
 *  @note ind array is not checked for duplicates, problems may appear if indeces are added more than once
 */
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** deletes all rows in the given range from LP-Block */
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** deletes LP rows from SCIP_SDPI; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of LP rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes lower and upper bounds of variables */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes right hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes a single coefficient in LP constraint matrix */
SCIP_RETCODE SCIPsdpiChgLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of LP-coefficient to change */
   int                   col,                /**< column number of LP-coefficient to change */
   SCIP_Real             newval              /**< new value of LP-coefficient */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes objective values of variables in the SDP */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of variables to change objective value for */
   int*                  ind,                /**< variable indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for variables */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes a single coefficient in constant matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPConstCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   rowind,             /**< row index */
   int                   colind,             /**< column index */
   const SCIP_Real*      val                 /**< new value of given entry of constant matrix in given SDP-Block */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** changes a single coefficient in a constraint matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index */
   int                   var,                /**< variable index */
   int                   rowind,             /**< row index */
   int                   colind,             /**< column index */
   const SCIP_Real*      val                 /**< new value of given entry of the give constraint matrix in specified SDP-Block */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/* @todo: # rows -> SDP constraints
 *         return dimension
 */
/** gets the number of LP-rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets the number of SDP-Blocks in the SDP */
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of rows */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets the number of variables in the SDP */
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets the number of nonzero elements in the SDP constraint matrices */
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP constraint matrcies */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets the number of nonzero elements in the LP Matrix */
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either sdpnnonz, sdpbegblock, sdprowind, sdpcolind and sdpval have to be NULL, or all of them have to be non-NULL, the same is true for the lp-part.
 */
SCIP_RETCODE SCIPsdpiGetVarInfos(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to extract information for */
   int                   lastvar,            /**< last variable to extract information for */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  sdpnnonz,           /**< pointer to store the number of nonzero elements of the sdp-constraints returned, or NULL */
   int*                  sdpbegblock,        /**< buffer to store start index of each block in sdpval-array, or NULL */
   int*                  sdprowind,          /**< buffer to store row indices of sdp constraint matrix entries, or NULL */
   int*                  sdpcolind,          /**< buffer to store column indices of sdp constraint matrix entries, or NULL */
   SCIP_Real*            sdpval,             /**< buffer to store values of sdp constraint matrix entries, or NULL */
   int*                  lpnnonz,            /**< pointer to store the number of nonzero elements of the lp-constraints returned, or NULL */
   int*                  lpbegrow,           /**< buffer to store start index of each row in lpval-array, or NULL */
   int*                  lpcolind,           /**< buffer to store column indices of lp constraint matrix entries, or NULL */
   SCIP_Real*            lpval               /**< buffer to store values of lp constraint matrix entries, or NULL */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets LP rows from SDP problem object; the arrays have to be large enough to store all values.
 *  rhs can be null,
 *  either nnonz, begrow, colind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPsdpiGetLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first LP row to get from SDP */
   int                   lastrow,            /**< last LP row to get from SDP */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begrow,             /**< buffer to store start index of each row in colind- and val-array, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets a number SDP blocks; the arrays have to be large enough to store all values.
 *  either constnnonz, constbegblock, constrow, constcolind, and constval have to be NULL, or all of them have to be non-NULL, same for the non-constant parts.
 */
SCIP_RETCODE SCIPsdpiGetSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstblock,         /**< first SDP block to get from SDP */
   int                   lastblock,          /**< last SDP block to get from SDP */
   int*                  constnnonz,         /**< pointer to store the number of nonzero elements returned for the constant part, or NULL */
   int*                  constbegblock,      /**< buffer to store the start indices of the different blocks in the constval-array, or NULL */
   int*                  constrowind,        /**< buffer to store row indices of entries of constant matrices, or NULL */
   int*                  constcolind,        /**< buffer to store column indices of entries of constant matrices, or NULL */
   SCIP_Real*            constval,           /**< buffer to store values of entries of constant matrices, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begblock,           /**< buffer to store the start indices of the different blocks in the val-array, or NULL */
   int*                  varind,             /**< buffer to store variable indices of constraint matrix entries, meaning the i in \f A_i^j \f, in val-array, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for */
   int                   lastvar,            /**< last variable to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets current variable bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get bounds for */
   int                   lastvar,            /**< last variable to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets current right hand sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets a single coefficient of LP constraint matrix */
SCIP_RETCODE SCIPsdpiGetLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets a single coefficient of constant SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPConstCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets a single coefficient of SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient */
   int                   var,                /**< variable index of coefficient, meaning the i in \f A_i^j \f, in val-array, or NULL */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}


/**@} */




/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   int* dsdpconstind;         /* indices for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpconstval;      /* non-zero values for constant SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdpind;              /* indices for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdpval;          /* non-zero values for SDP-constraint-matrices, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplpbegcol;         /* starting-indices for all columns in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int* dsdplprowind;         /* row indices in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   double* dsdplpval;         /* nonzeroes in LP, needs to be stored for DSDP during solving and be freed only afterwards */
   int i;
   int ind;

#ifdef SCIP_DEBUG
   DSDPTerminationReason* reason; /* this will later be used to check if DSDP converged */
   BMSallocBlockMemory(sdpi->blkmem, &reason);
#endif

   /* insert data */

   SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", sdpi->sdpid);

   DSDP_CALLM(DSDPCreate(sdpi->nvars, &(sdpi->dsdp)));
   DSDP_CALLM(DSDPCreateSDPCone(sdpi->dsdp, sdpi->nsdpblocks, &(sdpi->sdpcone)));
   DSDP_CALLM(DSDPCreateLPCone(sdpi->dsdp, &(sdpi->lpcone)));
   DSDP_CALLM(DSDPCreateBCone(sdpi->dsdp, &(sdpi->bcone)));

   DSDP_CALLM(BConeAllocateBounds(sdpi->bcone,2*sdpi->nvars)); /*allocate memory for lower and upper bounds */

   for (i = 0; i < sdpi->nvars; i++)
   {
      DSDP_CALL(DSDPSetDualObjective(sdpi->dsdp, i+1, -1 * sdpi->obj[i])); /*insert objective value, DSDP counts from 1 to n instead of 0 to n-1,
                                                                                               * *(-1) because DSDP maximizes instead of minimizing */
      DSDP_CALL(BConeSetLowerBound(sdpi->bcone, i+1, sdpi->lb[i])); /*insert lower bound, DSDP counts from 1 to n instead of 0 to n-1 */
      DSDP_CALL(BConeSetUpperBound(sdpi->bcone, i+1, sdpi->ub[i])); /*insert upper bound, DSDP counts from 1 to n instead of 0 to n-1 */
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpi->bcone);
#endif


   for(i = 0; i < sdpi->nsdpblocks; i++)
   {
      DSDP_CALL(SDPConeSetBlockSize(sdpi->sdpcone, i, sdpi->sdpblocksizes[i])); /*set the blocksizes (blocks are counted from 0 to m-1) */
   }


   /*start inserting the constant matrix*/
   if ( sdpi->nsdpblocks > 0 && sdpi->sdpconstnnonz > 0 )
   {
      int block;
      /*allocate memory*/
      /*This needs to be one long array, because DSDP uses it for solving, so all nonzeros have to be in it, and it may not be freed before the problem is solved. */

      /*indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpconstind, sdpi->sdpconstnnonz);
      /*values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpconstval, sdpi->sdpconstnnonz);

      for(block = 0; block < sdpi->nsdpblocks; block++)
      {
         int blocknnonz;   /*number of nonzeroes in the constant matrix for the current block */

         if ( block == sdpi->nsdpblocks - 1 )
            blocknnonz = sdpi->sdpconstnnonz - sdpi->sdpconstbegblock[block];
         else
            blocknnonz = sdpi->sdpconstbegblock[block + 1] - sdpi->sdpconstbegblock[block]; /* difference between first index of next block and first index of current block
                                                                                             * gives the number of non-zeroes in the current block */

         for(i = 0; i < blocknnonz; i++)
         {
            ind = sdpi->sdpconstbegblock[block] + i; /* compute the current position in the nonzero-arrays */
            dsdpconstind[ind] = compLowerTriangPos(sdpi->sdpconstrowind[ind], sdpi->sdpconstcolind[ind]);
            dsdpconstval[ind] = -1 * sdpi->sdpconstval[ind]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j)
                                                                * should be positive semidefinite */
         }

         DSDP_CALL(SDPConeSetASparseVecMat(sdpi->sdpcone, block, 0, sdpi->sdpblocksizes[block], 1, 0, dsdpconstind + sdpi->sdpconstbegblock[block],
            dsdpconstval + sdpi->sdpconstbegblock[block], blocknnonz));   /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block
                                                                           * by adding sdpi->sdpconstbegblock[block] */
      }
   }


   /*start inserting the other SDP-Constraint-Matrices */
   if(sdpi->nsdpblocks > 0 && sdpi->sdpnnonz > 0)
   {
      int block;
      int var;
      int k;

      /*allocate memory */
      /*This needs to be one long array, because DSDP uses it for solving so all nonzeros have to be in it and it may not be freed before the problem is solved. The distinct blocks/variables
       *(for the i,j-parts) are then given by dsdpind + sdpbegvarblock[nvars * block + var], which gives a pointer to the first array-element belonging to this block and then the number of
       *elements in this block is given to DSDP for iterating over it */

      /*indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpind, sdpi->sdpnnonz);
      /*values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpval, sdpi->sdpnnonz);

      for(block = 0; block < sdpi->nsdpblocks; block++)
      {
         for(var = 0; var < sdpi->nvars; var++)
         {
            int aij_nnonz; /* number of nonzeroes in Matrix Aij (current block and current variable) */
            ind = sdpi->nvars * block + var; /* compute current position in the sdpbegvarblock-array */
            if ( ind == (sdpi->nsdpblocks * sdpi->nvars) - 1 )
               aij_nnonz = sdpi->sdpnnonz - sdpi->sdpbegvarblock[ind];
            else
               aij_nnonz = sdpi->sdpbegvarblock[ind + 1] - sdpi->sdpbegvarblock[ind];

            for (k = 0; k < aij_nnonz; k++)
            {
               ind = sdpi->sdpbegvarblock[sdpi->nvars * block +var] + k; /* compute current position in the nonzero-arrays */
               dsdpind[ind] = compLowerTriangPos(sdpi->sdprowind[ind], sdpi->sdpcolind[ind]);
               dsdpval[ind] = -1 * sdpi->sdpval[ind];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0)
                                                                                                                                                * should be positive semidefinite */
            }
            DSDP_CALL(SDPConeSetASparseVecMat(sdpi->sdpcone, block, var + 1, sdpi->sdpblocksizes[block], 1, 0, dsdpind + sdpi->sdpbegvarblock[sdpi->nvars * block + var],
               dsdpval + sdpi->sdpbegvarblock[sdpi->nvars * block + var], aij_nnonz)); /* var+1 is needed because DSDP indexes the vars from 1 to nvars (var 0 is the constant matrix), adding
                                                         * sdpi->sdpbegvarblock[sdpi->nvars * block + var] shifts the arrays to the first nonzero belonging to this block and this variable */
         }
      }
      #ifdef SCIP_DEBUG
      SDPConeView2(sdpi->sdpcone);
      #endif
   }


   /*start inserting the LP constraints */
   if(sdpi->nlpcons > 0 && sdpi->lpnnonz > 0)
   {
      int* sortedlpcolind;
      int column;
      int constraint;

      /*memory allocation */

      /*these arrays are needed in DSDP during solving, so they may only be freed afterwards */
      /*lpbegcol[i] gives the number of nonzeroes in column 0 (right hand side) till i-1 (i going from 1 till n, with extra entries 0 (always 0) and n+1 (always lpcons + lpnnonz)) */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdplpbegcol, sdpi->nvars + 2);
      /*dsdplprowind saves the row indices of the LP for DSDP */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdplprowind, sdpi->nlpcons + sdpi->lpnnonz); /*length is lpnnonz + nlpcons, because right hand sides are also included in the vector */
      /*values given to DSDP */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdplpval, sdpi->nlpcons + sdpi->lpnnonz); /*length is lpnnonz + nlpcons, because right hand sides are also included in the vector */

      /*compute lpbegcol */

      /*to compute lpbegcol the column indices need to be sorted, for this they need to be copied into a non-constant array */
      SCIP_ALLOC(BMSallocMemoryArray(&sortedlpcolind, sdpi->lpnnonz));

      for(i = 0; i < sdpi->lpnnonz; i++)
      {
         sortedlpcolind[i] = sdpi->lpcolind[i];
         dsdplprowind[sdpi->nlpcons + i] = sdpi->lprowind[i] - 1;  /*the first nlpcons entries will be used for the right hand sides, so the matrix-entries are copied in the later ones, rowindices in DSDP start at 0 instead of 1 */
         dsdplpval[sdpi->nlpcons + i] = -1 * sdpi->lpval[i];   /*the first nlpcons entries will be used for the right hand sides, so the matrix-entries are copied in the later ones, *(-1) is needed, because
                                                          *DSDP wants <= instead of >= */
      }

      SCIPsortIntIntReal(sortedlpcolind, dsdplprowind + sdpi->nlpcons, dsdplpval + sdpi->nlpcons, sdpi->lpnnonz); /* all three arrays should now be sorted by non-decreasing column-indizes, for dsdplprowind and dsdplpval
      the sorting starts at position nlpcons (the first index is shifted by nlpcons), because the earlier entries are still empty and will only later be used for the right hand sides */

      dsdplpbegcol[0]=0;
      dsdplpbegcol[1]=sdpi->nlpcons; /*the first nlpcons indices are used to save the right hand sides */
      ind=0; /* this will be used for traversing the sortedlpcolind-array */

      for(column = 1; column < sdpi->nvars + 1; column++) /*columns are indexed 1 to nvars */
      {
         dsdplpbegcol[column+1]=dsdplpbegcol[column];  /*each new column can't start before the last one */
         while(ind < sdpi->lpnnonz && sortedlpcolind[ind] == column) /* look at all indices whose column index matches the current column */
         {
            dsdplpbegcol[column+1]++;   /*for each element with (lpcolind = current column) an additional entry in dsdplpval is needed, so the next column starts one spot later */
            ind++;
         }
      }
      assert(dsdplpbegcol[sdpi->nvars+1] == sdpi->lpnnonz + sdpi->nlpcons);

      BMSfreeMemoryArray(&sortedlpcolind); /*this was only needed to sort the column indices and compute dsdplpbegcol */

      for(column = 1; column < sdpi->nvars + 1; column++)
      {
         SCIPsortIntReal(dsdplprowind + dsdplpbegcol[column], dsdplpval + dsdplpbegcol[column], dsdplpbegcol[column+1] - dsdplpbegcol[column]);
         /*sort all the entries belonging to the same column by their row numbers */
      }

      for(constraint = 0; constraint < sdpi->nlpcons; constraint++)
      {
         dsdplprowind[constraint] = constraint; /*the row index of each constraint is the index of the constraint */
         dsdplpval[constraint] = -1 * sdpi->lprhs[constraint]; /*insert rhs values, *(-1) is needed, because DSDP wants <= instead of >= */
      }

      DSDP_CALL(LPConeSetData(sdpi->lpcone, sdpi->nlpcons, dsdplpbegcol, dsdplprowind, dsdplpval));
      #ifdef SCIP_DEBUG
      LPConeView(sdpi->lpcone);
      #endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpi->sdpid);
   DSDP_CALLM(DSDPSetup(sdpi->dsdp));
   DSDP_CALL(DSDPSolve(sdpi->dsdp));

   DSDP_CALL(DSDPComputeX(sdpi->dsdp)); /*computes X and determines feasibility and unboundedness of the solution */

   /*these arrays were used to give information to DSDP and were needed during solving and for computing X, so they may only be freed now*/
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdpconstind, sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdpconstval, sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdpind, sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdpval, sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdplpbegcol, sdpi->nvars +2);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdplprowind, sdpi->nlpcons + sdpi->lpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &dsdplpval, sdpi->nlpcons + sdpi->lpnnonz);

#ifdef SCIP_DEBUG
   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   switch ( *reason )
   {
      case DSDP_CONVERGED:
         SCIPdebugMessage("DSDP converged!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_INFEASIBLE_START:
         SCIPdebugMessage("DSDP started with an infeasible point!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_SMALL_STEPS:
         SCIPdebugMessage("Short step lengths created by numerical difficulties prevented progress in DSDP!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_INDEFINITE_SCHUR_MATRIX:
         SCIPdebugMessage("Schur Matrix in DSDP was indefinite but should have been positive semidefinite!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_MAX_IT:
         SCIPdebugMessage("DSDP reached maximum number of iterations!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_NUMERICAL_ERROR:
         SCIPdebugMessage("A numerical error occured in DSDP!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_UPPERBOUND:
         SCIPdebugMessage("Dual objective value in DSDP reached upper bound.\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case DSDP_USER_TERMINATION:
         SCIPdebugMessage("DSDP didn't stop solving, did you?\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         break;

      case CONTINUE_ITERATING:
         SCIPdebugMessage("DSDP wants to continue iterating but somehow was stopped!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         SCIPABORT();
         return SCIP_ERROR;
         break;

      default:
         SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
         SCIPABORT();
         return SCIP_ERROR;
         break;
   }
#endif

   return SCIP_OKAY;
}
/**@} */




/*
 * Solution Information Methods
 */

/**@name Solution Information Methods */
/**@{ */

/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType* pdfeasible;

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));

   switch ( *pdfeasible)
   {
      case DSDP_PDFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      case DSDP_UNBOUNDED:
         *primalfeasible = FALSE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      case DSDP_INFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = FALSE;
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         break;

      default: /* should only include DSDP_PDUNKNOWN */
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
         SCIPABORT();
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to be primal unbounded */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;
   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return SCIP_ERROR;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal infeasible */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPSolutionType* pdfeasible;
   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return SCIP_ERROR;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal feasible */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return !SCIPsdpiIsPrimalInfeasible(sdpi);
}

/** returns TRUE iff SDP is proven to be dual unbounded */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return SCIPsdpiIsPrimalInfeasible(sdpi);
}

/** returns TRUE iff SDP is proven to be dual infeasible */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return SCIPsdpiIsPrimalUnbounded(sdpi);
}

/** returns TRUE iff SDP is proven to be dual feasible */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return !SCIPsdpiIsPrimalUnbounded(sdpi);
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;
   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_CONVERGED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;
   BMSallocBlockMemory(sdpi->blkmem, &reason);

   SCIP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_UPPERBOUND)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;
   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_MAX_IT)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return SCIP_ERROR;
}

/** returns the internal solution status of the solver */
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiIgnoreInstability(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, objval));
   *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int                   dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL */
   )
{
   if ( objval != NULL )
   {
      DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, objval));
      *objval *= -1; /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   }
   if (dualsollength > 0)
   {
      assert(dualsol != NULL);
      DSDP_CALL(DSDPGetY(sdpi->dsdp, dualsol, dualsollength)); /*last entry needs to be the number of variables */
   }
   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPsdpiGetRealSolQuality(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */




/*
 * SDPi State Methods
 */

/**@name SDPi State Methods */
/**@{ */

/** stores SDPi state (like solve status since last data manipulation) into sdpistate object */
SCIP_RETCODE SCIPsdpiGetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information (like solve status since last data manipulation) */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** loads SDPi state into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetState()
 */
SCIP_RETCODE SCIPsdpiSetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE*       sdpistate           /**< SDPi state information */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** clears current SDPi state (like basis information) of the solver */
SCIP_RETCODE SCIPsdpiClearState(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** frees SDPi state information */
SCIP_RETCODE SCIPsdpiFreeState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDPi state information */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** reads SDPi state from a file */
SCIP_RETCODE SCIPsdpiReadState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** writes SDPi state to a file */
SCIP_RETCODE SCIPsdpiWriteState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */


/*
 * Parameter Methods
 */

/**@name Parameter Methods */
/**@{ */

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** sets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */




/*
 * Numerical Methods
 */

/**@name Numerical Methods */
/**@{ */

/** returns value treated as infinity in the SDP solver */
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   )
{
   return IS_POSINF(val) || IS_NEGINF(val);
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */
