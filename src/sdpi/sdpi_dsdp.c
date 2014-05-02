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
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */
#include "blockmemshell/memory.h"            /* for memory allocation */

#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */
#include "dsdpmem.h"                         /* for DSDPCALLOC2, DSDPFREE */


/** calls a DSDP-Function and transforms the return-code to a SCIP_ERROR if needed */
#define DSDP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_ERROR;                                                               \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/** same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced */
#define DSDP_CALLM(x)   do                                                                                     \
                       {                                                                                      \
                          int _dsdperrorcode_;                                                                \
                          if ( (_dsdperrorcode_ = (x)) != 0 )                                                 \
                          {                                                                                   \
                             SCIPerrorMessage("DSDP-Error <%d> in function call\n", _dsdperrorcode_);         \
                             SCIPABORT();                                                                     \
                             return SCIP_NOMEMORY;                                                            \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

/** this will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem */
#define CHECK_IF_SOLVED(sdpi)  do                                                                             \
                        {                                                                                     \
                           if (!(sdpi->solved))                                                               \
                           {                                                                                  \
                              SCIPerrorMessage("Tried to access solution information ahead of solving! \n");  \
                              SCIPABORT();                                                                    \
                              return SCIP_ERROR;                                                              \
                           }                                                                                  \
                        }                                                                                     \
                        while( FALSE )

/** data for SDP interface */
struct SCIP_SDPi
{
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP for handling variable bounds */
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
   int*                  sdpconstbegblock;   /**< start index of each block in sdpconstval-array, indices starting at 0 */
   int*                  sdpconstrowind;     /**< row-index for each entry in sdpconstval-array */
   int*                  sdpconstcolind;     /**< column-index for each entry in sdpconstval-array */
   SCIP_Real*            sdpconstval;        /**< values of entries of constant matrices in SDP-Block */
   int                   sdpnnonz;           /**< number of nonzero elements in the SDP-constraint matrices */
   int*                  sdpbegvarblock;     /**< entry j*nvars + i is the start index of matrix \f A_i^j \f in sdpval,
                                              *   particularly entry j * nvars gives the starting point of block j */
   int*                  sdprowind;          /**< row-index for each entry in sdpval-array */
   int*                  sdpcolind;          /**< column-index for each entry in sdpval-array */
   SCIP_Real*            sdpval;             /**< values of SDP-constraint matrix entries */
   int                   nlpcons;            /**< number of LP-constraints */
   SCIP_Real*            lprhs;              /**< right hand sides of LP rows */
   int                   lpnnonz;            /**< number of nonzero elements in the LP-constraint matrix */
   int*                  lprowind;           /**< row-index for each entry in lpval-array (going from 1 to nlpcons) */
   int*                  lpcolind;           /**< column-index for each entry in lpval-array (going from 1 to nvars) */
   SCIP_Real*            lpval;               /**< values of LP-constraint matrix entries */
   SCIP_Bool             solved;             /**< was the SDP solved since the problem was last changed */
};

static int nextsdpid     =  1;               /**< used to give ids to the generated sdps for debugging messages */
static double epsilon    = 1e-6;             /**< this is used for checking if primal and dual objective are equal */
static double feastol    = 1e-4;             /**< this is used for checking if a solution is feasible */

/*
 * Local Functions
 */

/** for given row and column (i,j) (indices going from 1 to n) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( i >= 1 );
   assert( j >= 1 );
   assert( i >= j );

   return i*(i-1)/2 + j - 1;
}

/**
 * For given row and column (i,j) (indices going from 1 to n) checks if i >= j, so that i and j give a position in the lower
 * triangular part, otherwise i and j will be switched. This function will be called whenever a position in a symmetric matrix
 * is given, to prevent problems if position (i,j) is given but later (j,i) should be changed.
 */
static void checkIfLowerTriang(
  int*                   i,                  /**< row index */
  int*                   j                   /**< column index */
  )
{
   if ( *i < *j )
   {
      int temp;
      temp = *i;
      *i = *j;
      *j = temp;
   }
}

/** This function checks feasibility (currently only LP-inequalities and only dual solution) and will be called if DSDP returns "primal-dual-unknown" */
static SCIP_Bool checkFeasibility(
   SCIP_SDPI*            sdpi               /**< pointer to an SDP interface structure */
   )
{
   int i;
   int ind;
   SCIP_Real rhs;
   SCIP_Real lhs;
   int nnonz;
   int* rowind;
   int* colind;
   SCIP_Real* val;
   SCIP_Real* sol;

   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemoryArray(sdpi->blkmem, &rowind, sdpi->nvars);
   BMSallocBlockMemoryArray(sdpi->blkmem, &colind, sdpi->nvars);
   BMSallocBlockMemoryArray(sdpi->blkmem, &val, sdpi->nvars);
   BMSallocBlockMemoryArray(sdpi->blkmem, &sol, sdpi->nvars);

   DSDP_CALL(DSDPGetY(sdpi->dsdp, sol, sdpi->nvars)); /* get the optimal solution */

   for (i = 1; i <= sdpi->nlpcons; i++)
   {
      SCIPsdpiGetLPRows(sdpi, i, i, &rhs, sdpi->nvars, &nnonz, rowind, colind, val);   /* get the next LPRow */

      lhs = 0; /* reset the left hand side */

      for (ind = 0; ind < nnonz; ind++)
      {
         assert ( rowind != NULL );
         assert ( rowind[ind] == i );
         lhs = lhs + val[ind] * sol[colind[ind]]; /* multiply the LP-coefficient with the value of the corresponding variable and
                                                   * summarize these for the left hand side value */
      }

      if (lhs + feastol < rhs)   /* this LP-inequality is violated */
      {
         BMSfreeBlockMemoryArray(sdpi->blkmem, &rowind, sdpi->nvars);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &colind, sdpi->nvars);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &val, sdpi->nvars);
         BMSfreeBlockMemoryArray(sdpi->blkmem, &sol, sdpi->nvars);
         return FALSE;
      }
   }

   BMSfreeBlockMemoryArray(sdpi->blkmem, &rowind, sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &colind, sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &val, sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &sol, sdpi->nvars);

   return TRUE; /* all tests were passed, so this solution has to be feasible */

   /* TODO: also check bounds and SDP-block and primal solution (possibly splitting this into checkDualFeasibility and checkPrimalFeasibility) */
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
   assert ( sdpi != NULL );
   assert ( *sdpi != NULL );

   if (((*sdpi)->dsdp) != NULL)
   {
   DSDP_CALL(DSDPDestroy((*sdpi)->dsdp));
   }

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
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks */
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
   int col;
   int row;
   SCIPdebugMessage("Calling SCIPsdpiLoadSDP (%d)\n",sdpi->sdpid);

   /* copy all inputs into the corresponding sdpi-parameters to later put them into DSDP prior to solving when the final
    * number of blocks and variables are known
    */

   sdpi->nvars = nvars;

   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), nvars);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), nvars);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), nvars);
   for (i = 0; i < nvars; i++)
   {
      (sdpi->obj)[i] = obj[i];
      (sdpi->lb)[i] = lb[i];
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
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpconstnnonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpconstnnonz);
   for (i = 0; i < sdpconstnnonz; i++)
   {
      row = sdpconstrowind[i];
      col = sdpconstcolind[i];
      checkIfLowerTriang(&row, &col);
      (sdpi->sdpconstrowind)[i] = row;
      (sdpi->sdpconstcolind)[i] = col;
      (sdpi->sdpconstval)[i] = sdpconstval[i];
   }

   sdpi->sdpnnonz = sdpnnonz;

   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), nvars * nsdpblocks);
   for (i = 0; i < nvars * nsdpblocks; i++)
   {
      (sdpi->sdpbegvarblock)[i] = sdpbegvarblock[i];
   }

   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpnnonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpnnonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpnnonz);
   for (i = 0; i < sdpnnonz; i++)
   {
      row = sdprowind[i];
      col = sdpcolind[i];
      checkIfLowerTriang(&row, &col);
      (sdpi->sdprowind)[i] = row;
      (sdpi->sdpcolind)[i] = col;
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
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), lpnnonz);
   BMSallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), lpnnonz);
   for (i = 0; i < lpnnonz; i++)
   {
      (sdpi->lprowind)[i] = lprowind[i];
      (sdpi->lpcolind)[i] = lpcolind[i];
      (sdpi->lpval)[i] = lpval[i];
   }

   sdpi->solved = FALSE;

   return SCIP_OKAY;
}

/** adds another SDP-Block to the problem
 *
 *  @note as \f A_i^j \f is symmetric, only the lower triangular part of it must be specified
 */
SCIP_RETCODE SCIPsdpiAddSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   blocksize,          /**< dimension of the matrices \f A_i^j \f */
   int                   constnnonz,         /**< sum of non-zeroes in the lower triagonal parts of \f A_0^j \f */
   const int*            constrowind,        /**< the row-indices of the non-zero entries of \f A_i^0 \f */
   const int*            constcolind,        /**< the column-indices of the non-zero entries of \f A_i^0 \f */
   const SCIP_Real*      constval,           /**< the values of \f A_i^0 \f as specified by constbegrow and constcolind */
   int                   nnonz,              /**< sum of non-zeroes in the lower triagonal parts of the \f A_i^j \f */
   const int*            begvar,             /**< start index of the matrix \f A_i^j \f for each i */
   const int*            rowind,             /**< the row-indices of the non-zero entries of \f A_i^j \f */
   const int*            colind,             /**< the column-indices of the non-zero entries of \f A_i^j \f */
   const SCIP_Real*      val                 /**< the values of of \f A_i^j \f as specified by begvar, rowind and colind */
   )
{
   int i;
   int row;
   int col;
   SCIPdebugMessage("Adding a block to SDP %d\n",nextsdpid);

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks + 1);
   (sdpi->sdpblocksizes)[sdpi->nsdpblocks] = blocksize; /* new SDP-Block will be added as the last block of the new SDP */


   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks + 1);
   (sdpi->sdpconstbegblock)[sdpi->nsdpblocks] = sdpi->sdpconstnnonz; /* new SDP-Block starts after all the old ones in the arrays */

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + constnnonz);
   for (i = 0; i < constnnonz; i++)
   {
      assert ( constrowind[i] <= blocksize );
      assert ( constcolind[i] <= blocksize );
      row = constrowind[i];
      col = constcolind[i];
      checkIfLowerTriang(&row, &col);
      (sdpi->sdpconstrowind)[sdpi->sdpconstnnonz + i] = row;
      (sdpi->sdpconstcolind)[sdpi->sdpconstnnonz + i] = col;
      (sdpi->sdpconstval)[sdpi->sdpconstnnonz + i] = constval[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nsdpblocks * sdpi->nvars, (sdpi->nsdpblocks + 1) * sdpi->nvars);
   for (i = 0; i < sdpi->nvars; i++)
   {
      /* new SDP-Block starts after all the old ones in the arrays */
      (sdpi->sdpbegvarblock)[(sdpi->nsdpblocks * sdpi->nvars) + i] = sdpi->sdpnnonz + begvar[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + nnonz);
   for (i = 0; i < nnonz; i++)
   {
      assert ( rowind[i] <= blocksize );
      assert ( colind[i] <= blocksize );
      row = rowind[i];
      col = colind[i];
      checkIfLowerTriang(&row, &col);
      (sdpi->sdprowind)[sdpi->sdpnnonz + i] = row;
      (sdpi->sdpcolind)[sdpi->sdpnnonz + i] = col;
      (sdpi->sdpval)[sdpi->sdpnnonz + i] = val[i];
   }

   sdpi->nsdpblocks++;
   sdpi->sdpconstnnonz = sdpi->sdpconstnnonz + constnnonz;
   sdpi->sdpnnonz = sdpi->sdpnnonz + nnonz;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** deletes a SDP-Block */
SCIP_RETCODE SCIPsdpiDelSDPBlock(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block              /**< index of the SDP-Block (indexing starting at 1) that should be deleted */
   )
{
   int movingblock;
   int var;
   int i;
   int deletedconstnonz;
   int newsdpconstnnonz;
   int deletednonz;
   int newsdpnnonz;

   SCIPdebugMessage("Deleting block %d from SDP %d\n",block, nextsdpid);

   assert ( block > 0);
   assert ( block <= sdpi->nsdpblocks );

   if (block == sdpi->nsdpblocks) /* the block can simply be deleted */
   {
      deletedconstnonz = sdpi->sdpconstnnonz - sdpi->sdpconstbegblock[block - 1]; /* sdpconstbegblock[block] gives the first index belonging to the deleted block, all thereafter need to be
                                                                                   * deleted, the indexshift is because DSDP starts counting the blocks at 0 */
      newsdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnonz;
      deletednonz = sdpi->sdpnnonz - sdpi->sdpbegvarblock[(block-1) * sdpi->nvars]; /* sdpbegvarblock[block*nvars] gives the first index of the deleted block */
      newsdpnnonz = sdpi->sdpnnonz - deletednonz;

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, sdpi->nvars * (sdpi->nsdpblocks - 1));
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, newsdpnnonz);

      sdpi->nsdpblocks--;
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnonz;
      sdpi->sdpnnonz = sdpi->sdpnnonz - deletednonz;
   }
   else /* all blocks after the deleted block need to be shifted in the arrays */
   {
      /* compute the new numbers of nonzeroes, these need to be computed before the begvar-arrays are updated, but the old values are still needed for iterating */
      deletedconstnonz =  sdpi->sdpconstbegblock[block] - sdpi->sdpconstbegblock[block-1]; /* starting index of the next block minus starting index of the deleted block gives the number of
                                                                                         * nonzeroes of the deleted block, which is then substracted from the old value, indices are
                                                                                         * shifted from block + 1 and block to block and block-1 because in DSDP indexing starts at 0 */
      newsdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnonz;
      deletednonz = sdpi->sdpbegvarblock[block * sdpi->nvars] - sdpi->sdpbegvarblock[(block-1) * sdpi->nvars]; /* same as above, but because of the structure of the
                                                                                                             * sdpbegvarblock-arrays block*nvars gives the first index of
                                                                                                             * that block, again this is shifted because of DSDPs indexing */
      newsdpnnonz = sdpi->sdpnnonz - deletednonz;


      /* all later blocks (because DSDP starts counting the blocks at 0 this starts at block) need to be moved to the left in the arrays to fill the spot of the deleted block */
      for (movingblock = block; movingblock < sdpi->nsdpblocks; movingblock++)
      {
         sdpi->sdpblocksizes[movingblock - 1] = sdpi->sdpblocksizes[movingblock];
         sdpi->sdpconstbegblock[movingblock - 1] = sdpi->sdpconstbegblock[movingblock] - deletedconstnonz;
         for (var = 0; var < sdpi->nvars; var++)
         {
            sdpi->sdpbegvarblock[movingblock * sdpi->nvars + var - sdpi->nvars] =  sdpi->sdpbegvarblock[movingblock * sdpi->nvars + var] - deletednonz; /* these are shifted nvars spaces
                                                                                                                                        * to the left, because there are nvars entries
                                                                                                                                        * in sdpbegvarblock belonging to the
                                                                                                                                        * deleted block */
         }
      }
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks, sdpi->nsdpblocks - 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, sdpi->nvars * (sdpi->nsdpblocks - 1));

      /* shift all nonzeroes to the left by a number of spots equal to the number of nonzeroes in the deleted block */
      for (i = sdpi->sdpconstbegblock[block]; i<sdpi->sdpconstnnonz; i++)
      {
         sdpi->sdpconstrowind[i - deletedconstnonz] = sdpi->sdpconstrowind[i];
         sdpi->sdpconstcolind[i - deletedconstnonz] = sdpi->sdpconstcolind[i];
         sdpi->sdpconstval[i - deletedconstnonz] = sdpi->sdpconstval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, newsdpconstnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, newsdpconstnnonz);

      for (i = sdpi->sdpbegvarblock[(block) * sdpi->nvars]; i < sdpi->sdpnnonz; i++)
      {
         sdpi->sdprowind[i - deletednonz] = sdpi->sdprowind[i];
         sdpi->sdpcolind[i - deletednonz] = sdpi->sdpcolind[i];
         sdpi->sdpval[i - deletednonz] = sdpi->sdpval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, newsdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, newsdpnnonz);

      sdpi->nsdpblocks--;
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz - deletedconstnonz;
      sdpi->sdpnnonz = sdpi->sdpnnonz - deletednonz;
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** adds additional variables to the SDP
 *
 *  @note arrays are not checked for duplicates, problems may appear if indices are added more than once
 */
SCIP_RETCODE SCIPsdpiAddVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to be added */
   const SCIP_Real*      obj,                /**< objective function values of new variables */
   const SCIP_Real*      lb,                 /**< lower bounds of new variables */
   const SCIP_Real*      ub,                 /**< upper bounds of new variables */
   int                   sdpnnonz,           /**< number of nonzero elements to be added to the SDP constraint matrices */
   const int*            sdpbegvarblock,     /**< start index of each new variable for each block in sdpval-array, or NULL if sdpnnonz == 0 */
   const int*            sdprowind,          /**< row indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const int*            sdpcolind,          /**< col indices of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   const SCIP_Real*      sdpval,             /**< values of SDP constraint matrix entries, or NULL if sdpnnonz == 0 */
   int                   lpnnonz,            /**< number of nonzero elements to be added to the LP constraint matrices */
   const int*            lprowind,           /**< row indices of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   const int*            lpcolind,           /**< col indices of LP constraint matrix entries, indexed from 1 to the number of added vars, or NULL if lpnnonz == 0 */
   const SCIP_Real*      lpval               /**< values of LP constraint matrix entries, or NULL if lpnnonz == 0 */
   )
{
   int i;
   int block;
   int toInsert;

   SCIPdebugMessage("Adding %d variables to SDP %d.\n", nvars, nextsdpid);

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars, sdpi->nvars + nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars, sdpi->nvars + nvars);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars, sdpi->nvars + nvars);

   for (i = 0; i < nvars; i++)
   {
      sdpi->obj[sdpi->nvars + i] = obj[i];
      sdpi->lb[sdpi->nvars + i] = lb[i];
      sdpi->ub[sdpi->nvars + i] = ub[i];
   }

   if (sdpnnonz > 0)
   {
      int row;
      int col;
      int* begblockold; /* these are the old starting indices of the blocks (only for the first variable in that block), with extra entry equal to spdnnonz */

      BMSallocBlockMemoryArray(sdpi->blkmem, &begblockold, sdpi->nsdpblocks + 1);

      begblockold[sdpi->nsdpblocks] = sdpi->sdpnnonz; /* often begblockold[block+1] is needed, so this extra entry removes some additional if-clauses */

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, (sdpi->nvars + nvars) * sdpi->nsdpblocks);

      for (block = sdpi->nsdpblocks - 1; block > -1; block--)
      {
         begblockold[block] = sdpi->sdpbegvarblock[block * sdpi->nvars];

         for (i = sdpi->nvars - 1; i > -1; i--)
         {
            /* all the old entries have to be shifted to the right to get the needed space for inserting the new variables for all the blocks before it, also the
             * number of nonzeroes that are added for the new variables for the earlier blocks have to be added, for not overwriting needed entries this iteration
             *  has to go from right to left
             */
            sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) +i ] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] + sdpbegvarblock[block*nvars];
         }

         for (i = 0; i < nvars; i++)
         {
            /* insert the new values, add to them the start-index of the next block in the original problem (as all nonzeroes in earlier and this block of the
             * original problem as well as the new ones in the earlier blocks have to be added before the new nonzeroes in this block), for the last block the
             * number of nonzeroes is used, as sdpbegvarblock[nblocks*nvars] doesn't exist */
               sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars + i] = sdpbegvarblock[block * nvars + i] + begblockold[block + 1];
         }
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + sdpnnonz);

      /* now insert the nonzero-entries at the right positions of the arrays */
      for (block=sdpi->nsdpblocks - 1; block > -1; block--)
      {
         for (i = begblockold[block + 1]; i > begblockold[block] - 1; i--)
         {
            /* all the entries belonging to block j have to be shifted sdpbegvarblock[block*nvars] entries to the left, as this many entries belonging to
             * earlier blocks and new variables have to be inserted in front of them
             */
            sdpi->sdprowind[i + sdpbegvarblock[block * nvars]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i + sdpbegvarblock[block * nvars]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i + sdpbegvarblock[block * nvars]] = sdpi->sdpval[i];
         }

         /* compute the number of new nonzeroes to be inserted into this block */
         if (block == sdpi->nsdpblocks - 1)
         {
            toInsert = sdpnnonz - sdpbegvarblock[block * nvars];
         }
         else
         {
            toInsert = sdpbegvarblock[(block+1) * nvars] - sdpbegvarblock[block * nvars];
         }
         for (i = 0; i < toInsert; i++)
         {
            /* insert the new values for this block, they are inserted where the first new variable for the specific block starts */
            assert ( sdprowind[sdpbegvarblock[block * nvars]+i] <= sdpi->sdpblocksizes[block] );
            assert ( sdpcolind[sdpbegvarblock[block * nvars]+i] <= sdpi->sdpblocksizes[block] ); /* the row and column indices shouldn't exceed blocksizes */
            row = sdprowind[sdpbegvarblock[block * nvars]+i];
            col = sdpcolind[sdpbegvarblock[block * nvars]+i];
            checkIfLowerTriang(&row, &col);
            sdpi->sdprowind[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = row;
            sdpi->sdpcolind[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = col;
            sdpi->sdpval[sdpi->sdpbegvarblock[block * (sdpi->nvars + nvars) + sdpi->nvars]+i] = sdpval[sdpbegvarblock[block * nvars]+i];
         }
      }

      BMSfreeBlockMemoryArray(sdpi->blkmem, &begblockold, sdpi->nsdpblocks);

      sdpi->sdpnnonz = sdpi->sdpnnonz + sdpnnonz;
   }

   if (lpnnonz > 0)
   {
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + lpnnonz);

      for (i = 0; i < lpnnonz; i++)
      {
         assert ( lprowind[i] <= sdpi->nlpcons ); /* only insert into existing LP constraints */
         sdpi->lprowind[sdpi->lpnnonz + i] = lprowind[i]; /* just add these at the end, they will be sorted before solving */
         sdpi->lpcolind[sdpi->lpnnonz + i] = lpcolind[i] + sdpi->nvars; /* the columns are added to the right of the old ones, so the column indices have to be shifted
                                                                         * by the number of old variables */
         sdpi->lpval[sdpi->lpnnonz + i] = lpval[i];
      }

      sdpi->lpnnonz = sdpi->lpnnonz + lpnnonz;
   }

   sdpi->nvars = sdpi->nvars + nvars;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** deletes all variables in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to be deleted, indexing starts at 1 */
   int                   lastvar             /**< last variable to be deleted, indexing starts at 1 */
   )
{
   int block;
   int i;
   int* deletedsdpnonz; /* entry i will give the number of deleted sdp nonzeroes in block 1 till i+1 (with indexing starting at 1) */
   int deletedvars;
   int firstvarlpind;
   int lastvarlpind;
   int deletedlpnonz;
   int lastindexforshifting;

   SCIPdebugMessage("Deleting vars %d to %d from SDP %d.\n", firstvar, lastvar, nextsdpid);
   assert ( firstvar > 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar <= sdpi->nvars );

   deletedvars = lastvar - firstvar + 1;

   BMSallocBlockMemoryArray(sdpi->blkmem, &deletedsdpnonz, sdpi->nsdpblocks);
   deletedsdpnonz[0] = sdpi->sdpbegvarblock[lastvar] - sdpi->sdpbegvarblock[firstvar - 1]; /* indexshift because these are indexed starting at 0,
                                                                                             * begvarblock[lastvar] gives the first index of the first
                                                                                             * non-deleted block */
   for (block = 1; block < sdpi->nsdpblocks; block++)
   {
      if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars)
         {
         deletedsdpnonz[block] = deletedsdpnonz[block-1] + sdpi->sdpnnonz - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar -1];
         }
      else
      {
         deletedsdpnonz[block] = deletedsdpnonz[block-1] + sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar]
                                                                             - sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar -1];
      }
   }

   for (i=lastvar; i < sdpi->nvars; i++) /* index shift, because the arrays start counting at 0 */
   {
      sdpi->obj[i - deletedvars] = sdpi->obj[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars, sdpi->nvars - deletedvars);

   for (i=lastvar; i < sdpi->nvars; i++) /* index shift, because the arrays start counting at 0 */
   {
      sdpi->lb[i - deletedvars] = sdpi->lb[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars, sdpi->nvars - deletedvars);

   for (i=lastvar; i < sdpi->nvars; i++) /* index shift, because the arrays start counting at 0 */
   {
      sdpi->ub[i - deletedvars] = sdpi->ub[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars, sdpi->nvars - deletedvars);

   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      if (block > 0)
      {
         /* first look at all nonzeroes in this given block before the deleted vars (index shift because the arreas count from 0), for the first block there's
          * nothing to do, as no entries before those were deleted */
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars]; i < sdpi->sdpbegvarblock[block * sdpi->nvars + firstvar - 1]; i++)
         {
            sdpi->sdprowind[i - deletedsdpnonz[block - 1]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i - deletedsdpnonz[block - 1]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i - deletedsdpnonz[block - 1]] = sdpi->sdpval[i];
         }
      }
      /* then look at all nonzeroes in this given block after the deleted vars, here they are shifted by deletsdpnnonz[block] instead of block - 1 */
      if (lastvar < sdpi->nvars - 1) /* if the deleted var is the last one then there aren't any left after it in the same block, this extra if is needed,
                                      * because otherwise the start index of the for loop would be outside the bounds of sdpbegvarblock for the last block*/
      {
         if (block == sdpi->nsdpblocks - 1)
         {
            lastindexforshifting = sdpi->sdpnnonz;
         }
         else
         {
            lastindexforshifting = sdpi->sdpbegvarblock[(block+1) * sdpi->nvars];
         }
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar]; i < lastindexforshifting; i++)
         {
            sdpi->sdprowind[i - deletedsdpnonz[block]] = sdpi->sdprowind[i];
            sdpi->sdpcolind[i - deletedsdpnonz[block]] = sdpi->sdpcolind[i];
            sdpi->sdpval[i - deletedsdpnonz[block]] = sdpi->sdpval[i];
         }
      }
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1]);

   /* sdpbegvarblock should be updated last to still be able to find the variables which should be moved or deleted */
   for (block = 0; block < sdpi->nsdpblocks; block++)
   {
      for (i = 0; i < sdpi->nvars; i++)
      {
         if (i < firstvar - 1 && block > 0) /* index shift, because counting starts at 0, for block 0 nothing needs to be done prior to the first deleted var */
         {
            /* the entry will be moved to the corresponding position with the decreased number of variables and the number of deleted nonzeroes in earlier blocks
             * will be substracted */
            sdpi->sdpbegvarblock[block * (sdpi->nvars - deletedvars) + i] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] - deletedsdpnonz[block - 1];
         }
         else if (i > lastvar - 1) /* index shift, because counting starts at 0 */
         {
            sdpi->sdpbegvarblock[block * (sdpi->nvars - deletedvars) + i - deletedvars] = sdpi->sdpbegvarblock[block * sdpi->nvars + i] -
                  deletedsdpnonz[block]; /* this time it is moved even further left because in this block the variables were also deleted, and the entry also
                                           * also has to be decreased further because now also the deleted nonzeroes from this block must be substracted */
         }
      }
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks, (sdpi->nvars - deletedvars) * sdpi->nsdpblocks);

   sdpi->sdpnnonz = sdpi->sdpnnonz - deletedsdpnonz[sdpi->nsdpblocks - 1];

   BMSfreeBlockMemoryArray(sdpi->blkmem, &(deletedsdpnonz), sdpi->nsdpblocks);

   /* now the LP arrays are cleared of the deleted vars, for this they will have to be sorted first to have those belonging to the deleted vars together */
   SCIPsortIntIntReal(sdpi->lpcolind, sdpi->lprowind, sdpi->lpval, sdpi->lpnnonz); /* now the arrays should be sorted by nondecreasing column indices */

   /*iterate over the lpcolind array to find the first index belonging to a deleted var */
   firstvarlpind = -1; /* if this stays at -1 the deleted variables weren't part of the LP block */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lpcolind[i] >= firstvar && sdpi->lpcolind[i] <= lastvar)
      {
         firstvarlpind = i;
         lastvarlpind = i;
         i++;
         break;
      }
   }

   if (firstvarlpind > -1) /* if this is still -1 nothing has to be done for the LP part */
   {
   /* now find the last occurence of a deleted variable (as these are sorted all in between also belong to deleted vars and will be removed) */
      while (i < sdpi->lpnnonz && sdpi->lpcolind[i] <= lastvar)
      {
         lastvarlpind++;
         i++;
      }

      deletedlpnonz = lastvarlpind - firstvarlpind + 1;

      /* finally shift all LP-array-entries after the deleted variables */
      for (i = lastvarlpind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcolind[i - deletedlpnonz] = sdpi->lpcolind[i] - deletedvars; /* the column indices have to be decreased by the number of vars deleted
          * before that var */
         sdpi->lprowind[i - deletedlpnonz] = sdpi->lprowind[i];
         sdpi->lpval[i - deletedlpnonz] = sdpi->lpval[i];
      }

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletedlpnonz);

      sdpi->lpnnonz = sdpi->lpnnonz - deletedlpnonz;
   }
   sdpi->nvars = sdpi->nvars - deletedvars;

   sdpi->solved = FALSE;

   /* at this point there could be checked if any SDP-blocks or LP-rows have become empty (no variables left), but this isn't done,
    * because then the indices of alle blocks/rows behind it would change, possibly creating problems if the user wanted to insert
    * variables into them (or even the deleted block/row) afterwards
    */

   return SCIP_OKAY;
}

/** deletes variables from SCIP_SDPI; the new position of a variable must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelVarset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of variables
                                              *   input:  1 if variable should be deleted, 0 if not
                                              *   output: new position of variable, -1 if variable was deleted */
   )
{
   int i;
   int deletedvars;
   int oldnvars;
   SCIPdebugMessage("Calling SCIPsdpiDelColset for sdpi %d.\n", sdpi->sdpid);

   deletedvars = 0;
   oldnvars = sdpi->nvars;

   for (i=0; i < oldnvars; i++)
   {
      if (dstat[i] == 1)
      {
         SCIPsdpiDelVars(sdpi, i + 1 - deletedvars, i + 1 - deletedvars); /* delete this variable, the index-shift is needed as SCIPsdpiDelVars
                                                                           * asks for an index between 1 and nvars, but earliers deletions are
                                                                           * already applied to the problem, so the index also has to be lowered
                                                                           * by deletedvars */
         dstat[i] = -1;
         deletedvars++;
      }
      else
      {
         dstat[i] = i - deletedvars;
      }
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** adds rows to the LP-Block
 *
 *  @note arrays are not checked for duplicates, problems may appear if indices are added more than once
 */
SCIP_RETCODE SCIPsdpiAddLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   int                   nnonz,              /**< number of nonzero elements to be added to the LP constraint matrix */
   const int*            rowind,             /**< row indices of constraint matrix entries, going from 1 to nrows, these will be changed to nlpcons + i */
   const int*            colind,             /**< column indices of constraint matrix entries, going from 1 to nvars */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{
   int i;
   SCIPdebugMessage("Adding %d LP-Constraints to SDP %d.\n", nrows, sdpi->sdpid);

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons + nrows);
   for (i=0; i < nrows; i++)
   {
      sdpi->lprhs[sdpi->nlpcons + i] = rhs[i];
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + nnonz);

   for (i=0; i < nnonz; i++)
   {
      sdpi->lprowind[sdpi->lpnnonz + i] = rowind[i] + sdpi->nlpcons; /* the new rows are added at the end, so the row indices are increased by the old
                                                                               * number of LP-constraints */
      assert(colind[i] <= sdpi->nvars); /* only existing vars should be added to the LP-constraints */
      sdpi->lpcolind[sdpi->lpnnonz + i] = colind[i];
      sdpi->lpval[sdpi->lpnnonz + i] = val[i];
   }

   sdpi->nlpcons = sdpi->nlpcons + nrows;
   sdpi->lpnnonz = sdpi->lpnnonz + nnonz;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** deletes all rows in the given range from LP-Block */
SCIP_RETCODE SCIPsdpiDelLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted, index between 1 and nlpcons */
   int                   lastrow             /**< last row to be deleted, index between 1 and nlpcons */
   )
{
   int i;
   int deletedrows;
   int firstrowind;
   int lastrowind;
   int deletednonz;
   SCIPdebugMessage("Deleting rows %d to %d from SDP %d.\n", firstrow, lastrow, sdpi->sdpid);

   assert ( firstrow > 0 );
   assert ( firstrow <= lastrow );
   assert ( lastrow <= sdpi->nlpcons );

   deletedrows = lastrow - firstrow + 1;
   deletednonz = 0;

   /* first delete the right-hand-sides */
   for (i = lastrow; i < sdpi->nlpcons; i++) /* shift all rhs after the deleted rows (indexshift because of starting the arrays at 0) */
   {
      sdpi->lprhs[i - deletedrows] = sdpi->lprhs[i];
   }
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons, sdpi->nlpcons - deletedrows);

   /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
   SCIPsortIntIntReal(sdpi->lprowind, sdpi->lpcolind, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

   firstrowind = -1;
   /*iterate over the lprowind array to find the first index belonging to a row that should be deleted */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprowind[i] >= firstrow && sdpi->lprowind[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
      {
         firstrowind = i;
         lastrowind = i;
         i++;
         break;
      }
   }

   if (firstrowind > -1) /* if this is still 0 there are no nonzeroes for the given rows */
   {
      /* now find the last occurence of one of the rows (as these are sorted all in between also belong to deleted rows and will be removed) */
      while (i < sdpi->lpnnonz && sdpi->lprowind[i] <= lastrow)
      {
         lastrowind++;
         i++;
      }
      deletednonz = lastrowind - firstrowind + 1;

      /* finally shift all LP-array-entries after the deleted rows */
      for (i = lastrowind + 1; i < sdpi->lpnnonz; i++)
      {
         sdpi->lpcolind[i - deletednonz] = sdpi->lpcolind[i];
         sdpi->lprowind[i - deletednonz] = sdpi->lprowind[i] - deletedrows; /* all rowindices after the deleted ones have to be lowered to still have ongoing
                                                                             * indices from 1 to nlpcons */
         sdpi->lpval[i - deletednonz] = sdpi->lpval[i];
      }
   }

   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz - deletednonz);
   sdpi->nlpcons = sdpi->nlpcons - deletedrows;
   sdpi->lpnnonz = sdpi->lpnnonz - deletednonz;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** deletes LP rows from SCIP_SDPI; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelLPRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of LP rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int i;
   int oldnlpcons;
   int deletedrows;

   SCIPdebugMessage("Calling SCIPsdpiDelLPRowset for SDP %d.\n", sdpi->sdpid);

   oldnlpcons = sdpi->nlpcons;
   deletedrows = 0;

   for (i = 0; i < oldnlpcons; i++)
   {
      if (dstat[i] == 1)
      {
         SCIPsdpiDelLPRows(sdpi, i + 1 - deletedrows, i + 1 - deletedrows); /* delete this row, index shift because of indexing starting
                                                                             * at 1 from DelLPRows, and - deletedrows, because in this
                                                                             * problem the earlier rows have already been deleted */
         dstat[i] = -1;
         deletedrows++;
      }
      else
         dstat[i] = i - deletedrows;
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert ( sdpi != NULL );
   SCIPdebugMessage("Calling SCIPsdpiClear for SDPI (%d) \n", sdpi->sdpid);

   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->obj), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lb), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->ub), sdpi->nvars);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpblocksizes), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstbegblock), sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpbegvarblock), sdpi->nvars * sdpi->nsdpblocks);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprhs), sdpi->nlpcons);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz);
   BMSfreeBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz);

   sdpi->nvars = 0;
   sdpi->nsdpblocks = 0;
   sdpi->sdpconstnnonz = 0;
   sdpi->sdpnnonz = 0;
   sdpi->nlpcons = 0;
   sdpi->lpnnonz = 0;

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes lower and upper bounds of variables */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of variables to change bounds for */
   const int*            ind,                /**< variables indices (between 1 and nvars) */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   int i;
   SCIPdebugMessage("Changing %d variable bounds in SDP %d\n", nvars, sdpi->sdpid);

   for (i = 0; i < nvars; i++)
   {
      assert ( ind[i] > 0 );
      assert ( ind[i] <= sdpi->nvars );
      sdpi->lb[ind[i] - 1] = lb[i]; /* index shift because the arrays are indexed starting at 0 */
      sdpi->ub[ind[i] - 1] = ub[i];
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes right hand sides of LP rows */
SCIP_RETCODE SCIPsdpiChgLPRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of LP rows to change right hand sides for */
   const int*            ind,                /**< row indices between 1 and nlpcons */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   int i;

   SCIPdebugMessage("Changing %d right hand sides of SDP %d\n", nrows, sdpi->sdpid);

   for (i = 0; i < nrows; i++)
   {
      assert ( ind[i] > 0 );
      assert ( ind[i] <= sdpi->nlpcons );
      sdpi->lprhs[ind[i] - 1] = rhs[i]; /* index shift because the arrays are indexed starting at 0 */
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes a single coefficient in LP constraint matrix */
SCIP_RETCODE SCIPsdpiChgLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of LP-coefficient to change, index between 1 and nlpcons */
   int                   col,                /**< column number of LP-coefficient to change, index between 1 and nvars */
   SCIP_Real             newval              /**< new value of LP-coefficient */
   )
{
   int i;
   SCIP_Bool found;

   SCIPdebugMessage("Changed the LP Coefficient in row %d and colum %d of SDP %d.\n", row, col, sdpi->sdpid);

   assert ( row > 0 );
   assert ( row <= sdpi->nlpcons );
   assert ( col > 0 );
   assert ( col <= sdpi->nvars );

   /* check if that entry already exists */
   found = FALSE;
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lprowind[i] == row && sdpi->lpcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->lpval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("An LP Coefficient in row %d and colum %d of SDP %d didn't exist so far or was zero, it is now added.\n", row, col, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lprowind), sdpi->lpnnonz, sdpi->lpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpcolind), sdpi->lpnnonz, sdpi->lpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->lpval), sdpi->lpnnonz, sdpi->lpnnonz + 1);

      sdpi->lprowind[sdpi->lpnnonz] = row;
      sdpi->lpcolind[sdpi->lpnnonz] = col;
      sdpi->lpval[sdpi->lpnnonz] = newval;
      sdpi->lpnnonz = sdpi->lpnnonz + 1;
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes objective values of variables in the SDP */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of variables to change objective value for */
   int*                  ind,                /**< variable indices to change objective value for, indexing starts at 1 */
   SCIP_Real*            obj                 /**< new objective values for variables */
   )
{
   int i;

   SCIPdebugMessage("Changing %d objective values of SDP %d.\n", ncols, sdpi->sdpid);

   for (i=0; i < ncols; i++)
   {
      assert ( ind[i] > 0 );
      assert ( ind[i] <= sdpi->nvars );
      sdpi->obj[ind[i] - 1] = obj[ind[i]]; /* index shift because indexing in arrays starts at 0 */
   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes a single coefficient in constant matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPConstCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index between 1 and nsdpblocks */
   int                   rowind,             /**< row index between 1 and corresponding sdpblocksize */
   int                   colind,             /**< column index between 1 and corresponding sdpblocksize*/
   const SCIP_Real      newval               /**< new value of given entry of constant matrix in given SDP-Block */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;
   SCIP_Bool found;

   SCIPdebugMessage("Changing a SDP coefficient in block %d of SDP %d.\n", block, sdpi->sdpid);

   assert ( block > 0 );
   assert ( block <= sdpi->nsdpblocks );
   assert ( rowind > 0 );
   assert ( rowind <= sdpi->sdpblocksizes[block - 1] ); /* index shift because arrays start at 0 */
   assert ( colind > 0 );
   assert ( colind <= sdpi->sdpblocksizes[block - 1] );

   row = rowind;
   col = colind;
   checkIfLowerTriang(&row, &col); /* make sure that this is a lower triangular position, otherwise it could happen that an upper
                                  * triangular position is added if the same entry is already filled in the lower triangular part */

   /* check if that entry already exists */
   found = FALSE;
   if (block == sdpi->nsdpblocks)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpconstbegblock[block]; /* index shift beause arrays start at 0 */
   }
   for (i = sdpi->sdpconstbegblock[block - 1]; i < lastiterationindex; i++)
   {
      if (sdpi->sdpconstrowind[i] == row && sdpi->sdpconstcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->sdpconstval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("A constant SDP Coefficient in row %d and colum %d of block %d of SDP %d didn't exist so far or was zero, it is now added.\n",
            rowind, colind, block, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstrowind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstcolind), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpconstval), sdpi->sdpconstnnonz, sdpi->sdpconstnnonz + 1);

      /* move all sdpconstnonzeroes of later blocks one space in the arrays to be able to insert this one at the right position */
      for (i = sdpi->sdpconstnnonz - 1; i >= sdpi->sdpconstbegblock[block]; i--)
      {
         sdpi->sdpconstrowind[i + 1] = sdpi->sdpconstrowind[i];
         sdpi->sdpconstcolind[i + 1] = sdpi->sdpconstcolind[i];
         sdpi->sdpconstval[i + 1] = sdpi->sdpconstval[i];
      }

      /* insert the new entries at the right position (namely what was originally the first position of the next block [again there's an index shift]) */
      sdpi->sdpconstrowind[sdpi->sdpconstbegblock[block]] = row;
      sdpi->sdpconstcolind[sdpi->sdpconstbegblock[block]] = col;
      sdpi->sdpconstval[sdpi->sdpconstbegblock[block]] = newval;

      /* update other information */
      for (i = block; i < sdpi->nsdpblocks; i++)
         {
         sdpi->sdpconstbegblock[i]++; /* all later blocks start one spot later in the arrays */
         }
      sdpi->sdpconstnnonz = sdpi->sdpconstnnonz + 1;

   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}

/** changes a single coefficient in a constraint matrix of given SDP-Block */
SCIP_RETCODE SCIPsdpiChgSDPCoeff(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index between 1 and nsdpblocks */
   int                   var,                /**< variable index between 1 and nvars */
   int                   rowind,             /**< row index between 1 and corresponding blocksize */
   int                   colind,             /**< column index between 1 and corresponding blocksize */
   const SCIP_Real      newval              /**< new value of given entry of the give constraint matrix in specified SDP-Block */
   )
{
   int i;
   int lastiterationindex;
   SCIP_Bool found;
   int row;
   int col;

   SCIPdebugMessage("Changing a Coefficient of Matrix A_%d^%d in SDP %d\n", var, block, sdpi->sdpid);

   assert ( block > 0 );
   assert ( block <= sdpi->nsdpblocks );
   assert ( var > 0 );
   assert ( var <= sdpi->nvars );
   assert ( rowind > 0 );
   assert ( rowind <= sdpi->sdpblocksizes[block - 1] ); /* index shift because arrays start at 0 */
   assert ( colind > 0 );
   assert ( colind <= sdpi->sdpblocksizes[block - 1] );

   row = rowind;
   col = colind;
   checkIfLowerTriang(&row, &col); /* make sure that this is a lower triangular position, otherwise it could happen that an upper
                                  * triangular position is added if the same entry is already filled in the lower triangular part */

   /* check if that entry already exists */
   found = FALSE;
   if (block == sdpi->nsdpblocks && var == sdpi->nvars)
   {
      lastiterationindex = sdpi->sdpnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpbegvarblock[(block-1) * sdpi->nvars + var]; /* index shift beause arrays start at 0 */
   }
   for (i = sdpi->sdpbegvarblock[(block - 1) * sdpi->nvars + (var - 1)]; i < lastiterationindex; i++)
   {
      if (sdpi->sdprowind[i] == row && sdpi->sdpcolind[i] == col)
      {
         found = TRUE;
         break;
      }
   }

   if (found)
   {
      sdpi->sdpval[i] = newval;
   }
   else
   {
      SCIPdebugMessage("A SDP Coefficient in row %d and colum %d of of Matrix A_%d^%d of SDP %d didn't exist so far or was zero, it is now added.\n",
            rowind, colind, var, block, sdpi->sdpid);

      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdprowind), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpcolind), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);
      BMSreallocBlockMemoryArray(sdpi->blkmem, &(sdpi->sdpval), sdpi->sdpnnonz, sdpi->sdpnnonz + 1);

      /* move all sdpnonzeroes of later blocks and vars one space in the arrays to be able to insert this one at the right position */
      for (i = sdpi->sdpnnonz - 1; i >= sdpi->sdpbegvarblock[(block - 1) * sdpi->nvars + var]; i--)
      {
         sdpi->sdprowind[i + 1] = sdpi->sdprowind[i];
         sdpi->sdpcolind[i + 1] = sdpi->sdpcolind[i];
         sdpi->sdpval[i + 1] = sdpi->sdpval[i];
      }

      /* insert the new entries at the right position (namely what was originally the first position of the next block [again there's an index shift]) */
      sdpi->sdprowind[sdpi->sdpbegvarblock[(block - 1)*sdpi->nvars + var]] = row;
      sdpi->sdpcolind[sdpi->sdpbegvarblock[(block - 1)*sdpi->nvars + var]] = col;
      sdpi->sdpval[sdpi->sdpbegvarblock[(block - 1)*sdpi->nvars + var]] = newval;

      /* update other information */
      for (i = (block - 1)*sdpi->nvars + var; i < sdpi->nsdpblocks * sdpi->nvars; i++)
         {
         sdpi->sdpbegvarblock[i]++; /* all later blocks start one spot later in the arrays */
         }
      sdpi->sdpnnonz = sdpi->sdpnnonz + 1;

   }

   sdpi->solved = FALSE;
   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */

/**@name Data Accessing Methods */
/**@{ */

/** gets the number of LP-rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nlprows             /**< pointer to store the number of rows */
   )
{
   *nlprows = sdpi->nlpcons;
   return SCIP_OKAY;
}

/** gets the number of SDP-Blocks in the SDP */
SCIP_RETCODE SCIPsdpiGetNSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpblocks          /**< pointer to store the number of rows */
   )
{
   *nsdpblocks = sdpi->nsdpblocks;
   return SCIP_OKAY;
}

/** gets the number of variables in the SDP */
SCIP_RETCODE SCIPsdpiGetNVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nvars               /**< pointer to store the number of variables */
   )
{
   *nvars = sdpi->nvars;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP constraint matrices */
SCIP_RETCODE SCIPsdpiGetSDPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the SDP constraint matrcies */
   )
{
   *nnonz = sdpi->sdpnnonz;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the constant matrices of the SDP-Blocks */
SCIP_RETCODE SCIPsdpiGetConstNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the constant matrices of the SDP-Blocks */
   )
{
   *nnonz = sdpi->sdpconstnnonz;
   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the LP Matrix */
SCIP_RETCODE SCIPsdpiGetLPNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros in the LP Matrix */
   )
{
   *nnonz = sdpi->lpnnonz;
   return SCIP_OKAY;
}

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either sdparraylength, sdpnnonz, sdpbegblock, sdprowind, sdpcolind and sdpval have to be 0 or NULL,
 *  or all of them have to be bigger than zero or non-NULL, the same is true for the lp-part.
 */
SCIP_RETCODE SCIPsdpiGetVarInfos(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to extract information for (between 1 and numvars) */
   int                   lastvar,            /**< last variable to extract information for (between 1 and numvars) */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int                   sdparraylength,     /**< length of sdparrays, if this is less than sdpnnonz an error will be returned */
   int*                  sdpnnonz,           /**< pointer to store the number of nonzero elements for the given variables in the sdp-constraints returned, or NULL */
   int*                  sdpbegvarblock,     /**< buffer to store start index of each block in sdpval-array (length (lastvar-firstvar + 1) * nblocks), or NULL */
   int*                  sdprowind,          /**< buffer to store row indices of sdp constraint matrix entries, or NULL */
   int*                  sdpcolind,          /**< buffer to store column indices of sdp constraint matrix entries, or NULL */
   SCIP_Real*            sdpval,             /**< buffer to store values of sdp constraint matrix entries, or NULL */
   int                   lparraylength,      /**< length of lparrays, if this is less than lpnnonz an error will be returned */
   int*                  lpnnonz,            /**< pointer to store the number of nonzero elements of the lp-constraints returned, or NULL */
   int*                  lprowind,           /**< buffer to store row indices of lp constraint matrix entries, or NULL */
   int*                  lpcolind,           /**< buffer to store column indices of lp constraint matrix entries (these will refer to the original variable indices), or NULL */
   SCIP_Real*            lpval               /**< buffer to store values of lp constraint matrix entries, or NULL */
   )
{
   int i;
   int numvars;
   int block;
   int lastiterationindex;
   int ind;
   int firstvarlpind;
   int lastvarlpind;

   assert ( firstvar > 0 );
   assert ( firstvar <= sdpi->nvars );
   assert ( lastvar > 0 );
   assert ( lastvar <= sdpi->nvars );

   numvars = lastvar - firstvar + 1;

   if (lb != NULL)
   {
      assert ( ub != NULL );
      for (i = firstvar - 1; i < lastvar; i++) /* indexshift because the arrays start at 0 */
      {
         lb[i] = sdpi->lb[i];
         ub[i] = sdpi->ub[i];
      }
   }

   if (sdparraylength > 0)
   {
      assert ( sdpnnonz != NULL );
      assert ( sdpbegvarblock != NULL );
      assert ( sdprowind != NULL );
      assert ( sdpcolind != NULL );
      assert ( sdpval != NULL );

      /* count the number of nonzeroes for the given variables */
      for (block = 0; block < sdpi->nsdpblocks; block++)
      {
         if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars)
         {
            *sdpnnonz = *sdpnnonz + sdpi->sdpnnonz - sdpi->sdpbegvarblock[block * sdpi->nvars + (firstvar - 1)];
         }
         else
         {
            *sdpnnonz = *sdpnnonz + sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar]
                                                         - sdpi->sdpbegvarblock[block * sdpi->nvars + (firstvar - 1)];
         }
      }

      /* check if the arrays are long enough to store all information */
      if (*sdpnnonz > sdparraylength)
      {
         SCIPerrorMessage("In SCIPsdpiGetVarInfos the given SDP-arrays were too short to store all information");
         SCIPABORT();
         return SCIP_ERROR;
      }

      ind = 0;
      for (block = 0; block < sdpi->nsdpblocks; block ++)
      {
         /* compute sdpbegvarblock */
         sdpbegvarblock[0] = 0;
         for (i = 0; i < numvars; i++)
         {
            if (block < sdpi->nsdpblocks - 1 || i < numvars + 1) /* for the last block-var-combination nothing has to be done, as this would only result
                                                                  * in sdpnnonz */
            {
               sdpbegvarblock[block * numvars + i + 1] = sdpbegvarblock[block * numvars + i]
                                                                    + sdpi->sdpbegvarblock[block * sdpi->nvars + (firstvar - 1) + i + 1]
                                                                                           - sdpi->sdpbegvarblock[block * sdpi->nvars + (firstvar - 1) + i];
            }
         }

         /* copy the nonzeroes in the corresponding arrays */
         if (block == sdpi->nsdpblocks - 1 && lastvar == sdpi->nvars)
         {
            lastiterationindex = sdpi->sdpnnonz;
         }
         else
         {
            lastiterationindex = sdpi->sdpbegvarblock[block * sdpi->nvars + lastvar];
         }
         for (i = sdpi->sdpbegvarblock[block * sdpi->nvars + (firstvar - 1)]; i < lastiterationindex; i++)
         {
            sdprowind[ind] = sdpi->sdprowind[i];
            sdpcolind[ind] = sdpi->sdpcolind[i];
            sdpval[ind] = sdpi->sdpval[i];
            ind++;
         }
      }
   }

   if (lparraylength > 0)
   {
      assert ( lpnnonz != 0 );
      assert ( lprowind != 0 );
      assert ( lpcolind != 0 );
      assert ( lpval != 0 );

      /* for copying the lp-nonzeroes first sort the arrays by columns, so that the entries corresponding to the asked for variables are in one block */
      SCIPsortIntIntReal(sdpi->lpcolind, sdpi->lprowind, sdpi->lpval, sdpi->lpnnonz);

      /*iterate over the lpcolind array to find the first index belonging to a one of the variables */
      firstvarlpind = -1; /* if this stays at -1 the variables weren't part of the LP block */
      for (i = 0; i < sdpi->lpnnonz; i++)
      {
         if (sdpi->lpcolind[i] >= firstvar && sdpi->lpcolind[i] <= lastvar)
         {
            firstvarlpind = i;
            lastvarlpind = i;
            i++;
            break;
         }
      }

      if (firstvarlpind > -1) /* if this is still -1 there are no entries */
      {
      /* now find the last occurence of one of the variable (as these are sorted all in between also belong to these vars) */
         while (i < sdpi->lpnnonz && sdpi->lpcolind[i] <= lastvar)
         {
            lastvarlpind++;
            i++;
         }

         *lpnnonz = lastvarlpind - firstvarlpind + 1;

         /* check if the provided arrays are sufficient for copying everything to them */
         if (lparraylength < *lpnnonz)
         {
            SCIPerrorMessage("In SCIPsdpiGetVarInfos the given LP-arrays were too short to store all information");
            SCIPABORT();
            return SCIP_ERROR;
         }

         /* now copy the lpnonzeroes */
         ind = 0;
         for (i = firstvarlpind; i <= lastvarlpind; i++)
         {
            lprowind[ind] = sdpi->lprowind[i];
            lpcolind[ind] = sdpi->lpcolind[i];
            lpval[ind] = sdpi->lpval[i];
            ind++;
         }
      }
   }

   return SCIP_OKAY;
}

/** gets LP rows from SDP problem object; the arrays have to be large enough to store all values.
 *  rhs can be null,
 *  either nnonz, begrow, colind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPsdpiGetLPRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first LP row to get from SDP (index between 1 and nlpcons) */
   int                   lastrow,            /**< last LP row to get from SDP (index between 1 and nlpcons) */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int                   arraylength,        /**< length of the rowind and colind arrays, if this is less than nnonz an ERROR will be returned */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int nrows;
   int i;
   int firstrowind;
   int lastrowind;
   int ind;

   assert ( firstrow > 0 );
   assert ( lastrow >= firstrow );
   assert ( lastrow <= sdpi->nlpcons );

   nrows = lastrow - firstrow + 1;

   if (rhs != NULL)
   {
      for (i = 0; i < nrows; i++)
      {
         rhs[(firstrow - 1) + i] = sdpi->lprhs[i]; /* indexshift because arrays start at 0 */
      }
   }

   if (arraylength > 0)
   {

      /* for deleting and reordering the lpnonzeroes, the arrays first have to be sorted to have the rows to be deleted together */
      SCIPsortIntIntReal(sdpi->lprowind, sdpi->lpcolind, sdpi->lpval, sdpi->lpnnonz); /* sort all arrays by non-decreasing row indices */

      firstrowind = -1;
      /*iterate over the lprowind array to find the first index belonging to one of the rows */
      for (i = 0; i < sdpi->lpnnonz; i++)
      {
         if (sdpi->lprowind[i] >= firstrow && sdpi->lprowind[i] <= lastrow) /* the and part makes sure that there actually were some nonzeroes in these rows */
         {
            firstrowind = i;
            lastrowind = i;
            i++;
            break;
         }
      }

      if (firstrowind > -1) /* if this is still -1 there are no nonzeroes for the given rows */
      {
         /* now find the last occurence of one of the rows (as these are sorted all in between also belong to these rows) */
         while (i < sdpi->lpnnonz && sdpi->lprowind[i] <= lastrow)
         {
            lastrowind++;
            i++;
         }

         *nnonz = lastrowind - firstrowind + 1;

         /* check if given arrays are long enough */
         if (*nnonz > arraylength)
         {
            SCIPerrorMessage("In SCIPsdpiGetLPRows the given LP-arrays were too short to store all information");
            SCIPABORT();
            return SCIP_ERROR;
         }

         /* copy the entries */
         ind = 0;
         for (i = firstrowind; i <= lastrowind; i++)
         {
            rowind[ind] = sdpi->lprowind[i];
            colind[ind] = sdpi->lpcolind[i];
            val[ind] = sdpi->lpval[i];
            ind++;
         }
      }
   }

   return SCIP_OKAY;
}

/** gets a number SDP blocks; the arrays have to be large enough to store all values.
 *  either constnnonz, constbegblock, constrow, constcolind, and constval have to be NULL, or all of them have to be non-NULL, same for the non-constant parts.
 */
SCIP_RETCODE SCIPsdpiGetSDPBlocks(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstblock,         /**< first SDP block to get from SDP (index between 1 and nblocks) */
   int                   lastblock,          /**< last SDP block to get from SDP (index between 1 and nblocks) */
   int                   constarraylength,   /**< how many elements can be stored in the const arrays, if this is less than constnnonz an error will be thrown */
   int*                  constnnonz,         /**< pointer to store the number of nonzero elements returned for the constant part, or NULL */
   int*                  constbegblock,      /**< buffer to store the start indices of the different blocks in the constval-array, or NULL */
   int*                  constrowind,        /**< buffer to store row indices of entries of constant matrices, or NULL */
   int*                  constcolind,        /**< buffer to store column indices of entries of constant matrices, or NULL */
   SCIP_Real*            constval,           /**< buffer to store values of entries of constant matrices, or NULL */
   int                   arraylength,        /**< how many elements can be stored in the arrays, if this is less than nnonz an error will be thrown */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  begvarblock,        /**< buffer to store the start indices of the different block/var-combinations in the val-array, or NULL */
   int*                  rowind,             /**< buffer to store row indices of constraint matrix entries, or NULL */
   int*                  colind,             /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   int i;

   assert ( 0 < firstblock );
   assert ( firstblock <= lastblock );
   assert ( lastblock <= sdpi->nsdpblocks);

   if (constarraylength > 0)
   {

      assert ( constnnonz != NULL );
      assert ( constbegblock != NULL );
      assert ( constrowind != NULL );
      assert ( constcolind != NULL );
      assert ( constval != NULL );

      if (lastblock == sdpi->nsdpblocks)
      {
         *constnnonz = sdpi->sdpconstnnonz - - sdpi->sdpconstbegblock[firstblock - 1]; /* indexshift because arrays start at 0 */
      }
      else
      {
         *constnnonz = sdpi->sdpconstbegblock[lastblock] - sdpi->sdpconstbegblock[firstblock - 1]; /* again indexshift because arrays start at 0 */
      }

      /* check if given arrays are sufficiently long */
      if (*constnnonz > constarraylength)
      {
         SCIPerrorMessage("In SCIPsdpiGetSDPBlocks the given const-arrays were too short to store all information");
         SCIPABORT();
         return SCIP_ERROR;
      }

      /* compute constbegblock */
      for (i = 0; i <= lastblock - firstblock; i++)
      {
         constbegblock[i] = sdpi->sdpconstbegblock[(firstblock - 1) + i] - sdpi->sdpconstbegblock[(firstblock - 1)]; /* starting index of each block is the starting index in the
                                                                                                                   * original problem minus that of the first block taken */
      }

      /* copy nonzeroes */
      for (i = 0; i < *constnnonz; i++)
      {
         constrowind[i] = sdpi->sdpconstrowind[sdpi->sdpconstbegblock[(firstblock - 1)] + i];
         constcolind[i] = sdpi->sdpconstcolind[sdpi->sdpconstbegblock[(firstblock - 1)] + i];
         constval[i] = sdpi->sdpconstval[sdpi->sdpconstbegblock[(firstblock - 1)] + i];
      }
   }

   if (arraylength > 0)
   {
      assert ( nnonz != NULL );
      assert ( begvarblock != NULL );
      assert ( rowind != NULL );
      assert ( colind != NULL );
      assert ( val != NULL );

      if (lastblock == sdpi->nsdpblocks)
      {
         *nnonz = sdpi->sdpnnonz - sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars]; /* indexshift because arrays start at 0 */
      }
      else
      {
         *nnonz = sdpi->sdpbegvarblock[lastblock * sdpi->nvars] - sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars]; /* indexshift because arrays start at 0 as usual*/
      }

      /* check if given arrays are long enough */
      if (*nnonz > arraylength)
      {
         SCIPerrorMessage("In SCIPsdpiGetSDPBlocks the given arrays for the non-constant part were too short to store all information");
         SCIPABORT();
         return SCIP_ERROR;
      }

      /* compute begvarblock */
      for (i = 0; i < (lastblock - firstblock +1) * sdpi->nvars; i++)
      {
         begvarblock[i] = sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars + i] - sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars];
      }

      /* copy nonzeroes */
      for (i = 0; i < *nnonz; i++)
      {
         rowind[i] = sdpi->sdprowind[sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars] + i];
         colind[i] = sdpi->sdpcolind[sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars] + i];
         val[i] = sdpi->sdpval[sdpi->sdpbegvarblock[(firstblock - 1) * sdpi->nvars] + i];
      }
   }

   return SCIP_OKAY;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get objective coefficient for (index between 1 and nvars) */
   int                   lastvar,            /**< last variable to get objective coefficient for (index between 1 and nvars) */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   int i;

   assert ( firstvar > 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar <= sdpi->nvars);

   for (i = 0; i < lastvar - firstvar + 1; i++)
   {
      vals[i] = sdpi->obj[(firstvar - 1) + i]; /* index shift because arrays start at 0 */
   }
   return SCIP_OKAY;
}

/** gets current variable bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstvar,           /**< first variable to get bounds for (index between 1 and nvars) */
   int                   lastvar,            /**< last variable to get bounds for (index between 1 and nvars) */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   int i;

   assert ( firstvar > 0 );
   assert ( firstvar <= lastvar );
   assert ( lastvar <= sdpi->nvars);

   for (i = 0; i < lastvar - firstvar + 1; i++)
   {
      if (lbs != NULL)
      {
         lbs[i] = sdpi->lb[(firstvar - 1) + i]; /* index shift because arrays start at 0 */
      }
      if (ubs != NULL)
      {
         ubs[i] = sdpi->ub[(firstvar - 1) + i]; /* index shift because arrays start at 0 */
      }
   }
   return SCIP_OKAY;
}

/** gets current right hand sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetRhSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for (index between 1 and nlpcons) */
   int                   lastrow,            /**< last row to get sides for (index between 1 and nlpcons) */
   SCIP_Real*            rhss                /**< array to store right hand side values */
   )
{
   int i;

   assert ( firstrow > 0 );
   assert ( firstrow <= lastrow );
   assert ( lastrow <= sdpi->nlpcons);

   for (i = 0; i < lastrow - firstrow + 1; i++)
   {
      rhss[(firstrow - 1) + i] = sdpi->lprhs[i]; /* indexshift because arrays start at 0 */
   }

   return SCIP_OKAY;
}

/** gets a single coefficient of LP constraint matrix */
SCIP_RETCODE SCIPsdpiGetLPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient (index between 1 and nlpcons) */
   int                   col,                /**< column number of coefficient (index between 1 and nvars) */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;

   assert ( row > 0 );
   assert ( row <= sdpi->nlpcons );
   assert ( col > 0 );
   assert ( col <= sdpi->nvars);

   /* search for the entry */
   for (i = 0; i < sdpi->lpnnonz; i++)
   {
      if (sdpi->lpcolind[i] == col && sdpi->lprowind[i] == row)
      {
         *val = sdpi->lpval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetLPCoef no entry for the given column = %d and row =%d was found.", col, row);
   SCIPABORT();
   return SCIP_ERROR;
}

/** gets a single coefficient of constant SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPConstCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient (index between 1 and nsdpblocks) */
   int                   rowind,             /**< row number of coefficient (index between 1 and corresponding blocksize) */
   int                   colind,             /**< column number of coefficient (index between 1 and corresponding blocksize) */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;

   assert ( block > 0 );
   assert ( block <= sdpi->nsdpblocks );
   assert ( rowind > 0 );
   assert ( rowind <= sdpi->sdpblocksizes[block - 1] ); /* indexshift */
   assert ( colind > 0 );
   assert ( colind <= sdpi->sdpblocksizes[block - 1] ); /* indexshift again */

   row = rowind;
   col = colind;
   checkIfLowerTriang(&row, &col); /* Because the matrices are symmetric it doesn't matter if a position in the upper or lower triangular
                                  * was given, but only positions in the lower triangular path are saved in the corresponding arrays */

   /* search for the entry */
   if (block == sdpi->nsdpblocks)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpconstbegblock[block]; /* indexshift */
   }
   for (i = sdpi->sdpconstbegblock[block - 1]; i < lastiterationindex; i++) /* indexshift */
   {
      if (sdpi->sdpconstcolind[i] == col && sdpi->sdpconstrowind[i] == row)
      {
         *val = sdpi->sdpconstval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetSDPConstCoef no entry for the given column = %d and row = %d was found.", colind, rowind);
   SCIPABORT();
   return SCIP_ERROR;
}

/** gets a single coefficient of SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetSDPCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   block,              /**< block index of coefficient  (index between 1 and nsdpblocks) */
   int                   var,                /**< variable index of coefficient, meaning the i in \f A_i^j \f, in val-array, or NULL */
   int                   rowind,             /**< row number of coefficient (index between 1 and corresponding blocksize) */
   int                   colind,             /**< column number of coefficient (index between 1 and corresponding blocksize) */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   int i;
   int lastiterationindex;
   int row;
   int col;

   assert ( block > 0 );
   assert ( block <= sdpi->nsdpblocks );
   assert ( var > 0 );
   assert ( var <= sdpi->nvars );
   assert ( rowind > 0 );
   assert ( rowind <= sdpi->sdpblocksizes[block - 1] ); /* indexshift */
   assert ( colind > 0 );
   assert ( colind <= sdpi->sdpblocksizes[block - 1] ); /* indexshift again */

   row = rowind;
   col = colind;
   checkIfLowerTriang(&row, &col); /* Because the matrices are symmetric it doesn't matter if a position in the upper or lower triangular
                                  * was given, but only positions in the lower triangular path are saved in the corresponding arrays */

   /* search for the entry */
   if (block == sdpi->nsdpblocks && var == sdpi->nvars)
   {
      lastiterationindex = sdpi->sdpconstnnonz;
   }
   else
   {
      lastiterationindex = sdpi->sdpbegvarblock[(block-1) * sdpi->nvars + var]; /* indexshift */
   }
   for (i = sdpi->sdpbegvarblock[(block - 1) * sdpi->nvars + var - 1]; i < lastiterationindex; i++) /* indexshift */
   {
      if (sdpi->sdpcolind[i] == col && sdpi->sdprowind[i] == row)
      {
         *val = sdpi->sdpval[i];
         return SCIP_OKAY;
      }
   }

   /* if this is reached no corresponding entry in the LP-constraints was found */
   SCIPerrorMessage("In SCIPsdpiGetSDPConstCoef no entry for the given column=%d and row=%d was found.", colind, rowind);
   SCIPABORT();
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
   int block;

#ifdef SCIP_DEBUG
   DSDPTerminationReason* reason; /* this will later be used to check if DSDP converged */
   BMSallocBlockMemory(sdpi->blkmem, &reason);
#endif

   /* insert data */

   SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", sdpi->sdpid);

   if (sdpi->dsdp != NULL)
   {
      DSDP_CALL(DSDPDestroy(sdpi->dsdp)); /* if there already exists a DSDP-instance, destroy the old one */
   }

   DSDP_CALLM(DSDPCreate(sdpi->nvars, &(sdpi->dsdp)));
   DSDP_CALLM(DSDPCreateSDPCone(sdpi->dsdp, sdpi->nsdpblocks, &(sdpi->sdpcone)));
   DSDP_CALLM(DSDPCreateLPCone(sdpi->dsdp, &(sdpi->lpcone)));
   DSDP_CALLM(DSDPCreateBCone(sdpi->dsdp, &(sdpi->bcone)));

   DSDP_CALLM(BConeAllocateBounds(sdpi->bcone,2*sdpi->nvars)); /*allocate memory for lower and upper bounds */

   for (i = 0; i < sdpi->nvars; i++)
   {
      DSDP_CALL(DSDPSetDualObjective(sdpi->dsdp, i+1, -1 * sdpi->obj[i])); /*insert objective value, DSDP counts from 1 to n instead of 0 to n-1,
                                                                                               * *(-1) because DSDP maximizes instead of minimizing */
      if (!SCIPsdpiIsInfinity(sdpi, sdpi->lb[i]))
      {
         DSDP_CALL(BConeSetLowerBound(sdpi->bcone, i+1, sdpi->lb[i])); /*insert lower bound, DSDP counts from 1 to n instead of 0 to n-1 */
      }
      if (!SCIPsdpiIsInfinity(sdpi, sdpi->ub[i]))
      {
         DSDP_CALL(BConeSetUpperBound(sdpi->bcone, i+1, sdpi->ub[i])); /*insert upper bound, DSDP counts from 1 to n instead of 0 to n-1 */
      }
   }

#ifdef SCIP_DEBUG
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpi->bcone);
#endif

   /*start inserting the constant matrix*/
   if ( sdpi->nsdpblocks > 0 && sdpi->sdpconstnnonz > 0 )
   {

      /*allocate memory*/
      /*This needs to be one long array, because DSDP uses it for solving, so all nonzeros have to be in it, and it may not be freed before the problem is solved. */

      /*indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpconstind, sdpi->sdpconstnnonz);
      /*values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMSallocBlockMemoryArray(sdpi->blkmem, &dsdpconstval, sdpi->sdpconstnnonz);

      for(i = 0; i < sdpi->nsdpblocks; i++)
      {
         DSDP_CALL(SDPConeSetBlockSize(sdpi->sdpcone, i, sdpi->sdpblocksizes[i])); /*set the blocksizes (blocks are counted from 0 to m-1) */
      }

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

         /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
         SCIPsortIntReal(dsdpconstind + sdpi->sdpconstbegblock[block], dsdpconstval + sdpi->sdpconstbegblock[block], blocknnonz);

         DSDP_CALL(SDPConeSetASparseVecMat(sdpi->sdpcone, block, 0, sdpi->sdpblocksizes[block], 1, 0, dsdpconstind + sdpi->sdpconstbegblock[block],
            dsdpconstval + sdpi->sdpconstbegblock[block], blocknnonz));   /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block
                                                                           * by adding sdpi->sdpconstbegblock[block] */
      }
   }


   /*start inserting the other SDP-Constraint-Matrices */
   if(sdpi->nsdpblocks > 0 && sdpi->sdpnnonz > 0)
   {
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
               dsdpval[ind] = -1 * sdpi->sdpval[ind];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
            }

            /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
            SCIPsortIntReal(dsdpind + sdpi->sdpbegvarblock[sdpi->nvars * block + var], dsdpval + sdpi->sdpbegvarblock[sdpi->nvars * block + var], aij_nnonz);

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

      /*to compute lpbegcol the column indices need to be sorted, for this they are copied in an extra array */
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

   DSDP_CALL(DSDPSetGapTolerance(sdpi->dsdp, 1e-3)); /* set DSDP's tolerance */


   DSDP_CALLM(DSDPSetup(sdpi->dsdp));
   DSDP_CALL(DSDPSolve(sdpi->dsdp));

   DSDP_CALL(DSDPComputeX(sdpi->dsdp)); /*computes X and determines feasibility and unboundedness of the solution */
   sdpi->solved = TRUE;

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

   switch ( *reason ) /* TODO: perhaps also check for feasibility and call the penalty-method here in that case */
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
         printf("Numerical Trouble in DSDP! \n");
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
         break;

      default:
         SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
         BMSfreeBlockMemory(sdpi->blkmem, &reason);
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
   return sdpi->solved;
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType* pdfeasible;

   CHECK_IF_SOLVED(sdpi);

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
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
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

   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      if (!checkFeasibility(sdpi)){
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible, but the dual solution actually isn't feasible\n");
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         return TRUE;
      }
      else
      {
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible, but the dual solution actually is feasible\n");
         BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
         return FALSE;
      }
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

   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpi->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &pdfeasible);
      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
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
   CHECK_IF_SOLVED(sdpi);

   return !SCIPsdpiIsPrimalInfeasible(sdpi);
}

/** returns TRUE iff SDP is proven to be dual unbounded */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   CHECK_IF_SOLVED(sdpi);

   return SCIPsdpiIsPrimalInfeasible(sdpi);
}

/** returns TRUE iff SDP is proven to be dual infeasible */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   CHECK_IF_SOLVED(sdpi);

   return SCIPsdpiIsPrimalUnbounded(sdpi);
}

/** returns TRUE iff SDP is proven to be dual feasible */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   CHECK_IF_SOLVED(sdpi);

   return !SCIPsdpiIsPrimalUnbounded(sdpi);
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

   if (*reason == DSDP_CONVERGED)
   {
      BMSfreeBlockMemory(sdpi->blkmem, &reason);
      return TRUE;
   }
   else
   {
      double* pobj;
      double* dobj;

      BMSfreeBlockMemory(sdpi->blkmem, &reason);

      /* if it didn't converge check the optimality gap */
      BMSallocBlockMemory(sdpi->blkmem, &pobj);
      BMSallocBlockMemory(sdpi->blkmem, &dobj);

      DSDP_CALL(DSDPGetPObjective(sdpi->dsdp, pobj));
      DSDP_CALL(DSDPGetDObjective(sdpi->dsdp, dobj));

      if (((abs(*pobj - *dobj))/ *dobj) < epsilon)
      {
         BMSfreeBlockMemory(sdpi->blkmem, &pobj);
         BMSfreeBlockMemory(sdpi->blkmem, &dobj);
         return TRUE;
      }
      else
      {
         return FALSE;
      }
   }
/* TODO: also check for primal feasibility, as this is also needed for optimality */
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   DSDPTerminationReason* reason;

   CHECK_IF_SOLVED(sdpi);

   BMSallocBlockMemory(sdpi->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpi->dsdp, reason));

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

   CHECK_IF_SOLVED(sdpi);

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
   CHECK_IF_SOLVED(sdpi);

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
   CHECK_IF_SOLVED(sdpi);

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
   CHECK_IF_SOLVED(sdpi);
   DSDP_CALL(DSDPGetIts(sdpi->dsdp, iterations));
   return SCIP_OKAY;
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
   return 10000000; /* DSDP has implicit bounds of +- 10⁷ for y */
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   )
{
   return (val <= -10000000) || (val >= 10000000);
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
