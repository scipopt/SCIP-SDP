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

#define SCIP_DEBUG
#define SCIP_MORE_DEBUG

/**@file   sdpi_dsdp.c
 * @brief  interface for dsdp
 * @author Tristan Gally
 */

#include <assert.h>

#include "sdpi/sdpi.h"

#include "dsdp5.h"                           /* for DSDPUsePenalty, etc */

#include "blockmemshell/memory.h"            /* for memory allocation */
#include "scip/def.h"                        /* for SCIP_Real, _Bool, ... */
#include "scip/pub_misc.h"                   /* for sorting */

/* calls a DSDP-Function and transforms the return-code to a SCIP_ERROR if needed */
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

/* same as DSDP_CALL, but this will be used for initialization methods with memory allocation and return a SCIP_NOMEMORY if an error is produced */
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

/* this will be called in all functions that want to access solution information to check if the problem was solved since the last change of the problem */
#define CHECK_IF_SOLVED(sdpisolver)  do                                                                             \
                        {                                                                                     \
                           if (!(sdpisolver->solved))                                                               \
                           {                                                                                  \
                              SCIPerrorMessage("Tried to access solution information ahead of solving! \n");  \
                              SCIPABORT();                                                                    \
                              return SCIP_ERROR;                                                              \
                           }                                                                                  \
                        }                                                                                     \
                        while( FALSE )

struct SCIP_SDPiSolver
{
   DSDP                  dsdp;               /**< solver-object */
   SDPCone               sdpcone;            /**< sdpcone-object of DSDP for handling SDP-constraints */
   LPCone                lpcone;             /**< lpcone-object of DSDP for handling LP-constraints */
   BCone                 bcone;              /**< bcone-object of DSDP to add variable bounds to */
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehandler to printing messages, or NULL */
   BMS_BLKMEM*           blkmem;             /**< block memory */
   SCIP_Bool             solved;             /**< was the SDP solved since the problem was last changed */
   int                   sdpcounter;         /**< used for debug messages */
};

static double epsilon    = 1e-6;             /**< this is used for checking if primal and dual objective are equal */

/*
 * Local Functions
 */

/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
static int compLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}

/**
 * For given row and column (i,j) checks if i >= j, so that i and j give a position in the lower
 * triangular part, otherwise i and j will be switched. This function will be called whenever a position in a symmetric matrix
 * is given, to prevent problems if position (i,j) is given but later (j,i) should be changed.
 */
static void ensureLowerTriangular(
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

/*
 * Miscellaneous Methods
 */

/**@name Miscellaneous Methods */
/**@{ */


/** gets name of SDP solver, getting version doesn't seem to be supported by DSDP */
const char* SCIPsdpiSolverGetSolverName(
   void
   )
{
   return "DSDP";
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiSolverGetSolverDesc(
   void
   )
{
   return "Dual-Scaling Interior Point Solver for Semidefinite Programming developed by Steve Benson, Yinyu Ye, and Xiong Zhang (http://www.mcs.anl.gov/hs/software/DSDP/)";
}

/** Does the solver have a way to solve a penalty formulation on its own or must one be provided */
SCIP_Bool SCIPsdpiSolverKnowsPenalty(
   void
   )
{
   return TRUE;
}

/** gets pointer for SDP solver - use only with great care
 *
 *  The behavior of this function depends on the solver and its use is
 *  therefore only recommended if you really know what you are
 *  doing. In general, it returns a pointer to the SDP solver object.
 */
void* SCIPsdpiSolverGetSolverPointer(
   SCIP_SDPISOLVER*      sdpisolver           /**< pointer to an SDP interface solver structure */
   )
{
   assert( sdpisolver != NULL );
   return (void*) sdpisolver->dsdp;
}

/**@} */


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiSolverCreate(
   SCIP_SDPISOLVER**     sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   DSDP newdsdp;
   SDPCone newsdpcone;
   LPCone newlpcone;
   BCone newbcone;

   assert ( sdpisolver != NULL );
   assert ( blkmem != NULL );

   /* these will be properly initialized only immediatly prior to solving because DSDP and the SDPCone need information about the number
    * of variables and sdpblocks during creation */
   newdsdp = NULL;
   newsdpcone = NULL;
   newlpcone = NULL;
   newbcone = NULL;

   SCIPdebugMessage("Calling SCIPsdpiCreate (%d)\n",nextsdpid);

   BMSallocBlockMemory(blkmem, sdpisolver);

   (*sdpisolver)->messagehdlr = messagehdlr;
   (*sdpisolver)->blkmem = blkmem;
   (*sdpisolver)->dsdp = newdsdp;
   (*sdpisolver)->sdpcone = newsdpcone;
   (*sdpisolver)->lpcone = newlpcone;
   (*sdpisolver)->bcone = newbcone;
   (*sdpisolver)->solved = FALSE;
   (*sdpisolver)->sdpcounter = 0;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiSolverFree(
   SCIP_SDPISOLVER**     sdpisolver          /**< pointer to an SDP interface solver structure */
   )
{
   assert ( sdpisolver != NULL );
   assert ( *sdpisolver != NULL );

   if (((*sdpisolver)->dsdp) != NULL)
   {
   DSDP_CALL(DSDPDestroy((*sdpisolver)->dsdp));
   }

   BMSfreeBlockMemory((*sdpisolver)->blkmem, sdpisolver);

   return SCIP_OKAY;
}

/**@} */


/*
 * Solving Methods
 */

/**@name Solving Methods */
/**@{ */

/** solves the SDP */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolve(
   SCIP_SDPISOLVER*      sdpisolver,          /**< SDP interface solver structure */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                *  number of entries  of sdpconst row/col/val [i] */
   const int**          sdpconstrow,        /**< pointers to row-indices for each block */
   const int**          sdpconstcol,        /**< pointers to column-indices for each block */
   const SCIP_Real**     sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpnblockvarnonz,   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i * nvars + j] */
   const int**          sdprow,             /**< pointer to the row-indices for each block and variable */
   const int**          sdpcol,             /**< pointer to the column-indices for each block and variable */
   const SCIP_Real**     sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
   )
{
   return SCIP_CALL(SCIPsdpiSolverLoadAndSolveWithPenalty(sdpisolver, 0.0, TRUE, nvars, obj, lb, ub, nsdpblocks, sdpblocksizes, sdpconstnnonz, sdpconstnblocknonz, sdpconstrow,
         sdpconstcol, sdpconstval, sdpnnonz, sdpnblockvarnonz, sdprow, sdpcol, sdpval, nlpcons, lprhs, lpnnonz, lprowind, lpcolind, lpval));
}

/** solves the following penalty formulation of the SDP:
 *      \f{eqnarray*}{
 *      \min & & b^T y + \Gamma r \\
 *      \mbox{s.t.} & & \sum_{j=1}^n A_j^i y_j - A_0^i + r \cdot \mathbb{I} \succeq 0 \quad \forall i \leq m \\
 *      & & Dy \geq d \\
 *      & & l \leq y \leq u}
 *   \f
 *   alternatively withObj can be to false to set \f b \f to false and only check for feasibility (if the optimal
 *   objective value is bigger than 0 the problem is infeasible, otherwise it's feasible) */
SCIP_RETCODE SCIPsdpiSolverLoadAndSolveWithPenalty(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real             gamma,              /**< the penalty parameter above, needs to be >= 0 */
   SCIP_Bool             withObj,            /**< if this is false, the objective is set to 0 */
   int                   nvars,              /**< number of variables */
   const SCIP_Real*      obj,                /**< objective function values of variables */
   const SCIP_Real*      lb,                 /**< lower bounds of variables */
   const SCIP_Real*      ub,                 /**< upper bounds of variables */
   int                   nsdpblocks,         /**< number of SDP-blocks */
   const int*            sdpblocksizes,      /**< sizes of the SDP-blocks (may be NULL if nsdpblocks = sdpconstnnonz = sdpnnonz = 0) */
   int                   sdpconstnnonz,      /**< number of nonzero elements in the constant matrices of the SDP-Blocks */
   const int*            sdpconstnblocknonz, /**< number of nonzeros for each variable in the constant part, also the i-th entry gives the
                                                *  number of entries  of sdpconst row/col/val [i] */
   const int**          sdpconstrow,        /**< pointers to row-indices for each block */
   const int**          sdpconstcol,        /**< pointers to column-indices for each block */
   const SCIP_Real**     sdpconstval,        /**< pointers to the values of the nonzeros for each block */
   int                   sdpnnonz,           /**< number of nonzero elements in the SDP-constraint matrix */
   int*                  sdpnblockvarnonz,   /**< entry i * nvars + j gives the number of nonzeros for block i and variable j, this is exactly
                                               *  the number of entries of sdp row/col/val [i * nvars + j] */
   const int**          sdprow,             /**< pointer to the row-indices for each block and variable */
   const int**          sdpcol,             /**< pointer to the column-indices for each block and variable */
   const SCIP_Real**     sdpval,             /**< values of SDP-constraint matrix entries (may be NULL if sdpnnonz = 0) */
   int                   nlpcons,            /**< number of LP-constraints */
   const SCIP_Real*      lprhs,              /**< right hand sides of LP rows (may be NULL if nlpcons = 0) */
   int                   lpnnonz,            /**< number of nonzero elements in the LP-constraint matrix */
   const int*            lprowind,           /**< row-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const int*            lpcolind,           /**< column-index for each entry in lpval-array (may be NULL if lpnnonz = 0) */
   const SCIP_Real*      lpval               /**< values of LP-constraint matrix entries (may be NULL if lpnnonz = 0) */
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
   int pos;
   int ind;
   int block;
   int aij_nnonz;
   int blocknnonz;
   int startind;

#ifdef SCIP_DEBUG
   DSDPTerminationReason* reason; /* this will later be used to check if DSDP converged */
#endif

   assert ( sdpisolver != NULL );
   assert ( penaltyParam >= 0.0 );

   /* insert data */

   SCIPdebugMessage("Inserting Data into DSDP for SDP (%d) \n", sdpisolver->sdpcounter);

   if (sdpisolver->dsdp != NULL)
   {
      DSDP_CALL(DSDPDestroy(sdpisolver->dsdp)); /* if there already exists a DSDP-instance, destroy the old one */
   }

   DSDP_CALLM(DSDPCreate(nvars, &(sdpisolver->dsdp)));
   DSDP_CALLM(DSDPCreateSDPCone(sdpisolver->dsdp, nsdpblocks, &(sdpisolver->sdpcone)));
   DSDP_CALLM(DSDPCreateLPCone(sdpisolver->dsdp, &(sdpisolver->lpcone)));
   DSDP_CALLM(DSDPCreateBCone(sdpisolver->dsdp, &(sdpisolver->bcone)));

   for (i = 0; i < nvars; i++)
   {
      if (withObj)
      {
         DSDP_CALL(DSDPSetDualObjective(sdpisolver->dsdp, i+1, -1 * obj[i])); /* insert objective value, DSDP counts from 1 to n instead of 0 to n-1,
                                                                                  * *(-1) because DSDP maximizes instead of minimizing */
      }
      else
      {
         DSDP_CALL(DSDPSetDualObjective(sdpisolver->dsdp, i+1, 0));
      }
      if (!SCIPsdpiSolverIsInfinity(sdpisolver, lb[i]))
      {
         DSDP_CALL(BConeSetLowerBound(sdpisolver->bcone, i+1, lb[i])); /*insert lower bound, DSDP counts from 1 to n instead of 0 to n-1 */
      }
      if (!SCIPsdpiSolverIsInfinity(sdpisolver, sdpisolver->ub[i]))
      {
         DSDP_CALL(BConeSetUpperBound(sdpisolver->bcone, i+1, ub[i])); /*insert upper bound, DSDP counts from 1 to n instead of 0 to n-1 */
      }
   }

#ifdef SCIP_MORE_DEBUG
   SCIPdebugMessage("ATTENTION: BConeView shows the WRONG sign for the lower bound!\n");
   BConeView(sdpisolver->bcone);
#endif

   /* set blocksizes */
   for(i = 0; i < nsdpblocks; i++)
   {
      DSDP_CALL(SDPConeSetBlockSize(sdpisolver->sdpcone, i, sdpblocksizes[i])); /*set the blocksizes (blocks are counted from 0 to m-1) */
   }

   /* start inserting the constant matrix */
   if ( sdpconstnnonz > 0 )
   {
      assert ( nsdpblocks > 0 );
      assert ( sdpconstbegblock != NULL );
      assert ( sdpconstcolind != NULL );
      assert ( sdpconstrowind != NULL );
      assert ( sdpconstval != NULL );

      /*allocate memory*/
      /* DSDP uses these for solving, so they may not be freed before the problem is solved. */

      /* indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz);
      /* values given to DSDP, for this the original values are mutliplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz);

      ind = 0;

      for(block = 0; block < nsdpblocks; block++)
      {
         startind = ind; /* starting index of this block in the dsdpconst arrays */

         for(i = 0; i < sdpconstnblocknonz[block]; i++)
         {
            dsdpconstind[ind] = compLowerTriangPos(sdpconstrow[block][i], sdpconstcol[block][i]);
            dsdpconstval[ind] = -1 * sdpconstval[ind]; /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0^j)
                                                                * should be positive semidefinite */
            ind++;
         }

         /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
         SCIPsortIntReal(dsdpconstind + startind, dsdpconstval + startind, sdpconstnblocknonz[block]);

         DSDP_CALL(SDPConeSetASparseVecMat(sdpisolver->sdpcone, block, 0, sdpblocksizes[block], 1, 0, dsdpconstind + startind[block],
            dsdpconstval + startind, blocknnonz));   /* constant matrix is given as variable 0, the arrays are shifted to the first element of this block
                                                      * by adding sdpisolver->sdpconstbegblock[block] */
      }
   }


   /*start inserting the other SDP-Constraint-Matrices */
   if(sdpnnonz > 0)
   {
      int var;
      int k;

      assert ( nsdpblocks > 0 );
      assert ( sdpbegvarblock != NULL );
      assert ( sdpcolind != NULL );
      assert ( sdprowind != NULL );
      assert ( sdpval != NULL );

      /*allocate memory */
      /*This needs to be one long array, because DSDP uses it for solving so all nonzeros have to be in it and it may not be freed before the problem is solved. The distinct blocks/variables
       *(for the i,j-parts) are then given by dsdpind + sdpbegvarblock[nvars * block + var], which gives a pointer to the first array-element belonging to this block and then the number of
       *elements in this block is given to DSDP for iterating over it */

      /*indices given to DSDP, for this the elements in the lower triangular part of the matrix are labeled from 0 to n*(n+1)/2 -1 */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz);
      /*values given to DSDP, these will be multiplied by -1 because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz);

      ind = 0; /* this will be used for iterating over the nonzeroes */

      for(block = 0; block < nsdpblocks; block++)
      {
         for(var = 0; var < nvars; var++)
         {
            startind = ind;

            for (k = 0; k < sdpnblockvarnonz[block * nvars + var]; k++)
            {
               dsdpind[ind] = compLowerTriangPos(sdprow[k], sdpcol[k]);
               dsdpval[ind] = -1 * sdpval[k];  /* *(-1) because in DSDP -1* (sum A_i^j y_i - A_0) should be positive semidefinite */
               ind++;
            }

            /* sort the arrays for this Matrix (by non decreasing indices) as this might help the solving time of DSDP */
            SCIPsortIntReal(dsdpind + startind, dsdpval + startind, sdpnblockvarnonz[block * nvars + var]);

            DSDP_CALL(SDPConeSetASparseVecMat(sdpisolver->sdpcone, block, var + 1, sdpblocksizes[block], 1, 0, dsdpind + startind,
               dsdpval + startind, sdpnblockvarnonz[block * nvars + var])); /* var+1 is needed because DSDP indexes the vars from 1 to nvars (var 0 is the constant matrix), adding
                                                                             * sdpisolver->sdpbegvarblock[sdpisolver->nvars * block + var] shifts the arrays to the first nonzero belonging
                                                                             *  to this block and this variable */
         }
      }
      #ifdef SCIP_MORE_DEBUG
      SDPConeView2(sdpisolver->sdpcone);
      #endif
   }


   /*start inserting the LP constraints */
   if(nlpcons > 0 || lpnnonz > 0)
   {
      int* sortedlpcolind;
      int column;
      int constraint;

      assert ( nlpcons > 0 );
      assert ( lprhs != NULL );
      assert ( lpcolind != NULL );
      assert ( lprowind != NULL );
      assert ( lpval != NULL );

      /* memory allocation */

      /* these arrays are needed in DSDP during solving, so they may only be freed afterwards */
      /* dsdplpbegcol[i] gives the number of nonzeroes in column 0 (right hand side) till i-1 (i going from 1 till n, with extra entries 0 (always 0) and n+1 (always lpcons + lpnnonz)) */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, nvars + 2);
      /* dsdplprowind saves the row indices of the LP for DSDP */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplprowind, nlpcons + lpnnonz); /*length is lpnnonz + nlpcons, because right hand sides are also included in the vector */
      /* values given to DSDP */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz); /*length is lpnnonz + nlpcons, because right hand sides are also included in the vector */

      /* compute lpbegcol */

      /* to compute lpbegcol the column indices need to be sorted, for this they are copied in an extra array */
      BMSallocBlockMemoryArray(sdpisolver->blkmem, &sortedlpcolind, lpnnonz);

      for(i = 0; i < lpnnonz; i++)
      {
         assert ( lpcolind[i] >= 0 );
         assert ( lpcolind[i] < nvars );
         sortedlpcolind[i] = lpcolind[i];
         assert ( lprowind[i] >= 0 );
         assert ( lprowind[i] < nlpcons );
         dsdplprowind[nlpcons + i] = lprowind[i];  /* the first nlpcons entries will be used for the right hand sides, so the matrix-entries are copied in the later ones */
         dsdplpval[nlpcons + i] = -1 * lpval[i];   /* the first nlpcons entries will be used for the right hand sides, so the matrix-entries are copied in the later ones, *(-1) is needed, because
                                                                * DSDP wants <= instead of >= */
      }



      SCIPsortIntIntReal(sortedlpcolind, dsdplprowind + nlpcons, dsdplpval + nlpcons, lpnnonz); /* all three arrays should now be sorted by non-decreasing column-indizes, for dsdplprowind and dsdplpval
      the sorting starts at position nlpcons (the first index is shifted by nlpcons), because the earlier entries are still empty and will only later be used for the right hand sides */

      dsdplpbegcol[0] = 0;
      dsdplpbegcol[1] = nlpcons; /* the first nlpcons indices are used to save the right hand sides */
      ind = 0; /* this will be used for traversing the sortedlpcolind-array */

      for(column = 1; column < nvars + 1; column++) /*columns are indexed 1 to nvars in dsdplpbegcol */
      {
         dsdplpbegcol[column+1] = dsdplpbegcol[column];  /* each new column can't start before the last one */
         while(ind < lpnnonz && sortedlpcolind[ind] == column - 1) /* look at all indices whose column index matches the current column
                                                                          * in lpcolind the columns are indexed 0 to nvars - 1, while here the indexing starts at 1*/
         {
            dsdplpbegcol[column + 1]++;   /*for each element with (lpcolind = current column) an additional entry in dsdplpval is needed, so the next column starts one spot later */
            ind++;
         }
      }

      assert(dsdplpbegcol[nvars + 1] == lpnnonz + nlpcons);

      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &sortedlpcolind, lpnnonz); /*this was only needed to sort the column indices and compute dsdplpbegcol */

      for(column = 1; column < nvars + 1; column++)
      {
         SCIPsortIntReal(dsdplprowind + dsdplpbegcol[column], dsdplpval + dsdplpbegcol[column], dsdplpbegcol[column + 1] - dsdplpbegcol[column]);
         /*sort all the entries belonging to the same column by their row numbers */
      }

      /* insert the right-hand-side-values into the first entries of the dsdplp-arrays */
      for(constraint = 0; constraint < nlpcons; constraint++)
      {
         dsdplprowind[constraint] = constraint; /* the row index of each constraint is the index of the constraint */
         dsdplpval[constraint] = -1 * lprhs[constraint]; /* insert rhs values, *(-1) is needed, because DSDP wants <= instead of >= */
      }

      DSDP_CALL(LPConeSetData(sdpisolver->lpcone, sdpisolver->nlpcons, dsdplpbegcol, dsdplprowind, dsdplpval));
      #ifdef SCIP_MORE_DEBUG
      LPConeView(sdpisolver->lpcone);
      #endif
   }

   SCIPdebugMessage("Calling DSDP-Solve for SDP (%d) \n", sdpisolver->sdpcounter);

   DSDP_CALL(DSDPSetGapTolerance(sdpisolver->dsdp, 1e-3));  /* set DSDP's tolerance for duality gap */
   DSDP_CALL(DSDPSetRTolerance(sdpisolver->dsdp, 1e-6));    /* set DSDP's tolerance for the psd-constraints */


   /* set the penalty parameter */
   if (penaltyParam != 0.0) /* in sdpisolverSolve this is called with an exact 0 */
   {
      DSDPSetPenaltyParameter(sdpisolver->dsdp, penaltyParam);
      DSDPUsePenalty(sdpisolver->dsdp, 1);
   }

   DSDP_CALLM(DSDPSetup(sdpisolver->dsdp));
   DSDP_CALL(DSDPSolve(sdpisolver->dsdp));

   DSDP_CALL(DSDPComputeX(sdpisolver->dsdp)); /*computes X and determines feasibility and unboundedness of the solution */
   sdpisolver->solved = TRUE;

   /*these arrays were used to give information to DSDP and were needed during solving and for computing X, so they may only be freed now*/
   if ( sdpisolver->sdpconstnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstind, sdpconstnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpconstval, sdpconstnnonz);
   }

   if ( sdpisolver->sdpnnonz > 0 )
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpind, sdpnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdpval, sdpnnonz);
   }

   if(nlpcons > 0 || lpnnonz > 0)
   {
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpbegcol, nvars + 2);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplprowind, nlpcons + lpnnonz);
      BMSfreeBlockMemoryArray(sdpisolver->blkmem, &dsdplpval, nlpcons + lpnnonz);
   }

#ifdef SCIP_DEBUG
   BMSallocBlockMemory(sdpisolver->blkmem, &reason);
   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   switch ( *reason ) /* TODO: perhaps also check for feasibility and call the penalty-method here in that case */
   {
      case DSDP_CONVERGED:
         SCIPdebugMessage("DSDP converged!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_INFEASIBLE_START:
         SCIPdebugMessage("DSDP started with an infeasible point!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_SMALL_STEPS:
         SCIPdebugMessage("Short step lengths created by numerical difficulties prevented progress in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_INDEFINITE_SCHUR_MATRIX:
         SCIPdebugMessage("Schur Matrix in DSDP was indefinite but should have been positive semidefinite!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_MAX_IT:
         SCIPdebugMessage("DSDP reached maximum number of iterations!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_NUMERICAL_ERROR:
         SCIPdebugMessage("A numerical error occured in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_UPPERBOUND:
         SCIPdebugMessage("Dual objective value in DSDP reached upper bound.\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case DSDP_USER_TERMINATION:
         SCIPdebugMessage("DSDP didn't stop solving, did you?\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      case CONTINUE_ITERATING:
         SCIPdebugMessage("DSDP wants to continue iterating but somehow was stopped!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         break;

      default:
         SCIPdebugMessage("Unknown stopping reason in DSDP!\n");
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
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
SCIP_Bool SCIPsdpiSolverWasSolved(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   assert ( sdpisolver != NULL );
   return sdpisolver->solved;
}

/** returns true if the solver could determine whether or not the problem is feasible, so it returns true if the
 *  solver knows that the problem is feasible/infeasible/unbounded, it returns false if the solver doesn't know
 *  anything about the feasibility status and thus the functions IsPrimalFeasible etc. shouldn't be used */
SCIP_Bool SCIPsdpiSolverFeasibilityKnown(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiSolverGetSolFeasibility(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   assert ( primalfeasible != NULL );
   assert ( dualfeasible != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));

   switch ( *pdfeasible)
   {
      case DSDP_PDFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      case DSDP_UNBOUNDED:
         *primalfeasible = FALSE;
         *dualfeasible = TRUE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      case DSDP_INFEASIBLE:
         *primalfeasible = TRUE;
         *dualfeasible = FALSE;
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         break;

      default: /* should only include DSDP_PDUNKNOWN */
         BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
         SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible\n");
         SCIPABORT();
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverExistsPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}


/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 *  this is not implemented for all Solvers, always returns false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverHasPrimalRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
/*      SCIPerrorMessage("DSDP doesn't know if primal and dual solutions are feasible");
      SCIPABORT();
      return SCIP_ERROR;*/
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be primal feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsPrimalFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverExistsDualRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 *  this is not implemented for all Solvers, will always return false (and a debug message) if it isn't
 */
EXTERN
SCIP_Bool SCIPsdpiSolverHasDualRay(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual unbounded
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualUnbounded(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_UNBOUNDED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual infeasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualInfeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
}

/** returns TRUE iff SDP is proven to be dual feasible
 *  returns FALSE with a debug-message if the solver couldnot determine feasibility */
SCIP_Bool SCIPsdpiSolverIsDualFeasible(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPSolutionType* pdfeasible;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &pdfeasible);
   DSDP_CALL(DSDPGetSolutionType(sdpisolver->dsdp, pdfeasible));
   if (*pdfeasible == DSDP_PDUNKNOWN)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      SCIPdebugMessage("DSDP doesn't know if primal and dual solutions are feasible");
      return FALSE;
   }
   else if (*pdfeasible == DSDP_INFEASIBLE)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return FALSE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &pdfeasible);
      return TRUE;
   }
}

/** returns TRUE iff the solver converged */
SCIP_Bool SCIPsdpiSolverIsConverged(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_CONVERGED)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiSolverIsObjlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_UPPERBOUND)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiSolverIsIterlimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   if (*reason == DSDP_MAX_IT)
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return TRUE;
   }
   else
   {
      BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
      return FALSE;
   }
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiSolverIsTimelimExc(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   SCIPdebugMessage("Not implemented in DSDP!\n");
   return SCIP_ERROR;
}

/** returns the internal solution status of the solver */
int SCIPsdpiSolverGetInternalStatus(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   DSDPTerminationReason* reason;

   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   BMSallocBlockMemory(sdpisolver->blkmem, &reason);

   DSDP_CALL(DSDPStopReason(sdpisolver->dsdp, reason));

   switch ( *reason)
   {
      case DSDP_CONVERGED:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 0;
      }
      case DSDP_INFEASIBLE_START:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 1;
      }
      case DSDP_SMALL_STEPS:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_INDEFINITE_SCHUR_MATRIX:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_MAX_IT:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 4;
      }
      case DSDP_NUMERICAL_ERROR:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 2;
      }
      case DSDP_UPPERBOUND:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 3;
      }
      case DSDP_USER_TERMINATION:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 6;
      }
      default:
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &reason);
         return 7;
      }
   }
}

/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiSolverIsOptimal(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   return (SCIPsdpiIsConverged(sdpisolver) && SCIPsdpiIsPrimalFeasible(sdpisolver) && SCIPsdpiIsDualFeasible(sdpisolver));
}

/** returns TRUE iff SDP was solved to optimality or some other status was reached,
 * that is still acceptable inside a Branch & Bound framework */
SCIP_Bool SCIPsdpiSolverIsAcceptable(
   SCIP_SDPISOLVER*      sdpisolver          /**< pointer to SDP interface solver structure */
   )
{
   if (SCIPsdpiIsConverged(sdpisolver))
   {
      return TRUE;
   }
   else
   {
      double* pobj;
      double* dobj;
      double gap;

      printf("Numerical Trouble in DSDP!\n");

      /* if it didn't converge check the optimality gap */
      BMSallocBlockMemory(sdpisolver->blkmem, &pobj);
      BMSallocBlockMemory(sdpisolver->blkmem, &dobj);

      DSDP_CALL(DSDPGetPObjective(sdpisolver->dsdp, pobj));
      DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, dobj));

      gap = abs(*pobj - *dobj);

      if ((gap < epsilon) || ((gap / (0.5 * (abs(*pobj) + abs(*dobj)))) < epsilon)) /* this is the duality gap used in SDPA */
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &pobj);
         BMSfreeBlockMemory(sdpisolver->blkmem, &dobj);
         return TRUE;
      }
      else
      {
         BMSfreeBlockMemory(sdpisolver->blkmem, &pobj);
         BMSfreeBlockMemory(sdpisolver->blkmem, &dobj);
         return FALSE;
      }
   }
/* TODO: also check for primal feasibility, as this is also needed for optimality */
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiSolverIgnoreInstability(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiSolverGetObjval(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert ( sdpisolver != NULL );
   assert ( objval != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, objval));
   *objval = -1*(*objval); /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */

   return SCIP_OKAY;
}

/** gets dual solution vector for feasible SDPs */
SCIP_RETCODE SCIPsdpiSolverGetSol(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   int                   dualsollength       /**< length of the dual sol vector, must be 0 if dualsol is NULL */
   )
{
   assert ( sdpisolver != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   if ( objval != NULL )
   {
      DSDP_CALL(DSDPGetDObjective(sdpisolver->dsdp, objval));
      *objval *= -1; /*DSDP maximizes instead of minimizing, so the objective values were multiplied by -1 when inserted */
   }

   if (dualsollength > 0)
   {
      assert(dualsol != NULL);
      DSDP_CALL(DSDPGetY(sdpisolver->dsdp, dualsol, dualsollength)); /*last entry needs to be the number of variables, will return an error otherwise */
   }
   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiSolverGetIterations(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   assert ( sdpi != NULL );
   CHECK_IF_SOLVED(sdpisolver);

   DSDP_CALL(DSDPGetIts(sdpisolver->dsdp, iterations));
   return SCIP_OKAY;
}

/** gets information about the quality of an SDP solution
 *
 *  Such information is usually only available, if also a (maybe not optimal) solution is available.
 *  The SDPI should return SCIP_INVALID for *quality, if the requested quantity is not available.
 */
SCIP_RETCODE SCIPsdpiSolverGetRealSolQuality(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   SCIP_SDPSOLQUALITY    qualityindicator,   /**< indicates which quality should be returned */
   SCIP_Real*            quality             /**< pointer to store quality number */
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
SCIP_Real SCIPsdpiSolverInfinity(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   return 1E+20; /* default infinity from SCIP */
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiSolverIsInfinity(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be checked for infinity */
   )
{
   return ((val <= -SCIPsdpiSolverInfinity(sdpisolver)) || (val >= SCIPsdpiSolverInfinity(sdpisolver)));
}

/** returns highest penalty parameter to be used */
SCIP_Real SCIPsdpiSolverMaxPenParam(
   SCIP_SDPI*           sdpi                 /**< SDP interface structure */
   )
{
   return 1E+10;  /* DSDP will start with penalty param 10^10 if called normally */
}

/** checks if given value is greater or equal to the highest penalty parameter to be used */
SCIP_Bool SCIPsdpiSolverIsGEMaxPenParam(
   SCIP_SDPI*           sdpi,               /**< SDP interface structure */
   SCIP_Real            val                 /**< value to be compared to maximum penalty parameter */
   )
{
   return ((val <= -SCIPsdpiMaxPenParam(sdpisolver)) || (val >= SCIPsdpiMaxPenParam(sdpisolver)));
}

/**@} */




/*
 * File Interface Methods
 */

/**@name File Interface Methods */
/**@{ */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiSolverReadSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiSolverWriteSDP(
   SCIP_SDPISOLVER*      sdpisolver,         /**< pointer to an SDP interface solver structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("Not implemented yet\n");
   return SCIP_ERROR;
}

/**@} */
