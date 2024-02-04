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

/**@file   cons_sdp.c
 * @brief  Constraint handler for SDP-constraints
 * @author Sonja Mars
 * @author Lars Schewe
 * @author Tristan Gally
 * @author Frederic Matter
 * @author Marc Pfetsch
 *
 * Constraint handler for semidefinite constraints of the form \f$ \sum_{j=1}^n A_j y_j - A_0 \succeq 0 \f$,
 * where the matrices \f$A_j\f$ and \f$A_0\f$ need to be symmetric. Only the nonzero entries of the matrices
 * are stored.
 *
 * This file also contains a separate constraint handler for handling rank 1 SDP constraints. The callback functions are
 * essentially the same, but some quadratic constraints are added to enforce the rank 1 condition.
 *
 * Many techniques used in this constraint handler are described in the paper
 * Presolving for Mixed-Integer Semidefinite Optimization @p
 * Frederic Matter and Marc E. Pfetsch@p
 * Optimization Online.
 */

/* #define SCIP_DEBUG */
/* #define SCIP_MORE_DEBUG         /\* shows all cuts added and prints constraint after parsing *\/ */
/* #define PRINT_HUMAN_READABLE /\* change the output of PRINTCONS to a better readable format (dense instead of sparse), WHICH CAN NO LONGER BE PARSED *\/ */
/* #define PRINTMATRICES     /\* Should all matrices appearing in best rank-1 approximation heuristic be printed? *\/ */

#include "cons_sdp.h"

#include <assert.h>                     /*lint !e451*/
#include <string.h>                     /* for NULL, strcmp */
#include <ctype.h>                      /* for isspace */
#include <math.h>
#include "sdpi/lapack_interface.h"
#include "sdpi/solveonevarsdp.h"
#include "relax_sdp.h"

#include "scipsdp/SdpVarmapper.h"
#include "scipsdp/SdpVarfixer.h"

#include "scip/cons_linear.h"           /* for SCIPcreateConsLinear */
#include "scip/cons_nonlinear.h"        /* for newer SCIP versions */
#include "scip/cons_quadratic.h"        /* for SCIPcreateConsBasicQuadratic */
#include "scip/cons_soc.h"              /* for SCIPcreateConsSOC */
#include "scip/cons_linear.h"           /* for separateSol() */
#include "scip/scip_cons.h"             /* for SCIPgetConsVars */
#include "scip/scip.h"                  /* for SCIPallocBufferArray, etc */
#include "scip/def.h"
#if SCIP_VERSION >= 900
#include "scip/symmetry_graph.h"
#include "symmetry/struct_symmetry.h"
#endif

#ifdef OMP
#include "omp.h"                        /* for changing the number of threads */
#endif

/* turn off lint warnings for whole file: */
/*lint --e{788,818}*/

#define CONSHDLR_NAME          "SDP"
#define CONSHDLR_DESC          "SDP constraints of the form \\sum_{j} A_j y_j - A_0 psd"

#define CONSHDLRRANK1_NAME     "SDPrank1"
#define CONSHDLRRANK1_DESC     "rank 1 SDP constraints"

#define CONSHDLR_SEPAPRIORITY  +1000000 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY  -2000000 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -2000000 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_PROPFREQ             1 /**< priority of the constraint handler for propagation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PRESOLTIMING     SCIP_PRESOLTIMING_EXHAUSTIVE
#define CONSHDLR_PROPTIMING       SCIP_PROPTIMING_BEFORELP

#define PARSE_STARTSIZE               1 /**< initial size of the consdata-arrays when parsing a problem */
#define PARSE_SIZEFACTOR             10 /**< size of consdata-arrays is increased by this factor when parsing a problem */

#define DEFAULT_PROPUPPERBOUNDS    TRUE /**< Should upper bounds be propagated? */
#define DEFAULT_PROPUBPRESOL       TRUE /**< Should upper bounds be propagated in presolving? */
#define DEFAULT_PROP3MINORS        TRUE /**< Should 3x3 minors be propagated? */
#define DEFAULT_NONCONST3MINORS   FALSE /**< Should 3x3 minors be propagated if the diagonal is not constant? */
#define DEFAULT_PROP3MPROBING     FALSE /**< Should 3x3 minors be propagated in probing? */
#define DEFAULT_PROPTIGHTENBOUNDS  TRUE /**< Should tighten bounds be propagated? */
#define DEFAULT_PROPTBPROBING     FALSE /**< Should tighten bounds be propagated in probing? */
#define DEFAULT_TIGHTENBOUNDSCONT FALSE /**< Should only bounds be tightend for continuous variables? */
#define DEFAULT_TIGHTENMATRICES   FALSE /**< If all matrices are psd, should the matrices be tightened if possible? */
#define DEFAULT_TIGHTENBOUNDS      TRUE /**< If all matrices are psd, should the bounds be tightened if possible? */
#define DEFAULT_DIAGGEZEROCUTS    FALSE /**< Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added? */
#define DEFAULT_DIAGZEROIMPLCUTS   TRUE /**< Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added? */
#define DEFAULT_TWOMINORLINCONSS  FALSE /**< Should linear cuts corresponding to 2 by 2 minors be added? */
#define DEFAULT_TWOMINORPRODCONSS FALSE /**< Should linear cuts corresponding to products of 2 by 2 minors be added? */
#define DEFAULT_TWOMINORVARBOUNDS  TRUE /**< Should linear cuts corresponding to variable bounds for 2 by 2 minors be added? */
#define DEFAULT_QUADCONSRANK1      TRUE /**< Should quadratic cons for 2x2 minors be added in the rank-1 case? */
#define DEFAULT_UPGRADEQUADCONSS  FALSE /**< Should quadratic constraints be upgraded to a rank 1 SDP? */
#define DEFAULT_UPGRADEKEEPQUAD   FALSE /**< Should the quadratic constraints be kept in the problem after upgrading and the corresponding SDP constraint be added without the rank 1 constraint? */
#define DEFAULT_MAXNVARSQUADUPGD   1000 /**< maximal number of quadratic constraints and appearing variables so that the QUADCONSUPGD is performed */
#define DEFAULT_RANK1APPROXHEUR   FALSE /**< Should the heuristic that computes the best rank-1 approximation for a given solution be executed? */
#define DEFAULT_SEPARATEONECUT    FALSE /**< Should only one cut corresponding to the most negative eigenvalue be separated? */
#define DEFAULT_CUTSTOPOOL         TRUE /**< Should the cuts be added to the pool? */
#define DEFAULT_SPARSIFYCUT       FALSE /**< Should the eigenvector cuts be sparsified? */
#define DEFAULT_SPARSIFYFACTOR      0.1 /**< target size for sparsification in relation to number of variables */
#define DEFAULT_SPARSIFYTARGETSIZE   -1 /**< absolute target size for sparsification (-1: use sparsifyfactor instead) */
#define DEFAULT_MULTIPLESPARSECUTS FALSE /**< Should multiple sparsified eigenvector cuts be added? */
#define DEFAULT_MAXNSPARSECUTS        0 /**< maximal number of sparse eigenvector cuts that should be added (-1: no limit) */
#define DEFAULT_ENFORCESDP        FALSE /**< Solve SDP if we do lp-solving and have an integral solution in enforcing? */
#define DEFAULT_ONLYFIXEDINTSSDP  FALSE /**< Should solving an SDP only be applied if all integral variables are fixed (instead of having integral values)? */
#define DEFAULT_ADDSOCRELAX       FALSE /**< Should a relaxation of SOC constraints be added */
#define DEFAULT_USEDIMACSFEASTOL  FALSE /**< Should a feasibility tolerance based on the DIMACS be used for computing negative eigenvalues? */
#define DEFAULT_GENERATEROWS       TRUE /**< Should rows be generated (constraints otherwise)? */
#define DEFAULT_GENERATECMIR       TRUE /**< Should CMIR cuts be generated? */
#define DEFAULT_PRESOLLINCONSSPARAM   0 /**< Parameters for linear constraints added during presolving: (0) propagate, if solving LPs also separate (1) initial and propagate, if solving LPs also separate, enforce and check */
#define DEFAULT_ADDITIONALSTATS   FALSE /**< Should additional statistics be output at the end? */
#define DEFAULT_ENABLEPROPTIMING  FALSE /**< Should timing be activated for propagation routines? */
#define DEFAULT_REMOVESMALLVAL    FALSE /**< Should small values in the constraints be removed? */

#ifdef OMP
#define DEFAULT_NTHREADS              1 /**< number of threads used for OpenBLAS */
#endif

/* defines for sparsification of eigenvector cuts using TPower */
#define DEFAULT_RECOMPUTESPARSEEV FALSE /**< Should the sparse eigenvalue returned from TPower be recomputed exactly by using Lapack for the corresponding submatrix? */
#define DEFAULT_RECOMPUTEINITIAL  FALSE /**< Should the inital vector for TPower be computed each time before calling TPower (instead of using the original smallest eigenvector)? */
#define DEFAULT_EXACTTRANS        FALSE /**< Should the matrix be transformed with the exact maximal eigenvalue before calling TPower (instead of using estimate)? */

/* default values for CMIR generation */
#define BOUNDSWITCH                0.51 /**< threshold for bound switching - see cuts.c */
#define POSTPROCESS               FALSE /**< apply postprocessing to the cut - see cuts.c */
#define USEVBDS                   FALSE /**< use variable bounds - see cuts.c */
#define ALLOWLOCAL                FALSE /**< allow to generate local cuts - see cuts. */
#define MINFRAC                   0.05  /**< minimal fractionality of floor(rhs) - see cuts.c */
#define MAXFRAC                   0.999 /**< maximal fractionality of floor(rhs) - see cuts.c */

#define COEFZERO                  1e-12 /**< tolerance below which coefficients are eliminated from cuts */

/** constraint data for sdp constraints */
struct SCIP_ConsData
{
   int                   nvars;              /**< number of variables in this SDP constraint */
   int                   nnonz;              /**< number of nonzeros in this SDP constraint */
   int                   blocksize;          /**< size of this SDP-block */
   int*                  nvarnonz;           /**< length of the arrays pointed to by col/row/val, number of nonzeros for each variable */
   int**                 col;                /**< pointers to the column indices of the nonzeros for each variable */
   int**                 row;                /**< pointers to the row indices of the nonzeros for each variable */
   SCIP_Real**           val;                /**< pointers to the values of the nonzeros for each variable */
   SCIP_VAR**            vars;               /**< SCIP_VARiables present in this SDP constraint, ordered by their begvar-indices */
   int*                  locks;              /**< whether each variable is up-locked (1), down-locked (-1) or both (0); -2 if not locked (yet) */
   int                   constnnonz;         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol;           /**< column indices of the constant nonzeros */
   int*                  constrow;           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval;           /**< values of the constant nonzeros */
   SCIP_Bool*            issymunique;        /**< for each matrix marks whether it is unique and cannot be symmetric */
   SCIP_Real             maxrhsentry;        /**< maximum entry of constant matrix (needed for DIMACS error norm) */
   SCIP_Bool             rankone;            /**< Should matrix be rank one? */
   int*                  maxevsubmat;        /**< two row indices of 2x2 subdeterminant with maximal eigenvalue [or -1,-1 if not available] */
   SCIP_Bool             addedquadcons;      /**< Are the quadratic 2x2-minor constraints already added (in the rank1-case)?  */
   /* alternative view via matrix entries for propagation */
   SCIP_VAR**            matrixvar;          /**< pointer to variable if given position is uniquely covered, NULL otherwise */
   SCIP_Real*            matrixval;          /**< value at given position of unique covering variable */
   SCIP_Real*            matrixconst;        /**< value of constant matrix */
   int                   nsingle;            /**< number of matrix entries that depend on a single variable only */
   SCIP_Bool             propubpossible;     /**< whether the propagation of upper bounds is possible */
   SCIP_Bool             diagconstantone;    /**< true if all diagonal entries are fixed to be 1 (used for speeding-up propagate3minors() */
   SCIP_Real             tracebound;         /**< possible bound on the trace */
   SCIP_Bool             allmatricespsd;     /**< true if all variables are positive semidefinite (excluding the constant matrix) */
   SCIP_Bool             initallmatricespsd; /**< true if allmatricespsd has been initialized */
};

/** SDP constraint handler data */
struct SCIP_ConshdlrData
{
   int                   neigveccuts;        /**< this is used to give the eigenvector-cuts distinguishable names */
   int                   ncmir;              /**< number of cmir cuts generated */
   SCIP_Bool             diaggezerocuts;     /**< Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added? */
   int                   ndiaggezerocuts;    /**< this is used to give the diagGEzero-cuts distinguishable names */
   int                   n1x1blocks;         /**< this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */
   int                   symuniqueid;        /**< this variable is used to create unique ids within the symmetry callback */
   SCIP_Bool             propupperbounds;    /**< Should upper bounds be propagated? */
   SCIP_Bool             propubpresol;       /**< Should upper bounds be propagated in presolving? */
   SCIP_Bool             prop3minors;        /**< Should 3x3 minors be propagated? */
   SCIP_Bool             nonconst3minors;    /**< Should 3x3 minors be propagated if the diagonal is not constant? */
   SCIP_Bool             prop3mprobing;      /**< Should 3x3 minors be propagated in probing? */
   SCIP_Bool             tightenboundscont;  /**< Should only bounds be tightend for continuous variables? */
   SCIP_Bool             proptightenbounds;  /**< Should tighten bounds be propagated? */
   SCIP_Bool             proptbprobing;      /**< Should tighten bounds be propagated in probing? */
   SCIP_Bool             tightenmatrices;    /**< If all matrices are psd, should the matrices be tightened if possible? */
   SCIP_Bool             tightenbounds;      /**< If all matrices are psd, should the bounds be tightened if possible? */
   SCIP_Bool             diagzeroimplcuts;   /**< Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added? */
   SCIP_Bool             twominorlinconss;   /**< Should linear cuts corresponding to 2 by 2 minors be added? */
   SCIP_Bool             twominorprodconss;  /**< Should linear cuts corresponding to products of 2 by 2 minors be added? */
   SCIP_Bool             twominorvarbounds;  /**< Should linear cuts corresponding to variable bounds for 2 by 2 minors be added? */
   SCIP_Bool             quadconsrank1;      /**< Should quadratic cons for 2x2 minors be added in the rank-1 case? */
   SCIP_Bool             upgradequadconss;   /**< Should quadratic constraints be upgraded to a rank 1 SDP? */
   SCIP_Bool             upgradekeepquad;    /**< Should the quadratic constraints be kept in the problem after upgrading and the corresponding SDP constraint be added without the rank 1 constraint? */
   SCIP_Bool             separateonecut;     /**< Should only one cut corresponding to the most negative eigenvalue be separated? */
   SCIP_Bool             cutstopool;         /**< Should the cuts be added to the pool? */
   SCIP_Bool             sparsifycut;        /**< Should the eigenvector cuts be sparsified? */
   SCIP_Real             sparsifyfactor;     /**< target size for sparsification in relation to number of variables */
   int                   sparsifytargetsize; /**< absolute target size for sparsification (-1: use sparsifyfactor instead) */
   SCIP_Bool             multiplesparsecuts; /**< Should multiple sparsified eigenvector cuts be added? */
   int                   maxnsparsecuts;     /**< maximal number of sparse eigenvector cuts that should be added (-1: no limit) */
   SCIP_Bool             enforcesdp;         /**< Solve SDP if we do lp-solving and have an integral solution in enforcing? */
   SCIP_Bool             onlyfixedintssdp;   /**< Should solving an SDP only be applied if all integral variables are fixed (instead of having integral values)? */
   SCIP_Bool             addsocrelax;        /**< Should a relaxation of SOC constraints be added */
   int                   maxnvarsquadupgd;   /**< maximal number of quadratic constraints and appearing variables so that the QUADCONSUPGD is performed */
   SCIP_Bool             triedlinearconss;   /**< Have we tried to add linear constraints? */
   SCIP_Bool             triedvarbounds;     /**< Have we tried to add variable bounds based on 2x2 minors */
   SCIP_Bool             rank1approxheur;    /**< Should the heuristic that computes the best rank-1 approximation for a given solution be executed? */
   SCIP_Bool             generaterows;       /**< Should rows be generated (constraints otherwise)? */
   SCIP_Bool             generatecmir;       /**< Should CMIR cuts be generated? */
#ifdef OMP
   int                   nthreads;           /**< number of threads used for OpenBLAS */
#endif
   int*                  quadconsidx;        /**< store index of variables appearing in quadratic constraints for upgrading */
   SCIP_VAR**            quadconsvars;       /**< temporary array to store variables appearing in quadratic constraints for upgrading */
   int                   nquadconsidx;       /**< size of quadconsidx/quadconsvars arrays */
   SCIP_VAR***           X;                  /**< matrix variables added within upgrading */
   int                   nsdpvars;           /**< number of variables in SDP constraint for quadratic constraints */
   SCIP_CONS*            sdpcons;            /**< SDP rank 1 constraint for quadratic constraints */
   SCIP_CONSHDLRDATA*    sdpconshdlrdata;    /**< possibly store SDP constraint handler for retrieving parameters */
   SCIP_RANDNUMGEN*      randnumgen;         /**< random number generator (for sparsifyCut) */
   SCIP_RELAX*           relaxsdp;           /**< SDP relaxator */
   SCIP_Bool             usedimacsfeastol;   /**< Should a feasibility tolerance based on the DIMACS be used for computing negative eigenvalues? */
   SCIP_Real             dimacsfeastol;      /**< feasibility tolerance for computing negative eigenvalues based on the DIMACS error */
   SCIP_Bool             recomputesparseev;  /**< Should the sparse eigenvalue returned from TPower be recomputed exactly by using Lapack for the corresponding submatrix? */
   SCIP_Bool             recomputeinitial;   /**< Should the inital vector for TPower be computed each time before calling TPower (instead of using the original smallest eigenvector)? */
   SCIP_Bool             exacttrans;         /**< Should the matrix be transformed with the exact maximal eigenvalue before calling TPower (instead of using estimate)? */
   int                   presollinconssparam; /**< Parameters for linear constraints added during presolving: (0) propagate, if solving LPs also separate (1) initial and propagate, if solving LPs also separate, enforce and check */
   SCIP_Bool             additionalstats;    /**< Should additional statistics be output at the end? */
   SCIP_Bool             enableproptiming;   /**< Should timing be activated for propagation routines? */
   SCIP_Bool             removesmallval;     /**< Should small values in the constraints be removed? */

   int                   ncallspropub;       /**< Number of calls of propagateUpperBounds in propagation */
   int                   ncallsproptb;       /**< Number of calls of tightenBounds in propagation */
   int                   ncallsprop3minor;   /**< Number of calls of propagate3Minors in propagation */
   SCIP_CLOCK*           propubtime;         /**< Time for propagateUpperBounds in propagation */
   SCIP_CLOCK*           proptbtime;         /**< Time for tightenBounds in propagation */
   SCIP_CLOCK*           prop3minortime;     /**< Time for propagate3Minors in propagation */
   SCIP_Real             maxtimepropub;      /**< Maximal time spent for one round of propagateUpperBounds in propagation */
   SCIP_Real             maxtimeproptb;      /**< Maximal time spent for one round of tightenBounds in propagation */
   SCIP_Real             maxtimeprop3minor;  /**< Maximal time spent for one round of propagate3Minors in propagation */
   int                   npropub;            /**< Number of propagations through upper bounds */
   int                   nproptb;            /**< Number of tightened bounds in propagation */
   int                   nprop3minor;        /**< Number of propagations through 3x3 minors */
   int                   npropcutoffub;      /**< Number of cutoffs in propagation through upper bounds */
   int                   npropcutofftb;      /**< Number of cutoffs in propagation of bound tightening */
   int                   npropcutoff3m;      /**< Number of cutoffs in propagation through 3x3 minors */
   int                   npropintrndub;      /**< Number of rounded bounds of integer variables in propagation through upper bounds */
   int                   npropintrndtb;      /**< Number of rounded bounds of integer variables in propagation of bound tightening */
   int                   nproppreub;         /**< Number of propagations through upper bounds in presolving */
   int                   nproppretb;         /**< Number of tightened bounds in propagation in presolving */
   int                   nproppre3m;         /**< Number of propagations through 3x3 minors in presolving */
   int                   nproppreintrndub;   /**< Number of rounded bounds of integer variables in propagation through upper bounds in presolving */
   int                   nproppreintrndtb;   /**< Number of rounded bounds of integer variables in propagation of bound tightening in presolving */
   int                   nproppreintrnd3m;   /**< Number of rounded bounds of integer variables in propagation through 3x3 minors in presolving */
   int                   npropprobub;        /**< Number of propagations through upper bounds in probing */
   int                   npropprobtb;        /**< Number of tightened bounds in propagation in probing */
   int                   npropprob3minor;    /**< Number of propagations through 3x3 minor in probing */
};

/** generates matrix in colum-first format (needed by LAPACK) from matrix given in full row-first format (SCIP-SDP
 *  default)
 */
static
SCIP_RETCODE convertRowToColFormatFullMatrix(
   int                   rows,               /**< number of rows */
   int                   cols,               /**< number of columns */
   SCIP_Real*            rowmatrix,          /**< matrix entries given as rows*cols array */
   SCIP_Real*            colmatrix           /**< pointer to array of length blocksize^2 to store rowmatrix in column-first
                                              *   format */
   )
{
   int i;
   int j;
   SCIP_Real act;

   assert( rows > 0 );
   assert( cols > 0 );
   assert( rowmatrix != NULL );
   assert( colmatrix != NULL );

   for (i = 0; i < rows; ++i)
   {
      for (j = 0; j < cols; ++j)
      {
         act = rowmatrix[i*cols + j];
         colmatrix[j*rows + i] = act;
      }
   }

   return SCIP_OKAY;
}

/** multiplies all entries in the i-th row by scale[i] */
static
SCIP_RETCODE scaleRowsMatrix(
   int                   blocksize,          /* number of rows and columns */
   SCIP_Real*            matrix,             /* matrix entries given as blocksize^2 array */
   SCIP_Real*            scale               /* array of length blocksize to multiply the rows of matrix with */
   )
{
   int r;                       /* rows */
   int c;                       /* columns */

   assert( blocksize >= 0 );
   assert( matrix != NULL );
   assert( scale != NULL );

   for (r = 0; r < blocksize; r++)
   {
      for (c = 0; c < blocksize; c++)
      {
         /* row-first format! */
         matrix[r * blocksize + c] *= scale[r];  /*lint !e679*/
      }
   }

   return SCIP_OKAY;
}

/** takes a 0.5*n*(n+1) array of a symmetric matrix and expands it to an n*n array of the full matrix to input into LAPACK */
static
SCIP_RETCODE expandSymMatrix(
   int                   size,               /**< size of the matrix, named n above */
   SCIP_Real*            symMat,             /**< symmetric matrix indexed via SCIPconsSdpCompLowerTriangPos that should be expanded */
   SCIP_Real*            fullMat             /**< pointer to store the n*n matrix, that is the symmetric expansion of symMat */
   )
{
   int i;
   int j;
   int ind = 0;

   assert( size >= 0 );
   assert( symMat != NULL );
   assert( fullMat != NULL );

   /* traverse the lower triangular part in the order of the indices and copy the values to both lower and upper triangular part */
   for (i = 0; i < size; i++)
   {
      for (j = 0; j <= i; j++)
      {
         assert( ind == SCIPconsSdpCompLowerTriangPos(i,j) );
         fullMat[i*size + j] = symMat[ind]; /*lint !e679*/
         fullMat[j*size + i] = symMat[ind]; /*lint !e679*/
         ind++;
      }
   }

   return SCIP_OKAY;
}

/** For a vector \f$y\f$ given by @p sol, computes the (length of y) * (length of y + 1) /2 -long array of the lower-triangular part
 *  of the matrix \f$ \sum_{j=1}^m A_j y_j - A_0 \f$ for this SDP block, indexed by SCIPconsSdpCompLowerTriangPos().
 */
static
SCIP_RETCODE computeSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real*            matrix              /**< pointer to store the SDP-Matrix */
   )
{
   SCIP_Real yval;
   int blocksize;
   int nvars;
   int i;
   int j;

   assert( consdata != NULL );
   assert( matrix != NULL );

   nvars = consdata->nvars;
   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1))/2; ++i)
      matrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      yval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
      if ( ! SCIPisZero(scip, yval) )
      {
         for (j = 0; j < consdata->nvarnonz[i]; ++j)
            matrix[SCIPconsSdpCompLowerTriangPos(consdata->row[i][j], consdata->col[i][j])] += yval * consdata->val[i][j];
      }
   }

   /* substract the constant part */
   for (j = 0; j < consdata->constnnonz; ++j)
      matrix[SCIPconsSdpCompLowerTriangPos(consdata->constrow[j], consdata->constcol[j])] -= consdata->constval[j];

   return SCIP_OKAY;
}

/** for a vector \f$y\f$ given by @p sol, computes the full matrix \f$ \sum_{j=1}^m A_j y_j - A_0 \f$ */
static
SCIP_RETCODE computeFullSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_SOL*             sol,                /**< solution to separate */
   SCIP_Real*            fullmatrix          /**< array for full matrix */
   )
{
   SCIP_Real aval;
   int blocksize;
   int nvars;
   int r;
   int c;
   int i;
   int j;

   assert( scip != NULL );
   assert( consdata != NULL );
   assert( fullmatrix != NULL );

   nvars = consdata->nvars;
   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < blocksize * blocksize; ++i)
      fullmatrix[i] = 0.0;

   /* add the non-constant-part */
   for (i = 0; i < nvars; i++)
   {
      SCIP_Real yval;

      yval = SCIPgetSolVal(scip, sol, consdata->vars[i]);
      if ( ! SCIPisZero(scip, yval) )
      {
         for (j = 0; j < consdata->nvarnonz[i]; ++j)
         {
            r = consdata->row[i][j];
            c = consdata->col[i][j];
            aval = consdata->val[i][j];

            if ( r == c )
               fullmatrix[r * blocksize + c] += yval * aval;
            else
            {
               fullmatrix[r * blocksize + c] += yval * aval;
               fullmatrix[c * blocksize + r] += yval * aval;
            }
         }
      }
   }

   /* substract the constant part */
   for (j = 0; j < consdata->constnnonz; ++j)
   {
      r = consdata->constrow[j];
      c = consdata->constcol[j];
      aval = consdata->constval[j];

      if ( r == c )
         fullmatrix[r * blocksize + c] -= aval;
      else
      {
         fullmatrix[r * blocksize + c] -= aval;
         fullmatrix[c * blocksize + r] -= aval;
      }
   }

   return SCIP_OKAY;
}

/** check whether propagation of upper bounds can be applied */
static
SCIP_RETCODE checkPropagateUpperbounds(
   SCIP_CONS*            cons                /**< constraints to check */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   consdata->propubpossible = FALSE;

   if ( consdata->nsingle > 0 )
   {
      int blocksize;
      int s;
      int t;

      blocksize = consdata->blocksize;
      assert( consdata->matrixval != NULL );
      assert( consdata->matrixvar != NULL );

      /* check all off-diagonal positions */
      for (s = 1; s < blocksize; ++s)
      {
         int diags;

         diags = s * (s + 1)/2 + s;
         if ( consdata->matrixval[diags] == SCIP_INVALID ) /*lint !e777*/
            continue;

         for (t = 0; t < s; ++t)
         {
            SCIP_VAR* varst;
            int diagt;
            int pos;

            pos = s * (s + 1)/2 + t;
            varst = consdata->matrixvar[pos];
            if ( varst == NULL || ! SCIPvarIsActive(varst) )
               continue;

            diagt = t * (t + 1)/2 + t;
            if ( consdata->matrixval[diagt] == SCIP_INVALID ) /*lint !e777*/
               continue;

            /* at this place propagation of upper bounds would be possible */
            consdata->propubpossible = TRUE;
            break;
         }

         if ( consdata->propubpossible )
            break;
      }
   }

   return SCIP_OKAY;
}

/** build matrixvar data
 *
 *  We have:
 *  - matrixvar[i] is NULL if position i is not uniquely covered by a variable (either because there is no variable or at least two variables that cover the position).
 *  - matrixval[i] == 0.0 if position i is not covered by any variable.
 *  - matrixval[i] == SCIP_INVALID if position i is covered by at least two variables.
 */
static
SCIP_RETCODE constructMatrixvar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata            /**< constraint data */
   )
{
   SCIP_Real* constmatrix;
   SCIP_Real** matrices;
   int blocksize;
   int cnt = 0;
   int s;
   int t;
   int i;

   assert( scip != NULL );
   assert( consdata != NULL );

   if ( consdata->matrixvar != NULL )
      return SCIP_OKAY;

   consdata->nsingle = 0;
   blocksize = consdata->blocksize;

   /* allocate matrices */
   SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
   SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, constmatrix) );

   SCIP_CALL( SCIPallocBufferArray(scip, &matrices, consdata->nvars) );
   for (i = 0; i < consdata->nvars; ++i)
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &matrices[i], blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, i, matrices[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->matrixvar, blocksize * (blocksize+1)/2) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->matrixval, blocksize * (blocksize+1)/2) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->matrixconst, blocksize * (blocksize+1)/2) );

   consdata->diagconstantone = TRUE;
   for (s = 0; s < blocksize; ++s)
   {
      for (t = 0; t <= s; ++t)
      {
         SCIP_VAR* var = NULL;
         SCIP_Real val = 0.0;
         int pos;

         pos = s * blocksize + t;

         for (i = 0; i < consdata->nvars; ++i)
         {
            if ( ! SCIPisZero(scip, matrices[i][pos]) )
            {
               if ( var == NULL )
               {
                  var = consdata->vars[i];
                  val = matrices[i][pos];
               }
               else
                  break;
            }
         }

         /* if at most one entry was found */
         if ( i >= consdata->nvars )
         {
            consdata->matrixvar[cnt] = var;  /* note that var == NULL is possible */
            consdata->matrixval[cnt] = val;
            consdata->matrixconst[cnt] = constmatrix[pos];
            ++consdata->nsingle;
            if ( s == t && (var != NULL || ! SCIPisZero(scip, val) || ! SCIPisEQ(scip, constmatrix[pos], -1.0)) )
               consdata->diagconstantone = FALSE;
         }
         else
         {
            consdata->matrixvar[cnt] = NULL;
            consdata->matrixval[cnt] = SCIP_INVALID;
            consdata->matrixconst[cnt] = SCIP_INVALID;
            if ( s == t )
               consdata->diagconstantone = FALSE;
         }
         ++cnt;
      }
   }
   assert( cnt == blocksize * (blocksize + 1)/2 );

   SCIPfreeBufferArray(scip, &constmatrix);
   for (i = consdata->nvars - 1; i >= 0; --i)
      SCIPfreeBufferArray(scip, &matrices[i]);
   SCIPfreeBufferArray(scip, &matrices);

   if ( SCIPgetSubscipDepth(scip) == 0 )
      SCIPdebugMsg(scip, "Total number of matrix entries that only depend on a single variable: %d.\n", consdata->nsingle);

   /* determine whether propagation of upper bounds is possible */
   SCIP_CALL( checkPropagateUpperbounds(cons) );

   return SCIP_OKAY;
}

/** checks feasibility for a single SDP constraint */
static
SCIP_RETCODE SCIPconsSdpCheckSdpCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< the constraint the solution should be checked for */
   SCIP_SOL*             sol,                /**< the solution to check feasibility for */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_RESULT*          result              /**< pointer to store the result of the feasibility checking call */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real eigenvalue;
   SCIP_Real tol;
   int blocksize;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 );
   assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLRRANK1_NAME) == 0 );
   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize) ); /*lint !e647*/

   SCIP_CALL( computeFullSdpMatrix(scip, consdata, sol, fullmatrix) );

   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, fullmatrix, 1, &eigenvalue, NULL) );

   if ( conshdlrdata->sdpconshdlrdata->usedimacsfeastol )
   {
      assert( conshdlrdata->dimacsfeastol != SCIP_INVALID );
      tol = conshdlrdata->dimacsfeastol;
   }
   else
      tol = SCIPfeastol(scip);

   if ( eigenvalue >= -tol )
      *result = SCIP_FEASIBLE;
   else
   {
      *result = SCIP_INFEASIBLE;
      if ( printreason )
      {
         SCIPinfoMessage(scip, NULL, "SDP-constraint <%s> violated: non psd matrix (eigenvalue %f).\n", SCIPconsGetName(cons), eigenvalue);
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      }
   }

   if ( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, -eigenvalue, (-eigenvalue) / (1.0 + consdata->maxrhsentry));

   SCIPfreeBufferArray(scip, &fullmatrix);

   return SCIP_OKAY;
}

/** Check whether current matrix is rank one, if not so, sets maxevsubmat */
static
SCIP_RETCODE isMatrixRankOne(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the SDP constraint to check the rank for */
   SCIP_SOL*             sol,                /**< solution to check for rank one */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool*            isrankone           /**< pointer to return whether matrix is rank one */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* matrix = NULL;
   SCIP_Real* fullmatrix = NULL;
   SCIP_Real eigenvalue;
   int blocksize;
   int i;
   int j;
   int ind1 = 0;
   int ind2 = 0;
   SCIP_Real submatrix[4];
   SCIP_Real largestminev = 0.0;

   assert( cons != NULL );
   assert( isrankone != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;
   *isrankone = TRUE;

   /* allocate memory to store full matrix */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize+1))/2 ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ) );

   /* compute the matrix \f$ \sum_j A_j y_j - A_0 \f$ - we need an undestroyed version in matrix below */
   SCIP_CALL( computeSdpMatrix(scip, consdata, sol, matrix) );

   /* expand it because LAPACK wants the full matrix instead of the lower triangular part */
   SCIP_CALL( expandSymMatrix(blocksize, matrix, fullmatrix) );

   /* compute the second largest eigenvalue */
   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, fullmatrix, blocksize - 1, &eigenvalue, NULL) );

   /* the matrix is rank 1 iff the second largest eigenvalue is zero (since the matrix is symmetric and psd) */
   if ( SCIPisFeasEQ(scip, eigenvalue, 0.0) )
      *isrankone = TRUE;
   else
   {
      *isrankone = FALSE;
      if ( printreason )
      {
         SCIPinfoMessage(scip, NULL, "SDPrank1-constraint <%s> is not rank1 (second largest eigenvalue %f).\n", SCIPconsGetName(cons), eigenvalue);
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
      }

      /* if the matrix is not rank 1, compute minimal eigenvalues of 2x2 minors and save the largest of these minimal eigenvalues */
      for (i = 0; i < blocksize; ++i)
      {
         for (j = 0; j < i; ++j)
         {
            submatrix[0] = matrix[SCIPconsSdpCompLowerTriangPos(i,i)];
            submatrix[1] = matrix[SCIPconsSdpCompLowerTriangPos(i,j)];
            submatrix[2] = matrix[SCIPconsSdpCompLowerTriangPos(i,j)];
            submatrix[3] = matrix[SCIPconsSdpCompLowerTriangPos(j,j)];

            SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, 2, submatrix, 1, &eigenvalue, NULL) );

            if ( eigenvalue > largestminev )
            {
               largestminev = eigenvalue;
               ind1 = i;
               ind2 = j;
            }
         }
      }

      assert( ind1 > 0 || ind2 > 0 );
      assert( SCIPisFeasPositive(scip, largestminev) ); /* because the second largest eigenvalue satisfies SCIPisFeasPositive. */

      /* save indices for submatrix with largest minimal eigenvalue */
      consdata->maxevsubmat[0] = ind1;
      consdata->maxevsubmat[1] = ind2;
   }

   if ( sol != NULL )
      SCIPupdateSolConsViolation(scip, sol, largestminev, (largestminev) / (1.0 + consdata->maxrhsentry));

   SCIPfreeBufferArray(scip, &fullmatrix);
   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}

/** For a given variable-index j and a Vector v computes \f$ v^T A_j v \f$. */
static
SCIP_RETCODE multiplyConstraintMatrix(
   SCIP_CONS*            cons,               /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   int                   j,                  /**< variable-index of the matrix to multiply with */
   SCIP_Real*            v,                  /**< vector to multiply with */
   SCIP_Real*            vAv                 /**< pointer to store the the resulting scalar \f$ v^T A_j v \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real s = 0.0;
   int r;
   int c;
   int i;

   assert( cons != NULL );
   assert( j >= 0 );
   assert( vAv != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   assert( j < consdata->nvars );

   for (i = 0; i < consdata->nvarnonz[j]; i++)
   {
      r = consdata->row[j][i];
      c = consdata->col[j][i];
      if ( r == c )
         s += v[c] * consdata->val[j][i] * v[r];
      else
      {
         /* Multiply by 2, because the matrix is symmetric and there is one identical contribution each from lower and upper triangular part. */
         s += 2.0 * v[c] * consdata->val[j][i] * v[r];
      }
   }

   *vAv = s;

   return SCIP_OKAY;
}

/** Set the maximum absolute value of an entry of the constant matrix.
 *  This must be done before presolving, because otherwise this is influenced by variable fixings (which might lead
 *  to solutions being feasible in presolving no longer being feasible afterwards) */
static
SCIP_RETCODE setMaxRhsEntry(
   SCIP_CONS*            cons                /**< the SDP constraint that includes the Matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real max = 0.0;    /* initialize max with zero (this is used if there is no constant-matrix) */
   int i;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* iterate over the entries of the constant matrix, updating max if a higher absolute value is found */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      if ( REALABS(consdata->constval[i]) > max )
         max = REALABS(consdata->constval[i]);
   }

   consdata->maxrhsentry = max;

   return SCIP_OKAY;
}

/** produce cut from (possibly modified) eigenvector */
static
SCIP_RETCODE produceCutFromEigenvector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             enforce,            /**< whether we are in enforcing */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   blocksize,          /**< size of block */
   SCIP_Real*            fullconstmatrix,    /**< precomputed full constant matrix */
   SCIP_Real*            eigenvector,        /**< original eigenvector */
   SCIP_Real*            vector,             /**< temporary workspace (length blocksize) */
   SCIP_VAR**            vars,               /**< temporary workspace (length nvars) */
   SCIP_Real*            vals,               /**< temporary workspace (length nvars) */
   int*                  ngen,               /**< pointer to store the number of generated cuts/constraints */
   SCIP_Bool*            success,            /**< pointer to store whether we have produced a cut/constraint */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_Real lhs = 0.0;
   int cnt = 0;
   int j;

   assert( conshdlr != NULL );
   assert( conshdlrdata != NULL );
   assert( cons != NULL );
   assert( consdata != NULL );
   assert( fullconstmatrix != NULL );
   assert( eigenvector != NULL );
   assert( vector != NULL );
   assert( vars != NULL );
   assert( vals != NULL );
   assert( ngen != NULL );
   assert( success != NULL );
   assert( result != NULL );

   *success = TRUE;

   /* multiply eigenvector with constant matrix to get lhs (after multiplying again with eigenvector from the left) */
   SCIP_CALL( SCIPlapackMatrixVectorMult(blocksize, blocksize, fullconstmatrix, eigenvector, vector) );

   for (j = 0; j < blocksize; ++j)
      lhs += eigenvector[j] * vector[j];

   /* compute \f$ v^T A_j v \f$ for eigenvector v and each matrix \f$ A_j \f$ to get the coefficients of the LP cut */
   for (j = 0; j < consdata->nvars; ++j)
   {
      SCIP_Bool isfixed;
      SCIP_Real coef;
      SCIP_Real lb;
      SCIP_Real ub;

      /* compute coefficient by multiplying eigenvector with jth matrix */
      SCIP_CALL( multiplyConstraintMatrix(cons, j, eigenvector, &coef) );

      lb = SCIPvarGetLbGlobal(consdata->vars[j]);
      ub = SCIPvarGetUbGlobal(consdata->vars[j]);

      if ( ! ( SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub) ) && SCIPisEQ(scip, ub, lb) )
         isfixed = TRUE;
      else
         isfixed = FALSE;

      /* safely round coefficients to 0 */
      if ( ! enforce && (isfixed || SCIPisZero(scip, coef)) )
      {
         if ( REALABS(coef) > COEFZERO )
         {
            if ( coef > 0.0 )
            {
               if ( SCIPisInfinity(scip, ub) )
               {
                  *success = FALSE;
                  break;
               }
               else
                  lhs -= coef * ub;
            }
            else
            {
               if ( SCIPisInfinity(scip, -lb) )
               {
                  *success = FALSE;
                  break;
               }
               else
                  lhs -= coef *lb;
            }
         }
      }
      else
      {
         vars[cnt] = consdata->vars[j];
         vals[cnt++] = coef;
      }
   }

   if ( *success && cnt > 0 )
   {
      (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_%d", ++(conshdlrdata->neigveccuts));
      if ( conshdlrdata->sdpconshdlrdata->generaterows )
      {
         SCIP_Bool infeasible;
         SCIP_ROW* row;

         SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, cutname, lhs, SCIPinfinity(scip), FALSE, FALSE, TRUE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, row, cnt, vars, vals) );

         /* if we are enforcing, we take any of the cuts, otherwise only efficacious cuts */
         if ( enforce || SCIPisCutEfficacious(scip, sol, row) )
         {
#ifdef SCIP_MORE_DEBUG
            SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            if ( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SEPARATED;
            SCIP_CALL( SCIPresetConsAge(scip, cons) );

            if ( conshdlrdata->sdpconshdlrdata->cutstopool )
            {
               SCIP_CALL( SCIPaddPoolCut(scip, row) );
            }
            ++(*ngen);
         }

         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
      else
      {
         SCIP_CONS* newcons;

         SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, cutname, cnt, vars, vals, lhs, SCIPinfinity(scip),
               TRUE, TRUE, enforce, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, TRUE) );
         SCIP_CALL( SCIPaddCons(scip, newcons) );
         SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
         *result = SCIP_CONSADDED;
         ++(*ngen);
      }

      if ( conshdlrdata->sdpconshdlrdata->generatecmir )
      {
         SCIP_AGGRROW* aggrrow;
         SCIP_Real* cutcoefs;
         SCIP_Real cutrhs;
         SCIP_Real cutefficacy;
         SCIP_Bool cutislocal;
         int* cutinds;
         int cutnnz = 0;
         int cutrank;

         SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, consdata->nvars) );

         SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );

         /* switch row */
         for (j = 0; j < cnt; ++j)
         {
            vals[j] *= -1.0;
            cutinds[j] = SCIPvarGetProbindex(vars[j]);
         }

         /* add produced row as custom row: multiply with -1.0 to convert >= into <= row */
         SCIP_CALL( SCIPaggrRowAddCustomCons(scip, aggrrow, cutinds, vals, cnt, -lhs, 1.0, 0, FALSE) );

         /* try to generate CMIR inequality */
         cutefficacy = - SCIPinfinity(scip);
         SCIP_CALL( SCIPcutGenerationHeuristicCMIR(scip, NULL, POSTPROCESS, BOUNDSWITCH, USEVBDS, ALLOWLOCAL, INT_MAX, NULL, NULL,
               MINFRAC, MAXFRAC, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank, &cutislocal, success) );

         if ( *success && SCIPisEfficacious(scip, cutefficacy) )
         {
            SCIP_Bool infeasible;
            SCIP_ROW* row;
            SCIP_VAR** scipvars;

            /* set name and create empty row */
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "sepa_eig_sdp_CMIR_%d", ++(conshdlrdata->ncmir));
            SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &row, conshdlr, cutname, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE, TRUE) );

            /* get SCIP vars, because inds refers to these indices */
            scipvars = SCIPgetVars(scip);

            /* cache the row extension and only flush them if the cut gets added */
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            /* collect all non-zero coefficients */
            for (j = 0; j < cutnnz; ++j)
            {
               SCIP_CALL( SCIPaddVarToRow(scip, row, scipvars[cutinds[j]], cutcoefs[j]) );
            }

            /* flush all changes before adding the cut */
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            if ( SCIProwGetNNonz(row) == 0 )
            {
               assert( SCIPisFeasNegative(scip, cutrhs) );
               *result = SCIP_CUTOFF;
            }
            else
            {
#ifdef SCIP_MORE_DEBUG
               SCIP_CALL( SCIPprintRow(scip, row, NULL) );
#endif
               SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
               if ( infeasible )
                  *result = SCIP_CUTOFF;
               else
                  *result = SCIP_SEPARATED;
               SCIP_CALL( SCIPresetConsAge(scip, cons) );

               if ( conshdlrdata->sdpconshdlrdata->cutstopool )
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, row) );
               }
               ++(*ngen);
            }

            /* release the row */
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
         }

         SCIPfreeBufferArray(scip, &cutinds);
         SCIPfreeBufferArray(scip, &cutcoefs);

         SCIPaggrRowFree(scip, &aggrrow);
      }
   }

   return SCIP_OKAY;
}

/** compute largest sparse eigenvalue for a given sparsity level and corresponding eigenvector of a given matrix
 *
 *  The truncated power method works like the ordinary power method to compute the largest eigenvalue of a matrix, but
 *  truncates the iterates in each step to the k largest entries in absolute value, if k is the sparsity level. See
 *  [Yuan, Zhang: Truncated Power Method for Sparse Eigenvalue Problems].
 */
static
SCIP_RETCODE truncatedPowerMethod(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   blocksize,          /**< size of matrix */
   SCIP_Real*            fullmatrix,         /**< matrix for which the smallest sparse eigenvalue and eigenvector should
                                              *   be computed, need to be given with all blocksze * blocksize entries */
   SCIP_Real*            vector,             /**< initial vector for starting the truncated power method */
   int                   sparsity,           /**< sparsity level of eigenvalue and eigenvector */
   SCIP_Real             convergencetol,     /**< tolerance to be used to detect convergence */
   SCIP_Real*            eigenvector,        /**< pointer to store computed eigenvector */
   int*                  support,            /**< pointer to store support of the computed eigenvector */
   SCIP_Real*            eigenvalue          /**< pointer to store computed eigenvalue */
   )
{
   SCIP_Real* sortedev;
   SCIP_Real* oldev;
   SCIP_Real oldeig;
   SCIP_Real neweig;
   SCIP_Real norm;
   int* indices;
   int i;

   assert( scip != NULL );
   assert( blocksize > 0 );
   assert( fullmatrix != NULL );
   assert( vector != NULL );
   assert( sparsity > 0 && sparsity <= blocksize );
   assert( eigenvector != NULL );
   assert( eigenvalue != NULL );

   SCIPdebugMsg(scip, "Executing truncatedPowerMethod\n");

   SCIP_CALL( SCIPallocBufferArray(scip, &indices, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sortedev, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldev, blocksize) );

   neweig = -1.0;
   oldeig = -2.0;

   for (i = 0; i < blocksize; i++)
      eigenvector[i] = vector[i];

   while ( neweig - oldeig > convergencetol )
   {
      /* reset eigenvector and eigenvalue */
      oldeig = neweig;
      neweig = 0.0;

      for (i = 0; i < blocksize; i++)
         oldev[i] = eigenvector[i];

      SCIP_CALL( SCIPlapackMatrixVectorMult(blocksize, blocksize, fullmatrix, oldev, eigenvector) );

      /* get indices of largest absolute values */
      for (i = 0; i < blocksize; i++)
      {
         indices[i] = i;
         sortedev[i] = REALABS(eigenvector[i]);
      }

      SCIPsortDownRealInt(sortedev, indices, blocksize);

      /* truncate iterate to entries with largest absolute value and save support */
      for (i = sparsity; i < blocksize; i++)
      {
         eigenvector[indices[i]] = 0.0;
         support[indices[i]] = 0;
      }

      /* normalize iterate */
      norm = 0.0;
      for (i = 0; i < sparsity; i++)
      {
         norm += eigenvector[indices[i]] * eigenvector[indices[i]];
         support[indices[i]] = 1;
      }
      norm = sqrt(norm);

      for (i = 0; i < sparsity; i++)
         eigenvector[indices[i]] /= norm;

      /* compute iterate^T * fullmatrix * iterate */
      SCIP_CALL( SCIPlapackMatrixVectorMult(blocksize, blocksize, fullmatrix, eigenvector, oldev) );

      for (i = 0; i < blocksize; i++)
         neweig += eigenvector[i] * oldev[i];
   }

   SCIPfreeBufferArray(scip, &oldev);
   SCIPfreeBufferArray(scip, &sortedev);
   SCIPfreeBufferArray(scip, &indices);

   *eigenvalue = neweig;

   return SCIP_OKAY;
}

/** try to sparsify cut
 *
 *  We currently take a small subset of the components of a given eigenvector and check whether the cut is
 *  violated. Note that this does not necessarily yield a sparse cut because the resulting vector has to be multiplied
 *  with the matrix pencil.
 */
static
SCIP_RETCODE sparsifyCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             enforce,            /**< whether we are in enforcing */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   blocksize,          /**< size of block */
   SCIP_Real*            fullconstmatrix,    /**< precomputed full constant matrix */
   SCIP_Real*            eigenvector,        /**< original eigenvector */
   SCIP_Real*            vector,             /**< temporary workspace (length blocksize) */
   SCIP_VAR**            vars,               /**< temporary workspace */
   SCIP_Real*            vals,               /**< temporary workspace */
   int*                  ngen,               /**< pointer to store the number of generated cuts */
   SCIP_Bool*            success,            /**< pointer to store whether we have produced a cut/constraint */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_Real* ev;
   SCIP_Real norm = 0.0;
   int size;
   int* idx;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( fullconstmatrix != NULL );
   assert( eigenvector != NULL );
   assert( vector != NULL );
   assert( success != NULL );
   assert( result != NULL );
   assert( *result != SCIP_CUTOFF );

   *success = FALSE;

   /* produce random indices */
   SCIP_CALL( SCIPallocBufferArray(scip, &idx, blocksize) );
   for (j = 0; j < blocksize; ++j)
      idx[j] = j;

   /* the following should only be allocated once */
   assert( conshdlrdata->randnumgen != NULL );
   SCIPrandomPermuteIntArray(conshdlrdata->randnumgen, idx, 0, blocksize);

   /* compute target size, use sparsifytargetsize if specified, and sparsifyfactor else */
   if ( conshdlrdata->sdpconshdlrdata->sparsifytargetsize > 0 )
      size = conshdlrdata->sdpconshdlrdata->sparsifytargetsize;
   else
      size = MAX(10, (int) conshdlrdata->sdpconshdlrdata->sparsifyfactor * consdata->nvars);

   /* if size is larger than blocksize, trigger a debug message and set size to blocksize */
   if ( size > blocksize )
   {
      SCIPdebugMsg(scip, "Eigenvector cut is not sparsified since desired sparsity %d is larger than SDP blocksize %d.\n", size, blocksize);
      size = blocksize;
   }

   assert( size <= blocksize );

   /* take random subset of eigenvector - the remaining entries are 0 */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &ev, blocksize) );
   for (j = 0; j < size; ++j)
   {
      ev[idx[j]] = eigenvector[j];
      norm += eigenvector[j] * eigenvector[j];
   }

   /* if we by chance selected a zero subvector */
   if ( SCIPisFeasZero(scip, norm) )
   {
      SCIPfreeBufferArray(scip, &ev);
      SCIPfreeBufferArray(scip, &idx);
      return SCIP_OKAY;
   }

   /* normalize */
   norm = sqrt(norm);
   for (j = 0; j < size; ++j)
      ev[idx[j]] /= norm;

   /* produce cut/constraint */
   SCIP_CALL( produceCutFromEigenvector(scip, conshdlr, conshdlrdata, cons, consdata, enforce, sol,
         blocksize, fullconstmatrix, ev, vector, vars, vals, ngen, success, result) );

   SCIPfreeBufferArray(scip, &ev);
   SCIPfreeBufferArray(scip, &idx);

   return SCIP_OKAY;
}

/** add multiple sparse eigenvector cuts
 *
 *  We use Algorithm 1 from [Dey et al: Cutting Plan Generation Through Sparse Principal Component Analysis] to produce
 *  maxncuts many sparse eigenvector cuts.
 */
static
SCIP_RETCODE addMultipleSparseCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_Bool             enforce,            /**< whether we are in enforcing */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   int                   blocksize,          /**< size of block */
   SCIP_Real*            fullmatrix,         /**< precomputed full matrix \f$ \sum_j A_j y_j - A_0 \f$ */
   SCIP_Real*            fullconstmatrix,    /**< precomputed full constant matrix */
   SCIP_Real*            eigenvector,        /**< original eigenvector */
   SCIP_Real             tol,                /**< tolerance to be used when computing eigenvectors  */
   int                   maxncuts,           /**< maximal number of cuts to be produced */
   SCIP_Real*            vector,             /**< temporary workspace (length blocksize) */
   SCIP_VAR**            vars,               /**< temporary workspace */
   SCIP_Real*            vals,               /**< temporary workspace */
   int*                  ncuts,              /**< pointer to store the number of produced cuts/constrains */
   SCIP_Bool*            success,            /**< pointer to store whether we have produced a cut/constraint */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_Real* fullmatrixcopy;
   SCIP_Real* modmatrix;
   SCIP_Real* submatrix;
   SCIP_Real* sparseev;
   SCIP_Real* liftedev;
   SCIP_Real* minev;
   SCIP_Real convergencetol;
   SCIP_Real eigenvalue;
   SCIP_Real maxeig;
   SCIP_Real scalar;
   int* support;
   int cnt = 0;
   int size;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( fullmatrix != NULL );
   assert( fullconstmatrix != NULL );
   assert( eigenvector != NULL );
   assert( vector != NULL );
   assert( ncuts != NULL );
   assert( result != NULL );
   assert( *result != SCIP_CUTOFF );

   *ncuts = 0;
   *success = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPdebugMsg(scip, "Executing addMultipleSparseCuts.\n");

   /* compute target size, use sparsifytargetsize if specified, and sparsifyfactor else */
   if ( conshdlrdata->sdpconshdlrdata->sparsifytargetsize > 0 )
      size = conshdlrdata->sdpconshdlrdata->sparsifytargetsize;
   else
      size = MAX(10, (int) conshdlrdata->sdpconshdlrdata->sparsifyfactor * consdata->nvars);

   /* if size is larger than blocksize, trigger a debug message and set size to blocksize */
   if ( size > blocksize )
   {
      SCIPdebugMsg(scip, "No sparse eigenvector cuts produced since desired sparsity %d is larger than SDP blocksize %d.\n", size, blocksize);
      return SCIP_OKAY;
   }

   assert( size <= blocksize );

   /* copy fullmatrix, since Lapack destroys it when computing eigenvalues! */
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrixcopy, blocksize * blocksize) );
   for (i = 0; i < blocksize * blocksize; i++)
      fullmatrixcopy[i] = fullmatrix[i];

   /* compute the largest eigenvalue of A(y), the eigenvector is not needed */
   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, fullmatrixcopy, blocksize, &maxeig, NULL) );
   SCIPdebugMsg(scip, "Largest eigenvalue of A(y): %.15g\n", maxeig);

   /* compute the modified matrix \lambda_{max} I - A(y), where \lambda_{max} is the largest eigenvalue of A(y) */
   SCIP_CALL( SCIPallocBufferArray(scip, &modmatrix, blocksize * blocksize) );
   for (i = 0; i < blocksize; i++)
   {
      for (j = 0; j < blocksize; j++)
      {
         if ( i == j )
            modmatrix[i * blocksize + j] = maxeig - fullmatrix[i * blocksize + j];
         else
            modmatrix[i * blocksize + j] = -fullmatrix[i * blocksize + j];
      }
   }

   /* allocate necessary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &submatrix, size * size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &liftedev, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &support, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sparseev, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &minev, blocksize) );

   /* compute sparse eigenvector with the truncated power method for the modified matrix */
   convergencetol = 1e-6;
   SCIP_CALL( truncatedPowerMethod(scip, blocksize, modmatrix, eigenvector, size, convergencetol, liftedev, support, &eigenvalue) );

   SCIPdebugMsg(scip, "Smallest sparse eigenvalue computed by TPower: %.15g\n", maxeig - eigenvalue);

   /* instead of adding the sparse eigenvector, use it as initial point for adding many sparse eigenvectors */
   /* compute v^T A(y) v, which is equal to maxeig - eigenvalue */
   scalar = maxeig - eigenvalue;

   if ( ! SCIPisFeasNegative(scip, scalar) )
   {
      SCIPdebugMsg(scip, "Did not add any sparse eigenvector cuts, since sparsified initial eigenvector cut does not cut off the solution!\n");
   }

   if ( maxncuts < 0 )
      maxncuts = INT_MAX;

   while ( SCIPisFeasNegative(scip, scalar) && *ncuts < maxncuts )
   {
      if ( conshdlrdata->sdpconshdlrdata->recomputesparseev )
      {
         SCIP_Real norm;

         /* Recompute smallest eigenvalue \lambda_{min} and corresponding unit norm eigenvector w of principal submatrix
          * A(y)_S, lift it to a full vector by adding zeros, and normalize, instead of using the eigenvector and
          * eigenvalue returned from TPower.
          */

         cnt = 0;
         for (i = 0; i < blocksize; i++)
         {
            if ( support[i] == 1 )
            {
               for (j = 0; j < blocksize; j++)
               {
                  if ( support[j] == 1 )
                     submatrix[cnt++] = fullmatrix[i * blocksize + j];
               }
            }
         }
         assert( cnt == size * size );

         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), TRUE, size, submatrix, 1, &eigenvalue, sparseev) );

         assert( SCIPisFeasNegative(scip, eigenvalue) );

         if ( eigenvalue >= -tol )
            break;

         SCIPdebugMsg(scip, "Smallest eigenvalue of principal submatrix: %.15g\n", eigenvalue);

         /* normalize sparse eigenvector */
         norm = 0.0;
         for (j = 0; j < size; j++)
            norm += sparseev[j] * sparseev[j];

         norm = sqrt(norm);
         for (j = 0; j < size; j++)
            sparseev[j] /= norm;

         /* lift eigenvector by setting all entries not in the support to 0 */
         cnt = 0;
         for (i = 0; i < blocksize; i++)
         {
            if ( support[i] == 1 )
            {
               /* assert( SCIPisFeasEQ(scip, liftedev[i], sparseev[cnt]) ); */
               liftedev[i] = sparseev[cnt++];
            }
            else
               liftedev[i] = 0.0;
         }
         assert( cnt == size );
      }
      else
      {
         /* scalar is the smallest sparse eigenvalue computed with TPower */
         eigenvalue = scalar;
      }

      /* produce cut/constraint */
      SCIP_CALL( produceCutFromEigenvector(scip, conshdlr, conshdlrdata, cons, consdata, enforce, sol,
            blocksize, fullconstmatrix, liftedev, vector, vars, vals, ncuts, success, result) );

      /* compute A(y) = A(y) - \lambda_{min} w w^T */
      for (i = 0; i < blocksize; i++)
      {
         for (j = 0; j < blocksize; j++)
            fullmatrix[i * blocksize + j] -= eigenvalue * liftedev[i] * liftedev[j];
      }

      if ( conshdlrdata->sdpconshdlrdata->recomputeinitial )
      {
         SCIP_Real mineig;
         /* compute smallest eigenvalue and corresponding eigenvector of A(y) as initial vector for TPpower, instead of
            using the eigenvector to the smallest eigenvalue of the original matrix */

         /* copy fullmatrix, since Lapack destroys it when computing eigenvalues! */
         for (i = 0; i < blocksize * blocksize; i++)
            fullmatrixcopy[i] = fullmatrix[i];

         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), TRUE, blocksize, fullmatrixcopy, 1, &mineig, minev) );

         SCIPdebugMsg(scip, "Smallest eigenvalue: %.15g\n", mineig);
      }
      else
      {
         for (i = 0; i < blocksize; i++)
            minev[i] = eigenvector[i];
      }

      if ( conshdlrdata->sdpconshdlrdata->exacttrans )
      {
         /* we need to modify the matrix again in order to use the truncated power method */

         /* copy fullmatrix, since Lapack destroys it when computing eigenvalues! */
         for (i = 0; i < blocksize * blocksize; i++)
            fullmatrixcopy[i] = fullmatrix[i];

         /* compute the largest eigenvalue of A(y), the eigenvector is not needed */
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, fullmatrixcopy, blocksize, &maxeig, NULL) );
         SCIPdebugMsg(scip, "Largest eigenvalue: %.15g\n", maxeig);
      }
      else
      {
         /* use estimate for largest eigenvalue of modified matrix: \lambda_{max}(A+B) \leq \lambda_{max}(A) + \lambda_{max}(B) */
         maxeig -= eigenvalue;
      }

      /* compute the modified matrix \lambda_{max} I - A(y), where \lambda_{max} is the largest eigenvalue of A(y) */
      for (i = 0; i < blocksize; i++)
      {
         for (j = 0; j < blocksize; j++)
         {
            if ( i == j )
               modmatrix[i * blocksize + j] = maxeig - fullmatrix[i * blocksize +j];
            else
               modmatrix[i * blocksize + j] = -fullmatrix[i * blocksize +j];
         }
      }

      /* compute smallest sparse eigenvalue and corresponding eigenvector of A(y) with the truncated power method */
      SCIP_CALL( truncatedPowerMethod(scip, blocksize, modmatrix, minev, size, convergencetol, liftedev, support, &eigenvalue) );

      SCIPdebugMsg(scip, "Smallest sparse eigenvalue computed by TPower: %.15g\n", maxeig - eigenvalue);

      /* compute v^T A(y) v for the new sparse eigenvector, which is equal to maxeig - eigenvalue, as computed by TPower */
      scalar = maxeig - eigenvalue;
   }

   if ( *ncuts > 0 || *result == SCIP_CUTOFF )
   {
      *success = TRUE;
      SCIPdebugMsg(scip, "Added %d sparse eigenvector cuts\n", *ncuts);
   }

   /* free all memory */
   SCIPfreeBufferArray(scip, &minev);
   SCIPfreeBufferArray(scip, &sparseev);
   SCIPfreeBufferArray(scip, &support);
   SCIPfreeBufferArray(scip, &liftedev);
   SCIPfreeBufferArray(scip, &submatrix);
   SCIPfreeBufferArray(scip, &modmatrix);
   SCIPfreeBufferArray(scip, &fullmatrixcopy);

   return SCIP_OKAY;
}


/** separate current solution with a cut using the eigenvectors and -values of the solution matrix */
static
SCIP_RETCODE separateSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< the constraint handler itself */
   SCIP_CONS*            cons,               /**< constraint to process */
   SCIP_SOL*             sol,                /**< primal solution that should be separated */
   SCIP_Bool             enforce,            /**< whether we are enforcing cuts */
   SCIP_RESULT*          result              /**< pointer to store the result of the separation call */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_RETCODE retcode;
   SCIP_VAR** vars;
   SCIP_Real* vals;
   SCIP_Real* eigenvectors;
   SCIP_Real* vector;
   SCIP_Real* fullmatrix;
   SCIP_Real* fullmatrixcopy;
   SCIP_Real* fullconstmatrix = NULL;
   SCIP_Real* eigenvalues;
   SCIP_Real tol;
   int neigenvalues;
   int blocksize;
   int nvars;
   int ngen = 0;
   int i;
   int j;

   assert( cons != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* prepare computations */
   nvars = consdata->nvars;
   blocksize = consdata->blocksize;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );

   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrixcopy, blocksize * blocksize ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvectors, blocksize * blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vector, blocksize) );

   /* compute the matrix \f$ \sum_j A_j y_j - A_0 \f$ */
   SCIP_CALL( computeFullSdpMatrix(scip, consdata, sol, fullmatrix) );

   /* copy fullmatrix, since Lapack destroys it when computing eigenvalues! */
   for (i = 0; i < blocksize * blocksize; i++)
      fullmatrixcopy[i] = fullmatrix[i];

   /* determine tolerance */
   if ( conshdlrdata->sdpconshdlrdata->usedimacsfeastol )
   {
      assert( conshdlrdata->dimacsfeastol != SCIP_INVALID );
      tol = conshdlrdata->dimacsfeastol;
   }
   else
   {
      if ( enforce )
         tol = SCIPfeastol(scip);
      else
         tol = SCIPgetSepaMinEfficacy(scip);
   }

   /* compute eigenvector(s) */
   if ( conshdlrdata->sdpconshdlrdata->separateonecut || conshdlrdata->sdpconshdlrdata->multiplesparsecuts )
   {
      /* compute smallest eigenvalue */
      retcode = SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), TRUE, blocksize, fullmatrixcopy, 1, eigenvalues, eigenvectors);
      if ( retcode == SCIP_OKAY )
      {
         if ( eigenvalues[0] < -tol )
            neigenvalues = 1;
         else
            neigenvalues = 0;
      }
   }
   else
   {
      /* compute all eigenvectors for negative eigenvalues */
      retcode = SCIPlapackComputeEigenvectorsNegative(SCIPbuffer(scip), blocksize, fullmatrixcopy, tol, &neigenvalues, eigenvalues, eigenvectors);
   }

   /* treat possible error */
   if ( retcode != SCIP_OKAY )
   {
      SCIPfreeBufferArray(scip, &vector);
      SCIPfreeBufferArray(scip, &eigenvalues);
      SCIPfreeBufferArray(scip, &eigenvectors);
      SCIPfreeBufferArray(scip, &fullmatrixcopy);
      SCIPfreeBufferArray(scip, &fullmatrix);

      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);

      /* if we are enforcing, generate error */
      if ( enforce )
      {
         SCIP_CALL( retcode );
      }
      else
         return SCIP_OKAY;  /* in separation just exit */
   }

   if ( neigenvalues > 0 )
   {
      /* there are no variables, but the matrix is negative definite -> cutoff */
      if ( consdata->nvars == 0 )
         *result = SCIP_CUTOFF;
      else
      {
         /* get full constant matrix */
         SCIP_CALL( SCIPallocBufferArray(scip, &fullconstmatrix, blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, cons, fullconstmatrix) );
      }
   }

   /* loop through eigenvectors and add cuts as long as cut is violated */
   for (i = 0; i < neigenvalues && *result != SCIP_CUTOFF; ++i)
   {
      SCIP_Real* eigenvector;
      SCIP_Bool success;
      int ncuts;

      /* get pointer to current eigenvector */
      eigenvector = &(eigenvectors[i * blocksize]);

      /* if we want to sparsify the cut */
      if ( ! enforce && conshdlrdata->sdpconshdlrdata->sparsifycut )
      {
         SCIP_CALL( sparsifyCut(scip, conshdlr, conshdlrdata, cons, consdata, enforce, sol, blocksize, fullconstmatrix, eigenvector, vector, vars, vals, &ngen, &success, result) );

         if ( success )
         {
            ++ngen;
            continue;
         }
      }
      else if ( conshdlrdata->sdpconshdlrdata->multiplesparsecuts )
      {
         SCIPdebugMsg(scip, "Smallest eigenvalue: %.15g\n", eigenvalues[i]);
         SCIP_CALL( addMultipleSparseCuts(scip, conshdlr, conshdlrdata, cons, consdata, enforce, sol, blocksize, fullmatrix, fullconstmatrix,
               eigenvector, tol, conshdlrdata->sdpconshdlrdata->maxnsparsecuts, vector, vars, vals, &ncuts, &success, result) );

         if ( success )
         {
            SCIPdebugMsg(scip, "Successfully added %d sparse eigenvector cuts.\n", ncuts);
            ngen += ncuts;
            continue;
         }
      }
      else
      {
         /* to avoid numerical trouble, we eliminate small entries in absolute value */
         for (j = 0; j < consdata->blocksize; ++j)
         {
            if ( SCIPisFeasZero(scip, eigenvector[j]) )
               eigenvector[j] = 0.0;
         }
      }

      /* produce cut/constraint */
      SCIP_CALL( produceCutFromEigenvector(scip, conshdlr, conshdlrdata, cons, consdata, enforce, sol,
            blocksize, fullconstmatrix, eigenvector, vector, vars, vals, &ngen, &success, result) );
   }
   SCIPdebugMsg(scip, "<%s>: Separated cuts = %d.\n", SCIPconsGetName(cons), ngen);

   SCIPfreeBufferArrayNull(scip, &fullconstmatrix);
   SCIPfreeBufferArray(scip, &vector);
   SCIPfreeBufferArray(scip, &eigenvalues);
   SCIPfreeBufferArray(scip, &eigenvectors);
   SCIPfreeBufferArray(scip, &fullmatrixcopy);
   SCIPfreeBufferArray(scip, &fullmatrix);

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** check whether all matrices are psd
 *
 *  Only needs to be called by rank-1 constraints to check whether all matrices are psd.
 */
static
SCIP_RETCODE computeAllmatricespsd(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* Aj;
   int blocksize;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );
   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->rankone );

   /* exit if allmatricespsd has been initialized already */
   if ( consdata->initallmatricespsd )
      return SCIP_OKAY;

   assert( consdata->allmatricespsd == FALSE );

   blocksize = consdata->blocksize;
   SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) );

   SCIPdebugMsg(scip, "Computing allmatricespsd for constraint <%s>.\n", SCIPconsGetName(cons));

   consdata->allmatricespsd = TRUE;
   for (v = 0; v < consdata->nvars && consdata->allmatricespsd; ++v)
   {
      SCIP_Real eigenvalue;

      SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, v, Aj) );

      /* compute minimal eigenvalue */
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
         consdata->allmatricespsd = FALSE;
   }

   SCIPfreeBufferArray(scip, &Aj);

   consdata->initallmatricespsd = TRUE;

   return SCIP_OKAY;
}

/** try to tighten matrices if all matrices are psd
 *
 *  Try to scale matrices without changing feasible solutions. The details are explained in the presolving paper (see
 *  the top of the file).
 */
static
SCIP_RETCODE tightenMatrices(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints to add cuts for */
   int                   nconss,             /**< number of constraints to add cuts for */
   int*                  nchgcoefs           /**< pointer to store how many matrices were tightened */
   )
{
   int c;

   assert( scip != NULL );
   assert( nchgcoefs != NULL );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_Real* constmatrix;
      int blocksize;
      int nvars;
      int i;

      assert( conss != NULL );
      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLR_NAME) == 0 );
      assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLRRANK1_NAME) == 0 );

      /* compute allmatricespsd for rank-1 constraints */
      if ( consdata->rankone )
      {
         SCIP_CALL( computeAllmatricespsd(scip, conss[c]) );
      }

      /* skip constraints in which not all matrices are psd */
      if ( ! consdata->allmatricespsd )
         continue;

      /* make sure that all lower bounds are nonnegative */
      nvars = consdata->nvars;
      for (i = 0; i < nvars; ++i)
      {
         if ( SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->vars[i])) )
            break;
      }
      if ( i < nvars )
         continue;

      SCIPdebugMsg(scip, "Trying to tighten matrices for constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* get matrices */
      blocksize = consdata->blocksize;
      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );

      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

      for (i = 0; i < nvars; ++i)
      {
         SCIP_Real objval;
         SCIP_Real factor;
         SCIP_Real lb;
         SCIP_Real ub;

         /* only treat binary variables */
         if ( ! SCIPvarIsBinary(consdata->vars[i]) )
            continue;

         /* skip fixed variables (will be removed anyway */
         lb = SCIPvarGetLbLocal(consdata->vars[i]);
         ub = SCIPvarGetUbLocal(consdata->vars[i]);
         if ( SCIPisEQ(scip, lb, ub) )
            continue;

         assert( SCIPisEQ(scip, lb, 0.0) );
         assert( SCIPisEQ(scip, ub, 1.0) );

         /* solve 1d SDP */
         SCIP_CALL( SCIPsolveOneVarSDPDense(SCIPbuffer(scip), 1.0, 0.0, 1.0, blocksize, constmatrix, consdata->nvarnonz[i], consdata->row[i], consdata->col[i], consdata->val[i],
               SCIPinfinity(scip), SCIPfeastol(scip), &objval, &factor) );

         if ( SCIPisInfinity(scip, objval) )
            continue;

         if ( factor == SCIP_INVALID ) /*lint !e777*/
            continue;

         if ( ! SCIPisFeasEQ(scip, factor, 1.0) )
         {
            int j;

            SCIPdebugMsg(scip, "Tightened coefficent matrix of variable <%s> with tightening factor %g.\n", SCIPvarGetName(consdata->vars[i]), factor);

            /* tighten matrix */
            for (j = 0; j < consdata->nvarnonz[i]; j++)
               consdata->val[i][j] *= factor;

            ++(*nchgcoefs);
         }
      }

      SCIPfreeBufferArray(scip, &constmatrix);
   }

   return SCIP_OKAY;
}

/** try to tighten bounds if all matrices are psd
 *
 *  Try to compute tighter variable bounds by solving one-variable SDPs. The details are explained in the presolving
 *  paper (see the top of the file).
 */
static
SCIP_RETCODE tightenBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints to add cuts for */
   int                   nconss,             /**< number of constraints to add cuts for */
   SCIP_Bool             tightenboundscont,  /**< Should only continuous variables be tightened? */
   int*                  nchgbds,            /**< pointer to store how many bounds were tightened */
   int*                  nintrnd,            /**< pointer to store how many tightened bounds of integer variables were rounded to be integral */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected */
   )
{
   int c;

   assert( scip != NULL );
   assert( nchgbds != NULL );
   assert( nintrnd != NULL );
   assert( infeasible != NULL );

   *nchgbds = 0;
   *nintrnd = 0;
   *infeasible = FALSE;

   for (c = 0; c < nconss && !(*infeasible); ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_Real* matrix = NULL;
      SCIP_Real* constmatrix;
      SCIP_Real factor;
      SCIP_Bool havebinaryvar = FALSE;
      int blocksize;
      int nvars;
      int i;

      assert( conss != NULL );
      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLR_NAME) == 0 );
      assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLRRANK1_NAME) == 0 );

      /* compute allmatricespsd for rank-1 constraints */
      if ( consdata->rankone )
      {
         SCIP_CALL( computeAllmatricespsd(scip, conss[c]) );
      }

      /* skip constraints in which not all matrices are psd */
      if ( ! consdata->allmatricespsd )
         continue;

      /* make sure that all upper bounds finite */
      nvars = consdata->nvars;
      for (i = 0; i < nvars; ++i)
      {
         if ( SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->vars[i])) )
            break;

         if ( SCIPvarIsBinary(consdata->vars[i]) )
            havebinaryvar = TRUE;
      }
      if ( i < nvars )
         continue;

      SCIPdebugMsg(scip, "Trying to tighten bounds for constraint <%s>.\n", SCIPconsGetName(conss[c]));

      /* get matrices */
      blocksize = consdata->blocksize;
      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
      if ( havebinaryvar )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * blocksize) );
      }

      for (i = 0; i < nvars; ++i)
      {
         SCIP_Real* sdpval;
         SCIP_Real lb;
         SCIP_Real ub;
         SCIP_Real ubk;
         SCIP_Real objval;
         int* sdprow;
         int* sdpcol;
         int sdpnnonz;
         int row;
         int col;
         int k;
         int l;

         /* possibly restrict tightening to continuous variables */
         if ( tightenboundscont && SCIPvarIsIntegral(consdata->vars[i]) )
            continue;

         /* skip fixed variables */
         lb = SCIPvarGetLbLocal(consdata->vars[i]);
         ub = SCIPvarGetUbLocal(consdata->vars[i]);
         if ( SCIPisEQ(scip, lb, ub) )
            continue;

         /* skip infinite lower bound (cannot currently deal with this in SCIPsolveOneVarSDPDense()) */
         assert( ! SCIPisInfinity(scip, ub) );
         if ( SCIPisInfinity(scip, -lb) )
            continue;

         /* get fresh copy of the constant matrix */
         SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

         /* loop over other variables */
         for (k = 0; k < nvars; ++k)
         {
            if ( k == i )
               continue;

            ubk = SCIPvarGetUbLocal(consdata->vars[k]);

            /* subtract matrix times upper bound from constant matrix (because of minus const. matrix) */
            if ( ! SCIPisZero(scip, ubk) )
            {
               sdprow = consdata->row[k];
               sdpcol = consdata->col[k];
               sdpval = consdata->val[k];
               sdpnnonz = consdata->nvarnonz[k];
               for (l = 0; l < sdpnnonz; ++l)
               {
                  row = sdprow[l];
                  col = sdpcol[l];
                  constmatrix[row * blocksize + col] -= sdpval[l] * ubk;
                  if ( row != col )
                     constmatrix[col * blocksize + row] -= sdpval[l] * ubk;
               }
            }
         }

         /* special handling of binary variables */
         if ( SCIPvarIsBinary(consdata->vars[i]) )
         {
            SCIP_Real eigenvalue;

            assert( matrix != NULL );

            /* first copy constant matrix (possibly needed later) */
            for (l = 0; l < blocksize * blocksize; ++l)
               matrix[l] = - constmatrix[l];

            /* compute maximal eigenvalue; this corresponds to checking the lower bound because of minus sign; constmatrix is destroyed. */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, constmatrix, blocksize, &eigenvalue, NULL) );

            /* take minus sign into account */
            eigenvalue = -eigenvalue;

            /* if the negative of the constant matrix is psd, then the lower bound (= 0) is feasible -> cannot tighten */
            if ( SCIPisFeasGE(scip, eigenvalue, 0.0) )
               continue;
            assert( SCIPisFeasNegative(scip, eigenvalue) );

            /* otherwise, we check whether the upper bound is feasible */

            /* add matrix i for evaluating upper bound */
            sdprow = consdata->row[i];
            sdpcol = consdata->col[i];
            sdpval = consdata->val[i];
            sdpnnonz = consdata->nvarnonz[i];
            for (l = 0; l < sdpnnonz; ++l)
            {
               row = sdprow[l];
               col = sdpcol[l];
               matrix[row * blocksize + col] += sdpval[l];
               if ( row != col )
                  matrix[col * blocksize + row] += sdpval[l];
            }

            /* compute minimal eigenvalue; matrix is destroyed */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, matrix, 1, &eigenvalue, NULL) );

            /* if matrix is not psd, we are infeasible */
            if ( SCIPisFeasNegative(scip, eigenvalue) )
            {
               SCIPdebugMsg(scip, "Binary variable <%s> yields infeasibility.\n", SCIPvarGetName(consdata->vars[i]));
               *infeasible = TRUE;
               break;
            }

            /* otherwise, we need at least the matrix for variable i to become feasible -> can tighten variable to 1 */
            factor = 0.5;
            SCIPdebugMsg(scip, "Tighten lower bound for binary variable <%s> to 1.\n", SCIPvarGetName(consdata->vars[i]));
         }
         else
         {
            /* solve 1d SDP */
            SCIP_CALL( SCIPsolveOneVarSDPDense(SCIPbuffer(scip), 1.0, lb, ub, blocksize, constmatrix, consdata->nvarnonz[i], consdata->row[i], consdata->col[i], consdata->val[i],
                  SCIPinfinity(scip), SCIPfeastol(scip), &objval, &factor) );

            /* if problem is infeasible */
            if ( SCIPisInfinity(scip, objval) )
            {
               SCIPdebugMsg(scip, "Variable <%s> yields infeasibility.\n", SCIPvarGetName(consdata->vars[i]));
               *infeasible = TRUE;
               break;
            }

            if ( factor == SCIP_INVALID ) /*lint !e777*/
               continue;
         }

         if ( SCIPisGT(scip, factor, lb) )
         {
            SCIP_Bool tightened;

            SCIP_CALL( SCIPinferVarLbCons(scip, consdata->vars[i], factor, conss[c], -(i+1), FALSE, infeasible, &tightened) );

            if ( *infeasible )
               break;

            if ( tightened )
            {
               SCIPdebugMsg(scip, "%sTightened lower bound of variable <%s> from %g to %g.\n", SCIPinProbing(scip) ? "In probing: " : "", SCIPvarGetName(consdata->vars[i]), lb, factor);
               ++(*nchgbds);

               /*  if variable is integral, the bound change should automatically produce an integer bound */
               if ( SCIPvarIsIntegral(consdata->vars[i]) && ! SCIPisFeasIntegral(scip, factor) )
               {
                  assert( SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(consdata->vars[i])) );
                  ++(*nintrnd);
               }
            }
         }
      }

      SCIPfreeBufferArrayNull(scip, &matrix);
      SCIPfreeBufferArray(scip, &constmatrix);
   }

   return SCIP_OKAY;
}

/** approximates the sdpcone using the fact that every diagonal entry must be non-negative, so it adds the LP-cut
 *  \f$ \sum_{j = 1}^m (A_j)_{kk} y_j - (A_0)_{kk} \geq 0 \quad \forall k \leq n \f$
 *
 *  For a context, see also the presolving paper (mentioned at the top of the file).
 */
static
SCIP_RETCODE diagGEzero(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints to add cuts for */
   int                   nconss,             /**< number of constraints to add cuts for */
   SCIP_Bool             solvesdps,          /**< are we solving SDPs or LPs? */
   int*                  naddconss,          /**< pointer to store how many constraints were added */
   int*                  nchgbds,            /**< pointer to store how many bounds were changed */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real* matrix;
   SCIP_Real* consvals;
   SCIP_VAR** consvars;
   SCIP_Real* diagentries;
   int blocksize;
   int nvars;
   int i;
   int j;
   int k;
   int c;

   assert( scip != NULL );
   assert( naddconss != NULL );
   assert( nchgbds != NULL );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   for (c = 0; c < nconss; ++c)
   {
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLR_NAME) == 0 );
      assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLRRANK1_NAME) == 0 );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

      SCIP_CALL( SCIPallocBufferArray(scip, &matrix, (blocksize * (blocksize + 1)) / 2) ); /*lint !e647*/
      SCIP_CALL( SCIPconsSdpGetLowerTriangConstMatrix(scip, conss[c], matrix) );

      /* allocate diagonal entries and init to 0.0 */
      SCIP_CALL( SCIPallocClearBufferArray(scip, &diagentries, blocksize * nvars) ); /*lint !e647*/

      /* get the (k,k)-entry of every matrix A_j */
      for (j = 0; j < nvars; ++j)
      {
         /* go through the nonzeros of A_j and look for diagonal entries */
         for (i = 0; i < consdata->nvarnonz[j]; i++)
         {
            if ( consdata->col[j][i] == consdata->row[j][i] )
               diagentries[consdata->col[j][i] * nvars + j] = consdata->val[j][i]; /*lint !e679*/
         }
      }

      /* add the LP-cuts to SCIP */
      for (k = 0; k < blocksize && !(*infeasible); ++k)
      {
         SCIP_CONS* cons;
         SCIP_Real lhs;
         SCIP_Real activitylb = 0.0;
         int cnt = 0;

         for ( j = 0; j < nvars; ++j)
         {
            SCIP_Real val;
            SCIP_VAR* var;

            val = diagentries[k * nvars + j];
            if ( ! SCIPisZero(scip, val) )
            {
               var = consdata->vars[j];
               consvals[cnt] = val;
               consvars[cnt++] = var;

               /* compute lower bound on activity */
               if ( val > 0 )
                  activitylb += val * SCIPvarGetLbGlobal(var);
               else
                  activitylb += val * SCIPvarGetUbGlobal(var);
            }
         }

         /* the lhs is the (k,k)-entry of the constant matrix */
         lhs = matrix[SCIPconsSdpCompLowerTriangPos(k,k)];

         if ( cnt == 0 )
         {
            /* if there are no variables, but the lhs is positive, we are infeasible */
            if ( SCIPisPositive(scip, lhs) )
               *infeasible = TRUE;
         }
         else if ( cnt == 1 && SCIPvarGetStatus(consvars[0]) != SCIP_VARSTATUS_MULTAGGR )
         {
            SCIP_Bool tightened;
            SCIP_VAR* var;
            SCIP_Real val;

            var = consvars[0];
            val = consvals[0];
            assert( var != NULL );

            /* try to tighten bound */
            if ( SCIPisPositive(scip, val) )
            {
               SCIP_CALL( SCIPtightenVarLb(scip, var, lhs / val, FALSE, infeasible, &tightened) );
               if ( tightened )
               {
                  SCIPdebugMsg(scip, "Tightend lower bound of <%s> to %g because of diagonal values of SDP-constraint <%s>!\n",
                     SCIPvarGetName(var), SCIPvarGetLbGlobal(var), SCIPconsGetName(conss[c]));
                  ++(*nchgbds);
               }
            }
            else if ( SCIPisNegative(scip, val) )
            {
               SCIP_CALL( SCIPtightenVarUb(scip, var, lhs / val, FALSE, infeasible, &tightened) );
               if ( tightened )
               {
                  SCIPdebugMsg(scip, "Tightend upper bound of <%s> to %g because of diagonal values of SDP-constraint <%s>!\n",
                     SCIPvarGetName(var), SCIPvarGetUbGlobal(var), SCIPconsGetName(conss[c]));
                  ++(*nchgbds);
               }
            }
         }
         /* generate linear inequality if lower bound on activity is less than the lhs, so the cut is not redundant */
         else if ( SCIPisLT(scip, activitylb, lhs) )
         {
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "diag_ge_zero_%d", ++(conshdlrdata->ndiaggezerocuts));

            /* Only separate if solving LPs */
            if ( conshdlrdata->sdpconshdlrdata->presollinconssparam == 1 )
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, cnt, consvars, consvals, lhs, SCIPinfinity(scip),
                     TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/
            else
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, cnt, consvars, consvals, lhs, SCIPinfinity(scip),
                     FALSE, ! solvesdps, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/

            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added diagonal lp-constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            ++(*naddconss);
         }
      }

      SCIPfreeBufferArray(scip, &diagentries);
      SCIPfreeBufferArray(scip, &matrix);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** Presolve-routine that enforces implications of diagonal entries of zero in SDP-matrices, namely that if \f$X_{ij} > 0\f$,
 *  then also \f$X_{ii} > 0\f$ and \f$ X_{jj} > 0\f$.
 *
 *  More precisely, if \f$ (A_0)_{k\ell} \neq 0\f$, \f$ (A_0)_{kk} = 0\f$, \f$ (A_i)_{k\ell} = 0\f$ for all \f$ i \leq m\f$,
 *  \f$ (A_i)_{kk} = 0\f$ for all continuous variables and \f$ \ell_i \geq 0\f$ for all integer variables, we add the cut
 *  \f$ \sum_{\substack{i \in \mathcal{I}:\\ (A_i)_{kk} > 0}} y_i \geq 1.\f$
 *
 *  The details are explained in the presolving paper (see the top of the file).
 */
static
SCIP_RETCODE diagZeroImpl(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   int c;

   assert( scip != NULL );
   assert( naddconss != NULL );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_Shortbool* nonzeroentries;
      SCIP_Shortbool* diagnonzero;
      SCIP_VAR** vars;
      SCIP_Real* vals;
      SCIP_Bool negintvar =  FALSE;
      int** diagvars;
      int* ndiagvars;
      int ndiagnonzero = 0;
      int blocksize;
      int rowidx;
      int colidx;
      int nvars;
      int pos;
      int j;
      int v;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      /* allocate storage */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &nonzeroentries, blocksize * (blocksize+1) / 2) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &diagnonzero, blocksize) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &ndiagvars, blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &diagvars, blocksize) );
      for (j = 0; j < blocksize; ++j)
      {
         SCIP_CALL( SCIPallocClearBufferArray(scip, &diagvars[j], nvars) );
      }

      /* collect nonzero entries of matrices */
      for (v = 0; v < nvars; ++v)
      {
         if ( SCIPvarIsIntegral(consdata->vars[v]) && SCIPisNegative(scip, SCIPvarGetLbGlobal(consdata->vars[v])) )
         {
            negintvar = TRUE;
            break;
         }

         for (j = 0; j < consdata->nvarnonz[v]; j++)
         {
            rowidx = consdata->row[v][j];
            colidx = consdata->col[v][j];
            assert( 0 <= rowidx && rowidx < blocksize );
            assert( 0 <= colidx && colidx < blocksize );

            pos = rowidx * (rowidx + 1)/2 + colidx;
            nonzeroentries[pos] = TRUE;

            /* treat diagonal entries */
            if ( rowidx == colidx )
            {
               /* collect variables for positive diagonal entries */
               if ( SCIPisPositive(scip, consdata->val[v][j]) )
                  diagvars[rowidx][ndiagvars[rowidx]++] = v;

               /* mark nonzero entries for non-integral variables */
               if ( ! SCIPvarIsIntegral(consdata->vars[v]) )
               {
                  if ( ! diagnonzero[rowidx] )
                  {
                     diagnonzero[rowidx] = TRUE;
                     ++ndiagnonzero;
                  }
               }
            }
         }
      }

      /* early termination if there is an integral variable with a negative lower bound */
      if ( negintvar )
      {
         for (j = blocksize - 1; j >= 0; j--)
         {
            SCIPfreeBufferArray(scip, &diagvars[j]);
         }
         SCIPfreeBufferArray(scip, &diagvars);
         SCIPfreeBufferArray(scip, &ndiagvars);
         SCIPfreeBufferArray(scip, &diagnonzero);
         SCIPfreeBufferArray(scip, &nonzeroentries);
         SCIPfreeBufferArray(scip, &vals);
         SCIPfreeBufferArray(scip, &vars);
         continue;
      }

      /* add nonzero diagonal entries of constant matrix */
      for (j = 0; j < consdata->constnnonz; j++)
      {
         rowidx = consdata->constrow[j];
         colidx = consdata->constcol[j];
         assert( 0 <= colidx && colidx < blocksize );
         assert( 0 <= rowidx && rowidx < blocksize );
         assert( ! SCIPisZero(scip, consdata->constval[j]) );

         if ( rowidx == colidx )
         {
            if ( ! diagnonzero[rowidx] )
            {
               diagnonzero[rowidx] = TRUE;
               ++ndiagnonzero;
            }
         }
      }
      assert( 0 <= ndiagnonzero && ndiagnonzero <= blocksize );

      /* early termination if all diagonals are marked to be nonzero */
      if ( ndiagnonzero >= blocksize )
      {
         for (j = blocksize - 1; j >= 0; j--)
         {
            SCIPfreeBufferArray(scip, &diagvars[j]);
         }
         SCIPfreeBufferArray(scip, &diagvars);
         SCIPfreeBufferArray(scip, &ndiagvars);
         SCIPfreeBufferArray(scip, &diagnonzero);
         SCIPfreeBufferArray(scip, &nonzeroentries);
         SCIPfreeBufferArray(scip, &vals);
         SCIPfreeBufferArray(scip, &vars);
         continue;
      }

      /* iterate over all nonzeros of the constant matrix to produce cuts */
      for (j = 0; j < consdata->constnnonz; j++)
      {
         SCIP_CONS* cons;

         rowidx = consdata->constrow[j];
         colidx = consdata->constcol[j];
         assert( 0 <= colidx && colidx < blocksize );
         assert( 0 <= rowidx && rowidx < blocksize );
         assert( ! SCIPisZero(scip, consdata->constval[j]) );

         /* skip diagonal entries */
         if ( rowidx == colidx )
            continue;

         pos = rowidx * (rowidx + 1)/2 + colidx;

         /* skip entry if it is non-zero in some non-constant matrix */
         if ( nonzeroentries[pos] )
            continue;

         /* if all continuous variables have a zero diagonal entry and the constant matrix is 0 as well */
         if ( ! diagnonzero[rowidx] && ndiagvars[rowidx] > 0 )
         {
            /* get the corresponding SCIP variables and set all coefficients to 1 */
            for (v = 0; v < ndiagvars[rowidx]; ++v)
            {
               assert( SCIPvarIsIntegral(consdata->vars[diagvars[rowidx][v]]) );
               assert( 0 <= diagvars[rowidx][v] && diagvars[rowidx][v] < nvars );
               vars[v] = consdata->vars[diagvars[rowidx][v]];
               vals[v] = 1.0;
            }
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "diag_0_impl_row_%d_%d", rowidx, colidx);

            /* add the linear constraint sum_v 1.0 * diagvars[v] >= 1.0 */
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ndiagvars[rowidx], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            (*naddconss)++;

            /* mark diagonal entry as covered */
            diagnonzero[rowidx] = TRUE;
         }

         /* same possibility for column index */
         if ( ! diagnonzero[colidx] && ndiagvars[colidx] > 0 )
         {
            /* get the corresponding SCIP variables and set all coefficients to 1 */
            for (v = 0; v < ndiagvars[colidx]; ++v)
            {
               assert( SCIPvarIsIntegral(consdata->vars[diagvars[colidx][v]]) );
               assert( 0 <= diagvars[colidx][v] && diagvars[colidx][v] < nvars );
               vars[v] = consdata->vars[diagvars[colidx][v]];
               vals[v] = 1.0;
            }
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "diag_0_impl_col_%d_%d", rowidx, colidx);

            /* add the linear constraint sum_v 1.0 * diagvars[v] >= 1.0 */
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, ndiagvars[colidx], vars, vals, 1.0, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            (*naddconss)++;

            /* mark diagonal entry as covered */
            diagnonzero[colidx] = TRUE;
         }
      }

      /* free space */
      for (j = blocksize - 1; j >= 0; j--)
      {
         SCIPfreeBufferArray(scip, &diagvars[j]);
      }
      SCIPfreeBufferArray(scip, &diagvars);
      SCIPfreeBufferArray(scip, &ndiagvars);
      SCIPfreeBufferArray(scip, &diagnonzero);
      SCIPfreeBufferArray(scip, &nonzeroentries);
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &vars);
   }

   return SCIP_OKAY;
}

/** presolve-routine that adds linear constraints arising from 2 by 2 minor inequalities
 *
 *  For a positive semidefinite matrix \f$X\f$ the following two inequalities hold: \f$X_{ss} + X_{tt} - 2\, X_{st} \geq
 *  0\f$ and \f$X_{ss} + X_{tt} + 2\, X_{st} \geq 0\f$. This follows by using a 2 by 2 minor and multiplying from left
 *  and right by the all-ones vector and \f$[1,-1]\f$, respectively. We add the corresponding linear constraint only to
 *  be propagated.
 *
 *  Translated to the matrix pencil notation the cut looks as follows:
 *  \f[
 *  \sum_{i=1}^m (A_i)_{ss}\, y_i - (A_0)_{ss} + \sum_{i=1}^m (A_i)_{tt}\, y_i - (A_0)_{tt} - 2 \Big(\sum_{i=1}^m (A_i)_{st}\, y_i - (A_0)_{st}\Big) \geq 0
 *  \quad\Leftrightarrow\quad
 *  \sum_{i=1}^m \Big((A_i)_{ss} + (A_i)_{tt} - 2\, (A_i)_{st}\Big)\, y_i \geq (A_0)_{ss} + (A_0)_{tt} - 2 (A_0)_{st}.
 *  \f]
 *
 *  @todo The cut \f$X_{ss} + X_{tt} - 2\, X_{st} \geq 0\f$ is a special form of an Eigenvector cut. Try out other
 *  Eigenvector cuts such as \f$X_{ss} + X_{tt} + 2\, X_{st} \geq 0\f$.
 *
 *  For more details, see the presolving paper (mentioned at the top of the file).
 */
static
SCIP_RETCODE addTwoMinorLinConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             solvesdps,          /**< are we solving SDPs or LPs? */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int blocksize;
   int nvars;
   int c;
   int i;

   assert( scip != NULL );
   assert( naddconss != NULL );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_Real** matrices = NULL;
      SCIP_Real* constmatrix;
      int s;
      int t;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert( conshdlrdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

      /* get matrices */
      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

      SCIP_CALL( SCIPallocBufferArray(scip, &matrices, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &matrices[i], blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullAj(scip, conss[c], i, matrices[i]) );
      }

      /* loop over all possible entries */
      for (s = 0; s < blocksize; ++s)
      {
         for (t = 0; t < s; ++t)
         {
            SCIP_CONS* cons;
            SCIP_Real val;
            SCIP_Real lhs;
            SCIP_Real activitylb = 0.0;
            int nconsvars = 0;
            int cnt = 0;

            /* collect coefficients */
            for (i = 0; i < nvars; ++i)
            {
               SCIP_Real coef = 0.0;

               val = matrices[i][s * blocksize + t];
               if ( ! SCIPisZero(scip, val) )
               {
                  coef = -2.0 * val;
                  ++cnt;
               }

               /* add diagonal entries for s and t */
               val = matrices[i][s * blocksize + s];
               if ( ! SCIPisZero(scip, val) )
                  coef += val;

               val = matrices[i][t * blocksize + t];
               if ( ! SCIPisZero(scip, val) )
                  coef += val;

               if ( ! SCIPisZero(scip, coef) )
               {
                  consvals[nconsvars] = coef;
                  consvars[nconsvars] = consdata->vars[i];
                  ++nconsvars;

                  /* compute lower bound on activity */
                  if ( coef > 0 )
                     activitylb += coef * SCIPvarGetLbGlobal(consdata->vars[i]);
                  else
                     activitylb += coef * SCIPvarGetUbGlobal(consdata->vars[i]);
               }
            }

            /* only proceed if off-diagonal is nonzero and cut is nontrivial */
            if ( cnt <= 0 || nconsvars <= 0 )
               continue;

            /* compute rhs */
            lhs = constmatrix[s * blocksize + s] + constmatrix[t * blocksize + t] - 2.0 * constmatrix[s * blocksize + t];

            /* only proceed if constraint is not redundant */
            if ( SCIPisGE(scip, activitylb, lhs) )
               continue;

            /* add linear constraint (if not solving LPs only propagate) */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "2x2minorlin#%d#%d", s, t);
            if ( conshdlrdata->sdpconshdlrdata->presollinconssparam == 1 )
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, lhs, SCIPinfinity(scip),
                     TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/
            else
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, lhs, SCIPinfinity(scip),
                     FALSE, ! solvesdps, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added 2x2 minor linear constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            ++(*naddconss);
         }
      }

      for (i = nvars - 1; i >= 0; --i)
         SCIPfreeBufferArray(scip, &matrices[i]);
      SCIPfreeBufferArray(scip, &matrices);
      SCIPfreeBufferArray(scip, &constmatrix);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** presolve-routine that adds SOC constraints arising from 2 by 2 minor inequalities
 *
 *  For a positive semidefinite matrix \f$X\f$ the following quadratic inequality holds by using the 2 by 2 minor:
 *  \f$X_{st}^2 \leq X_{ss} X_{tt}\f$. This constraint is SOC representable:
 *  \f[
 *     \left\|\binom{2\,X_{st}}{X_{ss} - X_{tt}}\right\| \leq X_{ss} + X_{tt}.
 *     \quad\Leftrightarrow\quad
 *     4\, X_{st}^2 + (X_{ss}^2 - 2\,X_{ss}\, X_{tt} + X_{tt}^2) \leq X_{ss}^2 + 2\, X_{ss}\,X_{tt} + X_{tt}^2.
 *  \f]
 *
 *  Translated to the matrix pencil notation the cut looks as follows:
 *  \f[
 *  \left\|\binom{2\, \sum_{i=1}^m (A_i)_{st}\, y_i - 2\, (A_0)_{st}}{
 *          \sum_{i=1}^m (A_i)_{ss}\, y_i - (A_0)_{ss} - \sum_{i=1}^m (A_i)_{tt}\, y_i + (A_0)_{tt}}\right\| \leq
 *  \sum_{i=1}^m (A_i)_{ss}\, y_i - (A_0)_{ss} + \sum_{i=1}^m (A_i)_{tt}\, y_i - (A_0)_{tt}.
 *  \f]
 *  We add add new variables for \f$X_{st}\f$ and corresponding linear constraints. Then we add SOC-constraints.
 */
static
SCIP_RETCODE addTwoMinorSOCConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             solvesdps,          /**< are we solving SDPs or LPs? */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR*** matrixvars;
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int blocksize;
   int nvars;
   int c;
   int i;

   assert( scip != NULL );
   assert( naddconss != NULL );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_Real** matrices = NULL;
      SCIP_Real* constmatrix;
      int size;
      int s;
      int t;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      size = MAX(nvars + 1, 3);
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, size) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, size) );

      /* get matrices */
      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

      SCIP_CALL( SCIPallocBufferArray(scip, &matrices, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &matrices[i], blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullAj(scip, conss[c], i, matrices[i]) );
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &matrixvars, blocksize) );

      /* loop over all possible entries */
      for (s = 0; s < blocksize; ++s)
      {
         SCIP_CALL( SCIPallocClearBufferArray(scip, &matrixvars[s], blocksize) );

         for (t = s; t >= 0; --t)
         {
            SCIP_VAR* matrixvar = NULL;
            SCIP_VAR* matrixsumvar;
            SCIP_VAR* matrixdiffvar;
            SCIP_CONS* cons;
            SCIP_Real val;
            SCIP_Real lhs;
            int nconsvars = 0;

            /* create linear constraint defining matrix variable */
            for (i = 0; i < nvars; ++i)
            {
               val = matrices[i][s * blocksize + t];
               if ( ! SCIPisZero(scip, val) )
               {
                  consvals[nconsvars] = val;
                  consvars[nconsvars++] = consdata->vars[i];
               }
            }

            /* only proceed if entry is nontrivial */
            if ( nconsvars == 0 )
               continue;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "M%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateVarBasic(scip, &matrixvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, matrixvar) );
            SCIPdebugMsg(scip, "Added new matrix variable <%s>\n", name );

            consvals[nconsvars] = -1.0;
            consvars[nconsvars++] = matrixvar;

            /* compute lhs and rhs */
            lhs = constmatrix[s * blocksize + t];

            /* create linear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "C%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, lhs, lhs,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added linear constraint for coupling the new matrix variable: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            assert( matrixvar != NULL );

            /* store variable (variable is captured) */
            matrixvars[s][t] = matrixvar;

            /* skip diagonal entries */
            if ( s == t )
               continue;

            /* skip if any of the variables is trivial */
            if ( matrixvars[s][s] == NULL || matrixvars[t][t] == NULL || matrixvars[s][t] == NULL )
               continue;

            /* create variable for sum of two matrix variables (note they they are nonnegative) */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sum%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateVarBasic(scip, &matrixsumvar, name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, matrixsumvar) );
            SCIPdebugMsg(scip, "Added new variable <%s> for sum of two matrix variables.\n", name);

            consvars[0] = matrixvars[s][s];
            consvars[1] = matrixvars[t][t];
            consvars[2] = matrixsumvar;

            consvals[0] = 1.0;
            consvals[1] = 1.0;
            consvals[2] = -1.0;

            /* create linear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Sum%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 3, consvars, consvals, 0.0, 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added linear constraint for sum of two matrix variables: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            /* create variable for difference of two matrix variablees */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "diff%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateVarBasic(scip, &matrixdiffvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, matrixdiffvar) );
            SCIPdebugMsg(scip, "Added new variable <%s> for difference of two matrix variables.\n", name);

            consvars[0] = matrixvars[s][s];
            consvars[1] = matrixvars[t][t];
            consvars[2] = matrixdiffvar;

            consvals[0] = 1.0;
            consvals[1] = -1.0;
            consvals[2] = -1.0;

            /* create linear constraint */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "Diff%d#%d#%d", c, s, t);
            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, 3, consvars, consvals, 0.0, 0.0,
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added linear constraint for difference of two matrix variables: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif

            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            /* construct SOC constraint */
            consvars[0] = matrixvars[s][t];
            consvars[1] = matrixdiffvar;

            consvals[0] = 2.0;
            consvals[1] = 1.0;

            /* add SOC constraint (if not solving LPs only propagate) */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "2x2minorSOC#%d#%d#%d", c, s, t);

#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
            consvals[0] = 4.0; /* correct to 4, because coefficient is outside of square */
            consvars[2] = matrixsumvar;
            consvals[2] = -1.0;
            SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, &cons, name, 0, NULL, NULL, 3, consvars, consvars, consvals,
                  - SCIPinfinity(scip), 0.0, TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE) );
#else
            SCIP_CALL( SCIPcreateConsSOC(scip, &cons, name, 2, consvars, consvals, NULL, 0.0, matrixsumvar, 1.0, 0.0,
                  TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE) );
#endif
            SCIP_CALL( SCIPaddCons(scip, cons) );

#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added 2x2 minor linear constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            ++(*naddconss);

            SCIP_CALL( SCIPreleaseVar(scip, &matrixsumvar) );
            SCIP_CALL( SCIPreleaseVar(scip, &matrixdiffvar) );
         }
      }

      for (s = blocksize - 1; s >= 0; --s)
      {
         for (t = 0; t <= s; ++t)
         {
            if ( matrixvars[s][t] != NULL )
            {
               SCIP_CALL( SCIPreleaseVar(scip, &matrixvars[s][t]) );
            }
         }
         SCIPfreeBufferArray(scip, &matrixvars[s]);
      }
      SCIPfreeBufferArray(scip, &matrixvars);

      for (i = 0; i < nvars; ++i)
         SCIPfreeBufferArray(scip, &matrices[i]);
      SCIPfreeBufferArray(scip, &matrices);
      SCIPfreeBufferArray(scip, &constmatrix);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** Presolve-routine that adds a linear cut arising from 2 by 2 minor inequalities \f$\sum_{i=1}^m (A_i)_{st} y_i \geq (A_0)_{st} - \sqrt{(A_0)_{ss} (A_0)_{tt}}\f$,
 *  if \f$(A_i)_{ss} = (A_i)_{tt} = 0\f$ and \f$(A_0)_{ss} (A_0)_{tt} > 0\f$ for all \f$i,s,t\f$.
 *
 *  See the dissertation of T. Gally, page 150.
 */
static
SCIP_RETCODE addTwoMinorProdConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             solvesdps,          /**< are we solving SDPs or LPs? */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** consvars;
   SCIP_Real* consvals;
   int blocksize;
   int nvars;
   int c;
   int i;
   int j;

   assert( scip != NULL );
   assert( naddconss != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_Real** matrices = NULL;
      SCIP_Real* constmatrix;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

      /* loop over all nonzero matrix entries in the constant part */
      for (j = 0; j < consdata->constnnonz; ++j)
      {
         SCIP_Real val;
         SCIP_Real prod;
         int s;
         int t;

         s = consdata->constrow[j];
         t = consdata->constcol[j];

         /* skip diagonal entries */
         if ( s == t )
            continue;

         /* compute matrices if not yet done */
         if ( matrices == NULL )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &matrices, nvars) );
            for (i = 0; i < nvars; ++i)
            {
               SCIP_CALL( SCIPallocBufferArray(scip, &matrices[i], blocksize * blocksize) );
               SCIP_CALL( SCIPconsSdpGetFullAj(scip, conss[c], i, matrices[i]) );
            }
         }
         assert( matrices != NULL );

         /* check whether diagonal entries in the matrices are all 0 */
         for (i = 0; i < nvars; ++i)
         {
            if ( ! SCIPisZero(scip, matrices[i][s * blocksize + s]) )
               break;
            if ( ! SCIPisZero(scip, matrices[i][t * blocksize + t]) )
               break;
         }
         if ( i < nvars )
            continue;

         /* check if product is sufficiently positive */
         prod = constmatrix[s * blocksize + s] * constmatrix[t * blocksize + t];
         if ( SCIPisEfficacious(scip, prod) )
         {
            SCIP_CONS* cons;
            SCIP_Real activitylb = 0.0;
            SCIP_Real lhs;
            int nconsvars = 0;

            /* fill in constraint */
            for (i = 0; i < nvars; ++i)
            {
               val = matrices[i][s * blocksize + t];
               if ( ! SCIPisZero(scip, val) )
               {
                  consvals[nconsvars] = val;
                  consvars[nconsvars] = consdata->vars[i];
                  ++nconsvars;

                  /* compute lower bound on activity */
                  if ( val > 0 )
                     activitylb += val * SCIPvarGetLbGlobal(consdata->vars[i]);
                  else
                     activitylb += val * SCIPvarGetUbGlobal(consdata->vars[i]);
               }
            }

            lhs = consdata->constval[j] - sqrt(prod);

            /* only proceed if constraint is not redundant */
            if ( SCIPisGE(scip, activitylb, lhs) )
               continue;

            /* add linear constraint (if not solving LPs only propagate) */
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "2x2minorprod#%d#%d", s, t);

            if ( conshdlrdata->sdpconshdlrdata->presollinconssparam == 1 )
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, lhs, SCIPinfinity(scip),
                     TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/
            else
               SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, lhs, SCIPinfinity(scip),
                     FALSE, ! solvesdps, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added 2x2 minor product constraint: ");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );
            ++(*naddconss);
         }
      }

      if ( matrices != NULL )
      {
         for (i = nvars - 1; i >= 0; --i)
            SCIPfreeBufferArray(scip, &matrices[i]);
         SCIPfreeBufferArray(scip, &matrices);
      }
      SCIPfreeBufferArray(scip, &constmatrix);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** add variable bounds based on 2 by 2 minors
 *
 *  We generate inequalities of the following form:
 *  \f[
 *      2\, \tilde{U}_{st} A(y)_{st} - \tilde{U}_{tt} A(y)_{ss} \leq \tilde{U}_{st}^2.
 *  \f]
 *
 *  The details are explained in the presolving paper (see the top of the file).
 */
static
SCIP_RETCODE addTwoMinorVarBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             solvesdps,          /**< are we solving SDPs or LPs? */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( naddconss != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   for (c = 0; c < nconss; ++c)
   {
      char name[SCIP_MAXSTRLEN];
      SCIP_VAR** consvars;
      SCIP_Real* consvals;
      SCIP_CONSDATA* consdata;
      SCIP_Real* constmatrix;
      SCIP_Real** matrices;
      int blocksize;
      int nvars;
      int s;
      int t;
      int i;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      blocksize = consdata->blocksize;
      nvars = consdata->nvars;

      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2 * nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2 * nvars) );

      /* get matrices */
      SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );

      SCIP_CALL( SCIPallocBufferArray(scip, &matrices, nvars) );
      for (i = 0; i < nvars; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &matrices[i], blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullAj(scip, conss[c], i, matrices[i]) );
      }

      /* loop over all possible entries */
      for (s = 1; s < blocksize; ++s)
      {
         SCIP_Real ubs = 0.0;
         SCIP_Real val;
         SCIP_Real bound;
         int pos;

         /* compute upper bound */
         pos = s * blocksize + s;
         for (i = 0; i < nvars; ++i)
         {
            val = matrices[i][pos];
            if ( ! SCIPisZero(scip, val) )
            {
               if ( val > 0.0 )
                  bound = SCIPvarGetUbLocal(consdata->vars[i]);
               else
                  bound = SCIPvarGetLbLocal(consdata->vars[i]);

               if ( SCIPisInfinity(scip, REALABS(bound)) )
               {
                  ubs = SCIPinfinity(scip);
                  break;
               }

               ubs += val * bound;
            }
            assert( ! SCIPisInfinity(scip, ubs) );
         }
         if ( ! SCIPisInfinity(scip, ubs) )
            ubs -= constmatrix[pos];

         /* loop over all other entries */
         for (t = s-1; t >= 0; --t)
         {
            SCIP_CONS* cons;
            SCIP_Real ubt = 0.0;
            SCIP_Real ubst = 0.0;

            /* compute upper bound for off-diagonal */
            pos = s * blocksize + t;
            for (i = 0; i < nvars; ++i)
            {
               val = matrices[i][pos];
               if ( ! SCIPisZero(scip, val) )
               {
                  if ( val > 0.0 )
                     bound = SCIPvarGetUbLocal(consdata->vars[i]);
                  else
                     bound = SCIPvarGetLbLocal(consdata->vars[i]);

                  if ( SCIPisInfinity(scip, REALABS(bound)) )
                  {
                     ubst = SCIPinfinity(scip);
                     break;
                  }

                  ubst += val * bound;
               }
               assert( ! SCIPisInfinity(scip, ubst) );
            }
            if ( ! SCIPisInfinity(scip, ubst) )
               ubst -= constmatrix[pos];

            /* if the bound is zero, no interesting inequality arises */
            if ( SCIPisZero(scip, ubst) )
               continue;

            /* we need a finite bound on the off-diagonal */
            if ( SCIPisInfinity(scip, ubst) )
               continue;

            /* compute upper bound for diagonal */
            pos = t * blocksize + t;
            for (i = 0; i < nvars; ++i)
            {
               val = matrices[i][pos];
               if ( ! SCIPisZero(scip, val) )
               {
                  if ( val > 0.0 )
                     bound = SCIPvarGetUbLocal(consdata->vars[i]);
                  else
                     bound = SCIPvarGetLbLocal(consdata->vars[i]);

                  if ( SCIPisInfinity(scip, REALABS(bound)) )
                  {
                     ubt = SCIPinfinity(scip);
                     break;
                  }

                  ubt += val * bound;
               }
               assert( ! SCIPisInfinity(scip, ubt) );
            }
            if ( ! SCIPisInfinity(scip, ubt) )
               ubt -= constmatrix[pos];

            assert( ! SCIPisZero(scip, ubst) );
            assert( ! SCIPisInfinity(scip, ubst) );

            /* first type of constraint */
            if ( ! SCIPisInfinity(scip, ubt) )
            {
               SCIP_Real rhs;
               int nconsvars = 0;

               rhs = ubst * ubst;

               if ( ! SCIPisZero(scip, ubst) )
               {
                  for (i = 0; i < nvars; ++i)
                  {
                     val = matrices[i][s * blocksize + t];
                     if ( ! SCIPisZero(scip, val) )
                     {
                        consvals[nconsvars] = 2.0 * ubst * val;
                        consvars[nconsvars++] = consdata->vars[i];
                     }
                  }
                  rhs += 2.0 * ubst * constmatrix[s * blocksize + t];
               }

               if ( ! SCIPisZero(scip, ubt) )
               {
                  for (i = 0; i < nvars; ++i)
                  {
                     val = matrices[i][s * blocksize + s];
                     if ( ! SCIPisZero(scip, val) )
                     {
                        consvals[nconsvars] = - ubt * val;
                        consvars[nconsvars++] = consdata->vars[i];
                     }
                  }
                  rhs -= ubt * constmatrix[s * blocksize + s];
               }

               if ( nconsvars >= 1 )
               {
                  /* add linear constraint (if not solving LPs only propagate) */
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "twominorvarbounda#%d#%d", s, t);

                  if ( conshdlrdata->sdpconshdlrdata->presollinconssparam == 1 )
                     SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
                           TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/
                  else
                     SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
                           FALSE, ! solvesdps, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

                  SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "\n");
#endif
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                  ++(*naddconss);
               }
            }

            /* second type of constraint */
            if ( ! SCIPisInfinity(scip, ubs) )
            {
               SCIP_Real rhs;
               int nconsvars = 0;

               rhs = ubst * ubst;

               if ( ! SCIPisZero(scip, ubst) )
               {
                  for (i = 0; i < nvars; ++i)
                  {
                     val = matrices[i][s * blocksize + t];
                     if ( ! SCIPisZero(scip, val) )
                     {
                        consvals[nconsvars] = 2.0 * ubst * val;
                        consvars[nconsvars++] = consdata->vars[i];
                     }
                  }
                  rhs += 2.0 * ubst * constmatrix[s * blocksize + t];
               }

               if ( ! SCIPisZero(scip, ubs) )
               {
                  for (i = 0; i < nvars; ++i)
                  {
                     val = matrices[i][t * blocksize + t];
                     if ( ! SCIPisZero(scip, val) )
                     {
                        consvals[nconsvars] = - ubs * val;
                        consvars[nconsvars++] = consdata->vars[i];
                     }
                  }
                  rhs -= ubs * constmatrix[t * blocksize + t];
               }

               if ( nconsvars >= 1 )
               {
                  /* add linear constraint (if not solving LPs only propagate) */
                  (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "twominorvarboundb#%d#%d", s, t);

                  if ( conshdlrdata->sdpconshdlrdata->presollinconssparam == 1 )
                     SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
                           TRUE, ! solvesdps, ! solvesdps, ! solvesdps, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) ); /*lint !e679*/
                  else
                     SCIP_CALL( SCIPcreateConsLinear(scip, &cons, name, nconsvars, consvars, consvals, -SCIPinfinity(scip), rhs,
                           FALSE, ! solvesdps, FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

                  SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
                  SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
                  SCIPinfoMessage(scip, NULL, "\n");
#endif
                  SCIP_CALL( SCIPreleaseCons(scip, &cons) );
                  ++(*naddconss);
               }
            }
         }
      }

      for (i = nvars - 1; i >= 0; --i)
         SCIPfreeBufferArray(scip, &matrices[i]);
      SCIPfreeBufferArray(scip, &matrices);
      SCIPfreeBufferArray(scip, &constmatrix);
      SCIPfreeBufferArray(scip, &consvals);
      SCIPfreeBufferArray(scip, &consvars);
   }

   return SCIP_OKAY;
}

/** add quadratic constraints to enforce rank-1 condition */
static
SCIP_RETCODE addRank1QuadConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddconss           /**< pointer to store how many constraints were added */
   )
{
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( naddconss != NULL );

   if ( conss == NULL )
      return SCIP_OKAY;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      consdata->maxevsubmat[0] = -1;
      consdata->maxevsubmat[1] = -1;

      /* For each constraint, if it should be rank one, add all quadratic constraints given by the 2x2 principal
       * minors. */
      if ( consdata->rankone && ! consdata->addedquadcons )
      {
         SCIP_VAR** quadvars1;
         SCIP_VAR** quadvars2;
         SCIP_VAR** linvars;
         SCIP_CONS* quadcons;
         SCIP_Real* lincoefs;
         SCIP_Real* quadcoefs;
         SCIP_Real* constmatrix;
         SCIP_Real** matrixAk;
         SCIP_Real lhs;
         SCIP_Real aiik;
         SCIP_Real ajjk;
         SCIP_Real aijk;
         SCIP_Real aiil;
         SCIP_Real ajjl;
         SCIP_Real aijl;
         SCIP_Real cii;
         SCIP_Real cjj;
         SCIP_Real cij;
         char name[SCIP_MAXSTRLEN];
         int* nnonzvars;
         int** nonzvars;
         int i;
         int j;
         int k;
         int l;
         int blocksize;
         int varind1;
         int varind2;

         blocksize = consdata->blocksize;

         SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, (blocksize * (blocksize + 1)) / 2) ); /*lint !e647*/
         SCIP_CALL( SCIPconsSdpGetLowerTriangConstMatrix(scip, conss[c], constmatrix) );

         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars1, consdata->nvars * consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadvars2, consdata->nvars * consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &linvars, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &quadcoefs, consdata->nvars * consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &lincoefs, consdata->nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &matrixAk, consdata->nvars) );

         for (i = 0; i < consdata->nvars; ++i)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &matrixAk[i], blocksize * blocksize) );
            SCIP_CALL( SCIPconsSdpGetFullAj(scip, conss[c], i, matrixAk[i]) );
         }

         SCIP_CALL( SCIPallocBufferArray(scip, &nnonzvars, (blocksize * (blocksize + 1)) / 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &nonzvars, (blocksize * (blocksize + 1)) / 2) );

         for (i = 0; i < blocksize; ++i)
         {
            for (j = 0; j <= i; ++j)
            {
               int varcnt = 0;

               SCIP_CALL( SCIPallocBufferArray(scip, &nonzvars[SCIPconsSdpCompLowerTriangPos(i,j)], consdata->nvars) );

               for (k = 0; k < consdata->nvars; ++k)
               {
                  if ( ! SCIPisZero(scip, matrixAk[k][i * blocksize + j]) || ! SCIPisZero(scip, matrixAk[k][i * blocksize + i]) || ! SCIPisZero(scip, matrixAk[k][j * blocksize + j]) )
                  {
                     nonzvars[SCIPconsSdpCompLowerTriangPos(i,j)][varcnt] = k;
                     varcnt++;
                  }
               }
               nnonzvars[SCIPconsSdpCompLowerTriangPos(i,j)] = varcnt;
            }
         }

         for (i = 0; i < blocksize; ++i)
         {
            for (j = 0; j < i; ++j)
            {
               int lincnt = 0;
               int quadcnt = 0;

               cii = constmatrix[SCIPconsSdpCompLowerTriangPos(i,i)];
               cjj = constmatrix[SCIPconsSdpCompLowerTriangPos(j,j)];
               cij = constmatrix[SCIPconsSdpCompLowerTriangPos(i,j)];

               for (k = 0; k < nnonzvars[SCIPconsSdpCompLowerTriangPos(i,j)]; ++k)
               {
                  varind1 = nonzvars[SCIPconsSdpCompLowerTriangPos(i,j)][k];
                  ajjk = matrixAk[varind1][j * consdata->blocksize + j];
                  aiik = matrixAk[varind1][i * consdata->blocksize + i];
                  aijk = matrixAk[varind1][j * consdata->blocksize + i];

                  if ( ! SCIPisZero(scip, -cii * ajjk - cjj * aiik + cij * aijk) )
                  {
                     linvars[lincnt] = consdata->vars[varind1];
                     lincoefs[lincnt] = -cii * ajjk - cjj * aiik + cij * aijk;
                     ++lincnt;
                  }

                  for (l = 0; l < k; ++l)
                  {
                     varind2 = nonzvars[SCIPconsSdpCompLowerTriangPos(i,j)][l];
                     ajjl = matrixAk[varind2][j * consdata->blocksize + j];
                     aiil = matrixAk[varind2][i * consdata->blocksize + i];
                     aijl = matrixAk[varind2][j * consdata->blocksize + i];

                     if ( ! SCIPisZero(scip, aiik * ajjl + ajjk * aiil - 2 * aijk * aijl) )
                     {
                        quadvars1[quadcnt] = consdata->vars[varind1];
                        quadvars2[quadcnt] = consdata->vars[varind2];
                        quadcoefs[quadcnt] = aiik * ajjl + ajjk * aiil - 2 * aijk * aijl;
                        ++quadcnt;
                     }
                  }

                  /* case l == k needs special treatment */
                  if ( ! SCIPisZero(scip, aiik * ajjk - aijk * aijk) )
                  {
                     quadvars1[quadcnt] = consdata->vars[varind1];
                     quadvars2[quadcnt] = consdata->vars[varind1];
                     quadcoefs[quadcnt] = aiik * ajjk - aijk * aijk;
                     ++quadcnt;
                  }
               }
               assert( quadcnt <= consdata->nvars * consdata->nvars );
               assert( lincnt <= consdata->nvars );

               lhs = cij * cij - cii * cjj;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quadcons#%d#%d#%d", i, j, c);

               /* create quadratic constraint */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
               SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, &quadcons, name, lincnt, linvars, lincoefs, quadcnt, quadvars1, quadvars2, quadcoefs, lhs, lhs,
                     TRUE,      /* initial */
                     TRUE,      /* separate */
                     TRUE,      /* enforce */
                     TRUE,      /* check */
                     TRUE,      /* propagate */
                     FALSE,     /* local */
                     FALSE,     /* modifiable */
                     FALSE,     /* dynamic */
                     TRUE) );   /* removable */
#else
               SCIP_CALL( SCIPcreateConsQuadratic(scip, &quadcons, name, lincnt, linvars, lincoefs, quadcnt, quadvars1, quadvars2, quadcoefs, lhs, lhs,
                     TRUE,      /* initial */
                     TRUE,      /* separate */
                     TRUE,      /* enforce */
                     TRUE,      /* check */
                     TRUE,      /* propagate */
                     FALSE,     /* local */
                     FALSE,     /* modifiable */
                     FALSE,     /* dynamic */
                     TRUE) );   /* removable */
#endif

#ifdef SCIP_MORE_DEBUG
               SCIP_CALL( SCIPprintCons(scip, quadcons, NULL) );
               SCIPinfoMessage(scip, NULL, "\n");
#endif

               SCIP_CALL( SCIPaddCons(scip, quadcons) );
               SCIP_CALL( SCIPreleaseCons(scip, &quadcons) );
               ++(*naddconss);
            }
         }

         for (i = blocksize - 1; i >= 0; --i)
         {
            for (j = i; j >= 0; --j)
               SCIPfreeBufferArray(scip, &nonzvars[SCIPconsSdpCompLowerTriangPos(i,j)]);
         }

         SCIPfreeBufferArray(scip, &nonzvars);
         SCIPfreeBufferArray(scip, &nnonzvars);

         for (i = consdata->nvars - 1; i >= 0; --i)
            SCIPfreeBufferArray(scip, &matrixAk[i]);

         SCIPfreeBufferArray(scip, &matrixAk);
         SCIPfreeBufferArray(scip, &lincoefs);
         SCIPfreeBufferArray(scip, &quadcoefs);
         SCIPfreeBufferArray(scip, &linvars);
         SCIPfreeBufferArray(scip, &quadvars2);
         SCIPfreeBufferArray(scip, &quadvars1);
         SCIPfreeBufferArray(scip, &constmatrix);
      }
      consdata->addedquadcons = TRUE;
   }

   return SCIP_OKAY;
}

/** checks quadratic constraints that will be added for a single rank-1 constraint  */
static
SCIP_RETCODE checkRank1QuadConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler */
   SCIP_CONS*            cons,               /**< the SDP constraint to check the rank for */
   SCIP_SOL*             sol,                /**< solution to check for rank one */
   SCIP_Bool             printreason,        /**< should the reason for the violation be printed? */
   SCIP_Bool*            isrankone           /**< pointer to return whether matrix is rank one */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* matrix;
   SCIP_Real submatrix[3];
   SCIP_Real tol;
   int blocksize;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( isrankone != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->rankone );
   assert( ! consdata->addedquadcons );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLRRANK1_NAME) == 0 );
   assert( conshdlrdata->sdpconshdlrdata->quadconsrank1 );

   blocksize = consdata->blocksize;
   *isrankone = TRUE;

   if ( conshdlrdata->sdpconshdlrdata->usedimacsfeastol )
   {
      assert( conshdlrdata->dimacsfeastol != SCIP_INVALID );
      tol = conshdlrdata->dimacsfeastol;
   }
   else
      tol = SCIPfeastol(scip);

   /* check quadratic constraints that will be added later */

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix, blocksize * (blocksize + 1)/2) ); /*lint !e647*/
   SCIP_CALL( computeSdpMatrix(scip, consdata, sol, matrix) );

   for (i = 0; i < blocksize && *isrankone; ++i)
   {
      for (j = 0; j < i && *isrankone; ++j)
      {
         SCIP_Real minor;

         submatrix[0] = matrix[SCIPconsSdpCompLowerTriangPos(i,i)];
         submatrix[1] = matrix[SCIPconsSdpCompLowerTriangPos(i,j)];
         submatrix[2] = matrix[SCIPconsSdpCompLowerTriangPos(j,j)];

         minor = submatrix[0] * submatrix[2] - submatrix[1] * submatrix[1];

         if ( minor < -tol || minor > tol )
         {
            *isrankone = FALSE;
            if ( printreason )
            {
               SCIPinfoMessage(scip, NULL, "SDPrank1-constraint <%s> is not rank1 (quadratic 2x2 minor for (%d,%d): %f).\n", SCIPconsGetName(cons), i, j, minor);
               SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            }
         }
         if ( sol != NULL )
            SCIPupdateSolConsViolation(scip, sol, REALABS(minor), REALABS(minor) / (1.0 + consdata->maxrhsentry));
      }
   }

   SCIPfreeBufferArray(scip, &matrix);

   return SCIP_OKAY;
}


/** detects if there are blocks with size one and transforms them to lp-rows */
static
SCIP_RETCODE move_1x1_blocks_to_lp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< array of constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  naddconss,          /**< pointer to store how many constraints were added */
   int*                  ndelconss,          /**< pointer to store how many constraints were deleted */
   int*                  nchgbds,            /**< pointer to store how many bounds were changed */
   SCIP_Bool*            infeasible          /**< pointer to store whether infeasibility was detected */
   )
{
   char cutname[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   SCIP_Real rhs;
   int nnonz;
   int nvars;
   int i;
   int j;
   int v;
#ifndef NDEBUG
   int snprintfreturn; /* used to assert the return code of snprintf */
#endif

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 || strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLRRANK1_NAME) == 0 );
   assert( infeasible != NULL );

   *infeasible = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   for (i = 0; i < nconss && !(*infeasible); ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );

      /* if there is a 1x1 SDP-Block */
      if ( consdata->blocksize == 1 )
      {
         int cnt = 0;

         nvars = consdata->nvars;
         nnonz = consdata->nnonz;
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, nnonz) );

         /* get all lhs-entries */
         for (v = 0; v < nvars; v++)
         {
            assert( consdata->nvarnonz[v] <= 1 ); /* in a 1x1 block there may be at most one entry per variable */

            for (j = 0; j < consdata->nvarnonz[v]; j++)
            {
               assert( consdata->col[v][j] == 0 && consdata->row[v][j] == 0 ); /* if the block has size 1, all entries should have row and col equal to 0 */
               if ( ! SCIPisZero(scip, consdata->val[v][j]) )
               {
                  coeffs[cnt] = consdata->val[v][j];
                  vars[cnt++] = consdata->vars[v];
               }
            }
         }

         /* get rhs */
         assert( consdata->constnnonz <= 1 ); /* the 1x1 constant matrix may only have one entry */

         rhs = (consdata->constnnonz == 1) ? consdata->constval[0] : 0.0; /* if this one entry is not 0, than this is the rhs, otherwise it's 0 */

         /* if there is more than one nonzero left-hand-side-entry, add a linear constraint, otherwise update the variable bound */
         if ( cnt > 1 )
         {
            /* add new linear cons */
#ifndef NDEBUG
            snprintfreturn = SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
            assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether name fits into string */
#else
            (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "1x1block_%d", ++(conshdlrdata->n1x1blocks));
#endif

            SCIP_CALL( SCIPcreateConsLinear(scip, &cons, cutname, cnt, vars, coeffs, rhs, SCIPinfinity(scip),
                  TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE, FALSE) );

            SCIP_CALL( SCIPaddCons(scip, cons) );
#ifdef SCIP_MORE_DEBUG
            SCIPinfoMessage(scip, NULL, "Added lp-constraint:\n");
            SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
            SCIPinfoMessage(scip, NULL, "\n");
#endif
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            (*naddconss)++;
         }
         else if ( cnt == 1 )
         {
            SCIP_Bool tightened;

            /* try to tighten bound */
            if ( SCIPisPositive(scip, coeffs[0]) )
            {
               SCIP_CALL( SCIPtightenVarLb(scip, vars[0], rhs / coeffs[0], FALSE, infeasible, &tightened) );
               if ( tightened )
               {
                  SCIPdebugMsg(scip, "Tightend lower bound of <%s> to %g because of diagonal values of SDP-constraint %s!\n",
                     SCIPvarGetName(vars[0]), SCIPvarGetLbGlobal(vars[0]), SCIPconsGetName(conss[i]));
                  ++(*nchgbds);
               }
            }
            else
            {
               assert( SCIPisNegative(scip, coeffs[0]) );
               SCIP_CALL( SCIPtightenVarUb(scip, vars[0], rhs / coeffs[0], FALSE, infeasible, &tightened) );
               if ( tightened )
               {
                  SCIPdebugMsg(scip, "Tightend upper bound of <%s> to %g because of diagonal values of SDP-constraint %s!\n",
                     SCIPvarGetName(vars[0]), SCIPvarGetUbGlobal(vars[0]), SCIPconsGetName(conss[i]));
                  ++(*nchgbds);
               }
            }
         }
         else
         {
            assert( cnt == 0 );
            SCIPdebugMsg(scip, "Detected 1x1 SDP-block without any nonzero coefficients \n");
            if ( SCIPisFeasGT(scip, rhs, 0.0) )
            {
               SCIPdebugMsg(scip, "Detected infeasibility in 1x1 SDP-block without any nonzero coefficients but with strictly positive rhs\n");
               *infeasible = TRUE;
            }
         }

         /* delete old 1x1 sdpcone */
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
         (*ndelconss)++;

         SCIPfreeBufferArray(scip, &coeffs);
         SCIPfreeBufferArray(scip, &vars);
      }
   }
   return SCIP_OKAY;
}

/** unlock variable */
static
SCIP_RETCODE unlockVar(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSDATA*        consdata,           /**< data of constraint */
   int                   v                   /**< index of variable */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );
   assert( 0 <= v && v < consdata->nvars );

   if ( consdata->locks != NULL )
   {
      assert( consdata->locks[v] == -2 || consdata->locks[v] == -1 || consdata->locks[v] == 0 || consdata->locks[v] == 1 );

      if ( consdata->locks[v] == 1 )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, 0, -1) );
      }
      else if ( consdata->locks[v] == - 1 )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, -1, 0) );
      }
      else if ( consdata->locks[v] == 0 )
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, -1, -1) );
      }
   }

   return SCIP_OKAY;
}

/** update locks of variable after aggregation */
static
SCIP_RETCODE updateVarLocks(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint */
   int                   v                   /**< index of variable */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real eigenvalue;
   SCIP_Real* Aj;
   int blocksize;
   int newlock = -2;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->locks != NULL );
   assert( 0 <= v && v < consdata->nvars );

   /* rank-1 constraints are always up- and down-locked */
   if ( consdata->rankone )
   {
      SCIP_CALL( unlockVar(scip, consdata, v) );
      consdata->locks[v] = 0;
      SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, 1, 1) );
      consdata->initallmatricespsd = FALSE; /* needs to be recomputed */
      return SCIP_OKAY;
   }

   blocksize = consdata->blocksize;
   SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) );

   SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, v, Aj) );

   /* compute new lock as in consLockSdp */
   SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
   if ( SCIPisNegative(scip, eigenvalue) )
   {
      newlock = 1;  /* up-lock */
      consdata->allmatricespsd = FALSE;
   }

   /* @todo check whether one can set allmatrices to true */
   if ( SCIPisPositive(scip, eigenvalue) )
      newlock = -1; /* down-lock */
   else
   {
      consdata->allmatricespsd = FALSE;
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
      if ( SCIPisPositive(scip, eigenvalue) )
      {
         if ( newlock == 1 )
            newlock = 0; /* up- and down-lock */
         else
            newlock = -1; /* down-lock */
      }
   }

   SCIPfreeBufferArray(scip, &Aj);

   /* if new lock is not equal to the old one, unlock variable and add new locks */
   if ( newlock != consdata->locks[v] )
   {
      SCIP_CALL( unlockVar(scip, consdata, v) );
      if ( newlock == 1 )  /* up-lock */
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, 0, 1) );
      }
      else if ( newlock == -1 )  /* down-lock */
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, 1, 0) );
      }
      else if ( newlock == 0 )  /* up and down lock */
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], SCIP_LOCKTYPE_MODEL, 1, 1) );
      }
      else
         assert( newlock == -2 );

      consdata->locks[v] = newlock;
   }

   return SCIP_OKAY;
}

#ifndef NDEBUG
/** check whether variable locks are correctly set */
static
SCIP_RETCODE checkVarsLocks(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons                /**< constraint */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* Aj;
   int blocksize;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->locks != NULL );

   /* rank-1 constraints should always be up- and down-locked */
   if ( consdata->rankone )
   {
      for (v = 0; v < consdata->nvars; ++v)
         assert( consdata->locks[v] == 0 );
      return SCIP_OKAY;
   }

   blocksize = consdata->blocksize;
   SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) );

   for (v = 0; v < consdata->nvars; ++v)
   {
      SCIP_Real eigenvalue;
      int newlock = -2;

      SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, v, Aj) );

      /* compute new lock as in consLockSdp */
      SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
      if ( SCIPisNegative(scip, eigenvalue) )
         newlock = 1;  /* up-lock */

      if ( SCIPisPositive(scip, eigenvalue) )
         newlock = -1; /* down-lock */
      else
      {
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
         if ( SCIPisPositive(scip, eigenvalue) )
         {
            if ( newlock == 1 )
               newlock = 0; /* up- and down-lock */
            else
               newlock = -1; /* down-lock */
         }
      }

      assert( newlock == consdata->locks[v] );
   }

   SCIPfreeBufferArray(scip, &Aj);

   return SCIP_OKAY;
}
#endif


/* defines for lexicographic sorting macros */
#define SORTTPL_NAMEEXT     RealRealIntInt
#define SORTTPL_KEYTYPE     SCIP_Real
#define SORTTPL_FIELD1TYPE  int
#define SORTTPL_FIELD2TYPE  int
#include "sorttpllex.c" /*lint !e451*/

/** check which matrices are unique across all SDP constraints */
static
SCIP_RETCODE checkSymUniqueMatrices(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   nconss,             /**< number of SDP constraints */
   SCIP_CONS**           conss               /**< SDP constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* minvals;
   SCIP_Real* maxvals;
   int* nnonz;
   int* considx;
   int* varidx;
   int nvars;
   int nmatrices = 0;
   int lastidx;
   int currentidx;
   int currentnnonz;
   int c;
   int v;
   int i;
   int j;

   assert( scip != NULL );

   if ( conss == NULL || nconss <= 0 )
      return SCIP_OKAY;
   assert( conss != NULL );

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &nnonz, nconss * nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &considx, nconss * nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varidx, nconss * nvars) );

   /* collect number of nonzeros of all matrices */
   for (c = 0; c < nconss; ++c)
   {
      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* init array */
      if ( consdata->issymunique == NULL )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->issymunique, consdata->nvars) );
      }
      else
         continue; /* skip constraints that have been considered before */

      /* collect data */
      for (v = 0; v < consdata->nvars; ++v)
      {
         consdata->issymunique[v] = FALSE;

         nnonz[nmatrices] = consdata->nvarnonz[v];
         considx[nmatrices] = c;
         varidx[nmatrices] = v;
         ++nmatrices;
      }
   }
   assert( nmatrices <= nconss * nvars );

   /* possibly stop early */
   if ( nmatrices == 0 )
   {
      SCIPfreeBufferArray(scip, &varidx);
      SCIPfreeBufferArray(scip, &considx);
      SCIPfreeBufferArray(scip, &nnonz);

      return SCIP_OKAY;
   }

   /* sort according to number of nonzeros */
   SCIPsortIntIntInt(nnonz, considx, varidx, nmatrices);

   SCIP_CALL( SCIPallocBufferArray(scip, &minvals, nmatrices) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxvals, nmatrices) );

   /* search matrices of equal number of nonzeros and compare w.r.t. min/max values */
   currentidx = 0;
   while ( currentidx < nmatrices )
   {
      /* find next matrix with different number of nonzeros */
      lastidx = currentidx;
      currentnnonz = nnonz[lastidx];
      do
      {
         ++currentidx;
      }
      while ( currentidx < nmatrices && nnonz[currentidx] == currentnnonz );

      /* if there are at least two matrices of the same size */
      if ( currentidx > lastidx + 1 )
      {
         int nmat = 0;

         /* loop through matrices and compute min/max values */
         for (i = lastidx; i < currentidx; ++i)
         {
            SCIP_Real minval = SCIP_REAL_MAX;
            SCIP_Real maxval = SCIP_REAL_MIN;

            /* compute min/max values */
            c = considx[i];
            v = varidx[i];
            assert( 0 <= c && c < nconss );
            assert( conss[c] != NULL );
            consdata = SCIPconsGetData(conss[c]);
            assert( consdata != NULL );

            assert( 0 <= v && v < consdata->nvars );
            for (j = 0; j < consdata->nvarnonz[v]; ++j)
            {
               SCIP_Real val;

               val = consdata->val[v][j];
               if ( val < minval )
                  minval = val;
               if ( val > maxval )
                  maxval = val;
            }
            SCIPdebugMsg(scip, "Constraint %d, variable %d: nnonz = %d, minval = %g, maxval = %g\n", c, v, consdata->nvarnonz[v], minval, maxval);

            minvals[nmat] = minval;
            maxvals[nmat] = maxval;
            ++nmat;
         }
         assert( nmat >= 2 );

         /* lexicographically sort entries according to minvals and maxvals */
         SCIPlexSortRealRealIntInt(minvals, maxvals, &considx[lastidx], &varidx[lastidx], nmat);

         /* loop through entries and check for different entries */
         i = 0;
         while ( i < nmat )
         {
            int len = 0;

            while ( i < nmat - 1 && SCIPisEQ(scip, minvals[i], minvals[i+1]) && SCIPisEQ(scip, maxvals[i], maxvals[i+1]) )
            {
               ++i;
               ++len;
            }

            /* if the entry is unique */
            if ( len == 0 )
            {
               /* matrix is unique */
               c = considx[lastidx + i];
               v = varidx[lastidx + i];
               assert( 0 <= c && c < nconss );
               assert( conss[c] != NULL );
               consdata = SCIPconsGetData(conss[c]);
               assert( consdata != NULL );
               consdata->issymunique[v] = TRUE;
               SCIPdebugMsg(scip, "unique: constraint %d, variable %d\n", c, v);
            }
            ++i;
         }
      }
      else
      {
         /* matrix is unique */
         c = considx[lastidx];
         v = varidx[lastidx];
         assert( 0 <= c && c < nconss );
         assert( conss[c] != NULL );
         consdata = SCIPconsGetData(conss[c]);
         assert( consdata != NULL );
         consdata->issymunique[v] = TRUE;
         SCIPdebugMsg(scip, "unique: constraint %d, variable %d\n", c, v);
      }
   }
   SCIPfreeBufferArray(scip, &maxvals);
   SCIPfreeBufferArray(scip, &minvals);
   SCIPfreeBufferArray(scip, &varidx);
   SCIPfreeBufferArray(scip, &considx);
   SCIPfreeBufferArray(scip, &nnonz);

   return SCIP_OKAY;
}


/** local function to perform (parts of) multiaggregation of a single variable within fixAndAggrVars */
static
SCIP_RETCODE multiaggrVar(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONS*            cons,               /**< constraint to multiaggregate for */
   int                   v,                  /**< position of the variable that gets (multi-)aggregated */
   SCIP_VAR**            aggrvars,           /**< variables this has to be (multi-)aggregated to */
   SCIP_Real*            scalars,            /**< scalar parts to multiply with for each variable this is aggregated to */
   int                   naggrvars,          /**< number of variables this is (multi-)aggregated to */
   SCIP_Real             constant,           /**< the constant part for the (multi-)aggregation */
   int*                  savedcol,           /**< array of columns for nonzeros that need to be added to the constant matrix */
   int*                  savedrow,           /**< array of rows for nonzeros that need to be added to the constant matrix */
   SCIP_Real*            savedval,           /**< array of values for nonzeros that need to be added to the constant matrix */
   int*                  nfixednonz,         /**< length of the arrays of saved nonzeros for the constant matrix */
   int*                  vararraylength      /**< length of the variable array */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real* vals;
   int* rows;
   int* cols;
   int nvarnonz;
   int aggrind;
   int aggrtargetlength;
   int globalnvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( scalars != NULL );
   assert( naggrvars > 0 );
   assert( savedcol != NULL );
   assert( savedrow != NULL );
   assert( savedval != NULL );
   assert( nfixednonz != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->locks != NULL );
   assert( 0 <= v && v < consdata->nvars );

   /* unlock variable */
   SCIP_CALL( unlockVar(scip, consdata, v) );

   /* save matrix of variable v (will be freed later) */
   rows = consdata->row[v];
   cols = consdata->col[v];
   vals = consdata->val[v];
   nvarnonz = consdata->nvarnonz[v];

   /* Sort matrix of variable v by nondecreasing row and then col to make merging below faster. */
   SCIPsdpVarfixerSortRowCol(rows, cols, vals, nvarnonz);

   /* fill the empty spot of the (multi-)aggregated variable with the last variable of this constraint (since variables do not have to be sorted) */
   SCIP_CALL( SCIPreleaseVar(scip, &consdata->vars[v]) );
   consdata->col[v] = consdata->col[consdata->nvars - 1];
   consdata->row[v] = consdata->row[consdata->nvars - 1];
   consdata->val[v] = consdata->val[consdata->nvars - 1];
   consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
   consdata->vars[v] = consdata->vars[consdata->nvars - 1];
   consdata->locks[v] = consdata->locks[consdata->nvars - 1];

   /* free issymunique, because it is invalid */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->issymunique, consdata->nvars);

   (consdata->nvars)--;

   /* iterate over all variables that variable v was aggregated to and insert the corresponding nonzeros */
   for (aggrind = 0; aggrind < naggrvars; aggrind++)
   {
      int aggrconsind = -1;

      assert( ! SCIPisZero(scip, scalars[aggrind]) );

      /* check if the variable already exists in this block */
      for (i = 0; i < consdata->nvars; i++)
      {
         if ( consdata->vars[i] == aggrvars[aggrind] )
         {
            aggrconsind = i;
            break;
         }
      }

      if ( aggrconsind > -1 )
      {
         /* if the variable to aggregate to is already part of this sdp-constraint, just add the nonzeros of the old variable to it */

         /* resize the arrays to the maximally needed length */
         aggrtargetlength = consdata->nvarnonz[aggrconsind] + nvarnonz;
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), consdata->nvarnonz[aggrconsind], aggrtargetlength) );

         /* merge: add scalar times matrix for variable v to aggregated matrix */
         SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), rows, cols, vals, nvarnonz, TRUE,
               scalars[aggrind], consdata->row[aggrconsind], consdata->col[aggrconsind],
               consdata->val[aggrconsind], &(consdata->nvarnonz[aggrconsind]), aggrtargetlength) );

         /* shrink them again if nonzeros could be combined */
         assert( consdata->nvarnonz[aggrconsind] <= aggrtargetlength );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->row[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->col[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->val[aggrconsind]), aggrtargetlength, consdata->nvarnonz[aggrconsind]) );

         SCIP_CALL( updateVarLocks(scip, cons, aggrconsind) );
      }
      else
      {
         int cnt = 0;

         /* the variable has to be added to this constraint */
         SCIPdebugMsg(scip, "adding variable %s to SDP constraint %s because of (multi-)aggregation\n", SCIPvarGetName(aggrvars[aggrind]), SCIPconsGetName(cons));

         /* check if we have to enlarge the arrays */
         if ( consdata->nvars == *vararraylength )
         {
            globalnvars = SCIPgetNVars(scip);

            /* we don't want to enlarge this by one for every variable added, so we immediately set it to the maximum possible size */
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, *vararraylength, globalnvars) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->locks, *vararraylength, globalnvars) );
            *vararraylength = globalnvars;
         }

         /* we insert this variable at the last position, as the ordering doesn't matter */
         SCIP_CALL( SCIPcaptureVar(scip, aggrvars[aggrind]) );
         consdata->vars[consdata->nvars] = aggrvars[aggrind];

         /* as there were no nonzeros thus far, we can just duplicate the saved arrays to get the nonzeros for the new variable */
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), cols, nvarnonz) );
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), rows, nvarnonz) );

         /* we have to multiply all entries by scalar before inserting them */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), nvarnonz) );
         for (i = 0; i < nvarnonz; i++)
         {
            /* if both scalar and savedval are small this might become too small */
            if ( ! SCIPisZero(scip, scalars[aggrind] * vals[i]) )
               consdata->val[consdata->nvars][cnt++] = scalars[aggrind] * vals[i];
         }
         consdata->nvarnonz[consdata->nvars] = cnt;

         /* free issymunique, because it is invalid */
         SCIPfreeBlockMemoryArrayNull(scip, &consdata->issymunique, consdata->nvars);

         consdata->locks[consdata->nvars] = -2;
         consdata->nvars++;
         SCIP_CALL( updateVarLocks(scip, cons, consdata->nvars-1) );
      }
   }

   /* if constant is not 0, insert entries for variable v into a long array of entries that is merged with constant matrix */
   if ( ! SCIPisZero(scip, constant) )
   {
      for (i = 0; i < nvarnonz; i++)
      {
         savedcol[*nfixednonz] = cols[i];
         savedrow[*nfixednonz] = rows[i];
         savedval[*nfixednonz] = vals[i] * constant; /* multiply with constant, since this is added to the constant matrix */
         (*nfixednonz)++;
      }
   }

   /* free the memory for the entries of the aggregated variable */
   SCIPfreeBlockMemoryArray(scip, &vals, nvarnonz);
   SCIPfreeBlockMemoryArray(scip, &rows, nvarnonz);
   SCIPfreeBlockMemoryArray(scip, &cols, nvarnonz);

#ifndef NDEBUG
   SCIP_CALL( checkVarsLocks(scip, cons) );
#endif

   return SCIP_OKAY;
}


/** presolve routine that looks through the data and handles fixed, (multi-)aggregated and negated variables */
static
SCIP_RETCODE fixAndAggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< array with constraints to check */
   int                   nconss,             /**< number of constraints to check */
   SCIP_Bool             aggregate           /**< do we want to (mutli-)aggregate variables ? */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int* savedcol;
   int* savedrow;
   SCIP_Real* savedval;
   int c;
   int v;
   int arraylength;
   SCIP_VAR* var;
   SCIP_VAR** aggrvars;
   SCIP_Real scalar;
   SCIP_Real* scalars;
   int naggrvars;
   SCIP_Real constant;
   int requiredsize;
   int globalnvars;
   int vararraylength;

   /* Loop over all variables once and collect all matrix entries that should be added to the constant matrix in
    * savedrow/col/val; this can then be merged with the constant matrix. */

   assert( scip != NULL );
   assert( conss != NULL );
   assert( nconss >= 0 );

   SCIPdebugMsg(scip, "Calling fixAndAggrVars with aggregate = %u.\n", aggregate);

   for (c = 0; c < nconss; ++c)
   {
      int nfixednonz = 0;

      assert( conss[c] != NULL );

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->locks != NULL );
      assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLR_NAME) == 0 );
      assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c])), CONSHDLRRANK1_NAME) == 0 );

      /* allocate memory to save nonzeros that need to be fixed */
      SCIP_CALL( SCIPallocBufferArray(scip, &savedcol, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedrow, consdata->nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &savedval, consdata->nnonz) );

      vararraylength = consdata->nvars;
      globalnvars = SCIPgetNVars(scip);

      for (v = 0; v < consdata->nvars; v++)/*lint --e{850}*/
      {
         SCIP_Bool negated = FALSE;

         /* if the variable is negated, get the negation var */
         if ( SCIPvarIsBinary(consdata->vars[v]) && SCIPvarIsNegated(consdata->vars[v]) )
         {
            negated = TRUE;
            var = SCIPvarGetNegationVar(consdata->vars[v]);
         }
         else
            var = consdata->vars[v];

         /* check if the variable is fixed in SCIP */
         if ( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED || SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) )
         {
            assert( SCIPisEQ(scip, SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var)) );

#ifdef SCIP_MORE_DEBUG
            SCIPdebugMsg(scip, "<%s>: Treating globally fixed variable %s with value %f!\n", SCIPconsGetName(conss[c]), SCIPvarGetName(var), SCIPvarGetLbGlobal(var));
#endif

            if ( (! negated && ! SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0)) || (negated && SCIPisEQ(scip, SCIPvarGetLbGlobal(var), 0.0)) )
            {
               /* the nonzeros are saved to later be inserted into the constant part (this is only done after all nonzeros of fixed variables have been
                * assembled, because we need to sort the constant nonzeros and loop over them, which we only want to do once and not once for each fixed
                * variable) */
               for (i = 0; i < consdata->nvarnonz[v]; i++)
               {
                  savedcol[nfixednonz] = consdata->col[v][i];
                  savedrow[nfixednonz] = consdata->row[v][i];

                  /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */
                  if ( ! negated )
                     savedval[nfixednonz] = consdata->val[v][i] * SCIPvarGetLbGlobal(var);
                  else
                     savedval[nfixednonz] = consdata->val[v][i]; /* if it is the negation of a variable fixed to zero, this variable is fixed to one */

                  nfixednonz++;
                  consdata->nnonz--;
               }
            }
            else
            {
               /* if the variable is fixed to zero, the nonzeros will just vanish, so we only reduce the number of nonzeros */
               consdata->nnonz -= consdata->nvarnonz[v];
            }

            /* free the memory of the corresponding entries in col/row/val */
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->val[v]), consdata->nvarnonz[v]);
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->row[v]), consdata->nvarnonz[v]);
            SCIPfreeBlockMemoryArrayNull(scip, &(consdata->col[v]), consdata->nvarnonz[v]);

            /* unlock variable */
            SCIP_CALL( unlockVar(scip, consdata, v) );

            /* as the variables don't need to be sorted, we just put the last variable into the empty spot and decrease sizes by one (at the end) */
            SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[v])) );
            if ( v < consdata->nvars - 1 )
            {
               consdata->col[v] = consdata->col[consdata->nvars - 1];
               consdata->row[v] = consdata->row[consdata->nvars - 1];
               consdata->val[v] = consdata->val[consdata->nvars - 1];
               consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
               consdata->vars[v] = consdata->vars[consdata->nvars - 1];
               consdata->locks[v] = consdata->locks[consdata->nvars - 1];
            }

            /* free issymunique, because it is invalid */
            SCIPfreeBlockMemoryArrayNull(scip, &consdata->issymunique, consdata->nvars);

            consdata->nvars--;
            v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
         }
         else if ( aggregate && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR) )
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, globalnvars) );
            SCIP_CALL( SCIPallocBufferArray(scip, &scalars, globalnvars) );

            /* this is how they should be initialized before calling SCIPgetProbvarLinearSum */
            if ( ! negated )
            {
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 0.0;
               scalars[0] = 1.0;
            }
            else
            {
               /* if this variable is the negation of var, than we look for a representation of 1.0 - var */
               aggrvars[0] = consdata->vars[v];
               naggrvars = 1;
               constant = 1.0;
               scalars[0] = -1.0;
            }

            /* get the variables this var was aggregated to */
            SCIP_CALL( SCIPgetProbvarLinearSum(scip, aggrvars, scalars, &naggrvars, globalnvars, &constant, &requiredsize, TRUE) );
            assert( requiredsize <= globalnvars ); /* requiredsize is the number of empty spots in aggrvars needed, globalnvars is the number
                                                    * of spots we provided */

            /* in the unlikely event that the multiaggregation reduces to 0 variables, the original one is fixed */
            if ( naggrvars == 0 )
            {
               assert( ! negated );
               if ( ! SCIPisZero(scip, constant) )
               {
                  for (i = 0; i < consdata->nvarnonz[v]; i++)
                  {
                     savedcol[nfixednonz] = consdata->col[v][i];
                     savedrow[nfixednonz] = consdata->row[v][i];

                     /* this is the final value to add, we no longer have to remember from which variable this comes, minus because we have +A_i but -A_0 */
                     savedval[nfixednonz] = consdata->val[v][i] * constant;

                     nfixednonz++;
                     consdata->nnonz--;
                  }
               }
               else
               {
                  /* if the variable is fixed to zero, the nonzeros will just vanish, so we only reduce the number of nonzeros */
                  consdata->nnonz -= consdata->nvarnonz[v];
               }

               /* free the memory of the corresponding entries in col/row/val */
               SCIPfreeBlockMemoryArrayNull(scip, &(consdata->val[v]), consdata->nvarnonz[v]);
               SCIPfreeBlockMemoryArrayNull(scip, &(consdata->row[v]), consdata->nvarnonz[v]);
               SCIPfreeBlockMemoryArrayNull(scip, &(consdata->col[v]), consdata->nvarnonz[v]);

               /* unlock variable */
               SCIP_CALL( unlockVar(scip, consdata, v) );

               /* as the variables don't need to be sorted, we just put the last variable into the empty spot and decrease sizes by one (at the end) */
               SCIP_CALL( SCIPreleaseVar(scip, &(consdata->vars[v])) );
               if ( v < consdata->nvars -1 )
               {
                  consdata->col[v] = consdata->col[consdata->nvars - 1];
                  consdata->row[v] = consdata->row[consdata->nvars - 1];
                  consdata->val[v] = consdata->val[consdata->nvars - 1];
                  consdata->nvarnonz[v] = consdata->nvarnonz[consdata->nvars - 1];
                  consdata->vars[v] = consdata->vars[consdata->nvars - 1];
                  consdata->locks[v] = consdata->locks[consdata->nvars - 1];
               }

               /* free issymunique, because it is invalid */
               SCIPfreeBlockMemoryArrayNull(scip, &consdata->issymunique, consdata->nvars);

               consdata->nvars--;
               v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
            }
            else
            {
               /* Debugmessages for the (multi-)aggregation */
#ifdef SCIP_DEBUG
               if ( SCIPvarGetStatus(consdata->vars[v]) == SCIP_VARSTATUS_AGGREGATED )
                  SCIPdebugMsg(scip, "aggregating variable %s to ", SCIPvarGetName(var));
               else
                  SCIPdebugMsg(scip, "multiaggregating variable %s to ", SCIPvarGetName(var));
               for (i = 0; i < naggrvars; i++)
                  SCIPdebugMessagePrint(scip, "+ %g %s ", scalars[i], SCIPvarGetName(aggrvars[i]));
               SCIPdebugMessagePrint(scip, "+ %g.\n", constant);
#endif

               /* add the nonzeros to the saved-arrays for the constant part, remove the nonzeros for the old variables and add them to the variables this variable
                * was (multi-)aggregated to */
               SCIP_CALL( multiaggrVar(scip, conss[c], v, aggrvars, scalars, naggrvars, constant, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );
               v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
            }

            SCIPfreeBufferArray(scip, &aggrvars);
            SCIPfreeBufferArray(scip, &scalars);
         }
         else if ( negated && (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN) && aggregate )
         {
             /* if var1 is the negation of var2, then this is equivalent to it being aggregated to -var2 + 1 = 1 - var2 */

            SCIPdebugMsg(scip, "Changing variable %s to negation of variable <%s>!\n", SCIPvarGetName(consdata->vars[v]), SCIPvarGetName(var));

            scalar = -1.0;

            SCIP_CALL( multiaggrVar(scip, conss[c], v, &var, &scalar, 1, 1.0, savedcol, savedrow, savedval, &nfixednonz, &vararraylength) );
            v--; /* we need to check again if the variable we just shifted to this position also needs to be fixed */
         }
      }

      /* shrink the variable arrays if they were enlarged too much (or more vars were removed than added) */
      assert( consdata->nvars <= vararraylength );
      if ( consdata->nvars < vararraylength )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, vararraylength, consdata->nvars) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->locks, vararraylength, consdata->nvars) );
      }

      /* allocate the maximally needed memory for inserting the fixed variables into the constant part */
      arraylength = consdata->constnnonz + nfixednonz;
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), consdata->constnnonz, arraylength) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), consdata->constnnonz, arraylength) );

      /* insert the fixed variables into the constant arrays, as we have +A_i but -A_0 we mutliply them by -1 */
      SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), savedrow, savedcol, savedval, nfixednonz, FALSE, -1.0, consdata->constrow,
            consdata->constcol, consdata->constval, &(consdata->constnnonz), arraylength) );

      assert( consdata->constnnonz <= arraylength ); /* the allocated memory should always be sufficient */

      /* shrink the arrays if nonzeros could be combined */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constcol), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constrow), arraylength, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(consdata->constval), arraylength, consdata->constnnonz) );

      /* free the saved arrays */
      SCIPfreeBufferArray(scip, &savedval);
      SCIPfreeBufferArray(scip, &savedrow);
      SCIPfreeBufferArray(scip, &savedcol);

      /* recompute sdpnnonz */
      consdata->nnonz = 0;
      for (v = 0; v < consdata->nvars; v++)
         consdata->nnonz += consdata->nvarnonz[v];

      /* possibly update issymunique */
      if ( consdata->issymunique == NULL )
      {
         SCIP_CALL( checkSymUniqueMatrices(scip, 1, &conss[c]) );
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflicting assignment on given constraint, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint */
   int                   diags,              /**< index for diagonal entry corresponding to s */
   int                   diagt,              /**< index for diagonal entry corresponding to t */
   int                   pos,                /**< index for off-diagonal entry corresponding to (s,t) */
   SCIP_Bool             upperbound,         /**< whether upper bound on pos caused infeasibility */
   SCIP_Bool             usepos              /**< whether the off-diagonal entry for (s,t) is necessary for analysis */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "Analyzing a conflict during propagation\n");

   /* conflict analysis can only be applied in solving stage and if it is applicable */
   if ( (SCIPgetStage(scip) != SCIP_STAGE_SOLVING && ! SCIPinProbing(scip)) || ! SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   assert( consdata->matrixvar != NULL );
   assert( consdata->matrixval != NULL );
   if ( consdata->matrixvar[diags] != NULL )
   {
      assert( consdata->matrixval[diags] != SCIP_INVALID );

      if ( consdata->matrixval[diags] > 0.0 )
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diags], NULL) );
      else
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diags], NULL) );
   }

   if ( consdata->matrixvar[diagt] != NULL )
   {
      assert( consdata->matrixval[diagt] != SCIP_INVALID );

      if ( consdata->matrixval[diagt] > 0.0 )
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diagt], NULL) );
      else
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diagt], NULL) );
   }

   if ( usepos )
   {
      assert( consdata->matrixvar[pos] != NULL);
      assert( consdata->matrixval[pos] != SCIP_INVALID);

      if ( upperbound )
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[pos], NULL) );
      else
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[pos], NULL) );
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, &success) );

   if ( success )
      SCIPdebugMsg(scip, "Succesfully analyzed and resolved conflict!\n");

   return SCIP_OKAY;
}

/** propagates upper bounds
 *
 *  The details are explained in the presolving paper (see the top of the file).
 */
static
SCIP_RETCODE propagateUpperBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility was detected */
   int*                  nprop,              /**< pointer to store the number of propagations performed */
   int*                  nintrnd             /**< pointer to store how many propagations of integer variables were rounded to be integral */
   )
{
   int c;

   assert( infeasible != NULL );
   assert( nprop != NULL );
   assert( nintrnd != NULL );

   *infeasible = FALSE;
   *nprop = 0;
   *nintrnd = 0;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      int blocksize;
      int i;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      blocksize = consdata->blocksize;

      /* build matrix if not yet done */
      SCIP_CALL( constructMatrixvar(scip, conss[c], consdata) );
      assert( consdata->matrixvar != NULL );
      assert( consdata->matrixval != NULL );
      assert( consdata->matrixconst != NULL );

      /* search for trace constraint */
      if ( consdata->tracebound < -1.5 )
      {
         SCIP_CONSHDLR* linconshdlr;
         linconshdlr = SCIPfindConshdlr(scip, "linear");
         if ( linconshdlr != NULL )
         {
            SCIP_CONS** linconss;
            int nlinconss;

            linconss = SCIPconshdlrGetConss(linconshdlr);
            nlinconss = SCIPconshdlrGetNConss(linconshdlr);

            for (i = 0; i < nlinconss; ++i)
            {
               SCIP_Real* linvals;
               SCIP_VAR** linvars;
               SCIP_Bool coefok = TRUE;
               int nlinvars;
               int j;

               nlinvars = SCIPgetNVarsLinear(scip, linconss[i]);

               /* if the number of variables is not equal to the dimension, we do not have a trace bound */
               if ( nlinvars != consdata->blocksize )
                  continue;

               linvars = SCIPgetVarsLinear(scip, linconss[i]);
               linvals = SCIPgetValsLinear(scip, linconss[i]);

               /* check whether all variables are diagonal entries */
               for (j = 0; j < nlinvars; ++j)
               {
                  SCIP_VAR* var;
                  int k;

                  if ( ! SCIPisEQ(scip, linvals[j], 1.0) )
                  {
                     coefok = FALSE;
                     break;
                  }

                  var = linvars[j];
                  for (k = 0; k < consdata->blocksize; ++k)
                  {
                     if ( consdata->matrixvar[k * (k + 1)/2 + k] == var )
                        break;
                  }

                  if ( k >= consdata->blocksize )
                  {
                     SCIPdebugMsg(scip, "Could not find variable <%s>.\n", SCIPvarGetName(linvars[j]));
                     break;
                  }
               }

               /* did not find variables or coefficients != 1 -> consider next linear constraint */
               if ( j < nlinvars || ! coefok )
                  continue;

               consdata->tracebound = SCIPgetRhsLinear(scip, linconss[i]);
               SCIPdebugMsg(scip, "Found tracebound constraint with bound = %g.\n", consdata->tracebound);
               break;
            }
         }
         if ( consdata->tracebound < -1.5 )
            consdata->tracebound = -1.0;
      }

      /* if there is at least one entry that only depends on a single variable */
      if ( consdata->propubpossible )
      {
         int s;
         int t;

         assert( consdata->nsingle > 0 );

         /* check all off-diagonal positions */
         for (s = 1; s < blocksize; ++s)
         {
            SCIP_VAR* vars;
            SCIP_Real ubs = 0.0;
            int diags;

            diags = s * (s + 1)/2 + s;
            if ( consdata->matrixval[diags] == SCIP_INVALID ) /*lint !e777*/
               continue;

            vars = consdata->matrixvar[diags];
            if ( vars != NULL )
            {
               if ( consdata->matrixval[diags] > 0.0 )
                  ubs = SCIPvarGetUbLocal(vars);
               else
                  ubs = SCIPvarGetLbLocal(vars);

               if ( SCIPisInfinity(scip, REALABS(ubs)) )
                  continue;

               ubs *= consdata->matrixval[diags];
            }
            assert( consdata->matrixconst[diags] != SCIP_INVALID );
            assert( ! SCIPisInfinity(scip, ubs) );

            ubs -= consdata->matrixconst[diags];

            for (t = 0; t < s; ++t)
            {
               SCIP_Bool tightened;
               SCIP_VAR* vart;
               SCIP_VAR* varst;
               SCIP_Real bound;
               SCIP_Real ubt = 0.0;
               int diagt;
               int pos;

               pos = s * (s + 1)/2 + t;
               varst = consdata->matrixvar[pos];
               if ( varst == NULL || ! SCIPvarIsActive(varst) )
                  continue;

               diagt = t * (t + 1)/2 + t;
               if ( consdata->matrixval[diagt] == SCIP_INVALID ) /*lint !e777*/
                  continue;

               vart = consdata->matrixvar[diagt];
               if ( vart != NULL )
               {
                  if ( consdata->matrixval[diagt] > 0.0 )
                     ubt = SCIPvarGetUbLocal(vart);
                  else
                     ubt = SCIPvarGetLbLocal(vart);

                  if ( SCIPisInfinity(scip, REALABS(ubt)) )
                     continue;

                  ubt *= consdata->matrixval[diagt];
               }
               assert( consdata->matrixconst[diagt] != SCIP_INVALID );
               assert( ! SCIPisInfinity(scip, ubt) );

               ubt -= consdata->matrixconst[diagt];

               if ( SCIPisFeasLT(scip, ubs, 0.0) || SCIPisFeasLT(scip, ubt, 0.0) )
               {
                  *infeasible = TRUE;
                  SCIP_CALL( analyzeConflict(scip, conss[c], diags, diagt, pos, TRUE, FALSE) );
                  return SCIP_OKAY;
               }

               /* compute upper bound without trace bound */
               if ( consdata->matrixval[pos] > 0.0 )
                  bound = (sqrt(ubs * ubt) + consdata->matrixconst[pos]) /  consdata->matrixval[pos];
               else
                  bound = (- sqrt(ubs * ubt) + consdata->matrixconst[pos]) /  consdata->matrixval[pos];

               /* check for stronger bound with trace bound */
               if ( consdata->tracebound > 0.0 )
               {
                  /* Note that tracebound is only computed for a primal SDP, thus it does not need to be retransformed
                   * to matrix pencil notation. */
                  if ( consdata->tracebound/2.0 < bound )
                     bound = consdata->tracebound/2.0;
               }

               assert( varst != NULL );
               if ( SCIPisFeasLT(scip, bound, SCIPvarGetUbLocal(varst)) )
               {
                  SCIP_CALL( SCIPinferVarUbCons(scip, varst, bound, conss[c], s * blocksize + t, FALSE, infeasible, &tightened) );
                  if ( *infeasible )
                  {
                     SCIPdebugMsg(scip, "Upper bound propagation detected infeasibility, call analyzeConfilct.\n");
                     SCIP_CALL( analyzeConflict(scip, conss[c], diags, diagt, pos, TRUE, TRUE) );
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                  {
                     SCIPdebugMsg(scip, "Upper bound propagation tightened bound of <%s> to %g.\n", SCIPvarGetName(varst), bound);
                     ++(*nprop);

                     /*  if variable is integral, the bound change should automatically produce an integer bound */
                     if ( SCIPvarIsIntegral(varst) && ! SCIPisFeasIntegral(scip, bound) )
                     {
                        assert( SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(varst)) );
                        ++(*nintrnd);
                     }
                  }
               }

               /* compute lower bound without trace bound */
               if ( consdata->matrixval[pos] > 0.0 )
                  bound = (- sqrt(ubs * ubt) + consdata->matrixconst[pos]) /  consdata->matrixval[pos];
               else
                  bound = (sqrt(ubs * ubt) + consdata->matrixconst[pos]) /  consdata->matrixval[pos];

               /* check for stronger bound with trace bound */
               if ( consdata->tracebound > 0.0 )
               {
                  if ( -consdata->tracebound/2.0 > bound )
                     bound = -consdata->tracebound/2.0;
               }

               if ( SCIPisFeasGT(scip, bound, SCIPvarGetLbLocal(varst)) )
               {
                  SCIP_CALL( SCIPinferVarLbCons(scip, varst, bound, conss[c], s * blocksize + t, FALSE, infeasible, &tightened) );
                  if ( *infeasible )
                  {
                     SCIPdebugMsg(scip, "Upper bound propagation detected infeasibility, call analyzeConfilct.\n");
                     SCIP_CALL( analyzeConflict(scip, conss[c], diags, diagt, pos, FALSE, TRUE) );
                     return SCIP_OKAY;
                  }
                  if ( tightened )
                  {
                     SCIPdebugMsg(scip, "Upper bound propagation tightened bound of <%s> to %g.\n", SCIPvarGetName(varst), bound);
                     ++(*nprop);

                     /*  if variable is integral, the bound change should automatically produce an integer bound */
                     if ( SCIPvarIsIntegral(varst) && ! SCIPisFeasIntegral(scip, bound) )
                     {
                        assert( SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(varst)) );
                        ++(*nintrnd);
                     }
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** analyzes conflicting assignment on given constraint from 3x3 minor propagation, and adds conflict constraint to problem */
static
SCIP_RETCODE analyzeConflict3Minor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint */
   int                   diagr,              /**< index for diagonal entry corresponding to r (or -1) */
   int                   diags,              /**< index for diagonal entry corresponding to s (or -1) */
   int                   posrs,              /**< index for off-diagonal entry corresponding to (r,s) */
   int                   pos1,               /**< index for one off-diagonal entry corresponding to fixed variable (or -1) */
   int                   pos2                /**< index for one off-diagonal entry corresponding to fixed variable (or -1) */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( posrs >= 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* conflict analysis can only be applied in solving stage */
   if ( SCIPgetStage(scip) != SCIP_STAGE_SOLVING && ! SCIPinProbing(scip) )
      return SCIP_OKAY;

   if ( ! SCIPisConflictAnalysisApplicable(scip) )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Analyzing a conflict during propagation of 3x3 minors ...\n");

   /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
   SCIP_CALL( SCIPinitConflictAnalysis(scip, SCIP_CONFTYPE_PROPAGATION, FALSE) );

   assert( consdata->matrixvar != NULL );
   assert( consdata->matrixval != NULL );

   if ( diagr >= 0 && consdata->matrixvar[diagr] != NULL )
   {
      assert( SCIPisFeasEQ(scip, consdata->matrixval[diagr] * SCIPvarGetLbLocal(consdata->matrixvar[diagr]) - consdata->matrixconst[diagr], 1.0) );
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->matrixvar[diagr]), SCIPvarGetUbLocal(consdata->matrixvar[diagr])) );

      if ( SCIPvarIsBinary(consdata->matrixvar[diagr]) )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->matrixvar[diagr]) );
      }
      else
      {
         /* add both bounds, because we do not know which bound caused the fixing */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diagr], NULL) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diagr], NULL) );
      }
   }

   if ( diags >= 0 && consdata->matrixvar[diags] != NULL )
   {
      assert( SCIPisFeasEQ(scip, consdata->matrixval[diags] * SCIPvarGetLbLocal(consdata->matrixvar[diags]) - consdata->matrixconst[diags], 1.0) );
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->matrixvar[diags]), SCIPvarGetUbLocal(consdata->matrixvar[diags])) );

      if ( SCIPvarIsBinary(consdata->matrixvar[diags]) )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->matrixvar[diags]) );
      }
      else
      {
         /* add both bounds, because we do not know which bound caused the fixing */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diags], NULL) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diags], NULL) );
      }
   }

   if ( consdata->matrixvar[posrs] != NULL )
   {
      assert( SCIPisFeasEQ(scip, consdata->matrixval[posrs] * SCIPvarGetLbLocal(consdata->matrixvar[posrs]) - consdata->matrixconst[posrs], 1.0) );
      assert( SCIPisEQ(scip, SCIPvarGetLbLocal(consdata->matrixvar[posrs]), SCIPvarGetUbLocal(consdata->matrixvar[posrs])) );

      if ( SCIPvarIsBinary(consdata->matrixvar[posrs]) )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->matrixvar[posrs]) );
      }
      else
      {
         /* add both bounds, because we do not know which bound cause the fixing */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[posrs], NULL) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[posrs], NULL) );
      }
   }

   if ( pos1 >= 0 && consdata->matrixvar[pos1] != NULL )
   {
      if ( SCIPvarIsBinary(consdata->matrixvar[pos1]) )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->matrixvar[pos1]) );
      }
      else
      {
         /* add both bounds, because we do not know which bound cause the fixing */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[pos1], NULL) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[pos1], NULL) );
      }
   }

   if ( pos2 >= 0 && consdata->matrixvar[pos2] != NULL )
   {
      if ( SCIPvarIsBinary(consdata->matrixvar[pos2]) )
      {
         SCIP_CALL( SCIPaddConflictBinvar(scip, consdata->matrixvar[pos2]) );
      }
      else
      {
         /* add both bounds, because we do not know which bound cause the fixing */
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[pos2], NULL) );
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[pos2], NULL) );
      }
   }

   /* analyze the conflict */
   SCIP_CALL( SCIPanalyzeConflictCons(scip, cons, &success) );

   if ( success )
      SCIPdebugMsg(scip, "Succesfully analyzed and resolved conflict!\n");

   return SCIP_OKAY;
}


/** propagates 3x3 minors
 *
 *  The idea is the following. If the diagonal entries of a 3x3 minor are fixed to 1 and one further entry is fixed to
 *  1, then the other two entries must be the same. For instance, assume the current 3x3 minor looks as follows:
 *  [1 1 a]
 *  [1 1 b]
 *  [a b c]
 *  Its determinant is c + 2 ab - a^2 - b^2 - c = - (a - b)^2. Since this determinant must be nonnegative for the
 *  complete matrix to be positive semidefinite, a = b follows. This argument can be repeatedly applied to show that the
 *  rows/column corresponding to the first two rows/columns must be equal.
 *
 *  The idea is motivated by the Masters Thesis "Facial Reduction on Binary Semidefinite Programs" by Jeremy Jany, TU
 *  Darmstadt, 2021.
 */
static
SCIP_RETCODE propagate3Minors(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool             nonconst3minors,    /**< Should 3x3 minors be propagated if the diagonal is not constant? */
   SCIP_Bool*            infeasible,         /**< pointer to store whether infeasibility was detected */
   int*                  nprop               /**< pointer to store the number of propagations performed */
   )
{
   int c;

   assert( infeasible != NULL );
   assert( nprop != NULL );

   *infeasible = FALSE;
   *nprop = 0;

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      int blocksize;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      blocksize = consdata->blocksize;

      /* build matrix if not yet done */
      SCIP_CALL( constructMatrixvar(scip, conss[c], consdata) );
      assert( consdata->matrixvar != NULL );
      assert( consdata->matrixval != NULL );
      assert( consdata->matrixconst != NULL );

      /* if there is at least one entry that only depends on a single variable */
      if ( consdata->nsingle > 0 )
      {
         SCIP_VAR* var;
         SCIP_Real val;
         int r;
         int s;
         int t;

         /* streamlined version if we know that the diagonals are fixed to be 1 */
         if ( consdata->diagconstantone )
         {
            /* check rows */
            for (r = 1; r < blocksize; ++r)
            {
               /* check column */
               for (s = 0; s < r; ++s)
               {
                  int posrs;

                  /* check whether position (r,s) is 1 */
                  posrs = r * (r + 1)/2 + s;

                  /* skip positions covered by at least two variables */
                  if ( consdata->matrixval[posrs] == SCIP_INVALID ) /*lint !e777*/
                     continue;

                  var = consdata->matrixvar[posrs];
                  if ( var != NULL )
                  {
                     /* skip unfixed variable */
                     if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                        continue;

                     val = SCIPvarGetLbLocal(var); /* fixed value */
                  }
                  else
                     val = 0.0;

                  /* check whether off-diagonal (r,s) is 1 */
                  if ( ! SCIPisFeasEQ(scip, consdata->matrixval[posrs] * val - consdata->matrixconst[posrs], 1.0) )
                     continue;

                  /* at this place rows/columns r and s are equal */

                  /* try to fix variables */
                  for (t = 1; t < blocksize; ++t)
                  {
                     SCIP_VAR* var1;
                     SCIP_VAR* var2;
                     int pos1;
                     int pos2;

                     if ( t == r || t == s )
                        continue;

                     if ( t > s )
                        pos1 = t * (t + 1)/2 + s;
                     else
                        pos1 = s * (s + 1)/2 + t;
                     var1 = consdata->matrixvar[pos1];
                     if ( var1 == NULL )
                        continue;

                     if ( t > r )
                        pos2 = t * (t + 1)/2 + r;
                     else
                        pos2 = r * (r + 1)/2 + t;
                     var2 = consdata->matrixvar[pos2];
                     if ( var2 == NULL )
                        continue;

                     /* if var1 is fixed */
                     if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) )
                     {
                        /* if var2 is also fixed */
                        if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
                        {
                           /* if the variables are fixed to different values, we are infeasible */
                           if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)) )
                           {
                              SCIPdebugMsg(scip, "Detected infeasibility for (%d, %d, %d) <%s>, <%s>.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              *infeasible = TRUE;
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], -1, -1, posrs, pos1, pos2) );
                              return SCIP_OKAY;
                           }
                        }
                        else
                        {
                           SCIP_Bool tightened;

                           /* fix var2 to the same value of var1 */
                           /* currently reverse propagation does not work for this case: use INT_MAX as inferinfo */
                           SCIP_CALL( SCIPinferVarFixCons(scip, var2, SCIPvarGetLbLocal(var1), conss[c], INT_MAX, FALSE, infeasible, &tightened) );
                           if ( *infeasible )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) <%s>, <%s> detected infeasibility.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], -1, -1, posrs, pos1, -1) );
                              return SCIP_OKAY;
                           }
                           if ( tightened )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) successfully tightened a bound of <%s> to %f.\n",
                                 r, s, t, SCIPvarGetName(var2), SCIPvarGetLbLocal(var1));
                              ++(*nprop);
                           }
                        }
                     }
                     else
                     {
                        /* if var2 is fixed (var1 is not fixed) */
                        if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
                        {
                           SCIP_Bool tightened;

                           /* fix var1 to the same value of var2 */
                           /* currently reverse propagation does not work for this case: use INT_MAX as inferinfo */
                           SCIP_CALL( SCIPinferVarFixCons(scip, var1, SCIPvarGetLbLocal(var2), conss[c], INT_MAX, FALSE, infeasible, &tightened) );
                           if ( *infeasible )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) <%s>, <%s> detected infeasibility.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], -1, -1, posrs, -1, pos2) );
                              return SCIP_OKAY;
                           }
                           if ( tightened )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) successfully tightened a bound of <%s> to %f.\n",
                                 r, s, t, SCIPvarGetName(var1), SCIPvarGetLbLocal(var2));
                              ++(*nprop);
                           }
                        }
                     }
                  }
               }
            }
         }
         else if ( nonconst3minors )
         {
            /* extended version */

            /* check rows */
            for (r = 1; r < blocksize; ++r)
            {
               int diagr;

               /* make sure that we have 1s on the diagonal */
               diagr = r * (r + 1)/2 + r;

               /* skip positions covered by at least two variables */
               if ( consdata->matrixval[diagr] == SCIP_INVALID ) /*lint !e777*/
                  continue;

               var = consdata->matrixvar[diagr];
               if ( var != NULL )
               {
                  /* skip unfixed variable */
                  if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                     continue;

                  val = SCIPvarGetLbLocal(var); /* fixed value */
               }
               else
                  val = 0.0;

               /* the result should be equal to 1 */
               if ( ! SCIPisFeasEQ(scip, consdata->matrixval[diagr] * val - consdata->matrixconst[diagr], 1.0) )
                  continue;

               /* check column */
               for (s = 0; s < r; ++s)
               {
                  int diags;
                  int posrs;

                  /* make sure that we have 1s on the diagonal */
                  diags = s * (s + 1)/2 + s;

                  /* skip positions covered by at least two variables */
                  if ( consdata->matrixval[diags] == SCIP_INVALID ) /*lint !e777*/
                     continue;

                  var = consdata->matrixvar[diags];
                  if ( var != NULL )
                  {
                     /* skip unfixed variable */
                     if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                        continue;

                     val = SCIPvarGetLbLocal(var); /* fixed value */
                  }
                  else
                     val = 0.0;

                  /* the result should be equal to 1 */
                  if ( ! SCIPisFeasEQ(scip, consdata->matrixval[diags] * val - consdata->matrixconst[diags], 1.0) )
                     continue;

                  /* check whether position (r,s) is 1 */
                  posrs = r * (r + 1)/2 + s;

                  /* skip positions covered by at least two variables */
                  if ( consdata->matrixval[posrs] == SCIP_INVALID ) /*lint !e777*/
                     continue;

                  var = consdata->matrixvar[posrs];
                  if ( var != NULL )
                  {
                     /* skip unfixed variable */
                     if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
                        continue;

                     val = SCIPvarGetLbLocal(var); /* fixed value */
                  }
                  else
                     val = 0.0;

                  /* check whether off-diagonal (r,s) is 1 */
                  if ( ! SCIPisFeasEQ(scip, consdata->matrixval[posrs] * val - consdata->matrixconst[posrs], 1.0) )
                     continue;

                  /* at this place rows/columns r and s are equal */

                  /* try to fix variables */
                  for (t = 1; t < blocksize; ++t)
                  {
                     SCIP_VAR* var1;
                     SCIP_VAR* var2;
                     int pos1;
                     int pos2;

                     if ( t > s )
                        pos1 = t * (t + 1)/2 + s;
                     else
                        pos1 = s * (s + 1)/2 + t;
                     var1 = consdata->matrixvar[pos1];
                     if ( var1 == NULL )
                        continue;

                     if ( t > r )
                        pos2 = t * (t + 1)/2 + r;
                     else
                        pos2 = r * (r + 1)/2 + t;
                     var2 = consdata->matrixvar[pos2];
                     if ( var2 == NULL )
                        continue;

                     /* if var1 is fixed */
                     if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetUbLocal(var1)) )
                     {
                        /* if var2 is also fixed */
                        if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
                        {
                           /* if the variables are fixed to different values, we are infeasible */
                           if ( ! SCIPisEQ(scip, SCIPvarGetLbLocal(var1), SCIPvarGetLbLocal(var2)) )
                           {
                              SCIPdebugMsg(scip, "Detected infeasibility for (%d, %d, %d) <%s>, <%s>.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              *infeasible = TRUE;
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], diagr, diags, posrs, pos1, pos2) );
                              return SCIP_OKAY;
                           }
                        }
                        else
                        {
                           SCIP_Bool tightened;

                           /* fix var2 to the same value of var1 */
                           /* currently reverse propagation does not work for this case: use INT_MAX as inferinfo */
                           SCIP_CALL( SCIPinferVarFixCons(scip, var2, SCIPvarGetLbLocal(var1), conss[c], INT_MAX, FALSE, infeasible, &tightened) );
                           if ( *infeasible )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) <%s>, <%s> detected infeasibility.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], diagr, diags, posrs, pos1, -1) );
                              return SCIP_OKAY;
                           }
                           if ( tightened )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) successfully tightened a bound of <%s> to %f.\n",
                                 r, s, t, SCIPvarGetName(var2), SCIPvarGetLbLocal(var1));
                              ++(*nprop);
                           }
                        }
                     }
                     else
                     {
                        /* if var2 is fixed (var1 is not fixed) */
                        if ( SCIPisEQ(scip, SCIPvarGetLbLocal(var2), SCIPvarGetUbLocal(var2)) )
                        {
                           SCIP_Bool tightened;

                           /* fix var1 to the same value of var2 */
                           /* currently reverse propagation does not work for this case: use INT_MAX as inferinfo */
                           SCIP_CALL( SCIPinferVarFixCons(scip, var1, SCIPvarGetLbLocal(var2), conss[c], INT_MAX, FALSE, infeasible, &tightened) );
                           if ( *infeasible )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) <%s>, <%s> detected infeasibility.\n", r, s, t,
                                 SCIPvarGetName(var1), SCIPvarGetName(var2));
                              SCIP_CALL( analyzeConflict3Minor(scip, conss[c], diagr, diags, posrs, -1, pos2) );
                              return SCIP_OKAY;
                           }
                           if ( tightened )
                           {
                              SCIPdebugMsg(scip, "Propagation on minor (%d, %d, %d) successfully tightened a bound of <%s> to %f.\n",
                                 r, s, t, SCIPvarGetName(var1), SCIPvarGetLbLocal(var2));
                              ++(*nprop);
                           }
                        }
                     }
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )

/** upgrade quadratic constraints to an SDP constraint with rank 1 */
static
SCIP_DECL_NONLINCONSUPGD(consQuadConsUpgdSdp)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* lincons;
   SCIP_VAR** linconsvars;
   SCIP_Real* linconsvals;
   SCIP_Real* linvalsterms;
   int nlinvarterms;
   int nquadvarterms;
   int nbilinterms;
   int nlinconsterms;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nupgdconss != NULL );
   assert( upgdconss != NULL );

   *nupgdconss = 0;

   /* do not upgrade modifiable/sticking at node constraints */
   if ( SCIPconsIsModifiable(cons) || SCIPconsIsStickingAtNode(cons) )
      return SCIP_OKAY;

   /* do not run in sub-SCIPs to avoid recursive reformulations due to rank 1 constraints */
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* do not upgrade after a restart */
   if ( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* make sure there is enough space to store the replacing constraints */
   if ( upgdconsssize < 1 )
   {
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   conshdlr = SCIPfindConshdlr(scip, CONSHDLRRANK1_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("rank 1 SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check whether upgrading should be performed */
   if ( ! conshdlrdata->sdpconshdlrdata->upgradequadconss )
      return SCIP_OKAY;

   /* we have to collect all variables appearing in quadratic constraints first */
   if ( conshdlrdata->sdpconshdlrdata->quadconsvars == NULL )
   {
      SCIP_CONSHDLR* nonlinearconshdlr;
      SCIP_CONS** conss;
      int nconss;
      int nvars;
      int c;
      int i;
      int nsdpvars = 0;

      int** cols;
      int** rows;
      SCIP_Real** vals;
      SCIP_VAR** vars;
      int* nvarnonz;
      int nnonz;
      int nvarscnt;
      int constcol = 0;
      int constrow = 0;
      SCIP_Real constval = -1.0;
      int nquadconss = 0;

      /* todo: The arrays quadconsidx and quadconsvars are needed to check if variables have already been seen in a
         quadratic constraint. This could be replaced with a hashmap. */
      nvars = SCIPgetNTotalVars(scip);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars) );
      conshdlrdata->sdpconshdlrdata->nquadconsidx = nvars;
      for (j = 0; j < nvars; ++j)
         conshdlrdata->sdpconshdlrdata->quadconsidx[j] = -1;

      nonlinearconshdlr = SCIPfindConshdlr(scip, "nonlinear");
      if ( nonlinearconshdlr == NULL )
      {
         SCIPerrorMessage("Nonlinear constraint handler not found\n");
         return SCIP_PLUGINNOTFOUND;
      }
      assert( nonlinearconshdlr != NULL );

      conss = SCIPconshdlrGetConss(nonlinearconshdlr);
      nconss = SCIPconshdlrGetNConss(nonlinearconshdlr);

      for (c = 0; c < nconss; ++c)
      {
         SCIP_Bool isquadratic;

         assert( conss[c] != NULL );

         SCIP_CALL( SCIPcheckQuadraticNonlinear(scip, conss[c], &isquadratic) );
         if ( ! isquadratic )
            continue;

         if ( ! SCIPexprAreQuadraticExprsVariables(SCIPgetExprNonlinear(conss[c])) )
            continue;

         ++nquadconss;

         /* Do not perform upgrade, if there are too many quadratic constraints present. */
         if ( nquadconss > conshdlrdata->sdpconshdlrdata->maxnvarsquadupgd )
         {
            SCIPdebugMsg(scip, "There are %d many quadratic constraints present in the problem, thus do not upgrade quadratic constraints to an SDPrank1 constraint.\n", nquadconss);
            SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
            SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
            return SCIP_OKAY;
         }

#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(scip, NULL, "Found quadratic constraint to upgrade:\n");
         SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIPexprGetQuadraticData(SCIPgetExprNonlinear(conss[c]), NULL, NULL, NULL, NULL, &nquadvarterms, &nbilinterms, NULL, NULL);
         assert( nquadvarterms + nbilinterms > 0 );

         for (i = 0; i < nquadvarterms; ++i)
         {
            SCIP_VAR* var;
            SCIP_EXPR* expr;
            int idx;

            /* get quadratic expression */
            SCIPexprGetQuadraticQuadTerm(SCIPgetExprNonlinear(conss[c]), i, &expr, NULL, NULL, NULL, NULL, NULL);

            var = SCIPgetVarExprVar(expr);
            assert( var != NULL );
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }
         }

         for (i = 0; i < nbilinterms; ++i)
         {
            SCIP_VAR* var;
            SCIP_EXPR* expr1;
            SCIP_EXPR* expr2;
            int idx;

            /* get bilinear expression */
            SCIPexprGetQuadraticBilinTerm(SCIPgetExprNonlinear(conss[c]), i, &expr1, &expr2, NULL, NULL, NULL);

            var = SCIPgetVarExprVar(expr1);
            assert( var != NULL );
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }

            var = SCIPgetVarExprVar(expr2);
            assert( var != NULL );
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }
         }
      }

      /* do not perform upgrade, if no sdpvars have been added */
      if ( nsdpvars == 0 )
      {
         SCIPdebugMsg(scip, "No sdp variables have been added\n");
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
         return SCIP_OKAY;
      }

      /* do not perform upgrade, if there are too many variables in the quadratic constraints, since we need sdpvars *
         sdpvars many variables for the (dual) SDPrank1 constraint */
      if ( nsdpvars > conshdlrdata->sdpconshdlrdata->maxnvarsquadupgd )
      {
         SCIPdebugMsg(scip, "There are %d many variables present in the quadratic constraints, thus do not upgrade quadratic constraints to an SDPrank1 constraint\n", nsdpvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
         return SCIP_OKAY;
      }

      /* create bilinear variables */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X, nsdpvars) );
      conshdlrdata->sdpconshdlrdata->nsdpvars = nsdpvars;

      for (i = 0; i < nsdpvars; ++i)
      {
         SCIP_Real lb1;
         SCIP_Real ub1;
         SCIP_VAR* var1;

         var1 = conshdlrdata->sdpconshdlrdata->quadconsvars[i];
         assert( var1 != NULL );
         lb1 = SCIPvarGetLbGlobal(var1);
         ub1 = SCIPvarGetUbGlobal(var1);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X[i], nsdpvars) );

         for (j = 0; j <= i; ++j)
         {
            SCIP_VARTYPE vartype;
            SCIP_VAR* var2;
            SCIP_Real lb2;
            SCIP_Real ub2;
            SCIP_Real lb;
            SCIP_Real ub;

            var2 = conshdlrdata->sdpconshdlrdata->quadconsvars[j];
            assert( var2 != NULL );
            lb2 = SCIPvarGetLbGlobal(var2);
            ub2 = SCIPvarGetUbGlobal(var2);

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "X%d#%d", i, j);

            lb = MIN3(lb1 * lb2, lb1 * ub2, ub1 * lb2);
            lb = MIN(lb, ub1 * ub2);
            ub = MAX3(lb1 * lb2, lb1 * ub2, ub1 * lb2);
            ub = MAX(ub, ub1 * ub2);

            if ( SCIPvarIsBinary(var1) && SCIPvarIsBinary(var2) )
               vartype = SCIP_VARTYPE_BINARY;
            else if ( SCIPvarIsIntegral(var1) && SCIPvarIsIntegral(var2) )
               vartype = SCIP_VARTYPE_INTEGER;
            else
               vartype = SCIP_VARTYPE_CONTINUOUS;

            SCIP_CALL( SCIPcreateVarBasic(scip, &(conshdlrdata->sdpconshdlrdata->X[i][j]), name, lb, ub, 0.0, vartype) );
            SCIP_CALL( SCIPaddVar(scip, conshdlrdata->sdpconshdlrdata->X[i][j]) );
         }
      }

      /* fill SDP data */
      nnonz = nsdpvars + nsdpvars * (nsdpvars + 1) / 2;
      SCIP_CALL( SCIPallocBufferArray(scip, &cols, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rows, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nvarnonz, nnonz) );

      /* first the terms for the original variables */
      for (j = 0; j < nsdpvars; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &cols[j], 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &rows[j], 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals[j], 1) );
         nvarnonz[j] = 1;
         cols[j][0] = 0;
         rows[j][0] = 1 + j;
         vals[j][0] = 1.0;
         vars[j] = conshdlrdata->sdpconshdlrdata->quadconsvars[j];
      }

      /* now the terms for the bilinear terms */
      nvarscnt = nsdpvars;
      for (i = 0; i < nsdpvars; ++i)
      {
         for (j = 0; j <= i; ++j)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &cols[nvarscnt], 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &rows[nvarscnt], 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals[nvarscnt], 1) );
            nvarnonz[nvarscnt] = 1;
            cols[nvarscnt][0] = 1 + j;
            rows[nvarscnt][0] = 1 + i;
            vals[nvarscnt][0] = 1.0;
            vars[nvarscnt] = conshdlrdata->sdpconshdlrdata->X[i][j];
            ++nvarscnt;
         }
      }
      assert( nvarscnt == nsdpvars + nsdpvars * (nsdpvars + 1)/2 );

      /* create corresponding rank 1 SDP constraint */
      if ( conshdlrdata->sdpconshdlrdata->upgradekeepquad )
      {
         SCIP_CALL( SCIPcreateConsSdp(scip, &conshdlrdata->sdpconshdlrdata->sdpcons, "QuadraticSDPcons", nvarscnt, nvarscnt, 1 + nsdpvars, nvarnonz,
               cols, rows, vals, vars, 1, &constcol, &constrow, &constval, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsSdpRank1(scip, &conshdlrdata->sdpconshdlrdata->sdpcons, "QuadraticSDPrank1cons", nvarscnt, nvarscnt, 1 + nsdpvars, nvarnonz,
               cols, rows, vals, vars, 1, &constcol, &constrow, &constval, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons) );
      }

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "In upgrade of quadratic constraint the following SDPrank1 constraint has been added:\n");
      SCIP_CALL( SCIPprintCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* free local memory */
      for (j = nvarscnt - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &vals[j]);
         SCIPfreeBufferArray(scip, &rows[j]);
         SCIPfreeBufferArray(scip, &cols[j]);
      }
      SCIPfreeBufferArray(scip, &nvarnonz);
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &rows);
      SCIPfreeBufferArray(scip, &cols);
   }

   /* create linear constraint for quadratic constraint */
   if ( ! conshdlrdata->sdpconshdlrdata->upgradekeepquad )
   {
      SCIP_EXPR** linexprs;
      SCIP_EXPR* expr;
      int cnt = 0;

      /* get information */
      SCIPexprGetQuadraticData(SCIPgetExprNonlinear(cons), NULL, &nlinvarterms, &linexprs, &linvalsterms, &nquadvarterms, &nbilinterms, NULL, NULL);

      /* a quadvarterm consists of a variable x and two coefficients, one for the linear term x and one for the quadratic
         term x^2, where at least one of the two coefficients is nonzero  */
      nlinconsterms = nlinvarterms + 2 * nquadvarterms + nbilinterms;
      SCIP_CALL( SCIPallocBufferArray(scip, &linconsvars, nlinconsterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linconsvals, nlinconsterms) );

      /* fill in constraint */
      for (j = 0; j < nlinvarterms; ++j)
      {
         linconsvals[cnt] = linvalsterms[j];

         SCIPexprGetQuadraticQuadTerm(SCIPgetExprNonlinear(cons), j, &expr, NULL, NULL, NULL, NULL, NULL);
         linconsvars[cnt] = SCIPgetVarExprVar(linexprs[j]);;
         assert( linconsvars[cnt] != NULL );
         ++cnt;
      }
      assert( cnt == nlinvarterms );

      for (j = 0; j < nquadvarterms; ++j)
      {
         SCIP_Real lincoef;
         SCIP_Real sqrcoef;
         SCIP_VAR* var;
         int idx;

         /* get quadratic expression */
         SCIPexprGetQuadraticQuadTerm(SCIPgetExprNonlinear(cons), j, &expr, &lincoef, &sqrcoef, NULL, NULL, NULL);
         var = SCIPgetVarExprVar(expr);
         assert( var != NULL );

         idx = SCIPvarGetIndex(var);
         idx = conshdlrdata->sdpconshdlrdata->quadconsidx[idx];
         assert( 0 <= idx && idx < conshdlrdata->sdpconshdlrdata->nsdpvars );

         /* add coefficient for linear term corresponding to the current variable (may be zero) */
         if ( ! SCIPisZero(scip, lincoef) )
         {
            linconsvals[cnt] = lincoef;
            linconsvars[cnt] = var;
            assert( linconsvars[cnt] != NULL );
            ++cnt;
         }

         /* add coefficient for quadratic term corresponding to the current variable (may be zero) */
         if ( ! SCIPisZero(scip, sqrcoef) )
         {
            linconsvals[cnt] = sqrcoef;
            linconsvars[cnt] = conshdlrdata->sdpconshdlrdata->X[idx][idx];
            assert( linconsvars[cnt] != NULL );
            ++cnt;
         }

         SCIPdebugMsg(scip, "New variable %s corresponds to squared original variable %s\n",
            SCIPvarGetName(conshdlrdata->sdpconshdlrdata->X[idx][idx]), SCIPvarGetName(var));
      }
      assert( cnt <= nlinvarterms + 2 * nquadvarterms );

      for (j = 0; j < nbilinterms; ++j)
      {
         SCIP_EXPR* expr1;
         SCIP_EXPR* expr2;
         SCIP_Real coef;
         SCIP_VAR* var1;
         SCIP_VAR* var2;
         int idx1;
         int idx2;

         /* get bilinear expression */
         SCIPexprGetQuadraticBilinTerm(SCIPgetExprNonlinear(cons), j, &expr1, &expr2, &coef, NULL, NULL);

         var1 = SCIPgetVarExprVar(expr1);
         assert( var1 != NULL );
         idx1 = SCIPvarGetIndex(var1);
         idx1 = conshdlrdata->sdpconshdlrdata->quadconsidx[idx1];
         assert( 0 <= idx1 && idx1 < conshdlrdata->sdpconshdlrdata->nsdpvars );

         var2 = SCIPgetVarExprVar(expr2);
         assert( var2 != NULL );
         idx2 = SCIPvarGetIndex(var2);
         idx2 = conshdlrdata->sdpconshdlrdata->quadconsidx[idx2];
         assert( 0 <= idx2 && idx2 < conshdlrdata->sdpconshdlrdata->nsdpvars );

         if ( idx2 > idx1 )
            SCIPswapInts(&idx1, &idx2);

         linconsvals[cnt] = coef;
         linconsvars[cnt] = conshdlrdata->sdpconshdlrdata->X[idx1][idx2];
         assert( linconsvars[cnt] != NULL );
         ++cnt;

         SCIPdebugMsg(scip, "New variable %s corresponds to product of original variables %s and %s\n",
            SCIPvarGetName(conshdlrdata->sdpconshdlrdata->X[idx1][idx2]), SCIPvarGetName(var1), SCIPvarGetName(var2));
      }
      assert( cnt <= nlinvarterms + 2 * nquadvarterms + nbilinterms );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin_%s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, cnt, linconsvars, linconsvals, SCIPgetLhsNonlinear(cons), SCIPgetRhsNonlinear(cons),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            FALSE, SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), FALSE) );

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "In upgrade of quadratic constraint the following linear constraint has been added:\n");
      SCIP_CALL( SCIPprintCons(scip, lincons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* fill in upgdconss - do not mention SDP constraint, since this has been added already */
      upgdconss[0] = lincons;
      *nupgdconss = 1;

      SCIPfreeBufferArray(scip, &linconsvals);
      SCIPfreeBufferArray(scip, &linconsvars);
   }
   else
   {
      /* todo: Check whether adding the linear constraints helps */
      *nupgdconss = 0;          /* the original quadratic constraint should be kept in the problem */
   }

   /* turn off upgrading in order to avoid a possibly infinite loop */
   conshdlrdata->sdpconshdlrdata->upgradequadconss = FALSE;

   return SCIP_OKAY;
}


/* -----------------------------------------------------------------*/
#else
/* -----------------------------------------------------------------*/


/** upgrade quadratic constraints to an SDP constraint with rank 1 */
static
SCIP_DECL_QUADCONSUPGD(consQuadConsUpgdSdp)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* lincons;
   SCIP_VAR** linconsvars;
   SCIP_Real* linconsvals;
   SCIP_VAR** linvarsterms;
   SCIP_Real* linvalsterms;
   SCIP_QUADVARTERM* quadvarterms;
   SCIP_BILINTERM* bilinterms;
   int nlinvarterms;
   int nquadvarterms;
   int nbilinterms;
   int nlinconsterms;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nupgdconss != NULL );
   assert( upgdconss != NULL );

   *nupgdconss = 0;

   /* do not upgrade modifiable/sticking at node constraints */
   if ( SCIPconsIsModifiable(cons) || SCIPconsIsStickingAtNode(cons) )
      return SCIP_OKAY;

   /* do not run in sub-SCIPs to avoid recursive reformulations due to rank 1 constraints */
   if ( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* do not upgrade after a restart */
   if ( SCIPgetNRuns(scip) > 1 )
      return SCIP_OKAY;

   /* make sure there is enough space to store the replacing constraints */
   if ( upgdconsssize < 1 )
   {
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   conshdlr = SCIPfindConshdlr(scip, CONSHDLRRANK1_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("rank 1 SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* check whether upgrading should be performed */
   if ( ! conshdlrdata->sdpconshdlrdata->upgradequadconss )
      return SCIP_OKAY;

   /* we have to collect all variables appearing in quadratic constraints first */
   if ( conshdlrdata->sdpconshdlrdata->quadconsvars == NULL )
   {
      SCIP_CONSHDLR* quadconshdlr;
      SCIP_CONS** conss;
      int nconss;
      int nvars;
      int c;
      int i;
      int nsdpvars = 0;

      int** cols;
      int** rows;
      SCIP_Real** vals;
      SCIP_VAR** vars;
      int* nvarnonz;
      int nnonz;
      int nvarscnt;
      int constcol = 0;
      int constrow = 0;
      SCIP_Real constval = -1.0;

      /* todo: The arrays quadconsidx and quadconsvars are needed to check if variables have already been seen in a
         quadratic constraint. This could be replaced with a hashmap. */
      nvars = SCIPgetNTotalVars(scip);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars) );
      conshdlrdata->sdpconshdlrdata->nquadconsidx = nvars;
      for (j = 0; j < nvars; ++j)
         conshdlrdata->sdpconshdlrdata->quadconsidx[j] = -1;

      quadconshdlr = SCIPfindConshdlr(scip, "quadratic");
      if ( quadconshdlr == NULL )
      {
         SCIPerrorMessage("Quadratic constraint handler not found\n");
         return SCIP_PLUGINNOTFOUND;
      }
      assert( quadconshdlr != NULL );

      conss = SCIPconshdlrGetConss(quadconshdlr);
      nconss = SCIPconshdlrGetNConss(quadconshdlr);

      /* Do not perform upgrade, if there are too many quadratic constraints present. */
      if ( nconss > conshdlrdata->sdpconshdlrdata->maxnvarsquadupgd )
      {
         SCIPdebugMsg(scip, "There are %d many quadratic constraints present in the problem, thus do not upgrade quadratic constraints to an SDPrank1 constraint\n", nconss);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
         return SCIP_OKAY;
      }

      for (c = 0; c < nconss; ++c)
      {
         assert( conss[c] != NULL );
#ifdef SCIP_MORE_DEBUG
         SCIPinfoMessage(scip, NULL, "Found quadratic constraint to upgrade:\n");
         SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
#endif
         nquadvarterms = SCIPgetNQuadVarTermsQuadratic(scip, conss[c]);
         quadvarterms = SCIPgetQuadVarTermsQuadratic(scip, conss[c]);

         for (i = 0; i < nquadvarterms; ++i)
         {
            SCIP_VAR* var;
            int idx;

            assert( quadvarterms != NULL );
            var = quadvarterms[i].var;
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }
         }

         nbilinterms = SCIPgetNBilinTermsQuadratic(scip, conss[c]);
         bilinterms =  SCIPgetBilinTermsQuadratic(scip, conss[c]);

         for (i = 0; i < nbilinterms; ++i)
         {
            SCIP_VAR* var;
            int idx;

            assert( bilinterms != NULL );
            var = bilinterms[i].var1;
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }

            var = bilinterms[i].var2;
            idx = SCIPvarGetIndex(var);
            assert( 0 <= idx && idx < nvars );
            if ( conshdlrdata->sdpconshdlrdata->quadconsidx[idx] < 0 )
            {
               conshdlrdata->sdpconshdlrdata->quadconsvars[nsdpvars] = var;
               conshdlrdata->sdpconshdlrdata->quadconsidx[idx] = nsdpvars++;
            }
         }
      }

      /* do not perform upgrade, if no sdpvars have been added */
      if ( nsdpvars == 0 )
      {
         SCIPdebugMsg(scip, "No sdp variables have been added\n");
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
         return SCIP_OKAY;
      }

      /* do not perform upgrade, if there are too many variables in the quadratic constraints, since we need sdpvars *
         sdpvars many variables for the (dual) SDPrank1 constraint */
      if ( nsdpvars > conshdlrdata->sdpconshdlrdata->maxnvarsquadupgd )
      {
         SCIPdebugMsg(scip, "There are %d many variables present in the quadratic constraints, thus do not upgrade quadratic constraints to an SDPrank1 constraint\n", nsdpvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, nvars);
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, nvars);
         return SCIP_OKAY;
      }

      /* create bilinear variables */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X, nsdpvars) );
      conshdlrdata->sdpconshdlrdata->nsdpvars = nsdpvars;

      for (i = 0; i < nsdpvars; ++i)
      {
         SCIP_Real lb1;
         SCIP_Real ub1;
         SCIP_VAR* var1;

         var1 = conshdlrdata->sdpconshdlrdata->quadconsvars[i];
         assert( var1 != NULL );
         lb1 = SCIPvarGetLbGlobal(var1);
         ub1 = SCIPvarGetUbGlobal(var1);

         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X[i], nsdpvars) );

         for (j = 0; j <= i; ++j)
         {
            SCIP_VARTYPE vartype;
            SCIP_VAR* var2;
            SCIP_Real lb2;
            SCIP_Real ub2;
            SCIP_Real lb;
            SCIP_Real ub;

            var2 = conshdlrdata->sdpconshdlrdata->quadconsvars[j];
            assert( var2 != NULL );
            lb2 = SCIPvarGetLbGlobal(var2);
            ub2 = SCIPvarGetUbGlobal(var2);

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "X%d#%d", i, j);

            lb = MIN3(lb1 * lb2, lb1 * ub2, ub1 * lb2);
            lb = MIN(lb, ub1 * ub2);
            ub = MAX3(lb1 * lb2, lb1 * ub2, ub1 * lb2);
            ub = MAX(ub, ub1 * ub2);

            if ( SCIPvarIsBinary(var1) && SCIPvarIsBinary(var2) )
               vartype = SCIP_VARTYPE_BINARY;
            else if ( SCIPvarIsIntegral(var1) && SCIPvarIsIntegral(var2) )
               vartype = SCIP_VARTYPE_INTEGER;
            else
               vartype = SCIP_VARTYPE_CONTINUOUS;

            SCIP_CALL( SCIPcreateVarBasic(scip, &(conshdlrdata->sdpconshdlrdata->X[i][j]), name, lb, ub, 0.0, vartype) );
            SCIP_CALL( SCIPaddVar(scip, conshdlrdata->sdpconshdlrdata->X[i][j]) );
         }
      }

      /* fill SDP data */
      nnonz = nsdpvars + nsdpvars * (nsdpvars + 1) / 2;
      SCIP_CALL( SCIPallocBufferArray(scip, &cols, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rows, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nnonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &nvarnonz, nnonz) );

      /* first the terms for the original variables */
      for (j = 0; j < nsdpvars; ++j)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &cols[j], 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &rows[j], 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &vals[j], 1) );
         nvarnonz[j] = 1;
         cols[j][0] = 0;
         rows[j][0] = 1 + j;
         vals[j][0] = 1.0;
         vars[j] = conshdlrdata->sdpconshdlrdata->quadconsvars[j];
      }

      /* now the terms for the bilinear terms */
      nvarscnt = nsdpvars;
      for (i = 0; i < nsdpvars; ++i)
      {
         for (j = 0; j <= i; ++j)
         {
            SCIP_CALL( SCIPallocBufferArray(scip, &cols[nvarscnt], 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &rows[nvarscnt], 1) );
            SCIP_CALL( SCIPallocBufferArray(scip, &vals[nvarscnt], 1) );
            nvarnonz[nvarscnt] = 1;
            cols[nvarscnt][0] = 1 + j;
            rows[nvarscnt][0] = 1 + i;
            vals[nvarscnt][0] = 1.0;
            vars[nvarscnt] = conshdlrdata->sdpconshdlrdata->X[i][j];
            ++nvarscnt;
         }
      }
      assert( nvarscnt == nsdpvars + nsdpvars * (nsdpvars + 1)/2 );

      /* create corresponding rank 1 SDP constraint */
      if ( conshdlrdata->sdpconshdlrdata->upgradekeepquad )
      {
         SCIP_CALL( SCIPcreateConsSdp(scip, &conshdlrdata->sdpconshdlrdata->sdpcons, "QuadraticSDPcons", nvarscnt, nvarscnt, 1 + nsdpvars, nvarnonz,
               cols, rows, vals, vars, 1, &constcol, &constrow, &constval, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsSdpRank1(scip, &conshdlrdata->sdpconshdlrdata->sdpcons, "QuadraticSDPrank1cons", nvarscnt, nvarscnt, 1 + nsdpvars, nvarnonz,
               cols, rows, vals, vars, 1, &constcol, &constrow, &constval, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons) );
      }

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "In upgrade of quadratic constraint the following SDPrank1 constraint has been added:\n");
      SCIP_CALL( SCIPprintCons(scip, conshdlrdata->sdpconshdlrdata->sdpcons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* free local memory */
      for (j = nvarscnt - 1; j >= 0; --j)
      {
         SCIPfreeBufferArray(scip, &vals[j]);
         SCIPfreeBufferArray(scip, &rows[j]);
         SCIPfreeBufferArray(scip, &cols[j]);
      }
      SCIPfreeBufferArray(scip, &nvarnonz);
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals);
      SCIPfreeBufferArray(scip, &rows);
      SCIPfreeBufferArray(scip, &cols);
   }

   if ( ! conshdlrdata->sdpconshdlrdata->upgradekeepquad )
   {
      int cnt = 0;

      /* create linear constraint for quadratic constraint */
      nlinvarterms = SCIPgetNLinearVarsQuadratic(scip, cons);
      linvarsterms = SCIPgetLinearVarsQuadratic(scip, cons);
      linvalsterms = SCIPgetCoefsLinearVarsQuadratic(scip, cons);
      nquadvarterms = SCIPgetNQuadVarTermsQuadratic(scip, cons);
      quadvarterms = SCIPgetQuadVarTermsQuadratic(scip, cons);
      nbilinterms = SCIPgetNBilinTermsQuadratic(scip, cons);
      bilinterms =  SCIPgetBilinTermsQuadratic(scip, cons);

      /* a quadvarterm consists of a variable x and two coefficients, one for the linear term x and one for the quadratic
         term x^2, where at least one of the two coefficients is nonzero  */
      nlinconsterms = nlinvarterms + 2 * nquadvarterms + nbilinterms;
      SCIP_CALL( SCIPallocBufferArray(scip, &linconsvars, nlinconsterms) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linconsvals, nlinconsterms) );

      /* fill in constraint */
      for (j = 0; j < nlinvarterms; ++j)
      {
         linconsvals[cnt] = linvalsterms[j];
         linconsvars[cnt] = linvarsterms[j];
         assert( linconsvars[cnt] != NULL );
         ++cnt;
      }
      assert( cnt == nlinvarterms );
      for (j = 0; j < nquadvarterms; ++j)
      {
         int idx;

         idx = SCIPvarGetIndex(quadvarterms[j].var);
         idx = conshdlrdata->sdpconshdlrdata->quadconsidx[idx];
         assert( 0 <= idx && idx < conshdlrdata->sdpconshdlrdata->nsdpvars );

         /* add coefficient for linear term corresponding to the current variable (may be zero) */
         if ( ! SCIPisZero(scip, quadvarterms[j].lincoef) )
         {
            linconsvals[cnt] = quadvarterms[j].lincoef;
            linconsvars[cnt] = quadvarterms[j].var;
            assert( linconsvars[cnt] != NULL );
            ++cnt;
         }

         /* add coefficient for quadratic term corresponding to the current variable (may be zero) */
         if ( ! SCIPisZero(scip, quadvarterms[j].sqrcoef) )
         {
            linconsvals[cnt] = quadvarterms[j].sqrcoef;
            linconsvars[cnt] = conshdlrdata->sdpconshdlrdata->X[idx][idx];
            assert( linconsvars[cnt] != NULL );
            ++cnt;
         }

         SCIPdebugMsg(scip, "New variable %s corresponds to squared original variable %s\n",
            SCIPvarGetName(conshdlrdata->sdpconshdlrdata->X[idx][idx]), SCIPvarGetName(quadvarterms[j].var));
      }
      assert( cnt <= nlinvarterms + 2 * nquadvarterms );

      for (j = 0; j < nbilinterms; ++j)
      {
         int idx1;
         int idx2;

         idx1 = SCIPvarGetIndex(bilinterms[j].var1);
         idx1 = conshdlrdata->sdpconshdlrdata->quadconsidx[idx1];
         assert( 0 <= idx1 && idx1 < conshdlrdata->sdpconshdlrdata->nsdpvars );

         idx2 = SCIPvarGetIndex(bilinterms[j].var2);
         idx2 = conshdlrdata->sdpconshdlrdata->quadconsidx[idx2];
         assert( 0 <= idx2 && idx2 < conshdlrdata->sdpconshdlrdata->nsdpvars );

         if ( idx2 > idx1 )
            SCIPswapInts(&idx1, &idx2);

         linconsvals[cnt] = bilinterms[j].coef;
         linconsvars[cnt] = conshdlrdata->sdpconshdlrdata->X[idx1][idx2];
         assert( linconsvars[cnt] != NULL );
         ++cnt;

         SCIPdebugMsg(scip, "New variable %s corresponds to product of original variables %s and %s\n",
            SCIPvarGetName(conshdlrdata->sdpconshdlrdata->X[idx1][idx2]), SCIPvarGetName(bilinterms[j].var1), SCIPvarGetName(bilinterms[j].var2));
      }
      assert( cnt <= nlinvarterms + 2 * nquadvarterms + nbilinterms );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lin_%s", SCIPconsGetName(cons));
      SCIP_CALL( SCIPcreateConsLinear(scip, &lincons, name, cnt, linconsvars, linconsvals, SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons),
            SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons), SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
            FALSE, SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons), FALSE) );

#ifdef SCIP_MORE_DEBUG
      SCIPinfoMessage(scip, NULL, "In upgrade of quadratic constraint the following linear constraint has been added:\n");
      SCIP_CALL( SCIPprintCons(scip, lincons, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* fill in upgdconss - do not mention SDP constraint, since this has been added already */
      upgdconss[0] = lincons;
      *nupgdconss = 1;

      SCIPfreeBufferArray(scip, &linconsvals);
      SCIPfreeBufferArray(scip, &linconsvars);
   }
   else
   {
      /* todo: Check whether adding the linear constraints helps */
      *nupgdconss = 0;          /* the original quadratic constraint should be kept in the problem */
   }

   /* turn off upgrading in order to avoid a possibly infinite loop */
   conshdlrdata->sdpconshdlrdata->upgradequadconss = FALSE;

   return SCIP_OKAY;
}

#endif


#if SCIP_VERSION >= 900

/* define indices for operator nodes in symmetry graph */
#define OP_SDP_DIMENSION 1

/** adds symmetry information of constraint to a symmetry detection graph */
static
SCIP_RETCODE addSymmetryInformation(
   SCIP*                 scip,               /**< SCIP pointer */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SYM_SYMTYPE           symtype,            /**< type of symmetries that need to be added */
   SCIP_CONS*            cons,               /**< SDP constraint */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   SCIP_Bool*            success             /**< pointer to store whether symmetry information could be added */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int* dimnodeidx = NULL;
   int consnodeidx;
   int nunique = 0;
   int v;
   int i;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( graph != NULL );
   assert( success != NULL );

   *success = TRUE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->issymunique != NULL );

   /* add node initializing constraint (with artificial rhs) */
   SCIP_CALL( SCIPaddSymgraphConsnode(scip, graph, cons, 0.0, 0.0, &consnodeidx) );

   /* check whether all variable matrices are unqiue */
   for (v = 0; v < consdata->nvars; ++v)
   {
      if ( consdata->issymunique[v] )
         ++nunique;
   }

   /* if there are some variable matrices that are not unique */
   if ( nunique < consdata->nvars )
   {
      /* for each constraint, add nodes for the dimensions of the matrix */
      SCIP_CALL( SCIPallocBufferArray(scip, &dimnodeidx, consdata->blocksize) );
      for (i = 0; i < consdata->blocksize; ++i)
      {
         SCIP_CALL( SCIPaddSymgraphOpnode(scip, graph, OP_SDP_DIMENSION, &dimnodeidx[i]) );

         /* connect new node to constraint node */
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, consnodeidx, dimnodeidx[i], FALSE, 0.0) );
      }
   }
   SCIPdebugMsg(scip, "Constraint <%s> has %d unique variable matrices among %d.\n", SCIPconsGetName(cons), nunique, consdata->nvars);

   /* for each variable matrix entry, add two nodes corresponding to row/column index and connect it with variable and
    * dimension nodes */
   for (v = 0; v < consdata->nvars; ++v)
   {
      int varnodeidx;

      varnodeidx = SCIPgetSymgraphVarnodeidx(scip, graph, consdata->vars[v]);

      /* connect variables with unique matrices to constraint node with unique color to signify that they should be fixed */
      if ( consdata->issymunique[v] )
      {
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, varnodeidx, consnodeidx, TRUE, conshdlrdata->symuniqueid++) );
         continue;
      }
      assert( dimnodeidx != NULL );

      /* loop over all entries in the matrix */
      for (i = 0; i < consdata->nvarnonz[v]; ++i)
      {
         int nodeidx;
         int node;

         /* add node for each symmetry row/col pair - we use values since this seems to improve the speed */
         SCIP_CALL( SCIPaddSymgraphValnode(scip, graph, consdata->val[v][i], &nodeidx) );

         /* add edge to variable node */
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, nodeidx, varnodeidx, FALSE, 0.0) );

         /* add edges to dimension nodes */
         assert( 0 <= consdata->row[v][i] && consdata->row[v][i] < consdata->blocksize );
         node = dimnodeidx[consdata->row[v][i]];
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, nodeidx, node, TRUE, consdata->val[v][i]) );

         assert( 0 <= consdata->col[v][i] && consdata->col[v][i] < consdata->blocksize );
         node = dimnodeidx[consdata->col[v][i]];
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, nodeidx, node, TRUE, consdata->val[v][i]) );
      }
   }

   /* if not all variable matrices are unique, we treat the constant matrix */
   if ( nunique < consdata->nvars )
   {
      assert( dimnodeidx != NULL );

      /* for each constant matrix entry, add two nodes corresponding to row/column index and connect it with dimension
       * nodes */
      for (i = 0; i < consdata->constnnonz; ++i)
      {
         int nodeidx;
         int node;

         /* add node for each symmetry row/col pair - we use values since this seems to improve the speed */
         SCIP_CALL( SCIPaddSymgraphValnode(scip, graph, consdata->constval[i], &nodeidx) );

         /* add edges to dimension nodes */
         assert( 0 <= consdata->constrow[i] && consdata->constrow[i] < consdata->blocksize );
         node = dimnodeidx[consdata->constrow[i]];
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, nodeidx, node, TRUE, consdata->constval[i]) );

         assert( 0 <= consdata->constcol[i] && consdata->constcol[i] < consdata->blocksize );
         node = dimnodeidx[consdata->constcol[i]];
         SCIP_CALL( SCIPaddSymgraphEdge(scip, graph, nodeidx, node, TRUE, consdata->constval[i]) );
      }
   }

   SCIPfreeBufferArrayNull(scip, &dimnodeidx);

   return SCIP_OKAY;
}
#endif


/*
 * callbacks
 */

/** informs constraint handler that the presolving process is being started */
static
SCIP_DECL_CONSINITPRE(consInitpreSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* reset numbers; called for both the rank1 and ordinary constraint handler */
   conshdlrdata->neigveccuts = 0; /* this is used to give the eigenvector-cuts distinguishable names */
   conshdlrdata->ncmir = 0;
   conshdlrdata->ndiaggezerocuts = 0; /* this is used to give the diagGEzero-cuts distinguishable names */
   conshdlrdata->n1x1blocks = 0; /* this is used to give the lp constraints resulting from 1x1 sdp-blocks distinguishable names */
   conshdlrdata->symuniqueid = 0;
   conshdlrdata->ncallspropub = 0;
   conshdlrdata->ncallsproptb = 0;
   conshdlrdata->ncallsprop3minor = 0;
   conshdlrdata->maxtimepropub = 0.0;
   conshdlrdata->maxtimeproptb = 0.0;
   conshdlrdata->maxtimeprop3minor = 0.0;
   conshdlrdata->npropub = 0;
   conshdlrdata->nproptb = 0;
   conshdlrdata->nprop3minor = 0;
   conshdlrdata->npropcutoffub = 0;
   conshdlrdata->npropcutofftb = 0;
   conshdlrdata->npropcutoff3m = 0;
   conshdlrdata->npropintrndub = 0;
   conshdlrdata->npropintrndtb = 0;
   conshdlrdata->nproppreub = 0;
   conshdlrdata->nproppretb = 0;
   conshdlrdata->nproppre3m = 0;
   conshdlrdata->nproppreintrndub = 0;
   conshdlrdata->nproppreintrndtb = 0;
   conshdlrdata->nproppreintrnd3m = 0;
   conshdlrdata->npropprobub = 0;
   conshdlrdata->npropprobtb = 0;
   conshdlrdata->npropprob3minor = 0;

   /* create clocks */
   if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
   {
      if ( conshdlrdata->sdpconshdlrdata->propubtime == NULL )
         SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->propubtime) );
      if ( conshdlrdata->sdpconshdlrdata->proptbtime == NULL )
         SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->proptbtime) );
      if ( conshdlrdata->sdpconshdlrdata->prop3minortime == NULL )
         SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->prop3minortime) );
   }

   /* determine unique matrices */
   SCIP_CALL( checkSymUniqueMatrices(scip, nconss, conss) );

   return SCIP_OKAY;
}

/** locks a variable up if the corresponding constraint matrix is not positive semidefinite, locks it down if it is not negative semidefinite */
static
SCIP_DECL_CONSLOCK(consLockSdp)
{/*lint --e{715}*/
   SCIP_Real* Aj;
   SCIP_CONSDATA* consdata;
   int nvars;
   int v;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   nvars = consdata->nvars;

   SCIPdebugMsg(scip, "locking method of <%s>.\n", SCIPconsGetName(cons));

   /* rank-1 constraints are always up- and down-locked */
   if ( consdata->rankone )
   {
      if ( consdata->locks == NULL )
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->locks, nvars) );

      for (v = 0; v < consdata->nvars; ++v)
      {
         consdata->locks[v] = 0;
         SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
      }
      return SCIP_OKAY;
   }

   /* if locks have not yet been computed */
   if ( consdata->locks == NULL )
   {
      SCIP_Real mineigenvalue = SCIP_REAL_MAX;
      SCIP_Real eigenvalue;
      int blocksize;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->locks, nvars) );

      blocksize = consdata->blocksize;

      SCIP_CALL( SCIPallocBufferArray(scip, &Aj, blocksize * blocksize) ); /*lint !e647*/

      for (v = 0; v < nvars; v++)
      {
         SCIP_CALL( SCIPconsSdpGetFullAj(scip, cons, v, Aj) );
         consdata->locks[v] = -2;  /* unintitialized */

         /* compute the smallest eigenvalue */
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, 1, &eigenvalue, NULL) );
         if ( SCIPisNegative(scip, eigenvalue) )
         {
            /* as the lowest eigenvalue is negative, the matrix is not positive semidefinite, so adding more of it can remove positive
             * semidefiniteness of the SDP-matrix */
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlocksneg, nlockspos) );
            consdata->locks[v] = 1; /* up-lock */
         }
         if ( eigenvalue < mineigenvalue )
            mineigenvalue = eigenvalue;

         /* if the smallest eigenvalue is already positive, we don't need to compute the biggest one */
         if ( SCIPisPositive(scip, eigenvalue) )
         {
            /* as an eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
             * semidefiniteness of the SDP-matrix */
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlockspos, nlocksneg) );
            consdata->locks[v] = -1; /* down-lock */
         }
         else
         {
            /* compute the biggest eigenvalue */
            SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, Aj, blocksize, &eigenvalue, NULL) );
            if ( SCIPisPositive(scip, eigenvalue) )
            {
               /* as the biggest eigenvalue is positive, the matrix is not negative semidefinite, so substracting more of it can remove positive
                * semidefiniteness of the SDP-matrix */
               SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlockspos, nlocksneg) );
               if ( consdata->locks[v] == 1 )
               {
                  consdata->locks[v] = 0;  /* up- and down-lock */
               }
               else
                  consdata->locks[v] = -1; /* down-lock */
            }
         }
      }

      if ( SCIPisFeasGE(scip, mineigenvalue, 0.0) )
      {
         consdata->allmatricespsd = TRUE;
         if ( SCIPgetSubscipDepth(scip) == 0 )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "All matrices are positive semidefinite (minimial eigenvalue: %g).\n", mineigenvalue);
      }
      consdata->initallmatricespsd = TRUE;

      SCIPfreeBufferArray(scip, &Aj);
   }
   else
   {
#ifndef NDEBUG
      SCIP_CALL( checkVarsLocks(scip, cons) );
#endif
      for (v = 0; v < nvars; v++)
      {
         if ( consdata->locks[v] == 1 )  /* up-lock */
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlocksneg, nlockspos) );
         }
         else if ( consdata->locks[v] == -1 )  /* down-lock */
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlockspos, nlocksneg) );
         }
         else if ( consdata->locks[v] == 0 )  /* up and down lock */
         {
            SCIP_CALL( SCIPaddVarLocksType(scip, consdata->vars[v], locktype, nlockspos + nlocksneg, nlockspos + nlocksneg) );
         }
         else
            assert( consdata->locks[v] == -2 );
      }
   }

   return SCIP_OKAY;
}

/** deinitialization method of constraint handler (called before transformed problem is freed)
 *
 * At the end of the solution process, the parameter for adding linear constraints in presolving needs to be reset.
 */
static
SCIP_DECL_CONSEXIT(consExitSdp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* possibly output more statistics */
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 )
   {
      if ( conshdlrdata->sdpconshdlrdata->additionalstats )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of calls of propagateUpperBounds in propagation: %d\n", conshdlrdata->ncallspropub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of calls of tightenBounds in propagation: %d\n", conshdlrdata->ncallsproptb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of calls of propagate3Minors in propagation: %d\n", conshdlrdata->ncallsprop3minor);
         if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Overall time spent for propagateUpperBounds in propagation: %f\n", SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->propubtime) );
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Overall time spent for tightenBounds in propagation: %f\n", SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->proptbtime) );
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Overall time spent for propagate3Minors in propagation: %f\n", SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->prop3minortime) );
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Maximal time spent for one round of propagateUpperBounds in propagation: %f\n", conshdlrdata->maxtimepropub);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Maximal time spent for one round of tightenBounds in propagation: %f\n", conshdlrdata->maxtimeproptb);
            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Maximal time spent for one round of propagate3Minors in propagation: %f\n", conshdlrdata->maxtimeprop3minor);
         }
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of propagations through upper bounds: %d\n", conshdlrdata->npropub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of tightened bounds in propagation:   %d\n", conshdlrdata->nproptb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of propagations through 3x3 minors:   %d\n", conshdlrdata->nprop3minor);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of cutoffs in propagation through upper bounds: %d\n", conshdlrdata->npropcutoffub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of cutoffs in propagation of bound tightening:  %d\n", conshdlrdata->npropcutofftb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of cutoffs in propagation through 3x3 minors:   %d\n", conshdlrdata->npropcutoff3m);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of rounded bounds of integer variables in propagation through upper bounds: %d\n", conshdlrdata->npropintrndub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of rounded bounds of integer variables in propagation of bound tightening:  %d\n", conshdlrdata->npropintrndtb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of propagations through upper bounds in presolving: %d\n", conshdlrdata->nproppreub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of tightened bounds in propagation in presolving:   %d\n", conshdlrdata->nproppretb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of propagations through 3x3 minors in presolving:   %d\n", conshdlrdata->nproppre3m);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of rounded bounds of integer variables in propagation through upper bounds in presolving: %d\n", conshdlrdata->nproppreintrndub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of rounded bounds of integer variables in propagation of bound tightening in presolving:  %d\n", conshdlrdata->nproppreintrndtb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of rounded bounds of integer variables in propagation through 3x3 minors in presolving:   %d\n", conshdlrdata->nproppreintrnd3m);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number propagations through upper bounds in probing:  %d\n", conshdlrdata->npropprobub);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of tightened bounds in propagation in probing: %d\n", conshdlrdata->npropprobtb);
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, 0, "Number of propagation through 3x3 minors in probing: %d\n", conshdlrdata->npropprob3minor);
      }

      /* reset counters */
      conshdlrdata->ncallspropub = 0;
      conshdlrdata->ncallsproptb = 0;
      conshdlrdata->ncallsprop3minor = 0;
      conshdlrdata->maxtimepropub = 0.0;
      conshdlrdata->maxtimeproptb = 0.0;
      conshdlrdata->maxtimeprop3minor = 0.0;
      conshdlrdata->npropub = 0;
      conshdlrdata->nproptb = 0;
      conshdlrdata->nprop3minor = 0;
      conshdlrdata->npropcutoffub = 0;
      conshdlrdata->npropcutofftb = 0;
      conshdlrdata->npropcutoff3m = 0;
      conshdlrdata->npropintrndub = 0;
      conshdlrdata->npropintrndtb = 0;
      conshdlrdata->nproppreub = 0;
      conshdlrdata->nproppretb = 0;
      conshdlrdata->nproppre3m = 0;
      conshdlrdata->nproppreintrndub = 0;
      conshdlrdata->nproppreintrndtb = 0;
      conshdlrdata->nproppreintrnd3m = 0;
      conshdlrdata->npropprobub = 0;
      conshdlrdata->npropprobtb = 0;
      conshdlrdata->npropprob3minor = 0;
   }

   /* reset parameter triedlinearconss */
   conshdlrdata->sdpconshdlrdata->triedlinearconss = FALSE;
   conshdlrdata->sdpconshdlrdata->triedvarbounds = FALSE;

   /* free clocks */
   if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
   {
      if ( conshdlrdata->sdpconshdlrdata->propubtime != NULL )
         SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->sdpconshdlrdata->propubtime) );
      if ( conshdlrdata->sdpconshdlrdata->proptbtime != NULL )
         SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->sdpconshdlrdata->proptbtime) );
      if ( conshdlrdata->sdpconshdlrdata->prop3minortime != NULL )
         SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->sdpconshdlrdata->prop3minortime) );
   }

   return SCIP_OKAY;
}

/** after presolving variables are fixed and multiaggregated */
static
SCIP_DECL_CONSEXITPRE(consExitpreSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   if ( conss == NULL )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Exitpre method of conshdlr <%s>.\n", SCIPconshdlrGetName(conshdlr));

   if ( SCIPgetStatus(scip) != SCIP_STATUS_OPTIMAL && SCIPgetStatus(scip) != SCIP_STATUS_INFEASIBLE )
   {
      SCIP_CALL( fixAndAggrVars(scip, conss, nconss, TRUE) );
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->sdpconshdlrdata->quadconsidx, conshdlrdata->sdpconshdlrdata->nquadconsidx);
   SCIPfreeBlockMemoryArrayNull(scip, &conshdlrdata->sdpconshdlrdata->quadconsvars, conshdlrdata->sdpconshdlrdata->nquadconsidx);
   if ( conshdlrdata->sdpconshdlrdata->X != NULL )
   {
      SCIPdebugMsg(scip, "Releasing additional variables from upgrading method\n");
      for (i = 0; i < conshdlrdata->sdpconshdlrdata->nsdpvars; ++i)
      {
         for (j = 0; j <= i; ++j)
         {
            SCIP_CALL( SCIPreleaseVar(scip, &(conshdlrdata->sdpconshdlrdata->X[i][j])) );
         }
         SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X[i], conshdlrdata->sdpconshdlrdata->nsdpvars);
      }
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->sdpconshdlrdata->X, conshdlrdata->sdpconshdlrdata->nsdpvars);
   }

   if ( conshdlrdata->sdpconshdlrdata->sdpcons != NULL )
   {
      SCIPdebugMsg(scip, "Releasing constraint %s from upgrading method\n", SCIPconsGetName(conshdlrdata->sdpconshdlrdata->sdpcons) );
      SCIP_CALL( SCIPreleaseCons(scip, &conshdlrdata->sdpconshdlrdata->sdpcons) );
   }

   return SCIP_OKAY;
}

/** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
 *
 *  At the beginning of the solution process the stored rank one submatrix is reset.
 */
static
SCIP_DECL_CONSINITSOL(consInitsolSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != NULL );
   assert( conshdlr != NULL );

   if ( conss == NULL )
      return SCIP_OKAY;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( (conshdlrdata->sdpconshdlrdata->sparsifycut || conshdlrdata->sdpconshdlrdata->multiplesparsecuts) && conshdlrdata->randnumgen == NULL )
   {
      SCIP_CALL( SCIPcreateRandom(scip, &conshdlrdata->randnumgen, 64293, FALSE) );
   }

   conshdlrdata->relaxsdp = SCIPfindRelax(scip, "SDP");

   /* make sure that quadratic constraints are added */
   if ( SCIPgetSubscipDepth(scip) == 0 && conshdlrdata->sdpconshdlrdata->quadconsrank1 )
   {
      int naddconss = 0;
      SCIP_CALL( addRank1QuadConss(scip, conshdlr, conss, nconss, &naddconss) );
      SCIPdebugMsg(scip, "Added %d quadratic constraints for rank 1 constraints.\n", naddconss);

      /* turn off upgrading in order to avoid upgrading to a rank-1 constraint again */
      conshdlrdata->sdpconshdlrdata->upgradequadconss = FALSE;
   }

   return SCIP_OKAY;
}

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSdp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int nprop = 0;
   int nintrnd = 0;
   SCIP_Real oldtime = 0.0;
   SCIP_Real newtime = 0.0;

   assert( conshdlr != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* if we want to propagate upper bounds */
   if ( conshdlrdata->sdpconshdlrdata->propupperbounds )
   {
      *result = SCIP_DIDNOTFIND;

      SCIPdebugMsg(scip, "Propagate upper bounds of conshdlr <%s> %s...\n", SCIPconshdlrGetName(conshdlr), SCIPinProbing(scip) ? "(in probing) " : "");

      if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
      {
         oldtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->propubtime);
         SCIP_CALL( SCIPstartClock(scip, conshdlrdata->sdpconshdlrdata->propubtime) );
      }

      SCIP_CALL( propagateUpperBounds(scip, conss, nconss, &infeasible, &nprop, &nintrnd) );

      if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
      {
         SCIP_CALL( SCIPstopClock(scip, conshdlrdata->sdpconshdlrdata->propubtime) );
         newtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->propubtime);
         if ( SCIPisGT(scip, newtime - oldtime, conshdlrdata->sdpconshdlrdata->maxtimepropub) )
            conshdlrdata->sdpconshdlrdata->maxtimepropub = newtime - oldtime;
      }

      ++conshdlrdata->sdpconshdlrdata->ncallspropub;

      if ( infeasible )
      {
         ++conshdlrdata->sdpconshdlrdata->npropcutoffub;
         SCIPdebugMsg(scip, "Propagation of upper bounds detected cutoff.\n");
         *result = SCIP_CUTOFF;
      }
      else
      {
         if ( nprop > 0 )
         {
            if ( ! SCIPinProbing(scip) )
            {
               conshdlrdata->sdpconshdlrdata->npropub += nprop;
               conshdlrdata->sdpconshdlrdata->npropintrndub += nintrnd;
            }
            else
               conshdlrdata->sdpconshdlrdata->npropprobub += nprop;

            SCIPdebugMsg(scip, "Propagation of upper bounds tightened %d bounds.\n", nprop);
            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   /* if we want to propagate bound tightening */
   if ( conshdlrdata->sdpconshdlrdata->proptightenbounds )
   {
      /* possibly avoid propagation in probing */
      if ( conshdlrdata->sdpconshdlrdata->proptbprobing || ! SCIPinProbing(scip) )
      {

         if ( *result == SCIP_DIDNOTRUN )
            *result = SCIP_DIDNOTFIND;

         SCIPdebugMsg(scip, "Propagate tighten bounds of conshdlr <%s> ...\n", SCIPconshdlrGetName(conshdlr));

         nprop = 0;
         nintrnd = 0;
         oldtime = 0.0;
         newtime = 0.0;

         if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
         {
            oldtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->proptbtime);
            SCIP_CALL( SCIPstartClock(scip, conshdlrdata->sdpconshdlrdata->proptbtime) );
         }

         SCIP_CALL( tightenBounds(scip, conss, nconss, conshdlrdata->sdpconshdlrdata->tightenboundscont, &nprop, &nintrnd, &infeasible) );

         if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
         {
            SCIP_CALL( SCIPstopClock(scip, conshdlrdata->sdpconshdlrdata->proptbtime) );
            newtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->proptbtime);
            if ( SCIPisGT(scip, newtime - oldtime, conshdlrdata->sdpconshdlrdata->maxtimeproptb) )
               conshdlrdata->sdpconshdlrdata->maxtimeproptb = newtime - oldtime;
         }

         ++conshdlrdata->sdpconshdlrdata->ncallsproptb;


         if ( infeasible )
         {
            ++conshdlrdata->sdpconshdlrdata->npropcutofftb;
            SCIPdebugMsg(scip, "Propagation of bound tightening detected cutoff.\n");
            *result = SCIP_CUTOFF;
         }
         else
         {
            if ( nprop > 0 )
            {
               if ( ! SCIPinProbing(scip) )
               {
                  conshdlrdata->sdpconshdlrdata->nproptb += nprop;
                  conshdlrdata->sdpconshdlrdata->npropintrndtb += nintrnd;
               }
               else
                  conshdlrdata->sdpconshdlrdata->npropprobtb += nprop;

               SCIPdebugMsg(scip, "Propagation of bound tightening tightened %d bounds.\n", nprop);
               *result = SCIP_REDUCEDDOM;
            }
         }
      }
   }

   /* if we want to propagate 3x3 minors and we are not in probing */
   if ( conshdlrdata->sdpconshdlrdata->prop3minors )
   {
      if ( conshdlrdata->sdpconshdlrdata->prop3mprobing || ! SCIPinProbing(scip) )
      {
         if ( *result == SCIP_DIDNOTRUN )
            *result = SCIP_DIDNOTFIND;

         SCIPdebugMsg(scip, "Propagate 3x3 minors of conshdlr <%s> %s...\n", SCIPconshdlrGetName(conshdlr), SCIPinProbing(scip) ? "(in probing) " : "");

         nprop = 0;
         oldtime = 0.0;
         newtime = 0.0;

         if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
         {
            oldtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->prop3minortime);
            SCIP_CALL( SCIPstartClock(scip, conshdlrdata->sdpconshdlrdata->prop3minortime) );
         }

         SCIP_CALL( propagate3Minors(scip, conss, nconss, conshdlrdata->sdpconshdlrdata->nonconst3minors, &infeasible, &nprop) );

         if ( conshdlrdata->sdpconshdlrdata->enableproptiming )
         {
            SCIP_CALL( SCIPstopClock(scip, conshdlrdata->sdpconshdlrdata->prop3minortime) );
            newtime = SCIPgetClockTime(scip, conshdlrdata->sdpconshdlrdata->prop3minortime);
            if ( SCIPisGT(scip, newtime - oldtime, conshdlrdata->sdpconshdlrdata->maxtimeprop3minor) )
               conshdlrdata->sdpconshdlrdata->maxtimeprop3minor = newtime - oldtime;
         }

         ++conshdlrdata->sdpconshdlrdata->ncallsprop3minor;

         if ( infeasible )
         {
            ++conshdlrdata->sdpconshdlrdata->npropcutoff3m;
            SCIPdebugMsg(scip, "Propagation of 3x3 minors detected cutoff.\n");
            *result = SCIP_CUTOFF;
         }
         else
         {
            if ( nprop > 0 )
            {
               assert( ! SCIPinProbing(scip) || conshdlrdata->sdpconshdlrdata->prop3mprobing );
               if ( ! SCIPinProbing(scip) )
                  conshdlrdata->sdpconshdlrdata->nprop3minor += nprop;
               else
                  conshdlrdata->sdpconshdlrdata->npropprob3minor += nprop;

               SCIPdebugMsg(scip, "Propagation of 3x3 minors tightened %d bounds.\n", nprop);
               *result = SCIP_REDUCEDDOM;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** propagation conflict resolving method of constraint handler */
static
SCIP_DECL_CONSRESPROP(consRespropSdp)
{
   SCIP_CONSDATA* consdata;
   int diags;
   int diagt;
   int s;
   int t;

   assert( conshdlr != NULL );
   assert( cons != NULL );
   assert( infervar != NULL );
   assert( result != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 );
   assert( ! consdata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLRRANK1_NAME) == 0 );

   SCIPdebugMsg(scip, "Executing conflict resolving method of <%s> constraint handler.\n", SCIPconshdlrGetName(conshdlr));

   /* skip unresolvable cases */
   if ( inferinfo == INT_MAX )
   {
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /* if inferinfo is >=, the bound change came from propagateUpperBounds() */
   if ( inferinfo >= 0 )
   {
      s = inferinfo / consdata->blocksize;
      t = inferinfo % consdata->blocksize;
      assert( 0 <= s && s < consdata->blocksize );
      assert( 0 <= t && t < consdata->blocksize );
      assert( consdata->matrixvar[s * (s + 1)/2 + t] == infervar );

      diags = s * (s + 1)/2 + s;
      diagt = t * (t + 1)/2 + t;

      assert( consdata->matrixvar[diags] != NULL );
      assert( consdata->matrixvar[diagt] != NULL );
      assert( consdata->matrixval[diags] != SCIP_INVALID ); /*lint !e777*/
      assert( consdata->matrixval[diagt] != SCIP_INVALID ); /*lint !e777*/

      if ( consdata->matrixval[diags] > 0.0 )
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diags], bdchgidx) );
      else
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diags], bdchgidx) );

      if ( consdata->matrixval[diagt] > 0.0 )
         SCIP_CALL( SCIPaddConflictUb(scip, consdata->matrixvar[diagt], bdchgidx) );
      else
         SCIP_CALL( SCIPaddConflictLb(scip, consdata->matrixvar[diagt], bdchgidx) );

      *result = SCIP_SUCCESS;
   }
   else
   {
      int i;
      int k;

      /* otherwise the bound change came from tightenbounds() */
      i = -inferinfo - 1;
      assert( 0 <= i && i < consdata->nvars );

      /* the upper bounds of all other variables are responsible */
      for (k = 0; k < consdata->nvars; ++k)
      {
         if ( k == i )
            continue;

         SCIP_CALL( SCIPaddConflictUb(scip, consdata->vars[k], bdchgidx) );
      }
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}

/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int nprop = 0;
   int nintrnd = 0;
   int c;

   assert( conshdlr != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *result = SCIP_DIDNOTRUN;

   /* if some variables have been fixed or aggregated */
   if ( nnewfixedvars + nnewaggrvars > 0 )
   {
      /* check whether some fixed or aggregated variables can be removed from constraint */
      SCIP_CALL( fixAndAggrVars(scip, conss, nconss, TRUE) );
   }

   /* check for empty constraints */
   for (c = 0; c < nconss && *result != SCIP_CUTOFF; ++c)
   {
      SCIP_CONSDATA* consdata;

      assert( conss[c] != NULL );
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      /* for empty constraint check whether constant matrix is not infeasible */
      if ( consdata->nvars <= 0 )
      {
         SCIP_Real* constmatrix;
         SCIP_Real eigenvalue;
         int blocksize;

         blocksize = consdata->blocksize;
         SCIP_CALL( SCIPallocBufferArray(scip, &constmatrix, blocksize * blocksize) );
         SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, conss[c], constmatrix) );
         SCIP_CALL( SCIPlapackComputeIthEigenvalue(SCIPbuffer(scip), FALSE, blocksize, constmatrix, blocksize, &eigenvalue, NULL) );

         /* if largest eigenvalue is positive then minus the constant matrix is not psd and we are infeasible */
         if ( SCIPisFeasPositive(scip, eigenvalue) )
         {
            SCIPdebugMsg(scip, "Infeasible constraint <%s> containts no variable.\n", SCIPconsGetName(conss[c]));
            *result = SCIP_CUTOFF;
         }
         else
         {
            SCIPdebugMsg(scip, "Feasible constraint <%s> containts no variable, removing.\n", SCIPconsGetName(conss[c]));
            SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
            ++(ndelconss);
            *result = SCIP_SUCCESS;
         }
         SCIPfreeBufferArray(scip, &constmatrix);
      }
      else if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;
   }

   if ( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   /* call propagation */
   if ( conshdlrdata->sdpconshdlrdata->propubpresol )
   {
      SCIPdebugMsg(scip, "Propagate upper bounds of conshdlr <%s> ...\n", SCIPconshdlrGetName(conshdlr));

      SCIP_CALL( propagateUpperBounds(scip, conss, nconss, &infeasible, &nprop, &nintrnd) );

      if ( infeasible )
      {
         SCIPdebugMsg(scip, "Presolving detected cutoff.\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMsg(scip, "Presolving upper bounds: %d.\n", nprop);
         if ( nprop > 0 )
         {
            *nchgbds += nprop;
            conshdlrdata->sdpconshdlrdata->nproppreub += nprop;
            conshdlrdata->sdpconshdlrdata->nproppreintrndub += nintrnd;
            *result = SCIP_SUCCESS;
         }
      }
   }

   /* add constraints in initial round */
   if ( SCIPconshdlrGetNPresolCalls(conshdlr) == 0 )
   {
      int noldaddconss;
      int nolddelconss;
      int noldchgbds;
      int noldchgcoefs;

      if ( *result == SCIP_DIDNOTRUN )
         *result = SCIP_DIDNOTFIND;

      noldaddconss = *naddconss;
      nolddelconss = *ndelconss;
      noldchgbds = *nchgbds;
      noldchgcoefs = *nchgcoefs;

      SCIP_CALL( move_1x1_blocks_to_lp(scip, conshdlr, conss, nconss, naddconss, ndelconss, nchgbds, &infeasible) );
      if ( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
      if ( noldaddconss != *naddconss || nolddelconss != *ndelconss || noldchgbds != *nchgbds )
         *result = SCIP_SUCCESS;

      /* possibly compute tightening of matrices */
      if ( conshdlrdata->sdpconshdlrdata->tightenmatrices )
      {
         SCIP_CALL( tightenMatrices(scip, conss, nconss, nchgcoefs) );
         if ( noldchgcoefs != *nchgcoefs )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Tightened %d SDP coefficient matrices.\n", *nchgcoefs - noldchgcoefs);
            *result = SCIP_SUCCESS;
         }
      }

      /* possibly tighten bounds */
      if ( conshdlrdata->sdpconshdlrdata->tightenbounds )
      {
         nprop = 0;
         nintrnd = 0;
         SCIP_CALL( tightenBounds(scip, conss, nconss, conshdlrdata->sdpconshdlrdata->tightenboundscont, &nprop, &nintrnd, &infeasible) );
         if ( infeasible )
         {
            *result = SCIP_CUTOFF;
            return SCIP_OKAY;
         }

         if ( nprop > 0 )
         {
            *nchgbds += nprop;
            conshdlrdata->sdpconshdlrdata->nproppretb += nprop;
            conshdlrdata->sdpconshdlrdata->nproppreintrndtb += nintrnd;

            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Tightened %d bounds using SDP constraints.\n", nprop);
            *result = SCIP_SUCCESS;
         }

      }

      /* In the following, we add linear constraints. This is needed only once. We assume that this is only necessary in
       * the main SCIP instance. The diagzeroimpl-cuts are added as basic linear constraints which are initial, as well
       * as separated, enforced, checked and propagated. All other linear constraints that are added below are redundant
       * for the SDP-constraint, and thus they are neither checked nor enforced and also not initial, so that they do
       * not appear in the SDP or LP relaxation. If LPs are solved, these linear constraints are separated and
       * propagated; if SDPs are solved, they are only propagated. */
      if ( SCIPgetSubscipDepth(scip) == 0 && ! conshdlrdata->sdpconshdlrdata->triedlinearconss )
      {
         int solvesdpsparam;
         SCIP_Bool solvesdps;

         SCIP_CALL( SCIPgetIntParam(scip, "misc/solvesdps", &solvesdpsparam) );

         if ( solvesdpsparam == 1 )
            solvesdps = TRUE;
         else
            solvesdps = FALSE;

         conshdlrdata->sdpconshdlrdata->triedlinearconss = TRUE;
         if ( conshdlrdata->sdpconshdlrdata->diaggezerocuts )
         {
            noldaddconss = *naddconss;
            noldchgbds = *nchgbds;
            SCIP_CALL( diagGEzero(scip, conshdlr, conss, nconss, solvesdps, naddconss, nchgbds, &infeasible) );
            SCIPdebugMsg(scip, "Diagonal entries: added %d constraints and changed %d bounds.\n", *naddconss - noldaddconss, *nchgbds - noldchgbds);

            if ( infeasible )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            if ( noldaddconss != *naddconss || noldchgbds != *nchgbds )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added %d constraints for SDP diagonals to be nonnegative and changed %d bounds.\n", *naddconss - noldaddconss, *nchgbds - noldchgbds);
               *result = SCIP_SUCCESS;
            }
         }

         if ( *result != SCIP_CUTOFF && conshdlrdata->sdpconshdlrdata->diagzeroimplcuts )
         {
            noldaddconss = *naddconss;
            SCIP_CALL( diagZeroImpl(scip, conss, nconss, naddconss) );
            SCIPdebugMsg(scip, "Added %d constraints for implication from 0 diagonal.\n", *naddconss - noldaddconss);
            if ( noldaddconss != *naddconss )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added %d constraints for implications on SDP diagonals.\n", *naddconss - noldaddconss);
               *result = SCIP_SUCCESS;
            }
         }

         if ( *result != SCIP_CUTOFF && conshdlrdata->sdpconshdlrdata->twominorlinconss )
         {
            noldaddconss = *naddconss;
            SCIP_CALL( addTwoMinorLinConstraints(scip, conshdlr, conss, nconss, solvesdps, naddconss) );
            SCIPdebugMsg(scip, "Added %d linear constraints for 2 by 2 minors.\n", *naddconss - noldaddconss);
            if ( noldaddconss != *naddconss )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added %d linear constraints based on 2 x 2 SDP-minors.\n", *naddconss - noldaddconss);
               *result = SCIP_SUCCESS;
            }
         }

         if ( *result != SCIP_CUTOFF && conshdlrdata->sdpconshdlrdata->twominorprodconss )
         {
            noldaddconss = *naddconss;
            SCIP_CALL( addTwoMinorProdConstraints(scip, conshdlr, conss, nconss, solvesdps, naddconss) );
            SCIPdebugMsg(scip, "Added %d linear constraints for products of 2 by 2 minors.\n", *naddconss - noldaddconss);
            if ( noldaddconss != *naddconss )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added %d linear constraints based on products of 2 x 2 SDP-minors.\n", *naddconss - noldaddconss);
               *result = SCIP_SUCCESS;
            }
         }

         /* add SOCP-approximation if required */
         if ( *result != SCIP_CUTOFF && conshdlrdata->sdpconshdlrdata->addsocrelax )
         {
            noldaddconss = *naddconss;
            SCIP_CALL( addTwoMinorSOCConstraints(scip, conss, nconss, solvesdps, naddconss) );
            SCIPdebugMsg(scip, "Added %d SOC constraints for 2 by 2 minors.\n", *naddconss - noldaddconss);
            if ( noldaddconss != *naddconss )
               *result = SCIP_SUCCESS;
         }
      }

      if ( SCIPgetSubscipDepth(scip) == 0 && *result != SCIP_CUTOFF && conshdlrdata->sdpconshdlrdata->quadconsrank1 )
      {
         noldaddconss = *naddconss;
         SCIP_CALL( addRank1QuadConss(scip, conshdlr, conss, nconss, naddconss) );
         SCIPdebugMsg(scip, "Added %d quadratic constraints for rank 1 constraints.\n", *naddconss - noldaddconss);
         if ( noldaddconss != *naddconss )
            *result = SCIP_SUCCESS;

         /* turn off upgrading in order to avoid upgrading to a rank-1 constraint again */
         conshdlrdata->sdpconshdlrdata->upgradequadconss = FALSE;
      }
   }

   /* add variable bounds based on 2x2 minors in final round */
   if ( SCIPisPresolveFinished(scip) && conshdlrdata->sdpconshdlrdata->twominorvarbounds && ! conshdlrdata->sdpconshdlrdata->triedvarbounds )
   {
      if ( SCIPgetSubscipDepth(scip) == 0 && *result != SCIP_CUTOFF )
      {
         int noldaddconss;
         int solvesdpsparam;
         SCIP_Bool solvesdps;

         SCIP_CALL( SCIPgetIntParam(scip, "misc/solvesdps", &solvesdpsparam) );

         if ( solvesdpsparam == 1 )
            solvesdps = TRUE;
         else
            solvesdps = FALSE;

         noldaddconss = *naddconss;
         SCIP_CALL( addTwoMinorVarBounds(scip, conshdlr, conss, nconss, solvesdps, naddconss) );
         SCIPdebugMsg(scip, "Added %d linear constraints for variables bounds from 2 by 2 minors.\n", *naddconss - noldaddconss);
         if ( noldaddconss != *naddconss )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Added %d linear constraints based on variable bounds from 2 x 2 SDP-minors.\n", *naddconss - noldaddconss);
            *result = SCIP_SUCCESS;
         }
         conshdlrdata->sdpconshdlrdata->triedvarbounds = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** creates transformed constraint */
static
SCIP_DECL_CONSTRANS(consTransSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif
   int i;
   char transname[SCIP_MAXSTRLEN];

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );

   SCIPdebugMsg(scip, "Transforming constraint <%s>\n", SCIPconsGetName(sourcecons));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

#ifdef OMP
   SCIPdebugMsg(scip, "Setting number of threads to %d via OpenMP in Openblas.\n", conshdlrdata->sdpconshdlrdata->nthreads);
   omp_set_num_threads(conshdlrdata->sdpconshdlrdata->nthreads);
#endif

   SCIP_CALL( SCIPallocBlockMemory(scip, &targetdata) );

   /* copy some general data */
   targetdata->nvars = sourcedata->nvars;
   targetdata->nnonz = sourcedata->nnonz;
   targetdata->blocksize = sourcedata->blocksize;
   targetdata->matrixvar = NULL;
   targetdata->matrixval = NULL;
   targetdata->matrixconst = NULL;
   targetdata->nsingle = 0;
   targetdata->propubpossible = TRUE;
   targetdata->tracebound = -2.0;
   targetdata->allmatricespsd = sourcedata->allmatricespsd;
   targetdata->initallmatricespsd = sourcedata->initallmatricespsd;

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->nvarnonz), sourcedata->nvarnonz, sourcedata->nvars) );

   /* copy the non-constant nonzeros */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->col), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->row), sourcedata->nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->val), sourcedata->nvars) );

   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->col[i]), sourcedata->col[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->row[i]), sourcedata->row[i], sourcedata->nvarnonz[i]) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->val[i]), sourcedata->val[i], sourcedata->nvarnonz[i]) );
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(targetdata->vars), sourcedata->nvars) );
   if ( sourcedata->locks != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->locks), sourcedata->locks, sourcedata->nvars) );
   }
   else
      targetdata->locks = NULL;

   /* copy & transform the vars array */
   for (i = 0; i < sourcedata->nvars; i++)
   {
      targetdata->vars[i] = SCIPvarGetTransVar(sourcedata->vars[i]);
      SCIP_CALL( SCIPcaptureVar(scip, targetdata->vars[i]) );
   }

   /* copy the constant nonzeros */
   targetdata->constnnonz = sourcedata->constnnonz;

   if ( sourcedata->constnnonz > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constcol), sourcedata->constcol, sourcedata->constnnonz));
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constrow), sourcedata->constrow, sourcedata->constnnonz));
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->constval), sourcedata->constval, sourcedata->constnnonz));
   }
   else
   {
      targetdata->constcol = NULL;
      targetdata->constrow = NULL;
      targetdata->constval = NULL;
   }

   /* copy the maxrhsentry */
   targetdata->maxrhsentry = sourcedata->maxrhsentry;

   /* copy uniqueness information */
   if ( sourcedata->issymunique != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->issymunique), sourcedata->issymunique, sourcedata->nvars) );
   }
   else
      targetdata->issymunique = NULL;

   /* copy rankone */
   targetdata->rankone = sourcedata->rankone;

   /* copy maxevsubmat */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->maxevsubmat), sourcedata->maxevsubmat, 2) );

   /* copy addedquadcons */
   targetdata->addedquadcons = sourcedata->addedquadcons;

   /* name the transformed constraint */
#ifndef NDEBUG
   snprintfreturn = SCIPsnprintf(transname, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
   assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
   (void) SCIPsnprintf(transname, SCIP_MAXSTRLEN, "t_%s", SCIPconsGetName(sourcecons));
#endif

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, transname, conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),  SCIPconsIsLocal(sourcecons),
         SCIPconsIsModifiable(sourcecons), SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons),
         SCIPconsIsStickingAtNode(sourcecons)) );

   /* we need to compute the DIMACS tolerance (if required) at this point, because it is needed in CONSCHECK */
   if ( conshdlrdata->sdpconshdlrdata->usedimacsfeastol )
   {
      SCIP_VAR** vars;
      SCIP_Real sum = 0.0;
      int nvars;
      int v;

      nvars = SCIPgetNOrigVars(scip);
      vars = SCIPgetOrigVars(scip);
      for ( v = 0; v < nvars; v++ )
         sum += REALABS( SCIPvarGetObj(vars[v]) );
      conshdlrdata->dimacsfeastol = 1e-5 * (1 + sum);
   }

   return SCIP_OKAY;
}

/** checks feasiblity of constraint, e.g., positive semidefiniteness */
static
SCIP_DECL_CONSCHECK(consCheckSdp)
{/*lint --e{715}*/
   SCIP_CONS* violcons;
   SCIP_CONSDATA* consdata;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** linvars;
   SCIP_VAR*  var;
   SCIP_SOL* bestrank1approx;
   SCIP_Real* fullmatrix;
   SCIP_Real* eigenvalues;
   SCIP_Real* eigenvectors;
   SCIP_Real* scaledeigenvectors;
   SCIP_Real* matrixC;
   SCIP_Real* matrixAj;
   SCIP_Real* linmatrix;
   SCIP_Real* rhsmatrix;
   SCIP_Real* lssolu;
   SCIP_Real* linvals;
   SCIP_Real* colmatrix;
   SCIP_Bool rank1result;
   SCIP_Bool stored;

   int c;
   int i;
   int j;
   int k;
   int l;
#ifdef PRINTMATRICES
   int r;
   int s;
#endif
   int idx;
   int blocksize;
   int nviolrank1 = 0;
   int nvars;
   int nrank1vars = 0;
   int linrows = 0; /* the number of rows of the linear equation system is given by the total number of entries
      (lower-triangular) in all violated rank-1 constraints altogether */
   int lincnt = 0;
   int nsdpvars;
   int* rank1considx;
   int* indviolrank1conss;
   SCIP_VAR** rank1consvars;

   assert( scip != NULL );
   assert( result != NULL );
   assert( conss != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *result = SCIP_FEASIBLE;

   /* early termination */
   if ( nconss == 0 )
      return SCIP_OKAY;

#ifdef PRINTMATRICES
   SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
#endif

   /* check positive semidefiniteness */
   for (i = 0; i < nconss && *result == SCIP_FEASIBLE; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conshdlrdata, conss[i], sol, printreason, result) );
#ifdef PRINTMATRICES
      SCIPinfoMessage(scip, NULL, "Solution is %d for constraint %s.\n", *result, SCIPconsGetName(conss[i]) );
#endif
   }

   /* if the SDP-constraint handler is the current constraint handler or we are already infeasible, we do not need to
    * check the rank1-part (if present) */
   if ( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 || *result == SCIP_INFEASIBLE )
      return SCIP_OKAY;

   /* otherwise, we need to check for rank1 */
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLRRANK1_NAME) == 0 && *result == SCIP_FEASIBLE );

   SCIP_CALL( SCIPallocBufferArray(scip, &indviolrank1conss, nconss) );

   for (i = 0; i < nconss; ++i)
   {
      consdata = SCIPconsGetData(conss[i]);
      assert( consdata != NULL );
      assert( consdata->rankone );

      /* If the quadratic constraints are not yet added, we need to check them manually. Otherwise, they are checked
       * by the quadratic constraint handler, and we do not need to check for rank-1 separately. */
      if ( conshdlrdata->sdpconshdlrdata->quadconsrank1 && ! consdata->addedquadcons )
      {
         /* check quadratic constraints manually */
         SCIP_CALL( checkRank1QuadConss(scip, conshdlrdata, conss[i], sol, printreason, &rank1result) );

#ifdef PRINTMATRICES
         SCIPinfoMessage(scip, NULL, "Solution is %d for rank-1 part of constraint %s.\n", rank1result, SCIPconsGetName(conss[i]) );
#endif
      }
      else if ( ! conshdlrdata->sdpconshdlrdata->quadconsrank1 )
      {
         /* We need to check for rank-1. */
         SCIP_CALL( isMatrixRankOne(scip, conss[i], sol, printreason, &rank1result) );
#ifdef PRINTMATRICES
         SCIPinfoMessage(scip, NULL, "Solution is %d for rank-1 part of constraint %s.\n", rank1result, SCIPconsGetName(conss[i]) );
#endif
      }
      else
         rank1result = TRUE;

      if ( ! rank1result )
      {
         /* save index of violated rank-1 constraint */
         indviolrank1conss[nviolrank1] = i;
         ++nviolrank1;
      }
   }

   /* If there are no (violated) rank-1 constraints, we are finished. Otherwise, try to compute a feasible primal
    * solution by computing the best rank-1 approximation for each violated rank-1 constraint and solve an LP to find a
    * solution for the appearing variables. */
   if ( nviolrank1 == 0 )
   {
      assert( *result == SCIP_FEASIBLE );
      SCIPdebugMsg(scip, "Found no violated rank-1 constraints.\n");
      SCIPfreeBufferArray(scip, &indviolrank1conss);
      return SCIP_OKAY;
   }
   else
   {
      assert( nviolrank1 > 0 );
      *result = SCIP_INFEASIBLE;
   }

   /* check if heuristic should be executed and do not run in sub-SCIPs to avoid recursive reformulations due to rank 1
    * constraints */
   if ( ! conshdlrdata->sdpconshdlrdata->rank1approxheur || SCIPgetSubscipDepth(scip) > 0 )
   {
      assert( nviolrank1 > 0 && *result == SCIP_INFEASIBLE );
      SCIPfreeBufferArray(scip, &indviolrank1conss);
      return SCIP_OKAY;
   }

   SCIPdebugMsg(scip, "Found %d violated rank-1 constraints, thus apply rank-1 approximation heuristic!\n", nviolrank1);

   /* we have to collect all variables appearing in violated SDPrank1-constraints first */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &rank1considx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rank1consvars, nvars) );

   for (j = 0; j < nvars; ++j)
      rank1considx[j] = -1;

   for (c = 0; c < nviolrank1; ++c)
   {
      assert( conss[indviolrank1conss[c]] != NULL );

      /* todo: write function to only get number of variables and variables of an SDP constraint */
      consdata = SCIPconsGetData(conss[indviolrank1conss[c]]);
      assert( consdata != NULL );

      nsdpvars = consdata->nvars;
      linrows += consdata->blocksize * (consdata->blocksize + 1) / 2;

      for (i = 0; i < nsdpvars; ++i)
      {
         var = consdata->vars[i];
         idx = SCIPvarGetProbindex(var);
         assert( 0 <= idx && idx < nvars );
         if ( rank1considx[idx] < 0 )
         {
            rank1consvars[nrank1vars] = var;
            rank1considx[idx] = nrank1vars++;
         }
      }
   }

   /* initialize matrix of linear equation system */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &linmatrix, linrows * nrank1vars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &rhsmatrix, MAX(linrows,nrank1vars)) );

   for (i = 0; i < nviolrank1; ++i)
   {
      /* get violated rank-1 constraint */
      violcons = conss[indviolrank1conss[i]];
      consdata = SCIPconsGetData(violcons);

      assert( consdata != NULL );

      blocksize = consdata->blocksize;

      SCIPdebugMsg(scip, "\n Start with violated rank-1 constraint %s, is %d out of %d violated rank-1 constraints.\n\n", SCIPconsGetName(violcons), i + 1, nviolrank1);

      /* allocate memory to store full matrix */
      SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, blocksize * blocksize ) );
      SCIP_CALL( SCIPallocBufferArray(scip, &eigenvalues, blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &eigenvectors, blocksize * blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &matrixC, blocksize * blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &matrixAj, blocksize * blocksize) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvars, consdata->nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linvals, consdata->nvars) );

      /* compute the matrix \f$ \sum_j A_j y_j - A_0 \f$ */
      SCIP_CALL( computeFullSdpMatrix(scip, consdata, sol, fullmatrix) );

#ifdef PRINTMATRICES
      /* SCIPSDP uses row-first format! */
      printf("Full SDP-constraint matrix Z: \n");
      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k < blocksize; ++k)
            printf("%.5f  ", fullmatrix[j*blocksize + k]);
         printf("\n");
      }

      /* Double-check that fullmatrix is in row-first format */
      printf("Full SDP-constraint matrix Z in row-first format: \n");
      for (j = 0; j < blocksize * blocksize; ++j)
         printf("%.5f  ", fullmatrix[j]);
      printf("\n");
#endif

      /* compute EVD */
      SCIP_CALL( SCIPlapackComputeEigenvectorDecomposition(SCIPbuffer(scip), blocksize, fullmatrix, eigenvalues, eigenvectors) );

#ifdef PRINTMATRICES
      /* caution: LAPACK uses column-first format! */
      printf("Eigenvectors of Z: \n");
      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k < blocksize; ++k)
            printf("%.5f  ", eigenvectors[k*blocksize + j]);
         printf("\n");
      }

      /* Double-check that eigenvectors is in column-first format */
      printf("Eigenvectors of Z in column-first format: \n");
      for (j = 0; j < blocksize * blocksize; ++j)
         printf("%.5f  ", eigenvectors[j]);
      printf("\n");

      printf("Eigenvalues of Z: \n");
      for (j = 0; j < blocksize; ++j)
         printf("%.5f ", eigenvalues[j]);
      printf("\n");
#endif

      /* duplicate memory of eigenvectors to compute diag(0,...,0,lambda_max) * U^T */
      SCIP_CALL( SCIPduplicateBufferArray(scip, &scaledeigenvectors, eigenvectors, blocksize*blocksize) );

      /* set all eigenvalues except the largest to zero (using the property that LAPACK returns them in ascending
         order) */
      for (j = 0; j < blocksize-1; ++j)
      {
         assert( ! SCIPisFeasNegative(scip, eigenvalues[j]) ); /* otherwise constraint is not psd */
         eigenvalues[j] = 0.0;
      }

      /* compute diag(0,...,0,lambda_max) * U^T. Note that scaledeigenvectors on entry is U in column-first format,
         i.e., U^T in row-first format. Thus, on exit, scaledeigenvectors contains diag(0,...,0,lambda_max) * U^T in
         row-first format. */
      SCIP_CALL( scaleRowsMatrix(blocksize, scaledeigenvectors, eigenvalues) );

#ifdef PRINTMATRICES
      printf("Scaled eigenvectors of Z (only keep largest eigenvalue and corresponding eigenvector) : \n");
      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k < blocksize; ++k)
         {
            printf("%.5f  ", scaledeigenvectors[j*blocksize + k]);
         }
         printf("\n");
      }

      /* Double-check that scaledeigenvectors is in row-first format */
      printf("Scaled eigenvectors of Z in row-first format: \n");
      for (j = 0; j < blocksize * blocksize; ++j)
         printf("%.5f  ", scaledeigenvectors[j]);
      printf("\n");
#endif

      /* compute U * [diag(0,...,0,lambda_max) * U^T]: Since eigenvectors, which contains U, already comes in LAPACK's
         column-first format, eigenvectors does not need to be transposed! scaledeigenvectors models
         [diag(0,...,0,lambda_max) * U^T in SCIPSDP's row-first format, so that scaledeigenvectors needs to be
         transposed for LAPACK! */
      SCIP_CALL( SCIPlapackMatrixMatrixMult(blocksize, blocksize, eigenvectors, FALSE, blocksize, blocksize, scaledeigenvectors,
            TRUE, fullmatrix) );

#ifdef PRINTMATRICES
      printf("Best rank-1 approximation of Z: \n");
      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k < blocksize; ++k)
            printf("%.5f  ", fullmatrix[j*blocksize + k]);
         printf("\n");
      }

      /* Double-check that fullmatrix (best rank-1 approximation) is in row-first format */
      printf("Best rank-1 approximation of Z in row-first format: \n");
      for (j = 0; j < blocksize * blocksize; ++j)
         printf("%.5f  ", fullmatrix[j]);
      printf("\n");
#endif

      /* update linear equation system */

      /* compute constant matrix A_0 in row-first format*/
      SCIP_CALL( SCIPconsSdpGetFullConstMatrix(scip, violcons, matrixC) );

#ifdef PRINTMATRICES
      printf("Constant matrix A_0 of SDP-constraint: \n");
      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k < blocksize; ++k)
            printf("%.5f  ", matrixC[j*blocksize + k]);
         printf("\n");
      }

      /* Double-check that matrixA0 is in row-first format */
      printf("Constant matrix A_0 of SDP-constraint in row-first format: \n");
      for (j = 0; j < blocksize * blocksize; ++j)
         printf("%.5f  ", matrixC[j]);
      printf("\n");
#endif

      for (j = 0; j < blocksize; ++j)
      {
         for (k = 0; k <= j; ++k)
         {
            for (l = 0; l < consdata->nvars; ++l)
            {
               /* compute matrix A_j in row-first format */
               SCIP_CALL( SCIPconsSdpGetFullAj(scip, violcons, l, matrixAj) );

#ifdef PRINTMATRICES
               printf("Coefficient matrix A_%d of SDP-constraint: \n", l+1);
               for (r = 0; r < blocksize; ++r)
               {
                  for (s = 0; s < blocksize; ++s)
                     printf("%.5f  ", matrixAj[r*blocksize + s]);
                  printf("\n");
               }

               /* Double-check that matrixAj is in row-first format */
               printf("Constant matrix A_0 of SDP-constraint in row-first format: \n");
               for (r = 0; r < blocksize * blocksize; ++r)
                  printf("%.5f  ", matrixAj[r]);
               printf("\n");
#endif

               idx = SCIPvarGetProbindex(consdata->vars[l]);
               idx = rank1considx[idx];
               assert( 0 <= idx && idx < nrank1vars );
               assert( lincnt <= linrows );
               linmatrix[lincnt * nrank1vars + idx] = matrixAj[j * blocksize + k];
            }
            rhsmatrix[lincnt] = matrixC[j * blocksize + k] + fullmatrix[j * blocksize + k];
            ++lincnt;
         }
      }

      /* free memory for full matrix, eigenvalues and eigenvectors */
      SCIPfreeBufferArray(scip, &linvals);
      SCIPfreeBufferArray(scip, &linvars);
      SCIPfreeBufferArray(scip, &matrixAj);
      SCIPfreeBufferArray(scip, &matrixC);
      SCIPfreeBufferArray(scip, &scaledeigenvectors);
      SCIPfreeBufferArray(scip, &eigenvectors);
      SCIPfreeBufferArray(scip, &eigenvalues);
      SCIPfreeBufferArray(scip, &fullmatrix);
   }

   assert( lincnt == linrows );

#ifdef PRINTMATRICES
   printf("Matrix for linear equation system, in row-first format:\n");
   for (j = 0; j < linrows; ++j)
   {
      for (k = 0; k < nrank1vars; ++k)
      {
         printf("%.5f  ", linmatrix[j * nrank1vars + k]);
      }
      printf("\n");
   }

   /* Double-check that linmatrix is in row-first format */
   printf("Matrix for linear equation system in row-first format: \n");
   for (r = 0; r < linrows * nrank1vars; ++r)
      printf("%.5f  ", linmatrix[r]);
   printf("\n");

   printf("Right-hand for linear equation system:\n");
   for (j = 0; j < nrank1vars; ++j)
   {
      printf("%.5f  ", rhsmatrix[j]);
   }
   printf("\n");
#endif

   /* solve linear equation system with LAPACK */
   SCIP_CALL( SCIPallocBufferArray(scip, &lssolu, nrank1vars) );

   /* caution: LAPACK wants matrices in columns-first format, but SCIPSDP represents matrices in row-first format */
   SCIP_CALL( SCIPallocBufferArray(scip, &colmatrix, linrows * nrank1vars ) );

   SCIP_CALL( convertRowToColFormatFullMatrix(linrows, nrank1vars, linmatrix, colmatrix) );

#ifdef PRINTMATRICES
   printf("Matrix for linear equation system, in col-first format:\n");
   for (j = 0; j < linrows; ++j)
   {
      for (l = 0; l < nrank1vars; ++l)
      {
         printf("%.5f  ", colmatrix[l * linrows + j]);
      }
      printf("\n");
   }

   /* Double-check that colmatrix is in col-first format */
   printf("Matrix for linear equation system in col-first format: \n");
   for (r = 0; r < linrows * nrank1vars; ++r)
      printf("%.5f  ", colmatrix[r]);
   printf("\n");
#endif

   SCIP_CALL( SCIPlapackLinearSolve( SCIPbuffer(scip), linrows, nrank1vars, colmatrix, rhsmatrix, lssolu) );

   /* copy current solution */
   SCIP_CALL( SCIPcreateSolCopy(scip, &bestrank1approx, sol) );

   /* update solution with values from linear equations system solution from LAPACK */
   for (i = 0; i < nrank1vars; ++i)
   {
      var = rank1consvars[i];
      SCIP_CALL( SCIPsetSolVal(scip, bestrank1approx, var, lssolu[i]) );
   }

   SCIP_CALL( SCIPtrySolFree(scip, &bestrank1approx, FALSE, TRUE, TRUE, TRUE, TRUE, &stored) );
   if ( stored )
      SCIPdebugMsg(scip, "Best Rank-1 Approximation Heuristic found feasible primal solution\n");
   else
      SCIPdebugMsg(scip, "Primal solution found by Best Rank-1 Approximation Heuristic is not feasible!\n");

   SCIPfreeBufferArray(scip, &colmatrix);
   SCIPfreeBufferArray(scip, &lssolu);
   SCIPfreeBufferArray(scip, &rhsmatrix);
   SCIPfreeBufferArray(scip, &linmatrix);
   SCIPfreeBufferArray(scip, &rank1consvars);
   SCIPfreeBufferArray(scip, &rank1considx);
   SCIPfreeBufferArray(scip, &indviolrank1conss);

   return SCIP_OKAY;
}

/** enforce pseudo solution method
 *
 *  Returns didnotrun if objinfeasible, computes feasibility otherwise.
 */
static
SCIP_DECL_CONSENFOPS(consEnfopsSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   assert( scip != NULL );
   assert( result != NULL );
   assert( conss != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( objinfeasible )
   {
      SCIPdebugMsg(scip, "-> pseudo solution is objective infeasible, return.\n");
      return SCIP_OKAY;
   }

   for (i = 0; i < nconss; ++i)
   {
      SCIP_CALL( SCIPconsSdpCheckSdpCons(scip, conshdlrdata, conss[i], NULL, FALSE, result) );

      if (*result == SCIP_INFEASIBLE)
      {
         /* if it is infeasible for one SDP constraint, it is infeasible for the whole problem */
         SCIPdebugMsg(scip, "-> pseudo solution infeasible for SDP-constraint %s, return.\n", SCIPconsGetName(conss[i]));
         return SCIP_OKAY;
      }
   }

   *result = SCIP_FEASIBLE;

   SCIPdebugMsg(scip, "-> pseudo solution feasible for all SDP-constraints.\n");

   return SCIP_OKAY;
}


/** Enforce lp solution; if some block is not psd, an eigenvector cut is added.
 */
static
SCIP_DECL_CONSENFOLP(consEnfolpSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   *result = SCIP_FEASIBLE;

   /* do not run if another constraint handler has declared the solution to be infeasible */
   if ( solinfeasible )
      return SCIP_OKAY;

   /* we first check whether the LP solution if feasible */
   for (c = 0; c < nconss && *result != SCIP_CUTOFF; ++c)
   {
      SCIP_RESULT separesult = SCIP_FEASIBLE;

      SCIP_CALL( separateSol(scip, conshdlr, conss[c], NULL, TRUE, &separesult) );
      assert( separesult == SCIP_FEASIBLE || separesult == SCIP_CUTOFF || separesult == SCIP_SEPARATED || separesult == SCIP_CONSADDED );

      if ( separesult == SCIP_CUTOFF )
         *result = SCIP_CUTOFF;
      else if ( separesult == SCIP_CONSADDED )
         *result = SCIP_CONSADDED;
      else if ( separesult == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;
   }

   /* Below, we enforce integral solutions. If the LP is unbounded, this might not be guaranteed due to the integrality
    * constraint handler. In this case, we exit. The same happens if no relaxation is available or if we reached a cutoff. */
   if ( SCIPgetLPSolstat(scip) == SCIP_LPSOLSTAT_UNBOUNDEDRAY || conshdlrdata->relaxsdp == NULL || *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   /* if all integer variables have integral values, then possibly solve SDP in addtion to separation */
   if ( conshdlrdata->sdpconshdlrdata->enforcesdp && (*result == SCIP_SEPARATED || *result == SCIP_CONSADDED) )
   {
      SCIP_Bool cutoff;
      SCIP_VAR** vars;
      int nfixed = 0;
      int nintvars;
      int v;

      /* all integer variables should have integer values, because enforcing is called after the integrality constraint handler */
      vars = SCIPgetVars(scip);
      nintvars = SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip);

      /* check integer variables which are fixed or have integeral values (integer variables are first in list) */
      for (v = 0; v < nintvars; ++v)
      {
         SCIP_VAR* var;

         var = vars[v];
         assert( SCIPvarIsIntegral(var) );

         assert( SCIPisFeasIntegral(scip, SCIPgetSolVal(scip, NULL, var)) );
         assert( SCIPisFeasIntegral(scip, SCIPvarGetLbLocal(var)) );
         assert( SCIPisFeasIntegral(scip, SCIPvarGetUbLocal(var)) );

         if ( SCIPvarGetLbLocal(var) + 0.5 > SCIPvarGetUbLocal(var) )
            ++nfixed;
      }

      /* solve SPD if either all integer variables are fixed or if required */
      if ( ! conshdlrdata->sdpconshdlrdata->onlyfixedintssdp || nfixed == nintvars )
      {
         /* start probing */
         SCIP_CALL( SCIPstartProbing(scip) );

         /* apply domain propagation (use parameter settings for maximal number of rounds) */
         SCIP_CALL( SCIPpropagateProbing(scip, 0, &cutoff, NULL) );
         if ( ! cutoff )
         {
            int freq;

            SCIPdebugMsg(scip, "Solving relaxation because all integer variable have integral values.\n");

            /* temporarily change relaxator frequency, since otherwise relaxation will not be solved */
            freq = SCIPrelaxGetFreq(conshdlrdata->relaxsdp);
            SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );

            /* solve SDP */
            SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );

            /* reset frequency of relaxator */
            SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", freq) );

            /* if solving was successfull */
            if ( SCIPrelaxSdpSolvedProbing(conshdlrdata->relaxsdp) && SCIPisRelaxSolValid(scip) )
            {
               /* if we are infeasible, we can cut off the node */
               if ( ! SCIPrelaxSdpIsFeasible(conshdlrdata->relaxsdp) )
               {
                  SCIPdebugMsg(scip, "Cut off node in enforcing, because remaining SDP was infeasible.\n");
                  *result = SCIP_CUTOFF;
               }
               else
               {
                  SCIP_SOL* enfosol;
                  SCIP_Bool feasible;

                  /* if we are feasible, we check whether the solution is valid */
                  SCIP_CALL( SCIPcreateSol(scip, &enfosol, NULL) );
                  SCIP_CALL( SCIPlinkRelaxSol(scip, enfosol) );

                  /* check all constraints, including integrality. Since there is an integral constraint handler,
                   * integrality gets checked as well, so we don't need to do this manually. */
                  SCIP_CALL( SCIPcheckSol(scip, enfosol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );

                  if ( feasible )
                  {
                     SCIP_Bool stored;

                     /* pass solution to SCIP */
                     SCIP_CALL( SCIPaddSol(scip, enfosol, &stored) );

                     /* if all integer variables are fixed, we can cut off the node, since this should be the final solution of this node */
                     if ( stored && nintvars == nfixed )
                        *result = SCIP_CUTOFF;
                  }
                  else if ( nfixed < nintvars )
                  {
                     /* if we have not obtained a feasible solution, we try to round it */
                     assert( *result != SCIP_CUTOFF );

                     for (v = 0; v < nintvars; ++v)
                     {
                        SCIP_Real val;
                        SCIP_VAR* var;

                        var = vars[v];
                        assert( SCIPvarIsIntegral(var) );

                        /* fix all unfixed variables */
                        if ( SCIPvarGetLbLocal(var) + 0.5 < SCIPvarGetUbLocal(var) )
                        {
                           /* round relaxation value to next integer value */
                           val = SCIPfeasFloor(scip, SCIPgetRelaxSolVal(scip, var) + 0.5);

                           if ( ! SCIPisEQ(scip, val, SCIPvarGetUbLocal(var)) )
                           {
                              SCIP_CALL( SCIPchgVarUbProbing(scip, var, val) );
                           }

                           if ( ! SCIPisEQ(scip, val, SCIPvarGetLbLocal(var)) )
                           {
                              SCIP_CALL( SCIPchgVarLbProbing(scip, var, val) );
                           }
                        }
                     }

                     /* solve SDP again */
                     SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", 1) );
                     SCIP_CALL( SCIPsolveProbingRelax(scip, &cutoff) );
                     SCIP_CALL( SCIPsetIntParam(scip, "relaxing/SDP/freq", freq) );

                     /* if solving was successfull */
                     if ( SCIPrelaxSdpSolvedProbing(conshdlrdata->relaxsdp) && SCIPisRelaxSolValid(scip) )
                     {
                        if ( SCIPrelaxSdpIsFeasible(conshdlrdata->relaxsdp) )
                        {
                           /* if we are feasible, we check whether the solution is valid */
                           assert( enfosol != NULL );
                           SCIP_CALL( SCIPlinkRelaxSol(scip, enfosol) );

                           /* Pass solution to SCIP: check all constraints, including integrality. Since there is an integral
                            * constraint handler, integrality gets checked as well, so we don't need to do this manually. */
                           SCIP_CALL( SCIPtrySol(scip, enfosol, FALSE, TRUE, TRUE, TRUE, TRUE, &feasible) );
                        }
                     }
                  }

                  /* We have to invalidate the relaxation solution, because SCIP will otherwise not check the relaxation solution for feasibility. */
                  SCIP_CALL( SCIPmarkRelaxSolInvalid(scip) );
                  SCIP_CALL( SCIPfreeSol(scip, &enfosol) );
               }
            }
         }

         /* free local problem */
         SCIP_CALL( SCIPendProbing(scip) );
      }
   }

   return SCIP_OKAY;
}

/** Enforce relaxation solution; if some block is not psd, an eigenvector cut is added. */
static
SCIP_DECL_CONSENFORELAX(consEnforelaxSdp)
{/*lint --e{715}*/
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( conss != NULL );
   assert( result != NULL );

   *result = SCIP_FEASIBLE;

   if ( solinfeasible )
      return SCIP_OKAY;

   /* Rounding errors might lead to infeasible relaxation solutions. We therefore perform a separation round in the hope
    * that this can resolve the problem. */
   for (c = 0; c < nconss && *result != SCIP_CUTOFF; ++c)
   {
      SCIP_RESULT separesult = SCIP_FEASIBLE;

      SCIP_CALL( separateSol(scip, conshdlr, conss[c], sol, TRUE, &separesult) );
      assert( separesult == SCIP_FEASIBLE || separesult == SCIP_CUTOFF || separesult == SCIP_SEPARATED || separesult == SCIP_CONSADDED );

      if ( separesult == SCIP_CUTOFF )
         *result = SCIP_CUTOFF;
      else if ( separesult == SCIP_CONSADDED )
         *result = SCIP_CONSADDED;
      else if ( separesult == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}

/** separates a solution using constraint specific ideas, gives cuts to SCIP */
static
SCIP_DECL_CONSSEPASOL(consSepasolSdp)
{/*lint --e{715}*/
   int i;

   assert( result != NULL );
   *result = SCIP_DIDNOTFIND;

   for (i = 0; i < nusefulconss && *result != SCIP_CUTOFF; ++i)
   {
      SCIP_RESULT separesult = SCIP_DIDNOTFIND;

      SCIP_CALL( separateSol(scip, conshdlr, conss[i], sol, FALSE, &separesult) );
      assert( separesult == SCIP_DIDNOTFIND || separesult == SCIP_CUTOFF || separesult == SCIP_SEPARATED || separesult == SCIP_CONSADDED );

      if ( separesult == SCIP_CUTOFF )
         *result = SCIP_CUTOFF;
      else if ( separesult == SCIP_CONSADDED )
         *result = SCIP_CONSADDED;
      else if ( separesult == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}

/** separation method of constraint handler for LP solution */
static
SCIP_DECL_CONSSEPALP(consSepalpSdp)
{/*lint --e{715}*/
   int i;

   assert( result != NULL );
   *result = SCIP_DIDNOTFIND;

   for (i = 0; i < nusefulconss && *result != SCIP_CUTOFF; ++i)
   {
      SCIP_RESULT separesult = SCIP_DIDNOTFIND;

      SCIP_CALL( separateSol(scip, conshdlr, conss[i], NULL, FALSE, &separesult) );

      assert( separesult == SCIP_DIDNOTFIND || separesult == SCIP_CUTOFF || separesult == SCIP_SEPARATED || separesult == SCIP_CONSADDED );

      if ( separesult == SCIP_CUTOFF )
         *result = SCIP_CUTOFF;
      else if ( separesult == SCIP_CONSADDED )
         *result = SCIP_CONSADDED;
      else if ( separesult == SCIP_SEPARATED )
         *result = SCIP_SEPARATED;
   }

   return SCIP_OKAY;
}

/** delete method of SDP constrainthandler */
static
SCIP_DECL_CONSDELETE(consDeleteSdp)
{/*lint --e{715}*/
   int i;

   assert( cons != NULL );
   assert( consdata != NULL );

   SCIPdebugMsg(scip, "deleting SDP constraint <%s>.\n", SCIPconsGetName(cons));

   /* release memory for rank one constraint */
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->maxevsubmat, 2);

   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row[i], (*consdata)->nvarnonz[i]);
      SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col[i], (*consdata)->nvarnonz[i]);
   }

   /* release all variables */
   for (i = 0; i < (*consdata)->nvars; i++)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &((*consdata)->vars[i])) );
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->vars, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->locks, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constval, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constrow, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->constcol, (*consdata)->constnnonz);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->val, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->row, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->col, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->nvarnonz, (*consdata)->nvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->matrixval, (*consdata)->blocksize * ((*consdata)->blocksize + 1) / 2);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->matrixvar, (*consdata)->blocksize * ((*consdata)->blocksize + 1) / 2);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->matrixconst, (*consdata)->blocksize * ((*consdata)->blocksize + 1) / 2);
   SCIPfreeBlockMemoryArrayNull(scip, &(*consdata)->issymunique, (*consdata)->nvars);
   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/** free method of SDP constrainthandler */
static
SCIP_DECL_CONSFREE(consFreeSdp)
{/*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIPdebugMsg(scip, "Freeing constraint handler <%s>.\n", SCIPconshdlrGetName(conshdlr));

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   if ( conshdlrdata->randnumgen != NULL )
   {
      SCIPfreeRandom(scip, &conshdlrdata->randnumgen);
   }

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}

/** copy SDP constraint handler */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySdp)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConshdlrSdp(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** copy rank 1 SDP constraint handler */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopySdpRank1)
{
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLRRANK1_NAME) == 0);

   SCIP_CALL( SCIPincludeConshdlrSdpRank1(scip) );

   *valid = TRUE;

   return SCIP_OKAY;
}

/** copy an SDP constraint */
static
SCIP_DECL_CONSCOPY(consCopySdp)
{/*lint --e{715}*/
   char copyname[SCIP_MAXSTRLEN];
   SCIP_CONSDATA* sourcedata;
   SCIP_Bool success;
   SCIP_VAR** targetvars;
   SCIP_VAR* var;
   int i;
#ifndef NDEBUG
   int snprintfreturn; /* used to check the return code of snprintf */
#endif

   assert( scip != NULL );
   assert( sourcescip != NULL );
   assert( sourcecons != NULL );
   assert( valid != NULL );

   SCIPdebugMsg(scip, "Copying SDP constraint <%s>\n", SCIPconsGetName(sourcecons));

   *valid = TRUE;

   /* As we can only map active variables, we have to make sure, that the constraint contains no fixed or
    * (multi-)aggregated vars, before presolving and after presolving this should always be the case,
    * earlier than that we need to call fixAndAggrVars. */
   if ( SCIPgetStage(sourcescip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(sourcescip) <= SCIP_STAGE_EXITPRESOLVE )
   {
      SCIP_CALL( fixAndAggrVars(sourcescip, &sourcecons, 1, TRUE) );
   }

   sourcedata = SCIPconsGetData(sourcecons);
   assert( sourcedata != NULL );
   assert( sourcedata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLR_NAME) == 0 );
   assert( ! sourcedata->rankone || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(sourcecons)), CONSHDLRRANK1_NAME) == 0 );

   SCIP_CALL( SCIPallocBufferArray(scip, &targetvars, sourcedata->nvars) );

   /* map all variables in the constraint */
   for (i = 0; i < sourcedata->nvars; i++)
   {
      SCIP_CALL( SCIPgetVarCopy(sourcescip, scip, sourcedata->vars[i], &var, varmap, consmap, global, &success) );
      if ( success )
         targetvars[i] = var;
      else
         *valid = FALSE;
   }

   /* name the copied constraint */
#ifndef NDEBUG
   snprintfreturn = SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", name == NULL ? SCIPconsGetName(sourcecons) : name);
   assert( snprintfreturn < SCIP_MAXSTRLEN ); /* check whether the name fits into the string */
#else
   (void) SCIPsnprintf(copyname, SCIP_MAXSTRLEN, "c_%s", name == NULL ? SCIPconsGetName(sourcecons) : name);
#endif

   /* create the new constraint */
   if ( ! sourcedata->rankone )
   {
      SCIP_CALL( SCIPcreateConsSdp(scip, cons, copyname, sourcedata->nvars, sourcedata->nnonz, sourcedata->blocksize, sourcedata->nvarnonz,
            sourcedata->col, sourcedata->row, sourcedata->val, targetvars, sourcedata->constnnonz,
            sourcedata->constcol, sourcedata->constrow, sourcedata->constval, FALSE) );
   }
   else
   {
      SCIP_CALL( SCIPcreateConsSdpRank1(scip, cons, copyname, sourcedata->nvars, sourcedata->nnonz, sourcedata->blocksize, sourcedata->nvarnonz,
            sourcedata->col, sourcedata->row, sourcedata->val, targetvars, sourcedata->constnnonz,
            sourcedata->constcol, sourcedata->constrow, sourcedata->constval, FALSE) );
   }

   /* copy lock information if available */
   if ( sourcedata->locks != NULL )
   {
      SCIP_CONSDATA* targetdata;

      targetdata = SCIPconsGetData(*cons);
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->locks), sourcedata->locks, sourcedata->nvars) );
      targetdata->allmatricespsd = sourcedata->allmatricespsd;
      targetdata->initallmatricespsd = sourcedata->initallmatricespsd;
   }

   /* copy issymunique information if available */
   if ( sourcedata->issymunique != NULL )
   {
      SCIP_CONSDATA* targetdata;

      targetdata = SCIPconsGetData(*cons);
      assert( targetdata->issymunique == NULL );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(targetdata->issymunique), sourcedata->issymunique, sourcedata->nvars) );
   }

   SCIPfreeBufferArray(scip, &targetvars);

   return SCIP_OKAY;
}

/** print an SDP constraint */
static
SCIP_DECL_CONSPRINT(consPrintSdp)
{/*lint --e{715}*/
#ifdef PRINT_HUMAN_READABLE
   SCIP_CONSDATA* consdata;
   SCIP_Real* fullmatrix;
   int v;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &fullmatrix, consdata->blocksize * consdata->blocksize) );

   /* print rank1 information */
   SCIPinfoMessage(scip, file, "rank-1? %d\n", consdata->rankone);

   /* print the non-constant matrices, for this they first have to be assembled in fullmatrix */
   for (v = 0; v < consdata->nvars; v++)
   {
      /* assemble the matrix */

      /* first initialize it with zero */
      for (i = 0; i < consdata->blocksize; i++)
      {
         for (j = 0; j < consdata->blocksize; j++)
            fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
      }

      /* then add the nonzeros */
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         fullmatrix[consdata->row[v][i] * consdata->blocksize + consdata->col[v][i]] = consdata->val[v][i]; /* lower triangular entry */ /*lint !e679*/
         fullmatrix[consdata->col[v][i] * consdata->blocksize + consdata->row[v][i]] = consdata->val[v][i]; /* upper triangular entry */ /*lint !e679*/
      }

      /* print it */
      SCIPinfoMessage(scip, file, "+\n");
      for (i = 0; i < consdata->blocksize; i++)
      {
         SCIPinfoMessage(scip, file, "( ");
         for (j = 0; j < consdata->blocksize; j++)
            SCIPinfoMessage(scip, file, "%g ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
         SCIPinfoMessage(scip, file, ")\n");
      }
      SCIPinfoMessage(scip, file, "* %s\n", SCIPvarGetName(consdata->vars[v]));
   }

   /* print the constant-matrix */

   /* assemble the matrix */

   /* first initialize it with zero */
   for (i = 0; i < consdata->blocksize; i++)
   {
      for (j = 0; j < consdata->blocksize; j++)
         fullmatrix[i * consdata->blocksize + j] = 0.0; /*lint !e679*/
   }

   /* then add the nonzeros */
   for (i = 0; i < consdata->constnnonz; i++)
   {
      fullmatrix[consdata->constrow[i] * consdata->blocksize + consdata->constcol[i]] = consdata->constval[i]; /* lower triangular entry */ /*lint !e679*/
      fullmatrix[consdata->constcol[i] * consdata->blocksize + consdata->constrow[i]] = consdata->constval[i]; /* upper triangular entry */ /*lint !e679*/
   }

   /* print it */
   SCIPinfoMessage(scip, file, "-\n");
   for (i = 0; i < consdata->blocksize; i++)
   {
      SCIPinfoMessage(scip, file, "( ");
      for (j = 0; j < consdata->blocksize; j++)
         SCIPinfoMessage(scip, file, "%g ", fullmatrix[i * consdata->blocksize + j]); /*lint !e679*/
      SCIPinfoMessage(scip, file, ")\n");
   }
   SCIPinfoMessage(scip, file, ">= 0\n");

   /* print rank1 information */
   SCIPinfoMessage(scip, file, "rank-1? %d\n", consdata->rankone);

   SCIPfreeBufferArray(scip, &fullmatrix);

   return SCIP_OKAY;
#else
   SCIP_CONSDATA* consdata;
   int i;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);

   /* print blocksize */
   SCIPinfoMessage(scip, file, "%d\n", consdata->blocksize);

   /* print rank1 information */
   SCIPinfoMessage(scip, file, "    rank-1? %u\n", consdata->rankone);

   /* print A_0 if it exists */
   if ( consdata->constnnonz > 0 )
   {
      SCIPinfoMessage(scip, file, "    A_0: ");

      for (i = 0; i < consdata->constnnonz; i++)
      {
         if ( i < consdata->constnnonz - 1 )
            SCIPinfoMessage(scip, file, "(%d,%d):%.15g, ", consdata->constrow[i], consdata->constcol[i], consdata->constval[i]);
         else
            SCIPinfoMessage(scip, file, "(%d,%d):%.15g", consdata->constrow[i], consdata->constcol[i], consdata->constval[i]);
      }
      SCIPinfoMessage(scip, file, "\n");
   }

   /* print other matrices */
   for (v = 0; v < consdata->nvars; v++)
   {
      SCIPinfoMessage(scip, file, "    <%s>: ", SCIPvarGetName(consdata->vars[v]));
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         if ( i < consdata->nvarnonz[v] - 1 || v < consdata->nvars - 1 )
            SCIPinfoMessage(scip, file, "(%d,%d):%.15g, ", consdata->row[v][i], consdata->col[v][i], consdata->val[v][i]);
         else
            SCIPinfoMessage(scip, file, "(%d,%d):%.15g", consdata->row[v][i], consdata->col[v][i], consdata->val[v][i]);
      }
      /* if this is not the last variable, add a newline */
      if (v < consdata->nvars - 1)
      {
         SCIPinfoMessage(scip, file, "\n");
      }
   }

   return SCIP_OKAY;
#endif
}

/** parse an SDP constraint */
static
SCIP_DECL_CONSPARSE(consParseSdp)
{  /*lint --e{715}*/
   SCIP_Bool parsesuccess;
   SCIP_CONSDATA* consdata = NULL;
   char* pos;
   int currentsize;
   int nvars;
   int i;
   int v;
   int rankoneint;

   assert( scip != NULL );
   assert( str != NULL );

   nvars = SCIPgetNVars(scip);

   assert( success != NULL );
   *success = TRUE;

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   consdata->nvars = 0;
   consdata->nnonz = 0;
   consdata->constnnonz = 0;
   consdata->rankone = 0;
   consdata->addedquadcons = FALSE;
   consdata->locks = NULL;
   consdata->matrixvar = NULL;
   consdata->matrixval = NULL;
   consdata->matrixconst = NULL;
   consdata->nsingle = 0;
   consdata->propubpossible = TRUE;
   consdata->diagconstantone = FALSE;
   consdata->tracebound = -2.0;
   consdata->allmatricespsd = FALSE;
   consdata->initallmatricespsd = FALSE;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars));

   consdata->constcol = NULL;
   consdata->constrow = NULL;
   consdata->constval = NULL;

   /* parse the blocksize */
   parsesuccess = SCIPstrToIntValue(str, &(consdata->blocksize), &pos);
   *success = *success && parsesuccess;

   /* skip whitespace */
   while ( isspace((unsigned char)*pos) )
      pos++;

   /* parse the rank1-information */
   if ( pos[0] == 'r' && pos[1] == 'a' && pos[2] == 'n' && pos[3] == 'k' && pos[4] == '-' && pos[5] == '1' && pos[6] == '?' )
   {
      pos += 8;                 /* we skip "rank1-? " */
      parsesuccess = SCIPstrToIntValue(pos, &rankoneint, &pos);
      consdata->rankone = (SCIP_Bool) rankoneint;
      *success = *success && parsesuccess;
   }

   /* skip whitespace */
   while( isspace((unsigned char)*pos) )
      pos++;

   /* check if there is a constant part */
   if ( pos[0] == 'A' && pos[1] == '_' && pos[2] == '0' )
   {
      pos += 5; /* we skip "A_0: " */

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, PARSE_STARTSIZE) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, PARSE_STARTSIZE) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, PARSE_STARTSIZE) );

      currentsize = PARSE_STARTSIZE;

      /* as long as there is another entry for the constant part, parse it */
      while (pos[0] == '(')
      {
         pos++; /* remove the '(' */

         /* check if we need to enlarge the arrays */
         if ( consdata->constnnonz == currentsize )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, currentsize, PARSE_SIZEFACTOR * currentsize) );
            currentsize *= PARSE_SIZEFACTOR;
         }

         parsesuccess = SCIPstrToIntValue(pos, &(consdata->constrow[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->constrow[consdata->constnnonz] < consdata->blocksize );
         pos++; /* remove the ',' */
         parsesuccess = SCIPstrToIntValue(pos, &(consdata->constcol[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->constcol[consdata->constnnonz] < consdata->blocksize );
         pos += 2; /* remove the "):" */
         parsesuccess = SCIPstrToRealValue(pos, &(consdata->constval[consdata->constnnonz]), &pos);
         *success = *success && parsesuccess;
         pos ++; /* remove the "," */

         /* if we got an entry in the upper triangular part, switch the entries for lower triangular */
         if ( consdata->constcol[consdata->constnnonz] > consdata->constrow[consdata->constnnonz] )
         {
            i = consdata->constcol[consdata->constnnonz];
            consdata->constcol[consdata->constnnonz] = consdata->constrow[consdata->constnnonz];
            consdata->constrow[consdata->constnnonz] = i;
         }

         consdata->constnnonz++;

         /* skip whitespace */
         while( isspace((unsigned char)*pos) )
            pos++;
      }

      /* resize the arrays to their final size */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, currentsize, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, currentsize, consdata->constnnonz) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, currentsize, consdata->constnnonz) );
   }

   /* skip whitespace */
   while( isspace((unsigned char)*pos) )
      pos++;

   /* parse the non-constant part */

   /* while there is another variable, parse it */
   while ( pos[0] == '<' )
   {
      /* add the variable to consdata->vars and create the corresponding nonzero arrays */
      SCIP_CALL( SCIPparseVarName(scip, pos, &(consdata->vars[consdata->nvars]), &pos) );
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[consdata->nvars]) );

      consdata->nvarnonz[consdata->nvars] = 0;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->col[consdata->nvars]), PARSE_STARTSIZE));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->row[consdata->nvars]), PARSE_STARTSIZE));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(consdata->val[consdata->nvars]), PARSE_STARTSIZE));
      consdata->nvars++;
      currentsize = PARSE_STARTSIZE;

      pos += 2; /* remove the ": " */

      /* while there is another entry, parse it */
      while (pos[0] == '(')
      {
         pos++; /* remove the '(' */

         /* check if we need to enlarge the arrays */
         if ( consdata->nvarnonz[consdata->nvars - 1] == currentsize )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[consdata->nvars - 1], currentsize, PARSE_SIZEFACTOR * currentsize) );
            currentsize *= PARSE_SIZEFACTOR;
         }

         parsesuccess = SCIPstrToIntValue(pos, &(consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] < consdata->blocksize );
         pos++; /* remove the ',' */
         parsesuccess = SCIPstrToIntValue(pos, &(consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         assert( consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] < consdata->blocksize );
         pos += 2; /* remove the "):" */
         parsesuccess = SCIPstrToRealValue(pos, &(consdata->val[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]]), &pos);
         *success = *success && parsesuccess;
         if ( *pos != '\0' )
            pos ++; /* remove the "," */

         /* if we got an entry in the upper triangular part, switch the entries for lower triangular */
         if ( consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] >
               consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] )
         {
            i = consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]];
            consdata->col[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] =
                  consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]];
            consdata->row[consdata->nvars - 1][consdata->nvarnonz[consdata->nvars - 1]] = i;
         }

         consdata->nvarnonz[consdata->nvars - 1]++;

         /* skip whitespace */
         while( isspace((unsigned char)*pos) )
            pos++;
      }

      /* resize the arrays to their final size */
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[consdata->nvars - 1], currentsize, consdata->nvarnonz[consdata->nvars - 1]) );

      /* skip whitespace */
      while ( isspace((unsigned char)*pos) )
         pos++;
   }

   /* resize the arrays to their final size */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars, consdata->nvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col, nvars, consdata->nvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row, nvars, consdata->nvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val, nvars, consdata->nvars) );
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->vars, nvars, consdata->nvars));

   /* compute sdpnnonz */
   for (v = 0; v < consdata->nvars; v++)
      consdata->nnonz += consdata->nvarnonz[v];

   /* set maxevsubmat */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->maxevsubmat, 2) );
   consdata->maxevsubmat[0] = -1;
   consdata->maxevsubmat[1] = -1;

   /* create the constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate, local, modifiable,
         dynamic, removable, stickingatnode) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

#ifdef SCIP_MORE_DEBUG
   SCIP_CALL( SCIPprintCons(scip, *cons, NULL) );
#endif

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the variables */
static
SCIP_DECL_CONSGETVARS(consGetVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int nvars;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( vars != NULL );
   assert( success != NULL );
   assert( varssize >= 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   nvars = consdata->nvars;

   if ( nvars > varssize )
   {
      SCIPdebugMsg(scip, "consGetVarsIndicator called for array of size %d, needed size %d.\n", varssize, nvars);
      *success = FALSE;
      return SCIP_OKAY;
   }

   for (i = 0; i < nvars; i++)
      vars[i] = consdata->vars[i];

   *success = TRUE;

   return SCIP_OKAY;
}

/** constraint method of constraint handler which returns the number of variables */
static
SCIP_DECL_CONSGETNVARS(consGetNVarsSdp)
{/*lint --e{715}*/
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( success != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   *nvars = consdata->nvars;
   *success = TRUE;

   return SCIP_OKAY;
}

#if SCIP_VERSION >= 900
/** constraint handler method which returns the permutation symmetry detection graph of a constraint */
static
SCIP_DECL_CONSGETPERMSYMGRAPH(consGetPermsymGraphSdp)
{  /*lint --e{715}*/
   SCIP_CALL( addSymmetryInformation(scip, conshdlr, SYM_SYMTYPE_PERM, cons, graph, success) );

   return SCIP_OKAY;
}

/** constraint handler method which returns the signed permutation symmetry detection graph of a constraint */
static
SCIP_DECL_CONSGETSIGNEDPERMSYMGRAPH(consGetSignedPermsymGraphSdp)
{
   *success = FALSE;
   return SCIP_OKAY;
}
#endif

/** creates the handler for SDP constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != NULL );

   /* allocate memory for the conshdlrdata */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   conshdlrdata->quadconsidx = NULL;
   conshdlrdata->quadconsvars = NULL;
   conshdlrdata->nquadconsidx = 0;
   conshdlrdata->X = NULL;
   conshdlrdata->nsdpvars = 0;
   conshdlrdata->sdpcons = NULL;
   conshdlrdata->triedlinearconss = FALSE;
   conshdlrdata->triedvarbounds = FALSE;
   conshdlrdata->randnumgen = NULL;
   conshdlrdata->relaxsdp = NULL;
   conshdlrdata->sdpconshdlrdata = conshdlrdata;  /* set this to itself to simplify access of parameters */
   conshdlrdata->dimacsfeastol = SCIP_INVALID;
   conshdlrdata->ncallspropub = 0;
   conshdlrdata->ncallsproptb = 0;
   conshdlrdata->ncallsprop3minor = 0;
   conshdlrdata->propubtime = NULL;
   conshdlrdata->proptbtime = NULL;
   conshdlrdata->prop3minortime = NULL;
   conshdlrdata->maxtimepropub = 0.0;
   conshdlrdata->maxtimeproptb = 0.0;
   conshdlrdata->maxtimeprop3minor = 0.0;
   conshdlrdata->npropub = 0;
   conshdlrdata->nproptb = 0;
   conshdlrdata->nprop3minor = 0;
   conshdlrdata->npropcutoffub = 0;
   conshdlrdata->npropcutofftb = 0;
   conshdlrdata->npropcutoff3m = 0;
   conshdlrdata->npropintrndub = 0;
   conshdlrdata->npropintrndtb = 0;
   conshdlrdata->nproppreub = 0;
   conshdlrdata->nproppretb = 0;
   conshdlrdata->nproppre3m = 0;
   conshdlrdata->nproppreintrndub = 0;
   conshdlrdata->nproppreintrndtb = 0;
   conshdlrdata->nproppreintrnd3m = 0;
   conshdlrdata->npropprobub = 0;
   conshdlrdata->npropprobtb = 0;


   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSdp, consEnfopsSdp, consCheckSdp, consLockSdp, conshdlrdata) );

   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSdp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSdp) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySdp, consCopySdp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitSdp) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSdp, CONSHDLR_PROPFREQ, FALSE, CONSHDLR_PROPTIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSdp) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSdp) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSdp) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSdp) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSdp) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSdp) );
#if SCIP_VERSION >= 900
   SCIP_CALL( SCIPsetConshdlrGetPermsymGraph(scip, conshdlr, consGetPermsymGraphSdp) );
   SCIP_CALL( SCIPsetConshdlrGetSignedPermsymGraph(scip, conshdlr, consGetSignedPermsymGraphSdp) );
#endif

   /* add parameter */
#ifdef OMP
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/threads", "number of threads used for OpenBLAS",
         &(conshdlrdata->nthreads), TRUE, DEFAULT_NTHREADS, 1, INT_MAX, NULL, NULL) );
#endif

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/propupperbounds",
         "Should upper bounds be propagated?",
         &(conshdlrdata->propupperbounds), TRUE, DEFAULT_PROPUPPERBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/propubpresol",
         "Should upper bounds be propagated in presolving?",
         &(conshdlrdata->propubpresol), TRUE, DEFAULT_PROPUBPRESOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/prop3minors",
         "Should 3x3 minors be propagated?",
         &(conshdlrdata->prop3minors), TRUE, DEFAULT_PROP3MINORS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/nonconst3minors",
         "Should 3x3 minors be propagated if the diagonal is not constant?",
         &(conshdlrdata->nonconst3minors), TRUE, DEFAULT_NONCONST3MINORS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/prop3mprobing",
         "Should 3x3 minors be propagated in probing?",
         &(conshdlrdata->prop3mprobing), TRUE, DEFAULT_PROP3MPROBING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/proptightenbounds",
         "Should tighten bounds be propagated?",
         &(conshdlrdata->proptightenbounds), TRUE, DEFAULT_PROPTIGHTENBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/proptbprobing",
         "Should tighten bounds be propagated in probing?",
         &(conshdlrdata->proptbprobing), TRUE, DEFAULT_PROPTBPROBING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/tightenboundscont",
         "Should only bounds be tightend for continuous variables?",
         &(conshdlrdata->tightenboundscont), TRUE, DEFAULT_TIGHTENBOUNDSCONT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/tightenmatrices",
         "If all matrices are psd, should the matrices be tightened if possible?",
         &(conshdlrdata->tightenmatrices), TRUE, DEFAULT_TIGHTENMATRICES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/tightenbounds",
         "If all matrices are psd, should the bounds be tightened if possible?",
         &(conshdlrdata->tightenbounds), TRUE, DEFAULT_TIGHTENBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/diaggezerocuts",
         "Should linear cuts enforcing the non-negativity of diagonal entries of SDP-matrices be added?",
         &(conshdlrdata->diaggezerocuts), TRUE, DEFAULT_DIAGGEZEROCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/diagzeroimplcuts",
         "Should linear cuts enforcing the implications of diagonal entries of zero in SDP-matrices be added?",
         &(conshdlrdata->diagzeroimplcuts), TRUE, DEFAULT_DIAGZEROIMPLCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/twominorlinconss",
         "Should linear cuts corresponding to 2 by 2 minors be added?",
         &(conshdlrdata->twominorlinconss), TRUE, DEFAULT_TWOMINORLINCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/twominorprodconss",
         "Should linear cuts corresponding to products of 2 by 2 minors be added?",
         &(conshdlrdata->twominorprodconss), TRUE, DEFAULT_TWOMINORPRODCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/twominorvarbounds",
         "Should linear cuts corresponding to variable bounds for 2 by 2 minors be added?",
         &(conshdlrdata->twominorvarbounds), TRUE, DEFAULT_TWOMINORVARBOUNDS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/quadconsrank1",
         "Should quadratic cons for 2x2 minors be added in the rank-1 case?",
         &(conshdlrdata->quadconsrank1), TRUE, DEFAULT_QUADCONSRANK1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/upgradequadconss",
         "Should quadratic constraints be upgraded to a rank 1 SDP?",
         &(conshdlrdata->upgradequadconss), TRUE, DEFAULT_UPGRADEQUADCONSS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/upgradekeepquad",
         "Should the quadratic constraints be kept in the problem after upgrading and the corresponding SDP constraint be added without the rank 1 constraint?",
         &(conshdlrdata->upgradekeepquad), TRUE, DEFAULT_UPGRADEKEEPQUAD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/separateonecut",
         "Should only one cut corresponding to the most negative eigenvalue be separated?",
         &(conshdlrdata->separateonecut), TRUE, DEFAULT_SEPARATEONECUT, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/cutstopool",
         "Should the cuts be added to the pool?",
         &(conshdlrdata->cutstopool), TRUE, DEFAULT_CUTSTOPOOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/sparsifycut",
         "Should the eigenvector cuts be sparsified?",
         &(conshdlrdata->sparsifycut), TRUE, DEFAULT_SPARSIFYCUT, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/SDP/sparsifyfactor",
         "target size for sparsification in relation to number of variables",
         &(conshdlrdata->sparsifyfactor), TRUE, DEFAULT_SPARSIFYFACTOR, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/sparsifytargetsize",
         "absolute target size for sparsification (-1: use sparsifyfactor instead)",
         &(conshdlrdata->sparsifytargetsize), TRUE, DEFAULT_SPARSIFYTARGETSIZE, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/multiplesparsecuts",
         "Should multiple sparsified eigenvector cuts be added?",
         &(conshdlrdata->multiplesparsecuts), TRUE, DEFAULT_MULTIPLESPARSECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/maxnsparsecuts",
         "maximal number of sparse eigenvector cuts that should be added (-1: no limit)",
         &(conshdlrdata->maxnsparsecuts), TRUE, DEFAULT_MAXNSPARSECUTS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/enforcesdp",
         "Solve SDP if we do lp-solving and have an integral solution in enforcing?",
         &(conshdlrdata->enforcesdp), TRUE, DEFAULT_ENFORCESDP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/onlyfixedintssdp",
         "Should solving an SDP only be applied if all integral variables are fixed (instead of having integral values)?",
         &(conshdlrdata->onlyfixedintssdp), TRUE, DEFAULT_ONLYFIXEDINTSSDP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/addsocrelax",
         "Should a relaxation of SOC constraints be added?",
         &(conshdlrdata->addsocrelax), TRUE, DEFAULT_ADDSOCRELAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/maxnvarsquadupgd",
         "maximal number of quadratic constraints and appearing variables so that the QUADCONSUPGD is performed",
         &(conshdlrdata->maxnvarsquadupgd), TRUE, DEFAULT_MAXNVARSQUADUPGD, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/rank1approxheur",
         "Should the heuristic that computes the best rank-1 approximation for a given solution be executed?",
         &(conshdlrdata->rank1approxheur), TRUE, DEFAULT_RANK1APPROXHEUR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/usedimacsfeastol",
         "Should a feasibility tolerance based on the DIMACS be used for computing negative eigenvalues?",
         &(conshdlrdata->usedimacsfeastol), TRUE, DEFAULT_USEDIMACSFEASTOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/generaterows",
         "Should rows be generated (constraints otherwise)?",
         &(conshdlrdata->generaterows), TRUE, DEFAULT_GENERATEROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/generatecmir",
         "Should CMIR cuts be generated?",
         &(conshdlrdata->generatecmir), TRUE, DEFAULT_GENERATECMIR, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/recomputesparseev",
         "Should the sparse eigenvalue returned from TPower be recomputed exactly by using Lapack for the corresponding submatrix?",
         &(conshdlrdata->recomputesparseev), TRUE, DEFAULT_RECOMPUTESPARSEEV, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/recomputeinitial",
         "Should the inital vector for TPower be computed each time before calling TPower (instead of using the original smallest eigenvector)?",
         &(conshdlrdata->recomputeinitial), TRUE, DEFAULT_RECOMPUTEINITIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/exacttrans",
         "Should the matrix be transformed with the exact maximal eigenvalue before calling TPower (instead of using estimate)?",
         &(conshdlrdata->exacttrans), TRUE, DEFAULT_EXACTTRANS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/SDP/presollinconssparam",
         "Parameters for linear constraints added during presolving: (0) propagate, if solving LPs also separate (1) initial and propagate, if solving LPs also separate, enforce and check",
         &(conshdlrdata->presollinconssparam), TRUE, DEFAULT_PRESOLLINCONSSPARAM, 0, 1, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/additionalstats",
         "Should additional statistics be output at the end?",
         &(conshdlrdata->additionalstats), TRUE, DEFAULT_ADDITIONALSTATS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/enableproptiming",
         "Should timing be activated for propagation routines?",
         &(conshdlrdata->enableproptiming), TRUE, DEFAULT_ENABLEPROPTIMING, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/SDP/removesmallval",
         "Should small values in the constraints be removed?",
         &(conshdlrdata->removesmallval), TRUE, DEFAULT_REMOVESMALLVAL, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates the handler for rank 1 SDP constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSdpRank1(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLR* sdpconshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   assert( scip != NULL );

   /* allocate memory for the conshdlrdata */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );

   /* only use one parameter */
   conshdlrdata->diaggezerocuts = FALSE;
   conshdlrdata->propupperbounds = FALSE;
   conshdlrdata->propubpresol = FALSE;
   conshdlrdata->prop3minors = FALSE;
   conshdlrdata->nonconst3minors = FALSE;
   conshdlrdata->prop3mprobing = FALSE;
   conshdlrdata->proptightenbounds = FALSE;
   conshdlrdata->proptbprobing = FALSE;
   conshdlrdata->tightenboundscont = FALSE;
   conshdlrdata->tightenmatrices = FALSE;
   conshdlrdata->tightenbounds = FALSE;
   conshdlrdata->diagzeroimplcuts = FALSE;
   conshdlrdata->twominorlinconss = FALSE;
   conshdlrdata->twominorprodconss = FALSE;
   conshdlrdata->twominorvarbounds = FALSE;
   conshdlrdata->quadconsrank1 = FALSE;
   conshdlrdata->upgradequadconss = FALSE;
   conshdlrdata->upgradekeepquad = FALSE;
   conshdlrdata->separateonecut = FALSE;
   conshdlrdata->cutstopool = FALSE;
   conshdlrdata->sparsifycut = FALSE;
   conshdlrdata->sparsifyfactor = SCIP_INVALID;
   conshdlrdata->sparsifytargetsize = -1;
   conshdlrdata->multiplesparsecuts = FALSE;
   conshdlrdata->maxnsparsecuts = 0;
   conshdlrdata->enforcesdp = FALSE;
   conshdlrdata->onlyfixedintssdp = FALSE;
   conshdlrdata->addsocrelax = FALSE;
   conshdlrdata->maxnvarsquadupgd = 0;
   conshdlrdata->triedlinearconss = FALSE;
   conshdlrdata->triedvarbounds = FALSE;
   conshdlrdata->rank1approxheur = FALSE;
   conshdlrdata->generaterows = FALSE;
   conshdlrdata->generatecmir = FALSE;
#ifdef OMP
   conshdlrdata->nthreads = 0;
#endif
   conshdlrdata->usedimacsfeastol = FALSE;
   conshdlrdata->recomputesparseev = FALSE;
   conshdlrdata->recomputeinitial = FALSE;
   conshdlrdata->exacttrans = FALSE;
   conshdlrdata->presollinconssparam = 0;
   conshdlrdata->additionalstats = FALSE;
   conshdlrdata->enableproptiming = FALSE;

   /* parameters are retrieved through the SDP constraint handler */
   sdpconshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( sdpconshdlr == NULL )
   {
      SCIPerrorMessage("Needs constraint handler <%s> to work.\n", CONSHDLR_NAME);
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata->sdpconshdlrdata = SCIPconshdlrGetData(sdpconshdlr);
   assert( conshdlrdata->sdpconshdlrdata != NULL );

   conshdlrdata->quadconsidx = NULL;
   conshdlrdata->quadconsvars = NULL;
   conshdlrdata->nquadconsidx = 0;
   conshdlrdata->X = NULL;
   conshdlrdata->nsdpvars = 0;
   conshdlrdata->sdpcons = NULL;
   conshdlrdata->randnumgen = NULL;
   conshdlrdata->relaxsdp = NULL;
   conshdlrdata->dimacsfeastol = SCIP_INVALID;
   conshdlrdata->ncallspropub = 0;
   conshdlrdata->ncallsproptb = 0;
   conshdlrdata->ncallsprop3minor = 0;
   conshdlrdata->propubtime = NULL;
   conshdlrdata->proptbtime = NULL;
   conshdlrdata->prop3minortime = NULL;
   conshdlrdata->maxtimepropub = 0.0;
   conshdlrdata->maxtimeproptb = 0.0;
   conshdlrdata->maxtimeprop3minor = 0.0;
   conshdlrdata->npropub = 0;
   conshdlrdata->nproptb = 0;
   conshdlrdata->nprop3minor = 0;
   conshdlrdata->npropcutoffub = 0;
   conshdlrdata->npropcutofftb = 0;
   conshdlrdata->npropcutoff3m = 0;
   conshdlrdata->npropintrndub = 0;
   conshdlrdata->npropintrndtb = 0;
   conshdlrdata->nproppreub = 0;
   conshdlrdata->nproppretb = 0;
   conshdlrdata->nproppre3m = 0;
   conshdlrdata->nproppreintrndub = 0;
   conshdlrdata->nproppreintrndtb = 0;
   conshdlrdata->nproppreintrnd3m = 0;
   conshdlrdata->npropprobub = 0;
   conshdlrdata->npropprobtb = 0;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLRRANK1_NAME, CONSHDLRRANK1_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSdp, consEnfopsSdp, consCheckSdp, consLockSdp, conshdlrdata) );

   assert( conshdlr != NULL );

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSdp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeSdp) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopySdpRank1, consCopySdp) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr,consInitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitSdp) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreSdp) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolSdp) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolSdp, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSdp, CONSHDLR_PROPFREQ, FALSE, CONSHDLR_PROPTIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropSdp) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpSdp, consSepasolSdp, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforelaxSdp) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSdp) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintSdp) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseSdp) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsSdp) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsSdp) );

   /* include upgrading function for quadratic constraints */
#if ( SCIP_VERSION >= 800 || ( SCIP_VERSION < 800 && SCIP_APIVERSION >= 100 ) )
   SCIP_CALL( SCIPincludeConsUpgradeNonlinear(scip, consQuadConsUpgdSdp, 0, TRUE, CONSHDLRRANK1_NAME) );
#else
   SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, consQuadConsUpgdSdp, 0, TRUE, CONSHDLRRANK1_NAME) );
#endif

   return SCIP_OKAY;
}

/** for given row and column (i,j) computes the position in the lower triangular part, if
 *  these positions are numbered from 0 to n(n+1)/2 - 1, this needs to be called for i >= j
 */
int SCIPconsSdpCompLowerTriangPos(
   int                   i,                  /**< row index */
   int                   j                   /**< column index */
   )
{
   assert( j >= 0 );
   assert( i >= j );

   return i*(i+1)/2 + j;
}

/** get the data belonging to a single SDP-constraint
 *
 *  In arraylength the length of the nvarnonz, col, row and val arrays has to be given, if it is not sufficient to store all block-pointers that
 *  need to be inserted, a debug message will be thrown and this variable will be set to the needed length.
 *  constnnonz should give the length of the const arrays, if it is too short it will also give the needed number and a debug message is thrown.
 *  rankone and maxevsubmat can be NULL-pointers, if the corresponding information is not needed.
 */
SCIP_RETCODE SCIPconsSdpGetData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   int*                  nvars,              /**< pointer to store the number of variables in this SDP constraint */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  blocksize,          /**< pointer to store the size of this SDP-block */
   int*                  arraylength,        /**< length of the given nvarnonz, col, row and val arrays; if this is too short, this will return the needed length */
   int*                  nvarnonz,           /**< pointer to store the number of nonzeros for each variable, also length of the arrays col/row/val are
                                              *   pointing to */
   int**                 col,                /**< pointer to store the column indices of the nonzeros for each variable */
   int**                 row,                /**< pointer to store the row indices of the nonzeros for each variable */
   SCIP_Real**           val,                /**< pointer to store the values of the nonzeros for each variable */
   SCIP_VAR**            vars,               /**< pointer to store the SCIP variables present in this constraint that correspond to the indices in col/row/val */
   int*                  constnnonz,         /**< pointer to store the number of nonzeros in the constant part of this SDP constraint, also length of
                                              *   the const arrays */
   int*                  constcol,           /**< pointer to store the column indices of the constant nonzeros */
   int*                  constrow,           /**< pointer to store the row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< pointer to store the values of the constant nonzeros */
   SCIP_Bool*            rankone,            /**< pointer to store if matrix should be rank one (or NULL, if information not necessary) */
   int**                 maxevsubmat,        /**< pointer to store two row indices of 2x2 subdeterminant with maximal eigenvalue [-1,-1 if not yet computed] (or NULL, if information not necessary) */
   SCIP_Bool*            addedquadcons       /**< pointer to store if the quadratic 2x2-minor constraints already added (in the rank1-case) (or NULL, if information not necessary) */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( nvars != NULL );
   assert( nnonz != NULL );
   assert( blocksize != NULL );
   assert( arraylength != NULL );
   assert( nvarnonz != NULL );
   assert( col != NULL );
   assert( row != NULL );
   assert( val != NULL );
   assert( vars != NULL );
   assert( constnnonz != NULL );

   consdata = SCIPconsGetData(cons);

   assert( consdata->constnnonz == 0 || ( constcol != NULL && constrow != NULL && constval != NULL ) );

   *nvars = consdata->nvars;
   *nnonz = consdata->nnonz;
   *blocksize = consdata->blocksize;

   for (i = 0; i < consdata->nvars; i++)
      vars[i] = consdata->vars[i];

   /* check that the sdp-arrays are long enough to store the information */
   if ( *arraylength < consdata->nvars )
   {
      SCIPdebugMsg(scip, "nvarnonz, col, row and val arrays were not long enough to store the information for cons %s, they need to be at least"
         "size %d, given was only length %d!\n", SCIPconsGetName(cons), consdata->nvars, *arraylength);
      *arraylength = consdata->nvars;
   }
   else
   {
      for (i = 0; i < consdata->nvars; i++)
      {
         nvarnonz[i] = consdata->nvarnonz[i];
         /* set the pointers for each variable */
         col[i] = consdata->col[i];
         row[i] = consdata->row[i];
         val[i] = consdata->val[i];
      }
   }

   /* set the constant pointers (if a constant part exists) */
   if ( consdata->constnnonz > 0 )
   {
      if ( consdata->constnnonz > *constnnonz )
      {
         SCIPdebugMsg(scip, "The constant nonzeros arrays were not long enough to store the information for cons %s, they need to be at least"
            "size %d, given was only length %d! \n", SCIPconsGetName(cons), consdata->constnnonz, *constnnonz);
      }
      else
      {
         for (i = 0; i < consdata->constnnonz; i++)
         {
            constcol[i] = consdata->constcol[i];
            constrow[i] = consdata->constrow[i];
            constval[i] = consdata->constval[i];
         }
      }
   }

   *constnnonz = consdata->constnnonz;

   /* set the information about rankone, the current submatrix with largest minimal eigenvalue ([-1,-1] if not yet
      computed), and the quadratic 2x2-minor constraints if desired */
   if ( rankone != NULL && maxevsubmat != NULL )
   {
      *rankone = consdata->rankone;
      *maxevsubmat[0] = consdata->maxevsubmat[0];
      *maxevsubmat[1] = consdata->maxevsubmat[1];
      *addedquadcons = consdata->addedquadcons;
   }

   return SCIP_OKAY;
}

/** gets the number of nonzeros and constant nonzeros for this SDP constraint
 *
 *  Either nnonz or constnnonz may be NULL if only the other one is needed.
 */
SCIP_RETCODE SCIPconsSdpGetNNonz(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get number of nonzeros for */
   int*                  nnonz,              /**< pointer to store the number of nonzeros in this SDP constraint */
   int*                  constnnonz          /**< pointer to store the number of nonzeros in the constant part of this SDP constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   if ( nnonz != NULL )
      *nnonz = consdata->nnonz;

   if ( constnnonz != NULL )
      *constnnonz = consdata->constnnonz;

   return SCIP_OKAY;
}

/** gets the number of variables of the SDP constraint */
int SCIPconsSdpGetNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get number of variables for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->nvars;
}

/** gets the variables of the SDP constraint */
SCIP_VAR** SCIPconsSdpGetVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get variables for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->vars;
}

/** gets the blocksize of the SDP constraint */
int SCIPconsSdpGetBlocksize(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< SDP constraint to get blocksize for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->blocksize;
}

/** gets the full constraint Matrix \f$ A_j \f$ for a given variable j */
SCIP_RETCODE SCIPconsSdpGetFullAj(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   int                   j,                  /**< the variable j to get the corresponding matrix \f$ A_j \f$ for */
   SCIP_Real*            Aj                  /**< pointer to store the full matrix \f$ A_j \f$ */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( j >= 0 );
   assert( Aj != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   blocksize = consdata->blocksize;

   assert( j < consdata->nvars );

   for (i = 0; i < blocksize * blocksize; i++)
      Aj[i] = 0;

   for (i = 0; i < consdata->nvarnonz[j]; i++)
   {
      Aj[consdata->col[j][i] * blocksize + consdata->row[j][i]] = consdata->val[j][i]; /*lint !e679*/
      Aj[consdata->row[j][i] * blocksize + consdata->col[j][i]] = consdata->val[j][i]; /*lint !e679*/
   }

   return SCIP_OKAY;
}

/** gives an n*n-long array with the full constant matrix */
SCIP_RETCODE SCIPconsSdpGetFullConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get matrix for */
   SCIP_Real*            mat                 /**< pointer to store the full constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   blocksize = consdata->blocksize;

   for (i = 0; i < blocksize; i++)
   {
      for (j = 0; j < blocksize; j++)
         mat[i * blocksize + j] = 0.0; /*lint !e679*/
   }

   for (i = 0; i < consdata->constnnonz; i++)
   {
      mat[consdata->constcol[i] * blocksize + consdata->constrow[i]] = consdata->constval[i]; /*lint !e679*/
      mat[consdata->constrow[i] * blocksize + consdata->constcol[i]] = consdata->constval[i]; /*lint !e679*/
   }

   return SCIP_OKAY;
}

/** gives a 0.5*n*(n+1)-long array with the lower triangular part of the constant matrix indexed by SCIPconsSdpCompLowerTriangPos */
SCIP_RETCODE SCIPconsSdpGetLowerTriangConstMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< SDP constraint to get data of */
   SCIP_Real*            mat                 /**< pointer to store the lower triangular part of the constant matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int blocksize;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( mat != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   blocksize = consdata->blocksize;

   /* initialize the matrix with 0 */
   for (i = 0; i < (blocksize * (blocksize + 1)) / 2; i++)
      mat[i] = 0.0;

   for (i = 0; i < consdata->constnnonz; i++)
      mat[SCIPconsSdpCompLowerTriangPos(consdata->constrow[i], consdata->constcol[i])] = consdata->constval[i];

   return SCIP_OKAY;
}

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
SCIP_RETCODE SCIPconsSdpGuessInitialPoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the initial point should be constructed */
   SCIP_Real*            lambdastar          /**< pointer to store the guess for the initial point */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real sparsity;
   SCIP_Real maxinfnorm;
   SCIP_Real maxconst;
   SCIP_Real mininfnorm;
   SCIP_Real maxobj;
   SCIP_Real maxbound;
   SCIP_Real primalguess;
   SCIP_Real dualguess;
   SCIP_Real compval;
   int blocksize;
   int i;
   int v;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( lambdastar != NULL );
   assert( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) == 0 || strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLRRANK1_NAME) == 0 );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* If there are no nonzeros, we cannot use the usual formula, since it divides through the number of nonzeros. In this case,
    * however, we will not solve an SDP anyways but at most an LP (more likely the problem will be solved in local presolving,
    * if all variables are fixed and not only those in the SDP-part), so we just take the default value of SDPA.
    */
   if ( consdata->nnonz == 0 )
   {
      *lambdastar = 100.0;
      return SCIP_OKAY;
   }

   blocksize = consdata->blocksize;

   sparsity = consdata->nnonz / (0.5 * blocksize * (blocksize + 1));

   /* compute the maximum entry of the A_i */
   maxinfnorm = 0.0;
   mininfnorm = SCIPinfinity(scip);
   for (v = 0; v < consdata->nvars; v++)
   {
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         if ( SCIPisGT(scip, REALABS(consdata->val[v][i]), maxinfnorm ) )
            maxinfnorm = REALABS(consdata->val[v][i]);
         if ( SCIPisLT(scip, REALABS(consdata->val[v][i]), mininfnorm) )
            mininfnorm = REALABS(consdata->val[v][i]);
      }
   }
   maxconst = 0.0;
   for (i = 0; i < consdata->constnnonz; i++)
   {
      if ( SCIPisGT(scip, REALABS(consdata->constval[i]), maxconst ) )
         maxconst = REALABS(consdata->constval[i]);
   }

   assert( SCIPisGT(scip, mininfnorm, 0.0) );

   /* compute maximum b_i and bound */
   maxobj = 0.0;
   maxbound = 0.0;
   for (v = 0; v < consdata->nvars; v++)
   {
      if ( SCIPisGT(scip, REALABS(SCIPvarGetObj(consdata->vars[v])), maxobj) )
         maxobj = REALABS(SCIPvarGetObj(consdata->vars[v]));
      compval = SCIPisInfinity(scip, REALABS(SCIPvarGetUbGlobal(consdata->vars[v]))) ? 1e+6 : REALABS(SCIPvarGetUbGlobal(consdata->vars[v]));
      if ( SCIPisGT(scip, compval, maxbound) )
         maxbound = compval;
      compval = SCIPisInfinity(scip, REALABS(SCIPvarGetLbGlobal(consdata->vars[v]))) ? 1e+6 : REALABS(SCIPvarGetUbGlobal(consdata->vars[v]));
      if ( SCIPisGT(scip, compval, maxbound) )
         maxbound = compval;
   }

   /* if all variables were unbounded, we set the value to 10^6 */
   if ( SCIPisEQ(scip, maxbound, 0.0) )
      maxbound = 1E+6;

   /* compute primal and dual guess */
   primalguess = maxobj / (sparsity * mininfnorm);
   dualguess = sparsity * maxinfnorm * maxbound + maxconst;

   if ( SCIPisGT(scip, primalguess, dualguess) )
      *lambdastar = primalguess;
   else
      *lambdastar = dualguess;

   return SCIP_OKAY;
}

/** Gets maximum absolute entry of constant matrix \f$ A_0 \f$ */
SCIP_Real SCIPconsSdpGetMaxConstEntry(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< the constraint to get the maximum constant matrix entry for */
   )
{
   SCIP_CONSDATA* consdata;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);

   return consdata->maxrhsentry;
}

/** Gets maximum absolute entry of all matrices \f$ A_i \f$ */
SCIP_Real SCIPconsSdpGetMaxSdpCoef(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< the constraint to get the maximum constant matrix entry for */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real maxcoef;
   int v;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);

   maxcoef = 0.0;

   for (v = 0; v < consdata->nvars; v++)
   {
      for (i = 0; i < consdata->nvarnonz[v]; i++)
      {
         if ( SCIPisGT(scip, REALABS(consdata->val[v][i]), maxcoef) )
            maxcoef = REALABS(consdata->val[v][i]);
      }
   }

   return maxcoef;
}

/** Computes an upper bound on the number of nonzeros of the (dual) SDP matrix \f$ Z = \sum_{j=1}^n A_j y_j - A_0 \f$,
 *  this should be used to allocate enough memory before calling SCIPconsSdpComputeSparseSdpMatrix.
 *
 *  Upper bound is computed as \f$ \min \{ \sum_{v \leq m} \text{nvarnonz}(v) + \text{constnnonz}, n \cdot (n+1) / 2 \} \f$.
 */
int SCIPconsSdpComputeUbSparseSdpMatrixLength(
   SCIP_CONS*            cons                /**< the constraint for which the matrix should be assembled */
   )
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   int v;
   int ub;
   int denselength;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   ub = consdata->constnnonz;

   for (v = 0; v < consdata->nvars; v++)
      ub += consdata->nvarnonz[v];

   denselength = consdata->blocksize * (consdata->blocksize + 1) / 2;

   return (ub <= denselength ? ub : denselength);
}

/** Computes (dual) SDP matrix \f$ Z = \sum_{j=1}^n A_j y_j - A_0 \f$ and returns it in sparse format
 *  @note row, col and val should have memory allocated equal to SCIPconsSdpComputeUbSparseSdpMatrixLength(),
 *        if the memory is not sufficient, length will be set to -1 and an error will be thrown
 */
SCIP_RETCODE SCIPconsSdpComputeSparseSdpMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< the constraint for which the Matrix should be assembled */
   SCIP_SOL*             sol,                /**< the solution to assemble the matrix for */
   int*                  length,             /**< input: allocated memory for row/col/val arrays
                                              *   output: number of nonzeros of the matrix / length of row/col/val arrays */
   int*                  row,                /**< pointer to store row indices of SDP-matrix */
   int*                  col,                /**< pointer to store column indices of SDP-matrix */
   SCIP_Real*            val                 /**< pointer to store values of SDP-matrix */
   )
{
   SCIP_CONSDATA* consdata;
   int i;
   int v;
   int nnonz;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( sol != NULL );
   assert( length != NULL );
   assert( row != NULL );
   assert( col != NULL );
   assert( val != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   /* initialize nnonz/row/col/val with constant arrays */
   nnonz = consdata->constnnonz;
   if ( *length < nnonz )
   {
      *length = -1;
      SCIPerrorMessage("Arrays not long enough in SCIPconsSdpComputeSparseSdpMatrix, length %d given, need at least %d (probably more)\n",
            *length, nnonz);
      return SCIP_ERROR;
   }

   for (i = 0; i < consdata->constnnonz; i++)
   {
      row[i] = consdata->constrow[i];
      col[i] = consdata->constcol[i];
      val[i] = -1 * consdata->constval[i];
   }

   /* add all variable arrays multiplied by corresponding solution value */
   for (v = 0; v < consdata->nvars; v++)
   {
      SCIP_CALL( SCIPsdpVarfixerMergeArrays(SCIPblkmem(scip), SCIPepsilon(scip), consdata->row[v], consdata->col[v], consdata->val[v], consdata->nvarnonz[v],
            FALSE, SCIPgetSolVal(scip, sol, consdata->vars[v]), row, col, val, &nnonz, *length) );
      if ( nnonz > *length )
      {
         *length = -1;
         SCIPerrorMessage("Arrays not long enough in SCIPconsSdpComputeSparseSdpMatrix, length %d given, need at least %d (probably more)\n",
               *length, nnonz);
         return SCIP_ERROR;
      }
   }

   /* update length pointer */
   *length = nnonz;

   return SCIP_OKAY;
}

/** returns whether the matrix should be rank one */
SCIP_Bool SCIPconsSdpShouldBeRankOne(
   SCIP_CONS*            cons                /**< the constraint for which the existence of a rank one constraint should be checked */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->rankone;
}

/** returns two row indices of 2x2 subdeterminant with maximal eigenvalue [or -1,-1 if not available] */
SCIP_RETCODE SCIPconsSdpGetMaxEVSubmat(
   SCIP_CONS*            cons,               /**< the constraint for which the existence of a rank one constraint should be checked */
   int**                 maxevsubmat         /**< pointer to store the two row indices of 2x2 subdeterminant with
                                              *   maximal eigenvalue [or -1,-1 if not available] */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( maxevsubmat != NULL );

   *maxevsubmat[0] = consdata->maxevsubmat[0];
   *maxevsubmat[1] = consdata->maxevsubmat[1];

   return SCIP_OKAY;
}

/** returns whether the quadratic 2x2-minor constraints are already added (in the rank1-case) */
SCIP_Bool SCIPconsSdpAddedQuadCons(
   SCIP_CONS*            cons                /**< the constraint for which it should be checked whether the quadratic 2x2-minor constraints are already added (in the rank1-case) */
   )
{
   SCIP_CONSDATA* consdata;

   assert( cons != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   return consdata->addedquadcons;
}

/** creates an SDP-constraint
 *
 *  The matrices should be lower triangular.
 */
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
   SCIP_VAR**            vars,               /**< variables present in this SDP constraint that correspond to the indices in col/row/val */
   int                   constnnonz,         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeros */
   int*                  constrow,           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< values of the constant nonzeros */
   SCIP_Bool             removeduplicates    /**< Should duplicate matrix entries be removed (then order of col/row/val might change)? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( name != NULL );
   assert( nvars >= 0 );
   assert( nnonz >= 0 );
   assert( blocksize >= 0 );
   assert( constnnonz >= 0 );
   assert( nvars == 0 || vars != NULL );
   assert( nnonz == 0 || (nvarnonz != NULL && col != NULL && row != NULL && val != NULL ));
   assert( constnnonz == 0 || (constcol != NULL && constrow != NULL && constval != NULL ));

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );

   for (i = 0; i < nvars; i++)
   {
      assert( nvarnonz[i] >= 0 );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i]));
   }

   consdata->nvars = nvars;
   consdata->nnonz = nnonz;
   consdata->constnnonz = constnnonz;
   consdata->blocksize = blocksize;
   consdata->locks = NULL;
   consdata->matrixvar = NULL;
   consdata->matrixval = NULL;
   consdata->matrixconst = NULL;
   consdata->nsingle = 0;
   consdata->propubpossible = TRUE;
   consdata->diagconstantone = FALSE;
   consdata->tracebound = -2.0;
   consdata->allmatricespsd = FALSE;
   consdata->initallmatricespsd = FALSE;
   consdata->issymunique = NULL;

   for (i = 0; i < nvars; i++)
   {
      consdata->nvarnonz[i] = nvarnonz[i];

      if ( nvarnonz[i] > 0 )
      {
         /* if duplicate matrix entries should be removed */
         if ( removeduplicates )
         {
            int cnt = 0;
            int c = 0;

            /* sort by rows, then columns */
            SCIPsdpVarfixerSortRowCol(row[i], col[i], val[i], nvarnonz[i]);

            while ( cnt < nvarnonz[i] )
            {
               assert( 0 <= row[i][cnt] && row[i][cnt] < blocksize );
               assert( 0 <= col[i][cnt] && col[i][cnt] < blocksize );
               assert( row[i][cnt] >= col[i][cnt] );

               while ( cnt < nvarnonz[i] - 1 && row[i][cnt] == row[i][cnt+1] && col[i][cnt] == col[i][cnt+1] )
               {
                  if ( ! SCIPisEQ(scip, val[i][cnt], val[i][cnt+1]) )
                  {
                     SCIPerrorMessage("Duplicate matrix entry (%d,%d) with different value (%g vs. %g).\n", row[i][cnt], col[i][cnt], val[i][cnt], val[i][cnt+1]);
                     return SCIP_INVALIDDATA;
                  }
                  ++cnt;
               }

               if ( ! conshdlrdata->sdpconshdlrdata->removesmallval || ! SCIPisZero(scip, val[i][cnt]) )
               {
                  consdata->row[i][c] = row[i][cnt];
                  consdata->col[i][c] = col[i][cnt];
                  consdata->val[i][c] = val[i][cnt];
                  ++c;
               }
               ++cnt;
            }

            /* possibly correct size; a reallocation should happen rarely */
            if ( c < nvarnonz[i] )
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i], c));
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i], c));
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i], c));
               consdata->nvarnonz[i] = c;
            }
         }
         else
         {
            for (j = 0; j < nvarnonz[i]; j++)
            {
               assert( 0 <= row[i][j] && row[i][j] < blocksize );
               assert( 0 <= col[i][j] && col[i][j] < blocksize );
               assert( row[i][j] >= col[i][j] );

               consdata->row[i][j] = row[i][j];
               consdata->col[i][j] = col[i][j];
               consdata->val[i][j] = val[i][j];
            }
         }
      }
   }

   if ( constnnonz > 0 )
   {
      /* if duplicate matrix entries should be removed */
      if ( removeduplicates )
      {
         int cnt = 0;
         int c = 0;

         /* sort by rows, then columns */
         SCIPsdpVarfixerSortRowCol(constrow, constcol, constval, constnnonz);

         while ( cnt < constnnonz )
         {
            while ( cnt < constnnonz - 1 && constrow[cnt] == constrow[cnt+1] && constcol[cnt] == constcol[cnt+1] )
            {
               if ( ! SCIPisEQ(scip, constval[cnt], constval[cnt+1]) )
               {
                  SCIPerrorMessage("Duplicate entry (%d,%d) with different value (%g vs. %g) in constant matrix.\n", constrow[cnt], constcol[cnt], constval[cnt], constval[cnt+1]);
                  return SCIP_INVALIDDATA;
               }
               ++cnt;
            }

            assert( 0 <= constrow[cnt] && constrow[cnt] < blocksize );
            assert( 0 <= constcol[cnt] && constcol[cnt] < blocksize );
            assert( constrow[cnt] >= constcol[cnt] );

            if ( ! conshdlrdata->sdpconshdlrdata->removesmallval || ! SCIPisZero(scip, constval[cnt]) )
            {
               consdata->constrow[c] = constrow[cnt];
               consdata->constcol[c] = constcol[cnt];
               consdata->constval[c] = constval[cnt];
               ++c;
            }
            ++cnt;
         }

         /* possibly correct size; a reallocation should happen rarely */
         if ( c < constnnonz )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, constnnonz, c) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, constnnonz, c) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, constnnonz, c) );
            consdata->constnnonz = c;
         }
      }
      else
      {
         for (j = 0; j < constnnonz; j++)
         {
            assert( 0 <= constrow[j] && constrow[j] < blocksize );
            assert( 0 <= constcol[j] && constcol[j] < blocksize );
            assert( constrow[j] >= constcol[j] );

            consdata->constrow[j] = constrow[j];
            consdata->constcol[j] = constcol[j];
            consdata->constval[j] = constval[j];
         }
      }
   }

   for (i = 0; i < nvars; i++)
   {
      consdata->vars[i] = vars[i];
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }
   SCIPdebugMsg(scip, "creating cons %s.\n", name);

   consdata->rankone = FALSE;

   /* allocate memory for rank one constraint */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->maxevsubmat, 2) );
   consdata->maxevsubmat[0] = -1;
   consdata->maxevsubmat[1] = -1;

   /* quadratic 2x2-minor constraints added? */
   consdata->addedquadcons = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

   return SCIP_OKAY;
}


/** creates a rank 1 SDP-constraint
 *
 *  The matrices should be lower triangular.
 */
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
   SCIP_VAR**            vars,               /**< variables present in this SDP constraint that correspond to the indices in col/row/val */
   int                   constnnonz,         /**< number of nonzeros in the constant part of this SDP constraint */
   int*                  constcol,           /**< column indices of the constant nonzeros */
   int*                  constrow,           /**< row indices of the constant nonzeros */
   SCIP_Real*            constval,           /**< values of the constant nonzeros */
   SCIP_Bool             removeduplicates    /**< Should duplicate matrix entries be removed (then order of col/row/val might change)? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( name != NULL );
   assert( nvars >= 0 );
   assert( nnonz >= 0 );
   assert( blocksize >= 0 );
   assert( constnnonz >= 0 );
   assert( nvars == 0 || vars != NULL );
   assert( nnonz == 0 || (nvarnonz != NULL && col != NULL && row != NULL && val != NULL ));
   assert( constnnonz == 0 || (constcol != NULL && constrow != NULL && constval != NULL ));

   conshdlr = SCIPfindConshdlr(scip, CONSHDLRRANK1_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("Rank 1 SDP constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   /* create constraint data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &consdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->nvarnonz, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constcol, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constrow, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->constval, constnnonz) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vars, nvars) );

   for (i = 0; i < nvars; i++)
   {
      assert( nvarnonz[i] >= 0 );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i]));
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i]));
   }

   consdata->nvars = nvars;
   consdata->nnonz = nnonz;
   consdata->constnnonz = constnnonz;
   consdata->blocksize = blocksize;
   consdata->locks = NULL;
   consdata->matrixvar = NULL;
   consdata->matrixval = NULL;
   consdata->matrixconst = NULL;
   consdata->nsingle = 0;
   consdata->propubpossible = TRUE;
   consdata->diagconstantone = FALSE;
   consdata->tracebound = -2.0;
   consdata->allmatricespsd = FALSE;
   consdata->initallmatricespsd = FALSE;
   consdata->issymunique = NULL;

   for (i = 0; i < nvars; i++)
   {
      consdata->nvarnonz[i] = nvarnonz[i];

      if ( nvarnonz[i] > 0 )
      {
         /* if duplicate matrix entries should be removed */
         if ( removeduplicates )
         {
            int cnt = 0;
            int c = 0;

            /* sort by rows, then columns */
            SCIPsdpVarfixerSortRowCol(row[i], col[i], val[i], nvarnonz[i]);

            while ( cnt < nvarnonz[i] )
            {
               assert( 0 <= row[i][cnt] && row[i][cnt] < blocksize );
               assert( 0 <= col[i][cnt] && col[i][cnt] < blocksize );
               assert( row[i][cnt] >= col[i][cnt] );

               while ( cnt < nvarnonz[i] - 1 && row[i][cnt] == row[i][cnt+1] && col[i][cnt] == col[i][cnt+1] )
               {
                  if ( ! SCIPisEQ(scip, val[i][cnt], val[i][cnt+1]) )
                  {
                     SCIPerrorMessage("Duplicate matrix entry (%d,%d) with different value (%g vs. %g).\n", row[i][cnt], col[i][cnt], val[i][cnt], val[i][cnt+1]);
                     return SCIP_INVALIDDATA;
                  }
                  ++cnt;
               }

               if ( ! conshdlrdata->sdpconshdlrdata->removesmallval || ! SCIPisZero(scip, val[i][cnt]) )
               {
                  consdata->row[i][c] = row[i][cnt];
                  consdata->col[i][c] = col[i][cnt];
                  consdata->val[i][c] = val[i][cnt];
                  ++c;
               }
               ++cnt;
            }

            /* possibly correct size; a reallocation should happen rarely */
            if ( c < nvarnonz[i] )
            {
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->col[i], nvarnonz[i], c));
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->row[i], nvarnonz[i], c));
               SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->val[i], nvarnonz[i], c));
               consdata->nvarnonz[i] = c;
            }
         }
         else
         {
            for (j = 0; j < nvarnonz[i]; j++)
            {
               assert( 0 <= row[i][j] && row[i][j] < blocksize );
               assert( 0 <= col[i][j] && col[i][j] < blocksize );
               assert( row[i][j] >= col[i][j] );

               consdata->row[i][j] = row[i][j];
               consdata->col[i][j] = col[i][j];
               consdata->val[i][j] = val[i][j];
            }
         }
      }
   }

   if ( constnnonz > 0 )
   {
      /* if duplicate matrix entries should be removed */
      if ( removeduplicates )
      {
         int cnt = 0;
         int c = 0;

         /* sort by rows, then columns */
         SCIPsdpVarfixerSortRowCol(constrow, constcol, constval, constnnonz);

         while ( cnt < constnnonz )
         {
            while ( cnt < constnnonz - 1 && constrow[cnt] == constrow[cnt+1] && constcol[cnt] == constcol[cnt+1] )
            {
               if ( ! SCIPisEQ(scip, constval[cnt], constval[cnt+1]) )
               {
                  SCIPerrorMessage("Duplicate entry (%d,%d) with different value (%g vs. %g) in constant matrix.\n", constrow[cnt], constcol[cnt], constval[cnt], constval[cnt+1]);
                  return SCIP_INVALIDDATA;
               }
               ++cnt;
            }

            assert( 0 <= constrow[cnt] && constrow[cnt] < blocksize );
            assert( 0 <= constcol[cnt] && constcol[cnt] < blocksize );
            assert( constrow[cnt] >= constcol[cnt] );

            if ( ! conshdlrdata->sdpconshdlrdata->removesmallval || ! SCIPisZero(scip, constval[cnt]) )
            {
               consdata->constrow[c] = constrow[cnt];
               consdata->constcol[c] = constcol[cnt];
               consdata->constval[c] = constval[cnt];
               ++c;
            }
            ++cnt;
         }

         /* possibly correct size; a reallocation should happen rarely */
         if ( c < constnnonz )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constcol, constnnonz, c) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constrow, constnnonz, c) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->constval, constnnonz, c) );
            consdata->constnnonz = c;
         }
      }
      else
      {
         for (j = 0; j < constnnonz; j++)
         {
            assert( 0 <= constrow[j] && constrow[j] < blocksize );
            assert( 0 <= constcol[j] && constcol[j] < blocksize );
            assert( constrow[j] >= constcol[j] );

            consdata->constrow[j] = constrow[j];
            consdata->constcol[j] = constcol[j];
            consdata->constval[j] = constval[j];
         }
      }
   }

   for (i = 0; i < nvars; i++)
   {
      consdata->vars[i] = vars[i];
      SCIP_CALL( SCIPcaptureVar(scip, consdata->vars[i]) );
   }
   SCIPdebugMsg(scip, "creating cons %s (rank 1).\n", name);
   consdata->rankone = TRUE;

   /* allocate memory for rank one constraint */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->maxevsubmat, 2) );
   consdata->maxevsubmat[0] = -1;
   consdata->maxevsubmat[1] = -1;

   /* quadratic 2x2-minor constraints added? */
   consdata->addedquadcons = FALSE;

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   /* compute maximum rhs entry for later use in the DIMACS Error Norm */
   SCIP_CALL( setMaxRhsEntry(*cons) );

   return SCIP_OKAY;
}
