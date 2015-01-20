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

/**@file   sdpi_msk.c
 * @ingroup SDPIS
 * @brief  SDP interface for MOSEK
 * @author Bo Jensen
 * @author Leif Naundorf
 */

#include <assert.h>

#define MSKCONST const
#include "mosek.h"

#include "sdpi/sdpi.h"
#include "scip/bitencode.h"
#include <string.h>

/* check version */
#if (MSK_VERSION_MAJOR < 7)
#error "MOSEK version 7 or higher is required for the SDP functionality"
#endif

/* do defines for windows directly her to make the sdpi more independent*/
#if defined(_WIN32) || defined(_WIN64)
#define snprintf _snprintf
#endif

#define scipmskobjsen MSKobjsensee
#define SENSE2MOSEK(objsen) (((objsen)==SCIP_OBJSEN_MINIMIZE)?(MSK_OBJECTIVE_SENSE_MINIMIZE):(MSK_OBJECTIVE_SENSE_MAXIMIZE))

#define MOSEK_CALL(x)  do                                                                                     \
                       {  /*lint --e{641}*/                                                                   \
                          MSKrescodee _restat_;                                                               \
                          _restat_ = (x);                                                                     \
                          if( (_restat_) != MSK_RES_OK && (_restat_ ) != MSK_RES_TRM_MAX_NUM_SETBACKS )       \
                          {                                                                                   \
                             SCIPerrorMessage("SDP Error: MOSEK returned %d\n", (int)_restat_);                \
                             return SCIP_ERROR;                                                             \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )

/* this macro is only called in functions returning SCIP_Bool; thus, we return FALSE if there is an error in optimized mode */
#define ABORT_FALSE(x) { int _restat_;                                  \
      if( (_restat_ = (x)) != 0 )                                       \
      {                                                                 \
         SCIPerrorMessage("SDP Error: MOSEK returned %d\n", (int)_restat_); \
         SCIPABORT();                                                   \
         return FALSE;                                                  \
      }                                                                 \
   }


#define IS_POSINF(x) ((x) >= SCIP_DEFAULT_INFINITY)
#define IS_NEGINF(x) ((x) <= -SCIP_DEFAULT_INFINITY)

static MSKenv_t MosekEnv =           NULL;
static int numsdp         =           0;

static int optimizecount            =  0;
static int nextsdpid                =  1;
static int numstrongbranchmaxiterup =  0;
static int numstrongbranchmaxiterdo =  0;
static int numprimalmaxiter         =  0;
static int numdualmaxiter           =  0;
static int numstrongbranchobjup     =  0;
static int numstrongbranchobjdo     =  0;
static int numprimalobj             =  0;
static int numdualobj               =  0;

#define DEBUG_PARAM_SETTING          0
#define DEBUG_PRINT_STAT             0
#define DEBUG_CHECK_DATA             0
#define DEBUG_EASY_REPRODUCE         0
#define DEBUG_DO_INTPNT_FEAS_CHECK   0
#define DEBUG_CHECK_STATE_TOL        1e-5
#define SHOW_ERRORS                  0
#define ASSERT_ON_NUMERICAL_TROUBLES 0
#define ASSERT_ON_WARNING            0
#define FORCE_MOSEK_LOG              0
#define FORCE_MOSEK_SUMMARY          0
#define FORCE_NO_MAXITER             0
#define FORCE_SILENCE                1
#define SETBACK_LIMIT                250
#define SCIP_CONTROLS_PRICING        1
#define SCIP_CONTROLS_TOLERANCES     1
#define STRONGBRANCH_PRICING         MSK_SIM_SELECTION_SE
#define SUPRESS_NAME_ERROR           1
#define WRITE_DUAL                   0
#define WRITE_PRIMAL                 0
#define WRITE_INTPNT                 0
#define WRITE_ABOVE                  0
#define DEGEN_LEVEL                  MSK_SIM_DEGEN_FREE
#define ALWAYS_SOLVE_PRIMAL          1

/** gives problem and solution status for a Mosek Task
 *
 * With Mosek 7.0, the routine MSK_getsolutionstatus was replaced by
 * MSK_getprosta and MSK_getsolsta.
 */
static
MSKrescodee MSK_getsolutionstatus(
   MSKtask_t             task,               /**< Mosek Task */
   MSKsoltypee           whichsol,           /**< for which type of solution a status is requested */
   MSKprostae*           prosta,             /**< buffer to store problem status, or NULL if not needed */
   MSKsolstae*           solsta              /**< buffer to store solution status, or NULL if not needed */
   )
{
   if( prosta != NULL )
   {
      MOSEK_CALL( MSK_getprosta(task, whichsol, prosta) );
   }
   if( solsta != NULL )
   {
      MOSEK_CALL( MSK_getsolsta(task, whichsol, solsta) );
   }

   return MSK_RES_OK;
}

/**********************************************/

struct SCIP_SDPi
{
   MSKtask_t             task;
   MSKrescodee           termcode;
   int                   itercount;
   SCIP_PRICING          pricing;            /**< SCIP pricing setting  */
   int                   sdpid;
   int                   skxsize;
   int                   skcsize;
   MSKstakeye*           skx;
   MSKstakeye*           skc;
   SCIP_MESSAGEHDLR*     messagehdlr;        /**< messagehdlr handler to printing messages, or NULL */
};

typedef SCIP_DUALPACKET COLPACKET;           /* each column needs two bits of information (basic/on_lower/on_upper) */
#define COLS_PER_PACKET SCIP_DUALPACKETSIZE
typedef SCIP_DUALPACKET ROWPACKET;           /* each row needs two bit of information (basic/on_lower/on_upper) */
#define ROWS_PER_PACKET SCIP_DUALPACKETSIZE

struct SCIP_SDPiState
{
   int                   num;
   MSKsolstae            solsta;
   int                   ncols;
   int                   nrows;
   COLPACKET*            skx;
   ROWPACKET*            skc;
};

/** returns the number of packets needed to store column packet information */
static
int colpacketNum(
   int                   ncols               /**< number of columns to store */
   )
{
   return (ncols+(int)COLS_PER_PACKET-1)/(int)COLS_PER_PACKET;
}

/** returns the number of packets needed to store row packet information */
static
int rowpacketNum(
   int                   nrows               /**< number of rows to store */
   )
{
   return (nrows+(int)ROWS_PER_PACKET-1)/(int)ROWS_PER_PACKET;
}

/** create error string */
static
void MSKAPI printstr(
   void*                 handle,             /**< error handle */
   const char*           str                 /**< string that contains string on output */
   )
{  /*lint --e{715}*/
#if SUPRESS_NAME_ERROR && !FORCE_SILENCE
   char errstr[32];
   snprintf(errstr,32,"MOSEK Error %d",MSK_RES_ERR_DUP_NAME);
   if (0 == strncmp(errstr,str,strlen(errstr)))
      return;
#endif

   SCIPdebugMessage("MOSEK: %s",str);
}

#if DEBUG_CHECK_DATA > 0
/** check data */
static SCIP_RETCODE scip_checkdata(
   SCIP_SDPI*            sdpi,               /**< pointer to an SDP interface structure */
   const char*           functionname        /**< function name */
   )
{
   int i;
   int numcon;
   int numvar;
   int gotbasicsol;
   MSKboundkeye* tbkc;
   MSKboundkeye* tbkx;
   MSKstakeye *tskc;
   MSKstakeye* tskx;
   double* tblc;
   double* tbuc;
   double* tblx;
   double* tbux;

   MOSEK_CALL( MSK_solutiondef(sdpi->task, MSK_SOL_BAS, &gotbasicsol) );

   MOSEK_CALL( MSK_getnumvar(sdpi->task,&numvar) );
   MOSEK_CALL( MSK_getnumcon(sdpi->task,&numcon) );

   /* allocate memory */
   SCIP_ALLOC( BMSallocMemoryArray( &tbkc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tskc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tblc, numcon) );
   SCIP_ALLOC( BMSallocMemoryArray( &tbuc, numcon) );

   SCIP_ALLOC( BMSallocMemoryArray( &tbkx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tskx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tblx, numvar) );
   SCIP_ALLOC( BMSallocMemoryArray( &tbux, numvar) );

   /* Check bounds */
   if( gotbasicsol )
   {
      MOSEK_CALL( MSK_getsolution(sdpi->task, MSK_SOL_BAS, NULL, NULL, tskc, tskx,
            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );
   }

   for( i = 0; i < numvar; i++ )
   {
      MOSEK_CALL( MSK_getbound(sdpi->task,MSK_ACC_VAR,i,&tbkx[i],&tblx[i],&tbux[i]) );
   }

   for( i = 0; i < numcon; i++ )
   {
      MOSEK_CALL( MSK_getbound(sdpi->task,MSK_ACC_CON,i,&tbkc[i],&tblc[i],&tbuc[i]) );
   }

   for( i = 0; i < numcon; ++i )
   {
      if( gotbasicsol )
      {
         if( ( tskc[i] == MSK_SK_FIX && tbkc[i] != MSK_BK_FX ) ||
            ( tskc[i] == MSK_SK_LOW && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_RA ) ) ||
            ( tskc[i] == MSK_SK_UPR && !(tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_RA ) ) )
         {
            SCIPerrorMessage("STATUS KEY ERROR i %d bkc %d skc %d %s\n", i, tbkc[i], tskc[i], functionname);
         }
      }

      if( tbkc[i] == MSK_BK_LO || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA )
      {
         if( isnan(tblc[i]) )
         {
            SCIPdebugMessage("nan in blc : %s\n", functionname);
         }
      }

      if( tbkc[i] == MSK_BK_UP || tbkc[i] == MSK_BK_FX || tbkc[i] == MSK_BK_RA )
      {
         if( isnan(tbuc[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
         }
      }
   }

   for( i = 0; i < numvar; ++i )
   {
      if( tbkx[i] == MSK_BK_LO || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA )
      {
         if( isnan(tblx[i]) )
         {
            SCIPdebugMessage("nan in blx : %s\n",functionname);
         }
      }

      if( tbkx[i] == MSK_BK_UP || tbkx[i] == MSK_BK_FX || tbkx[i] == MSK_BK_RA )
      {
         if( isnan(tbux[i]) )
         {
            SCIPdebugMessage("nan in bux : %s\n", functionname);
            getchar();
         }
      }
   }

   BMSfreeMemoryArray(&tbkc);
   BMSfreeMemoryArray(&tskc);
   BMSfreeMemoryArray(&tblc);
   BMSfreeMemoryArray(&tbuc);
   BMSfreeMemoryArray(&tbkx);
   BMSfreeMemoryArray(&tskx);
   BMSfreeMemoryArray(&tblx);
   BMSfreeMemoryArray(&tbux);

   return SCIP_OKAY;
}
#endif


/*
 * Local functions
 */

static
void generateMskBounds(
   int                   n,
   const double*         lb,
   const double*         ub,
   MSKboundkeye*         bk,
   double*               msklb,
   double*               mskub
   )
{
   int i;

   assert(lb != NULL);
   assert(ub != NULL);
   assert(bk != NULL);
   assert(msklb != NULL);
   assert(mskub != NULL);

   for( i = 0; i < n; i++ )
   {
      msklb[i] = lb[i];
      mskub[i] = ub[i];
      if (IS_NEGINF(lb[i]))
      {
         msklb[i] = -MSK_INFINITY;
         if (IS_POSINF(ub[i]))
         {
            mskub[i] = MSK_INFINITY;
            bk[i] = MSK_BK_FR;
         }
         else
         {
            assert(!IS_NEGINF(ub[i]));
            bk[i] = MSK_BK_UP;
         }
      }
      else
      {
         assert(!IS_POSINF(lb[i]));
         if (IS_POSINF(ub[i]))
         {
            mskub[i] = MSK_INFINITY;
            bk[i] = MSK_BK_LO;
         }
         else if (lb[i] == ub[i])  /**@todo is this good idea to compare the bound without any epsilontic? */
         {
            assert(lb[i]-ub[i]==0);
            assert(ub[i]-lb[i]==0);
            bk[i] = MSK_BK_FX;
         }
         else
         {
            assert(lb[i] < ub[i]);
            bk[i] = MSK_BK_RA;
         }
      }
   }
}

/** get end pointers of arrays */
static
SCIP_RETCODE getEndptrs(
   int                   n,                  /**< array size */
   const int*            beg,                /**< array of beginning indices */
   int                   nnonz,              /**< number of nonzeros */
   int**                 aptre               /**< pointer to store the result */
   )
{
   int i;

   assert(beg != NULL || nnonz == 0);

   SCIP_ALLOC( BMSallocMemoryArray( aptre, n) );

   /*    if (aptre == NULL)
         return NULL;
   */

   if (nnonz > 0)
   {
      assert(beg != NULL);
      for(i = 0; i < n-1; i++)
      {
         (*aptre)[i] = beg[i+1];
         assert((*aptre)[i] >= beg[i]);
      }

      (*aptre)[n-1] = nnonz;
      assert((*aptre)[n-1] >= beg[n-1]);
   }
   else
   {
      for( i = 0; i < n; i++ )
         (*aptre)[i] = 0;
   }

   return SCIP_OKAY;
}

/** compute indices from range */
static
SCIP_RETCODE getIndicesRange(
   int                   first,              /**< first index */
   int                   last,               /**< last index */
   int**                 sub                 /**< pointer to store the indices ranges */
   )
{
   int i;

   assert(first <= last);

   SCIP_ALLOC( BMSallocMemoryArray( sub, (last-first+1)) );

   for( i = first; i <= last; i++ )
   {
      (*sub)[i-first] = i;
   }

   return SCIP_OKAY;
}

/** compute indices from dense array */
static
SCIP_RETCODE getIndicesFromDense(
   int*                  dstat,              /**< array */
   int                   n,                  /**< size of array */
   int*                  count,              /**< array of counts (sizes) */
   int**                 sub                 /**< pointer to store array of indices */
   )
{
   int i;
   int j;

   assert(dstat != NULL);

   *count = 0;
   for( i = 0; i < n; i++ )
   {
      if (dstat[i] == 1)
      {
         (*count)++;
      }
   }

   if( (*count) > 0 )
   {
      SCIP_ALLOC( BMSallocMemoryArray( sub, (*count)) );
   }
   else
      return SCIP_OKAY;

   j = 0;
   for( i = 0; i < n; i++ )
   {
      if (dstat[i] == 1)
      {
         (*sub)[j++] = i;
      }
   }

   return SCIP_OKAY;
}

static
void scale_vec(
   int                   len,
   double*               vec,
   double                s
   )
{
   int i;
   for( i = 0; i < len; i++ )
   {
      vec[i] *= s;
   }
}

static
void scale_bound(
   MSKboundkeye*         bk,
   double*               bl,
   double*               bu,
   double                s
   )
{
   switch(*bk)
   {
   case MSK_BK_LO:
      *bl *= s;
      if (s < 0) *bk = MSK_BK_UP;
      break;
   case MSK_BK_UP:
      *bu *= s;
      if (s < 0) *bk = MSK_BK_LO;
      break;
   case MSK_BK_FX:
   case MSK_BK_RA:
      *bl *= s;
      *bu *= s;
      break;
   case MSK_BK_FR:
      break;
   default:
      assert(FALSE);
      break;
   }  /*lint !e788*/

   if (s < 0)
   {
      double tmp;
      tmp = *bl;
      *bl = *bu;
      *bu = tmp;
   }
}

static
SCIP_RETCODE ensureStateMem(
   SCIP_SDPI*            sdpi,
   int                   ncols,
   int                   nrows
   )
{
   if (sdpi->skxsize < ncols)
   {
      int newsize;
      newsize = MAX(2*sdpi->skxsize, ncols);

      SCIP_ALLOC( BMSreallocMemoryArray( &(sdpi->skx), newsize) );
      sdpi->skxsize = newsize;
   }

   if (sdpi->skcsize < nrows)
   {
      int newsize;
      newsize = MAX(2*sdpi->skcsize, nrows);

      SCIP_ALLOC( BMSreallocMemoryArray( &(sdpi->skc), newsize) );
      sdpi->skcsize = newsize;
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE getbase(
   SCIP_SDPI*            sdpi,
   int                   ncols,
   int                   nrows
   )
{
   SCIPdebugMessage("Calling getbase (%d)\n",sdpi->sdpid);

   SCIP_CALL( ensureStateMem(sdpi,ncols,nrows) );
   MOSEK_CALL( MSK_getsolution(sdpi->task, MSK_SOL_BAS, NULL, NULL, sdpi->skc, sdpi->skx,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE setbase(
   SCIP_SDPI*            sdpi                /**< pointer to an SDP interface structure */
   )
{
   SCIPdebugMessage("Calling setbase (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_putsolution(sdpi->task, MSK_SOL_BAS, sdpi->skc, sdpi->skx, NULL, NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}



/*
 * Miscellaneous Methods
 */

static char mskname[100];

/**@name Miscellaneous Methods */
/**@{ */

/** gets name and version of SDP solver */
const char* SCIPsdpiGetSolverName(
   void
   )
{
   sprintf(mskname, "MOSEK %.2f", (SCIP_Real)MSK_VERSION_MAJOR);
   return mskname;
}

/** gets description of SDP solver (developer, webpage, ...) */
const char* SCIPsdpiGetSolverDesc(
   void
   )
{
   return "Semi-definite Programming Solver developed by MOSEK Optimization Software (www.mosek.com)";
}

/** gets pointer for SDP solver - use only with great care */
void* SCIPsdpiGetSolverPointer(
   SCIP_SDPI*            sdpi                /**< pointer to an SDP interface structure */
   )
{
   return (void*) sdpi->task;
}


/*
 * SDPI Creation and Destruction Methods
 */

/**@name SDPI Creation and Destruction Methods */
/**@{ */

/** creates an SDP problem object */
SCIP_RETCODE SCIPsdpiCreate(
   SCIP_SDPI**           sdpi,               /**< pointer to an SDP interface structure */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler to use for printing messages, or NULL */
   const char*           name,               /**< problem name */
   SCIP_OBJSEN           objsen              /**< objective sense */
   )
{
   assert(sdpi != NULL);
   assert(numsdp >= 0);

   SCIPdebugMessage("Calling SCIPsdpiCreate\n");

   if (!MosekEnv)
   {
      MOSEK_CALL( MSK_makeenv(&MosekEnv, NULL) );
      MOSEK_CALL( MSK_linkfunctoenvstream(MosekEnv, MSK_STREAM_LOG, NULL, printstr) );
      MOSEK_CALL( MSK_initenv(MosekEnv) );
   }

   numsdp++;

   SCIP_ALLOC( BMSallocMemory(sdpi) );

   MOSEK_CALL( MSK_makeemptytask(MosekEnv, &((*sdpi)->task)) );

   MOSEK_CALL( MSK_linkfunctotaskstream((*sdpi)->task, MSK_STREAM_LOG, NULL, printstr) );

   MOSEK_CALL( MSK_putobjsense((*sdpi)->task, SENSE2MOSEK(objsen)) );
   /*
   MOSEK_CALL( MSK_putintparam((*sdpi)->task, MSK_IPAR_SIM_MAX_NUM_SETBACKS, SETBACK_LIMIT) );
   MOSEK_CALL( MSK_putintparam((*sdpi)->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_FREE_SIMPLEX) );
   MOSEK_CALL( MSK_putintparam((*sdpi)->task, MSK_IPAR_SIM_DEGEN, DEGEN_LEVEL) );
   MOSEK_CALL( MSK_putintparam((*sdpi)->task, MSK_IPAR_SIM_SWITCH_OPTIMIZER, MSK_ON) ); */
   /* We only have status keys (recalculate dual solution without dual superbasics) */
   /* MOSEK_CALL( MSK_putintparam((*sdpi)->task, MSK_IPAR_SIM_HOTSTART, MSK_SIM_HOTSTART_STATUS_KEYS) ); */
   MOSEK_CALL( MSK_puttaskname((*sdpi)->task, (char*) name) );

   (*sdpi)->termcode = MSK_RES_OK;
   (*sdpi)->itercount = 0;
   (*sdpi)->pricing = SCIP_PRICING_LPIDEFAULT;
   (*sdpi)->sdpid = nextsdpid++;
   (*sdpi)->skxsize = 0;
   (*sdpi)->skcsize = 0;
   (*sdpi)->skx = NULL;
   (*sdpi)->skc = NULL;
   (*sdpi)->messagehdlr = messagehdlr;

   return SCIP_OKAY;
}

/** deletes an SDP problem object */
SCIP_RETCODE SCIPsdpiFree(
   SCIP_SDPI**           sdpi                /**< pointer to an SDP interface structure */
   )
{
   assert(sdpi != NULL);
   assert(*sdpi != NULL);
   assert(numsdp > 0);

   SCIPdebugMessage("Calling SCIPsdpiFree (%d)\n",(*sdpi)->sdpid);

   MOSEK_CALL( MSK_deletetask(&(*sdpi)->task) );

   BMSfreeMemoryArrayNull(&(*sdpi)->skx);
   BMSfreeMemoryArrayNull(&(*sdpi)->skc);
   BMSfreeMemory(sdpi);

   numsdp--;
   if (numsdp == 0)
   {
      MOSEK_CALL( MSK_deleteenv(&MosekEnv) );
      MosekEnv = NULL;
   }

   return SCIP_OKAY;
}

/*
 * Modification Methods
 */


/** copies SDP data with column matrix into SDP solver */
SCIP_RETCODE SCIPsdpiLoadColSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN           objsen,             /**< objective sense */
   int                   ncols,              /**< number of columns */
   const SCIP_Real*      obj,                /**< objective function values of columns */
   const SCIP_Real*      lb,                 /**< lower bounds of columns */
   const SCIP_Real*      ub,                 /**< upper bounds of columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nrows,              /**< number of rows */
   const SCIP_Real*      lhs,                /**< left hand sides of rows */
   const SCIP_Real*      rhs,                /**< right hand sides of rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements in the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array */
   const int*            ind,                /**< row indices of constraint matrix entries */
   const SCIP_Real*      val                 /**< values of constraint matrix entries */
   )
{  /*lint --e{715}*/
   int* aptre;
   MSKboundkeye* bkc;
   MSKboundkeye* bkx;
   double* blc;
   double* buc;
   double* blx;
   double* bux;

   SCIPdebugMessage("Calling SCIPsdpiLoadColSDP (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   /* initialize all array with NULL */
   aptre = NULL;
   bkc = NULL;
   bkx = NULL;
   blc = NULL;
   buc = NULL;
   blx = NULL;
   bux = NULL;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiLoadColSDP") );
#endif

   if (nrows > 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &bkc, nrows) );
      SCIP_ALLOC( BMSallocMemoryArray( &blc, nrows) );
      SCIP_ALLOC( BMSallocMemoryArray( &buc, nrows) );

      generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);
   }

   if (ncols > 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray( &bkx, ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &blx, ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &bux, ncols) );

      generateMskBounds(ncols, lb, ub, bkx, blx, bux);

      SCIP_CALL( getEndptrs(ncols, beg, nnonz, &aptre) );
   }

   MOSEK_CALL( MSK_inputdata(sdpi->task, nrows, ncols, nrows, ncols, obj, 0.0, beg, aptre, ind, val,
         bkc, blc, buc, bkx, blx, bux) );

   MOSEK_CALL( MSK_putobjsense(sdpi->task, SENSE2MOSEK(objsen)) );


   if( ncols > 0 )
   {
      BMSfreeMemoryArray(&aptre);
      BMSfreeMemoryArray(&bux);
      BMSfreeMemoryArray(&blx);
      BMSfreeMemoryArray(&bkx);
   }

   if( nrows > 0 )
   {
      BMSfreeMemoryArray(&buc);
      BMSfreeMemoryArray(&blc);
      BMSfreeMemoryArray(&bkc);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiLoadColSDP") );
#endif

   return SCIP_OKAY;
}

/** adds a number of positive-semidefinite variables X_j to the SDP */
SCIP_RETCODE SCIPsdpiAddSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nvars,              /**< number of semidefinite variables to be added */
   const int*            dims                /**< dimensions of the semidefinite variables to be added */
   )
{
   assert(sdpi != NULL);
   MOSEK_CALL( MSK_appendbarvars(sdpi->task, nvars, dims) );

   return SCIP_OKAY;
}

/** deletes a number of positive-semidefinite variables X_j from the SDP */
SCIP_RETCODE SCIPsdpiDelSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of SDP variables
                                              *   input:  1 if SDP variable should be deleted, 0 if not
                                              *   output: new indices of SDP variable, -1 if variable was deleted */
   )
{
   int* sub;
   int count;
   int nbarvars;
   int barvar;
   int i;

   assert(sdpi != NULL);
   assert(dstat != NULL);
   MOSEK_CALL( MSK_getnumbarvar(sdpi->task, &nbarvars) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, nbarvars, &count, &sub) );

   barvar = 0;
   for( i = 0; i < nbarvars; i++)
   {
      if (dstat[i] == 1)
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = barvar;
         barvar++;
      }
   }
   SCIPdebugMessage("Deleting %d SDP vars %d, ...\n", count, sub[0]);
   MOSEK_CALL( MSK_removebarvars(sdpi->task, count, sub) );
   BMSfreeMemoryArray(&sub);

   return SCIP_OKAY;
}

/** adds a term of the form <C,X_j> to the objective funtion */
SCIP_RETCODE SCIPsdpiAddSDPTermObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix C */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of C */
   const int*            subi,               /**< the row-indices of the non-zero entries of C */
   const int*            subj,               /**< the column-indices of the non-zero entries of C */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   ind                 /**< index of the symmetric variable X */
   )
{
   SCIP_Longint idx;
   SCIP_Real falpha = 1.0;

   MOSEK_CALL( MSK_appendsparsesymmat(sdpi->task, dim, nnonz, subi, subj, valij, &idx) );
   MOSEK_CALL( MSK_putbarcj(sdpi->task, ind, 1, &idx, &falpha) );

   return SCIP_OKAY;
}


/** adds a term of the form <A,X_j> to a row of the SDP
 *
 * @note as A is symmetric, only the lower triangular part of it must be specified
 */
SCIP_RETCODE SCIPsdpiAddSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   dim,                /**< dimension of the matrix A */
   int                   nnonz,              /**< number of non-zeroes in the lower triagonal part of A */
   const int*            subi,               /**< the row-indices of the non-zero entries of A */
   const int*            subj,               /**< the column-indices of the non-zero entries of A */
   const SCIP_Real*      valij,              /**< the values at the indices specified by the subi and subj arrays */
   int                   row,                /**< row of the SDP where this term should be added */
   int                   ind                 /**< index of the symmetric variable X */
   )
{
   SCIP_Longint idx;
   SCIP_Real falpha = 1.0;
   int i;

   MOSEK_CALL( MSK_appendsparsesymmat(sdpi->task, dim, nnonz, subi, subj, valij, &idx) );
   MOSEK_CALL( MSK_putbaraij(sdpi->task, row, ind, 1, &idx, &falpha) );

   SCIPdebugMessage("Added a semidefinite constraint of dimension %d in row %d, nonzeroes: %d\n", dim, row, nnonz);

   for (i = 0; i < nnonz; i++)
   {
      SCIPdebugMessage("Entry (%d,%d) => %e\n", subi[i], subj[i], valij[i]);
   }

   return SCIP_OKAY;
}


/** deletes a term of the form <A,X_j> from a row of the SDP */
SCIP_RETCODE SCIPsdpiDelSDPTerm(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ind,                /**< index of the symmetric variable X whose term should be deleted */
   int                   row                 /**< row of the SDP where this term should be deleted */
   )
{
   /* Set the matrix <A,X_j> to a linear combination of zero matrices. This works since
    * the MOSEK documentation under http://docs.mosek.com/7.0/capi/MSK_putbaraij_.html states:
    * "Setting the same elements again will overwrite the earlier entry."
    */
   MOSEK_CALL( MSK_putbaraij(sdpi->task, row, ind, 0, NULL, NULL) );

   return SCIP_OKAY;
}


/** adds columns to the SDP */
SCIP_RETCODE SCIPsdpiAddCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to be added */
   const SCIP_Real*      obj,                /**< objective function values of new columns */
   const SCIP_Real*      lb,                 /**< lower bounds of new columns */
   const SCIP_Real*      ub,                 /**< upper bounds of  columns */
   char**                colnames,           /**< column names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each column in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< row indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{  /*lint --e{715}*/
   int* aptre;
   MSKboundkeye* bkx;
   double* blx;
   double* bux;
   int oldcols;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiAddCols (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiAddCols") );
#endif

   if (ncols == 0)
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&blx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&bux, ncols) );
   generateMskBounds(ncols, lb, ub, bkx, blx, bux);

   MOSEK_CALL( MSK_getnumvar(sdpi->task, &oldcols) );

   MOSEK_CALL( MSK_appendvars(sdpi->task, ncols) );
   MOSEK_CALL( MSK_putcslice(sdpi->task, oldcols, oldcols+ncols, obj) );
   MOSEK_CALL( MSK_putvarboundslice(sdpi->task, oldcols, oldcols+ncols, bkx, blx, bux) );

   if( nnonz > 0 )
   {
      SCIP_CALL( getEndptrs(ncols, beg, nnonz, &aptre) );
      MOSEK_CALL( MSK_putacolslice(sdpi->task, oldcols, oldcols+ncols, beg, aptre, ind, val) );
      BMSfreeMemoryArray(&aptre);
   }

   BMSfreeMemoryArray(&bux);
   BMSfreeMemoryArray(&blx);
   BMSfreeMemoryArray(&bkx);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiAddCols") );
#endif

   return SCIP_OKAY;
}

/** deletes all columns in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to be deleted */
   int                   lastcol             /**< last column to be deleted */
   )
{
   int* sub;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiDelCols (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelCols") );
#endif

   SCIP_CALL( getIndicesRange(firstcol, lastcol, &sub) );

   /*printf("Deleting vars %d to %d\n",firstcol,lastcol);*/
   MOSEK_CALL( MSK_removevars(sdpi->task, lastcol-firstcol+1, sub) );

   BMSfreeMemoryArray(&sub);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelCols") );
#endif

   return SCIP_OKAY;
}

/** deletes columns from SCIP_SDP; the new position of a column must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelColset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of columns
                                              *   input:  1 if column should be deleted, 0 if not
                                              *   output: new position of column, -1 if column was deleted */
   )
{
   int* sub;
   int count;
   int ncols;
   int col;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiDelColset (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelColset") );
#endif

   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, ncols, &count, &sub) );

   col = 0;
   for( i = 0; i < ncols; i++)
   {
      if (dstat[i] == 1)
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = col;
         col++;
      }
   }

   if (count > 0)
   {
      SCIPdebugMessage("Deleting %d vars %d,...\n", count, sub[0]);
      MOSEK_CALL( MSK_removevars(sdpi->task, count, sub) );
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelColset") );
#endif

   return SCIP_OKAY;
}

/** adds rows to the SDP */
SCIP_RETCODE SCIPsdpiAddRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to be added */
   const SCIP_Real*      lhs,                /**< left hand sides of new rows */
   const SCIP_Real*      rhs,                /**< right hand sides of new rows */
   char**                rownames,           /**< row names, or NULL */
   int                   nnonz,              /**< number of nonzero elements to be added to the constraint matrix */
   const int*            beg,                /**< start index of each row in ind- and val-array, or NULL if nnonz == 0 */
   const int*            ind,                /**< column indices of constraint matrix entries, or NULL if nnonz == 0 */
   const SCIP_Real*      val                 /**< values of constraint matrix entries, or NULL if nnonz == 0 */
   )
{  /*lint --e{715}*/
   int* aptre;
   MSKboundkeye* bkc;
   double* blc;
   double* buc;
   int oldrows;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiAddRows (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiAddRows") );
#endif

   if (nrows == 0)
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&blc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&buc, nrows) );

   generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &oldrows) );

   MOSEK_CALL( MSK_appendcons(sdpi->task, nrows) );
   MOSEK_CALL( MSK_putconboundslice(sdpi->task, oldrows, oldrows+nrows, bkc, blc, buc) );

   if( nnonz > 0 )
   {
      SCIP_CALL( getEndptrs(nrows, beg, nnonz, &aptre) );
      MOSEK_CALL( MSK_putarowslice(sdpi->task, oldrows, oldrows+nrows, beg, aptre, ind, val) );
      BMSfreeMemoryArray(&aptre);
   }

   BMSfreeMemoryArray(&buc);
   BMSfreeMemoryArray(&blc);
   BMSfreeMemoryArray(&bkc);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiAddRows") );
#endif

   return SCIP_OKAY;
}

/** deletes all rows in the given range from SDP */
SCIP_RETCODE SCIPsdpiDelRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to be deleted */
   int                   lastrow             /**< last row to be deleted */
   )
{
   int* sub;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiDelRows (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelRows") );
#endif

   SCIP_CALL( getIndicesRange(firstrow, lastrow, &sub) );

   SCIPdebugMessage("Deleting cons %d to %d\n",firstrow,lastrow);

   MOSEK_CALL( MSK_removecons(sdpi->task, lastrow-firstrow+1, sub) );

   BMSfreeMemoryArray(&sub);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelRows") );
#endif

   return SCIP_OKAY;
}

/** deletes rows from SCIP_SDP; the new position of a row must not be greater that its old position */
SCIP_RETCODE SCIPsdpiDelRowset(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  dstat               /**< deletion status of rows
                                              *   input:  1 if row should be deleted, 0 if not
                                              *   output: new position of row, -1 if row was deleted */
   )
{
   int* sub;
   int count;
   int nrows;
   int row;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiDelRowset (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelRowset") );
#endif

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );

   sub = NULL;
   SCIP_CALL( getIndicesFromDense(dstat, nrows, &count, &sub) );

   row = 0;
   for( i = 0; i < nrows; i++ )
   {
      if (dstat[i] == 1)
      {
         dstat[i] = -1;
      }
      else
      {
         dstat[i] = row;
         row++;
      }
   }

   if (count > 0)
   {
      SCIPdebugMessage("Deleting %d cons %d,...\n",count,sub[0]);
      MOSEK_CALL( MSK_removecons(sdpi->task, count, sub) );
      BMSfreeMemoryArray(&sub);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiDelRowset end") );
#endif

   return SCIP_OKAY;
}

/** clears the whole SDP */
SCIP_RETCODE SCIPsdpiClear(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   int nrows;
   int ncols;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiClear (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );

   SCIP_CALL( SCIPsdpiDelRows(sdpi, 0, nrows) );
   SCIP_CALL( SCIPsdpiDelCols(sdpi, 0, ncols) );

   return SCIP_OKAY;
}

/** changes lower and upper bounds of columns */
SCIP_RETCODE SCIPsdpiChgBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change bounds for */
   const int*            ind,                /**< column indices */
   const SCIP_Real*      lb,                 /**< values for the new lower bounds */
   const SCIP_Real*      ub                  /**< values for the new upper bounds */
   )
{
   MSKboundkeye* bkx;
   double* blx;
   double* bux;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiChgBounds (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgBounds") );
#endif

   if (ncols == 0)
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&blx, ncols) );
   SCIP_ALLOC( BMSallocMemoryArray(&bux, ncols) );

   generateMskBounds(ncols, lb, ub, bkx, blx, bux);
   MOSEK_CALL( MSK_putboundlist(sdpi->task, MSK_ACC_VAR, ncols, ind, bkx, blx, bux) );

   BMSfreeMemoryArray(&bux);
   BMSfreeMemoryArray(&blx);
   BMSfreeMemoryArray(&bkx);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgBounds") );
#endif

   return SCIP_OKAY;
}

/** changes left and right hand sides of rows */
SCIP_RETCODE SCIPsdpiChgSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   nrows,              /**< number of rows to change sides for */
   const int*            ind,                /**< row indices */
   const SCIP_Real*      lhs,                /**< new values for left hand sides */
   const SCIP_Real*      rhs                 /**< new values for right hand sides */
   )
{
   MSKboundkeye* bkc;
   double* blc;
   double* buc;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiChgSides (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgSides") );
#endif

   if (nrows == 0)
      return SCIP_OKAY;

   SCIP_ALLOC( BMSallocMemoryArray(&bkc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&blc, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&buc, nrows) );

   generateMskBounds(nrows, lhs, rhs, bkc, blc, buc);
   MOSEK_CALL( MSK_putboundlist(sdpi->task, MSK_ACC_CON, nrows, ind, bkc, blc, buc) );

   BMSfreeMemoryArray(&buc);
   BMSfreeMemoryArray(&blc);
   BMSfreeMemoryArray(&bkc);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgSides") );
#endif

   return SCIP_OKAY;
}

/** changes a single coefficient */
SCIP_RETCODE SCIPsdpiChgCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient to change */
   int                   col,                /**< column number of coefficient to change */
   SCIP_Real             newval              /**< new value of coefficient */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiChgCoef (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgCoef") );
#endif

   MOSEK_CALL( MSK_putaij(sdpi->task, row, col, newval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgCoef") );
#endif

   return SCIP_OKAY;
}

/** changes the objective sense */
SCIP_RETCODE SCIPsdpiChgObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN           objsen              /**< new objective sense */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiChgObjsen (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_putobjsense(sdpi->task, SENSE2MOSEK(objsen)) );

   return SCIP_OKAY;
}

/** changes objective values of columns in the SDP */
SCIP_RETCODE SCIPsdpiChgObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns to change objective value for */
   int*                  ind,                /**< column indices to change objective value for */
   SCIP_Real*            obj                 /**< new objective values for columns */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiChgObj (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiChgObj") );
#endif

   MOSEK_CALL( MSK_putclist(sdpi->task, ncols, ind, obj) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi,"SCIPsdpiChgObj") );
#endif

   return SCIP_OKAY;
}

/** multiplies a row with a non-zero scalar; for negative scalars, the row's sense is switched accordingly */
SCIP_RETCODE SCIPsdpiScaleRow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   int nnonz;
   int* sub;
   double* val;
   MSKboundkeye bkc;
   double blc;
   double buc;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiScaleRow (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiScaleRow") );
#endif

   assert(scaleval != 0);

   MOSEK_CALL( MSK_getarownumnz(sdpi->task, row, &nnonz) );

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

      MOSEK_CALL( MSK_getarow(sdpi->task, row, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putarow(sdpi->task, row, nnonz, sub, val) );

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

   MOSEK_CALL( MSK_getbound(sdpi->task, MSK_ACC_CON, row, &bkc, &blc, &buc) );
   scale_bound(&bkc, &blc, &buc, scaleval);
   MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_CON, row, bkc, blc, buc) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiScaleRow") );
#endif

   return SCIP_OKAY;
}

/** multiplies a column with a non-zero scalar; the objective value is multiplied with the scalar, and the bounds
 *  are divided by the scalar; for negative scalars, the column's bounds are switched
 */
SCIP_RETCODE SCIPsdpiScaleCol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column number to scale */
   SCIP_Real             scaleval            /**< scaling multiplier */
   )
{
   int nnonz;
   int *sub    = NULL;
   double *val = NULL;
   MSKboundkeye bkx;
   double blx, bux, c;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiScaleCol (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiScaleCol") );
#endif

   assert(scaleval != 0);
   MOSEK_CALL( MSK_getacolnumnz(sdpi->task, col, &nnonz) );

   if (nnonz != 0)
   {
      SCIP_ALLOC( BMSallocMemoryArray(&sub, nnonz) );
      SCIP_ALLOC( BMSallocMemoryArray(&val, nnonz) );

      MOSEK_CALL( MSK_getacol(sdpi->task, col, &nnonz, sub, val) );
      scale_vec(nnonz, val, scaleval);
      MOSEK_CALL( MSK_putacol(sdpi->task, col, nnonz, sub, val) );

      BMSfreeMemoryArray(&val);
      BMSfreeMemoryArray(&sub);
   }

   MOSEK_CALL( MSK_getbound(sdpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );
   scale_bound(&bkx, &blx, &bux, 1.0/scaleval);
   MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );

   MOSEK_CALL( MSK_getcslice(sdpi->task, col, col+1, &c) );
   MOSEK_CALL( MSK_putcj(sdpi->task, col, c*scaleval) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiScaleCol") );
#endif

   return SCIP_OKAY;
}


/*
 * Data Accessing Methods
 */

/** gets the number of rows in the SDP */
SCIP_RETCODE SCIPsdpiGetNRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nrows               /**< pointer to store the number of rows */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetNRows (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumcon(sdpi->task, nrows) );

   return SCIP_OKAY;
}

/** gets the number of columns in the SDP */
SCIP_RETCODE SCIPsdpiGetNCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  ncols               /**< pointer to store the number of cols */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetNCols (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumvar(sdpi->task, ncols) );

   return SCIP_OKAY;
}

/** gets the number of nonzero elements in the SDP constraint matrix */
SCIP_RETCODE SCIPsdpiGetNNonz(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nnonz               /**< pointer to store the number of nonzeros */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetNNonz (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumanz(sdpi->task, nnonz) );

   return SCIP_OKAY;
}

static
SCIP_RETCODE getASlice(
   SCIP_SDPI*            sdpi,
   MSKaccmodee           iscon,
   int                   first,
   int                   last,
   int*                  nnonz,
   int*                  beg,
   int*                  ind,
   double*               val
   )
{
   int* aptre;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);
   assert(first <= last);

   SCIPdebugMessage("Calling SCIPsdpiGetNNonz (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "getASlice") );
#endif

   if( nnonz != 0 )
   {
      int surplus;

      assert(beg != NULL);
      assert(ind != NULL);
      assert(val != NULL);

      SCIP_ALLOC( BMSallocMemoryArray(&aptre, last - first + 1) );

      MOSEK_CALL( MSK_getaslicenumnz(sdpi->task, iscon, first, last+1,nnonz) );
      surplus = *nnonz;
      MOSEK_CALL( MSK_getaslice(sdpi->task, iscon, first, last+1, *nnonz, &surplus, beg, aptre, ind, val) );

      assert(surplus == 0);

      BMSfreeMemoryArray(&aptre);
   }

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "getASlice") );
#endif

   return SCIP_OKAY;
}

/** gets the number of nonzero semidefinite terms of the form <A,X_j> in the SDP
  * @note Currently the returned number will be incorrect if we have deleted some of the terms, use with caution.
  */
SCIP_RETCODE SCIPsdpiGetNSDPTerms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpterms           /**< pointer to store the number of SDP terms */
   )
{
   MOSEK_CALL( MSK_getnumbaranz(sdpi->task, (SCIP_Longint*)nsdpterms) );

   return SCIP_OKAY;
}

/** gets the number of semidefinite variables */
SCIP_RETCODE SCIPsdpiGetNSDPVars(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  nsdpvars            /**< pointer to store the number of SDP variables */
   )
{
   MOSEK_CALL( MSK_getnumbarvar(sdpi->task, nsdpvars) );

   return SCIP_OKAY;
}

/** gets columns from SDP problem object; the arrays have to be large enough to store all values;
 *  Either both, lb and ub, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPsdpiGetCols(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get from SDP */
   int                   lastcol,            /**< last column to get from SDP */
   SCIP_Real*            lb,                 /**< buffer to store the lower bound vector, or NULL */
   SCIP_Real*            ub,                 /**< buffer to store the upper bound vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each column in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store column indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetCols (%d)\n",sdpi->sdpid);

   SCIP_CALL( SCIPsdpiGetBounds(sdpi, firstcol, lastcol, lb, ub) );
   SCIP_CALL( getASlice(sdpi, MSK_ACC_VAR, firstcol, lastcol, nnonz, beg, ind, val) );

   return SCIP_OKAY;
}

/** gets rows from SDP problem object; the arrays have to be large enough to store all values.
 *  Either both, lhs and rhs, have to be NULL, or both have to be non-NULL,
 *  either nnonz, beg, ind, and val have to be NULL, or all of them have to be non-NULL.
 */
SCIP_RETCODE SCIPsdpiGetRows(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get from SDP */
   int                   lastrow,            /**< last row to get from SDP */
   SCIP_Real*            lhs,                /**< buffer to store left hand side vector, or NULL */
   SCIP_Real*            rhs,                /**< buffer to store right hand side vector, or NULL */
   int*                  nnonz,              /**< pointer to store the number of nonzero elements returned, or NULL */
   int*                  beg,                /**< buffer to store start index of each row in ind- and val-array, or NULL */
   int*                  ind,                /**< buffer to store row indices of constraint matrix entries, or NULL */
   SCIP_Real*            val                 /**< buffer to store values of constraint matrix entries, or NULL */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetRows (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetRows") );
#endif


   SCIP_CALL( SCIPsdpiGetSides(sdpi, firstrow, lastrow, lhs, rhs) );
   SCIP_CALL( getASlice(sdpi, MSK_ACC_CON, firstrow, lastrow, nnonz, beg, ind, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetRows") );
#endif

   return SCIP_OKAY;
}

/** gets column names */
SCIP_RETCODE SCIPsdpiGetColNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get name from SDP */
   int                   lastcol,            /**< last column to get name from SDP */
   char**                colnames,           /**< pointers to column names (of size at least lastcol-firstcol+1) */
   char*                 namestorage,        /**< storage for col names */
   int                   namestoragesize,    /**< size of namestorage (if 0, storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPsdpiGetColNames() has not been implemented yet.\n");
   return SCIP_ERROR;
}

/** gets row names */
SCIP_RETCODE SCIPsdpiGetRowNames(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get name from SDP */
   int                   lastrow,            /**< last row to get name from SDP */
   char**                rownames,           /**< pointers to row names (of size at least lastrow-firstrow+1) */
   char*                 namestorage,        /**< storage for row names */
   int                   namestoragesize,    /**< size of namestorage (if 0, -storageleft returns the storage needed) */
   int*                  storageleft         /**< amount of storage left (if < 0 the namestorage was not big enough) */
   )
{
   SCIPerrorMessage("SCIPsdpiGetRowNames() has not been implemented yet.\n");
   return SCIP_ERROR;
}

/** gets the objective sense of the SDP */
SCIP_RETCODE SCIPsdpiGetObjsen(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_OBJSEN*          objsen              /**< pointer to store objective sense */
   )
{
   SCIPerrorMessage("SCIPsdpiGetObjsen() has not been implemented yet.\n");
   return SCIP_ERROR;
}

/** gets objective coefficients from SDP problem object */
SCIP_RETCODE SCIPsdpiGetObj(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get objective coefficient for */
   int                   lastcol,            /**< last column to get objective coefficient for */
   SCIP_Real*            vals                /**< array to store objective coefficients */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetObj (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getcslice(sdpi->task, firstcol, lastcol+1, vals) );

   return SCIP_OKAY;
}

/** gets current bounds from SDP problem object */
SCIP_RETCODE SCIPsdpiGetBounds(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstcol,           /**< first column to get bounds for */
   int                   lastcol,            /**< last column to get bounds for */
   SCIP_Real*            lbs,                /**< array to store lower bound values, or NULL */
   SCIP_Real*            ubs                 /**< array to store upper bound values, or NULL */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetBounds (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetBounds") );
#endif

   MOSEK_CALL( MSK_getboundslice(sdpi->task, MSK_ACC_VAR, firstcol, lastcol+1, NULL, lbs, ubs) );

   return SCIP_OKAY;
}

/** gets current row sides from SDP problem object */
SCIP_RETCODE SCIPsdpiGetSides(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   firstrow,           /**< first row to get sides for */
   int                   lastrow,            /**< last row to get sides for */
   SCIP_Real*            lhss,               /**< array to store left hand side values, or NULL */
   SCIP_Real*            rhss                /**< array to store right hand side values, or NULL */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetSides (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetSides") );
#endif

   MOSEK_CALL( MSK_getboundslice(sdpi->task, MSK_ACC_CON, firstrow, lastrow+1, NULL, lhss, rhss) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetSides") );
#endif

   return SCIP_OKAY;
}

/** gets a single coefficient */
SCIP_RETCODE SCIPsdpiGetCoef(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number of coefficient */
   int                   col,                /**< column number of coefficient */
   SCIP_Real*            val                 /**< pointer to store the value of the coefficient */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetCoef (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetCoef") );
#endif

   MOSEK_CALL( MSK_getaij(sdpi->task, row, col, val) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiGetCoef") );
#endif

   return SCIP_OKAY;
}

/** gets a symmetric matrix of the SDP */
SCIP_RETCODE SCIPsdpiGetSymmetricMatrix(
   SCIP_SDPI*           sdpi,                /**< SDP interface structure */
   int                  index,               /**< index of the semidefinite matrix */
   int                  length,              /**< length of the output arrays subi, subj and valij */
   int*                 subi,                /**< Row indices of nonzero entries */
   int*                 subj,                /**< Column indices of nonzero entries */
   SCIP_Real*           valij                /**< values of the nonzero entries defined by the arrays subi and subj */
   )
{
   MOSEK_CALL( MSK_getsparsesymmat(sdpi->task, index, length, subi, subj, valij) );
   return SCIP_OKAY;
}



/*
 * Solving Methods
 */

/** gets the internal solution status of the solver */
static
SCIP_RETCODE getSolutionStatus(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   MSKprostae*           prosta,             /**< pointer to store the problem status */
   MSKsolstae*           solsta              /**< pointer to store the solution status */
   )
{
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getsolutionstatus ( sdpi->task, MSK_SOL_ITR, prosta, solsta) );

   return SCIP_OKAY;
}


static
MSKrescodee filterTRMrescode(
   SCIP_MESSAGEHDLR*     messagehdlr,
   MSKrescodee*          termcode,
   MSKrescodee           res
   )
{
   if (   res == MSK_RES_TRM_MAX_ITERATIONS || res == MSK_RES_TRM_MAX_TIME
      || res == MSK_RES_TRM_OBJECTIVE_RANGE || res == MSK_RES_TRM_STALL
#if ASSERT_ON_NUMERICAL_TROUBLES > 0
      || res == MSK_RES_TRM_MAX_NUM_SETBACKS
      || res == MSK_RES_TRM_NUMERICAL_PROBLEM
#endif
          )
   {
      *termcode = res;
      if (res == MSK_RES_TRM_MAX_NUM_SETBACKS || res == MSK_RES_TRM_NUMERICAL_PROBLEM)
      {
         SCIPmessagePrintWarning(messagehdlr, "Return code %d in [%d]\n", res, optimizecount);

#if ASSERT_ON_WARNING
         assert(0);
#endif
      }

      return MSK_RES_OK;
   }
   else
   {
      *termcode = MSK_RES_OK;
      return res;
   }
}

static
SCIP_RETCODE SolveWSimplex(
   SCIP_SDPI*             sdpi                 /**< SDP interface structure */
   )
{
   int itercount_primal;
   int itercount_dual;
   int gotbasicsol;
   int presolve;
   int maxiter;
   MSKprostae prosta;
   MSKsolstae solsta;
   double pobj,dobj;

   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE, &presolve) );
   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &maxiter) );

#if DEBUG_EASY_REPRODUCE
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_AUTO_SORT_A_BEFORE_OPT, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_OFF) );
#else
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
#endif

   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_AUTO_UPDATE_SOL_INFO, MSK_OFF) );

#if FORCE_MOSEK_LOG

   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_ON) );
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG_SIM_FREQ, 1) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#else
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#endif

   MOSEK_CALL( MSK_solutiondef(sdpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if( gotbasicsol )
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_OFF) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE, MSK_PRESOLVE_MODE_ON) );
   }

#if ALWAYS_SOLVE_PRIMAL > 0
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_SOLVE_FORM, MSK_SOLVE_PRIMAL) );
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SolveWSimplex") );
#endif

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   if( gotbasicsol && maxiter < 20000 )
   {
      /* Since max iter often is set, we switch off restricted pricing */
      MOSEK_CALL( MSK_putintparam(sdpi->task,  MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, 0) );
   }

   if( FORCE_NO_MAXITER )
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, 2000000000) );
   }


#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "Begin optimize with simplex") );
#endif

#if FORCE_MOSEK_SUMMARY > 1
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(sdpi->task,MSK_STREAM_LOG) );
   }
#endif

#if !FORCE_SILENCE
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, 100) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG_SIM_FREQ, 100) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG_SIM, 100) );
#endif

   MOSEK_CALL( filterTRMrescode(sdpi->messagehdlr, &sdpi->termcode, MSK_optimize(sdpi->task)) );

   if( sdpi->termcode == MSK_RES_TRM_MAX_NUM_SETBACKS )
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_SCALING, MSK_SCALING_AGGRESSIVE) );

      MOSEK_CALL( filterTRMrescode(sdpi->messagehdlr, &sdpi->termcode, MSK_optimize(sdpi->task)) );
   }

#if FORCE_MOSEK_SUMMARY
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_solutionsummary(sdpi->task,MSK_STREAM_LOG) );
   }
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "End optimize with simplex") );
#endif

   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE, presolve) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

   MOSEK_CALL( MSK_getintinf(sdpi->task, MSK_IINF_SIM_PRIMAL_ITER, &itercount_primal) );
   MOSEK_CALL( MSK_getintinf(sdpi->task, MSK_IINF_SIM_DUAL_ITER, &itercount_dual) );

   sdpi->itercount = itercount_primal + itercount_dual;

   MOSEK_CALL( MSK_getprimalobj(sdpi->task, MSK_SOL_BAS, &pobj) );
   MOSEK_CALL( MSK_getdualobj(sdpi->task, MSK_SOL_BAS, &dobj) );
   MOSEK_CALL( MSK_getsolutionstatus(sdpi->task, MSK_SOL_BAS, &prosta, &solsta) );

#if DEBUG_PRINT_STAT
   SCIPdebugMessage("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, objval = %g : %g, iter = %d+%d\n",
      maxiter, sdpi->termcode, prosta, solsta, pobj, dobj, itercount_primal, itercount_dual);
#endif

   SCIPdebugMessage("maxiter = %d, termcode = %d, prosta = %d, solsta = %d, "
      "objval = %g : %g, iter = %d+%d\n",
      maxiter,sdpi->termcode,prosta,solsta,
      pobj,dobj,itercount_primal,itercount_dual);

   /*  SCIPdebugMessage("Iter dual %d primal %d\n",itercount_dual,itercount_primal); */
   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_FEAS:
   case MSK_SOL_STA_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
   case MSK_SOL_STA_UNKNOWN:
      break;
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
      SCIPmessagePrintWarning(sdpi->messagehdlr, "Simplex[%d] returned solsta = %d\n", optimizecount, solsta);

      if (sdpi->termcode == MSK_RES_OK)
         sdpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
   case MSK_SOL_STA_INTEGER_OPTIMAL:
   case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
   default:
#if SHOW_ERRORS && !FORCE_SILENCE
      SCIPerrorMessage("Simplex[%d] returned solsta = %d\n", optimizecount, solsta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_ERROR;
   }  /*lint !e788*/

   switch (prosta)
   {
   case MSK_PRO_STA_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_FEAS:
   case MSK_PRO_STA_DUAL_FEAS:
   case MSK_PRO_STA_PRIM_AND_DUAL_INFEAS:
   case MSK_PRO_STA_PRIM_INFEAS:
   case MSK_PRO_STA_DUAL_INFEAS:
   case MSK_PRO_STA_UNKNOWN:
      break;
   case MSK_PRO_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_PRO_STA_NEAR_PRIM_FEAS:
   case MSK_PRO_STA_NEAR_DUAL_FEAS:
   case MSK_PRO_STA_ILL_POSED:
   case MSK_PRO_STA_PRIM_INFEAS_OR_UNBOUNDED:
      SCIPmessagePrintWarning(sdpi->messagehdlr, "Simplex[%d] returned prosta = %d\n", optimizecount, prosta);

      if (sdpi->termcode == MSK_RES_OK)
         sdpi->termcode = MSK_RES_TRM_NUMERICAL_PROBLEM;

#if ASSERT_ON_WARNING
      assert(0);
#endif
      break;
   default:
#if SHOW_ERRORS && !FORCE_SILENCE
      SCIPerrorMessage("Simplex[%d] returned prosta = %d\n", optimizecount, prosta);
#endif

#if ASSERT_ON_WARNING
      assert(0);
#endif

      return SCIP_ERROR;
   }  /*lint !e788*/

   if( solsta == MSK_SOL_STA_OPTIMAL && fabs(dobj)+fabs(dobj) > 1.0e-6 && fabs(pobj-dobj)>0.0001*(fabs(pobj)+fabs(dobj)))
   {
      SCIPerrorMessage("Simplex[%d] returned optimal solution with different objvals %g != %g reldiff %.2g%%\n",
         optimizecount, pobj, dobj, 100*fabs(pobj-dobj)/ MAX(fabs(pobj),fabs(dobj)));
   }

   if (sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
   {
      if(solsta != MSK_SOL_STA_DUAL_FEAS && solsta != MSK_SOL_STA_OPTIMAL && solsta != MSK_SOL_STA_PRIM_AND_DUAL_FEAS)
      {
         SCIPerrorMessage("[%d] Terminated on objective range without dual feasible solsta.\n", optimizecount);

         SCIP_CALL( SCIPsdpiSolveBarrier(sdpi, TRUE) );
      }
      else
      {
         scipmskobjsen objsen;
         double bound;

         MOSEK_CALL( MSK_getobjsense(sdpi->task, &objsen) );

         if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
         {
            MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_UPPER_OBJ_CUT, &bound) );

            if (1.0e-6*(fabs(bound)+fabs(dobj)) < bound-dobj)
            {
               SCIPerrorMessage("[%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount, dobj, bound);

               SCIP_CALL( SCIPsdpiSolveBarrier(sdpi, TRUE) );
            }
         }
         else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
         {
            MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_LOWER_OBJ_CUT, &bound) );

            if (1.0e-6*(fabs(bound)+fabs(dobj)) < dobj-bound)
            {
               SCIPerrorMessage("[%d] Terminated on obj range, dobj = %g, bound = %g\n",
                  optimizecount, dobj, bound);

               SCIP_CALL( SCIPsdpiSolveBarrier(sdpi, TRUE) );
            }
         }
      }
   }

   if (maxiter >= 2000000000)
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, maxiter) );

      if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      {
         SCIPmessagePrintWarning(sdpi->messagehdlr, "Simplex[%d] failed to terminate in 10000 iterations, switching to interior point\n",
            optimizecount);

         SCIP_CALL( SCIPsdpiSolveBarrier(sdpi,TRUE) );
      }
   }

#if DEBUG_DO_INTPNT_FEAS_CHECK
   if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
   {
      SCIPdebugMessage("Checking infeasibility[%d]... ",optimizecount);

      SCIP_CALL( SCIPsdpiSolveBarrier(sdpi,true) );

      MOSEK_CALL(MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, &prosta, &solsta));

      if (solsta == MSK_SOL_STA_PRIM_INFEAS_CER || solsta == MSK_SOL_STA_DUAL_INFEAS_CER)
      {
         SCIPdebugPrintf("ok\n");
      }
      else
      {
         SCIPdebugPrintf("wrong [%d] prosta = %d, solsta = %d\n",optimizecount,prosta,solsta);
      }
   }
#endif


#if DEBUG_PRINT_STAT > 0
   SCIPdebugMessage("Max iter stat    : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      optimizecount, numstrongbranchmaxiterup, numstrongbranchmaxiterdo, numprimalmaxiter, numdualmaxiter);
   SCIPdebugMessage("Objcut iter stat : Count %d branchup = %d branchlo = %d primal %d dual %d\n",
      optimizecount, numstrongbranchobjup, numstrongbranchobjdo, numprimalobj, numdualobj);
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SolveWSimplex") );
#endif

   return SCIP_OKAY;
}

/** solves the SDP */
SCIP_RETCODE SCIPsdpiSolve(
   SCIP_SDPI*             sdpi                 /**< SDP interface structure */
   )
{
   MSKrescodee trmcode;
   MOSEK_CALL( MSK_optimizetrm(sdpi->task, &trmcode) );

   return SCIP_OKAY;
}



/** calls primal simplex to solve the SDP */
SCIP_RETCODE SCIPsdpiSolvePrimal(
   SCIP_SDPI*             sdpi                 /**< SDP interface structure */
   )
{
   optimizecount++;

   SCIPdebugMessage("Calling SCIPsdpiSolvePrimal[%d] (%d) ",optimizecount,sdpi->sdpid);

   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiSolvePrimal") );
#endif

   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_PRIMAL_SIMPLEX) );

#if WRITE_PRIMAL > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"primal_%d.sdp",optimizecount);
      SCIPdebugMessage("\nWriting sdp %s\n",fname);
      /*MOSEK_CALL( MSK_putintparam(sdpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
      MSK_writedata(sdpi->task,fname);
   }
#endif

   SCIP_CALL( SolveWSimplex(sdpi) );

   if ( sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, NULL, &solsta) );


      if( solsta != MSK_SOL_STA_PRIM_FEAS )
      {
         SCIP_CALL( SolveWSimplex(sdpi) );
      }
   }

   if (sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numprimalobj;

   if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numprimalmaxiter;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiSolvePrimal") );
#endif

   return SCIP_OKAY;
}

/** calls dual simplex to solve the SDP */
SCIP_RETCODE SCIPsdpiSolveDual(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   optimizecount++;

   SCIPdebugMessage("Calling SCIPsdpiSolveDual[%d] (%d)\n",optimizecount,sdpi->sdpid);


   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_INTEGER, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_DUAL_SIMPLEX) );

#if WRITE_DUAL > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"dual_%d.sdp",optimizecount);
      SCIPdebugMessage("\nWriting sdp %s\n",fname);
      MSK_writedata(sdpi->task,fname);
   }
#endif

#if !FORCE_SILENCE
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG_SIM_FREQ, 1) );
#endif

   SCIP_CALL( SolveWSimplex(sdpi) );

   if ( sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE )
   {
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, NULL, &solsta) );

      if( solsta != MSK_SOL_STA_DUAL_FEAS )
      {
         SCIP_CALL( SolveWSimplex(sdpi) );
      }
   }

   if (sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
      ++numdualobj;

   if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;

   return SCIP_OKAY;
}

/** calls barrier or interior point algorithm to solve the SDP with crossover to simplex basis */
SCIP_RETCODE SCIPsdpiSolveBarrier(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool             crossover           /**< perform crossover */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   optimizecount++;

#if FORCE_MOSEK_LOG
   if( optimizecount > WRITE_ABOVE )
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_ON) );
   }
   else
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#else
   {
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, MSK_OFF) );
   }
#endif


   SCIPdebugMessage("Calling SCIPsdpiSolveBarrier[%d] (%d) ",optimizecount,sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiSolveBarrier") );
#endif

   MOSEK_CALL( MSK_putintparam(sdpi->task,MSK_IPAR_INTPNT_BASIS, crossover ? MSK_BI_ALWAYS : MSK_BI_NEVER) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT) );


#if WRITE_INTPNT > 0
   if( optimizecount > WRITE_ABOVE )
   {
      char fname[40];
      snprintf(fname,40,"intpnt_%d.sdp",optimizecount);
      SCIPdebugMessage("\nWriting sdp %s\n",fname);
      /*MOSEK_CALL( MSK_putintparam(sdpi->task,MSK_IPAR_WRITE_GENERIC_NAMES,MSK_ON) );*/
      MSK_writedata(sdpi->task,fname);
   }
#endif

   MOSEK_CALL( filterTRMrescode(sdpi->messagehdlr, &sdpi->termcode, MSK_optimize(sdpi->task)) );

   if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
      ++numdualmaxiter;

   MOSEK_CALL( MSK_getintinf(sdpi->task, MSK_IINF_INTPNT_ITER, &sdpi->itercount) );

#ifdef SCIP_DEBUG
   {
      MSKprostae prosta;
      MSKsolstae solsta;

      MOSEK_CALL( MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, &prosta, &solsta) );
      SCIPdebugMessage("termcode = %d, prosta = %d, solsta = %d, iter = %d\n",
         sdpi->termcode, prosta, solsta, sdpi->itercount);
   }
#endif

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiSolveBarrier") );
#endif

   return SCIP_OKAY;
}

/** start strong branching - call before any strong branching */
SCIP_RETCODE SCIPsdpiStartStrongbranch(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** end strong branching - call after any strong branching */
SCIP_RETCODE SCIPsdpiEndStrongbranch(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   /* currently do nothing */
   return SCIP_OKAY;
}

/** performs strong branching iterations on all candidates */
static
SCIP_RETCODE SCIPsdpiStrongbranch(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   MSKobjsensee objsen;
   int olditerlim;
   int oldselection;
   int oldhotstart;

   double bound;
   int ncols;
   int nrows;
   MSKboundkeye bkx;
   double blx;
   double bux;
   double newub;
   double newlb;

   SCIPdebugMessage("Calling SCIPsdpiStrongbranch (%d)\n",sdpi->sdpid);

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiStrongbranch") );
#endif

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   if (sdpi->termcode != MSK_RES_OK)
   {
      SCIPmessagePrintWarning(sdpi->messagehdlr, "SB Warning: Previous termcode is %d\n",sdpi->termcode);
   }

   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );

   SCIP_CALL( getbase(sdpi, ncols, nrows) );

   MOSEK_CALL( MSK_getobjsense(sdpi->task, &objsen) );
   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, &olditerlim) );
   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_DUAL_SELECTION, &oldselection) );
   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART, &oldhotstart) );

   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, itlim) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_DUAL_SELECTION, STRONGBRANCH_PRICING) );

   if (objsen == MSK_OBJECTIVE_SENSE_MINIMIZE)
   {
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_UPPER_OBJ_CUT, &bound) );
   }
   else /* objsen == MSK_OBJECTIVE_SENSE_MAX */
   {
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_LOWER_OBJ_CUT, &bound) );
   }

   MOSEK_CALL( MSK_getbound(sdpi->task, MSK_ACC_VAR, col, &bkx, &blx, &bux) );

   *iter = 0;

   newub = EPSCEIL(psol-1.0, 1e-06);

   if (newub < blx - 0.5) /* infeasible */
   {
      *down = bound;
      *downvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_NEGINF(blx))
         newbk = MSK_BK_UP;
      else if (EPSEQ(blx,newub,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newub = blx;
      }
      else
         newbk = MSK_BK_RA;

      MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_VAR, col, newbk, blx, newub) );

      SCIP_CALL( SCIPsdpiSolveDual(sdpi) );

      *iter += sdpi->itercount;

      if (SCIPsdpiIsStable(sdpi))
         *downvalid = TRUE;
      else
         *downvalid = FALSE;

      if (SCIPsdpiExistsPrimalRay(sdpi))
      {
         SCIPmessagePrintWarning(sdpi->messagehdlr, "SB ERROR: Lp [%d] is dual infeasible\n",optimizecount);

         *down = -1e20;
         *downvalid = FALSE;
      }
      else if (SCIPsdpiExistsDualRay(sdpi))
      {
         *down = bound;
      }
      else
      {
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, NULL, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(sdpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n", optimizecount);

            *down = -1e20;
            *downvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(sdpi->task, MSK_SOL_BAS, down) );
         }
      }

      if (sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjup;

      if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterup;
   }

   /* Reset basis solution before doing the up branch */
   MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
   SCIP_CALL( setbase(sdpi) );

   newlb = EPSFLOOR(psol+1.0, 1e-06);
   if (newlb > bux + 0.5) /* infeasible */
   {
      *up = bound;
      *upvalid = TRUE;
   }
   else
   {
      MSKboundkeye newbk;

      if (IS_POSINF(bux))
         newbk = MSK_BK_LO;
      else if (EPSEQ(bux,newlb,1.0e-6))
      {
         newbk = MSK_BK_FX;
         newlb = bux;
      }
      else
         newbk = MSK_BK_RA;

      MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_VAR, col, newbk, newlb, bux) );
      SCIP_CALL( SCIPsdpiSolveDual(sdpi) );

      *iter += sdpi->itercount;

      if (SCIPsdpiIsStable(sdpi))
         *upvalid = TRUE;
      else
         *upvalid = FALSE;

      if (SCIPsdpiExistsPrimalRay(sdpi))
      {
         *up = -1e20;
         *upvalid = FALSE;
      }
      else if (SCIPsdpiExistsDualRay(sdpi))
      {
         *up = bound;
      }
      else
      {
         SCIP_Bool dfeas;

         SCIP_CALL( SCIPsdpiGetSolFeasibility(sdpi, NULL, &dfeas) );

         if (!dfeas)
         {
            SCIPmessagePrintWarning(sdpi->messagehdlr, "SB ERROR: Lp [%d] is not dual feasible\n",optimizecount);

            *up = -1e20;
            *upvalid = FALSE;
         }
         else
         {
            MOSEK_CALL( MSK_getdualobj(sdpi->task, MSK_SOL_BAS, up) );
         }
      }

      if (sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE)
         ++numstrongbranchobjdo;

      if (sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS)
         ++numstrongbranchmaxiterdo;
   }

   MOSEK_CALL( MSK_putbound(sdpi->task, MSK_ACC_VAR, col, bkx, blx, bux) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, olditerlim) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_DUAL_SELECTION, oldselection) );
   MOSEK_CALL( MSK_putintparam(sdpi->task,   MSK_IPAR_SIM_HOTSTART, oldhotstart) );

   SCIP_CALL( setbase(sdpi) );

   sdpi->termcode = MSK_RES_OK;
   sdpi->itercount = 0;

#if DEBUG_CHECK_DATA > 0
   SCIP_CALL( scip_checkdata(sdpi, "SCIPsdpiStrongbranch") );
#endif

   SCIPdebugMessage("End SCIPsdpiStrongbranch (%d)\n", sdpi->sdpid);

   return SCIP_OKAY;
}

/** performs strong branching iterations on one @b fractional candidate */
SCIP_RETCODE SCIPsdpiStrongbranchFrac(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< fractional current primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   /* pass call on to sdpiStrongbranch() */
   SCIP_CALL( SCIPsdpiStrongbranch(sdpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given @b fractional candidates */
SCIP_RETCODE SCIPsdpiStrongbranchesFrac(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< fractional current primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to sdpiStrongbranch() */
      SCIP_CALL( SCIPsdpiStrongbranch(sdpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}

/** performs strong branching iterations on one candidate with @b integral value */
SCIP_RETCODE SCIPsdpiStrongbranchInt(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   col,                /**< column to apply strong branching on */
   SCIP_Real             psol,               /**< current integral primal solution value of column */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bound after branching column down */
   SCIP_Real*            up,                 /**< stores dual bound after branching column up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up value is a valid dual bound;
                                              *   otherwise, it can only be used as an estimate value */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   /* pass call on to sdpiStrongbranch() */
   SCIP_CALL( SCIPsdpiStrongbranch(sdpi, col, psol, itlim, down, up, downvalid, upvalid, iter) );

   return SCIP_OKAY;
}

/** performs strong branching iterations on given candidates with @b integral values */
SCIP_RETCODE SCIPsdpiStrongbranchesInt(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cols,               /**< columns to apply strong branching on */
   int                   ncols,              /**< number of columns */
   SCIP_Real*            psols,              /**< current integral primal solution values of columns */
   int                   itlim,              /**< iteration limit for strong branchings */
   SCIP_Real*            down,               /**< stores dual bounds after branching columns down */
   SCIP_Real*            up,                 /**< stores dual bounds after branching columns up */
   SCIP_Bool*            downvalid,          /**< stores whether the returned down values are valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   SCIP_Bool*            upvalid,            /**< stores whether the returned up values are a valid dual bounds;
                                              *   otherwise, they can only be used as an estimate values */
   int*                  iter                /**< stores total number of strong branching iterations, or -1; may be NULL */
   )
{
   int j;

   assert( iter != NULL );
   assert( cols != NULL );
   assert( psols != NULL );
   assert( down != NULL );
   assert( up != NULL );
   assert( downvalid != NULL );
   assert( upvalid != NULL );
   assert( down != NULL );

   if ( iter != NULL )
      *iter = 0;

   for (j = 0; j < ncols; ++j)
   {
      /* pass call on to sdpiStrongbranch() */
      SCIP_CALL( SCIPsdpiStrongbranch(sdpi, cols[j], psols[j], itlim, &(down[j]), &(up[j]), &(downvalid[j]), &(upvalid[j]), iter) );
   }
   return SCIP_OKAY;
}


/*
 * Solution Information Methods
 */


/** returns whether a solve method was called after the last modification of the SDP */
SCIP_Bool SCIPsdpiWasSolved(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiWasSolved (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus (sdpi, &prosta, &solsta) );

   return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** gets information about primal and dual feasibility of the current SDP solution */
SCIP_RETCODE SCIPsdpiGetSolFeasibility(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            primalfeasible,     /**< stores primal feasibility status */
   SCIP_Bool*            dualfeasible        /**< stores dual feasibility status */
   )
{
   MSKsolstae solsta;
   SCIP_Bool pfeas;
   SCIP_Bool dfeas;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetSolFeasibility (%d)\n",sdpi->sdpid);

   pfeas = FALSE;
   dfeas = FALSE;

   MOSEK_CALL( MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, NULL, &solsta) );

   switch (solsta)
   {
   case MSK_SOL_STA_OPTIMAL:
   case MSK_SOL_STA_PRIM_AND_DUAL_FEAS:
      pfeas = TRUE;
      dfeas = TRUE;
      break;
   case MSK_SOL_STA_PRIM_FEAS:
      pfeas = TRUE;
      break;
   case MSK_SOL_STA_DUAL_FEAS:
      dfeas = TRUE;
      break;
   case MSK_SOL_STA_UNKNOWN:
   case MSK_SOL_STA_NEAR_OPTIMAL:
   case MSK_SOL_STA_NEAR_PRIM_FEAS:
   case MSK_SOL_STA_NEAR_DUAL_FEAS:
   case MSK_SOL_STA_NEAR_PRIM_AND_DUAL_FEAS:
   case MSK_SOL_STA_PRIM_INFEAS_CER:
   case MSK_SOL_STA_DUAL_INFEAS_CER:
   case MSK_SOL_STA_NEAR_PRIM_INFEAS_CER:
   case MSK_SOL_STA_NEAR_DUAL_INFEAS_CER:
   case MSK_SOL_STA_INTEGER_OPTIMAL:
   case MSK_SOL_STA_NEAR_INTEGER_OPTIMAL:
      break;
   default:
      return SCIP_ERROR;
   }  /*lint !e788*/

   if( primalfeasible != NULL )
      *primalfeasible = pfeas;

   if( dualfeasible != NULL )
      *dualfeasible = dfeas;

   return SCIP_OKAY;
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point);
 *  this does not necessarily mean, that the solver knows and can return the primal ray
 */
SCIP_Bool SCIPsdpiExistsPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiExistsPrimalRay (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus (sdpi, &prosta, &solsta));

   return (   solsta == MSK_SOL_STA_DUAL_INFEAS_CER
      || prosta == MSK_PRO_STA_DUAL_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff SDP is proven to have a primal unbounded ray (but not necessary a primal feasible point),
 *  and the solver knows and can return the primal ray
 */
SCIP_Bool SCIPsdpiHasPrimalRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiHasPrimalRay (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus (sdpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_DUAL_INFEAS_CER);
}

/** returns TRUE iff SDP is proven to be primal unbounded */
SCIP_Bool SCIPsdpiIsPrimalUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be primal infeasible */
SCIP_Bool SCIPsdpiIsPrimalInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return SCIPsdpiExistsDualRay(sdpi);
}

/** returns TRUE iff SDP is proven to be primal feasible */
SCIP_Bool SCIPsdpiIsPrimalFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKprostae prosta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiIsPrimalFeasible (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus (sdpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_PRIM_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point);
 *  this does not necessarily mean, that the solver knows and can return the dual ray
 */
SCIP_Bool SCIPsdpiExistsDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKprostae prosta;
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiExistsDualRay (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus(sdpi, &prosta, &solsta) );

   return (   solsta == MSK_SOL_STA_PRIM_INFEAS_CER
      || prosta == MSK_PRO_STA_PRIM_INFEAS
      || prosta == MSK_PRO_STA_PRIM_AND_DUAL_INFEAS);
}

/** returns TRUE iff SDP is proven to have a dual unbounded ray (but not necessary a dual feasible point),
 *  and the solver knows and can return the dual ray
 */
SCIP_Bool SCIPsdpiHasDualRay(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiHasDualRay (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus (sdpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_PRIM_INFEAS_CER);
}

/** returns TRUE iff SDP is proven to be dual unbounded */
SCIP_Bool SCIPsdpiIsDualUnbounded(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   return FALSE;
}

/** returns TRUE iff SDP is proven to be dual infeasible */
SCIP_Bool SCIPsdpiIsDualInfeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   return SCIPsdpiExistsPrimalRay(sdpi);
}

/** returns TRUE iff SDP is proven to be dual feasible */
SCIP_Bool SCIPsdpiIsDualFeasible(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKprostae prosta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiIsDualFeasible (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus(sdpi, &prosta, NULL) );

   return (prosta == MSK_PRO_STA_DUAL_FEAS || prosta == MSK_PRO_STA_PRIM_AND_DUAL_FEAS);
}


/** returns TRUE iff SDP was solved to optimality */
SCIP_Bool SCIPsdpiIsOptimal(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKsolstae solsta;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiIsOptimal (%d)\n",sdpi->sdpid);

   SCIP_CALL_ABORT( getSolutionStatus(sdpi, NULL, &solsta) );

   return (solsta == MSK_SOL_STA_OPTIMAL);
}

/** returns TRUE iff current SDP basis is stable */
SCIP_Bool SCIPsdpiIsStable(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return (   sdpi->termcode == MSK_RES_OK
      || sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS
      || sdpi->termcode == MSK_RES_TRM_MAX_TIME
      || sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE);
}

/** returns TRUE iff the objective limit was reached */
SCIP_Bool SCIPsdpiIsObjlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return sdpi->termcode == MSK_RES_TRM_OBJECTIVE_RANGE;
}

/** returns TRUE iff the iteration limit was reached */
SCIP_Bool SCIPsdpiIsIterlimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return sdpi->termcode == MSK_RES_TRM_MAX_ITERATIONS;
}

/** returns TRUE iff the time limit was reached */
SCIP_Bool SCIPsdpiIsTimelimExc(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return sdpi->termcode == MSK_RES_TRM_MAX_TIME;
}

/** returns the internal solution status of the solver */
int SCIPsdpiGetInternalStatus(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   MSKsolstae solsta;
   SCIP_RETCODE retcode;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetInternalStatus (%d)\n", sdpi->sdpid);

   retcode = getSolutionStatus(sdpi, NULL, &solsta);
   if ( retcode != SCIP_OKAY )
      return 0;

   return solsta; /*lint !e641*/
}

/** tries to reset the internal status of the SDP solver in order to ignore an instability of the last solving call */
SCIP_RETCODE SCIPsdpiIgnoreInstability(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Bool*            success             /**< pointer to store, whether the instability could be ignored */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiIgnoreInstability (%d)\n",sdpi->sdpid);

   *success = FALSE;

   return SCIP_OKAY;
}

/** gets objective value of solution */
SCIP_RETCODE SCIPsdpiGetObjval(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval              /**< stores the objective value */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetObjval (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getprimalobj(sdpi->task, MSK_SOL_ITR, objval) );

   /* TODO: tjek lighed med dual objektiv i de fleste tilfaelde. */

   return SCIP_OKAY;
}

/** gets primal and dual solution vectors */
SCIP_RETCODE SCIPsdpiGetSol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            objval,             /**< stores the objective value, may be NULL if not needed */
   SCIP_Real*            primsol,            /**< primal solution vector, may be NULL if not needed */
   SCIP_Real*            dualsol,            /**< dual solution vector, may be NULL if not needed */
   SCIP_Real*            activity,           /**< row activity vector, may be NULL if not needed */
   SCIP_Real*            redcost             /**< reduced cost vector, may be NULL if not needed */
   )
{  /*lint --e{715}*/
   double* sux;
   int ncols;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetSol (%d)\n",sdpi->sdpid);

   sux = NULL;
   ncols = 0;

   if( objval )
   {
      MOSEK_CALL( MSK_getprimalobj(sdpi->task, MSK_SOL_ITR, objval) );
   }

   if( redcost )
   {
      MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );
      SCIP_ALLOC( BMSallocMemoryArray( &sux, ncols) );
   }

   MOSEK_CALL( MSK_getsolution(sdpi->task, MSK_SOL_ITR, NULL, NULL, NULL, NULL, NULL, activity,
         primsol, dualsol, NULL, NULL, redcost, sux, NULL) );

   if( redcost )
   {
      for( i = 0; i < ncols; i++ )
      {
         assert(sux != NULL);
         redcost[i] -= sux[i];
      }
   }
   if( sux )
   {
      BMSfreeMemoryArray(&sux);
   }

   return SCIP_OKAY;
}

/** gets the primal solution for a semidefinite variable */
SCIP_RETCODE SCIPsdpiGetSolSDPVar(
   SCIP_SDPI*           sdpi,                /**< SDP interface structure */
   int                  ind,                 /**< index of semidefinite variable */
   SCIP_Real*           barxj                /**< primal semidefinite solution vector, must have size dim*(dim+1)/2 */
   )
{
   MOSEK_CALL(MSK_getbarxj(sdpi->task, MSK_SOL_ITR, ind, barxj));
   return SCIP_OKAY;
}

/** gets primal ray for unbounded SDPs */
SCIP_RETCODE SCIPsdpiGetPrimalRay(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            ray                 /**< primal ray */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetPrimalRay (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getsolution(sdpi->task, MSK_SOL_BAS, NULL, NULL, NULL, NULL, NULL, NULL, ray,
         NULL, NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets dual Farkas proof for infeasibility */
SCIP_RETCODE SCIPsdpiGetDualfarkas(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real*            dualfarkas          /**< dual Farkas row multipliers */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetDualfarkas (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getsolution(sdpi->task, MSK_SOL_BAS, NULL, NULL, NULL, NULL, NULL, NULL, NULL, dualfarkas,
         NULL, NULL, NULL, NULL, NULL) );

   return SCIP_OKAY;
}

/** gets the number of SDP iterations of the last solve call */
SCIP_RETCODE SCIPsdpiGetIterations(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  iterations          /**< pointer to store the number of iterations of the last solve call */
   )
{
   SCIPdebugMessage("Calling SCIPsdpiGetIterations (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   *iterations = sdpi->itercount;

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
   assert(sdpi != NULL);
   assert(quality != NULL);

   *quality = SCIP_INVALID;

   return SCIP_OKAY;
}

/** handle singular basis */
static
SCIP_RETCODE handle_singular(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  basis,              /**< array of basis indices */
   MSKrescodee           res                 /**< result */
   )
{
   if (res == MSK_RES_ERR_BASIS_SINGULAR)
   {
      SCIP_CALL( SCIPsdpiSolvePrimal(sdpi) );

      MOSEK_CALL( MSK_initbasissolve(sdpi->task, basis) );
   }
   else
   {
      MOSEK_CALL( res );
   }

   return SCIP_OKAY;
}


/*
 * SDP Basis Methods
 */

/** convert Mosek status to SCIP status */
static
SCIP_RETCODE convertstat_mosek2scip(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   MSKaccmodee           acc,                /**< ??? */
   MSKstakeye*           sk,                 /**< ??? */
   int                   n,                  /**< size */
   int*                  stat                /**< status array */
   )
{
   int i;

   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;

      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
         MOSEK_CALL( MSK_getsolutioni(sdpi->task, acc, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_INF:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_LOW:
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UPR:
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
      case MSK_SK_END:
         break;
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert Mosek to SCIP status - slack variables */
static
SCIP_RETCODE convertstat_mosek2scip_slack(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   MSKaccmodee           acc,                /**< ??? */
   MSKstakeye*           sk,                 /**< ??? */
   int                   n,                  /**< size */
   int*                  stat                /**< status array */
   )
{
   int i;

   /* slacks are stored as -1 in Mosek, i.e., bounds are reversed compared to SCIP  */

   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_BAS:
         stat[i] = (int)SCIP_BASESTAT_BASIC;
         break;
      case MSK_SK_SUPBAS:
         stat[i] = (int)SCIP_BASESTAT_ZERO;
         break;
      case MSK_SK_FIX:
         MOSEK_CALL( MSK_getsolutioni(sdpi->task, acc, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );

         if (sl < su) /* Negative reduced cost */
            stat[i] = (int)SCIP_BASESTAT_UPPER;
         else
            stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_UNK:
      case MSK_SK_INF:
      case MSK_SK_UPR: /* Reversed */
         stat[i] = (int)SCIP_BASESTAT_LOWER;
         break;
      case MSK_SK_LOW: /* Reversed */
         stat[i] = (int)SCIP_BASESTAT_UPPER;
         break;
      case MSK_SK_END:
         break;
      default:
         SCIPABORT();
         return SCIP_INVALIDDATA; /*lint !e527*/
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** convert SCIP to Mosek status */
static
void convertstat_scip2mosek(
   int*                  stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   int i;
   for( i = 0; i < n; i++ )
   {
      switch (stat[i])
      {
      case SCIP_BASESTAT_LOWER:
         resstat[i] = MSK_SK_LOW;
         break;
      case SCIP_BASESTAT_BASIC:
         resstat[i] = MSK_SK_BAS;
         break;
      case SCIP_BASESTAT_UPPER:
         resstat[i] = MSK_SK_UPR;
         break;
      case SCIP_BASESTAT_ZERO:
         resstat[i] = MSK_SK_SUPBAS;
         break;
      default:
         SCIPABORT();
      }
   }
}

/** convert SCIP to Mosek status - slack variables */
static
void convertstat_scip2mosek_slack(
   int*                  stat,               /**< SCIP status array */
   int                   n,                  /**< size of array */
   MSKstakeye*           resstat             /**< resulting Mosek status array */
   )
{
   /* slacks are stored as -1 in Mosek, i.e., bounds are reversed compared to SCIP  */
   int i;

   for( i = 0; i < n; i++ )
   {
      switch (stat[i])
      {
      case SCIP_BASESTAT_LOWER:
         resstat[i] = MSK_SK_UPR;/* Reversed */
         break;
      case SCIP_BASESTAT_BASIC:
         resstat[i] = MSK_SK_BAS;
         break;
      case SCIP_BASESTAT_UPPER:
         resstat[i] = MSK_SK_LOW; /* Reversed */
         break;
      case SCIP_BASESTAT_ZERO:
         resstat[i] = MSK_SK_SUPBAS;
         break;
      default:
         SCIPABORT();
      }
   }
}

/** gets current basis status for columns and rows; arrays must be large enough to store the basis status */
SCIP_RETCODE SCIPsdpiGetBase(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cstat,              /**< array to store column basis status, or NULL */
   int*                  rstat               /**< array to store row basis status, or NULL */
   )
{
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPsdpiGetBase (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );

   SCIP_CALL( getbase(sdpi, ncols, nrows) );

   if (cstat)
   {
      SCIP_CALL( convertstat_mosek2scip(sdpi, MSK_ACC_VAR, sdpi->skx, ncols, cstat) );
   }

   if (rstat)
   {
      SCIP_CALL( convertstat_mosek2scip_slack(sdpi, MSK_ACC_CON, sdpi->skc, nrows, rstat) );
   }

   return SCIP_OKAY;
}

/** sets current basis status for columns and rows */
SCIP_RETCODE SCIPsdpiSetBase(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  cstat,              /**< array with column basis status */
   int*                  rstat               /**< array with row basis status */
   )
{
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPsdpiSetBase (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );
   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );

   SCIP_CALL( ensureStateMem(sdpi, ncols, nrows) );

   convertstat_scip2mosek(cstat, ncols, sdpi->skx);
   convertstat_scip2mosek_slack(rstat, nrows, sdpi->skc);

   SCIP_CALL( setbase(sdpi) );

   return SCIP_OKAY;
}

/** returns the indices of the basic columns and rows; basic column n gives value n, basic row m gives value -1-m */
extern
SCIP_RETCODE SCIPsdpiGetBasisInd(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int*                  bind                /**< pointer to store basis indices ready to keep number of rows entries */
   )
{
   int nrows;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetBasisInd (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );

#if 0
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_OFF) );
#endif

   SCIP_CALL( handle_singular(sdpi,bind,MSK_initbasissolve(sdpi->task, bind)) );

#if 0
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );
#endif

   for (i = 0; i < nrows; i++ )
   {
      if (bind[i] < nrows) /* row bind[i] is basic */
         bind[i] = -1 - bind[i];
      else                 /* column bind[i]-nrows is basic */
         bind[i] = bind[i] - nrows;
   }

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPsdpiGetBInvCol(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   c,                  /**< column number of B^-1; this is NOT the number of the column in the SDP;
                                              *   you have to call SCIPsdpiGetBasisInd() to get the array which links the
                                              *   B^-1 column numbers to the row and column numbers of the SDP!
                                              *   c must be between 0 and nrows-1, since the basis has the size
                                              *   nrows * nrows */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the column */
   )
{
   int* sub;
   int nrows;
   int numnz;
   int i;

   SCIPdebugMessage("Calling SCIPsdpiGetBInvCol (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getnumcon(sdpi->task,&nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );

   for (i=0; i<nrows; i++)
      coef[i] = 0;

   numnz = 1;
   sub[0]= c;
   coef[c] = 1; /* Unit vector e_col */

   MOSEK_CALL( MSK_putnaintparam(sdpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_OFF) );
   MOSEK_CALL( MSK_solvewithbasis(sdpi->task, 0, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}

/** get dense column of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPsdpiGetBInvACol(
   SCIP_SDPI*            sdpi,                /**< SDP interface structure */
   int                   c,                  /**< column number */
   SCIP_Real*            coef                /**< vector to return coefficients */
   )
{  /*lint --e{715}*/
   SCIP_Real* val;
   int* sub;
   int nrows;
   int numnz;
   int i;

   SCIPdebugMessage("Calling SCIPsdpiGetBInvACol (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getnumcon(sdpi->task,&nrows) );
   MOSEK_CALL( MSK_getacolnumnz(sdpi->task,c,&numnz) );
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &val, numnz+1) );

   for (i=0; i<nrows; i++)
      coef[i] = 0;

   MOSEK_CALL( MSK_getacol(sdpi->task, c, &numnz, sub, val) );

   for (i=0; i<numnz; i++)
      coef[sub[i]] = val[i];

   MOSEK_CALL( MSK_putnaintparam(sdpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_OFF) );
   MOSEK_CALL( MSK_solvewithbasis(sdpi->task, 0, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);
   BMSfreeMemoryArray(&val);
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART_LU, MSK_ON) );

   return SCIP_OKAY;
}


/** get dense row of inverse basis matrix B^-1 */
SCIP_RETCODE SCIPsdpiGetBInvRow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number */
   SCIP_Real*            coef                /**< pointer to store the coefficients of the row */
   )
{
   int* sub;
   int nrows;
   int numnz;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetBInvRow (%d)\n",sdpi->sdpid);

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );
   SCIP_ALLOC( BMSallocMemoryArray( &sub, nrows) );

   for (i=0; i<nrows; i++)
      coef[i] = 0;

   numnz = 1;
   sub[0] = row;
   coef[row] = 1; /* Unit vector e_row */

   MOSEK_CALL( MSK_putnaintparam(sdpi->task, MSK_IPAR_BASIS_SOLVE_USE_PLUS_ONE_, MSK_ON) );
   MOSEK_CALL( MSK_solvewithbasis(sdpi->task, 1, &numnz, sub, coef) );

   BMSfreeMemoryArray(&sub);

   SCIPdebugMessage("End SCIPsdpiGetBInvRow (%d)\n",sdpi->sdpid);

   return SCIP_OKAY;
}

/** get dense row of inverse basis matrix times constraint matrix B^-1 * A */
SCIP_RETCODE SCIPsdpiGetBInvARow(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   row,                /**< row number */
   const SCIP_Real*      binvrow,            /**< row in (A_B)^-1 from prior call to SCIPsdpiGetBInvRow(), or NULL */
   SCIP_Real*            val                 /**< vector to return coefficients */
   )
{
   int nrows;
   int ncols;
   int numnz;
   int* csub;
   int didalloc;
   double* cval;
   double* binv;
   int i;
   int k;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiGetBInvARow (%d)\n",sdpi->sdpid);

   didalloc = 0;

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );

   SCIP_ALLOC( BMSallocMemoryArray(&csub, nrows) );
   SCIP_ALLOC( BMSallocMemoryArray(&cval, nrows) );

   for( i = 0; i < ncols; i++ )
      val[i] = 0;

   if( binvrow == NULL )
   {
      didalloc = 1;

      SCIP_ALLOC( BMSallocMemoryArray( &binv, nrows) );
      SCIP_CALL( SCIPsdpiGetBInvRow(sdpi, row, binv) );
   }
   else
      binv = (SCIP_Real*)binvrow;

   /* binvrow*A */
   for( i = 0; i < ncols; i++)
   {
      val[i] = 0;

      MOSEK_CALL( MSK_getacol(sdpi->task, i, &numnz, csub, cval) );

      for( k = 0; k < numnz; ++k )
         val[i] += binv[csub[k]] * cval[k];
   }

   /* free memory arrays */
   BMSfreeMemoryArray(&cval);
   BMSfreeMemoryArray(&csub);

   if( didalloc > 0 )
   {
      BMSfreeMemoryArray(&binv);
   }

   return SCIP_OKAY;
}

/*
 * SDP State Methods
 */

/** creates SDPi state information object */
static
SCIP_RETCODE sdpistateCreate(
   SCIP_SDPISTATE**      sdpistate,          /**< pointer to SDPi state */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   ncols,              /**< number of columns to store */
   int                   nrows               /**< number of rows to store */
   )
{
   assert(sdpistate != NULL);
   assert(blkmem != NULL);
   assert(ncols >= 0);
   assert(nrows >= 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, sdpistate) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*sdpistate)->skx, colpacketNum(ncols)) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*sdpistate)->skc, rowpacketNum(nrows)) );

   if( sdpistate[0] )
   {
      sdpistate[0]->solsta = MSK_SOL_STA_UNKNOWN;
      sdpistate[0]->num    = -1;
      sdpistate[0]->ncols  = ncols;
      sdpistate[0]->nrows  = nrows;
   }

   return SCIP_OKAY;
}

/** frees SDPi state information */
static
void sdpistateFree(
   SCIP_SDPISTATE**      sdpistate,          /**< pointer to SDPi state information (like basis information) */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(blkmem != NULL);
   assert(sdpistate != NULL);
   assert(*sdpistate != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*sdpistate)->skx, colpacketNum((*sdpistate)->ncols));
   BMSfreeBlockMemoryArray(blkmem, &(*sdpistate)->skc, rowpacketNum((*sdpistate)->nrows));
   BMSfreeBlockMemory(blkmem, sdpistate);
}

#ifndef NDEBUG
/** check state */
static
SCIP_RETCODE checkState1(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   n,                  /**< number of rows or columns */
   MSKstakeye*           sk,                 /**< ??? */
   MSKaccmodee           accmode,            /**< ??? */
   char                  xc                  /**< ??? */
   )
{
   int i;

   /* printout for all except LOW, UPR, FIX and BAS with sl[xc]==su[xc] */
   for( i = 0; i < n; i++ )
   {
      double sl;
      double su;
      switch (sk[i])
      {
      case MSK_SK_UNK:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unk\n", optimizecount, xc, i);
         break;
      case MSK_SK_BAS:
         MOSEK_CALL( MSK_getsolutioni(sdpi->task, accmode, i, MSK_SOL_BAS, NULL, NULL, &sl, &su, NULL) );
         if (fabs(sl-su) > DEBUG_CHECK_STATE_TOL)
            SCIPdebugMessage("STATE[%d]: %c[%d] = bas, sl%c = %g, su%c = %g\n", optimizecount, xc, i, xc, sl, xc, su);
         break;
      case MSK_SK_SUPBAS:
         SCIPdebugMessage("STATE[%d]: %c[%d] = supbas\n", optimizecount, xc, i);
         break;
      case MSK_SK_LOW:
      case MSK_SK_UPR:
      case MSK_SK_FIX:
         break;
      case MSK_SK_INF:
         SCIPdebugMessage("STATE[%d]: %c[%d] = inf\n", optimizecount, xc, i);
         break;
      default:
         SCIPdebugMessage("STATE[%d]: %c[%d] = unknown status <%d>\n", optimizecount, xc, i, sk[i]);
         break;
      }  /*lint !e788*/
   }

   return SCIP_OKAY;
}

/** check state */
static
SCIP_RETCODE checkState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   int                   ncols,              /**< number of columns */
   int                   nrows               /**< number of rows */
   )
{
   SCIP_CALL( checkState1(sdpi, ncols, sdpi->skx, MSK_ACC_VAR, 'x') );
   SCIP_CALL( checkState1(sdpi, nrows, sdpi->skc, MSK_ACC_CON, 'c') );

   return SCIP_OKAY;
 }
#endif

/** store row and column basis status in a packed SDPi state object */
static
SCIP_RETCODE sdpistatePack(
   SCIP_SDPI*             sdpi,                /**< SDP interface structure */
   SCIP_SDPISTATE*        sdpistate            /**< pointer to SDPi state data */
   )
{
   int *skxi = (int *) sdpi->skx; /* Used as temp. buffer */
   int *skci = (int *) sdpi->skc; /* Used as temp. buffer */

   assert(sizeof(int) == sizeof(MSKstakeye));

   SCIP_CALL( convertstat_mosek2scip(sdpi, MSK_ACC_VAR, sdpi->skx, sdpistate->ncols, skxi) );
   SCIP_CALL( convertstat_mosek2scip_slack(sdpi, MSK_ACC_CON, sdpi->skc, sdpistate->nrows, skci) );

   SCIPencodeDualBit(skxi, sdpistate->skx, sdpistate->ncols);
   SCIPencodeDualBit(skci, sdpistate->skc, sdpistate->nrows);

   return SCIP_OKAY;
}

/** unpacks row and column basis status from a packed SDPi state object */
static
void sdpistateUnpack(
   SCIP_SDPISTATE*       sdpistate,          /**< pointer to SDPi state data */
   MSKstakeye*           skx,                /**< ??? */
   MSKstakeye*           skc                 /**< ??? */
   )
{
   assert(sizeof(int) == sizeof(MSKstakeye));

   SCIPdecodeDualBit(sdpistate->skx, (int*) skx, sdpistate->ncols);
   SCIPdecodeDualBit(sdpistate->skc, (int*) skc, sdpistate->nrows);

   convertstat_scip2mosek((int*) skx, sdpistate->ncols, skx);
   convertstat_scip2mosek_slack((int*) skc, sdpistate->nrows, skc);
}

/** stores SDP state (like basis information) into sdpistate object */
SCIP_RETCODE SCIPsdpiGetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDP state information (like basis information) */
   )
{
   int gotbasicsol;
   int nrows;
   int ncols;

   SCIPdebugMessage("Calling SCIPsdpiGetState (%d)\n",sdpi->sdpid);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);
   assert(sdpistate != NULL);

   *sdpistate = NULL;

   MOSEK_CALL( MSK_solutiondef(sdpi->task, MSK_SOL_BAS, &gotbasicsol) );

   if ( gotbasicsol == 0 || SCIPsdpiExistsDualRay(sdpi) )
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );

   /* allocate sdpistate data */
   SCIP_CALL( sdpistateCreate(sdpistate, blkmem, ncols, nrows) );

   sdpistate[0]->num = optimizecount;

   MOSEK_CALL(MSK_getsolutionstatus ( sdpi->task, MSK_SOL_BAS, NULL, &sdpistate[0]->solsta));

   SCIP_CALL( getbase(sdpi, ncols, nrows) );

#ifndef NDEBUG
   SCIP_CALL( checkState(sdpi, ncols, nrows) );
#endif

   SCIP_CALL( sdpistatePack(sdpi, sdpistate[0]) );

   SCIPdebugMessage("Store into state from iter : %d\n",optimizecount);

   /*    if (r != SCIP_OKAY)
    *    sdpistateFree(sdpistate, blkmem );
    */

   return SCIP_OKAY;
}

/** loads SDPi state (like basis information) into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetState()
 */
SCIP_RETCODE SCIPsdpiSetState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE*       sdpistate           /**< SDP state information (like basis information) */
   )
{  /*lint --e{715}*/
   int nrows;
   int ncols;
   int i;

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   if (sdpistate == NULL)
   {
      SCIPdebugMessage("Setting NULL state\n");
      return SCIP_OKAY;
   }

   if (sdpistate->nrows == 0 || sdpistate->ncols == 0)
      return SCIP_OKAY;

   MOSEK_CALL( MSK_getnumcon(sdpi->task, &nrows) );
   MOSEK_CALL( MSK_getnumvar(sdpi->task, &ncols) );
   assert(sdpistate->nrows <= nrows);
   assert(sdpistate->ncols <= ncols);

   SCIP_CALL( ensureStateMem(sdpi, ncols, nrows) );
   SCIP_CALL( getbase(sdpi, ncols, nrows) );

   sdpistateUnpack(sdpistate, sdpi->skx, sdpi->skc);

   /* extend the basis to the current SDP beyond the previously existing columns */
   for (i = sdpistate->ncols; i < ncols; ++i)
   {
      SCIP_Real lb;
      SCIP_Real ub;
      MOSEK_CALL( MSK_getboundslice(sdpi->task, MSK_ACC_VAR, i, i, NULL, &lb, &ub) );
      if ( SCIPsdpiIsInfinity(sdpi, REALABS(lb)) )
      {
         /* if lower bound is +/- infinity -> try upper bound */
         if ( SCIPsdpiIsInfinity(sdpi, REALABS(ub)) )
            sdpi->skx[i] = MSK_SK_SUPBAS;  /* variable is free (super basic) */
         else
            sdpi->skx[i] = MSK_SK_UPR;     /* use finite upper bound */
      }
      else
         sdpi->skx[i] = MSK_SK_LOW;        /* use finite lower bound */
   }
   for (i = sdpistate->nrows; i < nrows; ++i)
      sdpi->skc[i] = MSK_SK_BAS;

   /* load basis information into MOSEK */
   SCIP_CALL( setbase(sdpi) );

   SCIPdebugMessage("Store from state into task iter : %d with solsta : %d\n", sdpistate->num, sdpistate->solsta);

   return SCIP_OKAY;
}

/** clears current SDPi state (like basis information) of the solver */
SCIP_RETCODE SCIPsdpiClearState(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{
   assert(sdpi != NULL);

   /**@todo implement SCIPsdpiClearState() for MOSEK */
   SCIPmessagePrintWarning(sdpi->messagehdlr, "MOSEK interface does not implement SCIPsdpiClearState()\n");

   return SCIP_OKAY;
}

/** frees SDP state information */
SCIP_RETCODE SCIPsdpiFreeState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPISTATE**      sdpistate           /**< pointer to SDP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiFreeState (%d)\n",sdpi->sdpid);

   if( *sdpistate != NULL )
   {
      sdpistateFree(sdpistate, blkmem);
   }

   return SCIP_OKAY;
}

/** checks, whether the given SDP state contains simplex basis information */
SCIP_Bool SCIPsdpiHasStateBasis(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPISTATE*       sdpistate           /**< SDP state information (like basis information) */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("Calling SCIPsdpiHasStateBasis (%d)\n",sdpi->sdpid);

   return ( sdpistate != NULL && sdpistate->num >= 0);
}

/** reads SDP state (like basis information from a file */
SCIP_RETCODE SCIPsdpiReadState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   SCIPdebugMessage("reading SDP state from file <%s>\n", fname);

   MOSEK_CALL( MSK_readsolution(sdpi->task, MSK_SOL_BAS, fname) );

   return SCIP_OKAY;
}

/** writes SDP state (like basis information) to a file */
SCIP_RETCODE SCIPsdpiWriteState(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   SCIPdebugMessage("writing SDP state to file <%s>\n", fname);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   /* set parameter to be able to write */
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_WRITE_SOL_HEAD, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_WRITE_SOL_VARIABLES, MSK_ON) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_WRITE_SOL_CONSTRAINTS, MSK_ON) );

   MOSEK_CALL( MSK_writesolution(sdpi->task, MSK_SOL_BAS, fname) );

   return SCIP_OKAY;
}




/*
 * SDP Pricing Norms Methods
 */

/**@name SDP Pricing Norms Methods */
/**@{ */

/** stores SDPi pricing norms information
 *  @todo should we store norm information?
 */
SCIP_RETCODE SCIPsdpiGetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
   )
{
   assert(sdpinorms != NULL);

   (*sdpinorms) = NULL;

   return SCIP_OKAY;
}

/** loads SDPi pricing norms into solver; note that the SDP might have been extended with additional
 *  columns and rows since the state was stored with SCIPsdpiGetNorms()
 */
SCIP_RETCODE SCIPsdpiSetNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS*       sdpinorms           /**< SDPi pricing norms information */
   )
{
   assert(sdpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/** frees pricing norms information */
SCIP_RETCODE SCIPsdpiFreeNorms(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SDPINORMS**      sdpinorms           /**< pointer to SDPi pricing norms information */
   )
{
   assert(sdpinorms == NULL);

   /* no work necessary */
   return SCIP_OKAY;
}

/**@} */

/*
 * Parameter Methods
 */

/** constant array containing the parameter names */
static const char* paramname[] = {
   "SCIP_SDPPAR_FROMSCRATCH",                 /**< solver should start from scratch at next call? */
   "SCIP_SDPPAR_FASTMIP",                     /**< fast mip setting of SDP solver */
   "SCIP_SDPPAR_SCALING",                     /**< should SDP solver use scaling? */
   "SCIP_SDPPAR_PRESOLVING",                  /**< should SDP solver use presolving? */
   "SCIP_SDPPAR_PRICING",                     /**< pricing strategy */
   "SCIP_SDPPAR_SDPINFO",                      /**< should SDP solver output information to the screen? */
   "SCIP_SDPPAR_FEASTOL",                     /**< feasibility tolerance for primal variables and slacks */
   "SCIP_SDPPAR_DUALFEASTOL",                 /**< feasibility tolerance for dual variables and reduced costs */
   "SCIP_SDPPAR_BARRIERCONVTOL",              /**< convergence tolerance used in barrier algorithm */
   "SCIP_SDPPAR_LOBJLIM",                     /**< lower objective limit */
   "SCIP_SDPPAR_UOBJLIM",                     /**< upper objective limit */
   "SCIP_SDPPAR_SDPITLIM",                     /**< SDP iteration limit */
   "SCIP_SDPPAR_SDPTILIM",                     /**< SDP time limit */
   "SCIP_SDPPAR_MARKOWITZ",                   /**< Markowitz tolerance */
   "SCIP_SDPPAR_ROWREPSWITCH",                /**< simplex algorithm shall use row representation of the basis
                                               *  if number of rows divided by number of columns exceeds this value */
   "SCIP_SDPPAR_THREADS"                      /**< number of threads used to solve the SDP */
};

/** method mapping parameter index to parameter name */
static
const char* paramty2str(
   SCIP_SDPPARAM          type
   )
{  /*lint --e{641}*/
   /* check if the parameters in this order */
   assert(SCIP_SDPPAR_FROMSCRATCH == 0);      /**< solver should start from scratch at next call? */
   assert(SCIP_SDPPAR_FASTMIP == 1);          /**< fast mip setting of SDP solver */
   assert(SCIP_SDPPAR_SCALING == 2);          /**< should SDP solver use scaling? */
   assert(SCIP_SDPPAR_PRESOLVING == 3);       /**< should SDP solver use presolving? */
   assert(SCIP_SDPPAR_PRICING == 4);          /**< pricing strategy */
   assert(SCIP_SDPPAR_SDPINFO == 5);           /**< should SDP solver output information to the screen? */
   assert(SCIP_SDPPAR_FEASTOL == 6);          /**< feasibility tolerance for primal variables and slacks */
   assert(SCIP_SDPPAR_DUALFEASTOL == 7);      /**< feasibility tolerance for dual variables and reduced costs */
   assert(SCIP_SDPPAR_BARRIERCONVTOL == 8);   /**< convergence tolerance used in barrier algorithm */
   assert(SCIP_SDPPAR_LOBJLIM == 9);          /**< lower objective limit */
   assert(SCIP_SDPPAR_UOBJLIM == 10);         /**< upper objective limit */
   assert(SCIP_SDPPAR_SDPITLIM == 11);         /**< SDP iteration limit */
   assert(SCIP_SDPPAR_SDPTILIM == 12);         /**< SDP time limit */
   assert(SCIP_SDPPAR_MARKOWITZ == 13);       /**< Markowitz tolerance */
   assert(SCIP_SDPPAR_ROWREPSWITCH == 14);    /**< row representation switch */
   assert(SCIP_SDPPAR_THREADS == 15);         /**< number of threads used to solve the SDP */

   return paramname[type];
}

/** gets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiGetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int*                  ival                /**< buffer to store the parameter value */
                              )
{  /*lint --e{641}*/
   SCIPdebugMessage("getting int parameter %s\n", paramty2str(type));

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   switch (type)
   {
   case SCIP_SDPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART, ival) );
      *ival = (*ival == MSK_SIM_HOTSTART_NONE);
      break;
   case SCIP_SDPPAR_FASTMIP:                   /* fast mip setting of SDP solver */
      return  SCIP_PARAMETERUNKNOWN;
   case SCIP_SDPPAR_SCALING:                   /* should SDP solver use scaling? */
      MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_SCALING, ival) );
      *ival = (*ival != MSK_SCALING_NONE);
      break;
   case SCIP_SDPPAR_PRESOLVING:                /* should SDP solver use presolving? */
      MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE, ival) );
      *ival = (*ival != MSK_PRESOLVE_MODE_OFF);
      break;
   case SCIP_SDPPAR_PRICING:                   /* pricing strategy */
      *ival = sdpi->pricing;
      break;
   case SCIP_SDPPAR_SDPINFO:                    /* should SDP solver output information to the screen? */
      MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_LOG, ival) );
      *ival = (*ival == MSK_ON);
      break;
   case SCIP_SDPPAR_SDPITLIM:                   /* SDP iteration limit */
      MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** sets integer parameter of SDP */
SCIP_RETCODE SCIPsdpiSetIntpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   int                   ival                /**< parameter value */
   )
{
   int scaling;

#if SCIP_CONTROLS_PRICING
   /*lint --e{641}*/
   static int pricing[7] = {
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_FULL,
      MSK_SIM_SELECTION_PARTIAL,
      MSK_SIM_SELECTION_SE,
      MSK_SIM_SELECTION_ASE,
      MSK_SIM_SELECTION_DEVEX,
   };
#endif

   SCIPdebugMessage("Calling SCIPsdpiSetIntpar (%d) Parameter=<%s>  Value=<%d>\n", sdpi->sdpid, paramty2str(type), ival);

   assert(SCIP_PRICING_LPIDEFAULT == 0);
   assert(SCIP_PRICING_AUTO == 1);
   assert(SCIP_PRICING_FULL == 2);
   assert(SCIP_PRICING_PARTIAL == 3);
   assert(SCIP_PRICING_STEEP == 4);
   assert(SCIP_PRICING_STEEPQSTART == 5);
   assert(SCIP_PRICING_DEVEX == 6);

   SCIPdebugMessage("Calling SCIPsdpiSetIntpar (%d) %s = %d\n",sdpi->sdpid,paramty2str(type),ival);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   switch (type)
   {
   case SCIP_SDPPAR_FROMSCRATCH:               /* solver should start from scratch at next call? */
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_HOTSTART,
            ival ? MSK_SIM_HOTSTART_NONE : MSK_SIM_HOTSTART_STATUS_KEYS ) );
      break;
   case SCIP_SDPPAR_FASTMIP:                   /* fast mip setting of SDP solver */
      return SCIP_PARAMETERUNKNOWN;
   case SCIP_SDPPAR_SCALING:                   /* should SDP solver use scaling? */
      scaling = (ival ? MSK_SCALING_FREE : MSK_SCALING_NONE);
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_SCALING, scaling) );
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_INTPNT_SCALING, scaling) );
      break;
   case SCIP_SDPPAR_PRESOLVING:                /* should SDP solver use presolving? */
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_PRESOLVE_USE,
            ival ? MSK_PRESOLVE_MODE_FREE : MSK_PRESOLVE_MODE_OFF) );

#ifdef SCIP_DEBUG
      if( ival )
      {
         SCIPdebugMessage("Setting presolve to on\n");
      }
#endif
      break;
   case SCIP_SDPPAR_PRICING:                   /* pricing strategy */
      assert(ival >= 0 && ival <= SCIP_PRICING_DEVEX);
      sdpi->pricing = (SCIP_PRICING)ival;

#ifdef SCIP_DEBUG
      switch( (SCIP_PRICING)ival )
      {
      case SCIP_PRICING_AUTO:
         SCIPdebugMessage("Setting pricing to auto\n");
         break;
      case SCIP_PRICING_FULL:
         SCIPdebugMessage("Setting pricing to full\n");
         break;
      case SCIP_PRICING_PARTIAL:
         SCIPdebugMessage("Setting pricing to partial\n");
         break;
      case SCIP_PRICING_LPIDEFAULT:
         SCIPdebugMessage("Setting pricing to lpi default\n");
         break;
      case SCIP_PRICING_STEEP:
         SCIPdebugMessage("Setting pricing to steep\n");
         break;
      case SCIP_PRICING_STEEPQSTART:
         SCIPdebugMessage("Setting pricing to steep quick start\n");
         break;
      case SCIP_PRICING_DEVEX:
         SCIPdebugMessage("Setting pricing to devex\n");
         break;
      }
#endif

#if SCIP_CONTROLS_PRICING
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_PRIMAL_SELECTION, pricing[ival]) );

      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_DUAL_SELECTION, pricing[ival]) );

      if( !(sdpi->pricing == SCIP_PRICING_PARTIAL || sdpi->pricing == SCIP_PRICING_AUTO ) )
      {
         /* No restrict */
         MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_DUAL_RESTRICT_SELECTION, 0) );

         MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_PRIMAL_RESTRICT_SELECTION, 0) );
      }
#else
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_PRIMAL_SELECTION, MSK_SIM_SELECTION_FREE) );

      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_DUAL_SELECTION, MSK_SIM_SELECTION_FREE) );
#endif
      break;
   case SCIP_SDPPAR_SDPINFO:
      /* should SDP solver output information to the screen? */
#if FORCE_MOSEK_LOG
      SCIPdebugMessage("Ignoring log setting!\n");
#else
      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_LOG, ival ? MSK_ON : MSK_OFF) );
#endif
      break;
   case SCIP_SDPPAR_SDPITLIM:                   /* SDP iteration limit */
#if DEBUG_PARAM_SETTING
      if( ival )
      {
         SCIPdebugMessage("Setting max iter to : %d\n",ival);
      }
#endif

      MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_SIM_MAX_ITERATIONS, ival) );
      break;
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** gets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiGetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real*            dval                /**< buffer to store the parameter value */
   )
{
   SCIPdebugMessage("getting real parameter %s\n", paramty2str(type));

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   switch (type)
   {
#if SCIP_CONTROLS_TOLERANCES
   case SCIP_SDPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_SDPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_SDPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
#endif
   case SCIP_SDPPAR_LOBJLIM:                   /* lower objective limit */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      break;
   case SCIP_SDPPAR_UOBJLIM:                   /* upper objective limit */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      break;
   case SCIP_SDPPAR_SDPTILIM:                   /* SDP time limit */
      MOSEK_CALL( MSK_getdouparam(sdpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_SDPPAR_MARKOWITZ:                 /* Markowitz tolerance */
   default:
      return SCIP_PARAMETERUNKNOWN;
   } /*lint !e788*/

   return SCIP_OKAY;
}

/** sets floating point parameter of SDP */
SCIP_RETCODE SCIPsdpiSetRealpar(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_SDPPARAM         type,               /**< parameter number */
   SCIP_Real             dval                /**< parameter value */
   )
{
   SCIPdebugMessage("setting real parameter %s to %g\n", paramty2str(type), dval);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   /**@todo Limits shouldn't be hardcoded */

   switch (type)
   {
#if SCIP_CONTROLS_TOLERANCES
   case SCIP_SDPPAR_FEASTOL:                   /* feasibility tolerance for primal variables and slacks */
      if (dval < 1e-9)
         dval = 1e-9;

      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_BASIS_TOL_X, dval) );
      break;
   case SCIP_SDPPAR_DUALFEASTOL:               /* feasibility tolerance for dual variables and reduced costs */
      if (dval < 1e-9)
         return SCIP_PARAMETERUNKNOWN;
      /*         dval = 1e-9; */

      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_BASIS_TOL_S, dval) );
      break;
   case SCIP_SDPPAR_BARRIERCONVTOL:            /* convergence tolerance used in barrier algorithm */
      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_INTPNT_TOL_REL_GAP, dval) );
      break;
#endif
   case SCIP_SDPPAR_LOBJLIM:                   /* lower objective limit */
      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_LOWER_OBJ_CUT, dval) );
      break;
   case SCIP_SDPPAR_UOBJLIM:                   /* upper objective limit */
      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_UPPER_OBJ_CUT, dval) );
      break;
   case SCIP_SDPPAR_SDPTILIM:                   /* SDP time limit */
      MOSEK_CALL( MSK_putdouparam(sdpi->task, MSK_DPAR_OPTIMIZER_MAX_TIME, dval) );
      break;
   case SCIP_SDPPAR_MARKOWITZ:                 /* Markowitz tolerance */
   default:
      return SCIP_PARAMETERUNKNOWN;
   }  /*lint !e788*/

   return SCIP_OKAY;
}


/*
 * Numerical Methods
 */


/** returns value treated as infinity in the SDP solver */
SCIP_Real SCIPsdpiInfinity(
   SCIP_SDPI*            sdpi                /**< SDP interface structure */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return MSK_INFINITY;
}

/** checks if given value is treated as infinity in the SDP solver */
SCIP_Bool SCIPsdpiIsInfinity(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   SCIP_Real             val                 /**< value to be checked for infinity */
   )
{  /*lint --e{715}*/
   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   return IS_POSINF(val);
}


/*
 * File Interface Methods
 */

/**@todo read/write are still from the LPI */

/** reads SDP from a file */
SCIP_RETCODE SCIPsdpiReadSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   int olddataformat;

   SCIPdebugMessage("Calling SCIPsdpiReadSDP (%d), filename <%s>\n", sdpi->sdpid, fname);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_READ_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_READ_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_readdata(sdpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_READ_DATA_FORMAT, olddataformat) );

   return SCIP_OKAY;
}

/** writes SDP to a file */
SCIP_RETCODE SCIPsdpiWriteSDP(
   SCIP_SDPI*            sdpi,               /**< SDP interface structure */
   const char*           fname               /**< file name */
   )
{
   int olddataformat;

   SCIPdebugMessage("Calling SCIPsdpiReadSDP (%d), filename <%s>\n", sdpi->sdpid, fname);

   assert(MosekEnv != NULL);
   assert(sdpi != NULL);
   assert(sdpi->task != NULL);

   MOSEK_CALL( MSK_getintparam(sdpi->task, MSK_IPAR_WRITE_DATA_FORMAT, &olddataformat) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_WRITE_DATA_FORMAT, MSK_DATA_FORMAT_LP) );
   MOSEK_CALL( MSK_writedata(sdpi->task, fname) );
   MOSEK_CALL( MSK_putintparam(sdpi->task, MSK_IPAR_WRITE_DATA_FORMAT, olddataformat) );

   return SCIP_OKAY;
}
