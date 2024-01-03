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

/* This file is copied from SCIP. */

#include "scip/scip.h"
#include "locale.h"

#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wshadow"
#pragma GCC diagnostic ignored "-Wredundant-decls"
#pragma GCC diagnostic ignored "-Wstrict-prototypes"
#pragma GCC diagnostic ignored "-Wdeclaration-after-statement"
#endif

#include <criterion/criterion.h>
#include <criterion/redirect.h>
#include <criterion/parameterized.h>
#include <criterion/theories.h>

#ifdef __GNUC__
#pragma GCC diagnostic warning "-Wredundant-decls"
#pragma GCC diagnostic warning "-Wstrict-prototypes"
#pragma GCC diagnostic warning "-Wshadow"
#endif

#undef SCIP_CALL
#define SCIP_CALL(x)   do                                                                                     \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             cr_assert(FALSE, "Error <%d> in function call\n", _restat_);                     \
                          }                                                                                   \
                       }                                                                                      \
                       while( FALSE )


CR_API int main(int argc, char *argv[])
{
   struct criterion_test_set* tests = criterion_initialize();
   int result = 0;

   if ( criterion_handle_args(argc, argv, true) )
   {
      setlocale(LC_ALL, "C");
      result = !criterion_run_all_tests(tests);
   }

   criterion_finalize(tests);
   return result;
}
