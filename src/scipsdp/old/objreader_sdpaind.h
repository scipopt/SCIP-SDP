/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/* This file is part of SCIPSDP - a solving framework for mixed-integer      */
/* semidefinite programs based on SCIP.                                      */
/*                                                                           */
/* Copyright (C) 2011-2013 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
/*               2014-2020 Discrete Optimization, TU Darmstadt               */
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
/* Copyright (C) 2002-2020 Zuse Institute Berlin                             */
/* SCIP is distributed under the terms of the SCIP Academic Licence,         */
/* see file COPYING in the SCIP distribution.                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objreader_sdpaind.h
 * @brief  Reader for SDPA-Files with indicator constraints (-var in linear constraint => indicator constraint with var as indicator variable)
 * @author Jakob Schelbert
 * @author Sonja Mars
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJREADER_SDPAIND_H__
#define __SCIP_OBJREADER_SDPAIND_H__

#include <utility>                      // for pair
#include <vector>                       // for vector

#include "objscip/objreader.h"          // for ObjReader
#include "objreader_sdpa.h"             // for LProw, LPblock, SDPblock
#include "scip/scip.h"

namespace scip
{
   /** C++ wrapper object for file readers */
   class ObjReaderSDPAind : public ObjReader
   {
   public:

      /** default constructor */
      ObjReaderSDPAind(SCIP* scip)
      : ObjReader(scip, "sdpaindreader", "file reader for SDPA files with indicator constraints", "dat-s-ind")
      {}

      /** destructor */
      virtual ~ObjReaderSDPAind()
      {}

      /** problem reading method of reader
       *
       *  possible return values for *result:
       *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
       *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
       *
       *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
       */
      virtual SCIP_DECL_READERREAD(scip_read);

   };

} /* namespace scip */


#endif
