/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the                                 */
/*      SDP-Package for SCIP: a solving framework for                        */
/*                            mixed-integer semidefinite programms           */
/*                                                                           */
/* Copyright (C) 2011-2012 Discrete Optimization, TU Darmstadt               */
/*                         EDOM, FAU Erlangen-Nürnberg                       */
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
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objreader_sdpa.h
 * @brief  Reader for SDPA-Files
 * @author Jakob Schelbert, Sonja Mars
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJREADER_SDPA_H__
#define __SCIP_OBJREADER_SDPA_H__

#include <utility>                      // for pair
#include <vector>                       // for vector
#include "objscip/objreader.h"          // for ObjReader
#include "scip/type_reader.h"           // for SCIP_READER
#include "scip/type_result.h"           // for SCIP_RESULT
#include "scip/type_retcode.h"          // for SCIP_RETCODE
#include "scip/type_scip.h"             // for SCIP

class SdpProblem;
class SdpVarMapper;

namespace scip
{
   /**struct with the lp-rows-data*/
   struct LProw
   {
      std::vector< std::pair<int, double> > data;
   };
   /**class of one sdpblock, very similar to sdpcone and the sdpa-format*/
   class SDPBlock
   {
   public:
      SDPBlock (int size) : blocksize(size), num_nonzeros(0), constnum_nonzeros(0) { }
      ~SDPBlock() {}

      int blocksize;
      int num_nonzeros;
      std::vector<int> variables;
      std::vector<int> columns;
      std::vector<int> rows;
      std::vector<double> values;

      std::vector<int> constcolumns;
      std::vector<int> constrows;
      std::vector<double> constvalues;
      int constnum_nonzeros;
   };

   /**class for the lp-blocks*/
   class LPBlock
   {
   public:
      LPBlock () { }
      ~LPBlock() {}

      int numrows;
      std::vector<LProw> rows;
   };

   /** C++ wrapper object for file readers */
   class ObjReaderSDPA : public ObjReader
   {
   public:

      /** default constructor */
   ObjReaderSDPA(SCIP* scip)
      : ObjReader(scip, "sdpareader", "file reader for SDPA files", "dat-s")
      {}

      /** destructor */
      virtual ~ObjReaderSDPA()
      {}

      /** problem reading method of reader
       *
       *  possible return values for *result:
       *  - SCIP_SUCCESS    : the reader read the file correctly and created an appropritate problem
       *  - SCIP_DIDNOTRUN  : the reader is not responsible for given input file
       *
       *  If the reader detected an error in the input file, it should return with RETCODE SCIP_READERR or SCIP_NOFILE.
       */
      virtual SCIP_RETCODE scip_read(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_READER*       reader,             /**< the file reader itself */
         const char*        filename,           /**< full path and name of file to read, or NULL if stdin should be used */
         SCIP_RESULT*       result              /**< pointer to store the result of the file reading call */
         );

   };
   
   /** writes the current problem in sdpa-format to a file*/
   SCIP_RETCODE write_sdpafile(
      SCIP*             scip,                /**< SCIP data structure */
      SdpProblem*       problemdata,         /**< node specific problem data */
      SdpVarMapper*     varmapper            /**<varmapper class object*/
      );
} /* namespace scip */


#endif
