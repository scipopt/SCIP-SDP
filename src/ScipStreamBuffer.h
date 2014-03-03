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

/**@file   ScipStreamBuffer.cpp
 * @brief  An std::streambuf that uses SCIP I/O routines (suitable for reading)
 * @author Lars Schewe
 */

#ifndef SCIPSTREAM_H
#define SCIPSTREAM_H

#include <streambuf>                    // for streamsize, streambuf

#include <cstddef>                      // for size_t
#include "scip/pub_fileio.h"            // for SCIP_FILE
#include "scip/type_scip.h"             // for SCIP

class ScipStreamBuffer : public std::streambuf
{
 public:
   ScipStreamBuffer(SCIP* scip, SCIP_FILE* file, bool close_on_exit);

   ~ScipStreamBuffer();
   
 private:
   /// the underflow function is responsible for the refilling of the buffer
   virtual int underflow();

   virtual std::streamsize xsgetn(char *dest, std::streamsize request);
   
   SCIP* scip_;
   SCIP_FILE* file_;
   char* g_buffer_; /// pointer to the get-buffer
   size_t g_buffer_size_; /// size of the get-buffer
   bool close_on_exit_;
};
#endif
