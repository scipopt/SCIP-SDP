# - Try to find SoPlex LP solver library
# See http://soplex.zib.de/ for more information on SoPlex
#
# Variables used by this module:
#
# SOPLEX_ROOT_DIR
#   Set this variable to the SoPlex source or install path.
#   Otherwise, default paths are searched, e.g. /usr/local/
# SOPLEX_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds.
#   The default is "opt".
# SOPLEX_OSTYPE
#   Set this variable to any of SoPlex's OS types, e.g. "linux", "win", etc.
#   The default is determined from `uname -s`.
# SOPLEX_ARCH
#   Set this variable to any of SoPlex's architecture types, e.g. "x86", "x86_64", etc.
#   The default is determined from `uname -m`.
# SOPLEX_COMP
#   Set this variable to any of SoPlex's compiler types, e.g. "gnu", "intel", etc.
#   The default is "gnu".
#
# Once done, this will define
#
#  SOPLEX_INCLUDE_DIRS   - where to find soplex.h, etc.
#  SOPLEX_LIBRARIES      - list of libraries when using SoPlex.
#  SOPLEX_FOUND          - true if SoPlex was found.
#
# Author:
# Matthias Walter <matthias.walter@ovgu.de>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

find_package(ZLIB REQUIRED)

# Hints and paths for the search
set(_SOPLEX_ROOT_HINTS $ENV{SOPLEX_ROOT_DIR} ${SOPLEX_ROOT_DIR})
set(_SOPLEX_ROOT_PATHS $ENV{SOPLEX_ROOT_DIR} ${SOPLEX_ROOT_DIR})

# Read SOPLEX_BUILD from (environment) variable.
if (${SOPLEX_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SOPLEX_BUILD ${SOPLEX_BUILD})
elseif ($ENV{SOPLEX_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SOPLEX_BUILD $ENV{SOPLEX_BUILD})
else()
  set(_SOPLEX_BUILD "opt")
endif()

# Note: To see how SoPlex determines OSTYPE and ARCH, look at soplex/make/make.detecthost.

# Read SOPLEX_OSTYPE from (environment) variable or from `uname -s`
if (${SOPLEX_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SOPLEX_OSTYPE ${SOPLEX_OSTYPE})
elseif ($ENV{SOPLEX_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SOPLEX_OSTYPE $ENV{SOPLEX_OSTYPE})
else()
  execute_process(COMMAND uname -s OUTPUT_VARIABLE _SOPLEX_OSTYPE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(TOLOWER ${_SOPLEX_OSTYPE} _SOPLEX_OSTYPE)
  string(REGEX REPLACE "cygwin.*" "cygwin" _SOPLEX_OSTYPE ${_SOPLEX_OSTYPE})
  string(REGEX REPLACE "irix.." "irix" _SOPLEX_OSTYPE ${_SOPLEX_OSTYPE})
  string(REGEX REPLACE "windows.*" "windows" _SOPLEX_OSTYPE ${_SOPLEX_OSTYPE})
  string(REGEX REPLACE "mingw.*" "mingw" _SOPLEX_OSTYPE ${_SOPLEX_OSTYPE})
endif()

# Read SOPLEX_ARCH from (environment) variable or from `uname -m`
if (${SOPLEX_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SOPLEX_ARCH ${SOPLEX_ARCH})
elseif ($ENV{SOPLEX_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SOPLEX_ARCH $ENV{SOPLEX_ARCH})
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE _SOPLEX_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "sun.." "sparc" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "i.86" "x86" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "i86pc" "x86" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "[0-9]86" "x86" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "amd64" "x86_64" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "IP.." "mips" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "9000...." "hppa" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "Power\ Macintosh" "ppc" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "00.........." "pwr4" _SOPLEX_ARCH ${_SOPLEX_ARCH})
  string(REGEX REPLACE "arm.*" "arm" _SOPLEX_ARCH ${_SOPLEX_ARCH})
endif()

# Read SOPLEX_COMP from (environment) variable.
if (${SOPLEX_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SOPLEX_COMP ${SOPLEX_COMP})
elseif ($ENV{SOPLEX_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SOPLEX_COMP $ENV{SOPLEX_COMP})
else()
  set(_SOPLEX_COMP "gnu")
endif()

# Root

find_path(SOPLEX_ROOT_DIR NAMES include/soplex.h soplex.h HINTS ${_SOPLEX_ROOT_HINTS} PATHS ${_SOPLEX_ROOT_PATHS})

# Includes

find_path(_SOPLEX_INCLUDE NAMES soplex.h PATHS ${SOPLEX_ROOT_DIR} PATH_SUFFIXES include src)
if (_SOPLEX_INCLUDE)
  set(SOPLEX_INCLUDE_DIRS ${_SOPLEX_INCLUDE})
  
  # Extract version from spxdefines.h
  file(STRINGS "${_SOPLEX_INCLUDE}/spxdefines.h" _SOPLEX_VERSION_STR REGEX "^#define[\t ]+SOPLEX_VERSION[\t ]+[0-9][0-9][0-9].*")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+([0-9]).*$" "\\1" SOPLEX_VERSION_MAJOR "${_SOPLEX_VERSION_STR}")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+[0-9]([0-9]).*$" "\\1" SOPLEX_VERSION_MINOR "${_SOPLEX_VERSION_STR}")
  string(REGEX REPLACE "^.*SOPLEX_VERSION[\t ]+[0-9][0-9]([0-9]).*$" "\\1" SOPLEX_VERSION_PATCH "${_SOPLEX_VERSION_STR}")
  file(STRINGS "${_SOPLEX_INCLUDE}/spxdefines.h" _SOPLEX_SUBVERSION_STR REGEX "^#define[\t ]+SOPLEX_SUBVERSION[\t ]+[0-9].*")
  string(REGEX REPLACE "^.*SOPLEX_SUBVERSION[\t ]+([0-9]).*$" "\\1" SOPLEX_VERSION_SUBVERSION "${_SOPLEX_SUBVERSION_STR}")
  set(SOPLEX_VERSION_STRING "${SOPLEX_VERSION_MAJOR}.${SOPLEX_VERSION_MINOR}.${SOPLEX_VERSION_PATCH}.${SOPLEX_VERSION_SUBVERSION}.${_SOPLEX_OSTYPE}.${_SOPLEX_ARCH}.${_SOPLEX_COMP}.${_SOPLEX_BUILD}")
 
  # Search for library corresponding to version.
  find_library(_SOPLEX_LIB NAMES "soplex-${SOPLEX_VERSION_STRING}" PATHS ${_SOPLEX_ROOT_PATHS} PATH_SUFFIXES lib)
  if (_SOPLEX_LIB)
    set(SOPLEX_LIBRARIES ${_SOPLEX_LIB} ${ZLIB_LIBRARIES})
  endif()
endif()

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SoPlex REQUIRED_VARS SOPLEX_ROOT_DIR SOPLEX_INCLUDE_DIRS SOPLEX_LIBRARIES VERSION_VAR SOPLEX_VERSION_STRING)

