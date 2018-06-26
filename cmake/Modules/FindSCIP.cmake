# - Try to find the SCIP library
# See http://scip.zib.de/ for more information on SCIP
#
# Variables used by this module:
#
# SCIP_ROOT_DIR
#   Set this variable to the SCIP source or install path.
#   Otherwise, default paths are searched, e.g. /usr/local/
# SCIP_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds.
#   The default is "opt".
# SCIP_LPS
#   Set this variable to "spx" (="spx2") or "spx1" (others not implemented, yet).
#   The default is "spx".
# SCIP_LPS_BUILD
#   Set this variable to "opt", "dbg" or "prf" for corresponding builds of the LP solver.
#   The default is "opt".
#   Note that this value is ignored if the LP solver was searched for already,
#   e.g., if find_package(SoPlex) was called before calling find_package(SCIP).
# SCIP_OSTYPE
#   Set this variable to any of SCIP's OS types, e.g. "linux", "win", etc.
#   The default is determined from `uname -s`.
# SCIP_ARCH
#   Set this variable to any of SCIP's architecture types, e.g. "x86", "x86_64", etc.
#   The default is determined from `uname -m`.
# SCIP_COMP
#   Set this variable to any of SCIP's compiler types, e.g. "gnu", "intel", etc.
#   The default is "gnu".
#
# Once done, this will define
#
#  SCIP_INCLUDE_DIRS   - where to find scip.h, etc.
#  SCIP_LIBRARIES      - list of libraries when using SCIP.
#  SCIP_FOUND          - true if SCIP was found.
#
# Author:
# Matthias Walter <matthias.walter@ovgu.de>
#
# Distributed under the Boost Software License, Version 1.0.
# (See http://www.boost.org/LICENSE_1_0.txt)

find_package(ZLIB REQUIRED)
find_package(Readline REQUIRED)
find_package(GMP REQUIRED)

find_package(ZIMPL QUIET)

# Hints and paths for the search
set(_SCIP_ROOT_HINTS $ENV{SCIP_ROOT_DIR} ${SCIP_ROOT_DIR})
set(_SCIP_ROOT_PATHS $ENV{SCIP_ROOT_DIR} ${SCIP_ROOT_DIR})

# Read SCIP_BUILD from (environment) variable.
if (${SCIP_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_BUILD ${SCIP_BUILD})
elseif ($ENV{SCIP_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_BUILD $ENV{SCIP_BUILD})
else()
  set(_SCIP_BUILD "opt")
endif()

# Read SCIP_LPS from (environment) variable.
if (${SCIP_LPS} MATCHES "^(spx|spx2|spx1|cpx|msk)$")
  set(_SCIP_LPS ${SCIP_LPS})
elseif ($ENV{SCIP_LPS} MATCHES "^(spx|spx2|spx1|cpx|msk)$")
  set(_SCIP_LPS $ENV{SCIP_LPS})
else()
  set(_SCIP_LPS "none")
endif()

# Read SCIP_LPS_BUILD from (environment) variable.
if (${SCIP_LPS_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_LPS_BUILD ${SCIP_LPS_BUILD})
elseif ($ENV{SCIP_LPS_BUILD} MATCHES "^(opt|dbg|prf)$")
  set(_SCIP_LPS_BUILD $ENV{SCIP_LPS_BUILD})
else()
  set(_SCIP_LPS_BUILD "opt")
endif()

# Note: To see how SCIP determines OSTYPE and ARCH, look at scip/make/make.detecthost.

# Read SCIP_OSTYPE from (environment) variable or from `uname -s`
if (${SCIP_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SCIP_OSTYPE ${SCIP_OSTYPE})
elseif ($ENV{SCIP_OSTYPE} MATCHES "^(aix|cygwin|darwin|freebsd|hp-ux|irix|linux|mingw|osf1|sunos|win)$")
  set(_SCIP_OSTYPE $ENV{SCIP_OSTYPE})
else()
  execute_process(COMMAND uname -s OUTPUT_VARIABLE _SCIP_OSTYPE OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(TOLOWER ${_SCIP_OSTYPE} _SCIP_OSTYPE)
  string(REGEX REPLACE "cygwin.*" "cygwin" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "irix.." "irix" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "windows.*" "windows" _SCIP_OSTYPE ${_SCIP_OSTYPE})
  string(REGEX REPLACE "mingw.*" "mingw" _SCIP_OSTYPE ${_SCIP_OSTYPE})
endif()

# Read SCIP_ARCH from (environment) variable or from `uname -m`
if (${SCIP_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SCIP_ARCH ${SCIP_ARCH})
elseif ($ENV{SCIP_ARCH} MATCHES "^(alpha|arm|clang|gnu|hppa|intel|mips|ppc|pwr4|sparc|x86|x86_64)$")
  set(_SCIP_ARCH $ENV{SCIP_ARCH})
else()
  execute_process(COMMAND uname -m OUTPUT_VARIABLE _SCIP_ARCH OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE "sun.." "sparc" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "i.86" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "i86pc" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "[0-9]86" "x86" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "amd64" "x86_64" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "IP.." "mips" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "9000...." "hppa" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "Power\ Macintosh" "ppc" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "00.........." "pwr4" _SCIP_ARCH ${_SCIP_ARCH})
  string(REGEX REPLACE "arm.*" "arm" _SCIP_ARCH ${_SCIP_ARCH})
endif()

# Read SCIP_COMP from (environment) variable.
if (${SCIP_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SCIP_COMP ${SCIP_COMP})
elseif ($ENV{SCIP_COMP} MATCHES "^(clang|compaq|gnu|hp|ibm|insure|intel|msv|purify|sgi|sun)$")
  set(_SCIP_COMP $ENV{SCIP_COMP})
else()
  set(_SCIP_COMP "gnu")
endif()

# Root

find_path(SCIP_ROOT_DIR NAMES include/scip/solve.h src/scip/solve.h HINTS ${_SCIP_ROOT_HINTS} PATHS ${_SCIP_ROOT_HINTS})

# Includes

find_path(_SCIP_INCLUDE NAMES scip/solve.h PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES include src)
if (_SCIP_INCLUDE)
  set(SCIP_INCLUDE_DIRS ${_SCIP_INCLUDE})

  # Extract version from scip/def.h
  file(STRINGS "${_SCIP_INCLUDE}/scip/def.h" _SCIP_VERSION_STR REGEX "^#define[\t ]+SCIP_VERSION[\t ]+[0-9][0-9][0-9].*")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+([0-9]).*$" "\\1" SCIP_VERSION_MAJOR "${_SCIP_VERSION_STR}")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+[0-9]([0-9]).*$" "\\1" SCIP_VERSION_MINOR "${_SCIP_VERSION_STR}")
  string(REGEX REPLACE "^.*SCIP_VERSION[\t ]+[0-9][0-9]([0-9]).*$" "\\1" SCIP_VERSION_PATCH "${_SCIP_VERSION_STR}")
  file(STRINGS "${_SCIP_INCLUDE}/scip/def.h" _SCIP_SUBVERSION_STR REGEX "^#define[\t ]+SCIP_SUBVERSION[\t ]+[0-9].*")
  string(REGEX REPLACE "^.*SCIP_SUBVERSION[\t ]+([0-9]).*$" "\\1" SCIP_VERSION_SUBVERSION "${_SCIP_SUBVERSION_STR}")
  set(SCIP_VERSION_STRING "${SCIP_VERSION_MAJOR}.${SCIP_VERSION_MINOR}.${SCIP_VERSION_PATCH}.${_SCIP_OSTYPE}.${_SCIP_ARCH}.${_SCIP_COMP}.${_SCIP_BUILD}")

  set(_SCIP_FOUND_ALL TRUE)

  # Check for ZIMPL.
  #if (NOT ZIMPL_FOUND)
  #  set(_SCIP_FOUND_ALL FALSE)
  #  if (NOT SCIP_FIND_QUIETLY)
  #    message(STATUS "SCIP dependency ZIMPL was not found.")
  #  endif()
  #endif()

  if(SHARED)
    set(SCIP_PATH_SUFFIX lib/shared)
  else()
    set(SCIP_PATH_SUFFIX lib/static)
  endif()

  # Search for libscip corresponding to version.
  find_library(_SCIP_LIB_SCIP NAMES "scip-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})

  if (NOT ${_SCIP_LIB_SCIP} MATCHES "scip")
    set(_SCIP_FOUND_ALL FALSE)
    if (NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libscip-${SCIP_VERSION_STRING} was not found. Search paths: ${_SCIP_ROOT_PATHS}")

      # Check if some other version was found and report to user.
      find_library(_SCIP_LIB_SCIP_TEST NAMES "scip" ${SCIP_ROOT_DIR} PATH_SUFFIXES lib)
      if (${_SCIP_LIB_SCIP_TEST} MATCHES "scip")
        message(STATUS "Found a SCIP library different from the one promised by ${SCIP_INCLUDE_DIRS}/scip/def.h")
      endif()
    endif()
  endif()

  # Search for libobjscip
  find_library(_SCIP_LIB_OBJSCIP NAMES "objscip-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})
  if (NOT ${_SCIP_LIB_OBJSCIP} MATCHES "objscip")
    set(_SCIP_FOUND_ALL FALSE)
    if (NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libobjscip-${SCIP_VERSION_STRING} was not found.")
    endif()
  endif()

  # Search for nlpi. TODO: cppad is currently hard-coded, while ipopt is not recognized.
  find_library(_SCIP_LIB_NLPI NAMES "nlpi.cppad-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})
  if (NOT ${_SCIP_LIB_NLPI} MATCHES "nlpi")
    set(_SCIP_FOUND_ALL FALSE)
    if (NOT SCIP_FIND_QUIETLY)
      message(STATUS "SCIP library libnlpi.cppad-${SCIP_VERSION_STRING} was not found.")
    endif()
  endif()

  # Search for the LP solver: spx
  if (${_SCIP_LPS} STREQUAL "spx")
    # Search for liblpispx
    find_library(_SCIP_LIB_LPI NAMES "lpispx-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})
    if (NOT ${_SCIP_LIB_LPI} MATCHES "spx")
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()

    # Search for SoPlex.
    if (NOT SOPLEX_FOUND)
      set(SOPLEX_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SOPLEX_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()

  # Search for the LP solver: spx2
  if (${_SCIP_LPS} STREQUAL "spx2")
    # Search for liblpispx2
    find_library(_SCIP_LIB_LPI NAMES "lpispx2-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})
    if (NOT ${_SCIP_LIB_LPI} MATCHES "spx2")
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx2-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()

    # Search for SoPlex.
    if (NOT SOPLEX_FOUND)
      set(SOPLEX_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SOPLEX_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()

  # Search for the LP solver: spx1
  if (${_SCIP_LPS} STREQUAL "spx1")
    # Search for liblpispx1
    find_library(_SCIP_LIB_LPI NAMES "lpispx1-${SCIP_VERSION_STRING}" PATHS ${SCIP_ROOT_DIR} PATH_SUFFIXES ${SCIP_PATH_SUFFIX})
    if (NOT ${_SCIP_LIB_LPI} MATCHES "liblpispx1")
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP library liblpispx1-${SCIP_VERSION_STRING} was not found.")
      endif()
    endif()

    # Search for SoPlex.
    if (NOT SOPLEX_FOUND)
      set(SOPLEX_BUILD ${SCIP_LPS_BUILD})
      find_package(SoPlex REQUIRED)
    endif()
    if (SOPLEX_FOUND)
      set(_SCIP_LIB_LPSOLVER ${SOPLEX_LIBRARIES})
    else()
      set(_SCIP_FOUND_ALL FALSE)
      if (NOT SCIP_FIND_QUIETLY)
        message(STATUS "SCIP dependency SoPlex was not found.")
      endif()
    endif()
  endif()

  if (_SCIP_FOUND_ALL)
    set(SCIP_LIBRARIES
      ${_SCIP_LIB_OBJSCIP} ${_SCIP_LIB_SCIP} ${_SCIP_LIB_LPI} ${_SCIP_LIB_NLPI} ${_SCIP_LIB_LPSOLVER}
      ${ZLIB_LIBRARIES} ${Readline_LIBRARY} ${GMP_LIBRARIES}
    )
    if (ZIMPL_FOUND)
      set(SCIP_LIBRARIES ${SCIP_LIBRARIES} ${ZIMPL_LIBRARIES})
    endif()
  endif()
endif()

# Let cmake process everything.
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCIP REQUIRED_VARS SCIP_ROOT_DIR SCIP_INCLUDE_DIRS SCIP_LIBRARIES VERSION_VAR SCIP_VERSION_STRING)
