cmake_minimum_required(VERSION 3.0)

set(PROJECT_NAME SCIPSDP)
project(SCIP-SDP)

# set variables
set(EXECUTABLE_NAME "scipsdp")
set(SCIPSDP_VERSION_MAJOR 4)
set(SCIPSDP_VERSION_MINOR 2)
set(SCIPSDP_VERSION_PATCH 0)

set(SRC src)
set(DOC doc)
set(CMAKE cmake)

option(SHARED "Build shared libraries" on)
set(BUILD_SHARED_LIBS ${SHARED})
message(STATUS "Build shared libraries: " ${SHARED})

set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -pedantic -Wno-unused-but-set-variable -Wno-unused-variable")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

find_package(SCIP REQUIRED VARIABLES SHARED)
include_directories(${SCIP_INCLUDE_DIRS})

set(SDPS none CACHE STRING "options for SDP solver")          # create the variable
set_property(CACHE SDPS PROPERTY STRINGS msk sdpa dsdp none)  # define list of values GUI will offer for the variable

set(SYM bliss CACHE STRING "options for symmetry computation")  #create the variable
set_property(CACHE SYM PROPERTY STRINGS bliss none )  #define list of values GUI will offer for the variable

# if we want long long ints for Lapack, e.g., for Matlab
option(LAPACKLONG "use long long int for Lapack" OFF)

#set the correct rpath for OS X
set(CMAKE_MACOSX_RPATH ON)

if(LAPACKLONG)
    add_definitions(-DLAPACKLONG)
endif()

#set defines for Windows
if(WIN32)
    add_definitions(-DNO_SIGACTION)
    add_definitions(-DNO_STRTOK_R)
endif()
if(MSVC)
#    add_definitions(/W4)
    add_definitions(/wd4100)
    add_definitions(/wd4244)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif()

# Visual Studio compiler with static runtime libraries
if(MSVC AND MT)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MT")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MTd")
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} /MT")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} /MTd")
endif()

# create a target for updating the current git hash
file(WRITE ${CMAKE_BINARY_DIR}/scip_update_githash.cmake "
find_program(GIT git)
if(EXISTS \${DST})
   file(STRINGS \${DST} GITHASH_OLD)
   string(REGEX REPLACE \"#define SCIPSDP_GITHASH \\\"(.*)\\\"\" \"\\\\1\" GITHASH_OLD \${GITHASH_OLD})
endif()
if((GIT) AND (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git))
   execute_process(
      COMMAND \${GIT} describe --always --dirty
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      OUTPUT_VARIABLE GITHASH OUTPUT_STRIP_TRAILING_WHITESPACE)
   string(REGEX REPLACE \"^.*-g\" \"\" GITHASH \${GITHASH})
   if(NOT \${GITHASH} STREQUAL \"\${GITHASH_OLD}\")
      file(WRITE \${DST} \"#define SCIPSDP_GITHASH \\\"\${GITHASH}\\\"\n\")
   endif()
else()
   set(GITHASH \${GITHASH_OLD})
endif()
message(STATUS \"Git hash: \" \${GITHASH})
")
add_custom_target(scip_update_githash
                  COMMAND ${CMAKE_COMMAND} -DDST=${PROJECT_SOURCE_DIR}/scipsdpgithash.c
                                           -P ${CMAKE_BINARY_DIR}/scip_update_githash.cmake)

# find lapack and blas
find_package(LAPACK REQUIRED)
if(LAPACK_FOUND)
  message(STATUS "Found lapack library: " ${LAPACK_LIBRARIES})
else()
  message(FATAL_ERROR "Lapack not found")
endif()
find_package(BLAS REQUIRED)
if(BLAS_FOUND)
  message(STATUS "Found blas library: " ${BLAS_LIBRARIES})
else()
  message(FATAL_ERROR "Blas not found")
endif()

# search the selected symmetry computation program
message(STATUS "Finding symmetry computation program \"${SYM}\"")
if(SYM STREQUAL "bliss")
    message(STATUS "Finding BLISS")
    # search for the bliss package delivered with SCIP
    find_library(BLISS_LIBRARY
      NAMES bliss
      HINTS ${SCIP_DIR}
      PATH_SUFFIXES lib)
    if(BLISS_LIBRARY)
        message(STATUS "Finding BLISS - found")
	find_path(BLISS_INCLUDE_DIRS
	  NAMES bliss/abstractgraph.hh
	  HINTS ${SCIP_DIR}
	  PATH_SUFFIXES ../src/bliss/include)
	include_directories(${BLISS_INCLUDE_DIRS})
        set(sym symmetry/compute_symmetry_bliss.cpp)
        set(SYM_LIBRARIES ${BLISS_LIBRARY})
        set(SYM_PIC_LIBRARIES ${BLISS_LIBRARIES})
        if(BLISS_DEFINITIONS STREQUAL "")
            add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:${BLISS_DEFINITIONS}>")
        endif()
    else()
        message(STATUS "Finding BLISS - not found")
        set(sym symmetry/compute_symmetry_none.cpp)
    endif()

    # compute_symmetry_bliss.cpp needsC++-17
    set(CMAKE_CXX_STANDARD 17)
elseif(SYM STREQUAL "none")
    message(STATUS "Support SYM: OFF")
    set(sym symmetry/compute_symmetry_none.cpp)
else()
    message(FATAL_ERROR "option SYM has wrong value")
endif()

# search for the selected SDP solver library
if(SDPS STREQUAL "msk")
    find_package(MOSEK REQUIRED)
elseif(SDPS STREQUAL "sdpa")
    find_package(SDPA REQUIRED)
elseif(SDPS STREQUAL "dsdp")
    find_package(DSDP REQUIRED)
elseif(SDPS STREQUAL "none")
    set(sdpi sdpi/sdpisolver_none.c)
else()
    message(FATAL_ERROR "option SDPS has wrong value")
endif()

if(MOSEK_FOUND)
    include_directories(${MOSEK_INCLUDE_DIRS})
    set(sdpi sdpi/sdpisolver_mosek.c)
    set(SDPS_LIBRARIES ${MOSEK_LIBRARIES})
    set(SDPS_PIC_LIBRARIES ${SDPS_LIBRARIES})
endif()

if(SDPA_FOUND)
    include_directories(${SDPA_INCLUDE_DIRS})
    set(sdpi sdpi/sdpisolver_sdpa.cpp)
    set(SDPS_LIBRARIES ${SDPA_LIBRARIES})
    set(SDPS_PIC_LIBRARIES ${SDPS_LIBRARIES})
endif()

if(DSDP_FOUND)
    include_directories(${DSDP_INCLUDE_DIRS})
    set(sdpi sdpi/sdpisolver_dsdp.c)
    set(SDPS_LIBRARIES ${DSDP_LIBRARIES})
    set(SDPS_PIC_LIBRARIES ${SDPS_LIBRARIES})
endif()

# go to src/ and compile the code
add_subdirectory(src)

include(CTest)

#
# add tests to build the application and run on some easy instances
#
add_test(NAME applications-${EXECUTABLE_NAME}-build
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target ${EXECUTABLE_NAME}
        )

#
# avoid that several build jobs try to concurrently build the SCIP library
# note that this ressource lock name is not the actual libscip target
#
set_tests_properties(applications-${EXECUTABLE_NAME}-build
                    PROPERTIES
                    RESOURCE_LOCK libscipsdp
                    )

#
# instances in the data subdirectory
#
set(instances
    example_small.dat-s
    example_small_cbf.cbf
    example_inf.dat-s
    example_TT.dat-s.gz
    example_CLS.dat-s.gz
    example_MkP.dat-s.gz
)

#
# loop over instances and define a test
#
foreach(instance ${instances})
  add_test(NAME ${EXECUTABLE_NAME}-${instance}
    COMMAND $<TARGET_FILE:${EXECUTABLE_NAME}> -f ${CMAKE_CURRENT_SOURCE_DIR}/instances/${instance}
    )
  set_tests_properties(${EXECUTABLE_NAME}-${instance}
    PROPERTIES
    PASS_REGULAR_EXPRESSION "SCIP Status        : problem is solved"
    DEPENDS applications-${EXECUTABLE_NAME}-build
    )
endforeach()
