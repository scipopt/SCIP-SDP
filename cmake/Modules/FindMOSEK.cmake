find_path(MOSEK_INCLUDE_DIRS
    NAMES mosek.h
    HINTS ${MOSEK_DIR} $ENV{MOSEK_DIR}
    PATH_SUFFIXES h)

find_library(MOSEK_LIBRARY
    NAMES mosek64
    HINTS ${MOSEK_DIR} $ENV{MOSEK_DIR}
    PATH_SUFFIXES bin)

# this library is not used anymore starting with Mosek 9
find_library(IOMP5_LIBRARY
    NAMES iomp5
    HINTS ${MOSEK_DIR} $ENV{MOSEK_DIR}
    PATH_SUFFIXES bin)

if(IOMPS_LIBRARY)
  set(MOSEK_LIBRARIES ${MOSEK_LIBRARY} ${IOMP5_LIBRARY} -pthread)
else()
  # if libiomps is not available, we skip it
  set(MOSEK_LIBRARIES ${MOSEK_LIBRARY} -pthread)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MOSEK DEFAULT_MSG MOSEK_INCLUDE_DIRS MOSEK_LIBRARIES)
