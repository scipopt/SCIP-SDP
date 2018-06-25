find_path(SCIP_INCLUDE_DIRS
    NAMES scip/solve.h
    HINTS ${SCIP_DIR} $ENV{SCIP_DIR}
    PATH_SUFFIXES src)

find_library(SCIP_LIBRARY
    NAMES scip
    HINTS ${SCIP_DIR} $ENV{SCIP_DIR}
    PATH_SUFFIXES lib/static)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCIP DEFAULT_MSG SCIP_INCLUDE_DIRS SCIP_LIBRARY)
