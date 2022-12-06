find_path(DSDP_INCLUDE_DIRS
    NAMES dsdp5.h
    HINTS ${DSDP_DIR} $ENV{DSDP_DIR}
    PATH_SUFFIXES include)

find_library(DSDP_LIBRARY
    NAMES dsdp
    HINTS ${DSDP_DIR} $ENV{DSDP_DIR}
    PATH_SUFFIXES lib)

set(DSDP_LIBRARIES ${DSDP_LIBRARY} -llapack -lblas)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(DSDP DEFAULT_MSG DSDP_INCLUDE_DIRS DSDP_LIBRARIES)
