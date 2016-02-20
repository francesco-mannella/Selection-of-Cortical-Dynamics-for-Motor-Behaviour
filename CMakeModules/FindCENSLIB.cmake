# - Try to find CENSLIB
# Once done this will define
#  LIBCENS_FOUND - System has CENSLIB
#  LIBCENS_INCLUDE_DIRS - The CENSLIB include directories
#  LIBCENS_LIBRARIES - The libraries needed to use CENSLIB
#  LIBCENS_DEFINITIONS - Compiler switches required for using CENSLIB

set(LIBCENS_DEFINITIONS "-std=c++11 ")

find_path(LIBCENS_INCLUDE_DIR CENS/cens_physics.h
    PATH_SUFFIXES libCENS )

find_library(LIBCENS_LIBRARY NAMES CENS libCENS )

set(LIBCENS_LIBRARIES ${LIBCENS_LIBRARY} )
set(LIBCENS_INCLUDES ${LIBCENS_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LIBCENS_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CENSLIB  DEFAULT_MSG
                                  LIBCENS_LIBRARY LIBCENS_INCLUDE_DIR)

mark_as_advanced(LIBCENS_INCLUDE_DIR LIBCENS_LIBRARY )
