# - Try to find cpplex lib
#
# Once done this will define
#
#  CPPLEX_FOUND - system has cpplex lib 
#  CPPLEX_INCLUDE_DIR - the cpplex include directory

# Copyright (c) 2014 Trevor Campbell <tdjc@mit.edu>
# Redistribution and use is allowed 

find_library(CPPLEX_LIBS NAMES cpplex-pilal cpplex-simplex
	PATHS
    "/usr/local/lib"
    "/usr/lib"
    ${CMAKE_INSTALL_PREFIX}/lib
    PATH_SUFFIXES cpplex
  )

find_path(CPPLEX_INCLUDE_DIR NAMES simplex.h pilal.h
    PATHS
    "/usr/local/include"
    "/usr/include"
    ${CMAKE_INSTALL_PREFIX}/include
    PATH_SUFFIXES cpplex
  )

if(NOT CPPLEX_INCLUDE_DIR OR NOT CPPLEX_LIBS)
	message("")
	message("----------------------------------------------------------------------")
	message("You have to build/install cpplex before trying to build/install dmeans")
	message("The cpplex library is provided with this project, located in:")
	message("${PROJECT_SOURCE_DIR}/cpplex")
	message("----------------------------------------------------------------------")
	message("")
endif(NOT CPPLEX_INCLUDE_DIR OR NOT CPPLEX_LIBS)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Cpplex DEFAULT_MSG CPPLEX_INCLUDE_DIR CPPLEX_LIBS)

mark_as_advanced(CPPLEX_INCLUDE_DIR)
mark_as_advanced(CPPLEX_LIBS)


