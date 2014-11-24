cmake_minimum_required(VERSION 2.8)
project(cpplex)

message("Running cmake for cpplex subdirectory")
message("cpplex: Copyright (c) 2010-2013 Tommaso Urli")
set(CMAKE_COLOR_MAKEFILE ON)

# add flags depending on build type
set(CMAKE_BUILD_TYPE Release) 
#set(CMAKE_BUILD_TYPE Debug) 
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++0x -Wall -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++0x -Wall -fopenmp")

file(GLOB SPLX_SRC 
	"simplex/src/*.cc")
file(GLOB SPLX_HDR
	"simplex/include/*.h")
file(GLOB PILAL_SRC 
	"pilal/src/*.cc")
file(GLOB PILAL_HDR
	"pilal/include/*.h")

include_directories(
  ./pilal/include
  ./simplex/include
  /usr/include
  /usr/local/include
  )
link_directories(
  /usr/lib
  /usr/local/lib/
  )

add_library(cpplex-simplex SHARED 
	${SPLX_SRC})
add_library(cpplex-pilal SHARED 
	${PILAL_SRC})


install(TARGETS cpplex-simplex cpplex-pilal
  LIBRARY DESTINATION lib
  )
install(FILES ${PILAL_HDR} ${SPLX_HDR}
	DESTINATION include/cpplex)
