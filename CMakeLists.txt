cmake_minimum_required(VERSION 2.8)
project(dmeans)
#include(ExternalProject)

set(CMAKE_COLOR_MAKEFILE ON)

#set build type
#set(CMAKE_BUILD_TYPE Debug) 
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++0x -Wall -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++0x -Wall -fopenmp")

#first compile cpplex for the examples
#ExternalProject_Add(cpplex,
#	PREFIX ${PROJECT_SOURCE_DIR}/cpplex
#	SOURCE_DIR ${PROJECT_SOURCE_DIR}/cpplex
#	)


#find eigen using a module in cmake/
find_package (Eigen3 REQUIRED)
find_package (Cpplex REQUIRED)
message("Eigen3 found: ${EIGEN3_INCLUDE_DIR}")
message("Cpplex found: ${CPPLEX_INCLUDE_DIR}")

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/cpplex/pilal/include
  ${PROJECT_SOURCE_DIR}/cpplex/simplex/include
  /usr/include
  /usr/local/include/
  )

link_directories(
  /usr/lib
  /usr/local/lib/
  )

add_executable(dmeans_example example/maindm.cpp)
add_executable(sdmeans_example example/mainsdm.cpp)
add_executable(kdmeans_example example/mainkdm.cpp)
add_executable(all_example example/mainall.cpp)

target_link_libraries(dmeans_example
	cpplex-pilal
	cpplex-simplex
  )
target_link_libraries(sdmeans_example
	cpplex-pilal
	cpplex-simplex
  )
target_link_libraries(kdmeans_example
	cpplex-pilal
	cpplex-simplex
  )
target_link_libraries(all_example
   cpplex-pilal
   cpplex-simplex
)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/dmeans
	DESTINATION include)

