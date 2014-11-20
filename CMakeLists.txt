cmake_minimum_required(VERSION 2.6)

project(dmeans)

set(CMAKE_COLOR_MAKEFILE ON)

#set build type
#set(CMAKE_BUILD_TYPE Debug) 
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++0x -Wall -fopenmp")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++0x -Wall -fopenmp")

#first compile cpplex for the examples
add_subdirectory(cpplex)

#install the header library
file(GLOB_RECURSE HDRS "dmeans/" 
	"*.hpp")

#find eigen using a module in cmake/
find_package (Eigen3 REQUIRED)
message(${EIGEN3_INCLUDE_DIR})

include_directories(
  ${EIGEN3_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}
  /usr/include
  /usr/local/include/
  /usr/local/include/cpplex
  /usr/local/include/dmeans
  )

link_directories(
  /usr/lib
  /usr/local/lib/
  )

add_executable(dmeans_example example/maindm.cpp)

add_dependencies(dmeans_example headers)

target_link_libraries(dmeans_example
	cpplex-pilal
	cpplex-simplex
  )

install(FILES ${HDRS}
	DESTINATION include/dmeans)

