cmake_minimum_required(VERSION 2.6)

project(dmeans)

#first compile cpplex for the examples
add_subdirectory(cpplex)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Debug) # Release; debug; ...
#set(CMAKE_BUILD_TYPE Release) # Release; debug; ...

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++0x -Wall -fopenmp")
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++0x -Wall -fopenmp")


#find eigen using a module in cmake/
find_package (Eigen3 REQUIRED)
message(${EIGEN3_INCLUDE_DIR})


include_directories(
  ./include
  ${EIGEN3_INCLUDE_DIR}
  /usr/include
  /usr/local/include/
  /usr/local/include/eigen3
  example/cpplex/include
  )
link_directories(
  /usr/lib
  /usr/local/lib/
  )


add_executable(DMeansExample example/maindm.cpp)

target_link_libraries(testVb 
  ampsBnps 
  gsl gslcblas)

INSTALL(TARGETS ampsBnps
  #  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  # ARCHIVE DESTINATION lib
  COMPONENT compLib
  )

## libraries
#add_library(ampsBnps SHARED
#  ./src/vbfuncs_gaussian.c
#  ./src/vb.c
#  ./src/updates.c
#  ./src/costfcn.c
#  ./src/amps.c
#  ./src/matching.c
#  )



install(DIRECTORY include/ DESTINATION include/
  FILES_MATCHING PATTERN "*.h")
install(DIRECTORY cuda/ DESTINATION include/
  FILES_MATCHING PATTERN "*.h")
