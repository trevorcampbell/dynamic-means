project(DMeans)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_COLOR_MAKEFILE ON)
set(CMAKE_BUILD_TYPE Debug) # Release; debug; ...
#set(CMAKE_BUILD_TYPE Release) # Release; debug; ...

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_RELEASE "-g -std=c++0x -Wall -fopenmp -Wno-deprecated -finline-functions -DNDEBUG -DCUDA") # -DNDEBUG
set(CMAKE_CXX_FLAGS_DEBUG "-g -std=c++0x -Wall -Wno-deprecated  -DCUDA")


#find eigen
find_package (Eigen3 REQUIRED)
message(${EIGEN3_INCLUDE_DIR})


include_directories(
  ./include
  ./cuda
  ${EIGEN3_INCLUDE_DIR}
  #${EIGEN_INCLUDE_DIR}
  ${GSL_INCLUDE_DIRS}
  /usr/include
  /usr/local/include/eigen3
  )
link_directories(
  $ENV{WORKSPACE_HOME}/3rdparty/lib
  ${GSL_LIBRARY_DIRS}
  /usr/lib
  /usr/local/cuda-5.5/targets/x86_64-linux/lib/
  ./cuda/
  /usr/local/lib/
  )

# libraries
add_library(ampsBnps SHARED
  ./src/vbfuncs_gaussian.c
  ./src/vb.c
  ./src/updates.c
  ./src/costfcn.c
  ./src/amps.c
  ./src/matching.c
  )


add_executable(DMeans example/maindm.cpp)
target_link_libraries(testVb 
  ampsBnps 
  gsl gslcblas)

INSTALL(TARGETS ampsBnps
  #  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  # ARCHIVE DESTINATION lib
  COMPONENT compLib
  )



install(DIRECTORY include/ DESTINATION include/                                 
  FILES_MATCHING PATTERN "*.h")                                         
install(DIRECTORY cuda/ DESTINATION include/                                 
  FILES_MATCHING PATTERN "*.h")                                         
