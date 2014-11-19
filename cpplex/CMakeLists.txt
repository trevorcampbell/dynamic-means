project(DMeans)
cmake_minimum_required(VERSION 2.6)

message("Running cmake for cpplex subdirectory")
message("cpplex: Copyright (c) 2010-2013 Tommaso Urli")

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
#add_executable(solver 
#	solver/solver.cc)
#target_link_libraries(solver 
#	simplex 
#	pilal)

install(TARGETS cpplex-simplex cpplex-pilal
  LIBRARY DESTINATION lib
  )
install(FILES ${PILAL_HDR} ${SPLX_HDR}
	DESTINATION include/cpplex)
