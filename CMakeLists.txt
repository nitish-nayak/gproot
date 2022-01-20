# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(gproot CXX C)

# source code
include_directories("${PROJECT_SOURCE_DIR}/")
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cxx)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.h)

# ROOT
find_program(ROOTCONF "root-config")
if(ROOTCONF)
  message("-- Detecting ROOT:    found at ${ROOTCONF}")
else()
  message(FATAL_ERROR "-- Detecting ROOT:    not found")
endif()
execute_process(COMMAND root-config --prefix OUTPUT_VARIABLE ROOT_PREFIX  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CFLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND root-config --libs   OUTPUT_VARIABLE ROOT_LINK    OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x ${ROOT_CFLAGS}")
add_library(gproot SHARED ${sources})
target_link_libraries(gproot ${ROOT_LINK})

message(${CMAKE_PROJECT_NAME} " will be installed to " ${CMAKE_INSTALL_PREFIX})

install(TARGETS gproot DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

file(GLOB dist_headers ${PROJECT_SOURCE_DIR}/*.h)
install(FILES ${dist_headers} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${CMAKE_PROJECT_NAME}/)

