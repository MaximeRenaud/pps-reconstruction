cmake_minimum_required(VERSION 2.6)
project(PPS_TB_RECO)

set(GCC_COMPILE_FLAGS "-Wall -fPIC -O2 -g -lsqlite3")
set(GCC_LINK_FLAGS "-lsqlite3")
add_definitions(${GCC_COMPILE_FLAGS})

# File reader
file(GLOB reader_sources ${PROJECT_SOURCE_DIR}/src/FileReader.cpp)
add_library(reader_lib OBJECT ${reader_sources})

#file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
#list(REMOVE_ITEM sources ${reader_sources})
#add_library(src_lib OBJECT ${sources})
#set_property(TARGET src_lib PROPERTY LINK_FLAGS "-lsqlite3")

include_directories("${PROJECT_SOURCE_DIR}/include")

# Here have tests
add_subdirectory(test EXCLUDE_FROM_ALL)

