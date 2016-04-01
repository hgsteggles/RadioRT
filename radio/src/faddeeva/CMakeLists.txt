cmake_minimum_required(VERSION 2.6)
project (faddeeva)

set (PROJECT_INCLUDE_DIR ../../..)
set (PROJECT_SOURCE_DIR .)

set(FADDEEVA_SRCS
        ${PROJECT_SOURCE_DIR}/Faddeeva.cpp
        )

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_INCLUDE_DIR}")

add_library(${PROJECT_NAME} STATIC ${FADDEEVA_SRCS})