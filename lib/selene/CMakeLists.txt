cmake_minimum_required(VERSION 3.0)
project (selene)

find_package(Lua REQUIRED)

message(WARNING ${LUA_VERSION_MAJOR} "." ${LUA_VERSION_MINOR} " " ${LUA_INCLUDE_DIR})
include_directories(${LUA_INCLUDE_DIR})

set (PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR})
set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_INCLUDE_DIR}")