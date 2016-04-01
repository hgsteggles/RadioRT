cmake_minimum_required(VERSION 2.6)
project (radio)

set (PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/src)
set (PROJECT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)

set(RADIO_SRCS
        ${PROJECT_SOURCE_DIR}/main.cpp
        ${PROJECT_SOURCE_DIR}/Bremsstrahlung.cpp
        ${PROJECT_SOURCE_DIR}/Constants.cpp
        ${PROJECT_SOURCE_DIR}/Converter.cpp
        ${PROJECT_SOURCE_DIR}/Fluid.cpp
        ${PROJECT_SOURCE_DIR}/Gaunt.cpp
        ${PROJECT_SOURCE_DIR}/GaussianIntegrator.cpp
        ${PROJECT_SOURCE_DIR}/Logger.cpp
        ${PROJECT_SOURCE_DIR}/ProgressBar.cpp
        ${PROJECT_SOURCE_DIR}/RadioRT.cpp
        ${PROJECT_SOURCE_DIR}/Ray.cpp
        ${PROJECT_SOURCE_DIR}/RayTracer.cpp
        ${PROJECT_SOURCE_DIR}/RecombinationLine.cpp
        ${PROJECT_SOURCE_DIR}/SplineData.cpp
        ${PROJECT_SOURCE_DIR}/Timer.cpp
        ${PROJECT_SOURCE_DIR}/WriteFITS2D.cpp
        ${PROJECT_SOURCE_DIR}/faddeeva/Faddeeva.cpp)

include_directories("${PROJECT_BINARY_DIR}")
include_directories("${PROJECT_INCLUDE_DIR}")

add_executable(${PROJECT_NAME} ${RADIO_SRCS})
target_link_libraries(${PROJECT_NAME} ${LUA_LIBRARIES} cfitsio)