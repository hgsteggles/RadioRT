cmake_minimum_required(VERSION 2.6)

set(RADIO_SRCS
        ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Bremsstrahlung.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Constants.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Converter.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Fluid.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Gaunt.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/GaussianIntegrator.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Logger.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ProgressBar.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/RadioRT.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Ray.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/RayTracer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/RecombinationLine.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/SplineData.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/Timer.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/WriteFITS2D.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/faddeeva/Faddeeva.cpp)

include_directories("${RADIO_SOURCE_DIR}/lib")
include_directories("${RADIO_SOURCE_DIR}/lib/lua-5.2.3")

add_executable(radio ${RADIO_SRCS})
#target_link_libraries(radio ${LUA_LIBRARIES} cfitsio)
target_link_libraries(radio ${RADIO_SOURCE_DIR}/lib/liblua.a dl cfitsio)
