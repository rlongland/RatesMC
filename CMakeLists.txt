cmake_minimum_required(VERSION 2.4)
project(RatesMC2)

file(GLOB CPP_FILES src/*.cpp)
file(GLOB H_FILES include/*.h)
set(SOURCE_FILES ${CPP_FILES} ${H_FILES})
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)

add_executable(RatesMC2 RatesMC2.cpp ${SOURCE_FILES}) # be sure there's exactly one main() function in the source files
