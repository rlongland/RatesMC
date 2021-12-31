cmake_minimum_required(VERSION 2.4)
project(RatesMC2)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

file(GLOB CPP_FILES src/*.cpp)
file(GLOB H_FILES include/*.h)
set(SOURCE_FILES ${CPP_FILES} ${H_FILES})
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)

set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
link_libraries(gsl)


add_executable(RatesMC2 RatesMC2.cpp ${SOURCE_FILES}) # be sure there's exactly one main() function in the source files

## Copy the RatesMC.in file over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/RatesMC.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Copy the deploy file over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/deploy
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Copy the R analysis scripts over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/analysis
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
