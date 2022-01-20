cmake_minimum_required(VERSION 3.0)
project(RatesMC)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 11)

file(GLOB CPP_FILES src/*.cpp)
file(GLOB H_FILES include/*.h)
set(SOURCE_FILES ${CPP_FILES} ${H_FILES})
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/include)

##set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
find_package(GSL REQUIRED)
link_libraries(GSL::gsl)


add_executable(RatesMC RatesMC.cpp ${SOURCE_FILES}) # be sure there's exactly one main() function in the source files

## Copy the RatesMC.in file over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/RatesMC.in
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Copy the R analysis scripts over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/analysis/
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Copy the AME Mass file over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mass_1.mas20
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

## Copy the NuBase file over if needed
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/nubase_3.mas20
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
