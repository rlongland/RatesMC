cmake_minimum_required(VERSION 3.5)
project(RatesMC)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -Wall")

set(CMAKE_CXX_STANDARD 14)

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

#ADD_CUSTOM_TARGET(
#	do_always ALL
#   	COMMAND git describe --abbrev=10 --always
#	OUTPUT_VARIABLE GIT_COMMIT_HASH
#)
#execute_process(COMMAND git describe --tags --always
#    OUTPUT_VARIABLE GIT_COMMIT_HASH)
#add_definitions(-DGIT_COMMIT_HASH=$(GIT_COMMIT_HASH))
#shell git describe --abbrev=10 --always))

execute_process(
    COMMAND git describe --abbrev=10 --always
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE)

target_compile_definitions(RatesMC PRIVATE
    "-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
