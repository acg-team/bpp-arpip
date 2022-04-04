# Licensed Materials - Property of Gholamhossein Jowkar
# Copyright (C) 2019-2023 by Gholamhossein Jowkar
#
# This work was supported by the Swiss National Science Foundation (SNF) grants $31003A\_176316$ to Dr. M. Anisimova.
# The funding body did not play any role in the design of the study and collection, analysis, and interpretation of
# data and in writing the code.
# -------------------------------------------------------------------------
#
# This file is part of ARPIP project
#
# ARPIP: Ancestral Sequence Reconstruction with insertions and deletions under the Poisson Indel Process
# ARPIP is a joint maximum likelihood approach for phylogenetic ancestral sequence reconstruction, capable of modeling
# indels both biological and mathematically.
#
#
# This software is based and extends the following libraries:
#
# - the Bio++ libraries
#   developed by the Bio++ Development Team <http://biopp.univ-montp2.fr>
#
# ARPIP is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.
#
# ARPIP is a free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU Lesser General Public
# License along with ARPIP. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 3.16)

set(ARPIP_SOFTWARENAME "ARPIP")
set(ARPIP_DESCRIPTION "Ancestral Sequence Reconstruction under PIP")

project(bpp_arpip ${ARPIP_VERSION}
        DESCRIPTION ${ARPIP_DESCRIPTION}
        LANGUAGES CXX)

set(CMAKE_CXX_FLAGS "-Wall -Wshadow -Wextra -Weffc++")


set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -fopenmp")

## Store the git hash of the current head
if(EXISTS "${PROJECT_SOURCE_DIR}/.git/HEAD")
    file(READ "${PROJECT_SOURCE_DIR}/.git/HEAD"
            PROJECT_SOURCE_VERSION)
    if("${PROJECT_SOURCE_VERSION}" MATCHES "^ref:")
        string(REGEX REPLACE "^ref: *([^ \n\r]*).*" "\\1"
                PROJECT_GIT_REF "${PROJECT_SOURCE_VERSION}")
        file(READ "${PROJECT_SOURCE_DIR}/.git/${PROJECT_GIT_REF}"
                PROJECT_SOURCE_VERSION)
    endif()
    string(STRIP "${PROJECT_SOURCE_VERSION}"
            PROJECT_SOURCE_VERSION)
endif()

# Store the build date
if(WIN32)
    execute_process(COMMAND "cmd" " /c date /t"
            OUTPUT_VARIABLE DATE)
    string(REGEX REPLACE "[^0-9]*(..).*" "\\1" MONTH "${DATE}")
    set(MONTHS ""
            "Jan" "Feb" "Mar" "Apr" "May" "Jun"
            "Jul" "Aug" "Sep" "Oct" "Nov" "Dec")
    list(GET MONTHS "${MONTH}" MONTH)
    string(REGEX REPLACE "[^/]*/(..)/(....).*" "\\1 ${MONTH} \\2"
            PROJECT_BUILD_DATE "${DATE}")
    execute_process(COMMAND "cmd" " /c echo %TIME%"
            OUTPUT_VARIABLE TIME)
    string(REGEX REPLACE "[^0-9]*(..:..:..).*" "\\1"
            PROJECT_BUILD_TIME "${TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
else()
    execute_process(COMMAND "date" "+%d %b %Y/%H:%M:%S"
            OUTPUT_VARIABLE DATE_TIME)
    string(REGEX REPLACE "([^/]*)/.*" "\\1"
            PROJECT_BUILD_DATE "${DATE_TIME}")
    string(REGEX REPLACE "[^/]*/([0-9:]*).*" "\\1"
            PROJECT_BUILD_TIME "${DATE_TIME}")
    execute_process(COMMAND git describe --abbrev=0 --tags
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            OUTPUT_VARIABLE GIT_TAG_VERSION)
    string(REGEX REPLACE "\n" "" GIT_TAG_VERSION_STRIPPED "${GIT_TAG_VERSION}")
endif()

set(ARPIP_RELTIME ${PROJECT_BUILD_TIME})
set(ARPIP_RELDATE ${PROJECT_BUILD_DATE})
set(ARPIP_VERSION ${GIT_TAG_VERSION_STRIPPED})

project(${ARPIP_SOFTWARENAME} )

add_definitions(-DPRJ_GITBRANCH=\"${PROJECT_GIT_REF}\"
        -DPRJ_VERSION=\"${ARPIP_VERSION}\"
        -DPRJ_GITREF=\"${PROJECT_SOURCE_VERSION}\"
        -DPRJ_DESC=\"${ARPIP_DESCRIPTION}\"
        -DPRJ_NAME=\"${ARPIP_SOFTWARENAME}\"
        -DPRJ_DATE=\"${ARPIP_RELDATE}\"
        -DPRJ_TIME=\"${ARPIP_RELTIME}\" )

message("-- Compilation will be performed with the following release of the software:
\tbranch   ${PROJECT_GIT_REF}
\tref      ${PROJECT_SOURCE_VERSION}
\ttime     ${ARPIP_RELDATE} ${ARPIP_RELTIME}
\tcurrent  ${ARPIP_VERSION} (latest version)")



## --- the IMPROVED_FIND_LIBRARY macro (taken from the Bio++ CMakeLists.txt) ---
macro(IMPROVED_FIND_LIBRARY OUTPUT_LIBS lib_name include_to_find)
    FIND_PATH(${lib_name}_INCLUDE_DIR ${include_to_find})
    SET(${lib_name}_NAMES ${lib_name} ${lib_name}lib ${lib_name}dll)
    FIND_LIBRARY(${lib_name}_LIBRARY NAMES ${${lib_name}_NAMES})

    IF(${lib_name}_LIBRARY)
        MESSAGE("-- Library ${lib_name} found here:")
        MESSAGE("    ${${lib_name}_INCLUDE_DIR}/")
        MESSAGE("    ${${lib_name}_LIBRARY}")
    ELSE(${lib_name}_LIBRARY)
        MESSAGE(FATAL_ERROR "${lib_name} required but not found.")
    ENDIF(${lib_name}_LIBRARY)

    #add the dependency:
    INCLUDE_DIRECTORIES(${${lib_name}_INCLUDE_DIR})
    SET(${OUTPUT_LIBS} ${${OUTPUT_LIBS}} ${${lib_name}_LIBRARY})
ENDMACRO(IMPROVED_FIND_LIBRARY)



#--- Libraries & includes ---
# Dependencies not covered by find package should be found in the following directories
if (${CMAKE_PREFIX_PATH})
    include_directories("${CMAKE_PREFIX_PATH}/include")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}}/lib")
    LINK_DIRECTORIES("${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib")
    message(STATUS "Looking for libraries in the following directory: ${CMAKE_PREFIX_PATH}/lib64")
    message(STATUS "Looking for headers in the following directory: ${CMAKE_PREFIX_PATH}/include")
endif ()


#include_directories("include/")
MESSAGE(STATUS "Loading libraries & includes:")

find_package (glog 0.5.0 REQUIRED)


# Find Boost Libraries
set(Boost_USE_MULTITHREADED      ON)
set(Boost_USE_STATIC_RUNTIME    OFF)
find_package(Boost REQUIRED COMPONENTS system filesystem)

##set(CMAKE_FIND_LIBRARY_SUFFIXES ".a ${CMAKE_FIND_LIBRARY_SUFFIXES}") # uncomment for static compilation
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-core "Bpp/Clonable.h")
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-seq "Bpp/Seq/Sequence.h")
IMPROVED_FIND_LIBRARY(LIBS_BPP bpp-phyl "Bpp/Phyl/Tree.h")

MESSAGE(STATUS "Done!")


#--- Targets ---
set(EXE ${PROJECT_NAME})

file(GLOB_RECURSE HEADERS include/*.h include/*.hpp)
file(GLOB_RECURSE SOURCES src/*.cpp src/*.hpp)
add_executable(${EXE} ${SOURCES} ${HEADERS})
target_link_libraries(${EXE} ${LIBS_BPP})
target_link_libraries(${EXE} glog::glog)

MESSAGE(STATUS "Targets OK")

add_subdirectory(src)
#include_directories(src)

#add_subdirectory(test)
#add_subdirectory(lib/googletest)

#add_subdirectory("/usr/src/googletest" ${CMAKE_BINARY_DIR}/gtest)
#add_subdirectory(/usr/src/googletest)

# Tests
enable_testing ()
include_directories(${GTEST_INCLUDE_DIRS})
#add_executable(gtest_test test/main.cpp)
add_test(AllTests test)

include (CTest)
if (BUILD_TESTING)
#    add_subdirectory(src)
    include_directories(test)
    add_subdirectory(test)
endif (BUILD_TESTING)
MESSAGE(STATUS "Test OK")

