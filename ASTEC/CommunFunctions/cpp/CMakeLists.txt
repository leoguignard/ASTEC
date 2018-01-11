############################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 

## #################################################################
## Header
## #################################################################
# Specify the version being used
cmake_minimum_required(VERSION 2.8)

# Project name
project(project_vt)

## #################################################################
## Options
## #################################################################
# Specify the version being used
option(BUILD_SHARED_LIBS      " Build shared libraries"  true)
option(BASIC_BUILD_TEMPLATES  " Regenerate templated files in basic lib" false)
option(BUILD_TESTING	       " Build testing" false)
option(vt_USE_OPENMP	       " Use OpenMP" true) 
# SET(FILTERS_BUILD_TEST CACHE BOOL ON)
# SET(FILTERS_BUILD_INSIGHT CACHE BOOL ON)
# SET(FILTERS_BUILD_PURIFY CACHE BOOL ON)
# SET(FILTERS_BUILD_QUANTIFY CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_COMPARE CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_TIME CACHE BOOL ON)
# SET(TOPOEXEC_BUILD_STATS CACHE BOOL ON)


## #################################################################
## Config file
## #################################################################
# Cmake config file directory path
SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Call include vtBuild.cmake file. This file configure builds ouput (build/bin, build/lib)
include(vtBuild)

# Call include vtDependencies.cmake file. This file configure the externals dependencies of vt (zlib, openmp, m)
include(vtDependencies)

include(vtPack)

# To improve
SET(vt_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/src/baladin
  ${PROJECT_SOURCE_DIR}/src/blockmatching
  ${PROJECT_SOURCE_DIR}/src/diameter
  ${PROJECT_SOURCE_DIR}/src/libbarbier
  ${PROJECT_SOURCE_DIR}/src/libbasic
  ${PROJECT_SOURCE_DIR}/src/libconfocal
  ${PROJECT_SOURCE_DIR}/src/libdavid
  ${PROJECT_SOURCE_DIR}/src/libelf
  ${PROJECT_SOURCE_DIR}/src/libfilters
  ${PROJECT_SOURCE_DIR}/src/libgems
  ${PROJECT_SOURCE_DIR}/src/libhisto
  ${PROJECT_SOURCE_DIR}/src/libio
  ${PROJECT_SOURCE_DIR}/src/libipmc
  ${PROJECT_SOURCE_DIR}/src/libjacobien
  ${PROJECT_SOURCE_DIR}/src/libmedtronic
  ${PROJECT_SOURCE_DIR}/src/libmosaic
  ${PROJECT_SOURCE_DIR}/src/libnoesis
  ${PROJECT_SOURCE_DIR}/src/librealviz
  ${PROJECT_SOURCE_DIR}/src/libspect
  ${PROJECT_SOURCE_DIR}/src/libtopo
  ${PROJECT_SOURCE_DIR}/src/libtransfo
  ${PROJECT_SOURCE_DIR}/src/libvp
  ${PROJECT_SOURCE_DIR}/src/libvt
  ${PROJECT_SOURCE_DIR}/src/src-4D
  ${PROJECT_SOURCE_DIR}/src/src-test-openmp
  ${PROJECT_SOURCE_DIR}/src/src-vt-test
  ${PROJECT_SOURCE_DIR}/src/topo-exec
  ${PROJECT_SOURCE_DIR}/src/vt-exec
  )

set(vt_LIBRARY_DIRS ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY})
set(vt_CMAKE_DIR    ${CMAKE_MODULE_PATH})
set(vt_USE_FILE     ${PROJECT_BINARY_DIR}/vtUse.cmake)

## Setup Config File
if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/vtConfig.cmake.in)
  configure_file( ## Build tree configure file
    ${PROJECT_SOURCE_DIR}/cmake/vtConfig.cmake.in
    ${PROJECT_BINARY_DIR}/vtConfig.cmake
    @ONLY IMMEDIATE)
endif()

## Setup use file
if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/vtUse.cmake.in)
configure_file( ## Common use file
  ${PROJECT_SOURCE_DIR}/cmake/vtUse.cmake.in
  ${PROJECT_BINARY_DIR}/vtUse.cmake
  @ONLY IMMEDIATE)
endif()


## #################################################################
## ADD_SUBIDIRECTORY
## #################################################################

# Source directory
add_subdirectory(src)
