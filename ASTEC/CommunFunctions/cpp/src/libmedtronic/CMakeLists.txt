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
#

# Name
project(libmedtronic)
set(LIB_NAME medtronic)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(C_FILES YelTr3D.c)
SET(EXE_NAME evalTr3D)

## #################################################################
## Build Directories 
## #################################################################
# Build library
add_library(${LIB_NAME} ${C_FILES})

# Build executable and link
add_executable(${EXE_NAME} ${EXE_NAME}.c)
target_link_libraries(${EXE_NAME} ${LIB_NAME})

