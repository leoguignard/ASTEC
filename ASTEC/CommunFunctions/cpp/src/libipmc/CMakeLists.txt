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
project(libipmc)
set(LIB_NAME ipmc)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(FILES fit-distribution-tools.c
	  tracking-tools.c
)

# Exe names
SET(EXE_NAMES 
    compute-speed
    fit-distribution
    remove-spurious
    speed-from-area
    tracking-colocalization
    tracking-pits
    tracking-profiles
)


## #################################################################
## Build 
## #################################################################
# Build binaries and link
add_library(${LIB_NAME} ${FILES})
target_link_libraries(${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME}  ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)





