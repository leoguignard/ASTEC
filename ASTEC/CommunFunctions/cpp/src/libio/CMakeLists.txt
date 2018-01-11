####################################################################
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
project(libio)
set(LIB_NAME io)

## #################################################################
## Input Files 
## #################################################################
# Create file list in "${SRC_FILES}" var with "${CMAKE_CURRENT_SOURCE_DIR}/*.c" regex
# file(GLOB SRC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.c)

SET(SRC_FILES  ImageIO.c
	       analyze.c
	       bmp.c
	       bmpendian.c
	       bmpread.c
	       gif.c
	       gis.c
	       inr.c
	       iris.c
	       pnm.c )

# Exe names
SET(EXE_NAMES test-libio
)
  

## #################################################################
## Build
## #################################################################
# Build library
add_library( ${LIB_NAME} ${SRC_FILES} )

# Add dependency to zlib
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} )

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} io)
endforeach(E)

