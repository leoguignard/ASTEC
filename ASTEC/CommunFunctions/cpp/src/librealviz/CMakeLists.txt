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
project(librealviz)
set(LIB_NAME realviz)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  is_distance.c
	rl_interpol.c
	vt_interpol.c
	vt_recUS2D.c
	_predef_morpho.c
	)

# Exe names
SET(EXE_NAMES test-unit
	test-inter-pnm
	test-dist3D-point
	im2dist
	dists2im
	dists2im-inc
	im2iso
	recUS2D
	gener-synth
	test-inter-inr
  )
  
## #################################################################
## Build
## #################################################################

# Build library
add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

