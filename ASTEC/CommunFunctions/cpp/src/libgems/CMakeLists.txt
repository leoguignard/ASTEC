####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2012
#
# AUTHOR:
# Etienne Delclaux (etienne.delclaux@inria.fr)
# From Gregoire Malandain (gregoire.malandain@inria.fr)
# 
#

# Name
project(libgems)
set(LIB_NAME gems)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  vt_gaussienne.c
	vt_geambdd.c 
	vt_geambdd2.c
	vt_geamdist.c
	vt_geline.c
	vt_geline2D.c
	vt_gelrec3D.c
	vt_gemask.c
	vt_gemedian.c
	vt_geremove.c
	vt_gesmooth.c
	vt_isolated.c
	vt_gerecfilters.c 
	vt_geline3D.c
#	vt_geseuil1.c
)


SET(EXE_NAMES ge_anev
	ge_demo
	ge_median
	ge_smooth
	ge_testeval
	ge_thin
	ge_recfilters
	ge_contours
	ge_line 
#	ge_seuil
#	ge_seuil1
#	ge_seuil2 
#	ge_seuil3
#	ge_seuil3_3mip
#	ge_seuil_test_1 
#	ge_seuil_test_2
#	ge_testseuil
  )
  
## #################################################################
## Build
## #################################################################
# Build library
add_library(${LIB_NAME} ${SRC_FILES})
# Dependancies
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )
foreach(E ${EXE_NAMES})
  # Build execs with dependancies
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

