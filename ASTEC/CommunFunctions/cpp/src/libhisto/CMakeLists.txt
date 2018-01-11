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
project(libhisto)
set(LIB_NAME histoLIB)

message(STATUS "Libhisto : LIB suffix added to the lib") 
# message(WARNING "Libhisto : LIB suffix added to the lib") 

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# SRC files	
SET(SRC_FILES  vt_histoLung.c
	vt_sliceHisto.c
	powell.c
#	vt_geseuil1.c
)

# Exe names
SET(EXE_NAMES evalBleaching
	      fitOnHisto
	fitOnHistoLung
	fitGaussOnHisto
	slicesHisto
	slicesHisto-2
	slicesCorrect
	test-parzen
	imagesHisto
	Robust_affine_equalization
	
)
  

## #################################################################
## Build
## #################################################################

add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)
