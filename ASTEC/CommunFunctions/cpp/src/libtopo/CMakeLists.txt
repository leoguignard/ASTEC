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
project(libtopo)
SET(LIB_NAME topo)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Src files
SET(SRC_FILES vt_bdd.c
	vt_skeleton.c
	vt_endpoints.c
)

# Exe names
SET(EXE_NAMES 
	thickLabels
	backLabels 
	countBackLabels
	forwLabels
	test-bdd
	test-endpoints
)

# For exe name with suffixes needed: presents in another project
SET(EXE_NAMES_DUPLICATED 
	skeleton  # vt-exec
	upperView # libconfocal
)

message(STATUS "bin names duplicated : skeleton(vt-exec) and upperView(libconfocal)") 
# message(WARNING "bin names duplicated : skeleton(vt-exec) and upperView(libconfocal)") 

## #################################################################
## Build 
## #################################################################
# Build lib and link
add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt )

# Build binaries and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

# Build binaries and link
foreach(E ${EXE_NAMES_DUPLICATED})
  add_executable(${LIB_NAME}_${E} ${E}.c)
  target_link_libraries(${LIB_NAME}_${E} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)





