############################################################
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
project(libconfocal)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Exe names
SET(EXE_NAMES
	best-translation	
)

# For exe name with suffixes needed: presents in another project
SET(EXE_NAMES_DUPLICATED 
	upperView # libtopo
)

message(STATUS "Bin names duplicated : upperview (libtopo)") 
# message(WARNING "Bin names duplicated : upperview (libtopo)") 

## #################################################################
## Build 
## #################################################################

# Build library and declare dependacies
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

# Build binaries and link
foreach(E ${EXE_NAMES_DUPLICATED})
  add_executable(libtopo_${E} ${E}.c)
  target_link_libraries(libtopo_${E} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)







