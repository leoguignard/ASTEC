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
project(baladin)
set(LIB_NAME baladinLIB)

message(STATUS "Baladin : LIB suffix added to the lib") 
# message(WARNING "Baladin : LIB suffix added to the lib") 

## #################################################################
## Input Files 
## #################################################################

# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})

# Src files
SET(SRC_FILES  
	balimage.c
	estimateur.c
	initialisation.c
	matrix.c
	mesure.c
	misc.c
	behavior.c
	py_image.c
	vecteur.c
	)

# Exe names
SET(EXE_NAME baladin
  )
  
  
## #################################################################
## Build
## #################################################################
# Build library and link to ZLIB, basic and io
add_library(${LIB_NAME} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} ${ZLIB_LIBRARIES} basic io)

# Build executable, and link to libbaladin ZLIB, basic and io
add_executable(${EXE_NAME} ${EXE_NAME}.c)
target_link_libraries(${EXE_NAME} ${LIB_NAME} ${ZLIB_LIBRARIES} basic io)

