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
project(diameter)
SET(LIB_NAME ${PROJECT_NAME})

## #################################################################
## Input Files 
## #################################################################

# Src files
SET(SRC_FILES alloc.c
	pick.c
	print.c
	rand.c
	read.c
	util.c
	exact-diam.c
	apprx-diam.c
	peled-like.c
	)

# Exe names
SET(EXE_NAMES test-diameter 
  diameter-test 
  gener-points
)


## #################################################################
## Build 
## #################################################################
add_library(${LIB_NAME} ${SRC_FILES})


foreach(E ${EXE_NAMES})

  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} )

endforeach(E)





