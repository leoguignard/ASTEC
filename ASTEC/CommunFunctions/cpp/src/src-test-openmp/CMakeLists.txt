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
project(src-test-openmp)

## #################################################################
## Input Files 
## #################################################################
# Include dirs
SET(EXE_NAMES hello
	for1
	for2
	test-schedule
)

## #################################################################
## Build 
## #################################################################

foreach(E ${EXE_NAMES})

  add_executable(${E} ${E}.c)

endforeach(E)





