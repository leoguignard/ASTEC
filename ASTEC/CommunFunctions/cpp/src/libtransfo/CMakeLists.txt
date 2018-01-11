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
#

# Project Name
project(libtransfo)

## #################################################################
## Input Files 
## #################################################################

# Src files
SET(EXE_NAMES CompoTrsfs
   SmoothTrsfs
   MultTrsfs  
   InvTrsf
   Real2VoxelTrsf
   TestTrsf 
)

## #################################################################
## Build
## #################################################################
foreach(E ${EXE_NAMES})

  add_executable(${E} ${E}.c)
  
endforeach(E) 

