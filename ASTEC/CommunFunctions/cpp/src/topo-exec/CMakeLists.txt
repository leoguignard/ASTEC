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
project(topo-exec)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# SRC FILES
SET(SRC_FILES tube3D.c
	tube3Dutil.c
	vt_tube3D.c
	vt_tube3Delf.c
	vt_tube2D.c	
	vt_tube2Dmatlab.c
	vt_contoursMatlab.c
	vt_tubeutils.c
	vt_link.c
	vt_contours.c
	vt_isocontours.c
	)

SET(EXE_NAMES test-distfield
	test-all-distfield
	)
		


## #################################################################
## Build 
## #################################################################

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${ZLIB_LIBRARIES} basic io vt)
endforeach(E)

## #################################################################
## COMPARE. history 
## #################################################################
#if(${TOPOEXEC_BUILD_COMPARE})
#
#  SET(COMPARE_NAMES compare_t06 
#	compare_t26
#	compare_bdd_simple 
#	compare_simple
# )
#  foreach(C ${COMPARE_NAMES})
#    add_executable(${C} ${C}.c)
#    target_link_libraries(${C} ${ZLIB_LIBRARIES} basic io vt)
#  endforeach(C)
#endif(${TOPOEXEC_BUILD_COMPARE})
#
#
## #################################################################
## TIME. history 
## #################################################################
#
#
#if(${TOPOEXEC_BUILD_TIME})
#SET(TIME_NAMES time_t06 
#  time_t26 
#  time_t06t26 
#  time_bdd_simple
#)
#  foreach(T ${TIME_NAMES})
#    add_executable(${NAME}_${T} ${T}.c)
#    target_link_libraries(${NAME}_${T} ${ZLIB_LIBRARIES} basic io vt)
#    set_target_properties(${NAME}_${T} PROPERTIES OUTPUT_NAME ${C})
#  endforeach(T)
#endif(${TOPOEXEC_BUILD_TIME})
#
## #################################################################
## STATS. history 
## #################################################################
#
#if(${TOPOEXEC_BUILD_STATS})
#SET(STATS_NAMES stats_t06 
#  stats_t06noopt 
#  stats_t26
#  stats_t26noopt
#  stats_t06t26 
#  stats_t06t26_simple
#	stats_bdd_simple
#	)
#  foreach(S ${STATS_NAMES})
#    add_executable(${NAME}_${S} ${S}.c)
#    target_link_libraries(${NAME}_${S} ${ZLIB_LIBRARIES} basic io vt)
#    set_target_properties(${NAME}_${S} PROPERTIES OUTPUT_NAME ${S})
#  endforeach(S)
#endif(${TOPOEXEC_BUILD_STATS})

	

	






