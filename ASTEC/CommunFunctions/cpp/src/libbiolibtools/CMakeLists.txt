####################################################################
#
# $Id$
#
# Copyright (c) INRIA 2013
#
# AUTHOR:
# Gregoire Malandain (gregoire.malandain@inria.fr)
# From Etienne Delclaux (etienne.delclaux@inria.fr)
#

# Name
project(libbiolibtools)
set(LIB_NAME biolibtools)




## #################################################################
## Input Files 
## #################################################################
# Includes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${BLOCKMATCHING_INCLUDE_DIRS})

# Src files
SET(SRC_FILES bal-biolib-tools.c
	      bal-biolib-tracker.cpp
)

# Exe names
SET(EXE_NAMES blApplyTrsfToDetectionResults
	      blBuildTracks
	      blCropDetectionResults
	      blCropTrackResults
	      blDiffTrackResults
	      blDriftQualityMeasure
	      blMaskDetectionResults
	      blMaskTrackResults
	      blPlotDetectionResults
	      blPlotTrackResults
	      blSelectTrackResults
	      blStatTrackResults
)
  
## #################################################################
## Build
## #################################################################

# build library and link
add_library(${LIB_NAME} ${LIB_STATUS} ${SRC_FILES})
target_link_libraries( ${LIB_NAME} 
		       blockmatchingLIB  
		       ${ZLIB_LIBRARIES} basic io
		       ${LEMON_LIBRARIES} )

# build execs and link
foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} ${LIB_NAME} blockmatchingLIB basic io)
endforeach(E)

