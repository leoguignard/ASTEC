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
project(vt-exec)

## #################################################################
## Input Files 
## #################################################################
# Inludes
include_directories(${BASIC_INCLUDE_DIRS})
include_directories(${IO_INCLUDE_DIRS})
include_directories(${VT_INCLUDE_DIRS})

# Exe names
SET(EXE_NAMES Arit  
	Convert  
	Create  
	Ctable  
	Ext  
	Ifmt  
	Logic  
	Par  
	Tran  
	amincir  
	axeMedianOfComponents  
	bbox  
	binary_noise  
	caracterisation  
	chamfer  
	changevals  
	connexe  
	contours  
	convolution  
	copy  
	createImage  
	detrend  
	distance  
	drawShapes
        extImage  
	filters  
	getVoxelValue
	gradient_unitaire  
	grey2color  
	histo  
	histoCopy
	histoJoint
	histoMerge
	histoStats
	laplacien  
	linearFilter  
	localFilter
	meanImages  
	melImages
	mergeSlices
	mescorr  
	minimum  
	mip_project  
	morpho  
	newTran  
	new_skiz  
	normalisation  
	norme_gradient  
	parcellisation  
	recfilters  
	reech3d  
	regionalext  
	regionalmax  
	relabel  
	sbinter  
	setVoxelSize  
	setVoxelValue  
	seuillage    
	slow_connexe  
	smooth  
	splitImage
	surface06  
	tgrowing
	thinning  
	watershed  
	zpar
)

# For exe name with suffixes needed: presents in another project
SET(EXE_NAMES_DUPLICATED 
	skeleton # libtopo
	skiz # libelf
)

message(STATUS "Bin names duplicated : skeleton (libtopo), skiz (libelf)") 
# message(WARNING "Bin names duplicated : skeleton (libtopo), skiz (libelf)") 


## #################################################################
## Build Directories 
## #################################################################

foreach(E ${EXE_NAMES})
  add_executable(${E} ${E}.c)
  target_link_libraries(${E} basic io vt)
endforeach(E)

foreach(E ${EXE_NAMES_DUPLICATED})
  add_executable(vt-exec_${E} ${E}.c)
  target_link_libraries(vt-exec_${E} basic io vt)
endforeach(E)



