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
SET( BASIC_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/libbasic )
SET( IO_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libio )
SET( VT_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libvt ) 
SET( FILTERS_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/libfilters ) 
SET( BLOCKMATCHING_INCLUDE_DIRS      ${CMAKE_CURRENT_SOURCE_DIR}/blockmatching ) 
SET( NOESIS_INCLUDE_DIRS  ${CMAKE_CURRENT_SOURCE_DIR}/libnoesis )
SET( DAVID_INCLUDE_DIRS   ${CMAKE_CURRENT_SOURCE_DIR}/libdavid )

# Independants
add_subdirectory(libtransfo)
add_subdirectory(diameter)
add_subdirectory(src-test-openmp)
add_subdirectory(libio)
add_subdirectory(libbasic)

# Depend on BASIC and IO
add_subdirectory(baladin)
add_subdirectory(blockmatching)
add_subdirectory(libvt)

# Depend on BASIC, IO and BLOCKMATCHING
add_subdirectory(libmatchingtools)
if (  NOT LEMON_FOUND ) 
  message( WARNING, " lemon was NOT found" )
  message( WARNING, " libbiolibtools will NOT be built" )
else( NOT LEMON_FOUND )
  add_subdirectory(libbiolibtools)
endif( NOT LEMON_FOUND )

# Depend on BASIC, IO, and VT
add_subdirectory(libbarbier)
add_subdirectory(libmedtronic)
add_subdirectory(libtopo)
add_subdirectory(libdavid)
add_subdirectory(libipmc)
add_subdirectory(libjacobien)
add_subdirectory(libmembrane)
add_subdirectory(libnoesis)
add_subdirectory(libconfocal)
add_subdirectory(libfilters)
add_subdirectory(libgems)
add_subdirectory(libhisto)
add_subdirectory(librealviz)
add_subdirectory(libspect)
add_subdirectory(libmosaic)
add_subdirectory(vt-exec)
add_subdirectory(src-vt-test)
add_subdirectory(src-4D)
add_subdirectory(topo-exec)

# Depend on BASIC, IO, VT, and NOESIS
add_subdirectory(libelf)
add_subdirectory(libtools)

# Depend on BASIC, IO, VT, and DAVID
add_subdirectory(libvp)







