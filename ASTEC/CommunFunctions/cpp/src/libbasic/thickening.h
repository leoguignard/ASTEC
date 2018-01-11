/****************************************************
 * thickening.h - 
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 20 jui 2014 07:19:01 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _thickening_h_
#define _thickening_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <chamferdistance.h>



extern void setVerboseInThickening( int v );
extern void incrementVerboseInThickening();
extern void decrementVerboseInThickening();



typedef enum {
  _NO_SORTING_,
  _ITERATION_SORTING_,
  _DISTANCE_SORTING_
} enumTypeSort;



typedef struct {
  int maxIteration;
  int connectivity;
  typeChamferMask *theMask;
  enumTypeSort additionalSorting;
} typeThickeningParameters;

extern void initTypeThickeningParameters( typeThickeningParameters *p );




extern int InitializeThickeningImage( unsigned char *resBuf,
				      unsigned short *theDistance,
				      int *theDim );

extern int ThickeningImage( unsigned char *resBuf,
			    unsigned short *theDistance,
			    unsigned short *thePropagation,
			    int *theDim,
			    typeThickeningParameters *par );


#ifdef __cplusplus

}
#endif

#endif
