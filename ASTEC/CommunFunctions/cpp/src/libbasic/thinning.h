/****************************************************
 * thinning.h - 
 *
 * $Id: thinning.h,v 1.4 2000/10/20 13:40:22 greg Exp $
 *
 * Copyright©INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Mon Aug  7 14:51:59 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#ifndef _thinning_h_
#define _thinning_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <string.h>

#include <typedefs.h>

#include <chamferdistance.h>
#include <t04t08.h>
#include <t06t26.h>



typedef enum {
  _04THICKNESS_ = 4,
  _08THICKNESS_ = 8,
  _06THICKNESS_ = 6,
  _18THICKNESS_ = 18,
  _26THICKNESS_ = 26
} enumThickness;
    
typedef enum {
  _SURFACE_,
  _PURE_SURFACE_,
  _CURVE_,
  _PURE_CURVE_,
  _NO_END_POINT_
} enumTypeEndPoint;




typedef struct {
  int anchorValue;
  int connectivity;
  int cyclesBeforeEnding;
  int distanceBeforeEnding;
  enumThickness typeThickness;
  enumTypeEndPoint typeEndPoint;
  typeChamferMask *theMask;
} typeThinningParameters;



extern void initTypeThinningParameters( typeThinningParameters *p );




extern int ThinImageWithAllParams( unsigned char *theBuf,
				   unsigned char *resBuf,
				   int *theDim,
				   const int chamfer,
				   typeThinningParameters *p );

extern int ThinImageWithDistanceAndAllParams( unsigned char *theBuf,
					      unsigned char *resBuf,
					      unsigned short *theDistance,
					      int *theDim,
					      typeThinningParameters *p );


extern void Thinning_verbose();
extern void Thinning_noverbose();
extern void Thinning_timeVerbose();
extern void Thinning_timeNoverbose();
#ifdef __cplusplus

}
#endif

#endif
