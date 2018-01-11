/*************************************************************************
 * vt_elfskiz.h - extraction de bulles dans des images de mousse
 *
 * $Id: vt_elfskiz.h,v 1.4 2000/04/07 07:52:00 greg Exp $
 *
 * DESCRIPTION: 
 *
 * Outils d'extraction de parties connexes a l'aide d'un squelette par zone
 * d'influence et reconstruction
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sat May 12 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 * - Mon Mar 27 17:57:41 MET DST 2000, Gregoire Malandain
 *   Ajout de l'enumeration enumElfDistance, afin de pouvoir
 *   calculer une distance du chamfrein 5x5x5
 *
 *
 */


#ifndef _vt_elfskiz_h_
#define _vt_elfskiz_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_connexe.h>
#include <vt_distance.h>

#include <chamfer.h>


typedef enum {
  _EUCLIDIENNE_,
  _CHAMFER_3_,
  _CHAMFER_5_
} enumElfDistance;


extern int VT_SkizFromBorder( vt_image *image, 
			      vt_image *labels, 
			      vt_image *dist,
			      vt_distance *par,
			      int height,
			      double divisor );



extern int VT_SkizFromBorderWithoutSort( vt_image *image, 
					 vt_image *labels, 
					 vt_image *dist,
					 enumElfDistance typeDistance,
					 int height,
					 double divisor );


extern int VT_MaximaRegionauxWithList( vt_image *theIm /* input image */, 
				       vt_image *resIm /* output image */,
				       int hauteur,
				       double divisor );


extern int VT_MaximaRegionaux( vt_image *theIm /* input image */, 
			       vt_image *resIm /* output image */,
			       int hauteur,
			       double divisor );





#ifdef __cplusplus
}
#endif

#endif /* _vt_elfskiz_h_ */
