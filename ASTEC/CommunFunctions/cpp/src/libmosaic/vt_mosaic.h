/*************************************************************************
 * vt_mosaic.h
 *
 * $Id: vt_mosaic.h,v 1.1 2005/07/20 15:05:03 greg Exp $
 *
 * Copyright (c) INRIA 1999-2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Wed Jul 20 11:50:30 MEST 2005
  *
 * ADDITIONS, CHANGES
 *
 *
 */


#ifndef _vt_mosaic_h_
#define _vt_mosaic_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <vt_common.h>


typedef enum {
  AVERAGE,
  WEIGHTING,
  MAXIMUM
} enumFusionMode;


extern void _set_fusion_mode( enumFusionMode f );

typedef struct {
  int neighbor;
  vt_ipt origin;
  vt_ipt dim;
} typeIntersection;


#define _MAX_NEIGHBORS_ 8



typedef struct {
  char imagename[STRINGLENGTH];
  
  /* position of the image in the mosaic
     (in voxels)
   */
  vt_ipt offset_init; /* initial offset of the mosaic */
  vt_ipt offset;      /* offset, recomputed (if required) */

  /* "interior" part of the image, ie that does not intersect
     with other images
  */
  vt_ipt interior_corner;
  vt_ipt interior_dim;

  /* neighbors
   */
  typeIntersection intersection[_MAX_NEIGHBORS_];
  int n_neighbors;
  int n_max_neighbors;

  /* ... 
   */
  int flag;
  double sum;

  /* the image
   */
  vt_image *image;
} typeMosaicPart;





typedef struct {
  int n;
  int x;
  int y;
} point;


typedef struct {
  point p1;
  point p2;
  int flag;
} pair;





extern void print_mosaic_info( FILE* f, typeMosaicPart *mosaic, int nb );


/* structure du fichier:
%d: nombre d'images
%s %d %d: image offset en x offset en y
...
%s %d %d: image offset en x offset en y
*/
extern typeMosaicPart * _ReadParamFileWithOffsets( char *name, int *nb );
extern int _WriteParamFileWithOffsets( char *name, typeMosaicPart *mosaic, int nb );

/* structure du fichier:
%d: nombre d'images
%s: nom d'image
...
%s: nom d'image
%d (%d,%d) <-> %d (%d,%d) : #image1 point dans image 1 <-> id image #2
...
%d (%d,%d) <-> %d (%d,%d)
*/
  extern typeMosaicPart * _ReadParamFileWithPairs( char *name, int *nb, pair **pairs, int *np );
extern void _PreProcessMosaicWithPairs( typeMosaicPart *mosaic, int nb_images, pair *p, int np );


extern int _ReadMosaicImages( typeMosaicPart *mosaic, int nb );



extern int search_neighbors( typeMosaicPart *mosaic, int nb );



extern void build_mosaic( typeMosaicPart *mosaic, int nb_images, 
			  int txmin, int txmax, int tymin, int tymax );


extern vt_image * build_color_mosaic( char *name, typeMosaicPart *mosaic, 
				      int nb_images );
extern vt_image * build_large_image( char *name, typeMosaicPart *mosaic, 
				     int nb_images );



#ifdef __cplusplus
}
#endif

#endif /* _vt_mosaic_h_ */
