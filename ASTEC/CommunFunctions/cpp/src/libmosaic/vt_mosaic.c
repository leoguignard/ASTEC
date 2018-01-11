/*************************************************************************
 * vt_mosaic.c - 
 *
 * $Id: vt_mosaic.c,v 1.1 2005/07/20 15:05:03 greg Exp $
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

static int _debug_ = 1;
static int _verbose_ = 1;

#include <stdlib.h>
#include <stdio.h>

#include <chunks.h>

#include <vt_mosaic.h>


/* how to fuse multiple images
 */


static enumFusionMode fusionMode = AVERAGE;

void _set_fusion_mode( enumFusionMode f )
{
  fusionMode = f;
}






/***************************************************
 *
 * 
 *
 ***************************************************/

static void _init_typeIntersection( typeIntersection *i ) 
{
  i->neighbor = -1;
  i->origin.x = i->origin.y = i->origin.z = -1;
  i->dim.x = i->dim.y = i->dim.z = 0;
}

static void _init_typeMosaicPart( typeMosaicPart *mp )
{
  int n;

  mp->imagename[0] = '\0';
  
  mp->offset_init.x = 0;
  mp->offset_init.y = 0;
  mp->offset_init.z = 0;

  mp->offset.x = 0;
  mp->offset.y = 0;
  mp->offset.z = 0;

  mp->interior_corner.x = 0;
  mp->interior_corner.y = 0;
  mp->interior_corner.z = 0;

  mp->interior_dim.x = 0;
  mp->interior_dim.y = 0;
  mp->interior_dim.z = 0;

  for ( n=0; n<_MAX_NEIGHBORS_; n++ )
    _init_typeIntersection( &(mp->intersection[n]) );
  mp->n_neighbors = 0;
  mp->n_max_neighbors = _MAX_NEIGHBORS_;

  mp->flag = 0;
  mp->sum = 0.0;

  mp->image = NULL;
}





static void print_mosaic_part_info( FILE* f, typeMosaicPart *mosaic, int print_neighbors )
{
  int i;
  fprintf( f, " '%30s' init (+%6d +%6d +%4d) final (+%6d +%6d +%4d) \n", 
	   mosaic->imagename, 
	   mosaic->offset_init.x, mosaic->offset_init.y, mosaic->offset_init.z, 
	   mosaic->offset.x, mosaic->offset.y, mosaic->offset.z );
  if ( 0 ) {
    fprintf( f, "     interior = (%3d %3d %3d) + [%3d %3d %3d]\n",
	     mosaic->interior_corner.x, mosaic->interior_corner.y, mosaic->interior_corner.z, 
	     mosaic->interior_dim.x, mosaic->interior_dim.y, mosaic->interior_dim.z );
  }
  if ( print_neighbors ) {
    for ( i=0; i<mosaic->n_neighbors; i++ )
      fprintf( f, "     - neighbor #%d: image #%d, at (%3d %3d %3d) + [%3d %3d %3d]\n", i,
	       mosaic->intersection[i].neighbor,
	       mosaic->intersection[i].origin.x,
	       mosaic->intersection[i].origin.y,
	       mosaic->intersection[i].origin.z,
	       mosaic->intersection[i].dim.x,
	       mosaic->intersection[i].dim.y,
	       mosaic->intersection[i].dim.z );
  }
}


void print_mosaic_info( FILE* f, typeMosaicPart *mosaic, int nb )
{
  char *proc = "print_mosaic_info";
  int i;
  
  if ( search_neighbors( mosaic, nb ) != 1 ) {
    fprintf( stderr, "%s: unable to compute neighborhood relationships\n", proc );
    return;
  }

  for (i=0; i< nb; i++ ) {
    fprintf( f, "#%2d:", i );
    print_mosaic_part_info( f, &(mosaic[i]), 1 );
  }
}



/***************************************************
 *
 * I/O
 *
 ***************************************************/


typeMosaicPart * _ReadParamFileWithOffsets( char *name, int *nb )
{
  char *proc = "_ReadParamFileWithOffsets";
  FILE *f;
  typeMosaicPart *mosaic;
  int i =0;
  
  /* number of files
   */

  *nb = 0;

  f = fopen( name, "r" );
 if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening file '%s'\n", proc, name );
    return( 0 );
  }
 
  if ( fscanf( f, "%d", nb ) != 1 ) {
    fprintf( stderr, "%s: Error in reading image number at the first line\n", proc );
    return( NULL );
  }
  
  mosaic = malloc( (*nb)*sizeof(typeMosaicPart) );

  /* file names
   */

  for ( i=0; i<(*nb); i++ ) {

    _init_typeMosaicPart( &(mosaic[i]) );
      
    if ( fscanf( f, "%s %d %d", mosaic[i].imagename, 
		 &(mosaic[i].offset_init.x), &(mosaic[i].offset_init.y) ) != 3 ) {
      fprintf( stderr, "%s: Error in reading filename #%d\n", proc, i );
      free( mosaic );
      *nb = 0;
      return( NULL );
    }
    
    mosaic[i].offset = mosaic[i].offset_init;
    
  }
  
  fclose( f );
  
  return( mosaic );
}





int _WriteParamFileWithOffsets( char *name, typeMosaicPart *mosaic, int nb )
{
  char *proc = "_WriteParamFileWithOffsets";
  FILE *f;
  int i =0;
  
  /* number of files
   */


  f = fopen( name, "w" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening file '%s'\n", proc, name );
    return( 0 );
  }

  fprintf( f, "%d\n", nb );
  
  /* file names
   */

  for ( i=0; i<nb; i++ ) {
    fprintf( f, "%s %d %d\n", mosaic[i].imagename, mosaic[i].offset.x, mosaic[i].offset.y ); 
  }
    
  
  fclose( f );
  
  return( 1 );
}





typeMosaicPart * _ReadParamFileWithPairs( char *name, int *nb,  pair **pairs, int *npairs )
{
  char *proc = "_ReadParamFileWithPairs";
  FILE *f;
  typeMosaicPart *mosaic;
  pair *p;
  int np = 10, tp = 10;
  int i =0;
  
  /* number of files
   */

  *nb = 0;

  f = fopen( name, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening file '%s'\n", proc, name );
    return( 0 );
  }

  if ( fscanf( f, "%d", nb ) != 1 ) {
    fprintf( stderr, "%s: Error in reading image number at the first line\n", proc );
    return( NULL );
  }
  
  mosaic = malloc( (*nb)*sizeof(typeMosaicPart) );

  /* file names
   */
  
  for ( i=0; i<(*nb); i++ ) {

    _init_typeMosaicPart( &(mosaic[i]) );

    if ( fscanf( f, "%s", mosaic[i].imagename ) != 1 ) {
      fprintf( stderr, "%s: Error in reading filename #%d\n", proc, i );
      free( mosaic );
      *nb = 0;
      return( NULL );
    }

  }
  
  /* pairs
   */

  i = 0;
  p = malloc( np * sizeof(pair) );

  while ( fscanf( f, "%d (%d,%d) <-> %d (%d,%d)",
		  &(p[i].p1.n), &(p[i].p1.x), &(p[i].p1.y), 
		  &(p[i].p2.n), &(p[i].p2.x), &(p[i].p2.y) ) == 6 ) {
    i ++;
    if ( i == np ) {
      np += tp;
      p = realloc( p, np * sizeof(pair) );
    }
  }
  
  fclose( f );

  /*
   */
  *pairs = p;
  *npairs = np = i;
  if ( np == 0 ) {
    fprintf( stderr, "%s: read no pairings\n", proc );
    return( mosaic );
  }
  
  if ( 1 && _debug_ )
    fprintf( stderr, "%s: read %d pairings\n", proc, i );

  /* initializing pairs
   */

  for (i=0; i< np; i++ ) {
    p[i].p1.n -= 1;
    p[i].p2.n -= 1;
    p[i].flag = 0;
  }


  if ( 0 ) {
    for (i=0; i< (*nb); i++ ) {
      fprintf( stdout, "#%2d: '%40s' +%8d +%8d\n", i, mosaic[i].imagename, 
	       mosaic[i].offset_init.x, mosaic[i].offset_init.y );
    }
    for (i=0; i< np; i++ ) {
      fprintf( stdout, "#%2d: %d (%d,%d) <-> %d (%d,%d)\n", i,
	       p[i].p1.n, p[i].p1.x, p[i].p1.y, 
	       p[i].p2.n, p[i].p2.x, p[i].p2.y );
      
    }
  }



  return( mosaic );
}




void _PreProcessMosaicWithPairs( typeMosaicPart *mosaic, int nb_images, pair *p, int np )
{
  int i, c;
  
  if ( 0 && _debug_ ) {
    for (i=0; i< np; i++ ) {
      fprintf( stdout, "#%2d: %d (%d,%d) <-> %d (%d,%d)\n", i,
	       p[i].p1.n, p[i].p1.x, p[i].p1.y, 
	       p[i].p2.n, p[i].p2.x, p[i].p2.y );
    }
  }


  /* correction (fiji, analyze)
   */
  for ( i=0; i< np; i++ ) {
    p[i].p1.y = mosaic[ p[i].p1.n ].image->dim.y - p[i].p1.y;
    p[i].p2.y = mosaic[ p[i].p2.n ].image->dim.y - p[i].p2.y;
  }

  /* offset computation	
   */
  mosaic[ p[0].p1.n ].flag = mosaic[ p[0].p2.n ].flag = 1;
  mosaic[ p[0].p2.n ].offset_init.x = mosaic[ p[0].p1.n ].offset_init.x 
    + p[0].p1.x - p[0].p2.x;
  mosaic[ p[0].p2.n ].offset_init.y = mosaic[ p[0].p1.n ].offset_init.y 
    + p[0].p1.y - p[0].p2.y;
  p[0].flag = 1;

  do {
    c = 0;
    for ( i=0; i<np; i++ ) {
      if ( p[i].flag == 1 ) continue;
      if ( mosaic[ p[i].p1.n ].flag == 1 && mosaic[ p[i].p2.n ].flag == 1 ) {
	p[i].flag = 1;
	continue;
      }
      if ( mosaic[ p[i].p1.n ].flag == 0 && mosaic[ p[i].p2.n ].flag == 0 ) 
	continue;
      
      if ( mosaic[ p[i].p1.n ].flag == 1 ) {
	mosaic[ p[i].p2.n ].offset_init.x = mosaic[ p[i].p1.n ].offset_init.x 
	  + p[i].p1.x - p[i].p2.x;
	mosaic[ p[i].p2.n ].offset_init.y = mosaic[ p[i].p1.n ].offset_init.y 
	  + p[i].p1.y - p[i].p2.y;
	mosaic[ p[i].p2.n ].flag = 1;
      }
      else {
	mosaic[ p[i].p1.n ].offset_init.x = mosaic[ p[i].p2.n ].offset_init.x 
	  + p[i].p2.x - p[i].p1.x;
	mosaic[ p[i].p1.n ].offset_init.y = mosaic[ p[i].p2.n ].offset_init.y 
	  + p[i].p2.y - p[i].p1.y;
	mosaic[ p[i].p1.n ].flag = 1;
      }
      p[i].flag = 1;
      c = 1;
    }
    
  } while ( c != 0 );

  for ( i=0; i<nb_images; i++ ) {
    mosaic[i].offset = mosaic[i].offset_init;
  }
}






int _ReadMosaicImages( typeMosaicPart *mosaic, int nb )
{
  char *proc = "_ReadMosaicImages";
  int i;

  if ( _verbose_ ) {
    fprintf( stderr, "Reading images:" );
  }

  for ( i=0; i<nb; i++ ) {
    if ( _verbose_ ) {
      fprintf( stderr, " #%d", i+1 );
    }
    mosaic[i].image = _VT_Inrimage( mosaic[i].imagename );
    if ( mosaic[i].image == NULL ) {
      fprintf( stderr, "%s: unable to read image '%s'\n", proc, mosaic[i].imagename );
      return( 0 );
    }
    mosaic[i].interior_dim.x = mosaic[i].image->dim.x;
    mosaic[i].interior_dim.y = mosaic[i].image->dim.y;
    mosaic[i].interior_dim.z = mosaic[i].image->dim.z;
  }

  if ( _verbose_ ) {
    fprintf( stderr, "\n" );
  }


  return( 1 );
}





/***************************************************
 *
 * operations sur les intersections
 *
 ***************************************************/

int intersection_dimensions( vt_ipt *offset1, vt_ipt *dim1,
			     vt_ipt *offset2, vt_ipt *dim2,
			     vt_ipt *o,
			     vt_ipt *d )
{
  vt_ipt origin;
  vt_ipt dim;

  origin.x = origin.y = origin.z = -1;
  dim.x = dim.y = dim.z = 0;
  
  /* intersection along X
   */
  if ( offset1->x < offset2->x ) {
    if ( dim1->x+offset1->x > offset2->x ) {
      origin.x = offset2->x;
      dim.x = dim1->x+offset1->x - offset2->x;
      if ( dim.x > dim2->x ) dim.x = dim2->x;
    }
  }
  else {
    if ( dim2->x+offset2->x > offset1->x) {
      origin.x = offset1->x;
      dim.x = dim2->x+offset2->x - offset1->x;
      if ( dim.x > dim1->x ) dim.x = dim1->x;
    }
  }

  /* intersection along Y
   */
  if ( offset1->y < offset2->y ) {
    if ( dim1->y+offset1->y > offset2->y ) {
      origin.y = offset2->y;
      dim.y = dim1->y+offset1->y - offset2->y;
      if ( dim.y > dim2->y ) dim.y = dim2->y;
    }
  }
  else {
    if ( dim2->y+offset2->y > offset1->y) {
      origin.y = offset1->y;
      dim.y = dim2->y+offset2->y - offset1->y;
      if ( dim.y > dim1->y ) dim.y = dim1->y;
    }
  }

  /* intersection along Z
   */
  if ( offset1->z < offset2->z ) {
    if ( dim1->z+offset1->z > offset2->z ) {
      origin.z = offset2->z;
      dim.z = dim1->z+offset1->z - offset2->z;
      if ( dim.z > dim2->z ) dim.z = dim2->z;
    }
  }
  else {
    if ( dim2->z+offset2->z > offset1->z) {
      origin.z = offset1->z;
      dim.z = dim2->z+offset2->z - offset1->z;
      if ( dim.z > dim1->z ) dim.z = dim1->z;
    }
  }

  /* is there any intersection ?
     => ( dim.x > 0 && dim.y > 0 && dim.z > 0 )
   */

  *o = origin;
  *d = dim;

  return( dim.x * dim.y * dim.z );
}





int intersection_dimensions_images( typeMosaicPart *image1,
			     typeMosaicPart *image2,
			     vt_ipt *o,
			     vt_ipt *d )
{
  vt_ipt dim1, dim2;

  dim1.x = image1->image->dim.x;
  dim1.y = image1->image->dim.y;
  dim1.z = image1->image->dim.z;
  dim2.x = image2->image->dim.x;
  dim2.y = image2->image->dim.y;
  dim2.z = image2->image->dim.z;

  return( intersection_dimensions( &(image1->offset), &dim1,
				   &(image2->offset), &dim2, o, d ) );
}





int search_neighbors( typeMosaicPart *mosaic, int nb )
{
  char *proc = "search_neighbors";
  int n, i, j, v;
  vt_ipt origin;
  vt_ipt dim;

  for ( i=0; i<nb; i++ ) {
    for ( n=0; n<_MAX_NEIGHBORS_; n++ )
      _init_typeIntersection( &(mosaic[i].intersection[n]) );
    mosaic[i].n_neighbors = 0;
  }
  
  for ( i=0; i<nb-1; i++ )
  for ( j=i+1; j<nb; j++ ) {

    v = intersection_dimensions_images( &(mosaic[i]), &(mosaic[j]), &origin, &dim );

    if ( v < 0 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when calculating intersection between #%d and #%d\n",
		 proc, i, j );
      return( -1 );
    }

    if ( v == 0 ) continue;
    
    if ( 0 && _debug_ ) {
      fprintf( stderr, "%s: intersection between #%d and #%d at (%d,%d,%d) + [%d,%d,%d]\n",
	       proc, i , j, origin.x, origin.y, origin.z, dim.x, dim.y, dim.z );
    }
  
    if ( mosaic[i].n_neighbors >= mosaic[i].n_max_neighbors ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: can not add neighbor to image #%d\n", proc, i );
      return( -1 );
    }
    mosaic[i].intersection[mosaic[i].n_neighbors].neighbor = j;
    mosaic[i].intersection[mosaic[i].n_neighbors].origin = origin;
    mosaic[i].intersection[mosaic[i].n_neighbors].origin.x -= mosaic[i].offset.x;
    mosaic[i].intersection[mosaic[i].n_neighbors].origin.y -= mosaic[i].offset.y;
    mosaic[i].intersection[mosaic[i].n_neighbors].origin.z -= mosaic[i].offset.z;
    mosaic[i].intersection[mosaic[i].n_neighbors].dim = dim;
    mosaic[i].n_neighbors ++;

    if ( mosaic[j].n_neighbors >= mosaic[j].n_max_neighbors ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: can not add neighbor to image #%d\n", proc, j );
      return( -1 );
    }
    mosaic[j].intersection[mosaic[j].n_neighbors].neighbor = j;
    mosaic[j].intersection[mosaic[j].n_neighbors].origin = origin;
    mosaic[j].intersection[mosaic[j].n_neighbors].origin.x -= mosaic[j].offset.x;
    mosaic[j].intersection[mosaic[j].n_neighbors].origin.y -= mosaic[j].offset.y;
    mosaic[j].intersection[mosaic[j].n_neighbors].origin.z -= mosaic[j].offset.z;
    mosaic[j].intersection[mosaic[j].n_neighbors].dim = dim;
    mosaic[j].n_neighbors ++;

  }
  
  return( 1 );
}





int intersection_volume_with_margin( typeMosaicPart *image1,
				     typeMosaicPart *image2, int margin )
{
  vt_ipt origin;
  vt_ipt dim;
  vt_ipt offset1, offset2;
  vt_ipt dim1, dim2;

  offset1 = image1->offset;
  dim1.x = image1->image->dim.x;
  dim1.y = image1->image->dim.y;
  dim1.z = image1->image->dim.z;

  offset2 = image2->offset;
  dim2.x = image2->image->dim.x;
  dim2.y = image2->image->dim.y;
  dim2.z = image2->image->dim.z;

  if ( margin > 0 ) {
    offset1.x -= margin;
    offset1.y -= margin;
    offset1.z -= margin;
    dim1.x += margin;
    dim1.y += margin;
    dim1.z += margin;

    offset2.x -= margin;
    offset2.y -= margin;
    offset2.z -= margin;
    dim2.x += margin;
    dim2.y += margin;
    dim2.z += margin;
  }
  
  return( intersection_dimensions( &offset1, &dim1,
				   &offset2, &dim2, &origin, &dim ) );
}





int intersection_volume( typeMosaicPart *image1,
			 typeMosaicPart *image2 )
{
  vt_ipt origin;
  vt_ipt dim;
  
  return( intersection_dimensions_images( image1, image2, &origin, &dim ) );
}



/* (x,y,z) is in local coordinates w.r.t the image
 */
int _is_in_intersection( typeMosaicPart *image, int x, int y, int z )
{
  int n;
  
  for ( n=0; n<image->n_neighbors; n++ ) {
    if ( x >= image->intersection[n].origin.x 
	 && x < image->intersection[n].origin.x + image->intersection[n].dim.x
	 && y >= image->intersection[n].origin.y 
	 && y < image->intersection[n].origin.y + image->intersection[n].dim.y
	 && z >= image->intersection[n].origin.z 
	 && z < image->intersection[n].origin.z + image->intersection[n].dim.z )
      return( 1 );
  }
  return( 0 );
}


/* (x,y,z) is in global coordinates w.r.t the mosaic
 */
int _is_in_image( typeMosaicPart *image, int x, int y, int z )
{
  if ( x >= image->offset.x && x < image->offset.x + (int)image->image->dim.x
       && y >= image->offset.y && y < image->offset.y + (int)image->image->dim.y
       && z >= image->offset.z && z < image->offset.z + (int)image->image->dim.z )
    return( 1 );
  return( 0 );
}











/***************************************************
 *
 * image registration
 *
 ***************************************************/
    



typedef enum {
  _LINF_COLOR_,
  _L1_LUMIN_
} enumCriterion;


double criterion_value( typeMosaicPart *reference,
			typeMosaicPart *floating,
			int xoffset,
			int yoffset,
			int zoffset,
			int *nb )
{
  char *proc = "criterion_value";
  double error_value = 1000000.0;

  int n = 0;
  double c, d;

  int x, y, z, v;
  int dimv=reference->image->dim.v;

  vt_ipt refOffset = reference->offset;
  vt_ipt floOffset = floating->offset;
  vt_ipt refDim;
  vt_ipt floDim;
  vt_ipt intersectionOffset;
  vt_ipt intersectionDim;
  int intersectionVolume;
  
  *nb = 0;

  floOffset.x += xoffset;
  floOffset.y += yoffset;
  floOffset.z += zoffset;
  
  refDim.x = reference->image->dim.x;
  refDim.y = reference->image->dim.y;
  refDim.z = reference->image->dim.z;

  floDim.x = floating->image->dim.x;
  floDim.y = floating->image->dim.y;
  floDim.z = floating->image->dim.z;

  intersectionVolume = intersection_dimensions( &floOffset, &floDim, &refOffset, &refDim,
						&intersectionOffset, &intersectionDim );

  if ( 0 && _debug_ ) {
    fprintf( stderr, "%s: intersection = (%d %d %d) + [%d %d %d]\n", proc,
	     intersectionOffset.x, intersectionOffset.y, intersectionOffset.z, 
	     intersectionDim.x, intersectionDim.y, intersectionDim.z );
    fprintf( stderr, "\t floating  = (%d %d %d) + [%d %d %d]\n",
	     floOffset.x, floOffset.y, floOffset.z, floDim.x, floDim.y, floDim.z );
    fprintf( stderr, "\t reference  = (%d %d %d) + [%d %d %d]\n",
	     refOffset.x, refOffset.y, refOffset.z, refDim.x, refDim.y, refDim.z );
  }

  if ( intersectionVolume <= 0 ) {
    if ( _verbose_ && 0 )
      fprintf( stderr, "%s: no intersection", proc );
    return( error_value );
  }

  /* seems to be to ensure a large intersection
    if ( xmax-xmin > ymax-ymin ) {
      if ( xmax-xmin < 10 ) return( 1000000.0 );
    }
    else {
      if ( ymax-ymin < 10 ) return( 1000000.0 );
    }
   */

  
  if ( reference->image->dim.v != floating->image->dim.v ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different vectorial dimensions\n", proc );
    return( error_value );
  }

  switch( reference->image->type ) {
  default :
    return( error_value );
  case UCHAR :
    {
      unsigned char ***ref = (unsigned char ***)reference->image->array;

      switch( floating->image->type ) {
      default :
	return( error_value );
      case UCHAR :
	{
	  unsigned char ***flo = (unsigned char ***)floating->image->array;
	  c = 0.0;
	  n = 0;
	  for ( z=intersectionOffset.z; z<intersectionOffset.z+intersectionDim.z; z++ )
	  for ( y=intersectionOffset.y; y<intersectionOffset.y+intersectionDim.y; y++ )
	  for ( x=intersectionOffset.x; x<intersectionOffset.x+intersectionDim.x; x++ )
	  for ( v=0; v<dimv; v++ ) {
	    d = flo[z-floOffset.z][y-floOffset.y][(x-floOffset.x)*dimv + v];
	    d -= ref[z-refOffset.z][y-refOffset.y][(x-refOffset.x)*dimv + v];
	    c += d*d;
	    n ++;
	  }
	} /* UCHAR : floating->image->type */
	break;
      } /* switch( floating->image->type ) */

    } /* UCHAR : reference->image->type */
    break;
  } /* switch( reference->image->type ) */
  
  *nb = n;
  return( c );
}







typedef struct {
  int dx;
  int dy;
  int dz;
  double criterion;
  int vol;
} typeDisplacement;


typedef struct {
  int reference_index;
  typeMosaicPart *mosaic;
  int nb_images;
  typeDisplacement *displacements;
} _displacementParameters;


int _computeCriterionValues( void *parameter,
			     size_t first,
			     size_t last )
{
  _displacementParameters *p = (_displacementParameters *)parameter;
  int reference_index = p->reference_index;
  typeMosaicPart *mosaic = p->mosaic;
  int nb_images = p->nb_images;
  typeDisplacement *displacements = p->displacements;

  double criterion;
  int i, v, vol;
  int x, y, z;

  size_t n;

  for ( n=first; n<=last; n ++ ) {

    displacements[n].criterion = 0.0;
    displacements[n].vol = 0;
    
    x = displacements[n].dx;
    y = displacements[n].dy;
    z = displacements[n].dz;

    for ( vol=0, criterion=0.0, i=0; i<nb_images; i++ ) {
      if ( i == reference_index ) continue;
      if ( mosaic[i].flag == 0 ) continue;
      if ( intersection_volume( &mosaic[i], &mosaic[reference_index] ) == 0 ) continue;
      criterion += criterion_value( &mosaic[i], &mosaic[reference_index], x, y, z, &v );
      vol += v;
    }
   
    displacements[n].criterion = criterion;
    displacements[n].vol = vol;
    if ( displacements[n].vol > 0 ) displacements[n].criterion /= displacements[n].vol;
  }

  return( 1 );
}








void build_mosaic( typeMosaicPart *mosaic, int nb_images, 
		   int txmin, int txmax, int tymin, int tymax )
{
  char *proc = "build_mosaic";
  int i, j, c, nb;
  
  vt_ipt correction;
  int n;
  
  vt_ipt best_displacement;
  double criterion, best_criterion;
  int vol;

  int x, y, z;
  int tzmin = 0;
  int tzmax = 0;
  int ndisplacements;

  typeDisplacement *displacements = NULL;
  typeChunks chunks;
  _displacementParameters p;


  /* parallelism stuff 
   */

  ndisplacements = (tzmax-tzmin+1) * (tymax-tymin+1) * (txmax-txmin+1);
  if ( ndisplacements <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: no displacements to be tested ?!\n", proc );
    return;
  }
  displacements = (typeDisplacement *)malloc( ndisplacements * sizeof(typeDisplacement) );

  p.reference_index = -1;
  p.mosaic = mosaic;
  p.nb_images = nb_images;
  p.displacements = displacements;

  initChunks( &chunks );
  if ( buildChunks( &chunks, 0, ndisplacements-1, proc ) != 1 ) {
    free( displacements );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }
  
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);









  for (i=0; i<nb_images; i++ )
    mosaic[i].flag = 0;

  mosaic[0].flag = 1;




  for (nb=nb_images-1; nb>0; nb-- ) {
    
    for ( i=0; i<nb_images; i++ ) {
      if ( mosaic[i].flag != 0 ) continue;
      mosaic[i].sum = 0;
    }

    /* potential candidate is first unlabeled image
     */
    for ( c=-1, i=0; i<nb_images && c==-1; i++ ) {
      if ( mosaic[i].flag != 0 ) continue;
      c = i;
    }

    /* choice of a better candidate, the one that maximizes the intersection
       with already processed ones
    */
    for ( i=0; i<nb_images;i++ ) {
      if ( mosaic[i].flag != 0 ) continue;
      mosaic[i].sum = 0;
      
      /* compute volume of intersection of #i
	 with  labeled images
      */
      for (j=0; j<nb_images; j++ ) {
	if ( mosaic[j].flag == 0 ) continue;
	
	if ( 0 && _debug_ ) {
	  fprintf( stderr, "%s: intersection of %2d and %2d = %d\n", proc, i, j,
		   intersection_volume( &mosaic[i], &mosaic[j] ) );
	}
	mosaic[i].sum += intersection_volume( &mosaic[i], &mosaic[j] );
      }
      if ( mosaic[i].sum > mosaic[c].sum ) c = i;
    }
    
    
    if ( 1 && _debug_ ) {
      fprintf( stderr, "%s: image %3d/%3d = #%3d, intersection volume = %lf\n", proc, nb_images-nb, nb_images, c, mosaic[c].sum );
    }



    /* case where there is no intersection with labeled images
       do nothing
     */
    if (  mosaic[c].sum == 0 ) {
      mosaic[c].offset = mosaic[c].offset_init;
      mosaic[c].flag = 1;
      continue;
    }



    /* correction d'apres le deplacement des images
       deja repositionnees
       
    */
    correction.x = correction.y = correction.z = 0;

    for ( n=0, i=0; i<nb_images; i++ ) {
      if ( i == c ) continue;
      if ( mosaic[i].flag == 0 ) continue;
      if ( intersection_volume( &mosaic[i], &mosaic[c] ) == 0 ) continue;
      correction.x += mosaic[i].offset.x - mosaic[i].offset_init.x;
      correction.y += mosaic[i].offset.y - mosaic[i].offset_init.y;
      correction.z += mosaic[i].offset.z - mosaic[i].offset_init.z;
      n ++;
    }

    if ( n > 0 ) {
      correction.x /= n;
      correction.y /= n;
      correction.z /= n;
    }

     if ( 1 && _debug_ ) {
       fprintf( stderr, "\t initial correction for #%3d is +[%d %d %d] from %d images \n", 
		c, correction.x, correction.y, correction.z, n );
    }
   

     /* initial value
	valeur du critere pour un deplacement de (0,0,0)
      */
    best_displacement.x = best_displacement.y = best_displacement.z = 0;

    for ( vol=0, best_criterion=0.0, i=0; i<nb_images; i++ ) {
      if ( i == c ) continue;
      if ( mosaic[i].flag == 0 ) continue;
      if ( intersection_volume( &mosaic[i], &mosaic[c] ) == 0 ) continue;
      best_criterion += criterion_value( &mosaic[i], &mosaic[c], 
					 best_displacement.x, 
					 best_displacement.y,
					 best_displacement.z, &n );
      vol += n;
      if ( 0 )
	fprintf( stderr, "image #%d: intersection %3d/%3d, criterion=%f, vol=%d (n=%d)\n",
		 c, i, nb_images, best_criterion, vol, n );
    }
    if ( vol > 0 ) 
      best_criterion /= (double)vol;
    else 
      fprintf( stderr, "%s: unable to compute first criterion value \n", proc );
    
    if ( _debug_ ) {
       fprintf( stderr, "\t initial criterion value is %lf\n", best_criterion );
    }

    /* criterion computation
     */
    p.reference_index = c;
    for ( n=0, z=correction.z+tzmin; z<=correction.z+tzmax; z++ )
    for ( y=correction.y+tymin; y<=correction.y+tymax; y++ )
    for ( x=correction.x+txmin; x<=correction.x+txmax; x++, n++ ) {
      displacements[n].dx = x;
      displacements[n].dy = y;
      displacements[n].dz = z;
    }
    
    if ( processChunks( &_computeCriterionValues, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute criterions\n", proc );
      free( displacements );
      freeChunks( &chunks );
      return;
    }

    for ( n=0; n<ndisplacements; n++ ) {
      if ( displacements[n].vol > 0 ) {
	if ( 0 && _debug_ ) 
	  fprintf( stderr, "test %d %d %d: criterion = %f <=> %f\n",
		   x, y, z, criterion, best_criterion );
	if ( displacements[n].criterion < best_criterion ) {
	  best_criterion = displacements[n].criterion;
	  best_displacement.x = displacements[n].dx;
	  best_displacement.y = displacements[n].dy;
	  best_displacement.z = displacements[n].dz;
	}
      }
    }
    
    if ( _debug_ ) {
      fprintf( stderr, "\t final criterion value is %lf\n", best_criterion );
      fprintf( stderr, "\t final correction for #%3d is +[%d %d %d]\n", 
	       c, best_displacement.x, best_displacement.y, best_displacement.z );
    }

    mosaic[c].offset.x = mosaic[c].offset_init.x + best_displacement.x;
    mosaic[c].offset.y = mosaic[c].offset_init.y + best_displacement.y;
    mosaic[c].offset.z = mosaic[c].offset_init.z + best_displacement.z;
    mosaic[c].flag = 1;
    

  }


  free( displacements );
  freeChunks( &chunks );
}






double old_criterion_value( typeMosaicPart *reference,
			typeMosaicPart *floating,
			int xoffset,
			int yoffset,
			int zoffset,
			int *nb )
{
  double error_value = 1000000.0;
  int x, y;
  int n = 0;
  int k, l[3];
  double c = 0.0;;
  enumCriterion criterion = _L1_LUMIN_;

  int tx = floating->offset_init.x + xoffset - reference->offset.x;
  int ty = floating->offset_init.y + yoffset - reference->offset.y;
  int tz = floating->offset_init.z + zoffset - reference->offset.z;

  int xmin=tx;
  int xmax=tx+floating->image->dim.x-1;
  int ymin=ty;
  int ymax=ty+floating->image->dim.y-1;
  int zmin=tz;
  int zmax=tz+floating->image->dim.z-1;

  int lr, lf;

  *nb = 0;

  if ( xmin < 0 ) xmin = 0;
  if ( xmin >= reference->image->dim.x ) return( error_value );
  if ( xmax < 0 ) return( error_value );
  if ( xmax >= reference->image->dim.x ) xmax = reference->image->dim.x-1;

  if ( ymin < 0 ) ymin = 0;
  if ( ymin >= reference->image->dim.y ) return( error_value );
  if ( ymax < 0 ) return( error_value );
  if ( ymax >= reference->image->dim.y ) ymax = reference->image->dim.y-1;

  if ( zmin < 0 ) zmin = 0;
  if ( zmin >= reference->image->dim.z ) return( error_value );
  if ( zmax < 0 ) return( error_value );
  if ( zmax >= reference->image->dim.z ) zmax = reference->image->dim.z-1;

  /* seems to be to ensure a large intersection
   */
  if ( 0 ) {
    if ( xmax-xmin > ymax-ymin ) {
      if ( xmax-xmin < 10 ) return( 1000000.0 );
    }
    else {
      if ( ymax-ymin < 10 ) return( 1000000.0 );
    }
  }

  switch( reference->image->type ) {
  default :
    return( 1000000.0 );
  case UCHAR :
    {
      unsigned char ***ref = (unsigned char ***)reference->image->array;
      switch( floating->image->type ) {
      default :
	return( 1000000.0 );
      case UCHAR :
	{
	  unsigned char ***flo = (unsigned char ***)floating->image->array;

	    switch ( criterion ) {

	    case _L1_LUMIN_ :
	       if ( floating->image->dim.v == 3 && reference->image->dim.v == 3 ) {
		 for ( y=ymin; y<=ymax; y += 10 )
		 for ( x=xmin; x<=xmax; x += 10 ) {
		   lr = (int)(0.2125 * (double)ref[0][y][3*x] +
			      0.7154 * (double)ref[0][y][3*x+1] +
			      0.0721 * (double)ref[0][y][3*x+2] + 0.5 );
		   lf = (int)(0.2125 * (double)flo[0][y-ty][3*(x-tx)] +
			      0.7154 * (double)flo[0][y-ty][3*(x-tx)+1] +
			      0.0721 * (double)flo[0][y-ty][3*(x-tx)+2] + 0.5 );
		   if ( lr > lf ) c += lr-lf;
		   else           c += lf-lr;
		   n++;
		 }
	       }
	       break;
	       
	    case _LINF_COLOR_ :
	      if ( floating->image->dim.v == 3 && reference->image->dim.v == 3 ) {
		for ( y=ymin; y<=ymax; y += 5 )
		for ( x=xmin; x<=xmax; x += 5 ) {
		  
		  for (k=0;k<3;k++) {
		    l[k] = (int)flo[0][y-ty][3*(x-tx)+k] - (int)ref[0][y][3*x+k];
		    if ( l[k] < 0 ) l[k] = (-l[k]);
		  }
		  if (l[0] > l[1] && l[0] > l[2] ) {
		    c += l[0];
		  }
		  else {
		    if ( l[1] > l[2] ) c += l[1];
		    else c += l[2];
		  }
		  n++;
		}
	      }
	      break;

	    } /* end switch ( criterion ) */
	}
	break;
      }
    }
    break;
  }
  
  *nb = n;
  if ( n > 0 ) c /= n;
  return( c );
}











/***************************************************
 *
 * image composition
 *
 ***************************************************/

typedef enum {
  _UNDEF_ = 10,
  _RED_ = 0,
  _GREEN_ = 1,
  _BLUE_ = 2
} enumColor;


/* on cree une grande image ou chacune des imagettes
   a une couleur (R, G ou B) pour bien les diffferencier
*/

vt_image * build_color_mosaic( char *name, typeMosaicPart *mosaic, int nb_images )
{
  char *proc = "build_color_mosaic";
  vt_image *image;
  int i, j;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  int dimv, dimx, dimy, dimz;
  int v, lx, ly, lz;
  int gx, gy, gz;

  int margin;
  int inc_margin = 10;
  int max_margin = 100;
  int color = _RED_; 

  

  unsigned char ***res;
  int nb, c, weight[3];
  
  xmin = mosaic[0].offset.x;
  ymin = mosaic[0].offset.y;
  zmin = mosaic[0].offset.z;
  
  xmax = mosaic[0].offset.x + mosaic[0].image->dim.x - 1;
  ymax = mosaic[0].offset.y + mosaic[0].image->dim.y - 1;
  zmax = mosaic[0].offset.z + mosaic[0].image->dim.z - 1;

  for (i=1; i<nb_images; i++ ) {

    if ( mosaic[i].image->dim.v != mosaic[0].image->dim.v ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: image #%d has different vectorial dimension than #0\n", proc, i );
      return( (vt_image*)NULL );
    }
    if ( mosaic[i].image->type != mosaic[0].image->type ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: image #%d has different type than #0\n", proc, i );
      return( (vt_image*)NULL );
    }

    if ( xmin > mosaic[i].offset.x )  xmin =  mosaic[i].offset.x;
    if ( ymin > mosaic[i].offset.y )  ymin =  mosaic[i].offset.y;
    if ( zmin > mosaic[i].offset.z )  zmin =  mosaic[i].offset.z;
    if ( xmax < mosaic[i].offset.x + (int)mosaic[i].image->dim.x - 1 ) 
      xmax = mosaic[i].offset.x + (int)mosaic[i].image->dim.x - 1;
    if ( ymax < mosaic[i].offset.y + (int)mosaic[i].image->dim.y - 1 ) 
      ymax = mosaic[i].offset.y + (int)mosaic[i].image->dim.y - 1;
    if ( zmax < mosaic[i].offset.z + (int)mosaic[i].image->dim.z - 1 ) 
      zmax = mosaic[i].offset.z + (int)mosaic[i].image->dim.z - 1;
  }

  dimv = 3;
  dimx = xmax-xmin+1;
  dimy = ymax-ymin+1;
  dimz = zmax-zmin+1;

  if ( _verbose_ ) {
    fprintf( stdout, "%s: min=(%d,%d,%d), max=(%d,%d,%d) dim=[%d %d %d]\n", 
	     proc, xmin, ymin, zmin, xmax, ymax, zmax, dimx, dimy, dimz );
  }
  
  if ( mosaic[0].image->dim.v != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: vectorial images not handled yet\n", proc );
    return( (vt_image*)NULL );
  }
  

  image = (vt_image *)malloc( sizeof(vt_image) );

  VT_Image( image );
  VT_InitVImage( image, name, dimv, dimx, dimy, dimz, mosaic[0].image->type );
  image->siz = mosaic[0].image->siz;

  VT_AllocImage( image );
  res = (unsigned char ***)image->array;
  
  

  /* attribution des couleurs
   */
  for (i=0; i<nb_images; i++ ) {
    mosaic[i].flag = _UNDEF_;
    mosaic[i].sum = 0;
  }
  
  mosaic[0].flag = color;
  if ( color == _RED_ ) color = _GREEN_;
  else if ( color == _GREEN_ ) color = _BLUE_;
  else if ( color == _BLUE_ ) color = _RED_;



  for (nb=nb_images-1; nb>0; nb-- ) {

    /* initialization required because of loop on margin
     */
    for ( i=0; i<nb_images; i++ ) {
      if ( mosaic[i].flag != _UNDEF_ ) continue;
      mosaic[i].sum = 0;
    }

    /* potential candidate is first unlabeled image
     */
    for ( c=-1, i=0; i<nb_images && c==-1; i++ ) {
      if ( mosaic[i].flag != _UNDEF_ ) continue;
      c = i;
    }

    /* choice of a better candidate, the one that maximizes the intersection
       with already processed ones
     */
    for ( margin=0; margin<max_margin && mosaic[c].sum == 0; margin+=inc_margin ) {

      for ( i=0; i<nb_images;i++ ) {
	if ( mosaic[i].flag != _UNDEF_ ) continue;
	mosaic[i].sum = 0;

	/* compute volume of intersection of #i
	   with  labeled images
	*/
	for (j=0; j<nb_images; j++ ) {
	  if ( mosaic[j].flag == _UNDEF_ ) continue;
	  
	  if ( 0 && _debug_ ) {
	    fprintf( stderr, "%s: intersection of %2d and %2d with margin of %2d = %d\n", proc, i, j, margin,
		     intersection_volume_with_margin( &mosaic[i], &mosaic[j], margin ) );
	  }
	  mosaic[i].sum += intersection_volume_with_margin( &mosaic[i], &mosaic[j], margin );
	}
	if ( mosaic[i].sum > mosaic[c].sum ) c = i;
      }
    }
    
    /* no intersection with others ...
     */
    if ( mosaic[c].sum == 0 ) {
      mosaic[c].flag = color;
      if ( color == _RED_ ) color = _GREEN_;
      else if ( color == _GREEN_ ) color = _BLUE_;
      else if ( color == _BLUE_ ) color = _RED_;
      continue;
    }

    if ( 0 && _debug_ ) {
      fprintf( stderr, "%s: image de plus grande intersection = %d\n", proc, c );
    }
    
    for ( j=0; j<3; j++ ) weight[j] = 0;
    for ( i=0; i<nb_images;i++ ) {
      if ( mosaic[i].flag == _UNDEF_ ) continue;
      weight[ mosaic[i].flag ] += intersection_volume_with_margin( &mosaic[i], &mosaic[c], margin );
    }

    if ( weight[ _RED_ ] <= weight[ _GREEN_ ] ) {
      if ( weight[ _RED_ ] <= weight[ _BLUE_ ] ) {
	mosaic[c].flag = _RED_;
      }
      else {
	mosaic[c].flag = _BLUE_;
      }
    }
    else {
      if ( weight[ _GREEN_ ] <= weight[ _BLUE_ ] ) {
	mosaic[c].flag = _GREEN_;
      }
      else  {
	mosaic[c].flag = _BLUE_;
      }
    }

  }


  
  if ( 0 && _debug_ ) {
    for (i=0; i<nb_images; i++ ) {
      fprintf( stderr, "#%3d: ", i );
      switch( mosaic[i].flag ) {
      case _RED_ : fprintf( stderr, "red\n" ); break;
      case _GREEN_ : fprintf( stderr, "green\n" ); break;	
      case _BLUE_ : fprintf( stderr, "blue\n" ); break;
      case _UNDEF_ : fprintf( stderr, "undef\n" ); break;
      default : fprintf( stderr, "default\n" ); break;
      }
    }
  }


  switch ( mosaic[0].image->type) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    VT_FreeImage( image );
    free( image );
    return( (vt_image*)NULL );

  case UCHAR :
    {
      unsigned char ***buf = (unsigned char***)NULL;
      
      for ( gz=0; gz<dimz; gz++ )
      for ( gy=0; gy<dimy; gy++ )
      for ( gx=0; gx<dimx; gx++ )
      for ( v=0; v<dimv; v++ )
	res[gz][gy][gx*dimv+v] = 0;

      for ( i=0; i<nb_images; i++ ) {
	buf = (unsigned char***)mosaic[i].image->array;
	for ( lz=0, gz=mosaic[i].offset.z-zmin; lz<mosaic[i].image->dim.z; lz++, gz++ )
	for ( ly=0, gy=mosaic[i].offset.y-ymin; ly<mosaic[i].image->dim.y; ly++, gy++ )
	for ( lx=0, gx=mosaic[i].offset.x-xmin; lx<mosaic[i].image->dim.x; lx++, gx++ ) {
	  switch ( mosaic[i].flag ) {
	  default :
	  case _RED_ :
	    res[gz][gy][gx*dimv]   = buf[lz][ly][lx]; break;
	  case _GREEN_ :
	    res[gz][gy][gx*dimv+1] = buf[lz][ly][lx]; break;
	  case _BLUE_ :
	    res[gz][gy][gx*dimv+2] = buf[lz][ly][lx]; break;
	  }
	}
      }
    }
    break;
  }
  return( image );
}





vt_image * build_large_image( char *name, typeMosaicPart *mosaic, int nb_images )
{
  char *proc = "build_large_image";
  vt_image *image;
  int i, j;
  int xmin, xmax, ymin, ymax, zmin, zmax;
  int dimv, dimx, dimy, dimz;
  int v, lx, ly, lz;
  int gx, gy, gz;
  int xj, yj, zj, wj;

  double sumvalues;
  double sumweights;

  xmin = mosaic[0].offset.x;
  ymin = mosaic[0].offset.y;
  zmin = mosaic[0].offset.z;

  xmax = mosaic[0].offset.x + mosaic[0].image->dim.x - 1;
  ymax = mosaic[0].offset.y + mosaic[0].image->dim.y - 1;
  zmax = mosaic[0].offset.z + mosaic[0].image->dim.z - 1;

  for ( i=1; i<nb_images; i++ ) {
    
    if ( mosaic[i].image->dim.v != mosaic[0].image->dim.v ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: image #%d has different vectorial dimension than #0\n", proc, i );
      return( (vt_image*)NULL );
    }
    if ( mosaic[i].image->type != mosaic[0].image->type ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: image #%d has different type than #0\n", proc, i );
      return( (vt_image*)NULL );
    }

    if ( xmin > mosaic[i].offset.x )  xmin =  mosaic[i].offset.x;
    if ( ymin > mosaic[i].offset.y )  ymin =  mosaic[i].offset.y;
    if ( zmin > mosaic[i].offset.z )  zmin =  mosaic[i].offset.z;
    if ( xmax < mosaic[i].offset.x + (int)mosaic[i].image->dim.x - 1 ) 
      xmax = mosaic[i].offset.x + (int)mosaic[i].image->dim.x - 1;
    if ( ymax < mosaic[i].offset.y + (int)mosaic[i].image->dim.y - 1 ) 
      ymax = mosaic[i].offset.y + (int)mosaic[i].image->dim.y - 1;
    if ( zmax < mosaic[i].offset.z + (int)mosaic[i].image->dim.z - 1 ) 
      zmax = mosaic[i].offset.z + (int)mosaic[i].image->dim.z - 1;
  }




  dimv = mosaic[0].image->dim.v;
  dimx = xmax-xmin+1;
  dimy = ymax-ymin+1;
  dimz = zmax-zmin+1;


  if ( _verbose_ ) {
    fprintf( stdout, "\n" );
    for ( i=0; i<nb_images; i++ ) {
      fprintf( stdout, "#%03d ", i );
      print_mosaic_part_info( stdout, &(mosaic[i]), 0 );
    }
    fprintf( stdout, "\n" );
    fprintf( stdout, "%s: min=(%d,%d,%d), max=(%d,%d,%d) dim=[%d %d %d]\n", 
	     proc, xmin, ymin, zmin, xmax, ymax, zmax, dimx, dimy, dimz );
    fprintf( stdout, "\n" );
  }

  image = (vt_image *)malloc( sizeof(vt_image) );

  VT_Image( image );
  VT_InitVImage( image, name, mosaic[0].image->dim.v, dimx, dimy, dimz, mosaic[0].image->type );
  image->siz = mosaic[0].image->siz;

  VT_AllocImage( image );



  switch ( image->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    VT_FreeImage( image );
    free( image );
    return( (vt_image*)NULL );

  case UCHAR :
    {
      unsigned char ***res = (unsigned char***)image->array;
      unsigned char ***buf = (unsigned char***)NULL;
      unsigned char value;

      
      for ( gz=0; gz<dimz; gz++ )
      for ( gy=0; gy<dimy; gy++ )
      for ( gx=0; gx<dimx; gx++ )
      for ( v=0; v<dimv; v++ )
	res[gz][gy][gx*dimv+v] = 0;
      
      for ( i=0; i<nb_images; i++ ) {
	buf = (unsigned char***)mosaic[i].image->array;
	for ( lz=0, gz=mosaic[i].offset.z-zmin; lz<mosaic[i].image->dim.z; lz++, gz++ )
	for ( ly=0, gy=mosaic[i].offset.y-ymin; ly<mosaic[i].image->dim.y; ly++, gy++ )
        for ( lx=0, gx=mosaic[i].offset.x-xmin; lx<mosaic[i].image->dim.x; lx++, gx++ ) {
	  if ( _is_in_intersection( &(mosaic[i]), lx, ly, lz ) ) {
	    switch ( fusionMode ) {
	    default :
	      if ( _verbose_ )
		fprintf( stderr, "%s: such fusion mode not handled yet\n", proc );
	      VT_FreeImage( image );
	      free( image );
	      return( (vt_image*)NULL );

	    case MAXIMUM :
	      for ( v=0; v<dimv; v++ ) {
		res[gz][gy][gx*dimv+v] = buf[lz][ly][lx*dimv+v];
		for ( j=0; j<nb_images; j++ ) {
		  if ( _is_in_image( &(mosaic[j]), gx+xmin, gy+ymin, gz+zmin ) == 0 )
		    continue;
		  value = ((unsigned char***)mosaic[j].image->array)[gz+zmin-mosaic[j].offset.z][gy+ymin-mosaic[j].offset.y][(gx+xmin-mosaic[j].offset.x)*dimv+v];
		  if ( res[gz][gy][gx*dimv+v] < value ) res[gz][gy][gx*dimv+v] = value;
		}
	      }
	      break;

	    case AVERAGE :
	      for ( v=0; v<dimv; v++ ) {
		sumvalues = 0.0;
		sumweights = 0.0;
		for ( j=0; j<nb_images; j++ ) {
		  if ( _is_in_image( &(mosaic[j]), gx+xmin, gy+ymin, gz+zmin ) == 0 )
		    continue;
		  sumvalues += ((unsigned char***)mosaic[j].image->array)[gz+zmin-mosaic[j].offset.z][gy+ymin-mosaic[j].offset.y][(gx+xmin-mosaic[j].offset.x)*dimv+v];
		  sumweights += 1;
		}
		res[gz][gy][gx*dimv+v] = (int)( sumvalues / sumweights + 0.5 );
	      }
	      break;

	    case WEIGHTING :
	      for ( v=0; v<dimv; v++ ) {
		sumvalues = 0.0;
		sumweights = 0.0;
		for ( j=0; j<nb_images; j++ ) {
		  if ( _is_in_image( &(mosaic[j]), gx+xmin, gy+ymin, gz+zmin ) == 0 )
		    continue;
		  xj = gx+xmin-mosaic[j].offset.x;
		  yj = gy+ymin-mosaic[j].offset.y;
		  zj = gz+zmin-mosaic[j].offset.z;
		  wj = xj + 1;
		  if ( wj > mosaic[j].image->dim.x - xj ) wj = mosaic[j].image->dim.x - xj;
		  if ( wj > yj + 1 ) wj = yj + 1;
		  if ( wj > mosaic[j].image->dim.y - yj ) wj = mosaic[j].image->dim.y - yj;
		  if ( mosaic[j].image->dim.z > 1 ) {
		    if ( wj > zj + 1 ) wj = zj + 1;
		    if ( wj > mosaic[j].image->dim.z - zj ) wj = mosaic[j].image->dim.z - zj;
		  }
		  sumvalues += wj * ((unsigned char***)mosaic[j].image->array)[zj][yj][(xj)*dimv+v];
		  sumweights += wj;
		}
		res[gz][gy][gx*dimv+v] = (int)( sumvalues / sumweights + 0.5 );
	      }
	      break;

	    } /* end switch ( fusionMode ) */

	  }
	  /* the point does not belong to an intersection
	     -> just copy
	   */
	  else {
	    for ( v=0; v<dimv; v++ ) 
	      res[gz][gy][gx*dimv+v] = buf[lz][ly][lx*dimv+v];
	  }
	}
      } /* end for ( i=0; i<nb_images; i++ ) */

    }
    break;
    
  }
  

  return( image );
}

