/*************************************************************************
 * tracking-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 25 jui 2013 16:52:38 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include <chunks.h>
#include <histogram.h>

#include <vt_inrimage.h>

#include <tracking-tools.h>

static int _verbose_ = 1;
static int _debug_ = 0;





/************************************************************
 *
 * MANAGEMENT (valueType)
 *
 ************************************************************/

void initValueType( valueType *v )
{
  v->data = (int *)NULL;
  v->n_data = 0;
  v->center = 0;
  v->median = 0;
  v->min = 0;
  v->max = 0;
  v->mean = 0;
}



int allocValueType( valueType *v, int n ) 
{
  char *proc = "allocValueType";
  int i;
  
  if ( n <= 0 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: null or negative size\n", proc );
    return( -1 );
  }
  
  v->data = (int*)malloc( n * sizeof(int) );
  if ( v->data == (int*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }
  
  for ( i=0; i<n; i++ )
    v->data[i] = 0;

  return( 1 );
}



void freeValueType( valueType *v ) 
{
  if ( v->data != (int*)NULL )
    free( v->data );
  initValueType( v );
}



void computeValueType( valueType *v ) 
{
  int i, j, tmp;
  int right, last, left, med=(int)(v->n_data / 2);
  
  for ( v->min=v->data[0], i=1; i<v->n_data; i++ )
    if ( v->min > v->data[i] ) v->min = v->data[i];

  for ( v->max=v->data[0], i=1; i<v->n_data; i++ )
    if ( v->max < v->data[i] ) v->max = v->data[i];

  for ( v->mean=0.0, i=0; i<v->n_data; i++ ) v->mean += v->data[i];
  v->mean /= (float)v->n_data;

  left = 0; right = v->n_data-1;
  do {
    /* swap left et (left+right)/2 */
    j = (left+right)/2;						
    tmp = v->data[left];   v->data[left] = v->data[j];   v->data[j] = tmp;		
    /* cut v->data into two */							
    last = left;								
    for ( i = left+1; i <= right; i++ ) {					
      if ( v->data[i] < v->data[left] ) {						
	last ++;								
	tmp = v->data[i];   v->data[i] = v->data[last];   v->data[last] = tmp;		
      }									
    }									
    tmp = v->data[left];   v->data[left] = v->data[last]; v->data[last] = tmp;		
    if ( last >  med ) right = last - 1;				
    if ( last <  med ) left  = last + 1;
  } while ( last != med );
  v->median = v->data[med]; 
}





/************************************************************
 *
 * MANAGEMENT (circleType)
 *
 ************************************************************/

void initCircleType( circleType *c )
{
  c->x = 0;
  c->y = 0;
  c->r = 0;
  c->image = 0;

  initValueType( &(c->intensity) );

  c->index = 0;

  c->type = _UNKNOWN_TYPE_POINT_;

  c->prev = (circleType **)NULL;
  c->n_prev = 0;
  c->n_prev_allocated = 0;

  c->next = (circleType **)NULL;
  c->n_next = 0;
  c->n_next_allocated = 0;

  c->n_chains = 0;
  c->n_colocs = 0;
}



void freeCircleType( circleType *c )
{
  freeValueType( &(c->intensity) );

  if ( c->prev != (circleType **)NULL ) 
    free( c->prev );
  if ( c->next != (circleType **)NULL ) 
    free( c->next );

  initCircleType( c );
}



static void printCircleType( FILE *f, circleType *c )
{
  fprintf( f, "%-7d %-7d %-7d\n", c->x, c->y, c->r );
}










/************************************************************
 *
 * MANAGEMENT (circleList)
 *
 ************************************************************/

void initCircleList( circleList *l )
{
  l->data = (circleType*)NULL;
  l->n_data = 0;
  l->n_selected_data = 0;
  l->n_allocated_data = 0;
  l->rmax = 0;
} 



void freeCircleList( circleList *l )
{
  int i;
  
  if ( l->data != (circleType*)NULL ) {
    for ( i=0; i<l->n_data; i++ )
      freeCircleType( &(l->data[i]) );
    free( l->data );
  } 
  initCircleList( l );
} 



void printCircleList( FILE *f, circleList *l )
{
  int i;
  
  for ( i=0; i<l->n_selected_data; i++ )
    printCircleType( f, &(l->data[i]) );
}






static int _size_to_be_allocated_ = 100;

static int addCircleToCircleList( circleList *l, circleType *c )
{
  char *proc = "addCircleToCircleList";
  int s =  l->n_allocated_data;
  circleType *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (circleType*)malloc( s * sizeof(circleType) );
    if ( data == (circleType*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(circleType) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  initCircleType( &(l->data[l->n_data]) );
  l->data[l->n_data] = *c;
  l->n_data ++;
  l->n_selected_data ++;
  return( 1 );
} 





int readCircleList( circleList *l, char *filename )
{
  char *proc = "readCircleList";
  FILE *f;
  int i, nread;
  int x, y, r, v;
  circleType c;
  int rmax = 0;

  f = fopen( filename, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
    return( -1 );
  }

  i = 0;
  do {
    nread = fscanf( f, "%d %d %d %d", &x, &y, &r, &v );

    initCircleType( &c );
    switch ( nread ) {
    case 4 :
      c.intensity.mean = v;
    case 3 :
      c.x = x;
      c.y = y;
      c.r = r;
      break;
    default :
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: weird number of entries (%d)\n", proc, nread );
      return( -1 );
    }

    c.type = _READ_;
    if ( addCircleToCircleList( l, &c ) != 1 ) {
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add measure #%d (=%d %d %d) to list\n", proc, i, x, y, r );
      return( -1 );
    }
    i ++;
  } while ( nread == 3 || nread == 4 );
  
  if ( 0 && _verbose_ )
    fprintf( stderr, "%s: read %d values\n", proc, i );

  for ( i=0; i<l->n_data; i++ ) 
    if ( rmax < l->data[i].r )
      rmax = l->data[i].r;

  l->rmax = rmax;

  fclose( f );
  return( 1 );
}










/************************************************************
 *
 * MANAGEMENT (circleListList)
 *
 ************************************************************/

void initCircleListList( circleListList *l )
{
  l->data = (circleList*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
  l->rmax = 0;
} 



int allocCircleListList( circleListList *l, int n )
{
  char *proc = "allocLCircleistList";
  int i;
  
  if ( n <= 0 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: null or negative size\n", proc );
    return( -1 );
  }

  l->data = (circleList *)malloc( n * sizeof( circleList ) );
  if ( l->data == (circleList*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  } 

  l->n_data = l->n_allocated_data = n;
  
  for ( i=0; i<n; i++ )
    initCircleList( &(l->data[i]) );

  return( 1 );
} 



static circleList *getCircleListFromCircleListList( circleListList *l )
{
  char *proc = "getCircleListFromCircleListList";
  int s =  l->n_allocated_data;
  circleList *data;
  int i;
  
  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (circleList*)malloc( s * sizeof(circleList) );
    if ( data == (circleList*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( (circleList *)NULL );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(circleList) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
      initCircleList( &(l->data[i]) );
  }
  l->n_data ++;
  return( &(l->data[l->n_data-1]) );
} 



void freeCircleListList( circleListList *l )
{
  char *proc = "freeCircleListList";
  int i;

  if ( l->data != (circleList*)NULL ) {
    for ( i=0; i<l->n_data; i++ ) {
      if ( _debug_ )
	fprintf( stderr, "%s: freeing list #%5d/%5d\n", proc, i, l->n_data-1 );
      freeCircleList( &(l->data[i]) );
    } 
    free( l->data );
  } 
  initCircleListList( l );
} 





int readCircleListList( circleListList *readlist, char *format, int first, int last )
{
  char *proc = "readCircleListList";
  char filename[1024];
  int nlist;
  int i, j, k;

  /* allocation de la structure pour lecture
   */
  nlist = last - first + 1;
  if ( nlist < 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: first=%d > last=%d\n", proc, first, last );
    return( -1 );
  }

  if ( allocCircleListList( readlist, nlist ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  /* lecture des fichiers #i
   */
  for ( j=0, i=first; i<=last; i++, j++ ) {

    sprintf( filename, format, i );
    if ( _debug_ )
      fprintf( stderr, "#%4d: processing '%s'\n", i, filename );
    
    if ( readCircleList( &(readlist->data[j]), filename ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "error while reading file '%s'\n", filename );
      continue;
    }
    
    for ( k=0; k<readlist->data[j].n_data; k++ ) {
      readlist->data[j].data[k].index = j;
      readlist->data[j].data[k].image = i;
    }
  }


  return( 1 );

}




void printTrackedCircleList( FILE *f, circleListList *l )
{
  int i, j;
  for ( j=0; j<l->n_data; j++ ) {
    fprintf( f, "track %d\n", j );
    for ( i=0; i<l->data[j].n_data; i++ ) {
      fprintf( f, "%d\t%d\t%d %d", 
	       l->data[j].data[i].x,
	       l->data[j].data[i].y,
	       l->data[j].data[i].r,
	       l->data[j].data[i].image );
      if ( 0 ) {
	switch( l->data[j].data[i].type ) {
	default :
	  fprintf( f, " DEFAULT ?" ); break; 
	case _UNKNOWN_TYPE_POINT_ :
	  fprintf( f, " _UNKNOWN_TYPE_POINT_" ); break; 
	case _READ_ :
	  fprintf( f, " _READ_" ); break; 
	case _ADDED_ :
	  fprintf( f, " _ADDED_" ); break; 
	} 
      }
      if ( 1 ) {
	fprintf( f, ", center = %3d", l->data[j].data[i].intensity.center );
		 fprintf( f, ", %d < %d %f < %d", 
		 l->data[j].data[i].intensity.min,
		 l->data[j].data[i].intensity.median,
		 l->data[j].data[i].intensity.mean,
		 l->data[j].data[i].intensity.max );
      }
      fprintf( f, "\n" );
    }
  }
}





int readTrackedCircleList( circleListList *readlist, char *filename )
{
  char *proc = "readTrackedCircleList";
  FILE *f;
  int track;
  int i;
  circleList *list;
  circleType c;

  if ( 0 )
    initCircleListList( readlist );

  initCircleType( &c );
  c.type = _READ_;

  f = fopen( filename, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
    return( -1 );
  }

  while ( fscanf( f, "track %d", &track ) == 1 ) {
    
    list = getCircleListFromCircleListList( readlist );
    if ( list == (circleList *)NULL ) {
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to get list for track #%d\n", proc, track );
      return( -1 );
    }

    c.index = track;

    i = 0;
    while ( fscanf( f, "%d\t%d\t%d %d", &(c.x), &(c.y), &(c.r), &(c.image) ) == 4 ) {
      if ( addCircleToCircleList( list, &c ) != 1 ) {
	fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add point #%d to track #%d\n", proc, i, track );
      return( -1 );
      }
      i ++;
    }
    
    if ( 0 && _verbose_ )
      fprintf( stderr, "%s: read %d points for track #%d\n", proc, i, track );
 
  }
  
  if ( _verbose_ )
    fprintf( stderr, "%s: read %d track\n", proc, track+1 );

  fclose( f );
  return( 1 );
}





int fillTrackedCircleList( circleListList *filledlist, circleListList *readlist )
{
  char *proc = "fillTrackedCircleList";
  int track, i, j, k, l;
  circleList *list;
  circleType c;
  double a, b;
  
  initCircleType( &c );
  c.type = _ADDED_;

  for ( track=0; track<readlist->n_data; track++ ) {

    if ( _debug_ ) {
      fprintf( stderr, "%s: process track #%d\n", proc, track );
    }

    list = getCircleListFromCircleListList( filledlist );

    if ( list ==  (circleList *)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to get list for track #%d\n", proc, track );
      return( -1 );
    }
    
    for ( j=0, i=0; i<readlist->data[track].n_data; i++, j++ ) {
      
      /* add read circle
       */
      if ( addCircleToCircleList( list, &(readlist->data[track].data[i]) ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to add read point #%d to track #%d\n", proc, i, track );
	return( -1 );
      }
 
      /* last point, nothing to do
       */
      if ( i == readlist->data[track].n_data-1 ) continue;

      /* should we fill a gap ?
	 if ( l > 1 ) yes
	 interpolate with (l-k)/l * data[i] and k/l * data[i+1]
       */
      l = readlist->data[track].data[i+1].image - readlist->data[track].data[i].image;
      for ( k=1; k < l; k++, j++ ) {
	a = (double)(l-k)/(double)(l);
	b = (double)(k)/(double)(l);
	c.x = ( a * readlist->data[track].data[i].x +
		a * readlist->data[track].data[i+1].x + 0.5 );
	c.y = ( a * readlist->data[track].data[i].y +
		a * readlist->data[track].data[i+1].y + 0.5 );
	c.r = ( a * readlist->data[track].data[i].r +
		a * readlist->data[track].data[i+1].r + 0.5 );
	c.image = readlist->data[track].data[i].image + k;
	if ( addCircleToCircleList( list, &c ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to add interpolated point #%d to track #%d\n", proc, j, track );
	return( -1 );
	}
      }
      
    }
  
  }

  return( 1 );
}





typedef struct {
  int dx;
  int dy;
} displacementType;

typedef struct {
  displacementType *data;
  int n_data;
} neighborhoodType;

void initNeighborhoodType( neighborhoodType *n )
{
  n->data = (displacementType *)NULL;
  n->n_data = 0;
}

void freeNeighborhoodType( neighborhoodType *n )
{
  if ( n->data != (displacementType *)NULL )
    free( n->data );
  initNeighborhoodType( n );
}





int intensityTrackedCircleList( circleListList *filledlist, char *format )
{
  char *proc = "intensityTrackedCircleList";
  int i, j, k, l;

  int rmax, dim;
  int *distance = (int*)NULL;

  circleType *c;

  neighborhoodType *neighbor = (neighborhoodType *)NULL;

  int index, indexmin, indexmax;
  char imagename[256];
  vt_image *image;
  


  /* maximum radius
   */
  rmax = filledlist->data[0].data[0].r;
  for ( j=0; j<filledlist->n_data; j++ )
  for ( i=0; i<filledlist->data[j].n_data; i++ )
    if ( rmax < filledlist->data[j].data[i].r )
      rmax = filledlist->data[j].data[i].r;

  

  /* indexes
   */
  indexmin = indexmax = filledlist->data[0].data[0].image;
  for ( j=0; j<filledlist->n_data; j++ )
  for ( i=0; i<filledlist->data[j].n_data; i++ ) {
    if ( indexmax < filledlist->data[j].data[i].image )
      indexmax = filledlist->data[j].data[i].image;
    if ( indexmin > filledlist->data[j].data[i].image )
      indexmin = filledlist->data[j].data[i].image;
  }

  if ( 0 ) {
    fprintf( stderr, "%s: image indexes from %d to %d\n", proc, indexmin, indexmax );
  }



  /* distance map
   */
  dim = 2*rmax+1;
  distance = (int*)malloc( dim*dim * sizeof(int) );
  if ( distance == (int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  for ( j=-rmax; j<=rmax; j++ )
  for ( i=-rmax; i<=rmax; i++ )
    distance[ (j+rmax)*dim + (i+rmax) ] = i*i+j*j;



  /* neighborhood building
   */
  neighbor = (neighborhoodType *)malloc( (rmax+1)*sizeof(neighborhoodType) );
  if ( neighbor == (neighborhoodType *)NULL ) {
    free( distance );
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }
  
  for ( k=0; k<=rmax; k++ )
    initNeighborhoodType( &(neighbor[k]) );
  for ( j=-rmax; j<=rmax; j++ )
  for ( i=-rmax; i<=rmax; i++ ) {
    if ( distance[(j+rmax)*dim + (i+rmax)] > rmax*rmax ) continue;
    for ( k = 1; k <= rmax; k ++ )
      if ( distance[(j+rmax)*dim + (i+rmax)] <= k*k ) {
	neighbor[k].n_data ++;
      }
  }

  if ( 0 ) {
    for ( k=0; k<=rmax; k++ )
      fprintf( stderr, "%s: neighborhood[%d] has %d points\n", proc, k, neighbor[k].n_data );
  }

  for ( k=1; k<=rmax; k++ ) {
    /* allocation
     */
    neighbor[k].data = (displacementType*)malloc( neighbor[k].n_data * sizeof(displacementType) );
    if ( neighbor[k].data == (displacementType*)NULL ) {
      for ( l=1; l<=rmax; l++ ) {
	if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
	neighbor[l].data = (displacementType*)NULL;
      }
      free( neighbor );
      free( distance );
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error for neighborhood #%d\n", proc, k );
      return( -1 );
    }
    /* displacements
     */
    for ( l=0, j=-rmax; j<=rmax; j++ )
      for ( i=-rmax; i<=rmax; i++ ) {
      if ( distance[(j+rmax)*dim + (i+rmax)] <= k*k ) {
	neighbor[k].data[ l ].dx = i;
	neighbor[k].data[ l ].dy = j;
	l++;
      }
    }
  }

  if ( 0 ) {
    for ( k=0; k<=rmax; k++ ) {
      fprintf( stderr, "neighbor[%d] =", k );
      for ( l=0; l<neighbor[k].n_data; l++ )
	fprintf( stderr, " (%d,%d)", neighbor[k].data[ l ].dx, neighbor[k].data[ l ].dy );
      fprintf( stderr, "\n" );
    }
  }

  free( distance );

  

  /* intensity
   */

  if ( format == (char*)NULL ) {
    for ( l=1; l<=rmax; l++ ) {
      if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
	neighbor[l].data = (displacementType*)NULL;
    }
    free( neighbor );
    if ( _verbose_ )
      fprintf( stderr, "%s: no image format\n", proc );
    return( -1 );
  }



  for ( j=0; j<filledlist->n_data; j++ )
  for ( i=0; i<filledlist->data[j].n_data; i++ ) {
    c = &(filledlist->data[j].data[i]);
    c->intensity.data = (int*)malloc( neighbor[ c->r ].n_data * sizeof(int) );
    if ( c->intensity.data == (int*)NULL ) {
      for ( l=1; l<=rmax; l++ ) {
	if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
	neighbor[l].data = (displacementType*)NULL;
      }
      free( neighbor );
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error for intensities of pit #%d of track #%d\n", proc, j, i );
      return( -1 );
    }
    c->intensity.n_data = neighbor[c->r].n_data;
  }


  for ( index=indexmin; index<=indexmax; index++ ) {

    sprintf( imagename, format, index );
    image = _VT_Inrimage( imagename );
    if ( image == (vt_image*)NULL ) {
      for ( l=1; l<=rmax; l++ ) {
	if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
	neighbor[l].data = (displacementType*)NULL;
      }
      free( neighbor );
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error for image #%d\n", proc, index );
      return( -1 );
    }

    switch( image->type ) {
    default :
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      for ( l=1; l<=rmax; l++ ) {
	if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
	neighbor[l].data = (displacementType*)NULL;
      }
      free( neighbor );
      if ( _verbose_ )
	fprintf( stderr, "%s: such type not handled for image #%d\n", proc, index );
      return( -1 );
    case UCHAR :
      {
	u8 ***theArray = (u8***)image->array;
	for ( j=0; j<filledlist->n_data; j++ )
	for ( i=0; i<filledlist->data[j].n_data; i++ ) {
	  if ( filledlist->data[j].data[i].image != index ) continue;
	  c = &(filledlist->data[j].data[i]);
	  c->intensity.center = theArray[0][ c->y ][ c->x ];
	  for ( l=0; l<neighbor[ c->r ].n_data; l++ ) {
	    c->intensity.data[l] = theArray[0][ c->y + neighbor[c->r].data[l].dy ][ c->x + neighbor[c->r].data[l].dx ];
	  }
	}
      }
      break;
    case SSHORT :
      {
	s16 ***theArray = (s16***)image->array;
	for ( j=0; j<filledlist->n_data; j++ )
	for ( i=0; i<filledlist->data[j].n_data; i++ ) {
	  if ( filledlist->data[j].data[i].image != index ) continue;
	  c = &(filledlist->data[j].data[i]);
	  c->intensity.center = theArray[0][ c->y ][ c->x ];
	  for ( l=0; l<neighbor[ c->r ].n_data; l++ ) {
	    c->intensity.data[l] = theArray[0][ c->y + neighbor[c->r].data[l].dy ][ c->x + neighbor[c->r].data[l].dx ];
	  }
	}
      }
      break;    }
    
    VT_FreeImage( image );
    VT_Free( (void**)&image );
  }
  


  for ( l=1; l<=rmax; l++ ) {
    if ( neighbor[l].data != (displacementType*)NULL ) free( neighbor[l].data );
    neighbor[l].data = (displacementType*)NULL;
  }
  free( neighbor );


  for ( j=0; j<filledlist->n_data; j++ )
  for ( i=0; i<filledlist->data[j].n_data; i++ ) {
    computeValueType( &(filledlist->data[j].data[i].intensity) );
  }

  return( 1 );
}






/************************************************************
 *
 * MANAGEMENT (chainType)
 *
 ************************************************************/

void initChainType( chainType *l )
{
  l->data = (circleType**)NULL;

  l->xmin = 0;
  l->xmax = 0;
  l->ymin = 0;
  l->ymax = 0;

  l->n_data = 0;
  l->n_allocated_data = 0;
  l->rmax = 0;
} 



void freeChainType( chainType *l )
{
  if ( l->data != (circleType**)NULL ) {
    free( l->data );
  } 
  initChainType( l );
} 



static int addCircleToChain( chainType *l, circleType *c )
{
  char *proc = "addCircleToChain";
  int s =  l->n_allocated_data;
  circleType **data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (circleType**)malloc( s * sizeof(circleType*) );
    if ( data == (circleType**)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(circleType*) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  l->data[l->n_data] = c;
  l->n_data ++;

  c->n_chains ++;

  return( 1 );
} 






/************************************************************
 *
 * MANAGEMENT (chainList)
 *
 ************************************************************/


void initChainList( chainList *l )
{
  l->data = (chainType*)NULL;
  l->n_data = 0;
  l->n_selected_data = 0;
  l->n_allocated_data = 0;
  l->rmax = 0;
} 



void freeChainList( chainList *l )
{
  int i;

  if ( l->data != (chainType*)NULL ) {
    for ( i=0; i<l->n_data; i++ ) {
      freeChainType( &(l->data[i]) );
    } 
    free( l->data );
  } 
  initChainList( l );
} 



static chainType *getChainTypeFromChainList( chainList *l )
{
  char *proc = "getChainTypeFromChainList";
  int s =  l->n_allocated_data;
  chainType *data;
  int i;
  
  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (chainType*)malloc( s * sizeof(chainType) );
    if ( data == (chainType*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( (chainType *)NULL );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(chainType) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
    for ( i=l->n_data; i<l->n_allocated_data; i++ )
      initChainType( &(l->data[i]) );
  }
  l->n_data ++;
  l->n_selected_data ++;
  return( &(l->data[l->n_data-1]) );
} 



static void computeChainListBoundingBoxes( chainList *list )
{
  int l, i;

  for ( l=0; l<list->n_data; l++ ) {
    list->data[l].xmin = list->data[l].data[0]->x - list->data[l].data[0]->r;
    list->data[l].xmax = list->data[l].data[0]->x + list->data[l].data[0]->r;
    list->data[l].ymin = list->data[l].data[0]->y - list->data[l].data[0]->r;
    list->data[l].ymax = list->data[l].data[0]->y + list->data[l].data[0]->r;
    
    for ( i=1; i<list->data[l].n_data; i++ ) {
      if ( list->data[l].xmin > list->data[l].data[i]->x - list->data[l].data[i]->r )
	list->data[l].xmin = list->data[l].data[i]->x - list->data[l].data[i]->r;
      if ( list->data[l].xmax < list->data[l].data[i]->x + list->data[l].data[i]->r )
	list->data[l].xmax = list->data[l].data[i]->x + list->data[l].data[i]->r;
      if ( list->data[l].ymin > list->data[l].data[i]->y - list->data[l].data[i]->r )
	list->data[l].ymin = list->data[l].data[i]->y - list->data[l].data[i]->r;
      if ( list->data[l].ymax < list->data[l].data[i]->y + list->data[l].data[i]->r )
	list->data[l].ymax = list->data[l].data[i]->y + list->data[l].data[i]->r;
    }
  }
}









/************************************************************
 *
 * GENERIC OUTPUT
 *
 ************************************************************/



static void initPlotType( plotType *p )
{
  p->n_data = 0;
  p->data = (int*)NULL;
  p->dim = 0;
  p->x = (int*)NULL;
  p->y = (int*)NULL;
} 



static void freePlotType( plotType *p )
{
  if ( p->data != (int*)NULL )
    free( p->data );
  if ( p->x != (int*)NULL )
    free( p->x );
  if ( p->y != (int*)NULL )
    free( p->y );
  initPlotType( p );
} 



static int allocPlotType( plotType *p, int n )
{
  char *proc = "allocPlotType";
  p->data = (int*)malloc( n * sizeof( int ) );
  if ( p->data == (int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  } 
  p->n_data = n;
  return( 1 );
} 



static int plotFromData( plotType *p  )
{
  char *proc = "plotFromData";
  int i, j;
  
  p->x = (int*)malloc( p->n_data * sizeof( int ) );
  if ( p->x == (int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  } 
  p->y = (int*)malloc( p->n_data * sizeof( int ) );
  if ( p->y == (int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  } 
  
  for ( i=0, j=0; i<p->n_data; i++ ) {
    if ( p->data[i] == 0 ) continue;
    p->x[j] = i; 
    p->y[j] = p->data[i]; 
    j++;
  } 
  p->dim = j;
  return( 1 );
} 



void printPlotList( FILE *f, plotType *p )
{
  int i, j;
  for ( i=0; i<p->dim; i++ )  {
    for (j=0; j<p->y[i]; j++ )
      fprintf( f, "%d\n", p->x[i] );
  }
}



void printPlotHist( FILE *f, plotType *p )
{
  int i;
  for ( i=0; i<p->dim; i++ )  {
    fprintf( f, "%d %d\n", p->x[i], p->y[i] );
  }
}



static int writePlotTypeScilab( FILE *f, int fd, plotType *p, char *s )
{
  char *proc="writePlotTypeScilab";
  int i, xmax, ymax;

  xmax = p->x[0];
  ymax = p->y[0];
  for ( i=0; i<p->dim; i++ )  {
    if ( xmax < p->x[i] ) xmax = p->x[i];
    if ( ymax < p->y[i] ) ymax = p->y[i];
  }
  
  if ( write( fd, p->x, p->dim * sizeof(s32) ) == -1 ) {
    fprintf( stderr, "%s: error when writing x\n", proc );
  }
  if ( write( fd, p->y, p->dim * sizeof(s32) ) == -1 ) {
    fprintf( stderr, "%s: error when writing y\n", proc );
  }

  fprintf( f, "\n" );
  fprintf( f, "//\n" );
  fprintf( f, "// read data\n" );
  fprintf( f, "//\n" );
  fprintf( f, "\n" );
  
  if ( s == (char*)NULL )
    fprintf( f, "INDEX=mget( %d, 'i', f);\n", p->dim );
  else 
    fprintf( f, "INDEX_%s=mget( %d, 'i', f);\n", s, p->dim );

  if ( s == (char*)NULL )
    fprintf( f, "DATA=mget( %d, 'i', f);\n", p->dim );
  else 
    fprintf( f, "DATA_%s=mget( %d, 'i', f);\n", s, p->dim );

  fprintf( f, "\n" );
  fprintf( f, "//\n" );
  fprintf( f, "// figure\n" );
  fprintf( f, "//\n" );
  fprintf( f, "\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "set(gca(),\"auto_clear\",\"off\");\n" );
  if ( s == (char*)NULL )
    fprintf( f, "plot( INDEX, DATA, \"k-\", \"thickness\", 2 );\n" );
  else 
    fprintf( f, "plot( INDEX_%s, DATA_%s, \"k-\", \"thickness\", 2 );\n", s, s );
  
  fprintf( f, "\n" );

  fprintf( f, "// a=get(\"current_axes\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "// removing the trailing ';' allows to see all properties\n" );
  fprintf( f, "// a.data_bounds = [0,0;%d,%d];\n", xmax, ymax  );
  fprintf( f, "// a.font_size=3;\n" );
  fprintf( f, "// a.title.text = \"Title\";\n" );
  fprintf( f, "// a.x_label.text = \"X Label\";\n" );
  fprintf( f, "// a.y_label.text = \"Y Label\";\n" );
  fprintf( f, "// ou \n" );
  fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
  if ( s != (char*)NULL ) {
    fprintf( f, "xtitle( \"%s\" );\n", s );
    fprintf( f, "a.title.font_size = 2;\n" );
  }
  else  {
    fprintf( f, "// a.title.font_size = 2;\n" );
  } 
  fprintf( f, "// a.x_label.font_size = 2;\n" );
  fprintf( f, "// a.y_label.font_size = 2;\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "// xs2jpg(gcf(),'FIG%s.jpg');\n", s  );
  fprintf( f, "// xs2png(gcf(),'FIG%s.png');\n", s  );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  return( 1 );
} 


/************************************************************
 *
 * RESULTS
 *
 ************************************************************/





void initStatisticType( statisticType *s )
{
  initPlotType( &(s->ncircle) );
  initPlotType( &(s->startingIndex) );
  initPlotType( &(s->endingIndex) );
  initPlotType( &(s->length) );
} 



void freeStatisticType( statisticType *s )
{
  freePlotType( &(s->ncircle) );
  freePlotType( &(s->startingIndex) );
  freePlotType( &(s->endingIndex) );
  freePlotType( &(s->length) );
  initStatisticType( s );
} 


static int writeStatisticTypeScilab( FILE *f, int fd, statisticType *p, char *s )
{
  char *proc = "writeStatisticTypeScilab";
  char str[256];

  if ( s == (char*)NULL )
    sprintf( str, "nc" );
  else 
    sprintf( str, "nc_%s", s );
  fprintf( f, "//\n" );
  fprintf( f, "// detected circle per image\n" );
  fprintf( f, "\n" );
  if ( writePlotTypeScilab( f, fd, &(p->ncircle), str ) != 1 ) {
    fprintf( stderr, "%s: error when writing number of samples\n", proc );
    return( -1 );
  }

  if ( s == (char*)NULL )
    sprintf( str, "stc" );
  else 
    sprintf( str, "stc_%s", s );
  fprintf( f, "//\n" );
  fprintf( f, "// start time (ie image number) of chains\n" );
  fprintf( f, "\n" );
  if ( writePlotTypeScilab( f, fd, &(p->startingIndex), str ) != 1 ) {
    fprintf( stderr, "%s: error when writing starting index\n", proc );
    return( -1 );
  }

  if ( s == (char*)NULL )
    sprintf( str, "etc" );
  else 
    sprintf( str, "etc_%s", s );
  fprintf( f, "//\n" );
  fprintf( f, "// end time (ie image number) of chains\n" );
  fprintf( f, "\n" );
  if ( writePlotTypeScilab( f, fd, &(p->endingIndex), str ) != 1 ) {
    fprintf( stderr, "%s: error when writing ending index\n", proc );
    return( -1 );
  }

  if ( s == (char*)NULL )
    sprintf( str, "cl" );
  else 
    sprintf( str, "cl_%s", s );
  fprintf( f, "//\n" );
  fprintf( f, "// chain length\n" );
  fprintf( f, "\n" );
  if ( writePlotTypeScilab( f, fd, &(p->length), str ) != 1 ) {
    fprintf( stderr, "%s: error when writing chain length\n", proc );
    return( -1 );
  }

  return( 1 );
} 



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}








static void printResultXxxlab( statisticType *stats1, 
			       char *desc1,
			       statisticType *stats2, 
			       char *desc2,
			       statisticType *stats3, 
			       char *desc3,
			       char *name, 
			       enumHistogramFile xxxlab )
{
  char *proc = "printResultXxxlab";
  char *defaultname = "histo";
  char *template;
  char *filename = (char*)NULL;
  FILE *f;
  int fd;


  template = ( name != (char*)NULL ) ? name : defaultname;
  filename = (char*)malloc( strlen( template ) + 5 );
  if ( filename == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate file name\n", proc );
    return;
  }

  /* open files
   */
  sprintf( filename, "%s.raw", template );

  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    free( filename );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return;
  }

 
  switch( xxxlab ) {
  default :
    free( filename );
    close( fd );
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not known\n", proc );
    return;
  case _MATLAB_ :
    sprintf( filename, "%s.m", template );
    break;
  case _SCILAB_ :
    sprintf( filename, "%s.sce", template );
    break;
  }

  f = fopen( filename, "w" );

  if ( f == (FILE*)NULL ) {
    free( filename );
    close( fd );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return;
  }

  
  /* write data and script
   */
  switch( xxxlab ) {
  default :
    free( filename );
    close( fd );
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not known\n", proc );
    return;

  case _MATLAB_ :

    fprintf( f, "\n" );
    fprintf( f, "f=fopen('%s.raw','r');\n", _BaseName( template ) );
    fprintf( f, "\n" );

    
      ;

    
    fprintf( f, "\n" );
    fprintf( f, "fclose(f);\n" );
    fprintf( f, "\n" );
    break;

  case _SCILAB_ :
    
    fprintf( f, "\n" );
    fprintf( f, "f=mopen('%s.raw','r');\n", _BaseName( template ) );
    fprintf( f, "\n" );

    if ( stats1 != (statisticType *)NULL ) {
      if ( writeStatisticTypeScilab( f, fd, stats1, desc1 ) != 1 ) {
	free( filename );
	close( fd );
	fclose( f );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when writing statistics #1 for scilab\n", proc );
	return;
      }
    }

    if ( stats2 != (statisticType *)NULL ) {
      if ( writeStatisticTypeScilab( f, fd, stats2, desc2 ) != 1 ) {
	free( filename );
	close( fd );
	fclose( f );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when writing statistics #2 for scilab\n", proc );
	return;
      }
    }

    if ( stats3 != (statisticType *)NULL ) {
      if ( writeStatisticTypeScilab( f, fd, stats3, desc3 ) != 1 ) {
	free( filename );
	close( fd );
	fclose( f );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when writing statistics #3 for scilab\n", proc );
	return;
      }
    }



    fprintf( f, "\n" );
    fprintf( f, "mclose(f);\n" );
    fprintf( f, "\n" );

    break;
  }


  /* close files
   */
  close( fd );
  fclose( f );

  free( filename );
}





void printChainResultXxxlab( statisticType *stats, 
			     char *desc,
			     char *name, 
			     enumHistogramFile xxxlab )
{
  printResultXxxlab( stats, desc, 
		     (statisticType *)NULL, (char*)NULL, 
		     (statisticType *)NULL, (char*)NULL, name, xxxlab );
}

void printColocalizationResultXxxlab( statisticType *stats1, 
				      char *desc1,
				      statisticType *stats2, 
				      char *desc2,
				      statisticType *stats3, 
				      char *desc3,
				      char *name, 
				      enumHistogramFile xxxlab )
{
  printResultXxxlab( stats1, desc1, stats2, desc2, stats3, desc3, name, xxxlab );
}



/************************************************************
 *
 * CHAINS
 *
 ************************************************************/


void rmaxCircleListList( circleListList *list )
{
  int l, i;
  int rmax;
  
  for ( l=0; l<list->n_data; l++ ) {
    for ( rmax=0, i=0; i<list->data[l].n_selected_data; i++ )
      if ( rmax < list->data[l].data[i].r )
	rmax = list->data[l].data[i].r;
    list->data[l].rmax = rmax;
  } 
  
  for ( rmax=0, l=0; l<list->n_data; l++ ) {
    if ( rmax < list->data[l].rmax )
	rmax = list->data[l].rmax;
    list->rmax = rmax;
  } 
} 



static int _neighbors_to_be_allocated_ = 10;

static int addForwardNeighbor( circleType *c, circleType *n )
{
  char *proc = "addForwardNeighbor";
  int s =  c->n_next_allocated;
  circleType **next;	    
  
  if ( c->n_next == c->n_next_allocated ) {
    s += _neighbors_to_be_allocated_;
    next = (circleType**)malloc( s * sizeof(circleType*) );
    if ( next == (circleType**)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( c->n_next_allocated > 0 ) {
      (void)memcpy( next, c->next, c->n_next_allocated*sizeof(circleType*) );
      free( c->next );
    }
    c->n_next_allocated = s;
    c->next = next;
  }

  c->next[c->n_next] = n;
  c->n_next ++;
  return( 1 );
} 







static int searchForwardNeighbors( circleListList *list, int depth, int margin )
{
  char *proc = "searchForwardNeighbors";
  int l, i, j, n;
  int found;
  int dx, dy, r;
  int rmax;

  int totaladdedpoint = 0;
  int addedpoint;
  circleType *current;
  int k;
  circleType c;

  rmax = list->rmax + list->rmax + margin;

  for ( l=0; l<list->n_data; l++ ) {
    for ( i=0; i<list->data[l].n_data; i++ ) {
      list->data[l].data[i].n_next = 0;
    }
  }

  /* process each list (except the last one)
   */
  for ( l=0; l<list->n_data-1; l++ ) {
    if ( 0 ) fprintf( stderr, "%s: processing list #%3d\n", proc, l );
    
    /* point loop 
       ce n'est pas possible de paralleliser car il y a des points ajoutes 
     */
    for ( addedpoint=0, i=0; i<list->data[l].n_selected_data; i++ ) {

      for ( found=0, n=l+1; found==0 && n<=l+depth && n<list->n_data; n++ ) {
	for ( j=0; j<list->data[n].n_selected_data; j++ ) {

	  /* a neighbor is found if the center distance
	     is below the sum of the 2 radii + the margin
	  */
	  dx = list->data[l].data[i].x - list->data[n].data[j].x; 
	  if ( dx < -rmax || rmax < dx ) continue;
	  dy = list->data[l].data[i].y - list->data[n].data[j].y;
	  if ( dy < -rmax || rmax < dy ) continue;
	  r = list->data[l].data[i].r + list->data[n].data[j].r + margin;
	  if ( dx*dx+dy*dy <= r*r ) {
	    current = &(list->data[l].data[i]);
	    if ( n > l+1 ) {
	      for ( k = 1; k <= n-l-1 ; k ++ ) {
		addedpoint ++;
		initCircleType( &c );
		c.x = (int)( (float)list->data[n].data[j].x * (float)k/(float)(n-l) 
			   + (float)list->data[l].data[i].x * (float)(n-l-k)/(float)(n-l) + 0.5); 
		c.y = (int)( (float)list->data[n].data[j].y * (float)k/(float)(n-l) 
			   + (float)list->data[l].data[i].y * (float)(n-l-k)/(float)(n-l) + 0.5); 
		c.r = (int)( (float)list->data[n].data[j].r * (float)k/(float)(n-l) 
			   + (float)list->data[l].data[i].r * (float)(n-l-k)/(float)(n-l) + 0.5); 
		c.index = l+k;
		c.type = _ADDED_;
		/* add c to list data[l+k], it will be the last point
		   if the (data[l+k].n_data-1)th point
		*/
		if ( addCircleToCircleList( &(list->data[l+k]), &c ) != 1 ) {
		  if ( _verbose_ )
		    fprintf( stderr, "%s: unable to add built circle (=%d %d %d) to list\n", proc, c.x, c.y, c.r );
		  return( -1 );
		}
		if ( addForwardNeighbor( current, &(list->data[l+k].data[ list->data[l+k].n_data-1 ]) ) != 1 ) {
		  if ( _verbose_ )
		    fprintf( stderr, "%s: unable to add built neighbor\n", proc );
		  return( -1 );
		}
		current = &(list->data[l+k].data[ list->data[l+k].n_data-1 ]);
	      }
	    }
	    
	    if ( addForwardNeighbor( current, &(list->data[n].data[j]) ) != 1 ) {
	      if ( _verbose_ )
		fprintf( stderr, "%s: unable to add neighbor\n", proc );
	      return( -1 );
	    }
	    found ++;
	  }
	}
      }
    } 
    /* end of point loop */
    if ( _verbose_ > 1 && addedpoint > 0 )
      fprintf( stderr, "%s: list #%3d: add %5d points\n", proc, l, addedpoint );

    totaladdedpoint += addedpoint;
    
  }
  
  if ( _verbose_ && totaladdedpoint > 0 )
    fprintf( stderr, "%s: build %d points to close jumps\n", proc, totaladdedpoint );
  return( 1 );
}




static void statsForwardNeighbors( FILE *f, circleListList *list )
{
  char *proc = "statsForwardNeighbors";
  int l, i;
  int max = 0, t;
  int *count;

  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( max < list->data[l].data[i].n_next )
	max = list->data[l].data[i].n_next;

  count = (int*)malloc( (max+1)*sizeof(int) );
  if ( count == (int*)NULL ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    count[i] = 0;
  
  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      count[ list->data[l].data[i].n_next ] ++;
  
  for ( t=0, i=0; i<=max; i++ )
    t += count[i];

  fprintf( f, "\n" );
  fprintf( f, "--- Forward neighbors\n" );
  for ( i=0; i<=max; i++ )
    fprintf( f, "    %2d neighbors: %8d %5.2f%%\n", i, count[i], 100.0*(float)count[i]/(float)t );
  fprintf( f, "\n" );
  
  free( count );
}





static int addBackwardNeighbor( circleType *c, circleType *n )
{
  char *proc = "addBackwardNeighbor";
  int s =  c->n_prev_allocated;
  circleType **prev;	    
  
  if ( c->n_prev == c->n_prev_allocated ) {
    s += _neighbors_to_be_allocated_;
    prev = (circleType**)malloc( s * sizeof(circleType*) );
    if ( prev == (circleType**)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( c->n_prev_allocated > 0 ) {
      (void)memcpy( prev, c->prev, c->n_prev_allocated*sizeof(circleType*) );
      free( c->prev );
    }
    c->n_prev_allocated = s;
    c->prev = prev;
  }

  c->prev[c->n_prev] = n;
  c->n_prev ++;
  return( 1 );
} 



static int searchBackwardNeighbors( circleListList *list, int depth, int margin )
{
  char *proc = "searchBackwardNeighbors";
  int l, i, j, p;
  int found;
  int dx, dy, r;
  int rmax;

  int totaladdedpoint = 0;
  int addedpoint;
  circleType *current;
  int k;
  circleType c;

  rmax = list->rmax + list->rmax + margin;

  for ( l=0; l<list->n_data; l++ ) {
    for ( i=0; i<list->data[l].n_data; i++ ) {
      list->data[l].data[i].n_prev = 0;
    }
  }

  /* process each list (except the last one)
   */
  for ( l=1; l<list->n_data; l++ ) {
    if ( 0 ) fprintf( stderr, "%s: processing list #%3d\n", proc, l );
    
    /* point loop */
    for ( addedpoint=0, i=0; i<list->data[l].n_selected_data; i++ ) {

      for ( found=0, p=l-1; found==0 && p>=l-depth && p>=0; p-- ) {
	for ( j=0; j<list->data[p].n_selected_data; j++ ) {
	  /* a neighbor is found if the center distance
	     is below the sum of the 2 radii + the margin
	  */
	  dx = list->data[l].data[i].x - list->data[p].data[j].x; 
	  if ( dx < -rmax || rmax < dx ) continue;
	  dy = list->data[l].data[i].y - list->data[p].data[j].y;
	  if ( dy < -rmax || rmax < dy ) continue;
	  r = list->data[l].data[i].r + list->data[p].data[j].r + margin;
	  if ( dx*dx+dy*dy <= r*r ) {
	    current = &(list->data[l].data[i]);
	    if ( p < l-1 ) {
	      /* add points
	       */
	      for ( k = 1; k <= l-p-1 ; k ++ ) {
		addedpoint ++;
		initCircleType( &c );
		c.x = (int)( (float)list->data[p].data[j].x * (float)k/(float)(l-p) 
			   + (float)list->data[l].data[i].x * (float)(l-p-k)/(float)(l-p) + 0.5); 
		c.y = (int)( (float)list->data[p].data[j].y * (float)k/(float)(l-p) 
			     + (float)list->data[l].data[i].y * (float)(l-p-k)/(float)(l-p) + 0.5); 
		c.r = (int)( (float)list->data[p].data[j].r * (float)k/(float)(l-p) 
		           + (float)list->data[l].data[i].r * (float)(l-p-k)/(float)(l-p) + 0.5);
		c.index = l-k;
		c.type = _ADDED_;
		/* add c to list data[l-k], it will be the last point
		   if the (data[l-k].n_data-1)th point
		*/
		if ( addCircleToCircleList( &(list->data[l-k]), &c ) != 1 ) {
		  if ( _verbose_ )
		    fprintf( stderr, "%s: unable to add built circle (=%d %d %d) to list\n", proc, c.x, c.y, c.r );
		  return( -1 );
		}
		if ( addBackwardNeighbor( current, &(list->data[l-k].data[ list->data[l-k].n_data-1 ]) ) != 1 ) {
		  if ( _verbose_ )
		    fprintf( stderr, "%s: unable to add built neighbor\n", proc );
		  return( -1 );
		}
		current = &(list->data[l-k].data[ list->data[l-k].n_data-1 ]);
	      }
	    }
	    
	    if ( addBackwardNeighbor( current, &(list->data[p].data[j]) ) != 1 ) {
	      if ( _verbose_ )
		fprintf( stderr, "%s: unable to add neighbor\n", proc );
	      return( -1 );
	    }
	    found ++;
	  }
	}
      }
    }
    /* end of point loop */
    if ( _verbose_ > 1 && addedpoint > 0 )
      fprintf( stderr, "%s: list #%3d: add %5d points\n", proc, l, addedpoint );

    totaladdedpoint += addedpoint;

  }
  
  if ( _verbose_ && totaladdedpoint > 0 )
    fprintf( stderr, "%s: build %d points to close jumps\n", proc, totaladdedpoint );
  return( 1 );
}



static void statsBackwardNeighbors( FILE *f, circleListList *list )
{
  char *proc = "statsBackwardNeighbors";
  int l, i;
  int max = 0, t;
  int *count;

  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( max < list->data[l].data[i].n_prev )
	max = list->data[l].data[i].n_prev;

  count = (int*)malloc( (max+1)*sizeof(int) );
  if ( count == (int*)NULL ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    count[i] = 0;
  
  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      count[ list->data[l].data[i].n_prev ] ++;
  
  for ( t=0, i=0; i<=max; i++ )
    t += count[i];

  fprintf( f, "\n" );
  fprintf( f, "--- Backward neighbors\n" );
  for ( i=0; i<=max; i++ )
    fprintf( f, "    %2d neighbors: %8d %5.2f%%\n", i, count[i], 100.0*(float)count[i]/(float)t );
  fprintf( f, "\n" );
  
  free( count );
}





static void statsNeighbors( FILE *f, circleListList *list )
{
  char *proc = "statsNeighbors";
  int l, i, j;
  int maxforward = 0;
  int maxbackward = 0;
  int t;
  int *count;

  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( maxforward < list->data[l].data[i].n_next )
	maxforward = list->data[l].data[i].n_next;

  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( maxbackward < list->data[l].data[i].n_prev )
	maxbackward = list->data[l].data[i].n_prev;
  

  count = (int*)malloc( (maxforward+1)*(maxbackward+1)*sizeof(int) );
  if ( count == (int*)NULL ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<(maxforward+1)*(maxbackward+1); i++ )
    count[i] = 0;
  
  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ ) 
      count[ list->data[l].data[i].n_next + list->data[l].data[i].n_prev * (maxforward+1) ] ++;
  
  for ( t=0, i=0; i<(maxforward+1)*(maxbackward+1); i++ )
    t += count[i];

  

  fprintf( f, "\n" );
  fprintf( f, "--- Backward \\ Forward neighbors\n" );

  fprintf( f, "        ");
  for ( i=0; i<=maxforward; i++ ) 
    fprintf( f, " %7d", i );
  fprintf( f, "\n" );
  
  for ( j=0; j<=maxbackward; j++ ) {
    fprintf( f, " %7d", j );
    for ( i=0; i<=maxforward; i++ ) {
      fprintf( f, " %7d", count[ i + j * (maxforward+1) ]  );
    }
    fprintf( f, "\n" );
  }

  fprintf( f, "\n" );
  fprintf( f, "--- Backward \\ Forward neighbors\n" );

  fprintf( f, "        ");
  for ( i=0; i<=maxforward; i++ ) 
    fprintf( f, " %6d ", i );
  fprintf( f, "\n" );
  
  for ( j=0; j<=maxbackward; j++ ) {
    fprintf( f, " %7d", j );
    for ( i=0; i<=maxforward; i++ ) {
      fprintf( f, "  %5.2f%%", 100.0*count[ i + j * (maxforward+1) ]/(float)t  );
    }
    fprintf( f, "\n" );
  }

  

  

  fprintf( f, "\n" );
  
  free( count );
}



static void statsCircleRadius( FILE *f, circleListList *list )
{
  char *proc = "statsCircleRadius";
  int l, i;
  int max = 0, t;
  int *count;

  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( max < list->data[l].data[i].r )
	max = list->data[l].data[i].r;

  count = (int*)malloc( (max+1)*sizeof(int) );
  if ( count == (int*)NULL ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    count[i] = 0;
  
  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      count[ list->data[l].data[i].r ] ++;
  
  for ( t=0, i=0; i<=max; i++ )
    t += count[i];

  fprintf( f, "\n" );
  fprintf( f, "--- Radii\n" );
  for ( i=0; i<=max; i++ ) {
    if ( count[i] == 0 ) continue;
    fprintf( f, "    radius = %2d: %8d %5.2f%%\n", i, count[i], 100.0*(float)count[i]/(float)t );
  }
  fprintf( f, "\n" );
  
  free( count );
}



static void statsNumberCircles( plotType *p, circleListList *list )
{
  char *proc = "statsNumberCircles";
  int i;

  if ( allocPlotType( p, list->n_data ) != 1 ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<list->n_data; i++ )
    p->data[i] = list->data[i].n_data;
  
  if ( plotFromData( p ) != 1 ) {
    fprintf( stderr, "%s: error when computing plot data\n", proc );
    return;
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "--- Circle number \n" );
    for ( i=0; i<list->n_data; i++ )
      if ( p->data[i] > 0 )
	fprintf( stderr, "    circle number of image %3d: %8d \n", i, p->data[i] );
    fprintf( stderr, "\n" );
  }
}



#ifdef _UNUSED_
static void statsNumberSelectedCircles( plotType *p, circleListList *list )
{
  char *proc = "statsNumberSelectedCircles";
  int i;

  if ( allocPlotType( p, list->n_data ) != 1 ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<list->n_data; i++ )
    p->data[i] = list->data[i].n_selected_data;
  
  if ( plotFromData( p ) != 1 ) {
    fprintf( stderr, "%s: error when computing plot data\n", proc );
    return;
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "--- Selected ircle number \n" );
    for ( i=0; i<list->n_data; i++ )
      if ( p->data[i] > 0 )
	fprintf( stderr, "    selected circle number of image %3d: %8d \n", i, p->data[i] );
    fprintf( stderr, "\n" );
  }
}
#endif



static void statsMembershipChains( FILE *f, circleListList *list, chainList *chainlist )
{
  char *proc = "statsMembershipChains";
  int l, i;
  int max = 0, t;
  int *count;


  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      if ( max < list->data[l].data[i].n_chains)
	max = list->data[l].data[i].n_chains;

  count = (int*)malloc( (max+1)*sizeof(int) );
  if ( count == (int*)NULL ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    count[i] = 0;
  
  for ( l=0; l<list->n_data; l++ ) 
    for ( i=0; i<list->data[l].n_selected_data; i++ )
      count[ list->data[l].data[i].n_chains ] ++;
  
  for ( t=0, i=0; i<=max; i++ )
    t += count[i];

  fprintf( f, "\n" );
  fprintf( f, "--- Membership to chains\n" );
  for ( i=0; i<=max; i++ )
    if ( count[i] > 0 ) 
      fprintf( f, "    membership to %2d chain(s): %8d %5.2f%%\n", i, count[i], 100.0*(float)count[i]/(float)t );
  fprintf( f, "\n" );
  
  free( count );
}




static void statsLengthChain( plotType *p, chainList *list )
{
  char *proc = "statsLengthChain";
  int l, i;
  int max = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    if ( max < list->data[l].n_data )
      max = list->data[l].n_data;

  if ( allocPlotType( p, max+1 ) != 1 ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    p->data[i] = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    p->data[ list->data[l].n_data ] ++;
 
  if ( plotFromData( p ) != 1 ) {
    fprintf( stderr, "%s: error when computing plot data\n", proc );
    return;
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "--- Chain length\n" );
    for ( i=0; i<=max; i++ )
      if ( p->data[i] > 0 )
	fprintf( stderr, "    chain length of %3d: %8d \n", i, p->data[i] );
    fprintf( stderr, "\n" );
  }
}



static void statsStartChain( plotType *p, chainList *list )
{
  char *proc = "statsStartChain";
  int l, i;
  int max = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    if ( max < list->data[l].data[0]->index )
      max = list->data[l].data[0]->index;

  if ( allocPlotType( p, max+1 ) != 1 ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    p->data[i] = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    p->data[ list->data[l].data[0]->index ] ++;
   
  if ( plotFromData( p ) != 1 ) {
    fprintf( stderr, "%s: error when computing plot data\n", proc );
    return;
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "--- Starting index\n" );
    for ( i=0; i<=max; i++ )
      if ( p->data[i] > 0 )
	fprintf( stderr, "    starting index at %3d: %8d \n", i, p->data[i] );
    fprintf( stderr, "\n" );
  }
}



static void statsEndChain( plotType *p, chainList *list )
{
  char *proc = "statsEndChain";
  int l, i;
  int max = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    if ( max < list->data[l].data[ list->data[l].n_data-1 ]->index )
      max = list->data[l].data[ list->data[l].n_data-1 ]->index;

 if ( allocPlotType( p, max+1 ) != 1 ) {
    fprintf( stderr, "%s: allocation error\n", proc );
    return;
  }
  for ( i=0; i<=max; i++ )
    p->data[i] = 0;
  
  for ( l=0; l<list->n_selected_data; l++ ) 
    p->data[ list->data[l].data[ list->data[l].n_data-1 ]->index ] ++;
   
  if ( plotFromData( p ) != 1 ) {
    fprintf( stderr, "%s: error when computing plot data\n", proc );
    return;
  }

  if ( 0 ) {
    fprintf( stderr, "\n" );
    fprintf( stderr, "--- Ending index\n" );
    for ( i=0; i<=max; i++ )
      if ( p->data[i] > 0 )
	fprintf( stderr, "    ending index at %3d: %8d \n", i, p->data[i] );
    fprintf( stderr, "\n" );
  }
}



void statsChainList( statisticType *s, circleListList *readlist, chainList *chainlist )
{
  statsNumberCircles( &(s->ncircle), readlist );
  statsLengthChain( &(s->length), chainlist );
  statsStartChain( &(s->startingIndex), chainlist );
  statsEndChain( &(s->endingIndex), chainlist );
}





/************************************************************
 *
 * FILTERING
 *
 ************************************************************/



void selectOnForwardNeighbors( circleListList *list, int neighbors )
{
  circleType tmp;  
  int l, i;

  if ( neighbors <= 0 ) return;

  for ( l=0; l<list->n_data; l++ ) {
    for ( i=0; i<list->data[l].n_selected_data; ) {
      if ( list->data[l].data[i].n_next >= neighbors ) {
	tmp = list->data[l].data[ list->data[l].n_selected_data-1 ];
	list->data[l].data[ list->data[l].n_selected_data-1 ] = list->data[l].data[ i ];
	list->data[l].data[ i ] = tmp;
	list->data[l].n_selected_data --;
      }
      else
	i++;
    }
  }
}



void selectOnBackwardNeighbors( circleListList *list, int neighbors )
{
  circleType tmp;  
  int l, i;

  if ( neighbors <= 0 ) return;

  for ( l=0; l<list->n_data; l++ ) {
    for ( i=0; i<list->data[l].n_selected_data; ) {
      if ( list->data[l].data[i].n_prev >= neighbors ) {
	tmp = list->data[l].data[ list->data[l].n_selected_data-1 ];
	list->data[l].data[ list->data[l].n_selected_data-1 ] = list->data[l].data[ i ];
	list->data[l].data[ i ] = tmp;
	list->data[l].n_selected_data --;
      }
      else
	i++;
    }
  }
}



void removeStartBoundChain( chainList *list, int index_first_circle ) 
{
  chainType tmp;  
  int l;

  for ( l=0; l<list->n_selected_data; ) {
    if ( list->data[l].data[0]->index <= index_first_circle ) {
      tmp = list->data[ list->n_selected_data-1 ];
      list->data[ list->n_selected_data-1 ] = list->data[l];
      list->data[l] = tmp;
      list->n_selected_data --;
    }
    else
      l++;
  }
}



void removeEndBoundChain( chainList *list, int index_last_circle ) 
{
  chainType tmp;
  int l;

  for ( l=0; l<list->n_selected_data; ) {
    if ( list->data[l].data[ list->data[l].n_data-1 ]->index >= index_last_circle ) {
      tmp = list->data[ list->n_selected_data-1 ];
      list->data[ list->n_selected_data-1 ] = list->data[l];
      list->data[l] = tmp;
      list->n_selected_data --;
     }
    else
      l++;
  }
}




/************************************************************
 *
 * BUILDING CHAINS
 *
 ************************************************************/



static circleType* bestNext( circleType *ref, circleType *c ) 
{
  int i, f;
  int score, bestscore;
  circleType *best = (circleType *)NULL;


  /* on cherche un point dans les suivants
     qui n'appartient a aucune chaine
  */
  for ( i=0, f=-1; f==-1 && i<c->n_next; i++ ) {
    if ( c->next[i]->n_chains == 0 ) f = i;
  }

  if ( f == -1 )
    return( (circleType *)NULL );
  
  /* ce point donne un score de reference
   */
  best = c->next[f];
  bestscore = (ref->x - c->next[f]->x)*(ref->x - c->next[f]->x) 
    + (ref->y - c->next[f]->y)*(ref->y - c->next[f]->y); 
  
  /* on cherche le point n'appartenant a aucune chaine 
     ayant le meilleur score
  */
  for ( i=0; i<c->n_next; i++ ) {
    if ( c->next[i]->n_chains > 0 ) continue;
    score = (ref->x - c->next[i]->x)*(ref->x - c->next[i]->x) 
	  + (ref->y - c->next[i]->y)*(ref->y - c->next[i]->y); 
    if ( bestscore > score )  {
      best = c->next[i];
      bestscore = score;
    }
  }
  
  return( best );
}


static int buildChain( chainType *list, circleType *first )
{
  char *proc = "buildChain";
  circleType *c, *next, *ref;
  enumChainBuilding typeChainBuilding = _CLOSEST_TO_PREVIOUS_;

  
  if ( addCircleToChain( list, first ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add first point\n", proc );
    return( -1 );
  }

  ref = c = first;

  while ( (next = bestNext( ref, c )) != (circleType *)NULL ) {
    
    if ( addCircleToChain( list, next ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add next point\n", proc );
      return( -1 );
    }
    
    
    /* go to next point
     */
    c = next;

    /* reference
     */
    switch( typeChainBuilding ) {
    default :
    case _CLOSEST_TO_PREVIOUS_ :
      /* reference = precedent
       */
      ref = c;
      break;
    case _CLOSEST_TO_REF11_ :
      /* reference = precedent "bien forem",
	 i.e. avec un seul voisin avant et apres
      */
      if ( c->n_prev == 1 && c->n_next == 1 )
	ref = c;
      break;
    }


  }

  return( 1 );
}




static int buildChainList( circleListList *readlist, chainList *chainlist)
{
  char *proc = "buildChainList";
  int i, l, nchains;
  chainType *list;  



  for ( nchains=0, l=0; l<readlist->n_data; l++ ) {
    for ( i=0; i<readlist->data[l].n_selected_data; i++ ) {

      /* test sur le nombre de predecesseurs
       */
      /*
      if ( readlist->data[l].data[i].n_prev > 0 )
	continue;
      */

      
      /* le point appartient deja a une chaine
	 => on passe au suivant
       */

      if ( readlist->data[l].data[i].n_chains > 0 )
	continue;


      /* */
      list = getChainTypeFromChainList( chainlist );
      if ( list == (chainType *)NULL ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to get a new chain (chain #%d)\n", proc, nchains );
	return( -1 );
      }

      if ( buildChain( list, &(readlist->data[l].data[i]) ) != 1 )  {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to build the new chain (chain #%d)\n", proc, nchains );
	return( -1 );
      }
      
      nchains++;
    }
  }

  return( 1 );
}





/************************************************************
 *
 * READING CHAINS
 *
 ************************************************************/

int chainListFromCircleListList( circleListList *readlist, chainList *chainlist, 
		       int depth, int margin,
		       int maxForwardNeighbors, int maxBackwardNeighbors )
{
  char *proc = "chainListFromCircleListList";
  
  /* preprocessing
   */

  rmaxCircleListList( readlist );



  /* selection on neighbor numbers
   */
  
  if ( maxForwardNeighbors >= 2 || maxBackwardNeighbors >= 2 ) {
    fprintf( stderr, "... selection on neighbors\n" );
    fprintf( stderr, "    ... looking for  forward neighbors with depth=%d\n", 1 );
    if ( searchForwardNeighbors( readlist, 1, margin ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when searching neighbors\n", proc );
      return( -1 );
    }
    
    fprintf( stderr, "    ... looking for backward neighbors with depth=%d\n", 1 );
    if ( searchBackwardNeighbors( readlist, 1, margin ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when searching neighbors\n", proc );
      return( -1 );
    }

    if ( maxForwardNeighbors >= 2 ) 
      selectOnForwardNeighbors( readlist, maxForwardNeighbors );

    if ( maxBackwardNeighbors >= 2 ) 
      selectOnBackwardNeighbors( readlist, maxBackwardNeighbors );
    
    rmaxCircleListList( readlist );
  }


  /* looking for neighbors 
   */

  fprintf( stderr, "... looking for  forward neighbors with depth=%d\n", depth );
  if ( searchForwardNeighbors( readlist, depth, margin ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when searching neighbors\n", proc );
    return( -1 );
  }

  fprintf( stderr, "... looking for backward neighbors with depth=%d\n", depth );
  if ( searchBackwardNeighbors( readlist, depth, margin ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when searching neighbors\n", proc );
    return( -1 );
  }

  if ( depth > 1 ) {

    fprintf( stderr, "    ... updating  forward neighbors with depth=%d\n", 1 );
    if ( searchForwardNeighbors( readlist, 1, margin ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when searching neighbors\n", proc );
      return( -1 );
    }
    
    fprintf( stderr, "    ... updating backward neighbors with depth=%d\n", 1 );
    if ( searchBackwardNeighbors( readlist, 1, margin ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when searching neighbors\n", proc );
      return( -1 );
    }
    
    rmaxCircleListList( readlist );

  }



  /* construction des chaines
   */
  if ( buildChainList( readlist, chainlist ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "error while building chains\n" );
    return( -1 );
  }

  /* post-processing
   */
  computeChainListBoundingBoxes( chainlist );

  return( 1 );
}









void printStatsChainList( FILE *f, circleListList *readlist, chainList *chainlist, char *s )
{
  int l;
  int ncircle = 0;
  int nscircle = 0;

  for ( l=0; l<readlist->n_data; l++ ) 
    ncircle += readlist->data[l].n_data;

  for ( l=0; l<readlist->n_data; l++ ) 
    nscircle += readlist->data[l].n_selected_data;

  fprintf( f, "\n" );
  fprintf( f, "---- some statistics" );
  if ( s != (char*)NULL ) fprintf( f, " on %s", s ); 
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  
  fprintf( f, "    number of samples:      %d\n", ncircle );
  fprintf( f, "    number of selected samples:      %d\n", nscircle );
  fprintf( f, "    maximal radius:         %d\n", readlist->rmax );
  fprintf( f, "    number of built chains: %d\n", chainlist->n_data );
  fprintf( f, "    number of selected chains: %d\n", chainlist->n_selected_data );
  fprintf( f, "\n" );

  statsCircleRadius( f, readlist );
  statsForwardNeighbors( f, readlist );
  statsBackwardNeighbors( f, readlist );
  statsNeighbors( f, readlist );
  statsMembershipChains( f, readlist, chainlist );
}








/************************************************************
 *
 * COLOCALIZATION
 *
 ************************************************************/



/* L'idee est de regarder les colocalisations d'abord et 
   de construire les chaines
*/
int colocalizationFromLists( circleListList *colocalizationlist, 
			     chainList *chainlist, 
			     circleListList *readlist1, circleListList *readlist2, 
			     int depth, int margin,
			     int maxForwardNeighbors, int maxBackwardNeighbors )
{
  char *proc = "colocalizationFromLists";
  int l, i, j;
  circleType c;
  circleList *list2, *list1, *reslist;
  int rmax, dx, dy, r;
  
  /* allocation de la liste resultat
   */
  l = (readlist2->n_data > readlist1->n_data) ? readlist2->n_data : readlist1->n_data;
  if ( allocCircleListList( colocalizationlist, l ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }


  
  for ( l=0; l<readlist2->n_data && l<readlist1->n_data; l++ ) {
    list1 = &(readlist1->data[l]);
    list2 = &(readlist2->data[l]);
    reslist = &(colocalizationlist->data[l]);

    rmax = list1->rmax + list2->rmax + margin;

    for ( i=0; i<list1->n_data; i++ ) {
      for ( j=0; j<list2->n_data; j++ ) {

	dx = list1->data[i].x - list2->data[j].x;
	if ( dx < -rmax || rmax < dx ) continue;
	dy = list1->data[i].y - list2->data[j].y;
	if ( dy < -rmax || rmax < dy ) continue;
	
	r = list1->data[i].r + list2->data[j].r + margin;

	/* found a colocalization
	 */
	if ( dx*dx+dy*dy <= r*r ) {
	  initCircleType( &c );
	  c.x = list2->data[j].x;
	  c.y = list2->data[j].y;
	  c.r = list2->data[j].r;
	  c.index = list2->data[j].index;
	  c.type = _READ_;
	  if ( addCircleToCircleList( reslist, &c) != 1 ) {
	    if ( _verbose_ )
	      fprintf( stderr, "%s: unable to add measure #%d of list #%d to colocalization list\n", proc, j, l );
	    return( -1 );
	  }
	}

      }
    }
  }

  if ( chainListFromCircleListList( colocalizationlist, chainlist, 
				    depth, margin,
				    maxForwardNeighbors, maxBackwardNeighbors ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: chain construction error\n", proc );
    return( -1 );
  }

  return( 1 );
}





/* L'idee est de construire les chaines d'abord et 
   de regarder la colocalisation entre les chaines
*/
int colocalizationFromChains( chainList *coloclist, 
			      chainList *chainlist1, chainList *chainlist2, 
			      int margin )
{
  char *proc = "colocalizationFromChains";
  int l1, l2;
  int nc=0, ncandidates = 0;
  chainType **candidates = (chainType **)NULL;


  /* pre-processing
   */
   for ( l2=0; l2<chainlist2->n_selected_data; l2++ ) {
     for ( nc=0, l1=0; l1<chainlist1->n_selected_data; l1++ ) {
       if ( chainlist1->data[l1].xmax < chainlist2->data[l2].xmin - margin
	    || chainlist2->data[l2].xmax < chainlist1->data[l1].xmin - margin
	    || chainlist1->data[l1].ymax < chainlist2->data[l2].ymin - margin
	    || chainlist2->data[l2].ymax < chainlist1->data[l1].ymin - margin )
	 continue;
       if ( chainlist1->data[l1].data[ chainlist1->data[l1].n_data-1 ]->index < chainlist2->data[l2].data[ 0 ]->index
	    || chainlist2->data[l2].data[ chainlist2->data[l2].n_data-1 ]->index < chainlist1->data[l1].data[ 0 ]->index )
	 continue;
       nc ++;
     }
     if ( ncandidates < nc ) ncandidates = nc;
   }

   fprintf( stderr, "... maximal number of possible chain intersection = %d\n", ncandidates );
   
   if ( ncandidates <= 0 ) {
     if ( _verbose_ )
       fprintf( stderr, "%s: no possible colocalization\n", proc );
     return( -1 );
   }

   candidates = (chainType **)malloc( ncandidates * sizeof(chainType *) );
   if ( candidates == (chainType **)NULL ) {
     if ( _verbose_ )
       fprintf( stderr, "%s: allocation error\n", proc );
     return( -1 );
   }
   


   /* colocalization
    */

   for ( l2=0; l2<chainlist2->n_selected_data; l2++ ) {

     /* chain candidates for colocalization
      */
     for ( nc=0, l1=0; l1<chainlist1->n_selected_data; l1++ ) {
       if ( chainlist1->data[l1].xmax < chainlist2->data[l2].xmin - margin
	    || chainlist2->data[l2].xmax < chainlist1->data[l1].xmin - margin
	    || chainlist1->data[l1].ymax < chainlist2->data[l2].ymin - margin
	    || chainlist2->data[l2].ymax < chainlist1->data[l1].ymin - margin )
	 continue;
       if ( chainlist1->data[l1].data[ chainlist1->data[l1].n_data-1 ]->index < chainlist2->data[l2].data[ 0 ]->index
	    || chainlist2->data[l2].data[ chainlist2->data[l2].n_data-1 ]->index < chainlist1->data[l1].data[ 0 ]->index )
	 continue;
       candidates[ nc ] = &(chainlist1->data[l1]);
       nc ++;
     }
     
     /* comparing chainlist2->data[l2] with *candidates
	A finir ...
      */
     ;
   }

   return( 1 );
}
