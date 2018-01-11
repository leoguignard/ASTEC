/*************************************************************************
 * bal-biolib-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 17 jan 2014 11:07:58 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <fcntl.h>
#include <unistd.h>

#include <histogram.h>

#include <bal-biolib-tools.h>


static int _debug_ = 0;
static int _verbose_ = 1;


static int _size_to_be_allocated_ = 100;
#define _LINE_LENGTH_ 1024







/************************************************************
 *
 * Detection  structure
 *
 ************************************************************/


void BAL_InitBlDetection( bal_blDetection *d )
{
  d->voxelcenter.x = 0;
  d->voxelcenter.y = 0;
  d->voxelcenter.z = 0;

  d->unit = UNDEF_UNIT;

  d->voxelhalfradius1 = 0;
  d->voxelhalfradius2 = 0;

  d->arg1 = 0;
  d->arg2 = 0;
  d->arg3 = 0;
  d->arg4 = 0;
  d->arg5 = 0;
  
  /* frame index
   */
  d->imageindex = -1;

  /* detection list index
   */
  d->identifier = -1;

  /*lemon index
   */
  d->idnode = -1;
}















/************************************************************
 *
 * Pointer detection list
 *
 ************************************************************/

void BAL_InitBlPointerDetectionList( bal_blPointerDetectionList *l )
{
  l->data = (bal_blDetection **)NULL;
  l->n = 0;
  l->n_allocated = 0;

  l->unit = UNDEF_UNIT;
}

int BAL_AllocBlPointerDetectionList( bal_blPointerDetectionList *l, int n ) 
{
  char *proc = "BAL_AllocBlPointerDetectionList";
  int i;
  l->data = (bal_blDetection **)malloc( n * sizeof(bal_blDetection*) );
  if ( l->data == (bal_blDetection **)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  l->n_allocated = n;
  l->n = 0;
  for ( i=0; i<n; i++ )
    l->data[i] = (bal_blDetection *)NULL;
  return( 0 );
}

extern void BAL_FreeBlPointerDetectionList( bal_blPointerDetectionList *l )
{
  if ( l->data != (bal_blDetection **)NULL )
    free( l->data );
  BAL_InitBlPointerDetectionList( l );
}

static int _addBlDetectionToPointerList( bal_blDetection *d, bal_blPointerDetectionList *l )
{
  char *proc = "_addBlDetectionToPointerList";
  bal_blDetection **data;
  int s = l->n_allocated;

  if ( l->n == l->n_allocated ) {
    s += _size_to_be_allocated_;
    data = (bal_blDetection **)malloc( s * sizeof(bal_blDetection*) );
    if ( data == (bal_blDetection **)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    if ( l->n_allocated > 0 && l->data != (bal_blDetection **)NULL ) {
      (void)memcpy( data, l->data, l->n_allocated*sizeof(bal_blDetection*) );
      free( l->data );
    }
    l->n_allocated = s;
    l->data = data;
  }
  
  l->data[l->n] = d;
  
  l->n ++;

  return( 0 );

}

void BAL_PrintBlPointerDetectionList( FILE *f, bal_blPointerDetectionList *l )
{
  fprintf( f, "data length = %d, data allocated = %d\n",
	   l->n, l->n_allocated );
  fprintf( f, "   pointer=%p\n", (void*)l->data );
}










/************************************************************
 *
 * search dedicated structure
 *
 ************************************************************/



void BAL_InitBlNeighborsSearch(  bal_blNeighborsSearch *s )
{
  s->buf = (bal_blPointerDetectionList *)NULL;
  s->array = (bal_blPointerDetectionList ***)NULL;

  s->lcorner.x = 0;
  s->lcorner.y = 0;
  s->lcorner.z = 0;

  s->rcorner.x = 0;
  s->rcorner.y = 0;
  s->rcorner.z = 0;

  s->bucketsizex = 0.0;
  s->bucketsizey = 0.0;
  s->bucketsizez = 0.0;

  s->dimx = 0;
  s->dimy = 0;
  s->dimz = 0;
}



int BAL_AllocBlNeighborsSearch(  bal_blNeighborsSearch *s )
{
  char *proc = "BAL_AllocBlNeighborsSearch";
  int size = s->dimx * s->dimy * s->dimz;
  int asize = 0;
  int i, j;
  bal_blPointerDetectionList ***z;
  bal_blPointerDetectionList **zy;
  bal_blPointerDetectionList *zyx;

  if ( s->dimx <= 0 || s->dimy <= 0 || s->dimz <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird dimensions\n", proc );
    return( -1 );
  }

  s->buf = (bal_blPointerDetectionList *)malloc( size * sizeof( bal_blPointerDetectionList ) );
  if ( s->buf ==  (bal_blPointerDetectionList *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: buffer allocation error\n", proc );
    return( -1 );
  }

  asize += s->dimy * s->dimz * sizeof( bal_blPointerDetectionList* );
  asize += s->dimz      * sizeof( bal_blPointerDetectionList** );

  s->array = (bal_blPointerDetectionList ***)malloc( asize );
  if ( s->array ==  (bal_blPointerDetectionList ***)NULL ) {
    free( s->buf );
    s->buf = (bal_blPointerDetectionList *)NULL;
    if ( _verbose_ )
      fprintf( stderr, "%s: array allocation error\n", proc );
    return( -1 );
  }

  /* pointer addresses calculation
   */
  z = (bal_blPointerDetectionList ***)(s->array);
  zy = (bal_blPointerDetectionList **)(z + s->dimz);
  zyx = (bal_blPointerDetectionList *)(s->buf);
  for ( i = 0; i < s->dimz; i ++ ) {
    *z = zy;
    z ++;
    for ( j = 0; j < s->dimy; j ++ ) {
      *zy = zyx;
      zy ++;
      zyx += s->dimx;
    }
  }

  for ( i = 0; i < size; i ++ ) {
    BAL_InitBlPointerDetectionList( &(s->buf[i]) );
  }

  return( 0 );
}





void BAL_FreeBlNeighborsSearch(  bal_blNeighborsSearch *s )
{
  int size = s->dimx*s->dimy*s->dimz;
  int i;

  if ( s->dimx > 0 &&  s->dimy > 0 &&  s->dimz > 0 
       && s->buf != (bal_blPointerDetectionList *)NULL ) {
    for ( i = 0; i < size; i ++ ) {
      BAL_FreeBlPointerDetectionList( &(s->buf[i]) );
    }
  }
  
  if ( s->buf != (bal_blPointerDetectionList *)NULL )
    free( s->buf );

  if ( s->array != (bal_blPointerDetectionList ***)NULL )
    free( s->array );

  BAL_InitBlNeighborsSearch( s );
}



void BAL_PrintBlNeighborsSearch( FILE *f, bal_blNeighborsSearch *s )
{
  int i, j, k;
  int n;

  fprintf( f, "buckets dimension = [%d %d %d], bucket size = [%f %f %f]\n",
	   s->dimx, s->dimy, s->dimz, s->bucketsizex, s->bucketsizey, s->bucketsizez );

  n = 0;
  for ( k=0; k<s->dimz; k++ )
    for ( j=0; j<s->dimy; j++ )
      for ( i=0; i<s->dimx; i++ ) {
	n += s->array[k][j][i].n;
      }
  fprintf( f, "   #points in array=%d\n", n );

  n = 0;
  for ( k=0; k<s->dimz * s->dimy * s->dimx; k++ )
    n += s->buf[k].n;
  fprintf( f, "   #points in buf=%d\n", n );

}










/************************************************************
 *
 * Detection list
 *
 ************************************************************/



void BAL_InitBlDetectionList( bal_blDetectionList *l )
{
  l->data = (bal_blDetection *)NULL;
  l->n = 0;
  l->n_allocated = 0;

  l->unit = UNDEF_UNIT;

  l->lcorner.x = -1.0;
  l->lcorner.y = -1.0;
  l->lcorner.z = -1.0;
  l->rcorner.x = -1.0;
  l->rcorner.y = -1.0;
  l->rcorner.z = -1.0;

  BAL_InitBlNeighborsSearch( &(l->buckets) );
}



int BAL_AllocBlDetectionList( bal_blDetectionList *l, int n )
{
  char *proc = "BAL_AllocBlDetectionList";
  int i;
  l->data = (bal_blDetection *)malloc( n * sizeof(bal_blDetection) );
  if ( l->data == (bal_blDetection *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  l->n_allocated = n;
  l->n = 0;
  for ( i=0; i<n; i++ )
    BAL_InitBlDetection( &(l->data[i]) );
  return( 0 );
}



void BAL_FreeBlDetectionList( bal_blDetectionList *l )
{
  BAL_FreeBlNeighborsSearch( &(l->buckets) );
  if ( l->data != (bal_blDetection *)NULL )
    free( l->data );
  BAL_InitBlDetectionList( l );
}




static int _addBlDetectionToList( bal_blDetection *d, bal_blDetectionList *l )
{
  char *proc = "_addBlDetectionToList";
  bal_blDetection *data;
  int s = l->n_allocated;

  if ( l->n == l->n_allocated ) {
    s += _size_to_be_allocated_;
    data = (bal_blDetection *)malloc( s * sizeof(bal_blDetection) );
    if ( data == (bal_blDetection *)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    if ( l->n_allocated > 0 && l->data != (bal_blDetection *)NULL ) {
      (void)memcpy( data, l->data, l->n_allocated*sizeof(bal_blDetection) );
      free( l->data );
    }
    l->n_allocated = s;
    l->data = data;
  }
  
  l->data[l->n] = *d;


  /* information from the structure
     frame/image index
     unit
  */
  l->data[l->n].unit = l->unit;

  /* identifier = list index
   */
  l->data[l->n].identifier = l->n;


  l->n ++;

  return( 0 );

}



int BAL_InitNeighborsSearchBlDetectionList(  bal_blDetectionList *l, bal_doublePoint* bucketsize )
{
  char *proc = "BAL_InitNeighborsSearchBlDetectionList";
  int n, i, j, k;

  if ( bucketsize->x <= 0.0 || bucketsize->y <= 0.0 || bucketsize->z <= 0.0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: negative dimensions for bucket sizes\n", proc );
    return( -1 );
  }

  BAL_BlDetectionListBoundingBox( l );
  
  l->buckets.dimx = (int)( (l->rcorner.x - l->lcorner.x) / bucketsize->x );
  l->buckets.dimy = (int)( (l->rcorner.y - l->lcorner.y) / bucketsize->y );
  l->buckets.dimz = (int)( (l->rcorner.z - l->lcorner.z) / bucketsize->z );
  
  if ( l->buckets.dimx <= 0 ) l->buckets.dimx = 1;
  if ( l->buckets.dimy <= 0 ) l->buckets.dimy = 1;
  if ( l->buckets.dimz <= 0 ) l->buckets.dimz = 1;
  
  l->buckets.bucketsizex = ( l->buckets.dimx > 1 ) ? (l->rcorner.x - l->lcorner.x) / (double)l->buckets.dimx : 1.0;
  l->buckets.bucketsizey = ( l->buckets.dimy > 1 ) ? (l->rcorner.y - l->lcorner.y) / (double)l->buckets.dimy : 1.0;
  l->buckets.bucketsizez = ( l->buckets.dimz > 1 ) ? (l->rcorner.z - l->lcorner.z) / (double)l->buckets.dimz : 1.0;

  if ( BAL_AllocBlNeighborsSearch( &(l->buckets) ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }


  for ( n=0; n<l->n; n++ ) {
    i = (l->data[n].voxelcenter.x - l->lcorner.x) / l->buckets.bucketsizex;
    if ( i >= l->buckets.dimx ) i = l->buckets.dimx-1;
    j = (l->data[n].voxelcenter.y - l->lcorner.y) / l->buckets.bucketsizey;
    if ( j >= l->buckets.dimy ) j = l->buckets.dimy-1;
    k = (l->data[n].voxelcenter.z - l->lcorner.z) / l->buckets.bucketsizez;
    if ( k >= l->buckets.dimz ) k = l->buckets.dimz-1;
    
    if ( _addBlDetectionToPointerList( &(l->data[n]), &(l->buckets.array[k][j][i]) ) != 0 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add detection (%f,%f,%f) to bucket [%d,%d,%d]\n", proc,
		 l->data[n].voxelcenter.x, l->data[n].voxelcenter.y, l->data[n].voxelcenter.z,
		 i, j, k );
      return( -1 );
    }

  }


  return( 1 );
}










/************************************************************
 *
 * Detection list I/O
 *
 ************************************************************/



int BAL_ReadBlDetectionList( bal_blDetectionList *l, char *name )
{
  char *proc = "BAL_ReadBlDetectionList";
  FILE *f;
  bal_blDetection d;
  int nread, n = 0;
  char line[_LINE_LENGTH_];

  f = fopen( name, "r" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening %s\n", proc, name);
    return( -1 );
  }

  BAL_InitBlDetection( &d );
  
  while ( fgets( line, _LINE_LENGTH_, f ) != (char*)NULL ) {
    nread = sscanf( line ,"%lf\t%lf\t%lf\t%d\t%d\%lf\t%d\t%d\t%lf\t%lf",
		    &(d.voxelcenter.x), &(d.voxelcenter.y), &(d.voxelcenter.z), 
		    &(d.voxelhalfradius1), &(d.voxelhalfradius2), 
		    &(d.arg1), &(d.arg2), &(d.arg3), &(d.arg4), &(d.arg5 ) );
    switch ( nread ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: found '%d' args when reading line #%d\n", proc, nread, n );
      break;
    case -1 : /* the end */
      break;
    case 6 : /* older format */
    case 10 :
      if ( _addBlDetectionToList( &d, l ) != 0 ) {
	fclose( f );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when adding detection #%d\n", proc, n );
	return( -1 );
      }
      n ++;
      break;
    }
  }


  fclose( f );

  if ( _verbose_ >= 2 ) 
    fprintf( stderr, "%s: has read %d detections in '%s'\n", proc, n, name );

  return( 0 );
}



int BAL_WriteBlDetectionList( bal_blDetectionList *l, char *name )
{
  char *proc = "BAL_WriteBlDetectionList";
  FILE *f;
  int i;

  f = fopen( name, "w" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening %s\n", proc, name);
    return( -1 );
  }
  
  for (i=0; i<l->n; i++ ) {
    if ( fprintf( f, "%g\t%g\t%g\t%d\t%d\t%g\t%d\t%d\t%g\t%g\t\n",
		  l->data[i].voxelcenter.x, l->data[i].voxelcenter.y, l->data[i].voxelcenter.z, 
		  l->data[i].voxelhalfradius1, l->data[i].voxelhalfradius2, 
		  l->data[i].arg1, l->data[i].arg2, l->data[i].arg3, 
		  l->data[i].arg4, l->data[i].arg5  ) < 0 ) {
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when writing detection #%d\n", proc, i );
      return( -1 );
    }
  }
  
  fclose( f );

  return( 0 );
}










/************************************************************
 *
 * Detection list properties
 *
 ************************************************************/



void BAL_BlDetectionListBoundingBox( bal_blDetectionList *t )
{
  int i;
  t->lcorner.x = t->data[0].voxelcenter.x;
  t->lcorner.y = t->data[0].voxelcenter.y;
  t->lcorner.z = t->data[0].voxelcenter.z;

  t->rcorner.x = t->data[0].voxelcenter.x;
  t->rcorner.y = t->data[0].voxelcenter.y;
  t->rcorner.z = t->data[0].voxelcenter.z;

  for ( i=1; i<t->n; i++ ) {

    if ( t->lcorner.x > t->data[i].voxelcenter.x ) 
      t->lcorner.x = t->data[i].voxelcenter.x;
    if ( t->lcorner.y > t->data[i].voxelcenter.y ) 
      t->lcorner.y = t->data[i].voxelcenter.y;
    if ( t->lcorner.z > t->data[i].voxelcenter.z ) 
      t->lcorner.z = t->data[i].voxelcenter.z;
    if ( t->rcorner.x < t->data[i].voxelcenter.x ) 
      t->rcorner.x = t->data[i].voxelcenter.x;
    if ( t->rcorner.y < t->data[i].voxelcenter.y ) 
      t->rcorner.y = t->data[i].voxelcenter.y;
    if ( t->rcorner.z < t->data[i].voxelcenter.z ) 
      t->rcorner.z = t->data[i].voxelcenter.z;
 }

}










/************************************************************
 *
 * Detection neighbor search
 *
 ************************************************************/


int _compareBlDetectionWrtZ( const void * a, const void * b )
{
  bal_blDetection *data = (bal_blDetection *)a;
  bal_blDetection *datb = (bal_blDetection *)b;
  if ( data->voxelcenter.z < datb->voxelcenter.z ) return( -1 );
  if ( data->voxelcenter.z > datb->voxelcenter.z ) return( 1 );
  return( 0 );
}

int _compareBlDetectionWrtY( const void * a, const void * b )
{
  bal_blDetection *data = (bal_blDetection *)a;
  bal_blDetection *datb = (bal_blDetection *)b;
  if ( data->voxelcenter.y < datb->voxelcenter.y ) return( -1 );
  if ( data->voxelcenter.y > datb->voxelcenter.y ) return( 1 );
  return( 0 );
}

int _compareBlDetectionWrtX( const void * a, const void * b )
{
  bal_blDetection *data = (bal_blDetection *)a;
  bal_blDetection *datb = (bal_blDetection *)b;
  if ( data->voxelcenter.x < datb->voxelcenter.x ) return( -1 );
  if ( data->voxelcenter.x > datb->voxelcenter.x ) return( 1 );
  return( 0 );
}





int _bucketNeighborsSearch( bal_blDetection *d,
			    bal_blDetectionList *l,
			    bal_blPointerDetectionList *neighbors,
			    bal_doublePoint* halfneighborhood )
{
  char *proc = "_bucketNeighborsSearch";

  double x = d->voxelcenter.x;
  double y = d->voxelcenter.y;
  double z = d->voxelcenter.z;

  double wx = halfneighborhood->x;
  double wy = halfneighborhood->y;
  double wz = halfneighborhood->z;

  int imin, imax;
  int jmin, jmax;
  int kmin, kmax;
  int i, j, k, n;

  bal_blPointerDetectionList *b;

  neighbors->n = 0;

  imin = (x - wx - l->lcorner.x) / l->buckets.bucketsizex;
  if ( imin < 0 ) imin = 0;
  imax = (x + wx - l->lcorner.x) / l->buckets.bucketsizex;
  if ( imax >= l->buckets.dimx ) imax = l->buckets.dimx - 1;

  jmin = (y - wy - l->lcorner.y) / l->buckets.bucketsizey;
  if ( jmin < 0 ) jmin = 0;
  jmax = (y + wy - l->lcorner.y) / l->buckets.bucketsizey;
  if ( jmax >= l->buckets.dimy ) jmax = l->buckets.dimy - 1;

  kmin = (z - wz - l->lcorner.z) / l->buckets.bucketsizez;
  if ( kmin < 0 ) kmin = 0;
  kmax = (z + wz - l->lcorner.z) / l->buckets.bucketsizez;
  if ( kmax >= l->buckets.dimz ) kmax = l->buckets.dimz - 1;

  for ( k=kmin; k<=kmax; k++ )
  for ( j=jmin; j<=jmax; j++ )
  for ( i=imin; i<=imax; i++ ) {

    b = &(l->buckets.array[k][j][i]);

    if ( b->n == 0 ) continue;
    for ( n=0; n<b->n; n ++ ) {
      if ( x - wx > b->data[n]->voxelcenter.x ) continue;
      if ( x + wx < b->data[n]->voxelcenter.x ) continue;
      if ( y - wy > b->data[n]->voxelcenter.y ) continue;
      if ( y + wy < b->data[n]->voxelcenter.y ) continue;
      if ( z - wz > b->data[n]->voxelcenter.z ) continue;
      if ( z + wz < b->data[n]->voxelcenter.z ) continue;
      if ( _addBlDetectionToPointerList( b->data[n], neighbors ) != 0 ) {
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: unable to add point to neighbor list\n", proc );
	}
	return( -1 );
      }
    }

  }

  return( neighbors->n );
}








int _voxelBruteForceNeighborsSearch( bal_blDetection *d,
				     bal_blDetectionList *l,
				     bal_blPointerDetectionList *neighbors,
				     bal_doublePoint* halfneighborhood )
{
  char *proc = "_voxelBruteForceNeighborsSearch";

  double x = d->voxelcenter.x;
  double y = d->voxelcenter.y;
  double z = d->voxelcenter.z;

  double wx = halfneighborhood->x;
  double wy = halfneighborhood->y;
  double wz = halfneighborhood->z;

  int i;

  neighbors->n = 0;

  /*
  fprintf( stdout, "%s: compare point (%f %f %f) in window [%f %f %f]\n",
	   proc, x, y, z, wx, wy ,wz );
  */

  for ( i=0; i<l->n; i++ ) {
    if ( x - wx > l->data[i].voxelcenter.x ) continue;
    if ( x + wx < l->data[i].voxelcenter.x ) continue;
    if ( y - wy > l->data[i].voxelcenter.y ) continue;
    if ( y + wy < l->data[i].voxelcenter.y ) continue;
    if ( z - wz > l->data[i].voxelcenter.z ) continue;
    if ( z + wz < l->data[i].voxelcenter.z ) continue;
    if ( _addBlDetectionToPointerList( &(l->data[i]), neighbors ) != 0 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to add point to neighbor list\n", proc );
      }
      return( -1 );
    }
  }
  return( neighbors->n );
}




int BAL_SearchBlDetectionNeighbors( bal_blDetection *d,
				    bal_blDetectionList *l,
				    bal_blPointerDetectionList *neighbors,
				    bal_doublePoint* halfneighborhood,
				    enumSearchType searchType )
{
  char *proc = "BAL_SearchBlDetectionNeighbors";

  switch( searchType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such search type not handled yet\n", proc );
    return( -1 );
  case _VOXEL_BUCKETS_ :
     return( _bucketNeighborsSearch( d, l, neighbors, halfneighborhood ) );
  case _VOXEL_BRUTEFORCE_ :
    return( _voxelBruteForceNeighborsSearch( d, l, neighbors, halfneighborhood ) );
  }

  return( 1 );
}











/************************************************************
 *
 * Track selection
 *
 ************************************************************/

void BAL_InitBlTrackSelectionCriteria( bal_blTrackSelectionCriteria *c )
{
  c->min_index_start = -1;
  c->max_index_end = -1;
  c->min_index_range = -1;
  
  c->min_fly_distance_2D = -1.0;
  c->min_length_2D = -1.0;
}



/************************************************************
 *
 * Track properties related structures 
 *
 ************************************************************/



void BAL_InitBlTrackProperties( bal_blTrackProperties *p )
{
  p->index_start = -1;
  p->index_end = -1;

  p->time_start = -1;
  p->time_end = -1;
  p->lifetime = 0;

  p->fly_distance_2D = 0.0;
  p->length_2D = 0.0;
  p->fly_distance_3D = 0.0;
  p->length_3D = 0.0;
}



void BAL_ComputeTrackProperties( bal_blTrack *t, 
				 double vx, 
				 double vy,
				 double vz,
				 double dt )
{
  int i;
  double dx, dy, dz;
  bal_blDetection *f = &(t->detectionList.data[ 0 ]);
  bal_blDetection *l = &(t->detectionList.data[ t->detectionList.n-1 ]);

  t->properties.index_start = f->imageindex;
  t->properties.index_end = l->imageindex;
  
  t->properties.time_start = f->imageindex * dt;
  t->properties.time_end = l->imageindex * dt;
  t->properties.lifetime = (t->properties.index_end - t->properties.index_start + 1) * dt;

  dx = l->voxelcenter.x-f->voxelcenter.x;
  dy = l->voxelcenter.y-f->voxelcenter.y;
  dz = l->voxelcenter.z-f->voxelcenter.z;
  
  t->properties.fly_distance_2D = sqrt( dx*dx*vx*vx + dy*dy*vy*vy );
  t->properties.fly_distance_3D = sqrt( dx*dx*vx*vx + dy*dy*vy*vy + dz*dz*vz*vz );

  if ( t->detectionList.n <=1 ) 
    return;

  t->properties.length_2D = 0;
  t->properties.length_3D = 0;

  for ( i=0; i<t->detectionList.n-1; i ++ ) {
    dx = t->detectionList.data[ i+1 ].voxelcenter.x - t->detectionList.data[ i ].voxelcenter.x;
    dy = t->detectionList.data[ i+1 ].voxelcenter.y - t->detectionList.data[ i ].voxelcenter.y;
    dz = t->detectionList.data[ i+1 ].voxelcenter.z - t->detectionList.data[ i ].voxelcenter.z;
    t->properties.length_2D += sqrt( dx*dx*vx*vx + dy*dy*vy*vy );
    t->properties.length_3D += sqrt( dx*dx*vx*vx + dy*dy*vy*vy + dz*dz*vz*vz );
  }

}



void BAL_ComputeTrackListProperties( bal_blTrackList *l,
				     double vx, 
				     double vy,
				     double vz,
				     double dt )
{
  int i;
  for ( i=0; i<l->n; i++ ) {
    BAL_ComputeTrackProperties( &(l->data[i]), vx, vy, vz, dt );
  }
}









/************************************************************
 *
 * Track properties output
 *
 ************************************************************/






int BAL_WriteLifetimeHistogram( FILE *f, int fd,
			      bal_blTrackList *list )
{
  char *proc = "BAL_WriteLifetimeHistogram";
  int i;
  double *l = (double*)NULL;

  typeHistogram h;
  unionValues min, max;
  
  l = (double*)malloc( list->n * sizeof(double) );
  if ( l == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate double array\n", proc );
    return( -1 );
  }
  
  min.val_r64 = max.val_r64 = list->data[0].properties.lifetime;
  for ( i=0; i<list->n; i++ ) {
    l[i] = list->data[i].properties.lifetime;
    if ( min.val_r64 > l[i] ) min.val_r64 = l[i];
    else if ( max.val_r64 < l[i] ) max.val_r64 = l[i];
  }
  
  initHistogram( &h );
  if ( allocHistogramHeader( &(h.xaxis), &min, &max, -1.0, DOUBLE ) != 1 ) {
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram header\n", proc );
    return( -1 );
  }

  if ( allocHistogramData( &h, DOUBLE ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram\n", proc );
    return( -1 );
  }
  
  if ( fill1DHistogramFromBuffer( &h, (void*)l, DOUBLE, list->n ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to fill histogram\n", proc );
    return( -1 );
  }

  free( l );
  
  if ( fprintf1DHistogramScilab( f, fd, &h, "lifetime"  ) != 1 ) {
    freeHistogram( &h );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write histogram\n", proc );
    return( -1 );
  }

  freeHistogram( &h );
  return( 1 );
}








int BAL_WriteDistanceHistogram( FILE *f, int fd,
			      bal_blTrackList *list )
{
  char *proc = "BAL_WriteDistanceHistogram";
  int i;
  double *l = (double*)NULL;

  typeHistogram h;
  unionValues min, max;
  
  l = (double*)malloc( list->n * sizeof(double) );
  if ( l == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate double array\n", proc );
    return( -1 );
  }
  
  min.val_r64 = max.val_r64 = list->data[0].properties.fly_distance_2D;
  for ( i=0; i<list->n; i++ ) {
    l[i] = list->data[i].properties.fly_distance_2D;
    if ( min.val_r64 > l[i] ) min.val_r64 = l[i];
    else if ( max.val_r64 < l[i] ) max.val_r64 = l[i];
  }
  
  initHistogram( &h );
  if ( allocHistogramHeader( &(h.xaxis), &min, &max, -1.0, DOUBLE ) != 1 ) {
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram header\n", proc );
    return( -1 );
  }

  if ( allocHistogramData( &h, DOUBLE ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram\n", proc );
    return( -1 );
  }
  
  if ( fill1DHistogramFromBuffer( &h, (void*)l, DOUBLE, list->n ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to fill histogram\n", proc );
    return( -1 );
  }

  free( l );
  
  if ( fprintf1DHistogramScilab( f, fd, &h, "distance"  ) != 1 ) {
    freeHistogram( &h );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write histogram\n", proc );
    return( -1 );
  }

  freeHistogram( &h );
  return( 1 );
}








int BAL_WriteLengthHistogram( FILE *f, int fd,
			      bal_blTrackList *list )
{
  char *proc = "BAL_WriteLengthHistogram";
  int i;
  double *l = (double*)NULL;

  typeHistogram h;
  unionValues min, max;
  
  l = (double*)malloc( list->n * sizeof(double) );
  if ( l == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate double array\n", proc );
    return( -1 );
  }
  
  min.val_r64 = max.val_r64 = list->data[0].properties.length_2D;
  for ( i=0; i<list->n; i++ ) {
    l[i] = list->data[i].properties.length_2D;
    if ( min.val_r64 > l[i] ) min.val_r64 = l[i];
    else if ( max.val_r64 < l[i] ) max.val_r64 = l[i];
  }
  
  initHistogram( &h );
  if ( allocHistogramHeader( &(h.xaxis), &min, &max, -1.0, DOUBLE ) != 1 ) {
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram header\n", proc );
    return( -1 );
  }

  if ( allocHistogramData( &h, DOUBLE ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram\n", proc );
    return( -1 );
  }
  
  if ( fill1DHistogramFromBuffer( &h, (void*)l, DOUBLE, list->n ) != 1 ) {
    freeHistogram( &h );
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to fill histogram\n", proc );
    return( -1 );
  }

  free( l );
  
  if ( fprintf1DHistogramScilab( f, fd, &h, "length"  ) != 1 ) {
    freeHistogram( &h );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write histogram\n", proc );
    return( -1 );
  }

  freeHistogram( &h );
  return( 1 );
}





int BAL_WriteStartEndTimesScilab( FILE *f, int fd,
				   bal_blTrackList *list )
{
  char *proc = "BAL_WriteStartEndTimesScilab";
  int i;
  double *s = (double*)NULL;
  double *e = (double*)NULL;

  s = e = (double*)malloc( 2*list->n * sizeof(double) );
  if ( s == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate double array\n", proc );
    return( -1 );
  }
  e += list->n;
  for ( i=0; i<list->n; i++ ) {
    s[i] = list->data[i].properties.time_start;
    e[i] = list->data[i].properties.time_end;
  }
  
  if ( write( fd, s, 2*list->n * sizeof(double) ) == -1 ) {
    free( s );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing raw data\n", proc );
    return( -1 );
  }
  free( s );
  

  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "//////////////////////////////////////////////////////////\n" );
  fprintf( f, "//\n" );
  fprintf( f, "// distance (start to end point) versus curvilinear length\n" );
  fprintf( f, "//\n" );
  fprintf( f, "//////////////////////////////////////////////////////////\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "STARTTIME=mget( %d, 'd', f);\n", list->n );
  fprintf( f, "ENDTIME=mget( %d, 'd', f);\n", list->n );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background=color(\"white\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
  fprintf( f, "a.font_size = 4;\n" );
  fprintf( f, "a.font_style = 8;\n" );
  fprintf( f, "xlabel( 'start time', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "ylabel( 'end time', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );
  fprintf( f, "plot( STARTTIME, ENDTIME, 'o' );\n" );
  fprintf( f, "a.data_bounds(1,1) = -10;\n" );
  fprintf( f, "\n" );

  fprintf( f, "\n" );
  fprintf( f, "\n" );

  return( 1 );
}




int BAL_WriteLengthDistanceScilab( FILE *f, int fd,
				   bal_blTrackList *list )
{
  char *proc = "BAL_WriteLengthDistanceScilab";
  int i;
  double *l = (double*)NULL;
  double *d = (double*)NULL;

  d = l = (double*)malloc( 2*list->n * sizeof(double) );
  if ( l == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate double array\n", proc );
    return( -1 );
  }
  d += list->n;
  for ( i=0; i<list->n; i++ ) {
    l[i] = list->data[i].properties.length_2D;
    d[i] = list->data[i].properties.fly_distance_2D;
  }
  
  if ( write( fd, l, 2*list->n * sizeof(double) ) == -1 ) {
    free( l );
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing raw data\n", proc );
    return( -1 );
  }
  free( l );
  

  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "//////////////////////////////////////////////////////////\n" );
  fprintf( f, "//\n" );
  fprintf( f, "// distance (start to end point) versus curvilinear length\n" );
  fprintf( f, "//\n" );
  fprintf( f, "//////////////////////////////////////////////////////////\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "LENGTH2D=mget( %d, 'd', f);\n", list->n );
  fprintf( f, "DISTANCE2D=mget( %d, 'd', f);\n", list->n );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background=color(\"white\");\n" );
  fprintf( f, "a=gca();\n" );
  fprintf( f, "set(a,\"auto_clear\",\"off\");\n" );
  fprintf( f, "a.font_size = 4;\n" );
  fprintf( f, "a.font_style = 8;\n" );
  fprintf( f, "xlabel( 'shortest distance', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "ylabel( 'tracking length', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );
  fprintf( f, "plot( DISTANCE2D, LENGTH2D, 'o' );\n" );
  fprintf( f, "\n" );

  fprintf( f, "plot( [0 a.data_bounds(2,1)], [0 a.data_bounds(2,1)], 'r');\n" );
  fprintf( f, "xstring( a.data_bounds(2,1)-15, a.data_bounds(2,1), '$d/l = 1$', 0, 0);\n" );
  fprintf( f, "t=get(\"hdl\");\n" );
  fprintf( f, "t.font_foreground=1;\n" ); 
  fprintf( f, "t.font_size=4;\n" );
  fprintf( f, "t.font_style=8;\n" );
  fprintf( f, "\n" );

  fprintf( f, "plot( [0 a.data_bounds(2,1)], [0 a.data_bounds(2,1)/ 0.2], 'r');\n" );
  fprintf( f, "xstring( a.data_bounds(2,1)-15, a.data_bounds(2,1)/0.2, '$d/l = 0.2$', 0, 0);\n" );
  fprintf( f, "t=get(\"hdl\");\n" );
  fprintf( f, "t.font_foreground=1;\n" ); 
  fprintf( f, "t.font_size=4;\n" );
  fprintf( f, "t.font_style=8;\n" );
  fprintf( f, "\n" );

  fprintf( f, "// a.data_bounds(2,1) = a.data_bounds(2,2);\n" );
  fprintf( f, "\n" );

  fprintf( f, "\n" );
  fprintf( f, "\n" );

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





static int _PrintAllStatsXxxlab( char *name, bal_blTrackList *list, 
				 enumHistogramFile xxxlab )
{
  char *proc = "_PrintAllStatsXxxlab";
  char *defaultname = "statistics";
  char *template;
  char *filename = (char*)NULL;
  FILE *f;
  int fd;

  template = ( name != (char*)NULL ) ? name : defaultname;
  filename = (char*)malloc( strlen( template ) + 5 );
  if ( filename == (char*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate file name\n", proc );
    return( -1 );
  }

  /* open files
   */
  sprintf( filename, "%s.raw", template );

  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    free( filename );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  switch( xxxlab ) {
  default :
    free( filename );
    close( fd );
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not known\n", proc );
    return( -1 );
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
    return( -1 );
  }

  
  free( filename );

  
  /* write data and script
   */
  switch( xxxlab ) {
  default :
    close( fd );
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: such output type not known\n", proc );
    return( -1 );

  case _MATLAB_ :

    fprintf( f, "\n" );
    fprintf( f, "f=fopen('%s.raw','r');\n", _BaseName( template ) );
    fprintf( f, "\n" );

    /* stuff here
     */
    
    fprintf( f, "\n" );
    fprintf( f, "fclose(f);\n" );
    fprintf( f, "\n" );

    break;
    
  case _SCILAB_ :

    fprintf( f, "\n" );
    fprintf( f, "f=mopen('%s.raw','r');\n", _BaseName( template ) );
    fprintf( f, "\n" );
    
    if ( BAL_WriteLengthDistanceScilab( f, fd,list ) != 1 ) {
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write length/distance data for scilab\n", 
		 proc );
      return( -1 );
    }

    if ( BAL_WriteStartEndTimesScilab( f, fd,list ) != 1 ) {
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write start/end times data for scilab\n", 
		 proc );
      return( -1 );
    }

    if ( BAL_WriteLengthHistogram( f, fd,list ) != 1 ) {
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write length histogram for scilab\n", 
		 proc );
      return( -1 );
    }

    if ( BAL_WriteDistanceHistogram( f, fd,list ) != 1 ) {
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write distance histogram for scilab\n", 
		 proc );
      return( -1 );
    }

    if ( BAL_WriteLifetimeHistogram( f, fd,list ) != 1 ) {
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write distance histogram for scilab\n", 
		 proc );
      return( -1 );
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


  return( 1 );
}





int BAL_PrintAllStatsScilab( char *filename, bal_blTrackList *list )
{
  char *proc = "BAL_PrintAllStatsScilab";
  if ( _PrintAllStatsXxxlab( filename, list, _SCILAB_ ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when computing/writing statistics\n", proc );
    return( -1 );
  }
  return( 1 );
}









/************************************************************
 *
 * Track related structures and I/O stuff
 *
 ************************************************************/


void BAL_InitBlTrack( bal_blTrack *t )
{
  BAL_InitBlDetectionList( &(t->detectionList) );
  t->index = 0;
  t->lcorner.x = -1.0;
  t->lcorner.y = -1.0;
  t->lcorner.z = -1.0;
  t->rcorner.x = -1.0;
  t->rcorner.y = -1.0;
  t->rcorner.z = -1.0;
  BAL_InitBlTrackProperties( &(t->properties) );
}



void BAL_FreeBlTrack( bal_blTrack *t )
{
  BAL_FreeBlDetectionList( &(t->detectionList) );
  BAL_InitBlTrack( t );
}



int BAL_CopyBlTrack( bal_blTrack *thet,  bal_blTrack *rest )
{
  char *proc = "BAL_CopyBlTrack";
  int i;

  if ( _debug_ ) {
    fprintf( stderr, "%s: will copy track of %d elements\n", proc, thet->detectionList.n );
  }

  *rest = *thet;

  BAL_InitBlDetectionList( &(rest->detectionList) );
  if ( BAL_AllocBlDetectionList( &(rest->detectionList), thet->detectionList.n ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result detection list\n", proc );
    return( -1 );
  }

  for ( i=0; i<thet->detectionList.n; i++ ) {
    rest->detectionList.data[i] = thet->detectionList.data[i];
  }

  rest->detectionList.n = thet->detectionList.n;

  rest->detectionList.unit = thet->detectionList.unit;
  rest->detectionList.lcorner = thet->detectionList.lcorner;
  rest->detectionList.rcorner = thet->detectionList.rcorner;

  return( 1 );
}





 /************************************************************
 *
 * Track list related stuff
 *
 ************************************************************/


void BAL_InitBlTrackList( bal_blTrackList *l )
{
  l->data = (bal_blTrack *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}



int BAL_AllocBlTrackList( bal_blTrackList *l, int n )
{
  char *proc = "BAL_AllocBlTrackList";
  int i;
  l->data = (bal_blTrack *)malloc( n * sizeof(bal_blTrack) );
  if ( l->data == (bal_blTrack *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  l->n_allocated = n;
  for ( i=0; i<n; i++ )
    BAL_InitBlTrack( &(l->data[i]) );
  return( 0 );
}



void BAL_FreeBlTrackList( bal_blTrackList *l )
{
  int i;

  if ( l->data != (bal_blTrack *)NULL ) {
    for ( i=0; i<l->n_allocated; i++ ) {
      BAL_FreeBlTrack( &(l->data[i]) );
    }
    free( l->data );
  }
  BAL_InitBlTrackList( l );
}





int BAL_AddBlTrackToList( bal_blTrack *d, bal_blTrackList *l )
{
  char *proc = "BAL_AddBlTrackToList";
  bal_blTrack *data;
  int i, s = l->n_allocated;

  if ( l->n == l->n_allocated ) {
    s += _size_to_be_allocated_;
    data = (bal_blTrack *)malloc( s * sizeof(bal_blTrack) );
    if ( data == (bal_blTrack *)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    if ( l->n_allocated > 0 && l->data != (bal_blTrack *)NULL ) {
      (void)memcpy( data, l->data, l->n_allocated*sizeof(bal_blTrack) );
      free( l->data );
    }
    l->n_allocated = s;
    l->data = data;
    
    for ( i=l->n; i<l->n_allocated; i ++ )
      BAL_InitBlTrack( &(l->data[i]) );

  }
  
  l->data[l->n] = *d;
  l->n ++;

  return( 0 );

}





int BAL_CopyBlTrackToList( bal_blTrack *d, bal_blTrackList *l )
{
  char *proc = "BAL_CopyBlTrackToList";
  bal_blTrack *data;
  int i, s = l->n_allocated;

  if ( l->n == l->n_allocated ) {
    s += _size_to_be_allocated_;
    data = (bal_blTrack *)malloc( s * sizeof(bal_blTrack) );
    if ( data == (bal_blTrack *)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation failed\n", proc );
      return( -1 );
    }
    if ( l->n_allocated > 0 && l->data != (bal_blTrack *)NULL ) {
      (void)memcpy( data, l->data, l->n_allocated*sizeof(bal_blTrack) );
      free( l->data );
    }
    l->n_allocated = s;
    l->data = data;
    
    for ( i=l->n; i<l->n_allocated; i ++ )
      BAL_InitBlTrack( &(l->data[i]) );

  }
  

  if ( BAL_CopyBlTrack( d, &(l->data[l->n]) ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when copying track\n", proc );
    return( -1 );
  }

  l->n ++;

  return( 0 );

}






int BAL_ReadBlTrackList( bal_blTrackList *l, char *name )
{
  char *proc = "BAL_ReadBlTrackList";
  FILE *f;
  bal_blTrack t;
  bal_blDetection d;
  int ntrack=0, itrack, nread, ndetection = 0;
  char line[_LINE_LENGTH_];
  int endisreached = 0;
  int i;

  f = fopen( name, "r" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening %s\n", proc, name );
    return( -1 );
  }

  BAL_InitBlTrack( &t );

  
  if ( fgets( line, _LINE_LENGTH_, f ) == (char*)NULL ) {
    fclose( f );
    if ( _verbose_ )
      fprintf( stderr, "%s: empty file '%s'?\n", proc, name );
    return( -1 );
  }

  
  do {
    
    BAL_InitBlTrack( &t );

    /* get a track, read the track number
     */
    if ( sscanf( line, "track %d", &(itrack) ) != 1 ) {
      BAL_FreeBlTrackList( l );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading track number in '%s'n", proc, line );
      break;
    }
    t.index = itrack;

    /* read track
     */
    ndetection = 0;
    do {

      /* end of file 
       */
      if ( fgets( line, _LINE_LENGTH_, f ) == (char*)NULL ) {
	endisreached = 1;
	break;
      }

      /* next track
       */
      if ( sscanf( line, "track %d", &(itrack) ) == 1 ) {
	break;
      }

      BAL_InitBlDetection( &d );
      nread = sscanf( line ,"%lf\t%lf\t%lf\t%d\t%d\%lf\t%d",
		      &(d.voxelcenter.x), &(d.voxelcenter.y), &(d.voxelcenter.z), 
		      &(d.voxelhalfradius1), &(d.voxelhalfradius2), 
		      &(d.arg1), &(d.imageindex) );
      
      switch ( nread ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: found '%d' args when reading line #%d of track #%d\n", 
		   proc, nread, ndetection, itrack );
	break;
      case -1 : /* end of file */
	break;
      case 4 :
	/* nikita style
	 */
	d.imageindex = d.voxelhalfradius1;
	d.voxelhalfradius1 = 0;
      case 7 :
	if ( _addBlDetectionToList( &d, &(t.detectionList) ) != 0 ) {
	  BAL_FreeBlTrack( &t );
	  BAL_FreeBlTrackList( l );
	  fclose( f );
	  if ( _verbose_ )
	    fprintf( stderr, "%s: error when adding detection #%d of track #%d\n", 
		     proc, ndetection, itrack );
	  return( -1 );
	}
	ndetection ++;
	break;
      }

    } while ( endisreached != 1 );
    
    /* a track has been read
     */
    if ( BAL_AddBlTrackToList( &t, l ) != 0 ) {
      BAL_FreeBlTrack( &t );
      BAL_FreeBlTrackList( l );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when adding track #%d\n", 
		 proc, itrack );
      return( -1 );
    }
    ntrack++;
  } while ( endisreached != 1 );

  fclose( f );


  for ( i=0; i<l->n; i++ ) {
    BAL_BlTrackBoundingBox( &(l->data[i]) );
  }

  if ( _verbose_ >= 2 ) 
    fprintf( stderr, "%s: has read %d tracks in '%s'\n", proc, ntrack, name );

  return( 0 );
}



int BAL_WriteBlTrackList( bal_blTrackList *l, char *name )
{
  char *proc = "BAL_WriteBlTrackList";
  FILE *f;
  int i, j;

  f = fopen( name, "w" );
  if ( f == NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when opening %s\n", proc, name);
    return( -1 );
  }
  
  for (i=0; i<l->n; i++ ) {
    if ( fprintf( f, "track %d\r\n", i ) < 0 ) {
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when writing heading of track #%d\n", proc, i );
      return( -1 );
    }
    for ( j=0; j<l->data[i].detectionList.n; j++ ) {
      if ( fprintf( f, "%g\t%g\t%g\t%d\t%d\t%g\t%d\r\n",
		    l->data[i].detectionList.data[j].voxelcenter.x,
		    l->data[i].detectionList.data[j].voxelcenter.y,
		    l->data[i].detectionList.data[j].voxelcenter.z,
		    l->data[i].detectionList.data[j].voxelhalfradius1,
		    l->data[i].detectionList.data[j].voxelhalfradius2,
		    l->data[i].detectionList.data[j].arg1,
		    l->data[i].detectionList.data[j].imageindex ) < 0 ) {
	fclose( f );
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when writing detection #%d of track #%d\n", proc, j, i );
	return( -1 );
      }
    }
  }
  
  fclose( f );

  return( 0 );
}






/************************************************************
 *
 * Track related procedures
 *
 ************************************************************/


void BAL_BlTrackBoundingBox( bal_blTrack *t ) 
{
  int i;
  int maxhalfradius;

  if ( t->detectionList.n < 0 ) return;
  
  maxhalfradius = t->detectionList.data[0].voxelhalfradius1;
  if ( maxhalfradius < t->detectionList.data[0].voxelhalfradius2 )
    maxhalfradius = t->detectionList.data[0].voxelhalfradius2;
  
  t->lcorner.x = t->detectionList.data[0].voxelcenter.x - maxhalfradius;
  t->lcorner.y = t->detectionList.data[0].voxelcenter.y - maxhalfradius;
  t->lcorner.z = t->detectionList.data[0].voxelcenter.z;

  t->rcorner.x = t->detectionList.data[0].voxelcenter.x + maxhalfradius;
  t->rcorner.y = t->detectionList.data[0].voxelcenter.y + maxhalfradius;
  t->rcorner.z = t->detectionList.data[0].voxelcenter.z;

  for ( i=1; i<t->detectionList.n; i++ ) {
    maxhalfradius = t->detectionList.data[i].voxelhalfradius1;
    if ( maxhalfradius < t->detectionList.data[i].voxelhalfradius2 )
      maxhalfradius = t->detectionList.data[i].voxelhalfradius2;

    if ( t->lcorner.x > t->detectionList.data[i].voxelcenter.x - maxhalfradius ) 
      t->lcorner.x = t->detectionList.data[i].voxelcenter.x - maxhalfradius;
    if ( t->lcorner.y > t->detectionList.data[i].voxelcenter.y - maxhalfradius ) 
      t->lcorner.y = t->detectionList.data[i].voxelcenter.y - maxhalfradius;
    if ( t->lcorner.z > t->detectionList.data[i].voxelcenter.z ) 
      t->lcorner.z = t->detectionList.data[i].voxelcenter.z;
    if ( t->rcorner.x < t->detectionList.data[i].voxelcenter.x + maxhalfradius ) 
      t->rcorner.x = t->detectionList.data[i].voxelcenter.x + maxhalfradius;
    if ( t->rcorner.y < t->detectionList.data[i].voxelcenter.y + maxhalfradius ) 
      t->rcorner.y = t->detectionList.data[i].voxelcenter.y + maxhalfradius;
    if ( t->rcorner.z < t->detectionList.data[i].voxelcenter.z ) 
      t->rcorner.z = t->detectionList.data[i].voxelcenter.z;
  }
}


