/****************************************************
 * thickening.c - 
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 20 jui 2014 07:27:29 CEST
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <string.h>


#include <t04t08.h>
#include <t06t26.h>

#include <thickening.h>



static int _verbose_ = 1;

void setVerboseInThickening( int v )
{
  _verbose_ = v;
}

void incrementVerboseInThickening(  )
{
  _verbose_ ++;
}

void decrementVerboseInThickening(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}








/**************************************************
 *
 *
 *
 **************************************************/


void initTypeThickeningParameters( typeThickeningParameters *p )
{
  p->maxIteration = -1;
  p->connectivity = 26;
  p->theMask = (typeChamferMask*)NULL;
  p->additionalSorting = _NO_SORTING_;
}







/**************************************************
 *
 * offset
 *
 **************************************************/



typedef struct {
  int dx;
  int dy;
  int dz;
  int di;
} typeOffset;



typedef struct {
  int nneighbors;
  typeOffset neighbors[27];
} typeNeighborhood;



static void _defineNeighborsTobBeAdded( typeNeighborhood *n, 
					int *theDim, 
					int connectivity ) 
{
  int i = 0;

  for ( i=0; i<27; i++ ) {
    n->neighbors[i].dx = 0;
    n->neighbors[i].dy = 0;
    n->neighbors[i].dz = 0;
    n->neighbors[i].di = 0;
  }

  i = 0;

  /* 04-neighbors
   */
  n->neighbors[i].dy = -1;   i++;
  n->neighbors[i].dx = -1;   i++;
  n->neighbors[i].dx =  1;   i++;
  n->neighbors[i].dy =  1;   i++;

  switch ( connectivity ) {
  case 4 :
    break;
  case 6 :
    n->neighbors[i].dz = -1;   i++;
    n->neighbors[i].dz =  1;   i++;
    break;
  case 8 :
  case 18 :
  case 26 :
  default : 
    n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;
    break;
  }
  
  switch ( connectivity ) {
  case 4 :
  case 8 :
  case 6 :
    break;
  case 26 :
  default : 
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;
  case 18 :
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   i++;
  }

  n->nneighbors = i;
  
  for ( i=0; i<n->nneighbors; i++ ) {
    n->neighbors[i].di = n->neighbors[i].dz * theDim[0] * theDim[1]
      + n->neighbors[i].dy * theDim[0]
      + n->neighbors[i].dx;
  }
}



static void _defineNeighborsForSimplicity( typeNeighborhood *n, 
					   int *theDim, 
					   int connectivity ) 
{
  int i = 0;

  for ( i=0; i<27; i++ ) {
    n->neighbors[i].dx = 0;
    n->neighbors[i].dy = 0;
    n->neighbors[i].dz = 0;
    n->neighbors[i].di = 0;
  }

  i = 0;

  if ( connectivity == 6 || connectivity == 18 || connectivity == 26 ) {
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   i++;
    n->neighbors[i].dz = -1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;
  }

  n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
  n->neighbors[i].dy = -1;   i++;
  n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
  n->neighbors[i].dx = -1;   i++;
  i++;
  n->neighbors[i].dx =  1;   i++;
  n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
  n->neighbors[i].dy =  1;   i++;
  n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;

  if ( connectivity == 6 || connectivity == 18 || connectivity == 26 ) {
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy = -1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dx =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx = -1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   i++;
    n->neighbors[i].dz =  1;   n->neighbors[i].dy =  1;   n->neighbors[i].dx =  1;   i++;
  }

  n->nneighbors = i;
  
  for ( i=0; i<n->nneighbors; i++ ) {
    n->neighbors[i].di += n->neighbors[i].dz * theDim[0] * theDim[1]
      + n->neighbors[i].dy * theDim[0]
      + n->neighbors[i].dx;
  }
}






/**************************************************
 *
 * point management
 *
 **************************************************/



typedef enum {
  _BACKGROUND_ = 0,
  _INQUEUE_ = 150,      /* in the list, but not queued for addition */
  _WILLBEADDED_ = 200,  /* in the list, and queued for addition */
  _ADDED_ = 255         
} enumTypePoint;



typedef struct {
  int x;
  int y;
  int z;
  int i;
  
  int iteration;
  int distance;

  char isinside;
  enumTypePoint type;
} typePoint;



typedef struct {
  typePoint *data;
  int n_data;
  int n_allocated_data;
} pointList;



static void initPoint( typePoint *p )
{
  p->x = -1;
  p->y = -1;
  p->z = -1;
  p->i = -1;

  p->iteration = -1;
  p->distance = -1;

  p->isinside = -1;
  p->type = _BACKGROUND_;
}



static void initPointList( pointList *l )
{
  l->data = (typePoint*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
}



static void freePointList( pointList *l )
{
  if ( l->data != (typePoint*)NULL )
    free( l->data );
  initPointList( l );
}



static int _size_to_be_allocated_ = 100;

static int addPointToPointList( pointList *l, typePoint *p )
{
  char *proc = "addPointToPointList";
  int s =  l->n_allocated_data;
  typePoint *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (typePoint*)malloc( s * sizeof(typePoint) );
    if ( data == (typePoint*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(typePoint) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }

  initPoint( &(l->data[l->n_data]) );
  l->data[l->n_data] = *p;
  l->n_data ++;
  return( 1 );
}











static void sortPointListWrtIteration( pointList *l,
				       int left, 
				       int right )
{
  int i, last;
  typePoint tmp;
  
  if ( left >= right ) return;

  tmp = l->data[left];   l->data[left] = l->data[(left+right)/2];   l->data[(left+right)/2] = tmp;
  
  last = left;
  for ( i = left+1; i <= right; i++ )       
    if ( l->data[i].iteration < l->data[left].iteration ) {
      tmp = l->data[++last];   l->data[last] = l->data[i];   l->data[i] = tmp;
    }

  tmp = l->data[left];   l->data[left] = l->data[last];   l->data[last] = tmp;
  
  sortPointListWrtIteration( l, left, last-1 );
  sortPointListWrtIteration( l, last+1, right );
}



static void sortPointListWrtDistance( pointList *l,
				       int left, 
				       int right )
{
  int i, last;
  typePoint tmp;
  
  if ( left >= right ) return;

  tmp = l->data[left];   l->data[left] = l->data[(left+right)/2];   l->data[(left+right)/2] = tmp;
  
  last = left;
  for ( i = left+1; i <= right; i++ )       
    if ( l->data[i].distance < l->data[left].distance ) {
      tmp = l->data[++last];   l->data[last] = l->data[i];   l->data[i] = tmp;
    }

  tmp = l->data[left];   l->data[left] = l->data[last];   l->data[last] = tmp;
  
  sortPointListWrtIteration( l, left, last-1 );
  sortPointListWrtIteration( l, last+1, right );
}



static void sortPointList( pointList *l, enumTypeSort sortingCriterium )
{
  switch( sortingCriterium ) {
  default :
    break;
  case _ITERATION_SORTING_ :
    sortPointListWrtIteration( l, 0, l->n_data-1 );
    break;
  case _DISTANCE_SORTING_ :
    sortPointListWrtDistance( l, 0, l->n_data-1 );
    break;
  }
}




/**************************************************
 *
 * thickening tools
 *
 **************************************************/

#ifdef _UNUSED_
static void _print2DNeighborhood( FILE *f, int *n )
{
  fprintf( stderr, "%3d %3d %3d\n", n[0], n[1], n[2] );
  fprintf( stderr, "%3d %3d %3d\n", n[3], n[4], n[5] );
  fprintf( stderr, "%3d %3d %3d\n", n[6], n[7], n[8] );
}
#endif


typedef int (*typeIsInsideFunction)( int, int, int, int * );



static int _is2DInside( int x, int y, int z, int *theDim )
{
  if ( 0 < x && x < theDim[0]-1 &&
       0 < y && y < theDim[1]-1 )
    return( 1 );
  return( 0 );
}



static int _is3DInside( int x, int y, int z, int *theDim )
{
  if ( 0 < x && x < theDim[0]-1 &&
       0 < y && y < theDim[1]-1 &&
       0 < z && z < theDim[2]-1 )
    return( 1 );
  return( 0 );
}









static int _addNeighbors( unsigned char *resBuf,
			  unsigned short *theDistance,
			  int *theDim,
			  pointList *listOfPoints,
			  int maxPossibleDistance,
			  typeNeighborhood *neighbors,
			  typeIsInsideFunction _isInside,
			  int iteration )
{
  char *proc = "_addNeighbors";
  int d, p, n;
  typePoint point;
  int x, y, z, i;
  int xn, yn, zn, in;
  int naddedpoints = 0;

  for ( d = maxPossibleDistance; d >= 1; d-- ) {
    for ( p = 0; p < listOfPoints[d].n_data; p++ ) {

      if ( listOfPoints[d].data[p].type != _WILLBEADDED_ ) continue;
      
      /* the point has just been added
	 get its neighbors and put them in queue
       */
      x = listOfPoints[d].data[p].x;
      y = listOfPoints[d].data[p].y;
      z = listOfPoints[d].data[p].z;
      i = listOfPoints[d].data[p].i;
      
      if ( listOfPoints[d].data[p].isinside ) {
	
	for ( n = 0; n < neighbors->nneighbors; n++ ) {
	  in = i + neighbors->neighbors[n].di;
	  if ( theDistance[in] > 0 && resBuf[in] == 0 ) {
	    initPoint( &point );
	    point.x = x + neighbors->neighbors[n].dx;
	    point.y = y + neighbors->neighbors[n].dy;
	    point.z = z + neighbors->neighbors[n].dz;
	    point.i = in;
	    point.iteration = iteration;
	    point.distance = 0;
	    point.isinside = (*_isInside)( point.x, point.y, point.z, theDim );
	    resBuf[in] =  point.type = _INQUEUE_;
	    if ( addPointToPointList( &(listOfPoints[theDistance[in]]), &point ) != 1 ) {
	      if ( _verbose_ )
		fprintf( stderr, "%s: unable to add point to list\n", proc );
	      return( -1 );
	    }
	    naddedpoints ++;
	  }
	}
	
      }
      else {
	
	for ( n = 0; n < neighbors->nneighbors; n++ ) {
	  
	  xn = x + neighbors->neighbors[n].dx;
	  if ( xn < 0 || xn >= theDim[0] ) continue;
	  yn = y + neighbors->neighbors[n].dy;
	  if ( yn < 0 || yn >= theDim[1] ) continue;
	  zn = z + neighbors->neighbors[n].dz;
	  if ( zn < 0 || zn >= theDim[2] ) continue;
	  
	  in = i + neighbors->neighbors[n].di;
	  if ( theDistance[in] > 0 && resBuf[in] == 0 ) {
	    initPoint( &point );
	    point.x = xn;
	    point.y = yn;
	    point.z = zn;
	    point.i = in;
	    point.iteration = iteration;
	    point.distance = 0;
	    point.isinside = (*_isInside)( point.x, point.y, point.z, theDim );
	    resBuf[in] =  point.type = _INQUEUE_;
	    if ( addPointToPointList( &(listOfPoints[theDistance[in]]), &point ) != 1 ) {
	      if ( _verbose_ )
		fprintf( stderr, "%s: unable to add point to list\n", proc );
	      return( -1 );
	    }
	    naddedpoints ++;
	  }
	}
	
      }
      
      /* neighbors have been added
	 change the point type
	 remove it from the list
       */
      resBuf[ listOfPoints[d].data[p].i ] = listOfPoints[d].data[p].type = _ADDED_;

      point = listOfPoints[d].data[ listOfPoints[d].n_data-1 ];
      listOfPoints[d].data[ listOfPoints[d].n_data-1 ] = listOfPoints[d].data[p];
      listOfPoints[d].data[p] = point;
      listOfPoints[d].n_data --;
      p --;
    }
  }

  return( naddedpoints );
}





static void _extractNeighborhood( typePoint *pt,
				 unsigned char *resBuf,
				 int *theDim,
				 int *neighb,
				 typeNeighborhood *neighbors )
{
  int n;

  if ( pt->isinside ) {
    
    for ( n=0; n<neighbors->nneighbors; n++ ) {
      neighb[n] = resBuf[ pt->i + neighbors->neighbors[n].di ];
    }

  }
  else {
    
    for ( n=0; n<neighbors->nneighbors; n++ ) {
      if ( pt->x+neighbors->neighbors[n].dx < 0 || pt->x+neighbors->neighbors[n].dx >= theDim[0] ||
	   pt->y+neighbors->neighbors[n].dy < 0 || pt->y+neighbors->neighbors[n].dy >= theDim[1] ||
	   pt->z+neighbors->neighbors[n].dz < 0 || pt->z+neighbors->neighbors[n].dz >= theDim[2] )
	neighb[n] = 0;
      else
	neighb[n] = resBuf[ pt->i + neighbors->neighbors[n].di ];
    }

  }
}





typedef int (*typeIsPointSimple)( int *,
				  int *,
				  int * );



static int _is2DPoint08Simple( int *neighb, int *t04, int *t08 )
{
  int checkT04, checkT08;
  int n;
  
  for ( n=0; n<9; n++ )
    if ( neighb[n] == _INQUEUE_ )
      neighb[n] = _BACKGROUND_;

  Compute_T04_and_T08( neighb, t04, t08 );

  if ( *t04 != 1 || *t08 != 1 ) return( 0 );
  
  for ( n=0; n<9; n++ )
    if ( neighb[n] == _WILLBEADDED_ )
      neighb[n] = _BACKGROUND_;
  
  Compute_T04_and_T08( neighb, &checkT04, &checkT08 );

  if ( checkT04 != 1 || checkT08 != 1 ) return( 0 );

  return( 1 );
}



static int _is3DPoint26Simple( int *neighb, int *t06, int *t26 )
{
  int checkT06, checkT26;
  int n;
  
  for ( n=0; n<27; n++ )
    if ( neighb[n] == _INQUEUE_ )
      neighb[n] = _BACKGROUND_;

  Compute_T06_and_T26( neighb, t06, t26 );
	  
  if ( *t06 != 1 || *t26 != 1 ) return( 0 );

  for ( n=0; n<27; n++ )
    if ( neighb[n] == _WILLBEADDED_ )
      neighb[n] = _BACKGROUND_;
  
  Compute_T06_and_T26( neighb, &checkT06, &checkT26 );

  if ( checkT06 != 1 || checkT26 != 1 ) return( 0 );

  return( 1 );
}








/**************************************************
 *
 * thickening initialization
 *
 **************************************************/


int InitializeThickeningImage( unsigned char *resBuf,
			       unsigned short *theDistance,
			       int *theDim )
{
  int i, max;
  int v = theDim[0]*theDim[1]*theDim[2];
  

  max = theDistance[0];
  for ( i=1; i<v; i++ )
    if ( max < theDistance[i] )
      max = theDistance[i];

  for ( i=0; i<v; i++ )
    resBuf[i] = 0;

  for ( i=0; i<v; i++ )
    if ( max == theDistance[i] ) {
      resBuf[i] = 255;
      break;
    }
  
  return( 1 );
}
			       











/**************************************************
 *
 * thickening
 *
 **************************************************/


int ThickeningImage( unsigned char *resBuf,
		     unsigned short *theDistance,
		     unsigned short *thePropagation,
		     int *theDim,
		     typeThickeningParameters *par )
{
  char *proc = "ThickeningImage";

  int connectivity = par->connectivity;
  
  int maxPossibleDistance = 65535;
  pointList *listOfPoints = (pointList *)NULL;
  int x, y, z;
  int d, i, n, p;
  int v = theDim[0]*theDim[1]*theDim[2];

  typeNeighborhood offsetToAddPoints;
  typeNeighborhood offsetForSimplicity;
  int neighb[27];
  int tback, tfore;

  typePoint point;

  int iteration = 0;
  int nInitialPoints = 0;
  int nAddedPointsQueue, nAddedPointsObjects;

  typeIsInsideFunction _isInside = (typeIsInsideFunction)NULL;
  typeIsPointSimple _isPointSimple = (typeIsPointSimple)NULL;

  enumTypeSort sortingCriterium = par->additionalSorting;
  typePoint firstPoint;
  int apointhasbeenadded, stopparsingpoints;




  /*--------------------------------------------------
   *
   * 
   *
   --------------------------------------------------*/
  
  if ( sortingCriterium == _DISTANCE_SORTING_ ) {
    if ( par->theMask == (typeChamferMask *)NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: no chamfer mask, switch to no additional sorting\n", proc );
      }
      par->additionalSorting = sortingCriterium = _NO_SORTING_;
    }
    else {
      for ( n=0; n<par->theMask->nb; n++ ) {
	par->theMask->list[n].o = par->theMask->list[n].z * theDim[0]*theDim[1]
	  + par->theMask->list[n].y * theDim[0] 
	  + par->theMask->list[n].x;
      }
    }
  }

  if ( _verbose_ >= 3 ) { 
    fprintf( stderr, "%s: sorting criterium = ", proc );
    switch ( sortingCriterium ) {
    default :
      fprintf( stderr, "default \n" ); break;
    case _NO_SORTING_ :
      fprintf( stderr, "none \n" ); break;
    case _ITERATION_SORTING_ :
      fprintf( stderr, "iteration \n" ); break;
    case _DISTANCE_SORTING_ :
      fprintf( stderr, "distance \n" ); break;
    }
  }


  /*--------------------------------------------------
   *
   * connectivity based choices
   *
   --------------------------------------------------*/

  switch( connectivity ) {
  case 4 :
  case 8 :
  case 6 :
  case 18 :
  case 26 :
    break;
  default :
    connectivity = 26;
  }
    
  if ( theDim[2] == 1 ) {
    switch( connectivity ) {
    case 6 :
      connectivity = 4; break;
    case 18 :
    case 26 :
      connectivity = 8; break;
    default :
      break;
    }
  }

  if ( _verbose_ >= 2 )
    fprintf( stderr, "%s: connectivity set to %d\n", proc, connectivity );
    
  switch( connectivity ) {
  default :
    if ( _verbose_ ) 
      fprintf( stderr, "%s: weird situation\n", proc );
    return( -1 );
  case 4 :
    _isInside = &_is2DInside;
    _isPointSimple = &_is2DPoint08Simple;
    break;
  case 8 :
    _isInside = &_is2DInside;
    _isPointSimple = &_is2DPoint08Simple;
    break;
  case 6 :
    _isInside = &_is3DInside;
    _isPointSimple = &_is3DPoint26Simple;
    break;
  case 18 :
    _isInside = &_is3DInside;
    _isPointSimple = &_is3DPoint26Simple;
    break;
  case 26 :
    _isInside = &_is3DInside;
    _isPointSimple = &_is3DPoint26Simple;
    break;
  }


  /*--------------------------------------------------
   *
   * pre-compute offsets for point access speed-up
   *
   --------------------------------------------------*/
 
  _defineNeighborsTobBeAdded( &offsetToAddPoints, theDim, connectivity );
  _defineNeighborsForSimplicity( &offsetForSimplicity, theDim, connectivity );




  /*--------------------------------------------------
   *
   * build point list
   *
   --------------------------------------------------*/
  
  listOfPoints = (pointList*)malloc( (maxPossibleDistance+1)*sizeof(pointList) );
  if ( listOfPoints == (pointList*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate point list array\n", proc );
    return( -1 );
  }
  for ( d=0; d<=maxPossibleDistance; d++ ) 
    initPointList( &(listOfPoints[d]) );

  
  /* looking for start points
     since the value is unknown, be sure it is correctly set
  */
  for ( nInitialPoints=0, i=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ )
  for ( x=0; x<theDim[0]; x++, i++ ) {
    if ( resBuf[i] > 0 ) {
      initPoint( &point );
      point.x = x;
      point.y = y;
      point.z = z;
      point.i = i;
      point.iteration = iteration;
      point.distance = 0;
      point.isinside = (*_isInside)( point.x, point.y, point.z, theDim );
      resBuf[i] =  point.type = _WILLBEADDED_;
      if ( addPointToPointList( &(listOfPoints[theDistance[i]]), &point ) != 1 ) {
	for ( d=0; d<=maxPossibleDistance; d++ ) 
	  freePointList( &(listOfPoints[d]) );
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to add point (%d,%d,%d) to list\n", proc, x, y, z );
	return( -1 );
      }
      nInitialPoints ++;
    }
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... found %d starting points\n", nInitialPoints );
  }

  

  

  /*--------------------------------------------------
   *
   * initialization of propagation image
   *
   --------------------------------------------------*/

  if ( thePropagation != (unsigned short int*) NULL ) {
    for ( i=0; i<v; i++ ) 
      thePropagation[i] = 0;
  }





  /*--------------------------------------------------
   *
   * thickening
   *
   --------------------------------------------------*/

  do {

    iteration ++;
    
    /* update the propgation image with iteration number
     */
    switch( sortingCriterium ) {
    default :
    case _ITERATION_SORTING_ :
      for ( d=0; d<=maxPossibleDistance; d++ ) {
	for ( p=0; p<listOfPoints[d].n_data; p++ ) {
	  if ( listOfPoints[d].data[p].type == _WILLBEADDED_ ) {
	    thePropagation[ listOfPoints[d].data[p].i ] = iteration;
	  }
	}
      }
      break;
    case _DISTANCE_SORTING_ :
      for ( d=0; d<=maxPossibleDistance; d++ ) {
	for ( p=0; p<listOfPoints[d].n_data; p++ ) {
	  if ( listOfPoints[d].data[p].type == _WILLBEADDED_ ) {
	    thePropagation[ listOfPoints[d].data[p].i ] = listOfPoints[d].data[p].distance;
	  }
	}
      }
      break;
    }


    /* turn the _WILLBEADDED_ points into _ADDED_
       remove them from the list
       add their neighboring points to the list
     */
    nAddedPointsQueue = _addNeighbors( resBuf, theDistance, theDim,
				       listOfPoints, maxPossibleDistance,
				       &offsetToAddPoints, _isInside, iteration );
    if ( nAddedPointsQueue == -1 ) {
      for ( d=0; d<=maxPossibleDistance; d++ ) 
	  freePointList( &(listOfPoints[d]) );
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to add points at iteration %d\n", proc, iteration );
	return( -1 );
    }
    
    
    /* identify points that can be added to the image
       turn them from _INQUEUE_ into _WILLBEADDED_
     */
    for ( nAddedPointsObjects=0, d=maxPossibleDistance; d>=1 && nAddedPointsObjects==0 ; d-- ) {

      if ( listOfPoints[d].n_data == 0 ) continue;

      /* propagation distance for _INQUEUE_ points  
	 (only for those that may be added)
       */
      
      switch( sortingCriterium ) {
      default :
      case _ITERATION_SORTING_ :
	break;
      case _DISTANCE_SORTING_ :
	for ( p=0; p<listOfPoints[d].n_data; p++ ) { 
	  if ( listOfPoints[d].data[p].type != _INQUEUE_ ) continue;
	  for ( n=0; n < par->theMask->nb; n++ ) {
	    x = listOfPoints[d].data[p].x +  par->theMask->list[n].x;
	    if ( x < 0 || theDim[0] <= x ) continue;
	    y = listOfPoints[d].data[p].y +  par->theMask->list[n].y;
	    if ( y < 0 || theDim[1] <= y ) continue;
	    z = listOfPoints[d].data[p].z +  par->theMask->list[n].z;
	    if ( z < 0 || theDim[2] <= z ) continue;
	    if ( resBuf[ listOfPoints[d].data[p].i + par->theMask->list[n].o ] != _ADDED_ )
	      continue;
	    if ( listOfPoints[d].data[p].distance > thePropagation[ listOfPoints[d].data[p].i + par->theMask->list[n].o ] +  par->theMask->list[n].inc )
	      listOfPoints[d].data[p].distance = thePropagation[ listOfPoints[d].data[p].i + par->theMask->list[n].o ] +  par->theMask->list[n].inc;
	  }
	}
	break;
      }


      /* sort points
       */
      if ( listOfPoints[d].n_data > 1 ) {
	sortPointList(  &(listOfPoints[d]), sortingCriterium );
      }

     if ( _verbose_ >= 3 )
	fprintf( stderr, "    parsing %d points at distance %d\n", listOfPoints[d].n_data, d );
     
     

      /* process points 
       */

      for ( apointhasbeenadded=0, stopparsingpoints=0, p=0; 
	    p<listOfPoints[d].n_data && stopparsingpoints == 0; p++ ) {
	
	switch( sortingCriterium ) {
	default :
	  break;
	case _ITERATION_SORTING_ :
	  if ( apointhasbeenadded )
	    if ( listOfPoints[d].data[p].iteration != firstPoint.iteration )
	      stopparsingpoints = 1;
	  break;
	case _DISTANCE_SORTING_ :
	  if ( apointhasbeenadded )
	    if ( listOfPoints[d].data[p].distance != firstPoint.distance )
	      stopparsingpoints = 1;
	  break;
	}



	/* get the neighborhood
	 */

	_extractNeighborhood( &(listOfPoints[d].data[p]),
			      resBuf, theDim, neighb,
			      &offsetForSimplicity );

	/* test the simplicity
	 */

	if ( (*_isPointSimple)( neighb, &tback, &tfore ) == 1 ) {
	  /* if there is an end condition for simple points,
	     it is to be tested here
	  */
	  listOfPoints[d].data[p].type = resBuf[ listOfPoints[d].data[p].i ] = _WILLBEADDED_;
	  nAddedPointsObjects ++;

	  if ( apointhasbeenadded == 0 ) {
	    apointhasbeenadded = 1;
	    firstPoint = listOfPoints[d].data[p];
	  }
	}

      }
      
      if ( _verbose_ >= 4 ) {
	fprintf( stderr, "  distance=%3d - %4d points added to object\n", d, nAddedPointsObjects );
      }

    }
    
    if ( _verbose_ >= 1 ) {
      fprintf( stderr, " #%6d", iteration );
      fprintf( stderr, " Points to queue=%6d", nAddedPointsQueue );
      fprintf( stderr, " Points to object=%6d", nAddedPointsObjects );
      fprintf( stderr, " distance=%d", d+1 );
      if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );
      else                  fprintf( stderr, "\r" );
    }

  } while ( ((par->maxIteration >= 0 && iteration < par->maxIteration) || (par->maxIteration < 0))
	    && nAddedPointsObjects > 0 );


  /* transform output image
   */
  for ( i=0; i<v; i++ ) {
    if ( resBuf[i] ==  _ADDED_ )
      resBuf[i] = 255;
    else
      resBuf[i] = 0;
  }

  return( 1 );
}
