/****************************************************
 * regionalmax.c - 
 *
 * $Id$
 *
 * Copyright (c) INRIA 2008
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Thu May 22 07:58:29 CEST 2008
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */




#include <chunks.h>
#include <connexe.h>
#include <regionalmax.h>



/*
 * static global variables
 * verbose,
 * management of ambiguous cases
 * memory management
 * 
 */



static int _verbose_ = 1;

void setVerboseInRegionalMax( int v )
{
  _verbose_ = v;
}

void incrementVerboseInRegionalMax(  )
{
  _verbose_ ++;
}

void decrementVerboseInRegionalMax(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}


static int _component_wise_processing_ = 1;

void allowComponentWiseProcessingInRegionalMax(  )
{
  _component_wise_processing_ = 1;
}

void disallowComponentWiseProcessingInRegionalMax(  )
{
  _component_wise_processing_ = 0;
}




static int _nPointsToBeAllocated_ = 1000;

void regionalmax_setNumberOfPointsForAllocation( int n )
{
   if ( n > 0 ) _nPointsToBeAllocated_ = n;
}



/*************************************************************
 *
 * static structures and functions: list management
 *
 ************************************************************/



typedef struct {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
  int i; /* x+y*dimx+z*dimx*dimy */
  char t; /* if 1, all neighbors are in the image (allows to skip tests) */
} typePoint;

typedef struct {
  int nPoints;
  int nAllocatedPoints;
  typePoint *pt;
} typePointList;



static void init_PointList( typePointList *l );
static void free_PointList( typePointList *l );
static void stats_PointList( typePointList *l, int n );




/*************************************************************
 *
 * static structures and functions: tools
 *
 ************************************************************/



static int find_maximal_value( void *theImage, 
			       bufferType theImageType, int *theDim );

static int build_first_list ( typePointList *thePointList,
			      void *theInput, void *theOutput, bufferType theType, int *theDim, 
			      int maxValue );

static int build_list ( typePointList *thePointList,
			void *theInput, void *theOutput, bufferType theType, int *theDim, 
			int maxValue );

static int process_list ( typePointList *thePointList,
			  void *theInput, void *theOutput, bufferType theType, int *theDim, 
			  int maxValue );

static int process_component ( typePointList *thePointList,
			       void *theInput, void *theOutput, bufferType theType, int *theDim, 
			       int maxValue );












/* The regional maxima are defined as the difference between the input I and 
min( I-height, I*multiplier) dilated infinitely under I.
With multiplier=1.0, we get the H-maxima
*/ 
int regionalmax( void *theInput, void *theOutput, bufferType theType,
		 int *theDim, int height, double multiplier )
{
  char *proc = "regionalmax";
  int error_value = -1;

  int v = theDim[0]*theDim[1]*theDim[2];
  void *theTmp = NULL;
  int imin;
  int i, j;

  int maxValue;
  typePointList *thePointList = NULL;

  int changes, localchanges;
  int iteration = 0;


  /* allocating the auxiliary image
   */
  if ( theInput == theOutput ) {

    switch ( theType ) {

    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: such image type not handled yet\n", proc );
      }
      return( error_value );

    case UCHAR :
      theTmp = malloc( v * sizeof( unsigned char ) );
      break;

    case USHORT :
      theTmp = malloc( v * sizeof( unsigned short int ) );
      break;
    }
    
    if ( theTmp == NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to allocate auxiliary image\n", proc );
      }
      return( error_value );
    }

  }
  else {

    theTmp = theOutput;

  }
  


  /* initializing the auxiliary image
     typically, substracting h 
   */
  switch ( theType ) {

    default :
      if ( theInput == theOutput ) free( theTmp );
      if ( _verbose_ ) {
	fprintf( stderr, "%s: such image type not handled yet\n", proc );
      }
      return( error_value );

  case UCHAR : 
    {
      unsigned char *theBuf = (unsigned char*)theInput;
      unsigned char *resBuf = (unsigned char*)theTmp;
      for ( i=0; i<v; i++ ) {
	imin = theBuf[i] - height;
	if ( imin < 0 ) imin = 0;
	else {
	  if ( imin > (int)( theBuf[i]*multiplier + 0.5 ) )
	    imin = (int)( theBuf[i]*multiplier + 0.5 );
	}
	resBuf[i] = imin;
      }
    }
    break;

  case USHORT : 
    {
      unsigned short int *theBuf = (unsigned short int*)theInput;
      unsigned short int *resBuf = (unsigned short int*)theTmp;
      for ( i=0; i<v; i++ ) {
	imin = theBuf[i] - height;
	if ( imin < 0 ) imin = 0;
	else {
	  if ( imin > (int)( theBuf[i]*multiplier + 0.5 ) )
	    imin = (int)( theBuf[i]*multiplier + 0.5 );
	}
	resBuf[i] = imin;
      }
    }
    break;

  }
    

  


  /* structure allocation and initialization
     we get the maximal value of the transformed image
     it is not necessary to allocate a list for points with this value
     since it can not change
   */
  
  maxValue = find_maximal_value( theTmp, theType, theDim );
  thePointList = (typePointList *)malloc( (maxValue)*sizeof(typePointList) );
  if ( thePointList == (typePointList*)NULL ) {
    if ( theInput == theOutput ) free( theTmp );
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate auxiliary array\n", proc );
    return( -1 );
  }
    
  for ( i=0; i<maxValue; i++ )
    init_PointList( &thePointList[i] );



  /* main loop
   */
  do {
    changes = 0;


    if ( _verbose_ >= 2 )
      fprintf( stderr, "iteration #%d: entering list building\n", iteration );

    for ( i=0; i<maxValue; i++ ) 
      thePointList[i].nPoints = 0;

    /* build list of points
       points are different from ref (typically below) and from maxValue
     */
    if ( iteration == 0 ) {
      if ( build_first_list( thePointList, theInput, theTmp, theType, theDim, maxValue ) < 0 ) {
	for ( i=0; i<maxValue; i++ )
	  free_PointList( &thePointList[i] );
	free( thePointList );
	if ( theInput == theOutput ) free( theTmp );
	
      }
    }
    else {
      if ( build_list( thePointList, theInput, theTmp, theType, theDim, maxValue ) < 0 ) {
	for ( i=0; i<maxValue; i++ )
	  free_PointList( &thePointList[i] );
	free( thePointList );
	if ( theInput == theOutput ) free( theTmp );
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: unable to build list at iteration %d\n", proc, iteration );
	return( -1 );
      }
    }
    

    if ( _verbose_ >= 2 )
      fprintf( stderr, "iteration #%d: exiting list building\n", iteration );
    

    if ( _verbose_ >= 3 ) stats_PointList( thePointList, maxValue );

    
    for ( j=maxValue-1; j>=0; j-- ) {
      
      if ( thePointList[j].nPoints == 0 ) continue;

      if ( _verbose_ >= 2 )
	fprintf( stderr, " - iteration #%d: list #%5d, n=%d\n", iteration, j, thePointList[j].nPoints );

      localchanges = process_list ( &(thePointList[j]),  theInput, theTmp, theType, theDim, maxValue );

      if ( localchanges < 0 ) {
	for ( i=0; i<maxValue; i++ ) free_PointList( &thePointList[i] );
	free( thePointList );
	if ( theInput == theOutput ) free( theTmp );
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: error when processing list %d at iteration %d\n", proc, j, iteration );
	return( -1 );
      }

      if ( _verbose_ >= 2 ) {
	fprintf( stderr, "                 n=%8d, changes = %8d\n", thePointList[j].nPoints, localchanges );
      }
      
      changes += localchanges;
      
    }
    
    if ( _verbose_ >= 1 )
      fprintf( stderr,"iteration #%2d - changes = %8d\n", iteration, changes );

    iteration ++;

  } while ( changes > 0 );






  for ( i=0; i<maxValue; i++ ) free_PointList( &thePointList[i] );
  free( thePointList );



  
  
  /* computing the ouput
   */
  switch ( theType ) {

    default :
      if ( theInput == theOutput ) free( theTmp );
      if ( _verbose_ ) {
	fprintf( stderr, "%s: such image type not handled yet\n", proc );
      }
      return( error_value );

  case UCHAR : 
    {
      unsigned char *theBuf = (unsigned char*)theInput;
      unsigned char *tmpBuf = (unsigned char*)theTmp;
      unsigned char *resBuf = (unsigned char*)theOutput;
      for ( i=0; i<v; i++ ) resBuf[i] = theBuf[i] - tmpBuf[i];
    }
    break;

  case USHORT : 
    {
      unsigned short int *theBuf = (unsigned short int*)theInput;
      unsigned short int *tmpBuf = (unsigned short int*)theTmp;
      unsigned short int *resBuf = (unsigned short int*)theOutput;
      for ( i=0; i<v; i++ ) resBuf[i] = theBuf[i] - tmpBuf[i];
    }
    break;

  }
  
  if ( theInput == theOutput ) free( theTmp );

  return( 1 );
}











/*************************************************************
 *
 * static structures and functions: tools
 *
 ************************************************************/

#ifdef _UNUSED_
static void print_Point( FILE *f, typePoint *p )
{
  fprintf( f, "[%5d %5d %5d] = %8d, interior=%1d\n",
	   p->x, p->y, p->z, p->i, p->t );
}
#endif

static void init_PointList( typePointList *l )
{
  l->nPoints = 0;
  l->nAllocatedPoints = 0;
  l->pt = NULL;
}

static void free_PointList( typePointList *l )
{
  if ( l->pt != NULL ) free ( l-> pt );
  init_PointList( l );
}

static void stats_PointList( typePointList *l, int n )
{
  int i, s=0, a=0;
  
  for ( i=n-1; i>=0; i-- ) {
    if ( l[i].nPoints > 0 )
      fprintf( stderr, " list #%3d: %8d points / %8d allocated\n", i, l[i].nPoints, l[i].nAllocatedPoints );
    s += l[i].nPoints;
    a += l[i].nAllocatedPoints;
  }
  fprintf( stderr, " total:     %8d points / %8d allocated\n", s, a );
}

static int addPointToList( typePointList *l, 
			   int x, int y, int z, int i, int t )
{
  char *proc = "addPointToList";
  int n = l->nPoints;
  int newn;
  typePoint *pt = NULL;

  if ( n == l->nAllocatedPoints ) {

    newn = l->nAllocatedPoints + _nPointsToBeAllocated_;

    pt = (typePoint *)malloc( newn * sizeof( typePoint ) );
    if ( pt == NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: can not reallocate point list\n", proc );
	fprintf( stderr, "\t failed to add (%d,%d,%d) in list", x, y, z );
      }
      return( -1 );
    }

    if ( l->pt != NULL ) {
      (void)memcpy( (void*)pt, (void*)l->pt,  l->nAllocatedPoints * sizeof( typePoint ) );
      free( l->pt );
    }
    
    l->pt = pt;
    l->nAllocatedPoints = newn;
    
  }

  l->pt[ n ].x = x;
  l->pt[ n ].y = y;
  l->pt[ n ].z = z;
  l->pt[ n ].i = i;
  l->pt[ n ].t = t;

  l->nPoints ++;

  return( l->nPoints );
}








/************************************************************
 *
 * maximal value of the image
 */

static int find_maximal_value( void *theImage, 
			       bufferType theImageType, int *theDim )
{
  char *proc = "find_maximal_value";
  int error_value = -1;
  int i;
  int v = theDim[0]*theDim[1]*theDim[2];
  int maximal_value = 0;

  switch( theImageType ) {

  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such gradient type not handled yet\n", proc );
    }
    return( error_value );

  case UCHAR :
    {
      u8 *buf = theImage;
      maximal_value = buf[0];
      for (i=1;i<v;i++)
	if ( maximal_value < buf[i] )  maximal_value = buf[i];
    }
    break;

  case USHORT :
    {
      u16 *buf = theImage;
      maximal_value = buf[0];
      for (i=1;i<v;i++)
	if ( maximal_value < buf[i] )  maximal_value = buf[i];
    }
    break;
  }

  return( maximal_value );
}





/************************************************************
 *
 * build list of point from image
 */

static int build_first_list ( typePointList *thePointList,
			      void *theInput, void *theOutput, bufferType theType, int *theDim, 
			      int maxValue )
{
  char *proc = "build_first_list";
  int error_value = -1;
  int *histo = NULL;

  int i;
  int v =  theDim[0]*theDim[1]*theDim[2];

  int x, y, z, t;
  typePoint *pt;
  int dimx = theDim[0];
  int dimx1 = theDim[0]-1;
  int dimy = theDim[1];
  int dimy1 = theDim[1]-1;
  int dimz = theDim[2];
  int dimz1 = theDim[2]-1;

  /* histogram computation
   */

  histo = (int*)malloc( (maxValue+1)*sizeof(int) );
  if ( histo == (int*)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate histogram\n", proc );
    }
    return( error_value );
  }

  for (i=0;i<=maxValue;i++) histo[i] = 0;
  
  switch ( theType ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    }
    return( error_value );

  case UCHAR :
    {
      unsigned char *theBuf = (unsigned char *)theOutput;
      
      for ( i=0; i<v; i++ ) histo[ theBuf[i] ] ++;
    }
    break;

  case USHORT :
    {
      unsigned short int *theBuf = (unsigned short int *)theOutput;
      
      for ( i=0; i<v; i++ ) histo[ theBuf[i] ] ++;
    }
    break;

  }
  


  /* list allocation
   */

  for ( i=0; i<maxValue; i++ ) {

    if ( histo[i] == 0 ) continue;

    thePointList[i].pt = (typePoint *)malloc( histo[i] * sizeof( typePoint ) );

    if ( thePointList[i].pt == (typePoint *)NULL ) {
      free( histo );
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to allocate list #%d\n", proc, i );
      }
      return( error_value );
    }

    thePointList[i].nAllocatedPoints = histo[i];
 
  }

  free( histo );

  /* fill the list
   */
  
  switch ( theType ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    }
    return( error_value );

  case UCHAR :
    {
      unsigned char *theRef = (unsigned char *)theInput;
      unsigned char *theBuf = (unsigned char *)theOutput;

      for ( i=0, z=0; z<dimz; z++ ) 
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++, i++ ) {

	if ( theBuf[i] == theRef[i] ) continue;
	if ( theBuf[i] == maxValue  ) continue;

	pt = &(thePointList[ theBuf[i] ].pt[ thePointList[ theBuf[i] ].nPoints ]);
	
	pt->x = x;
	pt->y = y;
	pt->z = z;
	pt->i = i;
	
	t = 1;
	if ( x == 0 ) t = 0;
	else if ( x == dimx1 ) t = 0;
	else if ( y == 0 ) t = 0;
	else if ( y == dimy1 ) t = 0;
	else if ( z == 0 ) t = 0;
	else if ( z == dimz1 ) t = 0;
	
	pt->t = t;
	
	thePointList[ theBuf[i] ].nPoints ++;
      }

    }
    break;

  case USHORT :
    {
      unsigned short int *theRef = (unsigned short int *)theInput;
      unsigned short int *theBuf = (unsigned short int *)theOutput;

      for ( i=0, z=0; z<dimz; z++ ) 
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++, i++ ) {

	if ( theBuf[i] == theRef[i] ) continue;
	if ( theBuf[i] == maxValue  ) continue;

	pt = &(thePointList[ theBuf[i] ].pt[ thePointList[ theBuf[i] ].nPoints ]);
	
	pt->x = x;
	pt->y = y;
	pt->z = z;
	pt->i = i;
	
	t = 1;
	if ( x == 0 ) t = 0;
	else if ( x == dimx1 ) t = 0;
	else if ( y == 0 ) t = 0;
	else if ( y == dimy1 ) t = 0;
	else if ( z == 0 ) t = 0;
	else if ( z == dimz1 ) t = 0;
	
	pt->t = t;
	
	thePointList[ theBuf[i] ].nPoints ++;
      }

    }
    break;
  }
  
  return( 1 );

}






#define _ADD_POINT( TEST ) {  \
    if ( theBuf[i] == theRef[i] ) continue; \
    if ( theBuf[i] == maxValue  ) continue; \
    if ( addPointToList( &(thePointList[theBuf[i]]), x, y, z, i, TEST ) <= 0 ) { \
      return( error_value ); \
    } \
  }


static int build_list ( typePointList *thePointList,
			void *theInput, void *theOutput, bufferType theType, int *theDim, 
			int maxValue )
{
  char *proc = "build_list";
  int error_value = -1;
  int i, x, y, z, t;
  int dimx = theDim[0];
  int dimx1 = theDim[0]-1;
  int dimy = theDim[1];
  int dimy1 = theDim[1]-1;
  int dimz = theDim[2];
  int dimz1 = theDim[2]-1;

  switch ( theType ) {

  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    }
    return( error_value );

  case UCHAR :
    {
      unsigned char *theRef = (unsigned char *)theInput;
      unsigned char *theBuf = (unsigned char *)theOutput;
      
      for ( i=0, z=0; z<dimz; z++ ) 
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++, i++ ) {
	if ( theBuf[i] == theRef[i] ) continue;
	if ( theBuf[i] == maxValue  ) continue;

	t = 1;
	if ( x == 0 ) t = 0;
	else if ( x == dimx1 ) t = 0;
	else if ( y == 0 ) t = 0;
	else if ( y == dimy1 ) t = 0;
	else if ( z == 0 ) t = 0;
	else if ( z == dimz1 ) t = 0;
	
	if ( addPointToList( &(thePointList[theBuf[i]]), x, y, z, i, t ) <= 0 ) {
	  return( error_value );
	}
      }

    }
    break;

  case USHORT :
    {
      unsigned short int *theRef = (unsigned short int *)theInput;
      unsigned short int *theBuf = (unsigned short int *)theOutput;
      

      /* gaffe,
	 faut verifier que les dimensions sont plus grandes que 1
	 pour faire comme ci-dessous
      */
      /*
      z = 0;
      for ( i=0, y=0; y<dimy; y++ ) 
      for ( x=0; x<dimx; x++, i++ ) _ADD_POINT( 0 )

      for ( z=1; z<dimz1; z ++ ) {

	y = 0;
	for ( x=0; x<dimx; x++, i++ ) _ADD_POINT( 0 )
	for ( y=1; y<dimy1; y++ ) {
	  x = 0; _ADD_POINT( 0 ) i++;
	  for ( x=1; x<dimx1; x++, i++ ) _ADD_POINT( 1 )
	  x = dimx1; _ADD_POINT( 0 ) i++;
	}
	y = dimy1;
	for ( x=0; x<dimx; x++, i++ ) _ADD_POINT( 0 )
      }
      z = dimz1;
      for ( y=0; y<dimy; y++ ) 
      for ( x=0; x<dimx; x++, i++ ) _ADD_POINT( 0 )
      */

      for ( i=0, z=0; z<dimz; z++ ) 
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++, i++ ) {
	if ( theBuf[i] == theRef[i] ) continue;
	if ( theBuf[i] == maxValue  ) continue;

	t = 1;
	if ( x == 0 ) t = 0;
	else if ( x == dimx1 ) t = 0;
	else if ( y == 0 ) t = 0;
	else if ( y == dimy1 ) t = 0;
	else if ( z == 0 ) t = 0;
	else if ( z == dimz1 ) t = 0;
	
	if ( addPointToList( &(thePointList[theBuf[i]]), x, y, z, i, t ) <= 0 ) {
	  return( error_value );
	}
      }

    }
    break;
  }
  
  return( 1 );
}





typedef struct {
  typePointList *thePointList;
  void *theInput;
  void *theOutput; 
  int *theDim; 
  
  int changes;
} _processListParam;



static int _processListSubroutine_u8( void *parameter,
				      size_t first,
				      size_t last )
{
  typePointList *thePointList = ((_processListParam*)parameter)->thePointList;
  unsigned char *theRef = (unsigned char *)(((_processListParam*)parameter)->theInput);
  unsigned char *theBuf = (unsigned char *)(((_processListParam*)parameter)->theOutput);
  int dimx = ((_processListParam*)parameter)->theDim[0];
  int dimy = ((_processListParam*)parameter)->theDim[1];
  int dimz = ((_processListParam*)parameter)->theDim[2];
  int dimxy = dimx*dimy;
  
  typePoint *p;
  int x, y, z, i;
  int n, val;
  int j, k, l;
  int c = 0;

  for ( n=first; n<=last; n ++ ) {
    
    /* get the point
     */
    p = &(thePointList->pt[n]);
    x = p->x;
    y = p->y;
    z = p->z;
    i = p->i;
    val = theBuf[i];
	  
    /* dilate (26-connectivity)
       warning: in parallelism, buffer value are accessed outside the allocated range
       of each chunk -> we may have a different behavior than the sequential
       computation (but should be the same result)
       can be corrected if two buffers are used (one as input, the second one as output)
       warning: this may cause access to a data being written by an other chunk :(
     */
    
    if ( p->t ) {
      for ( l= -1; l<=1; l++ )
	for ( k= -1; k<=1; k++ )
	  for ( j= -1; j<=1; j++ ) {
	    if ( val < theBuf[i+l*dimxy+k*dimx+j] ) val = theBuf[i+l*dimxy+k*dimx+j];
	  }
    } 
    else {
      for ( l= -1; l<=1; l++ ) {
	if ( (z+l >= 0) && (z+l < dimz) ) {
	  for ( k= -1; k<=1; k++ ) {
	    if ( (y+k >= 0) && (y+k < dimy) ) {
	      for ( j= -1; j<=1; j++ ) {
		if ( (x+j >= 0) && (x+j < dimx) ) {
		  if ( val < theBuf[i+l*dimxy+k*dimx+j] ) val = theBuf[i+l*dimxy+k*dimx+j];
		}
	      }
	    }
	  }
	}
      }
    }
	  
    
    /* val is the value of the input image after dilation
       check the correspondant value of the reference image
    */
    if ( val > theRef[i] ) val = theRef[i];

    /* there is a change if the computed value,
       ie min( dilated value, reference value), is higher than the buffer value
    */
    if ( val > theBuf[i] ) {
      theBuf[i] = val;
      c ++;
    }
	  
  }
  
  ((_processListParam*)parameter)->changes = c;
  return( 1 );
}



static int _processListSubroutine_u16( void *parameter,
				      size_t first,
				      size_t last )
{
  typePointList *thePointList = ((_processListParam*)parameter)->thePointList;
  unsigned short int *theRef = (unsigned short int *)(((_processListParam*)parameter)->theInput);
  unsigned short int *theBuf = (unsigned short int *)(((_processListParam*)parameter)->theOutput);
  int dimx = ((_processListParam*)parameter)->theDim[0];
  int dimy = ((_processListParam*)parameter)->theDim[1];
  int dimz = ((_processListParam*)parameter)->theDim[2];
  int dimxy = dimx*dimy;
  
  typePoint *p;
  int x, y, z, i;
  int n, val;
  int j, k, l;
  int c = 0;

  for ( n=first; n<=last; n ++ ) {
    
    /* get the point
     */
    p = &(thePointList->pt[n]);
    x = p->x;
    y = p->y;
    z = p->z;
    i = p->i;
    val = theBuf[i];
	  
    /* dilate (26-connectivity)
       warning: in parallelism, buffer value are accessed outside the allocated range
       of each chunk -> we may have a different behavior than the sequential
       computation (but should be the same result)
       can be corrected if two buffers are used (one as input, the second one as output)
       warning: this may cause access to a data being written by an other chunk :(
     */
    
    if ( p->t ) {
      for ( l= -1; l<=1; l++ )
	for ( k= -1; k<=1; k++ )
	  for ( j= -1; j<=1; j++ ) {
	    if ( val < theBuf[i+l*dimxy+k*dimx+j] ) val = theBuf[i+l*dimxy+k*dimx+j];
	  }
    } 
    else {
      for ( l= -1; l<=1; l++ ) {
	if ( (z+l >= 0) && (z+l < dimz) ) {
	  for ( k= -1; k<=1; k++ ) {
	    if ( (y+k >= 0) && (y+k < dimy) ) {
	      for ( j= -1; j<=1; j++ ) {
		if ( (x+j >= 0) && (x+j < dimx) ) {
		  if ( val < theBuf[i+l*dimxy+k*dimx+j] ) val = theBuf[i+l*dimxy+k*dimx+j];
		}
	      }
	    }
	  }
	}
      }
    }
	  
    
    /* val is the value of the input image after dilation
       check the correspondant value of the reference image
    */
    if ( val > theRef[i] ) val = theRef[i];

    /* there is a change is the computed value,
       ie min( dilated value, reference value), is higher than the buffer value
    */
    if ( val > theBuf[i] ) {
      theBuf[i] = val;
      c ++;
    }
	  
  }
  
  ((_processListParam*)parameter)->changes = c;
  return( 1 );
}










static int process_list ( typePointList *thePointList,
			  void *theInput, 
			  void *theOutput, 
			  bufferType theType, 
			  int *theDim, 
			  int maxValue )
{
  char *proc = "process_list";

  int n;
  typeChunks chunks;
  size_t last;
  _processListParam aux;
  _processListParam *auxList = NULL;
  int use_parallelism = 0;
  int magic_number = 10;

  typePoint tmp;
  int i;
  int val;
  int c, changes = 0;
  int iteration = 0;


  if ( thePointList->nPoints <= 0 ) return( 0 );

  switch( theType ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
  case USHORT :
    break;
  }



  /* allocations
   */
  
  last = thePointList->nPoints-1;

  if ( thePointList->nPoints > getMaxChunks() ) {
    use_parallelism = 1;

    initChunks( &chunks );
    if ( buildChunks( &chunks, 0, last, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
    }
    
    auxList = (_processListParam*)malloc( chunks.n_allocated_chunks * sizeof(_processListParam) );
    if ( auxList == (_processListParam*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate auxiliary variables\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    
    for ( n=0; n<chunks.n_allocated_chunks; n++ ) {	
      auxList[n].thePointList = thePointList;
      auxList[n].theInput = theInput;
      auxList[n].theOutput = theOutput;
      auxList[n].theDim = theDim;
      
      chunks.data[n].parameters = (void*)(&auxList[n]);
    }
  }
    
  aux.thePointList = thePointList;
  aux.theInput = theInput;
  aux.theOutput = theOutput;
  aux.theDim = theDim;



  do {

    if ( use_parallelism == 1 ) {
      for ( n=0; n<chunks.n_allocated_chunks; n++ )
	auxList[n].changes = 0;
    }
    else {
      aux.changes = 0;
    }



    if ( _component_wise_processing_ == 1 && 
	 iteration == 0 && thePointList->nPoints > theDim[0]*theDim[1]*theDim[2]/10 ) {
      if ( _verbose_ >= 3 ) {
	fprintf( stderr, "    switch to component-wise propagation\n" );
      }
      c = process_component( thePointList, theInput, theOutput, theType, theDim, maxValue );
    }
    else {

      switch( theType ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such image type not handled yet (changes computation)\n", proc );
	return( -1 );
      case UCHAR :
	if ( use_parallelism == 1 ) {
	  if ( processChunks( &_processListSubroutine_u8, &chunks, proc ) != 1 ) {
	    if ( _verbose_ )
	      fprintf( stderr, "%s: unable to compute changes\n", proc );
	    free( auxList );
	    freeChunks( &chunks );
	    return( -1 );
	  }
	}
	else {
	(void)_processListSubroutine_u8( (void *)&aux, 0, last );
	}
	break;
      case USHORT :
	if ( use_parallelism == 1 ) {
	  if ( processChunks( &_processListSubroutine_u16, &chunks, proc ) != 1 ) {
	    if ( _verbose_ )
	      fprintf( stderr, "%s: unable to compute changes\n", proc );
	    free( auxList );
	    freeChunks( &chunks );
	    return( -1 );
	  }
	}
	else {
	  (void)_processListSubroutine_u16( (void *)&aux, 0, last );
	}
	break;
      }
      
      /* number of changes at this iteration
       */
      c = 0;
      if ( use_parallelism == 1 ) {
	for ( n=0; n<chunks.n_allocated_chunks; n++ )
	  c += auxList[n].changes;
      }
      else 
	c = aux.changes;

    }

    if ( _verbose_ >= 3 ) {
      fprintf( stderr,"      changes = %8d", c );
      if ( _verbose_ >= 4 )
	  fprintf( stderr, "\n" );
	else
	  fprintf( stderr, "\r" );
    }
    
    changes += c;
    
    /* prepare next iteration if necessary
     */
    
    if ( c > 0 ) {
      /* put points which won't change any more 
	 at the beginning of the list
      */
      switch( theType ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such image type not handled yet (list re-ordering)\n", proc );
	return( -1 );
      case UCHAR :
	{
	  unsigned char *theRef = (unsigned char *)theInput;
	  unsigned char *theBuf = (unsigned char *)theOutput;
	  for ( n=0; n<thePointList->nPoints;) {
	    i = thePointList->pt[n].i;
	    val = theBuf[i];
	    
	    if ( val == theRef[i] || val == maxValue ) {
	      tmp = thePointList->pt[ thePointList->nPoints-1 ];
	      thePointList->pt[ thePointList->nPoints-1 ] = thePointList->pt[n];
	      thePointList->pt[n] = tmp;
	      thePointList->nPoints --;
	    }
	    else {
	      n ++;
	    }
	  }
	}
	break;
      case USHORT :
	{
	  unsigned short int *theRef = (unsigned short int *)theInput;
	  unsigned short int *theBuf = (unsigned short int *)theOutput;
	  for ( n=0; n<thePointList->nPoints; ) {
	    i = thePointList->pt[n].i;
	    val = theBuf[i];
	    
	    if ( val == theRef[i] || val == maxValue ) {
	      tmp = thePointList->pt[ thePointList->nPoints-1 ];
	      thePointList->pt[ thePointList->nPoints-1 ] = thePointList->pt[n];
	      thePointList->pt[n] = tmp;
	      thePointList->nPoints --;
	    }
	    else {
	      n ++;
	    }
	  }
	}
	break;
      }

      if ( (size_t)thePointList->nPoints > magic_number * getMaxChunks() ) {
	if ( buildEqualChunks( &(chunks.data[0]), 0, (size_t)(thePointList->nPoints-1), 
			       chunks.n_allocated_chunks ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to re-compute chunks\n", proc );
	  free( auxList );
	  freeChunks( &chunks );
	  return( -1 );
	}
      }
      else {
	if ( use_parallelism == 1 ) {
	  free( auxList );
	  freeChunks( &chunks );
	}
	use_parallelism = 0;
      }

    }
    
    iteration ++;

  } while ( c > 0 && thePointList->nPoints > 0 );
  
  if ( use_parallelism == 1 ) {
    free( auxList );
    freeChunks( &chunks );
  }

  return( changes );

}










static int process_component( typePointList *thePointList,
			      void *theInput, 
			      void *theOutput, 
			      bufferType theType, 
			      int *theDim, 
			      int maxValue )
{
  char *proc = "process_component";

  unsigned short int *theCC = (unsigned short int *)NULL;
  unsigned short int _low_threshold_ = 100;
  unsigned short int _high_threshold_ = 200;

  int i, j, k, l, n;
  
  
  int dimx = theDim[0];
  int dimy = theDim[1];
  int dimz = theDim[2];
  int dimxy = dimx*dimy;
  int v = dimxy*dimz;

  typePoint *p;
  int x, y, z;
  int val, nval, maxHigh, minRef;
  int c = 0;

  /* initialization
   */
  theCC = (unsigned short int *)malloc( v * sizeof( unsigned short int ) );
  if ( theCC == (unsigned short int *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }
  
  for ( i=0; i<v; i++ ) theCC[i] = 0;
  



  
  /* 
     maxHigh = maximum value of high threshold points
     minRef  = maximum value of reference for low threshold points

     by dilation, all "low threshold" points can have the maxHigh value,
     however, the "minRef" points can make a barrier to the propagation,
     thus, the maximal elevation value (for the components)
     is min( maxHigh, minRef ) 

     this can be done for each individual component
  */

  /* low threshold = points of the list to be processed

     high threshold = neighbooring points with higher value
     all points could have been considered, but for dilation
     only values of neighboring points are considered
   */


  switch ( theType ) {
  default :
    free( theCC );
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      unsigned char *theRef = (unsigned char *)theInput;
      unsigned char *theBuf = (unsigned char *)theOutput;

      val = theBuf[ thePointList->pt[0].i ];
      maxHigh = val;
      minRef = maxValue;

      for ( n=0; n<thePointList->nPoints; n ++ ) {
	
	p = &(thePointList->pt[n]);
	x = p->x;
	y = p->y;
	z = p->z;
	i = p->i;
	
	theCC[ i ] = _low_threshold_;
	if ( minRef > theRef[ i ] ) minRef = theRef[ i ];
	
	if ( p->t ) {
	  for ( l= -1; l<=1; l++ )
	    for ( k= -1; k<=1; k++ )
	      for ( j= -1; j<=1; j++ ) {
		if ( theBuf[i+l*dimxy+k*dimx+j] > val ) {
		  theCC[ i+l*dimxy+k*dimx+j ] = _high_threshold_;
		  if ( maxHigh < theRef[i+l*dimxy+k*dimx+j] ) maxHigh = theRef[i+l*dimxy+k*dimx+j];
		}
	      }
	} 
	else {
	  for ( l= -1; l<=1; l++ ) {
	    if ( (z+l >= 0) && (z+l < dimz) ) {
	      for ( k= -1; k<=1; k++ ) {
		if ( (y+k >= 0) && (y+k < dimy) ) {
		  for ( j= -1; j<=1; j++ ) {
		    if ( (x+j >= 0) && (x+j < dimx) ) {
		      if ( theBuf[i+l*dimxy+k*dimx+j] > val ) {
			theCC[ i+l*dimxy+k*dimx+j ] = _high_threshold_;
			if ( maxHigh < theRef[i+l*dimxy+k*dimx+j] ) maxHigh = theRef[i+l*dimxy+k*dimx+j];
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

      }
    }
    break;
  case USHORT :
    {
      unsigned short int *theRef = (unsigned short int *)theInput;
      unsigned short int *theBuf = (unsigned short int *)theOutput;

      val = theBuf[ thePointList->pt[0].i ];
      maxHigh = val;
      minRef = maxValue;

      for ( n=0; n<thePointList->nPoints; n ++ ) {
	
	p = &(thePointList->pt[n]);
	x = p->x;
	y = p->y;
	z = p->z;
	i = p->i;
	
	theCC[ i ] = _low_threshold_;
	if ( minRef > theRef[ i ] ) minRef = theRef[ i ];
	
	if ( p->t ) {
	  for ( l= -1; l<=1; l++ )
	    for ( k= -1; k<=1; k++ )
	      for ( j= -1; j<=1; j++ ) {
		if ( theBuf[i+l*dimxy+k*dimx+j] > val ) {
		  theCC[ i+l*dimxy+k*dimx+j ] = _high_threshold_;
		  if ( maxHigh < theRef[i+l*dimxy+k*dimx+j] ) maxHigh = theRef[i+l*dimxy+k*dimx+j];
		}
	      }
	} 
	else {
	  for ( l= -1; l<=1; l++ ) {
	    if ( (z+l >= 0) && (z+l < dimz) ) {
	      for ( k= -1; k<=1; k++ ) {
		if ( (y+k >= 0) && (y+k < dimy) ) {
		  for ( j= -1; j<=1; j++ ) {
		    if ( (x+j >= 0) && (x+j < dimx) ) {
		      if ( theBuf[i+l*dimxy+k*dimx+j] > val ) {
			theCC[ i+l*dimxy+k*dimx+j ] = _high_threshold_;
			if ( maxHigh < theRef[i+l*dimxy+k*dimx+j] ) maxHigh = theRef[i+l*dimxy+k*dimx+j];
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}

      }
    }
    break;
  }



  if (  HysteresisThresholding( (void *)theCC, USHORT, (void *)theCC, USHORT, theDim,
				(double)_low_threshold_, (double)_high_threshold_ ) < 0 ) {
    free( theCC );
    if ( _verbose_ )
      fprintf( stderr, "%s: computation failed\n", proc );
    return( -1 );
  }

  nval = maxHigh;
  if ( nval > minRef ) nval = minRef;
  
  if ( _verbose_ >= 3 ) 
    fprintf( stderr, "    %s: changes intensity %d into %d\n", proc, val, nval );
  
  switch( theType ) {
  default :
    free( theCC );
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      unsigned char *theBuf = (unsigned char *)theOutput;

      for ( n=0; n<thePointList->nPoints; n ++ ) {
	i = thePointList->pt[n].i;
	if ( theCC[i] > 0 ) {
	  theBuf[i] = nval;
	  c ++;
	}
      }
    }
    break;
  case USHORT :
    {
      unsigned short int *theBuf = (unsigned short int *)theOutput;
      for ( n=0; n<thePointList->nPoints; n ++ ) {
	i = thePointList->pt[n].i;
	if ( theCC[i] > 0 ) {
	  theBuf[i] = nval;
	  c ++;
	}
      }
    }
    break;
   }

  free( theCC );
  return( c );
}
