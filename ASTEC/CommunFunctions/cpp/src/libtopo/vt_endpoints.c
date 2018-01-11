/*************************************************************************
 * vt_endpoints.c -
 *
 * $Id: vt_endpoints.c,v 1.3 2001/04/04 07:24:58 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Tue Mar 20 16:41:37 MET 2001
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 *
 */

#include <vt_endpoints.h>


static int _verbose_ = 1;



typedef struct {
  int x;
  int y;
  int z;
  int v;
  int n;
} typeMaxPoint;


typedef struct {
  int x;
  int y;
  int z;
  int v;
} typeQueuePoint;


typedef struct {
  int n;
  int nalloc;
  typeQueuePoint *thePts;
} typeQueue;


static int _alloc_in_queue_ = 1000;

static int _PutPointInQueue( int x, int y, int z, int v,  typeQueue *theQueue );

static int inc_chamfer_3x3x3[5] = { 16, 21, 16, 21, 26 };











int VT_ComputeBackDistanceInsideObjects( vt_image *theCC,
					 vt_image *theDist,
					 vt_image *theBackDist,
					 int *theIncrements )
{
  char *proc = "VT_ComputeBackDistanceInsideObjects";

  int nbConnectedComponents = 0;
  int i, j, k, v, d, p, q;
  int c, x, y, z;
  typeMaxPoint *theMaxPoint = NULL;
  int maxDistance;
  int mult=1, previousDistance;

  int _infinity_ = 0;
  int _already_in_queue_ = 0;
  int _background_ = 0;
  typeQueue *theQueue = NULL;

  int dimx, dimy, dimz;
  int inc[3][3][3];
  int *increments = NULL;

  /* tests
   */
  if ( theCC->dim.x != theDist->dim.x ||
       theCC->dim.y != theDist->dim.y ||
       theCC->dim.z != theDist->dim.z ||
       theCC->dim.x != theBackDist->dim.x ||
       theCC->dim.y != theBackDist->dim.y ||
       theCC->dim.z != theBackDist->dim.z ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: images have different dimensions\n", proc );
    }
    return( -1 );
  }

  dimx = theCC->dim.x;
  dimy = theCC->dim.y;
  dimz = theCC->dim.z;

  if ( theDist->type != theBackDist->type ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: both distance images should have the same type\n", proc );
    }
    return( -1 );
  }



  /* increments
   */
  if ( theIncrements != NULL )
    increments = theIncrements;
  else
    increments = inc_chamfer_3x3x3;

  inc[1][1][1] = 0;

  inc[1][1][0] = inc[1][1][2] = inc[1][0][1] = inc[1][2][1] = increments[0];
  inc[1][0][0] = inc[1][0][2] = inc[1][2][0] = inc[1][2][2] = increments[1];

  inc[0][1][1] = increments[2];
  inc[0][1][0] = inc[0][1][2] = inc[0][0][1] = inc[0][2][1] = increments[3];
  inc[0][0][0] = inc[0][0][2] = inc[0][2][0] = inc[0][2][2] = increments[4];
 
  inc[2][1][1] = increments[2];
  inc[2][1][0] = inc[2][1][2] = inc[2][0][1] = inc[2][2][1] = increments[3];
  inc[2][0][0] = inc[2][0][2] = inc[2][2][0] = inc[2][2][2] = increments[4];
 


  /* count the number of connected components
   */
  v = theCC->dim.x * theCC->dim.y * theCC->dim.z;
  switch( theCC->type ) {
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: connected components type not handled yet\n", proc );
    }
    return( -1 );
  case UCHAR : 
    {
      u8 *theBuf = (u8*)theCC->buf;
      for (i=0; i<v; i++)
	if ( nbConnectedComponents < theBuf[i] )
	  nbConnectedComponents = theBuf[i];
    }
    break;
   case USHORT : 
    {
      u16 *theBuf = (u16*)theCC->buf;
      for (i=0; i<v; i++)
	if ( nbConnectedComponents < theBuf[i] )
	  nbConnectedComponents = theBuf[i];
    }
    break;
  }


  /* allocate the array for the maximum points
   */
  
  theMaxPoint = (typeMaxPoint*)malloc( (nbConnectedComponents+1)*sizeof(typeMaxPoint) );
  if ( theMaxPoint == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate array of maximum points\n", proc );
    }
    return( -1 );
  }
  for (i=0; i<=nbConnectedComponents; i++ ) {
    theMaxPoint[i].x = -1;
    theMaxPoint[i].y = -1;
    theMaxPoint[i].z = -1;
    theMaxPoint[i].v = 0; 
    theMaxPoint[i].n = 0; 
  }

  /* find the maximum for each connected components
   */

  switch( theDist->type ) {

  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: distance type not handled yet\n", proc );
    }
    free( theMaxPoint );
    return( -1 );

  case USHORT :
    {
      u16 ***dstBuf = (u16***)theDist->array;
      
      switch( theCC->type ) {
      default :
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: connected components type not handled yet (2)\n", proc );
	}
	free( theMaxPoint );
	return( -1 );
      case UCHAR :
	{
	  u8 ***theBuf = (u8***)theCC->array;
	  for ( z=0; z<dimz; z++ )
	  for ( y=0; y<dimy; y++ )
	  for ( x=0; x<dimx; x++ ) {
	    if ( theBuf[z][y][x] == 0 ) continue;
	    i = theBuf[z][y][x];
	    theMaxPoint[i].n ++;
	    if ( theMaxPoint[i].v < dstBuf[z][y][x] ) {
	      theMaxPoint[i].v = dstBuf[z][y][x];
	      theMaxPoint[i].x = x;
	      theMaxPoint[i].y = y;
	      theMaxPoint[i].z = z;
	    }
	  }
	}
	break;
      case USHORT :
	{
	  u16 ***theBuf = (u16***)theCC->array;
	  for ( z=0; z<dimz; z++ )
	  for ( y=0; y<dimy; y++ )
	  for ( x=0; x<dimx; x++ ) {
	    if ( theBuf[z][y][x] == 0 ) continue;
	    i = theBuf[z][y][x];
	    theMaxPoint[i].n ++;
	    if ( theMaxPoint[i].v < dstBuf[z][y][x] ) {
	      theMaxPoint[i].v = dstBuf[z][y][x];
	      theMaxPoint[i].x = x;
	      theMaxPoint[i].y = y;
	      theMaxPoint[i].z = z;
	    }
	  }
	}
	break;
       }

    }
    break;
  }




  /* maximum of Distance
   */
  maxDistance = 0;
  for (i=1; i<=nbConnectedComponents; i++ ) {
    if ( maxDistance < theMaxPoint[i].v )
      maxDistance = theMaxPoint[i].v;
  }



  /* allocations of queues and initialization
   */
  theQueue = (typeQueue*)malloc( (maxDistance+1)*sizeof(typeQueue) );
  if ( theQueue == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate array of queues\n", proc );
    }
    free( theMaxPoint );
    return( -1 );
  }
  for ( i=0; i<=maxDistance; i++ ) {
    theQueue[i].n = 0;
    theQueue[i].nalloc = 0;
    theQueue[i].thePts = NULL;
  }

  
  

  




  /* initialization of back distance map
   */
  
  switch( theBackDist->type ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: back distance type not handled yet\n", proc );
    }
    free( theQueue );
    free( theMaxPoint );
    return( -1 );
    
  case USHORT :
    {
      u16 *dstBuf = (u16*)theBackDist->buf;
      _infinity_         = 65000;
      _already_in_queue_ = 65100;
      _background_       = 65200;
      for (i=0; i<v; i++) dstBuf[i] = _background_;
    }
    break;
  }




  for ( c=1; c<= nbConnectedComponents; c++ ) {

    /* component #c does not exist
     */
    if ( theMaxPoint[c].x < 0 ||
	 theMaxPoint[c].y < 0 ||
	 theMaxPoint[c].z < 0 ||
	 theMaxPoint[c].v == 0 ||
	 theMaxPoint[c].n == 0 )
      continue;


    fprintf( stderr, " processing component #%d/%d\n", c, nbConnectedComponents );


    /* re-initialization of queues
     */
    for ( i=0; i<=maxDistance; i++ )
      theQueue[i].n = 0;
    



    /* initialization of back distance map
       with component #c
    */
    
    switch( theBackDist->type ) {
    
    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: back distance type not handled yet (2)\n", proc );
      }
      free( theQueue );
      free( theMaxPoint );
      return( -1 );
      
    case USHORT :
      {
	u16 *dstBuf = (u16*)theBackDist->buf;

	switch( theCC->type ) {
	default :
	  if ( _verbose_ ) {
	    fprintf( stderr, "%s: connected components type not handled yet (4)\n", proc );
	  }
	  free( theQueue );
	  free( theMaxPoint );
	  return( -1 );
	case UCHAR :
	  {
	    u8 *theBuf = (u8*)theCC->buf;
	    for (i=0; i<v; i++) 
	      if ( theBuf[i] == c ) dstBuf[i] = _infinity_;
	  }
	  break;
	case USHORT :
	  {
	    u16 *theBuf = (u16*)theCC->buf;
	    for (i=0; i<v; i++) 
	      if ( theBuf[i] == c ) dstBuf[i] = _infinity_;
	  }
	  break;
	}

      }
      break;
      
    }
    



    


    /* proceed
     */
    switch( theBackDist->type ) {
    
    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: back distance type not handled yet (3)\n", proc );
      }
      free( theQueue );
      free( theMaxPoint );
      return( -1 );
      
    case USHORT :
      {
	u16 ***dstBuf = (u16***)theDist->array;
	u16 ***bckBuf = (u16***)theBackDist->array;

	
	/* put neighbors of the maximum into the queue
	   
	   distance value is pu to 1 not to 0,
	   to avoid further problems (misleading with the background)
	 */
	x = theMaxPoint[c].x;
	y = theMaxPoint[c].y;
	z = theMaxPoint[c].z;

	bckBuf[ z ][ y ][ x ] = 0;
	
	for ( k = -1; k <= 1; k++ ) {
	  if ( z+k < 0 || z+k >= dimz ) continue;
	  for ( j = -1; j <= 1; j++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -1; i <= 1; i++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( bckBuf[z+k][y+j][x+i] != _infinity_ ) continue;
	      if ( _PutPointInQueue( x+i, y+j, z+k, (int)dstBuf[z+k][y+j][x+i], theQueue ) == 0 ) {
		for ( q=0; q<=maxDistance; q++ )
		  if ( theQueue[q].thePts != NULL ) free( theQueue[q].thePts );
		free( theQueue );
		free( theMaxPoint );
		return( -1 );
	      }
	      bckBuf[z+k][y+j][x+i] = _already_in_queue_;
	    }
	  }
	}


	previousDistance = dstBuf[z][y][x];

	
	/* loop for component #c
	 */
	do {
	  
	  for ( q=0, d=maxDistance; d > 0 && q == 0; d-- )
	    if ( theQueue[d].n > 0 ) q = d;

	  if ( q == 0 ) continue;

	  mult = ( previousDistance <= q ) ? 2 : 1;
	  mult = 1;
	  fprintf( stderr, " cc #%3d : distance = %5d (prev=%5d m=%d)  points = %5d\r", 
		   c, q, previousDistance, mult, theQueue[q].n );
	  previousDistance = q;
	    
	  
	  /* copie de la queue a traiter en #0
	   */
	  for ( p=0; p<theQueue[q].n; p++ ) {
	    if ( _PutPointInQueue( theQueue[q].thePts[p].x,
				   theQueue[q].thePts[p].y,
				   theQueue[q].thePts[p].z, 0, theQueue ) == 0 ) {
	      for ( d=0; d<=maxDistance; d++ )
		if ( theQueue[d].thePts != NULL ) free( theQueue[d].thePts );
	      free( theQueue );
	      free( theMaxPoint );
	      return( -1 );
	    }
	  }
	  theQueue[q].n = 0;
	  
	  /* on traite les points de la queue #q (maintenant 0)
	     1. on calcule les distances
	     2. on ajoute ses voisins
	  ...
	   */
	  
	  for ( p=0; p<theQueue[0].n; p++ ) {
	    x = theQueue[0].thePts[p].x;
	    y = theQueue[0].thePts[p].y;
	    z = theQueue[0].thePts[p].z;
	    v = _infinity_;

	    for ( k = -1; k <= 1; k++ ) {
	      if ( z+k < 0 || z+k >= dimz ) continue;
	      for ( j = -1; j <= 1; j++ ) {
		if ( y+j < 0 || y+j >= dimy ) continue;
		for ( i = -1; i <= 1; i++ ) {
		  if ( x+i < 0 || x+i >= dimx ) continue;
		  if ( i == 0 && j == 0 && k == 0 ) continue;
		  if ( bckBuf[z+k][y+j][x+i] == _infinity_ ||
		       bckBuf[z+k][y+j][x+i] == _already_in_queue_ ||
		       bckBuf[z+k][y+j][x+i] == _background_ ) continue;
		  if ( v > bckBuf[z+k][y+j][x+i]+mult*inc[1+k][1+j][1+i] )
		    v = bckBuf[z+k][y+j][x+i]+mult*inc[1+k][1+j][1+i];
		}
	      }
	    }
	    theQueue[0].thePts[p].v = v;
	  }
	  
	  
	  for ( p=0; p<theQueue[0].n; p++ ) {
	    x = theQueue[0].thePts[p].x;
	    y = theQueue[0].thePts[p].y;
	    z = theQueue[0].thePts[p].z;
	    bckBuf[z][y][x] = theQueue[0].thePts[p].v;

	    for ( k = -1; k <= 1; k++ ) {
	      if ( z+k < 0 || z+k >= dimz ) continue;
	      for ( j = -1; j <= 1; j++ ) {
		if ( y+j < 0 || y+j >= dimy ) continue;
		for ( i = -1; i <= 1; i++ ) {
		  if ( x+i < 0 || x+i >= dimx ) continue;
		  if ( bckBuf[z+k][y+j][x+i] != _infinity_ ) continue;
		  if ( _PutPointInQueue( x+i, y+j, z+k, (int)dstBuf[z+k][y+j][x+i], theQueue ) == 0 ) {
		    for ( d=0; d<=maxDistance; d++ )
		      if ( theQueue[d].thePts != NULL ) free( theQueue[d].thePts );
		    free( theQueue );
		    free( theMaxPoint );
		    return( -1 );
		  }
		  bckBuf[z+k][y+j][x+i] = _already_in_queue_;
		}
	      }
	    }
	  }
	  theQueue[0].n = 0;
	  
	} while ( q > 0 );


      } /* switch( theBackDist->type ) case USHORT */
      break; 

    } /* switch( theBackDist->type ) */
      




    

  } /* for ( c=1; c<= nbConnectedComponents; c++ ) */
  fprintf( stderr, "\n" );



  for ( d=0; d<=maxDistance; d++ )
    if ( theQueue[d].thePts != NULL ) free( theQueue[d].thePts );
  free( theQueue );
  free( theMaxPoint );
 

  v = dimx * dimy * dimz;
  switch( theBackDist->type ) {
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: back distance type not handled yet\n", proc );
    }
    return( -1 );
  case USHORT :
    {
      u16 *dstBuf = (u16*)theBackDist->buf;
      for (i=0; i<v; i++) {
	if ( dstBuf[i] == 0 )
	  dstBuf[i] = 1;
	else if ( dstBuf[i] == _background_ )
	  dstBuf[i] = 0;
      }
    }
    break;
  }



  return( 1 );
}





static int _PutPointInQueue( int x, int y, int z, int v,  typeQueue *theQueue )
{
  int i, n;
  typeQueuePoint *tmp = NULL;
  
  if ( theQueue[v].n == theQueue[v].nalloc ) {

    n = theQueue[v].nalloc + _alloc_in_queue_;
    tmp = (typeQueuePoint*)malloc( n*sizeof(typeQueuePoint) );
    if ( tmp == NULL ) {
      return( 0 );
    }
    (void)memcpy( tmp, theQueue[v].thePts, theQueue[v].nalloc*sizeof(typeQueuePoint) );
    free( theQueue[v].thePts );
    theQueue[v].thePts = tmp;
    theQueue[v].nalloc = n;

  }

  i = theQueue[v].n;
  theQueue[v].thePts[ i ].x = x;
  theQueue[v].thePts[ i ].y = y;
  theQueue[v].thePts[ i ].z = z;
  theQueue[v].n ++;
  return( 1 );
}
  







int VT_InitialiseImageToBeThinned( vt_image *theBackDist,
				   vt_image *theImToBeThinned,
				   int halfWindowSize,
				   int *theIncrements )
{
  char *proc = "VT_InitialiseImageToBeThinned";
  int x, y, z;
  int i, j, k, chg;
  int mayBeMax, mayBeMin, n=0;

  int wdim, dimx, dimy, dimz;
  int inc[3][3][3];
  int *increments = NULL;

  int l=halfWindowSize;

  vt_image imTmp, imFlag;
  u8 ***theFlag;

  u8 ***theBuf = (u8***)theImToBeThinned->array;

  if ( theBackDist->dim.x != theImToBeThinned->dim.x ||
       theBackDist->dim.y != theImToBeThinned->dim.y ||
       theBackDist->dim.z != theImToBeThinned->dim.z ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: images have different dimensions\n", proc );
    }
    return( -1 );
  }
   
  if ( theImToBeThinned->type != UCHAR ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: output image should be of unsigned char type\n", proc );
    }
    return( -1 );
  }
  
  dimx = theImToBeThinned->dim.x;
  dimy = theImToBeThinned->dim.y;
  dimz = theImToBeThinned->dim.z;



  wdim = 2*halfWindowSize+1;
  VT_InitImage( &imTmp,  NULL, wdim, wdim, wdim, theBackDist->type );
  VT_InitImage( &imFlag, NULL, wdim, wdim, wdim, UCHAR );
  if ( dimz <3 ) {
    imTmp.dim.z = imFlag.dim.z = 1;
  }
  if ( VT_AllocImage( &imTmp ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate first auxiliary subimage\n", proc );
    }
    return( -1 );
  }
  if ( VT_AllocImage( &imFlag ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: unable to allocate second auxiliary subimage\n", proc );
    }
    VT_FreeImage( &imTmp );
    return( -1 );
  }
  theFlag = (u8***)imFlag.array;


  /* increments
   */
  if ( theIncrements != NULL )
    increments = theIncrements;
  else
    increments = inc_chamfer_3x3x3;

  inc[1][1][1] = 0;

  inc[1][1][0] = inc[1][1][2] = inc[1][0][1] = inc[1][2][1] = increments[0];
  inc[1][0][0] = inc[1][0][2] = inc[1][2][0] = inc[1][2][2] = increments[1];

  inc[0][1][1] = increments[2];
  inc[0][1][0] = inc[0][1][2] = inc[0][0][1] = inc[0][2][1] = increments[3];
  inc[0][0][0] = inc[0][0][2] = inc[0][2][0] = inc[0][2][2] = increments[4];
 
  inc[2][1][1] = increments[2];
  inc[2][1][0] = inc[2][1][2] = inc[2][0][1] = inc[2][2][1] = increments[3];
  inc[2][0][0] = inc[2][0][2] = inc[2][2][0] = inc[2][2][2] = increments[4];


 
  for ( z=0; z<dimz; z++ )
  for ( y=0; y<dimy; y++ )
  for ( x=0; x<dimx; x++ ) {
    theBuf[z][y][x] = 0;
  }



  switch( theBackDist->type ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: back distance type not handled yet\n", proc );
    }
    VT_FreeImage( &imFlag );
    VT_FreeImage( &imTmp );
    return( -1 );

  case USHORT :
    {
      u16 ***theDist = (u16***)theBackDist->array;
      u16 ***theTmp  = (u16***)imTmp.array;
      
      if ( dimz >= 3 ) {
	for ( z=0; z<dimz; z++ ) {
	  fprintf( stderr, "processing slice %3d/%d\r", z+1,  dimz );
	  for ( y=0; y<dimy; y++ )
	  for ( x=0; x<dimx; x++ ) {
	    if ( theDist[z][y][x] == 0 ) {
	      continue;
	    }
	    
	    /* init
	     */
	    for ( k = 0; k < wdim; k ++ ) 
	      for ( j = 0; j < wdim; j ++ ) 
		for ( i = 0; i < wdim; i ++ ) {
		  theTmp[k][j][i] = 0;
		  theFlag[k][j][i] = 0;
		}
	    theBuf[z][y][x] = 1;
	    
	    
	    /* on remplit les imagettes
	     */
	    if ( x==0 || x==dimx-1 || y==0 || y==dimy-1 || z==0 || z==dimz-1 ) {
	      for ( k = -l; k <= l; k++ ) {
		if ( z+k < 0 || z+k >= dimz ) continue;
		for ( j = -l; j <= l; j++ ) {
		  if ( y+j < 0 || y+j >= dimy ) continue;
		  for ( i = -l; i <= l; i++ ) {
		    if ( z+k != 0 && z+k != dimz-1 &&
			 y+j != 0 && y+j != dimy-1 &&
			 x+i != 0 && x+i != dimx-1 ) continue;
		    theTmp[l+k][l+j][l+i] = theDist[z+k][y+j][x+i];
		    theFlag[l+k][l+j][l+i] = 1;
		  }
		}
	      }
	      
	    } else {
	      for ( k = -l; k <= l; k++ ) {
		if ( z+k < 0 || z+k >= dimz ) continue;
		for ( j = -l; j <= l; j++ ) {
		  if ( y+j < 0 || y+j >= dimy ) continue;
		  for ( i = -l; i <= l; i++ ) {
		    if ( x+i < 0 || x+i >= dimx ) continue;
		    theTmp[l+k][l+j][l+i] = theDist[z+k][y+j][x+i];
		    theFlag[l+k][l+j][l+i] = 1;
		  }
		}
	      }
	    }
	    
	    
	    /* on extrait la 6-composante
	       contenant le point central
	    */
	    theFlag[l][l][l] = 255;
	    do {
	      for ( chg = 0, k = 0; k < wdim; k ++ ) 
	      for ( j = 0; j < wdim; j ++ ) 
	      for ( i = 0; i < wdim; i ++ ) {
		if ( theFlag[k][j][i] == 0 ) continue;
		if ( theFlag[k][j][i] == 255 ) continue;
		if ( (k > 0      && theFlag[k-1][j][i] == 255) ||
		     (k < wdim-1 && theFlag[k+1][j][i] == 255) ||
		     (j > 0      && theFlag[k][j-1][i] == 255) ||
		     (j < wdim-1 && theFlag[k][j+1][i] == 255) ||
		     (i > 0      && theFlag[k][j][i-1] == 255) ||
		     (i < wdim-1 && theFlag[k][j][i+1] == 255) ) {
		  theFlag[k][j][i] = 255;
		  chg++;
		}
	      }
	      for ( k = wdim-1; k >= 0; k -- ) 
	      for ( j = wdim-1; j >= 0; j -- ) 
	      for ( i = wdim-1; i >= 0; i -- ) {
		if ( theFlag[k][j][i] == 0 ) continue;
		if ( theFlag[k][j][i] == 255 ) continue;
		if ( (k > 0      && theFlag[k-1][j][i] == 255) ||
		     (k < wdim-1 && theFlag[k+1][j][i] == 255) ||
		     (j > 0      && theFlag[k][j-1][i] == 255) ||
		     (j < wdim-1 && theFlag[k][j+1][i] == 255) ||
		     (i > 0      && theFlag[k][j][i-1] == 255) ||
		     (i < wdim-1 && theFlag[k][j][i+1] == 255) ) {
		  theFlag[k][j][i] = 255;
		  chg++;
		}
	      }
	    } while ( chg > 0 );
	    
	    /* on repercute sur les distances
	     */ 
	    for ( k = 0; k < wdim; k ++ ) 
	    for ( j = 0; j < wdim; j ++ ) 
	    for ( i = 0; i < wdim; i ++ ) 
	      if ( theFlag[k][j][i] != 255 ) 
		theTmp[k][j][i] = 0;
	    

	    /* on recherche les extremum
	     */


	    
	    /* bord de l'image : 
	       on cherche un minimum
	    */
	    if ( x==0 || x==dimx-1 || y==0 || y==dimy-1 || z==0 || z==dimz-1 ) {
	      mayBeMin = 1;
	      for ( k = -l; k <= l && mayBeMin; k++ ) 
	      for ( j = -l; j <= l && mayBeMin; j++ )
	      for ( i = -l; i <= l && mayBeMin; i++ ) {

		if ( k == 0 && j == 0 && i == 0 ) continue;
		if ( theTmp[l+k][l+j][l+i] == 0 ) continue;
		if ( theTmp[l][l][l] > theTmp[l+k][l+j][l+i] ) {
		  mayBeMin = 0;
		  continue;
		}
		if ( theTmp[l][l][l] == theTmp[l+k][l+j][l+i] ) {
		  if ( theBuf[z+k][y+j][x+i] == 255 ) mayBeMin = 0;
		  continue;
		}
		if ( i >= -1 && i <= 1 && j >= -1 && j <= 1 && k >= -1 && k <= 1 )
		  if ( theTmp[l+k][l+j][l+i] > 0 &&
		       theTmp[l][l][l] < theTmp[l+k][l+j][l+i]
		       - 2*inc[1+k][1+j][1+i] ) {
		    mayBeMin = 0;
		    continue;
		  }
	      }
	      
	      if ( mayBeMin == 1 ) {
		theBuf[z][y][x] = 255;
		n ++;
		if ( _verbose_ ) {
		  fprintf( stdout, "found border min #%2d at (%3d %3d %3d) = %5d\n",
			   n, x, y, z, theDist[z][y][x]  );
		}
	      }
	      continue;
	    }
	    
	    /* interieur de l'image :
	       on cherche un maximum
	    */
	    
	    mayBeMax = 1;
	    for ( k = -l; k <= l && mayBeMax; k++ )
	    for ( j = -l; j <= l && mayBeMax; j++ )
	    for ( i = -l; i <= l && mayBeMax; i++ ) {
	      if ( k == 0 && j == 0 && i == 0 ) continue;
	      if ( theTmp[l+k][l+j][l+i] == 0 ) continue;
	      if ( theTmp[l][l][l] < theTmp[l+k][l+j][l+i] ) {
		mayBeMax = 0;
		continue;
	      }
	      if ( theTmp[l][l][l] == theTmp[l+k][l+j][l+i] ) {
		if ( theBuf[z+k][y+j][x+i] == 255 ) mayBeMax = 0;
		continue;
	      }
	    if ( i >= -1 && i <= 1 && j >= -1 && j <= 1 && k >= -1 && k <= 1 )
	      if ( theTmp[l][l][l] > theTmp[l+k][l+j][l+i]
		   + 2*inc[1+k][1+j][1+i] ) {
		mayBeMax = 0;
		continue;
	      }
	    }
	    
	    if ( mayBeMax == 1 ) {
	      theBuf[z][y][x] = 255;
	      n ++;
	      if ( _verbose_ ) {
		fprintf( stdout, "found inner max #%3d at (%3d %3d %3d)\n",
			 n, x, y, z );
	      }
	    }
	  }
	}

      } else {

	for ( z=0; z<dimz; z++ )
	for ( y=0; y<dimy; y++ )
	for ( x=0; x<dimx; x++ ) {
	  if ( theDist[z][y][x] == 0 ) {
	    continue;
	  }
	  
	  theBuf[z][y][x] = 1;
	  

	  /* bord de l'image : 
	     on cherche un minimum
	   */
	  if ( x==0 || x==dimx-1 || y==0 || y==dimy-1 ) {
	    mayBeMin = 1;
	    for ( j = -l; j <= l && mayBeMin; j++ ) {
	      if ( y+j < 0 || y+j >= dimy ) continue;
	      for ( i = -l; i <= l && mayBeMin; i++ ) {
		if ( y+j != 0 && y+j != dimy-1 &&
		     x+i != 0 && x+i != dimx-1 ) continue;
		if ( j == 0 && i == 0 ) continue;
		if ( theDist[z][y+j][x+i] == 0 ) continue;
		if ( theDist[z][y][x] > theDist[z][y+j][x+i] ) {
		  mayBeMin = 0;
		  continue;
		}
		if ( theDist[z][y][x] == theDist[z][y+j][x+i] ) {
		  if ( theBuf[z][y+j][x+i] == 255 ) mayBeMin = 0;
		  continue;
		}
		if ( i >= -1 && i <= 1 && j >= -1 && j <= 1 )
		  if ( theDist[z][y][x] < theDist[z][y+j][x+i] 
		       - 2*inc[1][1+j][1+i] ) {
		    mayBeMin = 0;
		    continue;
		  }
	      }
	    }
	    if ( mayBeMin == 1 ) {
	      theBuf[z][y][x] = 255;
	      n ++;
	      if ( _verbose_ ) {
		fprintf( stdout, "found border min #%2d at (%3d %3d %3d) = %5d\n",
			 n, x, y, z, theDist[z][y][x]  );
	      }
	    }
	    continue;
	  }

	  /* interieur de l'image :
	     on cherche un maximum
	   */

	  mayBeMax = 1;
	  for ( j = -l; j <= l && mayBeMax; j++ ) {
	    if ( y+j < 0 || y+j >= dimy ) continue;
	    for ( i = -l; i <= l && mayBeMax; i++ ) {
	      if ( x+i < 0 || x+i >= dimx ) continue;
	      if ( j == 0 && i == 0 ) continue;
	      if ( theDist[z][y+j][x+i] == 0 ) continue;
	      if ( theDist[z][y][x] < theDist[z][y+j][x+i] ) {
		mayBeMax = 0;
		continue;
	      }
	      if ( theDist[z][y][x] == theDist[z][y+j][x+i] ) {
		if ( theBuf[z][y+j][x+i] == 255 ) mayBeMax = 0;
		continue;
	      }
	      if ( i >= -1 && i <= 1 && j >= -1 && j <= 1 )
		if ( theDist[z][y][x] > theDist[z][y+j][x+i] 
		     + 2*inc[1][1+j][1+i] ) {
		  mayBeMax = 0;
		  continue;
		}
	    }
	  }
	  
	  if ( mayBeMax == 1 ) {
	    theBuf[z][y][x] = 255;
	    n ++;
	    if ( _verbose_ ) {
	      fprintf( stdout, "found inner max #%2d at (%3d %3d %3d) = %5d\n",
		       n, x, y, z, theDist[z][y][x]  );
	    }
	    if ( 0 ) 
	    {
	      for ( j = -3; j <= 3; j++ ) {
		for ( i = -3; i <= 3; i++ ) {
		  printf( " %5d ",theDist[z][y+j][x+i] );
		}
		printf( "   ---   " );
		for ( i = -3; i <= 3; i++ ) {
		  if ( theDist[z][y+j][x+i] > 0 )
		    printf( " %5d ",theDist[z][y][x]-theDist[z][y+j][x+i] );
		  else
		    printf( "     . " );
		}
			
		printf( "\n" );
	      }
	    }
	  }
	}

      }
      
    }
    break;
  } 
  
  fprintf( stdout, " number of found max = %d\n", n );

  VT_FreeImage( &imFlag );
  VT_FreeImage( &imTmp );
  return( n );
}
