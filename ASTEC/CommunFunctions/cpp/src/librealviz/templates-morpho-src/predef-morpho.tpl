#include <typedefs.h>
#include <string.h>
#include <stdlib.h>


/*
 * operation with a centered SE a length 3 along X
 */

static int _TYPE_predef_binary_X_OPERATION( TYPE* inputSlice,
				       TYPE* outputSlice,
				       unsigned int dimx,
				       unsigned int dimy )
{
  unsigned int lineLength = dimx * sizeof( TYPE );
  unsigned int ilength, i, iy, x, y;
  
  int *l1, *l2;
  int lx;
  
  if ( inputSlice == outputSlice ) return( -1 );

  /* 
     we have to compare input:line[y][0] with input:line[y][1]
     for a length of (dimx-1)
     then output:line[y][1] with input:line[y][0]
     for the same length
  */
  
  ilength = ( dimx-1 ) * sizeof( TYPE ) / sizeof( int );
  lx      = ilength * sizeof( int ) / sizeof( TYPE );
  
  memcpy( outputSlice, inputSlice, dimy * lineLength );

  for ( y = 0; y < dimy; y ++ ) {

    iy = y*dimx;

    l2 = (int*)( &outputSlice[iy] );
    l1 = (int*)( &inputSlice [iy+1] );
    for ( i=0; i<ilength; i++ ) *l2++ _BINOP_= *l1++;
    for ( x=lx; x<dimx-1; x++ ) {
      if ( inputSlice[iy+x+1] _TEST_ outputSlice[iy+x] )
	outputSlice[iy+x] = inputSlice[iy+x+1];
    }
    
    l2 = (int*)( &outputSlice[iy+1] );
    l1 = (int*)( &inputSlice [iy] );
    for ( i=0; i<ilength; i++ ) *l2++ _BINOP_= *l1++;
    for ( x=lx; x<dimx-1; x++ ) {
      if ( inputSlice[iy+x] _TEST_ outputSlice[iy+x+1] )
	outputSlice[iy+x+1] = inputSlice[iy+x];
    }
  }
  
  return( 1 );
}






/*
 * operation with a centered SE a length 3 along Y
 */

static int _TYPE_predef_binary_Y_OPERATION( TYPE* inputSlice,
					    TYPE* centerSlice,
					    TYPE* outputSlice,
					    unsigned int dimx,
					    unsigned int dimy )
{
  unsigned int sliceLength = dimx * dimy * sizeof( TYPE );
  unsigned int ilength, i, x;
  
  int *l1, *l2;
  int lx;
  
  if ( inputSlice == outputSlice ) return( -1 );

  /*
    we have to compare input:line[0][0] with input:line[1][0]
    for a length of (dimy-1)*dimx
    then output:line[1][0] with input:line[0][0]
    for the same length
  */
  
  ilength = ( dimy-1 ) * dimx * sizeof( TYPE ) / sizeof( int );
  lx      = ilength * sizeof( int ) / sizeof( TYPE );
  
  if ( outputSlice != centerSlice )
    memcpy( outputSlice, centerSlice, sliceLength );

  l2 = (int*)(  outputSlice );
  l1 = (int*)( &inputSlice [dimx] );
  for ( i=0; i<ilength; i++ ) *l2++ _BINOP_= *l1++;

  for ( x=lx; x<(dimy-1)*dimx; x++ ) {
    if ( inputSlice[dimx+x] _TEST_ outputSlice[x] )
      outputSlice[x] = inputSlice[dimx+x];
  }
  
  l2 = (int*)( &outputSlice[dimx] );
  l1 = (int*)(  inputSlice  );
  for ( i=0; i<ilength; i++ ) *l2++ _BINOP_= *l1++;

  for ( x=lx; x<(dimy-1)*dimx; x++ ) {
    if ( inputSlice[x] _TEST_ outputSlice[dimx+x] )
      outputSlice[dimx+x] = inputSlice[x];
  }
  
  return( 1 );
}













int _TYPE_predef_binary_OPERATION( TYPE* inputBuf, /* buffer to be resampled */
				   TYPE* resultBuf, /* result buffer */
				   unsigned int dimx,
				   unsigned int dimy,
				   unsigned int dimz,
				   int connectivity, /* connectivity to be used */ 
				   int iterations  /* number of iterations */ )
{
  int conn = 0; 
  unsigned int z;
  int sliceSize = dimx * dimy * sizeof( TYPE );

  TYPE *theBuf;
  TYPE *resBuf;
  TYPE *input = inputBuf;
  TYPE *localBuf = NULL;
  TYPE *localTab[2];


  int iter;

  if ( iterations < 0 ) return( -1 );
  if ( iterations == 0 ) {
    if ( inputBuf != resultBuf ) {
      memcpy( resultBuf, inputBuf, sliceSize*dimz );
    }
    return( 1 );
  }

  /* test on connectivity */
  conn = connectivity;
  switch ( conn ) {
  case 4 :
  case 8 :
    break ;
  default :
    conn = 8;
  }


  
  
  localBuf = (TYPE*)malloc( 2 * sliceSize );
  if ( localBuf == NULL ) return( -1 );
  for ( z = 0; z < 2; z ++ ) {
    localTab[z] = &localBuf[z*dimx*dimy];
  }



  for ( iter = 0; iter < iterations; iter ++ ) {
    theBuf = input;
    resBuf = resultBuf;
    switch ( conn ) {
    case 4 :
      /*
	pour chaque plan 
	- operation selon X -> tmp1
	- operation selon Y -> tmp2
	- copie dans resultat
      */
      for ( z=0; z<dimz; z++ ) {
	_TYPE_predef_binary_X_OPERATION( &theBuf[z*dimx*dimy], 
					 localTab[0], dimx, dimy );
	_TYPE_predef_binary_Y_OPERATION( &theBuf[z*dimx*dimy], localTab[0], 
					 localTab[1], 
					 dimx, dimy );
	memcpy( &resBuf[z*dimx*dimy], localTab[1], sliceSize );
      }
      break;
    case 8 :
      for ( z=0; z<dimz; z++ ) {
	_TYPE_predef_binary_X_OPERATION( &theBuf[z*dimx*dimy], 
					 localTab[0], dimx, dimy );
	_TYPE_predef_binary_Y_OPERATION( localTab[0], localTab[0], 
					 localTab[1], 
					 dimx, dimy );
	memcpy( &resBuf[z*dimx*dimy], localTab[1], sliceSize );
      }
      break;
    }
    input = resultBuf;
  }

  return( 1 );
}








