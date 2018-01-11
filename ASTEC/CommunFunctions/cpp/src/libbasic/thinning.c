/****************************************************
 * thinning.c - 
 *
 * $Id: thinning.c,v 1.7 2001/05/09 15:55:48 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Mon Aug  7 14:56:09 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <thinning.h>




static int _verbose_ = 1;
static int _time_verbose_ = 0;




typedef enum {
  _BACKGROUND_ = 0,
  _TOBEDELETED_ = 150,
  _CANBEDELETED_ = 200,
  _ENDPOINT_ = 225,
  _ANCHOR_ = 255
} enumTypePoint;



#ifdef _UNUSED_
static void _printTypePoint( FILE *f, int s )
{
  switch ( s ) {
  default : fprintf( f, "unknown" );
  case _BACKGROUND_   : fprintf( f, "BACKGROUND" ); break;
  case _TOBEDELETED_  : fprintf( f, "TOBEDELETED" ); break;
  case _CANBEDELETED_ : fprintf( f, "CANBEDELETED" ); break;
  case _ENDPOINT_     : fprintf( f, "ENDPOINT" ); break;
  case _ANCHOR_       : fprintf( f, "ANCHOR" ); break;
  }
}
#endif




typedef struct {
  int x;
  int y;
  int z;
  int i;
  int value;

  char isinside;
  enumTypePoint type;
} typePoint;







static void _ComputeTimeFrom2Clocks( float c1, float c2,
				     int *hours, int *minutes, float *seconds );
static void _PrintTimeFrom2Clocks( char *str, float c1, float c2 );






/**************************************************
 *
 *
 *
 **************************************************/


void initTypeThinningParameters( typeThinningParameters *p )
{
  p->anchorValue = -1;
  p->connectivity = 26;
  p->cyclesBeforeEnding = -1;
  p->distanceBeforeEnding = -1;
  p->typeThickness = _26THICKNESS_;
  p->typeEndPoint = _SURFACE_;
  p->theMask = (typeChamferMask*)NULL;
}







/**************************************************
 *
 * thinning tools
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




typedef struct {
  unsigned char *resBuf;
  int *theDim;
  int *neighb;
  int *offset1D;
  int offset3D[3][3][3];
} typeExtractNeighborhoodParam;



typedef void (*typeExtractNeighborhoodFunction)( typePoint *,
						 typeExtractNeighborhoodParam * );



static void _extract2DNeighborhood( typePoint *pt,
				    typeExtractNeighborhoodParam *p )
{
  int n, i, j;
  
  if ( pt->isinside ) {
    for ( n=0; n<9; n++ )
      p->neighb[n] = p->resBuf[ pt->i + p->offset1D[9+n] ];
  } else {
    for ( n=0, j=-1; j<=1; j++ )
    for ( i=-1; i<=1; i++, n++ ) {
      if ( (pt->y+j < 0) || (pt->y+j>=p->theDim[1]) ||
	   (pt->x+i < 0) || (pt->x+i>=p->theDim[0]) ) 
	p->neighb[n] = 0;
      else 
	p->neighb[n] = p->resBuf[pt->i + p->offset3D[1][1+j][1+i] ];
    }
  }
}



#ifdef _UNUSED_
static void _extract2DNeighborhoodDistance( typePoint *pt,
					    typeExtractNeighborhoodParam *p )
{
  int n, i, j;
  
  if ( pt->isinside ) {
    for ( n=0; n<9; n++ )
      p->neighb[n] = p->resBuf[ pt->i + p->offset1D[9+n] ];
  } else {
    for ( n=0, j=-1; j<=1; j++ )
    for ( i=-1; i<=1; i++, n++ ) {
      if ( (pt->y+j < 0) || (pt->y+j>=p->theDim[1]) ||
	   (pt->x+i < 0) || (pt->x+i>=p->theDim[0]) ) 
	p->neighb[n] = 0;
      else 
	p->neighb[n] = p->resBuf[pt->i + p->offset3D[1][1+j][1+i] ];
    }
  }
}
#endif



static void _extract3DNeighborhood( typePoint *pt,
				    typeExtractNeighborhoodParam *p )
{
  int n, i, j, k;
  
  if ( pt->isinside ) {
    for ( n=0; n<27; n++ )
      p->neighb[n] = p->resBuf[ pt->i + p->offset1D[n] ];
  } else {
    for ( n=0, k=-1; k<=1; k++ )
    for ( j=-1; j<=1; j++ )
    for ( i=-1; i<=1; i++, n++ ) {
      if ( (pt->z+k < 0) || (pt->z+k>=p->theDim[2]) ||
	   (pt->y+j < 0) || (pt->y+j>=p->theDim[1]) ||
	   (pt->x+i < 0) || (pt->x+i>=p->theDim[0]) ) 
	p->neighb[n] = 0;
      else 
	p->neighb[n] = p->resBuf[ pt->i + p->offset3D[1+k][1+j][1+i] ];
    }
  }
}





typedef int (*typeCheckThickness)( int *,
				    enumThickness,
				    int );



static int _check2DThickness( int *neighb, 
			      enumThickness thickness,
			      int direction )
{
  /* neighborhood 
     0 1 2 
     3 4 5
     6 7 8
  */
  /* est-ce une epaisseur ?
     0 = . 0 .   1 = . 1 .   2 = . . .  3 = . . .
         . x .       . x .       1 x 0      0 x 1
         . 1 .       . 0 .       . . .      . . .
  */

  switch ( thickness ) {
  case _04THICKNESS_ :
  case _06THICKNESS_ :
    switch( direction ) {
    default :
    case 0 : 
      if ( neighb[1]  > 0 )  return( 0 );
      if ( neighb[7] == 0 )  return( 0 );
      break;
    case 1 : 
      if ( neighb[7]  > 0 )  return( 0 );
      if ( neighb[1] == 0 )  return( 0 );
      break;
    case 2 : 
      if ( neighb[5]  > 0 )  return( 0 );
      if ( neighb[3] == 0 )  return( 0 );
      break;
    case 3 : 
      if ( neighb[3]  > 0 )  return( 0 );
      if ( neighb[5] == 0 )  return( 0 );
      break;
    }
  case _08THICKNESS_ :
  case _18THICKNESS_ :
  case _26THICKNESS_ :
  default :
    switch( direction ) {
    default :
    case 0 : 
      if ( neighb[1]  > 0 )  return( 0 );
      if ( neighb[6] == 0 && neighb[7] == 0 && neighb[8] == 0 )  return( 0 );
      break;
    case 1 : 
      if ( neighb[7]  > 0 )  return( 0 );
      if ( neighb[0] == 0 && neighb[1] == 0 && neighb[2] == 0 )  return( 0 );
      break;
    case 2 : 
      if ( neighb[5]  > 0 )  return( 0 );
      if ( neighb[0] == 0 && neighb[3] == 0 && neighb[6] == 0 )  return( 0 );
      break;
    case 3 : 
      if ( neighb[3]  > 0 )  return( 0 );
      if ( neighb[2] == 0 && neighb[5] == 0 && neighb[8] == 0 )  return( 0 );
      break;
    }
    break;
  }
  return( 1 );
}



static int _check3DThickness( int *neighb, 
			      enumThickness thickness,
			      int direction )
{
  /* neighborhood 
     0  1  2     9 10 11    18 19 20 
     3  4  5    12 13 14    21 22 23 
     6  7  8    15 16 17    24 25 26
  */
  /* check the thickness ?
   */
  switch ( thickness ) {
  case _04THICKNESS_ :
  case _06THICKNESS_ :
    
    switch( direction ) {
    default :
    case 0 : 
      if ( neighb[ 4]  > 0 )  return( 0 );
      if ( neighb[22] == 0 )  return( 0 );
	      break;
    case 1 : 
      if ( neighb[22]  > 0 )  return( 0 );
      if ( neighb[ 4] == 0 )  return( 0 );
      break;
    case 2 : 
      if ( neighb[10]  > 0 )  return( 0 );
      if ( neighb[16] == 0 )  return( 0 );
      break;
    case 3 : 
      if ( neighb[16]  > 0 )  return( 0 );
      if ( neighb[10] == 0 )  return( 0 );
      break;
    case 4 : 
      if ( neighb[12]  > 0 )  return( 0 );
      if ( neighb[14] == 0 )  return( 0 );
      break;
    case 5 : 
      if ( neighb[14]  > 0 )  return( 0 );
      if ( neighb[12] == 0 )  return( 0 );
      break;
    }
    break;
	    
  case _18THICKNESS_ :
    
    switch( direction ) {
    default :
    case 0 : 
      if ( neighb[ 4]  > 0 )  return( 0 );
      if ( neighb[22] == 0 &&
	   neighb[19] == 0 && neighb[21] == 0 && 
	   neighb[23] == 0 && neighb[25] == 0 )  return( 0 );
      break;
    case 1 : 
      if ( neighb[22]  > 0 )  return( 0 );
      if ( neighb[ 4] == 0 &&
	   neighb[ 1] == 0 && neighb[ 3] == 0 && 
	   neighb[ 5] == 0 && neighb[ 7] == 0 )  return( 0 );
      break;
    case 2 : 
      if ( neighb[10]  > 0 )  return( 0 );
      if ( neighb[16] == 0 &&
	   neighb[ 7] == 0 && neighb[15] == 0 && 
	   neighb[17] == 0 && neighb[25] == 0 )  return( 0 );
      break;
    case 3 : 
      if ( neighb[16]  > 0 )  return( 0 );
      if ( neighb[10] == 0 &&
	   neighb[ 1] == 0 && neighb[ 9] == 0 && 
	   neighb[11] == 0 && neighb[19] == 0 )  return( 0 );
      break;
    case 4 : 
      if ( neighb[12]  > 0 )  return( 0 );
      if ( neighb[14] == 0 &&
	   neighb[ 5] == 0 && neighb[11] == 0 && 
	   neighb[17] == 0 && neighb[23] == 0 )  return( 0 );
      break;
    case 5 : 
      if ( neighb[14]  > 0 )  return( 0 );
      if ( neighb[12] == 0 &&
	   neighb[ 5] == 0 && neighb[ 9] == 0 && 
	   neighb[15] == 0 && neighb[21] == 0 )  return( 0 );
      break;
    }
    break;
    
  case _08THICKNESS_ :
  case _26THICKNESS_ :
  default :
    
    switch( direction ) {
    default :
    case 0 : 
      if ( neighb[ 4]  > 0 )  return( 0 );
      if ( neighb[22] == 0 &&
	   neighb[19] == 0 && neighb[21] == 0 && 
	   neighb[23] == 0 && neighb[25] == 0 &&
	   neighb[18] == 0 && neighb[20] == 0 && 
	   neighb[24] == 0 && neighb[26] == 0 )  return( 0 );
      break;
    case 1 : 
      if ( neighb[22]  > 0 )  return( 0 );
      if ( neighb[ 4] == 0 &&
	   neighb[ 1] == 0 && neighb[ 3] == 0 && 
	   neighb[ 5] == 0 && neighb[ 7] == 0 &&
	   neighb[ 0] == 0 && neighb[ 2] == 0 && 
	   neighb[ 6] == 0 && neighb[ 8] == 0 )  return( 0 );
      break;
    case 2 : 
      if ( neighb[10]  > 0 )  return( 0 );
      if ( neighb[16] == 0 &&
	   neighb[ 7] == 0 && neighb[15] == 0 && 
	   neighb[17] == 0 && neighb[25] == 0 && 
	   neighb[ 6] == 0 && neighb[ 8] == 0 && 
	   neighb[24] == 0 && neighb[26] == 0 )  return( 0 );
      break;
    case 3 : 
      if ( neighb[16]  > 0 )  return( 0 );
      if ( neighb[10] == 0 &&
	   neighb[ 1] == 0 && neighb[ 9] == 0 && 
	   neighb[11] == 0 && neighb[19] == 0 &&
	   neighb[ 0] == 0 && neighb[ 2] == 0 &&
	   neighb[18] == 0 && neighb[20] == 0 )  return( 0 );
      break;
    case 4 : 
      if ( neighb[12]  > 0 )  return( 0 );
      if ( neighb[14] == 0 &&
	   neighb[ 5] == 0 && neighb[11] == 0 && 
	   neighb[17] == 0 && neighb[23] == 0 && 
	   neighb[ 2] == 0 && neighb[ 8] == 0 && 
	   neighb[20] == 0 && neighb[26] == 0 )  return( 0 );
      break;
    case 5 : 
      if ( neighb[14]  > 0 )  return( 0 );
      if ( neighb[12] == 0 &&
	   neighb[ 5] == 0 && neighb[ 9] == 0 && 
	   neighb[15] == 0 && neighb[21] == 0 &&
	   neighb[ 0] == 0 && neighb[ 6] == 0 &&
	   neighb[18] == 0 && neighb[24] == 0 )  return( 0 );
      break;
    }
    break;
    
  }
  return( 1 );
}





typedef int (*typeIsPointSimple)( int *,
				  int *,
				  int * );



static int _is2DPointSimple( int *neighb, int *t04, int *t08 )
{
  int checkT04, checkT08;
  int n;
  
  Compute_T04_and_T08( neighb, t04, t08 );
	  
  if ( *t04 != 1 || *t08 != 1 ) return( 0 );
  
  for ( n=0; n<9; n++ )
    if ( neighb[n] == _TOBEDELETED_ )
      neighb[n] = _BACKGROUND_;
  Compute_T04_and_T08( neighb, &checkT04, &checkT08 );
  if ( checkT04 != 1 || checkT08 != 1 ) return( 0 );

  return( 1 );
}



static int _is3DPointSimple( int *neighb, int *t06, int *t26 )
{
  int checkT06, checkT26;
  int n;
  
  Compute_T06_and_T26( neighb, t06, t26 );
	  
  if ( *t06 != 1 || *t26 != 1 ) return( 0 );

  for ( n=0; n<27; n++ )
    if ( neighb[n] == _TOBEDELETED_ )
      neighb[n] = _BACKGROUND_;
  Compute_T06_and_T26( neighb, &checkT06, &checkT26 );
  if ( checkT06 != 1 || checkT26 != 1 ) return( 0 );

  return( 1 );
}





typedef int (*typeEndConditionSimplePoint)( int * );



static int _defaultEndConditionCurveSimplePoint( int *neighb )
{
  return( 0 );
}

static int _endConditionCurve2DSimplePoint( int *neighb )
{
  int i, n;
  for ( i=0, n=0; n<9; n++ ) {
    if ( neighb[n] == _CANBEDELETED_ || 
	 neighb[n] == _ENDPOINT_ ||
	 neighb[n] == _ANCHOR_ )
      i ++;
  }
  if ( i == 2 ) return( 1 );
  return( 0 );
}

static int _endConditionCurve3DSimplePoint( int *neighb )
{
  int i, n;
  for ( i=0, n=0; n<27; n++ ) {
    if ( neighb[n] == _CANBEDELETED_ || 
	 neighb[n] == _ENDPOINT_ ||
	 neighb[n] == _ANCHOR_ )
      i ++;
  }
  if ( i == 2 ) return( 1 );
  return( 0 );
}




typedef int (*typeEndConditionNonSimplePoint)( int *,
					       int,
					       int,
					       enumTypeEndPoint );



static int _endCondition2DNonSimplePoint( int *neighb, 
					  int t04, int t08,
					  enumTypeEndPoint typeEndPoint )
{
  int i, n;
  
  switch( typeEndPoint ) {
  default :
  case _SURFACE_ :
    if ( t04 >= 2 ) 
      return( 1 );
    break;
  case _PURE_SURFACE_ :
    if ( t04 == 2 ) 
      return( 1 );
    break;
  case _CURVE_ :
    if ( t08 >= 2 ) 
      return( 1 );
    break;
  case _PURE_CURVE_ :
    if ( t08 > 2 ) break;
    for ( i=0, n=0; n<9; n++ ) {
      if ( neighb[n] == _CANBEDELETED_ || 
	   neighb[n] == _ENDPOINT_ ||
	   neighb[n] == _ANCHOR_ )
	i ++;
    }
    if ( i == 3 ) 
      return( 1 );
    break;
  case _NO_END_POINT_ :
    break;
  }
  return( 0 );
}



static int _endCondition3DNonSimplePoint( int *neighb, 
					  int t06, int t26,
					  enumTypeEndPoint typeEndPoint )
{
  int i, n;
  
  switch( typeEndPoint ) {
  default :
  case _SURFACE_ :
    if ( t06 >= 2 ) 
      return( 1 );
    break;
  case _PURE_SURFACE_ :
    if ( t06 == 2 ) 
      return( 1 );
    break;
  case _CURVE_ :
    if ( t26 >= 2 ) 
      return( 1 );
    break;
  case _PURE_CURVE_ :
    if ( t26 > 2 ) break;
    for ( i=0, n=0; n<27; n++ ) {
      if ( neighb[n] == _CANBEDELETED_ || 
	   neighb[n] == _ENDPOINT_ ||
	   neighb[n] == _ANCHOR_ )
	i ++;
    }
    if ( i == 3 ) 
      return( 1 );
    break;
  case _NO_END_POINT_ :
    break;
  }
  return( 0 );
}





/**************************************************
 *
 * thinning procedures
 *
 **************************************************/


int ThinImageWithAllParams( unsigned char *theBuf,
			    unsigned char *resBuf,
			    int *theDim,
			    const int chamfer,
			    typeThinningParameters *p )
{
  char *proc = "ThinImageWithAllParams";
  unsigned short int *theDistance = (unsigned short int *)NULL;
  int v = theDim[0]*theDim[1]*theDim[2];
  int i;
  float exectime[4];
  int tmeas=0;
  int ret;

  exectime[tmeas++] = (float)clock();

  /* distance map computation
     1. initialisation
     2. computation
  */
  theDistance = (unsigned short *)malloc( v * sizeof(unsigned short) );
  if ( theDistance == (unsigned short *)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate distance map\n", proc );
    return( -1 );
  }

  for ( i=0; i<v; i++ ) {
    theDistance[i] = (unsigned short)( ( theBuf[i] > 0 ) ? 0 : 65535 );
  }
  

  if ( _verbose_ >= 1 ) 
    fprintf( stderr, " ... distance compute" );



  switch( chamfer ) {
  default :
  case 3 :
    if ( Compute3DNormalizedChamfer3x3x3( (void*)theDistance, USHORT,
				(void*)theDistance, USHORT,
				theDim ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute distance map\n", proc );
      free( theDistance );
      return( -1 );
    }
    break;
  case 5 :
    if ( Compute3DNormalizedChamfer5x5x5( (void*)theDistance, USHORT,
				(void*)theDistance, USHORT,
				theDim ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to compute distance map\n", proc );
      free( theDistance );
      return( -1 );
    }
    break;
  }



  if ( _verbose_ >= 1 ) 
    fprintf( stderr, "d\n" );




  exectime[tmeas++] = (float)clock();
  if ( _time_verbose_ ) 
    _PrintTimeFrom2Clocks( "total   passed time", exectime[tmeas-2], exectime[tmeas-1] );
 
  ret = ThinImageWithDistanceAndAllParams( theBuf, resBuf, theDistance, theDim, p );

  free( theDistance );

  return( ret );
}










int ThinImageWithDistanceAndAllParams( unsigned char *theBuf,
				       unsigned char *resBuf,
				       unsigned short *theDistance,
				       int *theDim,
				       typeThinningParameters *par )
{
  char *proc = "ThinImageWithDistanceAndAllParams";
  int i, j, k, n;
  int v = theDim[0]*theDim[1]*theDim[2];
  int x, y, z;

  int maxPossibleDistance = 65535;
  int maxFoundDistance = 0;
  int *nbPtsByDistance = (int*)NULL;

  int nbPoints;
  
  typePoint **thePts = (typePoint **)NULL;
  typePoint *tmpPts;
  typePoint tmp;

  int offset3D[3][3][3];
  int offset1D[27];
  int neighb[27];

  int tback, tfore;

  int nbMarkedPts;
  int nbDelPtsDirection, nbEndPtsDirection;
  int nbDelPtsCycle, nbEndPtsCycle;
  int totDelPtsCycle, totEndPtsCycle;
  

  int direction, ndirection=0;
  int cycle, successfuliteration, iteration;
  int minDistance, distance;
  int p;
  int nbPtsDistance;
  
  int d;

  float exectime[4];
  int tmeas = 0;

  
  typeExtractNeighborhoodFunction _extractNeighborhood = (typeExtractNeighborhoodFunction)NULL;
  typeExtractNeighborhoodParam _extractNeighborhoodParam;
  typeCheckThickness _checkThickness = (typeCheckThickness)NULL;
  typeIsPointSimple _isPointSimple = (typeIsPointSimple)NULL;
  typeEndConditionSimplePoint _endConditionSimplePoint = &_defaultEndConditionCurveSimplePoint;
  typeEndConditionNonSimplePoint _endConditionNonSimplePoint = (typeEndConditionNonSimplePoint)NULL;




  /*--------------------------------------------------
   *
   * start
   *
   --------------------------------------------------*/


  exectime[tmeas++] = (float)clock();

  
  
  

  /*--------------------------------------------------
   *
   * build point list
   *
   --------------------------------------------------*/
  

  /* combien de points par valeur de distance ?
   */
  nbPtsByDistance = (int*)malloc( (maxPossibleDistance+1)*sizeof(int) );
  if ( nbPtsByDistance == (int*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate auxiliary array\n", proc );
    return( -1 );
  }
  for ( i=0; i<=maxPossibleDistance; i++ ) nbPtsByDistance[i] = 0;
  
  maxFoundDistance = nbPoints = 0;

  

  /* on ne considere que les points
     > 0 et < par->anchorValue
  */
  if ( par->anchorValue > 0 ) {

    for ( i=0; i<v; i++ ) {
      if ( theBuf[i] == 0 || theBuf[i] >= par->anchorValue ) continue;
      nbPoints ++;
      nbPtsByDistance[ theDistance[i] ] ++;
      if ( maxFoundDistance < theDistance[i] ) maxFoundDistance = theDistance[i];
    }

  } else {

    for ( i=0; i<v; i++ ) {
      if ( theBuf[i] == 0 ) continue;
      nbPoints ++;
      nbPtsByDistance[ theBuf[i] ] ++;
      if ( maxFoundDistance < theDistance[i] ) maxFoundDistance = theDistance[i];
    }

  }
  
  

  if ( nbPoints <= 0 ) {

    if ( resBuf != theBuf ) 
      for ( i=0; i<v; i++ )
	resBuf[i] = ( theBuf[i] > 0 ) ? 255 : 0;

    if ( _verbose_ ) 
      fprintf( stderr, "%s: no points to be thinned\n", proc );
    free( nbPtsByDistance );
    return( 0 );
  }

  

  /* allocation de la liste
   */
  thePts = (typePoint **)malloc( (maxFoundDistance+1)*sizeof(typePoint *) +
			   nbPoints*sizeof(typePoint) );
  if ( thePts == (typePoint **) NULL ) {
     if ( _verbose_ ) 
       fprintf( stderr, "%s: unable to allocate points list\n", proc );
     free( nbPtsByDistance );
     return( -1 );
  }
  tmpPts = (typePoint *)(thePts + (maxFoundDistance+1));
  thePts[0] = (typePoint *)NULL;
  for (i=0; i<=maxFoundDistance; i++ ) {
    if ( nbPtsByDistance[i] > 0 ) {
      thePts[i] = tmpPts;
      tmpPts += nbPtsByDistance[i];
    }
    else {
      thePts[i] = (typePoint *)NULL;
    }
  }

  

  if ( _verbose_ >= 1 ) 
    fprintf( stderr, " ... list construction" );



  /* construction de la liste
   */
  for ( i=0; i<=maxPossibleDistance; i++ ) nbPtsByDistance[i] = 0;

  for ( i=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ )
  for ( x=0; x<theDim[0]; x++, i++ ) {
    if ( theBuf[i] == 0 ) {
      resBuf[i] = 0;
      continue;
    } 
    if ( par->anchorValue > 0  && theBuf[i] >= par->anchorValue ) {
      resBuf[i] = _ANCHOR_;
      continue;
    }

    resBuf[i] = _CANBEDELETED_;

    tmpPts = &(thePts[ theDistance[i] ][ nbPtsByDistance[ theDistance[i] ] ]);
    tmpPts->x = x;
    tmpPts->y = y;
    tmpPts->z = z;
    tmpPts->i = i;
    tmpPts->value = theDistance[i];
    tmpPts->type = _CANBEDELETED_;
    tmpPts->isinside = 1;
    if ( x == 0 ) { tmpPts->isinside = 0; }
    else if ( x == theDim[0] - 1 ) { tmpPts->isinside = 0; }
    else if ( y == 0 )             { tmpPts->isinside = 0; }
    else if ( y == theDim[1] - 1 ) { tmpPts->isinside = 0; }
    else if ( theDim[2] > 1 ) {
      if ( z == 0 ) { tmpPts->isinside = 0; }
      else if ( z == theDim[2] - 1 ) { tmpPts->isinside = 0; }
    }
    nbPtsByDistance[ theDistance[i] ] ++;
  }



  if ( _verbose_ >= 1 ) 
    fprintf( stderr, "\n" );



  exectime[tmeas++] = (float)clock();
  if ( _time_verbose_ ) {
    _PrintTimeFrom2Clocks( "partial passed time", exectime[tmeas-2], exectime[tmeas-1] );
    _PrintTimeFrom2Clocks( "total   passed time", exectime[0], exectime[tmeas-1] );
  }



  if ( _verbose_ ) {
    fprintf( stderr, " ... there are %d points, distances <= %d \n",
	     nbPoints, maxFoundDistance );
  }


  
  
  

  /*--------------------------------------------------
   *
   * list has been built
   * pre-compute offsets for point access speed-up
   *
   --------------------------------------------------*/
 
  /* les offsets [1+dz][1+dy][1+dx]
   */
  offset3D[0][0][0] = -theDim[0]*theDim[1] -theDim[0] -1;
  offset3D[0][0][1] = -theDim[0]*theDim[1] -theDim[0];
  offset3D[0][0][2] = -theDim[0]*theDim[1] -theDim[0] +1;
  offset3D[0][1][0] = -theDim[0]*theDim[1] -1;
  offset3D[0][1][1] = -theDim[0]*theDim[1];
  offset3D[0][1][2] = -theDim[0]*theDim[1] +1;
  offset3D[0][2][0] = -theDim[0]*theDim[1] +theDim[0] -1;
  offset3D[0][2][1] = -theDim[0]*theDim[1] +theDim[0];
  offset3D[0][2][2] = -theDim[0]*theDim[1] +theDim[0] +1;

  offset3D[1][0][0] = -theDim[0] -1;
  offset3D[1][0][1] = -theDim[0];
  offset3D[1][0][2] = -theDim[0] +1;
  offset3D[1][1][0] = -1;
  offset3D[1][1][1] = 0;
  offset3D[1][1][2] = 1;
  offset3D[1][2][0] = theDim[0] -1;
  offset3D[1][2][1] = theDim[0];
  offset3D[1][2][2] = theDim[0] +1;

  offset3D[2][0][0] = theDim[0]*theDim[1] -theDim[0] -1;
  offset3D[2][0][1] = theDim[0]*theDim[1] -theDim[0];
  offset3D[2][0][2] = theDim[0]*theDim[1] -theDim[0] +1;
  offset3D[2][1][0] = theDim[0]*theDim[1] -1;
  offset3D[2][1][1] = theDim[0]*theDim[1];
  offset3D[2][1][2] = theDim[0]*theDim[1] +1;
  offset3D[2][2][0] = theDim[0]*theDim[1] +theDim[0] -1;
  offset3D[2][2][1] = theDim[0]*theDim[1] +theDim[0];
  offset3D[2][2][2] = theDim[0]*theDim[1] +theDim[0] +1;

  for ( n=0, k=0; k<3; k++ )
  for ( j=0; j<3; j++ )
  for ( i=0; i<3; i++, n++ )
    offset1D[n] = offset3D[k][j][i];
  


  
  
  

  /*--------------------------------------------------
   *
   * 
   *
   --------------------------------------------------*/
 
  if ( 1 ) {
    fprintf( stderr, " ... cyclesBeforeEnding   = %d\n", par->cyclesBeforeEnding );
    fprintf( stderr, " ... distanceBeforeEnding = %d\n", par->distanceBeforeEnding );
    fprintf( stderr, " ... typeEndPoint         = " );
    switch ( par->typeEndPoint ) {
    default : fprintf( stderr, "unknown" );
    case _SURFACE_      : fprintf( stderr, "SURFACE" ); break;
    case _PURE_SURFACE_ : fprintf( stderr, "SURFACE" ); break;
    case _CURVE_        : fprintf( stderr, "CURVE" ); break;
    case _PURE_CURVE_   : fprintf( stderr, "CURVE" ); break;
    case _NO_END_POINT_ : fprintf( stderr, "NO_END_POINT" ); break;
    }
    fprintf( stderr, "\n" );
    fprintf( stderr, " ... typeThickness        = " );
    switch ( par->typeThickness ) {
    default : fprintf( stderr, "unknown" );
    case _04THICKNESS_ : fprintf( stderr, "04THICKNESS" ); break;
    case _08THICKNESS_ : fprintf( stderr, "08THICKNESS" ); break;
    case _06THICKNESS_ : fprintf( stderr, "06THICKNESS" ); break;
    case _18THICKNESS_ : fprintf( stderr, "18THICKNESS" ); break;
    case _26THICKNESS_ : fprintf( stderr, "26THICKNESS" ); break;
    }
    fprintf( stderr, "\n" );
  }





  if ( theDim[2] == 1 ) {
    _extractNeighborhood = &_extract2DNeighborhood;
    _checkThickness = &_check2DThickness;
    _isPointSimple = &_is2DPointSimple;
    _endConditionNonSimplePoint = &_endCondition2DNonSimplePoint;
    ndirection = 4;
    switch ( par->typeEndPoint ) {
    default : 
      break;
    case _SURFACE_ :
    case _CURVE_ :
      _endConditionSimplePoint = &_endConditionCurve2DSimplePoint;
      break;
    }
  }
  else {
    _extractNeighborhood = &_extract3DNeighborhood;
    _checkThickness = &_check3DThickness;
    _isPointSimple = &_is3DPointSimple;
    _endConditionNonSimplePoint = &_endCondition3DNonSimplePoint;
    ndirection = 6;
    switch ( par->typeEndPoint ) {
    default : 
      break;
    case _CURVE_ :
      _endConditionSimplePoint = &_endConditionCurve3DSimplePoint;
      break;
    }
  }
  

  _extractNeighborhoodParam.resBuf = resBuf;
  _extractNeighborhoodParam.theDim = theDim;
  _extractNeighborhoodParam.neighb = neighb;
  _extractNeighborhoodParam.offset1D = offset1D;
  for ( k=0; k<3; k++ )
  for ( j=0; j<3; j++ )
  for ( i=0; i<3; i++ )
    _extractNeighborhoodParam.offset3D[k][j][i] = offset3D[k][j][i];





  /*--------------------------------------------------
   *
   * thinning
   *
   --------------------------------------------------*/


  /* general loop from first possible distance to last possible distance 
  */
  distance = minDistance = 1;
  successfuliteration = iteration = 0;
  do {

    /* advance to first non-empty list
     */
    if ( distance == minDistance ) {
      while ( nbPtsByDistance[minDistance] == 0 ) {
	distance ++;
	minDistance ++;
	if ( minDistance > maxFoundDistance ) break;
      }
    }
    if ( minDistance > maxFoundDistance ) break;
    
    /* advance to next non-empty list
     */
    while ( nbPtsByDistance[distance] == 0 ) {
      distance ++;
      if ( distance > maxFoundDistance ) break;
    }
    if ( distance > maxFoundDistance ) break;

    /* loop on cycle
       cycle = a complete loop on all the directions
       iteration =  cycle
    */
    
    cycle = 0;
    totDelPtsCycle = totEndPtsCycle = 0;

    do { 

      /* loop on directions
       */
      nbDelPtsCycle = nbEndPtsCycle = 0;

      for ( direction = 0; direction < ndirection; direction ++ ) {
	
	/* loop on points for one direction
	 */
	nbDelPtsDirection = nbEndPtsDirection = 0;
	nbMarkedPts = 0;
	nbPtsDistance = nbPtsByDistance[distance];
	
	for ( p = 0; p < nbPtsByDistance[distance]; p ++ ) {

	  /* get the neighborhood
	   */
	  (*_extractNeighborhood)( &(thePts[distance][p]), &_extractNeighborhoodParam );

	  /* check the thickness for the direction
	   */
	  if ( (*_checkThickness)( neighb, par->typeThickness, direction ) == 0 )
	    continue;

	  /* test the simplicity
	     if simple and not endpoint, it is marked for deletion
	   */
	  if ( (*_isPointSimple)( neighb, &tback, &tfore ) == 1 ) {
	    /* if there is an end condition for simple points,
	       it is to be tested here
	    */
	    if ( (*_endConditionSimplePoint)( neighb ) == 0 ) {
	      thePts[distance][p].type = resBuf[ thePts[distance][p].i ] = _TOBEDELETED_;
	      nbMarkedPts ++;
	    }
	    continue;
	  }

	   /* shrinking: no further test
	   */
	  if ( par->typeEndPoint == _NO_END_POINT_ )
	    continue;
	  
	  /* end point condition for non-simple points
	     - check whether it is deep enough (whether the distance is large enough)
	     and whether the number of cycles is sufficient
 	     - if yes, it can be considered as a end point candidate
	     else it will be deleted
	     - to be honest, I have no idea why this condition on the cycle numbers
	   */
	  if ( par->distanceBeforeEnding > 0 && d < par->distanceBeforeEnding )
	    continue;
	  if ( par->cyclesBeforeEnding > 0 && cycle < par->cyclesBeforeEnding )
	    continue;

	  /* check whether the non-simple point satisfies an end point condition
	   */
	  if ( (*_endConditionNonSimplePoint)( neighb, tback, tfore, par->typeEndPoint ) == 1 ) {
	    thePts[distance][p].type = resBuf[ thePts[distance][p].i ] = _ENDPOINT_;
	    nbMarkedPts ++;
	  }
	}
	/* end of loop on points for one direction
	 */

	/* update points and image for one direction
	 */
	if ( nbMarkedPts > 0 ) {
	  for ( p = 0; p < nbPtsByDistance[distance]; p ++ ) {
	    switch( thePts[distance][p].type ) {
	    default :
	    case _BACKGROUND_ :
	    case _ANCHOR_ :
	      fprintf( stderr, "%s: WARNING this case should not occur\n", proc );
	      break;
	    case _CANBEDELETED_ :
	      break;
	    case _TOBEDELETED_ :
	      thePts[distance][p].type = resBuf[ thePts[distance][p].i ] = _BACKGROUND_;
	      nbDelPtsDirection ++;
	    case _ENDPOINT_ :
	      if ( thePts[distance][p].type == _ENDPOINT_ ) nbEndPtsDirection ++;
	      tmp = thePts[distance][nbPtsByDistance[distance]-1];
	      thePts[distance][nbPtsByDistance[distance]-1] = thePts[distance][p];
	      thePts[distance][p] = tmp;
	      nbPtsByDistance[distance] --;
	      p --;
	      break;
	    }
	  }
	  nbDelPtsCycle += nbDelPtsDirection;
	  nbEndPtsCycle += nbEndPtsDirection;
	}
	/* end of update points and image for one direction
	 */

	if ( _verbose_ >= 1 ) {
	  fprintf( stderr, " #%6d", iteration );
	  fprintf( stderr, " Dir=%d/%d Cyc.=%2d[%d] Dist.=%3d[%d]/%d",
		   direction, ndirection, 
		   cycle, par->cyclesBeforeEnding, 
		   distance, par->distanceBeforeEnding,
		   maxFoundDistance );
	  fprintf( stderr, " Pts=%8d/%8d Del=%6d End=%5d",
		   nbPtsByDistance[distance], nbPtsDistance, 
		   nbDelPtsDirection, nbEndPtsDirection );
	  if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );
	  else                  fprintf( stderr, "\r" );
	}

      }
      /* end of loop on direction 
       */

      /* go to next cycle if something happens
       */
      if ( nbDelPtsCycle > 0 || nbEndPtsCycle > 0 ) {
	totDelPtsCycle += nbDelPtsCycle;
	totEndPtsCycle += nbEndPtsCycle;
	cycle ++;
	successfuliteration ++;
      }
      iteration++;
      
    } while ( nbDelPtsCycle > 0 || nbEndPtsCycle > 0 );

    /* if some points have been deleted, go back to minDistance since points with
       smallest distances can be simple now, else increment distance
    */
    if ( totDelPtsCycle > 0 && distance > minDistance ) {
      distance = minDistance;
    }
    else {
      distance ++;
    }
	
  } while ( distance <= maxFoundDistance );
  /* fin de la boucle generale
   */




  if ( _verbose_ == 1 ) fprintf( stderr, "\n" );


  exectime[tmeas++] = (float)clock();
  if ( _time_verbose_ ) {
    _PrintTimeFrom2Clocks( "partial passed time", exectime[tmeas-2], exectime[tmeas-1] );
    _PrintTimeFrom2Clocks( "total   passed time", exectime[0], exectime[tmeas-1] );
  }









  /* release memory
  */
  free( thePts );
  free( nbPtsByDistance );



  /* binarise result image
   */
  for ( i=0; i<v; i++ ) {
    if ( resBuf[i] > 0 ) resBuf[i] = 255;
  }


  return( 1 );
}

			    












void Thinning_verbose()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}

void Thinning_noverbose()
{
  _verbose_ = 0;
}


void Thinning_timeVerbose()
{
  if ( _time_verbose_ <= 0 ) _time_verbose_ = 1;
  else _time_verbose_ ++;
}

void Thinning_timeNoverbose()
{
  _time_verbose_ = 0;
}






static void _ComputeTimeFrom2Clocks( float c1, float c2,
				     int *hours, int *minutes, float *seconds )
{
  double d = ( (double)c2 / (double)CLOCKS_PER_SEC ) - 
    ( (double)c1 / (double)CLOCKS_PER_SEC );
  *hours = *minutes = 0;
  *seconds = 0.0;

  if ( d > 3600 ) {
    *hours = (int)(d / 3600);
    d -= *hours * 3600.0;
  }
  if ( d > 60 ) {
    *minutes = (int)(d / 60);
    d -= *minutes * 60.0;
  }
  *seconds = d;
}





static void _PrintTimeFrom2Clocks( char *str, float c1, float c2 )
{
  int h, m;
  float s;
  int l, b=10;
  char format[10] = "%s";
  
  _ComputeTimeFrom2Clocks( c1, c2, &h, &m, &s );

  if ( str != (char*)NULL ) {
    l = strlen(str);
    if ( l % b == 0 ) sprintf(format, "%%%ds", l / b );
    else              sprintf(format, "%%%ds", 1+(l / b) );
    fprintf(stderr, format, str);
  } else {
    fprintf(stderr, "          " );
  }
  if ( h > 0 ) {
    fprintf(stderr, " %2d h", h);
  } else {
    fprintf(stderr, "     ");
  }
  if ( m > 0 ) {
    fprintf(stderr, " %2d mn", m);
  } else {
    fprintf(stderr, "      ");
  }
  fprintf(stderr, " %9.6f s\n", s);
}
