/*************************************************************************
 * vt_jointhisto.c -
 *
 * $Id: vt_jointhisto.c,v 1.1 2000/03/23 08:46:46 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 22:25:15 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */


#include <vt_jointhisto.h>


static int _verbose_ = 0;



typedef struct {
  int x, y, z;
  int i;
  double c;
} typeWeight;


static double _minsigma_ = 0.01;
static double _pourcentage_ = 0.95;
static double _minvalue_ = 1e-6;
static int _rmax_ = 30;
static double _rstep_ = 0.5;




int ComputeJointHistoWithTrsfAndMask( vt_image *imageHmpao,
				      vt_image *imageXenon,
				      vt_image *maskHmpao,
				      vt_image *maskXenon,
				      vt_image *histo,
				      double *mat,
				      float sigmaHmpao,
				      float sigmaXenon,
				      int minHmpao,
				      int minXenon )
{
  char *proc = "ComputeJointHistoWithTrsfAndMask";

  int i, j, k;
  double x, y, z;
  int ix, iy, iz;
  double dx, dy, dz;

  int hmpao, xenon;
  float ***theHist = (float***)NULL;
  

  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;

  float radVoxel;
  double sum, sigma, sigmaSquare;
  double last, c, d;
  typeWeight *theWeights = (typeWeight *)NULL;
  int nmax, nc, n;

  double r;
  vt_image imMask;
  float ***theMask = (float***)NULL;


  if ( histo->type != FLOAT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such histogram image type\n", proc );
    return( 0 );
  }
  theHist = (float***)histo->array;


  if ( maskHmpao != NULL && maskHmpao->type != imageHmpao->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: mask image should be of the same type as the whole image\n", proc );
    return( 0 );
  }
  if ( maskXenon != NULL && maskXenon->type != imageXenon->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: mask image should be of the same type as the whole image\n", proc );
    return( 0 );
  }

  
  for ( hmpao=0; hmpao<histo->dim.y; hmpao++ )
  for ( xenon=0; xenon<histo->dim.x; xenon++ )
    theHist[0][hmpao][xenon] = 0;



  
  if ( sigmaXenon < _minsigma_ && sigmaHmpao < _minsigma_ ) {

    fprintf( stderr, "masque non pris en compte ...\n" );

    switch ( imageHmpao->type ) {
    case USHORT :
      {
	u16 *** theHmpao = (u16 ***)imageHmpao->array;
	
	switch ( imageXenon->type ) {
	case USHORT :
	  {
	    u16 *** theXenon = (u16 ***)imageXenon->array;
	    	
	    for ( k=0; k<imageXenon->dim.z; k++ ) {
	      fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, k, imageXenon->dim.z-1 );
	      for ( j=0; j<imageXenon->dim.y; j++ )
	      for ( i=0; i<imageXenon->dim.x; i++ ) {

		xenon = theXenon[k][j][i];
		if ( xenon < minXenon ) continue;
		if ( xenon >= histo->dim.x ) continue;
		
		x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
		ix = (int)x;
		if ( x < 0.0 || ix >= imageHmpao->dim.x-1 ) continue;
		
		y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
		iy = (int)y;
		if ( y < 0.0 || iy >= imageHmpao->dim.y-1 ) continue;
		
		z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
		iz = (int)z;
		if ( z < 0.0 || iz >= imageHmpao->dim.z-1 ) continue;
		
		dx = x - ix;
		dy = y - iy;
		dz = z - iz;
		

		/*
		  v = 0;
		  v += (1.0-dz) * (1.0-dy) * (1.0-dx) * (double)theHmpao[iz][iy][ix];
		  v += (1.0-dz) * (1.0-dy) *      dx  * (double)theHmpao[iz][iy][ix+1];
		  v += (1.0-dz) *      dy  * (1.0-dx) * (double)theHmpao[iz][iy+1][ix];
		  v += (1.0-dz) *      dy  *      dx  * (double)theHmpao[iz][iy+1][ix+1];
		  v +=      dz  * (1.0-dy) * (1.0-dx) * (double)theHmpao[iz+1][iy][ix];
		  v +=      dz  * (1.0-dy) *      dx  * (double)theHmpao[iz+1][iy][ix+1];
		  v +=      dz  *      dy  * (1.0-dx) * (double)theHmpao[iz+1][iy+1][ix];
		  v +=      dz  *      dy  *      dx  * (double)theHmpao[iz+1][iy+1][ix+1];
		*/
		
		if ( theHmpao[iz][iy][ix] >= minHmpao &&
		     theHmpao[iz][iy][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy][ix] ][xenon]     += (1.0-dz)*(1.0-dy)*(1.0-dx);
		if ( theHmpao[iz][iy][ix+1] >= minHmpao &&
		     theHmpao[iz][iy][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy][ix+1] ][xenon]   += (1.0-dz)*(1.0-dy)*     dx ;
		if ( theHmpao[iz][iy+1][ix] >= minHmpao &&
		     theHmpao[iz][iy+1][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy+1][ix] ][xenon]   += (1.0-dz)*     dy *(1.0-dx);
		if ( theHmpao[iz][iy+1][ix+1] >= minHmpao &&
		     theHmpao[iz][iy+1][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy+1][ix+1] ][xenon] += (1.0-dz)*     dy *     dx ;
		if ( theHmpao[iz+1][iy][ix] >= minHmpao &&
		     theHmpao[iz+1][iy][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy][ix] ][xenon]   +=      dz *(1.0-dy)*(1.0-dx);
		if ( theHmpao[iz+1][iy][ix+1] >= minHmpao &&
		     theHmpao[iz+1][iy][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy][ix+1] ][xenon] +=      dz *(1.0-dy)*     dx ;
		if ( theHmpao[iz+1][iy+1][ix] >= minHmpao &&
		     theHmpao[iz+1][iy+1][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy+1][ix] ][xenon]   +=    dz *     dy *(1.0-dx);
		if ( theHmpao[iz+1][iy+1][ix+1] >= minHmpao &&
		     theHmpao[iz+1][iy+1][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy+1][ix+1] ][xenon] +=    dz *     dy *     dx ;
		
	      }
	    }
	  }
	  break;
	default :
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	  return( 0 );
	}

      }
      break;
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
      return( 0 );
    }
    return( 1 );
  }






  /* ici sigmaXenon >= _minsigma_ || sigmaHmpao >= _minsigma_
   */
  VT_Image( &imMask );
  VT_InitImage( &imMask, "coefficients.inr", 2*_rmax_+1, 2*_rmax_+1, 2*_rmax_+1, FLOAT );
  (void)VT_AllocImage( &imMask );
  theMask = (float***)imMask.array;

  /* sigma equivalent au carre 
     en mm
     puis en pixels dans l'image HMPAO
     
     => on se ramene donc dans l'image HMPAO
  */
  sigmaSquare = sigmaHmpao*sigmaHmpao + sigmaXenon*sigmaXenon;
  sigmaSquare /= (imageHmpao->siz.x * imageHmpao->siz.x);
  sigma = sqrt( sigmaSquare );
  c = sqrt( 2.0 * 3.1415926536 ) * sigma;
  c = c*c*c;

  
  for ( k=0; k<imMask.dim.z; k++ )
  for ( j=0; j<imMask.dim.y; j++ )
  for ( i=0; i<imMask.dim.x; i++ ) {
    d = (k-_rmax_)*(k-_rmax_) + (j-_rmax_)*(j-_rmax_) + (i-_rmax_)*(i-_rmax_);
    theMask[k][j][i] = exp( - (double)d / (2.0 * sigmaSquare) ) / c;
  }
  
  /* (void)VT_WriteInrimage( &imMask ); */
  
  r = 0.0;
  do { 
    last = sum = 0.0;
    r += _rstep_;
    n = 0;

    
    for ( k=0; k<imMask.dim.z; k++ )
    for ( j=0; j<imMask.dim.y; j++ )
    for ( i=0; i<imMask.dim.x; i++ ) {
      d = (k-_rmax_)*(k-_rmax_) + (j-_rmax_)*(j-_rmax_) + (i-_rmax_)*(i-_rmax_);
      if ( d > r*r ) continue;
      sum += theMask[k][j][i];
      last = theMask[k][j][i];
      n ++;
    }
    
    if ( _verbose_ ) {
      fprintf( stderr, " %s: r=%5.2f, sum = %6f, last = %8.6f  #pts=%5d pourcentage = %5.3f\n", 
	       proc, r, sum, last, n, _pourcentage_ );
    }
    
    /* } while ( sum < _pourcentage_ > 0.01 && r < _rmax_ && last > _minvalue_ ); */
  } while ( sum < _pourcentage_ && _pourcentage_> 0.01 && r < _rmax_ && last > _minvalue_ );
  
  VT_FreeImage( &imMask );



  radVoxel = r;
  if ( radVoxel < 1 ) radVoxel = 1.0;

  i = (int)(2*(radVoxel+1.0)+1.5);
  nmax = i*i*i;

  theWeights = (typeWeight *)malloc( nmax * sizeof( typeWeight ) );
  if ( theWeights == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: can not allocate\n", proc );
    }
    return( 0 );
  }
  
  for (i=0; i<nmax; i++) {
    theWeights[i].x = theWeights[i].y = theWeights[i].z = 0;
    theWeights[i].i = 0;
    theWeights[i].c = 0.0;
  }



  if ( 1 ) {
    fprintf( stderr, " %s: rayon = %f pixel (%f mm)\n",
	     proc, radVoxel, radVoxel*imageHmpao->siz.x );
    fprintf( stderr, "\t nb attendus de coeffs = %d\n", n );
    fprintf( stderr, "\t nb max      de coeffs = %d\n", nmax );
    fprintf( stderr, "\t sigmas (hmpao,xenon) = (%f,%f) => %f pixel \n",
	     sigmaHmpao/imageHmpao->siz.x, sigmaXenon/imageXenon->siz.x, sigma );
    fprintf( stderr, "\t                      = (%f,%f) => %f mm \n",
	     sigmaHmpao, sigmaXenon, sigma*imageHmpao->siz.x  );
    fprintf( stderr, "\t coefficient de normalisation = %f pixel\n", 1.0/c );
  }
  





  switch ( imageHmpao->type ) {
  case USHORT :
    {
      u16 *** theHmpao = (u16 ***)imageHmpao->array;
      u16 *** mskHmpao = ( maskHmpao != NULL ) ? (u16 ***)maskHmpao->array : NULL;
	
      switch ( imageXenon->type ) {
      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	  u16 *** mskXenon = ( maskXenon != NULL ) ? (u16 ***)maskXenon->array : NULL;
	    
	  for ( k=0; k<imageXenon->dim.z; k++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, k, imageXenon->dim.z-1 );
	    for ( j=0; j<imageXenon->dim.y; j++ )
	    for ( i=0; i<imageXenon->dim.x; i++ ) {

	      xenon = theXenon[k][j][i];
	      if ( xenon < minXenon ) continue;
	      if ( xenon >= histo->dim.x ) continue;

	      
	      x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
	      ix = (int)x;
	      if ( x < 0.0 || ix >= imageHmpao->dim.x-1 ) continue;
	      
	      y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
	      iy = (int)y;
	      if ( y < 0.0 || iy >= imageHmpao->dim.y-1 ) continue;
	      
	      z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
	      iz = (int)z;
	      if ( z < 0.0 || iz >= imageHmpao->dim.z-1 ) continue;
	  
	      xmin = (int)( x - radVoxel + 0.5 );
	      if ( xmin < 0 ) xmin = 0;
	      xmax = (int)( x + radVoxel + 0.5 );
	      if ( xmax >= imageHmpao->dim.x ) xmax = imageHmpao->dim.x-1;

	      ymin = (int)( y - radVoxel + 0.5 );
	      if ( ymin < 0 ) ymin = 0;
	      ymax = (int)( y + radVoxel + 0.5 );
	      if ( ymax >= imageHmpao->dim.y ) ymax = imageHmpao->dim.y-1;

	      zmin = (int)( z - radVoxel + 0.5 );
	      if ( zmin < 0 ) zmin = 0;
	      zmax = (int)( z + radVoxel + 0.5 );
	      if ( zmax >= imageHmpao->dim.z ) zmax = imageHmpao->dim.z-1;

	      
	      sum = 0;
	      n = nc = 0;
	      for ( iz=zmin; iz<=zmax; iz++ )
	      for ( iy=ymin; iy<=ymax; iy++ )
	      for ( ix=xmin; ix<=xmax; ix++ ) {
		d = (x - (double)ix) * (x - (double)ix) + 
		  (y - (double)iy) * (y - (double)iy) + 
		  (z - (double)iz) * (z - (double)iz);
		if ( d > radVoxel*radVoxel ) continue;
		if ( nc == nmax ) {
		  fprintf( stderr, "%s : FATAL ERROR\n", proc );
		  exit(0);
		}
		/*
		theWeights[nc].x = ix;
		theWeights[nc].y = iy;
		theWeights[nc].z = iz;
		*/

		/* version 1
		   On ne considere un point que si 
		   l'un des deux appartient au cerveau
		*/
		if ( 0 ) {
		  if ( (maskHmpao == NULL && maskXenon == NULL) ||
		       (maskHmpao == NULL && maskXenon != NULL && mskXenon[k][j][i] > 0 ) ||
		       (maskHmpao != NULL && mskHmpao[iz][iy][ix] > 0 && maskXenon == NULL) ||
		       (maskHmpao != NULL && maskXenon != NULL && (mskHmpao[iz][iy][ix] > 0 || mskXenon[k][j][i] > 0)) ) {
		    theWeights[nc].i = theHmpao[iz][iy][ix];
		    theWeights[ nc ].c = exp( - d / (2.0 * sigma) ) / c;
		    sum += theWeights[ nc ].c;
		    nc ++;
		  } 
		  else {
		    sum += exp( - d / (2.0 * sigma) ) / c;
		  }
		}

		/* version 2
		   On ne considere un point que si 
		   le point hmpao appartient au cerveau
		*/
		if ( 1 ) {
		  if ( (maskHmpao == NULL && maskXenon == NULL) ||
		       (maskHmpao == NULL && maskXenon != NULL && mskXenon[k][j][i] > 0 ) ||
		       (maskHmpao != NULL && mskHmpao[iz][iy][ix] > 0 ) ) {
		    theWeights[nc].i = theHmpao[iz][iy][ix];
		    theWeights[ nc ].c = exp( - d / (2.0 * sigma) ) / c;
		    sum += theWeights[ nc ].c;
		    nc ++;
		  } 
		  else {
		    sum += exp( - d / (2.0 * sigma) ) / c;
		  }
		}


	      }

	      
	      /* if ( nc == 0 || sum < 1e-8 ) continue; */
	      if ( nc == 0 ) continue;
	      
	      for ( n=0; n<nc; n++ ) {
		if ( theWeights[n].i < minHmpao ) continue;
		if ( theWeights[n].i >= histo->dim.y ) continue;
		theHist[0][ theWeights[n].i ][ xenon ] += theWeights[n].c/sum;
	      }

	    }
	  }
	}
	break;
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	free( theWeights );
	return( 0 );
      }
      
    }
    break;
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
    free( theWeights );
    return( 0 );
  }



  free( theWeights );

  return( 1 );
}






























int ComputeJointHistoWithTrsf( vt_image *imageHmpao,
			       vt_image *imageXenon,
			       vt_image *histo,
			       double *mat,
			       float sigmaHmpao,
			       float sigmaXenon,
			       int minHmpao,
			       int minXenon )
{
  char *proc = "ComputeJointHistoWithTrsf";

  int i, j, k;
  double x, y, z;
  int ix, iy, iz;
  double dx, dy, dz;

  int hmpao, xenon;
  float ***theHist = (float***)NULL;

  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;

  float radVoxel;
  double sum, sigma, sigmaSquare;
  double last, c, d;
  typeWeight *theWeights = (typeWeight *)NULL;
  int nmax, nc, n;

  double r;
  vt_image imMask;
  float ***theMask = (float***)NULL;


  if ( histo->type != FLOAT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such histogram image type\n", proc );
    return( 0 );
  }
  theHist = (float***)histo->array;


  
  for ( hmpao=0; hmpao<histo->dim.y; hmpao++ )
  for ( xenon=0; xenon<histo->dim.x; xenon++ )
    theHist[0][hmpao][xenon] = 0;



  
  if ( sigmaXenon < _minsigma_ && sigmaHmpao < _minsigma_ ) {
    switch ( imageHmpao->type ) {
    case USHORT :
      {
	u16 *** theHmpao = (u16 ***)imageHmpao->array;
	
	switch ( imageXenon->type ) {
	case USHORT :
	  {
	    u16 *** theXenon = (u16 ***)imageXenon->array;
	    
	    for ( k=0; k<imageXenon->dim.z; k++ ) {
	      fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, k, imageXenon->dim.z-1 );
	      for ( j=0; j<imageXenon->dim.y; j++ )
	      for ( i=0; i<imageXenon->dim.x; i++ ) {

		xenon = theXenon[k][j][i];
		if ( xenon < minXenon ) continue;
		if ( xenon >= histo->dim.x ) continue;
		
		x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
		ix = (int)x;
		if ( x < 0.0 || ix >= imageHmpao->dim.x-1 ) continue;
		
		y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
		iy = (int)y;
		if ( y < 0.0 || iy >= imageHmpao->dim.y-1 ) continue;
		
		z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
		iz = (int)z;
		if ( z < 0.0 || iz >= imageHmpao->dim.z-1 ) continue;
		
		dx = x - ix;
		dy = y - iy;
		dz = z - iz;
		

		/*
		  v = 0;
		  v += (1.0-dz) * (1.0-dy) * (1.0-dx) * (double)theHmpao[iz][iy][ix];
		  v += (1.0-dz) * (1.0-dy) *      dx  * (double)theHmpao[iz][iy][ix+1];
		  v += (1.0-dz) *      dy  * (1.0-dx) * (double)theHmpao[iz][iy+1][ix];
		  v += (1.0-dz) *      dy  *      dx  * (double)theHmpao[iz][iy+1][ix+1];
		  v +=      dz  * (1.0-dy) * (1.0-dx) * (double)theHmpao[iz+1][iy][ix];
		  v +=      dz  * (1.0-dy) *      dx  * (double)theHmpao[iz+1][iy][ix+1];
		  v +=      dz  *      dy  * (1.0-dx) * (double)theHmpao[iz+1][iy+1][ix];
		  v +=      dz  *      dy  *      dx  * (double)theHmpao[iz+1][iy+1][ix+1];
		*/
		
		if ( theHmpao[iz][iy][ix] >= minHmpao &&
		     theHmpao[iz][iy][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy][ix] ][xenon]     += (1.0-dz)*(1.0-dy)*(1.0-dx);
		if ( theHmpao[iz][iy][ix+1] >= minHmpao &&
		     theHmpao[iz][iy][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy][ix+1] ][xenon]   += (1.0-dz)*(1.0-dy)*     dx ;
		if ( theHmpao[iz][iy+1][ix] >= minHmpao &&
		     theHmpao[iz][iy+1][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy+1][ix] ][xenon]   += (1.0-dz)*     dy *(1.0-dx);
		if ( theHmpao[iz][iy+1][ix+1] >= minHmpao &&
		     theHmpao[iz][iy+1][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz][iy+1][ix+1] ][xenon] += (1.0-dz)*     dy *     dx ;
		if ( theHmpao[iz+1][iy][ix] >= minHmpao &&
		     theHmpao[iz+1][iy][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy][ix] ][xenon]   +=      dz *(1.0-dy)*(1.0-dx);
		if ( theHmpao[iz+1][iy][ix+1] >= minHmpao &&
		     theHmpao[iz+1][iy][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy][ix+1] ][xenon] +=      dz *(1.0-dy)*     dx ;
		if ( theHmpao[iz+1][iy+1][ix] >= minHmpao &&
		     theHmpao[iz+1][iy+1][ix] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy+1][ix] ][xenon]   +=    dz *     dy *(1.0-dx);
		if ( theHmpao[iz+1][iy+1][ix+1] >= minHmpao &&
		     theHmpao[iz+1][iy+1][ix+1] < histo->dim.y )
		  theHist[0][ theHmpao[iz+1][iy+1][ix+1] ][xenon] +=    dz *     dy *     dx ;
		
	      }
	    }
	  }
	  break;
	default :
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	  return( 0 );
	}

      }
      break;
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
      return( 0 );
    }
    return( 1 );
  }






  /* ici sigmaXenon >= _minsigma_ || sigmaHmpao >= _minsigma_
   */
  VT_Image( &imMask );
  VT_InitImage( &imMask, "coefficients.inr", 2*_rmax_+1, 2*_rmax_+1, 2*_rmax_+1, FLOAT );
  (void)VT_AllocImage( &imMask );
  theMask = (float***)imMask.array;

  /* sigma equivalent au carre 
     en mm
     puis en pixels dans l'image HMPAO
     
     => on se ramene donc dans l'image HMPAO
  */
  sigmaSquare = sigmaHmpao*sigmaHmpao + sigmaXenon*sigmaXenon;
  sigmaSquare /= (imageHmpao->siz.x * imageHmpao->siz.x);
  sigma = sqrt( sigmaSquare );
  c = sqrt( 2.0 * 3.1415926536 ) * sigma;
  c = c*c*c;

  
  for ( k=0; k<imMask.dim.z; k++ )
  for ( j=0; j<imMask.dim.y; j++ )
  for ( i=0; i<imMask.dim.x; i++ ) {
    d = (k-_rmax_)*(k-_rmax_) + (j-_rmax_)*(j-_rmax_) + (i-_rmax_)*(i-_rmax_);
    theMask[k][j][i] = exp( - (double)d / (2.0 * sigmaSquare) ) / c;
  }
  
  /* (void)VT_WriteInrimage( &imMask ); */
  
  r = 0.0;
  do { 
    last = sum = 0.0;
    r += _rstep_;
    n = 0;

    
    for ( k=0; k<imMask.dim.z; k++ )
    for ( j=0; j<imMask.dim.y; j++ )
    for ( i=0; i<imMask.dim.x; i++ ) {
      d = (k-_rmax_)*(k-_rmax_) + (j-_rmax_)*(j-_rmax_) + (i-_rmax_)*(i-_rmax_);
      if ( d > r*r ) continue;
      sum += theMask[k][j][i];
      last = theMask[k][j][i];
      n ++;
    }
    
    if ( _verbose_ ) {
      fprintf( stderr, " %s: r=%5.2f, sum = %6f, last = %8.6f  #pts=%5d pourcentage = %5.3f\n", 
	       proc, r, sum, last, n, _pourcentage_ );
    }
    
    /* } while ( sum < _pourcentage_ > 0.01 && r < _rmax_ && last > _minvalue_ ); */
  } while ( sum < _pourcentage_ && _pourcentage_> 0.01 && r < _rmax_ && last > _minvalue_ );
  
  VT_FreeImage( &imMask );



  radVoxel = r;
  if ( radVoxel < 1 ) radVoxel = 1.0;

  i = (int)(2*(radVoxel+1.0)+1.5);
  nmax = i*i*i;

  theWeights = (typeWeight *)malloc( nmax * sizeof( typeWeight ) );
  if ( theWeights == NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: can not allocate\n", proc );
    }
    return( 0 );
  }
  
  for (i=0; i<nmax; i++) {
    theWeights[i].x = theWeights[i].y = theWeights[i].z = 0;
    theWeights[i].i = 0;
    theWeights[i].c = 0.0;
  }



  if ( 1 ) {
    fprintf( stderr, " %s: rayon = %f pixel (%f mm)\n",
	     proc, radVoxel, radVoxel*imageHmpao->siz.x );
    fprintf( stderr, "\t nb attendus de coeffs = %d\n", n );
    fprintf( stderr, "\t nb max      de coeffs = %d\n", nmax );
    fprintf( stderr, "\t sigmas (hmpao,xenon) = (%f,%f) => %f pixel \n",
	     sigmaHmpao/imageHmpao->siz.x, sigmaXenon/imageXenon->siz.x, sigma );
    fprintf( stderr, "\t                      = (%f,%f) => %f mm \n",
	     sigmaHmpao, sigmaXenon, sigma*imageHmpao->siz.x  );
    fprintf( stderr, "\t coefficient de normalisation = %f pixel\n", 1.0/c );
  }
  





  switch ( imageHmpao->type ) {
  case USHORT :
    {
      u16 *** theHmpao = (u16 ***)imageHmpao->array;
      
      switch ( imageXenon->type ) {
      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	    
	  for ( k=0; k<imageXenon->dim.z; k++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, k, imageXenon->dim.z-1 );
	    for ( j=0; j<imageXenon->dim.y; j++ )
	    for ( i=0; i<imageXenon->dim.x; i++ ) {

	      xenon = theXenon[k][j][i];
	      if ( xenon < minXenon ) continue;
	      if ( xenon >= histo->dim.x ) continue;

	      
	      x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
	      ix = (int)x;
	      if ( x < 0.0 || ix >= imageHmpao->dim.x-1 ) continue;
	      
	      y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
	      iy = (int)y;
	      if ( y < 0.0 || iy >= imageHmpao->dim.y-1 ) continue;
	      
	      z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
	      iz = (int)z;
	      if ( z < 0.0 || iz >= imageHmpao->dim.z-1 ) continue;
	  
	      xmin = (int)( x - radVoxel + 0.5 );
	      if ( xmin < 0 ) xmin = 0;
	      xmax = (int)( x + radVoxel + 0.5 );
	      if ( xmax >= imageHmpao->dim.x ) xmax = imageHmpao->dim.x-1;

	      ymin = (int)( y - radVoxel + 0.5 );
	      if ( ymin < 0 ) ymin = 0;
	      ymax = (int)( y + radVoxel + 0.5 );
	      if ( ymax >= imageHmpao->dim.y ) ymax = imageHmpao->dim.y-1;

	      zmin = (int)( z - radVoxel + 0.5 );
	      if ( zmin < 0 ) zmin = 0;
	      zmax = (int)( z + radVoxel + 0.5 );
	      if ( zmax >= imageHmpao->dim.z ) zmax = imageHmpao->dim.z-1;

	      
	      sum = 0;
	      n = nc = 0;
	      for ( iz=zmin; iz<=zmax; iz++ )
	      for ( iy=ymin; iy<=ymax; iy++ )
	      for ( ix=xmin; ix<=xmax; ix++ ) {
		d = (x - (double)ix) * (x - (double)ix) + 
		  (y - (double)iy) * (y - (double)iy) + 
		  (z - (double)iz) * (z - (double)iz);
		if ( d > radVoxel*radVoxel ) continue;
		if ( nc == nmax ) {
		  fprintf( stderr, "%s : FATAL ERROR\n", proc );
		  exit(0);
		}
		/*
		theWeights[nc].x = ix;
		theWeights[nc].y = iy;
		theWeights[nc].z = iz;
		*/
		theWeights[nc].i = theHmpao[iz][iy][ix];
		theWeights[ nc ].c = exp( - d / (2.0 * sigma) ) / c;
		sum += theWeights[ nc ].c;
		nc ++;
	      }

	      
	      /* if ( nc == 0 || sum < 1e-8 ) continue; */
	      if ( nc == 0 ) continue;
	      
	      for ( n=0; n<nc; n++ ) {
		if ( theWeights[n].i < minHmpao ) continue;
		if ( theWeights[n].i >= histo->dim.y ) continue;
		theHist[0][ theWeights[n].i ][ xenon ] += theWeights[n].c/sum;
	      }

	    }
	  }
	}
	break;
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	free( theWeights );
	return( 0 );
      }
      
    }
    break;
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
    free( theWeights );
    return( 0 );
  }



  free( theWeights );

  return( 1 );
}






























int ComputeJointHistoWithoutTrsf( vt_image *imageHmpao,
				  vt_image *imageXenon,
				  vt_image *histo,
				  float sigmaHmpao,
				  float sigmaXenon,
				  int minHmpao )
{
  char *proc = "ComputeJointHistoWithoutTrsf";
  int x, y, z;
  int hmpao, xenon;
  int dimx, dimy, dimz;
  float ***theHist = (float***)NULL;
  
  vt_image imMask;
  int d, i, n;
  double r;
  double c, sigma, sigmaSquare;
  double last, sum;
  float ***theMask = (float***)NULL;
  typeWeight *theWeights = (typeWeight *)NULL;

  dimz = imageHmpao->dim.z;
  dimy = imageHmpao->dim.y;
  dimx = imageHmpao->dim.x;
  if ( dimz > imageXenon->dim.z ) dimz = imageXenon->dim.z;
  if ( dimy > imageXenon->dim.y ) dimy = imageXenon->dim.y;
  if ( dimx > imageXenon->dim.x ) dimx = imageXenon->dim.x;


  if ( histo->type != FLOAT ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such histogram image type\n", proc );
    return( 0 );
  }
  theHist = (float***)histo->array;



  
  for ( hmpao=0; hmpao<histo->dim.y; hmpao++ )
  for ( xenon=0; xenon<histo->dim.x; xenon++ )
    theHist[0][hmpao][xenon] = 0;



  
  if ( sigmaXenon < _minsigma_ && sigmaHmpao < _minsigma_ ) {

    if ( _verbose_ ) {
      fprintf( stderr, "%s: pas de sigmas\n", proc );
    }

    switch ( imageHmpao->type ) {

    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
      return( 0 );

    case UCHAR :
      {
	s8 *** theHmpao = (s8 ***)imageHmpao->array;
	
	switch ( imageXenon->type ) {

	default :
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	  return( 0 );

	case UCHAR :
	  {
	    u8 *** theXenon = (u8 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case USHORT :
	  {
	    u16 *** theXenon = (u16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case SSHORT :
	  {
	    s16 *** theXenon = (s16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	}

      }
      break;

    case USHORT :
      {
	u16 *** theHmpao = (u16 ***)imageHmpao->array;
	
	switch ( imageXenon->type ) {

	default :
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	  return( 0 );

	case UCHAR :
	  {
	    u8 *** theXenon = (u8 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case USHORT :
	  {
	    u16 *** theXenon = (u16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case SSHORT :
	  {
	    s16 *** theXenon = (s16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	}

      }
      break;

    case SSHORT :
      {
	s16 *** theHmpao = (s16 ***)imageHmpao->array;

	switch ( imageXenon->type ) {

	default :
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	  return( 0 );

	case UCHAR :
	  {
	    u8 *** theXenon = (u8 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case USHORT :
	  {
	    u16 *** theXenon = (u16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao < 0 ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	case SSHORT :
	  {
	    s16 *** theXenon = (s16 ***)imageXenon->array;
	    for ( z=0; z<dimz; z++ )
	    for ( y=0; y<dimy; y++ )
	    for ( x=0; x<dimx; x++ ) {
	      xenon = theXenon[z][y][x];
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      if ( xenon >= histo->dim.x ) continue;
	      theHist[0][hmpao][xenon] += 1;
	    }
	  }
	  break;

	}

      }
      break;

    }

    return( 1 );
  }

  /* ici sigmaXenon >= 0.01 || sigmaHmpao >= 0.01
   */
  VT_Image( &imMask );
  VT_InitImage( &imMask, "coefficients.inr", 2*_rmax_+1, 2*_rmax_+1, 2*_rmax_+1, FLOAT );
  (void)VT_AllocImage( &imMask );
  theMask = (float***)imMask.array;

  sigmaSquare = sigmaHmpao*sigmaHmpao + sigmaXenon*sigmaXenon;
  sigma = sqrt( sigmaSquare );
  c = sqrt( 2.0 * 3.1415926536 ) * sigma;
  c = c*c*c;

  for ( z=0; z<imMask.dim.z; z++ )
  for ( y=0; y<imMask.dim.y; y++ )
  for ( x=0; x<imMask.dim.x; x++ ) {
    d = (z-_rmax_)*(z-_rmax_) + (y-_rmax_)*(y-_rmax_) + (x-_rmax_)*(x-_rmax_);
    theMask[z][y][x] = exp( - (double)d / (2.0 * sigmaSquare) ) / c;
  }
  
  /* (void)VT_WriteInrimage( &imMask ); */

  
  r = 0.0;
  do { 
    last = sum = 0.0;
    r += _rstep_;
    n = 0;

    
    for ( z=0; z<imMask.dim.z; z++ )
    for ( y=0; y<imMask.dim.y; y++ )
    for ( x=0; x<imMask.dim.x; x++ ) {
      d = (z-_rmax_)*(z-_rmax_) + (y-_rmax_)*(y-_rmax_) + (x-_rmax_)*(x-_rmax_);
      if ( d > r*r ) continue;
      sum += theMask[z][y][x];
      last = theMask[z][y][x];
      n ++;
    }
    
    if ( 1 ) {
      fprintf( stderr, " %s: r=%5.2f, sum = %6f, #pts=%5d pourcentage = %5.3f\n", 
	       proc, r, sum, n, _pourcentage_ );
    }
    
  } while ( sum < _pourcentage_ && r < _rmax_ && last > _minvalue_ );
  
  
  theWeights = (typeWeight *)malloc( n * sizeof( typeWeight ) );
  if ( theWeights == NULL ) {
    if ( 1 ) {
      fprintf( stderr, "%s: can not allocate\n", proc );
    }
    return( 0 );
  }

  n = 0; 
  for ( z=0; z<imMask.dim.z; z++ )
  for ( y=0; y<imMask.dim.y; y++ )
  for ( x=0; x<imMask.dim.x; x++ ) {
    d = (z-_rmax_)*(z-_rmax_) + (y-_rmax_)*(y-_rmax_) + (x-_rmax_)*(x-_rmax_);
    if ( d > r*r ) continue;
    theWeights[n].x = x-_rmax_;
    theWeights[n].y = y-_rmax_;
    theWeights[n].z = z-_rmax_;
    theWeights[n].c = theMask[z][y][x]/sum;
    n ++;
  }

  VT_FreeImage( &imMask );


  if ( 1 ) {
    fprintf( stderr, " %s: rayon = %f pixel (%f mm)\n",
	     proc, r, r*imageXenon->siz.x );
    fprintf( stderr, "\t nb de coeffs = %d\n", n );
    fprintf( stderr, "\t sigmas (hmpao,xenon) = (%f,%f) => %f pixel \n",
	     sigmaHmpao, sigmaXenon, sigma );
    fprintf( stderr, "\t                      = (%f,%f) => %f mm \n",
	     sigmaHmpao*imageHmpao->siz.x, sigmaXenon*imageXenon->siz.x, sigma*imageHmpao->siz.x );
    fprintf( stderr, "\t coefficient de normalisation = %f pixel\n", 1.0/c );
  }
  



  switch ( imageHmpao->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to deal with such hmpao image type\n", proc );
    free( theWeights );
    return( 0 );

  case UCHAR :
    {
      u8 *** theHmpao = (u8 ***)imageHmpao->array;
      
      switch ( imageXenon->type ) {

      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	free( theWeights );
	return( 0 );

      case UCHAR :
	{
	  u8 *** theXenon = (u8 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case SSHORT :
	{
	  s16 *** theXenon = (s16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      }
      
    }
    break;

  case USHORT :
    {
      u16 *** theHmpao = (u16 ***)imageHmpao->array;
      
      switch ( imageXenon->type ) {

      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	free( theWeights );
	return( 0 );

      case UCHAR :
	{
	  u8 *** theXenon = (u8 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case SSHORT :
	{
	  s16 *** theXenon = (s16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      }
      
    }
    break;

  case SSHORT :
    {
      s16 *** theHmpao = (s16 ***)imageHmpao->array;
      
      switch ( imageXenon->type ) {

      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to deal with such xenon image type\n", proc );
	free( theWeights );
	return( 0 );

      case UCHAR :
	{
	  u8 *** theXenon = (u8 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      case SSHORT :
	{
	  s16 *** theXenon = (s16 ***)imageXenon->array;
	  for ( z=0; z<imageHmpao->dim.z; z++ ) {
	    fprintf( stderr, "%s: processing slice %3d/%lu\r", proc, z, imageHmpao->dim.z-1 );
	    for ( y=0; y<imageHmpao->dim.y; y++ )
	    for ( x=0; x<imageHmpao->dim.x; x++ ) {
	      hmpao = theHmpao[z][y][x];
	      if ( hmpao < minHmpao ) continue;
	      if ( hmpao >= histo->dim.y ) continue;
	      for ( i=0; i<n; i++ ) {
		if ( x+theWeights[i].x < 0 || x+theWeights[i].x >= imageXenon->dim.x )
		  continue;
		if ( y+theWeights[i].y < 0 || y+theWeights[i].y >= imageXenon->dim.y )
		  continue;
		if ( z+theWeights[i].z < 0 || z+theWeights[i].z >= imageXenon->dim.z )
		  continue;
		xenon = theXenon[z+theWeights[i].z][y+theWeights[i].y][x+theWeights[i].x];
		if ( xenon >= histo->dim.x ) continue;
		theHist[0][hmpao][xenon] += theWeights[i].c;
	      }
	    }
	  }
	}
	break;

      }
      
    }
    break;

  }
  


  free( theWeights );
  return( 1 );


}




























