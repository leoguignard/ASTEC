
#include <vt_sliceHisto.h>
#include <levenberg.h>
#include <powell.h>

#include <math.h>
double erf(double x);













































int _GetReferenceSlice( int **theHisto, int dimz, int maxval )
{
  int *sum = (int*)malloc( maxval * sizeof(int) );
  int i, imax, max, jmax, j;
  int partial_sum;

  sum[0] = theHisto[0][0];
  for ( i=1; i<maxval; i++ )
    sum[i] = sum[i-1] + theHisto[0][i];
  partial_sum = (int)( 0.9 * sum[ maxval-1 ] + 0.5 );
  for ( imax=0, i=1; i<maxval; i++ ) {
    if ( sum[i-1] < partial_sum && partial_sum <= sum[i] )
      imax = i;
  }
  max = imax;
  jmax = 0;

  for ( j=1; j<dimz; j++ ) {
    sum[0] = theHisto[j][0];
    for ( i=1; i<maxval; i++ )
      sum[i] = sum[i-1] + theHisto[j][i];
    partial_sum = (int)( 0.9 * sum[ maxval-1 ] + 0.5 );
    for ( imax=0, i=1; i<maxval; i++ ) {
      if ( sum[i-1] < partial_sum && partial_sum <= sum[i] )
	imax = i;
    }
    if ( imax > max ) {
      max = imax;
      jmax = j;
    }
  }

  free( sum );
  return( jmax );
}




double ** _GetPDFFromHisto( int **theHisto, int n, int max )
{
  double **thePDF = NULL;
  double *onePDF = NULL;
  double sum;
  int i, j;
  
  thePDF = (double **)malloc( n * sizeof( double * ) + (max) * n * sizeof(double) );
  onePDF = (double * )(thePDF + n);
  for (i=0; i<n; i++ ) {
    thePDF[i] = onePDF;
    memset( onePDF, 0, max*sizeof( double ) );
    onePDF += max;
  }

  for (i=0; i<n; i++ ) {
    sum = 0;
    for ( j =0; j<max; j++ ) {
      thePDF[i][j] = theHisto[i][j];
      sum += thePDF[i][j];
    }
    /* printf( " sum( histo[%d] ) = %f\n", i, sum ); */
    for ( j =0; j<max; j++ ) thePDF[i][j]/= sum;
  }

  return( thePDF );
}




int ** _GetSlicesHisto( vt_image *theIm, int *max )
{
  int **theHisto = NULL;
  int *oneHisto = NULL;
  
  int dimz = theIm->dim.z;
  int dimy = theIm->dim.y;
  int dimx = theIm->dim.x;
  int z, i;
  int maxval = 256;

  *max = 0;

  switch ( theIm->type ) {
  default :
    return( NULL );
  case UCHAR :
    {
      u8 *buf = theIm->buf;
      maxval = buf[0];
      for (i=1; i<dimx*dimy*dimz;i++ )
	if ( maxval < buf[i] ) maxval = buf[i];
    }
    break;
  case USHORT :
    {
      u16 *buf = theIm->buf;
      maxval = buf[0];
      for (i=1; i<dimx*dimy*dimz;i++ )
	if ( maxval < buf[i] ) maxval = buf[i];
    }
    break;
  }
      

  *max = maxval;
  
  theHisto = (int **)malloc( dimz * sizeof( int * ) + (maxval+1) * dimz * sizeof(int) );
  oneHisto = (int * )(theHisto + dimz);
  for (i=0; i<dimz; i++ ) {
    theHisto[i] = oneHisto;
    memset( oneHisto, 0, (maxval+1)*sizeof( int ) );
    oneHisto += maxval+1;
  }

  switch( theIm->type ) {
  case USHORT :
    {
      u16 *buf = theIm->buf;
      for ( z=0; z<dimz;z++ ) {
	for (i=0; i<dimx*dimy;i++, buf++ )
	  theHisto[z][ (*buf) ] ++;
      }
    }
    break;
  case UCHAR :
    {
      u8 *buf = theIm->buf;
      for ( z=0; z<dimz;z++ ) {
	for (i=0; i<dimx*dimy;i++, buf++ )
	  theHisto[z][ (*buf) ] ++;
      }
    }
    break;
  default :
    return( NULL );
  }
  return( theHisto );
}








int * _GetImageHisto( vt_image *theIm, vt_image *immask, int *max )
{
  int *oneHisto = NULL;
  
  int dimz = theIm->dim.z;
  int dimy = theIm->dim.y;
  int dimx = theIm->dim.x;
  int z, i;
  int maxval = 256;

  unsigned char *mask = NULL;

  if ( immask != NULL && immask->type != UCHAR )
    return( NULL );
  if ( immask != NULL )
    mask = (unsigned char*)immask->buf;

  *max = 0;

  switch ( theIm->type ) {
  default :
    return( NULL );
  case UCHAR :
    {
      u8 *buf = theIm->buf;
      maxval = buf[0];
      for (i=1; i<dimx*dimy*dimz;i++ )
	if ( maxval < buf[i] ) maxval = buf[i];
    }
    break;
  case USHORT :
    {
      u16 *buf = theIm->buf;
      maxval = buf[0];
      for (i=1; i<dimx*dimy*dimz;i++ )
	if ( maxval < buf[i] ) maxval = buf[i];
    }
    break;
  case SSHORT :
    {
      s16 *buf = theIm->buf;
      maxval = buf[0];
      for (i=1; i<dimx*dimy*dimz;i++ )
	if ( maxval < buf[i] ) maxval = buf[i];
      maxval += 32768;
    }
    break;
  }
      

  *max = maxval;
  
  oneHisto = (int *)malloc( (maxval+1) * sizeof(int) );
  memset( oneHisto, 0, (maxval+1)*sizeof( int ) );
  
  switch( theIm->type ) {
  case SSHORT :
    {
      s16 *buf = theIm->buf;
      if ( mask == NULL ) {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    oneHisto[ 32768+(int)(*buf) ] ++;
	}
      }
      else {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    if ( mask[i] > 0 ) oneHisto[ 32768+(int)(*buf) ] ++;
	}
      }
    }
    break;
  case USHORT :
    {
      u16 *buf = theIm->buf;
      if ( mask == NULL ) {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    oneHisto[ (*buf) ] ++;
	}
      }
      else {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    if ( mask[i] > 0 ) oneHisto[ (*buf) ] ++;
	}
      }
    }
    break;
  case UCHAR :
    {
      u8 *buf = theIm->buf;
      if ( mask == NULL ) {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    oneHisto[ (*buf) ] ++;
	}
      }
      else {
	for ( z=0; z<dimz;z++ ) {
	  for (i=0; i<dimx*dimy;i++, buf++ )
	    if ( mask[i] > 0 ) oneHisto[ (*buf) ] ++;
	}
      }
    }
    break;
  default :
    return( NULL );
  }
  return( oneHisto );
}








void _PrintOneHistoForMatlab( int fd, FILE *fp, 
			      double *theCoeff, int *theHisto,
			      int length, int id, 
			      enumIntensityCompensation c ) 
{
  char *proc = "_PrintOneHistoForMatlab";
  int i, m;

  if ( fd == 0 || fp == NULL ) return;

  m = theHisto[0];
  for (i=1; i<length; i++ ) if ( m < theHisto[i] )  m = theHisto[i];
  
  if ( write( fd, theHisto, length * sizeof(int) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "[HISTO%d, NREAD%d] = fread( fid, [%d], 'int%lu' );\n", 
	   id, id, length, 8*sizeof(int) );
  fprintf( fp, "\n" );

  if ( theCoeff == NULL ) {

    fprintf( fp, "plot( [0:%d], HISTO%d );\n", length-1, id );
    fprintf( fp, "axis([0 %d 0 %f]);\n", length-1, (double)m );

  } else {
    switch ( c ) {
    default :
    case _AFFINE_ :
      fprintf( fp, "plot( %f + %f*[0:%d], HISTO%d );\n", 
	       theCoeff[0], theCoeff[1], length-1, id );
      fprintf( fp, "axis([0 %d 0 %f]);\n", 
	       length-1, theCoeff[0] + m*theCoeff[1] );
    }
  }

  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );
}








void _PrintOneDoublePDFForMatlab( int fd, FILE *fp, 
			    double *theCoeff, double *theHisto,
			    int length, int id, 
			    enumIntensityCompensation c ) 
{
  char *proc = "_PrintOneDoublePDFForMatlab";
  int i;
  double m;

  if ( fd == 0 || fp == NULL ) return;

  m = theHisto[0];
  for (i=1; i<length; i++ ) if ( m < theHisto[i] )  m = theHisto[i];
  
  if ( write( fd, theHisto, length * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "[PDF%d, NREAD%d] = fread( fid, [%d], 'float%lu' );\n", 
	   id, id, length, 8*sizeof(double) );
  fprintf( fp, "\n" );

  if ( theCoeff == NULL ) {

    fprintf( fp, "plot( [0:%d], PDF%d );\n", length-1, id );
    fprintf( fp, "axis([0 %d 0 %f]);\n", length-1, (double)m );

  } else {
    switch ( c ) {
    default :
    case _AFFINE_ :
      fprintf( fp, "plot( %f + %f*[0:%d], PDF%d );\n", 
	       theCoeff[0], theCoeff[1], length-1, id );
      fprintf( fp, "axis([0 %d 0 %f]);\n", 
	       length-1, theCoeff[0] + m*theCoeff[1] );
    }
  }

  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );
}








void _Print2DSlicesHistoForMatlab( int fd, FILE *fp, int **theHisto, int dimz, int maxval, int id )
{
  char *proc = "_Print2DSlicesHistoForMatlab";

  if ( fd == 0 || fp == NULL ) return;

  if ( write( fd, theHisto[0], maxval * dimz * sizeof(int) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }

  fprintf( fp, "\n" );
  fprintf( fp, "[HISTOREAD, NREAD] = fread( fid, [%d,%d], 'int%lu' );\n", maxval, dimz, 8*sizeof(int) );
  fprintf( fp, "HISTO%d = HISTOREAD';\n", id );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "surf(0:%d,0:%d,HISTO%d)\n", maxval-1, dimz-1, id );
  fprintf( fp, "shading interp;\n" );
  fprintf( fp, "axis([0 %d 0 %d 0 1000]);\n", maxval-1, dimz-1 );
  fprintf( fp, "caxis([0 250]);\n" );
  fprintf( fp, "view(2);\n" );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );

}

void _Print2DSlicesPDFForMatlab( int fd, FILE *fp, double **theHisto, int dimz, int maxval, int id )
{
  char *proc = "_Print2DSlicesHistoForMatlab";

  if ( write( fd, theHisto[0], maxval * dimz * sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }

  fprintf( fp, "\n" );
  fprintf( fp, "[HISTOREAD, NREAD] = fread( fid, [%d,%d], 'float%lu' );\n", maxval, dimz, 8*sizeof(double) );
  fprintf( fp, "HISTO%d = HISTOREAD';\n", id );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "surf(0:%d,0:%d,HISTO%d)\n", maxval-1, dimz-1, id );
  fprintf( fp, "shading interp;\n" );
  fprintf( fp, "axis([0 %d 0 %d 0 0.02]);\n", maxval-1, dimz-1 );
  fprintf( fp, "caxis([0 0.005]);\n" );
  fprintf( fp, "view(2);\n" );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );

}







void _PrintOneSliceHistoForMatlab( int fd, FILE *fp, 
				   double ** theCoeff, int ** theHisto, 
				   int slice, int length, int id )
{
  
  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  if ( theCoeff == NULL ) {
    fprintf( fp, "plot( [0:%d], HISTO%d(%d,:) );\n", length-1, id,slice );
  }
  else {
    fprintf( fp, "plot( %f+%f*[0:%d], HISTO%d(%d,:) );\n", 
	     theCoeff[slice+1][0], theCoeff[slice+1][1], 
	     length-1, id, slice );
  }
  fprintf( fp, "axis([0 %d 0 1500]);\n", length-1 );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );

}

void _PrintOneSlicePDFForMatlab( int fd, FILE *fp, 
				 double ** theCoeff, double **theHisto, 
				 int slice, int maxval, int id )
{
  
  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  if ( theCoeff == NULL ) {
    fprintf( fp, "plot( [0:%d], HISTO%d(%d,:) );\n", maxval-1, id,slice );
  }
  else {
    fprintf( fp, "plot( %f+%f*[0:%d], HISTO%d(%d,:) );\n", 
	     theCoeff[slice+1][0], theCoeff[slice+1][1], 
	     maxval-1, id, slice );
  }
  fprintf( fp, "axis([0 %d 0 0.02]);\n", maxval-1 );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );

}

void _PrintSlicesHistoForMatlab( int fd, FILE *fp, double ** theCoeff,
				 int **theHisto, int dimz, int maxval, int id )
{
  int i;

  if ( fd == 0 || fp == NULL ) return;

  if ( theCoeff == NULL ) {
    fprintf( fp, "\n" );
    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );
    fprintf( fp, "for i=1:%d\n", dimz );
    fprintf( fp, "    plot( [0:%d], HISTO%d(i,:) );\n", maxval-1, id );
    fprintf( fp, "end\n" );
    fprintf( fp, "axis([0 %d 0 1500]);\n", maxval-1 );
    fprintf( fp, "hold off;\n" );
    fprintf( fp, "\n" );
  }
  else {
    fprintf( fp, "\n" );
    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );
    for ( i=0; i<dimz; i++ ) {
      fprintf( fp, "    plot( %f+%f*[0:%d], HISTO%d(%d,:) );\n", 
	       theCoeff[i][0], theCoeff[i][1], maxval-1, id, i+1 );
    }
    fprintf( fp, "axis([0 %d 0 1500]);\n", maxval-1 );
    fprintf( fp, "hold off;\n" );
    fprintf( fp, "\n" );
  }
  
}


















/* faire l'addition de 2 histograme
   ajouter un champ mesure + erreur (pour tri)
   ajouter un champ ratio (pour robust)
*/

typedef struct {
  int    ind;
  double err;
} typeMes;




typedef struct {

  /* histogramme discret
     pdf discrete associee a l'histogramme
     pdf discrete recalculee (a partir d un autre histogramme)
   */

  int *histo;
  double *dpdf;
  double *pdf;
  
  int min;
  int max;
  double sigma;

  /* REFERENCE : histogramme discret
     pdf discrete associee a l'histogramme
     pdf discrete recalculee (a partir d un autre histogramme)
   */
  int *href;
  double *dref;
  double *pref;

  int minref;
  int maxref;
  double sref;


  

  /* ftn permet de passer de histo a href (href = ftn(histo))
   */

  typeIntensityTrsf ftn;
  typeIntensityTrsf invftn;


} typeAuxParam;
  

void _RAZTypeAuxParam( typeAuxParam *p, int hlength , int rlength )
{
  /*
  memset( p->partialHisto, 0, maxval*sizeof(double) );
  p->histo = NULL;
  memset( p->sum, 0, maxval*sizeof(double) );
  memset( p->mes, 0, maxval*sizeof(typeMes) );
  */

  p->histo = NULL;
  p->href = NULL;

  memset( p->dpdf, 0, hlength*sizeof(double) );
  memset( p->pdf,  0, hlength*sizeof(double) );

  memset( p->dref, 0, rlength*sizeof(double) );
  memset( p->pref, 0, rlength*sizeof(double) );  
  
  p->min = 0;
  p->max = hlength-1;
  p->minref = 0;
  p->maxref = rlength-1;
  
}

int _AllocTypeAuxParam( typeAuxParam *p, int hlength, int rlength, double sigma ) 
{
  /*
  p->partialHisto = (double*)malloc( maxval*sizeof(double) );
  p->sum          = (double*)malloc( maxval*sizeof(double) );
  p->mes          =( typeMes*)malloc( maxval*sizeof(typeMes) );
  */
  /*
  p->partialHisto = NULL;
  p->sum          = NULL;
  p->mes          = NULL;
  */

  p->dpdf = (double*)malloc( hlength*sizeof(double) );
  p->pdf  = (double*)malloc( hlength*sizeof(double) );

  p->dref = (double*)malloc( rlength*sizeof(double) );
  p->pref = (double*)malloc( rlength*sizeof(double) );

  if ( sigma > 0.0 ) {
    p->sigma = sigma;
    p->sref  = sigma;
  }
  else {
    p->sigma = 5.0;
    p->sref  = 5.0;
  }

  _RAZTypeAuxParam( p, hlength, rlength );
  return 1;
}

void _FreeTypeAuxParam( typeAuxParam *p )
{
  /*
  free( p->partialHisto );
  free( p->sum );
  free( p->mes );
  */

  p->histo = NULL;
  p->href = NULL;
  free( p->pdf );
  free( p->pref );

  p->min = 0;
  p->max = 0;
}






double _affine( double x, double *c )
{
  return( c[0] + c[1]*x );
}
double _inv_affine( double x, double *c )
{
  return( (x-c[0])/ c[1]);
}


double _constant( double x, double *c )
{
  return( c[0] + x );
}
double _inv_constant( double x, double *c )
{
  return( (x-c[0]) );
}


/*
  faire des sommes a l'infini aux bords
*/

static int _full_computation_ = 0;

void _ComputeDPDFfromHisto( int *h, double *p, 
			    double (*func)(double, double*),
			    double (*invfunc)(double, double*),
			    double *coeff, double sigma,
			    int hmin, int hmax, int pmin, int pmax )
{
  int i, j;
  double sum = 0;
  double c = 1.0 / ( sigma * sqrt( 2.0 ) );
  double d;
  int imin, imax;
  
  for ( i=pmin; i<=pmax; i++ ) p[i] = 0.0;
  for ( j=hmin; j<=hmax; j++ ) sum += h[j];
  d = 0.5 / sum;


  if ( _full_computation_ || (func != NULL && invfunc == NULL) ) {

    if ( func != NULL ) {

      for ( i=pmin; i<=pmax; i++ ) {
	for ( j=hmin; j<=hmax; j++ ) {
	  p[i] += h[j] * d * ( erf( c * ((*func)( (double)i+0.5, coeff )-j) )
			       - erf( c * ((*func)( (double)i-0.5, coeff )-j) ) );
	}
      }
    } 
    else {

      for ( i=pmin; i<=pmax; i++ ) {
	for ( j=hmin; j<=hmax; j++ ) {
	  p[i] += h[j] * d * ( erf( c * ((double)i+0.5 -j) )
			       - erf( c * ((double)i-0.5 -j) ) );
	}
      }
    }

  }
  else {

    if ( func != NULL ) {

      for ( j=hmin; j<=hmax; j++ ) {

	if ( h[j] == 0 ) continue;

	imin = (int)( (*invfunc)( (double)j - 6.0/c, coeff ) - 0.5 ) - 1;
	imax = (int)( (*invfunc)( (double)j + 6.0/c, coeff ) + 0.5 ) + 1;

	if ( imin < pmin ) imin = pmin;
	if ( imax > pmax ) imax = pmax;
	for ( i=imin; i<=imax; i++ ) {
	  p[i] += h[j] * d * ( erf( c * ((*func)( (double)i+0.5, coeff )-j) )
			       - erf( c * ((*func)( (double)i-0.5, coeff )-j) ) );
	}
	
      }
    } 
    else {

      for ( j=hmin; j<=hmax; j++ ) {
	
	if ( h[j] == 0 ) continue;

	imin = (int)( (double)j - 6.0/c - 0.5 ) - 1;
	imax = (int)( (double)j + 6.0/c + 0.5 ) + 1;
	if ( imin < pmin ) imin = pmin;
	if ( imax > pmax ) imax = pmax;
	
	for ( i=imin; i<=imax; i++ ) {
	  p[i] += h[j] * d * ( erf( c * ((double)i+0.5 -j) )
			       - erf( c * ((double)i-0.5 -j) ) );
	}
      }
    }

  }

}



/* on calcule la pdf discrete de la ref
   dans le bining de la coupe
   que l'on compare avec la pdf discrete de la coupe
*/
   
double _DirecteSSD( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i;
  double r=0.0;
  _ComputeDPDFfromHisto( auxPar->href, auxPar->pdf, auxPar->ftn,
			 auxPar->invftn, c, auxPar->sref, 
			 auxPar->min, auxPar->max, auxPar->min, auxPar->max );
  for ( i=auxPar->min; i<=auxPar->max; i++ )
    r += (auxPar->pdf[i] - auxPar->dpdf[i]) * (auxPar->pdf[i] - auxPar->dpdf[i]);
  return( r );
}

/* on calcule la pdf discrete de la coupe
   dans le bining de la ref
   que l'on compare avec la pdf discrete de la ref
*/
   
double _InverseSSD( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i;
  double r=0.0;
  _ComputeDPDFfromHisto( auxPar->histo, auxPar->pref, auxPar->invftn,
			 auxPar->ftn, c, auxPar->sigma, 
			 auxPar->min, auxPar->max, 
			 auxPar->minref, auxPar->maxref );
  for ( i=auxPar->minref; i<=auxPar->maxref; i++ )
    r += (auxPar->pref[i] - auxPar->dref[i]) * (auxPar->pref[i] - auxPar->dref[i]);
  return( r );
}


/* on calcule la pdf discrete de la coupe
   dans le bining de la ref
   que l'on compare avec la pdf discrete de la ref
*/
   
double _SymetrieSSD( double *c, void *par )
{
  return( 0.5 * ( _InverseSSD( c, par ) + _DirecteSSD( c, par ) ) );

}





double _DirecteCorrelation( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i;
  double r=0.0;
  _ComputeDPDFfromHisto( auxPar->href, auxPar->pdf, auxPar->ftn,
			 auxPar->invftn, c, auxPar->sref, 
			 auxPar->min, auxPar->max, auxPar->min, auxPar->max );
  for ( i=auxPar->min; i<=auxPar->max; i++ )
    r += (auxPar->pdf[i] * auxPar->dpdf[i]);
  return( -r );
}
   
double _InverseCorrelation( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i;
  double r=0.0;
  _ComputeDPDFfromHisto( auxPar->histo, auxPar->pref, auxPar->invftn,
			 auxPar->ftn, c, auxPar->sigma, 
			 auxPar->min, auxPar->max, 
			 auxPar->minref, auxPar->maxref );
  for ( i=auxPar->minref; i<=auxPar->maxref; i++ )
    r += (auxPar->pref[i] * auxPar->dref[i]);
  return( -r );
}
   
double _SymetrieCorrelation( double *c, void *par )
{
  return( 0.5 * ( _InverseCorrelation( c, par ) + _DirecteCorrelation( c, par ) ) );

}







double _DirecteVraisemblance( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i, n=0;
  double r=0.0;
  double lp, pmin=1.0;
  
  _ComputeDPDFfromHisto( auxPar->href, auxPar->pdf, auxPar->ftn,
			 auxPar->invftn, c, auxPar->sref, 
			 auxPar->min, auxPar->max, auxPar->min, auxPar->max );

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    if ( 0.0 < auxPar->pdf[i] && auxPar->pdf[i] < pmin ) pmin = auxPar->pdf[i];
    n += auxPar->histo[i];
  }
  lp = log( pmin );
  
  for ( i=auxPar->min; i<=auxPar->max; i++ )
    r += auxPar->dpdf[i] * ((auxPar->pdf[i] > 0.0 ) ? log( auxPar->pdf[i] ) : lp );
  return( - n * r );
}
   
double _InverseVraisemblance( double *c, void *par )
{
  typeAuxParam *auxPar = (typeAuxParam *)par;
  int i, n=0;
  double r=0.0;
  double lp, pmin=1.0;
  
  _ComputeDPDFfromHisto( auxPar->histo, auxPar->pref, auxPar->invftn,
			 auxPar->ftn, c, auxPar->sigma, 
			 auxPar->min, auxPar->max, 
			 auxPar->minref, auxPar->maxref );

  for ( i=auxPar->minref; i<=auxPar->maxref; i++ ) {
    if ( 0.0 < auxPar->pref[i] && auxPar->pref[i] < pmin ) pmin = auxPar->pref[i];
    n += auxPar->href[i];
  }
  lp = log( pmin );
  
  for ( i=auxPar->minref; i<=auxPar->maxref; i++ )
    r += auxPar->dref[i] * ((auxPar->pref[i] > 0.0 ) ? log( auxPar->pref[i] ) : lp );
  return(  - n * r );
}
   
double _SymetrieVraisemblance( double *c, void *par )
{
  return( 0.5 * ( _InverseVraisemblance( c, par ) + _DirecteVraisemblance( c, par ) ) );

}

double _MinimumVraisemblance( double *c, void *par )
{
  double dv = _DirecteVraisemblance( c, par );
  double iv = _InverseVraisemblance( c, par );
  return( (dv<iv) ? dv : iv );
}










void _computeMoments( int *theHisto,
		      int length,
		      double *m1,
		      double *m2 )
{
  int i;
  double s0, s1, s2;

  s0 = 0.0;
  s1 = 0.0;
  s2 = 0.0;
  for ( i=0; i<length; i++ ) {
    s0 += theHisto[i];
    s1 += theHisto[i] * i;
    s2 += theHisto[i] * i * i;
  }
  
  if ( s0 > 0.0 ) *m1 = s1 / s0;
  if ( s0 > 1.0 ) *m2 = (s2 / s0 - s1*s1) * s0 / (s0 - 1.0);
}







/* si q() est la meme distribution que p(),
   mais apres transformation par f (ie q( f(x) ) = p(x))
   on a, avec les moments d'ordre 1 (moyennes) et 2 (variances)
   les relations suivantes
   M_1(q) = a * M_1(p) + b
   M_2(q) - (M_1(q))^2 = a^2 * M_2(p) - (M_1(p))^2
   pour ramener q sur p, il faut donc appliquer f^(-1)()
*/
   
void _initAffineTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength )
{
  double m0, m1, n0, n1;
  /* 
     function / f( auxPar.histo ) = auxPar.href;
  */
  _computeMoments( theHisto, hlength, &n0, &n1 );
  _computeMoments( refHisto, rlength, &m0, &m1 );
  /*
    parametres de f^(-1)
   */
  theCoeff[1] = sqrt( n1 / m1 );
  theCoeff[0] = n0 - theCoeff[1] * m0;
  /* parametres de f()
   */
  theCoeff[1] = sqrt( m1 / n1 );
  theCoeff[0] = m0 - theCoeff[1] * n0;
}








void _initConstantTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength )
{
  double m0, m1, n0, n1;
  /* 
     function / f( auxPar.histo ) = auxPar.href;
  */
  _computeMoments( theHisto, hlength, &n0, &n1 );
  _computeMoments( refHisto, rlength, &m0, &m1 );
  /*
    parametres de f^(-1)
   */
  theCoeff[1] = 1.0;
  theCoeff[0] = n0 - m0;
  /* parametres de f()
   */
  theCoeff[1] = 1.0;
  theCoeff[0] = m0 - n0;
}








void _evalTrsfBetweenTwoHisto( double *theCoeff,
			       int *theHisto, 
			       int hlength,
			       int *refHisto,
			       int rlength,
			       typeIntensityTrsf trsf,
			       typeIntensityTrsf invtrsf,
			       int nparam,
			       double (*func)(double *, void *),
			       double sigma )
{
  typeAuxParam auxPar;
  double measure;
  int iter=0;
  
  _AllocTypeAuxParam( &auxPar, hlength, rlength, sigma );

  auxPar.ftn = trsf;
  auxPar.invftn = invtrsf;
  
  auxPar.histo = theHisto;
  auxPar.href  = refHisto;

  /* initialisation
   */
  if ( 0 ) {
    theCoeff[0] = 0.0;
    theCoeff[1] = 1.0;
  }

  /* calcul des pdf discretes "lissees"
     function / f( auxPar.histo ) = auxPar.href;
   */
  _ComputeDPDFfromHisto( auxPar.histo, auxPar.dpdf, NULL, NULL, NULL,
			 auxPar.sigma, 
			 auxPar.min, auxPar.max, 
			 auxPar.min, auxPar.max );

  _ComputeDPDFfromHisto( auxPar.href, auxPar.dref, NULL, NULL, NULL,
			 auxPar.sref, 
			 auxPar.minref, auxPar.maxref, 
			 auxPar.minref, auxPar.maxref );

  if ( 0 ) {
    int fd = 0;
    FILE *fp = NULL;
    int id = 0;
    
    if ( 1 )
      _ComputeDPDFfromHisto( auxPar.href, auxPar.pdf, auxPar.ftn,
			     auxPar.invftn, theCoeff,
			     auxPar.sref, 
			     auxPar.minref, auxPar.maxref, auxPar.min, auxPar.max );
    {
      _ComputeDPDFfromHisto( auxPar.histo, auxPar.pref, auxPar.ftn,
			     auxPar.invftn, theCoeff,
			     auxPar.sref, 
			     auxPar.minref, auxPar.maxref, auxPar.min, auxPar.max );
    }


    fd = creat( "fooeval.raw", S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    fp = fopen( "fooeval.m", "w" );

    fprintf( fp, "fid = fopen('fooeval.raw', 'r' );\n" );
    fprintf( fp, "\n" );
    _PrintOneDoublePDFForMatlab( fd, fp, NULL, auxPar.dpdf, hlength, id++, _AFFINE_ );
    _PrintOneDoublePDFForMatlab( fd, fp, NULL, auxPar.dref, rlength, id++, _AFFINE_ );
    _PrintOneDoublePDFForMatlab( fd, fp, theCoeff, auxPar.dref, rlength, id++, _AFFINE_ );
    _PrintOneDoublePDFForMatlab( fd, fp, theCoeff, auxPar.pdf, hlength, id++, _AFFINE_ );
    fprintf( fp, "\n" );
    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );

    fprintf( fp, "plot( [0:%d], PDF%d, 'b-' );\n", hlength-1, id-4 ); /* theHisto1 */
    fprintf( fp, "plot( [0:%d], PDF%d, 'm-' );\n", rlength-1, id-3 ); /* theHisto2 */
    fprintf( fp, "plot( %f + %f*[0:%d], PDF%d, 'r-' );\n", 
	     theCoeff[0], theCoeff[1], rlength-1, id-3 );
    fprintf( fp, "plot( [0:%d], PDF%d, 'c-' );\n", rlength-1, id-1 ); /* theHisto2 */
    fprintf( fp, "legend('image 1', 'image 2', 'image 2 transformed', 'image 2 transformed (2)' );\n" );
    fprintf( fp, "hold off;\n" );
    fprintf( fp, "\n" );
    fclose( fp );
    close( fd );
  }

  PowellVerboseOFF();

  printf( "      S#0 A=%f B=%f -> ", theCoeff[1], theCoeff[0] );
  Powell( theCoeff, nparam,
	  0.01, 0.001,
	  &iter, &measure, func, &auxPar );
  printf( "S#0 A=%f B=%f\n", theCoeff[1], theCoeff[0] );

  _FreeTypeAuxParam( &auxPar );
}
			       









int *_subsamplehisto( int *theHisto, int length, int *nlength, int scale )
{
  int *oneHisto;
  int i, j;
  int s = 1;


  oneHisto = (int *)malloc( (length) * sizeof(int) );
  memset( oneHisto, 0, (length)*sizeof( int ) );

  for ( i=0; i<scale; i++ )
    s *= 2;
  
  for (i=0; i<length; i++ ) {
    if ( theHisto[i] == 0 ) continue;
    j = (int)(((double)i)/s + (0.5/s - 0.5));
    oneHisto[j] += theHisto[i];
  }

  *nlength = 0;
  for (i=length-1; i>=0 && *nlength == 0; i-- ) {
    if ( oneHisto[i] > 0 ) *nlength = i;
  }

  if ( 0 ) printf( "      (s=%d) length %d -> %d\n", s, length, *nlength );
  return( oneHisto );
}





void _multiScaleEvalTrsfBetweenTwoHisto( double *theCoeff,
					 int *theHisto, 
					 int hlength,
					 int *refHisto,
					 int rlength,
					 typeIntensityTrsf trsf,
					 typeIntensityTrsf invtrsf,
					 int nparam,
					 double (*func)(double *, void *),
					 double sigma,
					 int lscale, int nscale )
{
  int i, s;
  int *theSubHisto = NULL;
  int hslength = 0;
  int *refSubHisto = NULL;
  int rslength = 0;
  double subCoeff[2];

  typeAuxParam auxPar;
  double measure;
  int iter=0;
  
  _AllocTypeAuxParam( &auxPar, hlength, rlength, sigma );

  auxPar.ftn = trsf;
  auxPar.invftn = invtrsf;
  
  auxPar.histo = theHisto;
  auxPar.href  = refHisto;

  /* initialisation
   */
  if ( 0 ) {
    theCoeff[0] = 0.0;
    theCoeff[1] = 1.0;
  }


  
  for (s=1, i=0; i<lscale; i++, s*=2 ) ;

  for (i=0; i<lscale && i<nscale; i++, s/=2 ) {
    printf( "Scale %d (s=%d)\n", lscale-i, s );
    
    theSubHisto = _subsamplehisto( theHisto, hlength, &hslength, lscale-i );
    refSubHisto = _subsamplehisto( refHisto, rlength, &rslength, lscale-i );

    subCoeff[0] = (1-theCoeff[1])*(0.5/(double)s - 0.5) + theCoeff[0]/(double)s;
    subCoeff[1] = theCoeff[1];
    
    printf( "      S#0 A=%f B=%f -> ", theCoeff[1], theCoeff[0] );
    printf( "S#%d A=%f B=%f\n", lscale-i, subCoeff[1], subCoeff[0] );
    
    _AllocTypeAuxParam( &auxPar, hslength, rslength, sigma/(double)s );
    
    auxPar.ftn = trsf;
    auxPar.invftn = invtrsf;
  
    auxPar.histo = theSubHisto;
    auxPar.href  = refSubHisto;

    /* calcul des pdf discretes "lissees"
       function / f( auxPar.histo ) = auxPar.href;
    */
    _ComputeDPDFfromHisto( auxPar.histo, auxPar.dpdf, NULL, NULL, NULL,
			   auxPar.sigma, 
			   auxPar.min, auxPar.max, 
			   auxPar.min, auxPar.max );
    
    _ComputeDPDFfromHisto( auxPar.href, auxPar.dref, NULL, NULL, NULL,
			   auxPar.sref, 
			   auxPar.minref, auxPar.maxref, 
			   auxPar.minref, auxPar.maxref );

    PowellVerboseOFF();

    Powell( subCoeff, nparam,
	    0.01, 0.001,
	    &iter, &measure, func, &auxPar );

    
    theCoeff[0] = (1-subCoeff[1])*(0.5*(double)s - 0.5) + subCoeff[0]*(double)s;
    theCoeff[1] = subCoeff[1];
    
    printf( "      S#0 A=%f B=%f <- ", theCoeff[1], theCoeff[0] );
    printf( "S#%d A=%f B=%f\n", lscale-i, subCoeff[1], subCoeff[0] );

    _FreeTypeAuxParam( &auxPar );
    free( theSubHisto );
    theSubHisto = NULL;
    free( refSubHisto );
    refSubHisto = NULL;
  }

  /* fait-on le calcul a pleine resolution ?
   */
  if ( nscale > lscale ) {

    printf( "Scale 0\n" );

    _AllocTypeAuxParam( &auxPar, hlength, rlength, sigma );
    
    auxPar.ftn = trsf;
    auxPar.invftn = invtrsf;
    
    auxPar.histo = theHisto;
    auxPar.href  = refHisto;

    /* calcul des pdf discretes "lissees"
       function / f( auxPar.histo ) = auxPar.href;
    */
    _ComputeDPDFfromHisto( auxPar.histo, auxPar.dpdf, NULL, NULL, NULL,
			   auxPar.sigma, 
			   auxPar.min, auxPar.max, 
			   auxPar.min, auxPar.max );
    
    _ComputeDPDFfromHisto( auxPar.href, auxPar.dref, NULL, NULL, NULL,
			   auxPar.sref, 
			   auxPar.minref, auxPar.maxref, 
			   auxPar.minref, auxPar.maxref );
    PowellVerboseOFF();
    
    printf( "      S#0 A=%f B=%f -> ", theCoeff[1], theCoeff[0] );
    Powell( theCoeff, nparam,
	    0.01, 0.001,
	    &iter, &measure, func, &auxPar );
    printf( "S#0 A=%f B=%f\n", theCoeff[1], theCoeff[0] );
    
    _FreeTypeAuxParam( &auxPar );
  }
}
			       















void _evalTrsfsIn3DHistos( double **theCoeff,
			   int **theHisto, 
			   int length,
			   int dimz,
			   int iref,
			   typeIntensityTrsf trsf,
			   typeIntensityTrsf invtrsf,
			   double (*func)(double *, void *),
			   double sigma )
{
  typeAuxParam auxPar;
  double measure;
  int iter=0;
  int z, zz, zf, zl, zs, zref;

  _AllocTypeAuxParam( &auxPar, length, length, sigma );

  auxPar.ftn = trsf;
  auxPar.invftn = invtrsf;

  
  if ( 0 ) {
    theCoeff[iref][0] = 0.0;
    theCoeff[iref][1] = 1.0;
  }

  for ( zz = 0; zz < 2; zz ++ ) {
    if ( zz == 0 ) {
      zf = iref + 1;
      zl = dimz;
      zs = 1;
    }
    else if ( zz == 1 ) {
      zf = iref - 1;
      zl = -1;
      zs = -1;
    } 
    else {
      zf = 1;
      zl = dimz;
      zs = 1;
    }

    for ( z = zf ; z != zl; z += zs ) {
      
      if ( zs > 0 ) zref = z - 1;
      else          zref = z + 1;

      if ( 0 ) {
	theCoeff[z][0] = 0.0;
	theCoeff[z][1] = 1.0;
      }
      
      auxPar.histo = theHisto[z];
      auxPar.href  = theHisto[zref];

      _ComputeDPDFfromHisto( auxPar.histo, auxPar.dpdf, NULL, NULL, NULL,
			     auxPar.sigma, 
			     auxPar.min, auxPar.max, 
			     auxPar.min, auxPar.max );
      
      _ComputeDPDFfromHisto( auxPar.href, auxPar.dref, NULL, NULL, NULL,
			     auxPar.sref, 
			     auxPar.minref, auxPar.maxref, 
			     auxPar.minref, auxPar.maxref );
      PowellVerboseOFF();

      printf( "compensation de %d avec %d : j = %f + %f * i",
	      z, zref, theCoeff[z][0], theCoeff[z][1] );

      Powell( theCoeff[z], 2,
	      0.01, 0.001,
	      &iter, &measure, func, &auxPar );
      
      printf( " => j = %f + %f * i\n",
	      theCoeff[z][0], theCoeff[z][1] );
      
      theCoeff[z][0] = theCoeff[z][0]*theCoeff[zref][1] + theCoeff[zref][0];
      theCoeff[z][1] = theCoeff[z][1]*theCoeff[zref][1];

      printf( " => j = %f + %f * i\n",
	      theCoeff[z][0], theCoeff[z][1] );
    }
  }
  
  _FreeTypeAuxParam( &auxPar );


}
			       















double _minimizeOneSlice( double *c, typeAuxParam *a,
			  double (*func)(double *, void *) )
{
  double e, measure;
  int iter=0;
  e = (*func)( c, a );

  Powell( c,
	  2, 
	  0.01, 0.001,
	  &iter, &measure, func, a );
	  
  return( measure );
}


/*
  maxval est intensity_max +1
*/
void _minimizeFunctionWRTReference( double **theCoeff, int iref,
				    int **theHisto, int dimz, int maxval,
				    double (*func)(double *, void *) )
{
  int i;
  typeAuxParam auxPar;
  double e;

  int z, zz, zf, zl, zs, zref;

  _AllocTypeAuxParam( &auxPar, maxval, maxval, 5.0 );

  theCoeff[iref][0] = 0.0;
  theCoeff[iref][1] = 1.0;

  for ( zz = 0; zz < 2; zz ++ ) {
    if ( zz == 0 ) {
      zf = iref + 1;
      zl = dimz;
      zs = 1;
    }
    else if ( zz == 1 ) {
      zf = iref - 1;
      zl = -1;
      zs = -1;
    } 
    else {
      zf = 1;
      zl = dimz;
      zs = 1;
    }

    for ( z = zf ; z != zl; z += zs ) {
      
      if ( zs > 0 ) zref = z - 1;
      else          zref = z + 1;

      theCoeff[z][0] = 0.0;
      theCoeff[z][1] = 1.0;

      for ( i=0; i<maxval; i++ ) {
	/* auxPar.partialHisto[i] = theHisto[zref][i]; */
      }
      
      auxPar.histo = theHisto[z];
      e = _minimizeOneSlice( theCoeff[z], &auxPar, func );

      theCoeff[z][0] = theCoeff[z][0]*theCoeff[zref][1] + theCoeff[zref][0];
      theCoeff[z][1] = theCoeff[z][1]*theCoeff[zref][1];
    }
  }
  
  _FreeTypeAuxParam( &auxPar );
}






void _PrintJointHistoOfTwoSlicesForMatlab( int fd, FILE *fp, vt_image *image,
					   int maxval,
					   int s1, int s2 )
{
  char *proc = "_PrintJointHistoOfTwoSlicesForMatlab";
  int *jh = calloc( (maxval+1)*(maxval+1) , sizeof(int) );
  u8 ***buf = (u8***)image->array;
  int x,y;

  for (y=0;y<image->dim.y; y++ )
  for (x=0;x<image->dim.x; x++ ) {
    if ( buf[s1][y][x] == 0 || buf[s2][y][x] == 0 ) continue;
    jh[ buf[s1][y][x] * (maxval+1) + buf[s2][y][x] ] ++;
  }
  
  if ( write( fd,  jh, (maxval+1)*(maxval+1) * sizeof(int) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  fprintf( fp, "\n" );
  fprintf( fp, "JOINTREAD = fread( fid, [%d,%d], 'int%lu' );\n",
	   (maxval+1), (maxval+1), 8*sizeof(int) );
  fprintf( fp, "JOINT = JOINTREAD';\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "surf([0:%d],[0:%d], JOINT );\n", maxval, maxval );
  fprintf( fp, "shading interp;\n" );
  fprintf( fp, "axis( [0 %d 0 %d 0 100] );\n", maxval, maxval );
  fprintf( fp, "caxis([0 50]);\n" );
  fprintf( fp, "\n" );

  free( jh );
  
}



void _PrintFunctionOfTwoSlicesForMatlab( int fd, FILE *fp, 
					int **theHisto, int dimz, int maxval,
					 double (*func)(double *, void *),
					int s1, int s2 )
{
  char *proc = "_PrintFunctionOfTwoSlicesForMatlab";
  double *c0, *c1, *t;
  double c[2];

  double is0 = -10;
  double ss0 = 0.5;
  int    ns0 = 40;

  double is1 = 0.1;
  double ss1 = 0.05;
  int    ns1 = 40;

  int i,j;

  typeAuxParam auxPar;
  
  _AllocTypeAuxParam( &auxPar, maxval, maxval, 5.0 );


  


  for ( i=0; i<maxval; i++ ) {
    /* auxPar.partialHisto[i] += theHisto[s1][i]; */
  }
  auxPar.histo = theHisto[s2];
  

  c0 = malloc ( (ns0+1)*sizeof( double) );
  c1 = malloc ( (ns1+1)*sizeof( double) );
  t  = malloc ( (ns0+1)*(ns1+1)*sizeof( double) );

  for ( i=0; i<=ns0; i++ ) {
    c0[ i ] = is0 + i * ss0;
  }
  for ( j=0; j<=ns1; j++ ) {
    c1[ j ] = is1 + j * ss1;
  }

  for ( i=0; i<=ns0; i++ )
  for ( j=0; j<=ns1; j++ ) {
    c[0] = c0[ i ];
    c[1] = c1[ j ];
    t [ i * (ns1+1) + j ] = (*func)( c, (void*)&auxPar );
  }

  if ( write( fd, c0, (ns0+1)*sizeof( double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  if ( write( fd, c1, (ns1+1)*sizeof( double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  if ( write( fd, t,  (ns0+1)*(ns1+1)*sizeof( double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  
  fprintf( fp, "C0 = fread( fid, %d, 'float%lu' );\n", ns0+1, 8*sizeof(double) );
  fprintf( fp, "C1 = fread( fid, %d, 'float%lu' );\n", ns1+1, 8*sizeof(double) );
  fprintf( fp, "ENTROREAD = fread( fid, [%d,%d], 'float%lu' );\n",
	   ns1+1, ns0+1, 8*sizeof(double) );
  fprintf( fp, "ENTROPY = ENTROREAD';\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "surf(C1,C0,ENTROPY );\n" );

  free( c0 );
  free( c1 );
  free( t );
  
  _FreeTypeAuxParam( &auxPar );

}
















typeProbabilite * _buildOneNewPDF( double binsize, int maxval, int *nb )
{
  int i, n;
  double x;
  
  typeProbabilite *pdf = NULL;
  n = (double)maxval/binsize;
  pdf = ( typeProbabilite * )malloc( (n+1)*sizeof( typeProbabilite ) );
  
  x = 1.0;
  for ( i=0; i<=n; i ++, x+=binsize ) {
    pdf[i].x = x;
    pdf[i].xmin = x - binsize / 2.0;
    pdf[i].xmax = x + binsize / 2.0;
    pdf[i].p = 0.0;
  }
  *nb = n+1;
  return( pdf );
}




void _newPDFFromPDF( typeProbabilite *pdf, int npdf,
		    double *thePdf, int n, double sigma )
{
  int i, j;
  double c, s, t;
  for ( j=0; j<npdf; j++ ) pdf[j].p = 0.0;

  /* pour appliquer une transformation de l'intensite,
     il faut changer i en f(i)
     dans erf() ou exp()
  */

  if ( 1 ) {
    /* on a une densite de probabilite,
       en ce sens ou la somme de toutes les probabilites
       va faire 1 (ou a peu pres)
    */ 
    c = 1.0 / ( sigma * sqrt( 2.0 ) );
    for ( i = 0; i <n; i++ ) {
      if ( thePdf[i] < 0.0000001 ) continue;
      for ( j=0; j<npdf; j++ ) {
	pdf[j].p += thePdf[i] * 0.5 * ( erf( c * (pdf[j].xmax-i) ) - erf( c * (pdf[j].xmin-i) ) );
	if ( fabs( pdf[j].x - i) < 0.1  ) {
	  printf( "contrib de %d en %f = %f\n", i, pdf[j].x,
		  0.5 * ( erf( c * (pdf[j].xmax-i) ) - erf( c * (pdf[j].xmin-i) ) ) ); 
	}
      }
    }
  }
  else {
    /* ce n'est plus une densite de probabilite,
       on interpole juste la fonction
       a l'aide d'un noyau gaussien
       la somme des contributions d'un point
       a l'ensemble de la nouvelle PDF n'est pas 1
       (et peut etre variable si une transformation affine de 
       l'intensite est applique)
    */
    c = 1.0 / ( sigma * sqrt( 2.0 * 3.1415926536 ) );
    s = 1.0 / ( 2.0 * sigma * sigma );
    for ( i = 0; i <n; i++ ) {
      if ( thePdf[i] < 0.0000001 ) continue;
      for ( j=0; j<npdf; j++ ) {
	t = pdf[j].x - i;
	pdf[j].p += thePdf[i] * c * exp( - t * t  * s );
      }
    }
  }


  for ( j=0; j<npdf && j <n ; j++ ) {
    printf( "[%d] =% f  [%f] = %f \n",
	    j, thePdf[j], pdf[j].x, pdf[j].p );
  }


  if ( 0 ) {
    for ( j=0; j<npdf; j++ ) {
      printf( "%d [%f %f %f] = %f \n", j, pdf[j].xmin,pdf[j].x, pdf[j].xmax, pdf[j].p);
    }
  }

}

void _PrintOnePDFForMatlab( int fd, FILE *fp, typeProbabilite *pdf, int n, int id )
{
  char *proc = "_PrintOnePDFForMatlab";
  double *tmp;
  int i;

  tmp = (double*)malloc( n*sizeof( double ) );

  fprintf( fp, "\n" );
  for ( i=0; i<n; i++ ) tmp[i] = pdf[i].x;
  if ( write( fd, tmp, n*sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  fprintf( fp, "[XPDF%d, XPDFREAD] = fread( fid, %d, 'float%lu' );\n", id, n, 8*sizeof(double) );
  for ( i=0; i<n; i++ ) tmp[i] = pdf[i].p;
  if ( write( fd, tmp, n*sizeof(double) ) == -1 ) {
    fprintf( stderr, "%s: error when writing\n", proc );
  }
  fprintf( fp, "[YPDF%d, YPDFREAD] = fread( fid, %d, 'float%lu' );\n", id, n, 8*sizeof(double) );
  fprintf( fp, "\n" );

  free( tmp );

  fprintf( fp, "\n" );
  fprintf( fp, "figure;\n" );
  fprintf( fp, "hold on;\n" );
  fprintf( fp, "plot( XPDF%d, YPDF%d, 'ro-' );\n", id, id );
  fprintf( fp, "axis([0 %f 0 0.02]);\n", pdf[n-1].xmax );
  fprintf( fp, "hold off;\n" );
  fprintf( fp, "\n" );

}












































typedef struct {
  double val_x;
  double val_y;
  double err;
} typePoint;


static  int errorcompare( const void *pt1, const void *pt2)
{
  if ( ((typePoint *)pt1)->err > ((typePoint *)pt2)->err )
    return (1);
  if  ( ((typePoint *)pt1)->err < ((typePoint *)pt2)->err )
    return (-1);
  return (0);
}

static  double dabs( double val)
{
  if (val >= 0)
    return val;
  else
    return -val;
}

void _ComputeBiasFromJointHisto( vt_image *image, double ** theCoeff,
				 int iref )
{
  typePoint *thePts;
  int zz, zf, zl, zs, zref;
  int i, n, x, y, z;
  double moy_x, moy_y, moy_xx, moy_xy, moy_yy, sum;
  double A, a, b;

  int iteration = 0;
  int n_best;
  double ratio = 0.2;
  double a_err = 0.000001;
  double b_err = 0.000001;
  double old_a, old_b;
  double threshold = 1000;
  
  static int max_iterations = 100;

  thePts = (typePoint *)malloc( image->dim.x * image->dim.y * sizeof( typePoint ) );
  
  theCoeff[iref][0] = 0.0;
  theCoeff[iref][1] = 1.0;

  for ( zz = 0; zz < 2; zz ++ ) {
    if ( zz == 0 ) {
      zf = iref + 1;
      zl = image->dim.z;
      zs = 1;
    }
    else if ( zz == 1 ) {
      zf = iref - 1;
      zl = -1;
      zs = -1;
    } 
    else {
      zf = 1;
      zl = image->dim.z;
      zs = 1;
    }

    for ( z = zf ; z != zl; z += zs ) {

      if ( zs > 0 ) zref = z - 1;
      else          zref = z + 1;

      theCoeff[z][0] = theCoeff[zref][0];
      theCoeff[z][1] = theCoeff[zref][1];

      if ( z == iref ) continue;
      
      for ( i=0; i < image->dim.x * image->dim.y; i++ ) {
	thePts[i].val_x = 0.0;
	thePts[i].val_y = 0.0;
	thePts[i].err   = 0.0;
      }
      
      moy_x  = moy_y  = 0.0;
      moy_xx = moy_xy = moy_yy = 0.0;
      sum = 0.0;
      
      switch (  image->type ) {
      default :
	exit(0);
      case UCHAR :
	{
	  u8 ***buf = (u8***)image->array;
	  for (n=0, y=0; y <image->dim.y; y++ )
	  for (x=0; x <image->dim.x; x++ ) {
	    if ( buf[zref][y][x] == 0 || buf[   z][y][x] == 0 )
	      continue;
	    /* on veut corriger buf[   z]
	       c'est-a-dire qu'on (a,b) /

	       thePts[n].val_y = b + a * thePts[n].val_x;
	    */
	    thePts[n].val_y = theCoeff[zref][0] + theCoeff[zref][1] * buf[zref][y][x];
	    thePts[n].val_x = buf[   z][y][x];
	    moy_x += thePts[n].val_x;
	    moy_y += thePts[n].val_y;
	    moy_xx += thePts[n].val_x * thePts[n].val_x;
	    moy_xy += thePts[n].val_x * thePts[n].val_y;
	    moy_yy += thePts[n].val_y * thePts[n].val_y;
	    n++;
	    sum  += 1.0;
	  }
	}
	break;
      case USHORT :
	{
	  u16 ***buf = (u16***)image->array;
	  for (n=0, y=0; y <image->dim.y; y++ )
	  for (x=0; x <image->dim.x; x++ ) {
	    if ( buf[zref][y][x] == 0 || buf[   z][y][x] == 0 )
	      continue;
	    /* on veut corriger buf[   z]
	       c'est-a-dire qu'on (a,b) /

	       thePts[n].val_y = b + a * thePts[n].val_x;
	    */
	    thePts[n].val_y = theCoeff[zref][0] + theCoeff[zref][1] * buf[zref][y][x];
	    thePts[n].val_x = buf[   z][y][x];
	    moy_x += thePts[n].val_x;
	    moy_y += thePts[n].val_y;
	    moy_xx += thePts[n].val_x * thePts[n].val_x;
	    moy_xy += thePts[n].val_x * thePts[n].val_y;
	    moy_yy += thePts[n].val_y * thePts[n].val_y;
	    n++;
	    sum  += 1.0;
	  }
	}
	break;
      }
      
      moy_x  /= sum;
      moy_y  /= sum;
      moy_xx /= sum;
      moy_xy /= sum;
      moy_yy /= sum;
      
      A = 0.5 * ( (moy_yy - moy_y*moy_y) - (moy_xx - moy_x*moy_x) ) / ( moy_x*moy_y - moy_xy );
      /* Calcul du coeff dir a */
      if (moy_xy - moy_x*moy_y > 0)
	a = -A + sqrt(A*A+1);
      else
	a = -A - sqrt(A*A+1);
      
      /* Calcul de l'ordonnee a l'origine 
	 On veut bien Y = a * X + b
       */
      b = moy_y - a*moy_x;
      
      theCoeff[z][0] = b;
      theCoeff[z][1] = a;

      /* initialisation
       */
      iteration = 0;
      n_best = (int)(n * (1.0 - ratio)) ;
      if ( n_best > n ) n_best = n;
      old_a = a - a_err - 1;
      old_b = b - b_err - 1;

      fprintf(stdout, "#%3d (ref=%3d): Initial affine bias (%5d/%5d points) : y = %9.5f * x + %9.5f \n", z, zref, n_best, n, a, b );



      if ( 0 )
	fprintf(stderr, "      use %d / %d points\n", n_best, n );

      while ( n_best < n && 
	      iteration < max_iterations &&
	      ((dabs(a-old_a) > a_err) || (dabs(b-old_b) > b_err)) ) {

	iteration++;

	old_a = a;
	old_b = b;

	for ( i = 0 ; i < n; i++) {
	  thePts[i].err = 255.0 * abs ( thePts[i].val_y -  ( a * thePts[i].val_x + b ) ) / (1 + b * b );
	}

	qsort( thePts, n, sizeof(typePoint), errorcompare );

	if ( n_best < n )
	  threshold = thePts[n_best+1].err;

	moy_x  = moy_y  = 0.0;
	moy_xx = moy_xy = moy_yy = 0.0;
	sum = 0.0;
	
	if ( 0 ) {
	  for ( i = 0 ; i < n_best; i++) {
	    moy_x += thePts[i].val_x;
	    moy_y += thePts[i].val_y;
	    moy_xx += thePts[i].val_x * thePts[i].val_x;
	    moy_xy += thePts[i].val_x * thePts[i].val_y;
	    moy_yy += thePts[i].val_y * thePts[i].val_y;
	    sum += 1.0;
	  }
	} 
	else {
	  for ( i = 0 ; i < n; i++) {
	    if ( thePts[i].err >= threshold ) continue;
	    moy_x += thePts[i].val_x;
	    moy_y += thePts[i].val_y;
	    moy_xx += thePts[i].val_x * thePts[i].val_x;
	    moy_xy += thePts[i].val_x * thePts[i].val_y;
	    moy_yy += thePts[i].val_y * thePts[i].val_y;
	    sum += 1.0;
	  }
	}
	


	moy_x  /= sum;
	moy_y  /= sum;
	moy_xx /= sum;
	moy_xy /= sum;
	moy_yy /= sum;
	
	A = 0.5 * ( (moy_yy - moy_y*moy_y) - (moy_xx - moy_x*moy_x) ) / ( moy_x*moy_y - moy_xy );
	/* Calcul du coeff dir a */
	if (moy_xy - moy_x*moy_y > 0)
	  a = -A + sqrt(A*A+1);
	else
	  a = -A - sqrt(A*A+1);
	
	/* Calcul de l'ordonnee a l'origine 
	   On veut bien Y = a * X + b
	*/
	b = moy_y - a*moy_x;

	printf("     Iteration #%2d - ", iteration);
	printf(" %d instead %d/%d ", (int)(sum+0.5), n_best, n );
	printf("| errors = (%11.8f %11.8f)", dabs(a-old_a), dabs(b-old_b));
	printf(" -- y = %9.5f * x + %9.5f \r", a, b);
	
      }
      if ( n_best < n ) printf("\n" );
      printf("\n" );

      theCoeff[z][0] = b;
      theCoeff[z][1] = a;
    }

  }
  free( thePts );
}

































/*********************************************************************
OBSOLETE
*/



double _ComputeValueOfNGaussians( double *gauss, int ng, double x )
{
  double r, y;
  int i;
  r = 0;
  for (i =0; i<ng; i++ ) {
    y = x - gauss[3*i+1];
    r += gauss[3*i] * exp ( - y*y/(2.0*gauss[3*i+2]*gauss[3*i+2]) );
  }
  return ( r );
}


void _ComputePeaks( double *gauss, int ng, double *p1, double *p2, int imax )
{
  double  x=0;
  double dx = 3;
  double p, v, n;
  double max[4] = { -1.0, -1.0, -1.0, -1.0 };
  double xmax[4];
  int i = 0;
  int jismax, j, k;

  *p1 = gauss[1];
  *p2 = gauss[4];

  p = _ComputeValueOfNGaussians( gauss, ng, x );
  x += dx;
  v = _ComputeValueOfNGaussians( gauss, ng, x );
  x += dx;

  for ( ; x <= imax; x += dx ) {
    n = _ComputeValueOfNGaussians( gauss, ng, x );
    if ( v > p && v > n ) {
      xmax[i] = x;
      max[i]  = v;
      i++;
    }
    p = v;
    v = n;
  }

  for ( j = 0; j < i; j ++ ) {
    jismax = 1;
    for ( k = 0; k < i; k ++ ) {
      if ( k == j ) continue;
      if ( max[k] > max[j] ) jismax = 0;
    }
    if ( jismax == 1 ) {
      *p1 = xmax[j];
      max[j] = 0.0;
      j = i+1;
    }
  }
  
  for ( j = 0; j < i; j ++ ) {
    jismax = 1;
    for ( k = 0; k < i; k ++ ) {
      if ( k == j ) continue;
      if ( max[k] > max[j] ) jismax = 0;
    }
    if ( jismax == 1 ) {
      *p2 = xmax[j];
      max[j] = 0.0;
      j = i+1;
    }
  }
  

  dx = 0.1;
  x = *p1;
  p = _ComputeValueOfNGaussians( gauss, ng, x-dx );
  v = _ComputeValueOfNGaussians( gauss, ng, x );
  n = _ComputeValueOfNGaussians( gauss, ng, x+dx );
  if ( p > v ) {
    while ( p > v ) {
      x -= dx;
      v = p;
      p = _ComputeValueOfNGaussians( gauss, ng, x-dx );
    }
  } else {
    while ( n > v ) {
      x += dx;
      v = n;
      n = _ComputeValueOfNGaussians( gauss, ng, x+dx );
    }
  }
  *p1 = x;
    
  dx = 0.1;
  x = *p2;
  p = _ComputeValueOfNGaussians( gauss, ng, x-dx );
  v = _ComputeValueOfNGaussians( gauss, ng, x );
  n = _ComputeValueOfNGaussians( gauss, ng, x+dx );
  if ( p > v ) {
    while ( p > v ) {
      x -= dx;
      v = p;
      p = _ComputeValueOfNGaussians( gauss, ng, x-dx );
    }
  } else {
    while ( n > v ) {
      x += dx;
      v = n;
      n = _ComputeValueOfNGaussians( gauss, ng, x+dx );
    }
  }
  *p2 = x;
    
}

void _InitFirstGaussian( double *gauss, int *histo, int imin, int imax )
{
  double t = 750;
  double sp, s;
  int i, n;

  s = sp = 0.0;
  n = 0;
  for ( i=imin; i<=imax; i++ ) {
    if (histo[i] > t ) {
      s += histo[i] * i;
      sp += histo[i];
      n ++;
    }
  }
  if ( sp > 0.0 ) {
    gauss[1] = s / sp;
    if ( histo[(int)(s/sp)] > 0.0 ) gauss[0] = histo[(int)(s/sp)];
    gauss[2] = n/2;
  }

}

void _InitNextGaussian( double *gauss, int ng,
			  double *next, int *histo, int imin, int imax )
{
  int i, j;
  double *buf = malloc( (imax+1) * sizeof( double ) );

  double m, s;
  int f;
  double maxarea, area;

  for (i = 0; i<=imax; i ++ ) {
    buf[i] = 0;
  }
  for (i = imin; i<=imax; i ++ ) {
    buf[i] = histo[i];
    for ( j=0; j<ng; j++ ) {
      buf[i] -= gauss[3*j] * exp( - (i-gauss[3*j+1])*(i-gauss[3*j+1]) /
				  ( 2.0 * gauss[3*j+2] * gauss[3*j+2] ) );
    }
  }
  f = 0;
  m = 0.0;
  s = 1.0;
  maxarea = area = -1.0;
  
  for (i = imin; i<=imax; i ++ ) {
    if ( buf[i] > 0.0  ) {
      if (i == imin || buf[i-1] <= 0.0 ) {
	f = i;
	area = buf[i];
      }
      else {
	area += buf[i];
      }
    }
    else {
      if ( i > imin && buf[i-1] > 0.0 ) {
	if ( area > maxarea ) {
	  maxarea = area;
	  m = (f + i-1)/2.0;
	  s = (f - i+1)/2.0;
	}
      }
    }
  }
  
  if ( m > 0.0 ) {
    next[1] = m;
    if ( buf[(int)m] > 0 ) next[0] = buf[(int)m];
    if ( s > 0.0 ) next[2] = s;
  }

  free( buf );
}









#define NBARGS 15
#define NBGAUSS 5

void _ComputeBiasFromGaussians( FILE *fp,
				double ** theCoeff,
				int **theHisto, int dimz, int maxval,
				int imin, int imax, int iref, int id )

{
  int i, j, z;
  int ng = NBGAUSS;
  double init_gauss[NBARGS] = { 500,  40,  5,   200, 100, 10,
				 0,   50,  5,     0,  75,  5,  0,  75,  5 };
  double gauss[NBARGS];
  double ref_gauss[NBARGS];
  double ref_p1, ref_p2, p1, p2;

  memcpy( ref_gauss, init_gauss, NBARGS * sizeof(double) );
  fprintf( fp,  "\n" ); 
  fprintf( fp,  "\n" ); 
  fprintf( fp,  "\n" ); 
  fprintf( fp,  "X = [0:255];\n" ); 
  fprintf( fp,  "\n" ); 
  fprintf( fp,  "%% reference slice #%3d :\n", iref );

  for ( j=1; j<=ng; j++ ) {
      fprintf( fp,  "%% INIT: [%f %f %f] -> ", 
	       ref_gauss[3*(j-1)], ref_gauss[3*(j-1)+1], ref_gauss[3*(j-1)+2] );
      if ( j == 1 ) {
	_InitFirstGaussian( ref_gauss, theHisto[iref], imin, imax );
      }
      else if ( j >= 2 ) {
	_InitNextGaussian( ref_gauss, j-1, &(ref_gauss[3*(j-1)]), 
			   theHisto[iref], imin, imax );
      }
      fprintf( fp,  "[%f %f %f]\n", 
	      ref_gauss[3*(j-1)], ref_gauss[3*(j-1)+1], ref_gauss[3*(j-1)+2] );

      _FitGaussiansOnSliceHisto( ref_gauss, theHisto, dimz, maxval,
			       imin, imax, iref, j );
     fprintf( fp,  "%%" );
     for ( i=0; i<j; i++ )
       fprintf( fp, " [%f %f %f]", ref_gauss[i*3], ref_gauss[i*3+1], ref_gauss[i*3+2] );
     fprintf( fp,  "\n" ); 

     if ( j == ng ) {
       fprintf( fp,  "figure;\nhold on\n" );
       fprintf( fp,  "plot( X, HISTO%d(%d,:) );\n", id ,iref+1 );
       fprintf( fp,  "title('histogram of slice #%d');\n", iref );
     }
     for ( i=1; i<=j; i++ ) {
       if ( j != ng ) fprintf( fp,  "%% " );
       fprintf( fp,  "Y%d = (X - %f ) .* (X - %f );\n", 
		i, ref_gauss[(i-1)*3+1], ref_gauss[(i-1)*3+1] );
       if ( j != ng ) fprintf( fp,  "%% " );
       fprintf( fp,  "G%d = %f * exp( - Y%d / (2 * %f * %f) );\n", i,
		ref_gauss[(i-1)*3], i, ref_gauss[(i-1)*3+2], ref_gauss[(i-1)*3+2] );
       if ( j != ng ) fprintf( fp,  "%% " );
       fprintf( fp,  "plot( X, G%d, 'g' );\n", i );
     }
     if ( j != ng ) fprintf( fp,  "%% " );
     fprintf( fp,  "plot( X, " );
     for ( i=1; i<=j; i++ ) {
       fprintf( fp,  "G%d", i );
       if ( i < j ) fprintf( fp,  "+" );
     }
     fprintf( fp,  " ,'r' );\n" );
     fprintf( fp,  "\n" );
  }

  _ComputePeaks( ref_gauss, ng, &ref_p1, &ref_p2, imax );
  theCoeff[iref][0] = 0.0;
  theCoeff[iref][1] = 1.0;

  fprintf( fp,  "%% peaks of #%3d : %f %f \n", iref, p1, p2 );
  fprintf( fp,  "\n" );
  fprintf( fp,  "\n" );
  



  for ( z=0; z< dimz; z++ ) {
    if ( z == iref ) continue;
    memcpy( gauss, init_gauss, NBARGS * sizeof(double) );
    fprintf( fp,  "\n" ); 
    fprintf( fp,  "%% slice #%3d :\n", z );

    for ( j=1; j<=ng; j++ ) {
      fprintf( fp,  "%% INIT: [%f %f %f] -> ", 
	       gauss[3*(j-1)], gauss[3*(j-1)+1], gauss[3*(j-1)+2] );
      if ( j == 1 ) {
	_InitFirstGaussian( gauss, theHisto[z], imin, imax );
      }
      else if ( j >= 2 ) {
	_InitNextGaussian( gauss, j-1, &(gauss[3*(j-1)]), 
			   theHisto[z], imin, imax );
      }
      fprintf( fp,  "[%f %f %f]\n", 
	      gauss[3*(j-1)], gauss[3*(j-1)+1], gauss[3*(j-1)+2] );

      _FitGaussiansOnSliceHisto( gauss, theHisto, dimz, maxval,
				 imin, imax, z, j );
      if ( 0 ) {
	if ( gauss[3*(j-1)] < 0.0 ) {
	  gauss[3*(j-1)] = 0.0;
	  _FitGaussiansOnSliceHisto( gauss, theHisto, dimz, maxval,
				     imin, imax, z, j-1 );
	  j = ng+1;
	}
      }

      fprintf( fp,  "%%" );
      for ( i=0; i<j; i++ )
	fprintf( fp, " [%f %f %f]", gauss[i*3], gauss[i*3+1], gauss[i*3+2] );
      fprintf( fp,  "\n" ); 
      
      if ( j == ng ) {
	fprintf( fp,  "figure;\nhold on\n" );
	fprintf( fp,  "plot( X, HISTO%d(%d,:) );\n", id, z+1 );
	fprintf( fp,  "title('histogram of slice #%d');\n", z );
      }
      for ( i=1; i<=j; i++ ) {
	if ( j != ng ) fprintf( fp,  "%% " );
	fprintf( fp,  "Y%d = (X - %f ) .* (X - %f );\n", 
		 i, gauss[(i-1)*3+1], gauss[(i-1)*3+1] );
	if ( j != ng ) fprintf( fp,  "%% " );
	fprintf( fp,  "G%d = %f * exp( - Y%d / (2 * %f * %f) );\n", i,
		 gauss[(i-1)*3], i, gauss[(i-1)*3+2], gauss[(i-1)*3+2] );
	if ( j != ng ) fprintf( fp,  "%% " );
	fprintf( fp,  "plot( X, G%d, 'g' );\n", i );
      }
      if ( j != ng ) fprintf( fp,  "%% " );
      fprintf( fp,  "plot( X, " );
      for ( i=1; i<=j; i++ ) {
	fprintf( fp,  "G%d", i );
	if ( i < j ) fprintf( fp,  "+" );
      }
      fprintf( fp,  " ,'r' );\n" );
      fprintf( fp,  "\n" );
    }
    
    for ( i=0; i< NBARGS; i++ ) {
      if ( gauss[i] < 0.0 ) {
	fprintf( stderr, "WARNING, #%3d gauss[%2d] = %f\n", z, i, gauss[i] );
      }
    }


    _ComputePeaks( gauss, ng, &p1, &p2, imax );
    fprintf( fp,  "%% peaks of #%3d : %f %f \n", z, p1, p2 );
    theCoeff[z][1] = ( ref_p2 - ref_p1 ) / ( p2 - p1 );
    theCoeff[z][0] = ( ref_p2 - theCoeff[z][1] * p2 );
  }
}




void _FitGaussiansOnSliceHisto( double *par,
				   int **theHisto, int dimz, int maxval,
				   int imin, int imax,
				int z, int n )
{
  double *theX, *theY, *theC, *theS;
  int i, j;
  
  theX =  (double*)malloc( 4*(imax-imin+1)*sizeof(double) );
  theY =  theX;
  theY += imax-imin+1;
  theC =  theY;
  theC += imax-imin+1;
  theS =  theC;
  theS += imax-imin+1;

  for ( i=0; i<imax-imin+1; i++ ) 
    theX[i] = theY[i] = theC[i] = theS[i] = 1.0;

  for ( i=0, j=imin; j<=imax; j++, i++ ) {
    theX[i] = j;
    theY[i] = theHisto[z][j];
  }

  /*  for ( i=50; i<100; i++ ) {
    theC[i] = 0.0;
  }
  */
  switch ( n ) {
  case 1 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       imax-imin+1,
				       par, n*3, _GaussianForLM );
    break;
  case 2 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       imax-imin+1,
				       par, n*3, _MixtureOf2GaussiansForLM );
    break;
  case 3 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       imax-imin+1,
				       par, n*3, _MixtureOf3GaussiansForLM );
    break;
  case 4 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       imax-imin+1,
				       par, n*3, _MixtureOf4GaussiansForLM );
    break;
  case 5 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       imax-imin+1,
				       par, n*3, _MixtureOf5GaussiansForLM );
    break;
  }
  free( theX );
}
















double _CorrelationOfOneHisto( double *c, void *par )
{
  /*
  typeAuxParam *auxPar = (typeAuxParam *)par;
  double sum_sum, sum_histo;
  int i, j;
  double s = 0.0;
  double v, e=0.0;
  */

  exit( 0 );
  /*
  sum_sum = 0.0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_sum += auxPar->partialHisto[i];
  }

  sum_histo = 0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_histo += auxPar->histo[i];
  }
  
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    auxPar->sum[i] = 0.0;
  }

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = c[0] + c[1] * i;
    j = (int)v;

    auxPar->mes[i].ind = i;
    
    if ( j < auxPar->min || j >= auxPar->max ) {
      continue;
    }
    e = (j+1-v)*auxPar->partialHisto[j] + (v-j) * auxPar->partialHisto[j+1];
    s += e * auxPar->histo[i] / (sum_sum * sum_histo );
  }

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = ( i - c[0] ) / c[1];
    j = (int)v;
    if ( j < auxPar->min || j >= auxPar->max ) {
      continue;
    }
    e = (j+1-v)*auxPar->histo[j] + (v-j) * auxPar->histo[j+1];
    s += e * auxPar->partialHisto[i] / (sum_sum * sum_histo );

  }

  fprintf( stderr, "[%7.3f %7.5f] %10.5f \r", c[0], c[1], e );

  return( -s );
  */
}



double _SADOfOneHisto( double *c, void *par )
{
  /*
  typeAuxParam *auxPar = (typeAuxParam *)par;
  double sum_sum, sum_histo;
  int i, j;
  double s = 0.0;
  double v, e=0.0;
  */

  exit( 0 );
  /*
  sum_sum = 0.0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_sum += auxPar->partialHisto[i];
  }

  sum_histo = 0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_histo += auxPar->histo[i];
  }
  

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = c[0] + c[1] * i;
    j = (int)v;

    auxPar->mes[i].ind = i;
    
    if ( j < auxPar->min || j >= auxPar->max ) {
      e  = auxPar->histo[i]/sum_histo;
      auxPar->mes[i].err = fabs( e );
      s += fabs(e);
      
      continue;
    }
    e  = ( (v-j)*auxPar->partialHisto[j+1] +
	   (j+1-v)*auxPar->partialHisto[j] ) / (sum_sum);
    e -= auxPar->histo[i] / sum_histo;
    auxPar->mes[i].err = fabs( e );
    s += fabs(e);
  }

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = ( i- c[0]) /c[1];
    j = (int)v;

    auxPar->mes[i].ind = i;
    
    if ( j < auxPar->min || j >= auxPar->max ) {
      e  = auxPar->partialHisto[i]/sum_sum;
      auxPar->mes[i].err = fabs( e );
      s += fabs(e);
      
      continue;
    }
    e  = ( (v-j)*auxPar->histo[j+1] +
	   (j+1-v)*auxPar->histo[j] ) / (sum_histo);
    e -= auxPar->partialHisto[i] / sum_sum;
    auxPar->mes[i].err = fabs( e );
    s += fabs(e);
  }

  fprintf( stderr, "[%7.3f %7.5f] %10.5f \r", c[0], c[1], e );

  return( s );
  */
}


double _SSDOfOneHisto( double *c, void *par )
{
  /*
  typeAuxParam *auxPar = (typeAuxParam *)par;
  double sum_sum, sum_histo;
  int i, j;
  double s = 0.0;
  double v, e=0.0;
  */
  
  exit(0);
  /*
  sum_sum = 0.0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_sum += auxPar->partialHisto[i];
  }

  sum_histo = 0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_histo += auxPar->histo[i];
  }
  
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = c[0] + c[1] * i;
    j = (int)v;

    auxPar->mes[i].ind = i;
    
    if ( j < auxPar->min || j >= auxPar->max ) {
      e  = auxPar->histo[i]/sum_histo;
      auxPar->mes[i].err = fabs( e );
      s += e*e;
      
      continue;
    }
    e  = ( (v-j)*auxPar->partialHisto[j+1] +
	   (j+1-v)*auxPar->partialHisto[j] ) / (sum_sum);
    e -= auxPar->histo[i] / sum_histo;
    auxPar->mes[i].err = fabs( e );
    s += e*e;
  }



  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = ( i- c[0]) /c[1];
    j = (int)v;

    auxPar->mes[i].ind = i;
    
    if ( j < auxPar->min || j >= auxPar->max ) {
      e  = auxPar->partialHisto[i]/sum_sum;
      auxPar->mes[i].err = fabs( e );
      s += e*e;
      
      continue;
    }
    e  = ( (v-j)*auxPar->histo[j+1] +
	   (j+1-v)*auxPar->histo[j] ) / (sum_histo);
    e -= auxPar->partialHisto[i] / sum_sum;
    auxPar->mes[i].err = fabs( e );
    s += e*e;
  }

  fprintf( stderr, "[%7.3f %7.5f] %10.5f \r", c[0], c[1], e );

  return( s );
  */
}


double _entropyOfOneHisto( double *c, void *par )
{
  /*
  typeAuxParam *auxPar = (typeAuxParam *)par;
  double sum_sum, sum_histo;
  int i, j;
  double v;
  double s = 0;
  */

  exit( 0 );
  /*
  sum_sum = 0.0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_sum += auxPar->partialHisto[i];
  }

  sum_histo = 0;
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    sum_histo += auxPar->histo[i];
  }
  */

  /* initialisation de la somme 
   */
  /*
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    auxPar->sum[i] = auxPar->partialHisto[i] / sum_sum;
  }
  */

  /* on ajoute l'histogramme
   */
  /*
  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    v = c[0] + c[1] * i;
    j = (int)v;

    if ( j < auxPar->min || j >= auxPar->max ) {
      if ( auxPar->histo[i] > 0.0 ) {
	s += - ( auxPar->histo[i]/sum_histo ) * log (  auxPar->histo[i]/sum_histo );
      }
      continue;
    }
    auxPar->sum[ j ] += (j+1-v) * auxPar->histo[i]/sum_histo;
    auxPar->sum[ j+1 ] += (v-j) * auxPar->histo[i]/sum_histo;
  }

  for ( i=auxPar->min; i<=auxPar->max; i++ ) {
    if ( auxPar->sum[i] > 0.0 ) {
      s += - auxPar->sum[i] * log ( auxPar->sum[i] );
    }
  }

  fprintf( stderr, "[%7.3f %7.5f] %10.5f \r", c[0], c[1], s );

  return( s );
  */
}
