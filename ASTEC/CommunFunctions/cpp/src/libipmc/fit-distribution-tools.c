/*************************************************************************
 * fit-distribution-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 19 jui 2013 22:49:08 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <fcntl.h>
#include <math.h>
#include <unistd.h>

#include <fit-distribution-tools.h>
#include <levenberg.h>

static int _verbose_ = 1;




/************************************************************
 *
 * MANAGEMENT
 *
 ************************************************************/



void initMeasureList( measureList *l )
{
  l->data = (int*)NULL;
  l->n_data = 0;
  l->n_allocated_data = 0;
} 



void freeMeasureList( measureList *l )
{
  if ( l->data != (int*)NULL ) free( l->data );
  initMeasureList( l );
} 



static int _size_to_be_allocated_ = 1000;

int addMeasureToList( measureList *l, int m )
{
  char *proc = "addMeasureToList";
  size_t s =  l->n_allocated_data;
  int *data;

  if ( l->n_data == l->n_allocated_data ) {
    s += _size_to_be_allocated_;
    data = (int*)malloc( s * sizeof(int) );
    if ( data == (int*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: allocation error\n", proc );
      return( -1 );
    } 
    if ( l->n_allocated_data > 0 ) {
      (void)memcpy( data, l->data, l->n_allocated_data*sizeof(int) );
      free( l->data );
    }
    l->n_allocated_data = s;
    l->data = data;
  }
  l->data[l->n_data] = m;
  l->n_data ++;
  return( 1 );
} 



int readMeasureList( measureList *l, char *filename )
{
  char *proc = "readMeasureList";
  FILE *f;
  int a, m;

  f = fopen( filename, "r" );
  if ( f == (FILE*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s'\n", proc, filename );
    return( -1 );
  }

  a = 0;
  while ( fscanf( f, "%d", &m ) == 1 ) {
    if ( addMeasureToList( l, m ) != 1 ) {
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to add measure #%d (=%d) to list\n", proc, a, m );
      return( -1 );
    }
    a ++;
  }
  
  if ( _verbose_ )
    fprintf( stderr, "%s: read %d values\n", proc, a );

  fclose( f );
  return( 1 );
}



int maxMeasureList( measureList *l )
{
  char *proc = "_MaxMeasureList";
  int i, max;
  
  if ( l->n_data == 0 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: empty list\n", proc );
    return( -1 );
  }

  max = l->data[0];
  for ( i=1; i<l->n_data; i++ )
    if ( max < l->data[i] ) max = l->data[i];
  
  return( max );
}


int minMeasureList( measureList *l )
{
  char *proc = "_MinMeasureList";
  int i, min;
  
  if ( l->n_data == 0 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: empty list\n", proc );
    return( -1 );
  }

  min = l->data[0];
  for ( i=1; i<l->n_data; i++ )
    if ( min > l->data[i] ) min = l->data[i];
  
  return( min );
}



/************************************************************
 *
 * 
 *
 ************************************************************/

int build1DHistogramFromMeasureList( typeHistogram *h, measureList *l )
{
  char *proc = "build1DHistogramFromMeasureList";
  unionValues min, max;

  min.val_s32 = 0;
  max.val_s32 = maxMeasureList( l );

  if ( allocHistogramHeader( &(h->xaxis), &min, &max, -1.0, SINT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram header\n", proc );
    return( -1 );
  }

  if ( allocHistogramData( h, SINT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate histogram\n", proc );
    return( -1 );
  }
  
  if ( fill1DHistogramFromBuffer( h, (void*)(l->data), SINT, l->n_data ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to fill histogram\n", proc );
    return( -1 );
  }

  return( 1 );
}




/************************************************************
 *
 * output
 *
 ************************************************************/

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

static void printWeibullXxxlab( char *name, 
				typeHistogram *histo,
				typeHistogram *pdf,
				typeHistogram *cumul,
				double *theParam,
				enumHistogramFile xxxlab )
{
  char *proc = "printWeibullXxxlab";
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

    if ( fprintf1DHistogramMatlab( f, fd, histo, "histogram" ) != 1 ) {
      free( filename );
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write data and matlab script in '%s.raw' and '%s.m'\n", 
		 proc, filename, filename );
      return;
    }
    
    fprintf( f, "\n" );
    fprintf( f, "fclose(f);\n" );
    fprintf( f, "\n" );

    break;
    
  case _SCILAB_ :


    fprintf( f, "\n" );
    fprintf( f, "f=mopen('%s.raw','r');\n", _BaseName( template ) );
    fprintf( f, "\n" );
    
    if ( fprintf1DHistogramScilab( f, fd, histo, "histogram"  ) != 1 ) {
      free( filename );
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write data and scilab script in '%s.raw' and '%s.m'\n", 
		 proc, filename, filename );
      return;
    }
    
    if ( fprintf1DHistogramScilab( f, fd, pdf, "pdf"  ) != 1 ) {
      free( filename );
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write data and scilab script in '%s.raw' and '%s.m'\n", 
		 proc, filename, filename );
      return;
    }

    fprintf( f, "w1 = %f * (%f/%f) * (INDEX_pdf/%f)^(%f-1) .* exp(-(INDEX_pdf/%f)^%f);\n", 
	     theParam[0], 2.0, theParam[1], theParam[1], 2.0, theParam[1], 2.0 );
    fprintf( f, "w2 = %f * (%f/%f) * (INDEX_pdf/%f)^(%f-1) .* exp(-(INDEX_pdf/%f)^%f);\n", 
	     theParam[2], 2.0, theParam[3], theParam[3], 2.0, theParam[3], 2.0 );
    fprintf( f, "w3 = %f * (%f/%f) * (INDEX_pdf/%f)^(%f-1) .* exp(-(INDEX_pdf/%f)^%f);\n", 
	     theParam[4], 1.0, theParam[5], theParam[5], 1.0, theParam[5], 1.0 ); 
    fprintf( f, "\n" );

    fprintf( f, "plot( INDEX_pdf, w1, \"r-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, w2, \"g-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, w3, \"b-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, w1+w2+w3, \"m-\", \"thickness\", 2 );\n" );
    fprintf( f, "\n" );





    if ( fprintf1DHistogramScilab( f, fd, cumul, "cumulative"  ) != 1 ) {
      free( filename );
      close( fd );
      fclose( f );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write data and scilab script in '%s.raw' and '%s.m'\n", 
		 proc, filename, filename );
      return;
    }

    fprintf( f, "cw1 = %f * (1.0 - exp(-(INDEX_pdf/%f)^%f));\n", 
	     theParam[0], theParam[1], 2.0 );
    fprintf( f, "cw2 = %f * (1.0 - exp(-(INDEX_pdf/%f)^%f));\n", 
	     theParam[2], theParam[3], 2.0 );
    fprintf( f, "cw3 = %f * (1.0 - exp(-(INDEX_pdf/%f)^%f));\n", 
	     theParam[4], theParam[5], 2.0 );
    fprintf( f, "\n" );

    fprintf( f, "plot( INDEX_pdf, cw1, \"r-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, cw2, \"g-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, cw3, \"b-\", \"thickness\", 2 );\n" );
    fprintf( f, "plot( INDEX_pdf, cw1+cw2+cw3, \"m-\", \"thickness\", 2 );\n" );
    fprintf( f, "\n" );

    
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





/************************************************************
 *
 * output
 *
 ************************************************************/

/* la distribution cumulative de Weibull a la forme
   PDF(x) = a k/t (x/t)^(k-1)  exp( -(x/t)^k )

   la distribution cumulative de Weibull a la forme
   CDF(x) = a ( 1.0 - exp( -(x/t)^k ) )

   on fait la somme de 3 distributions i=1..3 dont les parametres
   sont 
   thePar[0] = a1
   thePar[1] = t1
   thePar[2] = a2
   thePar[3] = t2
   thePar[4] = a3
   thePar[5] = t3

   on fixe k1=2, k2=2, k3=1 (cf Loerke et al, Plos Biology, 7(3):e1000057, 2009

   derPar contient les derivees par rapport a ces parametres

*/



double _MixtureOf3Weibull2Parameters( double x, 
				      double *thePar, 
				      double *derPar ) 
{
  double r = 0.0;
  double e, xt;
  double k;
  int i;

  for ( i=0; i<3; i++ ) {
    /* (x/t)^k */
    xt  = x / thePar[i*2+1];
    /* exp( -(x/t)^k ) */
    if ( i == 0 || i == 1 ) e = exp( - xt * xt );
    else e = exp( - xt  );	   
    /* k */
    k = ( i == 0 || i == 1 ) ? 2 : 1;
    
    /* d/da */
    derPar[i*2]    = ( i == 0 || i == 1 ) ? xt : 1.0;
    derPar[i*2]   *= (k / thePar[i*2+1]);
    derPar[i*2]   *= e;
    /* d/dt */
    derPar[i*2+1]  = (k / thePar[i*2+1]) * xt;
    if ( i == 0 || i == 1 ) derPar[i*2+1]  *= xt;
    derPar[i*2+1] -= k * xt;
    derPar[i*2+1] *= thePar[i*2] * (k / thePar[i*2+1]) * (1.0 / thePar[i*2+1]);
    derPar[i*2+1] /= ( i == 0 || i == 1 ) ? 1.0 : xt;
    derPar[i*2+1] *= e;
    
    if ( i == 0 || i == 1 )
      r +=  thePar[i*2] * (k / thePar[i*2+1]) * xt * e;
    else 
      r +=  thePar[i*2] * (k / thePar[i*2+1]) * e;
  }
  return( r );
}



double _MixtureOf3CumulativeWeibull2Parameters( double x, 
						double *thePar, 
						double *derPar ) 
{
  double r = 0.0;
  double e, xtk;
  double k;
  int i;

  for ( i=0; i<3; i++ ) {
    /* (x/t)^k */
    xtk = x / thePar[i*2+1];
    if ( i == 0 || i == 1 ) xtk *= xtk;
    /* exp( -(x/t)^k ) */
    e = exp( - xtk );
    /* k */
    k = ( i == 0 || i == 1 ) ? 2 : 1;
    /* d/da */
    derPar[i*2]   = (1.0 - e);
    /* d/dt */
    derPar[i*2+1] = -( thePar[i*2] * (k / thePar[i*2+1]) * xtk ) * e;
    r +=  thePar[i*2] * ( 1.0 - e );
  }
  return( r );
}





/* la distribution cumulative de Weibull a la forme
   CDF(x) = a ( 1.0 - exp( -(x/t)^k ) )
   on fait la somme de 3 distributions i=1..3 dont les parametres
   sont 
   thePar[0] = a1
   thePar[1] = t1
   thePar[2] = k1

   thePar[3] = a2
   thePar[4] = t2
   thePar[5] = k2

   thePar[6] = a3
   thePar[7] = t3
   thePar[8] = k3

   derPar contient les derivees par rapport a ces parametres

*/


double _MixtureOf3CumulativeWeibull3Parameters( double x, 
						double *thePar, 
						double *derPar ) 
{
  double r = 0.0;
  double e, xtk;
  int i;

  for ( i=0; i<3; i++ ) {
    /* (x/t)^k */
    xtk = x / thePar[i*3+1];
    xtk = exp( thePar[i*3+2] * log( xtk ) );
    /* exp( -(x/t)^k ) */
    e = exp( - xtk );
    /* k */
    derPar[i*3]   = (1.0 - e);
    derPar[i*3+1] = -( thePar[i*3] * (thePar[i*3+2] / thePar[i*3+1]) * xtk ) * e;
    derPar[i*3+2] = -( thePar[i*3] * log(x / thePar[i*3+1]) * xtk ) * e;
    r +=  thePar[i*3] * ( 1.0 - e );
  }
  return( r );
}



/************************************************************
 *
 * output
 *
 ************************************************************/


int fit3WeibullDistributionsOnCumulative( typeHistogram *h, 
					  double *theParamATK,
					  char *filename,
					  enumHistogramFile xxxlab,
					  int minlength,
					  int maxlength )
{
  char *proc = "fitWeibullDistribution";
  typeHistogram c;
  typeHistogram d;
  
  int n_parameters = 2;
  double *theX, *theY, *theC, *theS;
  double theParamAT[6];
  int i, j;
  int firsti, lasti, nsamples;



  /* first and last index
     the minlength and maxlength parameters allow to change the interval of computation
   */
  firsti = 0;
  lasti = c.xaxis.dim-1;

  switch( h->typeHisto ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: histogram type not handled yet (computation of first and last index)\n", proc );
    return( -1 );
  case SINT :
    {
      s32 *theHisto = (s32*)h->data;
      for ( firsti=0; theHisto[firsti] == 0 && firsti<h->xaxis.dim; firsti++ ) 
	;
      if ( firsti >= h->xaxis.dim ) firsti = h->xaxis.dim - 1;
      for ( lasti=h->xaxis.dim-1; theHisto[lasti] == 0 && lasti>=0; lasti++ ) 
	;
      if ( lasti < 0 ) lasti = 0;
    }
    break;
  }
  
  if ( _verbose_ ) {
    fprintf( stderr, "%s: non-null interval is [%d,%d] in [%d,%d]\n",
	     proc, firsti, lasti, 0, h->xaxis.dim-1 );
  }
  
  if ( minlength >= 0 && firsti < minlength ) firsti = minlength;
  if ( maxlength >= 0 && lasti > maxlength ) lasti = maxlength;

  nsamples = lasti - firsti + 1;
  
  if ( _verbose_ ) {
    fprintf( stderr, "%s: retained interval is [%d,%d] in [%d,%d]\n",
	     proc, firsti, lasti, 0, h->xaxis.dim-1 );
  }





  /* histogram
   */
  initHistogram( &c );
  initHistogram( &d );

  if ( allocHistogramFromHistogram( &d, h, FLOAT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate pdf histogram\n", proc );
    return( -1 );
  }
  if ( pdfHistogram( &d, h ) != 1 ) {
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute pdf histogram\n", proc );
    return( -1 );
  }

  if ( allocHistogramFromHistogram( &c, &d, FLOAT ) != 1 ) {
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate cumulative histogram\n", proc );
    return( -1 );
  }
  if ( cumulative1DHistogram( &c, &d ) != 1 ) {
    freeHistogram( &c );
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute cumulative histogram\n", proc );
    return( -1 );
  }

  /* ajustement de courbes
   */
  theX = (double*)malloc( 4* nsamples * sizeof(double) );
  if ( theX == (double*)NULL ) {
    freeHistogram( &c );
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary array for fitting\n", proc );
    return( -1 );
  }
  theY = theS = theC = theX;
  theY += nsamples;
  theS += 2*nsamples;
  theC += 3*nsamples;
  
  for ( j=0; j<nsamples; j++ ) 
    theX[j] = theY[j] = theC[j] = theS[j] = 1.0;

  switch ( c.xaxis.typeIndex ) {
  default :
    free( theX );
    freeHistogram( &c );
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: index type not handled yet\n", proc );
    return( -1 );
  case SINT :
    {
      s32 *theIndex = (s32*)c.xaxis.index;
      for ( j=0, i=firsti; i<=lasti; j++, i++ ) 
	theX[j] = theIndex[i];
    }
    break;
  }

  switch ( c.typeHisto ) {
  default :
    free( theX );
    freeHistogram( &c );
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: histogram type not handled yet\n", proc );
    return( -1 );
  case FLOAT :
    {
      r32 *theHisto = (r32*)c.data;
      for ( j=0, i=firsti; i<=lasti; j++, i++ ) 
	theY[j] = theHisto[i];
    }
    break;
  }



  switch( n_parameters ) {
  default :
    free( theX );
    freeHistogram( &c );
    freeHistogram( &d );
    if ( _verbose_ )
      fprintf( stderr, "%s: such parameter numbers not handled yet\n", proc );
    return( -1 );
  case 2 :
    theParamAT[0] = theParamATK[0];
    theParamAT[1] = theParamATK[1];
    theParamAT[2] = theParamATK[3];
    theParamAT[3] = theParamATK[4];
    theParamAT[4] = theParamATK[6];
    theParamAT[5] = theParamATK[7];
    if ( Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				      theC, DOUBLE, theS, DOUBLE, 
				      nsamples, theParamAT, 6,
				      _MixtureOf3CumulativeWeibull2Parameters ) != 1 ) {
      free( theX );
      freeHistogram( &c );
      freeHistogram( &d );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when fitting\n", proc );
      return( -1 );
    }
    theParamATK[0] = theParamAT[0];
    theParamATK[1] = theParamAT[1];
    theParamATK[3] = theParamAT[2];
    theParamATK[4] = theParamAT[3];
    theParamATK[6] = theParamAT[4];
    theParamATK[7] = theParamAT[5];
    break;
  case 3 :
    if ( Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				      theC, DOUBLE, theS, DOUBLE, 
				      nsamples, theParamATK, 9,
				      _MixtureOf3CumulativeWeibull3Parameters ) != 1 ) {
      free( theX );
      freeHistogram( &c );
      freeHistogram( &d );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when fitting\n", proc );
      return( -1 );
    }
  }

  fprintf( stderr, "%s: %f %f %f %f %f %f %f %f %f\n", proc, 
	   theParamATK[0], theParamATK[1], theParamATK[2],
	   theParamATK[3], theParamATK[4], theParamATK[5],
	   theParamATK[6], theParamATK[7], theParamATK[8] );

  free( theX );

  /* output
   */
  printWeibullXxxlab( filename, h, &d, &c, theParamAT, xxxlab );
  freeHistogram( &c );
  freeHistogram( &d );

  return( 1 );
}
