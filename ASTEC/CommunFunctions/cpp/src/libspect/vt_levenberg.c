/*************************************************************************
 * vt_levenberg.c -
 *
 * $Id: vt_levenberg.c,v 1.5 2000/06/23 15:59:55 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 09:23:04 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <vt_levenberg.h>


static int _verbose_ = 0;


/* y = thePar[0] * x / (thePar[1] + x)
   x -> xenon
   y -> hmpao
*/
double _LassenFunction( double x, double *thePar, double *derPar )
{
  if ( thePar[1] + x == 0 ) {
    derPar[0] = derPar[1] = 0.0;
    return( 0.0 );
  }
  derPar[0] = x / (thePar[1] + x);
  derPar[1] = - thePar[0] * derPar[0] / (thePar[1] + x);
  return( thePar[0] * derPar[0] );
}



/* amplitude = thePar[0]
   moyenne = thePar[1]
   sigma = thePar[2]

   gaussienne = thePar[0] * exp ( - (x-thePar[1])*(x-thePar[1])/(2*thePar[2]*thePar[2]) )
*/ 
double _GaussianForLM ( double x, double *thePar, double *derPar )
{
  double a, e;
  
  if ( thePar[2] == 0.0 ) {
    derPar[0] = derPar[1] = derPar[2] = 0.0;
    return( 0.0 );
  }

  a = (x-thePar[1])/thePar[2];
  e = exp ( - a*a/2.0 );
  
  derPar[0] = e;
  derPar[1] = thePar[0] * e * a / thePar[2];
  derPar[2] = thePar[0] * e * a * a / thePar[2];
  
  return( thePar[0] * e );
}

/* amplitude = thePar[0]
   moyenne = thePar[1]
   sigma = thePar[2]

   gaussienne = thePar[0] * exp ( - (x-thePar[1])*(x-thePar[1])/(2*thePar[2]*thePar[2]) )
*/ 
double _NonSymetricGaussianForLM ( double x, double *thePar, double *derPar )
{
  double a, e;


  if ( x < thePar[1] ) {

    if ( thePar[2] == 0.0 ) {
      derPar[0] = derPar[1] = derPar[2] = derPar[3] = 0.0;
      return( 0.0 );
    }

    a = (x-thePar[1])/thePar[2];
    e = exp ( - a*a/2.0 );
  
    derPar[0] = e;
    derPar[1] = thePar[0] * e * a / thePar[2];
    derPar[2] = thePar[0] * e * a * a / thePar[2];
    derPar[3] = 0.0;

    return( thePar[0] * e );
  }
  

  if ( thePar[3] == 0.0 ) {
    derPar[0] = derPar[1] = derPar[2] = derPar[3] = 0.0;
    return( 0.0 );
  }

  a = (x-thePar[1])/thePar[3];
  e = exp ( - a*a/2.0 );
  
  derPar[0] = e;
  derPar[1] = thePar[0] * e * a / thePar[3];
  derPar[2] = 0.0;
  derPar[3] = thePar[0] * e * a * a / thePar[3];
  
  return( thePar[0] * e );
}














/* on evalue le chi2 
   = sum(i=0..length) w[i] ( ((y[i] - f(x[i]))/s[i])^2 )
*/

static double _ComputeChi2AndDerivatives( void *theX, bufferType xType,
					  void *theY, bufferType yType,
					  void *theW, bufferType wType,
					  void *theS, bufferType sType,
					  int length,
					  double *theParams,
					  double *derParams,
					  int nbParams,
					  typeFuncEval funcEval,
					  double *alpha,
					  double *beta )
{
  char *proc = "_ComputeChi2AndDerivatives";
  double chi = 0.0;
  double x, y, dy, c;
  int i, j, k;
  
  if ( nbParams <= 0 ) {
    return( 0 );
  }

  for ( i=0; i<nbParams*nbParams; i++ ) alpha[i] = 0.0;
  for ( i=0; i<nbParams; i++ ) beta[i] = 0.0;
  
  if ( xType != yType || xType != wType || xType != sType ) {
    if ( 1 || _verbose_ )
      fprintf( stderr, "%s: unable to handle different types\n", proc );
    return( 0.0 );
  }
  
  switch( xType ) {
  case DOUBLE :
    {
      double *bufX = (double*)theX;
      double *bufY = (double*)theY;
      double *bufW = (double*)theW;
      double *bufS = (double*)theS;
      for ( i = 0; i < length; i++ ) {

	if ( bufS[i] == 0.0 ) continue;
	x = bufX[i];
	c = bufW[i] / (bufS[i]*bufS[i]);
	y = (*funcEval)( x, theParams, derParams );
	
	
	dy = bufY[i] - y;
	for ( j = 0; j < nbParams; j++ ) {
	  beta[j] += c * dy * derParams[j];
	  for ( k = 0; k < nbParams; k++ ) {
	    alpha[j*nbParams+k] += c * derParams[j] * derParams[k];
	  }
	}
	
	chi += c * dy * dy;
      }
    }
    break;
  default :
    if ( 1 || _verbose_ )
      fprintf( stderr, "%s: unable to handle such type\n", proc );
    return( 0.0 );
  }
  return( chi );
}
















int VT_Modeling1DDataWithLevenberg( void *theX, bufferType xType,
				    void *theY, bufferType yType,
				    void *theW, bufferType wType,
				    void *theS, bufferType sType,
				    int theLength,
				    double *theParams, int nbParams,
				    typeFuncEval funcEval
				    )
{
  char *proc = "VT_Modeling1DDataWithLevenberg";
  double *matAlpha  = (double*)NULL;
  double *matPrime  = (double*)NULL;
  double *derParams = (double*)NULL;
  double *vecBeta   = (double*)NULL;
  double *incParams = (double*)NULL;
  double *newParams = (double*)NULL;
  double *theAllocatedBuffer = (double*)NULL;
  
  double lambda = 0.001;
  double chi2, newchi2;
  int stop = 0;
  int i, j;

  if ( nbParams <= 0 ) {
    return( 0 );
  }

  /* il faut 2 matrices nbParams * nbParams
             4 vecteurs nbParams
  */
  theAllocatedBuffer = (double*)malloc( (2*(nbParams*nbParams)+4*nbParams)*sizeof(double) );
  if ( theAllocatedBuffer == (double*)NULL ) {
    fprintf( stderr, "%s: unable to allocate buffer\n", proc );
    return( 0 );
  }
  matPrime = matAlpha = theAllocatedBuffer;
  matPrime += nbParams*nbParams;
  derParams = matPrime;
  derParams += nbParams*nbParams;
  vecBeta  = derParams;
  vecBeta  += nbParams;
  incParams = vecBeta;
  incParams += nbParams;
  newParams = incParams;
  newParams += nbParams;
 


  /* calcul du chi 2 initial et des derivees
   */
  chi2 = _ComputeChi2AndDerivatives( theX, xType, theY, yType, theW, wType, theS, sType, 
				     theLength, theParams, derParams, nbParams,
				     funcEval, matAlpha, vecBeta );
  if ( 0 ) {
    printf( " ...... [a] = { %f %f }\n", theParams[0], theParams[1] );
    printf( " ...... [BETA] = { %f %f }\n", vecBeta[0], vecBeta[1] );
    printf( " ...... [ALPHA] = [ %f %f ]\n", matAlpha[0], matAlpha[1] );
    printf( "                  [ %f %f ]\n", matAlpha[2], matAlpha[3] );
  }
  
  do {
    /* on construit A' = A + Lambda Diag (A)
     */
    (void)memcpy( matPrime, matAlpha, nbParams*nbParams*sizeof(double) );
    for ( i = 0; i < nbParams; i++ )
      matPrime[ i*nbParams+i ] += lambda * matAlpha[ i*nbParams+i ];

    /* on resoud A' * incParams = Beta
     */
    if ( _SolveLinearSystem( matPrime, vecBeta, incParams, nbParams ) != 1 ) {
      free( theAllocatedBuffer );
      return( 0 );
    }
    
    if ( 0 ) {
      printf( " ...... [ALPHA'] = [ %f %f ]\n", matPrime[0], matPrime[1] );
      printf( "                   [ %f %f ]\n", matPrime[2], matPrime[3] );
      printf( " ...... [da] = [ %f %f ]\n", incParams[0], incParams[1] );
    }

    if ( _verbose_ >= 2 ) {
	
      for ( i = 0; i < nbParams; i++ ) {
	if ( i == 0 )
	  printf( "... ALPHA =" );
	else 
	  printf( "           " );
	for (j=0; j<nbParams; j++ )
	  printf( "  %f", matPrime[i*nbParams+j] );
	printf( "\n" );
      }
      printf( "... BETA  =" );
      for ( i = 0; i < nbParams; i++ )
	printf( "  %f", incParams[i] );
      printf( "\n" );
    }
    

    /* nouveaux parametres
     */
    for ( i = 0; i < nbParams; i++ )
      newParams[i] = theParams[i] + incParams[i];
    
    newchi2 = _ComputeChi2AndDerivatives( theX, xType, theY, yType, 
					  theW, wType, theS, sType, 
					  theLength, 
					  newParams, derParams, nbParams,
					  funcEval, matPrime, incParams );
    if ( _verbose_ >= 3 ) {
	
      for ( i = 0; i < nbParams; i++ ) {
	if ( i == 0 )
	  printf( "... ALPHA PRIME =" );
	else 
	  printf( "                 " );
	for (j=0; j<nbParams; j++ )
	  printf( "  %f", matPrime[i*nbParams+j] );
	printf( "\n" );
      }
      printf( "... BETA  PRIME =" );
      for ( i = 0; i < nbParams; i++ )
	printf( "  %f", incParams[i] );
      printf( "\n" );
    }


    if ( newchi2 >= chi2 ) {

      lambda *= 10;
      if ( _verbose_ ) {
	printf( "... Chi2DChi = %f > Chi2 = %f", newchi2, chi2 );
	printf( "   =>   lambda *= 10 -> %f \n", lambda );
      }
      if ( lambda > 1e+10 ) stop = 1;

    } else {

      lambda /= 10;
      if ( _verbose_ ) {	
	printf( "... Chi2DChi = %f < Chi2 = %f", newchi2, chi2 );
	printf( "   =>   lambda /= 10 -> %f \n", lambda );
      }
      if ( _verbose_ ) {
	printf( "... old params =" );
	for ( i = 0; i < nbParams; i++ )
	  printf( "  %f", theParams[i] );
	printf( "\n" );
	printf( "... new params =" );
	for ( i = 0; i < nbParams; i++ )
	  printf( "  %f", newParams[i] );
	printf( "\n" );
      }
	
      (void)memcpy( theParams, newParams, nbParams*sizeof(double) );
      (void)memcpy( matAlpha, matPrime, nbParams*nbParams*sizeof(double) );
      (void)memcpy( vecBeta, incParams, nbParams*sizeof(double) );


      if ( (chi2-newchi2)/chi2 < 1e-4 ) stop = 1;

      chi2 = newchi2;

    }

  } while ( stop == 0 );
  




  free( theAllocatedBuffer );
  return( 1 );
}



void VT_Levenberg_verbose ()
{
  if ( _verbose_ <= 0 ) _verbose_ = 1;
  else _verbose_ ++;
}

void VT_Levenberg_noverbose ()
{
  _verbose_ = 0;
}
