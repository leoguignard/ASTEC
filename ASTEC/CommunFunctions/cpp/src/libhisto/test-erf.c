#include <stdio.h>
#include <math.h>

double erf( double x );

int main( int argc, char *argv[] )
{
  double dx, x, lx = 20;
  double s, ls = 4.0;
  int i;

  /* l'integrale de 0 a x de G(x) = 1 / (sigma * sqrt( 2 * pi )) exp ( - x*x/(2*sigma*sigma))
     est 1/2 erf( x / (sigma * sqrt(2)) )
  */

  for ( s = 1.0; s <= ls ; s += 1.0 ) {
    printf( "sigma = %f\n", s );
    
    for ( x = 0.0, dx=0.5; x <= lx ; x += dx ) {
      if ( x >= 5.0 ) dx = 1.0;
      if ( x >= 10.0 ) dx = 2.0;
      printf( "   1/2 * erf( x / (sigma*sqrt(2)) ) = 1/2 * erf( %f / (%f * %f) ) = %f\n",
	      x, s, sqrt(2.0), 0.5 * erf( x / ( s * sqrt( 2.0 ) ) ) );
    }
  }

  printf( "\n" );
  printf( "erf( 999) =  %f\n", erf( 999.0 ) );
  printf( "erf(-999) = %f\n", erf(-999.0 ) );
  printf( "\n" );


  for ( i=-10; i<10; i++ ) {
    printf( "%2d : 1.0 - |%f| = %g\n", i, erf((double)i), 1.0 - fabs(erf((double)i)) );
  }


}
