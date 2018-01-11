#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mat-rot.c"

main( argc, argv )
int argc;
char *argv[];
{
    int i, j;
    int n, nb=1000000;
    double max;
    long int init = time(0);
    
    double e;
    double v0[3];
    double mat0[4][4];
    double v1[3];
    double mat1[4][4];

    max = RAND_MAX; /* (2^31)-1 = 2147483647*/
    
    init = 1034870980;
    init = time(0);
    fprintf( stderr, "init = %ld\n", init );
    srandom( init );
 
    for (n=0;n<nb;n++) {
      
      fprintf( stderr, "\r test %7d/%d", n+1, nb );
      
      do {
	v0[0] = random() / max;
	v0[1] = random() / max;
	v0[2] = random() / max;
	
	for ( e=0.0, i=0; i<3; i++ )
	  e += v0[i] * v0[i];
      }
      while ( e < 1e-6 );


      e = sqrt( e );
      for ( i=0; i<3; i++ ) v0[i] /= e;
      
      e = 3.141592653589793 * random() / max;
      for ( i=0; i<3; i++ ) v0[i] *= e;

      if ( 0 ) fprintf( stderr, " angle = %f\n", e );

      Rvector_To_Rmatrix( mat0, v0 );
      Rmatrix_To_Rvector( mat0, v1 );
      Rvector_To_Rmatrix( mat1, v1 );

      for ( e=0.0, i=0; i<3; i++ )
	e += (v0[i]-v1[i]) * (v0[i]-v1[i]);
      e /= 3;
      
      if ( e > 1e-10 ) {
	fprintf( stderr, "\n" );
	printf( "vecteurs #0 = %f %f %f\n", v0[0], v0[1], v0[2] );
	printf( "         #1 = %f %f %f\n", v1[0], v1[1], v1[2] );
      }

      for ( e=0.0, i=0; i<3; i++ )
	for ( j=0; j<3; j++ )
	  e += (mat0[i][j]-mat1[i][j]) * (mat0[i][j]-mat1[i][j]);
      e /= 9;

      if ( e > 1e-10 ) {
	fprintf( stderr, "\n" );
	printf( "matrices #0 = %f %f %f\n\
              %f %f %f\n\
              %f %f %f\n", mat0[0][0], mat0[0][1], mat0[0][2],
		mat0[1][0], mat0[1][1], mat0[1][2],
		mat0[2][0], mat0[2][1], mat0[2][2] );
	printf( "         #1 = %f %f %f\n\
              %f %f %f\n\
              %f %f %f\n", mat1[0][0], mat1[0][1], mat1[0][2],
		mat1[1][0], mat1[1][1], mat1[1][2],
		mat1[2][0], mat1[2][1], mat1[2][2] );

      }
    }
    fprintf( stderr, "\n" );
}
