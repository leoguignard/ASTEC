#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <vt_distfield.h>

int main( int argc, char *argv[] )
{
  distanceField *dfield = (distanceField *)NULL;
  int x, y, z;
  double max;
  float v[3];

  unsigned char *b, *buf;

  max = 2147483647; /* (2^31)-1 */
  (void)srandom(time(0));
  
  dfield = VT_FillDistanceField( "/u/cajal/1/greg/Jclombar/test1.x",
				 "/u/cajal/1/greg/Jclombar/test1.y",
				 "/u/cajal/1/greg/Jclombar/test1.z" );
  if ( dfield == (distanceField *)NULL) {
    fprintf( stderr, "Quelque chose ne tourne pas rond.\n" );
    exit( 1 );
  }

  b = buf = (unsigned char*)malloc( dfield->dimx * dfield->dimy *
			      dfield->dimz * sizeof( unsigned char ) );
  
  for (z=0; z<dfield->dimz; z ++)
  for (y=0; y<dfield->dimy; y ++)
  for (x=0; x<dfield->dimx; x ++) {
    getDistanceFieldValues( v, x,y,z, dfield );
    if ( ((v[0] < 0.5) && (v[0] > -.5)) &&
	 ((v[1] < 0.5) && (v[1] > -.5)) && 
	 ((v[2] < 0.5) && (v[2] > -.5)) ) 
      *b++ = 255;
    else
      *b++ = 0;
  }

  {
    int fd;
    fd = open( "test-all-disfield.output", O_WRONLY | O_CREAT | O_TRUNC, 0644 );
    if ( write( fd, buf, dfield->dimx * dfield->dimy * dfield->dimz ) == -1 ) {
      fprintf( stderr, "error when writing test-all-disfield.output\n" );
      close ( fd );
      exit( 1 );
    }
    close ( fd );
  }
  VT_FreeDistanceField( dfield );
  free( buf );

  return( 0 );

}
