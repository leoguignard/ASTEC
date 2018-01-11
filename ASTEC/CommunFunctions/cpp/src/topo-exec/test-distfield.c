#include <vt_distfield.h>
#include <time.h>

int main( int argc, char *argv[] )
{
  distanceField *dfield = (distanceField *)NULL;
  int x, y, z;
  int i;
  double max;
  float v[3];
  
  max = 2147483647; /* (2^31)-1 */
  (void)srandom(time(0));
  
  dfield = VT_FillDistanceField( "/u/cajal/1/greg/Jclombar/test1.x",
				 "/u/cajal/1/greg/Jclombar/test1.y",
				 "/u/cajal/1/greg/Jclombar/test1.z" );
  if ( dfield == (distanceField *)NULL) {
    fprintf( stderr, "Quelque chose ne tourne pas rond.\n" );
    exit( 1 );
  }

  for (i=0;i<10;i++) {
    /* random() return successive pseudo-random numbers 
       in the range from 0 to (2^31)-1
    */
    x = (int)( random()/max * dfield->dimx + 0.5 );
    y = (int)( random()/max * dfield->dimy + 0.5 );
    z = (int)( random()/max * dfield->dimz + 0.5 );
    getDistanceFieldValues( v, x,y,z, dfield );
    printf(" %3d : x=%3d y=%3d z=%3d, v = (%7.2f, %7.2f, %7.2f)\n",
	   i,x,y,z,v[0],v[1],v[2] );
  }

  VT_FreeDistanceField( dfield );

  return( 0 );

}
