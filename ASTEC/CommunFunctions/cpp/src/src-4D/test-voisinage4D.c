#include <stdio.h>

typedef struct _my_vois{
  int x;
  int y;
  int z;
  int t;
} _my_vois;

main( int argc, char *argv[] )
{
  int n=1, s;
  int nb=0;
  _my_vois v[40];
  int t, z, y, x;

  if ( (argc == 1) || (sscanf( argv[1], "%d", &n ) != 1 ) )
    n = 1;

  for ( t=-1; t<=0; t++ )
  for ( z=-1; z<=1; z++ )
  for ( y=-1; y<=1; y++ )
  for ( x=-1; x<=1; x++ ) {
    if ( (t == -1)
	 || ( (t == 0) && (z == -1) )
	 || ( (t == 0) && (z == 0) && (y == -1) )
	 || ( (t == 0) && (z == 0) && (y == 0) && (x == -1) ) ) {
      s = 0;
      if ( (x == 1) || (x == -1) ) s++;
      if ( (y == 1) || (y == -1) ) s++;
      if ( (z == 1) || (z == -1) ) s++;
      if ( (t == 1) || (t == -1) ) s++;
      if ( (s > 0) && (s <= n) ) {
	v[nb].x = x;
	v[nb].y = y;
	v[nb].z = z;
	v[nb].t = t;
	nb ++;
      }
    }
  }

  printf( " il y a %2d voisins :\n", nb );
  for ( t=0; t<nb; t++ ) 
    printf( "%3d %3d %3d %3d\n", v[t].t, v[t].z, v[t].y, v[t].x );
}


