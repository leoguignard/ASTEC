#include <iopnm.h>

int main( int argc, char* argv[] )
{
  int dx, dy, dz;
  int nx=100, ny=100;
  int ix = 168, iy =230;
  unsigned char *theBuf;
  unsigned char *resBuf;
  int i, j;
  int ds, ns;

  theBuf = _readPnmImage( argv[1], &dx, &dy, &dz );
  ds = dx*dy;
  ns = nx*ny;

  resBuf = (unsigned char *)malloc( nx*ny*dz * sizeof(unsigned char) );
  for (j=0; j<ny; j++ )
  for (i=0; i<nx; i++ ) {
    resBuf[ j*nx + i ] = theBuf[ (j+iy)*dx + (i+ix) ];
    
  }
  if ( dz == 3 ) {
    for (j=0; j<ny; j++ )
    for (i=0; i<nx; i++ ) {
      resBuf[ ns + j*nx + i ] = theBuf[ ds + (j+iy)*dx + (i+ix) ];
      resBuf[ 2*ns + j*nx + i ] = theBuf[ 2*ds + (j+iy)*dx + (i+ix) ];
    }
  }
  

  _writePnmImage( argv[2], nx, ny, dz, resBuf );

}
