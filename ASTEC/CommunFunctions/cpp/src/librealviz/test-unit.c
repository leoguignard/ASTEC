
#include <iopnm.h>
#include <rl_interpol.h>

int main( int argc, char* argv[] )
{
  unsigned char *buffer = (unsigned char *)NULL;
  
  int bufferDims[3] = {0,0,0};
  int z, nbytes;

  if ( argc == 1 ) {
    fprintf( stderr, "Usage: %s image\n", argv[0] );
    exit(0);
  }
  

  /* lecture de l'image d'entree
   */
  buffer = _readPnmImage( argv[1], &bufferDims[0], &bufferDims[1], &bufferDims[2], &nbytes );

  if ( nbytes != 1 ) {
    fprintf(stderr, "%s: unable to deal with such buffer type\n", argv[0] );
    free( buffer );
    exit( 0 );
  }
  
  for ( z=0; z<bufferDims[2]; z++ ) {
    fprintf(stderr, "%s: testing slice #%d of image %s\n", argv[0], z, argv[1] );
    if ( _TestConversion( &(buffer[z*bufferDims[0]*bufferDims[1]]), bufferDims[0], bufferDims[1] ) 
	 != 1 ) {
      fprintf(stderr, "%s: error in testing slice #%d of image %s\n", argv[0], z, argv[1] );
    }
  }

  free( buffer );
  exit( 0 );
}
