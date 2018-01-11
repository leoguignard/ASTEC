#include <stdlib.h>
#include <stdio.h>
#include <time.h>


#include <vt_jacobien.h>

int main( int argc, char *argv[] )
{
  JCB_InitRandom( -1 );

  fprintf( stdout, "\n\n*** Tests translations ***\n\n" );

  JCB_TestTranslation( 1 );

  fprintf( stdout, "\n\n*** Tests similitudes ***\n\n" );

  JCB_TestSimilitude( 1 );

  fprintf( stdout, "\n\n*** Tests affines ***\n\n" );

  JCB_TestAffine( 2 );

  return( 1 );
}
