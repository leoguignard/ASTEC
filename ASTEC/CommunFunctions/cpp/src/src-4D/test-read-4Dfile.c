#include <stdio.h>

main( int argc, char *argv[] )
{
  FILE *fichier, *fopen();
  int readingStandardInput = 0;
  char name[256];

  if ( argc == 1 ) {
    fprintf( stderr, "standard input\n");
    fichier = stdin;
    readingStandardInput = 1;

  } else {
    fprintf( stderr, "fichier = %s\n", argv[1] );
    fichier = fopen( argv[1] , "r" );
  }
  
  while ( fscanf( fichier, "%*s", name ) == 0 ) {
    fprintf( stderr, "%s\n", name );

  }

}
