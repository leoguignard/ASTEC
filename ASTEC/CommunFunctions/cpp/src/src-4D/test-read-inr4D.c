#include <stdio.h>

main( int argc, char *argv[] )
{
  FILE *fichier, *fopen();
  int readingStandardInput = 0;
  char name[256];
  int i, status, nb;

  if ( argc == 1 ) {
    fprintf( stderr, "standard input\n");
    fichier = stdin;
    readingStandardInput = 1;

  } else {
    fprintf( stderr, "fichier = %s\n", argv[1] );
    fichier = fopen( argv[1] , "r" );
  }
  
/*
# 4D Inrimage File version 1.000000 
CONTENT 
*/

  do {
    status = fscanf( fichier, "%s", name );
  } while ( (status == 1) && (strcmp( name, "CONTENT" ) != 0) );
/*  
  Number of 4D INRIMAGE = 24
*/
  do {
    status = fscanf( fichier, "%s", name );
  } while ( (status == 1) && (strcmp( name, "=" ) != 0) );
  (void)fscanf( fichier, "%d", &nb );
  fprintf( stderr, "contient %d 3D images\n", nb );

/*
  ASCII
  Big endian
END
*/
  do {
    status = fscanf( fichier, "%s", name );
  } while ( (status == 1) && (strcmp( name, "END" ) != 0) );
/*
4D INRIMAGE kikinis.inr4D
  FileName
END
......blabla.inr
*/
  i = 0;
  do {
    do {
      status = fscanf( fichier, "%s", name );
    } while ( (status == 1) && (strcmp( name, "END" ) != 0) );
    if ( (status = fscanf( fichier, "%s", name )) == 1 )
      fprintf( stderr, "%3d: %s\n", i++, name );
  } while ( status == 1 );

}
