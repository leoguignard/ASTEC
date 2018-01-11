#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--- aide ---*/
static char *program;
static char *usage = "matrice-in [matrice-out] [-help]";
static char *detail ="\t Inverse 4x4  transformation";

static void _errorMessage( char *str, int flag )
{
  fprintf( stderr," Usage: %s %s\n", program, usage );
  if ( flag == 1 )
    fprintf( stderr," %s", detail );
  if ( str != NULL )
    fprintf( stderr," Error: %s", str );
  exit( -1 );
}

/******************************************

*******************************************/

#include "mat-manip.c"




/* Programme principal */

int main( int argc, char *argv[] )
{
  char *input_file = NULL;
  char *output_file = NULL;

  double inv[4][4], mat[4][4];  
  FILE *file;

  int i, status;
  int nnames;
  char message[256];

  /* 
   * lecture des arguments
   */
  program = argv[0];
  if ( argc == 1 ) _errorMessage( "\n", 0 );
  
  for ( nnames = 0, i=1; i<argc; i++ ) {
    if ( argv[i][0] == '-' ) {

      if ( strcmp ( argv[i], "-help" ) == 0 ) {
	_errorMessage( "\n", 1 );
      }
      else {
	sprintf( message, "unknown option '%s'\n", argv[i] );
	_errorMessage( message, 0 );
      }
    }
    /*
     * names
     */
    else {
      if ( nnames == 0 ) {
	input_file= argv[i];
	nnames ++;
      }
      else {
	_errorMessage( "too much file names when parsing\n", 0 );
      }
    }
  }


  /*
   * lecture
   */

  file = fopen( input_file, "r" );

  if ( file == NULL ) {
    sprintf( message, "Erreur a l'ouverture de '%s'\n", input_file );
    _errorMessage( message, 0 );
  }
  
  status = fscanf(file , "(\nO8\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf", 
		  &mat[0][0], &mat[0][1], &mat[0][2], &mat[0][3], 
		  &mat[1][0], &mat[1][1], &mat[1][2], &mat[1][3], 
		  &mat[2][0], &mat[2][1], &mat[2][2], &mat[2][3], 
		  &mat[3][0], &mat[3][1], &mat[3][2], &mat[3][3]);
  if ( status != 16 ){
    sprintf( message, "Erreur a la lecture de '%s'\n", input_file );
    _errorMessage( message, 0 );
  }

  fclose( file );

    

  fprintf( stdout,"plot3([0 %f], [0 %f], [0 %f],",
	   mat[0][3], mat[1][3], mat[2][3] );

  if ( fabs( mat[2][3] ) >= fabs( mat[1][3] ) && fabs( mat[2][3] ) >= fabs( mat[0][3] ) ) {
    fprintf( stdout,"'r-', 'LineWidth',2 );\n" );
  }
  else if ( fabs( mat[1][3] ) >= fabs( mat[2][3] ) && fabs( mat[1][3] ) >= fabs( mat[0][3] ) ) {
    fprintf( stdout,"'b-', 'LineWidth',2 );\n" );
  }
  else if ( fabs( mat[0][3] ) >= fabs( mat[2][3] ) && fabs( mat[0][3] ) >= fabs( mat[1][3] ) ) {
    fprintf( stdout,"'g-', 'LineWidth',2 );\n" );
  }
  else {
    fprintf( stdout,"'k-', 'LineWidth',2 );\n" );
  }


    return 0;

}
