#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--- aide ---*/
static char *program;
static char *usage = "matrix-in [-right-diag %lf %lf %lf] [-left-diag %lf %lf %lf] [matrix-out] [-help]";
static char *detail ="\t Multiply the 4x4 matrix 'matrix-in'\n\
 at the right by the diagonal matrix 'R' specified by '-diag-right'\n\
 at the left by the diagonal matrix 'L' specified by '-diag-left'\n\
 The result if then: 'L * matrix-in * R'\n\
-right-diag | -right | -rd \n\
 specifies the 3 first elements of the right diagonal matrix\n\
 (the fourth one being 1)\n\
-left-diag | -left | -ld \n\
 specifies the 3 first elements of the left diagonal matrix\n\
 (the fourth one being 1)\n";
  


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




/* Programme principal */

int main( int argc, char *argv[] )
{
  char *input_file = NULL;
  char *output_file = NULL;

  double r[3] = { 1.0, 1.0, 1.0 };
  double l[3] = { 1.0, 1.0, 1.0 };
  double mat[4][4];  
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

      if ( strcmp ( argv[i], "-right-diag" ) == 0 
	   || strcmp ( argv[i], "-right" ) == 0 
	   || strcmp ( argv[i], "-rd" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -right-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(r[0]) );
	if ( status <= 0 ) _errorMessage( "parsing -right-diag #1... \n", 1 );
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -right-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(r[1]) );
	if ( status <= 0 ) _errorMessage( "parsing -right-diag #2... \n", 1 );
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -right-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(r[2]) );
	if ( status <= 0 ) _errorMessage( "parsing -right-diag #3... \n", 1 );
      }

      else if ( strcmp ( argv[i], "-left-diag" ) == 0 
	   || strcmp ( argv[i], "-left" ) == 0 
	   || strcmp ( argv[i], "-ld" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -left-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(l[0]) );
	if ( status <= 0 ) _errorMessage( "parsing -left-diag #1... \n", 1 );
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -left-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(l[1]) );
	if ( status <= 0 ) _errorMessage( "parsing -left-diag #2... \n", 1 );
	i += 1;
	if ( i >= argc)    _errorMessage( "parsing -left-diag ... \n", 1 );
	status = sscanf( argv[i], "%lf", &(l[2]) );
	if ( status <= 0 ) _errorMessage( "parsing -left-diag #3... \n", 1 );
      }

      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
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
      else if ( nnames == 1 ) {
	output_file= argv[i];
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

    
  for ( i=0; i<4; i++ ) {
    mat[i][0] *= r[0];
    mat[i][1] *= r[1];
    mat[i][2] *= r[2];
  }

  for ( i=0; i<4; i++ ) {
    mat[0][i] *= l[0];
    mat[1][i] *= l[1];
    mat[2][i] *= l[2];
  }

  
  
  if ( output_file != NULL ) {
    file = fopen( output_file, "w" );
    if ( file == NULL ) {
      sprintf( message, "Erreur a l'ouverture de '%s'\n", output_file );
      _errorMessage( message, 0 );
    }
  }
  else 
    file = stdout;
  
  fprintf(file, "(\nO8\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n)\n",
	  mat[0][0],  mat[0][1], mat[0][2], mat[0][3], 
	  mat[1][0],  mat[1][1], mat[1][2], mat[1][3], 
	  mat[2][0],  mat[2][1], mat[2][2], mat[2][3], 
	  mat[3][0],  mat[3][1], mat[3][2], mat[3][3]);
  
  
  if ( output_file != NULL )
    fclose( file );



    return 0;

}
