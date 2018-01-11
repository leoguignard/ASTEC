#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--- aide ---*/
static char *program;
static char *usage = "matrice-in [-help]";
static char *detail ="\t Test 4x4  transformation w.r.t. identity ";

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

  double mat[4][4];  
  FILE *file;

  int i, status;
  int nnames;
  char message[256];
  int j, t;

  const double err = 0.001;
  double er = 0.0;
  double et = 0.0;
  double o, e = 0.0;
  int estimate_err = 1;

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
      else if ( strcmp ( argv[i], "-eval" ) == 0 ) {
	estimate_err = 1;
      }      
      else if ( strcmp ( argv[i], "-not-eval" ) == 0 ) {
	estimate_err = 0;
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

  if ( estimate_err == 0 ) {

    t = 0;
    for ( i=0; i<4; i++ )
      for ( j=0; j<4; j++ ) {
	if ( i == j ) {
	  if ( mat[i][j] < 1.0-err || 1.0+err <  mat[i][j] ) t++;
	}
	else {
	  if ( mat[i][j] < (-err) || err <  mat[i][j] ) t++;
	}
      }
    
    if ( t > 0 ) 
      fprintf( stderr, "%s is not equal to identity (e=%g)\n", input_file, err );
    else 
      fprintf( stderr, "%s is equal to identity (e=%g)\n", input_file, err );

  }
  else {

    for ( i=0; i<3; i++ )
      for ( j=0; j<3; j++ ) {
	if ( i == j ) {
	  if ( fabs( 1.0 - mat[i][j] ) > er ) er = fabs( 1.0 - mat[i][j] );
	}
	else {
	  if ( fabs( mat[i][j] ) > er ) er = fabs( mat[i][j] );
	}
      }

    for ( i=0; i<3; i++ ) {
      if ( fabs( mat[i][3] ) > et ) et = fabs( mat[i][3] );
    }

    e = er;
    if ( et > e ) e = et;

    for ( o = 1000.0; o > e; o/=10.0 )
      ;
    o *= 10.0;


    fprintf( stderr, "%s, error max wrt identity = %f\n", input_file, o );
    fprintf( stderr, "     err max rotation = %g\n", er ); 
    fprintf( stderr, "     err max translation = %g\n", et );


  }



  return 0;

}
