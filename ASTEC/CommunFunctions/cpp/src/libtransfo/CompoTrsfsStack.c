#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--- aide ---*/
static char *program;
static char *usage = "[-fr %d] [-lr %d] [-n %d] [-ref %d]\n\
[-id %d] [-ip %s] [-isep %s] [-is %s] [-od %s] [-op %s] [-os %s]\n\
[-inv input] output";
static char *detail ="\n Compose successive transformations\n\n\
-first_to_be_read | -fr %d:\n\
-last_to_be_read  | -lr %d:\n\
   see -input_suffix\n\
-n %d: number of slices\n\
   equivalent to '-first_to_be_read 1 -last_to_be_read %d'\n\
   see -input_suffix\n\
-ref %s: index of the reference slice\n\
-input_dir    | -id %s:\n\
-input_prefix | -ip %s:\n\
-input_separator | -isep %s:\n\
-input_suffix | -is %s: input file names are composed by\n\
   'input_dir'/'input_prefix'N'input_separator'N-1'input_suffix'\n\
   or\n\
   'input_dir'/'input_prefix'N-1'input_separator'N'input_suffix'\n\
   with N from n to 2\n\
   or with N from first_to_be_read+1 to last_to_be_read\n\
-output_dir    | -od %s:\n\
-output_prefix | -op %s:\n\
-output_suffix | -os %s: \n\
   if 'output' is not given, output file names are composed by\n\
   'output_dir'/'output_prefix'N'input_separator'ref'output_suffix'\n\
   with N from 1 to n\n\
   if 'output' is given, all the results will have the same name\n\
   'output'\n\
";

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
typedef double mat4x4[4][4];

#include "mat-manip.c"



typedef enum {
  _INVERSE_,
  _COMPOSE_ 
} enumComputation;


/* Programme principal */
#define NBTRSF 1500

static int _verbose_ = 1;

int main(int argc, char *argv[])
{
  char *def_prefixe = "p_";
  char *def_separator = "-ON-";
  char *def_suffixe = ".trsf";

  char *prefixe = def_prefixe;
  char *separator = def_separator;
  char *suffixe = def_suffixe;
  char *input_dir = NULL;
  
  char *oprefixe = def_prefixe;
  char *osuffixe = def_suffixe;
  char *output_dir = NULL;

  char *input_file  = NULL;
  char *output_file = NULL;

  enumComputation typeComputation = _COMPOSE_;

  int first_to_be_read = 1;
  int last_to_be_read  = -1;

  int refslice = -1;

  int nnames, i;
  char message[256];
  int status;

  char mat1_name[256];
  char mat2_name[256];
  int mat1_read, mat2_read;

  double Trsfs[NBTRSF][4][4];
  double invTrsfs[NBTRSF][4][4];
  double resTrsfs[NBTRSF][4][4];

  
  double inv[4][4], mat[4][4];  


  
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


      else if ( strcmp ( argv[i], "-first_to_be_read" )  == 0 ||
		strcmp ( argv[i], "-fr" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -first_to_be_read ...\n", 0 );
	status = sscanf( argv[i],"%d",&first_to_be_read );
	if ( status != 1 ) _errorMessage( "parsing -first_to_be_read ...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-last_to_be_read" )  == 0 ||
		strcmp ( argv[i], "-lr" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -last_to_be_read ...\n", 0 );
	status = sscanf( argv[i],"%d",&last_to_be_read );
	if ( status != 1 ) _errorMessage( "parsing -last_to_be_read ...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-n" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -n ...\n", 0 );
	status = sscanf( argv[i],"%d",&last_to_be_read );
	if ( status != 1 ) _errorMessage( "parsing -n ...\n", 0 );
	first_to_be_read = 1;
      }
      
      else if ( strcmp ( argv[i], "-ref" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -ref ...\n", 0 );
	status = sscanf( argv[i],"%d",&refslice );
	if ( status != 1 ) _errorMessage( "parsing -ref ...\n", 0 );
      }
      
      else if ( strcmp ( argv[i], "-input_dir" )  == 0 ||
		strcmp ( argv[i], "-id" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -input_dir ...\n", 0 );
	input_dir = argv[i];
      }
      else if ( strcmp ( argv[i], "-input_suffix" )  == 0 ||
		strcmp ( argv[i], "-is" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -input_suffix ...\n", 0 );
	suffixe = argv[i];
      }
      else if ( strcmp ( argv[i], "-input_separator" )  == 0 ||
		strcmp ( argv[i], "-isep" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -input_separator ...\n", 0 );
	separator = argv[i];
      }
      else if ( strcmp ( argv[i], "-input_prefix" )  == 0 ||
		strcmp ( argv[i], "-ip" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -input_prefix ...\n", 0 );
	prefixe = argv[i];
      }
      else if ( strcmp ( argv[i], "-output_dir" )  == 0 ||
		strcmp ( argv[i], "-od" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -output_dir ...\n", 0 );
	output_dir = argv[i];
      }
      else if ( strcmp ( argv[i], "-output_suffix" )  == 0 ||
		strcmp ( argv[i], "-os" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -output_suffix ...\n", 0 );
	osuffixe = argv[i];
      }
      else if ( strcmp ( argv[i], "-output_prefix" )  == 0 ||
		strcmp ( argv[i], "-op" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -output_prefix ...\n", 0 );
	oprefixe = argv[i];
      }


      else if ( strcmp ( argv[i], "-inv" )  == 0 ) {
	typeComputation = _INVERSE_;
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







  if ( typeComputation == _INVERSE_ ) {

    if ( MatRead( input_file, mat ) != 1 ) {
      sprintf( message, "Error while opening/reading '%s'\n", input_file );
      _errorMessage( message, 0 );
    }
    
    MatInverse( mat, inv );
    
    if ( output_file != NULL ) {
      if ( MatWrite( output_file, inv ) != 1 ) 
	sprintf( message, "Error while opening/writing '%s'\n", output_file );
      _errorMessage( message, 0 );
    }
    else {
      MatPrintf( stdout, inv );
    }

    return 0;
  
  }
  











  if ( refslice < 1 ) refslice = last_to_be_read / 2;


    
  /* 
   * lecture des transformations
   * TRSF[i] contains t-[i]-ON-[i-1]
   */
  for ( i=first_to_be_read+1; i<=last_to_be_read; i++ ) {

    if ( input_dir == NULL ) {
      sprintf( mat1_name, "%s%d%s%d%s", prefixe, i, separator, i-1, suffixe );
      sprintf( mat2_name, "%s%d%s%d%s", prefixe, i-1, separator, i, suffixe );
    }
    else {
      sprintf( mat1_name, "%s/%s%d%s%d%s", input_dir, prefixe, i, separator, i-1, suffixe );
      sprintf( mat2_name, "%s/%s%d%s%d%s", input_dir, prefixe, i-1, separator, i, suffixe );
    }

    mat1_read = MatRead( mat1_name, Trsfs[i] );
    mat2_read = MatRead( mat2_name, invTrsfs[i-1] );
    
    if ( mat1_read == 1 && mat2_read == 1 ) {
      sprintf( message, "Both '%s' and '%s' are present ...\n", mat1_name, mat2_name );
      _errorMessage( message, 0 );
    }
    if ( mat1_read != 1 && mat2_read != 1 ) {
      sprintf( message, "Neither '%s' nor '%s' are present ...\n", mat1_name, mat2_name );
      _errorMessage( message, 0 );
    }

    if ( mat1_read == 1 ) {
      if ( _verbose_ ) fprintf(stderr, "... has read '%s'\n", mat1_name );
      fprintf(stderr, "rank = %d\n", MatInverse( Trsfs[i], invTrsfs[i-1] ) );
      MatPrintf( stdout, Trsfs[i] );
      MatPrintf( stdout, invTrsfs[i] );
    }
    else {
      if ( _verbose_ ) fprintf(stderr, "... has read '%s'\n", mat2_name );
      MatInverse( invTrsfs[i-1], Trsfs[i] );
    }
  }
  if ( _verbose_ ) fprintf(stderr, "\n" );
  

  /* 
   * TRSF[i] contains t-[i]-ON-[i-1]
   * INVTRSF[i] contains t-[i]-ON-[i+1]
   */



  MatIdentity( resTrsfs[refslice] );


  if ( 0 ) MatPrintf( stderr, resTrsfs[refslice] );

  if ( refslice > 1 ) {
    /* on met t-[refslice-1]-ON-[refslice] dans resTrsfs[ refslice-1 ]
     */
    MatCopy( invTrsfs[ refslice-1 ], resTrsfs[ refslice-1 ] );
    /* on met 
       t-[i+1]-ON-[refslice] * t-[i]-ON-[i+1] = t-[i]-ON-[refslice]
       dans resTrsfs[i]
    */
    for ( i=refslice-2; i>=1; i-- ) {
      MatMult( resTrsfs[i+1], invTrsfs[i], resTrsfs[i] );
    }
  }
  
  if ( refslice < last_to_be_read ) {
    /* on met t-[refslice+1]-ON-[refslice] dans resTrsfs[ refslice+1 ]
     */
    MatCopy( Trsfs[ refslice+1 ], resTrsfs[ refslice+1 ] );
    for ( i=refslice+2; i<= last_to_be_read; i++ ) {
      /* on met 
	 t-[i-1]-ON-[refslice] * t-[i]-ON-[i-1] = t-[i]-ON-[refslice]
	 dans resTrsfs[i]
      */
      MatMult( resTrsfs[i-1], Trsfs[i], resTrsfs[i] );
    }
  }

  if ( 0 ) MatPrintf( stderr, resTrsfs[refslice] );

  for ( i=first_to_be_read; i<=last_to_be_read; i++ ) {
    if ( output_dir == NULL )
      sprintf( mat1_name, "%s%d%s%d%s", oprefixe, i, separator, refslice, suffixe );
    else 
      sprintf( mat1_name, "%s/%s%d%s%d%s", output_dir, oprefixe, i, separator, refslice, osuffixe );
    
    if ( MatWrite( mat1_name, resTrsfs[i] ) != 1 ) {
      sprintf( message, "Error while opening/writing '%s'\n", output_file );
      _errorMessage( message, 0 );
    }
    if ( _verbose_ ) fprintf(stderr, "... has written '%s'\n", mat1_name );

  }
  if ( _verbose_ ) fprintf(stderr, "\n" );
  
  return 0;
}

