#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*--- aide ---*/
static char *program;
static char *usage = "[-fr %d] [-lr %d] [-fp %d] [-lp %d]\n\
[-id %d] [-ip %s] [-is %s] [-od %s] [-op %s] [-os %s] [output]\n\
\t[-it %d] [-error %lf] [-s %lf] [-affine | -rigid ] [-2D | -3D]";
static char *detail ="\n Smooth transformations\n\n\
-first_to_be_processed | -fp %d: slice number of the first slice\n\
   to be processed\n\
-last_to_be_processed  | -lp %d:\n\
-first_to_be_read | -fr %d:\n\
-last_to_be_read  | -lr %d:\n\
-input_dir    | -id %s:\n\
-input_prefix | -ip %s:\n\
-input_suffix | -is %s: input file names are composed by\n\
   'input_dir'/'input_prefix''slice number''input_suffix'\n\
-output_dir    | -od %s:\n\
-output_prefix | -op %s:\n\
-output_suffix | -os %s: \n\
   if 'output' is not given, output file names are composed by\n\
   'output_dir'/'output_prefix''slice number''output_suffix'\n\
   if 'output' is given, all the results will have the same name\n\
   'output'\n\
-iterations | -it %d: maximal number of iterations for convergence\n\
-error %lf: maximal error for convergence\n\
-sigma | -s %lf\n\
-rigid:\n\
-affine:\n\
-2D:\n\
-3D\n\
";


/******************************************

*******************************************/

#include "mat-manip.c"
#include "mat-rot.c"

typedef struct {
  double rot[3];
  double trs[3];
} typeTrsfRigid;

void Matrix_To_TrsfRigid( double R[4][4], typeTrsfRigid *trsf, int dim )
{
  double theta;

  if ( dim == 2 ) {

    theta = acos( R[0][0] );
    if ( R[1][0] < 0.0 ) theta *= -1.0;
    
    trsf->trs[0] = R[0][3];
    trsf->trs[1] = R[1][3];
    trsf->trs[2] = 0.0;
    trsf->rot[0] = 0.0;
    trsf->rot[1] = 0.0;
    trsf->rot[2] = theta;

    return;
  }

  trsf->trs[0] = R[0][3];
  trsf->trs[1] = R[1][3];
  trsf->trs[2] = R[2][3];
  Rmatrix_To_Rvector( R, trsf->rot );
}

void TrsfRigid_To_Matrix( double R[4][4], typeTrsfRigid *trsf )
{
  Rvector_To_Rmatrix( R, trsf->rot );
  R[0][3] = trsf->trs[0];
  R[1][3] = trsf->trs[1];
  R[2][3] = trsf->trs[2];
  R[3][0] = R[3][1] = R[3][2] = 0.0;
  R[3][3] = 1.0;
}






/******************************************

*******************************************/

static void _errorMessage( char *str, int flag )
{
  fprintf( stderr," Usage: %s %s\n", program, usage );
  if ( flag == 1 )
    fprintf( stderr," %s", detail );
  if ( str != NULL )
    fprintf( stderr," Error: %s", str );
  exit( -1 );
}





typedef enum {
  _AFFINE_,
  _RIGID_ 
} enumTransformationType;


/* Programme principal */
#define NBTRSF 1500

static int _verbose_ = 1;

int main(int argc, char *argv[])
{
  char *def_prefixe = "s";
  char *def_suffixe = ".rigid.trsf";
  char *prefixe = def_prefixe;
  char *suffixe = def_suffixe;
  char *input_dir = NULL;
  
  char *def_oprefixe = "s";
  char *def_osuffixe = ".smooth.rigid.trsf";
  char *oprefixe = def_oprefixe;
  char *osuffixe = def_osuffixe;
  char *output_dir = NULL;

  char *output_file = NULL;

  int first_to_be_read = 132;
  int last_to_be_read  = 416;
  int first_to_be_processed = 132;
  int last_to_be_processed = 416;
  double sigma = 1.0;

  enumTransformationType transformationType = _AFFINE_;
  int dimensionality = 3;

  int nnames, i, j;
  char message[256];
  int status;

  char mat_name[256];
  FILE *file;
  double Trsfs[NBTRSF][4][4];
  typeTrsfRigid TrsfsRigid[NBTRSF];
  typeTrsfRigid tmpRigid;

  double Resid[NBTRSF][4][4];

  double inv[4][4], tmp[4][4];  
  double oldres[4][4];
  double newres[4][4];

  int w, wf, wl;
  double c[NBTRSF], sum;
  int a, b;

  int iterations;
  int max_iterations = 10;

  double error_rot;
  double error_trs;
  double error_aff;
  double error;

  double error_norm_rot = 3.1415926535897932;
  double error_norm_aff = 3.1415926535897932;
  double error_norm_trs = 256;

  double max_error = 1e-7;


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
      else if ( strcmp ( argv[i], "first_to_be_processed" )  == 0 ||
		strcmp ( argv[i], "-fp" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -first_to_be_processed ...\n", 0 );
	status = sscanf( argv[i],"%d",&first_to_be_processed );
	if ( status != 1 ) _errorMessage( "parsing -first_to_be_processed ...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-last_to_be_processed" )  == 0 ||
		strcmp ( argv[i], "-lp" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -last_to_be_processed ...\n", 0 );
	status = sscanf( argv[i],"%d",&last_to_be_processed );
	if ( status != 1 ) _errorMessage( "parsing -last_to_be_processed ...\n", 0 );
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



      else if ( strcmp ( argv[i], "-iterations" )  == 0 ||
		strcmp ( argv[i], "-it" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -iterations ...\n", 0 );
	status = sscanf( argv[i], "%d", &max_iterations );
	if ( status != 1 ) _errorMessage( "parsing -iterations ...\n", 0 );
	if ( max_iterations < 0 ) max_iterations = 0;
      }
      


      else if ( strcmp ( argv[i], "-error" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -error ...\n", 0 );
	status = sscanf( argv[i], "%lf", &max_error );
	if ( status != 1 ) _errorMessage( "parsing -error ...\n", 0 );
	if ( max_iterations < 0 ) max_iterations = 0;
      }
      


      else if ( strcmp ( argv[i], "-sigma" )  == 0 ||
		strcmp ( argv[i], "-s" )  == 0 ) {
	i++;
	if ( i >= argc )   _errorMessage( "parsing -sigma ...\n", 0 );
	status = sscanf( argv[i],"%lf",&sigma );
	if ( status != 1 ) _errorMessage( "parsing -sigma ...\n", 0 );
      }
      
      else if ( strcmp ( argv[i], "-2D" )  == 0 ) {
	dimensionality = 2;
      }
      else if ( strcmp ( argv[i], "-3D" )  == 0 ) {
	dimensionality = 3;
      }
      else if ( strcmp ( argv[i], "-rigid" )  == 0 ) {
	transformationType = _RIGID_;
      }
      else if ( strcmp ( argv[i], "-affine" )  == 0 ) {
	transformationType = _AFFINE_;
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
	output_file= argv[i];
	nnames ++;
      }
      else {
	_errorMessage( "too much file names when parsing\n", 0 );
      }
    }
  }


    
  /* lecture des transformations
   */
  for ( i=first_to_be_read; i<=last_to_be_read; i++ ) {

    if ( input_dir == NULL )
      sprintf( mat_name, "%s%d%s", prefixe, i, suffixe );
    else 
      sprintf( mat_name, "%s/%s%d%s", input_dir, prefixe, i, suffixe );

    file = fopen(mat_name, "r");
    if ( file == NULL ) {
      sprintf( message, "Erreur a l'ouverture de '%s'\n", mat_name );
      _errorMessage( message, 0 );
    }

    if ( _verbose_ ) fprintf(stderr, "... Lecture en cours de  %s\r", mat_name );

    status = fscanf(file , "(\nO8\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf\n%lf %lf %lf %lf", 
		    &Trsfs[i][0][0], &Trsfs[i][0][1], &Trsfs[i][0][2], &Trsfs[i][0][3], 
		    &Trsfs[i][1][0], &Trsfs[i][1][1], &Trsfs[i][1][2], &Trsfs[i][1][3], 
		    &Trsfs[i][2][0], &Trsfs[i][2][1], &Trsfs[i][2][2], &Trsfs[i][2][3], 
		    &Trsfs[i][3][0], &Trsfs[i][3][1], &Trsfs[i][3][2], &Trsfs[i][3][3]);
    if ( status != 16 ){
      sprintf( message, "Erreur a la lecture de '%s'\n", mat_name );
      _errorMessage( message, 0 );
    }
    fclose( file );
  }
  if ( _verbose_ ) fprintf(stderr, "\n" );

  w = (int)(3 * sigma) + 1;
  if ( w < 0 ) w = 1;
  

  for ( i=first_to_be_processed; i<=last_to_be_processed; i++ ) {
    /*
     * calcul du masque
     */
    wf = i - w;   if ( wf < first_to_be_read ) wf = first_to_be_read;
    wl = i + w;   if ( wl > last_to_be_read  ) wl = last_to_be_read;
    /*
     * fprintf( stderr, " lissage de %d sur [ %d %d ]\n", i, wf, wl );
     */
    for ( sum=0.0, j=wf; j<=wl; j++ ) {
      c[j-wf] = exp( - (j-i) * (j-i) / ( 2 * sigma * sigma ) );
      sum += c[j-wf];
    }
    for ( j=wf; j<=wl; j++ ) c[j-wf] /= sum;



    /* 
     * estimation de la moyenne => la matrice courante
     */
    MatCopy( Trsfs[i], newres );

    iterations = 0;

    do {

      iterations ++;

      /* copie de l'ancienne estimation
       */
      MatCopy( newres, oldres );
      /* inversion de la transformation
       */
      MatInverse( oldres, inv );
      /* calcul des residus
       * multiplication a gauche par l'inverse
       */
      for ( j=wf; j<=wl; j++ ) MatMult( inv, Trsfs[j], Resid[j-wf] );

      /*
       * lissage des coefficients
       */
      if ( transformationType == _RIGID_ ) {

	for ( j=wf; j<=wl; j++ ) 
	  Matrix_To_TrsfRigid( Resid[j-wf], &TrsfsRigid[j-wf], dimensionality );

	for ( a=0; a<3; a++ ) {
	  tmpRigid.rot[a] = tmpRigid.trs[a] = 0;
	  for ( j=wf; j<=wl; j++ ) {
	    tmpRigid.rot[a] += c[j-wf] * TrsfsRigid[j-wf].rot[a];
	    tmpRigid.trs[a] += c[j-wf] * TrsfsRigid[j-wf].trs[a];
	  }
	}
	
	TrsfRigid_To_Matrix( tmp, &tmpRigid );
	
	for ( error_rot=0.0, error_trs=0.0, a=0; a<3; a++ ) {
	  error_rot += tmpRigid.rot[a] * tmpRigid.rot[a];
	  error_trs += tmpRigid.trs[a] * tmpRigid.trs[a];
	}
	error = sqrt( error_rot / ( error_norm_rot * error_norm_rot )
		      +  error_trs  / ( error_norm_trs * error_norm_trs ) );

      }
      else {
	
	for ( a=0; a<4; a++ )
	for ( b=0; b<4; b++ ) {
	  tmp[a][b] = 0.0;
	  for ( j=wf; j<=wl; j++ ) tmp[a][b] += c[j-wf] * Resid[j-wf][a][b];
	}
	
	for ( error_aff = 0.0, a=0; a<3; a++ )
	for ( b=0; b<3; b++ ) {
	  if ( a == b ) error_aff += ( tmp[a][b] - 1 ) * ( tmp[a][b] - 1 );
	  else          error_aff += tmp[a][b] * tmp[a][b];
	}
	for ( error_trs=0.0, a=0; a<3; a++ )
	  error_trs += tmp[a][3] * tmp[a][3];
	error = sqrt( error_aff / ( error_norm_aff * error_norm_aff )
		      +  error_trs  / ( error_norm_trs * error_norm_trs ) );
	
      }

      /*
       * calcul de la matrice resultat, multiplication a gauche par la transfo 
       */
      MatMult( oldres, tmp, newres );
      
    } while ( error > max_error && iterations < max_iterations );
    

    if ( _verbose_ ) {
      fprintf( stderr, "#%3d error = %f, iterations = %d / %d\n", 
	       i, error, iterations, max_iterations );
    }


    /*
     * ecriture du resultat
     */
    if ( output_file != NULL ) {
      sprintf( mat_name, "%s", output_file );
    }
    else {
      if ( output_dir == NULL )
	sprintf( mat_name, "%s%d%s", oprefixe, i, osuffixe );
      else 
	sprintf( mat_name, "%s/%s%d%s", output_dir, oprefixe, i, osuffixe );
    }
    
    file = fopen( mat_name, "w");
    if ( file == NULL ) {
      sprintf( message, "Erreur a l'ouverture de '%s'\n", mat_name );
      _errorMessage( message, 0 );
    }
    
    fprintf(file, "(\nO8\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n%f %f %f %f\n)\n",
	    newres[0][0],  newres[0][1], newres[0][2], newres[0][3], 
	    newres[1][0],  newres[1][1], newres[1][2], newres[1][3], 
	    newres[2][0],  newres[2][1], newres[2][2], newres[2][3], 
	    newres[3][0],  newres[3][1], newres[3][2], newres[3][3]);
    fclose(file);
  }



    return 0;
}

