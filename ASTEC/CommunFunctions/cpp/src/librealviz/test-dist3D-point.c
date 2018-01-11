/*************************************************************************
 * test-dist3D.c - programme de test pour le calcul de distance
 *
 * On considere un point (voxel) au centre de l'image et on calcule
 * la distance a la surface de ce voxel (non a son centre).
 *
 * $Id:
 * 
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sun May 30 09:52:41 MET DST 1999
 *
 * Copyright (c) INRIA 
 *
 *
 * ADDITIONS, CHANGES:
 *
 * - Sat Jun  5 23:06:57 MET DST 1999
 *   Prise en compte de la structure typeDistanceMap
 *
 */
#include <vt_common.h>
#include <is_distance.h>

typedef enum {
  _ORIGINAL_=0,
  _DISTANCES_ENTIERES_=1,
  _DISTANCES_FLOTTANTES_=2,
  _ERREURS_=3
} typeOutput;

typedef struct local_par {
  vt_names names;
  int dim[3];
  int centre[3];
  float size[3];
  double refine;
  typeOutput output;

  int type;
} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[image-out]\n\
\t [-x %d] [-y %d] [-z %d] [-cx %d] [-cy %d] [-cz %d]\n\
\t [-refine %f] [-ori|-edist|-fdist|-errs] [-vz %f %f %f]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image image, imres, imerr;
  int x, y, z;
  unsigned char ***theBuf;
  short int ***resBuf;
  float ***errBuf;
  typeDistanceMap theDist;
  
  double px=0.0, py=0.0, pz=0.0;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- creation de l'image d'entree ---*/
  VT_InitImage( &image, par.names.in, par.dim[0], par.dim[1], par.dim[2], UCHAR );
  image.siz.x = par.size[0];
  image.siz.y = par.size[1];
  image.siz.z = par.size[2];
  (void)VT_AllocImage( &image );
  theBuf = (unsigned char ***)image.array;
  for ( z=0; z<image.dim.z; z++ )
  for ( y=0; y<image.dim.y; y++ )
  for ( x=0; x<image.dim.x; x++ ) {
    theBuf[z][y][x] = 255;
  }
  if ( par.centre[0] < 0 || par.centre[0] >= par.dim[0] ||
       par.centre[1] < 0 || par.centre[1] >= par.dim[1] ||
       par.centre[2] < 0 || par.centre[2] >= par.dim[2] ) {
    VT_ErrorParse("centre en dehors de l'image\n", 0 );
  }
  
  
  theBuf[par.centre[2]][par.centre[1]][par.centre[0]] = 0;



  if ( par.output == _ORIGINAL_ ) {
      (void)VT_WriteInrimage( &image );
      exit( 0 );
  }


  printf( " dim = %d x %d x %d \n",
	  par.dim[0], par.dim[1], par.dim[2] );
  printf( " voxel = %f x %f x %f mm^3\n",
	  par.size[0], par.size[1], par.size[2] );
  printf( " centre = (%d , %d , %d )\n",
	  par.centre[0], par.centre[1], par.centre[2] );

  /*--- initialisation de l'image resultat ---*/
  VT_InitFromImage( &imres, &image, par.names.in, SSHORT );
  (void)VT_AllocImage( &imres );

  


  /*-- initialisation de la structure carte de distance
   */
  theDist.buf = imres.buf;
  theDist.dim[0] = imres.dim.x;
  theDist.dim[1] = imres.dim.y;
  theDist.dim[2] = imres.dim.z;
  theDist.voxelSize[0] = par.size[0];
  theDist.voxelSize[1] = par.size[1];
  theDist.voxelSize[2] = par.size[2];
  theDist.multiplicativeCoefficient = 0.0;

  /*--- mise a jour de l'interieur et de l'exterieur
        = initialisation de la carte de distance
   */
  _InitSignedDistanceMap( (unsigned char*)image.buf,
		    (short int*)imres.buf,
		    par.dim );

  
  /*-- calcul de la carte
   */
  _ComputeSignedDistanceMap( &theDist, par.refine );


  /*--- sorties
   */

  resBuf = (short int ***)imres.array;
  resBuf[par.centre[2]][par.centre[1]][par.centre[0]] = 0;

  if ( par.output == _DISTANCES_ENTIERES_ ) {
    printf( " Les distances doivent etre multipliees par %f\n", 
	    theDist.multiplicativeCoefficient );
    (void)VT_WriteInrimage( &imres );
    exit( 0 );
  }


  VT_InitFromImage( &imerr, &image, par.names.in, FLOAT );
  (void)VT_AllocImage( &imerr );


  

  errBuf = (float ***)imerr.array;

  switch ( par.output ) {
  default :
  case _ERREURS_ :
  
    for ( z=0; z<image.dim.z; z++ )
    for ( y=0; y<image.dim.y; y++ )
    for ( x=0; x<image.dim.x; x++ ) {
      if ( x < par.centre[0] ) px = (double)par.centre[0] - 0.5;
      else if ( x == par.centre[0] ) px = (double)par.centre[0];
      else px = (double)par.centre[0] + 0.5;
      if ( y < par.centre[1] ) py = (double)par.centre[1] - 0.5;
      else if ( y == par.centre[1] ) py = (double)par.centre[1];
      else py = (double)par.centre[1] + 0.5;
      if ( z < par.centre[2] ) pz = (double)par.centre[2] - 0.5;
      else if ( z == par.centre[2] ) pz = (double)par.centre[2];
      else pz = (double)par.centre[2] + 0.5;

      errBuf[z][y][x] = (double)resBuf[z][y][x] * theDist.multiplicativeCoefficient -
	sqrt( (x-px)*(x-px)*par.size[0]*par.size[0] +
	      (y-py)*(y-py)*par.size[1]*par.size[1] +
	      (z-pz)*(z-pz)*par.size[2]*par.size[2] );
    }
    errBuf[32][32][32] = 0.0;
    break;
    
  case _DISTANCES_FLOTTANTES_ :
    
    for ( z=0; z<image.dim.z; z++ )
    for ( y=0; y<image.dim.y; y++ )
    for ( x=0; x<image.dim.x; x++ ) {
      errBuf[z][y][x] = (double)resBuf[z][y][x] * theDist.multiplicativeCoefficient;
    }

  }

 (void)VT_WriteInrimage( &imerr );


  /*--- liberations memoires ---*/
  return( 0 );
}

#if defined(_ANSI_)
static void VT_Parse( int argc, char *argv[], local_par *par )
#else
  static void VT_Parse( argc, argv, par )
  int argc;
char *argv[];
local_par *par;
#endif
{
  int i, nb, status;
  int o=0, s=0, r=0;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	if ( nb == 0 ) {
	  /*--- standart input ---*/
	  strcpy( par->names.in, "<" );
	  nb += 1;
	}
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }



      else if ( strcmp ( argv[i], "-ori" ) == 0 ) {
	par->output = _ORIGINAL_;
      }
      else if ( strcmp ( argv[i], "-edist" ) == 0 ) {
	par->output = _DISTANCES_ENTIERES_;
      }
      else if ( strcmp ( argv[i], "-fdist" ) == 0 ) {
	par->output = _DISTANCES_FLOTTANTES_;
      }
      else if ( strcmp ( argv[i], "-errs" ) == 0 ) {
	par->output = _ERREURS_;
      }


      else if ( strcmp ( argv[i], "-refine" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -refine...\n", 0 );
	status = sscanf( argv[i],"%lf",&par->refine );
	if ( status <= 0 ) VT_ErrorParse( "parsing -refine...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	status = sscanf( argv[i],"%d",&par->dim[0] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	status = sscanf( argv[i],"%d",&par->dim[1] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	status = sscanf( argv[i],"%d",&par->dim[2] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-cx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cx...\n", 0 );
	status = sscanf( argv[i],"%d",&par->centre[0] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -cx...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-cy" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cy...\n", 0 );
	status = sscanf( argv[i],"%d",&par->centre[1] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -cy...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-cz" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cz...\n", 0 );
	status = sscanf( argv[i],"%d",&par->centre[2] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -cz...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-vz" ) == 0 ) {
	if ( i+3 >= argc)    VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+1], "%f", &par->size[0] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+2], "%f", &par->size[1] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+3], "%f", &par->size[2] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	i += 3;
      }


      /*--- lecture du type de l'image de sortie ---*/
      else if ( strcmp ( argv[i], "-r" ) == 0 ) {
	r = 1;
      }
      else if ( strcmp ( argv[i], "-s" ) == 0 ) {
	s = 1;
      }
      else if ( strcmp ( argv[i], "-o" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
	status = sscanf( argv[i],"%d",&o );
	if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
      }
      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) { 
	strncpy( par->names.in, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
}

#if defined(_ANSI_)
static void VT_ErrorParse( char *str, int flag )
#else
  static void VT_ErrorParse( str, flag )
  char *str;
int flag;
#endif
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}

#if defined(_ANSI_)
static void VT_InitParam( local_par *par )
#else
  static void VT_InitParam( par )
  local_par *par;
#endif
{
  VT_Names( &(par->names) );
  par->dim[0] = par->dim[1] = par->dim[2] = 65;
  par->centre[0] = par->centre[1] = par->centre[2] = 32;
  par->size[0] = par->size[1] = par->size[2] = 1.0;
  par->refine = -1.0;
  par->output = _DISTANCES_ENTIERES_;
  par->type = TYPE_UNKNOWN;
}
