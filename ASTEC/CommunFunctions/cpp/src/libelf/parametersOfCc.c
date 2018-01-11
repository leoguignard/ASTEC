/*************************************************************************
 * parametersOfCc.c - calcul de parametres de granulometrie sur des composantes
 *                    connexes numerotees
 *
 * $Id: parametersOfCc.c,v 1.8 2001/04/04 07:24:58 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * May 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <vt_common.h>
#include <vt_elfparam.h>

#include <ccparameters.h>



typedef enum {
  _ESTIMATE_,
  _COMPUTE_
} enumMaxFeret;


typedef struct local_par {
  vt_names names;
  int type;
  int intervalLength;
  enumMaxFeret typeComputation;
  int methode;
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

static char *usage = "[image-in] [fichier-out]\n\
\t [-exact|-estime] [-nbest %d]\n\
\t [-il %d]\n\
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
  vt_image *image;
  int i, n=0, l=0;
  typeComponentParameters *thePar = (typeComponentParameters *)NULL;
  typeParameter *thePar2 = NULL;
  FILE *fopen(),*f;


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  

  
  (void)srandom(time(0));
  
  switch ( par.typeComputation ) {
  default :
  case _COMPUTE_ :
    _SetFeretDiameterComputationToComputation();
    break;
  case _ESTIMATE_ :
    _SetFeretDiameterComputationToEstimation();
    break;
  }
  
  switch ( par.methode ) {
  default :
  case 1 :
    n =  _CreateArrayOfParametersFromImage( image, -1, &thePar );
    break;
  case 3 :
    n =  _CreateArrayOfParametersFromImage( image, -1, &thePar );
  case 2 :
    {
      int dim[3];
      l = _MaxValueInImage( image, -1 );
      dim[0] = image->dim.x;
      dim[1] = image->dim.y;
      dim[2] = image->dim.z;
      thePar2 = NULL;
      if ( image->type == USHORT ) {
	thePar2 = ComputeParameterFromLabels( (unsigned short int *)image->buf, 
					      dim, NULL, l );
      }
      break;
    }
  }

  VT_FreeImage( image );
  VT_Free( (void**)&image );


  if ( thePar != NULL && thePar2 != NULL ) {
    for ( i=1; i<=n; i++ ) {
      if ( thePar[i].volume != thePar2[i].volume ||
	   fabs( thePar[i].maxDiameter - thePar2[i].maxDiameter ) > 0.0001 ||
	   fabs( thePar[i].medDiameter - thePar2[i].medDiameter ) > 0.0001 ||
	   fabs( thePar[i].minDiameter - thePar2[i].minDiameter ) > 0.0001 ) {
	printf( "#%4d REF = %8d (%7.4f %7.4f %7.4f)\n", i,
		thePar[i].volume, thePar[i].maxDiameter,
		thePar[i].medDiameter, thePar[i].minDiameter );
	printf( "      2nd = %8d (%7.4f %7.4f %7.4f)\n",
		thePar2[i].volume, thePar2[i].maxDiameter,
		thePar2[i].medDiameter, thePar2[i].minDiameter );
	
      }
    }
    fprintf( stderr, "verification effectuee\n" );
  }


  switch ( par.methode ) {
  case 1 :
    if ( n <= 0 ) {
      VT_ErrorParse( "unable to compute components parameters", 0 );
    }
    
    if ( _VT_VERBOSE_ ) {
      fprintf( stdout, "%s : %d components\n", par.names.in, n );
      
      for ( i=1; i<=n; i++ ) {
	fprintf( stdout, "component #%5d: volume = %d\n", i, thePar[i].volume );
	fprintf( stdout, "                  border = %d\n", thePar[i].border );
	fprintf( stdout, "                  rayon equivalent (volume) = %f\n", 
		 thePar[i].equivSphereRadiusFromVolume );
	fprintf( stdout, "                  rayon equivalent (surface) = %f\n", 
		 thePar[i].equivSphereRadiusFromSurface );
	fprintf( stdout, "                  coin inf = (%d %d %d)\n",
		 thePar[i].ptmin[0], thePar[i].ptmin[1], thePar[i].ptmin[2] );
	fprintf( stdout, "                  coin sup = (%d %d %d)\n",
		 thePar[i].ptmax[0], thePar[i].ptmax[1], thePar[i].ptmax[2] );
	fprintf( stdout, "                  diametre Feret max = %f (rayon = %6.3f)\n", 
		 thePar[i].maxDiameter, thePar[i].maxDiameter/2.0 );
	fprintf( stdout, "                  direction max      = (%5.2f %5.2f %5.2f)\n", 
		 thePar[i].maxDirection[0], thePar[i].maxDirection[1], 
		 thePar[i].maxDirection[2] );
	fprintf( stdout, "                  diametre Feret med = %f (rayon = %6.3f)\n", 
		 thePar[i].medDiameter, thePar[i].medDiameter/2.0 );
	fprintf( stdout, "                  direction med      = (%5.2f %5.2f %5.2f)\n", 
		 thePar[i].medDirection[0], thePar[i].medDirection[1], 
		 thePar[i].medDirection[2] );
	fprintf( stdout, "                  diametre Feret min = %f (rayon = %6.3f)\n", 
		 thePar[i].minDiameter, thePar[i].minDiameter/2.0 );
	fprintf( stdout, "                  direction min      = (%5.2f %5.2f %5.2f)\n", 
		 thePar[i].minDirection[0], thePar[i].minDirection[1], 
		 thePar[i].minDirection[2] );
      }
    }
    
    if ( par.names.out[0] != '\0' && par.names.out[0] != '>') {
      
      f = fopen( par.names.out, "w" );
      fprintf( f, "#" );
      for ( i=0; i<argc; i++ )
	fprintf( f, " %s", argv[i] );
      fprintf( f, "\n" );
      fprintf( f, "# volume border rayon-equiv-sphere-vol rayon-equiv-sphere-surf max-Feret med-Feret min-Feret\n");
      for ( i=1; i<=n; i++ ) 
	fprintf( f, "%9.1f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n", (float)thePar[i].volume,
		 (float)thePar[i].border,
		 thePar[i].equivSphereRadiusFromVolume, thePar[i].equivSphereRadiusFromSurface, 
		 thePar[i].maxDiameter,
		 thePar[i].medDiameter, thePar[i].minDiameter );
      fclose( f );
      
      /*
	(void)sprintf( filename, "%s.volumes", par.names.out );
	f = fopen( filename, "w" );
	for ( i=1; i<=n; i++ ) 
	fprintf( f, "%d\n", thePar[i].volume );
	fclose( f );
	
	(void)sprintf( filename, "%s.rayons", par.names.out );
	f = fopen( filename, "w" );
	for ( i=1; i<=n; i++ ) 
	fprintf( f, "%f\n", thePar[i].equivSphereRadius );
	fclose( f );
      */
    }
    break;

  case 2 :
    break;

  }


  free( thePar );
  return( 1 );
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
  int nb_estimates;
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

      else if ( strcmp ( argv[i], "-il" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -il...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->intervalLength) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -il...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-m" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -m...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->methode) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -m...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-exact" ) == 0 ) {
	par->typeComputation = _COMPUTE_;
      }
      else if ( strcmp ( argv[i], "-estime" ) == 0 ) {
	par->typeComputation = _ESTIMATE_;
      }
      else if ( strcmp ( argv[i], "-nbest" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nbest...\n", 0 );
	status = sscanf( argv[i],"%d",&nb_estimates );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nbest...\n", 0 );
	_SetNbTestsForMaxFeretEstimation( nb_estimates );
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
      else if ( nb == 1 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
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
  par->type = TYPE_UNKNOWN;
  par->intervalLength = 100;
  par->typeComputation = _ESTIMATE_;
  par->methode = 1;
}
