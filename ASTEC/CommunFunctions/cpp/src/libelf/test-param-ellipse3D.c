/*************************************************************************
 * test-param-ellipse3D.c - test de calcul de parametres de granulometrie
 *                          sur des ellipses 3D
 *
 * $Id: test-param-ellipse3D.c,v 1.6 2001/02/13 17:50:52 greg Exp $
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
 * Jul 21 1999
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
#include <time.h>

static int _NO_VERBOSE_ = 0;

typedef enum {
  _ESTIMATE_,
  _COMPUTE_,
  _COMPARE_
} enumMaxFeret;


typedef struct local_par {
  vt_names names;
  int x;
  int y;
  int z;
  int nb;
  double rmin;
  double rmax;
  enumMaxFeret typeComputation;
  int init;
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

static char *usage = "[image-out] [-x %d] [-y %d] [-z %d] [-nb% d]\n\
\t [-exact|-estime|-compare] [-rmin %lf] [-rmax %lf] [-nbest %d]\n\
\t [-v|-nv] [-D] [-help]";

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
  vt_image imres;
  int z;
  double radius[3], angle;
  typeComponentParameters thePar, thePar2;
  int corner[3], size[3];
  
  int init;

  double sumErrorMaxExact = 0.0;
  double sumErrorMedExact = 0.0;
  double sumErrorMinExact = 0.0;

  double sumErrorMaxEstime = 0.0;
  double sumErrorMedEstime = 0.0;
  double sumErrorMinEstime = 0.0;



  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitImage( &imres, par.names.in, par.x, par.y, par.z, USHORT );
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  corner[0] = corner[1] = corner[2] = 1;
  size[0] = imres.dim.x-1;
  size[1] = imres.dim.y-1;
  size[2] = imres.dim.z-1;



  switch ( par.typeComputation ) {
  default :
    break;
  case _COMPUTE_ :
    _SetFeretDiameterComputationToComputation();
    break;
  case _ESTIMATE_ :
    _SetFeretDiameterComputationToEstimation();
    break;
  }
  
  init = par.init;

  for ( z=0; z<par.nb; z++ ) {

    if ( par.init > 0 ) {
      _VerboseInElfParam();
      (void)srandom( par.init );
      par.nb = 1;
    }
    else {
      init = time(0);
      (void)srandom( init );
    }

    _GenerateEllipse3D( &imres, 255, par.rmax, par.rmin, radius, &angle );

    switch ( par.typeComputation ) {
    default :
    case _COMPUTE_ :
    case _ESTIMATE_ :
    _ComputeComponentParameters( &imres, &thePar, (int)-1, (int)255, corner, size );
    if ( _NO_VERBOSE_ == 0 )
      printf("#nb=%3d, rayons reels = (%7.3f %7.3f %7.3f) - estimes = (%7.3f %7.3f %7.3f)\n",
	     z, radius[0], radius[1], radius[2], 
	     thePar.maxDiameter/2.0, thePar.medDiameter/2.0, thePar.minDiameter/2.0 );
    break;

    case _COMPARE_ :
      if ( _NO_VERBOSE_ == 0 )
	printf("#nb=%3d, rayons reels = (%7.3f %7.3f %7.3f) ",
	       z, radius[0], radius[1], radius[2] );
      
      _SetFeretDiameterComputationToComputation();
      _ComputeComponentParameters( &imres, &thePar, (int)-1, (int)255, corner, size );
      sumErrorMaxExact += fabs( thePar.maxDiameter/2.0 - radius[0] ) / radius[0];
      sumErrorMedExact += fabs( thePar.medDiameter/2.0 - radius[1] ) / radius[1];
      sumErrorMinExact += fabs( thePar.minDiameter/2.0 - radius[2] ) / radius[2];
      if ( _NO_VERBOSE_ == 0 )
	printf("- calcules = (%7.3f %7.3f %7.3f) ",
	       thePar.maxDiameter/2.0, thePar.medDiameter/2.0, thePar.minDiameter/2.0 );

      if (0) (void)srandom( 1 );
      _SetFeretDiameterComputationToEstimation();
      _ComputeComponentParameters( &imres, &thePar2, (int)-1, (int)255, corner, size );
      sumErrorMaxEstime += fabs( thePar2.maxDiameter/2.0 - radius[0] ) / radius[0];
      sumErrorMedEstime += fabs( thePar2.medDiameter/2.0 - radius[1] ) / radius[1];
      sumErrorMinEstime += fabs( thePar2.minDiameter/2.0 - radius[2] ) / radius[2];

      if ( _NO_VERBOSE_ == 0 )
	printf("- estimes = (%7.3f %7.3f %7.3f) ",
	       thePar2.maxDiameter/2.0, thePar2.medDiameter/2.0, thePar2.minDiameter/2.0 );

      if ( _NO_VERBOSE_ == 0 )
	printf("\n");

      if ( fabs( thePar2.maxDiameter - thePar.maxDiameter ) > 1e-10 ||
	   fabs( thePar2.medDiameter - thePar.medDiameter ) > 1e-10 ||
	   fabs( thePar2.minDiameter - thePar.minDiameter ) > 1e-10 ) {
	printf("ERROR with INIT=%d\n", init );	
      }

      
    }
    
  }

  switch ( par.typeComputation ) {
  default :
    break;
  case _COMPARE_ :
    fprintf( stderr, " Erreurs (%%) sur rayons calcules Max -> %6.3f   Med -> %6.3f   Min -> %6.3f\n",
	     100.0 * sumErrorMaxExact/par.nb, 100.0 * sumErrorMedExact/par.nb, 
	     100.0 * sumErrorMinExact/par.nb );
    fprintf( stderr, "                        estimes  Max -> %6.3f   Med -> %6.3f   Min -> %6.3f\n",
	     100.0 * sumErrorMaxEstime/par.nb, 100.0 * sumErrorMedEstime/par.nb, 
	     100.0 * sumErrorMinEstime/par.nb );
  }
  



  /*--- ecriture de l'image resultat ---*/
  if ( 0 ) {
    if ( VT_WriteInrimage( &imres ) == -1 ) {
      VT_FreeImage( &imres );
      VT_ErrorParse("unable to write output image\n", 0);
    }
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( &imres );
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
	_NO_VERBOSE_ = 0;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }
      /*--- traitement eventuel de l'image d'entree ---*/

      else if ( strcmp ( argv[i], "-nv" ) == 0 ) {
	_NO_VERBOSE_ = 1;
      }
      /*---  ---*/
      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->x) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->y) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->z) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nb...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nb) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nb...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-rmin" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -rmin...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->rmin) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -rmin...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-rmax" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -rmax...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->rmax) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -rmax...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-exact" ) == 0 ) {
	par->typeComputation = _COMPUTE_;
      }
      else if ( strcmp ( argv[i], "-estime" ) == 0 ) {
	par->typeComputation = _ESTIMATE_;
      }
      else if ( strcmp ( argv[i], "-compare" ) == 0 ) {
	par->typeComputation = _COMPARE_;
      }
      else if ( strcmp ( argv[i], "-nbest" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nbest...\n", 0 );
	status = sscanf( argv[i],"%d",&nb_estimates );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nbest...\n", 0 );
	_SetNbTestsForMaxFeretEstimation( nb_estimates );
      }



      else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->init) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init...\n", 0 );
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
    strcpy( par->names.in,  ">" );  /* standart output */
  }


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
  par->x = 100;
  par->y = 100;
  par->z = 100;
  par->nb = 1;
  par->rmax = 40;
  par->rmin = 10;
  par->typeComputation = _COMPUTE_;
  par->init = -1;
}
