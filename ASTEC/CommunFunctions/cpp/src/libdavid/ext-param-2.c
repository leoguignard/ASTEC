/*************************************************************************
 * minimum.c -
 *
 * $Id: minimum.c,v 1.5 2000/08/16 16:31:56 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_histo.h>
#include <ccparameters.h>

typedef struct local_par {
  vt_names names;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-1] [image-2] [image-out]\n\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.5 $ $Date: 2000/08/16 16:31:56 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image1, *image2;
  typeImageParameter thePar1;
  typeImageParameter thePar2;

  char fname[STRINGLENGTH];
  int stats1[_CONFIGURATIONS_];
  int stats2[_CONFIGURATIONS_];

  FILE *f;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #1\n", 0);
  
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #2\n", 0);
  


  thePar1.i_image = 1;
  thePar1.name = par.names.in;

  thePar2.i_image = 2;
  thePar2.name = par.names.ext;

  thePar1.thePar = ComputeParameterFromLabels( image1, &(thePar1.n_labels) );
  thePar2.thePar = ComputeParameterFromLabels( image2, &(thePar2.n_labels) );

  fprintf(stderr,"found %d connect components in '%s'\n",
	  thePar1.n_labels, thePar1.name );
  fprintf(stderr,"found %d connect components in '%s'\n",
	  thePar2.n_labels, thePar2.name );

  _FollowsComponents( image1, &thePar1, image2, &thePar2 );

  if ( par.names.out[0] != '\0' && par.names.out[0] != '>' ) {
    sprintf( fname, "%s.grf", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open grf output file\n", 0 );
    _print_AllComponents( f, &thePar1, &thePar2 );
    fclose( f );
  }

  _NameComponents( image1, &thePar1, image2, &thePar2 );
  _NameComponents( image2, &thePar2, image1, &thePar1 );
   
  if ( par.names.out[0] == '\0' || par.names.out[0] == '>' ) {
    
    _print_Statistics( stdout, &thePar1 );
    _print_Statistics( stdout, &thePar2 );
    _print_UnknownComponents( stdout, &thePar1, &thePar2 );
    
    return( 1 );
  }

  sprintf( fname, "%s.log", par.names.out );
  f = fopen( fname, "w" );
  if ( f == NULL ) VT_ErrorParse("unable to open log output file\n", 0 );
  _print_Statistics( f, &thePar1 );
  _print_Statistics( f, &thePar2 );
  fprintf( f, "\n\n\n\n\n" );
  _print_AllComponents( f, &thePar1, &thePar2 );
  fclose( f );
  
  _get_Statistics( stats1, &thePar1 );
  _get_Statistics( stats2, &thePar2 );
  
  if ( stats1[_ISOLE_] > 0 ) {
    sprintf( fname, "%s.01", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open 01 output file\n", 0 );
    _print_same_components( f, _ISOLE_, &thePar1, &thePar2 );
    fclose( f );
  }

  if ( stats2[_ISOLE_] > 0 ) {
    sprintf( fname, "%s.02", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open 02 output file\n", 0 );
    _print_same_components( f, _ISOLE_, &thePar2, &thePar1 );
    fclose( f );
  }
  
  if ( stats1[_UNIVOQUE_] > 0 ) {
    sprintf( fname, "%s.1", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open 1 output file\n", 0 );
    _print_same_components( f, _UNIVOQUE_, &thePar1, &thePar2 );
    fclose( f );
  }

  if ( stats1[_ZIGZAG_] > 0 ) {
    sprintf( fname, "%s.1z", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open 1z output file\n", 0 );
    _print_same_components( f, _ZIGZAG_, &thePar1, &thePar2 );
    fclose( f );
  }

  if ( stats1[_SPLIT_] > 0 ) {
    sprintf( fname, "%s.N1", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open N1 output file\n", 0 );
    _print_same_components( f, _SPLIT_, &thePar1, &thePar2 );
    fclose( f );
  }

  if ( stats2[_SPLIT_] > 0 ) {
    sprintf( fname, "%s.N2", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open N2 output file\n", 0 );
    _print_same_components( f, _SPLIT_, &thePar2, &thePar1 );
    fclose( f );
  }
  
  if ( stats1[_INCONNU_] > 0 ) {
    sprintf( fname, "%s.X1", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open X1 output file\n", 0 );
    _print_UnknownComponents( f, &thePar1, &thePar2 );
    fclose( f );
  }
  
  if ( stats1[_INCONNU_] > 0 ) {
    sprintf( fname, "%s.X2", par.names.out );
    f = fopen( fname, "w" );
    if ( f == NULL ) VT_ErrorParse("unable to open X2 output file\n", 0 );
    _print_UnknownComponents( f, &thePar2, &thePar1 );
    fclose( f );
  }
  


  /*--- liberations memoires ---*/
  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
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
      else if ( nb == 1 ) {
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else if ( nb == 2 ) {
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
  if (nb == 2)
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






static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;
}
