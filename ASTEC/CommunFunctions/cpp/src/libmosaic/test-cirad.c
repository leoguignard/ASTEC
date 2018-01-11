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

typedef struct local_par {
  int txmin;
  int txmax;
  int tymin;
  int tymax;
  vt_names names;
  int type;
} local_par;

typedef struct {
  char imagename[STRINGLENGTH];
  int xoffset_init;
  int yoffset_init;
  int xoffset;
  int yoffset;
  int flag;
  double sum;
  vt_image *image;
} typeMosaicPart;



#define NB 20






typeMosaicPart * _ReadParamFile( char *name, int *nb )
{
  FILE *f;
  typeMosaicPart *mosaic;
  int i =0;
  
  mosaic = malloc( 20*sizeof(typeMosaicPart) );
  
  f = fopen( name, "r" );
  while ( fscanf( f, "%s %d %d", 
		  mosaic[i].imagename, 
		  &(mosaic[i].xoffset_init), &(mosaic[i].yoffset_init) ) == 3 ) {
    i++;
  }
  fclose( f );
  
  *nb = i;

  if ( 1 ) {
    for (i=0; i< (*nb); i++ ) {
      fprintf( stdout, "#%2d: '%40s' +%8d +%8d\n", i, mosaic[i].imagename, 
	       mosaic[i].xoffset_init, mosaic[i].yoffset_init );
    }
  }

  return( mosaic );
}

/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-tx %d %d] [-ty %d %d]\n\
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

  typeMosaicPart *mosaic;
  int nb_images;
  int i;
  
  vt_image *image;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  mosaic = _ReadParamFile( par.names.in, &nb_images);

  for ( i=0; i<nb_images; i++ ) {
    fprintf( stderr, "lecture image %d\n", i );
    mosaic[i].image = _VT_Inrimage( mosaic[i].imagename );
    if ( 1 ) {
      fprintf( stdout, "#%2d: DIMS=[%d,%d,%d]\n", i,
	       mosaic[i].image->dim.x, mosaic[i].image->dim.y,
	       mosaic[i].image->dim.v );
    }
  }
  

  for ( i=0; i<nb_images; i++ ) {
    mosaic[i].xoffset = mosaic[i].xoffset_init;
    mosaic[i].yoffset = mosaic[i].yoffset_init;
  }

  if ( 1 ) {
    image = build_large_mosaic( "init.ppm", mosaic, nb_images );
    VT_WriteInrimage( image );
    
    
    build_mosaic( mosaic, nb_images, par.txmin, par.txmax, par.tymin, par.tymax );
    image = build_large_mosaic( "final.ppm", mosaic, nb_images );
    VT_WriteInrimage( image );
  }
  else {
    image = build_large_mosaic( "final.ppm", mosaic, nb_images );
    VT_WriteInrimage( image );
  }
  image = build_large_image( par.names.out, mosaic, nb_images );
  VT_WriteInrimage( image );

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

      else if ( strcmp ( argv[i], "-tx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tx...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->txmin) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tx...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tx...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->txmax) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tx...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-ty" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ty...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->tymin) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -ty...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ty...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->tymax) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -ty...\n", 0 );
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
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = INT;
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
  par->txmin = -40;
  par->txmax =  40;
  par->tymin = -40;
  par->tymax =  40;
  par->type = TYPE_UNKNOWN;
}
