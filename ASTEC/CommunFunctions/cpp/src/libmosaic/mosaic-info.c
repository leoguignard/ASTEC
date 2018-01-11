/*************************************************************************
 * mosaic-info.c -
 *
 * $Id: mosaic-info.c,v 1.1 2005/07/20 15:05:03 greg Exp $
 *
 * Copyright (c) INRIA 1999-2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Wed Jul 20 11:50:30 MEST 2005
 *
 * ADDITIONS, CHANGES
 *
 */


#include <vt_common.h>
#include <vt_mosaic.h>

typedef enum {
  _PAIR_,
  _OFFSET_
} enumParamFile;

typedef struct local_par {

  char color_mosaic[STRINGLENGTH];

  enumParamFile typeParamFile;

  vt_names names;

} local_par;








/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[mosaic-in] [image-out]\n\
 [-color-mosaic|-cm %s] [-pair|-offset]\n\
 [-fusion-mode average|weighted|max]\n\
 [-v] [-D] [-help]";

static char *detail = "\
if '-offset' is specified, the mosaic input file is a text file:\n\
  'number-of-images'\n\
  '#1-image-filename x-offset y-offset'\n\
  '...'\n\
  '#N-image-filename x-offset y-offset'\n\
[-color-mosaic|-cm %s] # color mosaic (one channel per image)\n\
\n\
 $Revision: 1.1 $ $Date: 2005/07/20 15:05:03 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;

  typeMosaicPart *mosaic;
  int nb_images;
  pair *thePairs = (pair*)NULL;
  int npairs;

  vt_image *image;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  
  if ( par.names.in[0] == '\0' ) {
    VT_ErrorParse( "no input mosaic\n", 0 );
  }


  switch ( par.typeParamFile ) {
  default :
  case _PAIR_ :
    mosaic = _ReadParamFileWithPairs( par.names.in, &nb_images, &thePairs, &npairs );
    break;
  case _OFFSET_ :
    mosaic = _ReadParamFileWithOffsets( par.names.in, &nb_images);
    break;
  }

  if ( mosaic == NULL ) {
    VT_ErrorParse( "error when reading mosaic\n", 0 );
  }

  if ( _ReadMosaicImages ( mosaic, nb_images ) != 1 ) {
    VT_ErrorParse( "error when reading images\n", 0 );
  }

  /* pre-process pairs 
   */
  switch ( par.typeParamFile ) {
  default :
  case _PAIR_ :
    _PreProcessMosaicWithPairs( mosaic, nb_images, thePairs, npairs );
    free( thePairs );
    break;
  case _OFFSET_ :
    break;
  }

  fprintf( stdout, "- full information about mosaic '%s'\n", par.names.in );
  print_mosaic_info( stdout, mosaic, nb_images );
 
  
  if ( par.color_mosaic[0] != '\0' ) {
    image = build_color_mosaic( par.color_mosaic, mosaic, nb_images );
    VT_WriteInrimage( image );
    VT_FreeImage( image );
    VT_Free( (void**)&image );  
  }
  
  if ( par.names.out[0] != '\0' ) {
    image = build_large_image( par.names.out, mosaic, nb_images );
    VT_WriteInrimage( image );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
  }

 return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb;
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




      else if ( strcmp ( argv[i], "-color-mosaic" ) == 0 ||
		(strcmp ( argv[i], "-cm" ) == 0 && argv[i][3] == '\0' ) ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -color-mosaic...\n", 0 );
	strncpy( par->color_mosaic, argv[i], STRINGLENGTH );  
      }


      else if ( strcmp ( argv[i], "-pair" ) == 0 ) {
	par->typeParamFile = _PAIR_;
      }
      else if ( strcmp ( argv[i], "-offset" ) == 0 ) {
	par->typeParamFile = _OFFSET_;
      }


      else if ( strcmp ( argv[i], "-fusion-mode" ) == 0  ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-fusion-mode", 0 );
	if ( strcmp ( argv[i], "average" ) == 0 ) {
	  _set_fusion_mode( AVERAGE );
	}
	else if ( strcmp ( argv[i], "weighted" ) == 0 ) {
	  _set_fusion_mode( WEIGHTING );
	}
	else if ( strcmp ( argv[i], "max" ) == 0 ) {
	  _set_fusion_mode( MAXIMUM );
	}
	else {
	  fprintf( stderr, "unknown fusion mode: '%s'\n", argv[i] );
	  VT_ErrorParse( "-fusion-mode", 0 );
	}
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
  par->color_mosaic[0] = '\0';

  par->typeParamFile = _OFFSET_;

  VT_Names( &(par->names) );
}
