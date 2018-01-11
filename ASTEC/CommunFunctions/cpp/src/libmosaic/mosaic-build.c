/*************************************************************************
 * minimum.c -
 *
 * $Id: mosaic-build.c,v 1.1 2005/07/20 15:05:03 greg Exp $
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

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>
#include <vt_common.h>
#include <vt_mosaic.h>

typedef enum {
  _PAIR_,
  _OFFSET_
} enumParamFile;

typedef struct local_par {

  char color_mosaic[STRINGLENGTH];
  char output_mosaic[STRINGLENGTH];

  int txmin;
  int txmax;
  int tymin;
  int tymax;

  enumParamFile typeParamFile;

  vt_names names;

  int print_time;

} local_par;








/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );



static char *usage = "[image-in] [image-out]\n\
 [-color-mosaic|-cm %s] [-tx %d %d] [-ty %d %d] [-pair|-offset]\n\
 [-output-mosaic|-om %s]\n\
 [-fusion-mode average|weighted|max]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
 [-v] [-D] [-help]";

static char *detail = "\
if '-offset' is specified, the mosaic input file is a text file:\n\
  'number-of-images'\n\
  '#1-image-filename x-offset y-offset'\n\
  '...'\n\
  '#N-image-filename x-offset y-offset'\n\
if '-pair' is specified, the mosaic input file is a text file:\n\
  'number-of-images'\n\
  '#1-image-filename\n\
  '...'\n\
  '#N-image-filename\n\
  then a list of\n\
  '%d (%d,%d) <-> %d (%d,%d)'\n\
  where '%d (%d,%d)' means '#image (xpos,ypos)'\n\
[-color-mosaic|-cm %s] # color mosaic (one channel per image)\n\
\n\
\
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

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /* read input file
   */
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

  /* read images 
   */
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
  if ( search_neighbors( mosaic, nb_images ) != 1 ) {
    VT_ErrorParse( "unable to compute neighborhood relationships\n", 0 );
  }


  if ( 0 ) {
    fprintf( stdout, "- initial information about mosaic '%s'\n", par.names.in );
    print_mosaic_info( stdout, mosaic, nb_images );
  }


  /*
   */

  build_mosaic( mosaic, nb_images, 
		par.txmin, par.txmax, par.tymin, par.tymax );

  if ( 0 ) {
    fprintf( stdout, "- information about mosaic '%s'\n", par.names.in );
    print_mosaic_info( stdout, mosaic, nb_images );
  }

  if ( par.color_mosaic[0] != '\0' ) {
    image = build_color_mosaic( par.color_mosaic, mosaic, nb_images );
    VT_WriteInrimage( image );
    VT_FreeImage( image );
    VT_Free( (void**)&image );  
  }

  if ( par.output_mosaic[0] != '\0' ) {
    if ( _WriteParamFileWithOffsets( par.output_mosaic, mosaic, nb_images ) != 1 ) {
      VT_ErrorParse( "error when writing output mosaic\n", 0 );
    }
  }

  if ( par.names.out != '\0' ) {
    image = build_large_image( par.names.out, mosaic, nb_images );
    VT_WriteInrimage( image );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
  }



  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( program ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


 return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  char text[STRINGLENGTH];
  int tmp;

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

      else if ( strcmp ( argv[i], "-output-mosaic" ) == 0 ||
		(strcmp ( argv[i], "-om" ) == 0 && argv[i][3] == '\0' ) ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -output-mosaic...\n", 0 );
	strncpy( par->output_mosaic, argv[i], STRINGLENGTH );  
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

     /* parallelism
       */
      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
	setMaxChunks( 100 );
      }
      
      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
	setMaxChunks( 1 );
      }
      
      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-max-chunks", 0 );
	status = sscanf( argv[i], "%d", &tmp );
	if ( status <= 0 ) VT_ErrorParse( "-max-chunks", 0 );
	if ( tmp >= 1 ) setMaxChunks( tmp );
      }
      
      else if ( strcmp ( argv[i], "-parallel-scheduling" ) == 0 || 
		( strcmp ( argv[i], "-ps" ) == 0 && argv[i][3] == '\0') ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-parallel-scheduling", 0 );
	if ( strcmp ( argv[i], "default" ) == 0 ) {
	  setOpenMPScheduling( _DEFAULT_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "static" ) == 0 ) {
	  setOpenMPScheduling( _STATIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_ONE_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "guided" ) == 0 ) {
	  setOpenMPScheduling( _GUIDED_SCHEDULING_ );
	}
	else {
	  fprintf( stderr, "unknown scheduling type: '%s'\n", argv[i] );
	  VT_ErrorParse( "-parallel-scheduling", 0 );
	}
      }

      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	par->print_time = 0;
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
  par->output_mosaic[0] = '\0';

  par->txmin = -40;
  par->txmax =  40;
  par->tymin = -40;
  par->tymax =  40;

  par->typeParamFile = _OFFSET_;

  VT_Names( &(par->names) );

  par->print_time = 1;
}




static double _GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}

static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}
