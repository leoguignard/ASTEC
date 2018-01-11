/*************************************************************************
 * evalBleaching.c -
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 18 jul 2014 12:04:08 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */

#include <histogram.h>
#include <pixel-operation.h>
#include <string-tools.h>

#include <vt_common.h>


static int _time_ = 1;
static int _clock_ = 1;
static int _verbose_ = 1;





typedef enum {
  _MEANS_,
  _MEANS_DIFFS_
} enumComputation;



typedef struct local_par {

  vt_names names;
  bufferType type;

  char *trsfFormat;

  int firstindex;
  int lastindex;

  int fbound;
  int lbound;

  float quantile;
  float dquantile;

  enumComputation typeComputation;

  char *description;

} local_par;

static float minquantile = 0.4;


/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double VT_GetTime();
static double VT_GetClock();


static int _StreamingComputationMeanDiff( stringList *imageFileList, 
					  stringList *maskFileList,
					  stringList *realTrsfFileList,
					  int firstindex, 
					  int lastindex,
					  int fbound,
					  int lbound,
					  float quantile,
					  float dquantile,
					  char *template,
					  char *description );

static int _StreamingComputationMeans( stringList *imageFileList, 
				       stringList *maskFileList,
				       stringList *realTrsfFileList,
				       int firstindex, 
				       int lastindex,
				       int fbound,
				       int lbound,
				       float quantile,
				       float dquantile,
				       char *template,
				       char *description );
 



static char *usage = "[format-in] -f[irst] %d -l[ast] %d\n\
 [-mask format-mask] [-trsf format-tsrf] [file-out]\n\
 [-fbound %d] [-lbound %d] [-quantile %f] [-dquantile %f]\n\
 [-means | -diffs]\n\
 [-desc %s]\n\
 [-time|-notime] [-clock|-noclock]\n\
 [-inv] [-swap] [-v|-nv] [-D] [-help] [encoding-type]";

static char *detail = "\
if 'image-in' is equal to '-', we consider stdin\n\
if 'image-out' is not specified, we consider stdout\n\
if both are not specified, we consider stdin and stdout\n\
###\n\
[-inv]      # inverse 'image-in'\n\
[-swap]     # swap bytes of 'image-in' (if encoded on 2 bytes)\n\
[-v]        # be more verbose\n\
[-D]        # some debug information (if any)\n\
[encoding-type] # for the ouput image\n\
  -o 1    : unsigned char\n\
  -o 2    : unsigned short int\n\
  -o 2 -s : short int\n\
  -o 4 -s : int\n\
  -r      : float\n\
  -type s8|u8|s16|u16|...\n\
  default is same type than 'image-in's\n\
\n";



static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  double time_init, time_exit;
  double clock_init, clock_exit;

  stringList imageFileList, maskFileList, trsfFileList;

  time_init = VT_GetTime();
  clock_init = VT_GetClock();


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  


  /* reading list of images
   */

  initStringList( &imageFileList );
  initStringList( &maskFileList );
  initStringList( &trsfFileList );

  if ( par.names.in[0] != '\0' ) {
    if ( buildStringListFromFormat( par.names.in, par.firstindex, par.lastindex, &imageFileList ) != 1 ) {
      VT_ErrorParse( "unable to build input image list\n", 0);
    }
  }

  if ( _VT_DEBUG_ ) printStringList( stderr, &imageFileList, "Input images" );

  if ( par.names.ext[0] != '\0' ) {
    if ( buildStringListFromFormat( par.names.ext, par.firstindex, par.lastindex, &maskFileList ) != 1 ) {
      freeStringList( &imageFileList );
      VT_ErrorParse( "unable to build input mask list\n", 0);
    }
  }

  if ( _VT_DEBUG_ ) printStringList( stderr, &maskFileList, "Mask  images" );

  if ( par.trsfFormat != (char*)NULL ) {
    if ( buildStringListFromFormat( par.trsfFormat, par.firstindex, par.lastindex, &trsfFileList ) != 1 ) {
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "unable to build input mask list\n", 0);
    }
  }

  if ( _VT_DEBUG_ ) printStringList( stderr, &trsfFileList, "Transformation files" );

  /*
   */
  
  switch ( par.typeComputation ) {
  default :
  case _MEANS_DIFFS_ :
    if ( _StreamingComputationMeanDiff( &imageFileList, 
					&maskFileList,
					&trsfFileList,
					par.firstindex, par.lastindex,
					par.fbound,
					par.lbound,
					par.quantile,
					par.dquantile,
					par.names.out,
					par.description ) != 1 ) {
      freeStringList( &trsfFileList );
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "computation error\n", 0);
    }
    break;
  case _MEANS_ :
    if ( _StreamingComputationMeans( &imageFileList, 
					&maskFileList,
					&trsfFileList,
				     par.firstindex, par.lastindex,
					par.fbound,
					par.lbound,
					par.quantile,
					par.dquantile,
					par.names.out,
					par.description ) != 1 ) {
      freeStringList( &trsfFileList );
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "computation error\n", 0);
    }    
    break;
  }

 


  freeStringList( &trsfFileList );
  freeStringList( &maskFileList );
  freeStringList( &imageFileList );

  
  time_exit = VT_GetTime();
  clock_exit = VT_GetClock();

  if ( _time_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );

  if ( _clock_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );

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
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0' ) {
	_VT_VERBOSE_ = 0;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
	_VT_DEBUG_ = 1;
      }
      else if ( strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0' ) {
	_time_ = 1;
      }
      else if ( strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0' ) {
	_time_ = 0;
      }
      else if ( strcmp ( argv[i], "-clock" ) == 0 && argv[i][6] == '\0' ) {
	_clock_ = 1;
      }
      else if ( strcmp ( argv[i], "-noclock" ) == 0 && argv[i][8] == '\0' ) {
	_clock_ = 0;
      }


      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }



      else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -first ...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->firstindex) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -first ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -last ...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->lastindex) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -last ...", 0 );
      }

      else if ( strcmp ( argv[i], "-mask" ) == 0  && argv[i][5] == '\0' ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -mask...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-trsf" ) == 0  && argv[i][5] == '\0' ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -trsf...\n", 0 );
	par->trsfFormat = argv[i];
      }

      else if ( (strcmp ( argv[i], "-q" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-quantile" ) == 0 && argv[i][9] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -quantile ...\n", 0 );
	status = sscanf( argv[i], "%f", &(par->quantile) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -quantile ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-dq" ) == 0 && argv[i][3] == '\0') 
		|| (strcmp ( argv[i], "-dquantile" ) == 0 && argv[i][10] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -dquantile ...\n", 0 );
	status = sscanf( argv[i], "%f", &(par->dquantile) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -dquantile ...", 0 );
      }

      else if ( strcmp ( argv[i], "-diffs" ) == 0 ) {
	par->typeComputation = _MEANS_DIFFS_;
      }
      else if ( strcmp ( argv[i], "-means" ) == 0 ) {
	par->typeComputation = _MEANS_;
      }


      else if ( strcmp ( argv[i], "-fbound" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -fbound...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->fbound) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -fbound...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-lbound" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -lbound...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->lbound) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -lbound...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-desc" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -desc...\n", 0 );
	par->description = argv[i];
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
      else if ( strcmp ( argv[i], "-type" ) == 0 && argv[i][5] == '\0' ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -type...\n", 0 );
	if ( strcmp ( argv[i], "s8" ) == 0 && argv[i][2] == '\0' ) {
	   par->type = SCHAR;
	}
	else if ( strcmp ( argv[i], "u8" ) == 0 && argv[i][2] == '\0' ) {
	   par->type = UCHAR;
	}
	else if ( strcmp ( argv[i], "s16" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = SSHORT;
	}
	else if ( strcmp ( argv[i], "u16" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = USHORT;
	}
	else if ( strcmp ( argv[i], "s32" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = SINT;
	}
	else if ( strcmp ( argv[i], "u32" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = UINT;
	}
	else if ( strcmp ( argv[i], "s64" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = SLINT;
	}
	else if ( strcmp ( argv[i], "u64" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = ULINT;
	}
	else if ( strcmp ( argv[i], "r32" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = FLOAT;
	}
	else if ( strcmp ( argv[i], "r64" ) == 0 && argv[i][3] == '\0' ) {
	  par->type = DOUBLE;
	}
	else {
	  VT_ErrorParse( "parsing -type: unknown type...\n", 0 );
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
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) ) par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 1) && (r == 0) ) par->type = SSHORT;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 4) && (s == 1) && (r == 0) ) par->type = SINT;
  if ( (o == 4) && (s == 0) && (r == 0) ) par->type = UINT;
  if ( (o == 0 || o == 4) && (s == 0) && (r == 1) ) par->type = FLOAT;
  if ( (o == 8) && (s == 0) && (r == 1) ) par->type = DOUBLE;

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

  par->trsfFormat = (char*)NULL;

  par->firstindex = -1;
  par->lastindex = -1;

  par->fbound = 0;
  par->lbound = 4094;

  par->quantile = 0.95;
  par->dquantile = 0.20;

  par->typeComputation = _MEANS_DIFFS_;

  par->description = (char*)NULL;
}


static double VT_GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double VT_GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}










/************************************************************
 *
 ************************************************************/

/* from reech3d.c
 */
static int _readMatrice( char *name, double *mat )
{
  FILE *fopen(), *fp;
  char text[STRINGLENGTH];
  int i, nbelts = 0;
  int status;
  
  /* lecture de 4 double par ligne
     On prevoit le cas ou la ligne commence par "O8 xxxxx ...
     */

  fp = fopen( name, "r" );
  if ( fp == NULL ) return( 0 );
  
  while ( (nbelts < 16) && (fgets( text, STRINGLENGTH, fp ) != NULL) ) {
    if ( (text[0] == 'O') && (text[1] == '8') ) {
      status = sscanf( &(text[2]), "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    } else {
      status = sscanf( text, "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    }
    if ( _VT_DEBUG_ ) {
      fprintf( stderr, "read %d elements:", status );
      for (i=0; i<status; i++ )
	fprintf( stderr, " %lf", mat[nbelts+i] );
      fprintf( stderr, "\n" );
    }
    if ( status == 4 ) nbelts += 4;
  }
  fclose( fp );

  if ( _VT_DEBUG_ ) {
    fprintf( stderr, " lecture de la matrice %s\n", name );
    fprintf( stderr, " %d elements lus\n", nbelts );
    fprintf( stderr,"   %f %f %f %f\n", mat[0], mat[1], mat[2], mat[3] );
    fprintf( stderr,"   %f %f %f %f\n", mat[4], mat[5], mat[6], mat[7] );
    fprintf( stderr,"   %f %f %f %f\n", mat[8], mat[9], mat[10], mat[11] );
    fprintf( stderr,"   %f %f %f %f\n", mat[12], mat[13], mat[14], mat[15] );
  }
  if ( nbelts == 16 ) return ( 1 );
  return( 0 );
}



/* theTrsf allows to resample 'orig' into 'dest' thus goes
   from 'dest' to 'orig'
*/
static void _changeMatFromRealUnitToVoxelUnit( double *origVoxelSize,
					       double *destVoxelSize,
					       double *mat )
{
  mat[ 0] = mat[ 0] * destVoxelSize[0] / origVoxelSize[0];
  mat[ 1] = mat[ 1] * destVoxelSize[1] / origVoxelSize[0];
  mat[ 2] = mat[ 2] * destVoxelSize[2] / origVoxelSize[0];
  mat[ 3] = mat[ 3]                    / origVoxelSize[0];
  
  mat[ 4] = mat[ 4] * destVoxelSize[0] / origVoxelSize[1];
  mat[ 5] = mat[ 5] * destVoxelSize[1] / origVoxelSize[1];
  mat[ 6] = mat[ 6] * destVoxelSize[2] / origVoxelSize[1];
  mat[ 7] = mat[ 7]                    / origVoxelSize[1];
  
  mat[ 8] = mat[ 8] * destVoxelSize[0] / origVoxelSize[2];
  mat[ 9] = mat[ 9] * destVoxelSize[1] / origVoxelSize[2];
  mat[10] = mat[10] * destVoxelSize[2] / origVoxelSize[2];
  mat[11] = mat[11]                    / origVoxelSize[2];
}










/************************************************************
 *
 ************************************************************/



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



#ifdef _UNUSED_
static int _maxint( int *t, int l )
{
  int m = t[0];
  int i;
  for ( i=1; i<l; i++ )
    if ( m < t[i] ) m = t[i];
  return( m );
}
#endif



static double _maxdouble( double *t, int l )
{
  double m = t[0];
  int i;
  for ( i=1; i<l; i++ )
    if ( m < t[i] ) m = t[i];
  return( m );
}



static double _mindouble( double *t, int l )
{
  double m = t[0];
  int i;
  for ( i=1; i<l; i++ )
    if ( m > t[i] ) m = t[i];
  return( m );
}



#ifdef _UNUSED_
static int _UpdateMask( vt_image *mask, vt_image *im, int fbound, int lbound, int *outofbounds )
{
  char *proc = "_UpdateMask";
  int i;
  int v = im->dim.v * im->dim.x * im->dim.y * im->dim.z;
  int out = 0;

  switch ( mask->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mask type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    switch ( im->type ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such mask type not handled yet\n", proc );
      return( -1 );
    case SSHORT :
      {
	u8 *theMask = (u8*)mask->buf;
	s16 *theBuf = (s16*)im->buf;
	for ( i=0; i<v; i++ ) {
	  if ( theBuf[i] < fbound || lbound < theBuf[i] ) {
	    theMask[i] = 0;
	    out ++;
	  }
	}
      }
      break;
    }
    break; /* mask->type = UCHAR */
  }

  *outofbounds = out;

  return( 1 );
}
#endif










/************************************************************
 *
 ************************************************************/



static void _fprintfBeginFigure( FILE *f )
{
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "//\n" );
  fprintf( f, "// figure\n" );
  fprintf( f, "//\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );

  fprintf( f, "figure;\n" );
  fprintf( f, "myfig=gcf();\n" );
  fprintf( f, "myfig.background = color(\"white\");\n" );
}


static void _fprintfBeginPlot( FILE *f )
{
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  
  fprintf( f, "// myaxes=get(\"current_axes\");\n" );
  fprintf( f, "myaxes=gca();\n" );
  fprintf( f, "set(myaxes,\"auto_clear\",\"off\");\n" );
  fprintf( f, "// removing the trailing ';' allows to see all properties\n" );
  fprintf( f, "\n" );
}



static void _fprintfEndPlot( FILE *f, char *title )
{
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "myaxes.font_size = 4;\n" );
  fprintf( f, "myaxes.font_style = 8;\n" );
  fprintf( f, "\n" );
  if ( title != (char*)NULL ) {
    fprintf( f, "myaxes.title.text = \"%s\";\n", title );
    fprintf( f, "myaxes.title.font_size = 3;\n" );
  } 
  else {
  fprintf( f, "// myaxes.title.text = \"Title\";\n" );
  fprintf( f, "// myaxes.title.font_size = 3;\n" );
  }
  fprintf( f, "// myaxes.x_label.text = \"X Label\";\n" );
  fprintf( f, "// myaxes.y_label.text = \"Y Label\";\n" );
  fprintf( f, "// myaxes.x_label.font_size = 3;\n" );
  fprintf( f, "// myaxes.y_label.font_size = 3;\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xtitle( \"Title\", \"X Label\", \"Y Label\" );\n" );
  fprintf( f, "\n" );
  fprintf( f, "// or \n" );
  fprintf( f, "// xlabel( 'X Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "// ylabel( 'Y Label', 'fontsize', 4, 'fontname', 8 );\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
  fprintf( f, "myentity = gce();\n" );
  fprintf( f, "// myentity.children(1).thickness = 3;\n" );
  fprintf( f, "// myentity.children(1).foreground = 2;\n" );
  fprintf( f, "\n" );
  fprintf( f, "\n" );
}





#ifdef _UNUSED_
static int _fprintfOutOfBounds( FILE *f, int fd, 
				int* outofbounds, int firstindex, int lastindex,
				char *description )
{
  char *proc = "_fprintfOutOfBounds";
  /* writing data
   */

  if ( write( fd, outofbounds, (lastindex-firstindex+1) * sizeof(i32) ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing\n", proc );
    return( -1 );
  }

  /* reading data
   */

  if ( description != (char*)NULL ) {
    fprintf( f, "OUTOFBOUNDS_%s=mget( %d, 'i', f);\n", description, (lastindex-firstindex+1) );
  }
  else {
    fprintf( f, "OUTOFBOUNDS=mget( %d, 'i', f);\n", (lastindex-firstindex+1) );
  }

  /* figure
   */

  _fprintfBeginFigure( f );
  _fprintfBeginPlot( f );
  
  if ( description != (char*)NULL ) {
    fprintf( f, "plot( [%d:%d], OUTOFBOUNDS_%s, \"k-\", \"thickness\", 3 );\n", 
	     firstindex, lastindex, description );
  }
  else {
    fprintf( f, "plot( [%d:%d], OUTOFBOUNDS, \"k-\", \"thickness\", 3 );\n", firstindex, lastindex );
  }

  fprintf( f, "myleg = legend( ['out of bounds']);\n" );
  fprintf( f, "myleg.background = color(\"white\");\n" );
  fprintf( f, "// [%d:%d];\n", firstindex, lastindex );

  fprintf( f, "\n" );
  fprintf( f, "// myaxes.data_bounds = [%d,0;%d,%d];\n", firstindex, lastindex, 
	   _maxint( outofbounds, lastindex ) );
  fprintf( f, "myaxes.data_bounds(1,2) = 0;\n" );

  _fprintfEndPlot( f, description );
  
  if ( description != (char*)NULL ) {
    fprintf( f, "// xs2png(gcf(),'FIG_OUTOFBOUNDS_%s.png');\n", description  );
  }
  else {
    fprintf( f, "// xs2png(gcf(),'FIG_OUTOFBOUNDS.png');\n"  );
  }
  fprintf( f, "\n" );

  return( 1 );
}
#endif





#ifdef _UNUSED_
static int _fprintfMean( FILE *f, int fd, 
			 double *mean, double *rangemean, int firstindex, int lastindex,
			 char *description )
{
  char *proc = "_fprintfMean";
  
  /* writing data
   */

  if ( write( fd, mean, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing\n", proc );
    return( -1 );
  }
  
  if ( write( fd, rangemean, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing\n", proc );
    return( -1 );
  }

  /* reading data
   */

  if ( description != (char*)NULL ) {
    fprintf( f, "MEAN_%s=mget( %d, 'd', f);\n", description, (lastindex-firstindex+1) );
    fprintf( f, "RANGEMEAN_%s=mget( %d, 'd', f);\n", description, (lastindex-firstindex+1) );
  }
  else {
    fprintf( f, "MEAN=mget( %d, 'd', f);\n", (lastindex-firstindex+1) );
    fprintf( f, "RANGEMEAN=mget( %d, 'd', f);\n", (lastindex-firstindex+1) );
  }

  /* figure 
   */

  _fprintfBeginFigure( f );
  _fprintfBeginPlot( f );
  
  if ( description != (char*)NULL ) {
    fprintf( f, "plot( [%d:%d], MEAN_%s, \"k-\", \"thickness\", 3 );\n", 
	     firstindex, lastindex, description );
    fprintf( f, "plot( [%d:%d], RANGEMEAN_%s, \"r-\", \"thickness\", 3 );\n", 
	     firstindex, lastindex, description );
  }
  else {
    fprintf( f, "plot( [%d:%d], MEAN, \"k-\", \"thickness\", 3 );\n", firstindex, lastindex );
    fprintf( f, "plot( [%d:%d], RANGEMEAN, \"r-\", \"thickness\", 3 );\n", firstindex, lastindex );
  }

  fprintf( f, "myleg = legend( ['mean';'range mean']);\n" );
  fprintf( f, "myleg.background = color(\"white\");\n" );
  fprintf( f, "// [%d:%d];\n", firstindex, lastindex );

  fprintf( f, "\n" );
  fprintf( f, "// myaxes.data_bounds = [%d,0;%d,%f];\n", firstindex, lastindex, 
	   _maxdouble( mean, lastindex ) );
  fprintf( f, "myaxes.data_bounds(1,2) = 0;\n" );

  _fprintfEndPlot( f, description );

  if ( description != (char*)NULL ) {
    fprintf( f, "// xs2png(gcf(),'FIG_MEAN_%s.png');\n", description  );
  }
  else {
    fprintf( f, "// xs2png(gcf(),'FIG_MEAN.png');\n"  );
  }
  fprintf( f, "\n" );

  return( 1 );
}
#endif





#ifdef _UNUSED_
static int _fprintfDiff( FILE *f, int fd, 
			 double *diff, double *rangediff, int firstindex, int lastindex,
			 char *description )
{
  char *proc = "_fprintfDiff";
  
  /* writing data
   */

  if ( write( fd, diff, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing\n", proc );
    return( -1 );
  }
  
  if ( write( fd, rangediff, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when writing\n", proc );
    return( -1 );
  }

  /* reading data
   */

  if ( description != (char*)NULL ) {
    fprintf( f, "DIFF_%s=mget( %d, 'd', f);\n", description, (lastindex-firstindex+1) );
    fprintf( f, "RANGEDIFF_%s=mget( %d, 'd', f);\n", description, (lastindex-firstindex+1) );
  }
  else {
    fprintf( f, "DIFF=mget( %d, 'd', f);\n", (lastindex-firstindex+1) );
    fprintf( f, "RANGEDIFF=mget( %d, 'd', f);\n", (lastindex-firstindex+1) );
  }

  /* figure 
   */

  _fprintfBeginFigure( f );
  _fprintfBeginPlot( f );
    
  if ( description != (char*)NULL ) {
    fprintf( f, "plot( [%d:%d], DIFF_%s, \"k-\", \"thickness\", 3 );\n", firstindex, lastindex, description );
    fprintf( f, "plot( [%d:%d], RANGEDIFF_%s, \"r-\", \"thickness\", 3 );\n", firstindex, lastindex, description );
  }
  else {
    fprintf( f, "plot( [%d:%d], DIFF, \"k-\", \"thickness\", 3 );\n", firstindex, lastindex );
    fprintf( f, "plot( [%d:%d], RANGEDIFF, \"r-\", \"thickness\", 3 );\n", firstindex, lastindex );
  }

  fprintf( f, "myleg = legend( ['difference (N+1 - N)';'range difference']);\n" );
  fprintf( f, "myleg.background = color(\"white\");\n" );
  fprintf( f, "// [%d:%d];\n", firstindex, lastindex );

  fprintf( f, "\n" );
  fprintf( f, "// myaxes.data_bounds = [%d,%f;%d,%f];\n",  
	   firstindex, _mindouble( diff, firstindex ), lastindex,
	   _maxdouble( diff, lastindex ) );

  _fprintfEndPlot( f, description );
  
  if ( description != (char*)NULL ) {
    fprintf( f, "// xs2png(gcf(),'FIG_DIFFS_%s.png');\n", description  );
  }
  else {
    fprintf( f, "// xs2png(gcf(),'FIG_DIFFS.png');\n"  );
  }
  fprintf( f, "\n" );

return( 1 );
}
#endif





static void _fprintfData(  FILE *f, int fd, 
			   double *data, char *name, 
			   int firstindex, int lastindex,
			   char *description )
{
  char *proc = "_fprintfData";

  if ( data != (double*)NULL ) {
    if ( write( fd, data, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when writing\n", proc );
      return;
    }

    fprintf( f, "%s", name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, " = mget( %d, 'd', myfile );\n", (lastindex-firstindex+1) );
  }

}





static void _fprintfStddev(  FILE *f, int fd, 
			     double *stddev, char *name, 
			     int firstindex, int lastindex,
			     char *description )
{
  char *proc = "_fprintfData";

  if ( stddev != (double*)NULL ) {
    if ( write( fd, stddev, (lastindex-firstindex+1) * sizeof(r64) ) == -1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error when writing\n", proc );
      return;
    }

    fprintf( f, "SDEV%s", name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, " = mget( %d, 'd', myfile );\n", (lastindex-firstindex+1) );
  }

}





static void _fprintFigureData( FILE *f, int fd, 
			       double *data1, double *stddev1, char *name1, 
			       double *data2, double *stddev2, char *name2, 
			       double *data3, double *stddev3, char *name3, 
			       int firstindex, int lastindex,
			       char *description )
{
  /* writing data
   */
  fprintf( f, "\n\n" );

  _fprintfData( f, fd, data1, name1, firstindex, lastindex, description );
  _fprintfStddev( f, fd, stddev1, name1, firstindex, lastindex, description );

  _fprintfData( f, fd, data2, name2, firstindex, lastindex, description );
  _fprintfStddev( f, fd, stddev2, name2, firstindex, lastindex, description );

  _fprintfData( f, fd, data3, name3, firstindex, lastindex, description );
  _fprintfStddev( f, fd, stddev3, name3, firstindex, lastindex, description );

  fprintf( f, "\n\n" );
}





static void _fprintFigurePlot( FILE *f, int fd, 
			       double *data1, double *stddev1, char *name1, char *leg1,
			       double *data2, double *stddev2, char *name2, char *leg2,
			       double *data3, double *stddev3, char *name3, char *leg3,
			       int firstindex, int lastindex,
			       char *description,
			       char *title )
{
  double min, max, stdddevmax = 0;

  /*
   */
  min = _mindouble( data1, (lastindex-firstindex+1) );
  max = _maxdouble( data1, (lastindex-firstindex+1) );
  if ( data2 != (double*)NULL ) {
    if ( min > _mindouble( data2, (lastindex-firstindex+1) ) ) min = _mindouble( data2, (lastindex-firstindex+1) );
    if ( max < _maxdouble( data2, (lastindex-firstindex+1) ) ) max = _maxdouble( data2, (lastindex-firstindex+1) );
  }
  if ( data3 != (double*)NULL ) {
    if ( min > _mindouble( data3, (lastindex-firstindex+1) ) ) min = _mindouble( data3, (lastindex-firstindex+1) );
    if ( max < _maxdouble( data3, (lastindex-firstindex+1) ) ) max = _maxdouble( data3, (lastindex-firstindex+1) );
  }
  if ( stddev1 != (double*)NULL ) {
    if ( stdddevmax < _maxdouble( stddev1, (lastindex-firstindex+1) ) ) stdddevmax =  _maxdouble( stddev1, (lastindex-firstindex+1) );
  }
  if ( stddev2 != (double*)NULL ) {
    if ( stdddevmax < _maxdouble( stddev2, (lastindex-firstindex+1) ) ) stdddevmax =  _maxdouble( stddev2, (lastindex-firstindex+1) );
  }
  if ( stddev3 != (double*)NULL ) {
    if ( stdddevmax < _maxdouble( stddev3, (lastindex-firstindex+1) ) ) stdddevmax =  _maxdouble( stddev3, (lastindex-firstindex+1) );
  }

  _fprintfBeginPlot( f );

  fprintf( f, "plot( [%d:%d], %s", firstindex, lastindex, name1 );
  if ( description != (char*)NULL ) fprintf( f, "_%s", description );
  fprintf( f, ",\"k-\", \"thickness\", 3 );\n" );

  if ( data2 != (double*)NULL ) {
    fprintf( f, "plot( [%d:%d], %s", firstindex, lastindex, name2 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"r-\", \"thickness\", 3 );\n" );
  }

  if ( data3 != (double*)NULL ) {
    fprintf( f, "plot( [%d:%d], %s", firstindex, lastindex, name3 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"b-\", \"thickness\", 3 );\n" );
  }

  fprintf( f, "\n" );
  
  fprintf( f, "myleg = legend( ['%s'", leg1 );
  if ( data2 != (double*)NULL ) fprintf( f, ";'%s'", leg2 );
  if ( data3 != (double*)NULL ) fprintf( f, ";'%s'", leg3 );
  fprintf( f, "]" );
  if ( data2 == (double*)NULL && data3 == (double*)NULL )
    fprintf( f, ",1" );
  else 
    fprintf( f, ",4" );
  fprintf( f, ");\n" );
  fprintf( f, "myleg.background = color(\"white\");\n" );

  fprintf( f, "\n" );

  if ( stddev1 != (double*)NULL ) {
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name1 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "+SDEV%s", name1 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"k--\", \"thickness\", 2 );\n" );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name1 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "-SDEV%s", name1 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"k--\", \"thickness\", 2 );\n" );
  }

  if ( stddev2 != (double*)NULL ) {
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name2 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "+SDEV%s", name2 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"r--\", \"thickness\", 2 );\n" );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name2 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "-SDEV%s", name2 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"r--\", \"thickness\", 2 );\n" );
  }

  if ( stddev3 != (double*)NULL ) {
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name3 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "+SDEV%s", name3 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"b--\", \"thickness\", 2 );\n" );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name3 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "-SDEV%s", name3 );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"b--\", \"thickness\", 2 );\n" );
  }

  fprintf( f, "\n" );

  fprintf( f, "// myaxes.data_bounds = [%d,%f;%d,%f];\n",  
	   firstindex, min, lastindex, max );

  _fprintfEndPlot( f, title );

  fprintf( f, "\n\n" );

}





static void _fprintFigure( FILE *f, int fd, 
			   double *data1, double *stddev1, char *name1, char *leg1,
			   double *data2, double *stddev2, char *name2, char *leg2,
			   double *data3, double *stddev3, char *name3, char *leg3,
			   int firstindex, int lastindex,
			   char *figname,
			   char *description )
{
  char *proc = "_fprintFigure";

  if ( data1 == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL pointer \n", proc );
    return;
  }

  _fprintFigureData( f, fd, 
		     data1, stddev1, name1,
		     data2, stddev2, name2,
		     data3, stddev3, name3,
		     firstindex, lastindex,
		     description );

  _fprintfBeginFigure( f );

  
  _fprintFigurePlot( f, fd, 
		     data1, stddev1, name1, leg1,
		     data2, stddev2, name2, leg2,
		     data3, stddev3, name3, leg3,
		     firstindex, lastindex,
		     description, description );
  

  fprintf( f, "// xs2png(gcf(), 'FIG" );
  if ( figname != (char*)NULL ) fprintf( f, "_%s", figname );
  if ( description != (char*)NULL ) fprintf( f, "_%s", description );
  fprintf( f, ".png');\n" );

  fprintf( f, "\n\n" );

}





static void _fprintFigureMeansData( FILE *f, int fd, 
				    double *data, double *stddev, 
				    double **qdata, double **qstddev, 
				    int firstindex, int lastindex,
				    int nquantile,
				    char *description )
{
  int j;
  char name[128];

  /* writing data
   */
  fprintf( f, "\n\n" );
  
  _fprintfData( f, fd, data, "MEAN", firstindex, lastindex, description );
  _fprintfStddev( f, fd, stddev, "MEAN", firstindex, lastindex, description );

  for ( j=0; j<nquantile; j++ ) {
    sprintf( name, "QMEAN%02d", j );
    _fprintfData( f, fd, qdata[j], name, firstindex, lastindex, description );
    _fprintfStddev( f, fd, qstddev[j], name, firstindex, lastindex, description );
  }

  fprintf( f, "\n\n" );
}





static void _fprintFigureMeansPlot( FILE *f, int fd, 
				    double *data, double *stddev, 
				    double **qdata, double **qstddev, 
				    int firstindex, int lastindex,
				    float quantile, float dquantile, int nquantile,
				    char *description,
				    char *title )
{
  double min, max, stdddevmax = 0;
  int j;
  char name[128];
  
  /*
   */
  min = _mindouble( data, (lastindex-firstindex+1) );
  max = _maxdouble( data, (lastindex-firstindex+1) );
  for ( j=0; j<nquantile; j++ ) {
    if ( min > _mindouble( qdata[j], (lastindex-firstindex+1) ) ) min = _mindouble( qdata[j], (lastindex-firstindex+1) );
    if ( max < _maxdouble( qdata[j], (lastindex-firstindex+1) ) ) max = _maxdouble( qdata[j], (lastindex-firstindex+1) );
  }

  if ( stddev != (double*)NULL ) {
    if ( stdddevmax < _maxdouble( stddev, (lastindex-firstindex+1) ) ) stdddevmax =  _maxdouble( stddev, (lastindex-firstindex+1) );
  }
  for ( j=0; j<nquantile; j++ ) {
    if ( stdddevmax < _maxdouble( qstddev[j], (lastindex-firstindex+1) ) ) stdddevmax =  _maxdouble( qstddev[j], (lastindex-firstindex+1) );
  }


  /* figure
   */


  _fprintfBeginPlot( f );

  fprintf( f, "plot( [%d:%d], %s", firstindex, lastindex, "MEAN" );
  if ( description != (char*)NULL ) fprintf( f, "_%s", description );
  fprintf( f, ",\"k-\", \"thickness\", 3 );\n" );

  for ( j=0; j<nquantile; j++ ) {
    sprintf( name, "QMEAN%02d", j );
    fprintf( f, "plot( [%d:%d], %s", firstindex, lastindex, name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"" );
    switch ( j ) {
    default :
    case 0  : fprintf( f, "r-" ); break;
    case 1  : fprintf( f, "g-" ); break;
    case 2  : fprintf( f, "b-" ); break;
    case 3  : fprintf( f, "c-" ); break;
    case 4  : fprintf( f, "m-" ); break;
    case 5  : fprintf( f, "k--" ); break;
    case 6  : fprintf( f, "r--" ); break;
    case 7  : fprintf( f, "g--" ); break;
    case 8  : fprintf( f, "b--" ); break;
    case 9  : fprintf( f, "c--" ); break;
    case 10 : fprintf( f, "m--" ); break;
    }
    fprintf( f, "\", \"thickness\", 3 );\n" );
  }
  
  fprintf( f, "\n" );
  
  fprintf( f, "myleg = legend( ['%s'", "mean" );
  for ( j=0; j<nquantile; j++ ) {
    sprintf( name, "mean in [%4.2f%%-%4.2f%%]", quantile-(j+1)*dquantile, quantile-j*dquantile );
    fprintf( f, ";'%s'", name );
  }
  fprintf( f, "], 4);\n" );
  fprintf( f, "myleg.background = color(\"white\");\n" );

  fprintf( f, "\n" );

  if ( stddev != (double*)NULL ) {
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, "MEAN" );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "+SDEV%s", "MEAN" );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"k:\", \"thickness\", 2 );\n" );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, "MEAN" );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "-SDEV%s", "MEAN" );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"k:\", \"thickness\", 2 );\n" );
  }

  for ( j=0; j<nquantile; j++ ) {
    sprintf( name, "QMEAN%02d", j );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "+SDEV%s", name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, ",\"" );
    switch ( j ) {
    default :
    case 0  : fprintf( f, "r:" ); break;
    case 1  : fprintf( f, "g:" ); break;
    case 2  : fprintf( f, "b:" ); break;
    case 3  : fprintf( f, "c:" ); break;
    case 4  : fprintf( f, "m:" ); break;
    case 5  : fprintf( f, "k:" ); break;
    case 6  : fprintf( f, "r:" ); break;
    case 7  : fprintf( f, "g:" ); break;
    case 8  : fprintf( f, "b:" ); break;
    case 9  : fprintf( f, "c:" ); break;
    case 10 : fprintf( f, "m:" ); break;
    }
    fprintf( f, "\", \"thickness\", 2 );\n" );
    fprintf( f, "// plot( [%d:%d], %s", firstindex, lastindex, name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
    fprintf( f, "-SDEV%s", name );
    if ( description != (char*)NULL ) fprintf( f, "_%s", description );
        fprintf( f, ",\"" );
    switch ( j ) {
    default :
    case 0  : fprintf( f, "r:" ); break;
    case 1  : fprintf( f, "g:" ); break;
    case 2  : fprintf( f, "b:" ); break;
    case 3  : fprintf( f, "c:" ); break;
    case 4  : fprintf( f, "m:" ); break;
    case 5  : fprintf( f, "k:" ); break;
    case 6  : fprintf( f, "r:" ); break;
    case 7  : fprintf( f, "g:" ); break;
    case 8  : fprintf( f, "b:" ); break;
    case 9  : fprintf( f, "c:" ); break;
    case 10 : fprintf( f, "m:" ); break;
    }
    fprintf( f, "\", \"thickness\", 2 );\n" );
  }

  fprintf( f, "\n" );

  fprintf( f, "// myaxes.data_bounds = [%d,%f;%d,%f];\n",  
	   firstindex, min, lastindex, max );

  _fprintfEndPlot( f, title );

  fprintf( f, "\n\n" );
}





static void _fprintFigureMeans( FILE *f, int fd, 
				double *data, double *stddev, 
				double **qdata, double **qstddev, 
				int firstindex, int lastindex,
				float quantile, float dquantile, int nquantile,
				char *figname,
				char *description )
{
  char *proc = "_fprintFigureMeans";

  if ( data == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL pointer \n", proc );
    return;
  }

  _fprintFigureMeansData( f, fd,
			  data, stddev,
			  qdata, qstddev,
			  firstindex, lastindex, nquantile,
			  description );

  _fprintfBeginFigure( f );

  _fprintFigureMeansPlot( f, fd,
			  data, stddev,
			  qdata, qstddev,
			  firstindex, lastindex, 
			  quantile, dquantile, nquantile,
			  description, description );

  fprintf( f, "// xs2png(gcf(), 'FIG" );
  if ( figname != (char*)NULL ) fprintf( f, "_%s", figname );
  if ( description != (char*)NULL ) fprintf( f, "_%s", description );
  fprintf( f, ".png');\n" );

  fprintf( f, "\n\n" );

}






/************************************************************
 *
 ************************************************************/




static void _ComputeOutOfBounds( double *nb, 
				 double *percentage,  
				 int fbound,
				 int lbound,
				 typeHistogram *histo1D )
{
  char *proc = "_ComputeOutOfBounds";
  int i;

#define _COMPUTEOUTOFBOUNDS( TYPEINDEX, TYPEHISTO, TYPESUM ) { \
  TYPEINDEX *theIndex = (TYPEINDEX*)histo1D->xaxis.index;      \
  TYPEHISTO *theHisto = (TYPEHISTO*)histo1D->data;             \
  TYPESUM sum, sumout;                                         \
  for ( sum=0, sumout=0, i=0; i<histo1D->xaxis.dim; i++ ) {    \
    sum += (TYPESUM)theHisto[i];                               \
    if ( theIndex[i] < fbound || lbound < theIndex[i] ) {      \
      sumout += (TYPESUM)theHisto[i];                          \
    }                                                          \
  }                                                            \
  *nb = (double)sumout;                                        \
  *percentage = 100.0 * (double)sumout / (double)sum;	       \
}

  switch ( histo1D->xaxis.typeIndex ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such index type not handled yet\n", proc );
    return;
  case SINT :
    switch ( histo1D->typeHisto ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such histogram type not handled yet\n", proc );
      return;
    case UINT :
      _COMPUTEOUTOFBOUNDS( s32, u32, r64 );
      break;
    case SINT :
      _COMPUTEOUTOFBOUNDS( s32, s32, r64 );
      break;
    case FLOAT :
      _COMPUTEOUTOFBOUNDS( s32, r32, r64 );
      break;
    }
    break;
  }
}





static void _ComputeMeans( double *mean,
			   double *ect,
			   double *meanOut,
			   double *ectOut,
			   double *meanQuant,
			   double *ectQuant,
			   int fbound,
			   int lbound,
			   float quantile,
			   typeHistogram *histo1D )
{
  char *proc = "_ComputeMeans";
  int i;
  
#define _COMPUTEMEANS( TYPEINDEX, TYPEHISTO, TYPESUM ) {                                    \
  TYPEINDEX *theIndex = (TYPEINDEX*)histo1D->xaxis.index;                                   \
  TYPEHISTO *theHisto = (TYPEHISTO*)histo1D->data;                                          \
  TYPESUM sum, sumout, sumq, n, nout, nq;                                                   \
  for ( sum=0, sumout=0, n=0, nout=0, i=0; i<histo1D->xaxis.dim; i++ ) {               	    \
    sum += (TYPESUM)theHisto[i] * (TYPESUM)theIndex[i];                                     \
    n   += (TYPESUM)theHisto[i];                                                            \
    if ( theIndex[i] >= fbound && lbound >= theIndex[i] ) {                                 \
      sumout += (TYPESUM)theHisto[i] * (TYPESUM)theIndex[i];                                \
      nout   += (TYPESUM)theHisto[i];                                                       \
    }                                                                                       \
  }                                                                                         \
  for ( sumq=0, nq=0, i=0; i<histo1D->xaxis.dim; i++ ) {                                    \
    if ( nq + (TYPESUM)theHisto[i] < quantile * n ) {                                       \
      sumq += (TYPESUM)theHisto[i] * (TYPESUM)theIndex[i];                                  \
    }                                                                                       \
    else if ( nq < quantile * n && nq + (TYPESUM)theHisto[i] >= quantile * n ) {            \
      sumq += ( quantile * n - nq ) * (TYPESUM)theIndex[i];                                 \
    }                                                                                       \
    nq += (TYPESUM)theHisto[i];                                                             \
  }                                                                                         \
  *mean      = (double)sum / (double)n;                                                     \
  *meanOut   = (double)sumout / (double)nout;                                               \
  *meanQuant = (double)sumq / (double)( quantile * n );                                     \
  for ( sum=0, sumout=0, sumq=0, nq=0, i=0; i<histo1D->xaxis.dim; i++ ) {                   \
    sum += theHisto[i] * (double)(theIndex[i] - (*mean) )*(double)(theIndex[i] - (*mean) ); \
    if ( theIndex[i] >= fbound && lbound >= theIndex[i] ) {                                 \
      sumout += theHisto[i] * (double)(theIndex[i] - (*meanOut) )*(double)(theIndex[i] - (*meanOut) ); \
    }                                                                                       \
    if ( nq + (TYPESUM)theHisto[i] < quantile * n ) {                                       \
      sumq += theHisto[i] * (double)(theIndex[i] - (*meanQuant) )*(double)(theIndex[i] - (*meanQuant) ); \
    }                                                                                       \
    else if ( nq < quantile * n && nq + (TYPESUM)theHisto[i] >= quantile * n ) {            \
      sumq += ( quantile * n - nq ) *                                                       \
	(double)(theIndex[i] - (*meanQuant) )*(double)(theIndex[i] - (*meanQuant) );        \
    }                                                                                       \
    nq    += (TYPESUM)theHisto[i];                                                          \
  }                                                                                         \
  *ect      = sqrt( (double)sum / (double)n );                                              \
  *ectOut   = sqrt( (double)sumout / (double)nout );                                        \
  *ectQuant = sqrt( (double)sumq / (double)( quantile * n ) );                              \
}
  
  switch ( histo1D->xaxis.typeIndex ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such index type not handled yet\n", proc );
    return;
  case SINT :
    switch ( histo1D->typeHisto ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such histogram type not handled yet\n", proc );
      return;
    case UINT :
      _COMPUTEMEANS( s32, u32, r64 );
      break;
    case SINT :
      _COMPUTEMEANS( s32, s32, r64 );
      break;
    case FLOAT :
      _COMPUTEMEANS( s32, r32, r64 );
      break;
    }
    break;
  }
}





static void _ComputeQuantileMeans( double *mean,
				   double *ect,
				   float fquantile,
				   float lquantile,
				   typeHistogram *histo1D )
{
  char *proc = "_ComputeQuantileMeans";
  int i;
  
#define _COMPUTEQUANTILEMEANS( TYPEINDEX, TYPEHISTO, TYPESUM ) {                            \
  TYPEINDEX *theIndex = (TYPEINDEX*)histo1D->xaxis.index;                                   \
  TYPEHISTO *theHisto = (TYPEHISTO*)histo1D->data;                                          \
  TYPESUM sumq, n, nq;                                                                      \
  for ( n=0, i=0; i<histo1D->xaxis.dim; i++ ) {               	                            \
    n   += (TYPESUM)theHisto[i];                                                            \
  }                                                                                         \
  for ( sumq=0, nq=0, i=0; i<histo1D->xaxis.dim; i++ ) {                                    \
    if ( nq + (TYPESUM)theHisto[i] < fquantile * n ) {                                      \
      ;                                                                                     \
    }                                                                                       \
    else if ( nq < fquantile * n && fquantile * n <= nq + (TYPESUM)theHisto[i] ) {          \
      sumq += ( nq + (TYPESUM)theHisto[i] - fquantile * n  ) * (TYPESUM)theIndex[i];        \
    }                                                                                       \
    else if ( nq + (TYPESUM)theHisto[i] < lquantile * n ) {                                 \
      sumq += (TYPESUM)theHisto[i] * (TYPESUM)theIndex[i];                                  \
    }                                                                                       \
    else if ( nq < lquantile * n && lquantile * n <= nq + (TYPESUM)theHisto[i] ) {          \
      sumq += ( lquantile * n - nq ) * (TYPESUM)theIndex[i];                                \
    }                                                                                       \
    nq += (TYPESUM)theHisto[i];                                                             \
  }                                                                                         \
  *mean      = (double)sumq / (double)( (lquantile - fquantile) * n );	                    \
  for ( sumq=0, nq=0, i=0; i<histo1D->xaxis.dim; i++ ) {                                    \
    if ( nq + (TYPESUM)theHisto[i] < fquantile * n ) {                                      \
      ;                                                                                     \
    }                                                                                       \
    else if ( nq < fquantile * n && fquantile * n <= nq + (TYPESUM)theHisto[i] ) {          \
      sumq += ( nq + (TYPESUM)theHisto[i] - fquantile * n  ) * (double)(theIndex[i] - (*mean) )*(double)(theIndex[i] - (*mean) ); \
    }                                                                                       \
    else if ( nq + (TYPESUM)theHisto[i] < lquantile * n ) {                                 \
      sumq += (TYPESUM)theHisto[i] * (double)(theIndex[i] - (*mean) )*(double)(theIndex[i] - (*mean) ); \
    }                                                                                       \
    else if ( nq < lquantile * n && lquantile * n <= nq + (TYPESUM)theHisto[i] ) {          \
      sumq += ( lquantile * n - nq ) * (double)(theIndex[i] - (*mean) )*(double)(theIndex[i] - (*mean) ); \
    }                                                                                       \
    nq += (TYPESUM)theHisto[i];                                                             \
  }                                                                                         \
  *ect      = sqrt( (double)sumq / (double)( (lquantile - fquantile) * n ) );               \
}

  switch ( histo1D->xaxis.typeIndex ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such index type not handled yet\n", proc );
    return;
  case SINT :
    switch ( histo1D->typeHisto ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such histogram type not handled yet\n", proc );
      return;
    case UINT :
      _COMPUTEQUANTILEMEANS( s32, u32, r64 );
      break;
    case SINT :
      _COMPUTEQUANTILEMEANS( s32, s32, r64 );
      break;
    case FLOAT :
      _COMPUTEQUANTILEMEANS( s32, r32, r64 );
      break;
    }
    break;
  }
}





static void _ComputeDiffs( double *mean,
			   double *ect,
			   double *meanOut,
			   double *ectOut,
			   double *meanQuant,
			   double *ectQuant,
			   int fbound,
			   int lbound,
			   float quantile,
			   typeHistogram *histo2D )
{
  char *proc = "_ComputeDiffs";
  int xmin, xmax;
  int ymin, ymax;
  int i, j;

#define _COMPUTEDIFFS( TYPEINDEXX, TYPEINDEXY, TYPEHISTO, TYPESUM ) {            \
  TYPEINDEXX *theX = histo2D->xaxis.index;                                       \
  TYPEINDEXY *theY = histo2D->yaxis.index;                                       \
  TYPEHISTO *theJoint = histo2D->data;                                           \
  TYPEHISTO *projX = (TYPEHISTO*)NULL;                                           \
  TYPEHISTO *projY = (TYPEHISTO*)NULL;                                           \
  TYPESUM sum, sumout, sumq, n, nout, niq, njq, ci, cj, d, dout, dq;             \
  projX = (TYPEHISTO*)malloc( (histo2D->xaxis.dim + histo2D->yaxis.dim) * sizeof(TYPEHISTO) ); \
  if ( projX == (TYPEHISTO*)NULL ) {                                             \
    if ( _verbose_ )                                                             \
      fprintf( stderr, "%s: allocation errot\n", proc );                         \
    return;                                                                      \
  }                                                                              \
  projY = projX;                                                                 \
  projY += histo2D->xaxis.dim;                                                   \
  for ( i=0; i<histo2D->xaxis.dim; i++ ) {                                       \
    for ( projX[i]=0, j=0; j<histo2D->yaxis.dim; j++ )                           \
      projX[i] += theJoint[ j*histo2D->xaxis.dim + i ];                          \
  }                                                                              \
  for ( j=0; j<histo2D->yaxis.dim; j++ ) {                                       \
    for ( projY[j]=0, i=0; i<histo2D->xaxis.dim; i++ )                           \
      projY[j] += theJoint[ j*histo2D->xaxis.dim + i ];                          \
  }                                                                              \
  for ( sum=0, sumout=0, n=0, nout=0, j=0; j<histo2D->yaxis.dim; j++ )           \
    for ( i=0; i<histo2D->xaxis.dim; i++ ) {                                     \
      sum += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * (TYPESUM)( theY[j] - theX[i] ); \
      n   += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ];                      \
      if ( theX[i] >= fbound && lbound >= theX[i]                                \
	   && theY[j] >= fbound && lbound >= theY[j] ) {                         \
	sumout += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * (TYPESUM)( theY[j] - theX[i] ); \
	nout   += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ];                 \
      }                                                                          \
    }                                                                            \
  for ( sumq=0, njq=0, j=0; j<histo2D->yaxis.dim; j++ ) {                        \
    if ( njq + (TYPESUM)projY[j] < quantile * n ) {                              \
      cj = 1.0;                                                                  \
    }                                                                            \
    else if ( njq < quantile * n && njq + (TYPESUM)projY[j] >= quantile * n ) {  \
      cj = ( quantile * n - njq ) / projY[j];                                    \
    }                                                                            \
    else {                                                                       \
      cj = 0.0;                                                                  \
    }                                                                            \
    for ( niq=0, i=0; i<histo2D->xaxis.dim; i++ ) {                              \
      if ( niq + (TYPESUM)projX[i] < quantile * n ) {                            \
	ci = 1;                                                                  \
      }                                                                          \
      else if ( niq < quantile * n && niq + (TYPESUM)projX[i] >= quantile * n ) { \
	ci = ( quantile * n - niq ) / projX[i];                                  \
      }                                                                          \
      else {                                                                     \
	ci = 0.0;                                                                \
      }                                                                          \
      sumq += cj * ci * (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * (TYPESUM)( theY[j] - theX[i] ); \
      niq  += (TYPESUM)projX[i];                                                 \
    }                                                                            \
    njq  += (TYPESUM)projY[j];                                                   \
  }                                                                              \
  *mean      = (double)sum / (double)n;                                          \
  *meanOut   = (double)sumout / (double)nout;                                    \
  *meanQuant = (double)sumq / (double)( quantile * quantile * n );               \
  for ( sum=0, sumout=0, sumq=0, njq=0, j=0; j<histo2D->yaxis.dim; j++ ) {       \
    if ( njq + (TYPESUM)projY[j] < quantile * n ) {                              \
      cj = 1.0;                                                                  \
    }                                                                            \
    else if ( njq < quantile * n && njq + (TYPESUM)projY[j] >= quantile * n ) {   \
      cj = ( quantile * n - njq ) / projY[j];                                    \
    }                                                                            \
    else {                                                                       \
      cj = 0.0;                                                                  \
    }                                                                            \
    for ( niq=0, i=0; i<histo2D->xaxis.dim; i++ ) {                              \
      if ( niq + (TYPESUM)projX[i] < quantile * n ) {                            \
	ci = 1;                                                                  \
      }                                                                          \
      else if ( niq < quantile * n && niq + (TYPESUM)projX[i] >= quantile * n ) { \
	ci = ( quantile * n - niq ) / projX[i];                                  \
      }                                                                          \
      else {                                                                     \
	ci = 0.0;                                                                \
      }                                                                          \
      d    = (TYPESUM)( theY[j] - theX[i] ) - (*mean);                           \
      dout = (TYPESUM)( theY[j] - theX[i] ) - (*meanOut);                        \
      dq   = (TYPESUM)( theY[j] - theX[i] ) - (*meanQuant);                      \
      sum += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * d * d;              \
      if ( theX[i] >= fbound && lbound >= theX[i]                                \
	   && theY[j] >= fbound && lbound >= theY[j] ) {                         \
	sumout += (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * dout * dout;   \
      }                                                                          \
      sumq += ci * cj * (TYPESUM)theJoint[ j*histo2D->xaxis.dim + i ] * dq * dq; \
      niq  += (TYPESUM)projX[i];                                                 \
    }                                                                            \
    njq  += (TYPESUM)projY[j];                                                   \
  }                                                                              \
  *ect      = sqrt( (double)sum / (double)n );                                   \
  *ectOut   = sqrt( (double)sumout / (double)nout );                             \
  *ectQuant = sqrt( (double)sumq / (double)( quantile * quantile * n ) );        \
  free( projX );                                                                 \
}
 
  switch ( histo2D->xaxis.typeIndex ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such x index type not handled yet\n", proc );
   return;
  case SINT :
    xmin = histo2D->xaxis.min.val_s32;
    xmax = histo2D->xaxis.max.val_s32;
    switch ( histo2D->yaxis.typeIndex ) {
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such y ndex type not handled yet\n", proc );
      return;
    case SINT :
      ymin = histo2D->yaxis.min.val_s32;
      ymax = histo2D->yaxis.max.val_s32;
      switch ( histo2D->typeHisto ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such histogram type not handled yet\n", proc );
	return;
      case UINT :
	_COMPUTEDIFFS( s32, s32, u32, r64 )
	break;
      case SINT :
	_COMPUTEDIFFS( s32, s32, s32, r64 )
	break;
      case FLOAT :
	_COMPUTEDIFFS( s32, s32, r32, r64 )
	break;
      }
      /* histo2D->yaxis.typeIndex == SINT */
      break;
    }
    break;
  }
}










/************************************************************
 *
 ************************************************************/





static int _StreamingComputationMeanDiff( stringList *imageFileList, 
					  stringList *maskFileList,
					  stringList *realTrsfFileList,
					  int firstindex, 
					  int lastindex,
					  int fbound,
					  int lbound,
					  float quantile,
					  float dquantile,
					  char *template,
					  char *description )
{
  char *proc = "_StreamingComputationMeanDiff";

  vt_image *nextim = (vt_image*)NULL;
  vt_image *currim = (vt_image*)NULL;
  vt_image *nextmask = (vt_image*)NULL;
  vt_image *currmask = (vt_image*)NULL;

  double nextReadMat[16];
  double currReadMat[16];
  double *nextMat = (double*)NULL;
  double *currMat = (double*)NULL;

  int nquantile = 0;
  float q;

  int nallocated = 14;
  double *allocatedBuffer = (double*)NULL;
  double **allocatedArray = (double**)NULL;

  double *nsaturated;
  double *psaturated;

  double *intensityMean;
  double *intensityStddev;

  double *intensityNsMean;
  double *intensityNsStddev;

  double *intensityQuMean;
  double *intensityQuStddev;

  double *diffMean;
  double *diffStddev;

  double *diffNsMean;
  double *diffNsStddev;

  double *diffQuMean;
  double *diffQuStddev;

  double **intensityQqMean;
  double **intensityQqStddev;

  typeHistogram *histo1D = (typeHistogram *)NULL;
  typeHistogram histo2D;
  
  int i, j;

  int theDim[3];
  double voxelsize[3];

  int nimages;

  char filename[512];
  FILE *f;
  int fd;
  char legend1[128];


  initHistogram( &histo2D );


#define _DEALLOCATIONS_MEANDIFF_ {              \
  if ( currim != (vt_image*)NULL ) {   \
    VT_FreeImage( currim );            \
    VT_Free( (void**)&currim );        \
  }                                    \
  if ( nextim != (vt_image*)NULL ) {   \
    VT_FreeImage( nextim );            \
    VT_Free( (void**)&nextim );        \
  }                                    \
  if ( currmask != (vt_image*)NULL ) { \
    VT_FreeImage( currmask );          \
    VT_Free( (void**)&currmask );      \
  }                                    \
  if ( nextmask != (vt_image*)NULL ) { \
    VT_FreeImage( nextmask );          \
    VT_Free( (void**)&nextmask );      \
  }                                    \
  if ( allocatedBuffer != (double*)NULL ) { \
    free( allocatedBuffer );           \
    allocatedBuffer = (double*)NULL;   \
  }                                    \
  if ( allocatedArray != (double**)NULL ) { \
    free( allocatedArray );            \
    allocatedArray = (double**)NULL;   \
  }                                    \
  if ( histo1D != (typeHistogram*)NULL ) { \
    for ( i=0; i<imageFileList->n; i++ )   \
      freeHistogram( &(histo1D[i]) );      \
     free( histo1D );                      \
     histo1D = (typeHistogram*)NULL;	   \
  }                                        \
  freeHistogram( &histo2D );               \
}

  

  /* allocations
   */
  
  histo1D = (typeHistogram*)malloc( imageFileList->n * sizeof(typeHistogram) );
  if ( histo1D == (typeHistogram*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  for ( i=0; i<imageFileList->n; i++ )
    initHistogram( &(histo1D[i]) );

  
  for ( nquantile=0, q=quantile; q>minquantile; q-=dquantile, nquantile++ )
    ;
  

  allocatedArray = (double**)malloc( (2*nquantile) * sizeof(double*) );
  if ( allocatedArray == (double**)NULL ) {
    free( histo1D );
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation array error\n", proc );
    return( -1 );
  }
  
  intensityQqMean   = allocatedArray;
  intensityQqStddev = allocatedArray;   intensityQqStddev += nquantile;

  allocatedBuffer = (double*)malloc( (nallocated+2*nquantile)*imageFileList->n * sizeof(double) );
  if ( allocatedBuffer == (double*)NULL ) {
    free( allocatedArray );
    free( histo1D );
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  i = 0;

  nsaturated = allocatedBuffer;   nsaturated += i * imageFileList->n;   i++;
  psaturated = allocatedBuffer;   psaturated += i * imageFileList->n;   i++;

  intensityMean = allocatedBuffer;   intensityMean += i * imageFileList->n;   i++;
  intensityStddev = allocatedBuffer;   intensityStddev += i * imageFileList->n;   i++;

  intensityNsMean = allocatedBuffer;   intensityNsMean += i * imageFileList->n;   i++;
  intensityNsStddev = allocatedBuffer;   intensityNsStddev += i * imageFileList->n;   i++;

  intensityQuMean = allocatedBuffer;   intensityQuMean += i * imageFileList->n;   i++;
  intensityQuStddev = allocatedBuffer;   intensityQuStddev += i * imageFileList->n;   i++;

  diffMean = allocatedBuffer;   diffMean += i * imageFileList->n;   i++;
  diffStddev = allocatedBuffer;   diffStddev += i * imageFileList->n;   i++;

  diffNsMean = allocatedBuffer;   diffNsMean += i * imageFileList->n;   i++;
  diffNsStddev = allocatedBuffer;   diffNsStddev += i * imageFileList->n;   i++;

  diffQuMean = allocatedBuffer;   diffQuMean += i * imageFileList->n;   i++;
  diffQuStddev = allocatedBuffer;   diffQuStddev += i * imageFileList->n;   i++;

  for ( j=0; j<2*nquantile; j++ ) {
    allocatedArray[j] = allocatedBuffer;
    allocatedArray[j] += (i+j) * imageFileList->n;
  }



  /*************************************************************
   * first image
   ************************************************************/
  
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... processing image #0 " );
  }
  else if ( _verbose_ >= 1 ) {
    fprintf( stderr, "." );
  }
  
  /* lecture des donnees
   */
  
  currim = _VT_Inrimage( imageFileList->data[0] );
  if ( currim == (vt_image*)NULL ) {
    _DEALLOCATIONS_MEANDIFF_;
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading %s\n", proc, imageFileList->data[0] );
    return( -1 );
  }
  
  theDim[0] = currim->dim.x;
  theDim[1] = currim->dim.y;
  theDim[2] = currim->dim.z;

  voxelsize[0] = currim->siz.x;
  voxelsize[1] = currim->siz.y;
  voxelsize[2] = currim->siz.z;

  if ( maskFileList->n > 0 ) {
    currmask = _VT_Inrimage( maskFileList->data[0] );
    if ( currmask == (vt_image*)NULL ) {
      _DEALLOCATIONS_MEANDIFF_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, maskFileList->data[0] );
      return( -1 );
    }
  }
  
  if ( realTrsfFileList->n > 0 ) {
    if ( _readMatrice( realTrsfFileList->data[0], currReadMat ) != 1 ) {
      _DEALLOCATIONS_MEANDIFF_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, realTrsfFileList->data[0] );
      return( -1 );
    }
    _changeMatFromRealUnitToVoxelUnit( voxelsize, voxelsize, currReadMat );
    currMat = currReadMat;
  }
  else {
    currMat = (double*)NULL;
  }

  

  /* histogram
   */

  if ( maskFileList->n > 0 ) {
    if ( alloc1DHistogramFromImage( &(histo1D[0]), currim->buf, currim->type,
				    currmask->buf, currmask->type,
				    currMat,
				    theDim ) != 1 ) {
      _DEALLOCATIONS_MEANDIFF_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
      return( -1 );
    }
  }
  else {
    if ( alloc1DHistogramFromImage( &(histo1D[0]), currim->buf, currim->type,
				    (void*)NULL, TYPE_UNKNOWN,
				    currMat,
				    theDim ) != 1 ) {
      _DEALLOCATIONS_MEANDIFF_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
      return( -1 );
    }
  }

  

  /* computations
   */

  _ComputeOutOfBounds( &(nsaturated[0]), &(psaturated[0]), 
		       fbound, lbound, &(histo1D[0]) );
  
  _ComputeMeans( &(intensityMean[0]), &(intensityStddev[0]),
		 &(intensityNsMean[0]), &(intensityNsStddev[0]), 
		 &(intensityQuMean[0]), &(intensityQuStddev[0]),
		 fbound, lbound, quantile, &(histo1D[0]) );
  
  for ( j=0; j<nquantile; j++ ) {
    _ComputeQuantileMeans( &(intensityQqMean[j][0]), &(intensityQqStddev[j][0]),
			   quantile-(j+1)*dquantile, quantile-j*dquantile, &(histo1D[0]) );
  }

  



  /*************************************************************
   * loop on images
   ************************************************************/
  
  
  for ( nimages = 1; nimages < imageFileList->n; nimages++ ) {

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "#%d ", nimages );
    }
    else if ( _verbose_ >= 1 ) {
      fprintf( stderr, "." );
    }
    
    /* lecture des donnees
     */
  
    nextim = _VT_Inrimage( imageFileList->data[nimages] );
    if ( nextim == (vt_image*)NULL ) {
      _DEALLOCATIONS_MEANDIFF_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, imageFileList->data[nimages] );
      return( -1 );
    }
    
    if ( maskFileList->n > 0 ) {
      nextmask = _VT_Inrimage( maskFileList->data[nimages] );
      if ( nextmask == (vt_image*)NULL ) {
	_DEALLOCATIONS_MEANDIFF_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when reading %s\n", proc, maskFileList->data[nimages] );
	return( -1 );
      }
    }

    if ( realTrsfFileList->n > 0 ) {
      if ( _readMatrice( realTrsfFileList->data[nimages], nextReadMat ) != 1 ) {
	_DEALLOCATIONS_MEANDIFF_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when reading %s\n", proc, realTrsfFileList->data[nimages] );
	return( -1 );
      }
      _changeMatFromRealUnitToVoxelUnit( voxelsize, voxelsize, nextReadMat );
      nextMat = nextReadMat;
    }
    else {
      nextMat = (double*)NULL;
    }

    

    /* histogram
     */
    
    if ( maskFileList->n > 0 ) {
      if ( alloc1DHistogramFromImage( &(histo1D[nimages]), nextim->buf, nextim->type,
				      nextmask->buf, nextmask->type,
				      nextMat,
				      theDim ) != 1 ) {
	_DEALLOCATIONS_MEANDIFF_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when filling histogram\n", proc );
	return( -1 );
      }
      
      if ( nimages == 1 ) {
	if ( alloc2DHistogramFromImages( &histo2D,
					 currim->buf, currim->type, currmask->buf, currmask->type, currMat,
					 nextim->buf, nextim->type, nextmask->buf, nextmask->type, nextMat, 
					 theDim ) != 1 ) {
	  _DEALLOCATIONS_MEANDIFF_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: error when filling joint histogram\n", proc );
	  return( -1 );
	}
      }
      else {
	if ( fill2DHistogramFromImages( &histo2D,
					currim->buf, currim->type, currmask->buf, currmask->type, currMat,
					nextim->buf, nextim->type, nextmask->buf, nextmask->type, nextMat, 
					theDim ) != 1 ) {
	  _DEALLOCATIONS_MEANDIFF_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: error when filling joint histogram\n", proc );
	  return( -1 );
	}
      }
      
    }
    else {
      if ( alloc1DHistogramFromImage( &(histo1D[nimages]), nextim->buf, nextim->type,
				      (void*)NULL, TYPE_UNKNOWN,
				      nextMat,
				      theDim ) != 1 ) {
	_DEALLOCATIONS_MEANDIFF_;
	if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
	return( -1 );
      }

      if ( nimages == 1 ) {
	if ( alloc2DHistogramFromImages( &histo2D,
					 currim->buf, currim->type, (void*)NULL, TYPE_UNKNOWN, currMat,
					 nextim->buf, nextim->type, (void*)NULL, TYPE_UNKNOWN, nextMat, 
					 theDim ) != 1 ) {
	  _DEALLOCATIONS_MEANDIFF_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: error when filling joint histogram\n", proc );
	  return( -1 );
	}
      }
      else {
	if ( fill2DHistogramFromImages( &histo2D,
					currim->buf, currim->type, (void*)NULL, TYPE_UNKNOWN, currMat,
					nextim->buf, nextim->type, (void*)NULL, TYPE_UNKNOWN, nextMat, 
					theDim ) != 1 ) {
	  _DEALLOCATIONS_MEANDIFF_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: error when filling joint histogram\n", proc );
	  return( -1 );
	}
      }

    }


    /* computations
     */
    
    _ComputeOutOfBounds( &(nsaturated[nimages]), &(psaturated[nimages]), 
			 fbound, lbound, &(histo1D[nimages]) );

    _ComputeMeans( &(intensityMean[nimages]), &(intensityStddev[nimages]),
		   &(intensityNsMean[nimages]), &(intensityNsStddev[nimages]), 
		   &(intensityQuMean[nimages]), &(intensityQuStddev[nimages]),
		   fbound, lbound, quantile, &(histo1D[nimages]) );
    
    for ( j=0; j<nquantile; j++ ) {
      _ComputeQuantileMeans( &(intensityQqMean[j][nimages]), &(intensityQqStddev[j][nimages]),
			     quantile-(j+1)*dquantile, quantile-j*dquantile, &(histo1D[nimages]) );
    }
    
    _ComputeDiffs( &(diffMean[nimages-1]), &(diffStddev[nimages-1]),
		   &(diffNsMean[nimages-1]), &(diffNsStddev[nimages-1]), 
		   &(diffQuMean[nimages-1]), &(diffQuStddev[nimages-1]),
		   fbound, lbound, quantile, &histo2D );


    VT_FreeImage( currim );
    VT_Free( (void**)&currim );
    currim = nextim;
    nextim = (vt_image*)NULL;
    
    VT_FreeImage( currmask );
    VT_Free( (void**)&currmask );
    currmask = nextmask;
    nextmask = (vt_image*)NULL;

    for ( i=0; i<16; i++ ) currReadMat[i] = nextReadMat[i];
    currMat =  ( nextMat == (double*)NULL ) ? (double*)NULL : currReadMat;
    
  }

  if ( _verbose_ )
    fprintf( stderr, "\n" );

  /*************************************************************
   * end of loop
   ************************************************************/



  


  /* computations on images is done here
   */
  
  /* desallocations
   */
  if ( currim != (vt_image*)NULL ) {
    VT_FreeImage( currim );       
    VT_Free( (void**)&currim );   
  }                               
  if ( nextim != (vt_image*)NULL ) {
    VT_FreeImage( nextim );         
    VT_Free( (void**)&nextim );     
  }                                 
  if ( currmask != (vt_image*)NULL ) {
    VT_FreeImage( currmask );         
    VT_Free( (void**)&currmask );     
  }                                   
  if ( nextmask != (vt_image*)NULL ) {
    VT_FreeImage( nextmask );         
    VT_Free( (void**)&nextmask );     
  }                                   

  if ( histo1D != (typeHistogram*)NULL ) { 
    for ( i=0; i<imageFileList->n; i++ )
      freeHistogram( &(histo1D[i]) );
     free( histo1D ); 
     histo1D = (typeHistogram*)NULL; 
  }
  freeHistogram( &histo2D );    



  /* some statistics on differences
   */
  /*

  fprintf( stderr, "      differences: mean = %f +/ %f\n", diffmoy, diffect );
  fprintf( stderr, "range differences: mean = %f +/ %f\n", rdiffmoy, rdiffect );
  */





  /*************************************************************
   * scilab files
   ************************************************************/

  /* open files
   */
  sprintf( filename, "%s.raw", template );
  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    free( allocatedBuffer );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  sprintf( filename, "%s.sce", template );
  f = fopen( filename, "w" );
  if ( f == (FILE*)NULL ) {
    close( fd );
    free( allocatedBuffer );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }
  

  fprintf( f, "\n" );
  fprintf( f, "myfile = mopen('%s.raw','r');\n", _BaseName( template ) );
  fprintf( f, "\n" );

 
  /* outputs 
   */
  if ( 0 ) {
    _fprintFigure( f, fd,
		   nsaturated, (double*)NULL, "NSATURATED", "#4095",
		   (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		   (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		   firstindex, lastindex,
		   "NSATURATED", description );
  }
  
  _fprintFigure( f, fd,
		 psaturated, (double*)NULL, "PSATURATED", "percent. #4095",
		 (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		 (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		 firstindex, lastindex,
		 "PSATURATED", description );

  sprintf( legend1, "mean [0.00-%4.1f]", quantile );
  _fprintFigure( f, fd,
		 intensityMean, intensityStddev, "MEAN", "mean",
		 intensityNsMean, intensityNsStddev, "RMEAN", "mean - 4095",
		 intensityQuMean, intensityQuStddev, "QMEAN", legend1,
		  firstindex, lastindex,
		 "MEAN", description );

  sprintf( legend1, "mean diff [00.0%%-%4.1f]", quantile );
  _fprintFigure( f, fd,
		 diffMean, diffStddev, "DIFF", "mean diff",
		 diffNsMean, diffNsStddev, "RDIFF", "mean diff - 4095",
		 diffQuMean, diffQuStddev, "QDIFF", legend1,
		 firstindex+1, lastindex,
		 "DIFF", description );

  _fprintFigureMeans( f, fd,
		      intensityMean, intensityStddev, 
		      intensityQqMean, intensityQqStddev, 
		      firstindex, lastindex,
		      quantile, dquantile, nquantile,
		      "MEAN2", description );


  fprintf( f, "\n" );
  fprintf( f, "mclose( myfile );\n" );
  fprintf( f, "\n" );



  /*
   */
   _fprintfBeginFigure( f );
   fprintf( f, "myfig.figure_size= [1000,1000];\n" );

   fprintf( f, "\n" );
   fprintf( f, "subplot(2,2,1)\n" );
   fprintf( f, "\n" );

   _fprintFigurePlot( f, fd,
		      psaturated, (double*)NULL, "PSATURATED", "percent. #4095",
		      (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		      (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		      firstindex, lastindex,
		      description, "saturation percentage" );
   
   fprintf( f, "\n" );
   fprintf( f, "subplot(2,2,2)\n" );
   fprintf( f, "\n" );

   _fprintFigureMeansPlot( f, fd,
			   intensityMean, intensityStddev, 
			   intensityQqMean, intensityQqStddev, 
			   firstindex, lastindex,
			   quantile, dquantile, nquantile,
			   description, "intensity averages" );


   fprintf( f, "\n" );
   fprintf( f, "subplot(2,2,3)\n" );
   fprintf( f, "\n" );
   
   sprintf( legend1, "mean [0.00-%4.2f]", quantile );
   _fprintFigurePlot( f, fd,
		      intensityMean, intensityStddev, "MEAN", "mean",
		      intensityNsMean, intensityNsStddev, "RMEAN", "mean - 4095",
		      intensityQuMean, intensityQuStddev, "QMEAN", legend1,
		      firstindex, lastindex,
		      description, "intensity averages" );
  
   fprintf( f, "\n" );
   fprintf( f, "subplot(2,2,4)\n" );
   fprintf( f, "\n" );

  sprintf( legend1, "mean diff [0.00%%-%4.2f]", quantile );
  _fprintFigurePlot( f, fd,
		     diffMean, diffStddev, "DIFF", "mean diff",
		     diffNsMean, diffNsStddev, "RDIFF", "mean diff - 4095",
		     diffQuMean, diffQuStddev, "QDIFF", legend1,
		     firstindex+1, lastindex,
		     description, "intensity difference averages" );
  
  fprintf( f, "xs2png(gcf(), 'FIGCOMPOSITE" );
  if ( description != (char*)NULL ) fprintf( f, "_%s", description );
  fprintf( f, ".png');\n" );

  fprintf( f, "\n\n" );

  fclose( f );
  close( fd );


  free( allocatedBuffer );  

  return( 1 );
}










static int _StreamingComputationMeans( stringList *imageFileList, 
				       stringList *maskFileList,
				       stringList *realTrsfFileList,
				       int firstindex, 
				       int lastindex,
				       int fbound,
				       int lbound,
				       float quantile,
				       float dquantile,
				       char *template,
				       char *description )
{
  char *proc = "_StreamingComputationMeans";

  vt_image *nextim = (vt_image*)NULL;
  vt_image *currim = (vt_image*)NULL;
  vt_image *nextmask = (vt_image*)NULL;
  vt_image *currmask = (vt_image*)NULL;

  double nextReadMat[16];
  double currReadMat[16];
  double *nextMat = (double*)NULL;
  double *currMat = (double*)NULL;

  int nquantile = 0;
  float q;

  int nallocated = 4;
  double *allocatedBuffer = (double*)NULL;
  double **allocatedArray = (double**)NULL;

  double *nsaturated;
  double *psaturated;

  double *intensityMean;
  double *intensityStddev;

  double meanOut, ectOut;
  double meanQuant, ectQuant;

  double **intensityQqMean;
  double **intensityQqStddev;




  typeHistogram *histo1D = (typeHistogram *)NULL;
  
  int i, j;

  int theDim[3];
  double voxelsize[3];

  int nimages;

  char filename[512];
  FILE *f;
  int fd;



#define _DEALLOCATIONS_MEAN_ {              \
  if ( currim != (vt_image*)NULL ) {   \
    VT_FreeImage( currim );            \
    VT_Free( (void**)&currim );        \
  }                                    \
  if ( nextim != (vt_image*)NULL ) {   \
    VT_FreeImage( nextim );            \
    VT_Free( (void**)&nextim );        \
  }                                    \
  if ( currmask != (vt_image*)NULL ) { \
    VT_FreeImage( currmask );          \
    VT_Free( (void**)&currmask );      \
  }                                    \
  if ( nextmask != (vt_image*)NULL ) { \
    VT_FreeImage( nextmask );          \
    VT_Free( (void**)&nextmask );      \
  }                                    \
  if ( allocatedBuffer != (double*)NULL ) { \
    free( allocatedBuffer );           \
    allocatedBuffer = (double*)NULL;   \
  }                                    \
  if ( allocatedArray != (double**)NULL ) { \
    free( allocatedArray );            \
    allocatedArray = (double**)NULL;   \
  }                                    \
  if ( histo1D != (typeHistogram*)NULL ) { \
    for ( i=0; i<imageFileList->n; i++ )   \
      freeHistogram( &(histo1D[i]) );      \
     free( histo1D );                      \
     histo1D = (typeHistogram*)NULL;	   \
  }                                        \
}

  

  /* allocations
   */
  
  histo1D = (typeHistogram*)malloc( imageFileList->n * sizeof(typeHistogram) );
  if ( histo1D == (typeHistogram*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation error\n", proc );
    return( -1 );
  }

  for ( i=0; i<imageFileList->n; i++ )
    initHistogram( &(histo1D[i]) );


  for ( nquantile=0, q=quantile; q>minquantile; q-=dquantile, nquantile++ )
    ;
  

  allocatedArray = (double**)malloc( (2*nquantile) * sizeof(double*) );
  if ( allocatedArray == (double**)NULL ) {
    free( histo1D );
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation array error\n", proc );
    return( -1 );
  }

  intensityQqMean   = allocatedArray;
  intensityQqStddev = allocatedArray;   intensityQqStddev += nquantile;


  allocatedBuffer = (double*)malloc( (nallocated+2*nquantile)*imageFileList->n * sizeof(double) );
  if ( allocatedBuffer == (double*)NULL ) {
    free( allocatedArray );
    free( histo1D );
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation buffer error\n", proc );
    return( -1 );
  }

  i = 0;

  nsaturated = allocatedBuffer;   nsaturated += i * imageFileList->n;   i++;
  psaturated = allocatedBuffer;   psaturated += i * imageFileList->n;   i++;

  intensityMean = allocatedBuffer;   intensityMean += i * imageFileList->n;   i++;
  intensityStddev = allocatedBuffer;   intensityStddev += i * imageFileList->n;   i++;

  for ( j=0; j<2*nquantile; j++ ) {
    allocatedArray[j] = allocatedBuffer;
    allocatedArray[j] += (i+j) * imageFileList->n;
  }




  /*************************************************************
   * first image
   ************************************************************/
  
  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... processing image #0 " );
  }
  else if ( _verbose_ >= 1 ) {
    fprintf( stderr, "." );
  }

  /* lecture des donnees
   */
  
  currim = _VT_Inrimage( imageFileList->data[0] );
  if ( currim == (vt_image*)NULL ) {
    _DEALLOCATIONS_MEAN_;
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading %s\n", proc, imageFileList->data[0] );
    return( -1 );
  }
  
  theDim[0] = currim->dim.x;
  theDim[1] = currim->dim.y;
  theDim[2] = currim->dim.z;

  voxelsize[0] = currim->siz.x;
  voxelsize[1] = currim->siz.y;
  voxelsize[2] = currim->siz.z;

  if ( maskFileList->n > 0 ) {
    currmask = _VT_Inrimage( maskFileList->data[0] );
    if ( currmask == (vt_image*)NULL ) {
      _DEALLOCATIONS_MEAN_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, maskFileList->data[0] );
      return( -1 );
    }
  }
  
  if ( realTrsfFileList->n > 0 ) {
    if ( _readMatrice( realTrsfFileList->data[0], currReadMat ) != 1 ) {
      _DEALLOCATIONS_MEAN_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, realTrsfFileList->data[0] );
      return( -1 );
    }
    _changeMatFromRealUnitToVoxelUnit( voxelsize, voxelsize, currReadMat );
    currMat = currReadMat;
  }
  else {
    currMat = (double*)NULL;
  }

  

  /* histogram
   */

  if ( maskFileList->n > 0 ) {
    if ( alloc1DHistogramFromImage( &(histo1D[0]), currim->buf, currim->type,
				    currmask->buf, currmask->type,
				    currMat,
				    theDim ) != 1 ) {
      _DEALLOCATIONS_MEAN_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
      return( -1 );
    }
  }
  else {
    if ( alloc1DHistogramFromImage( &(histo1D[0]), currim->buf, currim->type,
				    (void*)NULL, TYPE_UNKNOWN,
				    currMat,
				    theDim ) != 1 ) {
      _DEALLOCATIONS_MEAN_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
      return( -1 );
    }
  }

  

  /* computations
   */

  _ComputeOutOfBounds( &(nsaturated[0]), &(psaturated[0]), 
		       fbound, lbound, &(histo1D[0]) );
  
  _ComputeMeans( &(intensityMean[0]), &(intensityStddev[0]),
		 &meanOut, &ectOut,
		 &meanQuant, &ectQuant,
		 fbound, lbound, quantile, &(histo1D[0]) );
  

  for ( j=0; j<nquantile; j++ ) {
    _ComputeQuantileMeans( &(intensityQqMean[j][0]), &(intensityQqStddev[j][0]),
			   quantile-(j+1)*dquantile, quantile-j*dquantile, &(histo1D[0]) );
  }

  



  /*************************************************************
   * loop on images
   ************************************************************/
  
  
  for ( nimages = 1; nimages < imageFileList->n; nimages++ ) {

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, "#%d ", nimages );
    }
    else if ( _verbose_ >= 1 ) {
      fprintf( stderr, "." );
    }
    
    /* lecture des donnees
     */
  
    nextim = _VT_Inrimage( imageFileList->data[nimages] );
    if ( nextim == (vt_image*)NULL ) {
      _DEALLOCATIONS_MEAN_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading %s\n", proc, imageFileList->data[nimages] );
      return( -1 );
    }
    
    if ( maskFileList->n > 0 ) {
      nextmask = _VT_Inrimage( maskFileList->data[nimages] );
      if ( nextmask == (vt_image*)NULL ) {
	_DEALLOCATIONS_MEAN_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when reading %s\n", proc, maskFileList->data[nimages] );
	return( -1 );
      }
    }

    if ( realTrsfFileList->n > 0 ) {
      if ( _readMatrice( realTrsfFileList->data[nimages], nextReadMat ) != 1 ) {
	_DEALLOCATIONS_MEAN_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when reading %s\n", proc, realTrsfFileList->data[nimages] );
	return( -1 );
      }
      _changeMatFromRealUnitToVoxelUnit( voxelsize, voxelsize, nextReadMat );
      nextMat = nextReadMat;
    }
    else {
      nextMat = (double*)NULL;
    }

    

    /* histogram
     */
    
    if ( maskFileList->n > 0 ) {
      if ( alloc1DHistogramFromImage( &(histo1D[nimages]), nextim->buf, nextim->type,
				      nextmask->buf, nextmask->type,
				      nextMat,
				      theDim ) != 1 ) {
	_DEALLOCATIONS_MEAN_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when filling histogram\n", proc );
	return( -1 );
      }
    }
    else {
      if ( alloc1DHistogramFromImage( &(histo1D[nimages]), nextim->buf, nextim->type,
				      (void*)NULL, TYPE_UNKNOWN,
				      nextMat,
				      theDim ) != 1 ) {
	_DEALLOCATIONS_MEAN_;
	if ( _verbose_ )
	fprintf( stderr, "%s: error when filling histogram\n", proc );
	return( -1 );
      }
    }


    /* computations
     */
    
    _ComputeOutOfBounds( &(nsaturated[nimages]), &(psaturated[nimages]), 
			 fbound, lbound, &(histo1D[nimages]) );

    _ComputeMeans( &(intensityMean[nimages]), &(intensityStddev[nimages]),
		   &meanOut, &ectOut,
		   &meanQuant, &ectQuant,
		   fbound, lbound, quantile, &(histo1D[nimages]) );
    
    for ( j=0; j<nquantile; j++ ) {
      _ComputeQuantileMeans( &(intensityQqMean[j][nimages]), &(intensityQqStddev[j][nimages]),
			     quantile-(j+1)*dquantile, quantile-j*dquantile, &(histo1D[nimages]) );
    }

    VT_FreeImage( currim );
    VT_Free( (void**)&currim );
    currim = nextim;
    nextim = (vt_image*)NULL;
    
    VT_FreeImage( currmask );
    VT_Free( (void**)&currmask );
    currmask = nextmask;
    nextmask = (vt_image*)NULL;

    for ( i=0; i<16; i++ ) currReadMat[i] = nextReadMat[i];
    currMat =  ( nextMat == (double*)NULL ) ? (double*)NULL : currReadMat;
    
  }

  if ( _verbose_ )
    fprintf( stderr, "\n" );

  /*************************************************************
   * end of loop
   ************************************************************/



  


  /* computations on images is done here
   */
  
  /* desallocations
   */
  if ( currim != (vt_image*)NULL ) {
    VT_FreeImage( currim );       
    VT_Free( (void**)&currim );   
  }                               
  if ( nextim != (vt_image*)NULL ) {
    VT_FreeImage( nextim );         
    VT_Free( (void**)&nextim );     
  }                                 
  if ( currmask != (vt_image*)NULL ) {
    VT_FreeImage( currmask );         
    VT_Free( (void**)&currmask );     
  }                                   
  if ( nextmask != (vt_image*)NULL ) {
    VT_FreeImage( nextmask );         
    VT_Free( (void**)&nextmask );     
  }                                   

  if ( histo1D != (typeHistogram*)NULL ) { 
    for ( i=0; i<imageFileList->n; i++ )
      freeHistogram( &(histo1D[i]) );
     free( histo1D ); 
     histo1D = (typeHistogram*)NULL; 
  }



  /* some statistics on differences
   */
  /*

  fprintf( stderr, "      differences: mean = %f +/ %f\n", diffmoy, diffect );
  fprintf( stderr, "range differences: mean = %f +/ %f\n", rdiffmoy, rdiffect );
  */





  /*************************************************************
   * scilab files
   ************************************************************/

  /* open files
   */
  sprintf( filename, "%s.raw", template );
  fd = open( filename, O_CREAT | O_TRUNC | O_WRONLY, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
  if ( fd == -1 ) {
    free( allocatedBuffer  );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }

  sprintf( filename, "%s.sce", template );
  f = fopen( filename, "w" );
  if ( f == (FILE*)NULL ) {
    close( fd );
    free( allocatedBuffer );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to open '%s' for writing\n", proc, filename );
    return( -1 );
  }
  

  fprintf( f, "\n" );
  fprintf( f, "myfile = mopen('%s.raw','r');\n", _BaseName( template ) );
  fprintf( f, "\n" );

 
  /* outputs 
   */
  if ( 0 ) {
    _fprintFigure( f, fd,
		   nsaturated, (double*)NULL, "NSATURATED", "#4095",
		   (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		   (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		   firstindex, lastindex,
		   "NSATURATED", description );
  }
  
  _fprintFigure( f, fd,
		 psaturated, (double*)NULL, "PSATURATED", "percent. #4095",
		 (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		 (double*)NULL, (double*)NULL, (char*)NULL, (char*)NULL,
		 firstindex, lastindex,
		 "PSATURATED", description );

  
  _fprintFigureMeans( f, fd,
		       intensityMean, intensityStddev, 
		       intensityQqMean, intensityQqStddev, 
		       firstindex, lastindex,
		       quantile, dquantile, nquantile,
		       "MEAN", description );

  fprintf( f, "\n" );
  fprintf( f, "mclose( myfile );\n" );
  fprintf( f, "\n" );






  fclose( f );
  close( fd );


  free( allocatedArray );  
  free( allocatedBuffer );  

  return( 1 );
}
						 
