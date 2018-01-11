/*************************************************************************
 * meanImages.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */


#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */

#include <convert.h>
#include <pixel-operation.h>
#include <local-operation.h>
#include <string-tools.h>

#include <vt_common.h>

static int _verbose_ = 1;
static int _debug_ = 0;


typedef enum {
  _MAX_,
  _MEAN_,
  _MEDIAN_,
  _MIN_,
  _QUANTILE_,
  _ROBUST_MEAN_,
  _STDDEV_,
  _SUM_,
  _VAR_
} typeOperation;

typedef enum {
  _MEMORY_,
  _STREAMING_
} typeComputation;


typedef struct local_par {

  vt_names names;
  bufferType type;

  char *nameformat;
  char *maskformat;
  int firstindex;
  int lastindex;

  vt_ipt window;
  typeOperation operation;

  double quantile; /* 0: min, 0.5: median, 1:max */
  double lts_fraction; /* samples to be kept */

  typeComputation computation;

  int print_time;

} local_par;








/*-------  ----------*/
static int _InMemoryComputationWithoutMasks( stringList *imageFileList, 
					     local_par *par );
static int _InMemoryComputationWithMasks( stringList *imageFileList, 
					  stringList *maskFileList,
					  local_par *par );
static int _StreamingComputationWithoutMasks( stringList *imageFileList, 
					      local_par *par );
static int _StreamingComputationWithMasks( stringList *imageFileList, 
					   stringList *maskFileList,
					   local_par *par );



/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double VT_GetTime();
static double VT_GetClock();




static char *usage = "[[-image-list|-list|-refl] %s]\n\
 [-image-format|-format %s] [-mask-format %s] -f[irst] %d -l[ast] %d\n\
 [[-res] image-out]\n\
 [-mask-list|-maskl %s]\n\
 [-max|-min|-mean|-median|-quantile|-robust-mean|-sum|-var|-stddev]\n\
 [-window %d %d [%d]]\n\
 [-quantile-value|-q %lf]\n\
 [-lts-cut|-lts-fraction %lf]\n\
 [-streaming | -memory]\n\
 [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
 -image-list|-list|-refl %s # text file = list of images to be processed\n\
 -mask-list|-maskl %s       # text file = list of image masks\n\
 -res %s      # output image\n\
 -format %s   # format 'a la printf' of images to be processed, must contain a '%d'\n\
 -first %d    # first value of the index in the format\n\
 -last %d     # last value of the index in the format\n\
 -max         # \n\
 -min         # \n\
 -mean        # \n\
 -median      # \n\
 -quantile    # \n\
 -robust-mean # \n\
 -sum         # \n\
 -var         # \n\
 -stddev      # \n\
 -window %d %d [%d] # default is 1 1 1\n\
 -quantile-value|-q %lf] # quantile of the retained value\n\
   0:  minimum value, thus '-quantile -q 0'   <=> '-min'\n\
   0.5: median value, thus '-quantile -q 0.5' <=> '-median'\n\
   1:  maximum value, thus '-quantile -q 1'   <=> '-max'\n\
 -lts-cut|-lts-fraction %lf] # fraction of points to be kept for the\n\
   calculation of the robust mean (trimmed estimation)\n\
 -streaming # computation is done by reading one image after the other\n\
   using masks is allowed, but some operations may be not implemented\n\
 -memory    # computation is done by loading all images in memory\n\
   using masks is not implemented\n\
 -inv : inverse 'image-in'\n\
 -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
 -v : mode verbose\n\
 -D : mode debug\n\
 options-de-type : -o 1    : unsigned char\n\
                   -o 2    : unsigned short int\n\
                   -o 2 -s : short int\n\
                   -o 4 -s : int\n\
                   -r      : float\n\
 si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  double time_init, time_exit;
  double clock_init, clock_exit;

  stringList imageFileList, maskFileList;


  time_init = VT_GetTime();
  clock_init = VT_GetClock();


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  

  /* differents cas :
     sequential / memory
     fenetre > 1x1x1 dans image / point dans l'image
     masques / pas masques

     => faire des calculs avec des fenetres > 1x1x1 et des masques 
        demande un peu de travail, donc pas pour maintenant
  */



  /* reading list of images
   */

  initStringList( &imageFileList );
  initStringList( &maskFileList );

  if ( par.nameformat != (char*)NULL ) {
    if ( buildStringListFromFormat( par.nameformat, par.firstindex, par.lastindex, &imageFileList ) != 1 ) {
      VT_ErrorParse( "unable to build input image list\n", 0);
    }
  }
  else if ( buildStringListFromFile( par.names.in, &imageFileList ) != 1 ) {
    VT_ErrorParse( "unable to read input image list\n", 0);
  }

  if ( 0 ) printStringList( stderr, &imageFileList, "Input images" );

  if ( par.maskformat != (char*)NULL ) {
    if ( buildStringListFromFormat( par.maskformat, par.firstindex, par.lastindex, &maskFileList ) != 1 ) {
      freeStringList( &imageFileList );
      VT_ErrorParse( "unable to build mask image list\n", 0);
    }
    if ( maskFileList.n != imageFileList.n ) {
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "image and mask lists have different length\n", 0);
    }
  }
  else if ( par.names.ext[0] != '\0' ) {
    if ( buildStringListFromFile( par.names.ext, &maskFileList ) != 1 ) {
      freeStringList( &imageFileList );
      VT_ErrorParse( "unable to read input mask list\n", 0);
    }
    if ( maskFileList.n != imageFileList.n ) {
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "image and mask lists have different length\n", 0);
    }
  }

  

  if ( par.maskformat != (char*)NULL || par.names.ext[0] != '\0' ) {
    
    switch ( par.computation ) {
    default :
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      VT_ErrorParse( "such computation type not handled yet\n", 0);
    case  _MEMORY_ :
      if ( _InMemoryComputationWithMasks( &imageFileList, &maskFileList, &par ) != 1 ) {
	freeStringList( &maskFileList );
	freeStringList( &imageFileList );
	VT_ErrorParse( "error when computing (memory with masks case)\n", 0);
      }
      break;
    case _STREAMING_ :
      if ( _StreamingComputationWithMasks( &imageFileList, &maskFileList, &par ) != 1 ) {
	freeStringList( &maskFileList );
	freeStringList( &imageFileList );
	VT_ErrorParse( "error when computing (streaming with masks case)\n", 0);
      }
      break;
    }
    
  }
  else {
    
    switch ( par.computation ) {
    default :
        freeStringList( &imageFileList );
      VT_ErrorParse( "such computation type not handled yet\n", 0);
    case  _MEMORY_ :
       if ( _InMemoryComputationWithoutMasks( &imageFileList, &par ) != 1 ) {
	freeStringList( &imageFileList );
	VT_ErrorParse( "error when computing (memory without masks case)\n", 0);
      }
      break;
    case _STREAMING_ :
      if ( _StreamingComputationWithoutMasks( &imageFileList, &par ) != 1 ) {
	freeStringList( &imageFileList );
	VT_ErrorParse( "error when computing (streaming without masks case)\n", 0);
      }
      break;
    }
    
  }
  


  freeStringList( &maskFileList );
  freeStringList( &imageFileList );









  /* end
   */
  
  time_exit = VT_GetTime();
  clock_exit = VT_GetClock();

  if (  par.print_time ) {
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );
  }

  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, status;
  int inputisread = 0;
  int maskisread = 0;
  int outputisread = 0;
  int o=0, s=0, r=0;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	VT_ErrorParse( "'-' is not an argument\n", 0 );
      }

      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) {
	VT_ErrorParse( "\n", 0 );
      }
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	_VT_VERBOSE_ = 1;
	_verbose_ ++;
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
	else _VT_VERBOSE_ ++;
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
	if ( _VT_DEBUG_ <= 0 ) _VT_DEBUG_ = 1;
	else _VT_DEBUG_ ++;
	if ( _debug_ <= 0 ) _debug_ = 1;
	else _debug_ ++;
      }
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }


      /*---  images ---*/
      else if ( strcmp ( argv[i], "-image-list" ) == 0  
		|| (strcmp ( argv[i], "-list" ) == 0 && argv[i][5] == '\0') 
		|| (strcmp ( argv[i], "-refl" ) == 0 && argv[i][5] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -image-list...\n", 0 );
	if ( par->names.in[0] != '\0' ) 
	  VT_ErrorParse( "parsing -image-list: input has already been parsed ...\n", 0 );
	strncpy( par->names.in, argv[i], STRINGLENGTH );  
	inputisread ++;
      }
      else if ( strcmp ( argv[i], "-mask-list" ) == 0  
		|| (strcmp ( argv[i], "-maskl" ) == 0  && argv[i][6] == '\0') ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -mask-list...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
	maskisread ++;
      }
      else if ( strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0' ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -res...\n", 0 );
	if ( par->names.out[0] != '\0' ) 
	  VT_ErrorParse( "parsing -res: output has already been parsed ...\n", 0 );
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	outputisread ++;
      }


      else if ( strcmp ( argv[i], "-format" ) == 0 || strcmp ( argv[i], "-image-format" ) == 0) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -image-format...\n", 0 );
	par->nameformat = argv[i];
	inputisread ++;
      }
      else if ( strcmp ( argv[i], "-mask-format" ) == 0) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -mask-format...\n", 0 );
	par->maskformat = argv[i];
	maskisread ++;
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



      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
	par->operation = _MAX_;
      }
      else if ( strcmp ( argv[i], "-min" ) == 0 ) {
	par->operation = _MIN_;
      }
      else if ( strcmp ( argv[i], "-median" ) == 0 ) {
	par->operation = _MEDIAN_;
      }
      else if ( strcmp ( argv[i], "-mean" ) == 0 ) {
	par->operation = _MEAN_;
      }
      else if ( strcmp ( argv[i], "-quantile" ) == 0 ) {
	par->operation = _QUANTILE_;
      }
      else if ( strcmp ( argv[i], "-robust-mean" ) == 0 
		|| strcmp ( argv[i], "-rmean" ) == 0 ) {
	par->operation = _ROBUST_MEAN_;
      }
      else if ( strcmp ( argv[i], "-sum" ) == 0 ) {
	par->operation = _SUM_;
      }
      else if ( strcmp ( argv[i], "-stddev" ) == 0 ) {
	par->operation = _STDDEV_;
      }
      else if ( strcmp ( argv[i], "-var" ) == 0 ) {
	par->operation = _VAR_;
      }

      else if ( strcmp ( argv[i], "-lts-fraction" ) == 0 
		|| strcmp ( argv[i], "-lts-cut" ) == 0) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "parsing -lts-fraction ...", 0 );
	status = sscanf( argv[i], "%lf", &(par->lts_fraction) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -lts-fraction ...", 0 );
      }
      
      else if ( strcmp ( argv[i], "-quantile-value" ) == 0 
		|| (strcmp ( argv[i], "-q" ) == 0  && argv[i][2] == '\0')  
		|| (strcmp ( argv[i], "-qv" ) == 0  && argv[i][3] == '\0') ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "parsing -quantile-value ...", 0 );
	status = sscanf( argv[i], "%lf", &(par->quantile) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -quantile-value ...", 0 );
      }


      else if ( strcmp (argv[i], "-window" ) == 0 
		|| (strcmp (argv[i], "-w") == 0  && argv[i][2] == '\0') ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "parsing -window %d ...\n", 0  );
	status = sscanf( argv[i], "%d", &(par->window.x) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -window %d", 0 );
	i ++;
	if ( i >= argc)    VT_ErrorParse( "parsing -window %d %d", 0 );
	status = sscanf( argv[i], "%d", &(par->window.y) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -window %d %d", 0 );
	i ++;
	if ( i >= argc) par->window.z = 1;
	else {
	  status = sscanf( argv[i], "%d", &(par->window.z) );
	  if ( status <= 0 ) {
	    i--;
	    par->window.z = 1;
	  }
	}
      }


      else if ( (strcmp ( argv[i], "-memory" ) == 0 && argv[i][7] == '\0') ) {
	par->computation = _MEMORY_;
      }
      else if ( (strcmp ( argv[i], "-streaming" ) == 0 && argv[i][10] == '\0') ) {
	par->computation = _STREAMING_;
      }


      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0')  
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	par->print_time = 0;
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
      if ( inputisread == 0 ) {
	strncpy( par->names.in, argv[i], STRINGLENGTH );  
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	outputisread = 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if ( inputisread == 0 )
    VT_ErrorParse( "no input file/format\n", 0 );
  if ( outputisread ==  0 )
    VT_ErrorParse( "no output file\n", 0 );

  if ( inputisread > 1 ) 
    VT_ErrorParse( "too many input files/formast\n", 0 );
  if ( outputisread > 1 ) 
    VT_ErrorParse( "too many output files\n", 0 );
  if ( maskisread > 1 ) 
    VT_ErrorParse( "too many mask files/formats\n", 0 );


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

  par->nameformat = (char*)NULL;
  par->maskformat = (char*)NULL;
  par->firstindex = 0;
  par->lastindex = 0;

  par->window.x =  par->window.y = par->window.z = 1;
  par->operation = _MEAN_;

  par->quantile = 0.50;
  par->lts_fraction = 0.75;

  par->computation = _MEMORY_;

  par->print_time = 0;
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






























/**************************************************
 *
 * image array
 *
 **************************************************/

static void _FreeImageList( vt_image **array, stringList *list ) 
{
  int i;

  if ( array == (vt_image**)NULL ) return;
  if ( list->n <= 0 ) return;

  for ( i=0; i<list->n; i++ ) {
    if ( array[i] != (vt_image*)NULL ) {
      VT_FreeImage( array[i] );
      VT_Free( (void**)&(array[i]) ); 
    }
  }
}

vt_image ** _ReadImageList( stringList *list ) 
{
  char *proc = "_ReadImageList";
  vt_image **a = (vt_image**)NULL;
  int i;
  
  if ( list->n <= 0 ) return( (vt_image**)NULL );

  a = (vt_image**)malloc( list->n * sizeof( vt_image* ) );
  if ( a == (vt_image**)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate array\n", proc );
    return( (vt_image**)NULL );
  }
  
  for ( i=0; i<list->n; i++ ) a[i] = (vt_image*)NULL;

  for ( i=0; i<list->n; i++ ) {
     a[i] = _VT_Inrimage( list->data[i] );
     if ( a[i] ==  (vt_image*)NULL ) {
       _FreeImageList( a, list );
       free( a );
       if ( _verbose_ ) 
	 fprintf( stderr, "%s: unable to read image #%d '%s'\n", proc, i, list->data[i] );
       return( (vt_image**)NULL );
     }
  }
  
  for ( i=1; i<list->n; i++ ) {
    if ( a[i]->dim.v != a[0]->dim.v
	 || a[i]->dim.x != a[0]->dim.x
	 || a[i]->dim.y != a[0]->dim.y
	 || a[i]->dim.z != a[0]->dim.z ) {
      _FreeImageList( a, list );
      free( a );
       if ( _verbose_ ) 
	 fprintf( stderr, "%s: image #%d '%s' has different dimensions than image #0 '%s'\n", 
		  proc, i, list->data[i], list->data[0]);
       return( (vt_image**)NULL );
    }
    if ( a[i]->type != a[0]->type ) {
      _FreeImageList( a, list );
      free( a );
       if ( _verbose_ ) 
	 fprintf( stderr, "%s: image #%d '%s' has a different type than image #0 '%s'\n", 
		  proc, i, list->data[i], list->data[0]);
       return( (vt_image**)NULL );
    }
  } 

  return( a );
}







/**************************************************
 *
 * common procedures
 *
 **************************************************/





static int _AllocAuxiliaryImages( vt_image *imRes,
				  vt_image *imSum,
				  vt_image *imSumSqr,
				  vt_image *image,
				  char *name,
				  typeOperation operation )
{
  char * proc = "_AllocAuxiliaryImages";
  
  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  case _MAX_ :
  case _MIN_ :
    VT_InitFromImage( imRes, image, name, image->type );
    if ( VT_AllocImage( imRes ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate result image\n", proc);
      return( -1 );
    }
    break;
  
  case _STDDEV_ :
  case _VAR_ :
    VT_InitFromImage( imSumSqr, image, name, FLOAT );
    if ( VT_AllocImage( imSumSqr ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate sum of squares image\n", proc );
      return( -1 );
    }
    
  case _MEAN_ :
  case _SUM_ :
    VT_InitFromImage( imSum, image, name, FLOAT );
    if ( VT_AllocImage( imSum ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate sum image\n", proc );
      return( -1 );
    }
    break;
  }
  
  return( 1 );
}





static int _AllocInitAuxiliaryWithMasks( vt_image *imTotalRes,
					 vt_image *imTotalSum,
					 vt_image *imTotalSumSqr,
					 vt_image *imTotalMask,
					 vt_image *image,
					 vt_image *mask,
					 char *name,
					 typeOperation operation )
{
  char *proc = "_AllocInitAuxiliaryWithMasks";

  unsigned short int *bufTotalMask = NULL;

  int theDim[3] = {0,0,0};
  size_t v = 0;
  size_t i;



  /* allocations
   */
  if ( _AllocAuxiliaryImages( imTotalRes, imTotalSum, imTotalSumSqr,
			      image, name, operation ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary images\n", proc);
    return( -1 );
  }

  VT_InitFromImage( imTotalMask, image, name, USHORT );
  if ( VT_AllocImage( imTotalMask ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate sum of masks image\n", proc );
    return( -1 );
  }
  bufTotalMask = (unsigned short int *)imTotalMask->buf;



  theDim[0] = image->dim.v * image->dim.x;
  theDim[1] = image->dim.y;
  theDim[2] = image->dim.z;
  v = (size_t)theDim[0] * (size_t)theDim[1] * (size_t)theDim[2];



  /* initialisation of auxiliary images with the first image
   */
  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  case _MAX_ :
  case _MIN_ :
    if ( maskImage( image->buf, image->type,
		    mask->buf, mask->type,
		    imTotalRes->buf, imTotalRes->type,
		    theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during min/max init\n", proc );
      return( -1 );
    }
    break;
  case _STDDEV_ :
  case _VAR_ :
    if ( sqrImage( image->buf, image->type,
		   imTotalSumSqr->buf, imTotalSumSqr->type,
		   theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during stddev/var init (1)\n", proc );
      return( -1 );
    }
    if ( maskImage( imTotalSumSqr->buf, imTotalSumSqr->type,
		    mask->buf, mask->type,
		    imTotalSumSqr->buf, imTotalSumSqr->type,
		    theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during stddev/var init (2)\n", proc );
      return( -1 );
    }
  case _SUM_ :
  case _MEAN_ :
    if ( ConvertBuffer( image->buf, image->type,
			imTotalSum->buf, imTotalSum->type, v ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during conversion (2)\n", proc );
      return( -1 );
    }
    if ( maskImage( imTotalSum->buf, imTotalSum->type,
		    mask->buf, mask->type,
		    imTotalSum->buf, imTotalSum->type,
		    theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during sum/mean init (2)\n", proc );
      return( -1 );
    }
    break;
  }

  switch ( mask->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mask image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 *theMask = (u8*)mask->buf;
      for ( i=0; i<v; i++ )
	bufTotalMask[i] = (theMask[i] > 0) ? 1 : 0;
    }
    break;
  case USHORT :
    {
      u16 *theMask = (u16*)mask->buf;
      for ( i=0; i<v; i++ )
	bufTotalMask[i] = (theMask[i] > 0) ? 1 : 0;
    }
    break;
  }

  return( 1 );
}





static int _UpdateAuxiliaryWithMasks( vt_image *imTotalRes,
				   vt_image *imTotalSum,
				   vt_image *imTotalSumSqr,
				   vt_image *imTotalMask,
				   vt_image *image,
				   vt_image *mask,
				   typeOperation operation )
{
  char * proc = "_UpdateAuxiliaryWithMasks";

  float *bufTotalSum = NULL;
  float *bufTotalSumSqr = NULL;
  unsigned short int *bufTotalMask = (unsigned short int *)imTotalMask->buf;

  size_t v = 0;
  size_t i;


  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  case _MAX_ :
  case _MIN_ :
    break;
    
  case _STDDEV_ :
  case _VAR_ :
    bufTotalSumSqr = (float*)imTotalSumSqr->buf;
    
  case _MEAN_ :
  case _SUM_ :
    bufTotalSum = (float*)imTotalSum->buf;
    break;
  }



  v = (size_t)image->dim.v * (size_t)image->dim.x * (size_t)image->dim.y * (size_t)image->dim.z;
  


#define _MAX_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->buf;     \
  TYPE *resBuf = (TYPE*)imTotalRes->buf; \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    if ( bufTotalMask[i] == 0 ) resBuf[i] = theBuf[i]; \
    else if ( resBuf[i] < theBuf[i] ) resBuf[i] = theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}
  
#define _MIN_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->buf;     \
  TYPE *resBuf = (TYPE*)imTotalRes->buf; \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    if ( bufTotalMask[i] == 0 ) resBuf[i] = theBuf[i]; \
    else if ( resBuf[i] > theBuf[i] ) resBuf[i] = theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}

#define _SQR_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->buf;     \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    bufTotalSumSqr[i] += (float)theBuf[i] * (float)theBuf[i]; \
    bufTotalSum[i] += (float)theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}

#define _SUM_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->buf;     \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    bufTotalSum[i] += (float)theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}


#define _OPERATION_WITHMASKS_( TYPEM ) { \
  TYPEM *theMask = (TYPEM*)mask->buf; \
  switch ( operation ) {              \
  default :                           \
    if ( _verbose_ )                  \
      fprintf( stderr, "%s: such operation not handled yet\n", proc ); \
    return( -1 );                     \
  case _MAX_ :                        \
    switch( image->type ) {           \
    default :                         \
      if ( _verbose_ )                \
	fprintf( stderr, "%s: such image type handled yet\n", proc ); \
      return( -1 );                   \
    case UCHAR :                      \
      _MAX_WITHMASKS_( u8 );          \
      break;                          \
    case USHORT :                     \
      _MAX_WITHMASKS_( u16 );         \
      break;                          \
    case SSHORT :                     \
      _MAX_WITHMASKS_( s16 );         \
      break;                          \
    }                                 \
    break;                            \
  case _MIN_ :                        \
    switch( image->type ) {           \
    default :                         \
      if ( _verbose_ )                \
	fprintf( stderr, "%s: such image type handled yet\n", proc ); \
      return( -1 );                   \
    case UCHAR :                      \
      _MIN_WITHMASKS_( u8 );          \
      break;                          \
    case USHORT :                     \
      _MIN_WITHMASKS_( u16 );         \
      break;                          \
    case SSHORT :                     \
      _MIN_WITHMASKS_( s16 );         \
      break;                          \
    }                                 \
    break;                            \
  case _STDDEV_ :                     \
  case _VAR_ :                        \
    switch( image->type ) {           \
    default :                         \
      if ( _verbose_ )                \
	fprintf( stderr, "%s: such image type handled yet\n", proc ); \
      return( -1 );                   \
    case UCHAR :                      \
      _SQR_WITHMASKS_( u8 );          \
      break;                          \
    case USHORT :                     \
      _SQR_WITHMASKS_( u16 );         \
      break;                          \
    case SSHORT :                     \
      _SQR_WITHMASKS_( s16 );         \
      break;                          \
    }                                 \
    break;                            \
  case _SUM_ :                        \
  case _MEAN_ :                       \
    switch( image->type ) {           \
    default :                         \
      if ( _verbose_ )                \
	fprintf( stderr, "%s: such image type handled yet\n", proc ); \
      return( -1 );                   \
    case UCHAR :                      \
      _SUM_WITHMASKS_( u8 );          \
      break;                          \
    case USHORT :                     \
      _SUM_WITHMASKS_( u16 );         \
      break;                          \
    case SSHORT :                     \
      _SUM_WITHMASKS_( s16 );         \
      break;                          \
    }                                 \
    break;                            \
  }                                   \
  break;                              \
}


  switch ( mask->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such mask image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    _OPERATION_WITHMASKS_( u8 );
    break;
  case USHORT :
    _OPERATION_WITHMASKS_( u16 );
    break;
  }

  return( 1 );
}





static int _ResultFromAuxiliaryWithMasks( vt_image *imTotalRes,
				       vt_image *imTotalSum,
				       vt_image *imTotalSumSqr,
				       vt_image *imTotalMask,
				       char *name,
				       typeOperation operation,
				       bufferType type )
{
  char * proc = "_ResultFromAuxiliaryWithMasks";

  vt_image *imResult = (vt_image*)NULL;

  float *bufTotalSum = NULL;
  float *bufTotalSumSqr = NULL;
  unsigned short int *bufTotalMask = (unsigned short int *)imTotalMask->buf;

  size_t v = 0;
  size_t i;


  switch ( operation ) {
  default :
    VT_FreeImage ( imTotalSumSqr );
    VT_FreeImage ( imTotalMask );
    VT_FreeImage ( imTotalSum );         
    VT_FreeImage ( imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  case _MAX_ :
  case _MIN_ :
    imResult = imTotalRes;
    break;

  case _SUM_ :
    imResult = imTotalSum;
    break;

  case _VAR_ :
    imResult = imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum->buf;
    bufTotalSumSqr = (float*)imTotalSumSqr->buf;
    for ( i=0; i<v; i++ ) {
      if ( bufTotalMask[i] == 0 ) continue;
      bufTotalSum[i] /= (float)bufTotalMask[i];
      bufTotalSumSqr[i] /= (float)bufTotalMask[i];
      bufTotalSum[i] = bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i];
    }
    break;

  case _STDDEV_ :
    imResult = imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum->buf;
    bufTotalSumSqr = (float*)imTotalSumSqr->buf;
    for ( i=0; i<v; i++ ) {
      if ( bufTotalMask[i] == 0 ) continue;
      bufTotalSum[i] /= (float)bufTotalMask[i];
      bufTotalSumSqr[i] /= (float)bufTotalMask[i];
      bufTotalSum[i] = sqrt( bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i] );
    }
    break;
    
  case _MEAN_ :
    imResult = imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum->buf;
    for ( i=0; i<v; i++ ) {
      if ( bufTotalMask[i] == 0 ) continue;
      bufTotalSum[i] /= (float)bufTotalMask[i];
    }
    break;
  }



  VT_FreeImage ( imTotalSumSqr );
  VT_FreeImage ( imTotalMask );



  if ( type == TYPE_UNKNOWN || type == imResult->type ) {
    if ( VT_WriteInrimage( imResult ) == -1 ) {
      VT_FreeImage ( imTotalSum );         
      VT_FreeImage ( imTotalRes );	
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
  }
  else {
     VT_InitFromImage( imTotalMask, imResult, name, type );
    if ( VT_AllocImage( imTotalMask ) != 1 ) {
      VT_FreeImage ( imTotalSum );         
      VT_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate auxiliary result image\n", proc );
      return( -1 );
    }
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    if ( ConvertBuffer( imResult->buf, imResult->type,
			imTotalMask->buf, imTotalMask->type,
			v ) != 1 ) {
      VT_FreeImage ( imTotalMask );
      VT_FreeImage ( imTotalSum );         
      VT_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to convert image\n", proc );
      return( -1 );
    }
    if ( VT_WriteInrimage( imTotalMask ) == -1 ) {
      VT_FreeImage ( imTotalMask );
      VT_FreeImage ( imTotalSum );         
      VT_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
    VT_FreeImage ( imTotalMask );
  }

  VT_FreeImage ( imTotalSum );         
  VT_FreeImage ( imTotalRes );	

  return( 1 );
}





static int _addSqrImage( void *bufferIn1, bufferType typeIn1,
			void *bufferIn2, bufferType typeIn2,
			void *bufferOut, bufferType typeOut,
			int *bufferDims )
{
  char *proc = "_addSqrImage";
  size_t i, v;
  
  v = (size_t)bufferDims[0] * (size_t)bufferDims[1] * (size_t)bufferDims[2];

#define _ADDSQR_OPERATION( TYPE1, TYPE2 ) { \
  TYPE1 *theBuf1 = (TYPE1*)bufferIn1;       \
  TYPE2 *theBuf2 = (TYPE2*)bufferIn2;       \
  switch ( typeOut ) {                      \
  default :                                 \
    if ( _verbose_ )			    \
      fprintf( stderr, "%s: such output image type not handled yet\n", proc ); \
    return( -1 );                           \
  case FLOAT :                              \
    {                                       \
       r32 *resBuf = (r32*)bufferOut;       \
       for ( i=0; i<v; i++ )                \
          resBuf[i] = (r32)theBuf1[i] + (r32)theBuf2[i]; \
    }                                       \
    break;                                  \
  case DOUBLE :                             \
    {                                       \
       r64 *resBuf = (r64*)bufferOut;       \
       for ( i=0; i<v; i++ )                \
          resBuf[i] = (r64)theBuf1[i] + (r64)theBuf2[i]; \
    }                                       \
    break;                                  \
  }                                         \
}

  switch ( typeIn1 ) {
  default :
    if ( _verbose_ )		
      fprintf( stderr, "%s: such first input image type not handled yet\n", proc ); 
    return( -1 );   
  case UCHAR :
    switch ( typeIn2 ) {
    default :
      if ( _verbose_ )		
	fprintf( stderr, "%s: such second input image type not handled yet\n", proc ); 
      return( -1 ); 
    case FLOAT :
      _ADDSQR_OPERATION( u8, r32 );
      break;
    case DOUBLE :
      _ADDSQR_OPERATION( u8, r64 );
      break;
    }
    break;
  case USHORT :
    switch ( typeIn2 ) {
    default :
      if ( _verbose_ )		
	fprintf( stderr, "%s: such second input image type not handled yet\n", proc ); 
      return( -1 ); 
    case FLOAT :
      _ADDSQR_OPERATION( u16, r32 );
      break;
    case DOUBLE :
      _ADDSQR_OPERATION( u16, r64 );
      break;
    }
    break;
  case SSHORT :
    switch ( typeIn2 ) {
    default :
      if ( _verbose_ )		
	fprintf( stderr, "%s: such second input image type not handled yet\n", proc ); 
      return( -1 ); 
    case FLOAT :
      _ADDSQR_OPERATION( s16, r32 );
      break;
    case DOUBLE :
      _ADDSQR_OPERATION( s16, r64 );
      break;
    }
    break;
  case FLOAT :
    switch ( typeIn2 ) {
    default :
      if ( _verbose_ )		
	fprintf( stderr, "%s: such second input image type not handled yet\n", proc ); 
      return( -1 ); 
    case UCHAR :
      _ADDSQR_OPERATION( r32, u8 );
      break;
    case USHORT :
      _ADDSQR_OPERATION( r32, u16 );
      break;
    case SSHORT :
      _ADDSQR_OPERATION( r32, s16 );
      break;
    case FLOAT :
      _ADDSQR_OPERATION( r32, r32 );
      break;
    case DOUBLE :
      _ADDSQR_OPERATION( r32, r64 );
      break;
    }
    break;
  case DOUBLE :
    switch ( typeIn2 ) {
    default :
      if ( _verbose_ )		
	fprintf( stderr, "%s: such second input image type not handled yet\n", proc ); 
      return( -1 ); 
    case UCHAR :
      _ADDSQR_OPERATION( r64, u8 );
      break;
    case USHORT :
      _ADDSQR_OPERATION( r64, u16 );
      break;
    case SSHORT :
      _ADDSQR_OPERATION( r64, s16 );
      break;
    case FLOAT :
      _ADDSQR_OPERATION( r64, r32 );
      break;
    case DOUBLE :
      _ADDSQR_OPERATION( r64, r64 );
      break;
    }
    break;
  }

  return( 1 );
}






/**************************************************
 *
 * in memory computation
 *
 **************************************************/





static int _InMemoryComputationWithoutMasks( stringList *imageFileList, 
					     local_par *par )
{
  char *proc = "_InMemoryComputationWithoutMasks";

  vt_image **imageStructure = (vt_image**)NULL;

  int theDim[3] = {0,0,0};
  int winDim[3] = {0,0,0};

  bufferType resType;
  vt_image imTotalRes;
  int i, r;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  imageStructure = _ReadImageList( imageFileList );
  if ( imageStructure == (vt_image**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading images\n", proc );
    return( -1 );
  }

  if ( par->type != TYPE_UNKNOWN ) {
    resType = par->type;
  }
  else {
    switch ( par->operation ) {
    default :
      _FreeImageList( imageStructure, imageFileList );
      VT_Free( (void**)&imageStructure );
      if ( _verbose_ )
	fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
    case _MAX_ :
    case _MEDIAN_ :
    case _MEAN_ :
    case _MIN_ :
    case _QUANTILE_ :
    case _ROBUST_MEAN_ :
      resType = imageStructure[0]->type;
      break;
    case _STDDEV_ :
    case _SUM_ :
    case _VAR_ :
      resType = FLOAT;
      break;
    }
  }

  VT_InitFromImage( &imTotalRes, imageStructure[0], par->names.out, resType );
  if ( VT_AllocImage( &imTotalRes ) != 1 ) {
    _FreeImageList( imageStructure, imageFileList );
    VT_Free( (void**)&imageStructure );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result image\n", proc);
    return( -1 );
  }

  theDim[0] = imageStructure[0]->dim.v * imageStructure[0]->dim.x;
  theDim[1] = imageStructure[0]->dim.y;
  theDim[2] = imageStructure[0]->dim.z;
  winDim[0] = imageStructure[0]->dim.v * par->window.x;
  winDim[1] = par->window.y;
  winDim[2] = par->window.z;


#define _OPERATION_MEMORY_WITHOUTMASKS_( TYPE ) {                       \
  TYPE **buffers = (TYPE**)NULL;                                        \
  buffers = (TYPE**)malloc( imageFileList->n * sizeof( TYPE* ) );       \
  if ( buffers == (TYPE**)NULL ) {                                      \
    VT_FreeImage ( &imTotalRes );	                                \
    _FreeImageList( imageStructure, imageFileList );                    \
    VT_Free( (void**)&imageStructure );                                 \
    if ( _verbose_ )                                                    \
      fprintf( stderr, "%s: unable to allocate buffer array\n", proc);  \
    return( -1 );                                                       \
  }                                                                     \
  for ( i=0; i<imageFileList->n; i++ )                                  \
    buffers[i] = (TYPE*)(imageStructure[i]->buf);                       \
  switch ( par->operation ) {                                           \
  default : r = 0; break;                                               \
  case _MAX_ :                                                          \
    r = maxFilteringBuffers( (void**)buffers, imageStructure[0]->type,  \
			     imageFileList->n, imTotalRes.buf, imTotalRes.type, \
			     theDim, winDim );                          \
    break;                                                              \
  case _MEAN_ :                                                         \
    r = meanFilteringBuffers( (void**)buffers, imageStructure[0]->type, \
				imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				theDim, winDim );                       \
    break;                                                              \
  case _MEDIAN_ :                                                       \
    r = medianFilteringBuffers( (void**)buffers, imageStructure[0]->type, \
				imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				theDim, winDim );                       \
    break;                                                              \
  case _MIN_ :                                                          \
    r = minFilteringBuffers( (void**)buffers, imageStructure[0]->type,  \
				imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				theDim, winDim );                       \
    break;                                                              \
  case _QUANTILE_ :                                                     \
    r = quantileFilteringBuffers( (void**)buffers, imageStructure[0]->type, \
				  imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				  theDim, winDim, par->quantile );      \
    break;                                                              \
  case _ROBUST_MEAN_ :                                                  \
    r = robustMeanFilteringBuffers( (void**)buffers, imageStructure[0]->type, \
				    imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				    theDim, winDim, par->lts_fraction ); \
    break;                                                              \
  case _STDDEV_ :                                                       \
    r = stddevFilteringBuffers( (void**)buffers, imageStructure[0]->type, \
				imageFileList->n, imTotalRes.buf, imTotalRes.type, \
				theDim, winDim );                       \
    break;                                                              \
  case _SUM_ :                                                          \
    r = sumFilteringBuffers( (void**)buffers, imageStructure[0]->type,  \
			     imageFileList->n, imTotalRes.buf, imTotalRes.type, \
			     theDim, winDim );                          \
    break;                                                              \
  case _VAR_ :                                                          \
    r = varFilteringBuffers( (void**)buffers, imageStructure[0]->type,  \
			     imageFileList->n, imTotalRes.buf, imTotalRes.type, \
			     theDim, winDim );                          \
    break;                                                              \
  }                                                                     \
  free( buffers );                                                      \
  if ( r != 1 ) {                                                       \
    VT_FreeImage ( &imTotalRes );	                                \
    _FreeImageList( imageStructure, imageFileList );                    \
    VT_Free( (void**)&imageStructure );                                 \
    if ( _verbose_ )                                                    \
      fprintf( stderr, "%s: error when filtering buffer array\n", proc); \
    return( -1 );                                                       \
  }	                                                                \
} 

  switch( imageStructure[0]->type ) {
  default :
    VT_FreeImage ( &imTotalRes );	
    _FreeImageList( imageStructure, imageFileList );
    VT_Free( (void**)&imageStructure );
    if ( _verbose_ )
      fprintf( stderr, "%s: such input image type not handled yet\n", proc);
    return( -1 );
  case UCHAR :
    _OPERATION_MEMORY_WITHOUTMASKS_( u8 );
    break;
  case SCHAR :
    _OPERATION_MEMORY_WITHOUTMASKS_( s8 );
    break;
  case USHORT :
    _OPERATION_MEMORY_WITHOUTMASKS_( u16 );
    break;
  case SSHORT :
    _OPERATION_MEMORY_WITHOUTMASKS_( s16 );
    break;
  case FLOAT :
    _OPERATION_MEMORY_WITHOUTMASKS_( r32 );
    break;
  }



  /* processing is over
   */
  _FreeImageList( imageStructure, imageFileList );
  VT_Free( (void**)&imageStructure );    

  if ( VT_WriteInrimage( &imTotalRes ) == -1 ) {
    VT_FreeImage ( &imTotalRes );	
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write output image\n", proc );
    return( -1 );
  }

  VT_FreeImage ( &imTotalRes );	

  return( 1 );
}





static int _InMemoryComputationWithMasks( stringList *imageFileList, 
					  stringList *maskFileList,
					  local_par *par )
{
  char *proc = "_InMemoryComputationWithMasks";

  vt_image **imageStructure = (vt_image**)NULL;
  vt_image **maskStructure = (vt_image**)NULL;
  
  vt_image imTotalSum;
  vt_image imTotalSumSqr;
  vt_image imTotalRes;
  vt_image imTotalMask;

  int nimages;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  VT_Image ( &imTotalSum );
  VT_Image ( &imTotalSumSqr );
  VT_Image ( &imTotalRes );
  VT_Image ( &imTotalMask );
  
  imageStructure = _ReadImageList( imageFileList );
  if ( imageStructure == (vt_image**)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading images\n", proc );
    return( -1 );
  }

  if ( maskFileList->n > 0 ) {
    maskStructure = _ReadImageList( maskFileList );
    if ( maskStructure == (vt_image**)NULL ) {
      _FreeImageList( imageStructure, imageFileList );
      VT_Free( (void**)&imageStructure );      
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading masks\n", proc );
      return( -1 );
    }
  }
  else {
    _FreeImageList( imageStructure, imageFileList );
    VT_Free( (void**)&imageStructure );
    return ( _InMemoryComputationWithoutMasks( imageFileList, par ) );
  }

  if ( maskStructure[0]->dim.v != imageStructure[0]->dim.v 
       || maskStructure[0]->dim.x != imageStructure[0]->dim.x 
       || maskStructure[0]->dim.y != imageStructure[0]->dim.y 
       || maskStructure[0]->dim.z != imageStructure[0]->dim.z ) { 
    if ( maskFileList->n > 0 ) {          
      _FreeImageList( maskStructure, maskFileList );   
      VT_Free( (void**)&maskStructure );  
    }                                     
    _FreeImageList( imageStructure, imageFileList );     
    VT_Free( (void**)&imageStructure );	
    if ( _verbose_ )
      fprintf( stderr, "%s: images and masks have different dimensions\n", proc );
    return( -1 );
  }




  if ( par->window.x != 1 || par->window.y != 1 || par->window.z != 1 ) {
    if ( maskFileList->n > 0 ) {          
      _FreeImageList( maskStructure, maskFileList ); 
      VT_Free( (void**)&maskStructure );  
    }                                     
    _FreeImageList( imageStructure, imageFileList ); 
    VT_Free( (void**)&imageStructure );	
    if ( _verbose_ )
      fprintf( stderr, "%s: computation with window > 1x1x1 not handled yet\n", proc );
    return( -1 );
  }



  /* we can go
   */

#define _DEALLOCATIONS_MEMORY_WITHMASKS_ { \
  VT_FreeImage ( &imTotalSum );         \
  VT_FreeImage ( &imTotalSumSqr );      \
  VT_FreeImage ( &imTotalRes );		\
  VT_FreeImage ( &imTotalMask );        \
  if ( maskFileList->n > 0 ) {          \
    _FreeImageList( maskStructure, maskFileList ); \
    VT_Free( (void**)&maskStructure );  \
  }                                     \
  _FreeImageList( imageStructure, imageFileList ); \
  VT_Free( (void**)&imageStructure );	\
}


  
  /* allocation and initialisation of auxiliary images
   */
  if ( _AllocInitAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				  &imTotalSumSqr, &imTotalMask,
				  imageStructure[0], maskStructure[0],
				  par->names.out, par->operation ) != 1 ) {
    _DEALLOCATIONS_MEMORY_WITHMASKS_;
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate or initialize auxiliary images\n", proc );
    return( -1 );
  }




  /* process other images
   */
  for ( nimages=1; nimages<imageFileList->n; nimages++ ) {

    if ( _UpdateAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				 &imTotalSumSqr, &imTotalMask,
				 imageStructure[nimages], maskStructure[nimages],
				 par->operation ) != 1 ) {
      _DEALLOCATIONS_MEMORY_WITHMASKS_;
      if ( _verbose_ ) {
	fprintf( stderr, "%s: unable to update auxiliary images\n", proc );
	fprintf( stderr, "\t while processing image #%d '%s'\n", nimages, imageFileList->data[nimages] );
      }
      return( -1 );
    }

  }




  /* processing is over
   */

  if ( maskFileList->n > 0 ) {
    _FreeImageList( maskStructure, maskFileList );
    VT_Free( (void**)&maskStructure );
  }
  _FreeImageList( imageStructure, imageFileList );
  VT_Free( (void**)&imageStructure );


  
  /* last computations
   */
  if ( _ResultFromAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				   &imTotalSumSqr, &imTotalMask,
				   par->names.out, par->operation, par->type ) != 1 ) {
    VT_FreeImage ( &imTotalSumSqr );
    VT_FreeImage ( &imTotalMask );
    VT_FreeImage ( &imTotalSum );         
    VT_FreeImage ( &imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get result\n", proc );
    return( -1 );
  }

  return( 1 );
}







/**************************************************
 *
 * streaming computation
 *
 **************************************************/





static int _StreamingComputationWithoutMasks( stringList *imageFileList, 
					      local_par *par )
{
  char *proc = "_StreamingComputationWithoutMasks";

  vt_image *imLocal = (vt_image*)NULL;
  vt_image *imResult = (vt_image*)NULL;

  vt_image imTotalSum;
  vt_image imTotalSumSqr;
  vt_image imTotalRes;

  vt_image imLocalSum;
  vt_image imLocalSumSqr;
  vt_image imLocalRes;

  float *bufTotalSum = NULL;
  float *bufTotalSumSqr = NULL;

  int theDim[3] = {0,0,0};
  int winDim[3] = {0,0,0};
  int nOffsets[3]= {0,0,0};
  int pOffsets[3]= {0,0,0};

  int nimages;
  size_t v = 0;
  size_t i;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  /* some initializations
   */
  winDim[0] = par->window.x;
  winDim[1] = par->window.y;
  winDim[2] = par->window.z;
  nOffsets[0] = -(int)(winDim[0] / 2); pOffsets[0] = winDim[0] - 1 + nOffsets[0];
  nOffsets[1] = -(int)(winDim[1] / 2); pOffsets[1] = winDim[1] - 1 + nOffsets[1];
  nOffsets[2] = -(int)(winDim[2] / 2); pOffsets[2] = winDim[2] - 1 + nOffsets[2];

  VT_Image ( &imTotalSum );
  VT_Image ( &imTotalSumSqr );
  VT_Image ( &imTotalRes );

  VT_Image ( &imLocalSum );
  VT_Image ( &imLocalSumSqr );
  VT_Image ( &imLocalRes );

  /* we can go
     il y aurait sans doute des ameliorations a faire dans le cas d'une fenetre 1x1x1
   */

#define _DEALLOCATIONS_STREAMING_WITHOUTMASKS_ { \
  VT_FreeImage ( &imTotalSum );      \
  VT_FreeImage ( &imTotalSumSqr );   \
  VT_FreeImage ( &imTotalRes );	     \
                                     \
  VT_FreeImage ( &imLocalSum );      \
  VT_FreeImage ( &imLocalSumSqr );   \
  VT_FreeImage ( &imLocalRes );	     \
                                     \
  if ( imLocal != NULL ) {           \
    VT_FreeImage( imLocal );         \
    VT_Free( (void**)&imLocal );     \
  }                                  \
}


  /* loop over images
   */
  for ( nimages = 0; nimages < imageFileList->n; nimages++ ) {

    if ( _VT_VERBOSE_ || _verbose_ )
      fprintf( stderr, "%d: image='%s'\n", nimages, imageFileList->data[nimages] );
    
    /* read input image #nimages
     */
    imLocal = _VT_Inrimage( imageFileList->data[nimages] );
    if ( imLocal == (vt_image*)NULL ) {
      _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when opening %s\n", proc, imageFileList->data[nimages] );
      return( -1 );
    }



    /* first reading: allocations and initialisations of auxiliary images
     */
    if ( nimages == 0 ) {
      
      /* allocations
       */
      if ( _AllocAuxiliaryImages( &imTotalRes, &imTotalSum, &imTotalSumSqr,
				  imLocal, par->names.out, par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate global auxiliary images\n", proc);
	return( -1 );
      }

      if ( winDim[0] > 1 || winDim[1] > 1 || winDim[2] > 1 ) {
	if ( _AllocAuxiliaryImages( &imLocalRes, &imLocalSum, &imLocalSumSqr,
				    imLocal, par->names.out, par->operation ) != 1 ) {
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate local auxiliary images\n", proc);
	  return( -1 );
	}
      }



      theDim[0] = imLocal->dim.v * imLocal->dim.x;
      theDim[1] = imLocal->dim.y;
      theDim[2] = imLocal->dim.z;
      v = (size_t)theDim[0] * (size_t)theDim[1] * (size_t)theDim[2];



      /* initialisations
       */
      if ( winDim[0] == 1 && winDim[1] == 1 && winDim[2] == 1 ) {
	
	switch ( par->operation ) {
	default :
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: such operation not handled yet\n", proc );
	  return( -1 );
	case _MAX_ :
	case _MIN_ :
	  if ( ConvertBuffer( imLocal->buf, imLocal->type,
			      imTotalRes.buf, imTotalRes.type, v ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during conversion (1)\n", proc );
	    return( -1 );
	  }
	  break;
	case _STDDEV_ :
	case _VAR_ :
	  if ( sqrImage( imLocal->buf, imLocal->type,
			 imTotalSumSqr.buf, imTotalSumSqr.type,
			 theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during stddev/var init (1)\n", proc );
	    return( -1 );
	  }
	case _SUM_ :
	case _MEAN_ :
	  if ( ConvertBuffer( imLocal->buf, imLocal->type,
			      imTotalSum.buf, imTotalSum.type, v ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during conversion (1)\n", proc );
	    return( -1 );
	  }
	  break;
	}
	
      }
      else {

	switch ( par->operation ) {
	default :
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: such operation not handled yet\n", proc );
	  return( -1 );
	case _MAX_ :
	  if ( maxFiltering( imLocal->buf, imLocal->type,
			     imTotalRes.buf, imTotalRes.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during max filtering (init)\n", proc );
	    return( -1 );
	  }
	  break;
	case _MIN_ :
	  if ( minFiltering( imLocal->buf, imLocal->type,
			     imTotalRes.buf, imTotalRes.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during min filtering (init)\n", proc );
	    return( -1 );
	  }
	  break;
	case _STDDEV_ :
	case _VAR_ :
	  if ( sumSquaresFiltering( imLocal->buf, imLocal->type,
				    imTotalSumSqr.buf, imTotalSumSqr.type,
				    theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering (init)\n", proc );
	    return( -1 );
	  }
	  break;
	case _SUM_ :
	case _MEAN_ :
	  if ( sumFiltering( imLocal->buf, imLocal->type,
			     imTotalSum.buf, imTotalSum.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum filtering (init)\n", proc );
	    return( -1 );
	  }
	  break;
	}

      }
      
    }

    /* other images 
     */
    else {

      if ( winDim[0] == 1 && winDim[1] == 1 && winDim[2] == 1 ) {

	switch ( par->operation ) {
	default :
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: such operation not handled yet\n", proc );
	  return( -1 );
	case _MAX_ :
	  if ( maxImages( imLocal->buf, imLocal->type,
			  imTotalRes.buf, imTotalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during max filtering\n", proc );
	    return( -1 );
	  }
	  break;
	case _MIN_ :
	  if ( minImages( imLocal->buf, imLocal->type,
			  imTotalRes.buf, imTotalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during min filtering\n", proc );
	    return( -1 );
	  }
	  break;
	case _STDDEV_ :
	case _VAR_ :
	  if ( _addSqrImage( imTotalSumSqr.buf, imTotalSumSqr.type,
			     imLocal->buf, imLocal->type,
			     imTotalSumSqr.buf, imTotalSumSqr.type,
			     theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering\n", proc );
	    return( -1 );
	  }
	case _SUM_ :
	case _MEAN_ :
	  if ( addImages( imLocal->buf, imLocal->type,
			  imTotalSum.buf, imTotalSum.type,
			  imTotalSum.buf, imTotalSum.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum filtering\n", proc );
	    return( -1 );
	  }
	  break;
	}
	
      }
      else {
	
	switch ( par->operation ) {
	default :
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: such operation not handled yet\n", proc );
	  return( -1 );
	case _MAX_ :
	  if ( maxFiltering( imLocal->buf, imLocal->type,
			     imLocalRes.buf, imLocalRes.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local max filtering\n", proc );
	    return( -1 );
	  }
	  if ( maxImages( imLocalRes.buf, imLocalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during max filtering\n", proc );
	    return( -1 );
	  }
	  break;
	case _MIN_ :
	  if ( minFiltering( imLocal->buf, imLocal->type,
			     imLocalRes.buf, imLocalRes.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local min filtering\n", proc );
	    return( -1 );
	  }
	  if ( minImages( imLocalRes.buf, imLocalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  imTotalRes.buf, imTotalRes.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during min filtering\n", proc );
	    return( -1 );
	  }
	  break;
	case _STDDEV_ :
	case _VAR_ :
	  if ( sumSquaresFiltering( imLocal->buf, imLocal->type,
				    imLocalSumSqr.buf, imLocalSumSqr.type,
				    theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local sum of squares filtering\n", proc );
	    return( -1 );
	  }
	  if ( addImages( imLocalSumSqr.buf, imLocalSumSqr.type,
			  imTotalSumSqr.buf, imTotalSumSqr.type,
			  imTotalSumSqr.buf, imTotalSumSqr.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering\n", proc );
	    return( -1 );
	  }
	case _SUM_ :
	case _MEAN_ :
	  if ( sumFiltering( imLocal->buf, imLocal->type,
			     imLocalSum.buf, imLocalSum.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local sum filtering\n", proc );
	    return( -1 );
	  }
	  if ( addImages( imLocalSum.buf, imLocalSum.type,
			  imTotalSum.buf, imTotalSum.type,
			  imTotalSum.buf, imTotalSum.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum filtering\n", proc );
	    return( -1 );
	  }
	  break;
	}
	
      }

    }
    


    if ( imLocal != NULL ) {
      VT_FreeImage( imLocal );
      VT_Free( (void**)&imLocal );
    }                              

  }
  /* end of loop over images
   */



  /* processing is over
   */

  VT_FreeImage ( &imLocalSum ); 
  VT_FreeImage ( &imLocalSumSqr );
  VT_FreeImage ( &imLocalRes );	 



 /* last computations
   */
  switch ( par->operation ) {
  default :
    VT_FreeImage ( &imTotalSumSqr );
    VT_FreeImage ( &imTotalSum );         
    VT_FreeImage ( &imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  case _MAX_ :
  case _MIN_ :
    imResult = &imTotalRes;
    break;

  case _SUM_ :
    imResult = &imTotalSum;
    break;

  case _VAR_ :
    imResult = &imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum.buf;
    bufTotalSumSqr = (float*)imTotalSumSqr.buf;
    for ( i=0; i<v; i++ ) {
      bufTotalSum[i] /= (float)imageFileList->n;
      bufTotalSumSqr[i] /= (float)imageFileList->n;
      bufTotalSum[i] = bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i];
    }
    break;

  case _STDDEV_ :
    imResult = &imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum.buf;
    bufTotalSumSqr = (float*)imTotalSumSqr.buf;
    for ( i=0; i<v; i++ ) {
      bufTotalSum[i] /= (float)imageFileList->n;
      bufTotalSumSqr[i] /= (float)imageFileList->n;
      bufTotalSum[i] = sqrt( bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i] );
    }
    break;
    
  case _MEAN_ :
    imResult = &imTotalSum;
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    bufTotalSum = (float*)imTotalSum.buf;
    for ( i=0; i<v; i++ ) {
      bufTotalSum[i] /= (float)imageFileList->n;
    }
    break;
  }



  VT_FreeImage ( &imTotalSumSqr );



  if ( par->type == TYPE_UNKNOWN || par->type == imResult->type ) {
    if ( VT_WriteInrimage( imResult ) == -1 ) {
      VT_FreeImage ( &imTotalSum );         
      VT_FreeImage ( &imTotalRes );	
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
  }
  else {
     VT_InitFromImage( &imTotalSumSqr, imResult, par->names.out, par->type );
    if ( VT_AllocImage( &imTotalSumSqr ) != 1 ) {
      VT_FreeImage ( &imTotalSum );         
      VT_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate auxiliary result image\n", proc );
      return( -1 );
    }
    v = (size_t)imResult->dim.v * (size_t)imResult->dim.x * (size_t)imResult->dim.y * (size_t)imResult->dim.z;
    if ( ConvertBuffer( imResult->buf, imResult->type,
		   imTotalSumSqr.buf, imTotalSumSqr.type,
			v ) != 1 ) {
      VT_FreeImage ( &imTotalSumSqr );
      VT_FreeImage ( &imTotalSum );         
      VT_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to convert image\n", proc );
      return( -1 );
    }
    if ( VT_WriteInrimage( &imTotalSumSqr ) == -1 ) {
      VT_FreeImage ( &imTotalSumSqr );
      VT_FreeImage ( &imTotalSum );         
      VT_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
    VT_FreeImage ( &imTotalSumSqr );
  }

  VT_FreeImage ( &imTotalSum );         
  VT_FreeImage ( &imTotalRes );	

  return( 1 );
}





static int _StreamingComputationWithMasks( stringList *imageFileList, 
					   stringList *maskFileList,
					   local_par *par )
{
  char *proc = "_StreamingComputationWithMasks";

  vt_image *imLocal = (vt_image*)NULL;
  vt_image *imLocalMask = (vt_image*)NULL;

  vt_image imTotalSum;
  vt_image imTotalSumSqr;
  vt_image imTotalRes;
  vt_image imTotalMask;

  int nimages;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  VT_Image ( &imTotalSum );
  VT_Image ( &imTotalSumSqr );
  VT_Image ( &imTotalRes );
  VT_Image ( &imTotalMask );

  if ( maskFileList->n == 0 ) {
    return ( _StreamingComputationWithoutMasks( imageFileList, par ) );
  }



  if ( par->window.x != 1 || par->window.y != 1 || par->window.z != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: computation with window > 1x1x1 not handled yet\n", proc );
    return( -1 );
  }



  /* we can go
   */

#define _DEALLOCATIONS_STREAMING_WITHMASKS_ { \
  VT_FreeImage ( &imTotalSum );      \
  VT_FreeImage ( &imTotalSumSqr );   \
  VT_FreeImage ( &imTotalRes );	     \
  VT_FreeImage ( &imTotalMask );     \
                                     \
  if ( imLocal != NULL ) {           \
    VT_FreeImage( imLocal );         \
    VT_Free( (void**)&imLocal );     \
  }                                  \
  if ( imLocalMask != NULL ) {       \
    VT_FreeImage( imLocalMask );     \
    VT_Free( (void**)&imLocalMask ); \
  }                                  \
}



  /* loop over images
   */
  for ( nimages = 0; nimages < imageFileList->n; nimages++ ) {

    if ( _VT_VERBOSE_ || _verbose_ )
      fprintf( stderr, "%d: image='%s'\n", nimages, imageFileList->data[nimages] );
    
    /* read input image #nimages
     */
    imLocal = _VT_Inrimage( imageFileList->data[nimages] );
    if ( imLocal == (vt_image*)NULL ) {
      _DEALLOCATIONS_STREAMING_WITHMASKS_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when opening %s\n", proc, imageFileList->data[nimages] );
      return( -1 );
    }

    imLocalMask = _VT_Inrimage( maskFileList->data[nimages] );
    if ( imLocalMask == (vt_image*)NULL ) {
      _DEALLOCATIONS_STREAMING_WITHMASKS_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when opening %s\n", proc, maskFileList->data[nimages]  );
      return( -1 );
    }



    /* first reading : allocations of auxiliary images
     */
    if ( nimages == 0 ) {

      if ( _AllocInitAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				      &imTotalSumSqr, &imTotalMask,
				      imLocal, imLocalMask,
				      par->names.out, par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHMASKS_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate or initialize auxiliary images\n", proc );
	return( -1 );
      }

    }

    else {

      if ( _UpdateAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				   &imTotalSumSqr, &imTotalMask,
				   imLocal, imLocalMask,
				   par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHMASKS_;
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: unable to update auxiliary images\n", proc );
	  fprintf( stderr, "\t while processing image #%d '%s'\n", nimages, imageFileList->data[nimages] );
	}
	return( -1 );
      }

    }
    
    if ( imLocal != NULL ) {
      VT_FreeImage( imLocal );
      VT_Free( (void**)&imLocal );
    }                              
    if ( imLocalMask != NULL ) {       
      VT_FreeImage( imLocalMask );     
      VT_Free( (void**)&imLocalMask ); 
    }

  }
  /* end of loop over images
   */



  
  /* last computations
   */
  if ( _ResultFromAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				   &imTotalSumSqr, &imTotalMask,
				   par->names.out, par->operation, par->type ) != 1 ) {
    VT_FreeImage ( &imTotalSumSqr );
    VT_FreeImage ( &imTotalMask );
    VT_FreeImage ( &imTotalSum );         
    VT_FreeImage ( &imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get result\n", proc );
    return( -1 );
  }



  return( 1 );
}
