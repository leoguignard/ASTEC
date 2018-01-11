/*************************************************************************
 * blDriftQualityMeasure.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 14 mar 2014 19:37:28 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#include <math.h> /* clock() */
#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */

#include <convert.h>
#include <pixel-operation.h>
#include <local-operation.h>
#include <string-tools.h>

#include <bal-image.h>

static int _verbose_ = 1;
static int _debug_ = 0;


typedef enum {
  _STDDEV_,
  _VAR_
} typeOperation;

typedef enum {
  _MEMORY_,
  _STREAMING_
} typeComputation;


typedef struct local_par {

  char *theimagelist;
  char *resimagename;
  char *maskimagelist;

  bufferType type;

  char *nameformat;
  int firstindex;
  int lastindex;

  typeOperation operation;

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
static void _Parse( int argc, char *argv[], local_par *par );
static void _ErrorParse( char *str, int l );
static void _InitParam( local_par *par );
static double _GetTime();
static double _GetClock();




static char *usage = "[[-image-list|-list|-refl] %s]\n\
 [-format %s -f[irst] %d -l[ast] %d]\n\
 [[-res] image-out]\n\
 [-mask-list|-maskl %s]\n\
 [-var|-stddev]\n\
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

static char *program;






int main( int argc, char *argv[] )
{
  local_par par;
  double time_init, time_exit;
  double clock_init, clock_exit;

  stringList imageFileList, maskFileList;


  time_init = _GetTime();
  clock_init = _GetClock();


  /*--- initialisation des parametres ---*/
  _InitParam( &par );
  
  /*--- lecture des parametres ---*/
  _Parse( argc, argv, &par );
  

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
      _ErrorParse( "unable to read input image list\n", 0);
    }
  }
  else if ( buildStringListFromFile( par.theimagelist, &imageFileList ) != 1 ) {
    _ErrorParse( "unable to read input image list\n", 0);
  }

  if ( 0 ) printStringList( stderr, &imageFileList, "Input images" );

  if ( par.maskimagelist != '\0' ) {
    if ( buildStringListFromFile( par.maskimagelist, &maskFileList ) != 1 ) {
      freeStringList( &imageFileList );
      _ErrorParse( "unable to read input mask list\n", 0);
    }
    if ( maskFileList.n != imageFileList.n ) {
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      _ErrorParse( "image and mask lists have different length\n", 0);
    }
  }

  

  if ( par.maskimagelist != '\0' ) {
    
    switch ( par.computation ) {
    default :
      freeStringList( &maskFileList );
      freeStringList( &imageFileList );
      _ErrorParse( "such computation type not handled yet\n", 0);
    case  _MEMORY_ :
      if ( _InMemoryComputationWithMasks( &imageFileList, &maskFileList, &par ) != 1 ) {
	freeStringList( &maskFileList );
	freeStringList( &imageFileList );
	_ErrorParse( "error when computing (memory with masks case)\n", 0);
      }
      break;
    case _STREAMING_ :
      if ( _StreamingComputationWithMasks( &imageFileList, &maskFileList, &par ) != 1 ) {
	freeStringList( &maskFileList );
	freeStringList( &imageFileList );
	_ErrorParse( "error when computing (streaming with masks case)\n", 0);
      }
      break;
    }
    
  }
  else {
    
    switch ( par.computation ) {
    default :
        freeStringList( &imageFileList );
      _ErrorParse( "such computation type not handled yet\n", 0);
    case  _MEMORY_ :
       if ( _InMemoryComputationWithoutMasks( &imageFileList, &par ) != 1 ) {
	freeStringList( &imageFileList );
	_ErrorParse( "error when computing (memory without masks case)\n", 0);
      }
      break;
    case _STREAMING_ :
      if ( _StreamingComputationWithoutMasks( &imageFileList, &par ) != 1 ) {
	freeStringList( &imageFileList );
	_ErrorParse( "error when computing (streaming without masks case)\n", 0);
      }
      break;
    }
    
  }
  


  freeStringList( &maskFileList );
  freeStringList( &imageFileList );









  /* end
   */
  
  time_exit = _GetTime();
  clock_exit = _GetClock();

  if (  par.print_time ) {
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );
  }

  return( 1 );
}








static void _Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;
  int o=0, s=0, r=0;
  char text[128];
  
  program = argv[0];

  if ( argc == 1 ) _ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {

      /*--- arguments generaux ---*/
      if ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) {
	_ErrorParse( "\n", 0 );
      }
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	_verbose_ ++;
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
	if ( _debug_ <= 0 ) _debug_ = 1;
	else _debug_ ++;
      }


      /*---  images ---*/
      else if ( strcmp ( argv[i], "-image-list" ) == 0  
		|| (strcmp ( argv[i], "-list" ) == 0 && argv[i][5] == '\0') 
		|| (strcmp ( argv[i], "-refl" ) == 0 && argv[i][5] == '\0') ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -image-list...\n", 0 );
	if ( par->theimagelist[0] != '\0' ) 
	  _ErrorParse( "parsing -image-list: input has already been parsed ...\n", 0 );
	par->theimagelist = argv[i];
	inputisread = 1;
      }
      else if ( strcmp ( argv[i], "-mask-list" ) == 0  
		|| (strcmp ( argv[i], "-maskl" ) == 0  && argv[i][6] == '\0') ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -mask-list...\n", 0 );
	par->maskimagelist = argv[i];
      }
      else if ( strcmp ( argv[i], "-res" ) == 0  && argv[i][4] == '\0' ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -res...\n", 0 );
	if ( par->resimagename[0] != '\0' ) 
	  _ErrorParse( "parsing -res: output has already been parsed ...\n", 0 );
	par->resimagename = argv[i];
	outputisread = 1;
      }


      else if ( strcmp ( argv[i], "-format" ) == 0 ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -format...\n", 0 );
	par->nameformat = argv[i];
	inputisread = 1;
      }
      else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -first ...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->firstindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -first ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
	i += 1;
	if ( i >= argc) _ErrorParse( "parsing -last ...\n", 0 );
	status = sscanf( argv[i], "%d", &(par->lastindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -last ...", 0 );
      }



      else if ( strcmp ( argv[i], "-stddev" ) == 0 ) {
	par->operation = _STDDEV_;
      }
      else if ( strcmp ( argv[i], "-var" ) == 0 ) {
	par->operation = _VAR_;
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
	if ( i >= argc)    _ErrorParse( "parsing -o...\n", 0 );
	status = sscanf( argv[i],"%d",&o );
	if ( status <= 0 ) _ErrorParse( "parsing -o...\n", 0 );
      }
      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 ) {
	par->theimagelist = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	par->resimagename = argv[i];
	outputisread = 1;
      }
      else 
	_ErrorParse( "too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if ( inputisread == 0 )
    _ErrorParse( "no input file names when parsing\n", 0 );
  if ( outputisread ==  0 )
    _ErrorParse( "no output file names when parsing\n", 0 );
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) _Warning("no specified type", program); */
}






static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void _InitParam( local_par *par )
{
  par->theimagelist = (char*)NULL;
  par->resimagename = (char*)NULL;
  par->maskimagelist = (char*)NULL;

  par->type = TYPE_UNKNOWN;

  par->nameformat = (char*)NULL;
  par->firstindex = 0;
  par->lastindex = 0;

  par->operation = _STDDEV_;

  par->computation = _MEMORY_;

  par->print_time = 0;
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






























/**************************************************
 *
 * image array
 *
 **************************************************/

static void _FreeImageList( bal_image *array, stringList *list ) 
{
  int i;

  if ( array == (bal_image*)NULL ) return;
  if ( list->n <= 0 ) return;

  for ( i=0; i<list->n; i++ ) {
    BAL_FreeImage( &(array[i]) );
  }
}




bal_image * _ReadImageList( stringList *list ) 
{
  char *proc = "_ReadImageList";
  bal_image *a = (bal_image*)NULL;
  int i;
  
  if ( list->n <= 0 ) return( (bal_image*)NULL );

  a = (bal_image*)malloc( list->n * sizeof( bal_image ) );
  if ( a == (bal_image*)NULL ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate array\n", proc );
    return( (bal_image*)NULL );
  }
  
  for ( i=0; i<list->n; i++ ) 
    BAL_FreeImage( &(a[i]) );

  for ( i=0; i<list->n; i++ ) {
    if ( BAL_ReadImage( &(a[i]), list->data[i], 0 ) != 1 ) {
      _FreeImageList( a, list );
      free( a );
      if ( _verbose_ ) 
	fprintf( stderr, "%s: unable to read image #%d '%s'\n", proc, i, list->data[i] );
      return( (bal_image*)NULL );
     }
  }
  
  for ( i=1; i<list->n; i++ ) {
    if ( a[i].vdim != a[0].vdim
	 || a[i].ncols != a[0].ncols
	 || a[i].nrows != a[0].nrows
	 || a[i].nplanes != a[0].nplanes ) {
      _FreeImageList( a, list );
      free( a );
       if ( _verbose_ ) 
	 fprintf( stderr, "%s: image #%d '%s' has different dimensions than image #0 '%s'\n", 
		  proc, i, list->data[i], list->data[0]);
       return( (bal_image*)NULL );
    }
    if ( a[i].type != a[0].type ) {
      _FreeImageList( a, list );
      free( a );
       if ( _verbose_ ) 
	 fprintf( stderr, "%s: image #%d '%s' has a different type than image #0 '%s'\n", 
		  proc, i, list->data[i], list->data[0]);
       return( (bal_image*)NULL );
    }
  } 

  return( a );
}







/**************************************************
 *
 * common procedures
 *
 **************************************************/





static int _AllocAuxiliaryImages( bal_image *imRes,
				  bal_image *imSum,
				  bal_image *imSumSqr,
				  bal_image *image,
				  char *name,
				  typeOperation operation )
{
  char * proc = "_AllocAuxiliaryImages";
  
  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
  
  case _STDDEV_ :
  case _VAR_ :
    BAL_InitImage( imSumSqr, name,
		   image->ncols, image->nrows, image->nplanes, image->vdim, FLOAT );
    if ( BAL_AllocImage( imSumSqr ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate sum of squares image\n", proc );
      return( -1 );
    }

    BAL_InitImage( imSum, name,
		   image->ncols, image->nrows, image->nplanes, image->vdim, FLOAT );
    if ( BAL_AllocImage( imSum ) != 1 ) {
      BAL_FreeImage( imSumSqr );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate sum image\n", proc );
      return( -1 );
    }
    break;
  }
  
  return( 1 );
}





static int _AllocInitAuxiliaryWithMasks( bal_image *imTotalRes,
					 bal_image *imTotalSum,
					 bal_image *imTotalSumSqr,
					 bal_image *imTotalMask,
					 bal_image *image,
					 bal_image *mask,
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

  BAL_InitImage( imTotalMask, name,
		 image->ncols, image->nrows, image->nplanes, image->vdim, USHORT );
  if ( BAL_AllocImage( imTotalMask ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate sum of masks image\n", proc );
    return( -1 );
  }
  bufTotalMask = (unsigned short int *)imTotalMask->data;



  theDim[0] = image->vdim * image->ncols;
  theDim[1] = image->nrows;
  theDim[2] = image->nplanes;
  v = (size_t)theDim[0] * (size_t)theDim[1] * (size_t)theDim[2];



  /* initialisation of auxiliary images with the first image
   */
  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );

  case _STDDEV_ :
  case _VAR_ :
    if ( sqrImage( image->data, image->type,
		   imTotalSumSqr->data, imTotalSumSqr->type,
		   theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during stddev/var init (1)\n", proc );
      return( -1 );
    }
    if ( maskImage( imTotalSumSqr->data, imTotalSumSqr->type,
		    mask->data, mask->type,
		    imTotalSumSqr->data, imTotalSumSqr->type,
		    theDim ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during stddev/var init (2)\n", proc );
      return( -1 );
    }

    if ( ConvertBuffer( image->data, image->type,
			imTotalSum->data, imTotalSum->type, v ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: error during conversion (2)\n", proc );
      return( -1 );
    }
    if ( maskImage( imTotalSum->data, imTotalSum->type,
		    mask->data, mask->type,
		    imTotalSum->data, imTotalSum->type,
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
      u8 *theMask = (u8*)mask->data;
      for ( i=0; i<v; i++ )
	bufTotalMask[i] = (theMask[i] > 0) ? 1 : 0;
    }
    break;
  case USHORT :
    {
      u16 *theMask = (u16*)mask->data;
      for ( i=0; i<v; i++ )
	bufTotalMask[i] = (theMask[i] > 0) ? 1 : 0;
    }
    break;
  }

  return( 1 );
}





static int _UpdateAuxiliaryWithMasks( bal_image *imTotalRes,
				   bal_image *imTotalSum,
				   bal_image *imTotalSumSqr,
				   bal_image *imTotalMask,
				   bal_image *image,
				   bal_image *mask,
				   typeOperation operation )
{
  char * proc = "_UpdateAuxiliaryWithMasks";

  float *bufTotalSum = NULL;
  float *bufTotalSumSqr = NULL;
  unsigned short int *bufTotalMask = (unsigned short int *)imTotalMask->data;

  size_t v = 0;
  size_t i;


  switch ( operation ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
    
  case _STDDEV_ :
  case _VAR_ :
    bufTotalSumSqr = (float*)imTotalSumSqr->data;
    bufTotalSum = (float*)imTotalSum->data;
    break;
  }



  v = (size_t)image->vdim * (size_t)image->ncols * (size_t)image->nrows * (size_t)image->nplanes;
  


#define _SQR_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->data;     \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    bufTotalSumSqr[i] += (float)theBuf[i] * (float)theBuf[i]; \
    bufTotalSum[i] += (float)theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}

#define _SUM_WITHMASKS_( TYPE ) {       \
  TYPE *theBuf = (TYPE*)image->data;     \
  for ( i=0; i<v; i++ ) {               \
    if ( theMask[i] == 0 ) continue;    \
    bufTotalSum[i] += (float)theBuf[i]; \
    bufTotalMask[i] ++;                 \
  }                                     \
}


#define _OPERATION_WITHMASKS_( TYPEM ) { \
  TYPEM *theMask = (TYPEM*)mask->data; \
  switch ( operation ) {              \
  default :                           \
    if ( _verbose_ )                  \
      fprintf( stderr, "%s: such operation not handled yet\n", proc ); \
    return( -1 );                     \
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





static int _ResultFromAuxiliaryWithMasks( bal_image *imTotalRes,
				       bal_image *imTotalSum,
				       bal_image *imTotalSumSqr,
				       bal_image *imTotalMask,
				       char *name,
				       typeOperation operation,
				       bufferType type )
{
  char * proc = "_ResultFromAuxiliaryWithMasks";

  bal_image *imResult = (bal_image*)NULL;

  float *bufTotalSum = NULL;
  float *bufTotalSumSqr = NULL;
  unsigned short int *bufTotalMask = (unsigned short int *)imTotalMask->data ;

  size_t v = 0;
  size_t i;


  switch ( operation ) {
  default :
    BAL_FreeImage ( imTotalSumSqr );
    BAL_FreeImage ( imTotalMask );
    BAL_FreeImage ( imTotalSum );         
    BAL_FreeImage ( imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );

  case _VAR_ :
    imResult = imTotalSum;
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    bufTotalSum = (float*)imTotalSum->data;
    bufTotalSumSqr = (float*)imTotalSumSqr->data;
    for ( i=0; i<v; i++ ) {
      if ( bufTotalMask[i] == 0 ) continue;
      bufTotalSum[i] /= (float)bufTotalMask[i];
      bufTotalSumSqr[i] /= (float)bufTotalMask[i];
      bufTotalSum[i] = bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i];
    }
    break;

  case _STDDEV_ :
    imResult = imTotalSum;
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    bufTotalSum = (float*)imTotalSum->data;
    bufTotalSumSqr = (float*)imTotalSumSqr->data;
    for ( i=0; i<v; i++ ) {
      if ( bufTotalMask[i] == 0 ) continue;
      bufTotalSum[i] /= (float)bufTotalMask[i];
      bufTotalSumSqr[i] /= (float)bufTotalMask[i];
      bufTotalSum[i] = sqrt( bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i] );
    }
    break;
    
  }



  BAL_FreeImage ( imTotalSumSqr );
  BAL_FreeImage ( imTotalMask );



  if ( type == TYPE_UNKNOWN || type == imResult->type ) {
    if ( BAL_WriteImage( imResult, imResult->name ) == -1 ) {
      BAL_FreeImage ( imTotalSum );         
      BAL_FreeImage ( imTotalRes );	
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
  }
  else {
    BAL_InitImage( imTotalMask, name,
		    imResult->ncols, imResult->nrows, imResult->nplanes, imResult->vdim, type );
    if ( BAL_AllocImage( imTotalMask ) != 1 ) {
      BAL_FreeImage ( imTotalSum );         
      BAL_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate auxiliary result image\n", proc );
      return( -1 );
    }
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    if ( ConvertBuffer( imResult->data, imResult->type,
			imTotalMask->data, imTotalMask->type,
			v ) != 1 ) {
      BAL_FreeImage ( imTotalMask );
      BAL_FreeImage ( imTotalSum );         
      BAL_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to convert image\n", proc );
      return( -1 );
    }
    if ( BAL_WriteImage( imTotalMask, imTotalMask->name ) == -1 ) {
      BAL_FreeImage ( imTotalMask );
      BAL_FreeImage ( imTotalSum );         
      BAL_FreeImage ( imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
    BAL_FreeImage ( imTotalMask );
  }

  BAL_FreeImage ( imTotalSum );         
  BAL_FreeImage ( imTotalRes );	

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

  bal_image *imageStructure = (bal_image*)NULL;

  int theDim[3] = {0,0,0};
  int winDim[3] = {0,0,0};

  bufferType resType;
  bal_image imTotalRes;
  int i, r;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  imageStructure = _ReadImageList( imageFileList );
  if ( imageStructure == (bal_image*)NULL ) {
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
      free( imageStructure );
      imageStructure = (bal_image*)NULL;
      if ( _verbose_ )
	fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );
    case _STDDEV_ :
    case _VAR_ :
      resType = FLOAT;
      break;
    }
  }

  BAL_InitImage( &imTotalRes, par->resimagename,
		 imageStructure[0].ncols, imageStructure[0].nrows, 
		 imageStructure[0].nplanes, imageStructure[0].vdim, resType );
  if ( BAL_AllocImage( &imTotalRes ) != 1 ) {
    _FreeImageList( imageStructure, imageFileList );
    free( imageStructure );
      imageStructure = (bal_image*)NULL;
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate result image\n", proc);
    return( -1 );
  }

  theDim[0] = imageStructure[0].vdim * imageStructure[0].ncols;
  theDim[1] = imageStructure[0].nrows;
  theDim[2] = imageStructure[0].nplanes;
  winDim[0] = imageStructure[0].vdim;
  winDim[1] = 1;
  winDim[2] = 1;


#define _OPERATION_MEMORY_WITHOUTMASKS_( TYPE ) {                       \
  TYPE **buffers = (TYPE**)NULL;                                        \
  buffers = (TYPE**)malloc( imageFileList->n * sizeof( TYPE* ) );       \
  if ( buffers == (TYPE**)NULL ) {                                      \
    BAL_FreeImage ( &imTotalRes );	                                \
    _FreeImageList( imageStructure, imageFileList );                    \
    free( imageStructure );                                             \
    imageStructure = (bal_image*)NULL;                                  \
    if ( _verbose_ )                                                    \
      fprintf( stderr, "%s: unable to allocate buffer array\n", proc);  \
    return( -1 );                                                       \
  }                                                                     \
  for ( i=0; i<imageFileList->n; i++ )                                  \
    buffers[i] = (TYPE*)(imageStructure[i].data );                       \
  switch ( par->operation ) {                                           \
  default : r = 0; break;                                               \
  case _STDDEV_ :                                                       \
    r = stddevFilteringBuffers( (void**)buffers, imageStructure[0].type, \
				imageFileList->n, imTotalRes.data, imTotalRes.type, \
				theDim, winDim );                       \
    break;                                                              \
  case _VAR_ :                                                          \
    r = varFilteringBuffers( (void**)buffers, imageStructure[0].type,  \
			     imageFileList->n, imTotalRes.data, imTotalRes.type, \
			     theDim, winDim );                          \
    break;                                                              \
  }                                                                     \
  free( buffers );                                                      \
  if ( r != 1 ) {                                                       \
    BAL_FreeImage ( &imTotalRes );	                                \
    _FreeImageList( imageStructure, imageFileList );                    \
    free( imageStructure );                                 \
    if ( _verbose_ )                                                    \
      fprintf( stderr, "%s: error when filtering buffer array\n", proc); \
    return( -1 );                                                       \
  }	                                                                \
} 

  switch( imageStructure[0].type ) {
  default :
    BAL_FreeImage ( &imTotalRes );	
    _FreeImageList( imageStructure, imageFileList );
    free( imageStructure );
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
  free( imageStructure ); 

  if ( BAL_WriteImage( &imTotalRes, par->resimagename ) == -1 ) {
    BAL_FreeImage ( &imTotalRes );	
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to write output image\n", proc );
    return( -1 );
  }

  BAL_FreeImage ( &imTotalRes );	

  return( 1 );
}





static int _InMemoryComputationWithMasks( stringList *imageFileList, 
					  stringList *maskFileList,
					  local_par *par )
{
  char *proc = "_InMemoryComputationWithMasks";

  bal_image *imageStructure = (bal_image*)NULL;
  bal_image *maskStructure = (bal_image*)NULL;
  
  bal_image imTotalSum;
  bal_image imTotalSumSqr;
  bal_image imTotalRes;
  bal_image imTotalMask;

  int nimages;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  BAL_InitImage ( &imTotalSum, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalSumSqr, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalRes, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalMask, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  
  imageStructure = _ReadImageList( imageFileList );
  if ( imageStructure == (bal_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error when reading images\n", proc );
    return( -1 );
  }

  if ( maskFileList->n > 0 ) {
    maskStructure = _ReadImageList( maskFileList );
    if ( maskStructure == (bal_image*)NULL ) {
      _FreeImageList( imageStructure, imageFileList );
      free( imageStructure );
      if ( _verbose_ )
	fprintf( stderr, "%s: error when reading masks\n", proc );
      return( -1 );
    }
  }
  else {
    _FreeImageList( imageStructure, imageFileList );
    free( imageStructure );
    return ( _InMemoryComputationWithoutMasks( imageFileList, par ) );
  }

  if ( maskStructure[0].vdim != imageStructure[0].vdim 
       || maskStructure[0].ncols != imageStructure[0].ncols 
       || maskStructure[0].nrows != imageStructure[0].nrows 
       || maskStructure[0].nplanes != imageStructure[0].nplanes ) { 
    if ( maskFileList->n > 0 ) {          
      _FreeImageList( maskStructure, maskFileList );   
      free( maskStructure ); 
    }                                     
    _FreeImageList( imageStructure, imageFileList );     
    free( imageStructure );
    if ( _verbose_ )
      fprintf( stderr, "%s: images and masks have different dimensions\n", proc );
    return( -1 );
  }




  /* we can go
   */

#define _DEALLOCATIONS_MEMORY_WITHMASKS_ { \
  BAL_FreeImage ( &imTotalSum );         \
  BAL_FreeImage ( &imTotalSumSqr );      \
  BAL_FreeImage ( &imTotalRes );		\
  BAL_FreeImage ( &imTotalMask );        \
  if ( maskFileList->n > 0 ) {          \
    _FreeImageList( maskStructure, maskFileList ); \
    free( maskStructure );			   \
  }                                     \
  _FreeImageList( imageStructure, imageFileList ); \
  free( imageStructure );	\
}


  
  /* allocation and initialisation of auxiliary images
   */
  if ( _AllocInitAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				     &imTotalSumSqr, &imTotalMask,
				     &(imageStructure[0]), &(maskStructure[0]),
				     par->resimagename, par->operation ) != 1 ) {
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
				    &(imageStructure[nimages]), &(maskStructure[nimages]),
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
    free( maskStructure );
  }
  _FreeImageList( imageStructure, imageFileList );
  free( imageStructure );


  
  /* last computations
   */
  if ( _ResultFromAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				      &imTotalSumSqr, &imTotalMask,
				      par->resimagename, par->operation, par->type ) != 1 ) {
    BAL_FreeImage ( &imTotalSumSqr );
    BAL_FreeImage ( &imTotalMask );
    BAL_FreeImage ( &imTotalSum );         
    BAL_FreeImage ( &imTotalRes );
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

  bal_image imLocal;
  bal_image *imResult = (bal_image*)NULL;

  bal_image imTotalSum;
  bal_image imTotalSumSqr;
  bal_image imTotalRes;

  bal_image imLocalSum;
  bal_image imLocalSumSqr;
  bal_image imLocalRes;

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
  winDim[0] = 1;
  winDim[1] = 1;
  winDim[2] = 1;
  nOffsets[0] = -(int)(winDim[0] / 2); pOffsets[0] = winDim[0] - 1 + nOffsets[0];
  nOffsets[1] = -(int)(winDim[1] / 2); pOffsets[1] = winDim[1] - 1 + nOffsets[1];
  nOffsets[2] = -(int)(winDim[2] / 2); pOffsets[2] = winDim[2] - 1 + nOffsets[2];

  BAL_InitImage ( &imLocal, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

  BAL_InitImage ( &imTotalSum, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalSumSqr, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalRes, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

  BAL_InitImage ( &imLocalSum, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imLocalSumSqr, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imLocalRes, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

  /* we can go
     il y aurait sans doute des ameliorations a faire dans le cas d'une fenetre 1x1x1
   */

#define _DEALLOCATIONS_STREAMING_WITHOUTMASKS_ { \
  BAL_FreeImage ( &imTotalSum );      \
  BAL_FreeImage ( &imTotalSumSqr );   \
  BAL_FreeImage ( &imTotalRes );      \
                                      \
  BAL_FreeImage ( &imLocalSum );      \
  BAL_FreeImage ( &imLocalSumSqr );   \
  BAL_FreeImage ( &imLocalRes );      \
                                      \
  BAL_FreeImage( &imLocal );          \
}


  /* loop over images
   */
  for ( nimages = 0; nimages < imageFileList->n; nimages++ ) {

    if ( _verbose_ )
      fprintf( stderr, "%d: image='%s'\n", nimages, imageFileList->data[nimages] );
    
    /* read input image #nimages
     */
    
    if ( BAL_ReadImage( &imLocal, imageFileList->data[nimages], 0 ) != 1 ) {
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
				  &imLocal, par->resimagename, par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate global auxiliary images\n", proc);
	return( -1 );
      }

      if ( winDim[0] > 1 || winDim[1] > 1 || winDim[2] > 1 ) {
	if ( _AllocAuxiliaryImages( &imLocalRes, &imLocalSum, &imLocalSumSqr,
				    &imLocal, par->resimagename, par->operation ) != 1 ) {
	  _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate local auxiliary images\n", proc);
	  return( -1 );
	}
      }



      theDim[0] = imLocal.vdim * imLocal.ncols;
      theDim[1] = imLocal.nrows;
      theDim[2] = imLocal.nplanes;
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
	case _STDDEV_ :
	case _VAR_ :
	  if ( sqrImage( imLocal.data, imLocal.type,
			 imTotalSumSqr.data, imTotalSumSqr.type,
			 theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during stddev/var init (1)\n", proc );
	    return( -1 );
	  }
	  if ( ConvertBuffer( imLocal.data, imLocal.type,
			      imTotalSum.data, imTotalSum.type, v ) != 1 ) {
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
	case _STDDEV_ :
	case _VAR_ :
	  if ( sumSquaresFiltering( imLocal.data, imLocal.type,
				    imTotalSumSqr.data, imTotalSumSqr.type,
				    theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering (init)\n", proc );
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
	case _STDDEV_ :
	case _VAR_ :
	  if ( _addSqrImage( imTotalSumSqr.data, imTotalSumSqr.type,
			     imLocal.data, imLocal.type,
			     imTotalSumSqr.data, imTotalSumSqr.type,
			     theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering\n", proc );
	    return( -1 );
	  }
	  if ( addImages( imLocal.data, imLocal.type,
			  imTotalSum.data, imTotalSum.type,
			  imTotalSum.data, imTotalSum.type,
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
	case _STDDEV_ :
	case _VAR_ :
	  if ( sumSquaresFiltering( imLocal.data, imLocal.type,
				    imLocalSumSqr.data, imLocalSumSqr.type,
				    theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local sum of squares filtering\n", proc );
	    return( -1 );
	  }
	  if ( addImages( imLocalSumSqr.data, imLocalSumSqr.type,
			  imTotalSumSqr.data, imTotalSumSqr.type,
			  imTotalSumSqr.data, imTotalSumSqr.type,
			  theDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during sum of squares filtering\n", proc );
	    return( -1 );
	  }
	  if ( sumFiltering( imLocal.data, imLocal.type,
			     imLocalSum.data, imLocalSum.type,
			     theDim, winDim ) != 1 ) {
	    _DEALLOCATIONS_STREAMING_WITHOUTMASKS_;
	    if ( _verbose_ )
	      fprintf( stderr, "%s: error during local sum filtering\n", proc );
	    return( -1 );
	  }
	  if ( addImages( imLocalSum.data, imLocalSum.type,
			  imTotalSum.data, imTotalSum.type,
			  imTotalSum.data, imTotalSum.type,
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
    


    BAL_FreeImage( &imLocal );

  }
  /* end of loop over images
   */



  /* processing is over
   */

  BAL_FreeImage ( &imLocalSum ); 
  BAL_FreeImage ( &imLocalSumSqr );
  BAL_FreeImage ( &imLocalRes );	 



 /* last computations
   */
  switch ( par->operation ) {
  default :
    BAL_FreeImage ( &imTotalSumSqr );
    BAL_FreeImage ( &imTotalSum );         
    BAL_FreeImage ( &imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: such operation not handled yet\n", proc );
    return( -1 );

  case _VAR_ :
    imResult = &imTotalSum;
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    bufTotalSum = (float*)imTotalSum.data;
    bufTotalSumSqr = (float*)imTotalSumSqr.data;
    for ( i=0; i<v; i++ ) {
      bufTotalSum[i] /= (float)imageFileList->n;
      bufTotalSumSqr[i] /= (float)imageFileList->n;
      bufTotalSum[i] = bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i];
    }
    break;

  case _STDDEV_ :
    imResult = &imTotalSum;
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    bufTotalSum = (float*)imTotalSum.data;
    bufTotalSumSqr = (float*)imTotalSumSqr.data;
    for ( i=0; i<v; i++ ) {
      bufTotalSum[i] /= (float)imageFileList->n;
      bufTotalSumSqr[i] /= (float)imageFileList->n;
      bufTotalSum[i] = sqrt( bufTotalSumSqr[i] - bufTotalSum[i] * bufTotalSum[i] );
    }
    break;
    
  }



  BAL_FreeImage ( &imTotalSumSqr );



  if ( par->type == TYPE_UNKNOWN || par->type == imResult->type ) {
    if ( BAL_WriteImage( imResult, imResult->name ) == -1 ) {
      BAL_FreeImage ( &imTotalSum );         
      BAL_FreeImage ( &imTotalRes );	
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
  }
  else {
    BAL_InitImage( &imTotalSumSqr, par->resimagename,
		   imResult->ncols, imResult->nrows, imResult->nplanes, imResult->vdim, par->type );
    if ( BAL_AllocImage( &imTotalSumSqr ) != 1 ) {
      BAL_FreeImage ( &imTotalSum );         
      BAL_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate auxiliary result image\n", proc );
      return( -1 );
    }
    v = (size_t)imResult->vdim * (size_t)imResult->ncols * (size_t)imResult->nrows * (size_t)imResult->nplanes;
    if ( ConvertBuffer( imResult->data, imResult->type,
		   imTotalSumSqr.data, imTotalSumSqr.type,
			v ) != 1 ) {
      BAL_FreeImage ( &imTotalSumSqr );
      BAL_FreeImage ( &imTotalSum );         
      BAL_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to convert image\n", proc );
      return( -1 );
    }
    if ( BAL_WriteImage( &imTotalSumSqr, imTotalSumSqr.name ) == -1 ) {
      BAL_FreeImage ( &imTotalSumSqr );
      BAL_FreeImage ( &imTotalSum );         
      BAL_FreeImage ( &imTotalRes );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write output image\n", proc );
      return( -1 );
    }
    BAL_FreeImage ( &imTotalSumSqr );
  }

  BAL_FreeImage ( &imTotalSum );         
  BAL_FreeImage ( &imTotalRes );	

  return( 1 );
}





static int _StreamingComputationWithMasks( stringList *imageFileList, 
					   stringList *maskFileList,
					   local_par *par )
{
  char *proc = "_StreamingComputationWithMasks";

  bal_image imLocal;
  bal_image imLocalMask;

  bal_image imTotalSum;
  bal_image imTotalSumSqr;
  bal_image imTotalRes;
  bal_image imTotalMask;

  int nimages;

  if ( _debug_ )
    fprintf( stderr, " ... entering %s\n", proc );

  BAL_InitImage ( &imLocal, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imLocalMask, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

  BAL_InitImage ( &imTotalSum, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalSumSqr, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalRes, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );
  BAL_InitImage ( &imTotalMask, (char*)NULL, 0, 0, 0, 0, TYPE_UNKNOWN );

  if ( maskFileList->n == 0 ) {
    return ( _StreamingComputationWithoutMasks( imageFileList, par ) );
  }



  /* we can go
   */

#define _DEALLOCATIONS_STREAMING_WITHMASKS_ { \
  BAL_FreeImage ( &imTotalSum );      \
  BAL_FreeImage ( &imTotalSumSqr );   \
  BAL_FreeImage ( &imTotalRes );      \
  BAL_FreeImage ( &imTotalMask );     \
                                      \
  BAL_FreeImage( &imLocal );          \
  BAL_FreeImage( &imLocalMask );      \
}



  /* loop over images
   */
  for ( nimages = 0; nimages < imageFileList->n; nimages++ ) {

    if ( _verbose_ )
      fprintf( stderr, "%d: image='%s'\n", nimages, imageFileList->data[nimages] );
    
    /* read input image #nimages
     */
    if ( BAL_ReadImage( &imLocal, imageFileList->data[nimages], 0 ) != 1 ) {
      _DEALLOCATIONS_STREAMING_WITHMASKS_;
      if ( _verbose_ )
	fprintf( stderr, "%s: error when opening %s\n", proc, imageFileList->data[nimages] );
      return( -1 );
    }

    if ( BAL_ReadImage( &imLocalMask, maskFileList->data[nimages], 0 ) != 1 ) {
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
				      &imLocal, &imLocalMask,
				      par->resimagename, par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHMASKS_;
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate or initialize auxiliary images\n", proc );
	return( -1 );
      }

    }

    else {

      if ( _UpdateAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				   &imTotalSumSqr, &imTotalMask,
				   &imLocal, &imLocalMask,
				   par->operation ) != 1 ) {
	_DEALLOCATIONS_STREAMING_WITHMASKS_;
	if ( _verbose_ ) {
	  fprintf( stderr, "%s: unable to update auxiliary images\n", proc );
	  fprintf( stderr, "\t while processing image #%d '%s'\n", nimages, imageFileList->data[nimages] );
	}
	return( -1 );
      }

    }
    
    BAL_FreeImage( &imLocal );
    BAL_FreeImage( &imLocalMask );     

  }
  /* end of loop over images
   */



  
  /* last computations
   */
  if ( _ResultFromAuxiliaryWithMasks( &imTotalRes, &imTotalSum,
				   &imTotalSumSqr, &imTotalMask,
				   par->resimagename, par->operation, par->type ) != 1 ) {
    BAL_FreeImage ( &imTotalSumSqr );
    BAL_FreeImage ( &imTotalMask );
    BAL_FreeImage ( &imTotalSum );         
    BAL_FreeImage ( &imTotalRes );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to get result\n", proc );
    return( -1 );
  }



  return( 1 );
}
