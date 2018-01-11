/*************************************************************************
 * blPlotTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 18 fév 2014 11:28:35 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

#include <math.h>

static int _verbose_ = 1;
static int _debug_ = 0;

#include <extract.h>
#include <bal-image.h>
#include <bal-biolib-tools.h>


typedef enum {
  _GREY_,
  _COLOR_
} enumColorPlot;



typedef struct local_parameter {

  char *theimage_name;
  char *resimage_name;
  char *track_name;

  char *redimage_name;
  char *greenimage_name;
  char *blueimage_name;

  char *colortable_name;

  enumColorPlot color;
  int value;

  int trace_growing_tracks;
  int keep_old_tracks;

  int dimz;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );





static int _buildColorTable( bal_image *colortable,
			     bufferType type,
			     int value );

static int _plotTrack3DinColor( bal_image *imr,
				bal_image *img,
				bal_image *imb,
				bal_blTrackList *list,
				bal_image *color,
			       int traceold );

static int _plotTrack3DinGrey( bal_image *theim,
			       bal_blTrackList *list,
			       int value,
			       int traceold );




static char *program = NULL;

static char *usage = "%s %s [-track %s]\n\
 [-value|-val %d]\n\
 [-keep-old-tracks | -k] [-ignore-old-tracks | -nk | -i]\n\
 [-growing-tracks | -growing] [-all-tracks | -a | -ng]\n\
 [-table-color|-color-table-tc %s] [-color]\n\
 [-output-red|-red|-r %s] [-output-green|-green|-g %s] [-output-blue|-blur|-b %s]\n\
 [-v] [-help]";

static char *detail = "\
 -v : mode verbose\n\
\n";







int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_image theim;

  int dimz;
  bal_image tmpim;
  int theDim[4];
  int resDim[4];
  int winDim[3] = {-1, -1, -1};
  int theLeftCorner[3] = {0, 0, 0};
  int resLeftCorner[3] = {0, 0, 0};


  bal_image img;
  bal_image imb;

  bal_image *resim = (bal_image*)NULL;
  bal_image colortable;

  bal_blTrackList list;
  bal_blDetectionList *detectionList;
  int n;
  int start_index, end_index;

  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );


  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );
  



  
  /* reading input image
   */
  BAL_InitImage( &theim, NULL, 0, 0, 0, 0, UCHAR );
  if ( BAL_ReadImage( &theim, p.theimage_name, 0 ) != 1 ) {
    _ErrorParse( "unable to read input image\n", 0 );
  }





  /* reading tracks 
   */
  BAL_InitBlTrackList( &list );
  if ( p.track_name == (char*)NULL ) {
    BAL_FreeImage( &theim );
    _ErrorParse( "no track file\n", 0 );
  }

  if ( BAL_ReadBlTrackList( &list, p.track_name ) != 0 ) {
    BAL_FreeImage( &theim );
    _ErrorParse( "unable to read track file\n", 0 );
  }

  if ( list.n <= 0 ) {
    BAL_FreeBlTrackList( &list );
    BAL_FreeImage( &theim );
    _ErrorParse( "empty track file ?\n", 0 );
  }


  /* output image
     evaluate minimum number of images
  */
  BAL_InitImage( &tmpim, NULL, 0, 0, 0, 0, UCHAR );

  if ( p.trace_growing_tracks ) {
    
    start_index = list.data[0].detectionList.data[0].imageindex;
    end_index = list.data[0].detectionList.data[ list.data[0].detectionList.n-1 ].imageindex;
    for ( n=1; n<list.n; n++ ) {
      detectionList = &(list.data[n].detectionList);
      if ( start_index > detectionList->data[0].imageindex )
	start_index = detectionList->data[0].imageindex;
      if ( end_index < detectionList->data[ detectionList->n-1 ].imageindex )
	end_index = detectionList->data[ detectionList->n-1 ].imageindex;
    }
    if ( _verbose_ >= 2 )
      fprintf( stderr, "indexes range in [%d - %d]\n", start_index, end_index );
    
    /* it is assumed that image indexes start at 0
       and that drawing always start from 0 
     */
    dimz = p.dimz;
    if ( dimz < end_index + 1 ) dimz = end_index + 1;

    
    if ( theim.nplanes > 1 ) {

      if ( theim.nplanes < dimz ) {
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( &theim );
	_ErrorParse( "the 3D image does not have a sufficient number of slices\n", 0 );
      }
      resim = &theim;

    }
    else {
      
      if ( BAL_InitAllocImage( &tmpim, (char*)NULL, 
			       theim.ncols, theim.nrows, dimz, 
			       theim.vdim, theim.type ) != 1 ) {
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( &theim );
	_ErrorParse( "unable to allocate output image\n", 0 );
      }
      
      resim = &tmpim;
      
      theDim[0] = theim.ncols;
      theDim[1] = theim.nrows;
      theDim[2] = theim.nplanes;
      theDim[3] = theim.vdim;
      
      theLeftCorner[0] = 0;
      theLeftCorner[1] = 0;
      theLeftCorner[2] = 0;
      
      winDim[0] = theim.ncols;
      winDim[1] = theim.nrows;
      winDim[2] = 1;
      
      resDim[0] = resim->ncols;
      resDim[1] = resim->nrows;
      resDim[2] = resim->nplanes;
      resDim[3] = resim->vdim;
      
      resLeftCorner[0] = 0;
      resLeftCorner[1] = 0;
      resLeftCorner[2] = 0;
    
      for ( n=0; n<dimz; n++ ) {
	resLeftCorner[2] = n;
	ExtractFromBufferAndMelt( theim.data, theDim, resim->data, resDim, 
				  theLeftCorner, resLeftCorner, winDim, theim.type );
      }

      BAL_FreeImage( &theim );
      
    }

  }
  else {
    resim = &theim;
  }





  switch ( p.color ) {
  default :
  case _GREY_ :

    if ( _plotTrack3DinGrey(  resim, &list, p.value, p.keep_old_tracks ) != 1 ) {
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "error when plotting\n", 0 );
    }
    if ( BAL_WriteImage( resim, p.resimage_name ) != 1 ) {
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to write output image\n", 0 );
    }
    break;

  case _COLOR_ :
    
    if ( BAL_InitAllocImage( &img, (char*)NULL, 
                             resim->ncols, resim->nrows, resim->nplanes, 
                             resim->vdim, resim->type ) != 1 ) {
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to allocate green image\n", 0 );
    }
    if ( BAL_CopyImage( resim, &img ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to copy green image\n", 0 );
    }

    if ( BAL_InitAllocImage( &imb, (char*)NULL, 
                             resim->ncols, resim->nrows, resim->nplanes, 
                             resim->vdim, resim->type ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to allocate green image\n", 0 );
    }
    if ( BAL_CopyImage( resim, &imb ) != 1 ) {
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to copy green image\n", 0 );
    }

    /* get a colortable
     */
    BAL_InitImage( &colortable, NULL, 0, 0, 0, 0, UCHAR );
    if ( p.colortable_name != (char*)NULL ) {
      /* color plot 
     */
      if ( BAL_ReadImage( &colortable, p.colortable_name, 0 ) != 1 ) {
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( resim );
	_ErrorParse( "unable to read color table\n", 0 );
      }
      if ( colortable.type != resim->type ) {
	BAL_FreeImage( &colortable );
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( resim );
	_ErrorParse( "color table and input image have different types\n", 0 );
      }
      if ( colortable.nrows != 3 
	   || colortable.nplanes != 1 
	   ||colortable.vdim != 1 ) {
	BAL_FreeImage( &colortable );
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( resim );
	_ErrorParse( "weird dimensions for the color table\n", 0 );
      }
    }
    else {
      if ( _buildColorTable( &colortable, resim->type, p.value ) != 1 ) {
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlTrackList( &list );
	BAL_FreeImage( resim );
	_ErrorParse( "unable to build color table\n", 0 );
      }
    }
      
    if ( _plotTrack3DinColor( resim, &img, &imb, &list, &colortable, p.keep_old_tracks ) != 1 ) {
      BAL_FreeImage( &colortable );
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "error when color plotting\n", 0 );
    }

    BAL_FreeImage( &colortable );

    if ( BAL_WriteImage( &imb, p.blueimage_name ) != 1 ) {
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to write blue output image\n", 0 );
    }
    BAL_FreeImage( &imb );
    if ( BAL_WriteImage( &img, p.greenimage_name ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to write green output image\n", 0 );
    }
    BAL_FreeImage( &img );
    if ( BAL_WriteImage( resim, p.redimage_name ) != 1 ) {
      BAL_FreeBlTrackList( &list );
      BAL_FreeImage( resim );
      _ErrorParse( "unable to write red output image\n", 0 );
    }

    break;
  }



  BAL_FreeBlTrackList( &list );
  BAL_FreeImage( resim );  
    

  return( 1 );
}








static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* file names
       */

      if ( strcmp ( argv[i], "-track" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -track...\n", 0 );
	if ( p->track_name != (char*)NULL ) 
	  _ErrorParse( "parsing -track: input has already been parsed ...\n", 0 );
	p->track_name = argv[i];
      }									
      


      else if ( (strcmp ( argv[i], "-table-color" ) == 0)
		|| (strcmp ( argv[i], "-color-table" ) == 0)
		|| (strcmp ( argv[i], "-tc" ) == 0  && argv[i][3] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -color-table...\n", 0 );
	if ( p->colortable_name != (char*)NULL ) 
	  _ErrorParse( "parsing -color-table: input has already been parsed ...\n", 0 );
	p->colortable_name = argv[i];
	p->color = _COLOR_;
      }									

      else if ( (strcmp ( argv[i], "-output-red" ) == 0)
		|| (strcmp ( argv[i], "-red" ) == 0  && argv[i][4] == '\0') 
		|| (strcmp ( argv[i], "-r" ) == 0  && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -output-red...\n", 0 );
	if ( p->redimage_name != (char*)NULL ) 
	  _ErrorParse( "parsing -output-red: input has already been parsed ...\n", 0 );
	p->redimage_name = argv[i];
	p->color = _COLOR_;
      }									
 
      else if ( (strcmp ( argv[i], "-output-green" ) == 0)
		|| (strcmp ( argv[i], "-green" ) == 0  && argv[i][6] == '\0') 
		|| (strcmp ( argv[i], "-g" ) == 0  && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -output-green...\n", 0 );
	if ( p->greenimage_name != (char*)NULL ) 
	  _ErrorParse( "parsing -output-green: input has already been parsed ...\n", 0 );
	p->greenimage_name = argv[i];
	p->color = _COLOR_;
      }									

      else if ( (strcmp ( argv[i], "-output-blue" ) == 0)
		|| (strcmp ( argv[i], "-blue" ) == 0  && argv[i][5] == '\0') 
		|| (strcmp ( argv[i], "-b" ) == 0  && argv[i][2] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -output-blue...\n", 0 );
	if ( p->blueimage_name != (char*)NULL ) 
	  _ErrorParse( "parsing -output-blue: input has already been parsed ...\n", 0 );
	p->blueimage_name = argv[i];
	p->color = _COLOR_;
      }									
     


      /* tracing options
       */

      else if ( (strcmp (argv[i], "-value" ) == 0 && argv[i][6] == '\0') 
		|| (strcmp (argv[i], "-val" ) == 0 && argv[i][4] == '\0') ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -value %d", 0 );
	status = sscanf( argv[i], "%d", &(p->value) );
	if ( status <= 0 ) _ErrorParse( "parsing -value %d", 0 );
	p->color = _GREY_;
      }

      else if ( strcmp ( argv[i], "-color") == 0 ) {
	p->color = _COLOR_;
      }

      else if ( strcmp ( argv[i], "-keep-old-tracks") == 0
		|| (strcmp ( argv[i], "-k") == 0 && argv[i][2] == '\0') ) {
	p->keep_old_tracks = 1;
      }
      
      else if ( strcmp ( argv[i], "-ignore-old-tracks") == 0
		|| (strcmp ( argv[i], "-nk") == 0 && argv[i][3] == '\0') 
		|| (strcmp ( argv[i], "-i") == 0 && argv[i][2] == '\0') ) {
	p->keep_old_tracks = 0;
      }

      else if ( strcmp ( argv[i], "-growing-tracks") == 0
		|| (strcmp ( argv[i], "-growing") == 0 ) ) {
	p->trace_growing_tracks = 1;
      }
      
      else if ( strcmp ( argv[i], "-all-tracks") == 0
		|| (strcmp ( argv[i], "-ng") == 0 && argv[i][3] == '\0') 
		|| (strcmp ( argv[i], "-a") == 0 && argv[i][2] == '\0') ) {
	p->trace_growing_tracks = 0;
      }

      else if ( strcmp ( argv[i], "-z") == 0 && argv[i][2] == '\0' ) {
	i ++;
        if ( i >= argc)    _ErrorParse( "parsing -z", 0 );
        status = sscanf( argv[i], "%d", &(p->dimz) );
        if ( status <= 0 ) _ErrorParse( "parsing -z", 0 );
      }


      /* misc
       */

      else if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
      }
      
      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }
      
       else if ( (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
	_debug_ = 1;
      }

      /* unknown option
       */
      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 ) {
	p->theimage_name = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->resimage_name = argv[i];
	outputisread = 1;
      }
      else 
	fprintf(stderr,"too many file names: '%s'\n",argv[i]);
    }

  }
  
}





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  if ( str != (char*)NULL )
    (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}





static void _InitParam( local_parameter *p )
{
  p->theimage_name = (char*)NULL;
  p->resimage_name = (char*)NULL;
  p->track_name = (char*)NULL;

  p->redimage_name = (char*)NULL;
  p->blueimage_name = (char*)NULL;
  p->greenimage_name = (char*)NULL;

  p->colortable_name = (char*)NULL;

  p->color = _GREY_;
  p->value = -1;

  p->trace_growing_tracks = 1;
  p->keep_old_tracks = 1;

  p->dimz = -1;;
}










/************************************************************
 *
 *
 *
 ************************************************************/


typedef struct {
  int dx;
  int dy;
  int dz;
} typeOffset;

typedef struct {
  typeOffset *data;
  int n;
  int n_allocated;
} typeOffsetList;

typedef struct {
  typeOffsetList *data;
  int n;
  int n_allocated;
} typeOffsetListList;



#ifdef _UNUSED_
static void _initTypeOffset( typeOffset *o )
{
  o->dx = 0;
  o->dy = 0;
  o->dz = 0;
}
#endif



#ifdef _UNUSED_
static void _initTypeOffsetList( typeOffsetList *l )
{
  l->data = (typeOffset *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}
#endif



#ifdef _UNUSED_
static void _initTypeOffsetListList( typeOffsetListList *l )
{
  l->data = (typeOffsetList *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}
#endif



#ifdef _UNUSED_
static void _freeTypeOffsetList( typeOffsetList *l )
{
  if ( l->data != (typeOffset *)NULL )
    free( l->data );
  _initTypeOffsetList( l );
}
#endif



#ifdef _UNUSED_
static void _freeTypeOffsetListList( typeOffsetListList *l )
{
  int i;
  if ( l->data != (typeOffsetList *)NULL ) {
    for ( i=0; i<l->n_allocated; i++ )
      _freeTypeOffsetList( &(l->data[i]) );
  }
  free( l->data );
  _initTypeOffsetListList( l );
}
#endif



#ifdef _UNUSED_
static int _allocTypeOffsetList( typeOffsetList *l, int n )
{
  char *proc = "_allocTypeOffsetList";
  int i;
  
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: null or negative number\n", proc );
    return( -1 );
  }

  l->data = (typeOffset *)malloc( n * sizeof( typeOffset ) );
  if ( l->data == (typeOffset *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }

  for ( i=0; i<n; i++ )
    _initTypeOffset( &(l->data[i]) );
  l->n = 0;
  l->n_allocated = n;

  return( 1 );
}
#endif



#ifdef _UNUSED_
static int _allocTypeOffsetListList( typeOffsetListList *l, int n )
{
  char *proc = "_allocTypeOffsetListList";
  int i;
  
  if ( n <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: null or negative number\n", proc );
    return( -1 );
  }

  l->data = (typeOffsetList *)malloc( n * sizeof( typeOffsetList ) );
  if ( l->data == (typeOffsetList *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: allocation failed\n", proc );
    return( -1 );
  }

  for ( i=0; i<n; i++ )
    _initTypeOffsetList( &(l->data[i]) );
  l->n = 0;
  l->n_allocated = n;

  return( 1 );
}
#endif





#ifdef _UNUSED_
static int _build2DOffsetListList( typeOffsetListList *l, int rmax )
{
  char *proc = "_buildTypeOffsetListList";
  int *buf = (int*)NULL;
  int rdim = rmax + 1;
  int dim = 2*rdim+1;
  int x, y, r, n;

  if ( rmax <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: null or negative radius\n", proc );
    return( -1 );
  }
    
  buf = (int*)malloc( dim*dim * sizeof(int) );
  if ( buf == (int*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: buffed allocation failed\n", proc );
    return( -1 );
  }
  
  if ( _allocTypeOffsetListList( l, rmax+1 ) != 1 ) {
    free( buf );
    if ( _verbose_ )
      fprintf( stderr, "%s: offset list allocation failed\n", proc );
    return( -1 );
  }
  
  for ( y= -rdim; y<=rdim; y++ )
  for ( x= -rdim; x<=rdim; x++ )
    buf[ (y+rdim)*dim + (x+rdim) ] =
      (int)( sqrt( (double)(y*y + x*x) ) + 0.5 );
  
  
  for ( r=1; r<=rmax; r++ ) {

    for ( n=0, y= -rdim; y<=rdim; y++ )
    for ( x= -rdim; x<=rdim; x++ ) {
      if ( buf[ (y+rdim)*dim + (x+rdim) ] > 0 ) continue;
      if ( buf[ (y+rdim)*dim + (x+rdim+1) ] > 0 ) n++;
      if ( buf[ (y+rdim)*dim + (x+rdim-1) ] > 0 ) n++;
      if ( buf[ (y+rdim+1)*dim + (x+rdim) ] > 0 ) n++;
      if ( buf[ (y+rdim-1)*dim + (x+rdim) ] > 0 ) n++;
    }
    
    if ( _verbose_ )
      fprintf( stderr, "%s: found %d points for radius=%d\n", proc, n, r );
    
    if ( _allocTypeOffsetList( &(l->data[r]), n ) != 1 ) {
      _freeTypeOffsetListList( l );
      free( buf );
      if ( _verbose_ )
	fprintf( stderr, "%s: offset list allocation failed for radius %d\n", proc, r );
      return( -1 );
    }
    
    for ( n=0, y= -rdim; y<=rdim; y++ )
    for ( x= -rdim; x<=rdim; x++ ) {
      if ( buf[ (y+rdim)*dim + (x+rdim) ] > 0 ) continue;
      if ( buf[ (y+rdim)*dim + (x+rdim+1) ] > 0 ) {
	l->data[r].data[n].dx = x+1;
	l->data[r].data[n].dy = y;
	n++;
      }
      if ( buf[ (y+rdim)*dim + (x+rdim-1) ] > 0 ) {
	l->data[r].data[n].dx = x-1;
	l->data[r].data[n].dy = y;
	n++;
      }
      if ( buf[ (y+rdim+1)*dim + (x+rdim) ] > 0 ) {
	l->data[r].data[n].dx = x;
	l->data[r].data[n].dy = y+1;
	n++;
      }
      if ( buf[ (y+rdim-1)*dim + (x+rdim) ] > 0 ) {
	l->data[r].data[n].dx = x;
	l->data[r].data[n].dy = y-1;
	n++;
      }
    }
    l->data[r].n = n;
    
    for ( n=0, y= -rdim; y<=rdim; y++ )
    for ( x= -rdim; x<=rdim; x++ ) {
      if ( buf[ (y+rdim)*dim + (x+rdim) ] > r ) continue;
      buf[ (y+rdim)*dim + (x+rdim) ] = 0;
    }

  }

  free( buf );
  return( 1 );
}
#endif









/************************************************************
 *
 *
 *
 ************************************************************/

static int _buildColorTable( bal_image *colortable,
			     bufferType type,
			     int value )
{
  char *proc = "_buildColorTable";
  int i, v = value;

  switch( type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 ***theColor;

      if ( BAL_InitAllocImage( colortable, (char*)NULL, 15, 3, 1, 1, type ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when allocating image (UCHAR)\n", proc );
	return( -1 );
      }

      if ( v <= 0 || v > 255 ) v = 255;
      theColor = (u8***)(colortable->array);

      if ( _verbose_ >= 3 )
	fprintf( stderr, "%s: build u8 color table with v=%d\n", proc, v );

      i = 0;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;

      theColor[0][0][i] = v;   theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = v;   i++;

      theColor[0][0][i] = v;   theColor[0][1][i] = v/2; theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = v/2; theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v/2; theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v/2; theColor[0][2][i] = v;   i++;
    }
    break;
  case SSHORT :
    {
      s16 ***theColor;

      if ( BAL_InitAllocImage( colortable, (char*)NULL, 15, 3, 1, 1, type ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when allocating image (SSHORT)\n", proc );
	return( -1 );
      }

      if ( v <= 0 || v > 4095 ) v = 4095;
      theColor = (s16***)(colortable->array);

      if ( _verbose_ >= 3 )
	fprintf( stderr, "%s: build s16 color table with v=%d\n", proc, v );

      i=0;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;

      theColor[0][0][i] = v;   theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = v;   i++;

      theColor[0][0][i] = v;   theColor[0][1][i] = v/2; theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v;   theColor[0][1][i] = v/2; theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v;   theColor[0][2][i] = 0;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v;   theColor[0][2][i] = v/2; i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = 0;   theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = 0;   theColor[0][1][i] = v/2; theColor[0][2][i] = v;   i++;
      theColor[0][0][i] = v/2; theColor[0][1][i] = v/2; theColor[0][2][i] = v;   i++;
    }
    break;
  }
  return( 1 );
}









/************************************************************
 *
 *
 *
 ************************************************************/



/*
 * line3d was dervied from DigitalLine.c published as "Digital Line Drawing"
 * by Paul Heckbert from "Graphics Gems", Academic Press, 1990
 * 
 * 3D modifications by Bob Pendleton. The original source code was in the public
 * domain, the author of the 3D version places his modifications in the
 * public domain as well.
 * 
 * line3d uses Bresenham's algorithm to generate the 3 dimensional points on a
 * line from (x1, y1, z1) to (x2, y2, z2)
 * 
 */

/* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))

/* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)

static void _traceline( bal_image *theim,
			int value,
			int pt1x, int pt1y, int pt1z,
			int pt2x, int pt2y, int pt2z )
{
  char *proc = "_traceline";
  int x1, y1, z1;
  int x2, y2, z2;

  int v = value;

  int xd, yd, zd;
  int x, y, z;
  int ax, ay, az;
  int sx, sy, sz;
  int dx, dy, dz;
  
  if ( (pt1x < 0 || pt1x >= theim->ncols)
       || (pt1y < 0 || pt1y >= theim->nrows)
       || (pt1z < 0 || pt1z >= theim->nplanes) ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: first point (%d %d %d) out of image\n",
	       proc, pt1x, pt1y, pt1z );
    return;
  }

  if ( (pt2x < 0 || pt2x >= theim->ncols)
       || (pt2y < 0 || pt2y >= theim->nrows)
       || (pt2z < 0 || pt2z >= theim->nplanes) ) {
    if ( _verbose_ >= 2 )
      fprintf( stderr, "%s: second point (%d %d %d) out of image\n",
	       proc, pt2x, pt2y, pt2z );
    return;
  }

  x1 = pt1x;
  y1 = pt1y;
  z1 = pt1z;

  x2 = pt2x;
  y2 = pt2y;
  z2 = pt2z;

  switch( theim->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return;
  case UCHAR :
    if ( v < 0 || v > 255 ) v = 255;
    break;
  case SSHORT :
    if ( v < 0 || v > 32767 ) v = 4095;
    break;
  case USHORT :
    if ( v < 0 || v > 65535 ) v = 4095;
    break;
  }

  

  dx = x2 - x1;
  dy = y2 - y1;
  dz = z2 - z1;
  
  ax = ABS(dx) << 1;
  ay = ABS(dy) << 1;
  az = ABS(dz) << 1;

  sx = ZSGN(dx);
  sy = ZSGN(dy);
  sz = ZSGN(dz);
  
  x = x1;
  y = y1;
  z = z1;
  
  if (ax >= MAX(ay, az)) {
    /* x dominant */
    yd = ay - (ax >> 1);
    zd = az - (ax >> 1);
    for (;;) {

      switch ( theim->type ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such image type not handled yet\n", proc );
	return;
      case UCHAR :
	((u8***)(theim->array))[z][y][x] = v;
	break;
      case SSHORT :
	((s16***)(theim->array))[z][y][x] = v;
	break;
      case USHORT :
	((u16***)(theim->array))[z][y][x] = v;
	break;
      }

      if (x == x2)
	return;

      if (yd >= 0) {
	y += sy;
	yd -= ax;
      }

      if (zd >= 0) {
	z += sz;
	zd -= ax;
      }

      x += sx;
      yd += ay;
      zd += az;
    }
  }

  else if (ay >= MAX(ax, az)) {
    /* y dominant */
    xd = ax - (ay >> 1);
    zd = az - (ay >> 1);
    for (;;) {

      switch ( theim->type ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such image type not handled yet\n", proc );
	return;
      case UCHAR :
	((u8***)(theim->array))[z][y][x] = v;
	break;
      case SSHORT :
	((s16***)(theim->array))[z][y][x] = v;
	break;
      case USHORT :
	((u16***)(theim->array))[z][y][x] = v;
	break;
      }

      if (y == y2)
	return;

      if (xd >= 0) {
	x += sx;
	xd -= ay;
      }

      if (zd >= 0) {
	z += sz;
	zd -= ay;
      }
      
      y += sy;
      xd += ax;
      zd += az;
    }
  }
  else if (az >= MAX(ax, ay)) {
    /* z dominant */
    xd = ax - (az >> 1);
    yd = ay - (az >> 1);
    for (;;) {
      
      switch ( theim->type ) {
      default :
	if ( _verbose_ )
	  fprintf( stderr, "%s: such image type not handled yet\n", proc );
	return;
      case UCHAR :
	((u8***)(theim->array))[z][y][x] = v;
	break;
      case SSHORT :
	((s16***)(theim->array))[z][y][x] = v;
	break;
      case USHORT :
	((u16***)(theim->array))[z][y][x] = v;
	break;
      }
      
      if (z == z2)
	return;
      
      if (xd >= 0) {
	x += sx;
	xd -= az;
      }
      
      if (yd >= 0) {
	y += sy;
	yd -= az;
      }
      
      z += sz;
      xd += ax;
      yd += ay;
    }
  }
}



static int _plotTrack3DinColor( bal_image *imr,
				bal_image *img,
				bal_image *imb,
				bal_blTrackList *trackList,
				bal_image *colortable,
				int traceold )
{
  char *proc = "_plotTrack3DinColor";
  bal_blDetectionList *detectionList;
  int icolor;
  int n;
  int i;
  int z;
  int px, py, nx, ny;
  int red, green, blue;



 /* all tracks are traced into the same image
   */
  if ( imr->nplanes == 1 ) {

    for ( n=0; n<trackList->n; n++ ) {
      icolor = n % (int)colortable->ncols;
      detectionList = &(trackList->data[n].detectionList);

      if ( detectionList->n == 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: track #%d has %d elements, skip it\n",
		   proc, n, detectionList->n );
	continue;
      }

      switch( colortable->type ) {
      default :
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: such color table type not handled yet\n", proc );
	return( -1 );
      case UCHAR :
	{
	  u8 ***color = (u8***)colortable->array;
	  red   = color[0][0][icolor];
	  green = color[0][1][icolor];
	  blue  = color[0][2][icolor];
	}
	break;
       case SSHORT :
	{
	  s16 ***color = (s16***)colortable->array;
	  red   = color[0][0][icolor];
	  green = color[0][1][icolor];
	  blue  = color[0][2][icolor];
	}
	break;
     }

      px = nx = (int)(detectionList->data[0].voxelcenter.x + 0.5);
      py = ny = (int)(detectionList->data[0].voxelcenter.y + 0.5);
      for ( i=0; i<detectionList->n-1; i++ ) {
	nx = (int)(detectionList->data[i+1].voxelcenter.x + 0.5);
	ny = (int)(detectionList->data[i+1].voxelcenter.y + 0.5);
	_traceline( imr, red, px, py, 0, nx, ny, 0 );
	_traceline( img, green, px, py, 0, nx, ny, 0 );
	_traceline( imb, blue, px, py, 0, nx, ny, 0 );
	px = nx;
	py = ny;
      }
    }

  }


  /* z represents the time, ie the field imageindex of detected spots
     thus a track is traced only from the time it appears
   */
  else {
    for ( z=0; z<imr->nplanes; z++ ) {
      
      for ( n=0; n<trackList->n; n++ ) {
	icolor = n % (int)colortable->ncols;
	detectionList = &(trackList->data[n].detectionList);
	if ( detectionList->n == 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: track #%d has %d elements, skip it\n",
		     proc, n, detectionList->n );
	  continue;
	}
	if ( detectionList->data[0].imageindex > z )
	  continue;
	if ( traceold == 0 && detectionList->data[detectionList->n-1].imageindex < z )
	  continue;

	switch( colortable->type ) {
	default :
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: such color table type not handled yet\n", proc );
	  return( -1 );
	case UCHAR :
	  {
	    u8 ***color = (u8***)colortable->array;
	    red   = color[0][0][icolor];
	    green = color[0][1][icolor];
	    blue  = color[0][2][icolor];
	  }
	  break;
	case SSHORT :
	  {
	    s16 ***color = (s16***)colortable->array;
	    red   = color[0][0][icolor];
	    green = color[0][1][icolor];
	    blue  = color[0][2][icolor];
	  }
	  break;
	}

	if ( _verbose_ >= 3 )
	  fprintf( stderr, "%s: use color (%d %d %d)\n", proc, red, green, blue );

	px = nx = (int)(detectionList->data[0].voxelcenter.x + 0.5);
	py = ny = (int)(detectionList->data[0].voxelcenter.y + 0.5);

	if ( detectionList->data[0].imageindex ==  z ) {
	  _traceline( imr, red, px, py, z, nx, ny, z );
	  _traceline( img, green, px, py, z, nx, ny, z );
	  _traceline( imb, blue, px, py, z, nx, ny, z );
	}
	else {
	  for ( i=0; i<detectionList->n-1 
		  && detectionList->data[i+1].imageindex <= z; i++ ) {
	    nx = (int)(detectionList->data[i+1].voxelcenter.x + 0.5);
	    ny = (int)(detectionList->data[i+1].voxelcenter.y + 0.5);
	    _traceline( imr, red, px, py, z, nx, ny, z );
	    _traceline( img, green, px, py, z, nx, ny, z );
	    _traceline( imb, blue, px, py, z, nx, ny, z );
	    px = nx;
	    py = ny;
	  }
	}
      }
      
    }
	
  }


 
  return( 1 );

}





static int _plotTrack3DinGrey( bal_image *theim,
			       bal_blTrackList *trackList,
			       int value,
			       int traceold )
{
  char *proc = "_plotTrack3DinGrey";
  bal_blDetectionList *detectionList;
  int n;
  int i;
  int z;
  int px, py, nx, ny;



  /* all tracks are traced into the same image
   */
  if ( theim->nplanes == 1 ) {

    for ( n=0; n<trackList->n; n++ ) {
      detectionList = &(trackList->data[n].detectionList);
      if ( detectionList->n == 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: track #%d has %d elements, skip it\n",
		   proc, n, detectionList->n );
	continue;
      }

      px = nx = (int)(detectionList->data[0].voxelcenter.x + 0.5);
      py = ny = (int)(detectionList->data[0].voxelcenter.y + 0.5);
      for ( i=0; i<detectionList->n-1; i++ ) {
	nx = (int)(detectionList->data[i+1].voxelcenter.x + 0.5);
	ny = (int)(detectionList->data[i+1].voxelcenter.y + 0.5);
	_traceline( theim, value, px, py, 0, nx, ny, 0 );
	px = nx;
	py = ny;
      }
    }

  }


  /* z represents the time, ie the field imageindex of detected spots
     thus a track is traced only from the time it appears
   */
  else {
    
    for ( z=0; z<theim->nplanes; z++ ) {
      
      for ( n=0; n<trackList->n; n++ ) {
	detectionList = &(trackList->data[n].detectionList);
	if ( detectionList->n == 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: track #%d has %d elements, skip it\n",
		     proc, n, detectionList->n );
	  continue;
	}
	if ( detectionList->data[0].imageindex > z )
	  continue;
	if ( traceold == 0 && detectionList->data[detectionList->n-1].imageindex < z )
	  continue;

	px = nx = (int)(detectionList->data[0].voxelcenter.x + 0.5);
	py = ny = (int)(detectionList->data[0].voxelcenter.y + 0.5);

	if ( detectionList->data[0].imageindex ==  z ) {
	  _traceline( theim, value, px, py, z, nx, ny, z );
	}
	else {
	  for ( i=0; i<detectionList->n-1 
		  && detectionList->data[i+1].imageindex <= z; i++ ) {
	    nx = (int)(detectionList->data[i+1].voxelcenter.x + 0.5);
	    ny = (int)(detectionList->data[i+1].voxelcenter.y + 0.5);
	    _traceline( theim, value, px, py, z, nx, ny, z );
	    px = nx;
	    py = ny;
	  }
	}
      }
      
    }

  }

  return( 1 );

}
