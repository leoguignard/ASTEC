/*************************************************************************
 * blPlotDetectionResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 17 jan 2014 10:19:36 CET
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

#include <bal-image.h>
#include <bal-biolib-tools.h>


typedef enum {
  _GREY_,
  _COLOR_,
  _RED_
} enumColorPlot;



typedef struct local_parameter {

  char *theimage_name;
  char *resimage_name;
  char *detection_name;

  char *redimage_name;
  char *greenimage_name;
  char *blueimage_name;

  char *colortable_name;

  enumColorPlot color;
  int value;

  int inversey;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static int _buildRedColorTable( bal_image *colortable,
				bufferType type,
				int value );

static int _buildColorTable( bal_image *colortable,
			     bufferType type,
			     int value );

static int _plotDetection3DinColor( bal_image *imr,
				    bal_image *img,
				    bal_image *imb,
				    bal_blDetectionList *list,
				    bal_image *color,
				    int inversey );

static int _plotDetection3DinGrey( bal_image *theim,
				   bal_blDetectionList *list,
				   int value,
				   int inversey );




static char *program = NULL;

static char *usage = "%s %s [-detection %s]\n\
 [-value|-val %d]\n\
 [-table-color|-color-table-tc %s] [-color] [-red-color] \n\
 [-output-red|-red|-r %s] [-output-green|-green|-g %s] [-output-blue|-blur|-b %s]\n\
 [-inverse-y | -sx] \n\
 [-v] [-help]";

static char *detail = "\
 [-inverse-y | -sx]  # inverse l'axe des y  =  symetrie / x\n\
 -v : mode verbose\n\
\n";







int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_image theim;
  bal_image img;
  bal_image imb;
  bal_image colortable;
  bal_blDetectionList list;



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



  /* reading detections 
   */
  BAL_InitBlDetectionList( &list );
  if ( p.detection_name == (char*)NULL ) {
    BAL_FreeImage( &theim );
    _ErrorParse( "no detection file\n", 0 );
  }

  if ( BAL_ReadBlDetectionList( &list, p.detection_name ) != 0 ) {
    BAL_FreeImage( &theim );
    _ErrorParse( "unable to read detection file\n", 0 );
  }


  switch ( p.color ) {
  default :
  case _GREY_ :

    if ( _plotDetection3DinGrey(  &theim, &list, p.value, p.inversey ) != 1 ) {
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "error when plotting\n", 0 );
    }
    if ( BAL_WriteImage( &theim, p.resimage_name ) != 1 ) {
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to write output image\n", 0 );
    }
    break;

  case _RED_ :
  case _COLOR_ :
    
    if ( BAL_InitAllocImage( &img, (char*)NULL, 
                             theim.ncols, theim.nrows, theim.nplanes, 
                             theim.vdim, theim.type ) != 1 ) {
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to allocate green image\n", 0 );
    }
    if ( BAL_CopyImage( &theim, &img ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to copy green image\n", 0 );
    }

    if ( BAL_InitAllocImage( &imb, (char*)NULL, 
                             theim.ncols, theim.nrows, theim.nplanes, 
                             theim.vdim, theim.type ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to allocate green image\n", 0 );
    }
    if ( BAL_CopyImage( &theim, &imb ) != 1 ) {
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
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
	BAL_FreeBlDetectionList( &list );
	BAL_FreeImage( &theim );
	_ErrorParse( "unable to read color table\n", 0 );
      }
      if ( colortable.type != theim.type ) {
	BAL_FreeImage( &colortable );
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlDetectionList( &list );
	BAL_FreeImage( &theim );
	_ErrorParse( "color table and input image have different types\n", 0 );
      }
      if ( colortable.nrows != 3 
	   || colortable.nplanes != 1 
	   ||colortable.vdim != 1 ) {
	BAL_FreeImage( &colortable );
	BAL_FreeImage( &imb );
	BAL_FreeImage( &img );
	BAL_FreeBlDetectionList( &list );
	BAL_FreeImage( &theim );
	_ErrorParse( "weird dimensions for the color table\n", 0 );
      }
    }
    else {
      switch( p.color ) {
      default : 
      case _RED_ :
	if ( _buildRedColorTable( &colortable, theim.type, p.value ) != 1 ) {
	  BAL_FreeImage( &imb );
	  BAL_FreeImage( &img );
	  BAL_FreeBlDetectionList( &list );
	  BAL_FreeImage( &theim );
	  _ErrorParse( "unable to build redcolor table\n", 0 );
	}
	break;
      case _COLOR_ :
	if ( _buildColorTable( &colortable, theim.type, p.value ) != 1 ) {
	  BAL_FreeImage( &imb );
	  BAL_FreeImage( &img );
	  BAL_FreeBlDetectionList( &list );
	  BAL_FreeImage( &theim );
	  _ErrorParse( "unable to build color table\n", 0 );
	}
	break;
      }
    }
      
    if ( _plotDetection3DinColor( &theim, &img, &imb, &list, &colortable, p.inversey ) != 1 ) {
      BAL_FreeImage( &colortable );
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "error when color plotting\n", 0 );
    }

    BAL_FreeImage( &colortable );

    if ( BAL_WriteImage( &imb, p.blueimage_name ) != 1 ) {
      BAL_FreeImage( &imb );
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to write blue output image\n", 0 );
    }
    BAL_FreeImage( &imb );
    if ( BAL_WriteImage( &img, p.greenimage_name ) != 1 ) {
      BAL_FreeImage( &img );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to write green output image\n", 0 );
    }
    BAL_FreeImage( &img );
    if ( BAL_WriteImage( &theim, p.redimage_name ) != 1 ) {
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to write red output image\n", 0 );
    }

    break;
  }


  if ( p.colortable_name != (char*)NULL ) {
    /* color plot 
     */
    BAL_InitImage( &colortable, NULL, 0, 0, 0, 0, UCHAR );
    if ( BAL_ReadImage( &colortable, p.colortable_name, 0 ) != 1 ) {
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "unable to read color table\n", 0 );
    }
    if ( colortable.type != theim.type ) {
      BAL_FreeImage( &colortable );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "color table and input image have different types\n", 0 );
    }
    if ( colortable.nrows != 3 
	 || colortable.nplanes != 1 
	 ||colortable.vdim != 1 ) {
      BAL_FreeImage( &colortable );
      BAL_FreeBlDetectionList( &list );
      BAL_FreeImage( &theim );
      _ErrorParse( "weird dimensions for the color table\n", 0 );
    }

    

    /* ... */
   




  }
  else {
    /* grey-level plot
     */
    

  }

  BAL_FreeBlDetectionList( &list );
  BAL_FreeImage( &theim );  
    

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

      if ( strcmp ( argv[i], "-detection" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -detection...\n", 0 );
	if ( p->detection_name != (char*)NULL ) 
	  _ErrorParse( "parsing -detection: input has already been parsed ...\n", 0 );
	p->detection_name = argv[i];
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
     


      /* some options
       */

      else if ( (strcmp (argv[i], "-value" ) == 0 && argv[i][6] == '\0') 
		|| (strcmp (argv[i], "-val" ) == 0 && argv[i][4] == '\0') ) {
	i ++;
	if ( i >= argc)    _ErrorParse( "parsing -value %d", 0 );
	status = sscanf( argv[i], "%d", &(p->value) );
	if ( status <= 0 ) _ErrorParse( "parsing -value %d", 0 );
	p->color = _GREY_;
      }

      else if ( strcmp ( argv[i], "-color") == 0 && argv[i][6] == '\0' ) {
	p->color = _COLOR_;
      }

      else if ( (strcmp ( argv[i], "-red-color") == 0 && argv[i][10] == '\0')
		|| (strcmp ( argv[i], "-color-red") == 0 && argv[i][10] == '\0') ) {
	p->color = _RED_;
      }

      else if ( (strcmp ( argv[i], "-inverse-y") == 0)
		|| (strcmp ( argv[i], "-sx") == 0 && argv[i][3] == '\0') ) {
	p->inversey = 1;
      }



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
  p->detection_name = (char*)NULL;

  p->redimage_name = (char*)NULL;
  p->blueimage_name = (char*)NULL;
  p->greenimage_name = (char*)NULL;

  p->colortable_name = (char*)NULL;

  p->color = _GREY_;
  p->value = -1;

  p->inversey = 0;
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



static void _initTypeOffset( typeOffset *o )
{
  o->dx = 0;
  o->dy = 0;
  o->dz = 0;
}

static void _initTypeOffsetList( typeOffsetList *l )
{
  l->data = (typeOffset *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}

static void _initTypeOffsetListList( typeOffsetListList *l )
{
  l->data = (typeOffsetList *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}

static void _freeTypeOffsetList( typeOffsetList *l )
{
  if ( l->data != (typeOffset *)NULL )
    free( l->data );
  _initTypeOffsetList( l );
}

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
    
    if ( _verbose_ >= 2 )
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










/************************************************************
 *
 *
 *
 ************************************************************/


static int _buildRedColorTable( bal_image *colortable,
				bufferType type,
				int value )
{
  char *proc = "_buildRedColorTable";
  int i, v = value;

  switch( type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    {
      u8 ***theColor;
      if ( BAL_InitAllocImage( colortable, (char*)NULL, 1, 3, 1, 1, type ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when allocating image (UCHAR)\n", proc );
	return( -1 );
      }

      if ( v <= 0 || v > 255 ) v = 255;
      theColor = (u8***)(colortable->array);

      i = 0;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = 0;   i++;
    }
    break;
  case SSHORT :
    {
      s16 ***theColor;
      if ( BAL_InitAllocImage( colortable, (char*)NULL, 1, 3, 1, 1, type ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when allocating image (SSHORT)\n", proc );
	return( -1 );
      }

      if ( v <= 0 || v > 4095 ) v = 4095;
      theColor = (s16***)(colortable->array);

      i = 0;
      theColor[0][0][i] = v;   theColor[0][1][i] = 0;   theColor[0][2][i] = 0;   i++;
    }
    break;
  }
  return( 1 );
}





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
  }
  return( 1 );
}









/************************************************************
 *
 *
 *
 ************************************************************/





static int _plotDetection3DinColor( bal_image *imr,
				    bal_image *img,
				    bal_image *imb,
				    bal_blDetectionList *list,
				    bal_image *colortable,
				    int inversey )
{
  char *proc = "_plotDetection3DinColor";
  typeOffsetListList offsetList;
  int rmax;
  int n;
  int ix, iy, iz;
  int i;
  int r;


  for ( rmax=0, n=0; n<list->n; n++ ) {
    if ( rmax < list->data[n].voxelhalfradius1 ) rmax = list->data[n].voxelhalfradius1;
    if ( rmax < list->data[n].voxelhalfradius2 ) rmax = list->data[n].voxelhalfradius2;
  }

  _initTypeOffsetListList( &offsetList );
  if ( _build2DOffsetListList( &offsetList, rmax ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build offset list\n", proc );
    return( -1 ); 
  }
  

  switch( imr->type ) {
  default :
    _freeTypeOffsetListList( &offsetList  );
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 ); 
  case UCHAR :
    {
      u8 ***rdata = (u8***)imr->array;
      u8 ***gdata = (u8***)img->array;
      u8 ***bdata = (u8***)imb->array;
      u8 ***color = (u8***)colortable->array;

      for ( n=0; n<list->n; n++ ) {
	ix = (int)( list->data[n].voxelcenter.x + 0.5 );
	iy = (int)( list->data[n].voxelcenter.y + 0.5 );
	iz = (int)( list->data[n].voxelcenter.z + 0.5 );
	if ( inversey )
	  iy = imr->nrows - 1 - iy;
	if ( iz < 0 || iz >= imr->nplanes ) continue;

	r = list->data[n].voxelhalfradius1;
	if ( r < list->data[n].voxelhalfradius2 ) r = list->data[n].voxelhalfradius2;

	for ( i=0; i<offsetList.data[r].n; i++ ) {
	  if ( ix+offsetList.data[r].data[i].dx < 0 
	       || ix+offsetList.data[r].data[i].dx >= imr->ncols ) continue;
	  if ( iy+offsetList.data[r].data[i].dy < 0 
	       || iy+offsetList.data[r].data[i].dy >= imr->nrows ) continue;
	   rdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][0][n % (int)colortable->ncols];
	   gdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][1][n % (int)colortable->ncols];
	   bdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][2][n % (int)colortable->ncols];
	}
      }
    }
    break;
   case SSHORT :
    {
      s16 ***rdata = (s16***)imr->array;
      s16 ***gdata = (s16***)img->array;
      s16 ***bdata = (s16***)imb->array;
      s16 ***color = (s16***)colortable->array;

      for ( n=0; n<list->n; n++ ) {
	ix = (int)( list->data[n].voxelcenter.x + 0.5 );
	iy = (int)( list->data[n].voxelcenter.y + 0.5 );
	iz = (int)( list->data[n].voxelcenter.z + 0.5 );
	if ( inversey )
	  iy = imr->nrows - 1 - iy;
	if ( iz < 0 || iz >= imr->nplanes ) continue;

	r = list->data[n].voxelhalfradius1;
	if ( r < list->data[n].voxelhalfradius2 ) r = list->data[n].voxelhalfradius2;

	for ( i=0; i<offsetList.data[r].n; i++ ) {
	  if ( ix+offsetList.data[r].data[i].dx < 0 
	       || ix+offsetList.data[r].data[i].dx >= imr->ncols ) continue;
	  if ( iy+offsetList.data[r].data[i].dy < 0 
	       || iy+offsetList.data[r].data[i].dy >= imr->nrows ) continue;
	   rdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][0][n % (int)colortable->ncols];
	   gdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][1][n % (int)colortable->ncols];
	   bdata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] 
	     = color[0][2][n % (int)colortable->ncols];
	}
      }
    }
    break;

  }

  _freeTypeOffsetListList( &offsetList  );
  return( 1 );

}





static int _plotDetection3DinGrey( bal_image *theim,
				   bal_blDetectionList *list,
				   int value,
				   int inversey )
{
  char *proc = "_plotDetection3DinGrey";
  typeOffsetListList offsetList;
  int rmax;

  int n;
  int ix, iy, iz;
  int i;
  int r;

  int v = value;

  for ( rmax=0, n=0; n<list->n; n++ ) {
    if ( rmax < list->data[n].voxelhalfradius1 ) rmax = list->data[n].voxelhalfradius1;
    if ( rmax < list->data[n].voxelhalfradius2 ) rmax = list->data[n].voxelhalfradius2;
  }

  _initTypeOffsetListList( &offsetList );
  if ( _build2DOffsetListList( &offsetList, rmax ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to build offset list\n", proc );
    return( -1 ); 
  }
  

  switch( theim->type ) {
  default :
    _freeTypeOffsetListList( &offsetList  );
    if ( _verbose_ )
      fprintf( stderr, "%s: such image type not handled yet\n", proc );
    return( -1 ); 
  case UCHAR :
    {
      u8 ***thedata = (u8***)theim->array;

      if ( v < 0 ) v = 255;
      
      for ( n=0; n<list->n; n++ ) {
	ix = (int)( list->data[n].voxelcenter.x + 0.5 );
	iy = (int)( list->data[n].voxelcenter.y + 0.5 );
	iz = (int)( list->data[n].voxelcenter.z + 0.5 );
	if ( inversey )
	  iy = theim->nrows - 1 - iy;
	if ( iz < 0 || iz >= theim->nplanes ) continue;

	r = list->data[n].voxelhalfradius1;
	if ( r < list->data[n].voxelhalfradius2 ) r = list->data[n].voxelhalfradius2;

	for ( i=0; i<offsetList.data[r].n; i++ ) {
	  if ( ix+offsetList.data[r].data[i].dx < 0 
	       || ix+offsetList.data[r].data[i].dx >= theim->ncols ) continue;
	  if ( iy+offsetList.data[r].data[i].dy < 0 
	       || iy+offsetList.data[r].data[i].dy >= theim->nrows ) continue;
	   thedata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] = v;
	}
      }
    }
    break;
  case SSHORT :
    {
      s16 ***thedata = (s16***)theim->array;

      if ( v < 0 ) v = 4095;

      for ( n=0; n<list->n; n++ ) {
	ix = (int)( list->data[n].voxelcenter.x + 0.5 );
	iy = (int)( list->data[n].voxelcenter.y + 0.5 );
	iz = (int)( list->data[n].voxelcenter.z + 0.5 );
	if ( inversey )
	  iy = theim->nrows - 1 - iy;
	if ( iz < 0 || iz >= theim->nplanes ) continue;

	r = list->data[n].voxelhalfradius1;
	if ( r < list->data[n].voxelhalfradius2 ) r = list->data[n].voxelhalfradius2;

	for ( i=0; i<offsetList.data[r].n; i++ ) {
	  if ( ix+offsetList.data[r].data[i].dx < 0 
	       || ix+offsetList.data[r].data[i].dx >= theim->ncols ) continue;
	  if ( iy+offsetList.data[r].data[i].dy < 0 
	       || iy+offsetList.data[r].data[i].dy >= theim->nrows ) continue;
	   thedata[iz][iy+offsetList.data[r].data[i].dy][ix+offsetList.data[r].data[i].dx] = v;
	}
      }
    }
    break;
  }

  _freeTypeOffsetListList( &offsetList  );
  return( 1 );

}
