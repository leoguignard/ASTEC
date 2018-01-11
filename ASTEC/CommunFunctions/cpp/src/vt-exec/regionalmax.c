/*************************************************************************
 * regionalmax.c - extraction des maxima regionaux
 *
 * $Id$
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
 * Tue Jul 22 11:00:11 CEST 2008
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>
#include <connexe.h>
#include <regionalmax.h>

#include <vt_common.h>
#include <vt_histo.h>





typedef struct local_par {
  vt_names names;

  char *imdiff_name;

  double height;
  double heightMultiplier;
  
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
 [-difference-image | -diff %s]\n\
 [-h %d] [-hm %lf]\n\
 [-component-wise-processing | -cwp] [-no-component-wise-processing | -ncwp]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
 [-time] [-notime]\n\
 [-inv] [-swap] [-v] [-D] [-help]";

static char *detail = "\
if 'image-in' is equal to '-', we consider stdin\n\
if 'image-out' is not specified, we consider stdout\n\
if both are not specified, we consider stdin and stdout\n\
 -h %d    # maxima height\n\
 -hm %lf  # maxima fraction\n\
            only for UCHAR and USHORT types\n\
\n\
 maxima are searched by dilating MIN( image-in*hm , image-in-h )\n\
 'below' image-in, and by subtracting it from image-in afterwards.\n\
\n\
 [-inv]   # inverse 'image-in'\n\
 [-swap]  # swap bytes of 'image-in' (if encoded on 2 bytes)\n\
 [-v]     # mode verbose\n\
 [-D]     #  mode debug\n";



static char program[STRINGLENGTH];



int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, *imres, imtmp, imbin;
  int tmpIsAllocated = 0;
  u8 *resBuf = NULL;

  int dim[3];
  
  int i, v;
  vt_3m m;
  double t, a=1.0, b=0.0;
  double min, max, sum;
  
  int h;
  double hm;
  
  int r;

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;



  m.min = m.moy = m.max = m.ect = (double)0.0;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  
  
  dim[0] = image->dim.x;
  dim[1] = image->dim.y;
  dim[2] = image->dim.z;
  v = dim[0] * dim[1] * dim[2];
  
  
  /*--- allocations :
    on alloue pour les types autres que UCHAR et USHORT 
    ---*/
  VT_Image( &imtmp );

  switch( image->type ) {
    
  case UCHAR :
  case USHORT :
  case SSHORT :
    imres = image;
    break;
    
  case UINT :
  case FLOAT :
  case DOUBLE :
    VT_InitFromImage( &imtmp, image, par.imdiff_name, USHORT );
    if ( VT_AllocImage( &imtmp ) != 1 ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to allocate auxiliary image\n", 0);
    }
    tmpIsAllocated = 1;
    imres = &imtmp;
    break;
    
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("such image type not handled yet\n", 0);
  }
  
  
  
  /*--- pre-processing ---*/
  
  switch ( image->type ) {
    
  case UCHAR :
  case USHORT :
    break;
    
  case SSHORT :
    {
      s16 *buf = (s16*)imres->buf;
      u16 *auxBuf = (u16*)imres->buf;
      m.min = buf[0];
      for ( i=0; i<v; i++ ) {
	if ( m.min > buf[i] ) m.min = buf[i];
	auxBuf[i] = buf[i] + 32768; 
      }
    }
    break;
    
  case UINT :
  case FLOAT :
  case DOUBLE :
    {
      u16 *theBuf = (u16*)imres->buf;
      if ( VT_3m( image, (vt_image*)NULL, &m ) == -1 ) {
	if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute minimum, mean and maximum\n", 0);
      }
      if ( m.max == m.min ) {
	if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("uniform image\n", 0);
      }
      a = 65535.0 / ( m.max - m.min );
      b = - m.min;
      
      fprintf( stderr, " input image statistics: min=%lf - mean=%lf - max=%lf\n", m.min, m.moy, m.max );
      fprintf( stderr, " input image normalized by (int)((I + %lf) * %lf + 0.5)\n", b, a );
      
      
      switch ( image->type ) {
      case UINT :
	{
	  u32 *buf = (u32*)image->buf;
	  for ( i=0; i<v; i++ ) {
	    t = ((double)buf[i] + b) * a;
	    if ( t < 0 ) theBuf[i] = 0;
	    else if ( t > 65535 ) theBuf[i] = 65535;
	    else theBuf[i] = (int)( t + 0.5 );
	  }
	}
	break;
      case FLOAT :
	{
	  r32 *buf = (r32*)image->buf;
	  for ( i=0; i<v; i++ ) {
	    t = ((double)buf[i] + b) * a;
	    if ( t < 0 ) theBuf[i] = 0;
	    else if ( t > 65535 ) theBuf[i] = 65535;
	    else theBuf[i] = (int)( t + 0.5 );
	  }
	}
	break;
      case DOUBLE :
	{
	  r64 *buf = (r64*)image->buf;
	  for ( i=0; i<v; i++ ) {
	    t = ((double)buf[i] + b) * a;
	    if ( t < 0 ) theBuf[i] = 0;
	    else if ( t > 65535 ) theBuf[i] = 65535;
	    else theBuf[i] = (int)( t + 0.5 );
	  }
	}
	break;
      default :
	VT_ErrorParse("such image type not handled yet (but should here)\n", 0);
      }
      
      min = max = theBuf[0];
      sum = 0;
      for ( i=0; i<v; i++ ) {
	if ( min > theBuf[i] ) min = theBuf[i];
	if ( max < theBuf[i] ) max = theBuf[i];
	sum += theBuf[i];
      }
      
      fprintf( stderr, " auxi. image statistics: min=%lf - mean=%lf - max=%lf\n", min, sum/v, max );
    }
    break;
    
  default :
    if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("such image type not handled yet\n", 0);
  }
  
  

  /*--- parameters ---*/
  
  switch ( image->type ) {
  case UCHAR :
  case USHORT :
  case SSHORT :
    h = (int)( par.height + 0.5 );
    if ( h < 1 ) h = 1;
    break;
  case UINT :
  case FLOAT :
  case DOUBLE :
    h = par.height * a;
    break;
  default :
    h = 1;
  }
  
  switch ( image->type ) {
  case UCHAR :
  case USHORT :
    hm = par.heightMultiplier;
    if ( hm <= 0.0 || hm > 1.0 ) hm = 1.0;
    break;
  case SSHORT :
  case UINT :
  case FLOAT :
  case DOUBLE :
    if ( m.min < 0.0 ) {
      hm = 1.0;
      fprintf( stderr, "%s: height multiplier set to 1.0\n", program );
    }
    else
      /* hm = par.heightMultiplier * a;
	 GM: no idea why I've multiplied by 'a' here :(
	 Wed Jun 15 17:11:14 CEST 2011
      */
      hm = par.heightMultiplier;
    break;
  default :
    hm = 1.0;
  }
  

  
  
  
  /*--- processing ---*/
  
  switch ( image->type ) {
  case UCHAR :
  case USHORT :
    r = regionalmax( imres->buf, imres->buf, image->type, dim, h, hm );
    break;
  case SSHORT :
    r = regionalmax( imres->buf, imres->buf, USHORT, dim, h, hm );
    break;
  case UINT :
  case FLOAT :
  case DOUBLE :
    r = regionalmax( imres->buf, imres->buf, USHORT, dim, h, hm );
    break;
  default :
    if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("such image type not handled yet\n", 0);
  }
  
  
  /*--- difference ---*/
  if ( par.imdiff_name != (char*)NULL ) {
    if ( VT_CopyName( imres->name, par.imdiff_name ) != 1 ) {
      if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("error when copying name\n", 0);
    }
    if ( VT_WriteInrimage( imres ) == -1 ) {
      if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to write difference image\n", 0);
    }
  }
  
  /*--- producing the result ---*/
  VT_Image( &imbin );
  VT_InitFromImage( &imbin, image, par.names.out, UCHAR );
  if ( VT_AllocImage( &imbin ) != 1 ) {
    if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  resBuf = (u8*)imbin.buf;
  
  switch ( image->type ) {
  case UCHAR :
    {
      u8 *buf = (u8*)imres->buf;
      for (i=0; i<v; i++ ) 
	resBuf[i] = ( buf[i] > 0 ) ? 255 : 0 ;
    }
    break;
  case USHORT :
  case SSHORT :
  case UINT :
  case FLOAT :
  case DOUBLE :
    {
      u16 *buf = (u16*)imres->buf;
      for (i=0; i<v; i++ ) 
	resBuf[i] = ( buf[i] > 0 ) ? 255 : 0 ;
    }
    break;
  default :
    VT_FreeImage( &imbin );
    if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("such image type not handled yet\n", 0);
  }
  
  
  
  /*--- liberations memoires ---*/
  if ( tmpIsAllocated ) VT_FreeImage( &imtmp );
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imbin ) == -1 ) {
    VT_FreeImage( &imbin );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( &imbin );

  time_exit = _GetTime();
  clock_exit = _GetClock();
  
  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( program ), time_exit - time_init );
    fprintf( stderr, "\t     elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t     ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }
  return( 1 );
}





static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, nb, status;
  int tmp;
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
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	_VT_VERBOSE_ = 1;
	incrementVerboseInRegionalMax();
	if ( 0 ) incrementVerboseInConnexe();
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0' ) {
	_VT_DEBUG_ = 1;
      }
      
      else if ( strcmp ( argv[i], "-difference-image" ) == 0 
		|| strcmp ( argv[i], "-diff" ) == 0 ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-difference-image...\n", 0 );
	par->imdiff_name = argv[i];
      }

      /*--- seuil ---*/

      else if ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -h...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->height) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -h...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-hm" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -hm...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->heightMultiplier) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -hm...\n", 0 );
      }
      
      /*--- ---*/
      
      else if ( strcmp ( argv[i], "-component-wise-processing" ) == 0 
		|| strcmp ( argv[i], "-cwp" ) == 0  ) {
	allowComponentWiseProcessingInRegionalMax();
      }
      else if ( strcmp ( argv[i], "-no-component-wise-processing" ) == 0 
		|| strcmp ( argv[i], "-ncwp" ) == 0  ) {
	disallowComponentWiseProcessingInRegionalMax();
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

      /*--- unknown option ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
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
  VT_Names( &(par->names) );
  par->imdiff_name = (char*)NULL;
  par->height = 1.0;
  par->heightMultiplier = 1.0;
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
