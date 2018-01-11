/*************************************************************************
 * removeCcOnBorder.c - met a zero les composantes connexes numerotees qui
 *                      "touchent" un bord (en X ou en Y)
 *
 * $Id: removeCcOnBorder.c,v 1.3 1999/07/28 14:04:16 greg Exp $
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
 * May 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  int type;
  int x;
  int y;
  int z;
} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[image-in] [image-out] [-x|-y|-z|-all] \n\
\t  [-v] [-D] [-help]";

static char *detail = "\
'image-in' est une image de labels. Les composantes dont les\n\
labels touchent un des bords specifies sont mises a 0.\n\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image;
  int x, y, z, i, j, n;
  int *v = (int*)NULL;
  
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
  
  /*--- initialisation de l'image resultat ---*/
  if ( VT_CopyName( image->name, par.names.out ) == 0 ) {
    VT_ErrorParse( "unable to copy name into image header", 0 );
  }

  if ( par.x == 0 && par.y == 0 && par.z == 0 )
    par.x = par.y = par.z = 1;

  switch ( image->type ) {

  case UCHAR :
    {
      unsigned char *** theBuf = (unsigned char ***)image->array;
      
      n = 0;
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];

      v = (int*)malloc( (n+1) * sizeof(int) );
      if ( v == (int*)NULL ) {
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse( "unable to allocate auxiliary array\n", 0);
      }

      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	v[ theBuf[z][y][x] ] = theBuf[z][y][x];

      /*
      if ( image->dim.z >= 3 ) {
	for ( y=0; y<image->dim.y ; y++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  v[ (int)theBuf[0][y][x] ] = 0;
	  v[ (int)theBuf[image->dim.z-1][y][x] ] = 0;
	}
      }
      */
      if ( image->dim.y >= 3 ) {
	for ( z=0; z<image->dim.z ; z++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  v[ (int)theBuf[z][0][x] ] = 0;
	  v[ (int)theBuf[z][image->dim.y-1][x] ] = 0;
	}
      }
      if ( image->dim.x >= 3 ) {
	for ( z=0; z<image->dim.z ; z++ )
        for ( y=0; y<image->dim.y ; y++ ) {
	  v[ (int)theBuf[z][y][0] ] = 0;
	  v[ (int)theBuf[z][y][image->dim.x-1] ] = 0;
	}
      }
      
      j = 0;
      for ( i=1; i<=n; i++ )
	if ( v[i] > 0 ) v[i] = ++j;

      fprintf( stderr," %d valid components out of %d\n", j, n );

     for ( z=0; z<image->dim.z ; z++ )
     for ( y=0; y<image->dim.y ; y++ )
     for ( x=0; x<image->dim.x ; x++ ) {
       if ( theBuf[z][y][x] == 0 ) continue;
       theBuf[z][y][x] = (unsigned short int)v[ (int)theBuf[z][y][x] ];
     }
     
     free( v );
    }
    /* end case UCHAR */
    break;

  case USHORT :
    {
      unsigned short int *** theBuf = (unsigned short int ***)image->array;
      
      n = 0;
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];

      v = (int*)malloc( (n+1) * sizeof(int) );
      if ( v == (int*)NULL ) {
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse( "unable to allocate auxiliary array\n", 0);
      }

      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	v[ theBuf[z][y][x] ] = theBuf[z][y][x];

      /*
      if ( image->dim.z >= 3 ) {
	for ( y=0; y<image->dim.y ; y++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  v[ (int)theBuf[0][y][x] ] = 0;
	  v[ (int)theBuf[image->dim.z-1][y][x] ] = 0;
	}
      }
      */
      if ( image->dim.y >= 3 ) {
	for ( z=0; z<image->dim.z ; z++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  v[ (int)theBuf[z][0][x] ] = 0;
	  v[ (int)theBuf[z][image->dim.y-1][x] ] = 0;
	}
      }
      if ( image->dim.x >= 3 ) {
	for ( z=0; z<image->dim.z ; z++ )
        for ( y=0; y<image->dim.y ; y++ ) {
	  v[ (int)theBuf[z][y][0] ] = 0;
	  v[ (int)theBuf[z][y][image->dim.x-1] ] = 0;
	}
      }
      
      j = 0;
      for ( i=1; i<=n; i++ )
	if ( v[i] > 0 ) v[i] = ++j;

      fprintf( stderr," %d valid components out of %d\n", j, n );

     for ( z=0; z<image->dim.z ; z++ )
     for ( y=0; y<image->dim.y ; y++ )
     for ( x=0; x<image->dim.x ; x++ ) {
       if ( theBuf[z][y][x] == 0 ) continue;
       theBuf[z][y][x] = (unsigned short int)v[ (int)theBuf[z][y][x] ];
     }
     
     free( v );
    }
    /* end case USHORT */
    break;

  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "unable to deal with such image type\n", 0);
  }




  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( image ) == -1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
  return( 1 );
}

#if defined(_ANSI_)
static void VT_Parse( int argc, char *argv[], local_par *par )
#else
  static void VT_Parse( argc, argv, par )
  int argc;
char *argv[];
local_par *par;
#endif
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


      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	par->x = 1;
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	par->y = 1;
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	par->z = 1;
      }
      else if ( strcmp ( argv[i], "-all" ) == 0 ) {
	par->x = par->y = par->z = 1;
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
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
}

#if defined(_ANSI_)
static void VT_ErrorParse( char *str, int flag )
#else
  static void VT_ErrorParse( str, flag )
  char *str;
int flag;
#endif
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}

#if defined(_ANSI_)
static void VT_InitParam( local_par *par )
#else
  static void VT_InitParam( par )
  local_par *par;
#endif
{
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;
  par->x = par->y = par->z = 0;
}
