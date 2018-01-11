/*************************************************************************
 * removeCcOnBorder.c - met a zero les composantes connexes numerotees qui
 *                      "touchent" un bord (en X ou en Y)
 *
 * $Id: relabel.c,v 1.2 2000/03/28 07:28:05 greg Exp $
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
 * - Mon Mar 27 21:32:00 MEST 2000, Gregoire Malandain
 *   Les labels en-dehors d'un cercle sont enleves.
 *
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  
} local_par;


typedef struct {
  int x;
  int y;
  int z;
  int v;
} typePoint;

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

static char *usage = "[image-in] [image-out]\n\
\t [-dist im-dist] [-v][-D]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, *imdist;
  int x, y, z, i, j, n;
  typePoint *theMax = NULL;
  double d, r;

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
  


  /* get max
   */
  n = 0;
  switch ( image->type ) {
  case UCHAR :
    {
      u8 *** theBuf = (u8 ***)image->array;
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];
    }
    break;
  case USHORT :
    {
      unsigned short int *** theBuf = (unsigned short int ***)image->array;
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];
    }
    break;
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "unable to deal with such image type\n", 0);
  }


  /* allocation du tableau
   */
  theMax = (typePoint*)malloc( (n+1) * sizeof(typePoint) );
  if ( theMax == (int*)NULL ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "unable to allocate auxiliary array\n", 0);
  }
  for (x=0; x<=n; x++ ) {
    theMax[x].x = -1;
    theMax[x].y = -1;
    theMax[x].z = -1;
    theMax[x].v = 0;
  }
  


  /* capture du maximum des distances
   */
  imDist = _VT_Inrimage( par.names.ext );
  if ( imDist == (vt_image*)NULL ) {
    VT_Free( &theMax );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to read input image\n", 0);
  }
  if ( imDist->type != USHORT ) {
    VT_FreeImage( imDist );
    VT_Free( (void**)&imDist );
    VT_Free( &theMax );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to deal with such distance type\n", 0);
  }
  


  switch( imDist->type ) {
  default :
    VT_FreeImage( imDist );
    VT_Free( (void**)&imDist );
    VT_Free( &theMax );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to deal with such distance type\n", 0);
    
  case USHORT :
    unsigned short int *** theDist = (unsigned short int ***)imDist->array;

    switch ( image->type ) {
     
    default :
      VT_FreeImage( imDist );
      VT_Free( (void**)&imDist );
      VT_Free( &theMax );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse( "unable to deal with such image type\n", 0);

    case USHORT :
      {
	unsigned short int *** theBuf = (unsigned short int ***)image->array;
	
	for ( z=0; z<image->dim.z ; z++ )
	for ( y=0; y<image->dim.y ; y++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  if ( theBuf[z][y][x] > 0 && theDist[z][y][x] > 0 ) {
	    if ( theMax[ theBuf[z][y][x] ].v < theDist[z][y][x] ) {
	      theMax[ theBuf[z][y][x] ].v = theDist[z][y][x];
	      theMax[ theBuf[z][y][x] ].x = x;
	      theMax[ theBuf[z][y][x] ].y = y;
	      theMax[ theBuf[z][y][x] ].z = z;
	    }
	  }
	}
      }
      break;
    case UCHAR :
      {
	u8 *** theBuf = (u8 ***)image->array;
      
	
	for ( z=0; z<image->dim.z ; z++ )
	for ( y=0; y<image->dim.y ; y++ )
	for ( x=0; x<image->dim.x ; x++ ) {
	  if ( theBuf[z][y][x] > 0 && theDist[z][y][x] > 0 ) {
	    if ( theMax[ theBuf[z][y][x] ].v < theDist[z][y][x] ) {
	      theMax[ theBuf[z][y][x] ].v = theDist[z][y][x];
	      theMax[ theBuf[z][y][x] ].x = x;
	      theMax[ theBuf[z][y][x] ].y = y;
	      theMax[ theBuf[z][y][x] ].z = z;
	    }
	  }
	}
      }
      break;
    }
  }
  

  /* faire le mask a parti de image.... */

  
  switch ( image->type ) {
     
  default :
    VT_FreeImage( imDist );
    VT_Free( (void**)&imDist );
    VT_Free( &theMax );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse( "unable to deal with such image type\n", 0);

  case USHORT :
    {
      unsigned short int *** theBuf = (unsigned short int ***)image->array;
	
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) {
	theBuf[z][y][x] = 0;
      }
      for (x=1; x<=n; x++ ) {
	if ( theMax[x].x >= 0 ) {
	  theBuf[ theMax[x].z ][ theMax[x].y ][ theMax[x].x ] = theMax[x].v;
	}
      }
    }
    break;
  case UCHAR :
    {
      u8 *** theBuf = (u8 ***)image->array;
      
      for ( z=0; z<image->dim.z ; z++ )
      for ( y=0; y<image->dim.y ; y++ )
      for ( x=0; x<image->dim.x ; x++ ) {
	theBuf[z][y][x] = 0;
      }
      for (x=1; x<=n; x++ ) {
	if ( theMax[x].x >= 0 ) {
	  theBuf[ theMax[x].z ][ theMax[x].y ][ theMax[x].x ] = theMax[x].v;
	}
      }
    }
    break;
  }


  /* traitement
   */


  /*--- initialisation de l'image resultat ---*/
  if ( VT_CopyName( image->name, par.names.out ) == 0 ) {
    VT_ErrorParse( "unable to copy name into image header", 0 );
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
      
      else if ( strcmp ( argv[i], "-dist" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -dist...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
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

  par->removeLabelsOutsideXYCircle = 0;
  par->xycenter[0] = par->xycenter[1] = 0.0;
  par->xyradius = 0.0;

  par->removeLabelsOnXBorder = 0;
  par->removeLabelsOnYBorder = 0;
  par->removeLabelsOnZBorder = 0;
}
