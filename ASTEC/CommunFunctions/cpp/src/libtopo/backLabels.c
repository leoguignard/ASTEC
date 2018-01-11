/*************************************************************************
 * backLabels.c -
 *
 * $Id: backLabels.c,v 1.1 2000/07/26 12:23:43 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Jul 26 11:16:21 MET DST 2000
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  char xname[STRINGLENGTH];
  char yname[STRINGLENGTH];
  char zname[STRINGLENGTH];
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  vt_image *imx = (vt_image*)NULL;
  vt_image *imy = (vt_image*)NULL;
  vt_image *imz = (vt_image*)NULL;
  int x, y, z;
  
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
  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, image->type );

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  

  if ( par.xname[0] != '\0' ) {
    imx = _VT_Inrimage( par.xname );
    if ( imx == (vt_image*)NULL ) 
      VT_ErrorParse("unable to read x image\n", 0);
  }
  if ( par.yname[0] != '\0' ) {
    imy = _VT_Inrimage( par.yname );
    if ( imy == (vt_image*)NULL ) 
      VT_ErrorParse("unable to read y image\n", 0);
  }
  if ( par.zname[0] != '\0' ) {
    imz = _VT_Inrimage( par.zname );
    if ( imz == (vt_image*)NULL ) 
      VT_ErrorParse("unable to read z image\n", 0);
  }

  
  switch ( image->type ) {
  case UCHAR :
    {
      unsigned char ***theBuf = (unsigned char ***)image->array;
      unsigned char ***resBuf = (unsigned char ***)imres.array;

      for ( z=0; z<image->dim.z; z++ )
      for ( y=0; y<image->dim.y; y++ )
      for ( x=0; x<image->dim.x; x++ )
	resBuf[z][y][x] = theBuf[z][y][x];

      switch ( imx->type ) {
      case SSHORT :
	{
	  s16 ***tX = (s16***)imx->array;
	  s16 ***tY = (s16***)imy->array;
	  s16 ***tZ = (s16***)NULL;
	  
	  if ( imz == (vt_image*)NULL ) {
	    for ( z=0; z<image->dim.z; z++ )
	    for ( y=0; y<image->dim.y; y++ )
	    for ( x=0; x<image->dim.x; x++ ) {
	      if ( theBuf[z][y][x] == 0 ) continue;
	      if ( resBuf[z][y+tY[z][y][x]][x+tX[z][y][x]] == 0 )
		resBuf[z][y+tY[z][y][x]][x+tX[z][y][x]] = theBuf[z][y][x];
	    }
	  } else {
	    tZ = (s16***)imz->array;
	    for ( z=0; z<image->dim.z; z++ )
	    for ( y=0; y<image->dim.y; y++ )
	    for ( x=0; x<image->dim.x; x++ ) {
	      if ( theBuf[z][y][x] == 0 ) continue;
	      if ( resBuf[z+tZ[z][y][x]][y+tY[z][y][x]][x+tX[z][y][x]] == 0 )
		resBuf[z+tZ[z][y][x]][y+tY[z][y][x]][x+tX[z][y][x]] = theBuf[z][y][x];
	    }
	  }


	}
	break;
      default :
	VT_ErrorParse("unable to deal with such vector type\n", 0);
      }

    }
    break;
  default :
    VT_ErrorParse("unable to deal with such input type\n", 0);
  }

  





  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_FreeImage( &imres );
  VT_Free( (void**)&image );
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
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }

      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	strncpy( par->xname, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	strncpy( par->yname, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	strncpy( par->zname, argv[i], STRINGLENGTH );  
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
  par->xname[0] = '\0';
  par->yname[0] = '\0';
  par->zname[0] = '\0';
}
