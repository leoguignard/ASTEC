/*************************************************************************
 * minimum.c -
 *
 * $Id: upperView.c,v 1.1 2001/04/04 07:25:45 greg Exp $
 *
 * Copyrigh (c) INRIA 1999
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

#include <vt_common.h>

typedef enum {
  _Z_,
  _COLOR_,
  _ZMAX_, 
  _MAX_
} enumTypeView;

typedef struct local_par {
  vt_names names;
  enumTypeView typeView;
  int type;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-z | -color | -zmax | -max]\n\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.1 $ $Date: 2001/04/04 07:25:45 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  int x, y, z;
  int zmax, max;
  u8 ***resZBuf = NULL;

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

  switch ( par.typeView ) {
  default :
  case _Z_ :
    VT_InitFromImage( &imres, image, par.names.out, UCHAR);
    break;
  case _COLOR_ :
    VT_InitFromImage( &imres, image, par.names.out, image->type );
  }
  imres.dim.z = 1;

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  switch ( par.typeView ) {
  default :
  case _Z_ :
      resZBuf = (u8***)imres.array;
      break;
  case _COLOR_ :
    break;
  }

  switch ( image->type ) {
  default :
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to deal with such image type\n", 0 );
  case UCHAR :
    {
      u8 ***theBuf = (u8***)image->array;
      u8 ***resBuf = (u8***)imres.array;
      for ( y=0; y<image->dim.y; y++ )
      for ( x=0; x<image->dim.x; x++ )
	resBuf[0][y][x] = 0;
      switch ( par.typeView ) {
      default :
      case _Z_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  for ( z=0; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > 0 ) {
	      resZBuf[0][y][x] = z+1;
	      z = image->dim.z;
	    }
	  }
	}
	break;
      case _ZMAX_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  zmax = 0;
	  max = theBuf[0][y][x];
	  for ( z=1; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > max ) {
	      zmax = z;
	      max = theBuf[z][y][x];
	    }
	  }
	  resBuf[0][y][x] = zmax+1;
	}
	break;
      case _MAX_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  zmax = 0;
	  max = theBuf[0][y][x];
	  for ( z=1; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > max ) {
	      zmax = z;
	      max = theBuf[z][y][x];
	    }
	  }
	  resBuf[0][y][x] = max;
	}
	break;
      case _COLOR_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  for ( z=0; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > 0 ) {
	      resBuf[0][y][x] = theBuf[z][y][x];
	      z = image->dim.z;
	    }
	  }
	}
	break;
      }
    }
    break;
   case USHORT :
    {
      u16 ***theBuf = (u16***)image->array;
      u16 ***resBuf = (u16***)imres.array;
      for ( y=0; y<image->dim.y; y++ )
      for ( x=0; x<image->dim.x; x++ )
	resBuf[0][y][x] = 0;
      switch ( par.typeView ) {
      default :
      case _Z_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  for ( z=0; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > 0 ) {
	      resZBuf[0][y][x] = z+1;
	      z = image->dim.z;
	    }
	  }
	}
	break;
      case _ZMAX_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  zmax = 0;
	  max = theBuf[0][y][x];
	  for ( z=1; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > max ) {
	      zmax = z;
	      max = theBuf[z][y][x];
	    }
	  }
	  resBuf[0][y][x] = zmax+1;
	}
	break;
      case _MAX_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  zmax = 0;
	  max = theBuf[0][y][x];
	  for ( z=1; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > max ) {
	      zmax = z;
	      max = theBuf[z][y][x];
	    }
	  }
	  resBuf[0][y][x] = max;
	}
	break;
      case _COLOR_ :
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ ) {
	  for ( z=0; z<image->dim.z; z++ ) {
	    if ( theBuf[z][y][x] > 0 ) {
	      resBuf[0][y][x] = theBuf[z][y][x];
	      z = image->dim.z;
	    }
	  }
	}
	break;
      }
    }
    break;
    
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
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }

      else if ( strcmp ( argv[i], "-zmax" ) == 0 ) {
	par->typeView = _ZMAX_;
      }
      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
	par->typeView = _MAX_;
      }
      else if ( strcmp ( argv[i], "-color" ) == 0 ) {
	par->typeView = _COLOR_;
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 && argv[i][2] == '\0' ) {
	par->typeView = _Z_;
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
  par->typeView = _Z_;
  par->type = TYPE_UNKNOWN;
}
