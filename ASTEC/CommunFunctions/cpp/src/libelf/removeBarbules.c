/*****************************************************************************
 * removeBarbules.c -
 *
 * $Id: removeBarbules.c,v 1.5 2001/04/13 18:12:12 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Feb, 24 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 */

#include <vt_common.h>
#include <vt_elfrelabel.h>
#include <vt_elfbarbules.h>


typedef struct local_par {
  vt_names names;
  int length;
  int type;
  int nojct;
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

static char *usage = "image-in1 image-in2 [image-out]\n\
\t [-l %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t image-in1 : composantes connexes (USHORT) des courbes+frontieres\n\
\t image-in2 : composantes connexes (USHORT) des jonctions\n\
\t -l %d : maximal length of barbules\n\
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
 $Revision: 1.5 $ $Date: 2001/04/13 18:12:12 $ $Author: greg $\n";

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
  vt_image *image1, *image2;
  int labelmin1, labelmax1, labelmin2, labelmax2;
  u16 *theBuf1, *theBuf2;
  int i, v, offset;
  int theDim[3];

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  


  /*--- lecture de l'image d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image 1\n", 0);
  if ( image1->type != USHORT ) {
    VT_ErrorParse("bad type for input image 1\n", 0);
  }

  
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image 2\n", 0);
  if ( image2->type != USHORT ) {
    VT_ErrorParse("bad type for input image 2\n", 0);
  }

  if ( image1->dim.x != image2->dim.x ||
       image1->dim.y != image2->dim.y ||
       image1->dim.z != image2->dim.z ) {
    VT_ErrorParse("input images must have same dimensions\n", 0);
  }
  
  v = image1->dim.x * image1->dim.y * image1->dim.z;

  theDim[0] = image1->dim.x;
  theDim[1] = image1->dim.y;
  theDim[2] = image1->dim.z;
  if ( RelabelConnectedComponentsSortBySize( image1->buf,
					     image1->type,
					     theDim ) != 1 ) {
    VT_ErrorParse("unable to sort ccs of image 1\n", 0 );
  }

  theBuf1 = (u16*)image1->buf;
  labelmin1 = 65536;
  labelmax1 = 0;
  for ( i=0; i<v; i++ ) {
    if ( labelmax1 < theBuf1[i] ) labelmax1 = theBuf1[i];
    if ( theBuf1[i] > 0 ) {
      if ( labelmin1 > theBuf1[i] ) labelmin1 = theBuf1[i];
    }
  }
  
  theBuf2 = (u16*)image2->buf;
  labelmin2 = 65536;
  labelmax2 = 0;
  for ( i=0; i<v; i++ ) {
    if ( labelmax2 < theBuf2[i] ) labelmax2 = theBuf2[i];
    if ( theBuf2[i] > 0 ) {
      if ( labelmin2 > theBuf2[i] ) labelmin2 = theBuf2[i];
    }
  }
  
  printf(" image 1, labels %5d -> %5d \n", labelmin1, labelmax1 );
  printf(" image 2, labels %5d -> %5d \n", labelmin2, labelmax2 );
  
  offset = 100;
  if ( labelmax1+offset+labelmax2 > 65535 ) {
    offset = 0;
    if ( labelmax1+offset+labelmax2 > 65535 ) {
      VT_ErrorParse( "unable to fuse labels into a single image\n", 0 );
    }
  }
  
  for ( i=0; i<v; i++ ) {
    if ( theBuf2[i] > 0 ) {
      if ( theBuf1[i] > 0 ) {
	fprintf( stderr, "WARNING: point %d was not null in image 1 (%s)\n",
		 i, par.names.in );
      }
      theBuf1[i] = labelmax1+offset+theBuf2[i];
    }
  }
  labelmin2 += labelmax1+offset;
  labelmax2 += labelmax1+offset;

  printf(" image 1/1, labels %5d -> %5d \n", labelmin1, labelmax1 );
  printf(" image 2/1, labels %5d -> %5d \n", labelmin2, labelmax2 );

  
  VT_FreeImage( image2 );
  VT_Free( (void**)&image2 );




  if ( VT_RemoveBarbules( image1, labelmin1, labelmax1,
			  labelmin2, labelmax2, par.length ) != 1 ) {
    VT_FreeImage( image1 );
    VT_Free( (void**)&image1 );
    VT_ErrorParse("unable to compute\n", 0);
  }

  if ( par.nojct == 0 ) {
    for ( i=0; i<v; i++ ) {
      if( theBuf1[i] >= labelmin2 ) theBuf1[i] = 0;
    }
  }


  /*--- ecriture de l'image resultat ---*/
  (void)VT_CopyName( image1->name, par.names.out );
  /* sprintf( image1->name, "%s", par.names.out ); */
  if ( VT_WriteInrimage( image1 ) == -1 ) {
    VT_FreeImage( image1 );
    VT_Free( (void**)&image1 );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image1 );
  VT_Free( (void**)&image1 );
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
	_VT_DEBUG_ = 2;
      }
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }


      else if ( strcmp ( argv[i], "-l" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -l...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->length) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -l...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-nojct" ) == 0 ) {
	par->nojct = 0;
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
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else if ( nb == 2 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb != 3) {
    VT_ErrorParse("must specify 3 file names\n", 0 );
  }
  
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
  par->length = 5;
  par->type = TYPE_UNKNOWN;
  par->nojct = 1;
}
