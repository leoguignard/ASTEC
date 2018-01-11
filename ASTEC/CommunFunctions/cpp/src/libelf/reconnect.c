/*************************************************************************
 * skiz.c - extraction des centres des "bulles" dans des images de mousse
 *          et calcul du squelette par zone d'influence
 *
 * $Id: skiz.c,v 1.6 2000/04/07 07:51:59 greg Exp $
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
 * - Mon Mar 27 17:57:41 MET DST 2000, Gregoire Malandain
 *   Ajout de l'enumeration enumElfDistance, afin de pouvoir
 *   calculer une distance du chamfrein 5x5x5
 *
 */

#include <vt_common.h>
#include <vt_elfskiz.h>

typedef struct {
  int outside;
  int n;
  int n_frontier;
  int n_realborder;
  int to_be_filled;
} typeComponent;

typedef struct local_par {
  vt_names names;
  enumElfDistance typeDistance;
  int bool_vector;
  int type;
  int heightMaximum;
  float heightMultiplier;
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

#ifndef NO_PROTO

static char *usage = "[image-in] [image-out] [-dist %s]\n\
\t [-h %d] [-hm %f]\n\
\t [-chamfer3 | -chamfer5 | -eucli ] \n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -h %d  : hauteur du maximum\n\
\t -hm %f : coefficient multiplicateur\n\
\t la recherche du maximum se fait en dilatant MIN( distance*hm , distance-h )\n\
\t 'en-dessous' de l'image originale, puis en soutrayant.\n\
\t -chamfer3 : distance du chanfrein, masque 3x3x3, cf Borgefors\n\
\t -chamfer5 : distance du chanfrein, masque 5x5x5\n\
\t -eucli : distance euclidienne, cf Danielsson\n\
\t\n\
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

#else

static char *usage = "";
static char *detail = "";

#endif

static char program[STRINGLENGTH];

#ifndef NO_PROTO
int main( int argc, char *argv[] )
#else
int main( argc, argv )
int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image, slice, imres, imdist;
  int max, x, y, z;
  int i;
  typeComponent *component;

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
  
  if ( par.bool_vector == 1 ) {
    fprintf( stderr, "-vect doesn't do anything\n" );
  }

  VT_Image( &slice );
  VT_InitFromImage( &slice, image, NULL, image->type );
  slice.dim.z = 1;

  VT_Image( &imres );
  VT_InitFromImage( &imres, &slice, par.names.out, USHORT );

  VT_Image( &imdist );
  VT_InitFromImage( &imdist, &slice, par.names.out, USHORT );

  (void)VT_AllocImage( &slice );
  (void)VT_AllocImage( &imres );
  (void)VT_AllocImage( &imdist );


  if ( par.heightMaximum < 1 ) par.heightMaximum = 1;
  if ( (par.heightMultiplier <= 0.0) || (par.heightMultiplier > 1.0) ) 
    par.heightMultiplier = 1.0;
	

  switch ( image->type ) {
  default :
    VT_ErrorParse("image type not handled yet\n", 0);
  case UCHAR :
    {
      unsigned char ***volArray = (unsigned char***)image->array;
      unsigned char ***theArray = (unsigned char***)slice.array;
      unsigned short int ***theSkiz = (unsigned short int ***)imres.array;
      for ( z=0; z<image->dim.z; z++ ) {
	/* copy buf 
	 */
	memcpy( slice.buf, &(volArray[z][0][0]), image->dim.x * image->dim.y );
	max = 0;
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ )
	  if ( max < theArray[0][y][x] ) max = theArray[0][y][x];
	if ( max == 0 ) {
	  fprintf( stderr, "no contours in slice #%d\n", z );
	  continue;
	}
	/* perform skiz on slice
	 */
	if ( VT_SkizFromBorderWithoutSort( &slice, &imres, &imdist, 
					   par.typeDistance, 
					   par.heightMaximum,
					   par.heightMultiplier ) != 1 ) {
	  fprintf( stderr, "error in processing skiz of slice #%d\n", z );
	  continue;
	}
	/* analyse skiz
	 */
	max = 0;
	for ( y=0; y<image->dim.y; y++ )
	for ( x=0; x<image->dim.x; x++ )
	  if ( max < theSkiz[0][y][x] ) max = theSkiz[0][y][x];
	if ( max == 0 ) {
	  fprintf( stderr, "no components in slice #%d ???\n", z );
	  continue;
	}
	fprintf( stderr, "found %d components in slice #%d\n", max, z );
	
	component = (typeComponent*)malloc( (max+1)*sizeof(typeComponent) );
	for (i=0; i<=max; i++ ) {
	  component[i].outside = 0;
	  component[i].n = 0;
	  component[i].n_frontier = 0;
	  component[i].n_realborder = 0;
	  component[i].to_be_filled = 0;
	}

	for ( y=0; y<image->dim.y; y++ ) {
	  component[ theSkiz[0][y][0] ].outside = 1;
	  component[ theSkiz[0][y][image->dim.x-1] ].outside = 1;
	}
	for ( x=0; x<image->dim.x; x++ ) {
	  component[ theSkiz[0][0][x] ].outside = 1;
	  component[ theSkiz[0][image->dim.y-1][x] ].outside = 1;
	}
	
	for ( y=1; y<image->dim.y-1; y++ )
	for ( x=1; x<image->dim.x-1; x++ ) {
	  if ( theSkiz[0][y][x] == 0 ) continue;
	  if ( component[ theSkiz[0][y][x] ].outside == 1 ) continue;
	  component[ theSkiz[0][y][x] ].n ++;
	  if ( theSkiz[0][y][x] == theSkiz[0][y][x-1]
	       && theSkiz[0][y][x] == theSkiz[0][y][x+1]
	       && theSkiz[0][y][x] == theSkiz[0][y-1][x]
	       && theSkiz[0][y][x] == theSkiz[0][y+1][x] ) continue;
	  component[ theSkiz[0][y][x] ].n_frontier ++;
	  if ( (theSkiz[0][y][x-1] == 0 || theSkiz[0][y][x-1] == theSkiz[0][y][x])
	       && (theSkiz[0][y][x+1] == 0 || theSkiz[0][y][x+1] == theSkiz[0][y][x])
	       && (theSkiz[0][y-1][x] == 0 || theSkiz[0][y-1][x] == theSkiz[0][y][x])
	       && (theSkiz[0][y+1][x] == 0 || theSkiz[0][y+1][x] == theSkiz[0][y][x]) )
	    component[ theSkiz[0][y][x] ].n_realborder ++;
	}

	/* analyse components
	 */
	for ( i=1; i<=max; i++ ) {
	  if ( 0 )
	  fprintf( stderr, "   *** #%2d out = %d, frontier = %4d real = %4d\n",
		   i, component[i].outside, component[i].n_frontier,
		   component[i].n_realborder );
	  if ( component[i].outside == 1 ) continue;
	  if ( component[i].n_frontier == component[i].n_realborder ) continue;
	  if ( component[i].n_frontier > 2*component[i].n_realborder ) continue;
	  if ( 1 ) {
	    fprintf( stderr, "   component #%d in slice #%d seems open ",
		     i, z );
	    fprintf( stderr, "(vol=%3d, frontier=%2d, real_frontier=%2d)\n",
		     component[i].n, component[i].n_frontier,
		     component[i].n_realborder );
	  }
	  component[i].to_be_filled = 1;
	}

	/* fill components
	 */
	for ( y=1; y<image->dim.y-1; y++ )
	for ( x=1; x<image->dim.x-1; x++ ) {
	  if ( theSkiz[0][y][x] == 0 ) continue;
	  if ( component[ theSkiz[0][y][x] ].outside == 1 ) continue;
	  if ( component[ theSkiz[0][y][x] ].n_frontier == component[ theSkiz[0][y][x] ].n_realborder
	      || component[ theSkiz[0][y][x] ].to_be_filled == 1 )
	    volArray[z][y][x] = 255;
	}
	
	/* fill residual holes
	 */
	for ( y=1; y<image->dim.y-1; y++ )
	for ( x=1; x<image->dim.x-1; x++ ) {
	  if ( theSkiz[0][y][x] == 0 ) continue;
	  if ( component[ theSkiz[0][y][x] ].to_be_filled == 0 ) continue;
	  if ( theSkiz[0][y][x-1] != 0 && theSkiz[0][y][x-1] != theSkiz[0][y][x] ) {
	    i = 0;
	    if ( volArray[z][y  ][x-2] == 255 ) i++;
	    if ( volArray[z][y  ][x  ] == 255 ) i++;
	    if ( volArray[z][y-1][x-1] == 255 ) i++;
	    if ( volArray[z][y+1][x-1] == 255 ) i++;
	    if ( i > 2 ) {
	      volArray[z][y  ][x-1] = 100;
	      fprintf( stderr, "   fill (%d %d %d)\n", x-1, y, z );
	    }
	  }
	  if ( theSkiz[0][y][x+1] != 0 && theSkiz[0][y][x+1] != theSkiz[0][y][x] ) {
	    i = 0;
	    if ( volArray[z][y  ][x  ] == 255 ) i++;
	    if ( volArray[z][y  ][x+2] == 255 ) i++;
	    if ( volArray[z][y-1][x+1] == 255 ) i++;
	    if ( volArray[z][y+1][x+1] == 255 ) i++;
	    if ( i > 2 ) {
	      volArray[z][y  ][x+1] = 100;
	      fprintf( stderr, "   fill (%d %d %d)\n", x+1, y, z );
	    }
	  }
	  if ( theSkiz[0][y-1][x] != 0 && theSkiz[0][y-1][x] != theSkiz[0][y][x] ) {
	    i = 0;
	    if ( volArray[z][y-1][x-1] == 255 ) i++;
	    if ( volArray[z][y-1][x+1] == 255 ) i++;
	    if ( volArray[z][y-2][x  ] == 255 ) i++;
	    if ( volArray[z][y  ][x  ] == 255 ) i++;
	    if ( i > 2 ) {
	      volArray[z][y-1][x  ] = 100;
	      fprintf( stderr, "   fill (%d %d %d)\n", x, y-1, z );
	    }
	  }
	  if ( theSkiz[0][y+1][x] != 0 && theSkiz[0][y+1][x] != theSkiz[0][y][x] ) {
	    i = 0;
	    if ( volArray[z][y+1][x-1] == 255 ) i++;
	    if ( volArray[z][y+1][x+1] == 255 ) i++;
	    if ( volArray[z][y  ][x  ] == 255 ) i++;
	    if ( volArray[z][y+2][x  ] == 255 ) i++;
	    if ( i > 2 ) {
	      volArray[z][y+1][x  ] = 100;
	      fprintf( stderr, "   fill (%d %d %d)\n", x, y+1, z );
	    }
	  }

	}

      }
    } /* case UCHAR */

  } /* switch */

  /*--- ecriture de l'image resultat ---*/
  strcpy( image->name, par.names.out );
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

#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par )
#else
static void VT_Parse( argc, argv, par )
int argc;
char *argv[];
local_par *par;
#endif
{
  int i, nb, status, connexite = 0;
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
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }
      
      else if ( strcmp ( argv[i], "-dist" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -dist...\n", 0 );
	strcpy( par->names.ext, argv[i] );
      }
      

      /*--- seuil ---*/

      else if ( strcmp ( argv[i], "-h" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -h...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->heightMaximum) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -h...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-hm" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -hm...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->heightMultiplier) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -hm...\n", 0 );
      }


      /*--- type de distance ---*/
      else if ( strcmp ( argv[i], "-chamfer3" ) == 0 ) {
	par->typeDistance = _CHAMFER_3_;
      }
      else if ( strcmp ( argv[i], "-chamfer5" ) == 0 ) {
	par->typeDistance = _CHAMFER_5_;
      }
      else if ( strcmp ( argv[i], "-eucli" ) == 0 ) {
	par->typeDistance = _EUCLIDIENNE_;
      }

      /*--- lecture du type de l'image ---*/
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
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
  /*--- conversion du type pour le calcul de l'euclidean mapping ---*/
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = INT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program);
}

#ifndef NO_PROTO
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

#ifndef NO_PROTO
static void VT_InitParam( local_par *par )
#else
static void VT_InitParam( par )
local_par *par;
#endif
{
	VT_Names( &(par->names) );
	par->bool_vector = 0;
	par->type = TYPE_UNKNOWN;
	par->heightMaximum = 1;
	par->heightMultiplier = 1.0;
}
