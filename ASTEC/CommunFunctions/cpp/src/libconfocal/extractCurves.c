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

typedef struct local_par {
  vt_names names;
  int ox, oy, oz;
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



typedef struct {
  int x;
  int y;
  int z;
  int d;
} typePoint;



static char *usage = "image-in1 image-in2 [image-out]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t image-in1 : composantes connexes (USHORT) des courbes\n\
\t image-in2 : image de distance\n\
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
  int labelmin1, labelmax1;
  u16 *theBuf1, *theBuf2;
  int x, y, z, l, i, v;
  int a, b, c, n;
  int ix, iy, iz, ia, ib, ic;
  typePoint *thePt;
  int nbv, nbp;

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


  theBuf1 = (u16*)image1->buf;
  labelmin1 = 65536;
  labelmax1 = 0;
  for ( i=0; i<v; i++ ) {
    if ( labelmax1 < theBuf1[i] ) labelmax1 = theBuf1[i];
    if ( theBuf1[i] > 0 ) {
      if ( labelmin1 > theBuf1[i] ) labelmin1 = theBuf1[i];
    }
  }
  
  thePt = (typePoint *)malloc( 10000 * sizeof( typePoint ) );
  theBuf2 = (u16*)image2->buf;



  for ( l=labelmin1; l<=labelmax1; l++ ) {

    fprintf( stderr, "processing   cc #%4d\n", l );

    nbp = 0;
    nbv = -1;

    for ( i=0, z=0; nbv != 0 && nbv != 1 && z<image1->dim.z; z++ )
    for ( y=0;      nbv != 0 && nbv != 1 && y<image1->dim.y; y++ )
    for ( x=0;      nbv != 0 && nbv != 1 && x<image1->dim.x; x++, i++ ) {
      if ( theBuf1[i] != l ) continue;

      nbv = 0;
      for ( c = -1; c <= 1; c ++ )
      for ( b = -1; b <= 1; b ++ )
      for ( a = -1; a <= 1; a ++ ) {
	if ( z+c < 0 || z+c >= image1->dim.z ) continue;
	if ( y+b < 0 || y+b >= image1->dim.y ) continue;
	if ( x+a < 0 || x+a >= image1->dim.x ) continue;
	if ( theBuf1[i+c*(image1->dim.y*image1->dim.x)+b*image1->dim.x+a] == l )
	  nbv++;
      }
      nbv --;
      if ( nbv == 0 || nbv == 1 ) {
	ix = x;
	iy = y;
	iz = z;
      }

    }



    if ( nbv == 0 ) {
      fprintf( stderr, "isolated point for #%4d = (%3d %3d %3d)\n",
	       l, ix, iy, iz );
      fprintf( stdout, "%-10d %4d %4d %4d %d\n", 
	       l, par.ox+ix, par.oy+iy, par.oz+iz, 
	       theBuf2[iz*(image1->dim.y*image1->dim.x)
		      +iy*image1->dim.x+ix]*2 );
      fprintf( stdout, "\n" );
      continue;
    }
    
    if ( nbv == 1 ) {
      fprintf( stderr, "first point for #%4d = (%3d %3d %3d)\n",
	       l, ix, iy, iz );

      n = 0;
      thePt[n].x = ix;
      thePt[n].y = iy;
      thePt[n].z = iz;

      
      do {

	nbv = 0;
	for ( c = -1; nbv != 1 && c <= 1; c ++ )
	for ( b = -1; nbv != 1 && b <= 1; b ++ )
	for ( a = -1; nbv != 1 && a <= 1; a ++ ) {
	  if ( a == 0 && b == 0 && c == 0 ) continue;
	  if ( thePt[n].z+c < 0 || thePt[n].z+c >= image1->dim.z ) continue;
	  if ( thePt[n].y+b < 0 || thePt[n].y+b >= image1->dim.y ) continue;
	  if ( thePt[n].x+a < 0 || thePt[n].x+a >= image1->dim.x ) continue;
	  if ( theBuf1[(thePt[n].z+c)*(image1->dim.y*image1->dim.x)
		      +(thePt[n].y+b)*image1->dim.x
		      +thePt[n].x+a] != l )
	    continue;
	  if ( n == 0 ) {
	    nbv = 1;
	    ia = a;
	    ib = b;
	    ic = c;
	    continue;
	  } 
	  if ( thePt[n].z+c != thePt[n-1].z ||
	       thePt[n].y+b != thePt[n-1].y ||
	       thePt[n].x+a != thePt[n-1].x ) {
	    nbv = 1;
	    ia = a;
	    ib = b;
	    ic = c;
	  }
	}
	if ( nbv == 1 ) {
	  n++;
	  thePt[n].x = thePt[n-1].x+ia;
	  thePt[n].y = thePt[n-1].y+ib;
	  thePt[n].z = thePt[n-1].z+ic;
	}

      } while ( nbv != 0 );

    }
    
    for (i=0; i<=n; i++ ) {
      thePt[i].d = 2 * theBuf2[(thePt[i].z)*(image1->dim.y*image1->dim.x)
			      +(thePt[i].y)*image1->dim.x
			      +thePt[i].x];
      
      fprintf( stdout, "%-10d %4d %4d %4d %d\n", 
	       l, par.ox+thePt[i].x, par.oy+thePt[i].y, par.oz+thePt[i].z, thePt[i].d );
    }
    fprintf( stdout, "\n" );
    
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


      else if ( strcmp ( argv[i], "-x" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->ox) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-y" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->oy) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->oz) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
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
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = INT;
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
  par->ox = par->oy = par->oz = 0;
  par->length = 5;
  par->type = TYPE_UNKNOWN;
  par->nojct = 1;
}
