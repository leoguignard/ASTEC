/*************************************************************************
 * copy.c - copie d'images
 *
 * $Id: copy.c,v 1.3 2000/08/16 16:31:55 greg Exp $
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
 * 1996
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 * - Thu Jul 22 11:07:12 MET DST 1999 (GM)
 *   cas ou l'on voulait juste swapper ou inverser l'image
 *   on ne donne pas de type ou d'options de normalisation
 *
 *
 */


#include <vt_common.h>
#include <vt_morpho.h>

#define VT_F2I( F ) ( (F) >= 0.0 ? ((int)((F)+0.5)) : ((int)((F)-0.5)) )

typedef enum {
  VT_NONE = 0,
  _NORMA_ = 1,
  _MULT_NORMA_ = 2
} TypeOperation;

typedef struct local_par {
  vt_names names;
  int type;
  double seuil1, seuil2;
  TypeOperation type_computation;
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

static char *usage = "[image-in] [image-out] [-inv] [-swap] [-v] [-D] [-help]\n\
\t [-norma] [-mult_norma] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t -norma : normalisation de 'image-in' avant copie dans 'image-out'\n\
\t          (transformation lineaire (ax+b) de l'intensite\n\
\t -mult_norma : normalisation de 'image-in' avant copie dans 'image-out'\n\
\t          (transformation multiplicative (ax) de l'intensite\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.3 $ $Date: 2000/08/16 16:31:55 $ $Author: greg $\n";

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
  vt_image *image, imres;
  int i, s, z ;
  double min, max;
  short int *res;


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  

  s = image->dim.v * image->dim.x * image->dim.y;
  
  min = max = 0;


  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, SSHORT );
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  res = (short int*)imres.buf;

  switch ( image->type ) {
  default :
    break;
  case FLOAT :
    {
      float *buf = (float*)image->buf;
      for ( z=0; z<image->dim.z; z++ ) {
	fprintf( stderr, "%3d\r", z );
	for ( i = z*s; i < (z+1)*s; i++ ) {
	  if ( buf[i] <= par.seuil1 ) {
	    res[i] = -32767;
	  }
	  else if ( buf[i] >= par.seuil2 ) {
	  res[i] = 32767;
	  }
	  else {
	    res[i] = -32767 + (2 * 32767) / (par.seuil2-par.seuil1) * (buf[i] - par.seuil1);
	    if ( min > buf[i] ) min = buf[i];
	    if ( max < buf[i] ) max = buf[i];
	  }
	}
      }
    }
  }

  printf( " min = %f max = %f \n ", min, max );

  

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
		
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_FreeImage( &imres );
  VT_Free( (void**)&image );
  return(1);
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
	    else if ( strcmp ( argv[i], "-norma" ) == 0 ) {
                                par->type_computation = _NORMA_;
	    }
	    else if ( strcmp ( argv[i], "-mult_norma" ) == 0 ) {
	      par->type_computation = _MULT_NORMA_;
	    }

	    else if ( strcmp ( argv[i], "-sb" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse( "parsing -sb...\n", 0 );
	      status = sscanf( argv[i],"%lf",&(par->seuil1) );
	      if ( status <= 0 ) VT_ErrorParse( "parsing -sb...\n", 0 );
	    }
	    else if ( strcmp ( argv[i], "-sh" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse( "parsing -sh...\n", 0 );
	      status = sscanf( argv[i],"%lf",&(par->seuil2) );
	      if ( status <= 0 ) VT_ErrorParse( "parsing -sh...\n", 0 );
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
	/*--- noms des images ---*/
        if (nb == 0) {
                strcpy( par->names.in,  "<" );  /* standart input */
                strcpy( par->names.out, ">" );  /* standart output */
        }
        if (nb == 1)
                strcpy( par->names.out, ">" );  /* standart output */
	/*--- type de l'image resultat ---*/
	if ( (o == 1) && (s == 1) && (r == 0) ) par->type = SCHAR;
	if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
	if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
	if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
	if ( (o == 4) && (s == 1) && (r == 0) )  par->type = INT;
	if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
	if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program);
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
    par->type_computation = VT_NONE;
    par->seuil1 = -300;
    par->seuil2 =  300;
}
