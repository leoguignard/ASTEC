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
#include <regionalmax.h>

typedef struct local_par {
  vt_names names;
  enumElfDistance typeDistance;
  int bool_vector;
  int type;

  int heightMaximum;
  double heightMultiplier;
  
  int method;

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
\t [-h %d] [-hm %lf]\n\
\t [-chamfer3 | -chamfer5 | -eucli ] \n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -h %d  : hauteur du maximum\n\
\t -hm %lf : coefficient multiplicateur\n\
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
	vt_image *image, imres;
	int dim[3];
	
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

	/*--- initialisation de l'image resultat ---*/
	VT_Image( &imres );
	VT_InitFromImage( &imres, image, par.names.out, USHORT );
	if ( VT_AllocImage( &imres ) != 1 ) {
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_ErrorParse("unable to allocate output image\n", 0);
	}
	
	if ( par.heightMaximum < 1 ) par.heightMaximum = 1;
	if ( (par.heightMultiplier <= 0.0) || (par.heightMultiplier > 1.0) ) 
	  par.heightMultiplier = 1.0;


	switch ( par.method ) {

	case 1 :
	  
	  if ( VT_MaximaRegionaux( image, &imres, 
				   par.heightMaximum,
				   par.heightMultiplier ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to compute method 1\n", 0);
	  }
	  break;

	case 2 :
	
	  if ( VT_MaximaRegionauxWithList( image, &imres, 
					   par.heightMaximum,
					   par.heightMultiplier ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to compute res 2\n", 0);
	  }
	  break;

	case 3 :
	  
	  dim[0] = image->dim.x;
	  dim[1] = image->dim.y;
	  dim[2] = image->dim.z;
	  if ( regionalmax( image->buf, imres.buf, USHORT, dim, par.heightMaximum,
			    par.heightMultiplier ) < 0 ) {
	    VT_FreeImage( &imres );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to compute res 3\n", 0);
	  }
	  
	}
	
	
	/*--- ecriture de l'image resultat ---*/
	if ( VT_WriteInrimage( &imres ) == -1 ) {
	  VT_FreeImage( &imres );
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_ErrorParse("unable to write output image 1\n", 0);
	}
	
	/*--- liberations memoires ---*/
	VT_FreeImage( &imres );


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
	regionalmax_setverbose();
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }
      
      else if ( strcmp ( argv[i], "-m" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->method) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
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
	status = sscanf( argv[i],"%lf",&(par->heightMultiplier) );
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

	par->method = 1;
}
