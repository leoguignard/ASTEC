/*************************************************************************
 * remove-spurious.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * ?
 *
 * ADDITIONS, CHANGES
 *
 * - Tue Apr 11 19:07:13 MET DST 2000, Gregoire Malandain
 *   propagation de la taille du voxel
 */

#include <vt_common.h>

#include <vt_connexe.h>

#define _VT_CONNECTED  1
#define _VT_HYSTERESIS 2
#define _VT_SEEDPT     3
#define _VT_SEEDSIM    4

typedef struct local_par {
  vt_names names;
  int threshold;
  int min_size;
  vt_connexe cpar;
  int type;
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

static char *usage = "[image-in] [image-out] [-d %d]\n\
\t [-con %d] [-tcc %d] [-2D] \n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t if 'image-in' is equal to '-', we consider stdin\n\
\t if 'image-out' is not specified, we consider stdout\n\
\t if both are not specified, we consider stdin and stdout\n";

static char program[STRINGLENGTH];


typedef struct {
  int xmin, xmax;
  int ymin, ymax;
  int s;
  int v;
} bb2D;

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
	vt_image imtmp;
	int retour;

	/*--- initialisation des parametres ---*/
	VT_InitParam( &par );

	/* parametres */
	par.cpar.dim = VT_2D;
	par.cpar.type_connexite = N04;
	par.cpar.type_output = VT_GREY;
	  

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
	VT_InitFromImage( &imtmp, image, par.names.out, USHORT );
	if ( VT_AllocImage( &imtmp ) != 1 ) {
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_ErrorParse("unable to allocate auxiliary image\n", 0);
	}
	
	
	retour = VT_ConnectedComponents( image, &imtmp, 1, &(par.cpar) );
	VT_FreeImage( image );
	VT_Free( (void**)&image );

	
	VT_InitFromImage( &imres, &imtmp, par.names.out, UCHAR );
	if ( VT_AllocImage( &imres ) != 1 ) {
	  VT_FreeImage( &imtmp );
	  VT_ErrorParse("unable to allocate output image\n", 0);
	}
	
	
	{
	  int z, y, x, i, s, max, d;
	  unsigned short int ***buf = (unsigned short int ***)imtmp.array;
	  unsigned char ***res = (unsigned char ***)imres.array;
	  bb2D *boxes;
	  
	  s = imtmp.dim.x*imtmp.dim.y;

	  for ( z=0; z<imtmp.dim.z; z++ ) {
	    
	    max = 0;
	    for ( y=0; y<imtmp.dim.y; y++ ) 
	    for ( x=0; x<imtmp.dim.x; x++ ) {
	      if ( buf[z][y][x] > max ) max = buf[z][y][x];
	    }
	    
	    boxes = (bb2D*)malloc( (max+1)*sizeof(bb2D) );
	    for ( i=0;i<=max;i++ ) {
	      boxes[i].xmin = imtmp.dim.x - 1;
	      boxes[i].xmax = 0;
	      boxes[i].ymin = imtmp.dim.y - 1;
	      boxes[i].ymax = 0;
	      boxes[i].s = 0;
	      boxes[i].v = 1;
	    }
	    
	    for ( y=0; y<imtmp.dim.y; y++ ) 
	    for ( x=0; x<imtmp.dim.x; x++ ) {
	      i = buf[z][y][x];
	      if ( i == 0 ) continue;
	      if ( boxes[i].xmin > x ) boxes[i].xmin = x;
	      if ( boxes[i].xmax < x ) boxes[i].xmax = x;
	      if ( boxes[i].ymin > y ) boxes[i].ymin = y;
	      if ( boxes[i].ymax < y ) boxes[i].ymax = y;
	      boxes[i].s ++;
	    }
	    
	    for ( i=1;i<=max;i++ ) {
	      d = imtmp.dim.y - boxes[i].ymin;
	      if ( d > boxes[i].ymax + 1 )
		d = boxes[i].ymax + 1;
	      if ( d <= par.threshold && 
		   (par.min_size == 0 || boxes[i].s <= par.min_size) ) boxes[i].v = 0;
	    }

	    for ( y=0; y<imtmp.dim.y; y++ ) 
	    for ( x=0; x<imtmp.dim.x; x++ ) {
	      i = buf[z][y][x];
	      if ( i == 0 ) {
		res[z][y][x] = 0;
		continue;
	      }
	      if ( boxes[i].v == 0 ) res[z][y][x] = 0;
	      else res[z][y][x] = 255;
	    }
	    
	    free( boxes );
	  }
	}


	/*--- ecriture de l'image resultat ---*/
        if ( VT_WriteInrimage( &imres ) == -1 ) {
	  VT_FreeImage( &imres );
	  VT_ErrorParse("unable to write output image\n", 0);
        }
		
	/*--- liberations memoires ---*/
	VT_FreeImage( &imtmp );
	VT_FreeImage( &imres );
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
    int local_connexite = 0;
    
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
		par->cpar.verbose = 1;
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

	    /*--- connexite ---*/
	    else if ( strcmp( argv[i], "-con" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse(" parsing -con...\n", 0 );
	      status = sscanf( argv[i],"%d",&local_connexite );
	      if ( status <= 0 ) VT_ErrorParse(" parsing -con...\n", 0 );
	    }
	    /*--- taille minimum des composantes ---*/
	    else if ( strcmp( argv[i], "-tcc" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse(" parsing -tcc...\n", 0 );
	      status = sscanf( argv[i],"%d",&(par->min_size) );
	      if ( status <= 0 ) VT_ErrorParse(" parsing -tcc...\n", 0 );
	    }

	    /*--- seuil(s) ---*/
	    else if ( strcmp ( argv[i], "-d" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse( "parsing -d...\n", 0 );
	      status = sscanf( argv[i],"%d",&(par->threshold) );
	      if ( status <= 0 ) VT_ErrorParse( "parsing -d...\n", 0 );
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

    /*--- connexite ---*/
    switch ( local_connexite ) {
    case 4 :
      par->cpar.type_connexite = N04;   break;
    case 6 :
      par->cpar.type_connexite = N06;   break;
    case 8 :
      par->cpar.type_connexite = N08;   break;
    case 10 :
      par->cpar.type_connexite = N10;   break;
    case 18 :
      par->cpar.type_connexite = N18;   break;
    case 26 :
      par->cpar.type_connexite = N26;   break;
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
	par->threshold = 0;
	par->min_size = 0;
	VT_Connexe( &(par->cpar) );
	par->type = UCHAR;
}
