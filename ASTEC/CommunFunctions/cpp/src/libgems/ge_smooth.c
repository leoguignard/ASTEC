#include <vt_common.h>

#include <vt_gesmooth.c>

typedef struct local_par {
    float sb;
    int dim;
} local_par;

/*-------Definition des fonctions statiques----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], vt_names *n, local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[image-in] [image-out] [-sb %f] [-2D]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";
static char *detail = "\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -sb %f : seuil bas\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";
static char program[STRINGLENGTH];

#ifndef NO_PROTO
int main( int argc, char *argv[] )
#else
int main( argc, argv )
int argc;
char *argv[];
#endif
{
	vt_names names;
	local_par par;
	vt_image *image, imres;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &names, &par );

	/*--- lecture de l'image d'entree ---*/
	image = _VT_Inrimage( names.in );
	if ( image == (vt_image*)NULL ) 
		VT_ErrorParse("unable to read input image\n", 0);

	/*--- operations eventuelles sur l'image d'entree ---*/
	if ( names.inv == 1 )  VT_InverseImage( image );
	if ( names.swap == 1 ) VT_SwapImage( image );

	/*--- initialisation de l'image resultat ---*/
	VT_Image( &imres );
	VT_InitImage( &imres, names.out, image->dim.x, image->dim.y, image->dim.z, (int)UCHAR );
	if ( VT_AllocImage( &imres ) != 1 ) {
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("unable to allocate output image\n", 0);
	}

	if ( VT_FastSmoothUC( &imres, image, (int)(par.sb + 0.5), par.dim ) != 1 ) {
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_FreeImage( &imres );
	    VT_ErrorParse("unable to process input image\n", 0);
	}

	VT_FreeImage( image );
	VT_Free( (void**)&image );
	
	/*--- ecriture de l'image resultat ---*/
	if ( VT_WriteInrimage( &imres ) == -1 ) {
	    VT_FreeImage( &imres );
	    VT_ErrorParse("unable to write output image\n", 0);
	}
		
	/*--- liberations memoires ---*/
	VT_FreeImage( &imres );
	return(1);
}


#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], vt_names *n, local_par *par )
#else
static void VT_Parse( argc, argv, n, par )
int argc;
char *argv[];
vt_names *n;
local_par *par;
#endif
{
	int i, nb, status;
	char text[STRINGLENGTH];

	if ( VT_CopyName( program, argv[0] ) != 1 )
		VT_Error("Error while copying program name", (char*)NULL);
	if ( argc == 1 ) VT_ErrorParse("\n", 0 );

	/*--- initialisation des parametres ---*/
	VT_Names( n );
	VT_InitParam( par );
	/*--- lecture des parametres ---*/
	i = 1; nb = 0;
	while ( i < argc ) {
	    if ( argv[i][0] == '-' ) {
		if ( argv[i][1] == '\0' ) {
		    if ( nb == 0 ) {
				/*--- standart input ---*/
			strcpy( n->in, "<" );
			nb += 1;
		    }
		}
		else if ( strcmp ( argv[i], "-help" ) == 0 ) {
		    VT_ErrorParse("\n", 1);
		}
		else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
                                n->inv = 1;
		}
		else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
		    n->swap = 1;
		}
		else if ( strcmp ( argv[i], "-v" ) == 0 ) {
		    _VT_VERBOSE_ = 1;
		}
		else if ( strcmp ( argv[i], "-D" ) == 0 ) {
		    _VT_DEBUG_ = 1;
		}
		/*--- seuils ---*/
		else if ( strcmp ( argv[i], "-sb" ) == 0 ) {
		    i += 1;
		    if ( i >= argc)    VT_ErrorParse( "parsing -sb...\n", 0 );
		    status = sscanf( argv[i],"%f",&(par->sb) );
		    if ( status <= 0 ) VT_ErrorParse( "parsing -sb...\n", 0 );
		}
		/*--- dimension ---*/
		else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
		    par->dim = VT_2D;
		}
		/*--- option inconnue ---*/
		else {
		    sprintf(text,"unknown option %s\n",argv[i]);
		    VT_ErrorParse(text, 0);
		}
	    }
	    else if ( argv[i][0] != 0 ) {
		if ( nb == 0 ) { 
		    strncpy( n->in, argv[i], STRINGLENGTH );  
		    nb += 1;
		}
		else if ( nb == 1 ) {
		    strncpy( n->out, argv[i], STRINGLENGTH );  
		    nb += 1;
		}
		else 
		    VT_ErrorParse("too much file names when parsing\n", 0 );
	    }
	    i += 1;
	}
	/*--- noms des images ---*/
        if (nb == 0) {
                strcpy( n->in,  "<" );  /* standart input */
                strcpy( n->out, ">" );  /* standart output */
        }
        if (nb == 1)
                strcpy( n->out, ">" );  /* standart output */
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
	par->sb = (float)1.0;
	par->dim = VT_3D;
}
