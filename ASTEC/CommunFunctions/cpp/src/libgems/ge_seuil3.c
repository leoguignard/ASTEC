#include <vt_common.h>

#include <vt_mip.h>

#include <vt_copy.h>

#include <vt_histo.h>

#include <vt_recfilters.h>

#include <vt_geseuil1.h>
#include <vt_gaussienne.h>

typedef struct local_par {
    vt_names names;
    vt_seuil1 spar;
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

static char *usage = "[image-in] [image-out]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
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
	local_par par;
	char image_name[STRINGLENGTH];
	vt_image *image, imXY, imXZ, imZY;
	vt_image auxXY, auxXZ, auxZY;
	vt_histo histoXY, histoXZ, histoZY;
	vt_image imres;
	vt_image imfloat;
	vt_3m m;
	i32 *buf;
	int i;
	double M, K, S;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- lecture de l'image d'entree ---*/
	image = _VT_Inrimage( par.names.in );
	if ( image == (vt_image*)NULL ) 
		VT_ErrorParse("unable to read input image\n", 0);

	/*--- operations eventuelles sur l'image d'entree ---*/
	if ( par.names.inv == 1 )  VT_InverseImage( image );
	if ( par.names.swap == 1 ) VT_SwapImage( image );

	/*--- initialisation des images MIP ---*/
	VT_Image( &imXY );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".xy_mip" );
	VT_InitImage( &imXY, image_name, image->dim.x, image->dim.y, (int)1, image->type );

	VT_Image( &imXZ );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".xz_mip" );
	VT_InitImage( &imXZ, image_name, image->dim.x, image->dim.z, (int)1, image->type );

	VT_Image( &imZY );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".zy_mip" );
	VT_InitImage( &imZY, image_name, image->dim.z, image->dim.y, (int)1, image->type );

	if ( VT_AllocImage( &imXY ) != 1 ) {
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }
	if ( VT_AllocImage( &imXZ ) != 1 ) {
	    VT_FreeImage( &imXY );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }	
	if ( VT_AllocImage( &imZY ) != 1 ) {
	    VT_FreeImage( &imXY );
	    VT_FreeImage( &imXZ );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }
	VT_Image( &auxXY );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".xy_mip" );
	VT_InitImage( &auxXY, image_name, image->dim.x, image->dim.y, (int)1, (int)UCHAR );

	VT_Image( &auxXZ );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".xz_mip" );
	VT_InitImage( &auxXZ, image_name, image->dim.x, image->dim.z, (int)1, (int)UCHAR );

	VT_Image( &auxZY );
	VT_Strncpy( image_name, par.names.in, STRINGLENGTH );
	(void)strcat( image_name, ".zy_mip" );
	VT_InitImage( &auxZY, image_name, image->dim.z, image->dim.y, (int)1, (int)UCHAR );

	if ( VT_AllocImage( &auxXY ) != 1 ) {
	    VT_FreeImage( &imXY );
	    VT_FreeImage( &imXZ );
	    VT_FreeImage( &imZY );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }
	if ( VT_AllocImage( &auxXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );
	    VT_FreeImage( &imXY );
	    VT_FreeImage( &imXZ );
	    VT_FreeImage( &imZY );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }	
	if ( VT_AllocImage( &auxZY ) != 1 ) {
	    VT_FreeImage( &auxXY );
	    VT_FreeImage( &auxXZ );
	    VT_FreeImage( &imXY );
	    VT_FreeImage( &imXZ );
	    VT_FreeImage( &imZY );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }

	/*--- initialisation histogramme ---*/
	VT_Histo( &histoXY );   histoXY.type = (int)UCHAR;
	VT_Histo( &histoXZ );   histoXZ.type = (int)UCHAR;
	VT_Histo( &histoZY );   histoZY.type = (int)UCHAR;
	if ( VT_AllocHisto( &histoXY ) != 1 ) {
	    VT_FreeImage( &auxXY );   VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );    VT_FreeImage( &imXZ );    VT_FreeImage( &imZY );
	    VT_FreeImage( image );    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}
	if ( VT_AllocHisto( &histoXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );   VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );    VT_FreeImage( &imXZ );    VT_FreeImage( &imZY );
	    VT_FreeImage( image );    VT_Free( (void**)&image );
	    VT_FreeHisto( &histoXY );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}
	if ( VT_AllocHisto( &histoZY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );    VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}

	/*--- initialisation output ---*/
	VT_Image( &imres );
	VT_InitImage( &imres, par.names.out, histoXY.size, 1, 1, INT );
	if ( VT_AllocImage( &imres ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );     VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );      VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );   VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to allocate output image\n", 0);
	}
	VT_Image( &imfloat );
	VT_InitImage( &imfloat, par.names.out, histoXY.size, 1, 1, FLOAT );
	(void)VT_AllocImage( &imfloat );

	/*------------------------------------------------------------
	  ------------------------------------------------------------
	  ------------------------------------------------------------*/

	/*--- calcul des max, min et moy de l'image ---*/
	if ( VT_3m( image, &m ) == -1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("computation of minimum and maximum not possible", 0);
	}
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," %s : minimum = %d, maximum = %d", par.names.in, (int)(m.min), (int)(m.max) );
	    VT_Message( message, program );
	}
	
	/*--- calcul des images MIP ---*/
	if ( VT_MIP( image, &imXY, &imXZ, &imZY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("computation of mip projections not possible", 0);
	}

	/*--- normalisation des images MIP ---*/
	if ( VT_NormaImage( &imXY, &auxXY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   
	if ( VT_NormaImage( &imXZ, &auxXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   
	if ( VT_NormaImage( &imZY, &auxZY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   

	/*--- calcul histogramme ---*/
	for ( i = 0; i < histoXY.size; i++ )
	    histoXY.buf[i] = histoXZ.buf[i] = histoZY.buf[i] = 0;

	if ( VT_Image2Histo( &histoXY, &auxXY, &(par.spar) ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to fill histogram from image\n", 0);
	}
	if ( VT_Image2Histo( &histoXZ, &auxXZ, &(par.spar) ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to fill histogram from image\n", 0);
	}
	if ( VT_Image2Histo( &histoZY, &auxZY, &(par.spar) ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to fill histogram from image\n", 0);
	}

	/*--- copie buffer ---*/
	buf = (i32*)imres.buf;
	for ( i = 0; i < histoXY.size; i++ )
	    buf[i] = histoXY.buf[i] + histoXZ.buf[i] + histoZY.buf[i];
	
	/*--- pas de point a 0 ---*/
	buf[0] = 0;

	/*--- lissage de l'histogramme ---*/
	if ( VT_RecFilterOnImage( &imres, &imfloat, &(par.spar.par_filt) ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to smooth output image\n", 0);
	}
	
	/*--- on ajuste une gaussienne ---*/
	if ( VT_FitGaussienne( (r32*)imfloat.buf, histoXY.size, &K, &S, &M ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	    VT_FreeImage( image );      VT_Free( (void**)&image );   VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to fit gaussian\n", 0);
	}
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," gaussienne : m = %f k = %f s = %f", (float)M, (float)K, (float)S );
	    VT_Message( message, program );
	    sprintf( message," seuil = %f", (float)(M * m.max / 255.0 ) );
	    VT_Message( message, program );
	}


	VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	VT_FreeImage( &imXY );      VT_FreeImage( &imXZ );       VT_FreeImage( &imZY );
	VT_FreeImage( image );      VT_Free( (void**)&image );   
	VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 

	/*--- ecriture de l'image resultat ---*/
	if ( VT_WriteInrimage( &imres ) == -1 ) {
	    VT_FreeImage( &imres );
	    VT_ErrorParse("unable to write output image\n", 0);
	}
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
    char text[STRINGLENGTH];
    
    if ( VT_CopyName( program, argv[0] ) != 1 )
	VT_Error("Error while copying program name", (char*)NULL);
    if ( argc == 1 ) VT_ErrorParse("\n", 0 );
    
    /*--- initialisation des parametres ---*/
    VT_InitParam( par );

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
	VT_Seuil1( &(par->spar) );
}
