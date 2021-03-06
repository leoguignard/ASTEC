#include <vt_common.h>

#include <vt_mip.h>

#include <vt_copy.h>

#include <vt_histo.h>

#include <vt_recfilters.h>

#include <vt_geseuil1.h>
#include <vt_gaussienne.h>

typedef struct local_par {
    char name_mip_XY[STRINGLENGTH];
    char name_mip_XZ[STRINGLENGTH];
    char name_mip_ZY[STRINGLENGTH];
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

static char *usage = "mip_XY mip_XZ mip_ZY [image-out]\n\
\t [-eroc %d] [-eroi %d] [-csb %f] [-csh %f] [-dilc %d] [-dili %d]\n\
\t [-v] [-D] [-help]";

static char *detail = "\
\t mip_XY mip_XZ mip_ZY : les trois images MIP\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t -eroc %d : connexite pour l'erosion de l'image MIP (defaut = 26)\n\
\t -eroi %d : nombre d'iterations pour l'erosion de l'image MIP (defaut = 1)\n\
\t -csb %f  : coefficient multiplicateur du max pour le calcul\n\
\t            du seuil bas du seuillage par hysteresis (defaut = 0.1666)\n\
\t -csh %f  : coefficient multiplicateur du max pour le calcul\n\
\t            du seuil haut du seuillage par hysteresis (defaut = 0.3333)\n\
\t -dilc %d : connexite pour la dilatation de l'image des contours (defaut = 6)\n\
\t -dili %d : nombre d'iterations pour la dilatation de l'image des contours (defaut = 1)\n\
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
	vt_image *imXY, *imXZ, *imZY;
	vt_image auxXY, auxXZ, auxZY;
	vt_histo histoXY, histoXZ, histoZY;
	vt_image imres;
	vt_image imfloat;
	vt_3m mXY, mXZ, mZY;
	i32 *buf;
	int i;
	double M, K, S;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- lectures des images d'entrees ---*/
	imXY = _VT_Inrimage( par.name_mip_XY );
	if ( imXY == (vt_image*)NULL )  VT_ErrorParse("unable to read input image\n", 0);
	imXZ = _VT_Inrimage( par.name_mip_XZ );
	if ( imXZ == (vt_image*)NULL )  VT_ErrorParse("unable to read input image\n", 0);
	imZY = _VT_Inrimage( par.name_mip_ZY );
	if ( imZY == (vt_image*)NULL )  VT_ErrorParse("unable to read input image\n", 0);
	
	VT_Image( &auxXY );
	(void)strcpy( image_name, "ge_seuil3_3mip" );
	(void)strcat( image_name, ".xy_mip" );
	VT_InitImage( &auxXY, image_name, imXY->dim.x, imXY->dim.y, (int)1, (int)UCHAR );

	VT_Image( &auxXZ );
	(void)strcpy( image_name, "ge_seuil3_3mip" );
	(void)strcat( image_name, ".xz_mip" );
	VT_InitImage( &auxXZ, image_name, imXZ->dim.x, imXZ->dim.y, (int)1, (int)UCHAR );

	VT_Image( &auxZY );
	(void)strcpy( image_name, "ge_seuil3_3mip" );
	(void)strcat( image_name, ".zy_mip" );
	VT_InitImage( &auxZY, image_name, imZY->dim.x, imZY->dim.y, (int)1, (int)UCHAR );

	if ( VT_AllocImage( &auxXY ) != 1 ) {
	    VT_FreeImage( imXY );
	    VT_FreeImage( imXZ );
	    VT_FreeImage( imZY );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }
	if ( VT_AllocImage( &auxXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );
	    VT_FreeImage( imXY );
	    VT_FreeImage( imXZ );
	    VT_FreeImage( imZY );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }	
	if ( VT_AllocImage( &auxZY ) != 1 ) {
	    VT_FreeImage( &auxXY );
	    VT_FreeImage( &auxXZ );
	    VT_FreeImage( imXY );
	    VT_FreeImage( imXZ );
	    VT_FreeImage( imZY );
	    VT_ErrorParse("unable to allocate first output image\n", 0);
        }

	/*--- initialisation histogramme ---*/
	VT_Histo( &histoXY );   histoXY.type = (int)UCHAR;
	VT_Histo( &histoXZ );   histoXZ.type = (int)UCHAR;
	VT_Histo( &histoZY );   histoZY.type = (int)UCHAR;
	if ( VT_AllocHisto( &histoXY ) != 1 ) {
	    VT_FreeImage( &auxXY );   VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );    VT_FreeImage( imXZ );    VT_FreeImage( imZY );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}
	if ( VT_AllocHisto( &histoXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );   VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );    VT_FreeImage( imXZ );    VT_FreeImage( imZY );
	    VT_FreeHisto( &histoXY );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}
	if ( VT_AllocHisto( &histoZY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );   VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );    VT_FreeImage( imZY );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );
	    VT_ErrorParse("unable to allocate histogram\n", 0);
	}

	/*--- initialisation output ---*/
	VT_Image( &imres );
	VT_InitImage( &imres, par.names.out, histoXY.size, 1, 1, INT );
	if ( VT_AllocImage( &imres ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );     VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );      VT_FreeImage( imZY );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );   VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to allocate output image\n", 0);
	}
	VT_Image( &imfloat );
	VT_InitImage( &imfloat, par.names.out, histoXY.size, 1, 1, FLOAT );
	(void)VT_AllocImage( &imfloat );

	/*------------------------------------------------------------
	  ------------------------------------------------------------
	  ------------------------------------------------------------*/

	/*--- normalisation des images MIP ---*/
	if ( VT_NormaImage( imXY, &auxXY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
	    VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   
	if ( VT_NormaImage( imXZ, &auxXZ ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
	    VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   
	if ( VT_NormaImage( imZY, &auxZY ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
	    VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to normalize MIP image\n", 0);
	}   

	/*--- calcul histogramme ---*/
	for ( i = 0; i < histoXY.size; i++ )
	    histoXY.buf[i] = histoXZ.buf[i] = histoZY.buf[i] = 0;

	(void)VT_3m( imXY, &mXY );
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," %s : minimum = %d, maximum = %d", imXY->name, (int)(mXY.min), (int)(mXY.max) );
	    VT_Message( message, program );
	}
	if ( (int)(mXY.min) == 0 ) {
	    if ( VT_Image2Histo( &histoXY, &auxXY, &(par.spar) ) != 1 ) {
		VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
		VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
		VT_FreeImage( &imres );
		VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
		VT_ErrorParse("unable to fill histogram from image\n", 0);
	    }
	}
	(void)VT_3m( imXZ, &mXZ );
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," %s : minimum = %d, maximum = %d", imXZ->name, (int)(mXZ.min), (int)(mXZ.max) );
	    VT_Message( message, program );
	}
	if ( (int)(mXZ.min) == 0 ) {
	    if ( VT_Image2Histo( &histoXZ, &auxXZ, &(par.spar) ) != 1 ) {
		VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
		VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
		VT_FreeImage( &imres );
		VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
		VT_ErrorParse("unable to fill histogram from image\n", 0);
	    }
	}
	(void)VT_3m( imZY, &mZY );
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," %s : minimum = %d, maximum = %d", imZY->name, (int)(mZY.min), (int)(mZY.max) );
	    VT_Message( message, program );
	}
	if ( (int)(mZY.min) == 0 ) {
	    if ( VT_Image2Histo( &histoZY, &auxZY, &(par.spar) ) != 1 ) {
		VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
		VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
		VT_FreeImage( &imres );
		VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
		VT_ErrorParse("unable to fill histogram from image\n", 0);
	    }
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
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
	    VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to smooth output image\n", 0);
	}
	
	/*--- on ajuste une gaussienne ---*/
	if ( VT_FitGaussienne( (r32*)imfloat.buf, histoXY.size, &K, &S, &M ) != 1 ) {
	    VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	    VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
	    VT_FreeImage( &imres );
	    VT_FreeHisto( &histoXY );   VT_FreeHisto( &histoXZ );    VT_FreeHisto( &histoZY ); 
	    VT_ErrorParse("unable to fit gaussian\n", 0);
	}
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," gaussienne : m = %f k = %f s = %f", (float)M, (float)K, (float)S );
	    VT_Message( message, program );
	    sprintf( message," seuil = %f", (float)(M * mXY.max / 255.0 ) );
	    VT_Message( message, program );
	}


	VT_FreeImage( &auxXY );     VT_FreeImage( &auxXZ );      VT_FreeImage( &auxZY );
	VT_FreeImage( imXY );      VT_FreeImage( imXZ );       VT_FreeImage( imZY );
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
    int eroc=0;
    int dilc=0;
    char text[STRINGLENGTH];
    
    if ( VT_CopyName( program, argv[0] ) != 1 )
	VT_Error("Error while copying program name", (char*)NULL);
    if ( argc == 1 ) VT_ErrorParse("\n", 0 );
    
    /*--- initialisation des parametres ---*/
    VT_InitParam( par );

    /*--- les 3 premiers arguments sont les 3 images MIP ---*/
    strcpy( par->name_mip_XY, argv[1] );
    strcpy( par->name_mip_XZ, argv[2] );
    strcpy( par->name_mip_ZY, argv[3] );
    
    /*--- lecture des parametres ---*/
    i = 4; nb = 0;
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
	    /*--- erosion ---*/
            else if ( strcmp ( argv[i], "-eroi" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -eroi...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->spar.ero_iterations) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -eroi...\n", 0 );
            }
            else if ( strcmp ( argv[i], "-eroc" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -eroc...\n", 0 );
                status = sscanf( argv[i],"%d",&eroc );
                if ( status <= 0 ) VT_ErrorParse( "parsing -eroc...\n", 0 );
            }
	    /*--- coefficient du seuillage ---*/
	    else if ( strcmp ( argv[i], "-csb" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -csb...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->spar.mul_max_sb) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -csb...\n", 0 );
            }
	    else if ( strcmp ( argv[i], "-csh" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -csh...\n", 0 );
                status = sscanf( argv[i],"%f",&(par->spar.mul_max_sh) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -csh...\n", 0 );
            }            
	    /*--- dilatation ---*/
            else if ( strcmp ( argv[i], "-dili" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -dili...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->spar.dil_iterations) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -dili...\n", 0 );
            }
            else if ( strcmp ( argv[i], "-dilc" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -dilc...\n", 0 );
                status = sscanf( argv[i],"%d",&dilc );
                if ( status <= 0 ) VT_ErrorParse( "parsing -dilc...\n", 0 );
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
	strcpy( par->names.out, ">" );  /* standart output */
    }
    /*--- type de connexite pour l'erosion ---*/
    switch ( eroc ) {
    case 4 :
        par->spar.ero_connexity = N04;   break;
    case 6 :
        par->spar.ero_connexity = N06;   break;
    case 8 :
        par->spar.ero_connexity = N08;   break;
    case 10 :
        par->spar.ero_connexity = N10;   break;
    case 18 :
        par->spar.ero_connexity = N18;   break;
    case 26 :
        par->spar.ero_connexity = N26;   break;
    }

    if ( ( par->spar.mul_max_sh > 1.0 ) || ( par->spar.mul_max_sb >= 1.0 ) ||
	 ( par->spar.mul_max_sb > par->spar.mul_max_sh ) )
	VT_ErrorParse("bad coefficients for thresholding\n", 0 );

    /*--- type de connexite pour la dilatation ---*/
    switch ( dilc ) {
    case 4 :
        par->spar.dil_connexity = N04;   break;
    case 6 :
        par->spar.dil_connexity = N06;   break;
    case 8 :
        par->spar.dil_connexity = N08;   break;
    case 10 :
        par->spar.dil_connexity = N10;   break;
    case 18 :
        par->spar.dil_connexity = N18;   break;
    case 26 :
        par->spar.dil_connexity = N26;   break;
    }

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
