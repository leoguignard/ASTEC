#include <vt_common.h>

#include <vt_histo.h>

#include <vt_recfilters.h>

#include <vt_geseuil1.h>
#include <vt_gaussienne.h>

#include <epidaure.h>
#include <epidaure.ee>
#include <epidaure.ei>

int _WRITE_HISTO_=0;

typedef struct local_par {
    int nb_names;
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

static char *usage = "image-1 [image-2] [image-3]\n\
\t [-eroc %d] [-eroi %d] [-csb %f] [-csh %f] [-dilc %d] [-dili %d]\n\
\t [-histo] [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t -eroc %d : connexite pour l'erosion de l'image MIP (defaut = 26)\n\
\t -eroi %d : nombre d'iterations pour l'erosion de l'image MIP (defaut = 1)\n\
\t -csb %f  : coefficient multiplicateur du max pour le calcul\n\
\t            du seuil bas du seuillage par hysteresis (defaut = 0.1666)\n\
\t -csh %f  : coefficient multiplicateur du max pour le calcul\n\
\t            du seuil haut du seuillage par hysteresis (defaut = 0.3333)\n\
\t -dilc %d : connexite pour la dilatation de l'image des contours (defaut = 6)\n\
\t -dili %d : nombre d'iterations pour la dilatation de l'image des contours (defaut = 1)\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n";

/* static char program[STRINGLENGTH]; */

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
	vt_image *image[3];
	vt_histo histo[3];
	vt_image hres[3], hflo[3];
	vt_image imres, imfloat;
	vt_3m m[3];
	i32 *buf, *ibuf;
	int i, n, na, filled_histo, filled_histos[3];
	double M, K, S;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- lecture des images d'entree ---*/
	for ( n = 0; n < par.nb_names; n ++ ) {
	    switch ( n ) {
	    case 0 :
		image[0] = _VT_Inrimage( par.names.in );
		if ( image[0] == (vt_image*)NULL ) 
		    VT_ErrorParse("unable to read first input image\n", 0);
		break;
	    case 1 :
		image[1] = _VT_Inrimage( par.names.out );
		if ( image[1] == (vt_image*)NULL ) {
		    VT_FreeImage( image[0] );    VT_Free( (void**)&image[0] );
		    VT_ErrorParse("unable to read second input image\n", 0);
		}
		break;
	    case 2 :
		image[2] = _VT_Inrimage( par.names.ext );
		if ( image[2] == (vt_image*)NULL ) {
		    VT_FreeImage( image[0] );    VT_Free( (void**)&image[0] );
		    VT_FreeImage( image[1] );    VT_Free( (void**)&image[1] );
		    VT_ErrorParse("unable to read third input image\n", 0);
		}
		break;
	    }
	}

	/*--- test de types ---*/
 	for ( n = 1; n < par.nb_names; n ++ ) 
	    if ( image[n]->type != image[0]->type ) {
		for ( na = 0; na < par.nb_names; na ++ ) {
		    VT_FreeImage( image[na] );    
		    VT_Free( (void**)&image[na] );
		}
		VT_ErrorParse("images have different types\n", 0);
	    }
	
	/*--- operations eventuelles sur l'image d'entree ---*/
	for ( n = 0; n < par.nb_names; n ++ ) {
	    if ( par.names.inv == 1 )  VT_InverseImage( image[n] );
	    if ( par.names.swap == 1 ) VT_SwapImage( image[n] );
	}

	/*--- initialisation histogramme ---*/
	for ( n = 0; n < par.nb_names; n ++ ) {
	    VT_Histo( &(histo[n]) );   histo[n].type = image[n]->type;
	    if ( VT_AllocHisto( &(histo[n]) ) != 1 ) {
		for ( na = 0; na < par.nb_names; na ++ ) {
		    VT_FreeImage( image[na] );    
		    VT_Free( (void**)&image[na] );
		}
		for ( na = 0; na < n; na ++ ) 
		    VT_FreeHisto( &(histo[na]) );
		VT_ErrorParse("unable to allocate auxiliary histograms\n", 0);
	    }
	}

	/*--- initialisation outputs ---*/
	for ( n = 0; n < par.nb_names; n ++ ) {

	  VT_Image( &hres[n] );
	  VT_InitImage( &hres[n], "", histo[0].size, 1, 1, INT );
	  if ( VT_AllocImage( &hres[n] ) != 1 ) {
	    VT_ErrorParse("unable to allocate auxiliary integer images\n", 0);
	  }

	  VT_Image( &hflo[n] );
	  VT_InitImage( &hflo[n], "", histo[0].size, 1, 1, FLOAT );
	  if ( VT_AllocImage( &hflo[n] ) != 1 ) {
	    VT_ErrorParse("unable to allocate auxiliary float images\n", 0);
	  }
	}

	/*--- initialisation output ---*/
	VT_Image( &imres );
	VT_Strncpy( image_name, program, STRINGLENGTH );
	(void)strcat( image_name, ".histo" );
	VT_InitImage( &imres, image_name, histo[0].size, 1, 1, INT );
	if ( VT_AllocImage( &imres ) != 1 ) {
	    for ( na = 0; na < par.nb_names; na ++ ) {
		VT_FreeImage( image[na] );    
		VT_Free( (void**)&image[na] );
		VT_FreeHisto( &(histo[n]) );
	    }
	    VT_ErrorParse("unable to allocate first auxiliary image\n", 0);
	}
	VT_Image( &imfloat );
	VT_InitImage( &imfloat, image_name, histo[0].size, 1, 1, FLOAT );
	if ( VT_AllocImage( &imfloat ) != 1 ) {
	    for ( na = 0; na < par.nb_names; na ++ ) {
		VT_FreeImage( image[na] );    
		VT_Free( (void**)&image[na] );
		VT_FreeHisto( &(histo[n]) );
	    }
	    VT_FreeImage( &imres );
	    VT_ErrorParse("unable to allocate second auxiliary image\n", 0);
	}

	/*--- on y va ---*/
	filled_histo = 0;
	for ( i = 0; i < 3; i ++ ) filled_histos[i] = 0;
	buf = (i32*)imres.buf;
	for ( i = 0; i < histo[0].size; i++ )
	    buf[i] = 0;
	for ( n = 0; n < par.nb_names; n ++ ) {
	  ibuf = (i32*)hres[n].buf;
	  for ( i = 0; i < histo[0].size; i++ )
	    ibuf[i] = 0;
	}

	/*--- calcul histogramme ---*/
	for ( n = 0; n < par.nb_names; n ++ ) {
	    for ( i = 0; i < histo[n].size; i++ )
		histo[n].buf[i] = 0;
	    
	    (void)VT_3m( image[n], &(m[n]) );
	    
	    /* pourquoi ce test ? */
	    /* if ( (int)(m[n].min) == 0 ) */ {
	        i32 *ibuf;
		/*--- calcul histogramme ---*/
		if ( VT_Image2Histo( &histo[n], image[n], &(par.spar) ) != 1 ) {
		    for ( na = 0; na < par.nb_names; na ++ ) {
			VT_FreeImage( image[na] );    
			VT_Free( (void**)&image[na] );
			VT_FreeHisto( &(histo[n]) );
		    }
		    VT_FreeImage( &imres );
		    VT_FreeImage( &imfloat );
		    {
		      char text[STRINGLENGTH];
		      switch ( n ) {
		      case 0 :
			sprintf(text,"unable to fill histogram from image %s\n", par.names.in );
			break;
		      case 1 :
			sprintf(text,"unable to fill histogram from image %s\n", par.names.out );
			break;
		      case 2 :
			sprintf(text,"unable to fill histogram from image %s\n", par.names.ext );
			break;
		      default :
			  sprintf(text,"unable to fill histogram from image ?\n");
		      }
		      VT_ErrorParse(text, 0);
		    }
		}
		
		filled_histo = filled_histos[n] = 1;
		ibuf = (i32 *)hres[n].buf;
		/*--- copie buffer ---*/
		for ( i = 0; i < histo[n].size; i++ ) {
		    buf[i] += histo[n].buf[i];
		    ibuf[i] = histo[n].buf[i];
		}
	    }
	    
	}

	/*--- pas de point a 0 ---*/
	buf[0] = 0;
	/*--- pas de point isole ---*/
	for ( i = 1; i < histo[0].size - 1; i++ ) {
	    if ( buf[i] == 0 ) continue;
	    if ( ( buf[i-1] == 0 ) && ( buf[i+1] == 0 ) ) buf[i] = 0;
	}

	for ( n = 0; n < par.nb_names; n ++ ) {
	  ibuf = (i32 *)hres[n].buf;
	  /*--- pas de point a 0 ---*/
	  ibuf[0] = 0;
	  /*--- pas de point isole ---*/
	  for ( i = 1; i < histo[0].size - 1; i++ ) {
	    if ( ibuf[i] == 0 ) continue;
	    if ( ( ibuf[i-1] == 0 ) && ( ibuf[i+1] == 0 ) ) ibuf[i] = 0;
	  }
	}

	/*--- ecriture de l'histogramme ---*/
	if ( _WRITE_HISTO_ == 1 ) {
	  (void)VT_WriteInrimage( &imres );
	}

	if ( filled_histo == 1 ) {
	    /*--- lissage de l'histogramme ---*/
	    if ( VT_RecFilterOnImage( &imres, &imfloat, &(par.spar.par_filt) ) != 1 ) {
		for ( na = 0; na < par.nb_names; na ++ ) {
		    VT_FreeImage( image[na] );    
		    VT_Free( (void**)&image[na] );
		    VT_FreeHisto( &(histo[n]) );
		    }
		VT_FreeImage( &imres );
		VT_FreeImage( &imfloat );
		VT_ErrorParse("unable to smooth output image\n", 0);
	    }

	    /*--- on ajuste une gaussienne ---*/
	    /* M = 1000.0;  K = 1000.0;  S = 100.0; */
	    if ( VT_FitGaussienne( (r32*)imfloat.buf, histo[0].size, &K, &S, &M ) != 1 ) {
		for ( na = 0; na < par.nb_names; na ++ ) {
		    VT_FreeImage( image[na] );    
		    VT_Free( (void**)&image[na] );
		    VT_FreeHisto( &(histo[n]) );
		}
		VT_FreeImage( &imres );
		VT_FreeImage( &imfloat );
		VT_ErrorParse("unable to fit gaussian\n", 0);
	    }
		
	    {
		char message[256];
		sprintf( message," gaussienne : m = %f k = %f s = %f", (float)M, (float)K, (float)S );
		VT_Message( message, program );
		sprintf( message," seuil = %f", (float)M );
		/* VT_Message( message, program ); */
		if ( par.nb_names > 1 ) {
		    fprintf( stdout, " Pour les images ");
		} else {
		    fprintf( stdout, " Pour l'image ");
		}
		for ( na = 0; na < par.nb_names; na ++ ) 
		    fprintf( stdout, "%s ", image[na]->name );
		fprintf( stdout, "\n le seuil calcule est %f (sigma = %f).\n", (float)M, (float)S );
		/* (void)VT_WriteInrimage( &imres ); */
	    }
	}

	if ( par.nb_names > 1 ) {
	  for ( n = 0; n < par.nb_names; n ++ ) {
	    if ( filled_histos[n] == 1 ) {
	      if ( VT_RecFilterOnImage( &hres[n], &hflo[n], &(par.spar.par_filt) ) != 1 ) {
		VT_ErrorParse("unable to smooth float output image\n", 0);
	      }
	      /*--- on ajuste une gaussienne ---*/
	      if ( VT_FitGaussienne( (r32*)hflo[n].buf, histo[0].size, &K, &S, &M ) != 1 ) {
		VT_ErrorParse("unable to fit gaussian\n", 0);
	      }
	      
	      {
		char message[256];
		sprintf( message," gaussienne : m = %f k = %f s = %f", (float)M, (float)K, (float)S );
		VT_Message( message, program );
		sprintf( message," seuil = %f", (float)M );
		/* VT_Message( message, program ); */
		fprintf( stdout, " Pour l'image ");
		fprintf( stdout, "%s ", image[n]->name );
		fprintf( stdout, "\n le seuil calcule est %f (sigma = %f).\n", (float)M, (float)S );
		/* (void)VT_WriteInrimage( &imres ); */
	      }
	    }
	  }
	}

	/*--- fin des tests ---*/
	VT_FreeImage( &imres );
	VT_FreeImage( &imfloat );
	for ( na = 0; na < par.nb_names; na ++ ) {
	    VT_FreeImage( image[na] );    
	    VT_Free( (void**)&image[na] );
	    VT_FreeHisto( &(histo[na]) );
	}

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
    
    /*--- initialisation des parametres ---*/
    VT_InitParam( par );

    (void*)strcpy( program, argv[0] );
    if ( argc == 1 ) VT_ErrorParse("\n", 0 );
    
    /*--- lecture des parametres ---*/
    i = 1; nb = 0;
    while ( i < argc ) {
	if ( argv[i][0] == '-' ) {
	    /*--- arguments generaux ---*/
	    if ( strcmp ( argv[i], "-help" ) == 0 ) {
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
	    /* */
	    else if ( strcmp ( argv[i], "-histo" ) == 0 ) {
	      _WRITE_HISTO_ = 1;
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
	    else if ( nb == 2 ) {
		strncpy( par->names.ext, argv[i], STRINGLENGTH );  
		nb += 1;
	    }
	    else 
		VT_ErrorParse("too much file names when parsing\n", 0 );
	}
	i += 1;
    }
    
    /*--- s'il n'y a pas assez de noms ... ---*/
    if (nb == 0) 
	VT_ErrorParse("no file name when parsing\n", 0 );
    par->nb_names = nb;

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
    par->nb_names = 0;
    VT_Names( &(par->names) );
    VT_Seuil1( &(par->spar) );
}
