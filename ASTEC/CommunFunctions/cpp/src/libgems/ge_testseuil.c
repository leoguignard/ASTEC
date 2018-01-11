#include <vt_common.h>

#include <vt_mip.h>
#include <vt_extract.h>

#include <vt_histo.h>

#include <vt_recfilters.h>

#include <vt_geseuil1.h>
#include <vt_gaussienne.h>

#include <epidaure.h>
#include <epidaure.ee>
#include <epidaure.ei>

#define _WRITE_HISTO_ 1

typedef struct local_par {
  int seuil;
  int size;
  int dim;
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

static char *usage = "image-in image-out [-seuil %d] [-2D | -3D] [-size %d]\n\
\t [-eroc %d] [-eroi %d] [-csb %f] [-csh %f] [-dilc %d] [-dili %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

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
	vt_image *first_image;
	vt_image sous_image, imres;
	vt_image mip[3];
	int i, v, x, y, z;
	int half_size, size;
	u16 *usbuf;
	u16 ***theBuf, ***resBuf;
	vt_ipt corner;
	vt_histo histo3D, histo2D[3];
	vt_image hres3D, hres2D[3];
	vt_image hflo3D, hflo2D[3];
	int filled3D, filled2D[3];
	int calc3D, calc2D[3];
	double M3, K3, S3, M2[3], K2[3], S2[3];
	i32 *ibuf;
	vt_3m m;

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- 1ere image ---*/
	first_image = _VT_Inrimage( par.names.in );
	
	/*--- sous-images ---*/
	half_size = par.size / 2;
	size = 2 * half_size + 1;
	VT_Image( &sous_image );
	VT_InitImage( &sous_image, "", size, size, size, first_image->type );
	(void)VT_AllocImage( &sous_image );
	for ( i = 0; i < 3; i ++ ) {
	  VT_Image( &mip[i] );
	  VT_InitImage( &mip[i], "", size, size, (int)1, first_image->type );
	  (void)VT_AllocImage( &mip[i] );
	}
	VT_Image( &imres );
	VT_InitFromImage( &imres, first_image, par.names.out, first_image->type );
	(void)VT_AllocImage( &imres );
	v = first_image->dim.x * first_image->dim.y * first_image->dim.z;
	usbuf = (u16*)imres.buf;
	for ( i = 0; i < v; i ++ ) 
	  *usbuf++ = (u16)0;

	/*--- histogrammes ----*/
	VT_Histo( &histo3D );   histo3D.type = first_image->type;
	(void)VT_AllocHisto( &histo3D ); 
	for ( i = 0; i < 3; i ++ ) {
	  VT_Histo( &histo2D[i] );   histo2D[i].type = first_image->type;
	  (void)VT_AllocHisto( &histo2D[i] ); 
	}

	/*--- images auxiliaires ---*/
	VT_Image( &hres3D );
	VT_InitImage( &hres3D, "", histo3D.size, 1, 1, INT );
	(void)VT_AllocImage( &hres3D );
	VT_Image( &hflo3D );
	VT_InitImage( &hflo3D, "", histo3D.size, 1, 1, FLOAT );
	(void)VT_AllocImage( &hflo3D );
	for ( i = 0; i < 3; i ++ ) {
	  VT_Image( &hres2D[i] );
	  VT_InitImage( &hres2D[i], "", histo2D[i].size, 1, 1, INT );
	  (void)VT_AllocImage( &hres2D[i] );
	  VT_Image( &hflo2D[i] );
	  VT_InitImage( &hflo2D[i], "", histo2D[i].size, 1, 1, FLOAT );
	  (void)VT_AllocImage( &hflo2D[i] );
	}

	/*--- traitement ---*/
	theBuf = (u16***)first_image->array;
	resBuf = (u16***)imres.array;	
	for ( z = half_size; z < first_image->dim.z - half_size; z ++ ) {
	  fprintf( stderr, "\r plan %3d", z );
	  for ( y = half_size; y < first_image->dim.y - half_size; y ++ ) {
	    for ( x = half_size; x < first_image->dim.x - half_size; x ++ ) {
	      if ( theBuf[z][y][x] < (u16)par.seuil ) continue;
	      
	      /*--- calcul des sous images ---*/
	      corner.x = x - half_size;
	      corner.y = y - half_size;
	      corner.z = z - half_size;
	      (void)VT_Extract( &sous_image, first_image, &corner );
	      
	      if ( par.dim == VT_3D ) {
		/*--- initialisation ---*/
		filled3D = calc3D = 0;
		
		/*--- calcul histogramme ---*/
		(void)VT_3m( &sous_image, &m );
		if ( (int)m.min == 0 ) {
		  if ( VT_Image2Histo( &histo3D, &sous_image, &(par.spar) ) == 1 ) {
		    filled3D = 1;
		    /*--- copie ---*/
		    ibuf = (i32*)hres3D.buf;
		    for ( v = 0; v < histo3D.size; v ++ ) ibuf[v] = histo3D.buf[v];
		    /*--- pas de point a 0 ---*/
		    ibuf[0] = 0;
		    /*--- pas de point isole ---*/
		    for ( v = 1; v < histo3D.size - 1; v++ ) {
		      if ( ibuf[v] == 0 ) continue;
		      if ( ( ibuf[v-1] == 0 ) && ( ibuf[v+1] == 0 ) ) ibuf[v] = 0;
		    }
		    /*--- lissage de l'histogramme ---*/
		    if ( VT_RecFilterOnImage( &hres3D, &hflo3D, &(par.spar.par_filt) ) == 1 )
		      calc3D = VT_FitGaussienne( (r32*)hflo3D.buf, histo3D.size, &K3, &S3, &M3 );
		    
		    if ( calc3D == 1 )
		      resBuf[z][y][x] = (u16)( M3 + 0.5 );
		  }
		}
		
	      } else {
		/*--- cas 2D ---*/
		(void)VT_MIP( &sous_image, &mip[0], &mip[1], &mip[2] );
		/*--- initialisation ---*/
		for ( i = 0; i < 3; i ++ ) {
		  filled2D[i] = calc2D[i] = 0;
		}
		
		/*--- calcul des histogrammes ---*/
		for ( i = 0; i < 3; i ++ ) {
		  (void)VT_3m( &mip[i], &m );
		  if ( (int)m.min == 0 ) {
		    if ( VT_Image2Histo( &histo2D[i], &mip[i], &(par.spar) ) == 1 ) {
		      filled2D[i] = 1;
		      /*--- copie ---*/
		      ibuf = (i32*)hres2D[i].buf;
		      for ( v = 0; v < histo2D[i].size; v ++ ) ibuf[v] = histo2D[i].buf[v];
		      /*--- pas de point a 0 ---*/
		      ibuf[0] = 0;
		      /*--- pas de point isole ---*/
		      for ( v = 1; v < histo2D[i].size - 1; v++ ) {
			if ( ibuf[v] == 0 ) continue;
			if ( ( ibuf[v-1] == 0 ) && ( ibuf[v+1] == 0 ) ) ibuf[v] = 0;
		      }
		      /*--- lissage de l'histogramme ---*/
		      if ( VT_RecFilterOnImage( &hres2D[i], &hflo2D[i], &(par.spar.par_filt) ) == 1 )
			calc2D[i] = VT_FitGaussienne( (r32*)hflo2D[i].buf, histo2D[i].size, &K2[i], &S2[i], &M2[i] );
		    }
		  }
		}
		
		/*--- calcul du seuil ---*/
		if ( (calc2D[0] != 1) || (calc2D[1] != 1) || (calc2D[2] != 1) ) {
		  fprintf( stderr, "calcul rate pour point (%d,%d,%d) : ok pour", x, y, z );
		  for ( i = 0; i < 3; i ++ ) 
		    if ( calc2D[i] == 1 ) fprintf( stderr, " %d", i );
		  fprintf( stderr, "\n" );
		} else {
		  /*--- calcul de la valeur mediane --*/
		  if ( M2[1] >= M2[0] ) {
		    if ( M2[2] >= M2[1] )      M3 = M2[1];
		    else {
		      if ( M2[2] >= M2[0] )    M3 = M2[2];
		      else                     M3 = M2[0];
		    }
		  } else {
		    if ( M2[2] >= M2[0] )      M3 = M2[0];
		    else {
		      if ( M2[2] >= M2[1] )    M3 = M2[1];
		      else                     M3 = M2[2];
		    }
		  }
		  resBuf[z][y][x] = (u16)( M3 + 0.5 );
		}
		
	      }
	     
	    }
	  }
	}
	fprintf( stderr, "\n" );
	(void)VT_WriteInrimage( &imres );
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
	    /*--- dimension ---*/
	    else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
		par->dim = VT_2D;
	    }
	    else if ( strcmp ( argv[i], "-3D" ) == 0 ) {
		par->dim = VT_3D;
	    }
	    /*--- seuil ---*/
            else if ( strcmp ( argv[i], "-seuil" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -seuil...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->seuil) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -seuil...\n", 0 );
            }
	    /*--- size ---*/
            else if ( strcmp ( argv[i], "-size" ) == 0 ) {
                i += 1;
                if ( i >= argc)    VT_ErrorParse( "parsing -size...\n", 0 );
                status = sscanf( argv[i],"%d",&(par->size) );
                if ( status <= 0 ) VT_ErrorParse( "parsing -size...\n", 0 );
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
  par->size = 64;
  par->seuil = 500;
  par->dim = VT_2D;
  par->nb_names = 0;
  VT_Names( &(par->names) );
  VT_Seuil1( &(par->spar) );
}
