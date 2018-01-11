#include <vt_common.h>
#include <vt_neighborhood.h>

#include <vt_distance.h>

#include <vt_amincir.h>
#include <vt_amseuil.h>
#include <vt_amliste.h>
#include <vt_geambdd.h>
#include <vt_geremove.h>
#include <vt_gemask.h>

/*
  1. on seuille l'image originale 
     INPUT = image ; OUTPUT = imres
     on peut liberer image
  2. on construit une carte de distance
     INPUT = imres ; OUTPUT = imdist
  3. on construit une liste de points
     INPUT = imres, imdist ; OUTPUT = liste
     on peut liberer imdist
  4. on amincit 
     INPUT = imres, liste ; OUTPUT = imres, liste

  5. on supprime les composantes en trop
  6. on restreint l'espace de reconstruction
  7. on reconstruit

     
 */

typedef struct local_par {
    float sb;
    float sh;
    int steps;
    int type_distance;
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

static char *usage = "[image-in] [image-out] [-sb %f] [-sh %f] [-steps %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help]";
static char *detail = "\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -sb %f : seuil bas\n\
\t -sh %f : seuil haut \n\
\t Les points dont la valeur est superieure ou egale au seuil bas, et,\n\
\t si un seuil haut est donne, inferieure ou egale au seuil haut, sont mis\n\
\t a 255, 32765, 65535, 2147483647 ou 1.0; les autres a 0.0 (defaut);\n\
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
	vt_image *image, imres, imdist;
	int steps = 0;
	char message[256], *proc="main";
	/* distance */
	vt_distance dpar;
	/* amincissement */
	vt_vpt_amincir *liste = (vt_vpt_amincir*)NULL;
	int nb=0, nb_ori=0;
	vt_amincir apar;

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

	/*--- seuillage pour l'amincissement :
	  1. val >= seuil haut 
             => label = VT_UNDELETABLE
	  2. val >= seuil bas && val < seuil haut 
             => label = VT_DELETABLE
	  3. val < seuil bas
	     => label = 0
	  ---*/
	sprintf( message, "# %3d # seuillage de l'image originale.", ++steps );
	VT_Message( message, proc );
	if ( _VT_Thin2Thresholds( &imres, image, par.sb, par.sh ) != 1 ) {
	    VT_FreeImage( image );
	    VT_FreeImage( &imres );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("erreur pendant le seuillage pour l'amincissement\n", 0);
	}

	/*--- liberation de l'input ---*/
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	    
	par.steps --;
	if ( par.steps <= 0 ) {
	    if ( VT_WriteInrimage( &imres ) == -1 ) {
		VT_FreeImage( &imres );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imres );
	    exit( 0 );
	}

	/*--- initialisation de l'image distance ---*/
	VT_InitImage( &imdist, names.out, imres.dim.x, imres.dim.y, imres.dim.z, (int)USHORT );
	if ( VT_AllocImage( &imdist ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_ErrorParse("unable to allocate distance image\n", 0);
	}
	
	/*--- calcul de la distance ---*/
	sprintf( message, "# %3d # calcul de la distance.", ++steps );
	VT_Message( message, proc );
	VT_Distance( &dpar );
	dpar.type = VT_DIST_CHMFR2;
	dpar.seuil = 100.0; /* les valeurs sont >= 200 */
	{
	    int i;
	    u8 *in;
	    u16 *out;
	    in = (u8*)imres.buf;
	    out = (u16*)imdist.buf;
	    for ( i = imres.dim.x * imres.dim.y * imres.dim.z; i > 0; i-- )
		*out++ = ( *in++ > (u8)100 ) ? (u16)0 : (u16)255;
	}
	if ( VT_Dist( &imdist, &imdist, &dpar ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_FreeImage( &imdist );
	    VT_ErrorParse("unable to compute distance image\n", 0);
	}
	
	par.steps --;
	if ( par.steps <= 0 ) {
	    VT_FreeImage( &imres );
	    if ( VT_WriteInrimage( &imdist ) == -1 ) {
		VT_FreeImage( &imdist );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imdist );
	    exit( 0 );
	}
	
	/*--- amincissement ---*/
	
	/*--- construction de la liste :
	  les points VT_DELETABLE sont mis dans la liste ---*/
	liste = _VT_ThinVPtList( &imres, &imdist, (int)VT_3D, &nb );
	nb_ori = nb;
	if ( (liste == (vt_vpt_amincir*)NULL) || (nb <= 0) ) {
	    VT_FreeImage( &imres );
	    VT_FreeImage( &imdist );
	    VT_Free( (void**)&liste );
	    VT_ErrorParse("impossible de construire une liste de points\n", 0);
	}
	
	/*--- liberation de l'image distance ---*/
	VT_FreeImage( &imdist );
	
	VT_Amincir( &apar );
	apar.bool_shrink = 1;
	apar.bool_end_surfaces = apar.bool_end_curves = 0;
	
	/*--- 1er amincissement : 6-connexite :
	  il reste 
	  VT_DELETABLE      (on aurait pu l'effacer mais on ne l'a pas fait)
	  VT_HASBEENDELETED (a effacer)
	  VT_UNDELETABLE    (rien a faire)
	  ---*/
	sprintf( message, "# %3d # premier amincissement.", ++steps );
	VT_Message( message, proc );
	apar.connexite = N06;
	apar.epaisseur = N06;
	if ( _VT_GEBDD_THINNING( &imres, liste, &nb, &apar ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_Free( (void**)&liste );
	    VT_ErrorParse("1er amincissement rate\n", 0);
	}
	
	par.steps --;
	if ( par.steps <= 0 ) {
	    VT_Free( (void**)&liste );
	    if ( VT_WriteInrimage( &imres ) == -1 ) {
		VT_FreeImage( &imres );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imres );
	    exit( 0 );
	}
	
	sprintf( message, "# %3d # second amincissement.", ++steps );
	VT_Message( message, proc );
	if ( nb > 0 ) {
	    /*--- 2nd amincissement : 26-connexite 
	      il reste 
	      VT_DELETABLE      (on aurait pu l'effacer mais on ne l'a pas fait)
	      VT_HASBEENDELETED (a effacer)
	      VT_UNDELETABLE    (rien a faire)
	      ---*/
	    apar.connexite = N26;
	    apar.epaisseur = N26;
	    if ( _VT_GEBDD_THINNING( &imres, liste, &nb, &apar ) != 1 ) {
		VT_FreeImage( &imres );
		VT_Free( (void**)&liste );
		VT_ErrorParse("2nd amincissement rate\n", 0);
	    }
	}
	
	par.steps --;
	if ( par.steps <= 0 ) {
	    if ( VT_WriteInrimage( &imres ) == -1 ) {
		VT_FreeImage( &imres );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imres );
	    exit( 0 );
	}
	
	/*
	{
	    int i;
	    int n_has, n_tobe, n_will, n_dele, n_unde, n;

	    n_has = n_tobe = n_will = n_dele = n_unde = n = 0;
	    for ( i = 0; i < nb_ori; i++ ) {
		switch ( liste[i].status ) {
		case VT_HASBEENDELETED :
		    n_has ++;
		    break;
		case VT_TOBEDELETED :
		    n_tobe ++;
		    break;
		case VT_WILLBEDELETED :
		    n_will ++;
		    break;
		case VT_DELETABLE :
		    n_dele ++;
		    break;
		case VT_UNDELETABLE :
		    n_unde ++;
		    break;
		default :
		    n ++;
		}
	    }
	    
	    printf(" sur %d points, il en reste %d\n", nb_ori, nb );
	    printf(" --- %d VT_HASBEENDELETED\n", n_has );
	    printf(" --- %d VT_TOBEDELETED\n", n_tobe );
	    printf(" --- %d VT_WILLBEDELETED\n", n_will );
	    printf(" --- %d VT_DELETABLE\n", n_dele );
	    printf(" --- %d VT_UNDELETABLE\n", n_unde );
	    printf(" --- %d unknown\n", n);
	}
	*/

	/*--- suppression des points isoles ---*/
	/*    Remarque :
              Il faudrait mieux conserver le resultat de 
	      la reconstruction des points "bas" par rapport
	      aux points "hauts" car il peut rester autre 
	      chose que des points isoles : des courbes fermees, etc. 
              On peut le faire sur la liste */
	/* dans l'image il doit y avoir :
	   VT_HASBEENDELETED : points deja effaces
	   VT_UNDELETABLE    : points a conserver (seuil haut)
	   VT_DELETABLE      : points effacables
	   dans la liste il doit y avoir :
	   VT_TOBEDELETED : points effaces
	   VT_DELETABLE   : points effacables (en debut de liste)
	   */
	sprintf( message, "# %3d # suppression des points en trop.", ++steps );
	VT_Message( message, proc );
	if ( _VT_GEREMOVE( &imres, liste, &nb ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_Free( (void**)&liste );	
	    VT_ErrorParse("suppression des points isoles ratee\n", 0);
	}

	par.steps --;
	if ( par.steps <= 0 ) {
	    VT_Free( (void**)&liste );	
	    if ( VT_WriteInrimage( &imres ) == -1 ) {
		VT_FreeImage( &imres );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imres );
	    exit( 0 );
	}

	/*--- contraindre la reconstruction ---*/
	/* dans l'image il doit y avoir :
	   VT_HASBEENDELETED : points deja effaces
	   VT_UNDELETABLE    : points a conserver (seuil haut)
	   VT_U_DELETABLE    : points a conserver (seuil bas mais topologie)
	   dans la liste il doit y avoir :
	   VT_U_DELETABLE : points a conserver (seuil bas mais topologie) (en debut de liste)
	   VT_TOBEDELETED : points effaces
	   */
	sprintf( message, "# %3d # contrainte de la reconstruction.", ++steps );
	VT_Message( message, proc );
	if ( _VT_GEMASK( &imres, liste, nb, (int)VT_DIST_CHMFR2 ) != 1 ) {
	    VT_FreeImage( &imres );
	    VT_Free( (void**)&liste );	
	    VT_ErrorParse("contrainte de la reconstruction ratee\n", 0);
	}

	par.steps --;
	if ( par.steps <= 0 ) {
	    if ( VT_WriteInrimage( &imres ) == -1 ) {
		VT_FreeImage( &imres );
		VT_ErrorParse("unable to write output image\n", 0);
	    }
	    VT_FreeImage( &imres );
	    exit( 0 );
	}

	/*--- on libere la liste ---*/
	VT_Free( (void**)&liste );	

	/*--- reconstruction ---*/
	/* dans l'image il doit y avoir maintenant :
	   VT_HASBEENDELETED  : points deja effaces
	   VT_RECONSTRUCTABLE : points pouvant etre reconstruits
	   VT_UNDELETABLE     : points a conserver (seuil haut)
	   VT_U_DELETABLE     : points a conserver (seuil bas mais topologie)
	   dans la liste il doit y avoir maintenant :
	   VT_U_DELETABLE : points a conserver (seuil bas mais topologie) (en debut de liste)
	   VT_TOBEDELETED : points effaces
	   */
	
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
		else if ( strcmp ( argv[i], "-sh" ) == 0 ) {
		    i += 1;
		    if ( i >= argc)    VT_ErrorParse( "parsing -sh...\n", 0 );
		    status = sscanf( argv[i],"%f",&(par->sh) );
		    if ( status <= 0 ) VT_ErrorParse( "parsing -sh...\n", 0 );
		}
		/*--- nombre de pas ---*/
		else if ( strcmp ( argv[i], "-steps" ) == 0 ) {
		    i += 1;
		    if ( i >= argc)    VT_ErrorParse( "parsing -steps...\n", 0 );
		    status = sscanf( argv[i],"%d",&(par->steps) );
		    if ( status <= 0 ) VT_ErrorParse( "parsing -steps...\n", 0 );
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
	par->sb = (float)0.0;
	par->sh = (float)0.0;
	par->steps = 100;
	par->type_distance = VT_DIST_CHMFR2;
}
