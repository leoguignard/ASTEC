#include <vt_common.h>
#include <vt_interpol.h>

typedef struct local_par {
  vt_names names;
  int nb;
  int indice;
  int writeOriginaux;
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

static char *usage = "[image-1] [image-2] [image-out]\n\
\t [-nb %d] [-2D] [-i %d] [-noext]\n\
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

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image1, *image2, imres;
  vt_image imDist1, imDist2;
  typeDistanceMap theDist1, theDist2;
  unsigned char ***resBuf;
  int i;
  int ifirst, ilast;
  int dimzres;
  double a, b;


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture des images d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);


  /*--- lecture des images de sortie ---*/
  VT_Image( &imDist1 );
  VT_InitImage( &imDist1, "distance-1.inr", image1->dim.x, image1->dim.y, 256, SSHORT );
  if ( VT_AllocImage( &imDist1 ) != 1 ) {
    VT_FreeImage( image1 );
    VT_Free( (void**)&image1 );
    VT_FreeImage( image2 );
    VT_Free( (void**)&image2 );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  VT_Image( &imDist2 );
  VT_InitImage( &imDist2, "distance-2.inr", image2->dim.x, image2->dim.y, 256, SSHORT );
  if ( VT_AllocImage( &imDist2 ) != 1 ) {
    VT_FreeImage( image1 );
    VT_Free( (void**)&image1 );
    VT_FreeImage( image2 );
    VT_Free( (void**)&image2 );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }


  /*--- calcul des distances ---*/

  theDist1.buf = imDist1.buf;
  theDist1.dim[0] = imDist1.dim.x;
  theDist1.dim[1] = imDist1.dim.y;
  theDist1.dim[2] = imDist1.dim.z;
  theDist1.voxelSize[0] = 1.0;
  theDist1.voxelSize[1] = 1.0;
  theDist1.voxelSize[2] = 1.0;
  theDist1.multiplicativeCoefficient = 0.0;
  theDist1.intensityMax = VT_InitDistanceImageFromSlice( image1, &imDist1, 0 );


  theDist2.buf = imDist2.buf;
  theDist2.dim[0] = imDist2.dim.x;
  theDist2.dim[1] = imDist2.dim.y;
  theDist2.dim[2] = imDist2.dim.z;
  theDist2.voxelSize[0] = 1.0;
  theDist2.voxelSize[1] = 1.0;
  theDist2.voxelSize[2] = 1.0;
  theDist2.multiplicativeCoefficient = 0.0;
  theDist2.intensityMax = VT_InitDistanceImageFromSlice( image2, &imDist2, 0 );



  /*-- calcul de la carte
   */
  _DistanceSetNoVerbose();
  _ComputeSignedDistanceMap( &theDist1, -1.0 );
  _ComputeSignedDistanceMap( &theDist2, -1.0 );
  _CombineTwoDistanceMaps2D( &theDist1, &theDist2 );

  if ( 0 ) {
    VT_WriteInrimage( &imDist1 );
    VT_WriteInrimage( &imDist2 );
  }
  /*--- initialisation de l'image resultat ---*/
  ifirst = 1;
  ilast = par.nb;
  if ( par.indice >= 1 && par.indice <= par.nb ) {
    ifirst = ilast = par.indice;
    par.writeOriginaux = 0;
  }
  dimzres = ilast - ifirst + 1;
  if ( par.indice < 1 || par.indice > par.nb ) {
    if ( par.writeOriginaux != 0 ) dimzres += 2;
  }

  
  VT_Image( &imres );
  VT_InitImage( &imres, par.names.out, image1->dim.x, image1->dim.y, 
		dimzres, image1->type );

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  resBuf = (unsigned char ***)imres.array;
  if ( par.writeOriginaux != 0 ) {
    (void)memcpy( (void*)&(resBuf[0][0][0]), image1->buf, 
		  image1->dim.x*image1->dim.y );
    (void)memcpy( (void*)&(resBuf[par.nb+1][0][0]), image2->buf, 
		  image2->dim.x*image2->dim.y );
  }
  
  VT_FreeImage( image1 );
  VT_Free( (void**)&image1 );
  VT_FreeImage( image2 );
  VT_Free( (void**)&image2 );



  for ( i=ifirst; i<=ilast; i++ ) {
    a = (double)(par.nb+1-i)/(double)(par.nb+1);
    b = (double)(i)/(double)(par.nb+1);
    fprintf( stderr, " creating image %3d \n", i );
    if ( par.indice >= 1 && par.indice <= par.nb ) {
      VT_ComputeSliceFromDistances( &imDist1, &imDist2,
				    &imres, 0, a, b );
    } else {
      if ( par.writeOriginaux != 0 ) {
	VT_ComputeSliceFromDistances( &imDist1, &imDist2,
				      &imres, i, a, b );
      } else {
	VT_ComputeSliceFromDistances( &imDist1, &imDist2,
				      &imres, i-1, a, b );
      }
    }
  }



  VT_FreeImage( &imDist1 );
  VT_FreeImage( &imDist2 );
  
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( &imres );
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
	_VT_DEBUG_ = 1;
      }
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }



      else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nb...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nb) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nb...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-i" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -i...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->indice) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -i...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	_DistanceSetComputationTo2D();
      }
      else if ( strcmp ( argv[i], "-noext" ) == 0 ) {
	par->writeOriginaux = 0;
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
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  
  /*--- type de l'image resultat ---*/
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
  par->nb = 1;
  par->indice = -1;
  par->writeOriginaux = 1;
}
