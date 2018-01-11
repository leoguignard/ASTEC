#include <vt_common.h>
#include <vt_interpol.h>

typedef struct local_par {
  vt_names names;
  int type;

  float size[3];
  double refine;

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

static char *usage = "[image-in] [image-out] [-vz %f %f %f] [-2D] \n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
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
  vt_image *image, imres, imschar;
  typeDistanceMap theDist;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  if ( (image->dim.z != 1) || (image->type != UCHAR) ) 
    VT_ErrorParse("unable to deal with such image\n", 0);

  /*--- operations eventuelles sur l'image d'entree ---*/
  /*
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  */

  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitImage( &imres, par.names.out, image->dim.x, image->dim.y, 256, SSHORT );
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  VT_InitDistanceImageFromSlice( image, &imres, 0 );

  VT_FreeImage( image );
  VT_Free( (void**)&image );
  
  theDist.buf = imres.buf;
  theDist.dim[0] = imres.dim.x;
  theDist.dim[1] = imres.dim.y;
  theDist.dim[2] = imres.dim.z;
  theDist.voxelSize[0] = par.size[0];
  theDist.voxelSize[1] = par.size[1];
  theDist.voxelSize[2] = par.size[2];
  theDist.multiplicativeCoefficient = 0.0;

  /*-- calcul de la carte
   */
  _ComputeSignedDistanceMap( &theDist, par.refine );


  switch ( par.type ) {
  default :
    /*--- ecriture de l'image resultat ---*/
    if ( VT_WriteInrimage( &imres ) == -1 ) {
      VT_FreeImage( &imres );
      VT_ErrorParse("unable to write output image\n", 0);
    }
    break;

  case SCHAR :
    VT_InitFromImage( &imschar, &imres, par.names.out, SCHAR );
    (void)VT_AllocImage( &imschar );
    {
      char *theBuf = (char *)imschar.buf;
      short int *resBuf = (short int *)imres.buf;
      int i, v=imschar.dim.z*imschar.dim.y*imschar.dim.x;
      double d;

      for (i=0; i<v; i++ ) {
	d = (double)resBuf[i] * theDist.multiplicativeCoefficient * 2.0;
	if ( d < -128.0 ) {
	  theBuf[i] = -128;
	} else if ( d < 0.0 ) {
	  theBuf[i] = (int)(d - 0.5);
	} else if ( d > 127.0 ) {
	  theBuf[i] = 127;
	} else {
	  theBuf[i] = (int)(d + 0.5);
	}
      }

    }
    (void)VT_WriteInrimage( &imschar );
    VT_FreeImage( &imschar );
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

      else if ( strcmp ( argv[i], "-refine" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -refine...\n", 0 );
	status = sscanf( argv[i],"%lf",&par->refine );
	if ( status <= 0 ) VT_ErrorParse( "parsing -refine...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-vz" ) == 0 ) {
	if ( i+3 >= argc)    VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+1], "%f", &par->size[0] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+2], "%f", &par->size[1] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+3], "%f", &par->size[2] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	i += 3;
      }


      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	_DistanceSetComputationTo2D();
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
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
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

  par->size[0] = par->size[1] = par->size[2] = 1.0;
  par->refine = -1.0;
}
