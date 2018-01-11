#include <vt_common.h>


typedef struct local_par {
  vt_names names;
  vt_names refs;
  int seuilXenon;
  int seuilHmpao;
  int type;
  int x;
  int y;
  int z;
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

static char *usage = "[image-hmpao] [image-xenon] [image-out]\n\
\t [-xref %s] [-href %s] [-sh %d] [-sx %d]\n\
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
  vt_image *imageHmpao, *imageXenon, imres;
  vt_image *refHmpao, *refXenon;

  int x, y, z;
  int dimx, dimy, dimz;
  int xenon, hmpao;

  double xrel, hrel;

  int valRefXenon=1.0, valRefHmpao=1.0;
  float ***resBuf = (float***)NULL;






  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  if ( 1 ) {
    fprintf( stderr, " ... image hmpao : %s\n", par.names.in );
    fprintf( stderr, " ... image xenon : %s\n", par.names.ext );
    fprintf( stderr, "     les deux images sont sensees etre recalees.\n" );
  }

  /*--- lecture de l'image d'entree ---*/
  imageHmpao = _VT_Inrimage( par.names.in );
  if ( imageHmpao == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read hmpao image\n", 0);
  imageXenon = _VT_Inrimage( par.names.ext );
  if ( imageXenon == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read xenon image\n", 0);


  refHmpao = _VT_Inrimage( par.names.in );
  if ( refHmpao == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read ref hmpao image\n", 0);
  refXenon = _VT_Inrimage( par.names.ext );
  if ( refXenon == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read ref xenon image\n", 0);



  VT_Image( &imres );
  VT_InitFromImage( &imres, imageHmpao, par.names.out, FLOAT );
  (void)VT_AllocImage( &imres );
  resBuf = (float***)imres.array;
  for ( z=0; z<imres.dim.z; z++ )
  for ( y=0; y<imres.dim.y; y++ )
  for ( x=0; x<imres.dim.x; x++ )
    resBuf[z][y][x] = 0.0;





  if ( 1 ) {
    fprintf( stderr, " reference at %d %d %d\n",
	     par.x, par.y, par.z );
  }
  switch ( refHmpao->type ) {
  case USHORT :
    {
      u16 *** theHmpao = (u16 ***)refHmpao->array;
      
      switch ( refXenon->type ) {
      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)refXenon->array;
	  
	  valRefXenon = theXenon[par.z][par.y][par.x];
	  valRefHmpao = theHmpao[par.z][par.y][par.x];

	}
	break;
      default :
	VT_ErrorParse("unable to deal with such ref xenon image type\n", 0);
      }
      
    }
    break;
  default :
    VT_ErrorParse("unable to deal with such ref hmpao image type\n", 0);
  }





  
  dimz = imageHmpao->dim.z;
  dimy = imageHmpao->dim.y;
  dimx = imageHmpao->dim.x;
  if ( dimz > imageXenon->dim.z ) dimz = imageXenon->dim.z;
  if ( dimy > imageXenon->dim.y ) dimy = imageXenon->dim.y;
  if ( dimx > imageXenon->dim.x ) dimx = imageXenon->dim.x;



  switch ( imageHmpao->type ) {
  case USHORT :
    {
      u16 *** theHmpao = (u16 ***)imageHmpao->array;
      
      switch ( imageXenon->type ) {
      case USHORT :
	{
	  u16 *** theXenon = (u16 ***)imageXenon->array;
	  
	  for ( z=0; z<dimz; z++ )
	  for ( y=0; y<dimy; y++ )
	  for ( x=0; x<dimx; x++ ) {
	    xenon = theXenon[z][y][x];
	    hmpao = theHmpao[z][y][x];
	    if ( xenon < par.seuilXenon ) continue;
	    if ( hmpao < par.seuilHmpao ) continue;

	    xrel = (double)xenon / (double)valRefXenon;
	    hrel = (double)hmpao / (double)valRefHmpao;
	    if ( xrel - hrel > 1e-8 || xrel - hrel < -1e-8 )
	      resBuf[z][y][x] = (1.0 - hrel) * xrel / (xrel - hrel);
	  }

	}
	break;
      default :
	VT_ErrorParse("unable to deal with such xenon image type\n", 0);
      }

    }
    break;
  default :
    VT_ErrorParse("unable to deal with such hmpao image type\n", 0);
  }












  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_ErrorParse("unable to write output image\n", 0);
  }
  VT_FreeImage( &imres );




  /*--- liberations memoires ---*/


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




      else if ( strcmp ( argv[i], "-sh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sh...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->seuilHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sh...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-sx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sx...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->seuilXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sx...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-href" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -href...\n", 0 );
	strncpy( par->refs.in, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-xref" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -xref...\n", 0 );
	strncpy( par->refs.ext, argv[i], STRINGLENGTH );  
      }

      
      else if ( strcmp ( argv[i], "-pt" ) == 0 ) {
	i += 1;
	if ( i+2 >= argc)    VT_ErrorParse( "parsing -pt...\n", 0 );
	status = sscanf( argv[i],"%d",&par->x );
	i += 1;
	status += sscanf( argv[i],"%d",&par->y );
	i += 1;
	status += sscanf( argv[i],"%d",&par->z );
	if ( status != 3 ) VT_ErrorParse( "parsing -pt...\n", 0 );
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
  VT_Names( &(par->refs) );
  par->type = TYPE_UNKNOWN;
  par->seuilHmpao = 0;
  par->seuilXenon = 0;
  par->x = 0;
  par->y = 0;
  par->z = 0;
}




