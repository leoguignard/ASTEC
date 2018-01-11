#include <vt_common.h>
#include <chamfer.h>
#include <time.h>

typedef struct local_par {
  vt_names names;
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











static void _ComputeTimeFrom2Clocks( float c1, float c2,
				     int *hours, int *minutes, float *seconds )
{
  double d = ( (double)c2 / (double)CLOCKS_PER_SEC ) - 
    ( (double)c1 / (double)CLOCKS_PER_SEC );
  *hours = *minutes = 0;
  *seconds = 0.0;

  if ( d > 3600 ) {
    *hours = (int)(d / 3600);
    d -= *hours * 3600.0;
  }
  if ( d > 60 ) {
    *minutes = (int)(d / 60);
    d -= *minutes * 60.0;
  }
  *seconds = d;
}

static void _PrintTimeFrom2Clocks( float c1, float c2 )
{
  int h, m;
  float s;
  _ComputeTimeFrom2Clocks( c1, c2, &h, &m, &s );
  if ( h > 0 ) printf(" %d h", h);
  if ( m > 0 ) printf(" %d mn", m);
  printf(" %f s\n", s);
}








static char *usage = "[image-in] [image-out]\n\
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
  vt_image image, imres;
  int x, y, z;
  u8 ***theBuf;
  int theDim[3];
  float t0, t1;


  /*--- initialisation des parametres ---*/
  /* VT_InitParam( &par ); */
  
  /*--- lecture des parametres ---*/
  /* VT_Parse( argc, argv, &par ); */
  
  /*--- creation de l'image d'entree ---*/
  VT_Image( &image );
  VT_InitImage( &image, "image-test-2D.inr", 65, 65, 65, UCHAR );
  (void)VT_AllocImage( &image );
  
  VT_InitFromImage( &imres, &image, "image-chamfer3-2D.inr", image.type );
  (void)VT_AllocImage( &imres );



  theBuf = (u8***)image.array;

  theDim[0] = image.dim.x;
  theDim[1] = image.dim.y;
  theDim[2] = image.dim.z;



  for ( z=0; z<65; z++ )
  for ( y=0; y<65; y++ )
  for ( x=0; x<65; x++ ) {
    theBuf[z][y][x] = 0;
  }
  for ( z=0; z<65; z++ )
    theBuf[z][z][z] = 255;
  (void)VT_WriteInrimage( &image );




  sprintf( imres.name, "image-chamfer3-2D.inr" );
  fprintf( stderr," ... compute 2D chamfer 3 ..." );
  t0 = (float)clock();
  (void)Compute2DChamfer3x3( image.buf, image.type, imres.buf, imres.type, 
			   theDim );
  t1 = (float)clock();
  fprintf( stderr,"\n" );
  _PrintTimeFrom2Clocks( t0, t1 );
  (void)VT_WriteInrimage( &imres );



  sprintf( imres.name, "image-chamfer5-2D.inr" );
  fprintf( stderr," ... compute 2D chamfer 5 ..." );
  t0 = (float)clock();
  (void)Compute2DChamfer5x5( image.buf, image.type, imres.buf, imres.type, 
			   theDim );
  t1 = (float)clock();
  fprintf( stderr,"\n" );
  _PrintTimeFrom2Clocks( t0, t1 );
  (void)VT_WriteInrimage( &imres );



  for ( z=0; z<65; z++ )
  for ( y=0; y<65; y++ )
  for ( x=0; x<65; x++ ) {
    theBuf[z][y][x] = 0;
  }
  for ( z=32; z<=32; z++ )
    theBuf[z][z][z] = 255;
  sprintf( image.name, "image-test-3D.inr" );
  (void)VT_WriteInrimage( &image );



  sprintf( imres.name, "image-chamfer3-3D.inr" );
  fprintf( stderr," ... compute 3D chamfer 3 ..." );
  t0 = (float)clock();
  (void)Compute3DChamfer3x3x3( image.buf, image.type, imres.buf, imres.type, 
			   theDim );
  t1 = (float)clock();
  fprintf( stderr,"\n" );
  _PrintTimeFrom2Clocks( t0, t1 );
  (void)VT_WriteInrimage( &imres );



  sprintf( imres.name, "image-chamfer5-3D.inr" );
  fprintf( stderr," ... compute 3D chamfer 5 ..." );
  t0 = (float)clock();
  (void)Compute3DChamfer5x5x5( image.buf, image.type, imres.buf, imres.type, 
			   theDim );
  t1 = (float)clock();
  fprintf( stderr,"\n" );
  _PrintTimeFrom2Clocks( t0, t1 );
  (void)VT_WriteInrimage( &imres );

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
}
