#include <vt_common.h>
#include <vt_copy.h>
typedef struct local_par {
  vt_names names;
  int size;
  int seuil;
  int offset;
  int type;
  int H;
} local_par;

typedef struct {
  int x;
  int y;
} pt2D;

typedef struct {
  pt2D centre;
  int nb;
  pt2D *pts;
} structuringElement;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static void ConvexHullWithPerimeterWindow( vt_image *theIm, 
			vt_image *resIm,
			int size, int seuil );
static void ConvexHullWithVolumeWindow( vt_image *theIm, 
			vt_image *resIm,
			int size, int seuil );
static void ConvexHullWithVolumeDisk( vt_image *theIm, 
				vt_image *resIm,
				int size, int offset);
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

#define _WINDOW_PERIMETER_ 1
#define _WINDOW_VOLUME_ 2
#define _DISK_VOLUME_ 3
static int _methode_ = _WINDOW_PERIMETER_ ;
static int _onlyborderpoints_ = 1;
static int percent = 50;
static char *usage = "[image-in] [image-out]\n\
\t[-size %d] [-seuil %d] [-offset %d] [-border|-noborder]\n\
\t [-wp|-wv|-dv]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -size %d : fenetre de taille N=2*size+1\n\
\t -seuil %d : on change le point si le nombre de ses voisins\n\
\t             est strictement superieur au seuil\n\
\t             old : seuil = 4*size -1\n\
\t             new : seuil = 4*size\n\
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

#ifndef NO_PROTO
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image, imtest, imres, imaux;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  if ( par.H <= 0 ) {
    image = _VT_Inrimage( par.names.in );
    if ( image == (vt_image*)NULL ) 
      VT_ErrorParse("unable to read input image\n", 0);
  } else {
    
    int x,y;
    unsigned char *test;
    int dim, first, second, middle;
    image = &imtest;
    VT_Image( &imtest );
    dim = 201;
    VT_InitImage( &imtest, "image.test.H", dim, dim, 1, UCHAR );
    VT_AllocImage( &imtest );
    test = (unsigned char *)(imtest.buf);
    first = dim/3; second = dim - first; middle = dim/2;
    /* objet selon les axes */
    /*
    for (y=0;y<first;y++) for(x=0;x<dim;x++) test[x+y*dim] = 0;
    for (y=first;y<second;y++) {
      for(x=0;x<middle-par.H/2;x++) test[x+y*dim] = 255;
      for(x=middle-par.H/2; x<middle-par.H/2+par.H; x++ ) test[x+y*dim] = 0;
      for(x=middle-par.H/2+par.H;x<dim;x++) test[x+y*dim] = 255;
    }
    for (y=second;y<dim;y++) for(x=0;x<dim;x++) test[x+y*dim] = 255;
    */ 
    /* objet a 45 deg. */
    for (y=0;y<dim;y++) 
    for(x=0;x<dim;x++) {
      if ( x -y < first -second ) { test[x+y*dim] = 255; continue; }
      if ( x -y > second -first ) { test[x+y*dim] = 0; continue; }
      if ( x+y < first+second -0.707*par.H ) { test[x+y*dim] = 255; continue; }
      if ( x+y > first+second +0.707*par.H ) { test[x+y*dim] = 255; continue; }
      test[x+y*dim] = 0;
    }
  }

  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitImage( &imres, par.names.out, image->dim.x, image->dim.y, image->dim.z, image->type );
  if ( par.type != TYPE_UNKNOWN ) imres.type = par.type;
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  VT_Image( &imaux );
  VT_InitImage( &imaux, par.names.out, image->dim.x, image->dim.y, image->dim.z, image->type );
  if ( VT_AllocImage( &imaux ) != 1 ) {
    VT_FreeImage( image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  VT_CopyImage( image, &imaux );

  switch( _methode_ ) {
  case _DISK_VOLUME_ :
    ConvexHullWithVolumeDisk( &imaux, &imres, par.size, par.offset );
    break;
  case _WINDOW_VOLUME_ :
    ConvexHullWithVolumeWindow( &imaux, &imres, par.size, par.seuil );
    break;
  case  _WINDOW_PERIMETER_ :
  default :
    ConvexHullWithPerimeterWindow( &imaux, &imres, par.size, par.seuil );
  }
      
  {
    unsigned char* res = (unsigned char*)(imres.buf);
    unsigned char* ori = (unsigned char*)(image->buf);
    int i,v=image->dim.x * image->dim.y;
    for ( i=0; i<v; i++ ) {
      if ( ori[i] > 0 ) {
	res[i] = 255;
      } else {
	if ( res[i] > 0 ) res[i] = 200;
      }
    }
  }

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_FreeImage( &imaux );
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
      /*--- ---*/
      else if ( strcmp ( argv[i], "-H" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -H...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->H) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -H...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-size" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -size...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->size) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -size...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-seuil" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -seuil...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->seuil) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -seuil...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-offset" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -offset...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->offset) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -offset...\n", 0 );
	_methode_ = _DISK_VOLUME_;
      }
      else if ( strcmp ( argv[i], "-border" ) == 0 ) {
	_onlyborderpoints_ = 1;
      }
      else if ( strcmp ( argv[i], "-noborder" ) == 0 ) {
	_onlyborderpoints_ = 0;
      }
      else if ( strcmp ( argv[i], "-wp" ) == 0 ) { /* window perimeter */
	_methode_ = _WINDOW_PERIMETER_;
      }
      else if ( strcmp ( argv[i], "-wv" ) == 0 ) { /* window volume */
	_methode_ = _WINDOW_VOLUME_;
      }
      else if ( strcmp ( argv[i], "-dv" ) == 0 ) { /* window volume */
	_methode_ = _DISK_VOLUME_;
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
  par->size = 1;
  par->seuil = 3;
  par->offset = 0;
  par->type = TYPE_UNKNOWN;
  par->H = 0;
}

static void ConvexHullWithPerimeterWindow( vt_image *theIm, 
			vt_image *resIm,
			int size, int seuil )
{
  unsigned char *theBuf=(unsigned char *)theIm->buf;
  unsigned char *resBuf=(unsigned char *)resIm->buf;
  int n, v=theIm->dim.x * theIm->dim.y;
  int dimx=theIm->dim.x;
  int x,y,i,j,l;
  int iteration=0;
  
  if ( size <= 0 ) {
    size = 1;
    seuil = 3; 
    /* borgefors size =1 seuil = 3
       zimmer 1 : size, N= 2 * size + 1, 
                        seuil = 2 * N - 3 = 4 * size - 1
       zimmer 2 : size, N= 2 * size + 1, 
                        seuil = 2 * N - 2 = 4 * size
     */
  }

  seuil = 4*size;
  for (x=0;x<v;x++) resBuf[x]=theBuf[x];

  do {
    n=0;
    iteration++;

    for (y=size; y<(theIm->dim.y)-size; y++)
    for (x=size; x<(theIm->dim.x)-size; x++) {
      if ( theBuf[x + y*dimx] > 0 ) continue;

      l = 0;
      
      if ( _onlyborderpoints_ == 1 ) {
	if ( (theBuf[x-1 + y*dimx] == 0) &&
	     (theBuf[x+1 + y*dimx] == 0) &&
	     (theBuf[x + (y-1)*dimx] == 0) &&
	     (theBuf[x + (y+1)*dimx] == 0) )
	  continue;
      }
      
      for ( i= (-size); i<=size; i++ ) {
	if ( theBuf[x+i + (y-size)*dimx] > 0 ) l++;
	if ( theBuf[x+i + (y+size)*dimx] > 0 ) l++;
      }
      for ( j= (-size)+1; j<size; j++ ) {
	if ( theBuf[x+size + (y+j)*dimx] > 0 ) l++;
	if ( theBuf[x-size + (y+j)*dimx] > 0 ) l++;
      }
      
      if ( l > seuil ) {
	resBuf[x + y*dimx] = 255;
	n++;
      }
    }
    
    for (x=0;x<v;x++) theBuf[x]=resBuf[x];

    if ( iteration % percent == 0 )
      fprintf(stderr,"... iteration %3d -> %5d points modifies\n", iteration, n);

  } while ( n > 0 );
}

static void ConvexHullWithVolumeWindow( vt_image *theIm, 
			vt_image *resIm,
			int size, int seuil )
{
  unsigned char *theBuf=(unsigned char *)theIm->buf;
  unsigned char *resBuf=(unsigned char *)resIm->buf;
  int n, v=theIm->dim.x * theIm->dim.y;
  int dimx=theIm->dim.x;
  int x,y,i,j,l;
  int iteration=0;
  
  if ( size <= 0 ) {
    size = 1;
    seuil = 3; 
    /* borgefors size =1 seuil = 3
       zimmer 1 : size, N= 2 * size + 1, 
                        seuil = 2 * N - 3 = 4 * size - 1
       zimmer 2 : size, N= 2 * size + 1, 
                        seuil = 2 * N - 2 = 4 * size
     */
  }

  for (x=0;x<v;x++) resBuf[x]=theBuf[x];

  do {
    n=0;
    iteration++;

    for (y=size; y<(theIm->dim.y)-size; y++)
    for (x=size; x<(theIm->dim.x)-size; x++) {
      if ( theBuf[x + y*dimx] > 0 ) continue;

      l = 0;
      
      if ( _onlyborderpoints_ == 1 ) {
	if ( (theBuf[x-1 + y*dimx] == 0) &&
	     (theBuf[x+1 + y*dimx] == 0) &&
	     (theBuf[x + (y-1)*dimx] == 0) &&
	     (theBuf[x + (y+1)*dimx] == 0) )
	  continue;
      }
      
      for ( i= (-size); i<=size; i++ ) 
      for ( j= (-size); j<=size; j++ ) {
	if ( theBuf[x+i + (y+j)*dimx] > 0 ) l++;
      }
      
      if ( l > seuil ) {
	resBuf[x + y*dimx] = 255;
	n++;
      }
    }
    
    for (x=0;x<v;x++) theBuf[x]=resBuf[x];

    if ( iteration % percent == 0 ) 
      fprintf(stderr,"... iteration %3d -> %5d points modifies\n", iteration, n);

  } while ( n > 0 );
}

static void ConvexHullWithVolumeDisk( vt_image *theIm, 
				vt_image *resIm,
				int size, int offset)
{
  unsigned char *theBuf=(unsigned char *)theIm->buf;
  unsigned char *resBuf=(unsigned char *)resIm->buf;
  int nb, dimx=theIm->dim.x;
  int x,y, v=theIm->dim.x * theIm->dim.y;
  int i,j, l, n, seuil;
  int iteration=0;
  structuringElement se;
  
  /* distance euclidienne */
  se.centre.x = i = theIm->dim.x / 2;
  se.centre.y = j = theIm->dim.y / 2;
  nb = 0;
  for (y=0; y<theIm->dim.y; y++)
  for (x=0; x<theIm->dim.x; x++) {
    if ( (x-i)*(x-i)+(y-j)*(y-j) <= size*size ) {
      resBuf[x + y*dimx] = 255; nb ++;
    } else 
      resBuf[x + y*dimx] = 0;
  }
  /* elemenet structurant */
  fprintf(stderr," element structurant de %d points\n", nb );
  se.nb = nb;
  se.pts = (pt2D*)malloc( nb * sizeof( pt2D ));
  nb = 0;
  for (y=0; y<theIm->dim.y; y++)
  for (x=0; x<theIm->dim.x; x++) {
    if ( resBuf[x + y*dimx] == 255 ) {
      se.pts[nb].x = x - se.centre.x;
      se.pts[nb].y = y - se.centre.y;
      nb++;
    }
  }

  seuil = (se.nb-1)/2 - 1 + offset;
  fprintf(stderr," seuil = %d ( = %d + %d ) \n", seuil,(se.nb-1)/2 - 1, offset  );
  for (x=0;x<v;x++) resBuf[x]=theBuf[x];

  do {
    n=0;
    iteration++;

    for (y=size; y<(theIm->dim.y)-size; y++)
      for (x=size; x<(theIm->dim.x)-size; x++) {
      if ( theBuf[x + y*dimx] > 0 ) continue;

      l = 0;
      
      if ( _onlyborderpoints_ == 1 ) {
	if ( (theBuf[x-1 + y*dimx] == 0) &&
	     (theBuf[x+1 + y*dimx] == 0) &&
	     (theBuf[x + (y-1)*dimx] == 0) &&
	     (theBuf[x + (y+1)*dimx] == 0) )
	  continue;
      }
      
      for ( i=0; i < se.nb; i++ ) {
	if ( theBuf[x+se.pts[i].x + (y+se.pts[i].y)*dimx] > 0 ) l++;
      }

      if ( l > seuil ) {
	resBuf[x + y*dimx] = 255;
	n++;
      }
    }
    
    for (x=0;x<v;x++) theBuf[x]=resBuf[x];

    if ( iteration % percent == 0 )
      fprintf(stderr,"... iteration %3d -> %5d points modifies\n", iteration, n);

  } while ( n > 0 );
}

  
