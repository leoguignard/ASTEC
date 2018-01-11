#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  int size;
  int bool_xmgr;
} local_par;

typedef struct local_pt {
  double sum_percent;
  int nb;
} local_pt;

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

static char *usage = "image-in [-size %d] [-xmgr]";
static char *detail = "";

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
  vt_image *image;
  local_pt *pt;
  register int x, y, z, i, j, k;
  int s2, d2, k2, j2, indz, indy, dxy, dx, dy, dz;
  u16 ***theBuf;
  double percent;

  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- 1ere image ---*/
  image = _VT_Inrimage( par.names.in );
  
  /*--- tableau ---*/
  s2 = par.size * par.size;
  pt = (local_pt*)VT_Malloc( (s2+1) * sizeof(local_pt) );
  for ( i = 0; i <= s2; i ++ ) {
    pt[i].nb = 0;
    pt[i].sum_percent = 0.0;
  }

  dx = image->dim.x;
  dy = image->dim.y;
  dz = image->dim.z;
  dxy = dx * dy;
  theBuf = (u16***)image->array;

  for ( z = 0; z < dz; z ++ )
  for ( y = 0; y < dy; y ++ )
  for ( x = 0; x < dx; x ++ ) {
    if ( theBuf[z][y][x] == (u16)0 ) continue;

    for ( k = 0 ; k <= par.size ; k ++ ) {
      if ( (z+k < 0) || (z+k >= dz) ) continue;
      k2 = k * k;
      indz = k * dxy;
      for ( j = (-par.size) ; j <= par.size ; j ++ ) {
	if ( (y+j < 0) || (y+j >= dy) ) continue;
	j2 = j * j;
	if ( k2 + j2 > s2 ) continue;
	indy = indz + j * dx;
	if ( indy <= (-par.size) ) continue;
	for ( i = (-par.size) ; i <= par.size ; i ++ ) {
	  if ( (x+i < 0) || (x+i >= dx) ) continue;
	  d2 = k2 + j2 + i * i;
	  if ( d2 > s2 ) continue;
	  if ( indy + i <= 0 ) continue;
	  /*--- le point est a l'interieur, apres et a moins de size du point central ---*/
	  if ( theBuf[z+k][y+j][x+i] == (u16)0 ) continue;
	  if ( theBuf[z+k][y+j][x+i] > theBuf[z][y][x] )
	    percent = ( (double)theBuf[z+k][y+j][x+i] - (double)theBuf[z][y][x] );
	  else
	    percent = ( (double)theBuf[z][y][x] - (double)theBuf[z+k][y+j][x+i] );
	  percent /= (double)theBuf[z][y][x];
	  pt[d2].sum_percent += (percent * 100.0);
	  pt[d2].nb ++;
	}
      }
    }
    
  }

  if ( par.bool_xmgr != 1 ) {
    for ( i = 0; i <= s2; i ++ ) {
      if ( pt[i].nb > 0 ) {
	fprintf( stdout, " distance = %f ; pourcentage moyen = %f \n", (float)sqrt( (double)i ), 
		 pt[i].sum_percent / (double)pt[i].nb );
      }
    }
    fprintf( stdout, "\n\n" );
  } else {
    for ( i = 0; i <= s2; i ++ ) {
      if ( pt[i].nb > 0 ) {
	fprintf( stdout, "%f %f \n", (float)sqrt( (double)i ), 
		 pt[i].sum_percent / (double)pt[i].nb );
      }
    }
  }
  return( 0 );
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

  /*--- initialisation des parametres ---*/
  VT_InitParam( par );
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);

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
      /*--- size ---*/
      else if ( strcmp ( argv[i], "-size" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -size...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->size) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -size...\n", 0 );
      }      
      else if ( strcmp ( argv[i], "-xmgr" ) == 0 ) {
	par->bool_xmgr = 1;
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
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }

  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) 
    VT_ErrorParse("no file name when parsing\n", 0 );

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
  par->size = 10;
  VT_Names( &(par->names) );
  par->bool_xmgr = 0;
}
