#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  int type;
  double a;
  double b;
  char matrice[STRINGLENGTH];
} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static int  ReadMatrice( char *name, double *mat );
static int  InverseMat4x4( double *matrice, double *inv );

#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
static int  ReadMatrice();
static int  InverseMat4x4();

#endif

static char *usage = "[image-hmapo] [image-xenon] -mat %s [image-out]\n\
\t [-a %lf] [-b %lf]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t Applique la transformation: a*I/(b-I)\n\
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
  vt_image *imageHmpao, *imageXenon;
  vt_3m mXenon;

  int n;
  double mat[16];

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  


  /*--- lecture de l'image d'entree ---*/
  imageHmpao = _VT_Inrimage( par.names.in );
  if ( imageHmpao == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read hmpao image\n", 0);
  imageXenon = _VT_Inrimage( par.names.ext );
  if ( imageXenon == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read hmpao image\n", 0);
  
  if ( imageHmpao->type != imageXenon->type ) 
    VT_ErrorParse("unable to deal with images of differents type\n", 0);


  (void)VT_3m( imageXenon, (vt_image*)NULL, &mXenon );

  sprintf( imageHmpao->name, "%s", par.names.out );






  /*  0  1  2  3
      4  5  6  7
      8  9 10 11
     12 13 14 15 
     */
  for ( n=0; n<16; n++ ) mat[n] = 0.0;
  mat[0] = mat[5] = mat[10] = mat[15] = 1.0;

  if ( par.matrice[0] != '\0' ) {
    if ( ReadMatrice( par.matrice, mat ) != 1 ) {
      VT_ErrorParse("unable to read matrice\n", 0);
    }
    if ( par.names.inv != 0 ) {
      double tmp[16];
      int rang;
      if ( (rang=InverseMat4x4( mat, tmp )) != 4 ) {
	fprintf( stderr, "Warning: la matrice %s est de rang %d.\n",
		 par.matrice, rang );
      }
      for ( n=0; n<16; n++ ) mat[n] = tmp[n];
    } 
  }




  fprintf( stderr, "Recalibre en appliquant la transformation %f I / (%f - I)\n",
	   par.a, par.b );

  

  
  switch ( imageHmpao->type ) {
  case USHORT :
    {
      unsigned short int ***theHmpao = (unsigned short int ***)imageHmpao->array;
      unsigned short int ***theXenon = (unsigned short int ***)imageXenon->array;
      int i, j, k;
      double hmpao, xenon, x, y, z;
      double dx, dy, dz;
      double d, dmax;
      int ix, iy, iz;
      

      for ( k=0; k<imageHmpao->dim.z; k++ ) {
	fprintf( stderr, "... processing slice %3d/%lu\r", k, imageHmpao->dim.z-1 );
      for ( j=0; j<imageHmpao->dim.y; j++ )
      for ( i=0; i<imageHmpao->dim.x; i++ ) {
	hmpao = theHmpao[k][j][i];
	theHmpao[k][j][i] = 0;
	
	x = mat[0] * i +  mat[1] * j + mat[2] * k + mat[3];
	ix = (int)x;
	if ( x < 0.0 || ix >= imageXenon->dim.x-1 ) continue;
	y = mat[4] * i +  mat[5] * j + mat[6] * k + mat[7];
	iy = (int)y;
	if ( y < 0.0 || iy >= imageXenon->dim.y-1 ) continue;
	z = mat[8] * i +  mat[9] * j + mat[10] * k + mat[11];
	iz = (int)z;
	if ( z < 0.0 || iz >= imageXenon->dim.z-1 ) continue;
	
	dx = x - ix;
	dy = y - iy;
	dz = z - iz;

	xenon = 0;
	xenon += theXenon[iz][iy][ix]       * (1.0-dz)*(1.0-dy)*(1.0-dx);
	xenon += theXenon[iz][iy][ix+1]     * (1.0-dz)*(1.0-dy)*     dx ;
	xenon += theXenon[iz][iy+1][ix]     * (1.0-dz)*     dy *(1.0-dx);
	xenon += theXenon[iz][iy+1][ix+1]   * (1.0-dz)*     dy *     dx ;
	xenon += theXenon[iz+1][iy][ix]     *      dz *(1.0-dy)*(1.0-dx);
	xenon += theXenon[iz+1][iy][ix+1]   *      dz *(1.0-dy)*     dx ;
	xenon += theXenon[iz+1][iy+1][ix]   *      dz *     dy *(1.0-dx);
	xenon += theXenon[iz+1][iy+1][ix+1] *      dz *     dy *     dx ;
	

	dmax = hmpao*hmpao + xenon*xenon;
	iy = xenon;
	for ( ix = 1; ix <= mXenon.max+0.5; ix ++ ) {
	  d = (hmpao-par.b*ix/(par.a+ix))*(hmpao-par.b*ix/(par.a+ix))
	    + (xenon-ix)*(xenon-ix);
	  if ( dmax > d ) {
	    dmax = d;
	    iy = ix;
	  }
	}
	theHmpao[k][j][i] = iy;
      }
      }
      
    }
    break;
  default :
    VT_ErrorParse("unable to deal with such input image type\n", 0);
  }
    

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( imageHmpao ) == -1 ) {
    VT_FreeImage( imageHmpao );
    VT_Free( (void**)&imageHmpao );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( imageHmpao );
  VT_Free( (void**)&imageHmpao );
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



      else if ( strcmp ( argv[i], "-mat" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mat...\n", 0 );
	strncpy( par->matrice, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-a" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -a...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->a) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -a...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-b" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -b...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->b) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -b...\n", 0 );
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
  par->type = TYPE_UNKNOWN;
  par->a = 1.0;
  par->b = 0.0;
  par->matrice[0] = '\0';
}






















					
#if defined(_ANSI_)
static int ReadMatrice( char *name, double *mat )
#else
static int ReadMatrice( name, mat )
char *name;
double *mat;
#endif
{
  FILE *fopen(), *fp;
  char text[STRINGLENGTH];
  int nbelts = 0;
  int status;
  
  /* lecture de 4 double par ligne
     On prevoit le cas ou la ligne commence par "O8 xxxxx ...
     */

  fp = fopen( name, "r" );
  if ( fp == NULL ) return( 0 );
  
  while ( (nbelts < 16) && (fgets( text, STRINGLENGTH, fp ) != NULL) ) {
    if ( (text[0] == 'O') && (text[1] == '8') ) {
      status = sscanf( &(text[2]), "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    } else {
      status = sscanf( text, "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    }
    if ( status == 4 ) nbelts += 4;
  }
  fclose( fp );

  if ( _VT_DEBUG_ == 1 ) {
    fprintf( stderr, " lecture de la matrice %s\n", name );
    fprintf( stderr, " %d elements lus\n", nbelts );
    fprintf( stderr,"   %f %f %f %f\n", mat[0], mat[1], mat[2], mat[3] );
    fprintf( stderr,"   %f %f %f %f\n", mat[4], mat[5], mat[6], mat[7] );
    fprintf( stderr,"   %f %f %f %f\n", mat[8], mat[9], mat[10], mat[11] );
    fprintf( stderr,"   %f %f %f %f\n", mat[12], mat[13], mat[14], mat[15] );
  }
  if ( nbelts == 16 ) return ( 1 );
  return( 0 );
}


#define TINY 1e-12

#if defined(_ANSI_)
static int InverseMat4x4( double *matrice, double *inv )
#else
static int InverseMat4x4( matrice, inv )
double *matrice;
double *inv;
#endif
{
  register int i, j, k;
  int kmax, rang = 4;
  register double c, max;
  double mat [16];
  
  for (i=0; i<16; i++ ) {
    mat[i] = matrice[i] ;
    inv[i] = 0.0;
  }
  inv[0] = inv[5] = inv[10] = inv[15] = 1.0;
  
  for ( j=0; j<4; j++ ) {
    if ( (mat[j*4+j] > (-TINY)) && (mat[j*4+j] < TINY) ) {
      /* recherche du plus grand element non nul sur la colonne j */
      kmax = j;
      max = 0.0;
      for (k=j+1; k<4; k++ ) {
	c = ( mat[k*4+j] > 0.0 ) ? mat[k*4+j] : (-mat[k*4+j]) ;
	if ( (c > TINY) && (c > max) ) { max = c; kmax = k; }
      }
      if ( kmax == j ) {
	/* la ligne est nulle */
	rang --;
      } else {
	/* sinon, on additionne */
	for ( i=0; i<4; i++ ) {
	  mat[j*4+i] += mat[kmax*4+i];
	  inv[j*4+i] += inv[kmax*4+i];
	}
      }
    }
    if ( (mat[j*4+j] < (-TINY)) || (mat[j*4+j] > TINY) ) {
      /* les autres lignes */
      for (k=0; k<4; k++) {
	if ( k != j ) {
	  c = mat[k*4 + j] / mat[j*4 + j];
	  for ( i=0; i<4; i++ ) {
	    mat[k*4 + i] -= c * mat[j*4 + i];
	    inv[k*4 + i] -= c * inv[j*4 + i];
	  }
	}
      }
      /* la ligne */
      c = mat[j*4 + j];
      for ( i=0; i<4; i++ ) {
	mat[j*4 + i] /= c;
	inv[j*4 + i] /= c;
      }
    }
  }

  return( rang );
}
