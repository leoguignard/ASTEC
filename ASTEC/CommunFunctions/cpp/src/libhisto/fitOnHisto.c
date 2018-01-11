/*************************************************************************
 * fitOnHisto.c -
 *
 * $Id: fitOnHisto.c,v 1.4 2000/07/04 10:38:04 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <levenberg.h>

typedef struct local_par {
  vt_names names;
  int zeromin, zeromax;
  
  int type;
  int min, max;
  double moy;
  double ect;
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
  vt_image *image, imres;
  double *theX, *theY, *theC, *theS;
  int i, j;
  double gauss[3];
  int xmax;
  double ymax, ymin;
  

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  










  VT_Image( &imres );
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    char name[256];
    sprintf( name, "%s.diff.inr", par.names.out );
    VT_InitFromImage( &imres, image, name, FLOAT );
    (void)VT_AllocImage( &imres );
  }








  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  

  if ( image->dim.y != 1 || image->dim.z != 1 
       || image->dim.v != 1 ) {
    VT_ErrorParse("bad image dimensions\n", 0);
  }



  if ( par.min < 0 ) par.min = 0;
  if ( par.max < 0 ) par.max = image->dim.x-1;
  if ( par.max <= par.min ) {
    VT_ErrorParse("bad interval\n", 0);
  }








  if ( par.zeromin >= 0 && par.zeromax >= 0 ) {

    printf(" on annule de %d a %d\n", par.zeromin, par.zeromax );

    switch ( image->type ) {
    case SINT :
      {
	int *theBuf = (int*)image->buf;
	for ( i=par.zeromin; i<image->dim.x && i<=par.zeromax; i++ ) 
	  theBuf[i] = 0;
      }
      break;
    case FLOAT :
      {
	float *theBuf = (float*)image->buf;
	for ( i=par.zeromin; i<image->dim.x && i<=par.zeromax; i++ ) 
	  theBuf[i] = 0;
      }
      break;
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
    }
  }







  xmax = image->dim.x-1;
  switch ( image->type ) {
  case SINT :
    {
      int *theBuf = (int*)image->buf;
      while ( theBuf[xmax] == 0 )
	xmax --;
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)image->buf;
      while ( theBuf[xmax] == 0.0 )
	xmax --;
    }
    break;
  default :
    VT_ErrorParse("unable to deal with such image\n", 0);
  }

  if ( xmax < par.max ) par.max = xmax;
  




  if ( par.moy < 0 ) par.moy = (par.min+par.max)/2;
  if ( par.ect < 0 ) par.ect = (par.max-par.min)/4;









  theX =  (double*)malloc( 4*(par.max-par.min+1)*sizeof(double) );
  theY =  theX;
  theY += par.max-par.min+1;
  theC =  theY;
  theC += par.max-par.min+1;
  theS =  theC;
  theS += par.max-par.min+1;

  for ( i=0; i<par.max-par.min+1; i++ ) 
    theX[i] = theY[i] = theC[i] = theS[i] = 1.0;
  
  /* 3 parametres pour la gaussienne
     amplitude -> gauss[0]
     moyenne   -> gauss[1]
     ecarttype -> gauss[2]
  */
     
  gauss[1] = par.moy;
  gauss[2] = par.ect;

  switch ( image->type ) {
  case SINT :
    {
      int *theBuf = (int*)image->buf;
      for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      gauss[0] = theBuf[ (int)(gauss[1]) ];
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)image->buf;
      for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
      gauss[0] = theBuf[ (int)(gauss[1]) ];
    }
    break;
  default :
    VT_ErrorParse("unable to deal with such image\n", 0);
  }

  ymax = theY[0];
  for ( i=1, j=par.min+1; j<=par.max; j++, i++ ) {
    if ( ymax < theY[i] ) ymax = theY[i];
  }

  printf( "INITIALISATION amplitude = %f moyenne = %f sigma = %f\n",
	  gauss[0], gauss[1], gauss[2] );
  
  (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				     theC, DOUBLE, theS, DOUBLE, 
				     par.max-par.min+1,
				     gauss, 3, _GaussianForLM );


  printf( "RESULTAT       amplitude = %f moyenne = %f sigma = %f\n",
	  gauss[0], gauss[1], gauss[2] );
  








  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    float *resBuf = (float*)imres.buf;
    double x;
    switch ( image->type ) {
    case SINT :
      {
	int *theBuf = (int*)image->buf;
	for (i=0; i<image->dim.x; i++ )
	  resBuf[i] = theBuf[i];
      }
      break;
    case FLOAT :
      {
	float *theBuf = (float*)image->buf;
	for (i=0; i<image->dim.x; i++ )
	  resBuf[i] = theBuf[i];
      }
      break;
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
    }
    for (i=0; i<image->dim.x; i++ ) {
      x = (i - gauss[1])*(i - gauss[1])/( 2.0*gauss[2]*gauss[2] );
      resBuf[i] -= gauss[0] * exp( - x );
    }
    (void)VT_WriteInrimage( &imres );
  }




  
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    char name[256];
    FILE *fopen(), *fp;
    int fd;
    double x;

    sprintf( name, "%s.raw", par.names.out );
    fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    switch( image->type ) {
    case SINT :
      if ( write( fd, image->buf, image->dim.x*sizeof( int ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      break;
    case FLOAT :
      if ( write( fd, image->buf, image->dim.x*sizeof( float ) ) == -1 ) {
	fprintf( stderr, "%s: error when writing\n", program );
      }
      break;
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
    }

    ymin = ymax;
    switch ( image->type ) {
    case SINT :
      {
	int *theBuf = (int*)image->buf;
	for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	  x = (j - gauss[1])*(j - gauss[1])/( 2.0*gauss[2]*gauss[2] );
	  if ( ymin >  theBuf[j] - gauss[0] * exp( - x ) )
	    ymin = theBuf[j] - gauss[0] * exp( - x );
	}
      }
      break;
    case FLOAT :
      {
	float *theBuf = (float*)image->buf;
	for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	  x = (j - gauss[1])*(j - gauss[1])/( 2.0*gauss[2]*gauss[2] );
	  if ( ymin >  theBuf[j] - gauss[0] * exp( - x ) )
	    ymin = theBuf[j] - gauss[0] * exp( - x );
	}
      }
      break;
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
    }
    close( fd );

    
    sprintf( name, "%s.m", par.names.out );
    fp = fopen( name, "w" );
    
    fprintf( fp, "\n" );
    fprintf( fp, "%% image = %s\n", par.names.in );
    fprintf( fp, "\n" );
    fprintf( fp, " XHISTO = 0:%lu;\n", image->dim.x-1 );
    {
      int k;
      k = strlen( par.names.out );
      for ( i = k-1; i >= 0 && par.names.out[i] != '/' ; i-- )
	;
      fprintf( fp, " fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
    }
    switch ( image->type ) {
    case SINT :
      fprintf( fp, " [TMP, NHISTO] = fread( fid, %lu, 'int%lu' );\n", 
	       image->dim.x, 8*sizeof( int ) );
      break;
    case FLOAT :
      fprintf( fp, " [TMP, NHISTO] = fread( fid, %lu, 'float%lu' );\n", 
	       image->dim.x, 8*sizeof( float ) );
      break;
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
    }
    fprintf( fp, " HISTO = TMP';\n" );
    fprintf( fp, " fclose( fid );\n" );
    fprintf( fp, "\n" );

    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );
    fprintf( fp, "plot(XHISTO, HISTO, 'b-' );\n" );
    fprintf( fp, "G = %f * exp( -(XHISTO - %f).*(XHISTO - %f)/(2*%f*%f) );\n",
	     gauss[0], gauss[1], gauss[1], gauss[2], gauss[2] );
    fprintf( fp, "plot(XHISTO, G, 'g-' );\n" );
    fprintf( fp, "D = HISTO - G;\n" );
    fprintf( fp, "plot(XHISTO, D, 'r-' );\n" );
    fprintf( fp, "axis([0 %d %f %f]);\n", xmax, ymin, ymax );
    fprintf( fp, "hold off;\n" );
    fclose(fp);
  }





  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_Free( (void**)&image );
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
	setVerboseInLevenberg( _VT_VERBOSE_ );
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



      else if ( strcmp ( argv[i], "-zeros" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -zeros...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->zeromin) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -zeros...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -zeros...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->zeromax) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -zeros...\n", 0 );
      }




      else if ( strcmp ( argv[i], "-min" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -min...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->min) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -min...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -max...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->max) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -max...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-moy" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -moy...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->moy) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -moy...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-ect" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -ect...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->ect) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -ect...\n", 0 );
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

  par->zeromin = 0;
  par->zeromax = 0;

  par->type = TYPE_UNKNOWN;
  par->min = -1;
  par->max = -1;
  par->moy = -1;
  par->ect = -1;
}
