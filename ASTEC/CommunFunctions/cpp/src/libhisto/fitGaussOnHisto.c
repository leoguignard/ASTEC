/*************************************************************************
 * fitOnHisto.c -
 *
 * $Id: fitGaussOnHisto.c,v 1.3 2000/07/12 17:54:38 greg Exp $
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


#define NBMAXGAUSS 7
#define NBMAXPARAM 21




typedef struct local_par {
  vt_names names;
  int zeromin, zeromax;
  
  double amp, moy, ect;

  int type;
  int min, max;
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

static char *usage = "[image-in] [nom-generique]\n\
\t [-amp %lf] [-moy %lf] [-ect %lf] [-init %s]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
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
  int nbGauss = 0;
  double gauss[NBMAXPARAM];
  int xmax;
  double ymax, ymin;
  
  FILE *fopen(), *fp;


  





  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  



  if ( par.names.ext[0] != '\0' ) {
    char text[256];
    int s;
    fp = fopen( par.names.ext, "r" );
    if ( fp == NULL )
      VT_ErrorParse("unable to read param file\n", 0);
    while ( fgets( text, STRINGLENGTH, fp ) != NULL ) {
      if ( text[0] == '#' ) continue;
      s = sscanf( text, "%lf %lf %lf",
		  &(gauss[nbGauss*3+0]), &(gauss[nbGauss*3+1]), &(gauss[nbGauss*3+2]) );
      if ( s == 3 ) nbGauss ++;
      else {
	fprintf( stderr, "lecture de '%s', probleme d'interpretation de ...\n   '%s'\n",
		 par.names.ext, text);
      }
    } 
    fclose( fp );
  } else {
    if ( par.moy > 0 ) {
      gauss[1] = par.moy;
      if ( par.amp > 0 ) gauss[0] = par.amp;
      else               gauss[0] = 100;
      if ( par.ect > 0 ) gauss[2] = par.ect;
      else               gauss[2] = 10;
      nbGauss = 1;
    } else {
      VT_ErrorParse("pas de gaussiennes ?\n",0 );
    }
  }







  VT_Image( &imres );
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    char name[256];
    sprintf( name, "%s.diff.inr", par.names.out );
    VT_InitFromImage( &imres, image, name, FLOAT );
    (void)VT_AllocImage( &imres );
  }


  fp = NULL;
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    char name[256];
    sprintf( name, "%s.txt", par.names.out );
    fp = fopen( name, "w" );
    fprintf( fp, "#\n" );
    fprintf( fp, "# " );
    for (i=0; i<argc; i++ )
      fprintf( fp, " %s", argv[i] );
     fprintf( fp, "\n" );
     fprintf( fp, "#\n" );
     fprintf( fp, "#\n" );
     fprintf( fp, "# lecture de %d gaussiennes (amplitude, moyenne, sigma)\n", nbGauss );
     fprintf( fp, "#\n" );
     for ( j=0; j<nbGauss; j++ )
       fprintf( fp, "#   #%2d = (%f %f %f)\n",
		j, gauss[j*3+0], gauss[j*3+1], gauss[j*3+2] ); 
     fprintf( fp, "#\n" );
     fprintf( fp, "#\n" );
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
     

  switch ( image->type ) {
  case SINT :
    {
      int *theBuf = (int*)image->buf;
      for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  case FLOAT :
    {
      float *theBuf = (float*)image->buf;
      for ( i=0, j=par.min; j<=par.max; j++, i++ ) {
	theX[i] = j;
	theY[i] = theBuf[j];
      }
    }
    break;
  default :
    VT_ErrorParse("unable to deal with such image\n", 0);
  }

  ymax = theY[0];
  for ( i=1, j=par.min+1; j<=par.max; j++, i++ ) {
    if ( ymax < theY[i] ) ymax = theY[i];
  }





  
  for ( j=0; j<nbGauss; j++ )
    printf( "#%2d INITIALISATION amplitude = %f moyenne = %f sigma = %f\n",
	  j, gauss[j*3+0], gauss[j*3+1], gauss[j*3+2] );
  





  switch( nbGauss ) {
  default :
    VT_ErrorParse("unable to deal with such number of gaussians\n", 0);
  case 1 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 3, _GaussianForLM );
    break;
  case 2 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 6, _MixtureOf2GaussiansForLM );
    break;
  case 3 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 9, _MixtureOf3GaussiansForLM );
    break;
  case 4 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 12, _MixtureOf4GaussiansForLM );
    break;
  case 5 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 15, _MixtureOf5GaussiansForLM );
    break;
  case 6 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 18, _MixtureOf6GaussiansForLM );
    break;
  case 7 :
    (void)Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, 
				       theC, DOUBLE, theS, DOUBLE, 
				       par.max-par.min+1,
				       gauss, 21, _MixtureOf7GaussiansForLM );
    break;
  }
    

  
  for ( j=0; j<nbGauss; j++ )
    printf( "#%2d RESULTAT       amplitude = %f moyenne = %f sigma = %f\n",
	    j, gauss[j*3+0], gauss[j*3+1], gauss[j*3+2] );
  

  if ( fp != NULL ) {
     fprintf( fp, "# RESULTATS\n" );
     for ( j=0; j<nbGauss; j++ )
       fprintf( fp, " %f %f %f \n",
		gauss[j*3+0], gauss[j*3+1], gauss[j*3+2] ); 
     fclose( fp );
  }







  ymin = 0;
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
      for ( j=0; j<nbGauss; j++ ) {
	x = (i - gauss[3*j+1])*(i - gauss[3*j+1])/( 2.0*gauss[3*j+2]*gauss[3*j+2] );
	resBuf[i] -= gauss[3*j+0] * exp( - x );
      }
      if ( ymin > resBuf[i] ) ymin = resBuf[i];
    }
    (void)VT_WriteInrimage( &imres );
  }




  
  if ( par.names.out[0] != '>' && par.names.out[0] != '\0' ) {
    char name[256];
    FILE *fopen(), *fp;
    int fd;

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
    fprintf( fp, "DIFF = HISTO;\n" );
    for ( j=0; j<nbGauss; j++ ) {
      fprintf( fp, "G%d = %f * exp( -(XHISTO - %f).*(XHISTO - %f)/(2*%f*%f) );\n",
	       j, gauss[3*j+0], gauss[3*j+1], gauss[3*j+1], gauss[3*j+2], gauss[3*j+2] );
      fprintf( fp, "%% plot(XHISTO, G%d, 'g-' );\n", j );
      if ( j==0 ) fprintf( fp, "G = G0;\n;" );
      else        fprintf( fp, "G = G+ G%d;\n;", j );
      fprintf( fp, "DIFF = DIFF - G%d;\n", j );
    }
    fprintf( fp, "plot(XHISTO, G, 'g-' );\n" );
    
    fprintf( fp, "plot(XHISTO, DIFF, 'r-' );\n" );
    fprintf( fp, "axis([0 %d %f %f]);\n", xmax, ymin, ymax );
    fprintf( fp, "grid;\n" );
    fprintf( fp, "hold off;\n" );

    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );
    fprintf( fp, "plot(XHISTO, HISTO, 'b-' );\n" );
    for ( j=0; j<nbGauss; j++ ) {
      fprintf( fp, "plot(XHISTO, G%d, 'g-' );\n", j );
    }
    fprintf( fp, "plot(XHISTO, DIFF, 'r-' );\n" );
    fprintf( fp, "axis([0 %d %f %f]);\n", xmax, ymin, ymax );
    fprintf( fp, "grid;\n" );
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



      else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
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
      else if ( strcmp ( argv[i], "-amp" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -amp...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->amp) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -amp...\n", 0 );
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

  par->amp = -1;
  par->moy = -1;
  par->ect = -1;
}
