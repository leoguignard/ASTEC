/*************************************************************************
 * hmpao-xenon.c - 
 *
 * $Id: hmpao-xenon.c,v 1.5 2000/04/13 10:54:56 greg Exp $
 *
 * Copyright (c) INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Thu Mar  2 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

/*
../../HMPAO-XENON/HMPAO/INRIA/fessy.inr 
../../HMPAO-XENON/XENON/INRIA/FESSYXe.inr test
-mat ../../HMPAO-XENON/AROCHE/fessy.trsf -inv -sh 90 -xsigma 14.0

/1/greg/HMPAO-XENON/HMPAO/INRIA/fessy.inr /1/greg/HMPAO-XENON/XENON/INRIA/FESSYXe.inr -mat /1/greg/HMPAO-XENON/AROCHE/fessy.trsf test11 -sh 90 -xsigma 8 -hsigma 2

*/



#include <vt_common.h>
#include <vt_3m.h>
#include <math.h>
#include <vt_histo.h>

#include <vt_matrices.h>
#include <vt_jointhisto.h>

#include <string.h>







static char *usage = "[image-1] [image-2] [nom-generique]\n\
\t [-inv] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t Construit un histogramme conjoint entre 2 images\n\
\t Genere une image inrimage et un fichier de commandes matlab\n\
\t\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];










typedef struct local_par {
  vt_names names;
  vt_names masks;

  int seuilXenon;
  int seuilHmpao;
  int type;
  float sigmaXenon;
  float sigmaHmpao;
  float voxelXenon;
  float voxelHmpao;
  char matrice[STRINGLENGTH];
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








#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image1, *image2, imtmp;
  vt_3m mmm1, mmm2, mConjo;

  
  float ***theHist = (float***)NULL;

  int n;

  double mat[16];
  
  





  /*  0  1  2  3
      4  5  6  7
      8  9 10 11
     12 13 14 15 
     */
  for ( n=0; n<16; n++ ) mat[n] = 0.0;
  mat[0] = mat[5] = mat[10] = mat[15] = 1.0;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /*--- lecture de l'image d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read image #1\n", 0);
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read image #2\n", 0);
  
  (void)VT_3m( image1, (vt_image*)NULL, &mmm1 );
  (void)VT_3m( image2, (vt_image*)NULL, &mmm2 );

  VT_Image( &imtmp );
  VT_InitImage( &imtmp, "", (int)(mmm2.max+0.5)+1, (int)(mmm1.max+0.5)+1,
		1, FLOAT );
  sprintf( imtmp.name, "%s.inr", par.names.out );
  (void)VT_AllocImage( &imtmp );
  theHist = (float***)imtmp.array;


  
  fprintf( stdout, "'%s' in [%f %f]\n", image1->name, mmm1.min, mmm1.max );
  fprintf( stdout, "'%s' in [%f %f]\n", image2->name, mmm2.min, mmm2.max );

  /* calcul de l'histogramme conjoint
   */
  
  (void)ComputeJointHistoWithoutTrsf( image1, image2, &imtmp,
				      0.0, 0.0, mmm1.min - 1 );
  (void)VT_3m( &imtmp, (vt_image*)NULL, &mConjo );

  {
    int dimx = (int)(mmm2.max+0.5)+1;
    int dimy = (int)(mmm1.max+0.5)+1;
    double *sumx = (double*)malloc( dimx * sizeof(double) );
    double *sumy = (double*)malloc( dimy * sizeof(double) );
    double sum = 0;
    int x, y;
    float *jhisto = (float*)imtmp.buf;

    double im=0.0;

    for ( x=0; x<dimx; x++ ) sumx[x] = 0.0;
    for ( y=0; y<dimy; y++ ) sumy[y] = 0.0;

    for ( y=0; y<dimy; y++ )
    for ( x=0; x<dimx; x++ ) {
      sumx[x] += jhisto[y*dimx+x];
      sumy[y] += jhisto[y*dimx+x];
      sum += jhisto[y*dimx+x];
    }

    for ( x=0; x<dimx; x++ ) sumx[x] /= sum;
    for ( y=0; y<dimy; y++ ) sumy[y] /= sum;
    
    /* information mutuelle
     */
    for ( y=0; y<dimy; y++ )
    for ( x=0; x<dimx; x++ ) {
      if ( jhisto[y*dimx+x] <= 0 ) continue;
      im += jhisto[y*dimx+x] / sum * log ( jhisto[y*dimx+x] / (sumx[x] * sumy[y] * sum) );
    }
    printf( "information mutuelle = %g\n", im );

  }



  /* ecriture fichier matlab
   */
  if ( par.names.out[0] != '\0' && par.names.out[0] != '>' ) {
    char name[256];
    FILE *f, *fopen();
    int fd;

    
    sprintf( name, "%s.raw", par.names.out );
    fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    if ( write( fd, imtmp.buf, (imtmp.dim.x*imtmp.dim.y)*sizeof( float ) ) == -1 ) {
      fprintf( stderr, "%s: error when reading\n", program );
    }
    close( fd );

    sprintf( name, "%s.m", par.names.out );
    f = fopen( name, "w" );

    
    fprintf( f, "echo off\n" );
    fprintf( f, "\n" );
    fprintf( f, "%% ... image #1 : %s\n", par.names.in );
    fprintf( f, "%% ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mmm1.min, mmm1.moy, mmm1.max );
    fprintf( f, "\n" );
    fprintf( f, "%% ... image #2 : %s\n", par.names.ext );
    fprintf( f, "%% ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mmm2.min, mmm2.moy, mmm2.max );
    fprintf( f, "\n\n" );



    {
      int i, k;
      k = strlen( par.names.out );
      for ( i = k-1; i >= 0 && par.names.out[i] != '/' ; i-- )
	;
      fprintf( f, " fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
    }
    fprintf( f, " [HCREAD, HCNBELTS] = fread( fid, [%d,%d], 'float%lu' );\n", 
	     1+(int)(mmm2.max+0.5), 1+(int)(mmm1.max+0.5), 8*sizeof( float ) );
    fprintf( f, " fclose( fid );\n" );
    fprintf( f, "\n" );
    fprintf( f, "HC = HCREAD';\n\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "%% mesh(HC);\n " );
    fprintf( f, "surf(0:%d,0:%d,HC);\n ", (int)(mmm2.max+0.5), (int)(mmm1.max+0.5) );
    fprintf( f, "axis([0 %f 0 %f 0 %f])\n", mmm2.max, mmm1.max, mConjo.max/100);
    fprintf( f, "caxis([0 %f])\n", mConjo.max/200 );
    fprintf( f, "shading interp\n " );
    
    fprintf( f, "xlabel('X = %s\\newline');\n", 
	     image1->name );
    fprintf( f, "ylabel('Y = %s\\newline');\n", 
	     image2->name );
    
    fprintf( f, "grid;\n" );

    fprintf( f, "view(130,60);\n" );
    fprintf( f, "hold off;\n" );
  }


  /*--- ecriture de l'image resultat ---*/

  if ( VT_WriteInrimage( &imtmp ) == -1 ) {
    VT_ErrorParse("unable to write output image\n", 0);
  }
  




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




      else if ( strcmp ( argv[i], "-mat" ) == 0 ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -mat...\n", 0 );
	strncpy( par->matrice, argv[i], STRINGLENGTH );  
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




      else if ( strcmp ( argv[i], "-vh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -vh...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->voxelHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vh...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-vx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -vx...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->voxelXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vx...\n", 0 );
      }




      else if ( strcmp ( argv[i], "-xsigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -xsigma...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->sigmaXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -xsigma...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-hsigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -hsigma...\n", 0 );
	status = sscanf( argv[i],"%f",&(par->sigmaHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -hsigma...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-hmask" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -hmask...\n", 0 );
	strncpy( par->masks.in, argv[i], STRINGLENGTH );
      }
      else if ( strcmp ( argv[i], "-xmask" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -xmask...\n", 0 );
	strncpy( par->masks.ext, argv[i], STRINGLENGTH );
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
  VT_Names( &(par->masks) );
  par->type = TYPE_UNKNOWN;
  par->seuilHmpao = 0;
  par->seuilXenon = 0;
  par->sigmaHmpao = 0;
  par->sigmaXenon = 0;
  par->voxelHmpao = 0;
  par->voxelXenon = 0;
  par->matrice[0] ='\0';
}
