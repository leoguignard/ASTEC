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
#include <math.h>
#include <vt_histo.h>

#include <vt_matrices.h>
#include <vt_jointhisto.h>

#include <string.h>







static char *usage = "[image-hmpao] [image-xenon] [-mat %s] [nom-generique]\n\
\t [-hsigma %f] [-xsigma %f] [-hmask %s] [-xmask %s]\n\
\t [-sh %d] [-sx %d] \n\
\t [-vh %f] [-vx %f]\n\
\t [-inv] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t Construit un histogramme conjoint entre 2 images\n\
\t Genere une image inrimage et un fichier de commandes matlab\n\
\t\n\
\t Chaque point de l'image 2 (XENON) est transforme vers l'image 1 (HMPAO)\n\
\t grace a la matrice de transformation. L'histogramme conjoint est calcule\n\
\t soit par PV (partial volume), soit par des ponderations gaussiennes\n\
\t\n\
\t -mat %s : matrice de transformation entre les images XENON et HMPAO\n\
\t           A utiliser telle quelle si c'est du XENON vers le HMPAO\n\
\t           c'est-a-dire si les coefficients de la diagonale sont > 1\n\
\t           Sinon il faut specifier -inv)\n\
\t           S'il n'y a pas de matrice, c'est la transformation identite\n\
\n\
\t -hsigma %f : sigma de la gaussienne pour l'image HMPAO (en mm)\n\
\t -xsigma %f : sigma de la gaussienne pour l'image XENON (en mm)\n\
\t Ces gaussiennes corespondent a des densites de probabilites centrees\n\
\n\
\t -sh %d : seuil pour l'image HMPAO\n\
\t -sx %d : seuil pour l'image XENON\n\
\t -vh %f : taille du voxel pour l'image HMPAO\n\
\t -vx %f : taille du voxel pour l'image XENON\n\
\t\n\
\t -inv : inverse la matrice\n\
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
  vt_image *imageHmpao, *imageXenon, imtmp;
  vt_image *maskHmpao = NULL, *maskXenon = NULL;
  vt_3m mHmpao, mXenon, mConjo;

  
  float ***theHist = (float***)NULL;

  int o;
  int n;
  int xenon, hmpao;

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

  if ( par.names.out[0] == '\0' || par.names.out[0] == '>' ) {
    VT_ErrorParse("unable to write output on stdout, please specify a filename\n", 0);
  }

  /*--- lecture de l'image d'entree ---*/
  imageHmpao = _VT_Inrimage( par.names.in );
  if ( imageHmpao == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read hmpao image\n", 0);
  imageXenon = _VT_Inrimage( par.names.ext );
  if ( imageXenon == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read hmpao image\n", 0);
  
  (void)VT_3m( imageHmpao, (vt_image*)NULL, &mHmpao );
  (void)VT_3m( imageXenon, (vt_image*)NULL, &mXenon );


  maskHmpao = _VT_Inrimage( par.masks.in );
  maskXenon = _VT_Inrimage( par.masks.ext );


  if ( 1 ) {
    fprintf( stderr, " ... image HMPAO : %s\n", par.names.in );
    fprintf( stderr, " ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mHmpao.min, mHmpao.moy, mHmpao.max );
    fprintf( stderr, " ...       taille du point %f x %f x %f\n",
	     imageHmpao->siz.x, imageHmpao->siz.y, imageHmpao->siz.z );
    if ( par.voxelHmpao > 0.1 )
    fprintf( stderr, "           corrigee en %f x %f x %f\n", 
	     par.voxelHmpao, par.voxelHmpao, par.voxelHmpao );
    fprintf( stderr, " ...       sigma = %f\n", par.sigmaHmpao );


    fprintf( stderr, " ... image XENON : %s\n", par.names.ext );
    fprintf( stderr, " ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mXenon.min, mXenon.moy, mXenon.max );
    fprintf( stderr, " ...       taille du point %f x %f x %f\n",
	     imageXenon->siz.x, imageXenon->siz.y, imageXenon->siz.z );
    if ( par.voxelXenon > 0.1 )
    fprintf( stderr, "           corrigee en %f x %f x %f\n", 
	     par.voxelXenon, par.voxelXenon, par.voxelXenon );
    fprintf( stderr, " ...       sigma = %f\n", par.sigmaXenon );

    fprintf( stderr, "     les deux images sont sensees etre recalees" );
    if ( par.matrice[0] != '\0' ) {
      fprintf( stderr, " (matrice de transfo %s)", par.matrice );
    }
    fprintf( stderr, ".\n\n" );
  }

  if ( par.voxelHmpao > 0.1 )
    imageHmpao->siz.x = imageHmpao->siz.y = imageHmpao->siz.z = par.voxelHmpao;
  if ( par.voxelXenon > 0.1 )
    imageXenon->siz.x = imageXenon->siz.y = imageXenon->siz.z = par.voxelXenon;
  
  


  

  



  /* allocation des images de calculs
     - histogramme conjoint
     - image des moyennes (pour une intensite hmpao fixee)
     - image des medianes (pour une intensite hmpao fixee)
  */

  VT_Image( &imtmp );
  VT_InitImage( &imtmp, "", (int)(mXenon.max+0.5)+1, (int)(mHmpao.max+0.5)+1,
		1, FLOAT );
  sprintf( imtmp.name, "%s.inr", par.names.out );
  (void)VT_AllocImage( &imtmp );
  theHist = (float***)imtmp.array;


  
  /* calcul de l'histogramme conjoint
   */
  if ( par.matrice[0] == '\0' ) {

    if ( 1 ) {
      fprintf( stderr, " ... calcul sans transformation geometrique\n" );
      fprintf( stderr, "     les points ont meme dimension dans les 2 images\n" );
      fprintf( stderr, "     sigma HMPAO = %f mm / %f = %f pixel\n", par.sigmaHmpao,
	       imageHmpao->siz.x, par.sigmaHmpao/imageHmpao->siz.x );
      fprintf( stderr, "     sigma XENON = %f mm / %f = %f pixel\n\n", par.sigmaXenon,
	       imageXenon->siz.x, par.sigmaXenon/imageXenon->siz.x );
    }
    (void)ComputeJointHistoWithoutTrsf( imageHmpao, imageXenon, &imtmp,
			     par.sigmaHmpao/imageHmpao->siz.x, 
			     par.sigmaXenon/imageXenon->siz.x, par.seuilHmpao );

  } else {
    
    if ( Read4x4Matrix( par.matrice, mat ) != 1 ) {
      VT_ErrorParse("unable to read matrice\n", 0);
    }
    if ( par.names.inv != 0 ) {
      double tmp[16];
      int rang;
      if ( (rang=Inverse4x4Matrix( mat, tmp )) != 4 ) {
	fprintf( stderr, "Warning: la matrice %s est de rang %d.\n",
		 par.matrice, rang );
      }
      for ( n=0; n<16; n++ ) mat[n] = tmp[n];
    } 
    if ( 0 ) {
      (void)ComputeJointHistoWithTrsf( imageHmpao, imageXenon, &imtmp,
				       mat, par.sigmaHmpao, par.sigmaXenon,
				       par.seuilHmpao, par.seuilXenon );
    }
    else {
      (void)ComputeJointHistoWithTrsfAndMask( imageHmpao, imageXenon, maskHmpao, maskXenon, &imtmp,
				       mat, par.sigmaHmpao, par.sigmaXenon,
				       par.seuilHmpao, par.seuilXenon );
    }
  }


  











  

  



  
  /* seuillage additionnel
   */
  
  for ( hmpao=0; hmpao<imtmp.dim.y; hmpao++ )
  for ( xenon=0; xenon<imtmp.dim.x; xenon++ ) {
    if ( xenon < par.seuilXenon || hmpao < par.seuilHmpao) 
      theHist[0][hmpao][xenon] = 0;
  }

  (void)VT_3m( &imtmp, (vt_image*)NULL, &mConjo );
  if ( 1 ) {
    fprintf(stderr," ... %s (apres seuillage) : %17.6f %17.6f %17.6f\n", 
	    imtmp.name, mConjo.min, mConjo.moy, mConjo.max );
  }
  











  /* ecriture fichier matlab
   */
  {
    char name[256];
    FILE *f, *fopen();
    vt_histo histo;
    int fd, i;

    char strtitle[256];
    char strmatrx[256];
    char strsigma[256];

    
    sprintf( name, "%s.raw", par.names.out );
    fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    if( write( fd, imtmp.buf, (imtmp.dim.x*imtmp.dim.y)*sizeof( float ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }
    close( fd );

    sprintf( name, "%s.m", par.names.out );
    f = fopen( name, "w" );

    
    fprintf( f, "echo off\n" );
    fprintf( f, "\n" );
    fprintf( f, "%% ... image HMPAO : %s\n", par.names.in );
    fprintf( f, "%% ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mHmpao.min, mHmpao.moy, mHmpao.max );
    fprintf( f, "%% ...       taille du point %f x %f x %f\n",
	     imageHmpao->siz.x, imageHmpao->siz.y, imageHmpao->siz.z );
    fprintf( f, "%%           seuillee au-dessus (au sens large) de %d\n",
	     par.seuilHmpao );
    fprintf( f, "\n" );
    fprintf( f, "%% ... image XENON : %s\n", par.names.ext );
    fprintf( f, "%% ...       min, moy, max = %17.6f %17.6f %17.6f\n", 
	     mXenon.min, mXenon.moy, mXenon.max );
    fprintf( f, "%% ...       taille du point %f x %f x %f\n",
	     imageXenon->siz.x, imageXenon->siz.y, imageXenon->siz.z );
    fprintf( f, "%%           seuillee au-dessus (au sens large) de %d\n",
	     par.seuilXenon );
    fprintf( f, "\n\n" );



    VT_Histo( &histo );
    histo.type = imageHmpao->type;
    (void)VT_AllocHisto( &histo );
    (void)VT_ComputeHisto( &histo, imageHmpao );

    o = 0;
    switch( histo.type ) {
    default : break;
    case SSHORT : o = -32768; break;
    }

    fprintf( f, " HISTOHMPAO = [" );
    for ( i = 0; i<=(int)(mHmpao.max+0.5); i++ )
      fprintf( f, " %ld", histo.buf[i-o] );
    fprintf( f, " ];\n\n" );

    fprintf( f, "echo on\n\n" );
    fprintf( f, "%%%% - Pour voir l'histogramme de l'image HMPAO\n" );
    fprintf( f, "%%%% - taper les lignes suivantes:\n\n" );
    fprintf( f, "%% figure;\n" );
    fprintf( f, "%% hold on;\n" );
    fprintf( f, "%% plot(0:%f, HISTOHMPAO);\n", mHmpao.max );
    fprintf( f, "%% axis([0 %f 0 3000]);\n", mHmpao.max );
    fprintf( f, "%% xlabel('Intensite HMPAO');\n" );
    fprintf( f, "%% ylabel('#points');\n" );
    fprintf( f, "%% t=title('Histogramme image Hmpao = %s');\n", imageHmpao->name );
    fprintf( f, "%% set(t,'FontWeight','bold');\n" );
    fprintf( f, "%% hold off;\n" );
    fprintf( f, "echo off\n" );
    VT_FreeHisto( &histo );


    VT_Histo( &histo );
    histo.type = imageXenon->type;
    (void)VT_AllocHisto( &histo );
    (void)VT_ComputeHisto( &histo, imageXenon );

    o = 0;
    switch( histo.type ) {
    default : break;
    case SSHORT : o = -32768; break;
    }

    fprintf( f, " HISTOXENON = [" );
    for ( i = 0; i<=(int)(mXenon.max+0.5); i++ )
      fprintf( f, " %ld", histo.buf[i-o] );
    fprintf( f, " ];\n\n" );

    fprintf( f, "echo on\n\n" );
    fprintf( f, "%%%% - Pour voir l'histogramme de l'image XENON\n" );
    fprintf( f, "%%%% - taper les lignes suivantes:\n\n" );
    fprintf( f, "%% figure;\n" );
    fprintf( f, "%% hold on;\n" );
    fprintf( f, "%% plot(0:%f, HISTOXENON);\n", mXenon.max );
    fprintf( f, "%% axis([0 %f 0 200]);\n", mHmpao.max );
    fprintf( f, "%% xlabel('Intensite XENON');\n" );
    fprintf( f, "%% ylabel('#points');\n" );
    fprintf( f, "%% t=title('Histogramme image Xenon = %s');\n", imageXenon->name );
    fprintf( f, "%% set(t,'FontWeight','bold');\n" );
    fprintf( f, "%% hold off;\n" );
    fprintf( f, "\n" );
    fprintf( f, "echo off\n" );
    VT_FreeHisto( &histo );



    {
      int i, k;
      k = strlen( par.names.out );
      for ( i = k-1; i >= 0 && par.names.out[i] != '/' ; i-- )
	;
      fprintf( f, " fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
    }
    fprintf( f, " [HCREAD, HCNBELTS] = fread( fid, [%d,%d], 'float%lu' );\n", 
	     1+(int)(mXenon.max+0.5), 1+(int)(mHmpao.max+0.5), 8*sizeof( float ) );
    fprintf( f, " fclose( fid );\n" );
    fprintf( f, "\n" );
    fprintf( f, "HC = HCREAD';\n\n" );

    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "%% mesh(HC);\n " );
    fprintf( f, "surf(0:%d,0:%d,HC);\n ", (int)(mXenon.max+0.5), (int)(mHmpao.max+0.5) );
    fprintf( f, "axis([0 %f 0 %f 0 3.0])\n", mXenon.max, mHmpao.max );
    fprintf( f, "caxis([0 2.0])\n" );
    fprintf( f, "shading interp\n " );

    fprintf( f, "xlabel('X = Xenon\\newline%s\\newline');\n", 
	     imageXenon->name );
    fprintf( f, "ylabel('Y = Hmpao\\newline%s\\newline');\n", 
	     imageHmpao->name );
    
    fprintf( f, "grid;\n" );

    strtitle[0] = strmatrx[0] = strsigma[0] = '\0';
    sprintf( strtitle, "Histogramme conjoint" );
    if ( par.matrice[0] != '\0' )
      sprintf( strmatrx, "Matrice de transformation = %s", par.matrice );
    if ( par.sigmaXenon > 0.1 || par.sigmaHmpao > 0.1 )
      sprintf( strsigma, "Sigma Hmpao = %7.4f - Sigma Xenon = %7.4f", par.sigmaHmpao, par.sigmaXenon );
    
    if ( strmatrx[0] == '\0' ) {
      if ( strsigma[0] == '\0' ) {
	 fprintf( f, "t=title('%s');\n", strtitle );
      } else {
	fprintf( f, "t=title('%s\\newline%s\\newline');\n", strtitle, strsigma );
      }
    } else {
      if ( strsigma[0] == '\0' ) {
	fprintf( f, "t=title('%s\\newline%s\\newline');\n", strtitle, strmatrx );
      } else {
	fprintf( f, "t=title('%s\\newline%s\\newline%s\\newline');\n", strtitle, strmatrx ,strsigma );
      }
    }
    fprintf( f, "set(t,'FontWeight','bold');\n" );
    fprintf( f, "set(t,'HorizontalAlignment','center');\n" );
    fprintf( f, "set(t,'VerticalAlignment','middle');\n" );
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
