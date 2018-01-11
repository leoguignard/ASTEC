/*************************************************************************
 * hmpao-xenon-fit.c - 
 *
 * $Id: hmpao-xenon-fit.c,v 1.16 2000/06/29 16:48:21 greg Exp $
 *
 * Copyright (c) INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Mon Mar 13 18:42:33 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

/*

*/



#include <vt_common.h>
#include <math.h>
#include <vt_histo.h>

#include <recbuffer.h>
#include <vt_levenberg.h>
#include <vt_statsutil.h>

#include <string.h>




static char *usage = "[histo-conjoint] [nom-generique]\n\
\t [-sh %d] [-sx %d] [-mh %d] [-mx %d]\n\
\t [-psh %lf] [-psx %lf] [-pmh %lf] [-pmx %lf]\n\
\t [-inv] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t -sh %d : seuil pour l'image HMPAO\n\
\t -sx %d : seuil pour l'image XENON\n\
\t -mh %d : max pour l'image HMPAO\n\
\t -mx %d : max pour l'image XENON\n\
\t -psh %lf : (pourcentage) seuil pour l'image HMPAO\n\
\t -psx %lf : (pourcentage) seuil pour l'image XENON\n\
\t -pmh %lf : (pourcentage) max pour l'image HMPAO\n\
\t -pmx %lf : (pourcentage) max pour l'image XENON\n\
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




extern double exp    (double x);













typedef struct local_par {
  vt_names names;
  int seuilXenon;
  int seuilHmpao;
  int maxXenon;
  int maxHmpao;

  double pourSeuilXenon;
  double pourSeuilHmpao;
  double pourMaxXenon;
  double pourMaxHmpao;

  int type;
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

int _UpdateLowThresholdXenon( vt_image *imHisto,
			      int old,
			      double p );
int _UpdateHighThresholdXenon( vt_image *imHisto,
			      int old,
			      double p );

void _InitValues( double *par,
		  double xenon1, double hmpao1,
		  double xenon2, double hmpao2 );



typedef struct {
  double a;
  double b;
} typeFtnCalib;






static void _PrintFonctionCalib( typeFtnCalib *ftn, char *s );
static int _InitFonctionCalib( typePoint *theXenon, typePoint *theHmpao,
			       int minx,
			       int maxx,
			       int miny,
			       int maxy,
			       typeFtnCalib *ftnInit,
			       int type );

static int _EstimeFonctionCalib( vt_image *theIm,
				 typePoint *theXenon, typePoint *theHmpao,
				 int minx,
				 int maxx,
				 int miny,
				 int maxy,
				 typeFtnCalib *ftnFinal,
				 int type );

#ifdef _UNUSED_
static int _EstimeFonctionCalibGradient( vt_image *theIm,
					 typePoint *theXenon, typePoint *theHmpao,
					 int minx,
					 int maxx,
					 int miny,
					 int maxy,
					 typeFtnCalib *ftnFinal,
					 int type );
#endif

#ifdef _UNUSED_
static void _TraceFonctionCalib( vt_image *theIm,
				 typeFtnCalib *ftn, int z, float v );
#endif




#define _XENON_MED_ 2
#define _XENON_MOY_ 3
#define _HMPAO_MED_ 5
#define _HMPAO_MOY_ 6

/*
#define _INIT_XENON_MOY_ 7
#define _INIT_HMPAO_MOY_ 8
#define _OPT_XMOY_ 9
#define _OPT_HMOY_ 10
#define _OPT_XHMOY_ 11

#define _XENON_ONEDERIV_ 12
#define _XENON_TWODERIV_ 13
#define _XENON_CURVATURE_ 14
#define _XENON_MAXCURVATURE_ 15
*/

#define NBSLICES 16

#define _IMAGE_XENON_MOY_ -3






static int _BORDER_ = 10;

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *imHisto, imCurves;
  
  float ***theHisto = (float***)NULL;
  float ***theCurve = (float***)NULL;

  
  int xenon, hmpao;

  typePoint *theXenon = (typePoint*)NULL;
  typePoint *theHmpao = (typePoint*)NULL;

  typeFtnCalib ftnInitXenon;
  typeFtnCalib ftnFinalXenon;
  typeFtnCalib ftnFinal2Xenon;

  int minXenon, maxXenon;
  int minHmpao, maxHmpao;

  double *theX, *theY, *theC, *theS;
  double *theAllocatedBuffer = (double*)NULL;
  int i, j, k, length;
  double thePar[2];



  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  if ( par.names.out[0] == '\0' || par.names.out[0] == '>' ) {
    VT_ErrorParse("unable to write output on stdout, please specify a filename\n", 0);
  }

  /*--- lecture de l'image histogramme Conjoint
    X = intensite xenon
    Y = intensite hmpao
  */
  
  imHisto = _VT_Inrimage( par.names.in );
  if ( imHisto == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read histogramme image\n", 0);
  if ( imHisto->type != FLOAT ) 
    VT_ErrorParse("bad type for histogramme image\n", 0);
  theHisto = (float***)imHisto->array;






  if ( par.maxXenon > imHisto->dim.x-1 || par.maxXenon < 0 ) 
    par.maxXenon = imHisto->dim.x-1;
  
  if ( par.maxHmpao > imHisto->dim.y-1 || par.maxHmpao < 0 ) 
    par.maxHmpao = imHisto->dim.y-1;
  
  minXenon = ( par.seuilXenon - _BORDER_ > 0 ) ? par.seuilXenon - _BORDER_ : 0;
  maxXenon = ( par.maxXenon + _BORDER_ < imHisto->dim.x ) ? par.maxXenon + _BORDER_ : imHisto->dim.x-1;

  minHmpao = ( par.seuilHmpao - _BORDER_ > 0 ) ? par.seuilHmpao - _BORDER_ : 0;
  maxHmpao = ( par.maxHmpao + _BORDER_ < imHisto->dim.y ) ? par.maxHmpao + _BORDER_ : imHisto->dim.y-1;
  

  /* calcul des tableaux Xenon et Hmpao
   */

  theXenon = (typePoint*)malloc( imHisto->dim.x * sizeof(typePoint) );
  theHmpao = (typePoint*)malloc( imHisto->dim.y * sizeof(typePoint) );
  
  minXenon = _UpdateLowThresholdXenon( imHisto, minXenon, par.pourSeuilXenon );
  maxXenon = _UpdateHighThresholdXenon( imHisto, maxXenon, par.pourMaxXenon );

  printf( " ... Seuils Xenon = [%d %d] ... Seuils Hmpao = [%d %d]\n",
	  minXenon, maxXenon, minHmpao, maxHmpao );

  _ComputeStatsValues( imHisto, theXenon, theHmpao, minXenon, maxXenon, minHmpao, maxHmpao );

  
  /* allocation pour la minimisation
   */

  length = imHisto->dim.y;
  if ( length < imHisto->dim.x ) length = imHisto->dim.x;

  theAllocatedBuffer = (double*)malloc( 4 * length * sizeof(double) );
  if ( theAllocatedBuffer == (double*)NULL ) {
    VT_ErrorParse("unable to allocate buffer\n", 0);
  }
  
  theX = theY = theC = theS = theAllocatedBuffer;
  theY += length;
  theC += 2*length;
  theS += 3*length;



  /* allocation des images de calculs
     - histogramme conjoint
     - image des moyennes (pour une intensite hmpao fixee)
     - image des medianes (pour une intensite hmpao fixee)
  */

  VT_Image( &imCurves );
  VT_InitImage( &imCurves, "", imHisto->dim.x, imHisto->dim.y, NBSLICES, FLOAT );
  sprintf( imCurves.name, "%s.inr", par.names.out );
  (void)VT_AllocImage( &imCurves );
  theCurve = (float***)imCurves.array;


  for ( k=0; k<imCurves.dim.z; k++ )
  for ( j=0; j<imCurves.dim.y; j++ )
  for ( i=0; i<imCurves.dim.x; i++ )
    theCurve[k][j][i] = 0.0;


  



  











  _InitFonctionCalib( theXenon, theHmpao, par.seuilXenon, par.maxXenon,
		      par.seuilHmpao, par.maxHmpao,
		      &ftnInitXenon, _XENON_MOY_ );
  if ( 0 )
    _PrintFonctionCalib( &ftnInitXenon, 
			 "initialisation sur les moyennes des xenons" );
  
  ftnFinalXenon = ftnInitXenon;
  _EstimeFonctionCalib( imHisto, theXenon, theHmpao, par.seuilXenon, par.maxXenon,
			par.seuilHmpao, par.maxHmpao, &ftnFinalXenon, _XENON_MOY_ );
  if ( 0 )
    _PrintFonctionCalib( &ftnFinalXenon, " optimisation sur les moyennes des xenons" );



  thePar[0] = ftnInitXenon.a;
  thePar[1] = ftnInitXenon.b;

  /* for ( length=0, xenon=par.seuilXenon; xenon<=par.maxXenon; xenon++ ) { */

  fprintf( stdout, "\n" );

  for ( length=0, xenon=minXenon; xenon<= maxXenon; xenon++ ) {
    if ( theXenon[xenon].sum0 <= _MINVAL_ ) continue;
    /* pas d'amplitude */
    if ( theXenon[xenon].nsgauss[0] <= 0.0 ) continue;
    /* pas de mode */
    if ( theXenon[xenon].nsgauss[1] <= 1.0 ) continue;
    /* ecart-type non valide */
    if ( theXenon[xenon].nsgauss[2] < 0.0 || theXenon[xenon].nsgauss[2] > 1000.0 ) continue;
    if ( theXenon[xenon].nsgauss[3] < 0.0 || theXenon[xenon].nsgauss[3] > 1000.0 ) continue;

    fprintf( stdout, " XENON #%2d : X = %3d AMP = %f MOY = %f EC1 = %f EC2 = %f\n",
	     length, xenon, theXenon[xenon].nsgauss[0], theXenon[xenon].nsgauss[1],
	     theXenon[xenon].nsgauss[2], theXenon[xenon].nsgauss[3] );
    theX[length] = xenon;
    theY[length] = theXenon[xenon].nsgauss[1]; /* mode */

    /* avant */
    theC[length] = theXenon[xenon].sum0;
    theS[length] = theXenon[xenon].nsgauss[2]; /* premier ecart-type */
    if ( theS[length] < theXenon[xenon].nsgauss[3] ) theS[length] = theXenon[xenon].nsgauss[3];

    /* on donne le meme poids a tous les points */
    theC[length] = 1.0;
    theS[length] = 1.0;

    /* on garde la variance qd meme */
    theC[length] = 1.0;
    theS[length] = theXenon[xenon].nsgauss[2]; /* premier ecart-type */
    if ( theS[length] < theXenon[xenon].nsgauss[3] ) theS[length] = theXenon[xenon].nsgauss[3];

    

    length++;
  }
  
  {
    double x1, x2, h1, h2;
    int l1 = (int)(length/3.0);
    int l2 = (int)((2.0 * length)/3.0);
    
    x1 = theX[l1];
    h1 = theY[l1];

    x2 = theX[l2];
    h2 = theY[l2];


    fprintf( stdout, "\n" );
    fprintf( stdout, " points #%d [X=%f H=%f] et #%d [X=%f H=%f]\n", l1,  x1, h1, l2, x2 ,h2 );
    _InitValues( thePar, x1, h1, x2, h2 );

    ftnFinal2Xenon.a = thePar[0];
    ftnFinal2Xenon.b = thePar[1];
    _PrintFonctionCalib( &ftnFinal2Xenon, " initialisation" );
  }
  /*
  VT_Levenberg_verbose();
  VT_Levenberg_verbose();
  */

  /*
    le critere est sum ( theC / (theS * theS) [ theY - f( theX ) ]^2 )
   */

  if ( VT_Modeling1DDataWithLevenberg( theX, DOUBLE, theY, DOUBLE, theC, DOUBLE, theS, DOUBLE, length,
				       thePar, 2, _LassenFunction ) == 1 ) {
    ftnFinal2Xenon.a = thePar[0];
    ftnFinal2Xenon.b = thePar[1];
    _PrintFonctionCalib( &ftnFinal2Xenon, " optimisation sur les modes des xenons" );
  }
  


  
  


  




  /* ecriture fichier matlab
   */
  {
    char name[256];
    FILE *f, *fopen();
    int fd;


    sprintf( name, "%s.raw", par.names.out );
    fd = creat( name, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    if ( write( fd, imHisto->buf, (imHisto->dim.x*imHisto->dim.y)*sizeof( float ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }

    for ( hmpao=0; hmpao<imHisto->dim.y; hmpao++ )
    for ( xenon=0; xenon<imHisto->dim.x; xenon++ ) { 
      if ( xenon < par.seuilXenon || hmpao < par.seuilHmpao || 
	   xenon > par.maxXenon || hmpao > par.maxHmpao ) 
	theHisto[0][hmpao][xenon] = 0.0;
    }

    if ( write( fd, imHisto->buf, (imHisto->dim.x*imHisto->dim.y)*sizeof( float ) ) == -1 ) {
      fprintf( stderr, "%s: error when writing\n", program );
    }
    close( fd );



    sprintf( name, "%s.m", par.names.out );
    f = fopen( name, "w" );


    fprintf( f, "echo off\n" );
    fprintf( f, "\n" );
    fprintf( f, "%% ... histogramme conjoint : %s\n", imHisto->name );
    fprintf( f, "%%     seuil HMPAO au-dessus (au sens large) de %d\n",
	     par.seuilHmpao );
    fprintf( f, "%%     seuil XENON au-dessus (au sens large) de %d\n",
	     par.seuilXenon );
    fprintf( f, "\n" );

    fprintf( f, "\n\n" );


    
    fprintf( f, "%% cree la matrice [ 0 1 2 3 ..\n" );
    fprintf( f, "%%                   0 1 2 3 ..\n" );
    fprintf( f, "%%                   0 1 2 3 .. ]\n" );
    fprintf( f, "%% XX = repmat((0:%lu),%lu,1);\n", imHisto->dim.x-1, imHisto->dim.y );
    fprintf( f, "%% YY = repmat((0:%lu)',1,%lu);\n\n", imHisto->dim.y-1, imHisto->dim.x );


    fprintf( f, "X = 0:%lu;\n", imHisto->dim.x-1 );
    fprintf( f, "Y = 0:%lu;\n\n", imHisto->dim.y-1 );


    k = strlen( par.names.out );
    for ( i = k-1; i >= 0 && par.names.out[i] != '/' ; i-- )
      ;
    fprintf( f, " fid = fopen('%s.raw', 'r' );\n", &(par.names.out[i+1]) );
    fprintf( f, " [HCREAD, HCNBELTS] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	     imHisto->dim.x, imHisto->dim.y, 8*sizeof( float ) );
    fprintf( f, " [HCSREAD, HCSNBELTS] = fread( fid, [%lu,%lu], 'float%lu' );\n", 
	     imHisto->dim.x, imHisto->dim.y, 8*sizeof( float ) );
    fprintf( f, " fclose( fid );\n" );
    fprintf( f, "\n" );
    fprintf( f, "HC = HCREAD';\n\n" );
    fprintf( f, "HCS = HCSREAD';\n\n" );

    
    
    fprintf( f, " INDXENONMAX = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].indMax > 1 )
	fprintf( f, " %d", xenon );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALXENONMAX = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].indMax > 1 )
	fprintf( f, " %d", theXenon[xenon].indMax );
    }
    fprintf( f, " ];\n\n" );


    
    fprintf( f, " INDXENON = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].sum0 > _MINVAL_ )
	fprintf( f, " %d", xenon );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALXENONMOY = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].sum0 > _MINVAL_ )
	fprintf( f, " %.2f", theXenon[xenon].moy );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALXENONMED = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].sum0 > _MINVAL_ )
	fprintf( f, " %.2f", theXenon[xenon].med );
    }
    fprintf( f, " ];\n\n" );

    
    fprintf( f, " INDXENONMODE = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].nsgauss[1] > 1 )
	fprintf( f, " %d", xenon );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALXENONMODE = [" );
    for ( xenon = 0; xenon < imHisto->dim.x; xenon++ ) {
      if ( theXenon[xenon].nsgauss[1] > 1 )
	fprintf( f, " %.2f", theXenon[xenon].nsgauss[1] );
    }
    fprintf( f, " ];\n\n" );

    

    
    /*
    fprintf( f, " INDHMPAO = [" );
    for ( hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
      if ( theHmpao[hmpao].sum0 > _MINVAL_ )
	fprintf( f, " %d", hmpao );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALHMPAOMOY = [" );
    for ( hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
      if ( theHmpao[hmpao].sum0 > _MINVAL_ )
	fprintf( f, " %.2f", theHmpao[hmpao].moy );
    }
    fprintf( f, " ];\n" );
    fprintf( f, " VALHMPAOMED = [" );
    for ( hmpao = 0; hmpao < imHisto->dim.y; hmpao++ ) {
      if ( theHmpao[hmpao].sum0 > _MINVAL_ )
	fprintf( f, " %.2f", theHmpao[hmpao].med );
    }
    fprintf( f, " ];\n\n" );
    */

    


    

    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );

    fprintf( f, "pcolor(0:%lu,0:%lu,HC);\n", imHisto->dim.x-1, imHisto->dim.y-1 );
    fprintf( f, "shading interp\n" );
    fprintf( f, "caxis([0 1.5]);\n\n" );
    
    fprintf( f, "xlabel('XENON');\n" );
    fprintf( f, "ylabel('HMPAO');\n" );
    fprintf( f, "axis([0 %lu 0 %lu]);\n\n", imHisto->dim.x-1, imHisto->dim.y-1 );
    
    fprintf( f, "plot(INDXENON, VALXENONMOY, 'b-' );\n" );
    fprintf( f, "plot(INDXENON, VALXENONMED, 'b--' );\n\n" );

    /*
    fprintf( f, "plot(VALHMPAOMOY, INDHMPAO, 'k-'  );\n" );
    fprintf( f, "plot(VALHMPAOMED, INDHMPAO, 'k--'  );\n\n" );
    fprintf( f, "legend('Moyenne', 'Mediane', 'Moyenne', 'Mediane'" );
    */

    fprintf( f, "legend('Moyenne', 'Mediane'" );

    if ( minXenon > imHisto->dim.x-1 - par.maxXenon )
      fprintf( f, ", 2" );
    else 
      fprintf( f, ", 1" );
    fprintf( f, ");\n" );
    fprintf( f, "hold off;\n\n" );



    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );

    fprintf( f, "pcolor(0:%lu,0:%lu,HCS);\n", imHisto->dim.x-1, imHisto->dim.y-1 );
    fprintf( f, "shading interp\n" );
    fprintf( f, "caxis([0 1.5]);\n" );
    
    fprintf( f, "xlabel('XENON');\n" );
    fprintf( f, "ylabel('HMPAO');\n" );
    fprintf( f, "axis([0 %lu 0 %lu]);\n", imHisto->dim.x-1, imHisto->dim.y-1 );
    
    fprintf( f, "plot(INDXENON, VALXENONMOY, 'b-' );\n" );
    fprintf( f, "plot(INDXENON, VALXENONMED, 'b--' );\n\n" );
    fprintf( f, "plot(INDXENONMAX, VALXENONMAX, 'k--' );\n" );
    fprintf( f, "plot(INDXENONMODE, VALXENONMODE, 'k-' );\n\n" );
    
    fprintf( f, "X=0:%lu;\n", imHisto->dim.x-1 );
    fprintf( f, "hmInitXenon = (%f * X) ./ (%f + X);\n", 
	     ftnInitXenon.a, ftnInitXenon.b );
    fprintf( f, "plot(X, hmInitXenon, 'r:');\n" );

    fprintf( f, "hmFinalXenon = (%f * X) ./ (%f + X);\n", 
	     ftnFinalXenon.a, ftnFinalXenon.b );
    fprintf( f, "plot(X, hmFinalXenon, 'r--');\n" );

    fprintf( f, "hmFinal2Xenon = (%f * X) ./ (%f + X);\n", 
	     ftnFinal2Xenon.a, ftnFinal2Xenon.b );
    fprintf( f, "plot(X, hmFinal2Xenon, 'r-');\n" );


    fprintf( f, "legend('Moyenne', 'Mediane', 'Maximum', 'Mode', 'Init. Xenon', 'Opt. Xenon / moyenne', 'Opt. Xenon / mode'" );
    if ( minXenon > imHisto->dim.x-1 - par.maxXenon )
      fprintf( f, ", 2" );
    else 
      fprintf( f, ", 1" );
    fprintf( f, ");\n" );
    fprintf( f, "hold off;\n\n" );

    fprintf( f, "%% figure;\n" );
    fprintf( f, "%% hold on;\n" );
    fprintf( f, "%% ZZ = abs( (%f * YY) - (%f * XX) + XX .* YY );\n", ftnInitXenon.b, ftnInitXenon.a );
    fprintf( f, "%% xlabel('XENON');\n" );
    fprintf( f, "%% ylabel('HMPAO');\n" );
    fprintf( f, "%% surf(XX,YY,ZZ);\n" );
    fprintf( f, "%% shading interp\n" );
    fprintf( f, "%% grid;\n\n" );
    fprintf( f, "%% hold off;\n\n" );

    


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


      else if ( strcmp ( argv[i], "-mh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mh...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->maxHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -mh...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-mx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mx...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->maxXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -mx...\n", 0 );
      }





      else if ( strcmp ( argv[i], "-psh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -psh...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pourSeuilHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -psh...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-psx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -psx...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pourSeuilXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -psx...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-pmh" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -pmh...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pourMaxHmpao) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -pmh...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-pmx" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -pmx...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pourMaxXenon) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -pmx...\n", 0 );
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

  par->seuilHmpao = 0;
  par->seuilXenon = 0;
  par->maxHmpao = 65535;
  par->maxXenon = 65535;

  par->pourSeuilXenon = -1.0;
  par->pourSeuilHmpao = -1.0;
  par->pourMaxXenon = -1.0;
  par->pourMaxHmpao = -1.0;

  par->voxelHmpao = 0;
  par->voxelXenon = 0;
  par->matrice[0] ='\0';
}























static void _PrintFonctionCalib( typeFtnCalib *ftn, char *s )
{
  printf( "\n");
  if ( s != (char*)NULL ) printf( "%s\n", s );
  printf( "     Hmpao = %f * Xenon / ( %f + Xenon )\n", ftn->a, ftn->b );
  printf( "     Xenon = %f * Hmpao / ( %f - Hmpao )\n", ftn->b, ftn->a );
  printf( "             %f = (1 + alpha) Hr \n", ftn->a );
  printf( "             %f = alpha Xr \n", ftn->b );
  printf( "\n");
}




#ifdef _UNUSED_
static void _TraceFonctionCalib( vt_image *theIm,
				 typeFtnCalib *ftn, int z, float v )
{
  float ***theCurve = (float***)theIm->array;
  int hmpao, xenon;
  double val;
  
  for ( xenon=0; xenon < theIm->dim.x; xenon++ ) {
    val = ftn->a * xenon / ( ftn->b + xenon );
    if ( val < 0 ) continue;
    if ( val + 0.5 >= theIm->dim.y ) continue;
    theCurve[z][(int)(val+0.5)][xenon] = v;
  }
  for ( hmpao=0; hmpao < theIm->dim.y; hmpao++ ) {
    val = ftn->b * hmpao / ( ftn->a - hmpao );
    if ( val < 0 ) continue;
    if ( val + 0.5 >= theIm->dim.x ) continue;
    theCurve[z][hmpao][(int)(val+0.5)] = v;
  }
}
#endif





static int _InitFonctionCalib( typePoint *theXenon, typePoint *theHmpao,
			       int minx,
			       int maxx,
			       int miny,
			       int maxy,
			       typeFtnCalib *ftn,
			       int type )
{
  char *proc = "_InitFonctionCalib";
  int fx=minx;
  int lx=maxx;
  int fy=miny;
  int ly=maxy;
  int ix1, ix2;
  double x1, x2;
  int iy1, iy2;
  double y1, y2;
  int i;
  

  
  ftn->a = ftn->b = 0.0;

  if ( fx < 0 ) fx = 0;
  if ( lx - fx < 9 ) {
    fprintf( stderr, "%s: bornes X [%d,%d] (transformees en [%d,%d]) non valides\n",
	     proc, minx, maxx, fx, lx );
    return( 0 );
  }

  if ( fy < 0 ) fy = 0;
  if ( ly - fy < 9 ) {
    fprintf( stderr, "%s: bornes Y [%d,%d] (transformees en [%d,%d]) non valides\n",
	     proc, miny, maxy, fy, ly );
    return( 0 );
  }

  switch ( type ) {
  default :
  case _XENON_MED_ :
    x1 = fx + (double)(lx-fx) / (double)3.0;
    x2 = fx + 2.0 * (double)(lx-fx) / (double)3.0;
    
    ix1 = (int)(x1+0.5);
    ix2 = (int)(x2+0.5);

    y1 = (double)(theXenon[ix1-1].med + theXenon[ix1].med + theXenon[ix1+1].med) / 3.0;
    y2 = (double)(theXenon[ix2-1].med + theXenon[ix2].med + theXenon[ix2+1].med) / 3.0;

    ftn->a = y1*y2 * (x1-x2) / (x1*y2-x2*y1);
    ftn->b = x1*x2 * (y1-y2) / (x1*y2-x2*y1);
    break;

  case _XENON_MOY_ :
    x1 = fx + (double)(lx-fx) / (double)3.0;
    x2 = fx + 2.0 * (double)(lx-fx) / (double)3.0;
    
    ix1 = (int)(x1+0.5);
    ix2 = (int)(x2+0.5);


    if ( theXenon[ix1].sum0 <= _MINVAL_ ) {
      i = 0;
      while( i>= minx && i<=maxx 
	     && theXenon[ix1-i].sum0 <= _MINVAL_ 
	     && theXenon[ix1+i].sum0 <= _MINVAL_ ) {
	i++;
      }
      if (theXenon[ix1-i].sum0 > _MINVAL_) {
	x1 = ix1-i;
	y1 = theXenon[ix1-i].moy;
      } else {
	x1 = ix1+i;
	y1 = theXenon[ix1+i].moy;
      }
    } else {
      if ( theXenon[ix1-1].sum0 > _MINVAL_ && theXenon[ix1+1].sum0 > _MINVAL_ )
	y1 = (double)(theXenon[ix1-1].moy + theXenon[ix1].moy + theXenon[ix1+1].moy) / 3.0;
      else 
	y1 = (double)theXenon[ix1].moy;
    }
    
    
    if ( theXenon[ix2].sum0 <= _MINVAL_ ) {
      i = 0;
      while( i>= minx && i<=maxx 
	     && theXenon[ix2-i].sum0 <= _MINVAL_ 
	     && theXenon[ix2+i].sum0 <= _MINVAL_ ) {
	i++;
      }
      if (theXenon[ix2-i].sum0 > _MINVAL_) {
	x2 = ix2-i;
	y2 = theXenon[ix2-i].moy;
      } else {
	x2 = ix2+i;
	y2 = theXenon[ix2+i].moy;
      }
    } else {
      if ( theXenon[ix2-1].sum0 > _MINVAL_ && theXenon[ix2+1].sum0 > _MINVAL_ )
	y2 = (double)(theXenon[ix2-1].moy + theXenon[ix2].moy + theXenon[ix2+1].moy) / 3.0;
      else 
	y2 = (double)theXenon[ix2].moy;
    }


    ftn->a = y1*y2 * (x1-x2) / (x1*y2-x2*y1);
    ftn->b = x1*x2 * (y1-y2) / (x1*y2-x2*y1);
    break;

  case _HMPAO_MED_ :
    y1 = fy + (double)(ly-fy) / (double)3.0;
    y2 = fy + 2.0 * (double)(ly-fy) / (double)3.0;
    
    iy1 = (int)(y1+0.5);
    iy2 = (int)(y2+0.5);

    x1 = (double)(theHmpao[iy1-1].med + theHmpao[iy1].med + theHmpao[iy1+1].med) / 3.0;
    x2 = (double)(theHmpao[iy2-1].med + theHmpao[iy2].med + theHmpao[iy2+1].med) / 3.0;

    ftn->a = y1*y2 * (x1-x2) / (x1*y2-x2*y1);
    ftn->b = x1*x2 * (y1-y2) / (x1*y2-x2*y1);
    break;

  case _HMPAO_MOY_ :
    y1 = fy + (double)(ly-fy) / (double)3.0;
    y2 = fy + 2.0 * (double)(ly-fy) / (double)3.0;
    
    iy1 = (int)(y1+0.5);
    iy2 = (int)(y2+0.5);

    x1 = (double)(theHmpao[iy1-1].moy + theHmpao[iy1].moy + theHmpao[iy1+1].moy) / 3.0;
    x2 = (double)(theHmpao[iy2-1].moy + theHmpao[iy2].moy + theHmpao[iy2+1].moy) / 3.0;

    ftn->a = y1*y2 * (x1-x2) / (x1*y2-x2*y1);
    ftn->b = x1*x2 * (y1-y2) / (x1*y2-x2*y1);
    break;

  }
  return( 1 );
}
			



	     

void _InitValues( double *par,
		  double xenon1, double hmpao1,
		  double xenon2, double hmpao2 )
{
  double xen1 = xenon1;
  double hmp1 = hmpao1;
  double xen2 = xenon2;
  double hmp2 = hmpao2;

  if ( xen2 < xen1 ) {
    xen1 = xenon2;
    hmp1 = hmpao2;
    xen2 = xenon1;
    hmp2 = hmpao1;
  }
  while ( hmp1 + 5.0 >= hmp2 ) {
    hmp1 -= 2.0;
    hmp2 += 2.0;
  }
  fprintf( stdout, " points [X=%f H=%f] et [X=%f H=%f]\n",  xen1, hmp1, xen2 ,hmp2 );
  par[1] = (hmp2 - hmp1) * xen1 * xen2 / ( hmp1 * xen2 - hmp2 * xen1 );
  par[0] = (xen2 - xen1) * hmp1 * hmp2 / ( hmp1 * xen2 - hmp2 * xen1 );
}












double _ComputeLMValuesInLists( typePoint *theXenon, typePoint *theHmpao,
				int minx,
				int maxx,
				int miny,
				int maxy,
				typeFtnCalib *ftn,
				double *mat,
				double *vec,
				int type )
{
  double d[2];
  double val, m, c = 0.0;
  int xenon ,hmpao;

  double n;

  mat[0] = mat[1] = mat[2] = mat[3] = 0.0;
  vec[0] = vec[1] = 0.0;

  n = ftn->a * ftn->a + ftn->b * ftn->b;

  switch ( type ) {
  default :

  case _XENON_MOY_ :
    for ( xenon=minx; xenon<=maxx; xenon++ ) {
      if ( theXenon[xenon].sum0 <= _MINVAL_ ) continue;
      val = ftn->a * (double)xenon / ( ftn->b + (double)xenon );
      if ( val < miny || val > maxy ) continue;
      d[0] =   (double)xenon / ( ftn->b + (double)xenon );
      d[1] = - ftn->a * d[0] / ( ftn->b + (double)xenon );
      m = theXenon[xenon].sum0;
      if ( theXenon[xenon].variance > 1.0 ) m /= theXenon[xenon].variance;
      c      += m * (theXenon[xenon].moy-val) * (theXenon[xenon].moy-val);
      vec[0] += m * (theXenon[xenon].moy-val) * d[0];
      vec[1] += m * (theXenon[xenon].moy-val) * d[1];
      mat[0] += m * d[0] * d[0];
      mat[1] += m * d[0] * d[1];
      mat[2] += m * d[0] * d[1];
      mat[3] += m * d[1] * d[1];
    }

    break;

  case _HMPAO_MOY_ :
    for ( hmpao=miny; hmpao<=maxy; hmpao++ ) {
      if ( theHmpao[hmpao].sum0 <= _MINVAL_ ) continue;
      val = ftn->b * (double)hmpao / ( ftn->a - (double)hmpao );
      if ( val < minx || val > maxx ) continue;
      d[1] =   (double)hmpao / ( ftn->a - (double)hmpao );
      d[0] = - ftn->b * d[1] / ( ftn->a - (double)hmpao );
      m = theHmpao[hmpao].sum0;
      if ( theHmpao[hmpao].variance > 1.0 ) m /= theHmpao[hmpao].variance;
      c      += m * (theHmpao[hmpao].moy-val) * (theHmpao[hmpao].moy-val);
      vec[0] += m * (theHmpao[hmpao].moy-val) * d[0];
      vec[1] += m * (theHmpao[hmpao].moy-val) * d[1];
      mat[0] += m * d[0] * d[0];
      mat[1] += m * d[0] * d[1];
      mat[2] += m * d[0] * d[1];
      mat[3] += m * d[1] * d[1];
    }

    break;

    /*
  case _XENON_HMPAO_MOY_ :

    for ( xenon=minx; xenon<=maxx; xenon++ ) {
      if ( theXenon[xenon].sum0 <= _MINVAL_ ) continue;
      val = ftn->a * (double)xenon / ( ftn->b + (double)xenon );
      if ( val < miny || val > maxy ) continue;
      d[0] = (double)xenon / ( ftn->b + (double)xenon );
      d[1] = - ftn->a * d[0] / ( ftn->b + (double)xenon );
      c += theXenon[xenon].sum0 * (theXenon[xenon].med-val) * (theXenon[xenon].med-val);
      vec[0] += theXenon[xenon].sum0 * (theXenon[xenon].med-val) * d[0];
      vec[1] += theXenon[xenon].sum0 * (theXenon[xenon].med-val) * d[1];
      mat[0] += theXenon[xenon].sum0 * d[0] * d[0];
      mat[1] += theXenon[xenon].sum0 * d[0] * d[1];
      mat[2] += theXenon[xenon].sum0 * d[0] * d[1];
      mat[3] += theXenon[xenon].sum0 * d[1] * d[1];
    }

    for ( hmpao=miny; hmpao<=maxy; hmpao++ ) {
      if ( theHmpao[hmpao].sum0 <= _MINVAL_ ) continue;
      val = ftn->b * (double)hmpao / ( ftn->a - (double)hmpao );
      if ( val < minx || val > maxx ) continue;
      d[1] = (double)hmpao / ( ftn->a - (double)hmpao );
      d[0] = - ftn->b * d[1] / ( ftn->a - (double)hmpao );
      c += theHmpao[hmpao].sum0 * (theHmpao[hmpao].med-val) * (theHmpao[hmpao].med-val);
      vec[0] += theHmpao[hmpao].sum0 * (theHmpao[hmpao].med-val) * d[0];
      vec[1] += theHmpao[hmpao].sum0 * (theHmpao[hmpao].med-val) * d[1];
      mat[0] += theHmpao[hmpao].sum0 * d[0] * d[0];
      mat[1] += theHmpao[hmpao].sum0 * d[0] * d[1];
      mat[2] += theHmpao[hmpao].sum0 * d[0] * d[1];
      mat[3] += theHmpao[hmpao].sum0 * d[1] * d[1];
    }

    break;
    */
  }

  return( c );
}








double _ComputeLMValuesInImage( vt_image *theIm,
				typePoint *theXenon, typePoint *theHmpao,
				 int minx,
				int maxx,
				int miny,
				int maxy,
				typeFtnCalib *ftn,
				double *mat,
				double *vec,
				int type )
{
  char *proc = "_ComputeLMValuesInImage";
  double d[2];
  double val, m, c = 0.0;
  int hmpao, xenon;
  
  mat[0] = mat[1] = mat[2] = mat[3] = 0.0;
  vec[0] = vec[1] = 0.0;

  switch ( theIm->type ) {
  case FLOAT :
    {
      float ***theHist = (float***)theIm->array;

      switch ( type ) {
      default :
      case _IMAGE_XENON_MOY_ :

	for ( xenon=minx; xenon<=maxx; xenon++ ) {

	  if ( theXenon[xenon].sum0 <= _MINVAL_ ) continue;
	  val = ftn->a * (double)xenon / ( ftn->b + (double)xenon );
	  if ( val < miny || val > maxy ) continue;
	  d[0] =   (double)xenon / ( ftn->b + (double)xenon );
	  d[1] = - ftn->a * d[0] / ( ftn->b + (double)xenon );
	  
	  for ( hmpao=miny; hmpao<=maxy; hmpao++ ) {
	    m = theHist[0][hmpao][xenon];
	    if ( theXenon[xenon].variance > 1.0 ) m /= theXenon[xenon].variance;
	    c +=      m * (hmpao-val) * (hmpao-val);
	    vec[0] += m * (hmpao-val) * d[0];
	    vec[1] += m * (hmpao-val) * d[1];
	    mat[0] += m * d[0] * d[0];
	    mat[1] += m * d[0] * d[1];
	    mat[2] += m * d[0] * d[1];
	    mat[3] += m * d[1] * d[1];
	  }
	}
	break;

      }

    }
    break;
    
  default :
    fprintf( stderr, "%s: can not deal with such image type\n", proc );
    return( 0.0 );
  }
  return( c );
}







static double _INFINI_ = 1000;



static int _EstimeFonctionCalib( vt_image *theIm,
				 typePoint *theXenon, typePoint *theHmpao,
				 int minx,
				 int maxx,
				 int miny,
				 int maxy,
				 typeFtnCalib *ftnFinal,
				 int type )
{
  int _verbose_ = 0;

  typeFtnCalib fdf, df;
  double matAlpha[4];
  double vecBeta[2];

  double matPrime[4];
  double incBeta[2];

  double lambda = 0.001;
  double chi2, chi2dchi;
  int stop = 0;

  










  /* on a       hmpao =  a xenon / ( b + xenon )
   */

  



  switch ( type ) {
  default :
    chi2 = _ComputeLMValuesInLists( theXenon, theHmpao, minx, maxx, miny, maxy,
				    ftnFinal, matAlpha, vecBeta, type );
    break;
  case _IMAGE_XENON_MOY_:
    chi2 = _ComputeLMValuesInImage( theIm, theXenon, theHmpao, minx, maxx, miny, maxy,
				    ftnFinal, matAlpha, vecBeta, type );
  }
			
  
  do {
    matPrime[0] = matAlpha[0] + lambda;
    matPrime[1] = matAlpha[1];
    matPrime[2] = matAlpha[2];
    matPrime[3] = matAlpha[3] + lambda;
    
    df.a = (matPrime[3] * vecBeta[0] - matPrime[1] * vecBeta[1] ) /
      (matPrime[0] * matPrime[3] - matPrime[1] * matPrime[2] );
    df.b = (matPrime[0] * vecBeta[1] - matPrime[2] * vecBeta[0] ) /
      (matPrime[0] * matPrime[3] - matPrime[1] * matPrime[2] );
    
    fdf.a = ftnFinal->a + df.a;
    fdf.b = ftnFinal->b + df.b;

    if ( 0 ) {
      printf( " ...... [a] = { %f %f } - [da] = { %f  %f }\n", ftnFinal->a, ftnFinal->b, df.a, df.b );
      printf( " ...... [BETA] = { %f %f }\n", vecBeta[0], vecBeta[1] );
      printf( " ...... [ALPHA] = [ %f %f ]\n", matAlpha[0], matAlpha[1] );
      printf( "                  [ %f %f ]\n", matAlpha[2], matAlpha[3] );
      printf( " ...... [ALPHA'] = [ %f %f ]\n", matPrime[0], matPrime[1] );
      printf( "                   [ %f %f ]\n", matPrime[2], matPrime[3] );
    }
    
    switch ( type ) {
    default :
      chi2dchi = _ComputeLMValuesInLists( theXenon, theHmpao, minx, maxx, miny, maxy,
					  &fdf, matPrime, incBeta, type );
      break;
    case _IMAGE_XENON_MOY_ :
      chi2dchi = _ComputeLMValuesInImage( theIm, theXenon, theHmpao, minx, maxx, miny, maxy,
					  &fdf, matPrime, incBeta, type );
    }

    
    if ( chi2dchi >= chi2 ) {
      lambda *= 10;
      if ( _verbose_ ) {
	printf( " ... Chi2DChi( %f, %f ) = %f > Chi2( %f, %f ) = %f",
		 fdf.a, fdf.b, chi2dchi, ftnFinal->a, ftnFinal->b, chi2 );
	printf( "   =>   lambda *= 10 -> %f \n", lambda );
      }
      if ( lambda > 1e+10 ) stop = 1;
    } else {
      if ( (chi2-chi2dchi)/chi2 < 1e-4 ) stop = 1;
      
      /* on tend vers le lineaire
	 => on change ?
      */
      if ( 0 ) {
	if ( ftnFinal->a > _INFINI_ && fdf.a > ftnFinal->a ) {
	  fdf.a = -fdf.a;
	  fdf.b = -fdf.b;
	} 
	if ( ftnFinal->a < -_INFINI_ && fdf.a < ftnFinal->a ) {
	  fdf.a = -fdf.a;
	  fdf.b = -fdf.b;
	} 
      }
      
      lambda /= 10;
      if ( _verbose_ ) {
	printf( " ... Chi2DChi( %f, %f ) = %f < Chi2( %f, %f ) = %f",
		 fdf.a, fdf.b, chi2dchi, ftnFinal->a, ftnFinal->b, chi2 );
	printf( "   =>   lambda /= 10 -> %g \n", lambda );
      }

      *ftnFinal = fdf;
      chi2 = chi2dchi;
      
      matAlpha[0] = matPrime[0];
      matAlpha[1] = matPrime[1];
      matAlpha[2] = matPrime[2];      
      matAlpha[3] = matPrime[3];

      vecBeta[0] = incBeta[0];
      vecBeta[1] = incBeta[1];
    }

  } while ( stop == 0 );





  return( 1 );
}
					
















#ifdef _UNUSED_
static int _EstimeFonctionCalibGradient( vt_image *theIm,
					 typePoint *theXenon, typePoint *theHmpao,
					 int minx,
					 int maxx,
					 int miny,
					 int maxy,
					 typeFtnCalib *ftnFinal,
					 int type )
{
  int _verbose_ = 0;
  
  typeFtnCalib fdf;
  double matAlpha[4];
  double vecBeta[2];
  double matPrime[4];
  double lambda = 0.1;
  double chi2, chi2dchi;
  int stop = 0;

  



  
  switch ( type ) {
  default :
    chi2 = _ComputeLMValuesInLists( theXenon, theHmpao, minx, maxx, miny, maxy,
				    ftnFinal, matAlpha, vecBeta, type );
    break;
  }
			 
  do {

    fdf.a = ftnFinal->a - lambda * vecBeta[0];
    fdf.b = ftnFinal->b - lambda * vecBeta[1];
    
    switch ( type ) {
    default :
      chi2dchi = _ComputeLMValuesInLists( theXenon, theHmpao, minx, maxx, miny, maxy,
					  &fdf, matPrime, matPrime, type );
      break;
    }

    
    if ( chi2dchi >= chi2 ) {
      lambda /= 10;
      if ( _verbose_ ) {
	printf( " ... Chi2( a+da ) = %f > Chi2( a ) = %f  -  G=(%f,%f) \n",
		 chi2dchi, chi2, vecBeta[0], vecBeta[1] );
	printf( "     =>   lambda /= 10 -> %f \n", lambda );
      }
      if ( lambda < 1e-12 ) stop = 1;
    } else {
      *ftnFinal = fdf;
      if ( _verbose_ ) {
	printf( " ... Chi2( a+da ) = %f < Chi2( a ) = %f -  G=(%f,%f) \n",
		 chi2dchi, chi2, vecBeta[0], vecBeta[1]  );
      }
      chi2 = chi2dchi;
    }

  } while ( stop == 0 );





  return( 1 );
}
#endif					








int _UpdateLowThresholdXenon( vt_image *imHisto,
			      int old,
			      double p )
{
  double *t;
  int x, y;
  double s, sum;
  

  if ( p < 0.0 || p >= 100 ) return( old );

  t = (double*)malloc( imHisto->dim.x * sizeof(double) );
  for ( x=0; x <imHisto->dim.x; x++ ) t[x] = 0;

  switch ( imHisto->type ) {
  default :
    free( t );
    return( old );
    break;
  case FLOAT :
    {
      float ***theHisto = (float***)imHisto->array;
      for ( x=0; x<imHisto->dim.x; x++ ) {
	for ( y=0; y<imHisto->dim.y; y++ )
	  t[x] += theHisto[0][y][x];
      }
    }
    break;
  }

  sum = 0.0;
  for ( x=0; x<imHisto->dim.x; x++ ) sum += t[x];
  
  for ( x=0 , s=0.0; s<=p*sum/100.0 && x<imHisto->dim.x; x++, s+= t[x] )
    ;
  return( x );

}


int _UpdateHighThresholdXenon( vt_image *imHisto,
			      int old,
			      double p )
{
  double *t;
  int x, y;
  double s, sum;
  

  if ( p < 0.0 || p >= 100 ) return( old );

  t = (double*)malloc( imHisto->dim.x * sizeof(double) );
  for ( x=0; x <imHisto->dim.x; x++ ) t[x] = 0;

  switch ( imHisto->type ) {
  default :
    free( t );
    return( old );
    break;
  case FLOAT :
    {
      float ***theHisto = (float***)imHisto->array;
      for ( x=0; x<imHisto->dim.x; x++ ) {
	for ( y=0; y<imHisto->dim.y; y++ )
	  t[x] += theHisto[0][y][x];
      }
    }
    break;
  }

  sum = 0.0;
  for ( x=0; x<imHisto->dim.x; x++ ) sum += t[x];
  
  for ( x=imHisto->dim.x-1 , s=0.0; s<=p*sum/100.0 && x >= 0; x--, s+= t[x] )
    ;
  return( x );

}

















