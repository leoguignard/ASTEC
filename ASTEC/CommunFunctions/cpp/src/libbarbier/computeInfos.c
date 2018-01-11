/*************************************************************************
 * minimum.c -
 *
 * $Id: minimum.c,v 1.5 2000/08/16 16:31:56 greg Exp $
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

typedef struct local_par {
  char imcc[STRINGLENGTH];
  char imdist[STRINGLENGTH];
  char immask[STRINGLENGTH];
  char base[STRINGLENGTH];
  double scale;

} local_par;



typedef struct {
  int x;
  int y;
  int z;
  int d;
} typePoint;

typedef struct {
  int n;
  
  int ne;
  typePoint ext[2];

  typePoint *pt;
  
  int l;

  int dmin;
  int dmax;
  double dmean;

} typeCC;


/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );


static void _printNeighborhood( unsigned short int ***buf,
				int dimx, int dimy, int dimz,
				int x, int y, int z );

static void _extractCC( unsigned short int ***buf,
			int dimx, int dimy, int dimz,
			typeCC *theCC, int l );

static void _lengthOfCC( unsigned char ***buf,
			 int dimx, int dimy, int dimz,
			 typeCC *theCC, int l );

static void _diameterOfCC( unsigned short int ***buf,
			   int dimx, int dimy, int dimz,
			   typeCC *theCC, int l );

static void _printStatsCC( FILE *fout, typeCC *theCC, int nCC );



static char *usage = "[-cc %s] [-dist %s] [-base %s] [-scale %lf]\n\
\t [-mask %s]";

static char *detail = "\
\n\
 $Revision: 1.5 $ $Date: 2000/08/16 16:31:56 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imcc, *imdi, *imma  ;
  
  int x, y, z;
  int i, j, k;
  unsigned short int ***ccArray, ***diArray;
  unsigned char *** maArray;

  int n, nc, nCC;
  typeCC *theCC;

  int l, nv;

  int nmin, nmax, nmean, dmin, dmax;
  double dmean;


  char nsummary[STRINGLENGTH];
  char nhisto[STRINGLENGTH];
  FILE *fopen(), *fsumma, *fhisto;

  int *histo;
  int hn;
  double hmean;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  imcc = _VT_Inrimage( par.imcc );
  if ( imcc == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read components image\n", 0);
  if ( imcc->type != USHORT )
    VT_ErrorParse("components image not of USHORT type\n", 0);
  
  ccArray = (unsigned short int ***)imcc->array;

  imdi = _VT_Inrimage( par.imdist );
  if ( imdi == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read distance image\n", 0);
  if ( imdi->type != USHORT )
    VT_ErrorParse("components image not of USHORT type\n", 0);
  
  diArray = (unsigned short int ***)imdi->array;

  imma = _VT_Inrimage( par.immask );
  if ( imma == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read mask image\n", 0);
  if ( imma->type != UCHAR )
    VT_ErrorParse("mask image not of UCHAR type\n", 0);
  
  maArray = (unsigned char ***)imma->array;



  sprintf( nsummary, "%s.summary", par.base );
  sprintf( nhisto,   "%s.histo",   par.base );

  fsumma = fopen( nsummary, "w" );
  fhisto = fopen( nhisto, "w" );


  nCC = 0;
  for ( z=0; z<imcc->dim.z; z++ )
  for ( y=0; y<imcc->dim.y; y++ )
  for ( x=0; x<imcc->dim.x; x++ ) {
    if ( ccArray[z][y][x] == 0 ) continue;
    if ( nCC < ccArray[z][y][x] ) nCC = ccArray[z][y][x];
  }
  fprintf( stdout, "found %d components\n", nCC );





  theCC = (typeCC*)malloc( (nCC+1)*sizeof(typeCC) );
  if ( theCC == NULL ) {
    VT_ErrorParse("allocation error", 0);
  }





  for ( n=0; n<=nCC; n++ ) {

    theCC[n].n = 0;

    theCC[n].ne = 0;

    theCC[n].ext[0].x = -1;
    theCC[n].ext[0].y = -1;
    theCC[n].ext[0].z = -1;

    theCC[n].ext[1].x = -1;
    theCC[n].ext[1].y = -1;
    theCC[n].ext[1].z = -1;

    theCC[n].pt = NULL;
  }

  for ( z=0; z<imcc->dim.z; z++ )
  for ( y=0; y<imcc->dim.y; y++ )
  for ( x=0; x<imcc->dim.x; x++ ) {

    if ( ccArray[z][y][x] == 0 ) continue;

    l = ccArray[z][y][x];
    theCC[ l ].n ++;
    
    nv = 0;
    for ( k=-1; k<=1; k++ ) {
      if ( z+k < 0 || z+k >= imcc->dim.z ) continue;
      for ( j=-1; j<=1; j++ ) {
	if ( y+j < 0 || y+j >= imcc->dim.y ) continue;
	for ( i=-1; i<=1; i++ ) {
	  if ( x+i < 0 || x+i >= imcc->dim.x ) continue;
	  if ( i == 0 && j == 0 && k == 0 ) continue;
	  if ( ccArray[z+k][y+j][x+i] == l )
	    nv ++;
	}
      }
    }

    switch( nv ) {
	default :
      if ( 1 ) {
	fprintf( stdout, "%d neighbors for point (%d,%d,%d) of component %d ?\n", nv, x, y, z, l );
	_printNeighborhood( ccArray, imcc->dim.x, imcc->dim.y, imcc->dim.z, x, y, z );
      }
      break;
    case 2 :
      break;
    case 1 :
      if ( theCC[ l ].ne >= 2 )
	fprintf( stdout, "more than 2 extremities for component %d ?\n", l );
      else {
	theCC[ l ].ext[theCC[ l ].ne].x = x;
	theCC[ l ].ext[theCC[ l ].ne].y = y;
	theCC[ l ].ext[theCC[ l ].ne].z = z;
	theCC[ l ].ne ++;
      }
    }

  }

  _printStatsCC( stdout, theCC, nCC );


  for ( n=1; n<=nCC; n++ ) {

    if ( theCC[n].ne != 2 ) {
      fprintf( stdout, "component #%d has %d extremities ?\n", n, theCC[n].ne );
      continue;
    }

    theCC[n].pt = (typePoint*)malloc( theCC[n].n * sizeof(typePoint) );

    if ( theCC[n].pt == NULL ) 
      VT_ErrorParse("allocation error (2)", 0);

    _extractCC( ccArray, imcc->dim.x, imcc->dim.y, imcc->dim.z, &(theCC[n]), n );

    if ( 0 ) {
      for ( i=0; i<theCC[n].n; i++ ) {
	fprintf( stdout, "CC[%d].PT[%d/%d] = %d %d %d\n ", 
		 n, i, theCC[n].n, 
		 theCC[n].pt[i].x, theCC[n].pt[i].y, theCC[n].pt[i].z );
      }
    }

    _lengthOfCC( maArray, imma->dim.x, imma->dim.y, imma->dim.z, &(theCC[n]), n );

    _diameterOfCC( diArray, imdi->dim.x, imdi->dim.y, imdi->dim.z, &(theCC[n]), n );

  }

  
  for ( n=1; theCC[n].ne != 2; n ++ )
    ;
  fprintf( stdout, "first component is %d\n", n );

  nmin = nmax = theCC[n].n;
  nmean = 0;
  dmin = theCC[n].dmin;
  dmax = theCC[n].dmax;
  dmean = 0;

  for ( nc=0, n=1; n<=nCC; n++ ) {

    if ( theCC[n].ne != 2 ) continue;

    if ( nmin > theCC[n].n ) nmin = theCC[n].n;
    if ( nmax < theCC[n].n ) nmax = theCC[n].n;
    nmean += theCC[n].n;

    if ( dmin > theCC[n].dmin ) dmin = theCC[n].dmin;
    if ( dmax < theCC[n].dmax ) dmax = theCC[n].dmax;
    dmean += theCC[n].dmean * theCC[n].n;
    nc ++;
  }

  fprintf( fsumma, "\n" );
  fprintf( fsumma, "nb total   de pts = %d\n", nmean );
  fprintf( fsumma, "nb minimum de pts / cc = %d\n", nmin );
  fprintf( fsumma, "nb maximum de pts / cc = %d\n", nmax );
  fprintf( fsumma, "nb moyen   de pts / cc = %f\n", nmean/(double)nc );
  fprintf( fsumma, "Rayon min = %6.3f\n", dmin/(double)par.scale );
  fprintf( fsumma, "Rayon moy = %6.3f\n", dmean/(double)nmean/(double)par.scale );
  fprintf( fsumma, "Rayon max = %6.3f\n", dmax/(double)par.scale );
  fprintf( fsumma, "\n" );
  
  fprintf( fsumma, "# num. nomb.  long.   R.min  R.moy  R.max   E1.x    E1.y    E1.z     E2.x    E2.y    E2.z \n" );
  for ( n=1; n<=nCC; n++ ) {
    if ( theCC[n].ne != 2 ) continue;
    fprintf( fsumma, "%5d ", n );
    fprintf( fsumma, "%5d ", theCC[n].n );
    fprintf( fsumma, "%8.3f ", theCC[n].l /(double)par.scale );
    fprintf( fsumma, "%6.3f ", theCC[n].dmin/(double)par.scale );
    fprintf( fsumma, "%6.3f ", theCC[n].dmean/(double)par.scale );
    fprintf( fsumma, "%6.3f ", theCC[n].dmax/(double)par.scale );

    fprintf( fsumma, " %7.3f ", theCC[n].ext[0].x*imcc->siz.x );
    fprintf( fsumma, "%7.3f ", theCC[n].ext[0].y*imcc->siz.y );
    fprintf( fsumma, "%7.3f ", theCC[n].ext[0].z*imcc->siz.z );
    fprintf( fsumma, " %7.3f ", theCC[n].ext[1].x*imcc->siz.x );
    fprintf( fsumma, "%7.3f ", theCC[n].ext[1].y*imcc->siz.y );
    fprintf( fsumma, "%7.3f ", theCC[n].ext[1].z*imcc->siz.z );
    
    fprintf( fsumma, "\n" );
  }


  fclose( fsumma );

  
  histo = (int*)malloc( 65536*sizeof(int));
  for (i=0; i<65536; i++ ) histo[i]=0;

  for ( n=1; n<=nCC; n++ ) {
    if ( theCC[n].ne != 2 ) continue;
    for ( i=0; i<theCC[n].n; i++ ) 
      histo[ theCC[n].pt[i].d ] ++;
  }
  hn = 0;
  hmean = 0;
  for (i=0; i<65536; i++ ) {
    if ( histo[i] == 0 ) continue;
    hn += histo[i];
    hmean += i*histo[i];
  }
  fprintf( stdout, "histogramme: %d points, moyenne = %f\n", 
	   hn, hmean/(double)hn/(double)par.scale );  
  
  
  fprintf( fhisto, "#     N     R\n" ); 
  for (i=0; i<65536; i++ ) {
    if ( histo[i] == 0 ) continue;
    fprintf( fhisto, " %6d", histo[i] );
    fprintf( fhisto, " %7.3f", i/(double)par.scale );
    fprintf( fhisto, "\n" );
  }

  fclose( fhisto );


  /*--- liberations memoires ---*/
  VT_FreeImage( imdi );
  VT_Free( (void**)&imdi );

  VT_FreeImage( imcc );
  VT_Free( (void**)&imcc );
  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
{
  int i, nb, status;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }

      /*---  ---*/

      else if ( strcmp ( argv[i], "-cc" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cc...\n", 0 );
	strncpy( par->imcc, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-dist" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -dist...\n", 0 );
	strncpy( par->imdist, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-mask" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mask...\n", 0 );
	strncpy( par->immask, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-base" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -base...\n", 0 );
	strncpy( par->base, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-scale" ) == 0  ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -scale...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->scale) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -scale...\n", 0 );
      }

      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else {
      sprintf(text,"unknown option %s\n",argv[i]);
      VT_ErrorParse(text, 0);
    }
    i += 1;
  }
  
}






static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void VT_InitParam( local_par *par )
{
  par->imcc[0] = '\0';
  par->imdist[0] = '\0';
  par->immask[0] = '\0';
  par->base[0] = '\0';

  par->scale = 1.0;
}





static void _printNeighborhood( unsigned short int ***buf,
				int dimx, int dimy, int dimz,
				int x, int y, int z )
{
  int i, j, k;
  
  for ( j=-1; j<=1; j++ ) {
    for ( k=-1; k<=1; k++ ) {
      for ( i=-1; i<=1; i++ ) {
	if ( z+k < 0 || z+k >= dimz ||
	     y+j < 0 || y+j >= dimy ||
	     x+i < 0 || x+i >= dimx )
	  fprintf( stdout, "   .  " );
	else
	  fprintf( stdout, " %5d", buf[z+k][y+j][x+i] );
      }
      if ( k < 1 )
	fprintf( stdout, "  --  " );
    }
    fprintf( stdout, "\n" );
  }

}





static void _extractCC( unsigned short int ***buf,
			int dimx, int dimy, int dimz,
			typeCC *theCC, int l )
{
  int i, j, k;
  int px, py, pz;
  int n=0;
  int f;

  px = theCC->pt[n].x = theCC->ext[0].x;
  py = theCC->pt[n].y = theCC->ext[0].y;
  pz = theCC->pt[n].z = theCC->ext[0].z;
  
  /* on cherche n+1 a prtir de n 
   */

  do {

    f = 0;
    
    for ( k=-1; k<=1; k++ ) {
      if ( pz+k < 0 || pz+k >= dimz ) continue;
      for ( j=-1; j<=1; j++ ) {
	if ( py+j < 0 || py+j >= dimy ) continue;
	for ( i=-1; i<=1; i++ ) {
	  if ( px+i < 0 || px >= dimx ) continue;

	  if ( i == 0 && j == 0 && k == 0 ) continue;

	  /* on a un voisin 
	   */
	  if ( buf[pz+k][py+j][px+i] == l
	       && ( n == 0 || 
		    (theCC->pt[n-1].x != px+i || 
		     theCC->pt[n-1].y != py+j || 
		     theCC->pt[n-1].z != pz+k ) ) ) {
	    if ( f == 0 ) {
	      if ( n < theCC->n-1 ) {
		if ( 0 )
		  fprintf( stdout, "next point (%d) is %d %d %d\n", n+1, px+i, py+j, pz+k );
		theCC->pt[n+1].x = px+i;
		theCC->pt[n+1].y = py+j;
		theCC->pt[n+1].z = pz+k;
		f = 1;
	      }
	      else {
		fprintf( stderr, "too many points for component %d ?\n", l );
	      }
	    }
	    else {
	      fprintf( stderr, "next neighbor already found for component %d ?\n", l );
	    }
	  }
	  
	      
	}
      }
    }

    if ( f == 1 ) {
      px = theCC->pt[n+1].x;
      py = theCC->pt[n+1].y;
      pz = theCC->pt[n+1].z;
    }
    n ++;
    
  } while ( f == 1 );

}


static void _lengthOfCC( unsigned char ***buf,
			 int dimx, int dimy, int dimz,
			 typeCC *theCC, int l )
{
  int n;
  int cx = dimx/2;
  int cy = dimy/2;
  int cz = dimz/2;
  int i, j, k;
  
  if ( 0 )
    fprintf( stderr, "center is %d %d %d\n", cx, cy ,cz );

  theCC->l = 0;
  for ( n=1; n<theCC->n; n ++ ) {
    i = theCC->pt[n].x - theCC->pt[n-1].x;
    j = theCC->pt[n].y - theCC->pt[n-1].y;
    k = theCC->pt[n].z - theCC->pt[n-1].z;
    theCC->l += buf[cz+k][cy+j][cx+i];
  }

}





static void _diameterOfCC( unsigned short int ***buf,
			   int dimx, int dimy, int dimz,
			   typeCC *theCC, int l )
{
  int n;
  
  
  for ( n=0; n<theCC->n; n ++ ) {
    theCC->pt[n].d = buf[theCC->pt[n].z][theCC->pt[n].y][theCC->pt[n].x];
  }

  theCC->dmin = theCC->pt[0].d;
  theCC->dmax = theCC->pt[0].d;
  theCC->dmean = theCC->pt[0].d;
    
  for ( n=1; n<theCC->n; n ++ ) {
    if ( theCC->dmin > theCC->pt[n].d ) theCC->dmin = theCC->pt[n].d;
    if ( theCC->dmax < theCC->pt[n].d ) theCC->dmax = theCC->pt[n].d;
    theCC->dmean += theCC->pt[n].d;
  }

  theCC->dmean /= (double)theCC->n;

}





static void _printStatsCC( FILE *fout, typeCC *theCC, int nCC )
{
  int n;
  int nP, minP, maxP;
  
  nP = 0;
  minP = theCC[1].n;
  maxP = theCC[1].n;
  
  for ( n=1; n <=nCC; n++ ) {
    if ( theCC[n].ne != 2 ) continue;
    if ( minP > theCC[n].n ) minP = theCC[n].n;
    if ( maxP < theCC[n].n ) maxP = theCC[n].n;
    nP += theCC[n].n;
  }

  fprintf( fout, "\n" );
  fprintf( fout, "nb total   de pts = %d\n", nP );
  fprintf( fout, "nb minimum de pts / cc = %d\n", minP );
  fprintf( fout, "nb maximum de pts / cc = %d\n", maxP );
  fprintf( fout, "nb moyen   de pts / cc = %f\n", nP/(double)nCC );
  fprintf( fout, "\n" );
  

}
