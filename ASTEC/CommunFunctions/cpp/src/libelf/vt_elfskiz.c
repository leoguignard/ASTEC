/*************************************************************************
 * vt_elfskiz.c - extraction de bulles dans des images de mousse
 *
 * $Id: vt_elfskiz.c,v 1.6 2000/04/07 07:52:00 greg Exp $
 *
 * DESCRIPTION: 
 *
 * Outils d'extraction de parties connexes a l'aide d'un squelette par zone
 * d'influence et reconstruction
 *
 *
 * fonction VT_SkizFromBorder()
 * ============================
 * premiere version
 * a ne pas utiliser puisqu'elle utilise un tri, ce qui n'est pas
 * necessaire
 * 1. calcul de la distance avec VT_Dist()
 * 2. extraction des maxima regionaux avec VT_MaximaRegionaux()
 *    cette fonction n'est pas non plus la meilleure
 * 3. numerotation des maxima regionaux avec VT_ConnectedComponents()
 * 4. construction de la liste des points a reconstruire
 *    et tri selon la distance
 * 5. on reconstruit
 *
 *
 *
 * fonction VT_SkizFromBorderWithoutSort()
 * =======================================
 * seconde version
 * a utiliser puisqu'elle n'utilise pas de tri
 * 1. calcul de la distance avec VT_Dist()
 * 2. extraction des maxima regionaux avec VT_MaximaRegionauxWithList()
 * 3. numerotation des maxima regionaux avec VT_ConnectedComponents()
 * 4. construction de la liste des points a reconstruire
 *    double liste, avec un des indices egal a la distance
 * 5. on reconstruit
 * 
 *
 *
 * fonction VT_MaximaRegionauxWithList()
 * =====================================
 * on transforme l'image d'entree theBuf
 * en resBuf = min( theBuf*multiplier , theBuf-hauteur )
 *
 * on construit une liste de points correspondant a resBuf
 * puis on dilate resBuf en 26-connexite en restant "sous" theBuf
 * 
 * Ceci est realise pour chaque valeur de distance par valeur
 * decroissante, on change de valeur de distance des que les
 * points a cette valeur ne bougent plus.
 * 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sat May 12 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */


#include <vt_elfskiz.h>

#define _OUTSIDE_ 0
#define _INSIDE_ 1
#define _PROCESSED_ 2

typedef struct {
  unsigned short int x;
  unsigned short int y;
  unsigned short int z;
  unsigned char status;
  unsigned short int valeur;
} typeSkizPoint;


#ifndef NO_PROTO
static void _QSortSkizPointsList( typeSkizPoint* list, int left, int right );
static void _MySortSkizPointsList( typeSkizPoint* list, int maxVal, int left, int right );
#else
static void _QSortSkizpointsList();
static void _MySortSkizpointsList();
#endif


static int nbMaxDistances = 65536;







/* l'image d'entree est sur 1 octet
   on calcule la distance a partir des points
   non nuls,
   on numerote les maxima regionaux et 
   on 
   on dilate les maxima regionaux obtenus
 */
int VT_SkizFromBorder( vt_image *image, 
		       vt_image *labels, 
		       vt_image *dist,
		       vt_distance *par,
		       int height,
		       double multiplier ) 
{
  char *proc = "VT_SkizFromBorder";
  vt_connexe connexePar;
  unsigned short int ***theLabels = (unsigned short int ***)NULL;
  unsigned short int ***theDist = (unsigned short int ***)NULL;
  int x=0, y=0, z=0;
  int n;
  int localHeight = height;
  double localMultiplier = multiplier;

  vt_image auxIm;
  vt_image *curIm, *nextIm, *tmpIm;
  unsigned short int ***theNextLabels = (unsigned short int ***)NULL;
  typeSkizPoint *thePts = (typeSkizPoint *)NULL;
  typeSkizPoint tmp;
  int max;
  int i, d;
  int first;
  int changes;
  int iter;
  int j, k, l;

  if ( image->type != UCHAR ) {
    VT_Error( "input image should be of type unsigned char", proc );
    return( -1 );
  }
  if ( labels->type != USHORT ) {
    VT_Error( "labels image should be of type unsigned short", proc );
    return( -1 );
  }
  if ( dist->type != USHORT ) {
    VT_Error( "distance image should be of type unsigned short", proc );
    return( -1 );
  }
  
  if ( VT_Test2Image( image, labels, proc ) == -1 ) return( -1 );
  if ( VT_Test2Image( image, dist, proc ) == -1 ) return( -1 );

  
 
  /* calcul de la distance
   */
  par->seuil = 1;
  if ( VT_Dist( dist, image, par ) != 1 ) {
    VT_Error( "unable to compute distance", proc );
    return( -1 );
  }
  theDist = (unsigned short int ***)(dist->array);



  /* extraction des maxima regionaux
   */
  if ( localHeight < 1 ) localHeight = 1;
  if ( (localMultiplier <= 0.0) || (localMultiplier > 1.0) ) 
    localMultiplier = 1.0;
  if ( VT_MaximaRegionaux( dist, labels, localHeight, localMultiplier ) != 1 ) {
    VT_Error( "unable to compute regional maxima", proc );
    return( -1 );
  }
  

  /* numerotation des maxima regionaux
   */
  VT_Connexe( &connexePar );
  if (_VT_VERBOSE_) connexePar.verbose = 1;

  if ( VT_ConnectedComponents( labels, labels, (float)1.0, &connexePar ) != 1 ) {
    VT_Error( "unable to label regional maxima", proc );
    return( -1 );
  }
  theLabels = (unsigned short int ***)(labels->array);






  /* construction d'une liste de points
     points de distance > 0 et non maxima 
  */
  n = 0;
  max = 0;
  for ( z=0; z<image->dim.z; z++ )
  for ( y=0; y<image->dim.y; y++ ) 
  for ( x=0; x<image->dim.x; x++ ) {
    if ( (theDist[z][y][x] > 0) && (theLabels[z][y][x] == 0) ) {
      if ( max < theDist[z][y][x] ) max = theDist[z][y][x];
      n++;
    }
  }
  
  if ( n==0 ) {
    VT_Error( "found 0 points", proc );
    return( -1 );
  }

  
  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"%s: found %d points (max=%d)", proc, n, max );
  }


  thePts = (typeSkizPoint*)VT_Malloc( n * sizeof(typeSkizPoint) );
  if ( thePts == (typeSkizPoint*)NULL ) {
    VT_Error( "unable to allocate points list", proc );
    return( -1 );
  }


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr," ... construct" );
  }

  i = 0;
  for ( z=0; z<image->dim.z; z++ )
  for ( y=0; y<image->dim.y; y++ ) 
  for ( x=0; x<image->dim.x; x++ ) {
    if ( (theDist[z][y][x] > 0) && (theLabels[z][y][x] == 0) ) {
      thePts[i].x = x;
      thePts[i].y = y;
      thePts[i].z = z;
      thePts[i].valeur = theDist[z][y][x];
      thePts[i].status = _INSIDE_;
      if ( x == 0 ) { thePts[i].status = _OUTSIDE_; }
      else if ( x == image->dim.x - 1 ) { thePts[i].status = _OUTSIDE_; }
      else if ( y == 0 ) { thePts[i].status = _OUTSIDE_; }
      else if ( y == image->dim.y - 1 ) { thePts[i].status = _OUTSIDE_; }
      else if ( z == 0 ) { thePts[i].status = _OUTSIDE_; }
      else if ( z == image->dim.z - 1 ) { thePts[i].status = _OUTSIDE_; }
      i++;
    }
  }

  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"ed ... sort\n"  );
  }

  if ( 1 )
    _QSortSkizPointsList( thePts, 0, n-1 );
  if ( 0 ) 
    _MySortSkizPointsList( thePts, max, 0, n-1 );

 
  /*
  for(i=0; i<n;i++ ) {
    fprintf(stdout, "%3d : %3d %3d %3d\n", thePts[i].valeur,
	    thePts[i].x, thePts[i].y,thePts[i].z );
  }
  */


  VT_Image( &auxIm );
  VT_InitImage( &auxIm, proc, image->dim.x, image->dim.y, image->dim.z, USHORT );
  if ( VT_AllocImage( &auxIm ) != 1 ) {
    VT_Free( (void**)&thePts );
    VT_Error( "can not allocate auxiliary image", proc );
    return( -1 );
  }


  curIm = labels;
  nextIm = &auxIm;


  iter = 0;
  first = 0;
  do {

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"%s: iteration #%5d ", proc, iter+1 );
    }
    
    theLabels = (unsigned short int ***)curIm->array;
    theNextLabels = (unsigned short int ***)nextIm->array;
    (void)memcpy( nextIm->buf, curIm->buf, image->dim.x*image->dim.y*image->dim.z*sizeof(unsigned short int) );

    
    changes = 0;

    /* on traite les points a la meme distance
     */
    d = thePts[first].valeur;


    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"first = %8d, distance = %5d ", first, thePts[first].valeur );
    }

      
    /* la 6-dilatation fonctionne
       si les maxima sont cherches en 6-connexite
       malheureusement, on peut trouver des maxima sur des
       "lignes minces"
       donc on cherche les maxima en 26-connexite
       alors il faut dilater en 26-connexite
       malheureusement 
       ca fait des effets bizarre, il faut donc dilater 
       en 6-connexite, 
       si rien passer a la 18-connexite,
       si rien passer a la 26-connexite
    */
    

    /*--- recherche en 6 -connexite ---*/
    
    for ( i=first; i<n && thePts[i].valeur==d; i++ ) {
      
      x = thePts[i].x;
      y = thePts[i].y;
      z = thePts[i].z;
      
      
      if ( thePts[i].status == _INSIDE_ ) {
	
	if ( theLabels[z-1][y][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z-1][y][x];
	}
	else if ( theLabels[z+1][y][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z+1][y][x];
	}
	else if ( theLabels[z][y-1][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y-1][x];
	}
	else if ( theLabels[z][y+1][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y+1][x];
	}
	else if ( theLabels[z][y][x-1] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y][x-1];
	}
	else if ( theLabels[z][y][x+1] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y][x+1];
	} else {
	  continue;
	}
	
	changes ++;
	tmp = thePts[first]; 
	thePts[first] = thePts[i];
	thePts[i] = tmp;
	first ++;
	
      } else {
	
	if ( z > 0 && theLabels[z-1][y][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z-1][y][x];
	}
	else if ( z < image->dim.z-1  && theLabels[z+1][y][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z+1][y][x];
	}
	else if ( y > 0  && theLabels[z][y-1][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y-1][x];
	}
	else if ( y < image->dim.y-1 && theLabels[z][y+1][x] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y+1][x];
	}
	else if ( x > 0 && theLabels[z][y][x-1] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y][x-1];
	}
	else if ( x < image->dim.x-1 && theLabels[z][y][x+1] > 0 ) {
	  theNextLabels[z][y][x] = theLabels[z][y][x+1];
	} else {
	  continue;
	}
	
	changes ++;
	tmp = thePts[first]; 
	thePts[first] = thePts[i];
	thePts[i] = tmp;
	first ++;
	
      }
      
    } /* fin de la boucle for i ... */
    






    /*--- si rien, on passe en 18-connexite ---*/
    
    if ( changes == 0 ) {
      for ( i=first; i<n && thePts[i].valeur==d; i++ ) {
	/* il faudrait traiter le cas ou le point
	   est voisin de plusieurs zones ...
	*/
	x = thePts[i].x;
	y = thePts[i].y;
	z = thePts[i].z;
	
	
	if ( thePts[i].status == _INSIDE_ ) {
	  
	  for ( l= -1; (l<=1) && (theNextLabels[z][y][x]==0) ; l++ )
	    for ( k= -1; (k<=1) && (theNextLabels[z][y][x]==0) ; k++ )
	      for ( j= -1; (j<=1) && (theNextLabels[z][y][x]==0) ; j++ ) {
		if ( ((l==0) && (k*j == 0)) ||
		     ( ((l==1) || (l==-1)) && ((k*j != 0) || ((k==0)&&(j==0))) ) )
		  continue;
		if ( theLabels[z+l][y+k][x+j] > 0 ) theNextLabels[z][y][x] = theLabels[z+l][y+k][x+j];
	      }
	  
	  if ( theNextLabels[z][y][x] > 0 ) {
	    changes ++;
	    tmp = thePts[first]; 
	    thePts[first] = thePts[i];
	    thePts[i] = tmp;
	    first ++;
	  }
	  
	} else {
	  
	  for ( l= -1; (l<=1) && (theNextLabels[z][y][x]==0) ; l++ ) {
	    if ( (z+l >= 0) && (z+l < image->dim.z) ) {
	      for ( k= -1; (k<=1) && (theNextLabels[z][y][x]==0) ; k++ ) {
		if ( (y+k >= 0) && (y+k < image->dim.y) ) {
		  for ( j= -1; (j<=1) && (theNextLabels[z][y][x]==0) ; j++ ) {
		    if ( (x+j >= 0) && (x+j < image->dim.x) ) {
		      if ( ((l==0) && (k*j == 0)) ||
			   ( ((l==1) || (l==-1)) && ((k*j != 0) || ((k==0)&&(j==0))) ) )
			continue;
		      if ( theLabels[z+l][y+k][x+j] > 0 ) theNextLabels[z][y][x] = theLabels[z+l][y+k][x+j];
		    }
		  }
		}
	      }
	    }
	  }
	  
	  if ( theNextLabels[z][y][x] > 0 ) {
	    changes ++;
	    tmp = thePts[first]; 
	    thePts[first] = thePts[i];
	    thePts[i] = tmp;
	    first ++;
	  }
	  
	}
      }
    } /* fin de la 18-connexite */
    
    
    
    
    
    
    /*--- si rien, on passe en 26-connexite ---*/
    
    if ( changes == 0 ) {
      for ( i=first; i<n && thePts[i].valeur==d; i++ ) {
	/* il faudrait traiter le cas ou le point
	   est voisin de plusieurs zones ...
	*/
	x = thePts[i].x;
	y = thePts[i].y;
	z = thePts[i].z;
	
	
	if ( thePts[i].status == _INSIDE_ ) {
	  
	  for ( l= -1; (l<=1) && (theNextLabels[z][y][x]==0) ; l++ )
	    for ( k= -1; (k<=1) && (theNextLabels[z][y][x]==0) ; k++ )
	      for ( j= -1; (j<=1) && (theNextLabels[z][y][x]==0) ; j++ ) {
		if (l*k*j == 0) continue;
		if ( theLabels[z+l][y+k][x+j] > 0 ) theNextLabels[z][y][x] = theLabels[z+l][y+k][x+j];
	      }
	  
	  if ( theNextLabels[z][y][x] > 0 ) {
	    changes ++;
	    tmp = thePts[first]; 
	    thePts[first] = thePts[i];
	    thePts[i] = tmp;
	    first ++;
	  }
	  
	} else {
	  
	  for ( l= -1; (l<=1) && (theNextLabels[z][y][x]==0) ; l++ ) {
	    if ( (z+l >= 0) && (z+l < image->dim.z) ) {
	      for ( k= -1; (k<=1) && (theNextLabels[z][y][x]==0) ; k++ ) {
		if ( (y+k >= 0) && (y+k < image->dim.y) ) {
		  for ( j= -1; (j<=1) && (theNextLabels[z][y][x]==0) ; j++ ) {
		    if ( (x+j >= 0) && (x+j < image->dim.x) ) {
		      if (l*k*j == 0) continue;
		      if ( theLabels[z+l][y+k][x+j] > 0 ) theNextLabels[z][y][x] = theLabels[z+l][y+k][x+j];
		    }
		  }
		}
	      }
	    }
	  }
	  
	  if ( theNextLabels[z][y][x] > 0 ) {
	    changes ++;
	    tmp = thePts[first]; 
	    thePts[first] = thePts[i];
	    thePts[i] = tmp;
	    first ++;
	  }
	  
	}
      }
    } /* fin de la 26-connexite */




    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"changes = %8d\r", changes );
    }
    

    tmpIm = curIm;
    curIm = nextIm;
    nextIm = tmpIm;

    iter ++;

  } while( (changes > 0) && (first < n) );


  

  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"\n" );
  }

  VT_Free( (void**)&thePts );
  return( 1 );
}





























/* l'image d'entree est sur 1 octet
   on calcule la distance a partir des points
   non nuls,
   on numerote les maxima regionaux et 
   on 
   on dilate les maxima regionaux obtenus
 */
int VT_SkizFromBorderWithoutSort( vt_image *image, 
				  vt_image *labels, 
				  vt_image *dist,
				  enumElfDistance typeDistance,
				  int height,
				  double multiplier ) 
{
  char *proc = "VT_SkizFromBorderWithoutSort";
  vt_connexe connexePar;
  unsigned short int ***theLabels = (unsigned short int ***)NULL;
  unsigned short int ***theDist = (unsigned short int ***)NULL;
  int x=0, y=0, z=0;
  int n;
  int localHeight = height;
  double localMultiplier = multiplier;

  int *nbPtsByDistance = (int*)NULL;

  typeSkizPoint **thePts = (typeSkizPoint **)NULL;
  typeSkizPoint *tmpPts;
  typeSkizPoint tmp;
  int max;
  int i;
  int first, oldfirst;
  int changes;
  int iter;
  int j, k, l;

  if ( image->type != UCHAR ) {
    VT_Error( "input image should be of type unsigned char", proc );
    return( -1 );
  }
  if ( labels->type != USHORT ) {
    VT_Error( "labels image should be of type unsigned short", proc );
    return( -1 );
  }
  if ( dist->type != USHORT ) {
    VT_Error( "distance image should be of type unsigned short", proc );
    return( -1 );
  }
  
  if ( VT_Test2Image( image, labels, proc ) == -1 ) return( -1 );
  if ( VT_Test2Image( image, dist, proc ) == -1 ) return( -1 );
 




  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"%s: compute distance\n", proc );
  }


  
  /* calcul de la distance
   */
  switch ( typeDistance ) {
  default :
  case _EUCLIDIENNE_ :
    {
      vt_distance par;
      VT_Distance( &par );
      par.seuil = 1;
      par.type = VT_DIST_EUCLI;
      if ( VT_Dist( dist, image, &par ) != 1 ) {
	if ( _VT_VERBOSE_ ) {
	  fprintf( stderr,"%s: unable to compute distance\n", proc );
	  return( -1 );
	}
      }
    }
    break;
  case _CHAMFER_3_ :
    {
      int theDim[3];
      theDim[0] = image->dim.x;
      theDim[1] = image->dim.y;
      theDim[2] = image->dim.z;
      Compute3DChamfer3x3x3( image->buf, image->type,
			     dist->buf, dist->type,
			     theDim );
    }
    break;
  case _CHAMFER_5_ :
    {
      int theDim[3];
      theDim[0] = image->dim.x;
      theDim[1] = image->dim.y;
      theDim[2] = image->dim.z;
      Compute3DChamfer5x5x5( image->buf, image->type,
			     dist->buf, dist->type,
			     theDim );
    }
  }



  theDist = (unsigned short int ***)(dist->array);



  /* extraction des maxima regionaux
   */
  if ( localHeight < 1 ) localHeight = 1;
  if ( (localMultiplier <= 0.0) || (localMultiplier > 1.0) ) 
    localMultiplier = 1.0;

  if ( VT_MaximaRegionauxWithList( dist, labels, localHeight, localMultiplier ) != 1 ) {
    VT_Error( "unable to compute regional maxima", proc );
    return( -1 );
  }
  


  /* numerotation des maxima regionaux
   */
  VT_Connexe( &connexePar );
  if (_VT_VERBOSE_) connexePar.verbose = 1;

  if ( VT_ConnectedComponents( labels, labels, (float)1.0, &connexePar ) != 1 ) {
    VT_Error( "unable to label regional maxima", proc );
    return( -1 );
  }
  theLabels = (unsigned short int ***)(labels->array);



  nbPtsByDistance = (int*)malloc( (nbMaxDistances+1) * sizeof(int) );
  if ( nbPtsByDistance == (int*)NULL ) {
    VT_Error( "unable to allocate auxiliary array", proc );
    return( -1 );
  }
  for ( i=0;  i<= nbMaxDistances; i++ )
    nbPtsByDistance[ i ] = 0;





  /* construction d'une liste de points
     points de distance > 0 et non maxima 
  */
  n = 0;
  max = 0;
  for ( z=0; z<image->dim.z; z++ )
  for ( y=0; y<image->dim.y; y++ ) 
  for ( x=0; x<image->dim.x; x++ ) {
    if ( (theDist[z][y][x] > 0) && (theLabels[z][y][x] == 0) ) {
      nbPtsByDistance[ theDist[z][y][x] ] ++;
      if ( max < theDist[z][y][x] ) max = theDist[z][y][x];
      n++;
    }
  }
  


  if ( n==0 ) {
    VT_Error( "found 0 points", proc );
    return( -1 );
  }
  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"%s: found %d points (max=%d)", proc, n, max );
  }



  thePts = (typeSkizPoint**)malloc( (max+1)*sizeof(typeSkizPoint*) +
				    n*sizeof(typeSkizPoint) );
  if ( thePts == (typeSkizPoint**)NULL ) {
    free( nbPtsByDistance );
    VT_Error( "unable to allocate points list", proc );
    return( -1 );
  }
  
  tmpPts = (typeSkizPoint*)(thePts + (max+1));
  thePts[0] = (typeSkizPoint*)NULL;
  for ( i=1; i<=max; i++ ) {
    thePts[i] = tmpPts;
    tmpPts += nbPtsByDistance[i];
  }






  if ( _VT_VERBOSE_ ) {
    fprintf( stderr," ... construct" );
  }


  for ( i=0;  i<= nbMaxDistances; i++ )
    nbPtsByDistance[ i ] = 0;

  i = 0;
  for ( z=0; z<image->dim.z; z++ )
  for ( y=0; y<image->dim.y; y++ ) 
  for ( x=0; x<image->dim.x; x++ ) {
    if ( (theDist[z][y][x] > 0) && (theLabels[z][y][x] == 0) ) {
      tmpPts = &(thePts[ (int)theDist[z][y][x] ][ nbPtsByDistance[(int)theDist[z][y][x]] ]);
      
      tmpPts->x = x;
      tmpPts->y = y;
      tmpPts->z = z;
      tmpPts->valeur = theDist[z][y][x];
      tmpPts->status = _INSIDE_;
      if ( x == 0 ) { tmpPts->status = _OUTSIDE_; }
      else if ( x == image->dim.x - 1 ) { tmpPts->status = _OUTSIDE_; }
      else if ( y == 0 ) { tmpPts->status = _OUTSIDE_; }
      else if ( y == image->dim.y - 1 ) { tmpPts->status = _OUTSIDE_; }
      else if ( z == 0 ) { tmpPts->status = _OUTSIDE_; }
      else if ( z == image->dim.z - 1 ) { tmpPts->status = _OUTSIDE_; }

      nbPtsByDistance[(int)theDist[z][y][x]] ++;
    }
  }

  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"ed\n"  );
  }
 








  /* mettre le nouveau label dans 
     le points et remettre a jour
     l'image a la fin de la boucle
     a partir des points traites
     
     en supprimant la distance et 
     en ajoutant le label (pas de changement)

     on economise l'image intermediaire

     
  */



  theLabels = (unsigned short int ***)labels->array;

  
  iter = 0;
  first = 0;
  tmpPts = thePts[max];

  do {

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"%s: iteration #%5d ", proc, iter+1 );
      fprintf( stderr,"(%8d ->%8d), distance = %5d ", 
	       first, nbPtsByDistance[ max ]-1, max );
    }
    

    
    changes = 0;

    
    /* la 6-dilatation fonctionne
       si les maxima sont cherches en 6-connexite
       malheureusement, on peut trouver des maxima sur des
       "lignes minces"
       donc on cherche les maxima en 26-connexite
       alors il faut dilater en 26-connexite
       malheureusement 
       ca fait des effets bizarre, il faut donc dilater 
       en 6-connexite, 
       si rien passer a la 18-connexite,
       si rien passer a la 26-connexite
    */
    
    
    /*--- recherche en 6 -connexite ---*/
    oldfirst = first;
    for ( i=first; i<nbPtsByDistance[ max ]; i++ ) {
      
      x = tmpPts[i].x;
      y = tmpPts[i].y;
      z = tmpPts[i].z;
      
      
      switch (  tmpPts[i].status ) {

      default : 
	break;

      case _INSIDE_ :
	
	if ( theLabels[z-1][y][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z-1][y][x];
	}
	else if ( theLabels[z+1][y][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z+1][y][x];
	}
	else if ( theLabels[z][y-1][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y-1][x];
	}
	else if ( theLabels[z][y+1][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y+1][x];
	}
	else if ( theLabels[z][y][x-1] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y][x-1];
	}
	else if ( theLabels[z][y][x+1] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y][x+1];
	} else {
	  continue;
	}
	
	changes ++;
	tmp = tmpPts[first]; 
	tmpPts[first] = tmpPts[i];
	tmpPts[i] = tmp;
	first ++;
	
	break;

      case _OUTSIDE_ :
	
	if ( z > 0 && theLabels[z-1][y][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z-1][y][x];
	}
	else if ( z < image->dim.z-1  && theLabels[z+1][y][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z+1][y][x];
	}
	else if ( y > 0  && theLabels[z][y-1][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y-1][x];
	}
	else if ( y < image->dim.y-1 && theLabels[z][y+1][x] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y+1][x];
	}
	else if ( x > 0 && theLabels[z][y][x-1] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y][x-1];
	}
	else if ( x < image->dim.x-1 && theLabels[z][y][x+1] > 0 ) {
	  tmpPts[i].valeur = theLabels[z][y][x+1];
	} else {
	  continue;
	}
	
	changes ++;
	tmp = tmpPts[first]; 
	tmpPts[first] = tmpPts[i];
	tmpPts[i] = tmp;
	first ++;
	
      }
      
    } /* fin de la boucle for i ... */
    
    
    
    
    
    
    
    /*--- si rien, on passe en 18-connexite ---*/
    
    if ( changes == 0 ) {
      for ( i=first; i<nbPtsByDistance[ max ]; i++ ) {
	/* il faudrait traiter le cas ou le point
	   est voisin de plusieurs zones ...
	*/
	x = tmpPts[i].x;
	y = tmpPts[i].y;
	z = tmpPts[i].z;
	
	
	switch (  tmpPts[i].status ) {

	default : 
	  break;

	case _INSIDE_ :
	  
	  for ( l= -1; (l<=1) && (tmpPts[i].status!=_PROCESSED_) ; l++ )
	  for ( k= -1; (k<=1) && (tmpPts[i].status!=_PROCESSED_) ; k++ )
	  for ( j= -1; (j<=1) && (tmpPts[i].status!=_PROCESSED_) ; j++ ) {
	    if ( ((l==0) && (k*j == 0)) ||
		 ( ((l==1) || (l==-1)) && ((k*j != 0) || ((k==0)&&(j==0))) ) )
	      continue;
	    if ( theLabels[z+l][y+k][x+j] > 0 ) {
	      tmpPts[i].valeur = theLabels[z+l][y+k][x+j];
	      tmpPts[i].status = _PROCESSED_;
	    }
	  }
	  if ( tmpPts[i].status == _PROCESSED_ ) {
	    changes ++;
	    tmp = tmpPts[first]; 
	    tmpPts[first] = tmpPts[i];
	    tmpPts[i] = tmp;
	    first ++;
	  }
	  
	  break;
	  
	case _OUTSIDE_ :
	  
	  for ( l= -1; (l<=1) && (tmpPts[i].status!=_PROCESSED_) ; l++ ) {
	    if ( (z+l >= 0) && (z+l < image->dim.z) ) {
	      for ( k= -1; (k<=1) && (tmpPts[i].status!=_PROCESSED_) ; k++ ) {
		if ( (y+k >= 0) && (y+k < image->dim.y) ) {
		  for ( j= -1; (j<=1) && (tmpPts[i].status!=_PROCESSED_) ; j++ ) {
		    if ( (x+j >= 0) && (x+j < image->dim.x) ) {
		      if ( ((l==0) && (k*j == 0)) ||
			   ( ((l==1) || (l==-1)) && ((k*j != 0) || ((k==0)&&(j==0))) ) )
			continue;
		      if ( theLabels[z+l][y+k][x+j] > 0 ) {
			tmpPts[i].valeur = theLabels[z+l][y+k][x+j];
			tmpPts[i].status = _PROCESSED_;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  if ( tmpPts[i].status == _PROCESSED_ ) {
	    changes ++;
	    tmp = tmpPts[first]; 
	    tmpPts[first] = tmpPts[i];
	    tmpPts[i] = tmp;
	    first ++;
	  }
	  
	}

      } /* boucle for i ... */
    } /* fin de la 18-connexite */
    
    
    
    
    
    
    /*--- si rien, on passe en 26-connexite ---*/
    
    if ( changes == 0 ) {
      for ( i=first; i<nbPtsByDistance[ max ]; i++ ) {
	/* il faudrait traiter le cas ou le point
	   est voisin de plusieurs zones ...
	*/
	x = tmpPts[i].x;
	y = tmpPts[i].y;
	z = tmpPts[i].z;
	
	
	switch (  tmpPts[i].status ) {

	default : 
	  break;
	  
	case _INSIDE_ :
	  
	  for ( l= -1; (l<=1) && (tmpPts[i].status!=_PROCESSED_) ; l++ )
	  for ( k= -1; (k<=1) && (tmpPts[i].status!=_PROCESSED_) ; k++ )
	  for ( j= -1; (j<=1) && (tmpPts[i].status!=_PROCESSED_) ; j++ ) {
	    if (l*k*j == 0) continue;
	    if ( theLabels[z+l][y+k][x+j] > 0 ) {
	      tmpPts[i].valeur = theLabels[z+l][y+k][x+j];
	      tmpPts[i].status = _PROCESSED_;
	    }
	  }

	  if ( tmpPts[i].status == _PROCESSED_ ) {
	    changes ++;
	    tmp = tmpPts[first]; 
	    tmpPts[first] = tmpPts[i];
	    tmpPts[i] = tmp;
	    first ++;
	  }
	  
	  break;

	case _OUTSIDE_ :
	  
	  for ( l= -1; (l<=1) && (tmpPts[i].status!=_PROCESSED_) ; l++ ) {
	    if ( (z+l >= 0) && (z+l < image->dim.z) ) {
	      for ( k= -1; (k<=1) && (tmpPts[i].status!=_PROCESSED_) ; k++ ) {
		if ( (y+k >= 0) && (y+k < image->dim.y) ) {
		  for ( j= -1; (j<=1) && (tmpPts[i].status!=_PROCESSED_) ; j++ ) {
		    if ( (x+j >= 0) && (x+j < image->dim.x) ) {
		      if (l*k*j == 0) continue;
		      if ( theLabels[z+l][y+k][x+j] > 0 ) {
			tmpPts[i].valeur = theLabels[z+l][y+k][x+j];
			tmpPts[i].status = _PROCESSED_;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	  
	  if ( tmpPts[i].status == _PROCESSED_ ) {
	    changes ++;
	    tmp = tmpPts[first]; 
	    tmpPts[first] = tmpPts[i];
	    tmpPts[i] = tmp;
	    first ++;
	  }
	  
	}

      } /* boucle for i ... */
    } /* fin de la 26-connexite */
    

    /* on met les valeurs dans l'image */
    for ( i=oldfirst; i<first; i++ ) {
      theLabels[ tmpPts[i].z ][ tmpPts[i].y ][ tmpPts[i].x ] = tmpPts[i].valeur;
    }


    if ( first >= nbPtsByDistance[ max ] ) {
      /* on passe a la distance inferieure */
      max --;
      tmpPts = thePts[max];
      first = 0;
    }
    

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"changes = %8d\r", changes );
    }
    
    iter ++;

  } while( max >= 1 );


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"\n" );
  }


  free( nbPtsByDistance );
  VT_Free( (void**)&thePts );
  return( 1 );
}



















/* faire un maxima ressemblant
   a la dilatation 
   conditionnelle
   avec une liste de points

   ne pas traiter les points a zero de l'image de depart

   si ca atteint le niveau de gris de l'image 
   de depart, on peut rejeter

*/



/* les maxima regionaux se  definissent 
   comme la difference entre 
   l'image originale et le resultat de la dilatation a l'infini de f-h par rapport a f
   g(n+1) = Inf( g(n) + B , f )
   avec g(0) = f-h
   (on dilate f-h sous f)
   maxima regionaux = f - g(infini)

   on peut aussi les voir comme 
   (f - g(1)) - reconstruction( f - g(1) par g(2) - g(1) )

   
   ici, on construit une image de la facon suivante :
   resBuf = MIN ( theBuf - hauteur , theBuf * multiplier )
   que l'on dilate "sous la courbe", 
   on recupere ensuite les extrema par soustraction.


*/ 

int VT_MaximaRegionauxWithList( vt_image *theIm /* input image = distance */, 
			vt_image *resIm /* output image = maxima */,
			int hauteur,
				double multiplier )
{
  char *proc="VT_MaximaRegionauxWithList";
  unsigned short int ***theBuf;
  unsigned short int ***resBuf;
  unsigned short val;
  int x, y, z, n;
  int iter;
  int i, j, k, l;

  int *nbPtsByDistance = (int*)NULL;

  typeSkizPoint **thePts = (typeSkizPoint **)NULL;
  typeSkizPoint *tmpPts;
  typeSkizPoint tmp;
  int max, min;
  int first, oldfirst;
  int changes;



  if ( VT_Test2Image( theIm, resIm, proc ) == -1 ) return( -1 );
  if ( theIm->type != USHORT ) {
    VT_Error( "input image type should be unsigned short", proc );
    return( -1 );
  }
  if ( resIm->type != USHORT ) {
    VT_Error( "result image type should be unsigned short", proc );
    return( -1 );
  }
  if ( theIm->array == resIm->array ) {
    VT_Error( "input and result images should be different", proc );
    return( -1 );
  }



  /* seuillage */
  theBuf = (unsigned short int ***)(theIm->array);
  resBuf = (unsigned short int ***)(resIm->array);
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {

    min = (int)theBuf[z][y][x] - hauteur;

    if ( min < 0 ) min = 0;
    else if (  min > (int)( (double)theBuf[z][y][x] * multiplier + 0.5 ) )
      min = (int)( (double)theBuf[z][y][x] * multiplier + 0.5 );

    resBuf[z][y][x] = (unsigned short int)min;
  }






  nbPtsByDistance = (int*)malloc( (nbMaxDistances+1) * sizeof(int) );
  if ( nbPtsByDistance == (int*)NULL ) {
    VT_Error( "unable to allocate auxiliary array", proc );
    return( -1 );
  }
  for ( i=0;  i<= nbMaxDistances; i++ )
    nbPtsByDistance[ i ] = 0;



  /* construction d'une liste de points
     points de distance > 0 dans l'image originale 
  */
  n = 0;
  max = 0;
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {
    if ( theBuf[z][y][x] == 0 ) continue;
    nbPtsByDistance[ (int)resBuf[z][y][x] ] ++;
    if ( max < resBuf[z][y][x] ) max = resBuf[z][y][x];
    n ++;
  }
  
  if ( (n==0) || (nbPtsByDistance[max] == n) ) {
    VT_Error( "found 0 points", proc );
    return( -1 );
  }
  
  
  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"%s: found %d points (max=%d, #points=%d)", proc, n, max, nbPtsByDistance[max]);
  }


  /* pas la peine de traiter les points a max 
   */
  n = 0;
  for (i=0; i<max; i++ ) n += nbPtsByDistance[i];



  thePts = (typeSkizPoint**)malloc( (max)*sizeof(typeSkizPoint*) +
				    n*sizeof(typeSkizPoint) );
  if ( thePts == (typeSkizPoint**)NULL ) {
    free( nbPtsByDistance );
    VT_Error( "unable to allocate points list", proc );
    return( -1 );
  }
  
  tmpPts = (typeSkizPoint*)(thePts + max);
  thePts[0] = (typeSkizPoint*)NULL;
  for ( i=0; i<max; i++ ) {
    thePts[i] = tmpPts;
    tmpPts += nbPtsByDistance[i];
  }


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr," ... construct" );
  }


  for ( i=0;  i<= nbMaxDistances; i++ )
    nbPtsByDistance[ i ] = 0;

  i = 0;
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {
    if ( theBuf[z][y][x] == 0 ) continue;
    if ( resBuf[z][y][x] == max ) continue;
    tmpPts = &(thePts[ (int)resBuf[z][y][x] ][ nbPtsByDistance[(int)resBuf[z][y][x]] ]);
    
    tmpPts->x = x;
    tmpPts->y = y;
    tmpPts->z = z;
    tmpPts->valeur = resBuf[z][y][x];
    tmpPts->status = _INSIDE_;
    if ( x == 0 ) { tmpPts->status = _OUTSIDE_; }
    else if ( x == theIm->dim.x - 1 ) { tmpPts->status = _OUTSIDE_; }
    else if ( y == 0 ) { tmpPts->status = _OUTSIDE_; }
    else if ( y == theIm->dim.y - 1 ) { tmpPts->status = _OUTSIDE_; }
    else if ( z == 0 ) { tmpPts->status = _OUTSIDE_; }
    else if ( z == theIm->dim.z - 1 ) { tmpPts->status = _OUTSIDE_; }

    nbPtsByDistance[ (int)resBuf[z][y][x] ] ++;
  }


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"ed\n"  );
  }
 

  iter = 0;
  first = 0;
  max --;
  tmpPts = thePts[max];

  do {

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"%s: iteration #%5d ", proc, iter+1 );
      fprintf( stderr,"(%8d ->%8d), distance = %5d ", 
	       first, nbPtsByDistance[ max ]-1, max );
    }
    
    changes = 0;



    
    /*--- recherche en 26 -connexite ---*/
    oldfirst = first;
    for ( i=first; i<nbPtsByDistance[ max ]; i++ ) {
      x = tmpPts[i].x;
      y = tmpPts[i].y;
      z = tmpPts[i].z;

      val = resBuf[z][y][x];
 
      switch (  tmpPts[i].status ) {

      default : 
	break;
	  
      case _INSIDE_ :
	  
	for ( l= -1; l<=1; l++ )
	for ( k= -1; k<=1; k++ )
	for ( j= -1; j<=1; j++ ) {
	  if ( val < resBuf[z+l][y+k][x+j] ) val = resBuf[z+l][y+k][x+j];
	}
	
	break;
	
      case _OUTSIDE_ :

	for ( l= -1; l<=1; l++ ) {
	  if ( (z+l >= 0) && (z+l < theIm->dim.z) ) {
	    for ( k= -1; k<=1; k++ ) {
	      if ( (y+k >= 0) && (y+k < theIm->dim.y) ) {
		for ( j= -1; j<=1; j++ ) {
		  if ( (x+j >= 0) && (x+j < theIm->dim.x) ) {
		    if ( val < resBuf[z+l][y+k][x+j] ) val = resBuf[z+l][y+k][x+j];
		  }
		}
	      }
	    }
	  }
	}
      } /* fin de switch */
	  
      if ( theBuf[z][y][x] <= val ) {
	/* on ne touchera plus le point */
	tmpPts[i].valeur =  theBuf[z][y][x];
	changes ++;
	tmp = tmpPts[first]; 
	tmpPts[first] = tmpPts[i];
	tmpPts[i] = tmp;
	first ++;
      } else {
	if ( val > resBuf[z][y][x] ) {
	  tmpPts[i].valeur = val;
	  changes++;
	}
      }

    } /* boucle for i ... */

    for ( i=oldfirst; i<nbPtsByDistance[ max ]; i++ ) {
      resBuf[ tmpPts[i].z ][ tmpPts[i].y ][ tmpPts[i].x ] = tmpPts[i].valeur;
    }

    if ( (first >= nbPtsByDistance[ max ]) || (changes==0) ) {
      /* on passe a la distance inferieure */
      max --;
      tmpPts = thePts[max];
      first = 0;
    }

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"changes = %8d\r", changes );
    }
    
    iter ++;


  } while( max >= 0 );


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"\n" );
  }


  /* on ne recupere que les maxima */
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {
    if ( theBuf[z][y][x] <= resBuf[z][y][x] ) resBuf[z][y][x] = 0;
    else resBuf[z][y][x] = (unsigned short int)((int)theBuf[z][y][x] - (int)resBuf[z][y][x]);
  }

  free( nbPtsByDistance );
  VT_Free( (void**)&thePts );
  return( 1 );
}


















/* les maxima regionaux se  definissent 
   comme la difference entre 
   l'image originale et le resultat de la dilatation a l'infini de f-h par rapport a f
   g(n+1) = Inf( g(n) + B , f )
   avec g(0) = f-h
   (on dilate f-h sous f)
   maxima regionaux = f - g(infini)

   on peut aussi les voir comme 
   (f - g(1)) - reconstruction( f - g(1) par g(2) - g(1) )

*/ 

int VT_MaximaRegionaux( vt_image *theIm /* input image = distance */, 
			vt_image *resIm /* output image = maxima */,
			int hauteur,
				double multiplier )
{
  char *proc="VT_MaximaRegionaux";
  unsigned short int ***theBuf;
  unsigned short int ***resBuf;
  unsigned short val;
  int x, y, z, n;
  int iter;
  int insideZ, insideY, inside;
  int i, j, k;
  int min;

  if ( VT_Test2Image( theIm, resIm, proc ) == -1 ) return( -1 );
  if ( theIm->type != USHORT ) {
    VT_Error( "input image type should be unsigned short", proc );
    return( -1 );
  }
  if ( resIm->type != USHORT ) {
    VT_Error( "result image type should be unsigned short", proc );
    return( -1 );
  }
  if ( theIm->array == resIm->array ) {
    VT_Error( "input and result images should be different", proc );
    return( -1 );
  }


  /* seuillage */
  theBuf = (unsigned short int ***)(theIm->array);
  resBuf = (unsigned short int ***)(resIm->array);
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {

    min = (int)theBuf[z][y][x] - hauteur;

    if ( min < 0 ) min = 0;
    else if (  min > (int)( (double)theBuf[z][y][x] * multiplier + 0.5 ) )
      min = (int)( (double)theBuf[z][y][x] * multiplier + 0.5 );

    resBuf[z][y][x] = (unsigned short int)min;
  }

  iter=0;
  do { 
    n = 0;

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"%s: iteration #%5d ", proc, iter+1 );
    }
    
    insideZ = insideY = inside = 0;

    for ( z=0; z<theIm->dim.z; z++ ) {
      if ( (z==0) || (z==theIm->dim.z-1) ) insideZ = 0;
      else insideZ = 1;

      for ( y=0; y<theIm->dim.y; y++ ) {
	if ( (insideZ == 0) || ((insideZ == 1) && ((y==0) || (y==theIm->dim.y-1))) ) insideY = 0;
	else insideY = 1;

	for ( x=0; x<theIm->dim.x; x++ ) {
	  if ( (insideY == 0) || ((insideY == 1) && ((x==0) || (x==theIm->dim.x-1))) ) inside = 0;
	  else inside = 1;
	  
	  val = resBuf[z][y][x];
	  
	  if ( val == theBuf[z][y][x] ) {
	    continue;
	  }


	  /* apres essai avec la 6-dilatation, 
	     il faut faire une 26-dilatation 
	  */

	  /* saisie du voisinage */
	  if ( inside == 1 ) {
	    for ( k= -1; k<=1; k++ )
	    for ( j= -1; j<=1; j++ )
	    for ( i= -1; i<=1; i++ ) 
	      if ( val < resBuf[z+k][y+j][x+i] ) val = resBuf[z+k][y+j][x+i];
	  } else {
	    
	    for ( k= -1; k<=1; k++ ) {
	      if ( (z+k >= 0) && (z+k < theIm->dim.z) ) {
		for ( j= -1; j<=1; j++ ) {
		  if ( (y+j >= 0) && (y+j < theIm->dim.y) ) {
		    for ( i= -1; i<=1; i++ ) {
		      if ( (x+i >= 0) && (x+i < theIm->dim.x) ) {
			if ( val < resBuf[z+k][y+j][x+i] ) val = resBuf[z+k][y+j][x+i];
		      }
		    } /* fin de boucle for i */
		  }
		} /* fin de boucle for j */
	      }
	    } /* fin de boucle for k */

	  } /* fin de inside == 0 */
	
	  
	  /*----
	    6-dilatation
	  val = resBuf[z][y][x];
	  if ( z > 0 ) { if ( val < resBuf[z-1][y][x] ) val = resBuf[z-1][y][x]; }
	  if ( z < theIm->dim.z-1 ) { if ( val < resBuf[z+1][y][x] ) val = resBuf[z+1][y][x]; }
	  if ( y > 0 ) { if ( val < resBuf[z][y-1][x] ) val = resBuf[z][y-1][x]; }
	  if ( y < theIm->dim.y-1 ) { if ( val < resBuf[z][y+1][x] ) val = resBuf[z][y+1][x]; }
	  if ( x > 0 ) { if ( val < resBuf[z][y][x-1] ) val = resBuf[z][y][x-1]; }
	  if ( x < theIm->dim.x-1 ) { if ( val < resBuf[z][y][x+1] ) val = resBuf[z][y][x+1]; }
	  ---*/
	  
	  if ( val > theBuf[z][y][x] ) val = theBuf[z][y][x];
	  if ( val > resBuf[z][y][x] ) n++;
	  resBuf[z][y][x] = val;

	}
      }
    }

    if ( _VT_VERBOSE_ ) {
      fprintf( stderr,"changes = %8d        \n", n );
    }

    iter++;

  } while( n > 0 ); 


  if ( _VT_VERBOSE_ ) {
    fprintf( stderr,"\n" );
  }


  /* le resultat est dans currentImage */
  
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ ) 
  for ( x=0; x<theIm->dim.x; x++ ) {
    
    resBuf[z][y][x] = (unsigned short int)((int)theBuf[z][y][x] - (int)resBuf[z][y][x]);
  }


  return( 1 );
}


















#ifndef NO_PROTO
static void _QSortSkizPointsList( typeSkizPoint* list, int left, int right )
#else
static void _QSortSkizpointsList( list, left, right )
typeSkizPoint* list;
int left;
int right;
#endif
{
    int i, last;
    typeSkizPoint aux;

    if ( left >= right ) return;

    /* on echange le premier point et le point du milieu
       pour obtenir un meilleur pivot
    */
    aux = list[left];   list[left] = list[(left+right)/2];   list[(left+right)/2] = aux;
  
    
    last = left;
    /* on teste les points par rapport au premier point
       a partir du second point 
       s'ils doivent etre classes avant, on les met
       a sa droite, puis on echangera le dernier 
       point deplace (last) avec le premier point
       ainsi on construit deux sous-listes */
    for ( i = left+1; i <= right; i++ )
      if ( list[i].valeur > list[left].valeur ) {
	last ++;
	aux = list[last];   list[last] = list[i];   list[i] = aux;
      }
    /*--- _VT_SwapVPts( &(list[left]), &(list[last]) ); ---*/
    aux = list[left];   list[left] = list[last];   list[last] = aux;

    _QSortSkizPointsList( list, left, last-1 );
    _QSortSkizPointsList( list, last+1, right );
}







#ifndef NO_PROTO
static void _MySortSkizPointsList( typeSkizPoint* list, int maxVal, int left, int right )
#else
static void _MySortSkizpointsList( list, maxVal, left, right )
typeSkizPoint* list;
int max;
int left;
int right;
#endif
{
  char *proc = "_MySortSkizpointsList";
    int last = right;
    int first = left;
    int i;
    int max = maxVal;
    int min = 1;
    typeSkizPoint aux;
    int iter;
    int changes = 0;

    if ( left >= right ) return;


    iter = 0;
    while ( max > min ) {

      changes = 0;
      for ( i=first; i<=last; i++ ) {

	if ( list[i].valeur == max ) {

	  aux = list[i];
	  list[i] = list[first];
	  list[first] = aux;
	  first++;
	  changes ++;

	} else if ( list[i].valeur == min ) {
	  aux = list[i];
	  list[i] = list[last];
	  list[last] = aux;
	  last--;
	  changes ++;
	}
      }

      if ( _VT_VERBOSE_ ) {
	fprintf( stderr,"%s : iteration #%5d  min=%5d max=%5d (%8d -> %8d) changes = %8d\r", 
		 proc, iter+1, min, max, first, last, changes );
      }


      if ( changes == 0 ) {
	max --;
	min ++;
      }

      iter ++;
    }
}

