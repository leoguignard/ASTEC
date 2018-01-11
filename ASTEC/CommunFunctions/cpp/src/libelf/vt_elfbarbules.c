/*****************************************************************************
 * vt_elfbarbules.h - enleve les barbules
 *
 * $Id: vt_elfbarbules.c,v 1.3 2000/09/08 17:53:11 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Feb, 24 2000
 *
 * ADDITIONS, CHANGES
 *
 * Fri Sep  8 19:49:47 MET DST 2000
 * correction d'un bug (ERROR(1) apparaissait)
 * Lors de la disparition d'une jonction (qui ne relie plus que 2 composantes #max et #min)
 * il fallait mettre a jour les voisins de la jonction (eventuelle) a 
 * l'autre 'bout' de #max.
 *
 *
 */

#include <vt_elfbarbules.h>


static int _verbose_ = 2;

typedef enum {
  UNKNOWN,
  COMPOSANTE,
  JONCTION
} enumType;

#define MAXNEIGHBORS 8

typedef struct {

  int label;
  int size;

  int xmin;
  int xmax;
  int ymin;
  int ymax;
  int zmin;
  int zmax;

  enumType type;

  int nbNeighbors;
  int labelsNeighbors[MAXNEIGHBORS];

} typeCC;




int VT_RemoveBarbules( vt_image *theIm,
		       int labelMinCc,
		       int labelMaxCc,
		       int labelMinJc,
		       int labelMaxJc,
		       int length )
{
  char *proc = "VT_RemoveBarbules";
  typeCC *theCC = (typeCC *)NULL;
  int i, j, k, l;
  int x, y, z;
  u16 *** theBuf = (u16 ***)NULL;
  int n, e, b, jct;
  int neighbors[27];
  int min, max;


  if ( labelMinCc > labelMaxCc ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: the smallest label for components (%d) is larger than the largest one (%d).\n", proc, labelMinCc, labelMaxCc );
    }
    return( -1 );
  }

  if ( labelMinJc > labelMaxJc ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: the smallest label for junction (%d) is larger than the largest one (%d).\n", proc, labelMinJc, labelMaxJc );
    }
    return( -1 );
  }

  if ( labelMinCc <= labelMaxJc && labelMinJc <= labelMaxCc ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: there is an overlap between labels of components [%d,%d] and labels of junctions [%d,%d].\n", 
	       proc, labelMinCc, labelMaxCc, labelMinJc, labelMaxJc );
    }
    return( -1 );
  }



  if ( theIm->type != USHORT ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: unable to deal with such image type\n", proc );
    }
    return( -1 );
  }


  n = labelMaxCc;
  if ( n < labelMaxJc ) n = labelMaxJc;
  

  /* theCC est une structure qui contient 
     les composantes et les jonctions
     
     pour chacune, on garde 
     - son label equivalent
     - sa taille
     - sa boite englobante
     - son type 
     - son nombre de voisins et les labels des voisins
  */

  theCC = (typeCC*)malloc( (n+1)*sizeof(typeCC) );
  if ( theCC == (typeCC *)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: can not allocate auxiliary array.\n", proc );
    }
    return( -1 );
  }
  
  for ( i=0; i<=n; i++ ) {

    theCC[i].label = i;
    theCC[i].size = 0;

    theCC[i].xmin = theIm->dim.x - 1;
    theCC[i].xmax = 0;
    theCC[i].ymin = theIm->dim.y - 1;
    theCC[i].ymax = 0;
    theCC[i].zmin = theIm->dim.z - 1;
    theCC[i].zmax = 0;

    theCC[i].type = UNKNOWN;

    theCC[i].nbNeighbors = 0;
    for ( j=0; j<MAXNEIGHBORS; j++ )
      theCC[i].labelsNeighbors[j] = 0;
  }



  if ( _verbose_ )
    fprintf( stderr, " %s: labels of components [%d,%d] and labels of junctions [%d,%d].\n", 
	     proc, labelMinCc, labelMaxCc, labelMinJc, labelMaxJc );



  theBuf = (u16 ***)theIm->array;



  /* on remplit le tableau theCC
   */


  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ )
  for ( x=0; x<theIm->dim.x; x++ ) {
    
    l = theBuf[z][y][x];
    if ( l == 0 ) continue;

    theCC[ l ].size ++;

    if ( theCC[ l ].xmin > x ) theCC[ l ].xmin = x;
    if ( theCC[ l ].xmax < x ) theCC[ l ].xmax = x;
    if ( theCC[ l ].ymin > y ) theCC[ l ].ymin = y;
    if ( theCC[ l ].ymax < y ) theCC[ l ].ymax = y;
    if ( theCC[ l ].zmin > z ) theCC[ l ].zmin = z;
    if ( theCC[ l ].zmax < z ) theCC[ l ].zmax = z;
    
    if ( l >= labelMinCc && l <= labelMaxCc )
      theCC[ l ].type = COMPOSANTE;
    if ( l >= labelMinJc && l <= labelMaxJc )
      theCC[ l ].type = JONCTION;
    
    for ( k=-1; k<=1; k++ ) {
      if ( z+k < 0 || z+k >= theIm->dim.z ) continue;
      for ( j=-1; j<=1; j++ ) {
	if ( y+j < 0 || y+j >= theIm->dim.y ) continue;
	for ( i=-1; i<=1; i++ ) {
	  if ( x+i < 0 || x+i >= theIm->dim.x ) continue;
	  if ( theBuf[z+k][y+j][x+i] == 0 ) continue;
	  if ( theBuf[z+k][y+j][x+i] == l ) continue;
	  /* on a un voisin d'un numero different
	   */
	  e = 0;
	  for ( n=0; n < theCC[ l ].nbNeighbors; n++ )
	    if ( theCC[ l ].labelsNeighbors[n] == theBuf[z+k][y+j][x+i] )
	      e = 1;
	  if ( e == 1 ) continue;
	  /* le voisin n'a pas ete reference
	   */
	  theCC[ l ].labelsNeighbors[ theCC[ l ].nbNeighbors ] = theBuf[z+k][y+j][x+i];
	  theCC[ l ].nbNeighbors ++;
	}
      }
    }

  }












  /* OK,
     ici on a tout pour traiter les barbules
  */
  if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );




  


  /* boucle sur les composantes
   */

  for ( b = labelMaxCc; (b >= labelMinCc) && (theCC[ b ].size <= length); b-- ) {
  
    /* on ne peut traiter que celles qui sont independantes
     */
    if ( theCC[ b ].label != b ) continue;

    /* on ne peut traiter que celles qui sont terminales
     */
    if ( theCC[ b ].nbNeighbors != 1 ) continue;



    /* bon on y va
     */
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, " %s: del cpt #%05d", proc, b );
    }

    theCC[ b ].label = 0;
    for ( z=theCC[ b ].zmin; z<=theCC[ b ].zmax; z++ )
    for ( y=theCC[ b ].ymin; y<=theCC[ b ].ymax; y++ )
    for ( x=theCC[ b ].xmin; x<=theCC[ b ].xmax; x++ )
      if ( theBuf[z][y][x] == b ) theBuf[z][y][x] = 0;
    





    /* et la jonction ?
       1. on supprime #b de ses voisins
     */
    jct = theCC[ b ].labelsNeighbors[0];
    e = -1;
    for ( n=0; n<theCC[ jct ].nbNeighbors; n++ )
      if ( theCC[ jct ].labelsNeighbors[n] == b )
	e = n;
    if ( e == -1 ) {
      if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );
      fprintf( stderr, " ... ERROR (1): %d was not found among the neigbors of junction %d\n",
	       b, jct );
    } else {
      for ( n=e; n<theCC[ jct ].nbNeighbors-1; n++ )
	theCC[ jct ].labelsNeighbors[n] = theCC[ jct ].labelsNeighbors[n+1];
      theCC[ jct ].nbNeighbors --;
    }





    /* et la jonction ?
       2. est-ce toujours une jonction ?
    */
    if ( theCC[ jct ].nbNeighbors != 2 ) {
      if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );
      continue; 
    }

    /* et la jonction ?
       3. il faut fusionner deux composantes
          3.1 on supprime les points simples de la jonction
    */
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, " .. del jct #%05d", jct );
    }

    e = 0;
    for ( z=theCC[ jct ].zmin; z<=theCC[ jct ].zmax; z++ )
    for ( y=theCC[ jct ].ymin; y<=theCC[ jct ].ymax; y++ )
    for ( x=theCC[ jct ].xmin; x<=theCC[ jct ].xmax; x++ ) {
      
      if ( theBuf[z][y][x] != jct ) continue;
      
      /* on extrait le voisinage
       */
      for ( n=0, k=-1; k<=1; k++ )
      for ( j=-1; j<=1; j++ )
      for ( i=-1; i<=1; i++, n++ ) {
	neighbors[n] = 0;
	if ( z+k < 0 || z+k >= theIm->dim.z ) continue;
	if ( y+j < 0 || y+j >= theIm->dim.y ) continue;
	if ( x+i < 0 || x+i >= theIm->dim.x ) continue;
	neighbors[n] = theBuf[z+k][y+j][x+i];
      }

      if ( IsA3DPointSimple( neighbors ) == 1 ) {
	theBuf[z][y][x] = 0;
	e ++;
      }
    }

    if ( e > 0 ) {
      if ( _verbose_ >= 2 ) {
	fprintf( stderr, " (delete %d/%d pts in jct #%05d)", 
		 e, theCC[ jct ].size, jct );
      }
    }
    theCC[ jct ].size -= e;
    
    /* bon, faut fusionner les 2 composantes et la jonction
       => on prend le plus petit label
     */
    if ( theCC[ theCC[ jct ].labelsNeighbors[0] ].label 
	 > theCC[ theCC[ jct ].labelsNeighbors[1] ].label ) {
      max = theCC[ theCC[ jct ].labelsNeighbors[0] ].label;
      min = theCC[ theCC[ jct ].labelsNeighbors[1] ].label;
    } else {
      min = theCC[ theCC[ jct ].labelsNeighbors[0] ].label;
      max = theCC[ theCC[ jct ].labelsNeighbors[1] ].label;
    }
    
    if ( _verbose_ >= 2 ) {
      fprintf( stderr, " ... merging cpts #%05d and #%05d by jct #%05d\n",
	       min, max, jct );
    }

    /* on lance l'equivalence dans les tables
       1. pour les composantes
       2. pour les jonctions

       Fri Sep  8 19:33:32 MET DST 2000, Gregoire Malandain
       je verifie les equivalences par rapport aux voisins
       des jonctions
     */
    for ( i=max; i<=labelMaxCc; i++ ) {
      if ( theCC[i].label == max ) theCC[i].label = min;
    }
    for ( i=labelMinJc; i<=labelMaxJc; i++ ) {
      if ( theCC[i].label == max ) theCC[i].label = min;
    }

    /* on ajoute la taille de la cc #max
     */
    theCC[ min ].size += theCC[ max ].size;
    theCC[ max ].size = 0;

    /* la cc #max avait deux voisins au plus,
       et un au moins (jonction #jct)
       - si c'est deux, le second voisin remplace #jct
         dans les voisins de #min
       - si c'est un, on enleve jct des voisins de #min
    */

    if ( theCC[ max ].nbNeighbors == 2 ) {
      if ( theCC[ max ].labelsNeighbors[0] == jct ) 
	e = theCC[ max ].labelsNeighbors[1];
      else
	e = theCC[ max ].labelsNeighbors[0];
      /* les 2 voisins de #max sont #jct et #e
	 #min a aussi 2 voisins, #jct et autre chose
	 on remplace #jct par #e
	 Modification Fri Sep  8 19:33:32 MET DST 2000
	 et on remplace #max dans les voisins de #e
	 par #min      
       */
      if ( theCC[ min ].labelsNeighbors[0] == jct )
	theCC[ min ].labelsNeighbors[0] = e;
      else
	theCC[ min ].labelsNeighbors[1] = e;
      for ( n=0; n<theCC[ e ].nbNeighbors; n++ ) {
	if ( theCC[ e ].labelsNeighbors[n] == max ) theCC[ e ].labelsNeighbors[n] = min;
      }
      
    } else {
      e = -1;
      for ( n=0; n<theCC[ min ].nbNeighbors; n++ )
	if ( theCC[ min ].labelsNeighbors[n] == jct )
	  e = n;
      if ( e == -1 ) {
	fprintf( stderr, " ... ERROR (2): %d was not found among the neigbors of component %d\n",
		 jct, min );
      } else {
	for ( n=e; n<theCC[ min ].nbNeighbors-1; n++ )
	  theCC[ min ].labelsNeighbors[n] = theCC[ min ].labelsNeighbors[n+1];
	theCC[ min ].nbNeighbors --;
      }
    }
    
    theCC[ jct ].label = min;
    theCC[ min ].size += theCC[ jct ].size;
    theCC[ jct ].size = 0;
    theCC[ jct ].type = COMPOSANTE;
    
  }

  if ( _verbose_ >= 2 ) fprintf( stderr, "\n" );





  /* la on pourrait re-numeroter
     et compter les composantes isolees
  */
  
  j = k = 0;
  if ( labelMinCc < labelMinJc ) {
    for ( i=labelMinCc; i<=labelMaxCc; i++ ) {
      if ( theCC[i].label == i ) 
	theCC[i].label = ++j;
      else 
	theCC[i].label = theCC[ theCC[i].label ].label;
    }
    for ( i=labelMinJc; i<=labelMaxJc; i++ ) {
      if ( theCC[i].label == i ) 
	theCC[i].label = ++k;
      else 
	theCC[i].label = theCC[ theCC[i].label ].label;
    }
  }
  
  if ( _verbose_ ) {
    fprintf( stderr, " %s: there are now %d components (instead of %d)\n", proc, j, labelMaxCc-labelMinCc+1 );
    fprintf( stderr, " \t there are now %d junctions (instead of %d)\n", k, labelMaxJc-labelMinJc+1 );
  }
  
  
    
  free( theCC );
  return( 1 );
}
