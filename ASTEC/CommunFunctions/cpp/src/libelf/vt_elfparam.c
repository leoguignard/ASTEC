/*************************************************************************
 * vt_elfparam.c - extraction de parametres sur des parties numerotees
 *
 * $Id: vt_elfparam.c,v 1.13 2001/02/13 17:50:52 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 * void _ComputeComponentParameters()
 *
 *
 *
 *
 *
 *
 * void _SetNbTestsForMaxFeretEstimation()
 * void _SetFeretDiameterComputationToComputation()
 * void _SetFeretDiameterComputationToEstimation()
 *
 *
 * Pour 10 essais avec test-param-ellipse3D / Version 1
 * ...  temps moyen pour le calcul exhaustif = 6061424.200000 (6.06142 secs)
 * ...  temps moyen pour le calcul estime    = 36665.200000 (0.0366652 secs)
 *  
 * Pour 100 essais avec test-param-ellipse3D / Version 1
 * ...  temps moyen pour le calcul exhaustif = 5600275.980000 (5.60028 secs)
 * ...  temps moyen pour le calcul estime    = 32165.380000 (0.0321654 secs)
 *
 * static double _ComputeMaxDifferenceInList()
 * static double _EstimateMaxSquareDistanceInList()
 *
 * static double _ComputeMaxDifferenceInList()
 *
 *
 *
 * static void _ComputeFeretDiametersIn2DList()
 * static void _ComputeFeretDiametersIn3DList()
 *
 *
 *
 * static void _FillListWithBorderPoints()
 *
 *
 *
 * static void _PrintTypeListe()
 * static void _InitTypeListe()
 * static void _FreeTypeListe()
 * static int _AllocTypeListe()
 *
 *
 *
 *
 * void _GenerateEllipse2D()
 * void _GenerateEllipse3D()
 *
 *
 * static void _MultiplyMatrixByVector()
 * static void _RotationMatrixFromRotationVector()
 *
 *
 *
 * int _MaxValueInImage()
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Jul 20 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef _LINUX_
extern long int random();
extern void srandom(unsigned int seed);
#endif

#include <vt_elfparam.h>


static int _verbose_ = 0;

void _VerboseInElfParam()
{
  _verbose_ = 1;
}
void _NoVerboseInElfParam()
{
  _verbose_ = 0;
}





typedef struct {
  double x;
  double y;
  double z;
} typePoint;

typedef struct {
  int n;
  int nbAllocs;
  typePoint *point; /* point "normal" */
  typePoint *projs; /* point "projetes" */
  double *abscisse;
} typeListe;

typedef enum {
  _ESTIMATE_,
  _COMPUTE_,
} enumMaxFeret;

static enumMaxFeret typeMaxFeret = _ESTIMATE_;
static int _MAXTESTS_FOR_MAXFERET_ESTIMATION_ = 5;
static double _EPSILON_FOR_NORMS_ = 0.0001;



static double _ComputeMaxSquareDistanceInList( typePoint *theList,
					    int n,
					    typePoint *first,
					    typePoint *last );
static double _EstimateMaxSquareDistanceInList( typePoint *theList,
					     int n,
					     typePoint *first,
					     typePoint *last  );
static double _ComputeMaxDifferenceInList( double *theList,
					   int n );



static void _ComputeFeretDiametersIn2DList( typeListe *theList,
					    typeComponentParameters *thePar );

static void _ComputeFeretDiametersIn3DList( typeListe *theList,
					    typeComponentParameters *thePar );



static void _FillListWithBorderPoints( vt_image *theIm,
				       typeListe *theList,
				       int slice,
				       int color,
				       int *corner,
				       int *windowSize );



static void _PrintTypeListe( typeListe *theListe );
static void _InitTypeListe( typeListe *theListe );
static void _FreeTypeListe( typeListe *theListe );
static int _AllocTypeListe( typeListe *theListe, int n );
static void _MultiplyMatrixByVector( const double *m,
				     const double *a,
				     double *x );
static void _RotationMatrixFromRotationVector( double *mat,
					       double *rot );













int _CreateArrayOfParametersFromImage( vt_image *theIm,
				       int slice,
				       typeComponentParameters **thePar )
{
  char *proc="_CreateArrayOfParametersFromImage";
  int i, n;
  int dim[3];
  int s=slice;
  typeComponentParameters *theComp = (typeComponentParameters *)NULL;
  typeListe theList;


  _InitTypeListe( &theList );


  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= theIm->dim.z) ) s = -1;


  /* c'est pas optimal. mais ca ira bien pour l'instant
     on alloue n+1 composantes pour adresser chaque
     element par son intensite
   */
  n = _MaxValueInImage( theIm, s );
  if ( n == 0 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: null image.\n", proc );
    }
    return( 0 );
  }

  theComp = (typeComponentParameters *)malloc( (n+1)*sizeof(typeComponentParameters) );
  if ( theComp == (typeComponentParameters *)NULL ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate components list.\n", proc );
    }
    return( -1 );
  }

  _InitArrayOfParametersFromImage( theIm, theComp, n );
  _FillArrayOfParametersFromImage( theIm, s, theComp );
  

  for ( i=1; i<=n; i++ ) {
    /* ca c'est le volume de la sphere
     */
    theComp[i].equivSphereRadiusFromVolume = cbrt( 3.0 * (double)theComp[i].volume
						   / (4.0 * 3.14159) );
    dim[0] = theComp[i].ptmax[0] - theComp[i].ptmin[0] + 1;
    dim[1] = theComp[i].ptmax[1] - theComp[i].ptmin[1] + 1;
    dim[2] = theComp[i].ptmax[2] - theComp[i].ptmin[2] + 1;


    _FillListWithBorderPoints( theIm, &theList, s, i, theComp[i].ptmin, dim );

    if ( theIm->dim.z == 1 || s >= 0 ) 
      _ComputeFeretDiametersIn2DList( &theList, &theComp[i] );
    else
      _ComputeFeretDiametersIn3DList( &theList, &theComp[i] );


    /* ca c'est la surface de la sphere
     */
    theComp[i].border = theList.n;
    theComp[i].equivSphereRadiusFromSurface = sqrt( (double)theList.n / (4.0 * 3.14159) );

  }


  
  _FreeTypeListe( &theList );

  *thePar = theComp;
  return( n );
}

















/* _ComputeComponentParameters()
 
   extrait les points frontieres en 2D (si la dimension en Z
   est egale a 1 ou si l'indice de coupe 'slice' est dans les
   bornes) ou en 3D,
   puis calcule les diametres de Feret en 2D ou en 3D
*/

void _ComputeComponentParameters( vt_image *theIm,
				  typeComponentParameters *thePar,
				  int slice,
				  int color,
				  int *corner,
				  int *windowSize )
{
  typeListe theList;
  int s=slice;

  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= theIm->dim.z) ) s = -1;
 

  _InitTypeListe( &theList );
  _FillListWithBorderPoints( theIm, &theList, s, color, corner, windowSize );

  if ( 0 ) _PrintTypeListe( &theList );
  

  /* cas 2D
   */
  if ( theIm->dim.z == 1 || s >= 0 ) {

    _ComputeFeretDiametersIn2DList( &theList, thePar );

  } else {

    _ComputeFeretDiametersIn3DList( &theList, thePar );

  }

  _FreeTypeListe( &theList );
}


















void _SetNbTestsForMaxFeretEstimation( int n ) 
{
  if ( n > 0 )
    _MAXTESTS_FOR_MAXFERET_ESTIMATION_ = n;
}
void _SetFeretDiameterComputationToComputation()
{
  typeMaxFeret = _COMPUTE_;
}
void _SetFeretDiameterComputationToEstimation()
{
  typeMaxFeret = _ESTIMATE_;
}



















static double _ComputeMaxSquareDistanceInList( typePoint *theList,
					    int n,
					    typePoint *first,
					    typePoint *last )
{
  int i, j;
  double max, d;

  if ( n == 1 ) {
    *first = *last = theList[0];
    return( (double)0.0 );
  }

  max = d = 0.0;
  for ( i=0; i<n-1; i++ ) 
  for ( j=i+1; j<n; j++ ) {
    d = (theList[i].z - theList[j].z) * (theList[i].z - theList[j].z) 
      + (theList[i].y - theList[j].y) * (theList[i].y - theList[j].y) 
      + (theList[i].x - theList[j].x) * (theList[i].x - theList[j].x);
    if ( d > max ) {
      max = d;
      *first = theList[i];
      *last  = theList[j];
    }
  }

  return( max );
}




















































static double _EstimateMaxSquareDistanceInList( typePoint *theList,
						int nbpts,
						typePoint *first,
						typePoint *last  )
{
  int i, i2, j, n;
  typePoint p1, p2;
  typePoint pf, pl;
  double max, old, m, d;
  int endIsReached = 0;
  int nb=nbpts;
  
  
  if ( nb == 1 ) {
    *first = *last = theList[0];
    return( (double)0.0 );
  } 

  max = 0;
  

  /* on tire un point au hasard dans les n premiers
     on cherche le plus eloigne
     on echange les points
     on itere jusqu'a stabilite
  */

  i2 = (int)( (nb - 1) * ((double)random() / (double)2147483647.0) + 0.5);

  do {

    m = 0;

    do {

      p1 = theList[i2];  
      theList[i2] = theList[nb-1]; theList[nb-1] = p1;
      nb --;
      
      old = m;
      for ( i=0; i<nb; i++ ) {
	d = (theList[i].z - p1.z) * (theList[i].z - p1.z) 
	  + (theList[i].y - p1.y) * (theList[i].y - p1.y) 
	  + (theList[i].x - p1.x) * (theList[i].x - p1.x);
	if ( d > m ) {
	  m = d;
	  p2 = theList[i];
	  i2 = i;
	  pf = p1;
	  pl = p2;
	}
      }

    } while ( m > old );
    
  
    if ( m > max ) {
      /* on a un maximum
       */
      max    = m;
      *first = pf;
      *last  = pl;
      

      /* on cherche le point le plus eloigne du cercle
       */
      i2 = -1;
      m = 0.0;
      for ( i=0; i<nb; i++ ) {
	d = (theList[i].z - pf.z) * (theList[i].z - pl.z) 
	  + (theList[i].y - pf.y) * (theList[i].y - pl.y) 
	  + (theList[i].x - pf.x) * (theList[i].x - pl.x);
	if ( m < d ) {
	  i2 = i;
	  m  = d;
	}
      }

      if ( i2 < 0 ) return( max );
      
    } else {

      endIsReached = 1;
    }

  } while( endIsReached == 0 );
  

  /* on met les points en dehors du cercle
     devant la liste
  */
  n = -1;
  for ( i=0; i<nb; i++ ) {
    d = (theList[i].z - pf.z) * (theList[i].z - pl.z) 
      + (theList[i].y - pf.y) * (theList[i].y - pl.y) 
      + (theList[i].x - pf.x) * (theList[i].x - pl.x);
    if ( d > 0 ) {
      n ++;
      p1 = theList[n]; theList[n] = theList[i]; theList[i] = p1;
    }
  }
  
  if (n <0) return( max );


  

  /* on finit par une recherche exhaustive
   */
  for ( i=0; i<=n; i++ ) 
  for ( j=i+1; j<nb; j++ ) {
    d = (theList[i].z - theList[j].z) * (theList[i].z - theList[j].z) 
      + (theList[i].y - theList[j].y) * (theList[i].y - theList[j].y) 
      + (theList[i].x - theList[j].x) * (theList[i].x - theList[j].x);
    if ( d > max ) {
      max = d;
      *first = theList[i];
      *last  = theList[j];
    }
  }

  return( max );
}






#ifdef _UNUSED_
static double _OLD_EstimateMaxSquareDistanceInList( typePoint *theList,
						int nb,
						typePoint *first,
						typePoint *last  )
{
  int i, j, n;
  typePoint p1, p2;
  typePoint c;
  typePoint pf, pl;
  double max, old, m, d;
  int endIsReached = 0;

  if ( nb == 1 ) {
    *first = *last = theList[0];
    return( (double)0.0 );
  } 

  max = 0;
  n = nb;

  do {

    /* on tire un point au hasard dans les n premiers
       on cherche le plus eloigne
       on echange les points
       on itere jusqu'a stabilite
    */
    i = (int)( (n - 1) * ((double)random() / (double)2147483647.0) + 0.5);
    p1 = theList[i];  

    m = 0;
    do {
      old = m;
      for ( i=0; i<nb; i++ ) {
	d = (theList[i].z - p1.z) * (theList[i].z - p1.z) 
	  + (theList[i].y - p1.y) * (theList[i].y - p1.y) 
	  + (theList[i].x - p1.x) * (theList[i].x - p1.x);
	if ( d > m ) {
	  m = d;
	  p2 = theList[i];
	  pf = p1;
	  pl = p2;
	}
      }
      p1 = p2;
    } while ( m > old );
    
    if ( m > max ) {
      /* on a un maximum
       */
      max    = m;
      *first = pf;
      *last  = pl;
      
      /* si PF-PL n'est pas le plus grand diametre,
	 alors une extremite de celui-ci 
	 est a l'exterieur du cercle de diametre PF-PL
      */
      c.x = (pf.x + pl.x) / 2.0;
      c.y = (pf.y + pl.y) / 2.0;
      c.z = (pf.z + pl.z) / 2.0;
      m /= 4.0;
      n = nb;
      for ( i=0; i<n; i++ ) {
	d = (theList[i].z - c.z) * (theList[i].z - c.z) 
	  + (theList[i].y - c.y) * (theList[i].y - c.y) 
	  + (theList[i].x - c.x) * (theList[i].x - c.x);
	if ( d <= m ) {
	  p1 = theList[i];
	  theList[i] = theList[n-1];
	  theList[n-1] = p1;
	  n--;
	  i--;
	}
      }

      /*
	fprintf( stderr, " ... max = %f (R=%f) ... reste %d points\n", max, sqrt(max)/2.0, n );
      */

      /* il y a n points en dehors du cercle 
	 on reitere, le prochain points sera tire au hasard parmi ces n.
       */
    } else {
      endIsReached = 1;
    }

  } while( (n > 0) && (endIsReached == 0) );
  
  /* on finit par une recherche exhaustive
   */
  if ( n > 0 ) {
    for ( i=0; i<n; i++ ) 
    for ( j=i+1; j<nb; j++ ) {
      d = (theList[i].z - theList[j].z) * (theList[i].z - theList[j].z) 
	+ (theList[i].y - theList[j].y) * (theList[i].y - theList[j].y) 
	+ (theList[i].x - theList[j].x) * (theList[i].x - theList[j].x);
      if ( d > max ) {
	max = d;
	*first = theList[i];
	*last  = theList[j];
      }
    }
  }
  

  return( max );
}
#endif





















static double _ComputeMaxDifferenceInList( double *theList,
					   int n )
{
  int i;
  double max, min;

  min = max = theList[0];
  for ( i=1; i<n; i++ ) {
    if ( theList[i] < min ) min = theList[i];
    else if ( theList[i] > max ) max = theList[i];
  }
  return( max - min );
}











/* calcul des diametres de Feret en 2D
 *
 * On suppose que la composante en Z des points
 * est la meme pour tous les points et donc on ne
 * regarde que les composantes en X et Y pour
 * le calcul de la distance
 *
 * On calcule le diametre maximum de Feret comme etant
 * la distance maximale entre 2 points de la liste
 * (qui represente la 6-frontiere de l'objet) etant
 * donne que l'on peut montrer que ce diametre est
 * effectivement une distance entre 2 points.
 *
 * Soit un diametre de Feret (distance entre 2 plans paralleles
 * tels que chacun des plans soit d'un cote de l'objet et touche
 * l'objet) : si les deux points de contact ne forme pas une orthogonale
 * aux plans, alors on peut trouver un diametre de Feret plus grand
 * (la distance entre les points de contact) et donc ce premier diametre
 * n'etait pas maximum.
 *
 * Le second diametre est le diametre orthogonal a ce diametre maximal
 * ce N'EST PAS le diametre minimum au sens de Feret, puisque les diametres
 * minimum et maximum n'ont pas de raisons d'etre orthogonaux.
 * 
 */
static void _ComputeFeretDiametersIn2DList( typeListe *theList,
					    typeComponentParameters *thePar )
{
  int i;
  typePoint pt1, pt2;
  double n, max=0;

  switch ( typeMaxFeret ) {
  default :
  case _COMPUTE_ :
    max = _ComputeMaxSquareDistanceInList( theList->point, theList->n,
					   &pt1, &pt2 );
    break;
  case _ESTIMATE_ :
    max = _EstimateMaxSquareDistanceInList( theList->point, theList->n,
					    &pt1, &pt2 );
  }

  /* directions
   */
  thePar->maxDirection[0] = (double)(pt2.x - pt1.x);
  thePar->maxDirection[1] = (double)(pt2.y - pt1.y);
  thePar->maxDirection[2] = 0.0;
  
  n = thePar->maxDirection[0]*thePar->maxDirection[0] 
    + thePar->maxDirection[1]*thePar->maxDirection[1];

  /* c'est un point isole
   */
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->maxDiameter = 0.0;
    thePar->maxDirection[0] = thePar->maxDirection[1] = thePar->maxDirection[2] = 0.0;
    thePar->medDiameter = 0.0;
    thePar->medDirection[0] = thePar->medDirection[1] = thePar->medDirection[2] = 0.0;
    thePar->minDiameter = 0.0;
    thePar->minDirection[0] = thePar->minDirection[1] = thePar->minDirection[2] = 0.0;
    return;
  }



  n = sqrt( n );

  thePar->maxDirection[0] /= n;
  thePar->maxDirection[1] /= n;
  
  thePar->maxDiameter = n;

  thePar->minDirection[0] = -(thePar->maxDirection[1]);
  thePar->minDirection[1] =   thePar->maxDirection[0];
  thePar->minDirection[2] = 0.0;

  /* abscisses curvilignes
   */
  for ( i=0; i<theList->n; i++ )
    theList->abscisse[i] = theList->point[i].x * thePar->minDirection[0] 
      + theList->point[i].y * thePar->minDirection[1];
  
  
  thePar->minDiameter = _ComputeMaxDifferenceInList( theList->abscisse, 
						   theList->n );

  /*
   */
  thePar->medDiameter = 0.0;
  thePar->medDirection[0] = thePar->medDirection[1] = thePar->medDirection[2] = 0.0;
}















static void _ComputeFeretDiametersIn3DList( typeListe *theList,
					    typeComponentParameters *thePar )
{
  typePoint pt1, pt2;
  double n, max=0;
  int i;
  double ps;
  
  switch ( typeMaxFeret ) {
  default :
  case _COMPUTE_ :
    max = _ComputeMaxSquareDistanceInList( theList->point, theList->n,
					   &pt1, &pt2 );
    break;
  case _ESTIMATE_ :
    max = _EstimateMaxSquareDistanceInList( theList->point, theList->n,
					    &pt1, &pt2 );
    break;
  }


  /* directions
   */
  thePar->maxDirection[0] = (double)(pt2.x - pt1.x);
  thePar->maxDirection[1] = (double)(pt2.y - pt1.y);
  thePar->maxDirection[2] = (double)(pt2.z - pt1.z);

  n = thePar->maxDirection[0]*thePar->maxDirection[0] 
    + thePar->maxDirection[1]*thePar->maxDirection[1] 
    + thePar->maxDirection[2]*thePar->maxDirection[2];

  if ( _verbose_ ) {
    printf("MAX = (%9.5g %9.5g %9.5g) - (%9.5g %9.5g %9.5g) = %g\n",
	   pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z, n );
  }


  /* c'est un point isole
   */
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->maxDiameter = 0.0;
    thePar->maxDirection[0] = thePar->maxDirection[1] = thePar->maxDirection[2] = 0.0;
    thePar->medDiameter = 0.0;
    thePar->medDirection[0] = thePar->medDirection[1] = thePar->medDirection[2] = 0.0;
    thePar->minDiameter = 0.0;
    thePar->minDirection[0] = thePar->minDirection[1] = thePar->minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );
  thePar->maxDirection[0] /= n;
  thePar->maxDirection[1] /= n;
  thePar->maxDirection[2] /= n;

  thePar->maxDiameter = n;

  if ( _verbose_ ) {
    printf( "   DIR = (%9.5g %9.5g %9.5g)\n", 
	    thePar->maxDirection[0], thePar->maxDirection[1], thePar->maxDirection[2] );
  }


  /* on enleve la projection du point sur la premiere direction
     ainsi tous les points seront coplanaires
  */
  for ( i=0; i<theList->n; i++ ) {
    ps = theList->point[i].x * thePar->maxDirection[0] +
         theList->point[i].y * thePar->maxDirection[1] +
         theList->point[i].z * thePar->maxDirection[2];
    theList->projs[i].x = theList->point[i].x - ps * thePar->maxDirection[0];
    theList->projs[i].y = theList->point[i].y - ps * thePar->maxDirection[1];
    theList->projs[i].z = theList->point[i].z - ps * thePar->maxDirection[2];
  }


  switch ( typeMaxFeret ) {
  default :
  case _COMPUTE_ :
    max = _ComputeMaxSquareDistanceInList( theList->projs, theList->n,
					   &pt1, &pt2 );
    break;
  case _ESTIMATE_ :
    max = _EstimateMaxSquareDistanceInList( theList->projs, theList->n,
					    &pt1, &pt2 );
    break;
  }

  thePar->medDirection[0] = (double)(pt2.x - pt1.x);
  thePar->medDirection[1] = (double)(pt2.y - pt1.y);
  thePar->medDirection[2] = (double)(pt2.z - pt1.z);

  n = thePar->medDirection[0]*thePar->medDirection[0] 
    + thePar->medDirection[1]*thePar->medDirection[1] 
    + thePar->medDirection[2]*thePar->medDirection[2];

  if ( _verbose_ ) {
    printf("MED = (%9.5g %9.5g %9.5g) - (%9.5g %9.5g %9.5g) = %g\n",
	   pt1.x, pt1.y, pt1.z, pt2.x, pt2.y, pt2.z, n );
  }

  /* c'est un segment de droite
   */
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->medDiameter = 0.0;
    thePar->medDirection[0] = thePar->medDirection[1] = thePar->medDirection[2] = 0.0;
    thePar->minDiameter = 0.0;
    thePar->minDirection[0] = thePar->minDirection[1] = thePar->minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );
  thePar->medDirection[0] /= n;
  thePar->medDirection[1] /= n;
  thePar->medDirection[2] /= n;

  thePar->medDiameter = n;

  if ( _verbose_ ) {
    printf( "   DIR = (%9.5g %9.5g %9.5g)\n", 
	    thePar->medDirection[0], thePar->medDirection[1], thePar->medDirection[2] );
  }


  /* 3eme direction
     c'est le produit vectoriel des deux autres
   */
  thePar->minDirection[0] = thePar->maxDirection[1] * thePar->medDirection[2]
                          - thePar->maxDirection[2] * thePar->medDirection[1];
  thePar->minDirection[1] = thePar->maxDirection[2] * thePar->medDirection[0]
                          - thePar->maxDirection[0] * thePar->medDirection[2];
  thePar->minDirection[2] = thePar->maxDirection[0] * thePar->medDirection[1]
                          - thePar->maxDirection[1] * thePar->medDirection[0];


  n = thePar->minDirection[0]*thePar->minDirection[0] 
    + thePar->minDirection[1]*thePar->minDirection[1] 
    + thePar->minDirection[2]*thePar->minDirection[2];
  
  if ( n <= _EPSILON_FOR_NORMS_ ) {
    thePar->minDiameter = 0.0;
    thePar->minDirection[0] = thePar->minDirection[1] = thePar->minDirection[2] = 0.0;
    return;
  }

  n = sqrt( n );

  thePar->minDirection[0] /= n;
  thePar->minDirection[1] /= n;
  thePar->minDirection[2] /= n;

  if ( _verbose_ ) {
    printf( "MINDIR = (%9.5g %9.5g %9.5g)\n", 
	    thePar->minDirection[0], thePar->minDirection[1], thePar->minDirection[2] );
  }

  /* abscisses curvilignes
   */
  for ( i=0; i<theList->n; i++ )
    theList->abscisse[i] = theList->point[i].x * thePar->minDirection[0] 
      + theList->point[i].y * thePar->minDirection[1] 
      + theList->point[i].z * thePar->minDirection[2];
  
  thePar->minDiameter = _ComputeMaxDifferenceInList( theList->abscisse, 
						     theList->n );

}





















static void _FillListWithBorderPoints( vt_image *theIm,
				       typeListe *theList,
				       int slice,
				       int color,
				       int *corner,
				       int *windowSize )
{
  char *proc ="_FillListWithBorderPoints";
  int x, y, z;
  int i, s=slice;
  typePoint theCorner;
  typePoint theSize;
  


  theList->n = 0;



  theCorner.x = corner[0];
  theCorner.y = corner[1];
  theCorner.z = corner[2];
  if ( theCorner.x < 0 ) theCorner.x = 0;
  if ( theCorner.y < 0 ) theCorner.y = 0;
  if ( theCorner.z < 0 ) theCorner.z = 0;
  if ( theCorner.x >= theIm->dim.x ) theCorner.x = theIm->dim.x - 1;
  if ( theCorner.y >= theIm->dim.y ) theCorner.y = theIm->dim.y - 1;
  if ( theCorner.z >= theIm->dim.z ) theCorner.z = theIm->dim.z - 1;

  theSize.x = windowSize[0];
  theSize.y = windowSize[1];
  theSize.z = windowSize[2];
  if ( theSize.x < 1 ) theSize.x = 1;
  if ( theSize.y < 1 ) theSize.y = 1;
  if ( theSize.z < 1 ) theSize.z = 1;
  if ( theCorner.x+theSize.x > theIm->dim.x ) theSize.x = theIm->dim.x - theCorner.x;
  if ( theCorner.y+theSize.y > theIm->dim.y ) theSize.y = theIm->dim.y - theCorner.y;
  if ( theCorner.z+theSize.z > theIm->dim.z ) theSize.z = theIm->dim.z - theCorner.z;

  
  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= theIm->dim.z) ) s = -1;
  /* restent les cas 
     (theIm->dim.z > 1) && (s >=0 0 && s < theIm->dim.z)
     qui sont les cas 2D
  */

  if ( 0 ) {
    if ( (theIm->dim.z == 1) || (s >= 0) ) {
      fprintf( stderr, "%s: 2D case\n", proc );
    } else {
      fprintf( stderr, "%s: 3D case\n", proc );
    }
  }

  if ( _AllocTypeListe( theList, theSize.x*theSize.y*theSize.z ) != 1 ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate points list.\n", proc );
    }
    return;
  }



  switch ( theIm->type ) {
  case USHORT :
    {
      unsigned short ***theBuf = (unsigned short ***)theIm->array;
      
      /* cas 2D
       */
      if ( theIm->dim.z == 1 || s >= 0 ) {

	if ( theCorner.x >= 1 && theCorner.y >= 1 && 
	     theCorner.y+theSize.y < theIm->dim.y &&
	     theCorner.x+theSize.x < theIm->dim.x ) {
	  /* 2D: ici on ne fait pas de tests
	   */

	  for ( y=theCorner.y; y<theCorner.y+theSize.y; y++ )
	  for ( x=theCorner.x; x<theCorner.x+theSize.x; x++ ) {
	    if ( theBuf[s][y][x] != color ) continue;
	    if ( theBuf[s][y][x-1] != color ||
		 theBuf[s][y][x+1] != color ||
		 theBuf[s][y-1][x] != color ||
		 theBuf[s][y+1][x] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->n ++;
	    }
	  }

	} else {
	  /* 2D: ici on fait des tests
	   */

	  for ( y=theCorner.y; y<theCorner.y+theSize.y; y++ )
	  for ( x=theCorner.x; x<theCorner.x+theSize.x; x++ ) {
	    if ( theBuf[s][y][x] != color ) continue;
	    if ( x-1 >= 0 ) {
	      if ( theBuf[s][y][x-1] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->n ++;
		continue;
	      }
	    }
	    if ( x+1 < theIm->dim.x ) {
	      if ( theBuf[s][y][x+1] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->n ++;
		continue;
	      }
	    }
	    if ( y-1 >= 0 ) {
	      if ( theBuf[s][y-1][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->n ++;
		continue;
	      }
	    }
	    if ( y+1 < theIm->dim.y ) {
	      if ( theBuf[s][y+1][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->n ++;
		continue;
	      }
	    }
	  }

	}

	for ( i=0; i<theList->n; i++ )
	  theList->point[i].z = 0;

	return;
      } /* fin du cas 2D */

      

      if ( 0 ) {
	fprintf( stderr, "%s: can not deal with 3D image.\n", proc );
	return;
      }

      /* cas 3D
       */

	if ( theCorner.x >= 1 && theCorner.y >= 1 && theCorner.z >= 1 && 
	     theCorner.z+theSize.z < theIm->dim.z && 
	     theCorner.y+theSize.y < theIm->dim.y &&
	     theCorner.x+theSize.x < theIm->dim.x ) {
	  /* 3D: ici on ne fait pas de tests
	   */
	  
	  for ( z=theCorner.z; z<theCorner.z+theSize.z; z++ )
	  for ( y=theCorner.y; y<theCorner.y+theSize.y; y++ )
	  for ( x=theCorner.x; x<theCorner.x+theSize.x; x++ ) {
	    if ( theBuf[z][y][x] != color ) continue;
	    if ( theBuf[z][y][x-1] != color ||
		 theBuf[z][y][x+1] != color ||
		 theBuf[z][y-1][x] != color ||
		 theBuf[z][y+1][x] != color ||
		 theBuf[z-1][y][x] != color ||
		 theBuf[z+1][y][x] != color ) {
	      theList->point[theList->n].x = x;
	      theList->point[theList->n].y = y;
	      theList->point[theList->n].z = z;
	      theList->n ++;
	    }
	  }

	} else {
	  /* 3D: ici on fait des tests
	   */

	  for ( z=theCorner.z; z<theCorner.z+theSize.z; z++ )
	  for ( y=theCorner.y; y<theCorner.y+theSize.y; y++ )
	  for ( x=theCorner.x; x<theCorner.x+theSize.x; x++ ) {
	    if ( theBuf[z][y][x] != color ) continue;
	    if ( x-1 >= 0 ) {
	      if ( theBuf[z][y][x-1] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	    if ( x+1 < theIm->dim.x ) {
	      if ( theBuf[z][y][x+1] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	    if ( y-1 >= 0 ) {
	      if ( theBuf[z][y-1][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	    if ( y+1 < theIm->dim.y ) {
	      if ( theBuf[z][y+1][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	    if ( z-1 >= 0 ) {
	      if ( theBuf[z-1][y][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	    if ( z+1 < theIm->dim.z ) {
	      if ( theBuf[z+1][y][x] != color ) {
		theList->point[theList->n].x = x;
		theList->point[theList->n].y = y;
		theList->point[theList->n].z = z;
		theList->n ++;
		continue;
	      }
	    }
	  }

	}
    }
    break;


  default :
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not deal with such image type.\n", proc );
    }
    return;
  }

  return;
  
}
















static void _PrintTypeListe( typeListe *theListe )
{
  int i;
  for ( i=0; i<theListe->n; i++ ) {
    printf("%5d : %5.2f %5.2f %5.2f\n", i, theListe->point[i].x, 
	   theListe->point[i].y, theListe->point[i].z );
  }
}







static void _InitTypeListe( typeListe *theListe )
{
  theListe->n = 0;
  theListe->nbAllocs = 0;
  theListe->point = (typePoint *)NULL;
  theListe->projs = (typePoint *)NULL;
  theListe->abscisse = (double*)NULL;
}








static void _FreeTypeListe( typeListe *theListe )
{
  if ( theListe->nbAllocs > 0 ) {
    if ( theListe->point != (typePoint *)NULL ) 
      free( theListe->point );
    if ( theListe->projs != (typePoint *)NULL ) 
      free( theListe->projs );
    if ( theListe->abscisse != (double*)NULL )
      free( theListe->abscisse );
  }
  theListe->n = 0;
  theListe->nbAllocs = 0;
  theListe->point = (typePoint *)NULL;
  theListe->projs = (typePoint *)NULL;
  theListe->abscisse = (double*)NULL;
}







static int _AllocTypeListe( typeListe *theListe, int n )
{
  char *proc = "_AllocTypeListe";
  typePoint *p = (typePoint *)NULL;
  double *d = (double*)NULL;

  
  if ( n <= theListe->nbAllocs ) return( 1 );
  if ( n <= 0 )  return( 1 );


  p = (typePoint *)malloc( n * sizeof( typePoint ) );
  if ( p == (typePoint *)NULL ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate 1st new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->point != (typePoint *)NULL ) 
    free( theListe->point );
  
  theListe->point = p;


  
  p = (typePoint *)malloc( n * sizeof( typePoint ) );
  if ( p == (typePoint *)NULL ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate 2nd new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->projs != (typePoint *)NULL ) 
    free( theListe->projs );
  
  theListe->projs = p;



  d = (double*)malloc( n * sizeof( double ) );
  if ( d == (double*)NULL ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not allocate 3rd new points list (size=%d).\n", 
	       proc, n );
    }
    return( 0 );
  }
  if ( theListe->nbAllocs > 0 && theListe->abscisse != (double *)NULL ) 
    free( theListe->abscisse );
  
  theListe->abscisse = d;

  theListe->nbAllocs = n;
  return( 1 );
}








void _GenerateEllipse2D( vt_image *theIm,
			 int slice,
			 int color,
			 double rmax,
			 double rmin,
			 double *radius,
			 double *angle )
{
  char *proc = "_GenerateEllipse2D";
  int ix, iy;
  double t;
  double cx, cy, tx, ty;
  double x, y;
  double mat[4];
  int i;
  unsigned short int ***theBuf = (unsigned short int ***)theIm->array;

  
  if ( theIm->type != USHORT ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not deal with such image type.\n", proc );
    }
    return;
  }


  /* le centre de l'ellipse c'est le centre de l'image
   */
  cx = (double)(theIm->dim.x - 1) / 2.0;
  cy = (double)(theIm->dim.y - 1) / 2.0;



  /* on tire les rayons de l'ellipse au hasard
     entre rmin et rmax 
  */
  for ( i=0; i<2; i++ ) 
    radius[i] = rmin + (rmax - rmin) * ((double)random() / (double)2147483647.0);
  if ( radius[1] > radius[0] ) {
    t = radius[0]; 
    radius[0] = radius[1];
    radius[1] = t;
  }



  /* on tire un angle de rotation au hasard
     c'est en degre
   */
  *angle = 360.0 * ((double)random() / (double)2147483647.0);
  t = (*angle) * (3.1415927 / 180.0);


  
  /* la matrice de rotation d'angle t
   */
  mat[0] = cos( t );   mat[1] = -sin( t );
  mat[2] = -mat[1];    mat[3] = mat[0];



  /* la translation qui conserve le centre / M*C + t = C
     => t = C - M*C
   */
  tx = cx - ( mat[0]*cx + mat[1]*cy );
  ty = cy - ( mat[2]*cx + mat[3]*cy );
  



  for ( iy=0; iy<theIm->dim.y; iy++ )
  for ( ix=0; ix<theIm->dim.x; ix++ ) {
    /* coordonnees par rapport au centre 
     */
    x = mat[0]*(double)ix + mat[1]*(double)iy + tx - cx;
    y = mat[2]*(double)ix + mat[3]*(double)iy + ty - cy;
    if ( (x*x)/(radius[0]*radius[0]) +
	 (y*y)/(radius[1]*radius[1]) <= 1.0 )
      theBuf[slice][iy][ix] = color;
    else 
      theBuf[slice][iy][ix] = 0;
  }

}




















void _GenerateEllipse3D( vt_image *theIm,
			 int color,
			 double rmax,
			 double rmin,
			 double *radius,
			 double *angle )
{
  char *proc = "_GenerateEllipse3D";
  int ix, iy, iz;
  double t;
  double cx, cy, cz, tx, ty, tz;
  double x, y, z;
  double r[3];
  double mat[9];
  int i;
  unsigned short int ***theBuf = (unsigned short int ***)theIm->array;

  
  if ( theIm->type != USHORT ) {
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) {
      fprintf( stderr, "%s: can not deal with such image type.\n", proc );
    }
    return;
  }


  /* le centre de l'ellipse c'est le centre de l'image
   */
  cx = (double)(theIm->dim.x - 1) / 2.0;
  cy = (double)(theIm->dim.y - 1) / 2.0;
  cz = (double)(theIm->dim.y - 1) / 2.0;



  /* on tire les rayons de l'ellipse au hasard
     entre rmin et rmax 
  */
  for ( i=0; i<3; i++ ) 
    radius[i] = rmin + (rmax - rmin) * ((double)random() / (double)2147483647.0);
  
  if ( radius[1] >= radius[0] && radius[1] >= radius[2] ) {
    t = radius[0]; 
    radius[0] = radius[1];
    radius[1] = t;
  } else if ( radius[2] >= radius[0] && radius[2] >= radius[1] ) {
    t = radius[0]; 
    radius[0] = radius[2];
    radius[2] = t;
  }

  if ( radius[2] >= radius[1] ) {
    t = radius[2]; 
    radius[2] = radius[1];
    radius[1] = t;
  }



  /* on tire un vecteur au hasard 
   */
  do {
    r[0] = 2.0 * (double)random() / (double)2147483647.0 - 1.0;
    r[1] = 2.0 * (double)random() / (double)2147483647.0 - 1.0;
    r[2] = 2.0 * (double)random() / (double)2147483647.0 - 1.0;
    t = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
  } while( t > 1 || t < 1e-6 );
  r[0] /= t;
  r[1] /= t;
  r[2] /= t;


  /* on tire un angle de rotation au hasard
     c'est en degre
   */
  *angle = 360.0 * ((double)random() / (double)2147483647.0);
  t = (*angle) * (3.1415927 / 180.0);
  r[0] *= t;
  r[1] *= t;
  r[2] *= t;



  
  /* la matrice de rotation d'angle t
   */
  _RotationMatrixFromRotationVector( mat, r );
  



  /* la translation qui conserve le centre  / M*C + t = C
     => t = C - M*C
   */
  r[0] = cx;
  r[1] = cy;
  r[2] = cz;
  _MultiplyMatrixByVector( mat, r, r );
  tx = cx - r[0];
  ty = cy - r[1];
  tz = cz - r[2];
  



  for ( iz=0; iz<theIm->dim.z; iz++ )
  for ( iy=0; iy<theIm->dim.y; iy++ )
  for ( ix=0; ix<theIm->dim.x; ix++ ) {
    /* coordonnees par rapport au centre 
     */
    r[0] = (double)ix;
    r[1] = (double)iy;
    r[2] = (double)iz;
    _MultiplyMatrixByVector( mat, r, r );
    x = r[0] + tx - cx;
    y = r[1] + ty - cy;
    z = r[2] + tz - cz;
    if ( (x*x)/(radius[0]*radius[0]) +
	 (y*y)/(radius[1]*radius[1]) +
	 (z*z)/(radius[2]*radius[2]) <= 1.0 )
      theBuf[iz][iy][ix] = color;
    else 
      theBuf[iz][iy][ix] = 0;
  }

}
















static void _MultiplyMatrixByVector( const double *m,
				     const double *a,
				     double *x )
{
  double t[3];
  
  t[ 0 ] = m[ 0 ] * a[ 0 ] + m[ 1 ] * a[ 1 ] + m[ 2 ] * a[ 2 ];
  t[ 1 ] = m[ 3 ] * a[ 0 ] + m[ 4 ] * a[ 1 ] + m[ 5 ] * a[ 2 ];
  t[ 2 ] = m[ 6 ] * a[ 0 ] + m[ 7 ] * a[ 1 ] + m[ 8 ] * a[ 2 ];
  x[ 0 ] = t[ 0 ];
  x[ 1 ] = t[ 1 ];
  x[ 2 ] = t[ 2 ];
}






static void _RotationMatrixFromRotationVector( double *mat,
					       double *rot )
{
  double f, g, theta, t2;
  
  t2 = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
  theta = sqrt( t2 );
  
  if ( theta > 1e-8 ) {
    f = sin( theta ) / theta;
    g = ( 1.0 - cos( theta ) ) / ( t2 );
    
    mat[0] = 1.0 - g * (rot[1]*rot[1] + rot[2]*rot[2]);
    mat[4] = 1.0 - g * (rot[2]*rot[2] + rot[0]*rot[0]);
    mat[8] = 1.0 - g * (rot[0]*rot[0] + rot[1]*rot[1]);
    
    mat[3] = mat[1] = g * rot[0] * rot[1];
    mat[6] = mat[2] = g * rot[0] * rot[2];
    mat[7] = mat[5] = g * rot[2] * rot[1];
    
    mat[1] -= f * rot[2];
    mat[2] += f * rot[1];
    mat[5] -= f * rot[0];
    
    mat[3] += f * rot[2];
    mat[6] -= f * rot[1];
    mat[7] += f * rot[0];
  }
  else {
    mat[0] = mat[4] = mat[8] = 1.0;
    mat[1] = mat[2] = mat[3] = 0.0;
    mat[5] = mat[6] = mat[7] = 0.0;
  }
}











int _MaxValueInImage( vt_image *theIm,
		      int slice )
{
  char *proc = "_MaxValueInImage";
  int v, offset=0;
  int volume = theIm->dim.x * theIm->dim.y * theIm->dim.z;
  int max=0, s = slice;

  
  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= theIm->dim.z) ) s = -1;
  /* cas 2D
   */
  if ( theIm->dim.z == 1 || s >= 0 ) {
    offset = s * theIm->dim.x * theIm->dim.y;
    volume = theIm->dim.x * theIm->dim.y;
  }

  switch ( theIm->type ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theIm->buf;
      max = theBuf[offset];
      for ( v=1; v<volume; v++ )
	if ( theBuf[offset+v] > max ) max = theBuf[offset+v];
    }
    break;
  case USHORT :
    {
      u16 *theBuf = (u16*)theIm->buf;
      max = theBuf[offset];
      for ( v=1; v<volume; v++ )
	if ( theBuf[offset+v] > max ) max = theBuf[offset+v];
    }
    break;
  default :
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) 
      fprintf( stderr, "%s: such image type is not handled yet.\n", proc );
    return( 0 );
  }
  
  return( max );
}

















void _InitArrayOfParametersFromImage( vt_image *theIm,
				     typeComponentParameters *thePar,
				     int n )
{
  int i;
  for ( i=0; i<=n; i++ ) {
    thePar[i].ptmin[0] = theIm->dim.x-1;
    thePar[i].ptmin[1] = theIm->dim.y-1;
    thePar[i].ptmin[2] = theIm->dim.z-1;
    thePar[i].ptmax[0] = 0;
    thePar[i].ptmax[1] = 0;
    thePar[i].ptmax[2] = 0;
    thePar[i].volume = 0;
  }
}










void _FillArrayOfParametersFromImage( vt_image *theIm,
				      int slice,
				      typeComponentParameters *thePar )
{
  char *proc = "_FillArrayOfParametersFromImage";
  int s=slice;
  int x, y, z, color;
  
  /* 2D ou 3D ?
   */
  if ( theIm->dim.z == 1 ) s = 0;
  if ( (theIm->dim.z > 1) && (s < 0 || s >= theIm->dim.z) ) s = -1;

  switch ( theIm->type ) {
  case USHORT :
    {
      u16 *** theBuf = (u16 ***)theIm->array;

      /* cas 2D
       */
      if ( theIm->dim.z == 1 || s >= 0 ) {
	for ( y=0; y<theIm->dim.y; y++ ) 
	for ( x=0; x<theIm->dim.x; x++ ) {

	  if ( theBuf[s][y][x] == 0 ) continue;

	  color = theBuf[s][y][x];

	  thePar[ color ].volume ++;
	  if ( x < thePar[ color ].ptmin[0] ) thePar[ color ].ptmin[0] = x;
	  if ( y < thePar[ color ].ptmin[1] ) thePar[ color ].ptmin[1] = y;
	  if ( x > thePar[ color ].ptmax[0] ) thePar[ color ].ptmax[0] = x;
	  if ( y > thePar[ color ].ptmax[1] ) thePar[ color ].ptmax[1] = y;
	  thePar[ color ].ptmin[2] = thePar[ color ].ptmax[2] = s;
	}
	break;
      }

      /* cas 3D 
       */
      for ( z=0; z<theIm->dim.z; z++ ) 
      for ( y=0; y<theIm->dim.y; y++ ) 
      for ( x=0; x<theIm->dim.x; x++ ) {

	if ( theBuf[z][y][x] == 0 ) continue;
	
	color = theBuf[z][y][x];
	
	thePar[ color ].volume ++;
	if ( x < thePar[ color ].ptmin[0] ) thePar[ color ].ptmin[0] = x;
	if ( y < thePar[ color ].ptmin[1] ) thePar[ color ].ptmin[1] = y;
	if ( z < thePar[ color ].ptmin[2] ) thePar[ color ].ptmin[2] = z;
	if ( x > thePar[ color ].ptmax[0] ) thePar[ color ].ptmax[0] = x;
	if ( y > thePar[ color ].ptmax[1] ) thePar[ color ].ptmax[1] = y;
	if ( z > thePar[ color ].ptmax[2] ) thePar[ color ].ptmax[2] = z;
	
      }

    }
    break;
  default :
    if ( _VT_VERBOSE_ || _VT_DEBUG_ ) 
      fprintf( stderr, "%s: such image type is not handled yet.\n", proc );
    return;
  }

}












void _TestsMaxDiameterInRandomList( int nbpoints,
				    int nbtests,
				    char *filename,
				    int testall )
{
  FILE *f, *fopen();

  typePoint *listOfPoints = (typePoint *)NULL;
  int t, i;
  double realmax, max;
  typePoint first, last;

  double clock1, clock2, duree;

  double totalTempsRecherche = 0.0;
  double minTempsRecherche = 0.0;
  double maxTempsRecherche = 0.0;

  double totalTempsMethodeMixte = 0.0;
  double minTempsMethodeMixte = 0.0;
  double maxTempsMethodeMixte = 0.0;

  if ( nbpoints <= 0 || nbtests <= 0 ) return;

  
  listOfPoints = (typePoint *)malloc( nbpoints * sizeof( typePoint ) );
  if ( listOfPoints == (typePoint *)NULL ) {
    fprintf( stderr, "unable to allocate points list\n" );
    return;
  }


  if ( filename[0] != '\0' && filename[0] != '>' && filename[0] != '<' )
    f = fopen( filename, "w" );
  else 
    f = stdout;


  (void)srandom(time(0));


  if ( testall != 0 ) {

    for ( t=0; t<nbtests; t++ ) {
      /* on tire au hasard des points de coordonnees entre 0 et 1 
       */
      for ( i=0; i<nbpoints; i++ ) {
	listOfPoints[i].x = (double)random() / (double)2147483647.0;
	listOfPoints[i].y = (double)random() / (double)2147483647.0;
	listOfPoints[i].z = (double)random() / (double)2147483647.0;
      }
      
      fprintf( f, "test #%4d:", t );
      
      clock1 = (double)clock();
      realmax = _ComputeMaxSquareDistanceInList( listOfPoints, nbpoints, &first, &last );
      clock2 = (double)clock();
      duree = clock2 - clock1;
      
      fprintf( f, " max = %9.7f", realmax );
      
      if ( t == 0 ) {
	minTempsRecherche = maxTempsRecherche = totalTempsRecherche = duree;
      } else {
	totalTempsRecherche += duree;
	if ( duree < minTempsRecherche )      minTempsRecherche = duree;
	else if ( duree > maxTempsRecherche ) maxTempsRecherche = duree;
      }
      
      clock1 = (double)clock();
      max = _EstimateMaxSquareDistanceInList( listOfPoints, nbpoints, &first, &last );
      clock2 = (double)clock();
      duree = clock2 - clock1;
      
      fprintf( f, "  - methode 2 = %9.7f (err=%7.4f%%)", max, 100.0*fabs( max - realmax ) / realmax );
      
      if ( t == 0 ) {
	minTempsMethodeMixte = maxTempsMethodeMixte = totalTempsMethodeMixte = duree;
      } else {
	totalTempsMethodeMixte += duree;
	if ( duree < minTempsMethodeMixte )      minTempsMethodeMixte = duree;
	else if ( duree > maxTempsMethodeMixte ) maxTempsMethodeMixte = duree;
    }
      
      fprintf( f, "\n" );
      
    }
    
    free( listOfPoints );
    
    fprintf( f, "\n\n" );
    
    fprintf( f, "... %d tests avec %d points\n", nbtests, nbpoints );
    fprintf( f, " TEMPS (secs)              MIN        MAX        AVERAGE\n" );
    fprintf( f, " recherche exhaustive  %9f  %9f  %9f\n",
	     minTempsRecherche/(double)CLOCKS_PER_SEC,
	     maxTempsRecherche/(double)CLOCKS_PER_SEC,
	     totalTempsRecherche/(nbtests * (double)CLOCKS_PER_SEC) );
    fprintf( f , " methode mixte         %9f  %9f  %9f\n",
	     minTempsMethodeMixte/(double)CLOCKS_PER_SEC,
	     maxTempsMethodeMixte/(double)CLOCKS_PER_SEC,
	     totalTempsMethodeMixte/(nbtests * (double)CLOCKS_PER_SEC) );
    
    if ( filename[0] != '\0' && filename[0] != '>' && filename[0] != '<' )
      fclose( f );
    
    return;
  }



  for ( t=0; t<nbtests; t++ ) {
    /* on tire au hasard des points de coordonnees entre 0 et 1 
     */
    for ( i=0; i<nbpoints; i++ ) {
      listOfPoints[i].x = (double)random() / (double)2147483647.0;
      listOfPoints[i].y = (double)random() / (double)2147483647.0;
      listOfPoints[i].z = (double)random() / (double)2147483647.0;
    }
    
    fprintf( f, "test #%4d:", t );
    
    clock1 = (double)clock();
    max = _EstimateMaxSquareDistanceInList( listOfPoints, nbpoints, &first, &last );
    clock2 = (double)clock();
    duree = clock2 - clock1;
      
    fprintf( f, "  methode 2 = %9.7f", max );
      
    if ( t == 0 ) {
      minTempsMethodeMixte = maxTempsMethodeMixte = totalTempsMethodeMixte = duree;
    } else {
      totalTempsMethodeMixte += duree;
      if ( duree < minTempsMethodeMixte )      minTempsMethodeMixte = duree;
      else if ( duree > maxTempsMethodeMixte ) maxTempsMethodeMixte = duree;
    }
      
    fprintf( f, "\n" );
      
  }
    
  free( listOfPoints );
  
  fprintf( f, "\n\n" );
    
  fprintf( f, "... %d tests avec %d points\n", nbtests, nbpoints );
  fprintf( f, " TEMPS (secs)              MIN        MAX      AVERAGE\n" );
  fprintf( f , " methode mixte         %9f  %9f  %9f\n",
	   minTempsMethodeMixte/(double)CLOCKS_PER_SEC,
	   maxTempsMethodeMixte/(double)CLOCKS_PER_SEC,
	   totalTempsMethodeMixte/(nbtests * (double)CLOCKS_PER_SEC) );
  
  if ( filename[0] != '\0' && filename[0] != '>' && filename[0] != '<' )
    fclose( f );
}
