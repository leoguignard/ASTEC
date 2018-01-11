/*************************************************************************
 * is_distance.h - calcul de distance pour le points/surface matching
 *
 * $Id: is_distance.h,v 1.5 1999/09/22 17:26:26 greg Exp $
 * 
  *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sat May 22 11:06:11 MET DST 1999
 *
 * Copyright Gregoire Malandain
 *
 *
 * ADDITIONS, CHANGES:
 *
 * - Sat Jun  5 10:05:18 MET DST 1999 (Gregoire Malandain)
 *   Ajout de la structure typeDistanceMap
 *   Modification des prototypes en consequence
 *
 *
 */

#ifndef _IS_DISTANCE_H_
#define _IS_DISTANCE_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>


extern void _DistanceSetVerbose ();
extern void _DistanceSetNoVerbose ();


extern void _DistanceSetComputationTo2D ();
extern void _DistanceSetComputationTo3D ();

typedef int v333[3][3][3];
typedef int v555[5][5][5];
typedef int v777[7][7][7];





/* structure des coeffients pour le calcul de distance
   
   Les valeurs sont entieres, 
   apres le calcul des distances par propagation,
   la distance estimee en millimetres s'obtiendra
   en multipliant la valeur obtenue par
   'multiplicativeCoefficient'.
*/
typedef struct {

  double multiplicativeCoefficient;

  /* masque 7x7x7 pour l'initialisation
     on calcule les distances a la surface 
     du voxel central (3,3,3) 
  */
  v777 init; 

  /* plus petite valeur d'initialisation
     a la surface du voisinage 7x7x7
  */
  int minInit;
  /* plus grande valeur d'initialisation
     a la surface du voisinage 7x7x7 
     => on sait que les voisins des points
     de valeur superieur a cette plus grande 
     valeur sont du meme cote de l'interface,
     ce n'est pas la peine de tester leur signe
     (ie leur appartenance au fond ou a l'objet)
  */
  int maxInit;

  /* coefficients pour une propagation
   */
  v555 incr;
   
} typeDistanceCoefficients;




  
typedef struct {
  short int *buf;
  int dim[3];
  float voxelSize[3];
  double multiplicativeCoefficient;
  int intensityMax;
} typeDistanceMap;



/* Procedure d'evaluation de la distance
   en un point reel
*/
extern int _InterpoleSignedDistance( double *d, 
				     const double x, const double y, 
				     const double z, 
				     const typeDistanceMap *theDist );


/* Procedures de calcul d'une carte de distance
 */

extern void _ComputeSignedDistanceMap( typeDistanceMap *theDist,
				       const double maxDistUpdate );


extern void _InitBorderSignedDistanceMap( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par );


extern void _ComputeSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
						       const typeDistanceCoefficients *par );
  

extern void _UpdateSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate );



/* Procedures d'initialisation
 */

extern void _InitSignedDistanceMap( const unsigned char *theBuf,
				    short int *resBuf,
				    const int *theBufferDimensions );




extern void _InitDistanceCoefficients( typeDistanceCoefficients *par,
				       float *voxel );


extern void _PrintDistanceCoefficients( typeDistanceCoefficients *par );


extern void _CombineTwoDistanceMaps2D( typeDistanceMap *theDist1,
				       typeDistanceMap *theDist2 );



#ifdef __cplusplus
}
#endif

#endif
