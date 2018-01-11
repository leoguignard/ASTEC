/*************************************************************************
 * bal-biolib-tracker.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Lun 24 mar 2014 20:37:39 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */


#ifndef BAL_BIOLIB_TRACKER_H
#define BAL_BIOLIB_TRACKER_H

#ifdef __cplusplus
extern "C" {
#endif




#include <bal-biolib-tools.h>
#include <bal-image.h>
#include <bal-transformation.h>


extern void BAL_SetVerboseInBalBiolibTracker( int v );
extern void BAL_IncrementVerboseInBalBiolibTracker(  );
extern void BAL_DecrementVerboseInBalBiolibTracker(  );


typedef struct {

  int index;

  /* les detections dans l'image originale
   */
  char *detectionName;
  bal_blDetectionList readDetectionList;
  bal_blDetectionList trsfDetectionList;

  bal_blDetectionList *detectionList;

  /* la transformation liee a l'image originale (en coordonnees reelles)
     permet de reechantillonner l'image originale dans
     l'image de reference donc va de la reference vers l'image originale

     il faut donc l'inverser pour l'appliquer aux detections

     apres lecture des images de vesselness, on multiplie ces transformations
     pour avoir des transformations de l'espace voxel des images originales
     ou sont les detections, a l'espace reel de la reference
   */
  char *transformationName;
  bal_transformation theTrsf;
  bal_transformation invTrsf;

  /* l'image de vesselness
   */
  char *vesselnessName;
  bal_image vesselnessImage;

} bal_blTrackerElement;





extern void BAL_InitBlTrackerElement( bal_blTrackerElement *e );
extern void BAL_FreeBlTrackerElement( bal_blTrackerElement *e );



typedef struct {
  bal_blTrackerElement *data;
  int n;          
  int n_allocated;
} bal_blTrackerElementList;

extern void BAL_InitBlTrackerElementList( bal_blTrackerElementList *l );
extern int BAL_AllocBlTrackerElementList( bal_blTrackerElementList *l, int n );
extern void BAL_FreeBlTrackerElementList( bal_blTrackerElementList *l );


extern int BAL_ReadDetectionsBlTrackerElementList( bal_blTrackerElementList *l );





/* cost = _FROM_NEXT_ corresponds to axonal growth (an axon appears)
   cost = _FROM_BOTH_ compares current and previous frames to evaluate
    retraction
*/

typedef enum {
  _FROM_PREV_,
  _FROM_NEXT_,
  _FROM_BOTH_
} enumCostLocation;

typedef enum {
  _EXACT_,
  _APPROXIMATE_
} enumCostCalculation;





typedef struct {
  /* search neighborhood
   */
  bal_doublePoint halfneighborhood_unit;
  bal_doublePoint halfneighborhood_voxel;

  /* minimal distance above which the cost is computed
     else it is set to the minimal value
  */
  double costMin2Ddistance_unit;
  double costMin2Ddistance_voxel;
  /* how to compute the cose
   */
  enumCostLocation costLocation;
  enumCostCalculation costCalculation;

  /* rejection parameter for paths
     - on the measure
     - on the length
     - on the number
  */
  double minVesselness;
  double maxVesselness;
  int minTrackLength;
  int maxTracks;

  enumSearchType searchType;

} bal_blTrackerParameter;

extern void BAL_InitBlTrackerParameter( bal_blTrackerParameter *p );





extern int BAL_ExtractBlTracks( bal_blTrackList *reslist,
			 bal_blTrackerElementList *thelist, 
			 bal_blTrackerParameter *par  );




#ifdef __cplusplus
}
#endif

#endif
