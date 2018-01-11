/*************************************************************************
 * bal-biolib-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 17 jan 2014 11:23:42 CET
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */

#ifndef BAL_BIOLIB_TOOLS_H
#define BAL_BIOLIB_TOOLS_H


#ifdef __cplusplus
extern "C" {
#endif



#include <bal-stddef.h>






/************************************************************
 *
 * Detection  structure
 *
 ************************************************************/

typedef struct {

  /* spot center (voxel coordinates)
   */
  bal_doublePoint voxelcenter;

  /*   VOXEL_UNIT,
       REAL_UNIT
  */
  enumUnitTransfo unit;

  /* spot radii (voxel coordinates)
   */
  int voxelhalfradius1;
  int voxelhalfradius2;

  double arg1;
  int    arg2;
  int    arg3;
  double arg4;
  double arg5;

  
  /* frameindex
   */
  int imageindex;

  /* pour retrouver le point si on le met dans une liste de pointeurs
   */
  int identifier;

  /* for lemon stuff
   */
  int idnode;

} bal_blDetection;

extern void BAL_InitBlDetection( bal_blDetection *d );





/************************************************************
 *
 * Pointer detection list
 *
 ************************************************************/

typedef struct {

  bal_blDetection **data;
  int n;          
  int n_allocated;
  
  /*   VOXEL_UNIT,
       REAL_UNIT
  */
  enumUnitTransfo unit;

} bal_blPointerDetectionList;

extern void BAL_InitBlPointerDetectionList( bal_blPointerDetectionList *l );
extern int BAL_AllocBlPointerDetectionList( bal_blPointerDetectionList *l, int n );
extern void BAL_FreeBlPointerDetectionList( bal_blPointerDetectionList *l );





/************************************************************
 *
 * search dedicated structure
 *
 ************************************************************/

typedef struct {

  bal_blPointerDetectionList *buf;
  bal_blPointerDetectionList ***array;

  bal_doublePoint lcorner;
  bal_doublePoint rcorner;

  double bucketsizex;
  double bucketsizey;
  double bucketsizez;

  int dimx;
  int dimy;
  int dimz;

} bal_blNeighborsSearch;





/************************************************************
 *
 * Detection list
 *
 ************************************************************/

typedef struct {

  bal_blDetection *data;
  int n;          
  int n_allocated;

  /*   VOXEL_UNIT,
       REAL_UNIT
  */
  enumUnitTransfo unit;


  /* bounding box of the detection list
     this only takes into account the centers
   */
  bal_doublePoint lcorner;
  bal_doublePoint rcorner;

  /* search structure
   */
  bal_blNeighborsSearch buckets;

} bal_blDetectionList;





extern void BAL_InitBlDetectionList( bal_blDetectionList *l );
extern int BAL_AllocBlDetectionList( bal_blDetectionList *l, int n );
extern void BAL_FreeBlDetectionList( bal_blDetectionList *l );
extern int BAL_InitNeighborsSearchBlDetectionList(  bal_blDetectionList *l, bal_doublePoint* bucketsize );

extern int BAL_ReadBlDetectionList( bal_blDetectionList *l, char *name );
extern int BAL_WriteBlDetectionList( bal_blDetectionList *l, char *name );

extern void BAL_BlDetectionListBoundingBox( bal_blDetectionList *t );








/************************************************************
 *
 *
 *
 ************************************************************/


typedef enum {
  _VOXEL_BRUTEFORCE_,
  _VOXEL_BUCKETS_
} enumSearchType;


extern int BAL_SearchBlDetectionNeighbors( bal_blDetection *d,
					   bal_blDetectionList *l,
					   bal_blPointerDetectionList *neighbors,
					   bal_doublePoint* halfneighborhood,
					   enumSearchType searchType );






/************************************************************
 *
 *
 *
 ************************************************************/

typedef struct {
  /* indexes of the first and the last images 
   */
  int index_start;
  int index_end;

  double time_start;
  double time_end;

  double lifetime;

  /* the euclidean distance from the first to the last point
     the total length (sum of elementary length)
   */
  double fly_distance_2D;
  double length_2D;
  double fly_distance_3D;
  double length_3D;
} bal_blTrackProperties;



typedef struct {
  bal_blDetectionList detectionList;
  /* track number when reading file
   */
  int index;

  /* bounding box of the track
     this includes the half radius of the detected spots (only in X and Y)
   */
  bal_doublePoint lcorner;
  bal_doublePoint rcorner;

  /* some properties
   */
  bal_blTrackProperties properties;

} bal_blTrack;



typedef struct {
  bal_blTrack *data;
  int n;          
  int n_allocated;
} bal_blTrackList;



typedef struct {
  int min_index_start;
  int max_index_end;
  int min_index_range;
  
  double min_fly_distance_2D;
  double min_length_2D;
} bal_blTrackSelectionCriteria;


extern void BAL_InitBlTrackSelectionCriteria( bal_blTrackSelectionCriteria *c );


extern void BAL_ComputeTrackListProperties( bal_blTrackList *l,
					    double vx, 
					    double vy,
					    double vz,
					    double dt );

extern int BAL_PrintAllStatsScilab( char *filename, bal_blTrackList *list );


extern void BAL_InitBlTrack( bal_blTrack *t );
extern void BAL_FreeBlTrack( bal_blTrack *t );
extern int BAL_CopyBlTrack( bal_blTrack *thet,  bal_blTrack *rest );

extern void BAL_InitBlTrackList( bal_blTrackList *l );
extern int BAL_AllocBlTrackList( bal_blTrackList *l, int n );
extern void BAL_FreeBlTrackList( bal_blTrackList *l );

extern int BAL_AddBlTrackToList( bal_blTrack *d, bal_blTrackList *l );
extern int BAL_CopyBlTrackToList( bal_blTrack *d, bal_blTrackList *l );

extern int BAL_ReadBlTrackList( bal_blTrackList *l, char *name );
extern int BAL_WriteBlTrackList( bal_blTrackList *l, char *name );



extern void BAL_BlTrackBoundingBox( bal_blTrack *t );





#ifdef __cplusplus
}
#endif

#endif
