/*************************************************************************
 * tracking-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 25 jui 2013 16:46:25 CEST
 *
 * ADDITIONS, CHANGES
 *
 */


#ifndef _tracking_tools_h_
#define _tracking_tools_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <histogram.h>


typedef enum {
  _UNKNOWN_TYPE_POINT_,
  _READ_,
  _ADDED_
} enumTypePoint;


typedef struct {
  int *data;
  int n_data;
  int center;
  int median;
  int min;
  int max;
  float mean;
} valueType;


struct circleType {
  /* cercle information :
     centre
     rayon
     index image
  */
  int x;
  int y;
  int r;
  int image;

  /* image intensities
   */
  valueType intensity;

  /* index de la liste de cercles
     contenant ce cercle
  */
  int index;
  
  /* type of the point
   */
  enumTypePoint type;

  /* neighbors in previous and following 
     circle list
     (if trying to build chains)
  */
  struct circleType **prev;
  int n_prev;
  int n_prev_allocated;

  struct circleType **next;
  int n_next;
  int n_next_allocated;

  /* number of chains the circle belongs to
     (if trying to build chains)
  */
  int n_chains;  

  /* not used 
   */
  int n_colocs;
};

typedef struct circleType circleType;



typedef struct {
  circleType *data;
  int n_data;
  int n_selected_data;
  int n_allocated_data;
  int rmax;
} circleList;



typedef struct {
  circleList *data;
  int n_data;
  int n_allocated_data;
  int rmax;
} circleListList;



typedef struct {
  circleType **data;

  int xmin;
  int xmax;
  int ymin;
  int ymax;

  int n_data;
  int n_allocated_data;
  int rmax;
} chainType;



typedef struct {
  chainType *data;
  int n_data;
  int n_selected_data;
  int n_allocated_data;
  int rmax;
} chainList;



typedef struct {
  int n_data;
  int *data;
  int dim;
  int *x;
  int *y;
} plotType;

typedef struct {
  plotType ncircle;
  plotType startingIndex;
  plotType endingIndex;
  plotType length;
} statisticType;



typedef enum {
  _CLOSEST_TO_PREVIOUS_,
  _CLOSEST_TO_REF11_
} enumChainBuilding;



extern void initCircleType( circleType *c );
extern void freeCircleType( circleType *c );

extern void initCircleList( circleList *l );
extern void freeCircleList( circleList *l );
extern void printCircleList( FILE *f, circleList *l );
extern int readCircleList( circleList *l, char *filename );

extern void initCircleListList( circleListList *l );
extern int allocCircleListList( circleListList *l, int n );
extern void freeCircleListList( circleListList *l );
extern int readCircleListList( circleListList *readlist, char *format, int first, int last );
extern void printTrackedCircleList( FILE *f, circleListList *l );
extern int readTrackedCircleList( circleListList *readlist, char *filename );
extern int fillTrackedCircleList( circleListList *filledlist, circleListList *readlist );
extern int intensityTrackedCircleList( circleListList *filledlist, char *format );

extern void initChainType( chainType *l );
extern void freeChainType( chainType *l );
extern void initChainList( chainList *l );
extern void freeChainList( chainList *l );




extern void initStatisticType( statisticType *s );
extern void freeStatisticType( statisticType *s );

extern void statsChainList( statisticType *s, circleListList *readlist, chainList *chainlist );


extern void printPlotList( FILE *f, plotType *p );
extern void printPlotHist( FILE *f, plotType *p );


extern void printChainResultXxxlab( statisticType *stats, 
			       char *desc,
			       char *name, 
			       enumHistogramFile xxxlab );
extern void printColocalizationResultXxxlab( statisticType *stats1, 
				      char *desc1,
				      statisticType *stats2, 
				      char *desc2,
					     statisticType *stats3, 
				      char *desc3,
				      char *name, 
				      enumHistogramFile xxxlab );



extern void removeStartBoundChain( chainList *list, int index_first_circle );
extern void removeEndBoundChain( chainList *list, int index_last_circle );

extern int chainListFromCircleListList( circleListList *readlist, chainList *chainlist, 
		       int depth, int margin,
				 int maxForwardNeighbors, int maxBackwardNeighbors );
extern void printStatsChainList( FILE *f, circleListList *readlist, chainList *chainlist, char *s );

extern int colocalizationFromLists( circleListList *colocalizationlist, 
			     chainList *chainlist, 
			     circleListList *readlist1, circleListList *readlist2, 
			     int depth, int margin,
			     int maxForwardNeighbors, int maxBackwardNeighbors );

extern int colocalizationFromChains( chainList *coloclist, 
			      chainList *chainlist1, chainList *chainlist2, 
				     int margin );

#ifdef __cplusplus
}
#endif

#endif 
