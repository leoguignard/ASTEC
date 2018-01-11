/*************************************************************************
 * fit-distribution-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 19 jui 2013 22:49:08 CEST
 *
 * ADDITIONS, CHANGES
 *
 */


#ifndef _fit_distribution_tools_h_
#define _fit_distribution_tools_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <histogram.h>



typedef struct {
  int *data;
  size_t n_data;
  size_t n_allocated_data;
} measureList;


extern void initMeasureList( measureList *l );
extern void freeMeasureList( measureList *l );
extern int addMeasureToList( measureList *l, int m );
  extern int readMeasureList( measureList *l, char *filename );


extern int maxMeasureList( measureList *l );
extern int minMeasureList( measureList *l );

extern int build1DHistogramFromMeasureList( typeHistogram *h, measureList *l );

extern int fit3WeibullDistributionsOnCumulative( typeHistogram *h,
						 double *initValues,
						 char *filename,
						 enumHistogramFile xxxlab,
						 int minlength,
						 int maxlength );

#ifdef __cplusplus
}
#endif

#endif 
