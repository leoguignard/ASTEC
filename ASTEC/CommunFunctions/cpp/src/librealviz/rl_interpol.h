/*************************************************************************
 * rl_interpol.h - procedures for image interpolation
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * July, 6 1999
 *
 * ADDITIONS, CHANGES
 *
 *
 */




#ifndef _RL_INTERPOL_H_
#define _RL_INTERPOL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <typedefs.h>
#include <is_distance.h>


extern int _InitDistanceImageFromSlice( const void *theSlice,
				 bufferType typeSlice,
				  const int dimx,
				  const int dimy,
				  short int *theDist,
				  const int dimz );


extern void _ComputeSliceFromDistances( void *theSlice,
				 bufferType typeSlice,
				 const int dimx,
				 const int dimy,
				 const short int *prevDist,
				 const short int *nextDist,
				 int dimz,
				 const double prevCoef,
				 const double nextCoef );

extern  int _TestConversion( const unsigned char *theSlice,
		     const int dimx,
		     const int dimy );
#ifdef __cplusplus
}
#endif

#endif
