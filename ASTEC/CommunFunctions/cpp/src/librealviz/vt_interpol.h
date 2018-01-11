
#ifndef _VT_INTERPOL_H_
#define _VT_INTERPOL_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <vt_common.h>
#include <rl_interpol.h>


extern int VT_InitDistanceImageFromSlice( vt_image *theIm,
					 vt_image *theDist,
					 int index );


extern void VT_ComputeSliceFromDistances( vt_image *thePrev,
					vt_image *theNext,
					vt_image *theRes,
					int index,
				 double prevCoef,
				 double nextCoef );

extern int *VT_BuildTranslationTable( vt_image *theIm, int *min, int *max, int *nb );

extern void VT_ApplyTranslationTable( vt_image *theIm, int *t );

#ifdef __cplusplus
}
#endif

#endif
