#ifndef _VT_RECUSDD_H_
#define _VT_RECUSDD_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include <vt_common.h>

typedef struct {
  double thetaMin;
  double thetaMax;
  double radiusMin;
  double radiusMax;
} typeGeometryUS2D;



extern int _ReconstructUS2D( vt_image *theUs,
		       vt_image *theRec,
		       typeGeometryUS2D *theGeom );



#ifdef __cplusplus
}
#endif

#endif
