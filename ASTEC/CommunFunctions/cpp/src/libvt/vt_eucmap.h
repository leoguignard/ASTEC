#ifndef _vt_eucmap_h_
#define _vt_eucmap_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>

#include <vt_eucmapsc.h>
#include <vt_eucmapss.h>

#include <vt_distance.h>

#ifndef NO_PROTO
extern int  VT_EucliDist( vt_image *resIm, vt_image *theIm, vt_distance *par );
#else /* NO_PROTO */
extern int  VT_EucliDist();
#endif /* NO_PROTO */


#ifdef __cplusplus
}
#endif

#endif
