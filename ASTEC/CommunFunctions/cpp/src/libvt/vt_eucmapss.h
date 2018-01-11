#ifndef _vt_eucmapss_h_
#define _vt_eucmapss_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_distance.h>

#ifndef NO_PROTO
extern int  VT_VecteurPPP_SS( vt_image *I, vt_image *X, vt_image *Y, vt_image *Z, vt_distance *p );
#else /* NO_PROTO */
extern int  VT_VecteurPPP_SS();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif
