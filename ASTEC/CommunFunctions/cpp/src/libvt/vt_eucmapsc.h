#ifndef _vt_eucmapsc_h_
#define _vt_eucmapsc_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_distance.h>

#ifndef NO_PROTO
extern int  VT_VecteurPPP_SC( vt_image *I, vt_image *X, vt_image *Y, vt_image *Z, vt_distance *p );
#else /* NO_PROTO */
extern int  VT_VecteurPPP_SC();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif
