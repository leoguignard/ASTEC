#ifndef _vt_daneucmapss_h_
#define _vt_daneucmapss_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_distance.h>

#ifndef NO_PROTO
extern int  VT_DanVecteurPPP_SS( vt_image *I, vt_image *X, vt_image *Y, vt_image *Z, vt_distance *p );
#else /* NO_PROTO */
extern int  VT_DanVecteurPPP_SS();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif
