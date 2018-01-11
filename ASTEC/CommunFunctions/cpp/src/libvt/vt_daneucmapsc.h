#ifndef _vt_daneucmapsc_h_
#define _vt_daneucmapsc_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_distance.h>

#ifndef NO_PROTO
extern int  VT_DanVecteurPPP_SC( vt_image *I, vt_image *X, vt_image *Y, vt_image *Z, vt_distance *p );
#else /* NO_PROTO */
extern int  VT_DanVecteurPPP_SC();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif
