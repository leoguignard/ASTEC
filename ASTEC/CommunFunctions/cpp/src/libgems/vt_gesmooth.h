#ifndef _vt_gesmooth_h_
#define _vt_gesmooth_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern int VT_FastSmoothUC( vt_image *resIm, vt_image *theIm, int seuil, int dim );
#else
extern int VT_FastSmoothUC();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_gesmooth_h_ */
