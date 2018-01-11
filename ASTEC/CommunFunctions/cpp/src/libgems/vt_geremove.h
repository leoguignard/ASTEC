#ifndef _vt_geremove_h_
#define _vt_geremove_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern int _VT_GEREMOVE( vt_image *im, vt_vpt_amincir *liste, int *lnb );
#else
extern int _VT_GEREMOVE();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_geremove_h_ */
