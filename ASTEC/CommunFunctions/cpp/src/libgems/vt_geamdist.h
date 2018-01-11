#ifndef _vt_geamdist_h_
#define _vt_geamdist_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern int _VT_GEDIST_THINNING( vt_image *im, vt_vpt_amincir *l, int *n, vt_amincir *p, int d );
#else
extern int _VT_GEDIST_THINNING();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_geamdist_h_ */
