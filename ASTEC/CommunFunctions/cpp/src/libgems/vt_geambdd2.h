#ifndef _vt_geambdd_bis_h_
#define _vt_geambdd_bis_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern int _VT_GEBDD2_THINNING( vt_image *im, vt_vpt_amincir *l, int *n, vt_amincir *p, int m );
#else
extern int _VT_GEBDD2_THINNING();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_geambdd_bis_h_ */
