#ifndef _vt_gemask_h_
#define _vt_gemask_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>
#include <vt_distance.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern int _VT_GEMASK( vt_image *im, vt_vpt_amincir *liste, int lnb, int type );
#else
extern int _VT_GEMASK();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_gemask_h_ */
