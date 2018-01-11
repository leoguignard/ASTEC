#ifndef _vt_isolated_h_
#define _vt_isolated_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_neighborhood.h>

#ifndef NO_PROTO
extern int VT_DeleteIsolatedPoints( vt_image *resIm, vt_image *theIm, int connexite );
#else
extern int VT_DeleteIsolatedPoints();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_isolated_h_ */
