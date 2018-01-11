#ifndef _vt_zerocrossings_h_
#define _vt_zerocrossings_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_seuil.h>

#ifndef NO_PROTO
extern int  VT_ZeroCrossings( vt_image *r,  vt_image *t, DimType d );
extern int _VT_ZeroCrossings( vt_image *r,  vt_image *t, DimType d, int ze, int zo );
#else /* NO_PROTO */
extern int  VT_ZeroCrossings();
extern int _VT_ZeroCrossings();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif
