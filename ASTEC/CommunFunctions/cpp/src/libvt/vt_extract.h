#ifndef _vt_extract_h_
#define _vt_extract_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>

#ifndef NO_PROTO
extern int VT_Extract( vt_image *im1, vt_image *im2, vt_ipt *corner );
#else
extern int VT_Extract();
#endif

#ifdef __cplusplus
}
#endif

#endif /* _vt_extract_h_ */
