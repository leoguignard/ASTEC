#ifndef _vt_gemedian_h_
#define _vt_gemedian_h_

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>
#include <vt_copy.h>

#ifndef NO_PROTO
extern int GE_MedianFilter( vt_image *theIm, vt_image *resIm, vt_ipt *window );
#else 
extern int GE_MedianFilter();
#endif

#ifdef __cplusplus
}
#endif

#endif /* _vt_gemedian_h_ */
