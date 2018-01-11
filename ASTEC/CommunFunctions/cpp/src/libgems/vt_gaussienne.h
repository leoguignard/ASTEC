#ifndef _vt_gaussienne_h_
#define _vt_gaussienne_h_


#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>

#ifndef NO_PROTO
extern int VT_FitGaussienne( r32 *tab, int l, double *K, double *sigma, double *moy );
#else
extern int VT_FitGaussienne();
#endif /* NO_PROTO */



#ifdef __cplusplus
}
#endif

#endif /* _vt_gaussienne_h_ */
