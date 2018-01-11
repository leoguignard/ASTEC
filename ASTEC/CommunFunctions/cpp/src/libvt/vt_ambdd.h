#ifndef _vt_ambdd_h_
#define _vt_ambdd_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>
#include <vt_bdd.h>

#ifndef NO_PROTO
extern int _VT_BDD_THINNING( vt_image *im, vt_pt_amincir *l, int *n, vt_amincir *p );
#else
extern int _VT_BDD_THINNING();
#endif /* NO_PROTO */



#ifdef __cplusplus
}
#endif

#endif /* _vt_ambdd_h_ */
