#ifndef _vt_slconnexe_h_
#define _vt_slconnexe_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_neighborhood.h>
#include <vt_connexe.h>


#ifndef NO_PROTO
extern int VT_SlowConnexe( vt_image *image, 
			   vt_connexe *par );
extern int VT_SlowConnexeWithThresholds( vt_image *i, 
					 vt_connexe *p, 
					 vt_ipt *s, 
					 vt_ipt *f );
#else
extern int VT_SlowConnexe();
extern int VT_SlowConnexeWithThresholds();
#endif

#ifdef __cplusplus
}
#endif

#endif
