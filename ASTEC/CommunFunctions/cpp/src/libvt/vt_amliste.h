#ifndef _vt_amliste_h_
#define _vt_amliste_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>

#ifndef NO_PROTO
extern vt_pt_amincir*  _VT_ThinPtList( vt_image *image, 
				       int dim, 
				       int *nb );
extern vt_vpt_amincir* _VT_ThinVPtList( vt_image *image, 
					vt_image *value, 
					int dim, 
					int *nb );
#else
extern vt_pt_amincir*  _VT_ThinPtList();
extern vt_vpt_amincir* _VT_ThinVPtList();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_amliste_h_ */
