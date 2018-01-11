#ifndef _vt_gelrec3D_h_ 
#define _vt_gelrec3D_h_ 

#ifdef __cplusplus
extern "C" {
#endif



#include <vt_common.h>

#include <vt_geline.h>

#ifndef NO_PROTO 
extern int VT_Reconstruct3D( vt_image *ext, vt_image *scl, vt_line *par );
extern int VT_ScaleReconstruct3D( vt_resline *res, vt_line *par );
extern int VT_CsteReconstruct3D( vt_resline *res, double r );
extern int VT_GreyReconstruct3D( vt_resline *res );
#else
extern int VT_Reconstruct3D();
extern int VT_ScaleReconstruct3D();
extern int VT_CsteReconstruct3D();
extern int VT_GreyReconstruct3D();
#endif /* NO_PROTO */




#ifdef __cplusplus
}
#endif

#endif /* _vt_gelrec3D_h_ */ 
