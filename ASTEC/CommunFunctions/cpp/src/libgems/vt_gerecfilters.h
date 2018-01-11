#ifndef _vt_gerecfilters_h_
#define _vt_gerecfilters_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_recfilters.h>
extern int GE_RecFilterOnImage( vt_image *theIm, vt_image *resIm, vt_recfilters *par );

#include <vt_contours.h>
extern int GE_MaximaGradient( vt_image *theIm, vt_image *resIm, vt_contours *par );


#ifdef __cplusplus
}
#endif

#endif /* _vt_gerecfilters_h_ */
