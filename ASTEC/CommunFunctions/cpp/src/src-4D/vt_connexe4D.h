#ifndef _vt_connexe_quatreD_h_
#define _vt_connexe_quatreD_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_image4D.h>


extern int VT_Hysteresis4D( vt_image4D *theIm, 
		     vt_image4D *resIm,
		     float seuilBas,
		     float seuilHaut );



#ifdef __cplusplus
}
#endif

#endif /* _vt_connexe_quatreD_h_ */
