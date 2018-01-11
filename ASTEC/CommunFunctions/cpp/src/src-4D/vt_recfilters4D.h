#ifndef _vt_recfilters_quatreD_h_
#define _vt_recfilters_quatreD_h_

#ifdef __cplusplus
extern "C" {
#endif


#include <vt_image4D.h>
#include <vt_recfilters.h>
#include <vt_contours.h>


extern int VT_MaximaGradient4D( vt_image4D *theIm, 
			 vt_image4D *resIm, 
			 vt_contours *spacePar,
			 vt_contours *timePar );


extern int VT_RecFilterTOnImage4D( vt_image4D *theIm, 
			     vt_image4D *resIm, 
			     vt_recfilters *par );

extern int VT_RecFilter3DOnImage4D( vt_image4D *theIm, 
			     vt_image4D *resIm, 
			     vt_recfilters *par );


extern int VT_NormeGradient4DImage4D( vt_image4D *theIm,
			       vt_image4D *resIm,
			       vt_contours *timePar,
			       vt_contours *spacePar,
			       int derivative );


extern int VT_NormeGradient4DImage4DWithDerivatives( vt_image4D *theX, 
			       vt_image4D *theY,
			       vt_image4D *theZ,
			       vt_image4D *theT,
			       vt_image4D *theNorme );


extern int VT_NormeGradient3DImage4DWithDerivatives( vt_image4D *theX, 
			       vt_image4D *theY,
			       vt_image4D *theZ,
			       vt_image4D *theNorme );



extern void VT_Recfilters4DVerbose( );
extern void VT_Recfilters4DNoVerbose( );


#ifdef __cplusplus
}
#endif

#endif /* _vt_recfilters_quatreD_h_ */


