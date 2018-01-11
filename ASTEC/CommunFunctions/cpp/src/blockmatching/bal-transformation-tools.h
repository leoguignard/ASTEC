/*************************************************************************
 * bal-transformation-tools.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_TRANSFORMATION_TOOLS_H
#define BAL_TRANSFORMATION_TOOLS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <string-tools.h>

#include <bal-image.h>
#include <bal-transformation.h>
#include <bal-field.h>
#include <bal-estimator.h>



typedef enum {
  NEAREST,
  LINEAR
} enumTransformationInterpolation;




extern void BAL_SetVerboseInBalTransformationTools( int v );
extern void BAL_IncrementVerboseInBalTransformationTools(  );
extern void BAL_DecrementVerboseInBalTransformationTools(  );




/***************************************************
 *
 * transformation from one image to an other
 *
 ***************************************************/

extern int BAL_ComputeImageToImageTransformation( bal_image *subsampled_image,
						  bal_image *image_to_be_subsampled,
						  bal_transformation *subsampling_trsf );



/***************************************************
 *
 * transformation from a pairing field
 *
 ***************************************************/

extern int BAL_ComputeIncrementalTransformation( bal_transformation *T,  
						 FIELD *field, 
						 bal_estimator *estimator );

extern int BAL_ComputeTransformationResiduals( bal_transformation *T,  
					       FIELD *field );

/***************************************************
 *
 * transformation composition
 * t_res = t1 o t2 
   t_{I3<-I1} = t_{I3<-I2} o t_{I2<-I1}
 *
 ***************************************************/

extern int BAL_TransformationComposition( bal_transformation *t_res,
					  bal_transformation *t1,  
					  bal_transformation *t2 );

extern enumTypeTransfo BAL_TypeTransformationComposition( bal_transformation *t1, 
							  bal_transformation *t2 );

extern int BAL_AllocTransformationComposition( bal_transformation *res, 
					       bal_transformation *t1, 
					       bal_transformation *t2,
					       bal_image *ref );



/***************************************************
 *
 * transformation use
 *
 ***************************************************/

/* point transformation
 */
extern int BAL_TransformPoint( bal_doublePoint *thePt, bal_doublePoint *resPt, bal_transformation *theTr );

/* image resampling
 */
extern int BAL_ResampleImage( bal_image *image, bal_image *resim, bal_transformation *theTr,
			      enumTransformationInterpolation interpolation );



/***************************************************
 *
 * transformation conversion
 *
 ***************************************************/
/* unit conversion

   'theTrsf' is the transformation that goes from 'resim' to 'image' 
   ie allows to resample 'image' into the geometry of 'resim'.
   When 'theTrsf' is the transformation issued from matching,
   'resim' is the reference image and 'image' the floating image.
   
*/

extern int BAL_ChangeTransformationToRealUnit( bal_image *image,
					       bal_image *resim, 
					       bal_transformation *theTrsf,
					       bal_transformation *resTrsf );

extern int BAL_ChangeTransformationToVoxelUnit( bal_image *image,
						bal_image *resim, 
						bal_transformation *theTrsf,
						bal_transformation *resTrsf );



/***************************************************
 *
 * transformation inversion
 *
 ***************************************************/

extern int BAL_InverseTransformation( bal_transformation *theTrsf,  
				      bal_transformation *invTrsf );



/***************************************************
 *
 * transformation construction
 *
 ***************************************************/

extern int BAL_Sinusoid3DVectorField( bal_transformation *theTrsf );
extern int BAL_Sinusoid2DVectorField( bal_transformation *theTrsf );

/* random transformation 
 */
extern int BAL_Random2DTranslationMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DTranslationMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DTranslationScalingMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DTranslationScalingMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DRigidMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DRigidMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DSimilitudeMatrix( bal_transformation *theTrsf ) ;
extern int BAL_Random3DSimilitudeMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DAffineMatrix( bal_transformation *theTrsf );
extern int BAL_Random3DAffineMatrix( bal_transformation *theTrsf );
extern int BAL_Random2DVectorField( bal_transformation *theTrsf );
extern int BAL_Random3DVectorField( bal_transformation *theTrsf );
extern int BAL_RandomTransformation( bal_transformation *t, enumTypeTransfo type, bal_image *ref );

#ifdef __cplusplus
}
#endif

#endif
