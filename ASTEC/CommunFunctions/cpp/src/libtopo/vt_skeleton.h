/*************************************************************************
 * vt_skeleton.h -
 *
 * $Id: vt_skeleton.h,v 1.1 2000/07/26 07:25:02 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue Jul 25 21:39:22 MEST 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_skeleton_h_
#define _vt_skeleton_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>




typedef enum {
  _ANGLE_,
  _DISTANCE_,
  _PRODUIT_,
  _PRODUIT_LOG_,
  _LOG_PRODUIT_,
  _PRODUIT_2_,
  _PRODUIT_3_,
  _PRODUIT_4_,
  _PRODUIT_5_,
  _PRODUIT_6_,
  _PRODUIT_7_
} enumCoefficient;

/* ZYX,
   B =-1
   H = 0
   F = 1
*/
typedef enum {
  _BBB_ =  1, /* -1 -1 -1 */
  _BBH_ =  2,
  _BBF_ =  3,
  _BHB_ =  4,
  _BHH_ =  5,
  _BHF_ =  6,
  _BFB_ =  7,
  _BFH_ =  8,
  _BFF_ =  9,
  _HBB_ = 10, /*  0 -1 -1 */
  _HBH_ = 11,
  _HBF_ = 12,
  _HHB_ = 13,
  _HHH_ =  0,
  _HHF_ = 14,
  _HFB_ = 15,
  _HFH_ = 16,
  _HFF_ = 17,
  _FBB_ = 18, /*  1 -1 -1 */
  _FBH_ = 19,
  _FBF_ = 20,
  _FHB_ = 21,
  _FHH_ = 22,
  _FHF_ = 23,
  _FFB_ = 24,
  _FFH_ = 25,
  _FFF_ = 26
} enumNeighbor;



extern int VT_Compute2DSkeletonCoefficient( vt_image *imx,
					    vt_image *imy,
					    vt_image *imc,
					    vt_image *imn,
					    enumCoefficient typeCoefficient );





#ifdef __cplusplus
}
#endif

#endif /* _vt_skeleton_h_ */

