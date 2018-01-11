#ifndef _vt_geseuil1_h_
#define _vt_geseuil1_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_neighborhood.h>

#include <vt_contours.h>
#include <vt_recfilters.h>

typedef struct vt_seuil1 {
    /*--- erosion de l'image initiale --*/
    int ero_connexity;
    int ero_iterations;
    /*--- detection de contours ---*/
    vt_contours par_cont;
    /*--- seuillage ---*/
    float mul_max_sh;
    float mul_max_sb;
    /*--- dilatation de l'image seuillee ---*/
    Neighborhood dil_connexity;
    int dil_iterations;
    /*--- lissage de l'histogramme ---*/
    vt_recfilters par_filt;
} vt_seuil1;



#include <vt_morpho.h>

#include <vt_seuil.h>

#include <vt_histo.h>

#ifndef NO_PROTO
extern void VT_Seuil1( vt_seuil1 *par );
extern int  VT_Image2Histo ( vt_histo *histo, vt_image *image, vt_seuil1 *par );
#else
extern void VT_Seuil1();
extern int  VT_Image2Histo();
#endif /* NO_PROTO */



#ifdef __cplusplus
}
#endif

#endif /* _vt_geseuil1_h_ */
