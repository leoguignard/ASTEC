#ifndef _vt_image_quatreD_h_
#define _vt_image_quatreD_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_image.h>
#include <vt_inrimage.h>
#include <stdio.h>
#include <vt_unix.h>


typedef struct vt_name4D {
  int nbAllocatedNames;
  int nbNames;
  char **names;
} vt_name4D;


extern void VT_Name4D( vt_name4D *name );
extern int VT_AllocateName4D( vt_name4D *name );
extern void VT_FreeName4D( vt_name4D *name );

extern int VT_ReadNameInr4D( FILE *f,
		  vt_name4D *name4D );

extern int VT_ReadName4D( FILE *f,
		  vt_name4D *name );

extern int VT_CreateName4D( vt_name4D *name4D,
		     char *base,
		     char *suffixe,
		     int n );


extern void VT_WriteNameInr4D( vt_name4D *name4D, char *nom );

extern void VT_WriteName4D( vt_name4D *name4D, char *nom );





typedef struct vt_image4D {
  vt_image *images;
  void ****array;
  int dimt;
  vt_ipt dim;
  ImageType type;
} vt_image4D;


    
extern int VT_AllocAndInitImage4D( vt_image4D *image,
		     vt_name4D *name4D,
		     int dimx,
		     int dimy,
		     int dimz,
		     int dimt,
		     int type );

extern int VT_AllocAndReadImage4D( vt_image4D *image,
		     vt_name4D *name4D );

extern void VT_FreeImage4D( vt_image4D *image );



#include <vt_copy.h>

extern int VT_CopyImage4D( vt_image4D *theIm,
		    vt_image4D *resIm );


extern int VT_WriteImage4D( vt_image4D *theIm );



#ifdef __cplusplus
}
#endif

#endif /* _vt_image_quatreD_h_ */
