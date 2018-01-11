#ifndef _PREDEF_MORPHO_H_
#define _PREDEF_MORPHO_H_

#ifdef __cplusplus
extern "C" {
#endif

extern int _u8_predef_binary_Dilation( u8* inputBuf, /* buffer to be resampled */
	    u8* resultBuf, /* result buffer */
	    unsigned int dimx,
	    unsigned int dimy,
	    unsigned int dimz,
	    int connectivity, /* connectivity to be used */ 
	    int iterations  /* number of iterations */ );

extern int _u16_predef_binary_Dilation( u16* inputBuf, /* buffer to be resampled */
	    u16* resultBuf, /* result buffer */
	    unsigned int dimx,
	    unsigned int dimy,
	    unsigned int dimz,
	    int connectivity, /* connectivity to be used */ 
	    int iterations  /* number of iterations */ );

extern int _u8_predef_binary_Erosion( u8* inputBuf, /* buffer to be resampled */
	    u8* resultBuf, /* result buffer */
	    unsigned int dimx,
	    unsigned int dimy,
	    unsigned int dimz,
	    int connectivity, /* connectivity to be used */ 
	    int iterations  /* number of iterations */ );

extern int _u16_predef_binary_Erosion( u16* inputBuf, /* buffer to be resampled */
	    u16* resultBuf, /* result buffer */
	    unsigned int dimx,
	    unsigned int dimy,
	    unsigned int dimz,
	    int connectivity, /* connectivity to be used */ 
	    int iterations  /* number of iterations */ );

#ifdef __cplusplus
}
#endif

#endif
