
#include <vt_geseuil1.h>

#include <epidaure.h>
#include <epidaure.ee>
#include <epidaure.ei>

#ifndef NO_PROTO
#else 
#endif

#ifndef NO_PROTO
void VT_Seuil1( vt_seuil1 *par )
#else
void VT_Seuil1( par )
vt_seuil1 *par;
#endif
{
    /*--- erosion de l'image initiale --*/
    par->ero_connexity = N26;
    par->ero_iterations = 1;
    /*--- detection de contours ---*/
    VT_Contours( &(par->par_cont) );
    par->par_cont.length_continue.x = par->par_cont.length_continue.y = par->par_cont.length_continue.z = 15;
    par->par_cont.value_coefficient.x = par->par_cont.value_coefficient.y = par->par_cont.value_coefficient.z = 1.5;
    par->par_cont.type_filter = VT_RECFILTERS_DERICHE;
    par->par_cont.type_contours = VT_MAXIMA_GRADIENT;
    /*--- seuillage ---*/
    par->mul_max_sh = 0.33333333333333333333; /* 1/3 */
    par->mul_max_sb = 0.16666666666666666666; /* 1/6 */
    /*--- dilatation de l'image seuillee ---*/
    par->dil_connexity = N06;
    par->dil_iterations = 1;
    /*--- lissage de l'histogramme ---*/
    VT_RecFilters( &(par->par_filt) );
    par->par_filt.type_filter = VT_RECFILTERS_DERICHE;
    par->par_filt.derivative.x = VT_SMOOTHING;
    par->par_filt.length_continue.x = 15;
    par->par_filt.value_coefficient.x = 0.75;
}

#ifndef NO_PROTO
int VT_Image2Histo( vt_histo *histo, vt_image *image, vt_seuil1 *par )
#else
int VT_Image2Histo( histo, image, par )
vt_histo *histo;
vt_image *image;
vt_seuil1 *par; 
#endif
{
    vt_image imMorpho, imExt;
    vt_3m m;
    char *proc="VT_Image2Histo";
    char *name_imMorpho="VT_Image2Histo.morpho";
    char *name_imExt="VT_Image2Histo.ext";
    char *name_imThres = "VT_Image2Histo.thres";
    char *name_imLogic = "VT_Image2Histo.logic";
    char *name_imHyst = "VT_Image2Histo.hyst";
    char *name_imMorpho2 ="VT_Image2Histo.morpho2";
    char *name_imLogic2 = "VT_Image2Histo.logic2";

    /*--- erosion de l'image initiale ---*/
    VT_Image( &imMorpho );
    VT_InitImage( &imMorpho, name_imMorpho, image->dim.x, image->dim.y, image->dim.z, image->type );
    if ( VT_AllocImage( &imMorpho ) != 1 ) {
	VT_Error( "unable to allocate morphology image", proc );
	return( -1 );
    }
    /*
    if ( VT_ErosionWithConnectivity( image, &imMorpho, par->ero_connexity, par->ero_iterations ) != 1 ) {
    */
    if ( GreyLevelErosion( image, &imMorpho, par->ero_connexity, par->ero_iterations ) != 1 ) {
	VT_Error( "unable to erode input image", proc );
	VT_FreeImage( &imMorpho );
	return( -1 );
    }
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_WriteInrimage( &imMorpho );
    }

    /*--- seuillage de l'image erodee ---*/
    if ( VT_Threshold( &imMorpho, &imMorpho, (float)1.0 ) != 1 ) {
	VT_Error( "unable to threshold eroded image", proc );
	VT_FreeImage( &imMorpho );
	return( -1 );
    }
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_CopyName( imMorpho.name, name_imThres );
      (void)VT_WriteInrimage( &imMorpho );
    }

    /*--- calcul des contours ---*/
    VT_Image( &imExt );
    VT_InitImage( &imExt, name_imExt, image->dim.x, image->dim.y, image->dim.z, image->type );
    if ( VT_AllocImage( &imExt ) != 1 ) {
	VT_Error( "unable to allocate contours image", proc );
	VT_FreeImage( &imMorpho );
	return( -1 );
    }
    if ( VT_ExtractEdges( image, &imExt, &(par->par_cont) ) != 1 ) {
	VT_Error( "unable to compute contours", proc );
	VT_FreeImage( &imMorpho );
	VT_FreeImage( &imExt );
	return( -1 );
    }
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_WriteInrimage( &imExt );
    }

    /*--- on enleve les contours proches des 0 ---*/
    VT_LogicEt( &imExt, &imMorpho, &imExt );
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_CopyName( imExt.name, name_imLogic );
      (void)VT_WriteInrimage( &imExt );
    }

    /*--- calcul des min, moy et max des extrema restants ---*/
    if ( VT_3m( &imExt, &m ) != 1 ) {
	VT_Error( "unable to compute maximum value of gradient extrema", proc );
	VT_FreeImage( &imMorpho );
	VT_FreeImage( &imExt );
	return( -1 );
    }
    if ( _VT_DEBUG_ == 1 ) {
      char message[256];
      sprintf( message,"extrema: min = %f moy = %f max = %f", m.min, m.moy, m.max );
      VT_Message( message, proc );
    }

    /*--- seuillage des contours :
          on passe dans EpidaureLib ---*/
    {
	e_ParamConnex cpar;
	t_Image theIm, resIm;

	E_InitConnexe( &cpar );
	cpar.lowThres = (float)( m.max * par->mul_max_sb );
	cpar.higThres = (float)( m.max * par->mul_max_sh );
	if ( _VT_VERBOSE_ == 1 ) {
	    char message[256];
	    sprintf( message," hysteresis : seuil bas = %f, seuil haut = %f", cpar.lowThres, cpar.higThres );
	    VT_Message( message, proc );
	}
	switch ( image->type ) {
	case UCHAR :
	    cpar.lowThres /= (float)( 255.0 );
	    cpar.higThres /= (float)( 255.0 );
	    T_Init3Image( &theIm, image->dim.x, image->dim.y, image->dim.z, TYPEIM_256, image->name );
	    T_Init3Image( &resIm, image->dim.x, image->dim.y, image->dim.z, TYPEIM_256, image->name );
	    break;
	case USHORT :
	    cpar.lowThres /= (float)( 65535.0 );
	    cpar.higThres /= (float)( 65535.0 );
	    T_Init3Image( &theIm, image->dim.x, image->dim.y, image->dim.z, TYPEIM_16B, image->name );
	    T_Init3Image( &resIm, image->dim.x, image->dim.y, image->dim.z, TYPEIM_16B, image->name );
	    break;
	default :
	    VT_Error( "unable to deal with such image type", proc );
	    VT_FreeImage( &imMorpho );
	    VT_FreeImage( &imExt );
	    return( -1 );
	}

	theIm.buf = (char*)(imExt.buf);
	resIm.buf = (char*)(imMorpho.buf);

	if ( E_Hysteresis( &theIm, &resIm, &cpar ) != TRUE ) {
	    VT_Error( "unable to threshold gradient extrema", proc );
	    VT_FreeImage( &imMorpho );
	    VT_FreeImage( &imExt );
	    return( -1 );
	}	
    }
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_CopyName( imMorpho.name, name_imHyst );
      (void)VT_WriteInrimage( &imMorpho );
    }
    
    /*--- dilatation de l'image des contours seuilles ---*/
    /*
    if ( VT_DilationWithConnectivity( &imMorpho, &imMorpho, par->dil_connexity, par->dil_iterations ) != 1 ) {
    */
    if ( BinaryDilation( &imMorpho, &imMorpho, par->dil_connexity, par->dil_iterations ) != 1 ) {
	VT_Error( "unable to dilate contours image", proc );
	VT_FreeImage( &imMorpho );
	VT_FreeImage( &imExt );
	return( -1 );
    }
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_CopyName( imMorpho.name, name_imMorpho2 );
      (void)VT_WriteInrimage( &imMorpho );
    }
    
    /*--- on recupere les valeurs sur les contours ---*/
    VT_LogicEt( &imMorpho, image, &imExt );
    VT_FreeImage( &imMorpho );
    if ( _VT_DEBUG_ == 1 ) {
      (void)VT_CopyName( imExt.name, name_imLogic2 );
      (void)VT_WriteInrimage( &imExt );
    }

    /*--- on calcule l'histogramme ---*/
    if ( VT_ComputeHisto( histo, &imExt ) != 1 ) {
	VT_Error( "unable to fill histogram from image", proc );
	VT_FreeImage( &imExt );
	return( -1 );
    }

    /*--- liberations ---*/
    VT_FreeImage( &imExt );

    return( 1 );
}



