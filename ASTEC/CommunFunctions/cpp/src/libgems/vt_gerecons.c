#include <vt_common.h>
#include <vt_common.ei>
#include <vt_common.ee>

#include <vt_amincir.h>

#include <vt_distance.h>
#include <vt_eucmapsc.ee>

#if defined(_ANSI_)
int _VT_GERECONSTRUCT( vt_image *resIm )
#else
int _VT_GERECONSTRUCT( resIm )
vt_image *resIm;
#endif
{
    char *proc="_VT_GERECONSTRUCT"
    vt_image auxX, auxY, auxZ;
    s8  ***scx, ***scy, ***scz, ***pb;
    register int x, y, z, n;
    vt_distance dpar;

    /* dans l'image il doit y avoir maintenant :
       VT_HASBEENDELETED  : points deja effaces
       VT_RECONSTRUCTABLE : points pouvant etre reconstruits
       VT_UNDELETABLE     : points a conserver (seuil haut)
       VT_U_DELETABLE     : points a conserver (seuil bas mais topologie)
       */

    if ( VT_Test1Image( im, proc ) == -1 ) return( -1 );
    if ( ( resIm->dim.x < 3 ) || ( resIm->dim.y < 3 ) ) {
	VT_Error( "images have bad dimensions", "VT_EucliDist" );
	return( -1 );
    }
    scx = scy = scz = (s8***)NULL;

    /*--- allocations ---*/
    VT_InitImage( &auxX, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)SCHAR );
    VT_InitImage( &auxY, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)SCHAR );
    VT_InitImage( &auxZ, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)SCHAR );
    if ( VT_AllocImage( &auxX ) != 1 ) {
	VT_Error( "unable to allocate auxiliary image", proc );
	return( -1 );
    }
    if ( VT_AllocImage( &auxY ) != 1 ) {
	VT_Error( "unable to allocate auxiliary image", proc );
	VT_FreeImage( &auxX );
	return( -1 );
    }
    if ( VT_AllocImage( &auxZ ) != 1 ) {
	VT_Error( "unable to allocate auxiliary image", proc );
	VT_FreeImage( &auxX );
	VT_FreeImage( &auxY );
	return( -1 );
    }
    
    /*--- calcul ---*/
    VT_Distance( &dpar );
    dpar.seuil = 235.0;
    if ( VT_VecteurPPP_SC( theIm, &auxX, &auxY, &auxZ, par ) != 1 ) {
	VT_FreeImage( &auxX );
	VT_FreeImage( &auxY );
	VT_FreeImage( &auxZ );
	return( -1 );
    }

    /*---  ecriture de l'image de sortie ---*/
    scx = (s8***)(auxX.array);
    scy = (s8***)(auxY.array);
    scz = (s8***)(auxZ.array);
    pb = (s8***)(resIm->array);
