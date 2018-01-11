
#include <vt_isolated.h>

/*-------Definition des fonctions statiques----------*/
#ifndef NO_PROTO
#else
#endif

#ifndef NO_PROTO
int VT_DeleteIsolatedPoints( vt_image *resIm /* result image */,
			     vt_image *theIm /* input image */,
			     int connexite   /* connectivity */ )
#else
int VT_DeleteIsolatedPoints( resIm, theIm, connexite )
vt_image *resIm; /* result image */
vt_image *theIm; /* input image */
int connexite;   /* connectivity */
#endif
{
    int xo[26], yo[26], zo[26], nb, dim;
    register int n, x, y, z, p;
    int bool_Z, bool_YZ, bool_XYZ;
    int dx1, dy1, dz1;
    char *local_name = "VT_DeleteIsolatedPoints";
    u8 ***resBuf, ***theBuf;

    if ( VT_Test2Image( resIm, theIm, local_name ) == -1 ) return( -1 );
    if ( (resIm->type != UCHAR) || (theIm->type != UCHAR) ) {
	VT_Error( "incorrect image type", local_name );
	return( -1 );
    }
    
    /*--- nombre de voisins et dimension ---*/
    nb = 26;
    dim = VT_3D;
    switch ( connexite ) {
    case N06 :    
	nb = 6;
	dim = VT_3D;
	break;
    case N26 :
    default :
	nb = 26;
	dim = VT_3D;
	break;
    case N18 :
	nb = 18;
	dim = VT_3D;
	break;
    case N10 :
	nb = 10;
	dim = VT_3D;
	break;
    case N08 :
	nb = 8;
	dim = VT_2D;
	break;
    case N04 :
	nb = 4;
	dim = VT_2D;
    }

    for ( x = 0; x < 26; x ++ )
	xo[x] = yo[x] = zo[x] = 0;

    switch ( connexite ) {
    case N06 :
	xo[ 5] =  0;   yo[ 5] =  0;  zo[ 5] = -1;
	xo[ 4] =  0;   yo[ 4] = -1;  zo[ 4] =  0;
	xo[ 3] = -1;   yo[ 3] =  0;  zo[ 3] =  0;
	xo[ 2] =  1;   yo[ 2] =  0;  zo[ 2] =  0;
	xo[ 1] =  0;   yo[ 1] =  1;  zo[ 1] =  0;
	xo[ 0] =  0;   yo[ 0] =  0;  zo[ 0] =  1;
	break;
    case N26 :
	xo[25] = -1;   yo[25] = -1;  zo[25] = -1;
	xo[24] =  1;   yo[24] = -1;  zo[24] = -1;
	xo[23] = -1;   yo[23] =  1;  zo[23] = -1;
	xo[22] =  1;   yo[22] =  1;  zo[22] = -1;
	xo[21] = -1;   yo[21] = -1;  zo[21] =  1;
	xo[20] =  1;   yo[20] = -1;  zo[20] =  1;
	xo[19] = -1;   yo[19] =  1;  zo[19] =  1;
	xo[18] =  1;   yo[18] =  1;  zo[18] =  1;
    case N18 :
	xo[17] =  0;   yo[17] = -1;  zo[17] = -1;
	xo[16] = -1;   yo[16] =  0;  zo[16] = -1;
	xo[15] =  1;   yo[15] =  0;  zo[15] = -1;
	xo[14] =  0;   yo[14] =  1;  zo[14] = -1;
	xo[13] =  0;   yo[13] = -1;  zo[13] =  1;
	xo[12] = -1;   yo[12] =  0;  zo[12] =  1;
	xo[11] =  1;   yo[11] =  0;  zo[11] =  1;
	xo[10] =  0;   yo[10] =  1;  zo[10] =  1;
    case N10 :
	xo[ 9] =  0;   yo[ 9] =  0;  zo[ 9] = -1;
	xo[ 8] =  0;   yo[ 8] =  0;  zo[ 8] =  1;
    case N08 :
	xo[ 7] = -1;   yo[ 7] = -1;  zo[ 7] =  0;
	xo[ 6] =  1;   yo[ 6] = -1;  zo[ 6] =  0;
	xo[ 5] = -1;   yo[ 5] =  1;  zo[ 5] =  0;
	xo[ 4] =  1;   yo[ 4] =  1;  zo[ 4] =  0;
    case N04 :
	xo[ 3] =  0;   yo[ 3] = -1;  zo[ 3] =  0;
	xo[ 2] = -1;   yo[ 2] =  0;  zo[ 2] =  0;
	xo[ 1] =  1;   yo[ 1] =  0;  zo[ 1] =  0;
	xo[ 0] =  0;   yo[ 0] =  1;  zo[ 0] =  0;
    }

    /*--- traitement ---*/
    bool_XYZ = bool_YZ = bool_Z = 0;
    dx1 = theIm->dim.x - 1;
    dy1 = theIm->dim.y - 1;
    dz1 = theIm->dim.z - 1;

    theBuf = (u8 ***)theIm->array;
    resBuf = (u8 ***)resIm->array;

    for ( z = 0; z < theIm->dim.z; z ++ ) {
	if ( (dim == VT_3D) && ((z == 0) || (z == dz1)) ) bool_Z = 0;
	else                                              bool_Z = 1;
	for ( y = 0; y < theIm->dim.y; y ++ ) {
	    if ( (bool_Z == 0) || (y == 0) || (y == dy1) ) bool_YZ = 0;
	    else                                                  bool_YZ = 1;
	    for ( x = 0; x < theIm->dim.x; x ++ ) {
		if ( (bool_YZ == 0) || (x == 0) || (x == dx1) ) bool_XYZ = 0;
		else                                                   bool_XYZ = 1;

		/*--- on compte les voisins ---*/
		p = 0;
		if ( bool_XYZ == 1 ) {
		    for ( n = 0; n < nb; n ++ )
			if ( theBuf[z + zo[n]][y + yo[n]][x + xo[n]] > 0 ) p++;
		} else {
		    for ( n = 0; n < nb; n ++ ) {
			if ( (z + zo[n] >= 0) && (z + zo[n] < theIm->dim.z) &&
			     (y + yo[n] >= 0) && (y + yo[n] < theIm->dim.y) &&
			     (x + xo[n] >= 0) && (x + xo[n] < theIm->dim.x) )
			    if ( theBuf[z + zo[n]][y + yo[n]][x + xo[n]] > 0 ) p++;
		    }
		}
		
		/*--- nouvelle valeur du point ---*/
		if ( p == 0 ) resBuf[ z ][ y ][ x ] = (u8)0;
		else          resBuf[ z ][ y ][ x ] = theBuf[ z ][ y ][ x ];
	    }
	}
    }

    return( 1 );
}
