
#include <vt_gesmooth.h>

#ifndef NO_PROTO
int VT_FastSmoothUC( vt_image *resIm, vt_image *theIm, int seuil, int dim )
#else
int VT_FastSmoothUC( resIm, theIm, seuil, dim )
vt_image *resIm;
vt_image *theIm;
int seuil;
int dim;
#endif
{
    register int x, y, z, new, old;
    int dx1, dy1, dz1, local_dim;
    register double mul;
    u8 ***theBuf, ***resBuf;

    if ( VT_Test2Image( resIm, theIm, "VT_FastSmoothUC" ) == -1 ) return( -1 );
    if ( (resIm->type != UCHAR) ) return( -1 );
    if ( (theIm->type != UCHAR) ) return( -1 );
    if ( (theIm->dim.x < 2) || (theIm->dim.y < 2) ) return( -1 );
    
    /*--- le masque est separable et est [1 2 1]
          au bord, on prend [3 1] ou [1 3]        ---*/

    local_dim = VT_3D;
    if ( (dim == VT_2D) || (theIm->dim.z < 2) ) local_dim = VT_2D;

    dx1 = theIm->dim.x - 1;
    dy1 = theIm->dim.y - 1;
    dz1 = theIm->dim.z - 1;

    theBuf = (u8 ***)theIm->array;
    resBuf = (u8 ***)resIm->array;

    /*--- seuillage ---*/
    for ( z = 0; z <= dz1; z ++ )
    for ( y = 0; y <= dy1; y ++ )
    for ( x = 0; x <= dx1; x ++ )
	resBuf[z][y][x] = ( theBuf[z][y][x] < seuil ) ? (u8)0 : (u8)1;
    
    /*--- traitement selon X ---*/
    for ( z = 0; z <= dz1; z ++ )
    for ( y = 0; y <= dy1; y ++ ) {
	/*--- premier point ---*/
	old = resBuf[z][y][0] + resBuf[z][y][0] + resBuf[z][y][0] + resBuf[z][y][1];
	/*--- points centraux ---*/
	for ( x = 1; x < dx1; x ++ ) {
	    new = resBuf[z][y][x-1] + resBuf[z][y][x] + resBuf[z][y][x] + resBuf[z][y][x+1];
	    resBuf[z][y][x-1] = old;
	    old = new;
	}
	/*--- dernier point ---*/
	new = resBuf[z][y][dx1-1] + resBuf[z][y][dx1] + resBuf[z][y][dx1] + resBuf[z][y][dx1];
	resBuf[z][y][dx1-1] = old;
	resBuf[z][y][dx1] = new;
    }

    /*--- traitement selon Y ---*/
    for ( z = 0; z <= dz1; z ++ )
    for ( x = 0; x <= dx1; x ++ ) {
	/*--- premier point ---*/
	old = resBuf[z][0][x] + resBuf[z][0][x] + resBuf[z][0][x] + resBuf[z][1][x];
	/*--- points centraux ---*/
	for ( y = 1; y < dy1; y ++ ) {
	    new = resBuf[z][y-1][x] + resBuf[z][y][x] + resBuf[z][y][x] + resBuf[z][y+1][x];
	    resBuf[z][y-1][x] = old;
	    old = new;
	}
	/*--- dernier point ---*/
	new = resBuf[z][dy1-1][x] + resBuf[z][dy1][x] + resBuf[z][dy1][x] + resBuf[z][dy1][x];
	resBuf[z][dy1-1][x] = old;
	resBuf[z][dy1][x] = new;
    }
  
    if ( local_dim == VT_3D ) {
	/*--- traitement selon Z ---*/
	for ( y = 0; y <= dy1; y ++ )
	for ( x = 0; x <= dx1; x ++ ) {
	    /*--- premier point ---*/
	    old = resBuf[0][y][x] + resBuf[0][y][x] + resBuf[0][y][x] + resBuf[1][y][x];
	    /*--- points centraux ---*/
	    for ( z = 1; z < dz1; z ++ ) {
		new = resBuf[z-1][y][x] + resBuf[z][y][x] + resBuf[z][y][x] + resBuf[z+1][y][x];
		resBuf[z-1][y][x] = old;
		old = new;
	    }
	    /*--- dernier point ---*/
	    new = resBuf[dz1-1][y][x] + resBuf[dz1][y][x] + resBuf[dz1][y][x] + resBuf[dz1][y][x];
	    resBuf[dz1-1][y][x] = old;
	    resBuf[dz1][y][x] = new;
	}
    }

    mul = (double)255.0;
    if ( local_dim == VT_3D ) mul /= 64.0;
    else                      mul /= 16.0;
    
    /*--- multiplication ---*/
    for ( z = 0; z <= dz1; z ++ )
    for ( y = 0; y <= dy1; y ++ )
    for ( x = 0; x <= dx1; x ++ )
	resBuf[z][y][x] = (int)( (double)resBuf[z][y][x] * mul + 0.5 );
    
    return( 1 );
}

 

