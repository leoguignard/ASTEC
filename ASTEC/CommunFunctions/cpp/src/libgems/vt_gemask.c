
#include <vt_gemask.h>

#ifndef NO_PROTO
int _VT_GEMASK( vt_image *im, vt_vpt_amincir *liste, int lnb, int type_distance )
#else
int _VT_GEMASK( im, liste, lnb, type_distance )
vt_image *im;
vt_vpt_amincir *liste;
int lnb;
int type_distance;
#endif
{
    register int i, j, k, l;
    int x, y, z, v, d;
    int max_dist;
    /* masque */
    vt_distance dpar;
    vt_image masque;
    u16 ***pm;
    u8 ***pb;
    
    /* dans l'image il doit y avoir :
       VT_HASBEENDELETED : points deja effaces
       VT_UNDELETABLE    : points a conserver (seuil haut)
       VT_U_DELETABLE    : points a conserver (seuil bas mais topologie)
       dans la liste il doit y avoir :
       VT_U_DELETABLE : points a conserver (seuil bas mais topologie) (en debut de liste)
       VT_TOBEDELETED : points effaces
       */

    /*--- tests ---*/
    if ( VT_Test1Image( im, "_VT_GEMASK" ) == -1 ) return( -1 );
    if ( (im->type != UCHAR) ) return( -1 );
    if ( (type_distance != VT_DIST_CHMFR2) ) return( -1 );
    if ( lnb <= 0 ) return( 0 );

    /*--- maximum de distance ---*/
    max_dist = 0;
    for ( l = 0; l < lnb; l++ ) 
	if ( liste[l].pt.v > max_dist ) max_dist = liste[l].pt.v;
    max_dist = (int)( (float)max_dist / 16.0 + 0.5 ) + 1;
    
    /*--- calcul du masque ---*/
    VT_Image( &masque );
    VT_InitImage( &masque, "", max_dist+2, max_dist+2, max_dist+2, (int)USHORT );

    if ( VT_AllocImage( &masque ) != 1 ) {
	return( -1 );
    }
    pm = (u16 ***)masque.array;
    for ( k = 0; k < masque.dim.z; k ++ )
    for ( j = 0; j < masque.dim.y; j ++ )
    for ( i = 0; i < masque.dim.x; i ++ )
	pm[k][j][i] = (u16)0;
    pm[0][0][0] = (u16)255;
    VT_Distance( &dpar );
    dpar.type = VT_DIST_CHMFR2;
    dpar.seuil = 100.0; /* les valeurs sont >= 200 */
    if ( VT_Dist( &masque, &masque, &dpar ) != 1 ) {
	VT_FreeImage( &masque );
	return( -1 );
    }
    
    pb = (u8 ***)im->array;

    /*--- parcours de la liste ---*/
    for ( l = 0; l < lnb; l++ ) {
	x = liste[l].pt.x;
	y = liste[l].pt.y;
	z = liste[l].pt.z;
	v = liste[l].pt.v;
	d = (int)( (float)v / 16.0 + 0.5 );

	if ( (x - d >= 0) && (x + d < im->dim.x) &&
	     (y - d >= 0) && (y + d < im->dim.y) &&
	     (z - d >= 0) && (z + d < im->dim.z) ) {

	    for ( k = 0; k < masque.dim.z; k ++ )
	    for ( j = 0; j < masque.dim.y; j ++ )
	    for ( i = 0; i < masque.dim.x; i ++ ) {
		if ( pm[k][j][i] <= v ) {
		    if ( pb[z-k][y-j][x-i] == VT_HASBEENDELETED ) pb[z-k][y-j][x-i] = VT_RECONSTRUCTABLE;
		    if ( pb[z-k][y-j][x+i] == VT_HASBEENDELETED ) pb[z-k][y-j][x+i] = VT_RECONSTRUCTABLE;
		    if ( pb[z-k][y+j][x-i] == VT_HASBEENDELETED ) pb[z-k][y+j][x-i] = VT_RECONSTRUCTABLE;
		    if ( pb[z-k][y+j][x+i] == VT_HASBEENDELETED ) pb[z-k][y+j][x+i] = VT_RECONSTRUCTABLE;
		    if ( pb[z+k][y-j][x-i] == VT_HASBEENDELETED ) pb[z+k][y-j][x-i] = VT_RECONSTRUCTABLE;
		    if ( pb[z+k][y-j][x+i] == VT_HASBEENDELETED ) pb[z+k][y-j][x+i] = VT_RECONSTRUCTABLE;
		    if ( pb[z+k][y+j][x-i] == VT_HASBEENDELETED ) pb[z+k][y+j][x-i] = VT_RECONSTRUCTABLE;
		    if ( pb[z+k][y+j][x+i] == VT_HASBEENDELETED ) pb[z+k][y+j][x+i] = VT_RECONSTRUCTABLE;
		}
	    }
	    
	} else {

	    for ( k = 0; k < masque.dim.z; k ++ )
	    for ( j = 0; j < masque.dim.y; j ++ )
	    for ( i = 0; i < masque.dim.x; i ++ ) {
		if ( pm[k][j][i] > v ) continue;
		if ( (z-k >= 0) && (z-k < im->dim.z) ) {
		    if ( (y-j >= 0) && (y-j < im->dim.y) ) {
			if ( (x-i >= 0) && (x-i < im->dim.x) )
			    if ( pb[z-k][y-j][x-i] == VT_HASBEENDELETED ) pb[z-k][y-j][x-i] = VT_RECONSTRUCTABLE;
			if ( (x+i >= 0) && (x+i < im->dim.x) )
			    if ( pb[z-k][y-j][x+i] == VT_HASBEENDELETED ) pb[z-k][y-j][x+i] = VT_RECONSTRUCTABLE;
		    }
		    if ( (y+j >= 0) && (y+j < im->dim.y) ) {
			if ( (x-i >= 0) && (x-i < im->dim.x) )
			    if ( pb[z-k][y+j][x-i] == VT_HASBEENDELETED ) pb[z-k][y+j][x-i] = VT_RECONSTRUCTABLE;
			if ( (x+i >= 0) && (x+i < im->dim.x) )
			    if ( pb[z-k][y+j][x+i] == VT_HASBEENDELETED ) pb[z-k][y+j][x+i] = VT_RECONSTRUCTABLE;
		    }
		}
		if ( (z+k >= 0) && (z+k < im->dim.z) ) {
		    if ( (y-j >= 0) && (y-j < im->dim.y) ) {
			if ( (x-i >= 0) && (x-i < im->dim.x) )
			    if ( pb[z+k][y-j][x-i] == VT_HASBEENDELETED ) pb[z+k][y-j][x-i] = VT_RECONSTRUCTABLE;
			if ( (x+i >= 0) && (x+i < im->dim.x) )
			    if ( pb[z+k][y-j][x+i] == VT_HASBEENDELETED ) pb[z+k][y-j][x+i] = VT_RECONSTRUCTABLE;
		    }
		    if ( (y+j >= 0) && (y+j < im->dim.y) ) {
			if ( (x-i >= 0) && (x-i < im->dim.x) )
			    if ( pb[z+k][y+j][x-i] == VT_HASBEENDELETED ) pb[z+k][y+j][x-i] = VT_RECONSTRUCTABLE;
			if ( (x+i >= 0) && (x+i < im->dim.x) )
			    if ( pb[z+k][y+j][x+i] == VT_HASBEENDELETED ) pb[z+k][y+j][x+i] = VT_RECONSTRUCTABLE;
		    }
		}
	    }
	    
	}
    }

    /* dans l'image il doit y avoir maintenant :
       VT_HASBEENDELETED  : points deja effaces
       VT_RECONSTRUCTABLE : points pouvant etre reconstruits
       VT_UNDELETABLE     : points a conserver (seuil haut)
       VT_U_DELETABLE     : points a conserver (seuil bas mais topologie)
       dans la liste il doit y avoir maintenant :
       VT_U_DELETABLE : points a conserver (seuil bas mais topologie) (en debut de liste)
       VT_TOBEDELETED : points effaces
       */

    VT_FreeImage( &masque );
    return( 1 );
}


