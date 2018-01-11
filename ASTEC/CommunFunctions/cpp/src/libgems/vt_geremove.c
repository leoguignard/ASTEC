
#include <vt_geremove.h>

#ifndef NO_PROTO
int _VT_GEREMOVE( vt_image *im, vt_vpt_amincir *liste, int *lnb )
#else
int _VT_GEREMOVE( im, liste, lnb )
vt_image *im;
vt_vpt_amincir *liste;
int *lnb;
#endif
{
    vt_vpt_amincir point;
    register int i, j, k, l, changes=0, bool=0;
    int n, nb, v[3][3][3];
    register u8 ***pb;
    
    /* dans l'image il doit y avoir :
       VT_HASBEENDELETED : points deja effaces
       VT_UNDELETABLE    : points a conserver (seuil haut)
       VT_DELETABLE      : points effacables
       dans la liste il doit y avoir :
       VT_TOBEDELETED : points effaces
       VT_DELETABLE   : points effacables (en debut de liste)
       */

    /*--- tests ---*/
    if ( VT_Test1Image( im, "_VT_GEREMOVE" ) == -1 ) return( -1 );
    if ( (im->type != UCHAR) ) return( -1 );
    if ( *lnb <= 0 ) return( 0 );

    pb = (unsigned char***)(im->array);

    n = nb = *lnb;
    do {
	changes = 0;
	/*--- parcours de la liste ---*/
	for ( l = 0; l < n; l ++ ) {

	    /*--- traitement d'un point : saisie du voisinage ---*/
	    if ( liste[l].inside == 1 ) {
		for ( k = -1; k < 2; k++ )
		for ( j = -1; j < 2; j++ )
		for ( i = -1; i < 2; i++ )
		    v[1+i][1+j][1+k] = (int)(pb[liste[l].pt.z + k][liste[l].pt.y + j][liste[l].pt.x + i]);
	    } else {
		for ( k = -1; k < 2; k++ )
		for ( j = -1; j < 2; j++ )
		for ( i = -1; i < 2; i++ ) {
		    if ( (liste[l].pt.x+i >= 0) && (liste[l].pt.x+i < im->dim.x) &&
			 (liste[l].pt.y+j >= 0) && (liste[l].pt.y+j < im->dim.y) &&
			 (liste[l].pt.z+k >= 0) && (liste[l].pt.z+k < im->dim.z) )
			v[1+i][1+j][1+k] = (int)(pb[liste[l].pt.z + k][liste[l].pt.y + j][liste[l].pt.x + i]);
		    else
			v[1+i][1+j][1+k] = (int)0;
		}
	    }

	    /*--- changement de status ---*/
	    bool = 0;
	    for ( k = 0; (k < 3) && (bool == 0); k++ )
	    for ( j = 0; (j < 3) && (bool == 0); j++ )
            for ( i = 0; (i < 3) && (bool == 0); i++ ) {
		if ( (v[i][j][k] == VT_UNDELETABLE) || (v[i][j][k] == VT_U_DELETABLE) ) {
		    bool = 1;
		    liste[l].status = VT_U_DELETABLE;
		    pb[ liste[l].pt.z ][ liste[l].pt.y ][ liste[l].pt.x ] = (u8)VT_U_DELETABLE;
		}
	    }
	}
	
	/*--- on rearrange la liste ---*/
	i = 0;
	for ( l = 0; l < n; l++ ) {
	    switch ( liste[l].status ) {
		case VT_DELETABLE :
		    /*--- on garde le point dans la liste :
		          on swappe les points pour tout
		          garder dans la liste              ---*/
		    point    = liste[i];
		    liste[i] = liste[l];
		    liste[l] = point;
		    i++;
		    break;
	    case VT_U_DELETABLE :
		/*--- on ne fait rien ---*/
		changes ++;
		break;
	    default :
		/*--- on s'inquiete ---*/
		VT_Error( "point with unexpected label", "_VT_GEREMOVE" );
	    }
	}

	/*-- nouveau nombre de points ---*/
	n = i;

    } while ( (n > 0) && (changes > 0) );

    /*--- on peut effacer tous les points effacable qui restent :
          ce sont les n premiers de la liste                      ---*/
    for ( l = 0; l < n; l++ ) {
	pb[ liste[l].pt.z ][ liste[l].pt.y ][ liste[l].pt.x ] = (u8)VT_HASBEENDELETED;
	liste[l].status = VT_TOBEDELETED;
    }
    
    /*--- on regroupe les VT_U_DELETABLE en debut de liste ---*/
    i = 0;
    for ( l = n; l < nb; l++ ) {
	if ( liste[l].status == VT_U_DELETABLE ) {
	    point    = liste[i];
	    liste[i] = liste[l];
	    liste[l] = point;
	    i++;
	}
    }
	
    *lnb = i;

    /* dans l'image il doit y avoir maintenant :
       VT_HASBEENDELETED : points deja effaces
       VT_UNDELETABLE    : points a conserver (seuil haut)
       VT_U_DELETABLE    : points a conserver (seuil bas mais topologie)
       dans la liste il doit y avoir maintenant :
       VT_U_DELETABLE : points a conserver (seuil bas mais topologie) (en debut de liste)
       VT_TOBEDELETED : points effaces
       */

    return( 1 );
}





