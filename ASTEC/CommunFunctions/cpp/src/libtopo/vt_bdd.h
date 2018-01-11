/*************************************************************************
 * vt_bdd.h -
 *
 * $Id: vt_bdd.h,v 1.1 2001/03/21 18:58:29 greg Exp $
 *
 * Copyright©INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */


#ifndef _vt_bdd_h_
#define _vt_bdd_h_

#ifdef __cplusplus
extern "C" {
#endif



/* Calcul si un point est simple ou non
 * renvoie 1 s'il est simple
 *         0 sinon
 *
 * la numerotation du voisinage est la suivante
 *
 *  0  1  2     9 10 11    17 18 19  
 *  3  4  5    12    13    20 21 22
 *  6  7  8    14 15 16    23 24 25
 *
 * les points de l'objet doivent avoir une valeur non nulle,
 * les points du fond doivent avoir une valeur nulle
 *
 * on utilise la 26-connexite pour l'objet
 * on utilise la  6-connexite pour le fond
 *
 * RETURN:
 *  Retourne 1 si le point est simple
 *  0 sinon
 */
extern int VT_IsSimple( int *V /* tableau de 26 valeurs (indices de 0 a 25) */ );



#ifdef __cplusplus
}
#endif

#endif /* _vt_bdd_h_ */
