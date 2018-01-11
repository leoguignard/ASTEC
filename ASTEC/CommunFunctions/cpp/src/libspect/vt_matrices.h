/*************************************************************************
 * vt_matrices.h -
 *
 * $Id: vt_matrices.h,v 1.1 2000/03/23 08:46:47 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Mar 22 22:07:09 MET 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _vt_matrices_h_
#define _vt_matrices_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
  
extern int Read4x4Matrix( char *name, double *mat );
extern int Inverse4x4Matrix( double *matrice, double *inv );


#ifdef __cplusplus
}
#endif

#endif /* _vt_matrices_h_ */
