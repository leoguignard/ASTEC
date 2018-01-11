/*************************************************************************
 * vt_meristemeFormationAxes.h -
 *
 * $Id: vt_meristemeFormationAxes.h,v 1.0 2014/06/16 16:36:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 * 
 * CREATION DATE: 
 * 2014/06/16 
 *
 * ADDITIONS, CHANGES
 *
 *
 */
 

#ifndef _vt_meristemeFormationAxes_h_
#define _vt_meristemeFormationAxes_h_


#ifdef __cplusplus
extern "C" {
#endif

/* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))

/* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))

/* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)

extern int addLineToBuf(void *buf, int *dim, bufferType t, int *pt, double tht, double phi, int r);

extern int addConeToBuf(void ***array, int *dim, bufferType t, int *pt, double tht, double phi, double r, double alpha);

extern int extractLabelFromLine(void *buf, int *dim, bufferType t, int *pt, double tht, double phi, int r, int *label, int *p1, int *p2);

extern int extractLabelFromCone(void ***array, int *dim, bufferType t, int *pt, double tht, double phi, double r, double alpha, int *label, int *p1, int *p2);

extern int setLineToLabel(void *buf, int *dim, bufferType t, int *p1, int *p2, int label);

extern int setConeToLabel(void ***array, int *dim, bufferType t, int *p1, int *p2, double tht, double phi, double alpha, int label);

extern int VT_conicVote(vt_image **imagesIn, vt_image *imres, double r, double alpha, int NangleIter );

extern void VT_ComputeConicField3D(unsigned char ***field, double *angles, int *dimField, double r, double alpha);

extern int VT_nearestAngle(double theta, double phi, double **angles, int Nangles);

extern void VT_AddConicField(float ***array, int *dimArray, unsigned char ***field, int *dimField, float coef, int x, int y, int z );

#ifdef __cplusplus
}
#endif

#endif /* _vt_meristemeFormationAxes_h_ */

