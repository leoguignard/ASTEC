/***************************************************************************
 * powell.h
 *
 * $Author: greg $
 * $Revision: 1.1 $
 * $Log: powell.h,v $
 * Revision 1.1  2002/10/18 18:02:07  greg
 * *** empty log message ***
 *
 * Revision 1.3  2001/12/11 09:13:47  greg
 * Simplification des interfaces
 *
 * Revision 1.2  2001/12/07 15:48:11  greg
 * En-tete
 *
 *
 *
 * $Id: powell.h,v 1.1 2002/10/18 18:02:07 greg Exp $
 ***************************************************************************/

#ifndef POWELL_H
#define POWELL_H


#ifdef __cplusplus
extern "C" {
#endif

#include <macro.h> 


/* --- Constantes utilisees par la fonction MNBRAK --- */
#define GLIMIT 100.0
#define TINY 1.0e-20 
#define GOLD 1.618034


/* --- Constantes utilisees par la fonction BRENT --- */
#define ITMAXBRENT 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-15 /*1.0e-10*/ 
/*#define TOL 1.0e-2   TOL 3.0e-8 */

/* --- Constantes utilisees par la fonction LINMIN --- */
#define INITAX -0.1
#define INITXX  0.0



/* --- Constantes utilisees par la fonction POWELL --- */
#define ITMAXPOWELL 200

#define ITMAXCAUCHY 200


/* --- Macros --- */
  /*#define SQR(a) ((a) * (a))
    #define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
    #define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d)
  */

/* --- SPECIFICATION FONCTIONS --- */


extern double Brent (double ax, 
		  double bx, 
		  double cx, 
		  double tol,
		  double *xmin,
		  double *pt, /* le point */
		  double *di, /* la direction */
		  double *aux, 
		  int n,
		  double (*func)(double *, void *),
		  void *par
		  );


extern void Mnbrak(double *ax, 
		double *bx, 
		double *cx, 
		double *pt, /* le point */
		double *di, /* la direction */
		double *aux,
		int n,
		double (*func)(double *, void *),
		void *par );

extern void Linmin (double *p, /* le point */
		 double *xi, /* la direction */
		 double *aux, 
		 int n, /* nombre de parametres */
		 double *fret, 
		 double tol,
		 double (*func)(double *, void *),
		 void *par );

extern void Powell (double *p, 
		 int n, 
		 double tol, 
		 double ftol, 
		 int *iter, 
		 double *fret, 
                 double (*func)(double *, void *),
		 void *par );

extern void SteepestDescent (double *p, 
	     int n, 
	     double tol, 
	     double ftol, 
	     int *iter, 
	     double *fret,
	     double (*func)(double *, void *),
	     void (*dfunc)(double *, double *, void *),
	     double * scl,
	     void *par );


extern void ConjugateGradient (double *p, 
	     int n, 
	     double tol, 
	     double ftol, 
	     int *iter, 
	     double *fret,
	     double (*func)(double *, void *),
	     void (*dfunc)(double *, double *, void *),
	     double * scl,
	     void *par );

void QuasiNewton (double *p, 
	     int n, 
	     double tol, 
	     double ftol, 
	     int *iter, 
	     double *fret,
	     double (*func)(double *, void *),
	     void (*dfunc)(double *, double *, void *),
	     double * scl,
	     void *par );

extern double ** _Powell_Matrix ( int p, int q );
extern void _Powell_FreeMatrix ( double ** H, int p, int q );
extern void _Powell_MultMatVect ( double ** A, double * x, double * ax, int p, int q );
extern void _Powell_OuterProduct ( double * x, double * y, double ** A, int p );


extern void PowellVerboseON();
extern void PowellVerboseOFF();
extern void PowellSetVerbose( int v );
extern void PowellVerboseOnStdout( );
extern void PowellVerboseOnStderr( );

#ifdef __cplusplus
}
#endif

#endif
