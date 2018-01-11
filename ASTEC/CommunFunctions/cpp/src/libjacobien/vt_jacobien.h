#ifndef _vt_jacobien_h_
#define _vt_jacobien_h_



#ifdef __cplusplus
extern "C" {
#endif

extern void JCB_InitRandom( int seed );
extern double JCB_Random();
extern void JCB_PrintMatrix( FILE *f, double *mat, char *str );
extern void JCB_RndTranslation( double *mat, double tmin, double tmax );
extern void JCB_RndTranslationMatrix( double *mat, double tmin, double tmax );
extern void JCB_RndScaleMatrix( double *mat, double smin, double smax, double tmin, double tmax );
extern void JCB_RndSimilitudeMatrix( double *mat, double smin, double smax, double tmin, double tmax );
extern void JCB_RndAffineMatrix( double *mat, double amin, double amax, double tmin, double tmax );


typedef struct {
  double x;
  double y;
  double z;
} typePoint3D;

typedef struct {
  int n;
  typePoint3D *pts;
} typeListPoint3D;

extern void JCB_initListPoint3D( typeListPoint3D *l );
extern int JCB_allocListPoint3D( typeListPoint3D *l, int n );
extern void JCB_freeListPoint3D( typeListPoint3D *l );
extern void JCB_PrintListPoint3D( FILE *f, typeListPoint3D *l , char *str );
extern int JCB_get6Neighborhood( typeListPoint3D *l );
extern void JCB_transformListPoint3D( typeListPoint3D *theList, typeListPoint3D *resList, double *mat );

extern void JCB_TestTranslation( int n );
extern void JCB_TestSimilitude( int n );
extern void JCB_TestAffine( int n );
extern double JCB_ComputeJacobien( typeListPoint3D *r, typeListPoint3D *l );

#ifdef __cplusplus
}
#endif

#endif /* _vt_jacobien_h_ */
