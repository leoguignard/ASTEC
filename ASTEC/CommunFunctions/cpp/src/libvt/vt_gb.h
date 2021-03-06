#ifndef _vt_gb_h_
#define _vt_gb_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_PROTO
extern int _VT_GB_IsSimple( int vois[3][3][3], int *up, int *no, int *ea, int *bo, int *so, int *we );
extern int _VT_GB_IsSimple0( int *vois );
extern int _VT_GB_IsSimple1( int *vois );
extern int _VT_GB_IsSimple2( int *vois );
extern int _VT_GB_IsSimple3( int *vois );
extern int _VT_GB_IsSimple4( int *vois );
extern int _VT_GB_IsSimple5( int *vois );
#else
extern int _VT_GB_IsSimple();
extern int _VT_GB_IsSimple0();
extern int _VT_GB_IsSimple1();
extern int _VT_GB_IsSimple2();
extern int _VT_GB_IsSimple3();
extern int _VT_GB_IsSimple4();
extern int _VT_GB_IsSimple5();
#endif

#ifdef __cplusplus
}
#endif

#endif /* _vt_gb_h_ */
