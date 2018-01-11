#ifndef _vt_cc_h_
#define _vt_cc_h_

#ifdef __cplusplus
extern "C" {
#endif

#ifndef NO_PROTO
extern int _VT_26CCin27V( int voisins[27] );
extern int _VT_26CCin26V( int voisins[27] );
extern int _VT_26CCin26VisOne( int voisins[27] );
extern int _VT_06CCin19V( int voisins[27] );
extern int _VT_06CCin18V( int voisins[27] );
extern int _VT_06CCin18VisOne( int voisins[27] );
#else
extern int _VT_26CCin27V();
extern int _VT_26CCin26V();
extern int _VT_26CCin26VisOne();
extern int _VT_06CCin19V();
extern int _VT_06CCin18V();
extern int _VT_06CCin18VisOne();
#endif /* NO_PROTO */

#ifdef __cplusplus
}
#endif

#endif /* _vt_cc_h_ */
