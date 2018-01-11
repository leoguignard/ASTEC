#ifndef _vt_etoilebarre_h_
#define _vt_etoilebarre_h_

#ifdef __cplusplus
extern "C" {
#endif



#ifndef NO_PROTO

extern int VT_ComputeCetoile( int voisins[27] );
extern int VT_ComputeCbarre( int voisins[27] );

#else 

extern int VT_ComputeCetoile();
extern int VT_ComputeCbarre();

#endif



#ifdef __cplusplus
}
#endif

#endif /* _vt_etoilebarre_h_ */
