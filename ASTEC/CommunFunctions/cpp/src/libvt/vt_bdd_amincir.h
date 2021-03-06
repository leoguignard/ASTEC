#ifndef _vt_bdd_amincir_h_
#define _vt_bdd_amincir_h_

#ifdef __cplusplus
extern "C" {
#endif

#include <vt_common.h>
#include <vt_amincir.h>
#include <vt_amliste.h>
#include <vt_bdd.h>
#include <vt_gb.h>

#ifndef NO_PROTO
extern int _VT_BDD_GONGBERTRAND( vt_image *im, vt_amincir *par );
extern int _VT_BDD_GREG_PLANES( vt_image *im, vt_amincir *par );
extern int _VT_BDD_GREG_CURVES( vt_image *im, vt_amincir *par );
#else
extern int _VT_BDD_GONGBERTRAND();
extern int _VT_BDD_GREG_PLANES();
extern int _VT_BDD_GREG_CURVES();
#endif

#ifdef __cplusplus
}
#endif

#endif /* _vt_bdd_amincir_h_ */
