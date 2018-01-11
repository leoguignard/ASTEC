/*************************************************************************
 * bal-matching.h -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */



#ifndef BAL_BLOCKMATCHING_H
#define BAL_BLOCKMATCHING_H

#ifdef __cplusplus
extern "C" {
#endif

#include <bal-blockmatching-param.h> 
#include <bal-image.h>
#include <bal-transformation.h>




extern void BAL_SetVerboseInBalMatching( int v );
extern void BAL_IncrementVerboseInBalMatching(  );
extern void BAL_DecrementVerboseInBalMatching(  );
extern void BAL_SetDebugInBalMatching( int v );
extern void BAL_IncrementDebugInBalMatching(  );
extern void BAL_DecrementDebugInBalMatching(  );
extern void BAL_SetTraceInBalMatching( int v );
extern void BAL_IncrementTraceInBalMatching(  );
extern void BAL_DecrementTraceInBalMatching(  );

extern int BAL_PyramidalBlockMatching( bal_image *theInrimage_ref, /* reference image 
							     */
				       bal_image *theInrimage_flo, /* floating image
								    */  
				       bal_transformation *theTr, /* allows to goes from Iref to Iflo 
								     (after initial transformation, if any)
								     ie to resample Iflo into Iref
								     that is T_{floating <- reference}
								  */
				       bal_transformation *theInit, /* initial transformation
								     */
				       bal_blockmatching_pyramidal_param *param );

extern int BAL_BlockMatching( bal_image *theInrimage_ref, /* reference image 
							     (may be subsampled) */
			      bal_image *theInrimage_flo, /* floating image
							     (in its original geometry,
							     may be filtered) */  
			      bal_transformation *theTr, /* allows to goes from Iref to Iflo 
							    (after initial transformation, if any)
							    ie to resample Iflo into Iref
							    that is T_{floating <- reference}
							 */
			      bal_transformation *theInit, /* initial transformation
							    */
			      bal_blockmatching_param *param );

#ifdef __cplusplus
}
#endif

#endif
