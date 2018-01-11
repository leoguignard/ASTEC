/*************************************************************************
 * vt_interpol.c - procedures for image interpolation
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * July, 6 1999
 *
 * ADDITIONS, CHANGES
 *
 * - Thu Aug 10 20:28:26 MEST 2000, G. Malandain
 *   Ajout du cas USHORT dans VT_CountIntensityLevels()
 * 
 * - Fri Aug  4 16:57:17 MET DST 2000, G. Malandain
 *   ajout de VT_CountIntensityLevels()
 *
 */

#include <vt_interpol.h>






int VT_InitDistanceImageFromSlice( vt_image *theIm,
				   vt_image *theDist,
				   int index )
{
  char *proc="VT_InitDistanceImageFromSlice";
  short int *resBuf = (short int *)theDist->buf;

  if ( theDist->type != SSHORT ) {
    VT_Error( "bad distance image type", proc );
    return( -1 );
  }
  if ( theIm->dim.x != theDist->dim.x || theIm->dim.y != theDist->dim.y ) {
    VT_Error( "images have different slice dimensions", proc );
    return( -1 );
  }
  if ( index < 0 || index >= theIm->dim.z ) {
    VT_Error( "bad z index", proc );
    return( -1 );
  }

  switch ( theIm->type ) {
  case UCHAR :
    return( _InitDistanceImageFromSlice( (void*)&(((u8***)(theIm->array))[index][0][0]),
					 theIm->type, theIm->dim.x, theIm->dim.y,
					 resBuf,
					 theDist->dim.z ) );
    break;
  case USHORT :
    return( _InitDistanceImageFromSlice( (void*)&(((u16***)(theIm->array))[index][0][0]),
					 theIm->type, theIm->dim.x, theIm->dim.y,
					 resBuf,
					 theDist->dim.z ) );
    break;
  default :
    VT_Error( "such type not handled in switch", proc );
    return( -1 );
  }
  return( -1 );
}














void VT_ComputeSliceFromDistances( vt_image *thePrev,
				 vt_image *theNext,
				 vt_image *theRes,
				 int index,
				 double prevCoef,
				 double nextCoef )
{
  char *proc = "_ComputeSliceFromDistances";

  if ( thePrev->type != SSHORT || theNext->type != SSHORT ) {
    VT_Error( "bad distance images type", proc );
    return;
  }
  if ( theRes->dim.x != thePrev->dim.x || theRes->dim.y != thePrev->dim.y ||
       theRes->dim.x != theNext->dim.x || theRes->dim.y != theNext->dim.y ) {
    VT_Error( "images have different slice dimensions", proc );
    return;
  }
  if ( thePrev->dim.z != theNext->dim.z ) {
    VT_Error( "distance images have different z dimensions", proc );
    return;
  }

  if ( index < 0 || index >= theRes->dim.z ) {
    VT_Error( "bad z index", proc );
    return;
  }

  switch ( theRes->type ) {
  case UCHAR :
    _ComputeSliceFromDistances( (void*)&(((u8***)(theRes->array))[index][0][0]),
				theRes->type, theRes->dim.x, theRes->dim.y, 
				thePrev->buf, theNext->buf,
				thePrev->dim.z,
				prevCoef, nextCoef );
    break;
  case USHORT :
    _ComputeSliceFromDistances( (void*)&(((u16***)(theRes->array))[index][0][0]),
				theRes->type, theRes->dim.x, theRes->dim.y, 
				thePrev->buf, theNext->buf,
				thePrev->dim.z,
				prevCoef, nextCoef );
    break;
  default :
    VT_Error( "such type not handled in switch", proc );
  }
}











/* renvoie le nombre de niveau d'intensites
   
   donne egalement l'intensite min et l'intensite max
 */
int *VT_BuildTranslationTable( vt_image *theIm, int *min, int *max, int *nb )
{
  char *proc = "VT_BuildTranslationTable";
  int *intensity = (int*)NULL;
  int nbintensities = 0;
  int maxIntensities = 0;
  int i, v;

  *min = *max = *nb = 0;


  v = theIm->dim.x * theIm->dim.y * theIm->dim.z;

  switch ( theIm->type ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theIm->buf;
      maxIntensities = 256;
      
      intensity = (int *)VT_Malloc( maxIntensities*sizeof(int) );
      if ( intensity == (int *)NULL ) {
	VT_Error( "unable to allocate auxiliary array", proc );
	return( (int*)NULL );
      }

      for ( i=0; i<maxIntensities; i++ ) intensity[i] = 0;
      
      *min = maxIntensities-1;
      *max = 0;

      for ( i=0; i<v; i++ ) {
	if ( intensity[ (int)theBuf[i] ] == 0 ) {
	  intensity[ (int)theBuf[i] ] = 1;
	  nbintensities ++;
	  if ( *min > (int)theBuf[i] ) *min = (int)theBuf[i];
	  if ( *max < (int)theBuf[i] ) *max = (int)theBuf[i];
	}
      }
      
    }
    break;

  case USHORT :
    {
      u16 *theBuf = (u16*)theIm->buf;
      maxIntensities = 65536;
      
     intensity = (int *)VT_Malloc( maxIntensities*sizeof(int) );
      if ( intensity == (int *)NULL ) {
	VT_Error( "unable to allocate auxiliary array", proc );
	return( (int*)NULL );
      }
      for ( i=0; i<maxIntensities; i++ ) intensity[i] = 0;
      
      *min = maxIntensities-1;
      *max = 0;

      for ( i=0; i<v; i++ ) {
	if ( intensity[ (int)theBuf[i] ] == 0 ) {
	  intensity[ (int)theBuf[i] ] = 1;
	  nbintensities ++;
	  if ( *min > (int)theBuf[i] ) *min = (int)theBuf[i];
	  if ( *max < (int)theBuf[i] ) *max = (int)theBuf[i];
	}
      }
      
    }
    break;

    
  default :
    VT_Error( "such image type not handled in switch", proc );
    return( (int*)NULL );
  }

  
  v = 0;
  for ( i=0; i<maxIntensities; i++ ) {
    if ( intensity[i] == 0 ) {
      intensity[i] = -1;
    } else {
      intensity[i] = v;
      v ++;
    }
  }

  if ( v != nbintensities ) {
    VT_Error( "erreur dans le denombrement", proc );
  }
  
  *nb = nbintensities;
  
  return( intensity );
}











void VT_ApplyTranslationTable( vt_image *theIm, int *t )
{
  char *proc = "VT_ApplyTranslationTable";
  int max;
  int i, v;

  v = theIm->dim.x * theIm->dim.y * theIm->dim.z;

  switch ( theIm->type ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theIm->buf;
      max = 255;
      for ( i=0; i<v; i++ ) {
	if ( t[ theBuf[i] ] < 0 ) {
	  theBuf[i] = 0;
	} else if ( t[ theBuf[i] ] > max ) {
	  theBuf[i] = max;
	} else {
	  theBuf[i] = (u8)t[ theBuf[i] ];
	}
      }
    }
    break;
  case USHORT :
    {
      u16 *theBuf = (u16*)theIm->buf;
      max = 65535;
      for ( i=0; i<v; i++ ) {
	if ( t[ theBuf[i] ] < 0 ) {
	  theBuf[i] = 0;
	} else if ( t[ theBuf[i] ] > max ) {
	  theBuf[i] = max;
	} else {
	  theBuf[i] = (u16)t[ theBuf[i] ];
	}
      }
    }
    break;
  default :
    VT_Error( "such image type not handled in switch", proc );
  }
}
