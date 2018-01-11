/*************************************************************************
 * reech-def.c -
 *
 * $Id: reech-def.begin,v 1.2 2000/10/18 08:45:56 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 *
 *
 * ADDITIONS, CHANGES
 *	
 *	
 *	
 *
 */


/* CAUTION
   DO NOT EDIT THIS FILE,
   UNLESS YOU HAVE A VERY GOOD REASON 
 */


#include <stdio.h>
#include <stdlib.h>

#include <typedefs.h>
#include <chunks.h>

#include <reech-def.h>




#define _CONVERTR_(R) ( R )
#define _CONVERTI_(R) ( (R) >= 0.0 ? ((int)((R)+0.5)) : ((int)((R)-0.5)) )
static int _verbose_ = 1;





void setVerboseInReechDef( int v )
{
  _verbose_ = v;
}

void incrementVerboseInReechDef(  )
{
  _verbose_ ++;
}

void decrementVerboseInReechDef(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}





typedef struct {
  void* theBuf; /* buffer to be resampled */
  int *theDim; /* dimensions of this buffer */
  void* resBuf; /* result buffer */
  int *resDim;  /* dimensions of this buffer */
  void** theDef;   /* deformation field */
  int *defDim; /* dimensions of these buffers */
  double* mat_aft;  /* transformation matrix */
  double* mat_bef;  /* transformation matrix */
  float gain;
  float bias;
} _VectorFieldResamplingParam;
