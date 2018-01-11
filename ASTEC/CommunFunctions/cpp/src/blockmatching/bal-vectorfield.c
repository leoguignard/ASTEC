/*************************************************************************
 * bal-vectorfield.c -
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



#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <chunks.h>

#include <bal-vectorfield.h>
#include <bal-field-tools.h>
#include <bal-behavior.h>



/************************************************************
 *
 * static variables
 *
 ************************************************************/


/* for some generic writing
 */

static int _verbose_ = 1;

void BAL_SetVerboseInBalVectorField( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalVectorField(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalVectorField(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}










/************************************************************
 *
 *
 *
 ************************************************************/



int BAL_SmoothVectorField( bal_transformation* theTrsf, bal_doublePoint *theSigma )
{
  char *proc="BAL_SmoothVectorField";

  switch ( theTrsf->type ) {
  default :
    if ( _verbose_ )
	fprintf( stderr, "%s: transformation is not a vector field\n", proc );
    return( -1 );
    
  case VECTORFIELD_3D :
    if ( BAL_SmoothImage( &(theTrsf->vz), theSigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth Z component\n", proc );
      return( -1 );
    }

  case VECTORFIELD_2D :
    if ( BAL_SmoothImage( &(theTrsf->vx), theSigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth X component\n", proc );
      return( -1 );
    }
    if ( BAL_SmoothImage( &(theTrsf->vy), theSigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth Y component\n", proc );
      return( -1 );
    }
  }
  return( 1 );
}





/************************************************************
 *
 * procedure to parallelize vector field initialization
 * and division by the weights
 *
 ************************************************************/

static double _min_weight_ = 0.00001;

typedef struct {
  float *bufVx;
  float *bufVy;
  float *bufVz;
  float *bufWeights;
  double v;
} _VectorFieldAuxiliaryParam;



static void _InitVectorFieldAuxiliaryParam( _VectorFieldAuxiliaryParam *p )
{
  p->bufVx = (float*)NULL;
  p->bufVy = (float*)NULL;
  p->bufVz = (float*)NULL;
  p->bufWeights = (float*)NULL;
  p->v = _min_weight_;
}



static int _Init2DVectorFieldSubroutine( void *parameter,
					size_t first,
					size_t last )
{
  float *bufVx = ((_VectorFieldAuxiliaryParam*)parameter)->bufVx;
  float *bufVy = ((_VectorFieldAuxiliaryParam*)parameter)->bufVy;
  float *bufWe = ((_VectorFieldAuxiliaryParam*)parameter)->bufWeights;
  size_t i;

  for ( i=first; i<=last; i++ ) {
    bufVx[i] = 0.0;
    bufVy[i] = 0.0;
    bufWe[i] = 0.0;
  }
  return( 1 );
}



static int _Init3DVectorFieldSubroutine( void *parameter,
					size_t first,
					size_t last )
{
  float *bufVx = ((_VectorFieldAuxiliaryParam*)parameter)->bufVx;
  float *bufVy = ((_VectorFieldAuxiliaryParam*)parameter)->bufVy;
  float *bufVz = ((_VectorFieldAuxiliaryParam*)parameter)->bufVz;
  float *bufWe = ((_VectorFieldAuxiliaryParam*)parameter)->bufWeights;
  size_t i;

  for ( i=first; i<=last; i++ ) {
    bufVx[i] = 0.0;
    bufVy[i] = 0.0;
    bufVz[i] = 0.0;
    bufWe[i] = 0.0;
  }
  return( 1 );
}



static int _Normalize2DVectorFieldSubroutine( void *parameter,
					      size_t first,
					      size_t last )
{
  float *bufVx = ((_VectorFieldAuxiliaryParam*)parameter)->bufVx;
  float *bufVy = ((_VectorFieldAuxiliaryParam*)parameter)->bufVy;
  float *bufWe = ((_VectorFieldAuxiliaryParam*)parameter)->bufWeights;
  double v = ((_VectorFieldAuxiliaryParam*)parameter)->v;
  size_t i;

  for ( i=first; i<=last; i++ ) {
    if ( bufWe[i] > v ) {
      bufVx[i] /= bufWe[i];
      bufVy[i] /= bufWe[i];
    }
    else {
      bufVx[i] /= v;
      bufVy[i] /= v;
    }
  }
  return( 1 );
}



static int _Normalize3DVectorFieldSubroutine( void *parameter,
					      size_t first,
					      size_t last )
{
  float *bufVx = ((_VectorFieldAuxiliaryParam*)parameter)->bufVx;
  float *bufVy = ((_VectorFieldAuxiliaryParam*)parameter)->bufVy;
  float *bufVz = ((_VectorFieldAuxiliaryParam*)parameter)->bufVz;
  float *bufWe = ((_VectorFieldAuxiliaryParam*)parameter)->bufWeights;
  double v = ((_VectorFieldAuxiliaryParam*)parameter)->v;
  size_t i;

  for ( i=first; i<=last; i++ ) {
    if ( bufWe[i] > v ) {
      bufVx[i] /= bufWe[i];
      bufVy[i] /= bufWe[i];
      bufVz[i] /= bufWe[i];
    }
    else {
      bufVx[i] /= v;
      bufVy[i] /= v;
      bufVz[i] /= v;
    }
  }
  return( 1 );
}





/************************************************************
 *
 *
 *
 ************************************************************/



static int _LS_VectorField_Estimation( bal_transformation* theTrsf, 
				       FIELD * field, 
				       bal_estimator *estimator )
{
  char *proc="_LS_VectorField_Estimation";
  bal_image theWeights;
  float ***bufWeights, ***bufVx, ***bufVy, ***bufVz;
  int i, j, k;
  int n;
  bal_doublePoint sigma = estimator->sigma;

  typeChunks chunks;
  size_t first, last;
  _VectorFieldAuxiliaryParam aux;



  /* initialization
   */
  if ( BAL_InitAllocImage( &theWeights, "weights.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate weight image\n", proc );
    return ( -1 );
  }
  bufWeights = (float***)theWeights.array;



  /* parallelization stuff 
   */
  first = 0;
  last = theTrsf->vx.ncols*theTrsf->vx.nrows*theTrsf->vx.nplanes - 1;
  _InitVectorFieldAuxiliaryParam( &aux );


  initChunks( &chunks );
  if ( 0 )
    fprintf( stderr, "\n%s: build chunks [%lu %lu]\n\n", proc, first, last );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    BAL_FreeImage( &theWeights );
    return( -1 );
  }
  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {	
    chunks.data[n].parameters = (void*)(&aux);
  }


  /* processing
   */
  switch ( theTrsf->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet (but weird)\n", proc );
    freeChunks( &chunks );
    BAL_FreeImage( &theWeights );
    return( -1 );
    
  case VECTORFIELD_2D :
  
    bufVx = (float***)(theTrsf->vx.array);
    bufVy = (float***)(theTrsf->vy.array);
    
    aux.bufVx = (float*)(theTrsf->vx.data);
    aux.bufVy = (float*)(theTrsf->vy.data);
    aux.bufWeights = (float*)(theWeights.data);

    if ( processChunks( &_Init2DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to initialize vectorfield (2D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    for ( n=0; n<field->n_selected_pairs; n++ ) {
      i = (int)(field->pointer[n]->origin.x + 0.5);
      j = (int)(field->pointer[n]->origin.y + 0.5);
      k = (int)(field->pointer[n]->origin.z + 0.5);
      bufVx[k][j][i] += field->pointer[n]->vector.x;
      bufVy[k][j][i] += field->pointer[n]->vector.y;
      bufWeights[k][j][i] ++;
    }

    sigma.z = -1.0;

    if ( BAL_SmoothVectorField( theTrsf, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth vectorfield (2D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    
    if ( BAL_SmoothImage( &theWeights, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth weights (2D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    
    if ( processChunks( &_Normalize2DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to normalize vectorfield (2D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    break;
    
  case VECTORFIELD_3D :
  
    bufVx = (float***)(theTrsf->vx.array);
    bufVy = (float***)(theTrsf->vy.array);
    bufVz = (float***)(theTrsf->vz.array);

    aux.bufVx = (float*)(theTrsf->vx.data);
    aux.bufVy = (float*)(theTrsf->vy.data);
    aux.bufVz = (float*)(theTrsf->vz.data);
    aux.bufWeights = (float*)(theWeights.data);

    if ( processChunks( &_Init3DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to initialize vectorfield (3D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    for ( n=0; n<field->n_selected_pairs; n++ ) {
      i = (int)(field->pointer[n]->origin.x + 0.5);
      j = (int)(field->pointer[n]->origin.y + 0.5);
      k = (int)(field->pointer[n]->origin.z + 0.5);
      bufVx[k][j][i] += field->pointer[n]->vector.x;
      bufVy[k][j][i] += field->pointer[n]->vector.y;
      bufVz[k][j][i] += field->pointer[n]->vector.z;
      bufWeights[k][j][i] ++;
    }

    if ( BAL_SmoothVectorField( theTrsf, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth vectorfield (3D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    if ( BAL_SmoothImage( &theWeights, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth weights (3D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    
    if ( processChunks( &_Normalize3DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to normalize vectorfield (3D LS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );
  BAL_FreeImage( &theWeights );

  return( 1 );

}




static int _WLS_VectorField_Estimation( bal_transformation* theTrsf, 
				       FIELD * field, 
				       bal_estimator *estimator )
{
  char *proc="_WLS_VectorField_Estimation";
  bal_image theWeights;
  float ***bufWeights, ***bufVx, ***bufVy, ***bufVz;
  int i, j, k;
  int n;
  bal_doublePoint sigma = estimator->sigma;

  typeChunks chunks;
  size_t first, last;
  _VectorFieldAuxiliaryParam aux;



  /* initialization
   */
  if ( BAL_InitAllocImage( &theWeights, "weights.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT) == -1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate weight image\n", proc );
    return ( -1 );
  }
  bufWeights = (float***)theWeights.array;



  /* parallelization stuff 
   */
  first = 0;
  last = theTrsf->vx.ncols*theTrsf->vx.nrows*theTrsf->vx.nplanes - 1;
  _InitVectorFieldAuxiliaryParam( &aux );


  initChunks( &chunks );
  if ( 0 )
    fprintf( stderr, "\n%s: build chunks [%lu %lu]\n\n", proc, first, last );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    BAL_FreeImage( &theWeights );
    return( -1 );
  }
  for ( n=0; n<chunks.n_allocated_chunks; n++ ) {	
    chunks.data[n].parameters = (void*)(&aux);
  }


  /* processing
   */
  switch ( theTrsf->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet (but weird)\n", proc );
    freeChunks( &chunks );
    BAL_FreeImage( &theWeights );
    return( -1 );
    
  case VECTORFIELD_2D :
  
    bufVx = (float***)(theTrsf->vx.array);
    bufVy = (float***)(theTrsf->vy.array);

    aux.bufVx = (float*)(theTrsf->vx.data);
    aux.bufVy = (float*)(theTrsf->vy.data);
    aux.bufWeights = (float*)(theWeights.data);

    if ( processChunks( &_Init2DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to initialize vectorfield (2D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    for ( n=0; n<field->n_selected_pairs; n++ ) {
      i = (int)(field->pointer[n]->origin.x + 0.5);
      j = (int)(field->pointer[n]->origin.y + 0.5);
      k = (int)(field->pointer[n]->origin.z + 0.5);
      bufVx[k][j][i] += field->pointer[n]->rho * field->pointer[n]->vector.x;
      bufVy[k][j][i] += field->pointer[n]->rho * field->pointer[n]->vector.y;
      bufWeights[k][j][i] += field->pointer[n]->rho;
    }

    sigma.z = -1.0;

    if ( BAL_SmoothVectorField( theTrsf, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth vectorfield (2D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    if ( BAL_SmoothImage( &theWeights, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth weights (2D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    
    if ( processChunks( &_Normalize2DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to normalize vectorfield (2D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    break;
    
  case VECTORFIELD_3D :
  
    bufVx = (float***)(theTrsf->vx.array);
    bufVy = (float***)(theTrsf->vy.array);
    bufVz = (float***)(theTrsf->vz.array);

    aux.bufVx = (float*)(theTrsf->vx.data);
    aux.bufVy = (float*)(theTrsf->vy.data);
    aux.bufVz = (float*)(theTrsf->vz.data);
    aux.bufWeights = (float*)(theWeights.data);

    if ( processChunks( &_Init3DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to initialize vectorfield (3D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    for ( n=0; n<field->n_selected_pairs; n++ ) {
      i = (int)(field->pointer[n]->origin.x + 0.5);
      j = (int)(field->pointer[n]->origin.y + 0.5);
      k = (int)(field->pointer[n]->origin.z + 0.5);
      bufVx[k][j][i] += field->pointer[n]->rho * field->pointer[n]->vector.x;
      bufVy[k][j][i] += field->pointer[n]->rho * field->pointer[n]->vector.y;
      bufVz[k][j][i] += field->pointer[n]->rho * field->pointer[n]->vector.z;
      bufWeights[k][j][i] += field->pointer[n]->rho;
    }

    if ( BAL_SmoothVectorField( theTrsf, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth vectorfield (3D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    if ( BAL_SmoothImage( &theWeights, &sigma ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to smooth weights (3D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }
    
    if ( processChunks( &_Normalize3DVectorFieldSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to normalize vectorfield (3D WLS case)\n", proc );
      freeChunks( &chunks );
      BAL_FreeImage( &theWeights );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );
  BAL_FreeImage( &theWeights );

  return( 1 );

}





/************************************************************
 *
 *
 *
 ************************************************************/



typedef struct {
  bal_transformation *theTrsf;
  FIELD *field;
} _VectorFieldResidualsParam;





static int _VectorField_Residuals_2DSubroutine( void *parameter,
						size_t first,
						size_t last )
{
  FIELD *field = ((_VectorFieldResidualsParam*)parameter)->field;
  float ***bufVx = (float***)(((_VectorFieldResidualsParam*)parameter)->theTrsf->vx.array);
  float ***bufVy = (float***)(((_VectorFieldResidualsParam*)parameter)->theTrsf->vy.array);

  size_t n;
  int i, j, k;
  double dx, dy;

  for ( n = first; n <= last; n ++ ) {
    i = (int)(field->pointer[n]->origin.x + 0.5);
    j = (int)(field->pointer[n]->origin.y + 0.5);
    k = (int)(field->pointer[n]->origin.z + 0.5);
    dx = bufVx[k][j][i] - field->pointer[n]->vector.x;
    dy = bufVy[k][j][i] - field->pointer[n]->vector.y;
    field->pointer[n]->error = dx*dx + dy*dy;
  }
  return( 1 );
}



static int _VectorField_Residuals_3DSubroutine( void *parameter,
						size_t first,
						size_t last )
{
  FIELD *field = ((_VectorFieldResidualsParam*)parameter)->field;
  float ***bufVx = (float***)(((_VectorFieldResidualsParam*)parameter)->theTrsf->vx.array);
  float ***bufVy = (float***)(((_VectorFieldResidualsParam*)parameter)->theTrsf->vy.array);
  float ***bufVz = (float***)(((_VectorFieldResidualsParam*)parameter)->theTrsf->vz.array);

  size_t n;
  int i, j, k;
  double dx, dy, dz;

  for ( n = first; n <= last; n ++ ) {
    i = (int)(field->pointer[n]->origin.x + 0.5);
    j = (int)(field->pointer[n]->origin.y + 0.5);
    k = (int)(field->pointer[n]->origin.z + 0.5);
    dx = bufVx[k][j][i] - field->pointer[n]->vector.x;
    dy = bufVy[k][j][i] - field->pointer[n]->vector.y;
    dz = bufVz[k][j][i] - field->pointer[n]->vector.z;
    field->pointer[n]->error = dx*dx + dy*dy + dz*dz; 
  }
  return( 1 );
}



int BAL_VectorField_Residuals( bal_transformation* theTrsf, 
				  FIELD * field )
{
  char *proc = "BAL_VectorField_Residuals";
  int n;

  _VectorFieldResidualsParam p;
  typeChunks chunks;
  size_t first, last;
  first = 0;
  last = field->n_computed_pairs-1;

  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  
  p.theTrsf = theTrsf;
  p.field = field;

  for ( n=0; n<chunks.n_allocated_chunks; n++ )	
    chunks.data[n].parameters = (void*)(&p);





  switch( theTrsf->type ) {
  default :
    if ( _verbose_ )
      fprintf (stderr, "%s: such transformation not handled in switch\n", proc );
    freeChunks( &chunks );
    return( -1 );
    
  case VECTORFIELD_2D :
    
    if ( processChunks( &_VectorField_Residuals_2DSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute residuals (2D vector field case)\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    break;
      
  case VECTORFIELD_3D :
      
    if ( processChunks( &_VectorField_Residuals_3DSubroutine, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute residuals (3D vector field case)\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    break;
    
  }
  
  freeChunks( &chunks );
  return( 1 );

}




static int _VectorField_Estimation( bal_transformation* theTrsf, 
				    FIELD * field, 
				    bal_estimator *estimator )
{
  char *proc="_VectorField_Estimation";
  switch( estimator->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: estimator type not handled in switch", proc );
    return( -1 );
  case TYPE_LS :
  case TYPE_LTS :
    return( _LS_VectorField_Estimation( theTrsf, field, estimator ) );
  case TYPE_WLS :
  case TYPE_WLTS :
    return( _WLS_VectorField_Estimation( theTrsf, field, estimator ) );
  }
  return( -1 );



}






static int _VectorField_Trimmed_Estimation( bal_transformation* theTrsf, 
					    FIELD * field, 
					    bal_estimator *estimator )
{
  char *proc="_VectorField_Trimmed_Estimation";

  int nretainedpairs;

  int iter;



  /* initial transformation estimation
   */
  field->n_selected_pairs = field->n_computed_pairs;
  if ( _VectorField_Estimation( theTrsf, field, estimator ) != 1 ) {
    if ( _verbose_ )
      fprintf ( stderr, "%s: error when estimating initial transformation\n", proc );
    return( -1 );
  }



  /* loop initialisation
   */
  for ( iter = 0; iter < estimator->max_iterations; iter ++ ) {

    /* calcul des residus
     */
    if ( BAL_VectorField_Residuals( theTrsf, field ) != 1 ) {
      if ( _verbose_ )
	fprintf (stderr, "%s: unable to compute residuals\n", proc );
      return( -1 );
    }


    /* sort residuals
     */
    nretainedpairs = BAL_SelectSmallestResiduals( field, estimator );
    if ( nretainedpairs <= 0 ) {
      if ( _verbose_ )
	fprintf (stderr, "%s: no retained residuals? Returned value was %d\n", proc, nretainedpairs );
      return( -1 );
    }
    field->n_selected_pairs = nretainedpairs;
    
    
    
    /* transformation estimation
     */
    if ( _VectorField_Estimation( theTrsf, field, estimator ) != 1 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: something goes wrong in the transformation estimation\n", proc );
	fprintf( stderr, "\t iteration #%d of the iterated (trimmed) estimation\n", iter );
      }
      return( -1 );
    }


    
    if ( _verbose_ >= 3 ) {
#ifdef _ORIGINAL_BALADIN_PRINTS_
      switch ( estimator->type ) {
      default: break;
      case TYPE_LS :
      case TYPE_LTS :
	fprintf( stderr, "      LTS Iteration n. %d \r", iter);
	break;
      case TYPE_WLS :
      case TYPE_WLTS :
	fprintf( stderr, "      WLTS Iteration n. %d \r", iter);
	break;
      }
#else
      switch ( estimator->type ) {
      default: break;
      case TYPE_LS :
      case TYPE_LTS :
	fprintf( stderr, "      LTS: iteration #%2d ... \r", iter ); 
	break;
      case TYPE_WLS :
      case TYPE_WLTS :
	fprintf( stderr, "      WLTS: iteration #%2d ... \r", iter ); 
	break;
      }
#endif
    }

  }
    
  return( 1 );
}




/* Compute the transformation from the reference image 
   towards the floating one, thus allows to resample the floating
   in the reference frame.
   
   Pairings are in voxel units.
*/
int BAL_ComputeVectorFieldTransformation( bal_transformation* theTrsf, 
					  FIELD * field, 
					  bal_estimator *estimator )
{
  char *proc= "BAL_ComputeVectorFieldTransformation";
  switch ( estimator->type ) {
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation estimation not handled yet\n", proc );
    return( -1 );
  case TYPE_LS :
  case TYPE_WLS :
    return( _VectorField_Estimation( theTrsf, field, estimator ) );
  case TYPE_LTS :
  case TYPE_WLTS :
    return( _VectorField_Trimmed_Estimation( theTrsf, field, estimator ) );
  }
  
  return( 1 );
}
