/*************************************************************************
 * bal-transformation-tools.c -
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

#include <bal-transformation-tools.h>
#include <bal-lineartrsf.h>
#include <bal-lineartrsf-tools.h>
#include <bal-vectorfield.h>

#include <chunks.h>
#include <reech4x4.h>
#include <reech-def.h>

static int _verbose_ = 1;
static int _debug_ = 0;



void BAL_SetVerboseInBalTransformationTools( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalTransformationTools(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalTransformationTools(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}





/*************************************************************
 *
 * static functions
 *
 ************************************************************/

static int _2DVectorField2DVectorFieldComposition( bal_transformation *t1, 
						bal_transformation *t2, 
						bal_transformation *res );
static int _3DVectorField3DVectorFieldComposition( bal_transformation *t1, 
						bal_transformation *t2, 
						bal_transformation *res );
static int _2DVectorFieldMatrixComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					   bal_transformation *res );
static int _3DVectorFieldMatrixComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					   bal_transformation *res );
static int _Matrix2DVectorFieldComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					    bal_transformation *res );
static int _Matrix3DVectorFieldComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					    bal_transformation *res );

static int BAL_Reech3DTriLin4x4( bal_image *theIm, bal_image *resIm,
				 bal_transformation *theTr );
static int BAL_Reech3DTriLin2DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr );
static int BAL_Reech3DTriLin3DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr );

static int BAL_Reech3DNearest4x4( bal_image *theIm, bal_image *resIm,
				 bal_transformation *theTr );
static int BAL_Reech3DNearest2DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr );
static int BAL_Reech3DNearest3DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr );

static int BAL_Inverse2DVectorField( bal_transformation *theTrsf,  
				     bal_transformation *invTrsf );
static int BAL_Inverse3DVectorField( bal_transformation *theTrsf,  
				     bal_transformation *invTrsf );





/*--------------------------------------------------
 *
 * TRANSFORMATION COMPUTATION
 *
 --------------------------------------------------*/




/* compute the subsampling matrix T (in 'real' coordinates)
   that allows to subsample the image 'image_to_be_subsampled'
   into the geometry of 'subsampled_image'

   It goes then from 'subsampled_image' to 'image_to_be_subsampled'
   and we have 'subsampled_image = image_to_be_subsampled o T'

   It puts the center of the field of view of 'subsampled_image' on the 
   center of the field of view of 'image_to_be_subsampled'.
   Assuming that the origin is at the center of the upper left voxel, 
   the coordinates of the center of the FOV are 
   ((dx-1)*vx/2, (dy-1)*vy/2, (dz-1)*vz/2 )
  
*/
int BAL_ComputeImageToImageTransformation( bal_image *subsampled_image,
					   bal_image *image_to_be_subsampled,
					   bal_transformation *subsampling_trsf )
{
  char *proc ="BAL_ComputeImageToImageTransformation";
  _MATRIX mat;
  size_t i, v;
  bal_doublePoint image1center;
  bal_doublePoint image2center;
  float *vx, *vy, *vz;

  if ( _alloc_mat( &mat, 4, 4) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: error when allocating matrix\n", proc );
    return( -1 );
  }
  for ( i=0; i<16; i++ ) mat.m[i] = 0.0;
  mat.m[0] = mat.m[5] = mat.m[10] = mat.m[15] = 1.0;

  image1center.x = image_to_be_subsampled->vx * (double)(image_to_be_subsampled->ncols - 1) / 2.0;
  image1center.y = image_to_be_subsampled->vy * (double)(image_to_be_subsampled->nrows - 1) / 2.0;
  image1center.z = image_to_be_subsampled->vz * (double)(image_to_be_subsampled->nplanes - 1) / 2.0;

  image2center.x = subsampled_image->vx * (double)(subsampled_image->ncols - 1) / 2.0;
  image2center.y = subsampled_image->vy * (double)(subsampled_image->nrows - 1) / 2.0;
  image2center.z = subsampled_image->vz * (double)(subsampled_image->nplanes - 1) / 2.0;

  mat.m[ 3] = image1center.x - image2center.x;
  mat.m[ 7] = image1center.y - image2center.y;
  mat.m[11] = image1center.z - image2center.z;

  switch ( subsampling_trsf->type ) {

  default :

    if ( _verbose_ ) 
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    if ( subsampling_trsf->mat.m ==  (double*)NULL ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: matrix is not allocated\n", proc );
      return( -1 );
    }

    for ( i=0; i<16; i++ ) 
      subsampling_trsf->mat.m[i] = mat.m[i];;
    break;

  case VECTORFIELD_2D :

    if ( subsampling_trsf->vx.ncols != subsampled_image->ncols
	 || subsampling_trsf->vx.nrows != subsampled_image->nrows
	 || subsampling_trsf->vx.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    if ( subsampling_trsf->vy.ncols != subsampled_image->ncols
	 || subsampling_trsf->vy.nrows != subsampled_image->nrows
	 || subsampling_trsf->vy.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: Y component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    
    v = subsampling_trsf->vx.ncols * subsampling_trsf->vx.nrows * subsampling_trsf->vx.nplanes;
    vx = (float*)(subsampling_trsf->vx.data);
    vy = (float*)(subsampling_trsf->vy.data);

    for ( i=0; i<v; i++ ) {
      vx[i] =  mat.m[ 3];
      vy[i] =  mat.m[ 7];
    }
    break;

  case VECTORFIELD_3D :

    if ( subsampling_trsf->vx.ncols != subsampled_image->ncols
	 || subsampling_trsf->vx.nrows != subsampled_image->nrows
	 || subsampling_trsf->vx.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    if ( subsampling_trsf->vy.ncols != subsampled_image->ncols
	 || subsampling_trsf->vy.nrows != subsampled_image->nrows
	 || subsampling_trsf->vy.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: Y component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }
    
    if ( subsampling_trsf->vz.ncols != subsampled_image->ncols
	 || subsampling_trsf->vz.nrows != subsampled_image->nrows
	 || subsampling_trsf->vz.nplanes != subsampled_image->nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: Z component of vector field has different dimensions than 'subsampled image'\n", proc );
      return( -1 );
    }

    v = subsampling_trsf->vx.ncols * subsampling_trsf->vx.nrows * subsampling_trsf->vx.nplanes;
    vx = (float*)(subsampling_trsf->vx.data);
    vy = (float*)(subsampling_trsf->vy.data);
    vz = (float*)(subsampling_trsf->vz.data);
    for ( i=0; i<v; i++ ) {
      vx[i] =  mat.m[ 3];
      vy[i] =  mat.m[ 7];
      vz[i] =  mat.m[11];
    }
    break;

  }
  _free_mat( &mat );

  /* the transformation has been calculated in real units
   */
  subsampling_trsf->transformation_unit = REAL_UNIT;

  return ( 1 );
}










/*--------------------------------------------------
 *
 * transformation from a pairing field
 *
 --------------------------------------------------*/



int BAL_ComputeIncrementalTransformation( bal_transformation *T,  FIELD *field, 
					  bal_estimator *estimator )

{
  char *proc = "BAL_ComputeIncrementalTransformation";
  
  switch( T->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    /* le champ est en voxel,
       avec le point de depart dans l'image de reference
       on le passe en reel (important pour la transformation rigide)
    */
    BAL_ChangeFieldToRealUnit( field );
    
    if ( BAL_ComputeLinearTransformation( T, field, estimator ) <= 0 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: incremental linear transformation computation failed\n", proc );
      return( -1 );
    }
    
    /* the transformation has been calculated in real units
     */
    T->transformation_unit = REAL_UNIT;

    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    /* le champ est en voxel,
       avec le point de depart dans l'image de reference
    */
    BAL_ChangeFieldToVoxelUnit( field );

    if ( BAL_ComputeVectorFieldTransformation( T, field, estimator ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: incremental vector field transformation computation failed\n", proc );
      return( -1 );
    }

    /* the transformation has been calculated in voxel units
     */
    T->transformation_unit = VOXEL_UNIT;

    /* change into real units
     */
    if ( BAL_ChangeTransformationToRealUnit( &(T->vx), &(T->vx), T, T ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error when translating transformation vector field from voxel to real\n", proc );
      return( -1 );
    }
    
  }

  return( 1 );

}





int BAL_ComputeTransformationResiduals( bal_transformation *T,  FIELD *field )
{ 
  char *proc = "BAL_ComputeTransformationResiduals";

  switch( T->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    /* check the units
    */
    switch ( T->transformation_unit ) {
    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: weird, transformation type is not defined\n", proc );
	fprintf( stderr, "\t set it to REAL_UNIT, but this has to be fixed\n" );
      }
      T->transformation_unit = REAL_UNIT;
    case REAL_UNIT :
      BAL_ChangeFieldToRealUnit( field );
      break;
    case VOXEL_UNIT :
      BAL_ChangeFieldToVoxelUnit( field );
      break;
    }

    if ( BAL_LinearTrsf_Residuals( T, field ) <= 0 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: linear transformation residual computation failed\n", proc );
      return( -1 );
    }

    break;
    
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :

    if ( T->transformation_unit != VOXEL_UNIT && field->unit != VOXEL_UNIT ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: units have to be checked\n", proc );
      return( -1 );
    }

    if ( BAL_VectorField_Residuals( T, field  ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: vector field residual transformation computation failed\n", proc );
      return( -1 );
    }
    
  }

  return( 1 );

}





/*--------------------------------------------------
 *
 * TRANSFORMATION COMPOSITION 
 *
 --------------------------------------------------*/



int BAL_TransformationComposition( bal_transformation *t_res, /* t_res = t1 o t2 */
				   bal_transformation *t1, 
				   bal_transformation *t2 )
{
  char *proc = "BAL_TransformationComposition";
  _MATRIX tmp;

  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  switch ( t1->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such first transformation type not handled yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    
    switch ( t2->type ) {
      
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (matrix o unknown)\n", proc );
      return( -1 );
      
    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      
      if ( _alloc_mat( &tmp, 4, 4) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate auxiliary matrix\n", proc );
	return( -1 );
      }
      _mult_mat( &(t1->mat), &(t2->mat), &tmp );
      _copy_mat( &tmp, &(t_res->mat) );
      _free_mat( &tmp );
      break;
    
    case VECTORFIELD_2D :

      if ( _Matrix2DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing 'matrix' o '2D vector field'\n", proc );
	return( -1 );
      }
      break;

    case VECTORFIELD_3D :

      if ( _Matrix3DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing 'matrix' o '3D vector field'\n", proc );
	return( -1 );
      }
      break;
      
    }

    break;
    /* end of matrix case 
     */
    
  case VECTORFIELD_2D :
    
    switch ( t2->type ) {
      
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (2D vector field o unknown)\n", proc );
      return( -1 );
      
    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      if ( _2DVectorFieldMatrixComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing '2D vector field' o 'matrix'\n", proc );
	return( -1 );
      }
      break;
      
    case VECTORFIELD_2D:
      if ( _2DVectorField2DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing '2D vector field' o '2D vector field'\n", proc );
	return( -1 );
      };
      break;
    }

    break;
    /* end of 2D vector field case
     */

  case VECTORFIELD_3D :
    
    switch ( t2->type ) {
      
    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (3D vector field o unknown)\n", proc );
      return( -1 );
      
    case TRANSLATION_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_2D :
    case TRANSLATION_SCALING_3D :
    case RIGID_2D :
    case RIGID_3D :
    case SIMILITUDE_2D :
    case SIMILITUDE_3D :
    case AFFINE_2D :
    case AFFINE_3D :
      if ( _3DVectorFieldMatrixComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing '3D vector field' o 'matrix'\n", proc );
	return( -1 );
      }
      break;
      
    case VECTORFIELD_3D:
      if ( _3DVectorField3DVectorFieldComposition( t1, t2, t_res ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: error when composing '3D vector field' o '3D vector field'\n", proc );
	return( -1 );
      };
      break;
    }

    break;
    /* end of 3D vector field case
     */

  }
  
  t_res->transformation_unit =  t1->transformation_unit;

  return( 1 );
}





enumTypeTransfo BAL_TypeTransformationComposition( bal_transformation *t1, 
						   bal_transformation *t2 )
{
  char *proc = "BAL_TypeTransformationComposition";
  
  switch ( t1->type ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such first transformation type not handled yet\n", proc );
    return( UNDEF_TRANSFORMATION );

  case TRANSLATION_2D :
  case TRANSLATION_SCALING_2D :
  case RIGID_2D :
  case SIMILITUDE_2D :
  case AFFINE_2D :

    /* 2D matrix case */
    switch ( t2->type ) {

    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (2D matrix o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
      return( AFFINE_2D );

    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D : 
      return( AFFINE_3D );

    case VECTORFIELD_2D :
      return( VECTORFIELD_2D );

    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 2D matrix case */

  case TRANSLATION_3D :
  case TRANSLATION_SCALING_3D :
  case RIGID_3D :
  case SIMILITUDE_3D :
  case AFFINE_3D :
    
    /* 3D matrix case */
    switch ( t2->type ) {

    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (3D matrix o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D : 
      return( AFFINE_3D );

    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 3D matrix case */

  case VECTORFIELD_2D :
    
    /* 2D vector field case */
    switch ( t2->type ) {

    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (2D vector field o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
      return( VECTORFIELD_2D );

    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D : 
      return( VECTORFIELD_3D );

    case VECTORFIELD_2D :
      return( VECTORFIELD_2D );

    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 2D vector field case */

  case VECTORFIELD_3D :
    
    /* 3D vector field case */
    switch ( t2->type ) {

    default :
      if ( _verbose_ )
	fprintf( stderr, "%s: such second transformation type not handled yet (3D vector field o unknown)\n", proc );
      return( UNDEF_TRANSFORMATION );
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D : 
    case VECTORFIELD_2D :
    case VECTORFIELD_3D :
      return( VECTORFIELD_3D );

    }
    /* 3D vector field case */
    
  }
  
  if ( _verbose_ )
    fprintf( stderr, "%s: should not reach this ...\n", proc );
  return( UNDEF_TRANSFORMATION );
}





int BAL_AllocTransformationComposition( bal_transformation *res, 
					bal_transformation *t1, 
					bal_transformation *t2,
					bal_image *ref )
{
  char *proc = "BAL_AllocTransformationComposition";
  enumTypeTransfo resType;

  resType = BAL_TypeTransformationComposition( t1, t2 );

   switch ( resType ) {

    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: such type '", proc );
	BAL_PrintTransformationType( stderr, resType );
	fprintf( stderr, "' not handled yet\n" );
      }
      return( -1 );
      
    case TRANSLATION_2D :
    case TRANSLATION_SCALING_2D :
    case RIGID_2D :
    case SIMILITUDE_2D :
    case AFFINE_2D :
    case TRANSLATION_3D :
    case TRANSLATION_SCALING_3D :
    case RIGID_3D :
    case SIMILITUDE_3D :
    case AFFINE_3D : 
      if ( BAL_AllocTransformation( res, resType, (bal_image *)NULL ) != 1 ) {
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to allocate transformation (linear case)", proc );
	return( -1 );
      }
      return( 1 );
      
    case VECTORFIELD_2D :
      if ( t2->type == VECTORFIELD_2D ) {
	if ( BAL_AllocTransformation( res, resType, &(t2->vx) ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (1))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      else if ( ref != (bal_image *)NULL ) {
	if ( BAL_AllocTransformation( res, resType, ref ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (2))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      else if ( t1->type == VECTORFIELD_2D ) {
	if ( BAL_AllocTransformation( res, resType, &(t1->vx) ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (2D vector field case (3))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate transformation (2D vector field case, weird case)", proc );
      return( -1 );

    case VECTORFIELD_3D :
      if ( t2->type == VECTORFIELD_3D ) {
	if ( BAL_AllocTransformation( res, resType, &(t2->vx) ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (1))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      else if ( ref != (bal_image *)NULL ) {
	if ( BAL_AllocTransformation( res, resType, ref ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (2))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      else if ( t1->type == VECTORFIELD_3D ) {
	if ( BAL_AllocTransformation( res, resType, &(t1->vx) ) != 1 ) {
	  if ( _verbose_ )
	    fprintf( stderr, "%s: unable to allocate transformation (3D vector field case (3))", proc );
	  return( -1 );
	}
	return( 1 );
      }
      
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate transformation (3D vector field case, weird case)", proc );
      return( -1 );

   }

  if ( _verbose_ )
    fprintf( stderr, "%s: should not reach this ...\n", proc );
   return( -1 );
}





/*--------------------------------------------------
 *
 * TRANSFORMATION COMPOSITION  (static operations)
 *
 --------------------------------------------------*/



typedef struct {
  bal_transformation *t1;
  bal_transformation *t2;
  bal_transformation *res;
} _TransformationCompositionParam;





static int _2DVectorField2DVectorFieldCompositionVoxelUnit( void *parameter,
							    size_t first,
							    size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x2, y2;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in real coordinates
    */
    x2 = i + vx[k][j][i];
    y2 = j + vy[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = x2 + BAL_GetXYKvalue( &(t1->vx), x2, y2, k );
    y1 = y2 + BAL_GetXYKvalue( &(t1->vy), x2, y2, k );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;

  } 
  return( 1 );
}





static int _2DVectorField2DVectorFieldCompositionRealUnit( void *parameter,
							    size_t first,
							    size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x, y;
  double x2, y2;
  double x2_1, y2_1;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    x = t2->vx.vx * i;
      y = t2->vx.vy * j;
      
      /* (x2, y2, z2) is the transformed point (i,j,k) by t2
	 in real coordinates
      */
      x2 = x + vx[k][j][i];
      y2 = y + vy[k][j][i];
      
      /* (x2_1, y2_1, z2_1) is the same point in voxel coordinates
	 in t1
      */
      x2_1 = x2 / t1->vx.vx;
      y2_1 = y2 / t1->vy.vy;
      
      /* (x1, y1, z1) is the transformed point (x2_1, y2_1, z2_1) by t1
	 in real coordinates
      */
      x1 = x2 + BAL_GetXYKvalue( &(t1->vx), x2_1, y2_1, k );
      y1 = y2 + BAL_GetXYKvalue( &(t1->vy), x2_1, y2_1, k );
      
      /* the displacement vector is calculated by substracting
	 the point (in real coordinates)
      */
      resx[k][j][i] = x1 - x;
      resy[k][j][i] = y1 - y;

  } 
  return( 1 );
}





/* t_res = t1 o t2
 */
static int _2DVectorField2DVectorFieldComposition( bal_transformation *t1, 
						   bal_transformation *t2,
						   bal_transformation *res ) 
{
  char *proc = "_2DVectorField2DVectorFieldComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;

  if ( t1 == res ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }
  
  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  /* the dimension of res and t2 should be the same
     as well as the voxel size
   */

  if ( t1->type != VECTORFIELD_2D || t2->type != VECTORFIELD_2D || res->type != VECTORFIELD_2D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types\n", proc );
    return( -1 );
  }
  

  
  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);



  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case VOXEL_UNIT :

    if ( processChunks( &_2DVectorField2DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :
    
    if ( processChunks( &_2DVectorField2DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  }

  freeChunks( &chunks );

  /* res has the same characteristics than t2
   */
  res->vx.vx = t2->vx.vx;
  res->vx.vy = t2->vx.vy;
  res->vx.vz = t2->vx.vz;

  res->vy.vx = t2->vy.vx;
  res->vy.vy = t2->vy.vy;
  res->vy.vz = t2->vy.vz;

  res->transformation_unit = t2->transformation_unit;

  return( 1 );
}





static int _3DVectorField3DVectorFieldCompositionVoxelUnit( void *parameter,
							    size_t first,
							    size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;
  float ***vz = (float ***)(t2->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x2, y2, z2;
  double x1, y1, z1;
  
  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in voxel coordinates
    */
    x2 = i + vx[k][j][i];
    y2 = j + vy[k][j][i];
    z2 = k + vz[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = x2 + BAL_GetXYZvalue( &(t1->vx), x2, y2, z2 );
    y1 = y2 + BAL_GetXYZvalue( &(t1->vy), x2, y2, z2 );
    z1 = z2 + BAL_GetXYZvalue( &(t1->vz), x2, y2, z2 );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;
    resz[k][j][i] = z1 - (double)k;
    
  }
  return( 1 );
}





static int _3DVectorField3DVectorFieldCompositionRealUnit( void *parameter,
							   size_t first,
							   size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;
  float ***vz = (float ***)(t2->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x, y, z;
  double x2, y2, z2;
  double x2_1, y2_1, z2_1;
  double x1, y1, z1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    x = t2->vx.vx * i;
    y = t2->vx.vy * j;
    z = t2->vx.vz * k;
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in real coordinates
    */
    x2 = x + vx[k][j][i];
    y2 = y + vy[k][j][i];
    z2 = z + vz[k][j][i];
    
    /* (x2_1, y2_1, z2_1) is the same point in voxel coordinates
       in t1
    */
    x2_1 = x2 / t1->vx.vx;
    y2_1 = y2 / t1->vy.vy;
    z2_1 = z2 / t1->vz.vz;
    
    /* (x1, y1, z1) is the transformed point (x2_1, y2_1, z2_1) by t1
       in real coordinates
    */
    x1 = x2 + BAL_GetXYZvalue( &(t1->vx), x2_1, y2_1, z2_1 );
    y1 = y2 + BAL_GetXYZvalue( &(t1->vy), x2_1, y2_1, z2_1 );
    z1 = z2 + BAL_GetXYZvalue( &(t1->vz), x2_1, y2_1, z2_1 );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - x;
    resy[k][j][i] = y1 - y;
    resz[k][j][i] = z1 - z;
    
  }
  return( 1 );
}





/* t_res = t1 o t2
 */
static int _3DVectorField3DVectorFieldComposition( bal_transformation *t1, 
						   bal_transformation *t2,
						   bal_transformation *res ) 
{
  char *proc = "_3DVectorField3DVectorFieldComposition";
  
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;
  
  if ( t1 == res ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }

  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  /* the dimension of res and t2 should be the same
     as well as the voxel size
   */

  if ( t1->type != VECTORFIELD_3D || t2->type != VECTORFIELD_3D || res->type != VECTORFIELD_3D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types\n", proc );
    return( -1 );
  }

  if ( t1->vx.nplanes <=1 || t2->vx.nplanes <=1 || res->vx.nplanes <=1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    return( -1 );
  }


  
  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);



  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case VOXEL_UNIT :

    if ( processChunks( &_3DVectorField3DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :

    if ( processChunks( &_3DVectorField3DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    
    break;

  }

  freeChunks( &chunks );

  /* res has the same characteristics than t2
   */
  res->vx.vx = t2->vx.vx;
  res->vx.vy = t2->vx.vy;
  res->vx.vz = t2->vx.vz;

  res->vy.vx = t2->vy.vx;
  res->vy.vy = t2->vy.vy;
  res->vy.vz = t2->vy.vz;

  res->vz.vx = t2->vz.vx;
  res->vz.vy = t2->vz.vy;
  res->vz.vz = t2->vz.vz;

  res->transformation_unit = t2->transformation_unit;

  return( 1 );
}





static int _2DVectorFieldMatrixCompositionVoxelUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;
  
  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  
  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x2, y2;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in voxel coordinates
    */
    x2 = t2->mat.m[ 0] * i + t2->mat.m[ 1] * j + t2->mat.m[ 3];
    y2 = t2->mat.m[ 4] * i + t2->mat.m[ 5] * j + t2->mat.m[ 7];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = x2 + BAL_GetXYKvalue( &(t1->vx), x2, y2, k );
    y1 = y2 + BAL_GetXYKvalue( &(t1->vy), x2, y2, k );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;
    
  }
  return( 1 );
}





static int _2DVectorFieldMatrixCompositionRealUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;
  
  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  
  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x, y;
  double x2, y2;
  double x2_1, y2_1;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    x = res->vx.vx * i;
    y = res->vx.vy * j;
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in real coordinates
    */
    x2 = t2->mat.m[ 0] * x + t2->mat.m[ 1] * y + t2->mat.m[ 3];
    y2 = t2->mat.m[ 4] * x + t2->mat.m[ 5] * y + t2->mat.m[ 7];
    
    /* (x2_1, y2_1, z2_1) is the same point in voxel coordinates
       in t1
    */
    x2_1 = x2 / t1->vx.vx;
    y2_1 = y2 / t1->vy.vy;
    
    /* (x1, y1, z1) is the transformed point (x2_1, y2_1, z2_1) by t1
       in real coordinates
    */
    x1 = x2 + BAL_GetXYKvalue( &(t1->vx), x2_1, y2_1, k );
    y1 = y2 + BAL_GetXYKvalue( &(t1->vy), x2_1, y2_1, k );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - x;
    resy[k][j][i] = y1 - y;

  }
  return( 1 );
}





/* t_res = t1 o t2
 */
static int _2DVectorFieldMatrixComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					   bal_transformation *res ) 
{
  char *proc = "_2DVectorFieldMatrixComposition";

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;

  if ( t1 == res ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }

  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( t1->type != VECTORFIELD_2D || res->type != VECTORFIELD_2D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 1st input and/or result\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( t2 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 2nd input\n", proc );
    return( -1 );
  }


 
  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);


  
  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :
    
    if ( processChunks( &_2DVectorFieldMatrixCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :
    
    if ( processChunks( &_2DVectorFieldMatrixCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    
    break;
    
  }

  freeChunks( &chunks );

  res->transformation_unit = t1->transformation_unit;

  return( 1 );
}





static int _3DVectorFieldMatrixCompositionVoxelUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x2, y2, z2;
  double x1, y1, z1;
  
  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in voxel coordinates
    */
    x2 = t2->mat.m[ 0] * i + t2->mat.m[ 1] * j + t2->mat.m[ 2] * k + t2->mat.m[ 3];
    y2 = t2->mat.m[ 4] * i + t2->mat.m[ 5] * j + t2->mat.m[ 6] * k + t2->mat.m[ 7];
    z2 = t2->mat.m[ 8] * i + t2->mat.m[ 9] * j + t2->mat.m[10] * k + t2->mat.m[11];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = x2 + BAL_GetXYZvalue( &(t1->vx), x2, y2, z2 );
    y1 = y2 + BAL_GetXYZvalue( &(t1->vy), x2, y2, z2 );
    z1 = z2 + BAL_GetXYZvalue( &(t1->vz), x2, y2, z2 );
    
    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;
    resz[k][j][i] = z1 - (double)k;

  }
  return( 1 );
}





static int _3DVectorFieldMatrixCompositionRealUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x, y, z;
  double x2, y2, z2;
  double x2_1, y2_1, z2_1;
  double x1, y1, z1;
  
  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    x = res->vx.vx * i;
    y = res->vx.vy * j;
    z = res->vx.vz * k;
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in real coordinates
    */
    x2 = t2->mat.m[ 0] * x + t2->mat.m[ 1] * y + t2->mat.m[ 2] * z + t2->mat.m[ 3];
    y2 = t2->mat.m[ 4] * x + t2->mat.m[ 5] * y + t2->mat.m[ 6] * z + t2->mat.m[ 7];
    z2 = t2->mat.m[ 8] * x + t2->mat.m[ 9] * y + t2->mat.m[10] * z + t2->mat.m[11];
    
    /* (x2_1, y2_1, z2_1) is the same point in voxel coordinates
       in t1
    */
    x2_1 = x2 / t1->vx.vx;
    y2_1 = y2 / t1->vy.vy;
    z2_1 = z2 / t1->vz.vz;
    
    /* (x1, y1, z1) is the transformed point (x2_1, y2_1, z2_1) by t1
       in real coordinates
    */
    x1 = x2 + BAL_GetXYZvalue( &(t1->vx), x2_1, y2_1, z2_1 );
    y1 = y2 + BAL_GetXYZvalue( &(t1->vy), x2_1, y2_1, z2_1 );
    z1 = z2 + BAL_GetXYZvalue( &(t1->vz), x2_1, y2_1, z2_1 );
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - x;
    resy[k][j][i] = y1 - y;
    resz[k][j][i] = z1 - z;
    
  }
  return( 1 );
}





/* t_res = t1 o t2
 */
static int _3DVectorFieldMatrixComposition( bal_transformation *t1, 
					    bal_transformation *t2,
					    bal_transformation *res ) 
{
  char *proc = "_3DVectorFieldMatrixComposition";
  
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;
  
  if ( t1 == res ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: output can not be equal to 1st input\n", proc );
    return( -1 );
  }

  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( t1->type != VECTORFIELD_3D || res->type != VECTORFIELD_3D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 1st input and/or result\n", proc );
    return( -1 );
  }

  if ( t1->vx.nplanes <=1 || res->vx.nplanes <=1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( t2 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 2nd input\n", proc );
    return( -1 );
  }


 
  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);



  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :
    
    if ( processChunks( &_3DVectorFieldMatrixCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :
    
    if ( processChunks( &_3DVectorFieldMatrixCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;
    
  }

  freeChunks( &chunks );

  res->transformation_unit = t1->transformation_unit;

  return( 1 );
}





static int _Matrix2DVectorFieldCompositionVoxelUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;
  
  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x2, y2;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in voxel coordinates
    */
    x2 = i + vx[k][j][i];
    y2 = j + vy[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = t1->mat.m[ 0] * x2 + t1->mat.m[ 1] * y2 + t1->mat.m[ 3];
    y1 = t1->mat.m[ 4] * x2 + t1->mat.m[ 5] * y2 + t1->mat.m[ 7];
    
    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;

  }
  return( 1 ); 
}





static int _Matrix2DVectorFieldCompositionRealUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;
  
  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;
  
  double x, y;
  double x2, y2;
  double x1, y1;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    x = t2->vx.vx * i;
    y = t2->vx.vy * j;
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
	 in real coordinates
    */
    x2 = x + vx[k][j][i];
    y2 = y + vy[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in real coordinates
    */
    x1 = t1->mat.m[ 0] * x2 + t1->mat.m[ 1] * y2 + t1->mat.m[ 3];
    y1 = t1->mat.m[ 4] * x2 + t1->mat.m[ 5] * y2 + t1->mat.m[ 7];
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - x;
    resy[k][j][i] = y1 - y;

  }
  return( 1 ); 
}





/* t_res = t1 o t2
 */
static int _Matrix2DVectorFieldComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					   bal_transformation *res ) 
{
  char *proc = "_Matrix2DVectorFieldComposition";
  
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;
  
  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( t2->type != VECTORFIELD_2D || res->type != VECTORFIELD_2D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 2nd input and/or result\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( t1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 1st input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);



  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :
    
    if ( processChunks( &_Matrix2DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
    
    break;

  case  REAL_UNIT :
    
    if ( processChunks( &_Matrix2DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;
    
  }

  freeChunks( &chunks );

  /* res has the same characteristics than t2
   */
  res->vx.vx = t2->vx.vx;
  res->vx.vy = t2->vx.vy;
  res->vx.vz = t2->vx.vz;

  res->vy.vx = t2->vy.vx;
  res->vy.vy = t2->vy.vy;
  res->vy.vz = t2->vy.vz;

  res->transformation_unit = t2->transformation_unit;

  return( 1 );
}





static int _Matrix3DVectorFieldCompositionVoxelUnit( void *parameter,
						     size_t first,
						     size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;
  float ***vz = (float ***)(t2->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x2, y2, z2;
  double x1, y1, z1;
  
  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {

    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in voxel coordinates
    */
    x2 = i + vx[k][j][i];
    y2 = j + vy[k][j][i];
    z2 = k + vz[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in voxel coordinates
    */
    x1 = t1->mat.m[ 0] * x2 + t1->mat.m[ 1] * y2 + t1->mat.m[ 2] * z2 + t1->mat.m[ 3];
    y1 = t1->mat.m[ 4] * x2 + t1->mat.m[ 5] * y2 + t1->mat.m[ 6] * z2 + t1->mat.m[ 7];
    z1 = t1->mat.m[ 8] * x2 + t1->mat.m[ 9] * y2 + t1->mat.m[10] * z2 + t1->mat.m[11];
    
    /* the displacement vector is calculated by substracting
       the point (in voxel coordinates)
    */
    resx[k][j][i] = x1 - (double)i;
    resy[k][j][i] = y1 - (double)j;
    resz[k][j][i] = z1 - (double)k;

  }
  return( 1 );
}





static int _Matrix3DVectorFieldCompositionRealUnit( void *parameter,
						    size_t first,
						    size_t last )
{
  _TransformationCompositionParam *p = (_TransformationCompositionParam *)parameter;

  bal_transformation *t1 = p->t1;
  bal_transformation *t2 = p->t2;
  bal_transformation *res = p->res;
  float ***resx = (float ***)(res->vx).array;
  float ***resy = (float ***)(res->vy).array;
  float ***resz = (float ***)(res->vz).array;
  float ***vx = (float ***)(t2->vx).array;
  float ***vy = (float ***)(t2->vy).array;
  float ***vz = (float ***)(t2->vz).array;

  size_t ncols = (res->vx).ncols; /* dimx */
  size_t nrows = (res->vx).nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  double x, y, z;
  double x2, y2, z2;
  double x1, y1, z1;
  
  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    
    x = t2->vx.vx * i;
    y = t2->vx.vy * j;
    z = t2->vx.vz * k;
    
    /* (x2, y2, z2) is the transformed point (i,j,k) by t2
       in real coordinates
    */
    x2 = x + vx[k][j][i];
    y2 = y + vy[k][j][i];
    z2 = z + vz[k][j][i];
    
    /* (x1, y1, z1) is the transformed point (x2, y2, z2) by t1
       in real coordinates
    */
    x1 = t1->mat.m[ 0] * x2 + t1->mat.m[ 1] * y2 + t1->mat.m[ 2] * z2 + t1->mat.m[ 3];
    y1 = t1->mat.m[ 4] * x2 + t1->mat.m[ 5] * y2 + t1->mat.m[ 6] * z2 + t1->mat.m[ 7];
    z1 = t1->mat.m[ 8] * x2 + t1->mat.m[ 9] * y2 + t1->mat.m[10] * z2 + t1->mat.m[11];
    
    /* the displacement vector is calculated by substracting
       the point (in real coordinates)
    */
    resx[k][j][i] = x1 - x;
    resy[k][j][i] = y1 - y;
    resz[k][j][i] = z1 - z;

  }
  return( 1 );
}





/* t_res = t1 o t2
 */
static int _Matrix3DVectorFieldComposition( bal_transformation *t1, 
					   bal_transformation *t2,
					   bal_transformation *res ) 
{
  char *proc = "_Matrix3DVectorFieldComposition";
   
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _TransformationCompositionParam p;
  
  if ( t1->transformation_unit != t2->transformation_unit ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations are not with the same unit\n", proc );
    return( -1 );
  }

  if ( t2->type != VECTORFIELD_3D || res->type != VECTORFIELD_3D ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation types for 2nd input and/or result\n", proc );
    return( -1 );
  }

  if ( t2->vx.nplanes <=1 || res->vx.nplanes <=1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations seem to be 2D\n", proc );
    return( -1 );
  }

  if ( BAL_IsTransformationLinear( t1 ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: bad transformation for 1st input\n", proc );
    return( -1 );
  }



  /* preparing parallelism
   */
  first = 0;
  last = (res->vx).nplanes * (res->vx).nrows * (res->vx).ncols - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return( -1 );
  }
  p.t1 = t1;
  p.t2 = t2;
  p.res = res;
  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);



  switch ( t1->transformation_unit ) {

  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such unit not handled yet\n", proc );
    return( -1 );

  case  VOXEL_UNIT :
    
    if ( processChunks( &_Matrix3DVectorFieldCompositionVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute voxel unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case  REAL_UNIT :
    
    if ( processChunks( &_Matrix3DVectorFieldCompositionRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute real unit composition\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;
    
  }

  freeChunks( &chunks );

  /* res has the same characteristics than t2
   */
  res->vx.vx = t2->vx.vx;
  res->vx.vy = t2->vx.vy;
  res->vx.vz = t2->vx.vz;

  res->vy.vx = t2->vy.vx;
  res->vy.vy = t2->vy.vy;
  res->vy.vz = t2->vy.vz;

  res->vz.vx = t2->vz.vx;
  res->vz.vy = t2->vz.vy;
  res->vz.vz = t2->vz.vz;

  res->transformation_unit = t2->transformation_unit;

  return( 1 );
}










/*--------------------------------------------------
 *
 * transformation use
 *
 --------------------------------------------------*/

int BAL_TransformPoint( bal_doublePoint *thePt, bal_doublePoint *resPt, bal_transformation *theTr )
{
  char *proc = "BAL_TransformPoint";
  double *mat;
  bal_doublePoint tmpPt;


  switch ( theTr->type ) {

  default : 
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    mat = theTr->mat.m;
    tmpPt.x = mat[ 0] * thePt->x + mat[ 1] * thePt->y + mat[ 2] * thePt->z + mat[ 3];
    tmpPt.y = mat[ 4] * thePt->x + mat[ 5] * thePt->y + mat[ 6] * thePt->z + mat[ 7];
    tmpPt.z = mat[ 8] * thePt->x + mat[ 9] * thePt->y + mat[10] * thePt->z + mat[11];

    resPt->x = tmpPt.x;
    resPt->y = tmpPt.y;
    resPt->z = tmpPt.z;

    break;

  }

  return( 1 );
}









/*--------------------------------------------------
 *
 * IMAGE RESAMPLING 
 *
 --------------------------------------------------*/

/* Resample 'image' into the geometry of 'resim'
   theTr is then the transformation that goes from 'resim' to 'image'
*/

int BAL_ResampleImage( bal_image *image, bal_image *resim, 
		       bal_transformation *theTr,
		       enumTransformationInterpolation interpolation )
{
  char *proc = "BAL_ResampleImage";
  bal_transformation voxTr;
  enumUnitTransfo unit = theTr->transformation_unit;

  BAL_InitTransformation( &voxTr );

  switch ( theTr->type ) {

  default : 
    if ( _verbose_ ) 
      fprintf( stderr, "%s: such a transformation type is not implemented yet\n", proc );
    return( -1 );

  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    
    if ( BAL_AllocTransformation( &voxTr, theTr->type, (bal_image *)NULL ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error when allocating voxel transformation matrix\n", proc );
      return( -1 );
    }
    if ( BAL_ChangeTransformationToVoxelUnit( image, resim, theTr, &voxTr ) != 1 ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: error when translating transformation matrix from real to voxel\n", proc );
      return( -1 );
    }

    switch ( interpolation ) {
    default :
      if ( _verbose_ ) 
	fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
      BAL_FreeTransformation( &voxTr );
      return( -1 );
    case NEAREST :
      if ( BAL_Reech3DNearest4x4( image, resim, &voxTr ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
	BAL_FreeTransformation( &voxTr );
	return( -1 );
      }
      break;    
    case LINEAR :
      if ( BAL_Reech3DTriLin4x4( image, resim, &voxTr ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: error when resampling image with a linear transformation\n", proc );
	BAL_FreeTransformation( &voxTr );
	return( -1 );
      }
      break;
    }

    BAL_FreeTransformation( &voxTr );
    break;

  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    
    if ( resim->ncols != theTr->vx.ncols 
       || resim->nrows != theTr->vx.nrows 
       || resim->nplanes != theTr->vx.nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: output image should have the same dimensions than the transformation\n", proc );
      return( -1 );
    }

    if ( fabs( resim->vx - theTr->vx.vx) > theTr->vx.vx / 1000.0 
	 || fabs( resim->vy - theTr->vx.vy) > theTr->vx.vy / 1000.0 
	 || fabs( resim->vz - theTr->vx.vz) > theTr->vx.vz / 1000.0 ) {
	   if ( _verbose_ ) 
	fprintf( stderr, "%s: output image should have the same voxel sizes than the transformation\n", proc );
      return( -1 );
    }

    if ( unit != VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToVoxelUnit( image, resim, theTr, theTr ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: error when translating transformation vector field from real to voxel\n", proc );
	return( -1 );
      }
    }



    switch ( theTr->type ) {
    default :
      break;
    case VECTORFIELD_2D :
      switch ( interpolation ) {
      default :
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
	return( -1 );
      case NEAREST :
	if ( BAL_Reech3DNearest2DVectorField( image, resim, theTr ) != 1 ) {
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
	  return( -1 );
	}
	break;
      case LINEAR :
	if ( BAL_Reech3DTriLin2DVectorField( image, resim, theTr ) != 1 ) {
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
	  return( -1 );
	}
	break;
      }
      break;
    case VECTORFIELD_3D :
      switch ( interpolation ) {
      default :
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: such interpolation type not handled yet\n", proc );
	return( -1 );
      case NEAREST :
	if ( BAL_Reech3DNearest3DVectorField( image, resim, theTr ) != 1 ) {
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
	  return( -1 );
	}
	break;
      case LINEAR :
	if ( BAL_Reech3DTriLin3DVectorField( image, resim, theTr ) != 1 ) {
	  if ( _verbose_ ) 
	    fprintf( stderr, "%s: error when resampling image with a vector field transformation\n", proc );
	  return( -1 );
	}
	break;
      }
      break;
    }

    if ( unit != VOXEL_UNIT ) {
      if ( BAL_ChangeTransformationToRealUnit( image, resim, theTr, theTr ) != 1 ) {
	if ( _verbose_ ) 
	  fprintf( stderr, "%s: error when translating transformation vector field from voxel to real\n", proc );
	return( -1 );
      }
    }
    
    break;

  }
  return( 1 );
}





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Reech3DTriLin4x4( bal_image *theIm, bal_image *resIm,
				bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin4x4";
  int theDim[3];
  int resDim[3];
  double *m = (double*)theTr->mat.m;

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case SSHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case USHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case FLOAT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DTriLin4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DTriLin4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DTriLin2DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin2DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols 
       || resIm->nrows != theTr->vx.nrows 
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;  

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)NULL;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech2DTriLinVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim, 
				     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech2DTriLinVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech2DTriLinVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech2DTriLinVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DTriLin3DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DTriLin3DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols 
       || resIm->nrows != theTr->vx.nrows 
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;  

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)(theTr->vz.data);

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech3DTriLinVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,  
				     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech3DTriLinVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech3DTriLinVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech3DTriLinVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





/* here, the transformation is assumed to be in voxel coordinates
   ie from image frame to image frame
 */
static int BAL_Reech3DNearest4x4( bal_image *theIm, bal_image *resIm,
				bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest4x4";
  int theDim[3];
  int resDim[3];
  double *m = (double*)theTr->mat.m;

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_u8( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case SSHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_s16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case USHORT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_u16( theIm->data, theDim, resIm->data, resDim, m );
    break;
  case FLOAT :
    if ( theIm->nplanes == 1 && resIm->nplanes == 1 )
      Reech2DNearest4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    else
      Reech3DNearest4x4_r32( theIm->data, theDim, resIm->data, resDim, m );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DNearest2DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest2DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols 
       || resIm->nrows != theTr->vx.nrows 
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;  

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)NULL;

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech2DNearestVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim, 
				     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech2DNearestVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech2DNearestVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech2DNearestVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}





static int BAL_Reech3DNearest3DVectorField( bal_image *theIm, bal_image *resIm,
					 bal_transformation *theTr )
{
  char *proc = "BAL_Reech3DNearest3DVectorField";
  int theDim[3];
  int resDim[3];
  int vecDim[3];
  float *vecBuf[3];

  if ( theIm->type != resIm->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: images have different types\n", proc );
    return( -1 );
  }

  if ( resIm->ncols != theTr->vx.ncols 
       || resIm->nrows != theTr->vx.nrows 
       || resIm->nplanes != theTr->vx.nplanes ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: result and X component images have different dimensions\n", proc );
      fprintf( stderr, "\t result image      is %lu x %lu x %lu\n", resIm->ncols, resIm->nrows, resIm->nplanes );
      fprintf( stderr, "\t X component image is %lu x %lu x %lu\n", theTr->vx.ncols, theTr->vx.nrows, theTr->vx.nplanes );
    }
    return( -1 );
  }

  theDim[0] = theIm->ncols;
  theDim[1] = theIm->nrows;
  theDim[2] = theIm->nplanes;
  
  resDim[0] = resIm->ncols;
  resDim[1] = resIm->nrows;
  resDim[2] = resIm->nplanes;  

  vecDim[0] = theTr->vx.ncols;
  vecDim[1] = theTr->vx.nrows;
  vecDim[2] = theTr->vx.nplanes;

  vecBuf[0] = (float*)(theTr->vx.data);
  vecBuf[1] = (float*)(theTr->vy.data);
  vecBuf[2] = (float*)(theTr->vz.data);

  switch( theIm->type ) {
  default :
    if ( _verbose_ )
      fprintf(stderr, "%s: such type not handled yet\n", proc );
    return( -1 );
  case UCHAR :
    Reech3DNearestVectorField_r32_u8( theIm->data, theDim, resIm->data, resDim,  
				     vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case SSHORT :
    Reech3DNearestVectorField_r32_s16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case USHORT :
    Reech3DNearestVectorField_r32_u16( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  case FLOAT :
    Reech3DNearestVectorField_r32_r32( theIm->data, theDim, resIm->data, resDim,  
				      vecBuf, vecDim, (double*)NULL, (double*)NULL );
    break;
  }
  return( 1 );
}










/*--------------------------------------------------
 *
 * transformation conversion
 *
 --------------------------------------------------*/



typedef struct {
  bal_image *image;
  bal_image *resim; 
  bal_transformation *theTrsf;
  bal_transformation *resTrsf;
} _UnitTransformationChangeParam;



/* 
   theTr is the transformation that goes from 'resim' to 'image'

   - linear case: 
     we have M_image = T * M_resim, with M_image and M_resim in voxels
     we go into real coordinates with
     H_image * M_image = H_image * T * M_resim
                       = H_image * T * H^(-1)_resim * H_resim * M_resim
     thus H_image * T * H^(-1)_resim is the transformation in real coordinates

   - vector field case:
     we have M_Image = M_resim + V(M_resim)
     we go into real coordinates with
     H_image * M_image = H_image * M_resim + H_image * V(M_resim)
                       = H_resim * M_resim + (H_image - H_resim) * M_resim + H_image * V(M_resim)
     thus, the displacement vector in real coordinates is
     (H_image - H_resim) * M_resim + H_image * V(M_resim)
*/



static int _Change2DVectorFieldToRealUnit( void *parameter,
					   size_t first,
					   size_t last )
{
  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *image = p->image;
  bal_image *resim = p->resim;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);

  double vx_from = resim->vx;
  double vy_from = resim->vy;
  double vx_to = image->vx;
  double vy_to = image->vy;
  double dx = vx_to - vx_from;
  double dy = vy_to - vy_from;
  
  size_t ncols = theTrsf->vx.ncols; /* dimx */
  size_t nrows = theTrsf->vx.nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    resVecX[k][j][i] = dx * i + vx_to * theVecX[k][j][i]; 
    resVecY[k][j][i] = dy * j + vy_to * theVecY[k][j][i]; 
  }
  return( 1 );
}





static int _Change3DVectorFieldToRealUnit( void *parameter,
					   size_t first,
					   size_t last )
{
  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;

  bal_image *image = p->image;
  bal_image *resim = p->resim;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***theVecZ = (float***)(theTrsf->vz.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float ***resVecZ = (float***)(resTrsf->vz.array);

  double vx_from = resim->vx;
  double vy_from = resim->vy;
  double vz_from = resim->vz;
  double vx_to = image->vx;
  double vy_to = image->vy;
  double vz_to = image->vz;
  double dx = vx_to - vx_from;
  double dy = vy_to - vy_from;
  double dz = vz_to - vz_from;
  
  size_t ncols = theTrsf->vx.ncols; /* dimx */
  size_t nrows = theTrsf->vx.nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    resVecX[k][j][i] = dx * i + vx_to * theVecX[k][j][i]; 
    resVecY[k][j][i] = dy * j + vy_to * theVecY[k][j][i]; 
    resVecZ[k][j][i] = dz * k + vz_to * theVecZ[k][j][i]; 
  }
  return( 1 );
}




int BAL_ChangeTransformationToRealUnit( bal_image *image,
					bal_image *resim, 
					bal_transformation *theTrsf,
					bal_transformation *resTrsf )
{
  char *proc = "BAL_ChangeTransformationToRealUnit";
  double vx_from;
  double vy_from;
  double vz_from;
  double vx_to;
  double vy_to;
  double vz_to;

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _UnitTransformationChangeParam p;
  
  if ( image == (bal_image*)NULL || resim == (bal_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL template images\n", proc );
    return( -1 );
  }
  
  vx_from = resim->vx;
  vy_from = resim->vy;
  vz_from = resim->vz;
  vx_to = image->vx;
  vy_to = image->vy;
  vz_to = image->vz;

  if ( theTrsf->type != resTrsf->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations have different types\n", proc );
    return( -1 );
  }
  
  /* nothing to do, but copy the transformation
   */
  if ( theTrsf->transformation_unit == REAL_UNIT ) {
    if ( BAL_CopyTransformation( theTrsf, resTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to copy transformation\n", proc );
      return( -1 );
    }
    return( 1 );
  }



  /* preparing parallelism
   */
  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    first = 0;
    last = theTrsf->vx.nplanes * theTrsf->vx.nrows * theTrsf->vx.ncols - 1;
    initChunks( &chunks );
    if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
    }
    p.image = image;
    p.resim = resim;
    p.theTrsf = theTrsf;
    p.resTrsf = resTrsf;
    for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
      chunks.data[i].parameters = (void*)(&p);
    break;
  }

  

  switch ( theTrsf->type ) {
  default :
    fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    resTrsf->mat.m[ 0] = theTrsf->mat.m[ 0] * vx_to / vx_from;
    resTrsf->mat.m[ 1] = theTrsf->mat.m[ 1] * vx_to / vy_from;
    resTrsf->mat.m[ 2] = theTrsf->mat.m[ 2] * vx_to / vz_from;
    resTrsf->mat.m[ 3] = theTrsf->mat.m[ 3] * vx_to;
    
    resTrsf->mat.m[ 4] = theTrsf->mat.m[ 4] * vy_to / vx_from;
    resTrsf->mat.m[ 5] = theTrsf->mat.m[ 5] * vy_to / vy_from;
    resTrsf->mat.m[ 6] = theTrsf->mat.m[ 6] * vy_to / vz_from;
    resTrsf->mat.m[ 7] = theTrsf->mat.m[ 7] * vy_to;
    
    resTrsf->mat.m[ 8] = theTrsf->mat.m[ 8] * vz_to / vx_from;
    resTrsf->mat.m[ 9] = theTrsf->mat.m[ 9] * vz_to / vy_from;
    resTrsf->mat.m[10] = theTrsf->mat.m[10] * vz_to / vz_from;
    resTrsf->mat.m[11] = theTrsf->mat.m[11] * vz_to;

    break;

  case VECTORFIELD_2D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols 
	 || resTrsf->vx.nrows != theTrsf->vx.nrows 
	 || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
		 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
		 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      return( -1 );
    }

    if ( processChunks( &_Change2DVectorFieldToRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to change 2D vector field into real unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }

    break;

  case VECTORFIELD_3D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols 
	 || resTrsf->vx.nrows != theTrsf->vx.nrows 
	 || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
		 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
		 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      return( -1 );
    }

    if ( processChunks( &_Change3DVectorFieldToRealUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to change 3D vector field into real unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
      
    break;

  }
  
  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    freeChunks( &chunks );
    break;
  }
  
  resTrsf->transformation_unit = REAL_UNIT;

  return( 1 );
}






/* 
   theTr is the transformation that goes from 'resim' to 'image'

   - linear case: 
     we have H_image * M_image = T * H_resim * M_resim, 
        with M_image and M_resim in voxel coordinates
     we go into voxel coordinates with
     M_image = H^(-1)_image * T * H_resim * M_resim, 
     thus H^(-1)_image * T * H_resim is the transformation in voxel coordinates

   - vector field case:
     we have H_image * M_Image = H_resim * M_resim + V(M_resim)
     we go into voxel coordinates with
     M_image = H^(-1)_image * H_resim * M_resim + H^(-1)_image * V(M_resim)
             = M_resim + ( H^(-1)_image * H_resim - Id ) * M_resim + H^(-1)_image * V(M_resim)
     thus, the displacement vector in real coordinates is
     ( H^(-1)_image * H_resim - Id ) * M_resim + H^(-1)_image * V(M_resim)
*/



static int _Change2DVectorFieldToVoxelUnit( void *parameter,
					   size_t first,
					   size_t last )
{
  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;
  
  bal_image *image = p->image;
  bal_image *resim = p->resim;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);

  double vx_from = resim->vx;
  double vy_from = resim->vy;
  double vx_to = image->vx;
  double vy_to = image->vy;
  double px = vx_from / vx_to - 1.0;
  double py = vy_from / vy_to - 1.0;

  size_t ncols = theTrsf->vx.ncols; /* dimx */
  size_t nrows = theTrsf->vx.nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    resVecX[k][j][i] = px * i + theVecX[k][j][i] / vx_to; 
    resVecY[k][j][i] = py * j + theVecY[k][j][i] / vy_to; 
  }
  return( 1 );
}





static int _Change3DVectorFieldToVoxelUnit( void *parameter,
					   size_t first,
					   size_t last )
{
  _UnitTransformationChangeParam *p = (_UnitTransformationChangeParam *)parameter;
  
  bal_image *image = p->image;
  bal_image *resim = p->resim;
  bal_transformation *theTrsf = p->theTrsf;
  bal_transformation *resTrsf = p->resTrsf;
  float ***theVecX = (float***)(theTrsf->vx.array);
  float ***theVecY = (float***)(theTrsf->vy.array);
  float ***theVecZ = (float***)(theTrsf->vz.array);
  float ***resVecX = (float***)(resTrsf->vx.array);
  float ***resVecY = (float***)(resTrsf->vy.array);
  float ***resVecZ = (float***)(resTrsf->vz.array);

  double vx_from = resim->vx;
  double vy_from = resim->vy;
  double vz_from = resim->vz;
  double vx_to = image->vx;
  double vy_to = image->vy;
  double vz_to = image->vz;
  double px = vx_from / vx_to - 1.0;
  double py = vy_from / vy_to - 1.0;
  double pz = vz_from / vz_to - 1.0;

  size_t ncols = theTrsf->vx.ncols; /* dimx */
  size_t nrows = theTrsf->vx.nrows; /* dimy */
  size_t dimxy = ncols*nrows;
  int i, j, k;
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / dimxy;
  j = jfirst = (first - kfirst*dimxy) / ncols;
  i = ifirst = (first - kfirst*dimxy - jfirst*ncols);

  klast = last / dimxy;
  jlast = (last - klast*dimxy) / ncols;
  ilast = (last - klast*dimxy - jlast*ncols);

  for ( ; k<=klast; k++, j=0 )
  for ( ; (j<nrows && k<klast) || (j<=jlast && k==klast); j++, i=0 )
  for ( ; (i<ncols && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++ ) {
    resVecX[k][j][i] = px * i + theVecX[k][j][i] / vx_to; 
    resVecY[k][j][i] = py * j + theVecY[k][j][i] / vy_to; 
    resVecZ[k][j][i] = pz * k + theVecZ[k][j][i] / vz_to; 
  }
  return( 1 );
}






/* theTrsf allows to resample 'image' into 'resim' thus goes
   from 'resim' to 'image'
*/
int BAL_ChangeTransformationToVoxelUnit( bal_image *image,
					 bal_image *resim, 
					 bal_transformation *theTrsf,
					 bal_transformation *resTrsf )
{
  char *proc = "BAL_ChangeTransformationToVoxelUnit";
  double vx_from;
  double vy_from;
  double vz_from;
  double vx_to;
  double vy_to;
  double vz_to;

  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _UnitTransformationChangeParam p;

  if ( image == (bal_image*)NULL || resim == (bal_image*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: NULL template images\n", proc );
    return( -1 );
  }
  
  vx_from = resim->vx;
  vy_from = resim->vy;
  vz_from = resim->vz;
  vx_to = image->vx;
  vy_to = image->vy;
  vz_to = image->vz;

  if ( theTrsf->type != resTrsf->type ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: transformations have different types\n", proc );
    return( -1 );
  }

  /* nothing to do, but copy the transformation
   */
  if ( theTrsf->transformation_unit == VOXEL_UNIT ) {
    if ( BAL_CopyTransformation( theTrsf, resTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to copy transformation\n", proc );
      return( -1 );
    }
    return( 1 );
  }



  /* preparing parallelism
   */
  switch ( theTrsf->type ) {
  default :
    break;
  case VECTORFIELD_2D :
  case VECTORFIELD_3D :
    first = 0;
    last = theTrsf->vx.nplanes * theTrsf->vx.nrows * theTrsf->vx.ncols - 1;
    initChunks( &chunks );
    if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compute chunks\n", proc );
      return( -1 );
    }
    p.image = image;
    p.resim = resim;
    p.theTrsf = theTrsf;
    p.resTrsf = resTrsf;
    for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
      chunks.data[i].parameters = (void*)(&p);
    break;
  }



  switch ( theTrsf->type ) {
  default :
    fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
    return( -1 );
    
  case TRANSLATION_2D :
  case TRANSLATION_3D :
  case TRANSLATION_SCALING_2D :
  case TRANSLATION_SCALING_3D :
  case RIGID_2D :
  case RIGID_3D :
  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :

    resTrsf->mat.m[ 0] = theTrsf->mat.m[ 0] * vx_from / vx_to;
    resTrsf->mat.m[ 1] = theTrsf->mat.m[ 1] * vy_from / vx_to;
    resTrsf->mat.m[ 2] = theTrsf->mat.m[ 2] * vz_from / vx_to;
    resTrsf->mat.m[ 3] = theTrsf->mat.m[ 3]           / vx_to;
    
    resTrsf->mat.m[ 4] = theTrsf->mat.m[ 4] * vx_from / vy_to;
    resTrsf->mat.m[ 5] = theTrsf->mat.m[ 5] * vy_from / vy_to;
    resTrsf->mat.m[ 6] = theTrsf->mat.m[ 6] * vz_from / vy_to;
    resTrsf->mat.m[ 7] = theTrsf->mat.m[ 7]           / vy_to;
    
    resTrsf->mat.m[ 8] = theTrsf->mat.m[ 8] * vx_from / vz_to;
    resTrsf->mat.m[ 9] = theTrsf->mat.m[ 9] * vy_from / vz_to;
    resTrsf->mat.m[10] = theTrsf->mat.m[10] * vz_from / vz_to;
    resTrsf->mat.m[11] = theTrsf->mat.m[11]           / vz_to;

    break;

  case VECTORFIELD_2D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols 
	 || resTrsf->vx.nrows != theTrsf->vx.nrows 
	 || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
		 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
		 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      return( -1 );
    }

     if ( processChunks( &_Change2DVectorFieldToVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to change 2D vector field into voxel unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
      
    break;

  case VECTORFIELD_3D :

    if ( resTrsf->vx.ncols != theTrsf->vx.ncols 
	 || resTrsf->vx.nrows != theTrsf->vx.nrows 
	 || resTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
      if ( _verbose_ ) 
	fprintf( stderr, "%s: X-components have different dimensions, input=(%lu,%lu,%lu) and output=(%lu,%lu,%lu)\n",
		 proc, theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes,
		 resTrsf->vx.ncols, resTrsf->vx.nrows, resTrsf->vx.nplanes );
      return( -1 );
    }

    if ( processChunks( &_Change3DVectorFieldToVoxelUnit, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to change3D vector field into voxel unit\n", proc );
      freeChunks( &chunks );
      return( -1 );
    }
      
    break;

  }

  resTrsf->transformation_unit = VOXEL_UNIT;

  return( 1 );

}








/*--------------------------------------------------
 *
 * TRANSFORMATION INVERSION 
 *
 --------------------------------------------------*/

int BAL_InverseTransformation( bal_transformation *theTrsf,  
			       bal_transformation *invTrsf )
{
  char *proc = "BAL_InverseTransformation";
  double *theMat = theTrsf->mat.m;
  double *invMat = invTrsf->mat.m;

  if ( theTrsf == invTrsf ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: output transformation has to be different from the input one\n", proc );
    }
    return( -1 );
  }

  

  if ( BAL_DoesTransformationExist( invTrsf ) != 1 ) {
    if ( _verbose_ ) {
      fprintf( stderr, "%s: inverse transformation does not exist or is not allocated\n", proc );
    }
    return( -1 );
  }
  


  switch ( theTrsf->type ) {
    
  default :

    if ( _verbose_ ) {
      fprintf( stderr, "%s: such transformation type not handled yet\n", proc );
      BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
    }
    return( -1 );
    
  case TRANSLATION_2D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 3] = (- theMat[ 3]);
    invMat[ 7] = (- theMat[ 7]);
    break;

  case TRANSLATION_3D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 3] = (- theMat[ 3]);
    invMat[ 7] = (- theMat[ 7]);
    invMat[11] = (- theMat[11]);
    break;

  case TRANSLATION_SCALING_2D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 0] = 1.0/theMat[ 0];
    invMat[ 5] = 1.0/theMat[ 5];
    invMat[ 3] = (- theMat[ 3])/theMat[ 0];
    invMat[ 7] = (- theMat[ 7])/theMat[ 5];
    break;

  case TRANSLATION_SCALING_3D :
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 0] = 1.0/theMat[ 0];
    invMat[ 5] = 1.0/theMat[ 5];
    invMat[10] = 1.0/theMat[10];
    invMat[ 3] = (- theMat[ 3])/theMat[ 0];
    invMat[ 7] = (- theMat[ 7])/theMat[ 5];
    invMat[11] = (- theMat[11])/theMat[10];
    break;

  case RIGID_2D :
    /* l'inverse d'une matrice de rotation R est sa transposee R^t
       la translation inverse est -R^t * T
     */
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 1] = theMat[ 4];
    invMat[ 4] = theMat[ 1];
    invMat[ 3] = ( - (theMat[ 0]*theMat[ 3] + theMat[ 4]*theMat[ 7]) );
    invMat[ 7] = ( - (theMat[ 1]*theMat[ 3] + theMat[ 5]*theMat[ 7]) );
    break;

  case RIGID_3D :
    /* l'inverse d'une matrice de rotation R est sa transposee R^t
       la translation inverse est -R^t * T
     */
    BAL_SetTransformationToIdentity( invTrsf );
    invMat[ 1] = theMat[ 4];
    invMat[ 2] = theMat[ 8];
    invMat[ 4] = theMat[ 1];
    invMat[ 6] = theMat[ 9];
    invMat[ 8] = theMat[ 2];
    invMat[ 9] = theMat[ 6];
    invMat[ 3] = ( - (theMat[ 0]*theMat[ 3] + theMat[ 4]*theMat[ 7] + theMat[ 8]*theMat[11]) );
    invMat[ 7] = ( - (theMat[ 1]*theMat[ 3] + theMat[ 5]*theMat[ 7] + theMat[ 9]*theMat[11]) );
    invMat[11] = ( - (theMat[ 2]*theMat[ 3] + theMat[ 6]*theMat[ 7] + theMat[10]*theMat[11]) );
    break;

  case SIMILITUDE_2D :
  case SIMILITUDE_3D :
  case AFFINE_2D :
  case AFFINE_3D :
    if ( InverseMat4x4( theMat, invMat ) != 4 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: input transformation is not invertible (linear case)\n", proc );
	BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }
    break;
    
  case VECTORFIELD_2D :
    if ( BAL_Inverse2DVectorField( theTrsf, invTrsf ) != 1 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: input transformation is not invertible (2D vector field)\n", proc );
	BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }
    break;

  case VECTORFIELD_3D :

    if ( BAL_Inverse3DVectorField( theTrsf, invTrsf ) != 1 ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: input transformation is not invertible (3D vector field)\n", proc );
	BAL_PrintTransformation( stderr, theTrsf, "input transformation" );
      }
      return( -1 );
    }
    break;
  
  }

  return( 1 );

}








/*--------------------------------------------------
 *
 * inversion of vector fields (static functions)
 *
 --------------------------------------------------*/

static void _get2DMatrixVector( double x, double y, int k, 
				float ***xdx, float ***xdy, float ***ydx, float ***ydy, 
				float ***vecx, float*** vecy,
				int *dim,
				double *mat, double *vec )
{
  int ix, iy;
  double dx, dy;
  double a00, a01, a10, a11;
  
  ix = (int)x;
  iy = (int)y;
  dx = x - ix;
  dy = y - iy;
	
  if ( ix < 0 ) {
    ix = 0;
    dx = 0; 
  }
  if ( iy < 0 ) { 
    iy = 0;
    dy = 0;
  }

  if ( ix >= dim[0] - 1 ) {
    ix = dim[0]-2;
    dx = 1;
  }
  if ( iy >= dim[1] - 1 ) {
    iy = dim[1]-2;
    dy = 1;
  }
  
  /* res = (1-dx) * ( (1-dy) * buf[k][iy][ix] 
	              + (dy) * buf[k][iy+1][ix] )
	   + (dx) * ( (1-dy) * buf[k][iy][ix+1] 
       	              + (dy) * buf[k][iy+1][ix+1] );
	 =   a00 * buf[k][iy][ix] 
	   + a10 * buf[k][iy+1][ix]
	   + a01 * buf[k][iy][ix+1] 
       	   + a11 * buf[k][iy+1][ix+1];
     a00 = (1-dx) * (1-dy);
     a10 = (1-dx) *    dy;
     a01 =    dx  * (1-dy);
     a11 =    dx  *    dy	   
  */
  a11 = dx * dy;
  a01 = dx - a11;
  a10 = dy - a11;
  a00 = 1 - dx - a10;

  mat[0] = a00 * xdx[k][iy][ix] + a10 * xdx[k][iy+1][ix] + a01 * xdx[k][iy][ix+1]  + a11 * xdx[k][iy+1][ix+1];
  mat[1] = a00 * xdy[k][iy][ix] + a10 * xdy[k][iy+1][ix] + a01 * xdy[k][iy][ix+1]  + a11 * xdy[k][iy+1][ix+1];
  mat[2] = a00 * ydx[k][iy][ix] + a10 * ydx[k][iy+1][ix] + a01 * ydx[k][iy][ix+1]  + a11 * ydx[k][iy+1][ix+1];
  mat[3] = a00 * ydy[k][iy][ix] + a10 * ydy[k][iy+1][ix] + a01 * ydy[k][iy][ix+1]  + a11 * ydy[k][iy+1][ix+1];

  vec[0] = a00 * vecx[k][iy][ix] + a10 * vecx[k][iy+1][ix] + a01 * vecx[k][iy][ix+1]  + a11 * vecx[k][iy+1][ix+1];
  vec[1] = a00 * vecy[k][iy][ix] + a10 * vecy[k][iy+1][ix] + a01 * vecy[k][iy][ix+1]  + a11 * vecy[k][iy+1][ix+1];
}
	



static void _get3DMatrixVector( double x, double y, double z, 
				float ***xdx, float ***xdy, float ***xdz, 
				float ***ydx, float ***ydy, float ***ydz,
				float ***zdx, float ***zdy, float ***zdz,
				float ***vecx, float*** vecy, float*** vecz,
				int *dim,
				double *mat, double *vec )
{
  int ix, iy, iz;
  double dx, dy, dz;
  double a000, a001, a010, a011, a100, a101, a110, a111;
  double dxdy, dxdz, dydz;

  ix = (int)x;
  iy = (int)y;
  iz = (int)z;
  dx = x - ix;
  dy = y - iy;
  dz = z - iz;
	
  if ( ix < 0 ) {
    ix = 0;
    dx = 0; 
  }
  if ( iy < 0 ) { 
    iy = 0;
    dy = 0;
  }
  if ( iz < 0 ) { 
    iz = 0;
    dz = 0;
  }

  if ( ix >= dim[0] - 1 ) {
    ix = dim[0]-2;
    dx = 1;
  }
  if ( iy >= dim[1] - 1 ) {
    iy = dim[1]-2;
    dy = 1;
  }
  if ( iz >= dim[2] - 1 ) {
    iz = dim[2]-2;
    dz = 1;
  }
  
  /* res = (1-dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix] 
				  + (dz)*buf[iz+1][iy][ix] ) 
		       + (dy) * ( (1-dz)*buf[iz][iy+1][ix] 
				  + (dz)*buf[iz+1][iy+1][ix] ) )
	    + (dx) * ( (1-dy) * ( (1-dz)*buf[iz][iy][ix+1] 
				  + (dz)*buf[iz+1][iy][ix+1] )
		       + (dy) * ( (1-dz)*buf[iz][iy+1][ix+1]
				  + (dz)*buf[iz+1][iy+1][ix+1] ) );

	 =   a000 * buf[iz][iy][ix] 
	   + a100 * buf[iz+1][iy][ix]
	   + a010 * buf[iz][iy+1][ix] 
	   + a110 * buf[iz+1][iy+1][ix]
	   + a001 * buf[iz][iy][ix+1] 
	   + a101 * buf[iz+1][iy][ix+1] 
	   + a011 * buf[iz][iy+1][ix+1] 
	   + a111 * buf[iz+1][iy+1][ix+1] 

     a000 = (1-dx) * (1-dy) * (1-dz);
     a100 = (1-dx) * (1-dy) *    dz;
     a010 = (1-dx) *    dy  * (1-dz);
     a110 = (1-dx) *    dy  *    dz;
     a001 =    dx  * (1-dy) * (1-dz);
     a101 =    dx  * (1-dy) *    dz;
     a011 =    dx  *    dy  * (1-dz);
     a111 =    dx  *    dy  *    dz;

  */
  
  dxdy = dx*dy;
  dxdz = dx*dz;
  dydz = dy*dz;

  a111 = dxdy * dz;
  a011 = dxdy - a111;
  a101 = dxdz - a111;
  a110 = dydz - a111;
  a001 = dx - dxdy - a101;
  a010 = dy - dydz - a011;
  a100 = dz - dxdz - a110;
  a000 = 1 - dy -dz + dydz - a001;


  mat[0] = a000 * xdx[iz][iy][ix] + a100 * xdx[iz+1][iy][ix] + a010 * xdx[iz][iy+1][ix]  + a110 * xdx[iz+1][iy+1][ix]
    + a001 * xdx[iz][iy][ix+1] + a101 * xdx[iz+1][iy][ix+1] + a011 * xdx[iz][iy+1][ix+1]  + a111 * xdx[iz+1][iy+1][ix+1];
  mat[1] = a000 * xdy[iz][iy][ix] + a100 * xdy[iz+1][iy][ix] + a010 * xdy[iz][iy+1][ix]  + a110 * xdy[iz+1][iy+1][ix]
    + a001 * xdy[iz][iy][ix+1] + a101 * xdy[iz+1][iy][ix+1] + a011 * xdy[iz][iy+1][ix+1]  + a111 * xdy[iz+1][iy+1][ix+1];
  mat[2] = a000 * xdz[iz][iy][ix] + a100 * xdz[iz+1][iy][ix] + a010 * xdz[iz][iy+1][ix]  + a110 * xdz[iz+1][iy+1][ix]
    + a001 * xdz[iz][iy][ix+1] + a101 * xdz[iz+1][iy][ix+1] + a011 * xdz[iz][iy+1][ix+1]  + a111 * xdz[iz+1][iy+1][ix+1];

  mat[3] = a000 * ydx[iz][iy][ix] + a100 * ydx[iz+1][iy][ix] + a010 * ydx[iz][iy+1][ix]  + a110 * ydx[iz+1][iy+1][ix]
    + a001 * ydx[iz][iy][ix+1] + a101 * ydx[iz+1][iy][ix+1] + a011 * ydx[iz][iy+1][ix+1]  + a111 * ydx[iz+1][iy+1][ix+1];
  mat[4] = a000 * ydy[iz][iy][ix] + a100 * ydy[iz+1][iy][ix] + a010 * ydy[iz][iy+1][ix]  + a110 * ydy[iz+1][iy+1][ix]
    + a001 * ydy[iz][iy][ix+1] + a101 * ydy[iz+1][iy][ix+1] + a011 * ydy[iz][iy+1][ix+1]  + a111 * ydy[iz+1][iy+1][ix+1];
  mat[5] = a000 * ydz[iz][iy][ix] + a100 * ydz[iz+1][iy][ix] + a010 * ydz[iz][iy+1][ix]  + a110 * ydz[iz+1][iy+1][ix]
    + a001 * ydz[iz][iy][ix+1] + a101 * ydz[iz+1][iy][ix+1] + a011 * ydz[iz][iy+1][ix+1]  + a111 * ydz[iz+1][iy+1][ix+1];

  mat[6] = a000 * zdx[iz][iy][ix] + a100 * zdx[iz+1][iy][ix] + a010 * zdx[iz][iy+1][ix]  + a110 * zdx[iz+1][iy+1][ix]
    + a001 * zdx[iz][iy][ix+1] + a101 * zdx[iz+1][iy][ix+1] + a011 * zdx[iz][iy+1][ix+1]  + a111 * zdx[iz+1][iy+1][ix+1];
  mat[7] = a000 * zdy[iz][iy][ix] + a100 * zdy[iz+1][iy][ix] + a010 * zdy[iz][iy+1][ix]  + a110 * zdy[iz+1][iy+1][ix]
    + a001 * zdy[iz][iy][ix+1] + a101 * zdy[iz+1][iy][ix+1] + a011 * zdy[iz][iy+1][ix+1]  + a111 * zdy[iz+1][iy+1][ix+1];
  mat[8] = a000 * zdz[iz][iy][ix] + a100 * zdz[iz+1][iy][ix] + a010 * zdz[iz][iy+1][ix]  + a110 * zdz[iz+1][iy+1][ix]
    + a001 * zdz[iz][iy][ix+1] + a101 * zdz[iz+1][iy][ix+1] + a011 * zdz[iz][iy+1][ix+1]  + a111 * zdz[iz+1][iy+1][ix+1];

  vec[0] = a000 * vecx[iz][iy][ix] + a100 * vecx[iz+1][iy][ix] + a010 * vecx[iz][iy+1][ix]  + a110 * vecx[iz+1][iy+1][ix]
    + a001 * vecx[iz][iy][ix+1] + a101 * vecx[iz+1][iy][ix+1] + a011 * vecx[iz][iy+1][ix+1]  + a111 * vecx[iz+1][iy+1][ix+1];
  vec[1] = a000 * vecy[iz][iy][ix] + a100 * vecy[iz+1][iy][ix] + a010 * vecy[iz][iy+1][ix]  + a110 * vecy[iz+1][iy+1][ix]
    + a001 * vecy[iz][iy][ix+1] + a101 * vecy[iz+1][iy][ix+1] + a011 * vecy[iz][iy+1][ix+1]  + a111 * vecy[iz+1][iy+1][ix+1];
  vec[2] = a000 * vecz[iz][iy][ix] + a100 * vecz[iz+1][iy][ix] + a010 * vecz[iz][iy+1][ix]  + a110 * vecz[iz+1][iy+1][ix]
    + a001 * vecz[iz][iy][ix+1] + a101 * vecz[iz+1][iy][ix+1] + a011 * vecz[iz][iy+1][ix+1]  + a111 * vecz[iz+1][iy+1][ix+1];

}
	



static double alpha = 0.5;
static double sigma = 2.0;
static double ERRMAX = 0.05;
static int ITERMAX = 10;


/* On utilise l'algorithme de Newton, comme P. Cachier
   
   principe: soit I le champ inverse et V le champ direct,
   on a M' = M + V(M)
   M = M' + I(M') 
   soit M' = M' + I(M') + V(M' +I(M'))
   donc  I(M) + V(M +I(M)) = 0

   On cherche a minimiser f =  I(M) + V(M +I(M))

   on ajoute une variation delta (d) a I(M) on a   
   f = I(M) + d + V(M + I(M) + d) 
   V(M'+d) = V(M') + (V.N^t)(M') d avec N l'operateur nabla
   donc f = I(M) + d + V(M + I(M)) +  (V.N^t)(M + I(M)) d
   qui s'annule pour 
   d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )

   Question: faut-il faire une adaptation selon les tailles de pixel/voxel ?

*/

static int BAL_Inverse2DVectorField( bal_transformation *theTrsf,  
				     bal_transformation *invTrsf )
{
  char *proc = "BAL_Inverse2DVectorField";
  
  bal_image theXdX, theXdY, theYdX, theYdY;
  bal_doublePoint theSigma;

  int i, j, k;
  float ***arrXdX, ***arrXdY, ***arrYdX, ***arrYdY;  
  float ***arrInvX, ***arrInvY;
  float ***arrTrsfX, ***arrTrsfY;

  int theDim[3] = { theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes };

  double tmpMat[4], tmpVec[2], det;
  double x, y;
  double dx, dy;

  int iterations;
  int niterationsmax = 0;


		    
  theSigma.x = theSigma.y = theSigma.z = sigma;



  if ( invTrsf->vx.nplanes != theTrsf->vx.nplanes ) {
    fprintf( stderr, "%s: input and inverse transformation should have the same z dimension\n", proc );
    return( -1 );
  }


  /* allocation of derivative images
   */
  
  if ( BAL_InitAllocImage( &theXdX, "XderivativeOfX.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of X component\n", proc );
    return( -1 );
  }
  theXdX.vx = theTrsf->vx.vx;   theXdX.vy = theTrsf->vx.vy;   theXdX.vz = theTrsf->vx.vz;

  if ( BAL_InitAllocImage( &theXdY, "YderivativeOfX.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of X component\n", proc );
    return( -1 );
  }
  theXdY.vx = theTrsf->vx.vx;   theXdY.vy = theTrsf->vx.vy;   theXdY.vz = theTrsf->vx.vz;

  if ( BAL_InitAllocImage( &theYdX, "XderivativeOfY.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Y component\n", proc );
    return( -1 );
  }
  theYdX.vx = theTrsf->vx.vx;   theYdX.vy = theTrsf->vx.vy;   theYdX.vz = theTrsf->vx.vz;

  if ( BAL_InitAllocImage( &theYdY, "YderivativeOfY.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Y component\n", proc );
    return( -1 );
  }
  theYdY.vx = theTrsf->vx.vx;   theYdY.vy = theTrsf->vx.vy;   theYdY.vz = theTrsf->vx.vz;




  /* calculation of derivatives
   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... %s: derivatives computation\n", proc );
  }


  if ( BAL_2DDerivativesOfImage( &(theTrsf->vx), 
				 &theXdX, &theXdY,  
				 &theSigma ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of X component\n", proc );
    return( -1 );
  }

  if ( BAL_2DDerivativesOfImage( &(theTrsf->vy), 
				 &theYdX, &theYdY,  
				 &theSigma ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Y component\n", proc );
    return( -1 );
  }



  /* add Identity to V.Nabla^t
     and invert (Id + V.Nabla^t)
     = ( xdx + 1     xdy     )
       ( ydx         ydy + 1 )
       
     The inverse is computed with
     | a11 a12 |-1             |  a22 -a12 |
     | a21 a22 |    =  1/DET * | -a21  a11 |
       
     with DET  =  a11 a22- a12 a21

   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... %s: matrices computation\n", proc );
  }

  arrXdX = (float***)(theXdX.array);
  arrXdY = (float***)(theXdY.array);
  arrYdX = (float***)(theYdX.array);
  arrYdY = (float***)(theYdY.array);

  for ( k=0; k<invTrsf->vx.nplanes; k++ )
  for ( j=0; j<invTrsf->vx.nrows; j++ )
  for ( i=0; i<invTrsf->vx.ncols; i++ ) {
    arrXdX[k][j][i] += 1.0;
    arrYdY[k][j][i] += 1.0;
  }

  for ( k=0; k<invTrsf->vx.nplanes; k++ )
  for ( j=0; j<invTrsf->vx.nrows; j++ )
  for ( i=0; i<invTrsf->vx.ncols; i++ ) {
    det = arrXdX[k][j][i]*arrYdY[k][j][i] - arrXdY[k][j][i] * arrYdX[k][j][i];
    tmpMat[0] = arrYdY[k][j][i] / det;
    tmpMat[1] = (- arrXdY[k][j][i]) / det;
    tmpMat[2] = (- arrYdX[k][j][i]) / det;
    tmpMat[3] = arrXdX[k][j][i] / det;
    arrXdX[k][j][i] = tmpMat[0];
    arrXdY[k][j][i] = tmpMat[1];
    arrYdX[k][j][i] = tmpMat[2];
    arrYdY[k][j][i] = tmpMat[3];
  }


  /*
   */
  arrTrsfX = (float***)(theTrsf->vx.array);
  arrTrsfY = (float***)(theTrsf->vy.array);
  arrInvX = (float***)(invTrsf->vx.array);
  arrInvY = (float***)(invTrsf->vy.array);

  for ( k=0; k<invTrsf->vx.nplanes; k++ ) {

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, " ... processing slice #%3d/%3lu in %s\n", k, invTrsf->vx.nplanes, proc );
    }

    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ ) {

      for ( iterations = 0; iterations <= ITERMAX; iterations++ ) {
	
	/* M + I(M) in millimeter
	 */
	x = i*invTrsf->vx.vx + arrInvX[k][j][i];
	y = j*invTrsf->vx.vy + arrInvY[k][j][i];
	
	/* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
	   and the vector v at M + I(M)
	*/
	
	_get2DMatrixVector( x/theTrsf->vx.vx, y/theTrsf->vx.vy, k,
			    arrXdX, arrXdY, arrYdX, arrYdY,
			    arrTrsfX, arrTrsfY,
			    theDim, tmpMat, tmpVec );
	/* I(M) + v(M+I(M))
       */
	x = arrInvX[k][j][i] + tmpVec[0];
	y = arrInvY[k][j][i] + tmpVec[1];
	/* test on I(M) + v(M+I(M))
	 */
	if ( fabs( x ) + fabs( y ) <= ERRMAX ) 
	  break;
	
	/* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
	   is the optimal variation of I(M)
	*/
	dx = tmpMat[0] * x + tmpMat[1] * y;
	dy = tmpMat[2] * x + tmpMat[3] * y;
	arrInvX[k][j][i] -= alpha * dx;
	arrInvY[k][j][i] -= alpha * dy;
	
      }
      
      if ( iterations > ITERMAX ) niterationsmax ++;
    
    }

  }

  if ( _verbose_ && niterationsmax > 0 ) {
    fprintf( stderr, "%s: non-convergence: %d/%lu\n", proc, niterationsmax, 
	     invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols );
  }
    
  BAL_FreeImage( &theYdY );
  BAL_FreeImage( &theYdX );
  BAL_FreeImage( &theXdY );
  BAL_FreeImage( &theXdX );

  return( 1 );

}











static int BAL_Inverse3DVectorField( bal_transformation *theTrsf,  
				     bal_transformation *invTrsf )
{
  char *proc = "BAL_Inverse3DVectorField";
  
  bal_image theXdX, theXdY, theXdZ, theYdX, theYdY, theYdZ, theZdX, theZdY, theZdZ;
  bal_doublePoint theSigma;
  
  int i, j, k;
  float ***arrXdX, ***arrXdY, ***arrXdZ, ***arrYdX, ***arrYdY, ***arrYdZ, ***arrZdX, ***arrZdY, ***arrZdZ;  
  float ***arrInvX, ***arrInvY, ***arrInvZ;
  float ***arrTrsfX, ***arrTrsfY, ***arrTrsfZ;

  int theDim[3] = { theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes };

  double tmpMat[9], tmpVec[3], det;
  double x, y, z;
  double dx, dy, dz;

  int iterations;
  int niterationsmax = 0;

  

  if ( _debug_ ) fprintf( stderr, "%s: entry\n", proc );



  theSigma.x = theSigma.y = theSigma.z = sigma;



  /* allocation of derivative images
   */

  if ( BAL_InitAllocImage( &theXdX, "XderivativeOfX.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of X component\n", proc );
    return( -1 );
  }
  theXdX.vx = theTrsf->vx.vx;   theXdX.vy = theTrsf->vx.vy;   theXdX.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of X dX\n", proc );


  if ( BAL_InitAllocImage( &theXdY, "YderivativeOfX.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of X component\n", proc );
    return( -1 );
  }
  theXdY.vx = theTrsf->vx.vx;   theXdY.vy = theTrsf->vx.vy;   theXdY.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of X dY\n", proc );


  if ( BAL_InitAllocImage( &theXdZ, "ZderivativeOfX.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of X component\n", proc );
    return( -1 );
  }
  theXdZ.vx = theTrsf->vx.vx;   theXdZ.vy = theTrsf->vx.vy;   theXdZ.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of X dZ\n", proc );


  if ( BAL_InitAllocImage( &theYdX, "XderivativeOfY.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Y component\n", proc );
    return( -1 );
  }
  theYdX.vx = theTrsf->vx.vx;   theYdX.vy = theTrsf->vx.vy;   theYdX.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Y dX\n", proc );


  if ( BAL_InitAllocImage( &theYdY, "YderivativeOfY.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Y component\n", proc );
    return( -1 );
  }
  theYdY.vx = theTrsf->vx.vx;   theYdY.vy = theTrsf->vx.vy;   theYdY.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Y dY\n", proc );


  if ( BAL_InitAllocImage( &theYdZ, "ZderivativeOfY.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of Y component\n", proc );
    return( -1 );
  }
  theYdZ.vx = theTrsf->vx.vx;   theYdZ.vy = theTrsf->vx.vy;   theYdZ.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Y dZ\n", proc );


  if ( BAL_InitAllocImage( &theZdX, "XderivativeOfZ.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate X derivative of Z component\n", proc );
    return( -1 );
  }
  theZdX.vx = theTrsf->vx.vx;   theZdX.vy = theTrsf->vx.vy;   theZdX.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Z dX\n", proc );


  if ( BAL_InitAllocImage( &theZdY, "YderivativeOfZ.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theZdX );   
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Y derivative of Z component\n", proc );
    return( -1 );
  }
  theZdY.vx = theTrsf->vx.vx;   theZdY.vy = theTrsf->vx.vy;   theZdY.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Z dY\n", proc );


  if ( BAL_InitAllocImage( &theZdZ, "ZderivativeOfZ.inr", 
			   theTrsf->vx.ncols, theTrsf->vx.nrows, theTrsf->vx.nplanes, 1, FLOAT ) != 1 ) {
    BAL_FreeImage( &theZdY );   
    BAL_FreeImage( &theZdX );   
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate Z derivative of Z component\n", proc );
    return( -1 );
  }
  theZdZ.vx = theTrsf->vx.vx;   theZdZ.vy = theTrsf->vx.vy;   theZdZ.vz = theTrsf->vx.vz;
  if ( _debug_ ) fprintf( stderr, "%s: allocation of Z dZ\n", proc );






  /* calculation of derivatives
   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... %s: derivatives computation\n", proc );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vx), 
				 &theXdX, &theXdY, &theXdZ,  
				 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );   
    BAL_FreeImage( &theZdY );   
    BAL_FreeImage( &theZdX );   
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of X component\n", proc );
    return( -1 );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vy), 
				 &theYdX, &theYdY, &theYdZ,  
				 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );   
    BAL_FreeImage( &theZdY );   
    BAL_FreeImage( &theZdX );   
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Y component\n", proc );
    return( -1 );
  }

  if ( BAL_3DDerivativesOfImage( &(theTrsf->vz), 
				 &theZdX, &theZdY, &theZdZ,  
				 &theSigma ) != 1 ) {
    BAL_FreeImage( &theZdZ );   
    BAL_FreeImage( &theZdY );   
    BAL_FreeImage( &theZdX );   
    BAL_FreeImage( &theYdZ );
    BAL_FreeImage( &theYdY );
    BAL_FreeImage( &theYdX );
    BAL_FreeImage( &theXdZ );
    BAL_FreeImage( &theXdY );
    BAL_FreeImage( &theXdX );    
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute derivatives of Z component\n", proc );
    return( -1 );
  }
  


  /* add Identity to V.Nabla^t
     and invert (Id + V.Nabla^t)
     = ( xdx + 1     xdy         xdz     )
       ( ydx         ydy + 1     ydz     )
       ( zdx         zdy         zdz + 1 )

     The inverse is computed with
     | a11 a12 a13 |-1             |   a33a22-a32a23  -(a33a12-a32a13)   a23a12-a22a13  |
     | a21 a22 a23 |    =  1/DET * | -(a33a21-a31a23)   a33a11-a31a13  -(a23a11-a21a13) |
     | a31 a32 a33 |               |   a32a21-a31a22  -(a32a11-a31a12)   a22a11-a21a12  |

     with DET  =  a11(a33a22-a32a23)-a21(a33a12-a32a13)+a31(a23a12-a22a13)

   */

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, " ... %s: matrices computation\n", proc );
  }

  arrXdX = (float***)(theXdX.array);
  arrXdY = (float***)(theXdY.array);
  arrXdZ = (float***)(theXdZ.array);
  arrYdX = (float***)(theYdX.array);
  arrYdY = (float***)(theYdY.array);
  arrYdZ = (float***)(theYdZ.array);
  arrZdX = (float***)(theZdX.array);
  arrZdY = (float***)(theZdY.array);
  arrZdZ = (float***)(theZdZ.array);

  for ( k=0; k<invTrsf->vx.nplanes; k++ )
  for ( j=0; j<invTrsf->vx.nrows; j++ )
  for ( i=0; i<invTrsf->vx.ncols; i++ ) {
    arrXdX[k][j][i] += 1.0;
    arrYdY[k][j][i] += 1.0;
    arrZdZ[k][j][i] += 1.0;
  }

  for ( k=0; k<invTrsf->vx.nplanes; k++ )
  for ( j=0; j<invTrsf->vx.nrows; j++ )
  for ( i=0; i<invTrsf->vx.ncols; i++ ) {
    det = arrXdX[k][j][i] * ( arrYdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrYdZ[k][j][i] )
      - arrYdX[k][j][i] * ( arrXdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrXdZ[k][j][i] )
      + arrZdX[k][j][i] * ( arrXdY[k][j][i]*arrYdZ[k][j][i] - arrYdY[k][j][i]*arrXdZ[k][j][i] );
    tmpMat[0] = ( arrYdY[k][j][i]*arrZdZ[k][j][i] - arrZdY[k][j][i]*arrYdZ[k][j][i] ) / det;
    tmpMat[1] = ( - arrXdY[k][j][i]*arrZdZ[k][j][i] + arrZdY[k][j][i]*arrXdZ[k][j][i] ) / det;
    tmpMat[2] = ( arrXdY[k][j][i]*arrYdZ[k][j][i] - arrYdY[k][j][i]*arrXdZ[k][j][i] ) / det;
    tmpMat[3] = ( - arrYdX[k][j][i]*arrZdZ[k][j][i] + arrZdX[k][j][i]*arrYdZ[k][j][i] ) / det;
    tmpMat[4] = ( arrXdX[k][j][i]*arrZdZ[k][j][i] - arrZdX[k][j][i]*arrXdZ[k][j][i] ) / det;
    tmpMat[5] = ( - arrXdX[k][j][i]*arrYdZ[k][j][i] + arrYdX[k][j][i]*arrXdZ[k][j][i] ) / det;
    tmpMat[6] = ( arrYdX[k][j][i]*arrZdY[k][j][i] - arrZdX[k][j][i]*arrYdY[k][j][i] ) / det;
    tmpMat[7] = ( - arrXdX[k][j][i]*arrZdY[k][j][i] + arrZdX[k][j][i]*arrXdY[k][j][i] ) / det;
    tmpMat[8] = ( arrXdX[k][j][i]*arrYdY[k][j][i] - arrYdX[k][j][i]*arrXdY[k][j][i] ) / det;

    arrXdX[k][j][i] = tmpMat[0];
    arrXdY[k][j][i] = tmpMat[1];
    arrXdZ[k][j][i] = tmpMat[2];
    arrYdX[k][j][i] = tmpMat[3];
    arrYdY[k][j][i] = tmpMat[4];
    arrYdZ[k][j][i] = tmpMat[5];
    arrZdX[k][j][i] = tmpMat[6];
    arrZdY[k][j][i] = tmpMat[7];
    arrZdZ[k][j][i] = tmpMat[8];
  }


  /*
   */
  arrTrsfX = (float***)(theTrsf->vx.array);
  arrTrsfY = (float***)(theTrsf->vy.array);
  arrTrsfZ = (float***)(theTrsf->vz.array);
  arrInvX = (float***)(invTrsf->vx.array);
  arrInvY = (float***)(invTrsf->vy.array);
  arrInvZ = (float***)(invTrsf->vz.array);

  for ( k=0; k<invTrsf->vx.nplanes; k++ ) {

    if ( _verbose_ >= 2 ) {
      fprintf( stderr, " ... processing slice #%3d/%3lu in %s\n", k, invTrsf->vx.nplanes, proc );
    }

    for ( j=0; j<invTrsf->vx.nrows; j++ )
    for ( i=0; i<invTrsf->vx.ncols; i++ ) {

      for ( iterations = 0; iterations<= ITERMAX; iterations++ ) {
	
	/* M + I(M) in millimeter
	 */
	x = i*invTrsf->vx.vx + arrInvX[k][j][i];
	y = j*invTrsf->vx.vy + arrInvY[k][j][i];
	z = k*invTrsf->vx.vz + arrInvZ[k][j][i];
	
	/* get the matrix  ( Id + (V.N^t) )^(-1) at M + I(M)
	   and the vector v at M + I(M)
	*/
	
	_get3DMatrixVector( x/theTrsf->vx.vx, y/theTrsf->vx.vy, z/theTrsf->vx.vz,
			    arrXdX, arrXdY, arrXdZ, 
			    arrYdX, arrYdY, arrYdZ,
			    arrZdX, arrZdY, arrZdZ,
			    arrTrsfX, arrTrsfY, arrTrsfZ,
			    theDim, tmpMat, tmpVec );
	/* I(M) + v(M+I(M))
	 */
	x = arrInvX[k][j][i] + tmpVec[0];
	y = arrInvY[k][j][i] + tmpVec[1];
	z = arrInvZ[k][j][i] + tmpVec[2];
	/* test on I(M) + v(M+I(M))
	 */
	if ( fabs( x ) + fabs( y ) + fabs( z ) <= ERRMAX ) 
	  break;
	
	/* d = - ( Id + (V.N^t)(M+I(M)) )^(-1) ( I(M)+v(M+I(M)) )
	   is the optimal variation of I(M)
	*/
	dx = tmpMat[0] * x + tmpMat[1] * y + tmpMat[2] * z;
	dy = tmpMat[3] * x + tmpMat[4] * y + tmpMat[5] * z;
	dz = tmpMat[6] * x + tmpMat[7] * y + tmpMat[8] * z;
	arrInvX[k][j][i] -= alpha * dx;
	arrInvY[k][j][i] -= alpha * dy;
	arrInvZ[k][j][i] -= alpha * dz;
	
      }
      
      if ( iterations > ITERMAX ) niterationsmax ++;
    
    }
    
  }

  if ( _verbose_ && niterationsmax > 0 ) {
    fprintf( stderr, "%s: non-convergence: %d/%lu\n", proc, niterationsmax, 
	     invTrsf->vx.nplanes * invTrsf->vx.nrows * invTrsf->vx.ncols );
  }
  
  BAL_FreeImage( &theZdZ );   
  BAL_FreeImage( &theZdY );   
  BAL_FreeImage( &theZdX );   
  BAL_FreeImage( &theYdZ );
  BAL_FreeImage( &theYdY );
  BAL_FreeImage( &theYdX );
  BAL_FreeImage( &theXdZ );
  BAL_FreeImage( &theXdY );
  BAL_FreeImage( &theXdX );

  return( 1 );

}










/*--------------------------------------------------
 *
 * Non-Linear Transformations Construction
 *
 --------------------------------------------------*/



static double sinusoid_amplitude[3] = {2.0, 2.0, 2.0};
static double sinusoid_period[3] = {50.0, 50.0, 50.0};

static int BAL_SinusoidVectorField( bal_transformation *theTrsf, int drawx, int drawy, int drawz )
{
  int i, x, y, z;
  float *xbuf, *ybuf, *zbuf;
  double wx, wy, wz;

  if ( drawz && drawy && drawx ) {
    xbuf = (float*)(theTrsf->vx.data);
    ybuf = (float*)(theTrsf->vy.data);
    zbuf = (float*)(theTrsf->vz.data);
    wx = 2 * 3.14159265 / sinusoid_period[0];
    wy = 2 * 3.14159265 / sinusoid_period[1];
    wz = 2 * 3.14159265 / sinusoid_period[2];
    for ( i=0, z=0; z<theTrsf->vx.nplanes; z++ )
    for ( y=0; y<theTrsf->vx.nrows; y++ )
    for ( x=0; x<theTrsf->vx.ncols; x++, i++ )
      xbuf[i] = ybuf[i] = zbuf[i] = sinusoid_amplitude[2] * sin( wz * (z-theTrsf->vz.nplanes/2) )
	* sinusoid_amplitude[1] * sin( wy * (y-theTrsf->vy.nrows/2) )
	* sinusoid_amplitude[0] * sin( wx * (x-theTrsf->vx.ncols/2) );
  }

  if ( !drawz && drawy && drawx ) {
    xbuf = (float*)(theTrsf->vx.data);
    ybuf = (float*)(theTrsf->vy.data);
    wx = 2 * 3.14159265 / sinusoid_period[0];
    wy = 2 * 3.14159265 / sinusoid_period[1];
    for ( i=0, z=0; z<theTrsf->vx.nplanes; z++ )
    for ( y=0; y<theTrsf->vx.nrows; y++ )
    for ( x=0; x<theTrsf->vx.ncols; x++, i++ )
      xbuf[i] = ybuf[i] = sinusoid_amplitude[1] * sin( wx * (y-theTrsf->vy.nrows/2) )
	* sinusoid_amplitude[0] * sin( wy * (x-theTrsf->vx.ncols/2) );
  }

  return( 1 );

}


int BAL_Sinusoid3DVectorField( bal_transformation *theTrsf )
{
  sinusoid_amplitude[0] = 1.3;
  sinusoid_amplitude[1] = 1.3;
  sinusoid_amplitude[2] = 1.3;
  return( BAL_SinusoidVectorField( theTrsf, 1, 1, 1 ) );
}

int BAL_Sinusoid2DVectorField( bal_transformation *theTrsf )
{
  sinusoid_amplitude[0] = 1.75;
  sinusoid_amplitude[1] = 1.75;
  sinusoid_amplitude[2] = 1.75;
  return( BAL_SinusoidVectorField( theTrsf, 1, 1, 0 ) );
}



/*--------------------------------------------------
 *
 * Random Transformations
 *
 --------------------------------------------------*/

static double angle_interval[2] = {0.0, 1.57079}; /* 0, pi/2 */
static double scale_interval[2] = {0.7, 1.4};
static double shear_interval[2] = {0.0, 0.3};
static double translation_interval[2] = {-10, 10};



int BAL_Random2DTranslationMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _IdentityMatrix( m );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DTranslationMatrix( bal_transformation *theTrsf )
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _IdentityMatrix( m );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DTranslationScalingMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DScaleMatrix( m, scale_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DTranslationScalingMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DScaleMatrix( m, scale_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DRigidMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DRotationMatrix( m, angle_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DRigidMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DRotationMatrix( m, angle_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DSimilitudeMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DSimilitudeMatrix( m, angle_interval, scale_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DSimilitudeMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DSimilitudeMatrix( m, angle_interval, scale_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DAffineMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random2DAffineMatrix( m, angle_interval, scale_interval, shear_interval );
  _Random2DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random3DAffineMatrix( bal_transformation *theTrsf ) 
{
  double m[9];
  double v[3] = {0.0, 0.0, 0.0};

  BAL_SetTransformationToIdentity( theTrsf );
  
  _Random3DAffineMatrix( m, angle_interval, scale_interval, shear_interval );
  _Random3DTranslationVector( v, translation_interval );
  
  theTrsf->mat.m[ 0] = m[0];   theTrsf->mat.m[ 1] = m[1];   theTrsf->mat.m[ 2] = m[2];
  theTrsf->mat.m[ 4] = m[3];   theTrsf->mat.m[ 5] = m[4];   theTrsf->mat.m[ 6] = m[5];
  theTrsf->mat.m[ 8] = m[6];   theTrsf->mat.m[ 9] = m[7];   theTrsf->mat.m[10] = m[8];

  theTrsf->mat.m[ 3] = v[0];   theTrsf->mat.m[ 7] = v[1];   theTrsf->mat.m[11] = v[2]; 

  return( 1 );
}



int BAL_Random2DVectorField( bal_transformation *theTrsf ) 
{
  return( -1 );
}



int BAL_Random3DVectorField( bal_transformation *theTrsf ) 
{
  return( -1 );
}



int BAL_RandomTransformation( bal_transformation *theTrsf, enumTypeTransfo type, bal_image *ref )
{
  char *proc = "BAL_RandomTransformation";
  
  if ( BAL_AllocTransformation( theTrsf, type, ref ) != 1 ) {
    if ( _verbose_ ) 
      fprintf( stderr, "%s: unable to allocate transformation\n", proc );
    return( -1 );
  }

  switch( theTrsf->type ) {
    
  default :
    if ( _verbose_ )
      fprintf( stderr, "%s: such type not handled yet\n", proc );
    return( -1 );

  case VECTORFIELD_2D :
    if (  BAL_Random2DVectorField( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D vector field\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;
    
  case VECTORFIELD_3D :
    if (  BAL_Random3DVectorField( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D vector field\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;
    
  case TRANSLATION_2D :
    if ( BAL_Random2DTranslationMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D translation matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_3D :
    if ( BAL_Random3DTranslationMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D translation matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_SCALING_2D :
    if ( BAL_Random2DTranslationScalingMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D translation and scaling matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case TRANSLATION_SCALING_3D :
    if ( BAL_Random3DTranslationScalingMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D translation and scaling matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case RIGID_2D :
    if ( BAL_Random2DRigidMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D rigid matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case RIGID_3D :
    if ( BAL_Random3DRigidMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D rigid matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case SIMILITUDE_2D :
    if ( BAL_Random2DSimilitudeMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D similitude matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case SIMILITUDE_3D :
    if ( BAL_Random3DSimilitudeMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D similitude matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case AFFINE_2D :
    if ( BAL_Random2DAffineMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 2D affine matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  case AFFINE_3D :
    if ( BAL_Random3DAffineMatrix( theTrsf ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to generate 3D affine matrix\n", proc );
      BAL_FreeTransformation( theTrsf );
      return( -1 );
    }
    break;

  }

  return( 1 );
}
