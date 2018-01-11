/*************************************************************************
 * bal-transformation-list-tools.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mar 21 jan 2014 18:00:31 CET
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
#include <bal-transformation-list-tools.h>

static int _verbose_ = 1;
static int _debug_ = 0;






static void _computeDim( int *ldim, int *rdim, double leftb, double rightb, double vsize )
{
  char *proc = "_computeDim";
  int l=0;
  int r=0;

  if ( rightb <= 0.0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird value of right boundary (%f)\n", proc, rightb );
  }

  if ( leftb > 0.0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: weird value of left boundary (%f)\n", proc, leftb );
  }

  r = (int)( rightb / vsize );
  if ( rightb - r * vsize > 0.0 
       && rightb - r * vsize > 0.0001 * vsize ) r++;

  if ( _debug_ >= 2  )
    fprintf( stderr, "%f / %f = %f -> %d\n", rightb, vsize, rightb / vsize, r );
  
  if ( leftb < 0.0 ) {
    l = (int)( (-leftb) / vsize );
    if ( (-leftb) - l * vsize > 0.0 
	 && (-leftb) - l * vsize >  0.0001 * vsize ) l++;
    if ( _debug_ >= 2  )
      fprintf( stderr, "%f / %f = %f -> %d\n", (-leftb), vsize, (-leftb) / vsize, l );
  }
  
  if ( _debug_ >= 2 ) {
    fprintf( stderr, "%s: interval of [%f %f]\n", proc, leftb, rightb );
    fprintf( stderr, "\t -> dimensions of [%d + %d] with voxel size of %f\n", l, r, vsize );
  }
  
  *ldim = l;
  *rdim = r;
}





/* les transformations de depart sont du type T_{i<-ref}
   on se donne une reference T_{r<-ref}

   1. avec les T_{r<-ref} o T_{i<-ref}^{-1} = T_{r<-i}
      on connait la transformation de la boite englobante 
      des images I_i et on deduit une boite englobante globale,
   2. on calcule la translation T_{r<-newref} qui permet de mettre
      toutes les images dans un seul template, avec I_r a une 
      frontiere de pixel (translation entiere)
   3. on calcule les transformations resultats
      T_{i<-ref} o T_{r<-ref}^{-1} o T_{r<-newref} 
*/


int BAL_ChangeTransformationList( bal_transformationList *theList,
				  bal_image *theIm,
				  bal_transformationList *resList,
				  bal_image *resIm,
				  int refTrsf,
				  int margin,
				  int isotropic )
{
  char *proc = "BAL_ChangeTransformationList";
  bal_transformation tmpTrsf, invTrsf;
  int i, n;

  bal_doublePoint corner[8];
  bal_doublePoint tmpPt;
  bal_floatPoint mincorner;
  bal_floatPoint maxcorner;

  bal_integerPoint minindex;
  bal_integerPoint maxindex;

  float vx = theIm->vx;   
  float vy = theIm->vy;   
  float vz = theIm->vz;   

  bal_integerPoint negativedim;
  bal_integerPoint positivedim;

  if ( refTrsf < 0 || refTrsf >= theList->n_trsfs ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: index of reference transformation out of range\n", proc );
    return( -1 );
  }

  if ( isotropic ) {
    if ( vx < vy ) {
      if ( vx < vz ) {
	vy = vz = vx;
      }
      else {
	vx = vy = vz;
      }
    }
    else {
      if ( vy < vz ) {
	vx = vz = vy;
      }
      else {
	vx = vy = vz;
      }
    }
    if ( _debug_ >= 2 ) {
      fprintf( stderr, "%s: new voxel sizes = [%f %f %f]\n", proc, vx, vy, vz );
    }
  }



  /* boite englobante
     les 8 coins de l'image
   */

  mincorner.x = mincorner.y = mincorner.z = 0.0;
  maxcorner.x = (theIm->ncols - 1) * theIm->vx;
  maxcorner.y = (theIm->nrows - 1) * theIm->vy;
  maxcorner.z = (theIm->nplanes - 1) * theIm->vz;

  minindex.x = minindex.y = minindex.z = refTrsf;
  maxindex.x = maxindex.y = maxindex.z = refTrsf;

  corner[0].z = mincorner.z;
  corner[1].z = mincorner.z;
  corner[2].z = mincorner.z;
  corner[3].z = mincorner.z;

  corner[4].z = maxcorner.z;
  corner[5].z = maxcorner.z;
  corner[6].z = maxcorner.z;
  corner[7].z = maxcorner.z;

  corner[0].y = mincorner.y;
  corner[1].y = mincorner.y;
  corner[2].y = maxcorner.y;
  corner[3].y = maxcorner.y;

  corner[4].y = mincorner.y;
  corner[5].y = mincorner.y;
  corner[6].y = maxcorner.y;
  corner[7].y = maxcorner.y;

  corner[0].x = mincorner.x;
  corner[1].x = maxcorner.x;
  corner[2].x = mincorner.x;
  corner[3].x = maxcorner.x;

  corner[4].x = mincorner.x;
  corner[5].x = maxcorner.x;
  corner[6].x = mincorner.x;
  corner[7].x = maxcorner.x;

  if ( _debug_ >= 2 ) {
    fprintf( stderr, "%s: Image corners\n", proc );
    if ( 0 ) {
      for ( i=0; i<8; i++ )
	fprintf( stderr, "   #%d = [%9.3f %9.3f %9.3f]\n", 
		 i, corner[i].x, corner[i].y, corner[i].z ); 
    }
    fprintf( stderr, "   left corner  = [%9.3f %9.3f %9.3f]\n", 
	     mincorner.x, mincorner.y, mincorner.z ); 
    fprintf( stderr, "   right corner = [%9.3f %9.3f %9.3f]\n", 
	     maxcorner.x, maxcorner.y, maxcorner.z ); 
  }
  
  

  BAL_InitTransformation( &invTrsf );
  if ( BAL_AllocTransformation( &invTrsf, theList->data[0].type, (bal_image *)NULL ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate inverse transformation\n", proc );
    return( -1 );
  }

  for ( n=0; n<theList->n_trsfs; n++ ) {

    if ( n == refTrsf ) continue;

    if ( BAL_InverseTransformation(  &(theList->data[n]), &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to invert transformation #%d\n", proc, n );
      return( -1 );
    }

    if ( BAL_TransformationComposition( &invTrsf, 
					&(theList->data[refTrsf]), &invTrsf ) != 1 ) {
      BAL_FreeTransformation( &invTrsf );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compose with inverse transformation #%d\n", proc, n );
      return( -1 );
    }
    
    for ( i=0; i<8; i++ ) {
      if ( BAL_TransformPoint( &(corner[i]), &tmpPt, &invTrsf ) != 1 ) {
	BAL_FreeTransformation( &invTrsf );
	if ( _verbose_ )
	  fprintf( stderr, "%s: unable to transform corner #%d with transformation #%d\n", proc, i, n );
	return( -1 );
      }
      if ( mincorner.x > tmpPt.x ) {
	mincorner.x = tmpPt.x;
	minindex.x  = n;
      }
      if ( maxcorner.x < tmpPt.x ) {
	maxcorner.x = tmpPt.x;
	maxindex.x  = n;
      }
      if ( mincorner.y > tmpPt.y ) {
	mincorner.y = tmpPt.y;
	minindex.y  = n;
      }
      if ( maxcorner.y < tmpPt.y ) {
	maxcorner.y = tmpPt.y;
	maxindex.y  = n;
      }
      if ( mincorner.z > tmpPt.z ) {
	mincorner.z = tmpPt.z;
	minindex.z  = n;
      }
      if ( maxcorner.z < tmpPt.z ) {
	maxcorner.z = tmpPt.z;
	maxindex.z  = n;
      }
    }

  }
  
  if ( _debug_ >= 2 ) {
    fprintf( stderr, "%s: Computed corners\n", proc );
    fprintf( stderr, "   left corner  = [%9.3f %9.3f %9.3f]\n", 
	     mincorner.x, mincorner.y, mincorner.z ); 
    fprintf( stderr, "   right corner = [%9.3f %9.3f %9.3f]\n", 
	     maxcorner.x, maxcorner.y, maxcorner.z ); 
    fprintf( stderr, "%s: Computed extremums\n", proc );
    fprintf( stderr, "   left indices  = [%3d %3d %3d]\n", 
	     minindex.x, minindex.y, minindex.z ); 
    fprintf( stderr, "   right indices = [%3d %3d %3d]\n", 
	     maxindex.x, maxindex.y, maxindex.z ); 
  }



  /* dimensions of the large image
   */

  _computeDim( &(negativedim.x), &(positivedim.x), mincorner.x, maxcorner.x, vx );
  _computeDim( &(negativedim.y), &(positivedim.y), mincorner.y, maxcorner.y, vy );
  _computeDim( &(negativedim.z), &(positivedim.z), mincorner.z, maxcorner.z, vz );

  if ( margin > 0 ) {
    negativedim.x += margin;
    negativedim.y += margin;
    negativedim.z += margin;
    positivedim.x += margin;
    positivedim.y += margin;
    positivedim.z += margin;
  }

  if ( BAL_InitImage( resIm, (char *)NULL,
		      negativedim.x + 1 + positivedim.x,
		      negativedim.y + 1 + positivedim.y,
		      negativedim.z + 1 + positivedim.z, 1, theIm->type ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to initialize result template\n", proc );
    return( -1 );
  }
  resIm->vx = vx;
  resIm->vy = vy;
  resIm->vz = vz;




  /* construction of T_{r<-newref}
   */

  BAL_InitTransformation( &tmpTrsf );
  if ( BAL_AllocTransformation( &tmpTrsf, TRANSLATION_3D, (bal_image *)NULL ) != 1 ) {
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate auxiliary transformation\n", proc );
    return( -1 );
  }
  BAL_SetTransformationToIdentity( &tmpTrsf );
  tmpTrsf.mat.m[ 3] = - negativedim.x * vx;
  tmpTrsf.mat.m[ 7] = - negativedim.y * vy;
  tmpTrsf.mat.m[11] = - negativedim.z * vz;

  if ( BAL_InverseTransformation(  &(theList->data[refTrsf]), &invTrsf ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to invert reference transformation #%d\n", proc, refTrsf );
    return( -1 );
  }

  /* compose transformations
     int BAL_TransformationComposition( bal_transformation *t_res,
     bal_transformation *t1, 
     bal_transformation *t2 ) 
     with t_res = t1 o t2 
  */
  if ( BAL_TransformationComposition( &tmpTrsf, &invTrsf, &tmpTrsf ) != 1 ) {
    BAL_FreeTransformation( &tmpTrsf );
    BAL_FreeTransformation( &invTrsf );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compose transformations (preparation)\n", proc );
    return( -1 );
  }

  BAL_FreeTransformation( &invTrsf );


  for ( n=0; n<theList->n_trsfs; n++ ) {
    if ( BAL_TransformationComposition( &(resList->data[n]),
					&(theList->data[n]), &tmpTrsf ) != 1 ) {
      BAL_FreeTransformation( &tmpTrsf );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to compose transformations (#%d)\n", proc, n );
      return( -1 );
    }
  }

  BAL_FreeTransformation( &tmpTrsf );
  return( 0 );
}


