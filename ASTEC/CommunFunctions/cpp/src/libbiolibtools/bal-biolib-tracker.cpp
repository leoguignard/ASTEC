/*************************************************************************
 * bal-biolib-tracker.cpp -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 21 mar 2014 17:41:11 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <iostream>
#include <string>


#include <math.h>

#include <lemon/list_graph.h>
#include <lemon/bellman_ford.h>

#include <bal-biolib-tools.h>
#include <bal-biolib-tracker.h>
#include <bal-transformation-tools.h>




typedef struct {
  lemon::ListDigraph graph;
  lemon::ListDigraph::Node source;
  lemon::ListDigraph::Node target;
} typeGraph;



static int _debug_ = 0;
static int _verbose_ = 2;





void BAL_SetVerboseInBalBiolibTracker( int v )
{
  _verbose_ = v;
}

void BAL_IncrementVerboseInBalBiolibTracker(  )
{
  _verbose_ ++;
}

void BAL_DecrementVerboseInBalBiolibTracker(  )
{
  _verbose_ --;
  if ( _verbose_ < 0 ) _verbose_ = 0;
}










/************************************************************
 *
 * bal_blTrackerElement related stuff
 *
 ************************************************************/


void BAL_InitBlTrackerElement( bal_blTrackerElement *e )
{

  e->index = -1;
  
  e->detectionName = (char*)NULL;

  BAL_InitBlDetectionList( &(e->readDetectionList) );
  BAL_InitBlDetectionList( &(e->trsfDetectionList) );
  e->detectionList = (bal_blDetectionList *)NULL;

  e->transformationName =  (char*)NULL;
  BAL_InitTransformation( &(e->theTrsf) );
  BAL_InitTransformation( &(e->invTrsf) );

  e->vesselnessName = (char*)NULL;
  BAL_InitImage ( &(e->vesselnessImage), NULL, 0, 0, 0, 0, UCHAR );
}



void BAL_FreeBlTrackerElement( bal_blTrackerElement *e )
{

  BAL_FreeBlDetectionList( &(e->readDetectionList) );
  BAL_FreeBlDetectionList( &(e->trsfDetectionList) );

  if ( e->transformationName != (char*)NULL ) {
    BAL_FreeTransformation( &(e->theTrsf) );
    BAL_FreeTransformation( &(e->invTrsf) );
  }

  if ( e->vesselnessName != (char*) NULL ) {
    BAL_FreeImage ( &(e->vesselnessImage) );
  }

  BAL_InitBlTrackerElement( e );
}







/* lit tous les elements
 */

int BAL_ReadDetectionsBlTrackerElement( bal_blTrackerElement *e )
{
  std::string proc = "BAL_ReadDetectionsBlTrackerElement";

  
  if ( BAL_ReadBlDetectionList( &(e->readDetectionList), e->detectionName ) != 0 ) {
    if ( _verbose_ ) 
      std::cerr << proc << ": unable to read detection file '" << e->detectionName << "'" << std::endl;
    return( -1 );
  }
  
  if ( e->transformationName != (char*)NULL ) {
    if ( BAL_ReadTransformation( &(e->theTrsf), e->transformationName  ) != 1 ) {
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to read transformation file '" << e->transformationName << "'" << std::endl;
      return( -1 );
    }
    if ( BAL_AllocTransformation( &(e->invTrsf), e->theTrsf.type, (bal_image*)NULL ) != 1 ) {
      BAL_FreeTransformation( &(e->theTrsf) );
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ )
	std::cerr << proc << ": unable to  allocate result transformation" << std::endl;
      return( -1 );
    }
    if ( BAL_InverseTransformation( &(e->theTrsf), &(e->invTrsf) ) != 1 ) {
      BAL_FreeTransformation( &(e->theTrsf) );
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to invert transformation from file '" << e->transformationName << "'" << std::endl;
      return( -1 );
    }
  }
  
  if ( e->vesselnessName != (char*)NULL ) {
    if ( BAL_ReadImage( &(e->vesselnessImage), e->vesselnessName, 0 ) != 1 ) {
      if ( e->transformationName != (char*)NULL ) {
	BAL_FreeTransformation( &(e->invTrsf) );
	BAL_FreeTransformation( &(e->theTrsf) );
      }
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ )
	std::cerr << proc << ": unable to read image '" << e->vesselnessName << "'" << std::endl;
      return( -1 );
    }
  }


  return( 1 );
}





/* lit tous les elements
   sauf la liste des detections
 */
int BAL_PartialReadDetectionsBlTrackerElement( bal_blTrackerElement *e )
{
  std::string proc = "BAL_PartialReadDetectionsBlTrackerElement";

  if ( e->transformationName != (char*)NULL ) {
    if ( BAL_ReadTransformation( &(e->theTrsf), e->transformationName  ) != 1 ) {
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to read transformation file '" << e->transformationName << "'" << std::endl;
      return( -1 );
    }
    if ( BAL_AllocTransformation( &(e->invTrsf), e->theTrsf.type, (bal_image*)NULL ) != 1 ) {
      BAL_FreeTransformation( &(e->theTrsf) );
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ )
	std::cerr << proc << ": unable to  allocate result transformation" << std::endl;
      return( -1 );
    }
    if ( BAL_InverseTransformation( &(e->theTrsf), &(e->invTrsf) ) != 1 ) {
      BAL_FreeTransformation( &(e->theTrsf) );
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to invert transformation from file '" << e->transformationName << "'" << std::endl;
      return( -1 );
    }
  }
  
  if ( e->vesselnessName != (char*)NULL ) {
    if ( BAL_ReadImage( &(e->vesselnessImage), e->vesselnessName, 0 ) != 1 ) {
      if ( e->transformationName != (char*)NULL ) {
	BAL_FreeTransformation( &(e->invTrsf) );
	BAL_FreeTransformation( &(e->theTrsf) );
      }
      BAL_FreeBlDetectionList( &(e->readDetectionList) );
      if ( _verbose_ )
	std::cerr << proc << ": unable to read image '" << e->vesselnessName << "'" << std::endl;
      return( -1 );
    }
  }


  return( 1 );
}



void BAL_PartialFreeBlTrackerElement( bal_blTrackerElement *e )
{

  BAL_FreeBlDetectionList( &(e->trsfDetectionList) );

  if ( e->vesselnessName != (char*) NULL ) {
    BAL_FreeImage ( &(e->vesselnessImage) );
  }

}





int BAL_UpdateDetectionsBlTrackerElement( bal_blTrackerElement *e )
{
  std::string proc = "BAL_UpdateDetectionsBlTrackerElement";
  bal_transformation theTrsf;
  int n;

  BAL_InitTransformation( &theTrsf );

  /* pas de transformation
     si pas d'image on fait rien
     si image mais deja en reel, on fait rien
   */
  if ( e->transformationName == (char*)NULL ) {
    
    if ( e->vesselnessName == (char*)NULL 
	 || e->readDetectionList.unit == REAL_UNIT ) {
      e->detectionList = &(e->readDetectionList);
      return( 1 );
    }
  }

  
  
  /* matrice homothetie pour le passage des voxels aux reels
   */

  if ( BAL_AllocTransformation( &theTrsf, AFFINE_3D, (bal_image*)NULL ) != 1 ) {
    if ( _verbose_ )
      std::cerr << proc << ": unable to allocate transformation" << std::endl;
    return( -1 );
  }
  
  BAL_SetTransformationToIdentity( &theTrsf );

  if ( e->readDetectionList.unit == VOXEL_UNIT ) {
    if ( e->vesselnessName == (char*)NULL ) {
      if ( _verbose_ )
	std::cerr << proc << ":weird, no image is given" << std::endl;
    }
    else {
      theTrsf.mat.m[0]  = e->vesselnessImage.vx;
      theTrsf.mat.m[5]  = e->vesselnessImage.vy;
      theTrsf.mat.m[10] = e->vesselnessImage.vz;
      theTrsf.mat.m[15] = 1;
    }
  }
  


  /* matrice transformation finale
     pour le passage des voxels a des reels 
     de la frame i dans une reference commune
     T_{ref <-i} = T_{ref <-i} o H_{i}
     e->invTrsf = e->invTrsf o theTrsf
   */

  if ( e->transformationName != (char*)NULL ) {
    if ( BAL_TransformationComposition( &(e->invTrsf), &(e->invTrsf), &theTrsf ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
	std::cerr << proc << ": unable to compose transformation (1)" << std::endl;
      return( -1 );
    }
  } else {
    /* should allocate and copy
     */
    if ( _verbose_ )
      std::cerr << proc << ":weird, no transformation is given" << std::endl;
    return( -1 );
  }
  

  theTrsf.mat.m[0]  = 1.0 / theTrsf.mat.m[0];
  theTrsf.mat.m[5]  = 1.0 / theTrsf.mat.m[5];
  theTrsf.mat.m[10] = 1.0 / theTrsf.mat.m[10];


  /* matrice transformation finale
     pour le passage des reels a des voxels 
     d'une reference commune a la fram i
     T_{i <- ref} = H_{i}^(-1) o  T_{i <- ref}
     theTrsf = e->invTrsf o theTrsf
   */

  if ( e->transformationName != (char*)NULL ) {
    if ( BAL_TransformationComposition( &(e->theTrsf), &theTrsf, &(e->theTrsf) ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
	std::cerr << proc << ": unable to compose transformation (2)" << std::endl;
      return( -1 );
    }
  } else {
    /* should allocate and copy
     */
     if ( _verbose_ )
       std::cerr << proc << ":weird, no transformation is given" << std::endl;
     return( -1 );
  }

  
  /* transformed list allocation
   */
  if ( e->trsfDetectionList.n_allocated < e->readDetectionList.n ) {
    BAL_FreeBlDetectionList( &(e->trsfDetectionList) );
    if ( BAL_AllocBlDetectionList( &(e->trsfDetectionList), e->readDetectionList.n ) != 0 ) {
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
	std::cerr << proc << ": unable to allocate detection list" << std::endl;
      return( -1 );
    }
  }

  for ( n=0; n<e->readDetectionList.n; n++ ) {
    e->trsfDetectionList.data[n] = e->readDetectionList.data[n];
    if ( BAL_TransformPoint( &(e->readDetectionList.data[n].voxelcenter), 
			     &(e->trsfDetectionList.data[n].voxelcenter), &(e->invTrsf) ) != 1 ) {
      BAL_FreeTransformation( &theTrsf );
      if ( _verbose_ )
	std::cerr << proc << ": unable to transform center #" << n << std::endl;
      return( -1 );
    }

  }
  e->trsfDetectionList.n = e->readDetectionList.n;
  
  e->trsfDetectionList.unit = REAL_UNIT;
  e->detectionList = &(e->trsfDetectionList);

  BAL_FreeTransformation( &theTrsf );

  return( 1 );
}










/************************************************************
 *
 * bal_blTrackerElementList related stuff
 *
 ************************************************************/

void BAL_InitBlTrackerElementList( bal_blTrackerElementList *l )
{
  l->data = (bal_blTrackerElement *)NULL;
  l->n = 0;
  l->n_allocated = 0;
}


int BAL_AllocBlTrackerElementList( bal_blTrackerElementList *l, int n )
{
  std::string proc = "BAL_AllocBlTrackerElementList";
  int i;
  l->data = (bal_blTrackerElement *)malloc( n * sizeof(bal_blTrackerElement) );
  if ( l->data == (bal_blTrackerElement *)NULL ) {
    if ( _verbose_ )
      std::cerr << proc << ": allocation failed" << std::endl;
    return( -1 );
  }
  l->n_allocated = n;
  l->n = 0;
  for ( i=0; i<n; i++ )
    BAL_InitBlTrackerElement( &(l->data[i]) );
  return( 1 );
}



void BAL_FreeBlTrackerElementList( bal_blTrackerElementList *l )
{
  int i;
  if ( l->data != (bal_blTrackerElement *)NULL ) {
    for ( i=0; i<l->n_allocated; i++ )
      BAL_InitBlTrackerElement( &(l->data[i]) );
    free( l->data );
  }
  BAL_InitBlTrackerElementList( l );
}





int BAL_ReadDetectionsBlTrackerElementList( bal_blTrackerElementList *l )
{
  std::string proc = "BAL_ReadDetectionsBlTrackerElementList";
  int i, j;
  
  for ( i=0; i<l->n; i++ ) {

    if ( BAL_ReadBlDetectionList( &(l->data[i].readDetectionList), l->data[i].detectionName ) != 0 ) {
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to read file '" << l->data[i].detectionName << "' for detection element #" << i << std::endl;
      for ( j=0; j<i; i++ )
	BAL_FreeBlDetectionList( &(l->data[i].readDetectionList) );
      return( -1 );
    }
    
    for ( j=0; j<l->data[i].readDetectionList.n; j++ ) 
      l->data[i].readDetectionList.data[j].imageindex =  l->data[i].index;

  }

  return( 1 );
}










/************************************************************
 *
 * bal_blTrackerParameter related stuff
 *
 ************************************************************/

void BAL_InitBlTrackerParameter( bal_blTrackerParameter *p )
{
  /* where to look for neighbors
   */
  p->halfneighborhood_unit.x = -1;
  p->halfneighborhood_unit.y = -1;
  p->halfneighborhood_unit.z = -1;

  p->halfneighborhood_voxel.x = -1;
  p->halfneighborhood_voxel.y = -1;
  p->halfneighborhood_voxel.z = -1;

  /* where to compute cost
   */
  p->costMin2Ddistance_unit = -1;
  p->costMin2Ddistance_voxel = -1;
  p->costLocation = _FROM_NEXT_;
  p->costCalculation = _EXACT_;

  
  p->minVesselness = 0.1;
  p->maxVesselness = 1.0;

  p->minTrackLength = 1;
  p->maxTracks = -1;


  /* how to search neighbors
   */
  p->searchType = _VOXEL_BRUTEFORCE_;
}












/************************************************************
 *
 * 
 *
 ************************************************************/

double _computeGeneralLineMeanIntensity( bal_doublePoint *pt1,
					 bal_doublePoint *pt2,
					 bal_image *theIm ) 
{
  std::string proc = "_computeGeneralLineMeanIntensity";
  double length, dl;
  int i, l;
  double v = 0.0;
  double dx, dy ,dz;

  /* points are supposed to be in voxel units
   */
  dz = pt2->z - pt1->z;
  dy = pt2->y - pt1->y;
  dx = pt2->x - pt1->x;

  length = sqrt( dx*dx + dy*dy + dz*dz );
  l = (int)length;


  if ( l < 1 ) {
    return( BAL_GetXYZvalue( theIm, pt1->x, pt1->y, pt1->z ) );
  }
  else if ( l < 2 ) {
    return( (BAL_GetXYZvalue( theIm, pt1->x, pt1->y, pt1->z ) + BAL_GetXYZvalue( theIm, pt2->x, pt2->y, pt2->z ))/2.0 );
  }
  else {
    dl = length / (double)l;
    dx /= length;
    dy /= length;
    dz /= length;
    v +=  (dl / 2.0) * BAL_GetXYZvalue( theIm, pt1->x, pt1->y, pt1->z );
    for ( i=1; i<l; i++ ) {
      v += dl * BAL_GetXYZvalue( theIm, pt1->x + i*dl*dx, pt1->y + i*dl*dy, pt1->z + i*dl*dz );
    }
    v +=  (dl / 2.0) * BAL_GetXYZvalue( theIm, pt2->x, pt2->y, pt2->z );
    return( v /length );
  }
}




double _computeApproximateLineMeanIntensity( bal_doublePoint *pt1,
					     bal_doublePoint *pt2,
					     bal_image *theIm ) 
{
  std::string proc = "_computeApproximateLineMeanIntensity";
  bal_integerPoint p1;
  bal_integerPoint p2;
  bal_integerPoint *line;
  int i, l;
  double s, v = 0.0;

  int xd, yd, zd;
  int x, y, z;
  int ax, ay, az;
  int sx, sy, sz;
  int dx, dy, dz;
  

  /* pick the closest points
   */
  x = p1.x = ( pt1->x < 0 ) ? (int)(pt1->x - 0.5) : (int)(pt1->x + 0.5);
  y = p1.y = ( pt1->y < 0 ) ? (int)(pt1->y - 0.5) : (int)(pt1->y + 0.5);
  z = p1.z = ( pt1->z < 0 ) ? (int)(pt1->z - 0.5) : (int)(pt1->z + 0.5);

  p2.x = ( pt2->x < 0 ) ? (int)(pt2->x - 0.5) : (int)(pt2->x + 0.5);
  p2.y = ( pt2->y < 0 ) ? (int)(pt2->y - 0.5) : (int)(pt2->y + 0.5);
  p2.z = ( pt2->z < 0 ) ? (int)(pt2->z - 0.5) : (int)(pt2->z + 0.5);


  /* find maximum of a and b */
#define MAX(a,b) (((a)>(b))?(a):(b))

  /* absolute value of a */
#define ABS(a) (((a)<0) ? -(a) : (a))

  /* take sign of a, either -1, 0, or 1 */
#define ZSGN(a) (((a)<0) ? -1 : (a)>0 ? 1 : 0)

  dx = p2.x - p1.x;
  dy = p2.y - p1.y;
  dz = p2.z - p1.z;
  
  ax = ABS(dx) << 1;
  ay = ABS(dy) << 1;
  az = ABS(dz) << 1;

  sx = ZSGN(dx);
  sy = ZSGN(dy);
  sz = ZSGN(dz);
  
  l = MAX( ABS(dx), MAX( ABS(dy), ABS(dz) ) );

  if ( l < 2 ) {

    v = 0.0;
    switch( theIm->type ) {
    default :
      if ( _verbose_ > 1 )
	std::cerr << proc << ": such image type not handled yet" << std::endl;
      return( -1 );
    case UCHAR : 
      {
	u8 ***buf = (u8***)theIm->array;
	if ( p1.x >= 0 && p1.x < (int)theIm->ncols
	     && p1.y >= 0 && p1.y < (int)theIm->nrows
	     && p1.z >= 0 && p1.z < (int)theIm->nplanes )
	  v += buf[p1.z][p1.y][p1.x];
	if ( l == 1 ) {
	  if ( p2.x >= 0 && p2.x < (int)theIm->ncols
	     && p2.y >= 0 && p2.y < (int)theIm->nrows
	     && p2.z >= 0 && p2.z < (int)theIm->nplanes )
	    v += buf[p2.z][p2.y][p2.x];
	  v /= 2.0;
	}
      }
      break;
    case FLOAT : 
      {
	r32 ***buf = (r32***)theIm->array;
	if ( p1.x >= 0 && p1.x < (int)theIm->ncols
	     && p1.y >= 0 && p1.y < (int)theIm->nrows
	     && p1.z >= 0 && p1.z < (int)theIm->nplanes )
	  v += buf[p1.z][p1.y][p1.x];
	if ( l == 1 ) {
	  if ( p2.x >= 0 && p2.x < (int)theIm->ncols
	     && p2.y >= 0 && p2.y < (int)theIm->nrows
	     && p2.z >= 0 && p2.z < (int)theIm->nplanes )
	    v += buf[p2.z][p2.y][p2.x];
	  v /= 2.0;
	}
      }
      break;
    }

  }
  else {

    line = new bal_integerPoint[l+1];

    if (ax >= MAX(ay, az)) {
      /* x dominant */
      yd = ay - (ax >> 1);
      zd = az - (ax >> 1);
      for ( i=0; x!=p2.x; i++ ) {
	line[i].x = x;
	line[i].y = y;
	line[i].z = z;
	if (yd >= 0) {
	  y += sy;
	  yd -= ax;
	}
	if (zd >= 0) {
	  z += sz;
	  zd -= ax;
	}
	x += sx;
	yd += ay;
	zd += az;
      }
    }
    else if (ay >= MAX(ax, az)) {
      /* y dominant */
      xd = ax - (ay >> 1);
      zd = az - (ay >> 1);
      for ( i=0; y!=p2.y; i++ ) {
	line[i].x = x;
	line[i].y = y;
	line[i].z = z;
	if (xd >= 0) {
	  x += sx;
	  xd -= ay;
	}
	if (zd >= 0) {
	  z += sz;
	  zd -= ay;
	}
	y += sy;
	xd += ax;
	zd += az;
      }
    }
    else if (az >= MAX(ax, ay)) {
      /* z dominant */
      xd = ax - (az >> 1);
      yd = ay - (az >> 1);
      for ( i=0; z!=p2.z; i++ ) {
	line[i].x = x;
	line[i].y = y;
	line[i].z = z;
	if (xd >= 0) {
	  x += sx;
	  xd -= az;
	}
	if (yd >= 0) {
	  y += sy;
	  yd -= az;
	}
	z += sz;
	xd += ax;
	yd += ay;
      }
    }
    else {
      delete[] line;
      if ( _verbose_ > 1 )
	std::cerr << proc << ": this should not occur" << std::endl;
      return( -1 );
    }
    
    line[l].x = p2.x;
    line[l].y = p2.y;
    line[l].z = p2.z;
    
    //std::cerr << proc << ": interval [(" << p1.x << "," << p1.y << "," << p1.z << "),(" << p2.x << "," << p2.y << "," << p2.z << ")]" << std::endl;
    for ( i=0; i<=l; i++ ) {
      //std::cerr << "   point #" << i << " = [" << line[i].x << "," << line[i].y << "," <<  line[i].z << "]" << std::endl;
    }
    

    s = v = 0.0;
    switch( theIm->type ) {
    default :
      if ( _verbose_ > 1 )
	std::cerr << proc << ": such image type not handled yet" << std::endl;
      return( -1 );
    case UCHAR : 
      {
	u8 ***buf = (u8***)theIm->array;
	if ( line[0].x >= 0 && line[0].x < (int)theIm->ncols
	     && line[0].y >= 0 && line[0].y < (int)theIm->nrows
	     && line[0].z >= 0 && line[0].z < (int)theIm->nplanes ) {
	  //std::cerr << " #0 :";
	  //std::cerr << " + 0.5 * " << buf[line[0].z][line[0].y][line[0].x] << " = ";
	  v += 0.5 * buf[line[0].z][line[0].y][line[0].x];
	  s += 0.5;
	  //std::cerr << v << std::endl;
	}
	for ( i=1; i<l; i++ ) {
	  if ( line[i].x >= 0 && line[i].x < (int)theIm->ncols
	       && line[i].y >= 0 && line[i].y < (int)theIm->nrows
	       && line[i].z >= 0 && line[i].z < (int)theIm->nplanes ) {
	    //std::cerr << " #" << i << " :";
	    //std::cerr << " + 1.0 * " << buf[line[i].z][line[i].y][line[i].x] << " = ";
	    v += buf[line[i].z][line[i].y][line[i].x];
	    s += 1.0;
	    //std::cerr << v << std::endl;
	  }
	}
	if ( line[l].x >= 0 && line[l].x < (int)theIm->ncols
	     && line[l].y >= 0 && line[l].y < (int)theIm->nrows
	     && line[l].z >= 0 && line[l].z < (int)theIm->nplanes ) {
	  //std::cerr << " #" << l << " :";
	  //std::cerr << " + 0.5 * " << buf[line[l].z][line[l].y][line[l].x] << " = ";
	  v += 0.5 * buf[line[l].z][line[l].y][line[0].x];
	  s += 0.5;
	  //std::cerr << v << std::endl;
	}
      }
      break;
    case FLOAT : 
      {
	r32 ***buf = (r32***)theIm->array;
	if ( line[0].x >= 0 && line[0].x < (int)theIm->ncols
	     && line[0].y >= 0 && line[0].y < (int)theIm->nrows
	     && line[0].z >= 0 && line[0].z < (int)theIm->nplanes ) {
	  //std::cerr << " #0 :";
	  //std::cerr << " + 0.5 * " << buf[line[0].z][line[0].y][line[0].x] << " = ";
	  v += 0.5 * buf[line[0].z][line[0].y][line[0].x];
	  s += 0.5;
	  //std::cerr << v << std::endl;
	}
	for ( i=1; i<l; i++ ) {
	  if ( line[i].x >= 0 && line[i].x < (int)theIm->ncols
	       && line[i].y >= 0 && line[i].y < (int)theIm->nrows
	       && line[i].z >= 0 && line[i].z < (int)theIm->nplanes ) {
	    //std::cerr << " #" << i << " :";
	    //std::cerr << " + 1.0 * " << buf[line[i].z][line[i].y][line[i].x] << " = ";
	    v += buf[line[i].z][line[i].y][line[i].x];
	    s += 1.0;
	    //std::cerr << v << std::endl;
	  }
	}
	if ( line[l].x >= 0 && line[l].x < (int)theIm->ncols
	     && line[l].y >= 0 && line[l].y < (int)theIm->nrows
	     && line[l].z >= 0 && line[l].z < (int)theIm->nplanes ) {
	  //std::cerr << " #" << l << " :";
	  //std::cerr << " + 0.5 * " << buf[line[l].z][line[l].y][line[l].x] << " = ";
	  v += 0.5 * buf[line[l].z][line[l].y][line[0].x];
	  s += 0.5;
	  //std::cerr << v << std::endl;
	}
      }
      break;
    }
    delete[] line;
    //std::cerr << " v / s = " << v << " / " << s << " = " << v/s << std::endl;
    if ( s > 0.0 ) v /= s;
  }

  return( v );
 
}










/* xxxxVoxelPoint corresponds to the point in voxel coordinates in xxxxImage,
   xxxxRealPoint corresponds to the point in real coordinates w.r.t a reference
   xxxxTrsf is the matrix that transform a point in real coordinates w.r.t a reference
       in a point in voxel coordinates

 * le cout de la connexion est calcule dans  blAssociationTracking/blAtCostAxons.cpp

   float blAtCostAxonVesselness::CalculateCost(blAtConnection* connection)

   si le deplacement est trop grand, on renvoie 1.0
   sinon on calcule le cout (vesselness moyen) dans l'image suivante et l'image precedente
   s'il n'y a pas de retractation, le cout est celui de l'image suivante
     (il y a une structure lineique qui est apparue)
   sinon c'est la valeur absolue de la difference (on aurait du prendre le max par coherence ?!)
   on divise par 255 pour normaliser entre 0 et 1 
   on renvoie 1.0 - la valeur normalisee

   float blAtCostAxonVesselness2::CalculateCost(blAtConnection* connection)

*/ 
double _computeCost( bal_blDetection *prevVoxelPoint,
		     bal_blDetection *prevRealPoint,
		     bal_transformation *prevTrsf,
		     bal_image *prevIm,
		     bal_blDetection *nextVoxelPoint,
		     bal_blDetection *nextRealPoint,
		     bal_transformation *nextTrsf,
		     bal_image *nextIm,
		     bal_blTrackerParameter *par )
{
  std::string proc = "_computeCost";
  double average = 0.0;
  double prevAverage = 0.0;
  double nextAverage = 0.0;
  bal_doublePoint pt;
  
  if ( par->costLocation == _FROM_PREV_ 
       || par->costLocation == _FROM_BOTH_ ) {
    
    if ( BAL_TransformPoint( &(nextRealPoint->voxelcenter), &pt, prevTrsf ) != 1 ) {
      if ( _verbose_ > 1 )
	std::cerr << proc << ": unable to transform point" << std::endl;
      return( -1 );
    }
    switch( par->costCalculation ) {
    default :
      if ( _verbose_ > 1 )
	std::cerr << proc << ": weird, this should not occur" << std::endl;
    case _EXACT_ :
      prevAverage = _computeGeneralLineMeanIntensity( &(prevVoxelPoint->voxelcenter), &pt, prevIm );
      break;
    case _APPROXIMATE_ :
      prevAverage = _computeApproximateLineMeanIntensity( &(prevVoxelPoint->voxelcenter), &pt, prevIm );
      break;
    }

  }
    
  if ( par->costLocation == _FROM_NEXT_ 
       || par->costLocation == _FROM_BOTH_ ) {
    
    if ( BAL_TransformPoint( &(prevRealPoint->voxelcenter), &pt, nextTrsf ) != 1 ) {
      if ( _verbose_ > 1 )
	std::cerr << proc << ": unable to transform point" << std::endl;
      return( -1 );
    }
    switch( par->costCalculation ) {
    default :
      if ( _verbose_ > 1 )
	std::cerr << proc << ": weird, this should not occur" << std::endl;
    case _EXACT_ :
      nextAverage = _computeGeneralLineMeanIntensity( &(nextVoxelPoint->voxelcenter), &pt, nextIm );
      break;
    case _APPROXIMATE_ :
      nextAverage = _computeApproximateLineMeanIntensity( &(nextVoxelPoint->voxelcenter), &pt, nextIm );
      break;
    }

  }

  switch( par->costLocation ) {
  case _FROM_PREV_ :
    average = prevAverage;
    break;
  case _FROM_NEXT_ :
    average = nextAverage;
    break;
  default :
  case _FROM_BOTH_ :
    average = ( prevAverage > nextAverage ) ? prevAverage : nextAverage;
    break;
  }
  
  if ( average > 1.0 ) {
    if ( _verbose_ > 1 )
      std::cerr << proc << ": weird, this should not occur, since vesselness is assumed to be normalized" << std::endl;
  }
  
  return( average );
}











/************************************************************
 *
 * building graph
 * blAssociationTracking/blAtTrackerGraphDetections.cpp
   optimisations
   void blAtTrackerGraphDetections::updateShortestPath()
   void blAtTrackerGraphDetections::updateMinCostFlow()
   construction du graphe

   void blAtTrackerGraphDetections::buildGraphNegative( lemon::ListDigraph &g, 
                                                        lemon::ListDigraph::ArcMap<int> &cap, 
                                                        lemon::ListDigraph::ArcMap<int> &cost, 
                                                        lemon::ListDigraph::Node &source, 
                                                        lemon::ListDigraph::Node &target, 
                                                        vector<lemon::ListDigraph::Arc> &arcs )

   void blAtTrackerGraphDetections::buildGraphNegativeOptimized( lemon::ListDigraph &g, 
                                                                 lemon::ListDigraph::ArcMap<int> &cap, 
                                                                 lemon::ListDigraph::ArcMap<int> &cost, 
                                                                 lemon::ListDigraph::Node &source, 
                                                                 lemon::ListDigraph::Node &target, 
                                                                 vector<lemon::ListDigraph::Arc> &arcs, 
                                                                 vector<vector<int> > &nodesConnectedArcsIds)

   void blAtTrackerGraphDetections::buildGraphNegative( lemon::ListDigraph &g, 
                                                        lemon::ListDigraph::ArcMap<int> &cap, 
                                                        lemon::ListDigraph::ArcMap<int> &cost, 
                                                        lemon::ListDigraph::Node &source, 
                                                        lemon::ListDigraph::Node &target, 
                                                        vector<lemon::ListDigraph::Arc> &arcs, 
                                                        vector< vector<int> > &arcsInOutMap )


   
 * blMpp/blMppOptimizerGraphcut.cpp
   - void blMppOptimizerGraphcut::Cut(int n)
     construit le graphe et optimise par graph cut (detection)
 * blMpp/blMppOptimizerMOGraphcut.cpp
 *
 ************************************************************/


// brute force
// blBuildTracks -detection DETECTIONS/axons_init%04d_detection.txt -vesselness VESSELNESS/axons_init_vesselness%04d.hdr -trsf TRSFS/res%04d.trsf -first 1 -last 29 -hnv 100 100 5
// blBuildTracks: elapsed time = 1.743614
// blBuildTracks -detection DETECTIONS/axons_init%04d_detection.txt -vesselness VESSELNESS/axons_init_vesselness%04d.hdr -trsf TRSFS/res%04d.trsf -first 1 -last 29 -hnv 10 10 5
// blBuildTracks: elapsed time = 1.592032

static int _buildLinksBetweenPair( typeGraph *mygraph,
				   lemon::ListDigraph::ArcMap<int> &cost,
				   bal_blTrackerElement *prev,
				   bal_blTrackerElement *next,
				   bal_blTrackerParameter *par )
{
  std::string proc = "_buildLinksBetweenPair";
  int i, n;
  bal_blPointerDetectionList neighbors;
  double dx, dy, d;
  double vesselness;
  double multiplierCoef = 1000.0;

  lemon::ListDigraph::Arc myarc;
  lemon::ListDigraph::Node mynode;

  BAL_InitBlPointerDetectionList( &neighbors );

  if ( par->searchType == _VOXEL_BUCKETS_ ) {
    if ( _debug_ >= 1 )
      std::cerr << proc << ": initialize search neighborhood" << std::endl;
    if (  BAL_InitNeighborsSearchBlDetectionList( next->detectionList, 
						  &(par->halfneighborhood_unit) ) != 1 ) {
      if ( _verbose_ >= 3 )
	std::cerr << proc << ": unable to initialize search neighborhood" << std::endl;
      return( -1 );
    }
  }


  /* loop on the detection on the previous frame
     detectionList gives point coordinates in real units 
     and in 'prev' geometry after registration

     to compute vesselness, points have to be transformed into prev image geometry
   */


  for ( i=0; i<prev->detectionList->n; i++ ) {

    neighbors.n = 0;

    if ( BAL_SearchBlDetectionNeighbors( &(prev->detectionList->data[i]), 
					 next->detectionList,
					 &neighbors, &(par->halfneighborhood_unit), 
					 par->searchType ) < 0 ) {
      BAL_FreeBlPointerDetectionList( &neighbors );
      if ( _verbose_ > 1 )
	std::cerr << proc << ": unable to search neighbors for detection element #" << i << std::endl;
      return( -1 );
    }
    
    /* loop on the neighbors
     */
    for ( n=0; n<neighbors.n; n++ ) {


      /* le cout (ie vesselness) est entre 0 (mauvais) et 1 (bon)
	 - a priori Sylvain renvoyait 1 - vesselness : entre 0 (bon) et 1 (mauvais)
	 - mais utilisait (1 - vesselness) - 1.0, soit (- vesselness) : entre -1 (bon) et 0 (mauvais)

       */
      dx = prev->detectionList->data[i].voxelcenter.x - neighbors.data[n]->voxelcenter.x;
      dy = prev->detectionList->data[i].voxelcenter.y - neighbors.data[n]->voxelcenter.y;
      d = sqrt( dx*dx + dy+dy );
      
      if ( d <= par->costMin2Ddistance_unit ) {
	vesselness = 1.0;
      }
      else {
	vesselness = _computeCost( &(prev->readDetectionList.data[i]),
				   &(prev->detectionList->data[i]),
				   &(prev->theTrsf),
				   &(prev->vesselnessImage),
				   &(next->readDetectionList.data[ neighbors.data[n]->identifier ]),
				   neighbors.data[n],
				   &(next->theTrsf),
				   &(next->vesselnessImage),
				   par );
      }

      if ( par->minVesselness < vesselness &&  vesselness < par->maxVesselness ) {
	myarc = mygraph->graph.addArc( mygraph->graph.nodeFromId( prev->readDetectionList.data[i].idnode ),
				       mygraph->graph.nodeFromId( neighbors.data[n]->idnode ) );
	/* le vesselness est dans [0,1], pour le traduire en int,
	   on le multiplie par 1000
	*/
	cost[ myarc ] = (int)( (-1.0) * vesselness * multiplierCoef );
	// std::cout << " " << cost[ myarc ] << " ";

      }

    }

  }

  

  BAL_FreeBlPointerDetectionList( &neighbors );


  return( 1 );
}















int BAL_BuildBlTrackerGraph( typeGraph *mygraph, 
			     lemon::ListDigraph::ArcMap<int> &cost,
			     bal_blTrackerElementList *theList, 
			     bal_blTrackerParameter *par  )
{
  std::string proc = "BAL_BuildBlTrackerGraph";
  int i, j;
  int nnodes;

  lemon::ListDigraph::Node mynode;
  lemon::ListDigraph::Arc myarc;

  bal_blTrackerElement *previous = (bal_blTrackerElement *)NULL;
  bal_blTrackerElement *next = (bal_blTrackerElement *)NULL;


  
  /****************************************
   *
   * NODES
   * the detection files are read
   *
   ****************************************/

  if ( BAL_ReadDetectionsBlTrackerElementList( theList ) != 1 ) {
    if ( _verbose_ ) 
      std::cerr << proc << ": unable to read detection lists" << std::endl;
    return( -1 );
  }

  /****************************************
   *
   * NODES
   *
   ****************************************/
  /* nodes are
     - source and sink
     - all detections
  */
  for ( i=0, nnodes=0; i<theList->n; i++ )
    nnodes += theList->data[i].readDetectionList.n;
  if ( _verbose_ >= 2 )
    std::cout << " ... there are " << nnodes << " detections" << std::endl;

  mygraph->graph.reserveNode(nnodes+2);

  mygraph->source = mygraph->graph.addNode();
  mygraph->target = mygraph->graph.addNode();

  for ( i=0, nnodes=0; i<theList->n; i++ ) {
    for ( j=0; j<theList->data[i].readDetectionList.n; j++ ) {
      mynode = mygraph->graph.addNode();
      theList->data[i].readDetectionList.data[j].idnode = mygraph->graph.id( mynode );
      // std::cout << " id=" << theList->data[i].readDetectionList.data[j].idnode;
      myarc = mygraph->graph.addArc( mygraph->source, mynode );
      cost[ myarc ] = 0;
      myarc = mygraph->graph.addArc( mynode, mygraph->target );
      cost[ myarc ] = 0;
    }
  }

  if ( _verbose_ >= 2 ) {
    std::cout << " ... max node id is " << mygraph->graph.maxNodeId() << std::endl;
    std::cout << "     suggesting " << mygraph->graph.maxNodeId() + 1 << " nodes" << std::endl;
    std::cout << "     #nodes is  " << lemon::countNodes( mygraph->graph ) << std::endl;
    std::cout << "     #arcs  is  " << lemon::countArcs( mygraph->graph ) << std::endl;
  }

  /****************************************
   *
   * LINKS
   *
   ****************************************/

  /* the detection list has already been read
     - read now the other element : transformation, vesselness
  */
  next = &(theList->data[0]);
  if ( BAL_PartialReadDetectionsBlTrackerElement( next ) != 1 ) {
    if ( _verbose_ ) 
      std::cerr << proc << ": unable to read element tracker #0" << std::endl;
    return( -1 );
  }
  /* transform detection list into real units
   */
  if ( BAL_UpdateDetectionsBlTrackerElement( next ) != 1 ) {
    if ( _verbose_ ) 
      std::cerr << proc << ": unable to update element tracker #" << i-1 << std::endl;
    return( -1 );
  }
  


  /* detections are in voxel units
     update search neighborhood
     need at least one image to do that
   */
  if ( par->halfneighborhood_voxel.x > 0.0 && par->halfneighborhood_voxel.y > 0.0 && par->halfneighborhood_voxel.z > 0.0 
       && par->halfneighborhood_unit.x <=  0.0 && par->halfneighborhood_unit.y <= 0.0 && par->halfneighborhood_unit.z <= 0.0) {
    if ( next->vesselnessImage.vx > 0.0 
	 && next->vesselnessImage.vy > 0.0 
	 && next->vesselnessImage.vz > 0.0 ) {
      par->halfneighborhood_unit.x = par->halfneighborhood_voxel.x * next->vesselnessImage.vx;
      par->halfneighborhood_unit.y = par->halfneighborhood_voxel.y * next->vesselnessImage.vy;
      par->halfneighborhood_unit.z = par->halfneighborhood_voxel.z * next->vesselnessImage.vz;
    }
    else {
      par->halfneighborhood_unit.x = par->halfneighborhood_voxel.x;
      par->halfneighborhood_unit.y = par->halfneighborhood_voxel.y;
      par->halfneighborhood_unit.z = par->halfneighborhood_voxel.z;
    }
  }
  else if ( par->halfneighborhood_unit.x > 0.0 && par->halfneighborhood_unit.y > 0.0 && par->halfneighborhood_unit.z > 0.0 
	    && par->halfneighborhood_voxel.x <= 0.0 && par->halfneighborhood_voxel.y <= 0.0 && par->halfneighborhood_voxel.z <= 0.0 ) {
    if ( next->vesselnessImage.vx > 0.0 
	 && next->vesselnessImage.vy > 0.0 
	 && next->vesselnessImage.vz > 0.0 ) {
      par->halfneighborhood_voxel.x = par->halfneighborhood_unit.x / next->vesselnessImage.vx;
      par->halfneighborhood_voxel.y = par->halfneighborhood_unit.y / next->vesselnessImage.vy;
      par->halfneighborhood_voxel.z = par->halfneighborhood_unit.z / next->vesselnessImage.vz;
    }
    else {
      par->halfneighborhood_voxel.x = par->halfneighborhood_unit.x;
      par->halfneighborhood_voxel.y = par->halfneighborhood_unit.y;
      par->halfneighborhood_voxel.z = par->halfneighborhood_unit.z;
    }
  }
  else {
    BAL_FreeBlTrackerElement( next );
    if ( _verbose_ ) 
      std::cerr << proc << ": weird, both neighborhoods (unit and voxel) are either initialized or not" << std::endl;
    return( -1 );
  }


  if ( par->costMin2Ddistance_voxel <= 0.0 && par->costMin2Ddistance_unit <= 0.0 ) {
    par->costMin2Ddistance_voxel = 1.0;
  }

  if ( par->costMin2Ddistance_voxel > 0.0 && par->costMin2Ddistance_unit <= 0.0 ) {
    if ( next->vesselnessImage.vx > 0.0 ) {
      par->costMin2Ddistance_unit = par->costMin2Ddistance_voxel * next->vesselnessImage.vx;
    }
    else {
      par->costMin2Ddistance_unit = par->costMin2Ddistance_voxel;
    }
  }
  else if ( par->costMin2Ddistance_unit > 0.0 && par->costMin2Ddistance_voxel <= 0.0 ) {
    if ( next->vesselnessImage.vx > 0.0 ) {
      par->costMin2Ddistance_voxel = par->costMin2Ddistance_unit / next->vesselnessImage.vx;
    }
    else {
      par->costMin2Ddistance_voxel = par->costMin2Ddistance_unit;
    }
  }





  /* loop on pairs of consecutive frames
   */
  
  for ( i=1; i<theList->n; i++ ) {

    if ( _verbose_ >= 2 ) {
      std::cout << proc << ": process element tracker pair [" << i-1 << "," << i << "]" << std::endl;
    }
    else if ( _verbose_ >= 1 ) {
      std::cerr << ".";
    }

    
    previous = next;

    /* the detection list has already been read
       - read now the other element : transformation, vesselness
       - transform the detections into real unit, and in frame #i-1 coordinates
    */
    next = &(theList->data[i]);
    if ( BAL_PartialReadDetectionsBlTrackerElement( next ) != 1 ) {
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to read element tracker #" << i << std::endl;
      return( -1 );
    }
    if ( BAL_UpdateDetectionsBlTrackerElement( next ) != 1 ) {
      if ( _verbose_ ) 
	std::cerr << proc << ": unable to update element tracker #" << i << std::endl;
      return( -1 );
    }




    /* here we can connect 'previous' and 'next'
     */
    /* build pairings between #i and #(i+1)
     */
    if (  _buildLinksBetweenPair( mygraph, cost, previous, next, par  ) != 1 ) {
      if ( _verbose_ )
	std::cerr << proc << ": unable to build links between sets #" << i << " and #" << i+1 << std::endl;
      return( -1 );
    }


    if ( _verbose_ >= 3 ) {
      std::cerr << "     #nodes is  " << lemon::countNodes( mygraph->graph ) << std::endl;
      std::cerr << "     #arcs  is  " << lemon::countArcs( mygraph->graph ) << std::endl;
    }
  
    if ( previous->vesselnessName != (char*) NULL ) {
      BAL_FreeImage ( &(previous->vesselnessImage) );
    }
  }

  if ( next->vesselnessName != (char*) NULL ) {
    BAL_FreeImage ( &(next->vesselnessImage) );
  }

  
  /****************************************
   *
   * SOME FIGURES
   *
   ****************************************/
  if ( _verbose_ < 2 && _verbose_ >= 1 ) {
      std::cerr <<  std::endl;
    }

  if ( _verbose_ >=  2 ) {
    std::cout << std::endl;
    std::cout << " ... max node id is " << mygraph->graph.maxNodeId() << std::endl;
    std::cout << "     max arc  id is " << mygraph->graph.maxArcId() << std::endl;
  }
  else if ( _verbose_ >=  1 ) {
    std::cout << "     #nodes is  " << lemon::countNodes( mygraph->graph ) << std::endl;
    std::cout << "     #arcs  is  " << lemon::countArcs( mygraph->graph ) << std::endl;
  }



  return( 1 );
}








int _addBlDetectionToList( int idnode,
			   bal_blDetection *resDetection,
			   bal_blTrackerElementList *thelist )
{
  std::string proc = "_addBlDetectionToList";
  int firstlistid;
  int lastlistid;
  int l;
  bal_blDetectionList *trsfDetectionList;

  firstlistid = 0;
  lastlistid = 1;
  for ( l=0; l<thelist->n && idnode > lastlistid; l++ ) {
    firstlistid = lastlistid + 1;
    lastlistid = firstlistid + thelist->data[l].readDetectionList.n - 1;
  }
  l--;

  trsfDetectionList = &(thelist->data[l].trsfDetectionList);
  
  if ( trsfDetectionList->data[ idnode - firstlistid ].idnode != idnode ) {
    if ( _verbose_ ) {
      std::cerr << proc << ": weird, found point has id=" << trsfDetectionList->data[ idnode - firstlistid ].idnode << std::endl;
      std::cerr << "      while searched id was " << idnode << std::endl;
      return( -1 );
    }
  }

  /* point is transformed into voxel coordinates with respect to the
     first image
  */

  *resDetection = trsfDetectionList->data[ idnode - firstlistid ];
  if ( BAL_TransformPoint( &(trsfDetectionList->data[ idnode - firstlistid ].voxelcenter), 
			   &(resDetection->voxelcenter), 
			   &(thelist->data[0].theTrsf) ) != 1 ) {
    if ( _verbose_ ) {
      std::cerr << proc << ": unable to transform point" << std::endl;
      return( -1 );
    }
  }

  return( 1 );
}










int _addBlTrack( typeGraph *mygraph, 
		 lemon::Path<lemon::ListDigraph> &path,
		 bal_blTrackList *reslist,
		 bal_blTrackerElementList *thelist )
{
  std::string proc = "_addBlTrack";
  int narcs;
  bal_blTrack t;
  lemon::ListDigraph::Node fnode, lnode;  


  BAL_InitBlTrack( &t );
  
  narcs = 0;
  for ( lemon::Path<lemon::ListDigraph>::ArcIt it( path ); it != lemon::INVALID; ++it ) {
    narcs ++;
  }

  if ( _verbose_ >= 2 )
    std::cout << proc << ": add path #" << reslist->n << " of length " << narcs << std::endl;
  
  if ( narcs <= 2 ) {
    if ( _verbose_ >= 3 )
      std::cout << proc << ": too small path, #arcs=" << narcs << std::endl;
    return( -1 );
  }
  
  if ( BAL_AllocBlDetectionList( &(t.detectionList), narcs-1 ) != 0 ) {
    if ( _verbose_ )
      std::cerr << proc << ": error when allocating track" << std::endl;
    return( -1 );
  }


  /* add the arc to the tracks
   */
  for ( lemon::Path<lemon::ListDigraph>::ArcIt it( path ); it != lemon::INVALID; ++it ) {
    /* arc extremities
     */
    fnode = mygraph->graph.source( it );
    lnode = mygraph->graph.target( it );
    if ( mygraph->graph.id( fnode ) == mygraph->graph.id( mygraph->source ) )
      continue;
    if ( _addBlDetectionToList( mygraph->graph.id( fnode ), &(t.detectionList.data[ t.detectionList.n ]), thelist ) != 1 ) {
      BAL_FreeBlTrack( &t );
      if ( _verbose_ )
	std::cerr << proc << ": error when adding detection to list" << std::endl;
      return( -1 );
    }
    t.detectionList.n ++;
  }
  
  if ( BAL_AddBlTrackToList( &t, reslist ) != 0 ) {
    BAL_FreeBlTrack( &t );
    if ( _verbose_ )
	std::cerr << proc << ": error when adding track to list" << std::endl;
    return( -1 );
  }
  
 

  return( 1 );
}










int _ExtractBlTracks( typeGraph *mygraph, 
		      lemon::ListDigraph::ArcMap<int> &cost,
		      bal_blTrackList *reslist,
		      bal_blTrackerElementList *thelist, 
		      bal_blTrackerParameter *par )
{
  std::string proc = "_ExtractBlTracks";
  lemon::BellmanFord<lemon::ListDigraph, lemon::ListDigraph::ArcMap<int> > bf( mygraph->graph , cost );
  lemon::Path<lemon::ListDigraph> path;
  lemon::ListDigraph::Node fnode, lnode;  
  
  int narcs = 0;
  int stop = 0;
  int iterations = 0;
  int ind=0;
  int indicator = 10;
  int indicatorlength = 60;  


  do {

    if (  _verbose_ < 2  && _verbose_ >= 1 ) {
      if ( iterations % indicator == 0 ) {
	std::cerr << ".";
	if ( narcs > 0 && narcs < 20 && indicator < 20 ) 
	  indicator *= 2;
	else if ( narcs > 0 && narcs < 10 && indicator < 50 ) 
	  indicator *= 2;
	else if ( narcs > 0 && narcs < 5 && indicator < 100 ) 
	  indicator *= 2;
	else if ( narcs > 0 && narcs < 5 && indicator < 500 ) 
	  indicator *= 2;
	ind ++;
	if ( ind == indicatorlength ) {
	  std::cerr << std::endl;
	  std::cerr << "     #nodes = " << lemon::countNodes( mygraph->graph );
	  std::cerr << "     #arcs = " << lemon::countArcs( mygraph->graph ) << std::endl;
          ind = 0;
	}
      }
    }

    
    bf.run( mygraph->source );
    path = bf.path( mygraph->target );

    narcs = 0;
    for ( lemon::Path<lemon::ListDigraph>::ArcIt it( path ); it != lemon::INVALID; ++it ) {
      narcs ++;
    }

    /* for a path with narcs arcs, there are narcs+1 nodes, including 
       source and target, then narcs-1 from the data
    */
    if ( narcs >=3 && narcs >= par->minTrackLength+1 ) {

      if ( 0 && _verbose_ >= 2 ) {
	std::cout << "add track #" << reslist->n << " of length " << narcs-2 << std::endl;
      }
      
      if ( _addBlTrack( mygraph, path, reslist, thelist ) != 1 ) {
	if ( _verbose_ )
	  std::cerr << proc << ": error when adding path" << std::endl;
	stop = 1;
      }
    }
      
    /* remove the nodes
     */
    for ( lemon::Path<lemon::ListDigraph>::ArcIt it( path ); it != lemon::INVALID; ++it ) {
      /* arc extremities
       */
      fnode = mygraph->graph.source( it );
      lnode = mygraph->graph.target( it );
      if ( mygraph->graph.id( fnode ) == mygraph->graph.id( mygraph->source ) )
	continue;
      mygraph->graph.erase( fnode );
    }
    
    if ( _verbose_ >= 3 ) {
      std::cerr << "     #nodes = " << lemon::countNodes( mygraph->graph );
      std::cerr << "     #arcs = " << lemon::countArcs( mygraph->graph ) << std::endl;
    }

    /* empty graph with only source and target
     */
    if ( lemon::countNodes( mygraph->graph ) <= 2 || lemon::countArcs( mygraph->graph ) == 0 ) 
      stop = 1;
    
    iterations ++;

  } while ( stop != 1 
	    && ( reslist->n < par->maxTracks || par->maxTracks < 0 ) );
  
  if (  _verbose_ < 2  && _verbose_ >= 1 )
    std::cerr << std::endl;


  return( 1 );
}










int BAL_ExtractBlTracks( bal_blTrackList *reslist,
			 bal_blTrackerElementList *thelist, 
			 bal_blTrackerParameter *par  )
{
  std::string proc = "BAL_ExtractBlTracks";
  
  typeGraph mygraph;
  lemon::ListDigraph::ArcMap<int> cost( mygraph.graph );

  

  if ( _verbose_ ) {
     std::cout << proc << ": build the graph" << std::endl;
  }

  if ( BAL_BuildBlTrackerGraph( &mygraph, cost, thelist, par ) != 1 ) {
    if ( _verbose_ )
      std::cerr << proc << ": unable to build graph" << std::endl;
    return( -1 );
  }
  
  if ( _verbose_ ) {
     std::cout << proc << ": extract the tracks" << std::endl;
  }


  if ( _ExtractBlTracks(  &mygraph, cost, reslist, thelist, par ) != 1 ) {
    if ( _verbose_ )
      std::cerr << proc << ": unable to extract tracks from the graph" << std::endl;
    return( -1 );
  }
  
  if ( _verbose_ ) {
    std::cout << proc << ": found " << reslist->n << " tracks" << std::endl;
  }

		
  return( 1 );
}
