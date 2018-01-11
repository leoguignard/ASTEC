/****************************************************
 * parcelling.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2008
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Tue Apr 15 17:41:49 CEST 2008
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <time.h>
#include <stdlib.h>

#include <convert.h>
#include <parcelling.h>


/*
 * static global variables
 * verbose,
 * management of ambiguous cases
 * memory management
 * 
 */





static int _verbose_ = 1;

void parcelling_setnoverbose()
{
  _verbose_ = 0;
}

void parcelling_setverbose()
{
  if ( _verbose_ < 1 ) _verbose_ = 1;
  else _verbose_ ++;
}





static long int seedRandom = -1;

void parcelling_setRandomSeed( long int seed )
{
  seedRandom = seed;
}





static int _nPointsToBeAllocated_ = 100;

void parcelling_setNumberOfPointsForAllocation( int n )
{
   if ( n > 0 ) _nPointsToBeAllocated_ = n;
}





static int _max_iterations_ = -1;

void parcelling_setNumberOfIterations( int n )
{
  _max_iterations_ = n;
}





static int _force_exact_center_calculation_ = 0;

void parcelling_ForceExactCenterCalculation()
{
  _force_exact_center_calculation_ = 1;
}

void parcelling_DoNotForceExactCenterCalculation()
{
  _force_exact_center_calculation_ = 0;
}





/*
 *
 *
 *
 */

typedef struct {
  int x;
  int y;
  int z;
} typeParcelCenter;


typedef struct {
  int xmin, xmax;
  int ymin, ymax;
  int zmin, zmax;
  int n;
  double sx, sy, sz;
  int x, y, z;
} typeParcel;


typedef struct {
  int x;
  int y;
  int z;
  int i;
  int l;
  int d;
  int inside;
} typePoint;

typedef struct {
  int nPoints;
  int nAllocatedPoints;
  typePoint *pt;
} typePointList;





/*
 *
 * some static functions
 *
 */

static int _get_centers_from_parcels( u16* theDist, u16* theLabels, int *theDim, 
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] );

static int _get_centers_from_parcels_u32( u32* theDist, u32* theLabels, int *theDim,
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] );

static int _get_parcels_from_centers( u16* theDist, u16* theLabels, int *theDim, 
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] );

static int _get_parcels_from_centers_u32( u32* theDist, u32* theLabels, int *theDim,
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] );

static int _remove_first_points_from_list( typePointList *thePointList, int npts );

static int _get_parcel_neighbors( typePointList *thePointList, int npts,
                                  u16* theLabels, int *theDim,
                                  int INLIST, int NONINLIST );

static int _get_parcel_neighbors_u32( typePointList *thePointList, int npts,
                                  u32* theLabels, int *theDim,
                                  int INLIST, int NONINLIST );

static void _update_distances_in_list( typePointList *thePointList, int inc[3][3][3],
                                       int MAXDIST,
                                       u16* theDist, u16* theLabels, int *theDim,
                                       int nparcels );

static void _update_distances_in_list_u32( typePointList *thePointList, int inc[3][3][3],
                                       int MAXDIST,
                                       u32* theDist, u32* theLabels, int *theDim,
                                       int nparcels );

static int _get_min_distance( typePointList *thePointList, int *nb );

static int _get_initial_parcel_centers( u16* theLabels, int *theDim, 
					typeParcelCenter *center, int nparcels );

static int _get_initial_parcel_centers_u32( u32* theLabels, int *theDim,
                                        typeParcelCenter *center, int nparcels );

static int _count_points( void *theBuf, bufferType theBufType, int *theDim );
static int _check_free_borders( void *theBuf, bufferType theBufType, 
				int *theDim );



static void _init_point_list( typePointList *l );
static void _free_point_list( typePointList *l );
static int _add_point_to_list( typePointList *l, 
			       int x, int y, int z, int i, 
			       int label, int dist, int *theDim );




 
/************************************************************
 *
 *
 *
 ************************************************************/
int parcelling( void *theBuf, bufferType theBufType, 
		int **theSeeds, int nparcels,
		void *theOutputLabels, bufferType theOutputLabelsType, 
		void *theOutputDistance, bufferType theOutputDistanceType,
		int *theDim, int inSeeds )
{
  char *proc = "parcelling";
  typeParcelCenter *centers;
  typeParcelCenter *theCenter, *oldCenter;
  u16 *theLabels, *theDist;
  u32 *theLabels_u32, *theDist_u32;
  int flag_32 = (nparcels >= 65535) ? 1 : 0;
  int i, v=theDim[2]*theDim[1]*theDim[0];
  int iteration, change=0;

  int inc[3][3][3] = {{{5,4,5},{4,3,4},{5,4,5}},{{4,3,4},{3,0,3},{4,3,4}},{{5,4,5},{4,3,4},{5,4,5}}};

  if ( 0 && _check_free_borders( theBuf, theBufType, theDim ) != 1 ) {
    /* il faudrait en fait encapsuler cette fonction et l'appeler
       avec des images plus grandes si necessaire
    */
    if ( _verbose_ )
      fprintf( stderr, "%s: non-empty borders not handled yet\n", proc );
    return( -1 );
  }



  /* some allocations
     - 2 arrays of parcel centers
     - 2 data arrays
   */

  centers = (typeParcelCenter *)malloc( 2*nparcels*sizeof( typeParcelCenter ) );
  if ( centers == (typeParcelCenter *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate parcel center list\n", proc );
    return( -1 );
  }
  theCenter = centers;
  oldCenter = centers + nparcels;
  
  if (flag_32)
  {
    theLabels_u32 = (u32*)malloc( 2*theDim[2]*theDim[1]*theDim[0]*sizeof(u32) );
    if ( theLabels_u32 == (u32*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary array\n", proc );
      free( centers );
      return( -1 );
    }
    theDist_u32 = theLabels_u32 + theDim[2]*theDim[1]*theDim[0];
  }
  else
  {
    theLabels = (u16*)malloc( 2*theDim[2]*theDim[1]*theDim[0]*sizeof(u16) );
    if ( theLabels == (u16*)NULL ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to allocate auxiliary array\n", proc );
      free( centers );
      return( -1 );
    }
    theDist = theLabels + theDim[2]*theDim[1]*theDim[0];
  }




  /* put initial data into data array
  */
  if (flag_32)
  {
    if ( ConvertBuffer( theBuf, theBufType, theLabels_u32, UINT,  theDim[2]*theDim[1]*theDim[0] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to convert buffer\n", proc );
      free( theLabels_u32 );
      free( centers );
      return( -1 );
    }
  }
  else
  {
    if ( ConvertBuffer( theBuf, theBufType, theLabels, USHORT,  theDim[2]*theDim[1]*theDim[0] ) != 1 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to convert buffer\n", proc );
      free( theLabels );
      free( centers );
      return( -1 );
    }
  }

  
  /* draw random initial seeds (if necessary)
     (assume that foreground points are non-null points)
  */

  if ( theSeeds != (int**) NULL && inSeeds == 1 ) {
    for ( i=0; i<nparcels; i++ ) {
      theCenter[i].x = theSeeds[i][0];
      theCenter[i].y = theSeeds[i][1];
      if ( theDim[2] > 1 )
	theCenter[i].z = theSeeds[i][2];
      else 
	theCenter[i].z = 0;
    }
  }
  else {
    if (flag_32)
    {
      if ( _get_initial_parcel_centers_u32( theLabels_u32, theDim, theCenter, nparcels ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get initial centers\n", proc );
        free( theLabels_u32 );
        free( centers );
        return( -1 );
      }
    }
    else
    {
      if ( _get_initial_parcel_centers( theLabels, theDim, theCenter, nparcels ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get initial centers\n", proc );
        free( theLabels );
        free( centers );
        return( -1 );
      }
    }
  }

 
  


  iteration = 0;
  do {

    /* compute geodesic influence areas of parcels
     */
    if (flag_32)
    {
      if ( _get_parcels_from_centers_u32( theDist_u32, theLabels_u32, theDim,
                                      theCenter, nparcels, inc ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get parcels from centers\n", proc );
        free( theLabels_u32 );
        free( centers );
        return( -1 );
      }

    }
    else
    {
      if ( _get_parcels_from_centers( theDist, theLabels, theDim,
                                      theCenter, nparcels, inc ) != 1 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: unable to get parcels from centers\n", proc );
        free( theLabels );
        free( centers );
        return( -1 );
      }
    }



    /* if some iterations remains to do 
       -> update centers
    */

    if ( (_max_iterations_ < 0) ||
	 ( _max_iterations_ > 0 && iteration < _max_iterations_ ) ) {

      /* copy of parcel centers
       */
      memcpy( (void*)oldCenter, (void*)theCenter, nparcels*sizeof( typeParcelCenter ) );
      
      if (flag_32)
      {
        if ( _get_centers_from_parcels_u32( theDist_u32, theLabels_u32, theDim,
                                        theCenter, nparcels, inc ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get centers from parcels\n", proc );
          free( theLabels_u32 );
          free( centers );
          return( -1 );
        }
      }
      else
      {
        if ( _get_centers_from_parcels( theDist, theLabels, theDim,
                                        theCenter, nparcels, inc ) != 1 ) {
          if ( _verbose_ )
            fprintf( stderr, "%s: unable to get centers from parcels\n", proc );
          free( theLabels );
          free( centers );
          return( -1 );
        }
      }
      
      change = 0;
      for ( i=0; i<nparcels; i++ ) {
	if ( theCenter[i].x != oldCenter[i].x 
	     || theCenter[i].y != oldCenter[i].y 
	   || theCenter[i].z != oldCenter[i].z )
	change ++;
      }
    
      if ( _verbose_ ) {
	fprintf( stderr, "%s: iteration #%2d, %3d changes\n", proc, iteration, change );
      }
    }


    iteration ++;


    /* stop if 
       _max_iterations_ == 0 (we just want geodesic influence areas)
	OR (maximal number of iterations is positive and is reached)
	
	OR (there are no changes)
    */
  } while ( !( (_max_iterations_ == 0)
	       || (_max_iterations_ > 0 && iteration > _max_iterations_) 
	       || (change == 0) ) );
    

  
  /* export results :
     - parcel centers
     - labels 
     - distance
  */
  

  if ( theSeeds != (int**) NULL ) {
    for ( i=0; i<nparcels; i++ ) {
      theSeeds[i][0] = theCenter[i].x;
      theSeeds[i][1] = theCenter[i].y;
      if ( theDim[2] > 1 )
	theSeeds[i][2] = theCenter[i].z;
      else
	theSeeds[i][2] = 0;
    }
  }

  free( centers );


  
  /* copy auxiliary label array
   */
  if ( theOutputLabels != NULL ) {
    switch( theOutputLabelsType ) {
    default :
      if ( _verbose_ ) {
        fprintf( stderr, "%s: such label type not handled yet\n", proc );
      }
      free( theLabels );
      return( -1 );
    case UCHAR :
      {
        u8 *buf = theOutputLabels;
        for ( i=0; i<v; i++ )
	  buf[i] = theLabels[i];
      }
      break;
    case USHORT :
      {
        u16 *buf = theOutputLabels;
        for ( i=0; i<v; i++ )
          buf[i] = theLabels[i];
      }
      break;
    case UINT :
      {
        u32 *buf = theOutputLabels;
        for ( i=0; i<v; i++ )
          buf[i] = theLabels[i];
      }
      break;
    }
  }
  

  if ( theOutputDistance != NULL ) {
    switch( theOutputDistanceType ) {
    default :
      if ( _verbose_ ) {
	fprintf( stderr, "%s: such distance type not handled yet\n", proc );
      }
      free( theLabels );
      return( -1 );
    case UCHAR :
      {
	u8 *buf = theOutputDistance;
	for ( i=0; i<v; i++ )
	  buf[i] = theDist[i];
      }
      break;
    case USHORT :
      {
        u16 *buf = theOutputDistance;
        for ( i=0; i<v; i++ )
          buf[i] = theDist[i];
      }
      break;
    case UINT :
      {
        u32 *buf = theOutputDistance;
        for ( i=0; i<v; i++ )
          buf[i] = theDist[i];
      }
      break;
    }
  }


  if (flag_32)
    free( theLabels_u32 );
  else
    free( theLabels );

  return( 1 );
}










/*************************************************************
 *
 * some static functions
 *
 *************************************************************/

static void _print_neighorhood( int x, int y, int z, void *buf, bufferType type, int *theDim, char *str )
{
  int i, j, k;

  fprintf( stdout, "Voisinage:" );
  if ( str != (char*)NULL ) fprintf( stdout, "%s", str );
  fprintf( stdout, "\n" );

  switch ( type ) {
  default :
    fprintf( stdout, "type not handled yet\n" );
    break;

  case DOUBLE :
    {
      double *theBuf = (double*)buf;
      for ( j=-1; j<=1; j++ ) {
	for ( k=-1; k<=1; k++ ) { 
	  for ( i=-1; i<=1; i++ ) {
	    if ( x+i < 0 || x+i >= theDim[0] 
		 || y+j < 0 || y+j >= theDim[1] 
		 || z+k < 0 || z+k >= theDim[2] )
	      fprintf( stdout, " ....." );
	    else 
	      fprintf( stdout, " %9f", theBuf[ ((z+k)*theDim[1]+(y+j))*theDim[0]+(x+i) ] );
	}
	  if ( k < 1 ) fprintf( stdout, "  - " );
	}
	fprintf( stdout, "\n" );
      }
      fprintf( stdout, "\n" );
    }
    break;

  case USHORT :
    {
      u16 *theBuf = (u16*)buf;
      for ( j=-1; j<=1; j++ ) {
        for ( k=-1; k<=1; k++ ) {
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0]
                 || y+j < 0 || y+j >= theDim[1]
                 || z+k < 0 || z+k >= theDim[2] )
              fprintf( stdout, " ....." );
          else
            fprintf( stdout, " %5d", theBuf[ ((z+k)*theDim[1]+(y+j))*theDim[0]+(x+i) ] );
          }
          if ( k < 1 ) fprintf( stdout, "  - " );
        }
        fprintf( stdout, "\n" );
      }
      fprintf( stdout, "\n" );
    }
    break;

  case UINT :
    {
      u32 *theBuf = (u32*)buf;
      for ( j=-1; j<=1; j++ ) {
        for ( k=-1; k<=1; k++ ) {
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0]
                 || y+j < 0 || y+j >= theDim[1]
                 || z+k < 0 || z+k >= theDim[2] )
              fprintf( stdout, " ....." );
          else
            fprintf( stdout, " %9d", theBuf[ ((z+k)*theDim[1]+(y+j))*theDim[0]+(x+i) ] );
          }
          if ( k < 1 ) fprintf( stdout, "  - " );
        }
        fprintf( stdout, "\n" );
      }
      fprintf( stdout, "\n" );
    }
    break;

  }
}


static int _get_parcel_center( typeParcel *parcel, typeParcelCenter *center,
                               int label,
                               u16* theDist, u16* theLabels, int *theDim, int inc[3][3][3] )
{
  char *proc = "_get_parcel_center";
  u16 *theLocalLabels = NULL;
  u16 *theLocalDist = NULL;
  double *theLocalCriteria = NULL;
  int theLocalDim[3];

  typeParcelCenter theBestCenter, theCenter;
  int indexBestCenter, indexCenter;

  int x, y, z;
  int lx, ly, lz, i, j, k, v;
  int flag;
  
  int localverbose = 0;


  theLocalDim[0] = parcel[0].xmax - parcel[0].xmin + 1;
  theLocalDim[1] = parcel[0].ymax - parcel[0].ymin + 1;
  theLocalDim[2] = parcel[0].zmax - parcel[0].zmin + 1;

  theLocalLabels = (u16*)malloc( 2*theLocalDim[0]*theLocalDim[1]*theLocalDim[2]*sizeof(u16) );
  if ( theLocalLabels == (u16*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 1st auxiliary array\n", proc );
    return( -1 );
  }
  
  theLocalDist = theLocalLabels + theLocalDim[0]*theLocalDim[1]*theLocalDim[2];

  theLocalCriteria = (double*)malloc( theLocalDim[0]*theLocalDim[1]*theLocalDim[2]*sizeof(double) );
  if ( theLocalCriteria == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 2nd auxiliary array\n", proc );
    free( theLocalLabels );
    return( -1 );
  }

  /* get the parcel and the distances associated to the actual center
   */
  for ( i=0, lz=0; lz<theLocalDim[2]; lz++ ) {
    z = lz + parcel[0].zmin;
    for ( ly=0; ly<theLocalDim[1]; ly++ ) {
      y = ly + parcel[0].ymin;
      j = (z * theDim[1] + y) * theDim[0] + parcel[0].xmin;
      for ( lx=0; lx<theLocalDim[0]; lx++, i++, j++ ) {
        theLocalCriteria[i] = 0.0;
        if ( theLabels[j] == label ) {
          theLocalLabels[i] = label;
          theLocalDist[i] = theDist[j];
        }
        else {
          theLocalLabels[i] = 0;
          theLocalDist[i] = 0;
        }
      }
    }
  }

  /* computed criterium for the actual center
   */
  theBestCenter.z = z = center[0].z - parcel[0].zmin;
  theBestCenter.y = y = center[0].y - parcel[0].ymin;
  theBestCenter.x = x = center[0].x - parcel[0].xmin;

  indexBestCenter = (z*theLocalDim[1] + y)*theLocalDim[0] + x;

  theLocalCriteria[ indexBestCenter ] = 0.0;
  for ( v=0, lz=0; lz<theLocalDim[2]; lz++ )
  for ( ly=0; ly<theLocalDim[1]; ly++ )
  for ( lx=0; lx<theLocalDim[0]; lx++, v++ ) {
    theLocalCriteria[ indexBestCenter ] += theLocalDist[v]*theLocalDist[v];
  }

  if ( localverbose )
    fprintf( stdout, "   1st Criterium (%d %d %d)= %f\n",
             theBestCenter.x, theBestCenter.y, theBestCenter.z,
             theLocalCriteria[ indexBestCenter ] );

  do {

    flag = 0;
    for ( k=-1; !flag && k<=1; k++ ) {

      theCenter.z = theBestCenter.z+k;
      if ( theCenter.z < 0 || theCenter.z >= theLocalDim[2] ) continue;

      for ( j=-1; !flag && j<=1; j++ ) {

        theCenter.y = theBestCenter.y+j;
        if ( theCenter.y < 0 || theCenter.y >= theLocalDim[1] ) continue;

        for ( i=-1; !flag && i<=1; i++ ) {

          theCenter.x = theBestCenter.x+i;
          if ( theCenter.x < 0 || theCenter.x >= theLocalDim[0] ) continue;

          indexCenter = (theCenter.z*theLocalDim[1] + theCenter.y)*theLocalDim[0] + theCenter.x;

          if ( theLocalCriteria[ indexCenter ] > 0.0 ) continue;

          /* not in the parcel
           */
          if ( theLocalLabels[ indexCenter ] == 0 ) continue;

          if ( _get_parcels_from_centers( theLocalDist, theLocalLabels, theLocalDim,
                                          &theCenter, 1, inc ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to compute criteria\n", proc );
            free( theLocalCriteria );
            free( theLocalLabels );
            return( -1 );
          }

          theLocalCriteria[ indexCenter ] = 0.0;
          for ( v=0, lz=0; lz<theLocalDim[2]; lz++ )
            for ( ly=0; ly<theLocalDim[1]; ly++ )
              for ( lx=0; lx<theLocalDim[0]; lx++, v++ ) {
                theLocalCriteria[ indexCenter ] += theLocalDist[v]*theLocalDist[v];
              }

          if ( localverbose )
            fprintf( stdout, "      test  (%d %d %d)= %f",
                     theCenter.x, theCenter.y, theCenter.z,
                     theLocalCriteria[ indexCenter ] );

          if ( theLocalCriteria[ indexCenter ] < theLocalCriteria[ indexBestCenter ] ) {
            flag = 1;
            theBestCenter = theCenter;
            indexBestCenter = indexCenter;
            if ( localverbose )
              fprintf( stdout, " => SUCCESS\n");
          }
          else {
            if ( localverbose )
              fprintf( stdout, " -> failure\n");
          }
          if ( localverbose )
            _print_neighorhood( theCenter.x, theCenter.y, theCenter.z, theLocalCriteria, DOUBLE, theLocalDim, NULL );

        }
      }
    }


  } while ( flag );


  center[0].x = theBestCenter.x + parcel[0].xmin;
  center[0].y = theBestCenter.y + parcel[0].ymin;
  center[0].z = theBestCenter.z + parcel[0].zmin;

  free( theLocalCriteria );
  free( theLocalLabels );
  return( 1 );
}

static int _get_parcel_center_u32( typeParcel *parcel, typeParcelCenter *center,
                               int label,
                               u32* theDist, u32* theLabels, int *theDim, int inc[3][3][3] )
{
  char *proc = "_get_parcel_center";
  u32 *theLocalLabels = NULL;
  u32 *theLocalDist = NULL;
  double *theLocalCriteria = NULL;
  int theLocalDim[3];

  typeParcelCenter theBestCenter, theCenter;
  int indexBestCenter, indexCenter;

  int x, y, z;
  int lx, ly, lz, i, j, k, v;
  int flag;

  int localverbose = 0;


  theLocalDim[0] = parcel[0].xmax - parcel[0].xmin + 1;
  theLocalDim[1] = parcel[0].ymax - parcel[0].ymin + 1;
  theLocalDim[2] = parcel[0].zmax - parcel[0].zmin + 1;

  theLocalLabels = (u32*)malloc( 2*theLocalDim[0]*theLocalDim[1]*theLocalDim[2]*sizeof(u32) );
  if ( theLocalLabels == (u32*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 1st auxiliary array\n", proc );
    return( -1 );
  }

  theLocalDist = theLocalLabels + theLocalDim[0]*theLocalDim[1]*theLocalDim[2];

  theLocalCriteria = (double*)malloc( theLocalDim[0]*theLocalDim[1]*theLocalDim[2]*sizeof(double) );
  if ( theLocalCriteria == (double*)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate 2nd auxiliary array\n", proc );
    free( theLocalLabels );
    return( -1 );
  }

  /* get the parcel and the distances associated to the actual center
   */
  for ( i=0, lz=0; lz<theLocalDim[2]; lz++ ) {
    z = lz + parcel[0].zmin;
    for ( ly=0; ly<theLocalDim[1]; ly++ ) {
      y = ly + parcel[0].ymin;
      j = (z * theDim[1] + y) * theDim[0] + parcel[0].xmin;
      for ( lx=0; lx<theLocalDim[0]; lx++, i++, j++ ) {
        theLocalCriteria[i] = 0.0;
        if ( theLabels[j] == label ) {
          theLocalLabels[i] = label;
          theLocalDist[i] = theDist[j];
        }
        else {
          theLocalLabels[i] = 0;
          theLocalDist[i] = 0;
        }
      }
    }
  }

  /* computed criterium for the actual center
   */
  theBestCenter.z = z = center[0].z - parcel[0].zmin;
  theBestCenter.y = y = center[0].y - parcel[0].ymin;
  theBestCenter.x = x = center[0].x - parcel[0].xmin;
  
  indexBestCenter = (z*theLocalDim[1] + y)*theLocalDim[0] + x;
  
  theLocalCriteria[ indexBestCenter ] = 0.0;
  for ( v=0, lz=0; lz<theLocalDim[2]; lz++ ) 
  for ( ly=0; ly<theLocalDim[1]; ly++ )
  for ( lx=0; lx<theLocalDim[0]; lx++, v++ ) {
    theLocalCriteria[ indexBestCenter ] += theLocalDist[v]*theLocalDist[v];
  }

  if ( localverbose )
    fprintf( stdout, "   1st Criterium (%d %d %d)= %f\n", 
             theBestCenter.x, theBestCenter.y, theBestCenter.z,
             theLocalCriteria[ indexBestCenter ] );
  
  do {

    flag = 0;
    for ( k=-1; !flag && k<=1; k++ ) {

      theCenter.z = theBestCenter.z+k;
      if ( theCenter.z < 0 || theCenter.z >= theLocalDim[2] ) continue;

      for ( j=-1; !flag && j<=1; j++ ) {

        theCenter.y = theBestCenter.y+j;
        if ( theCenter.y < 0 || theCenter.y >= theLocalDim[1] ) continue;

        for ( i=-1; !flag && i<=1; i++ ) {

          theCenter.x = theBestCenter.x+i;
          if ( theCenter.x < 0 || theCenter.x >= theLocalDim[0] ) continue;

          indexCenter = (theCenter.z*theLocalDim[1] + theCenter.y)*theLocalDim[0] + theCenter.x;

          if ( theLocalCriteria[ indexCenter ] > 0.0 ) continue;

          /* not in the parcel
           */
          if ( theLocalLabels[ indexCenter ] == 0 ) continue;

          if ( _get_parcels_from_centers_u32( theLocalDist, theLocalLabels, theLocalDim,
                                          &theCenter, 1, inc ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to compute criteria\n", proc );
            free( theLocalCriteria );
            free( theLocalLabels );
            return( -1 );
          }

          theLocalCriteria[ indexCenter ] = 0.0;
          for ( v=0, lz=0; lz<theLocalDim[2]; lz++ )
            for ( ly=0; ly<theLocalDim[1]; ly++ )
              for ( lx=0; lx<theLocalDim[0]; lx++, v++ ) {
                theLocalCriteria[ indexCenter ] += theLocalDist[v]*theLocalDist[v];
              }

          if ( localverbose )
            fprintf( stdout, "      test  (%d %d %d)= %f",
                     theCenter.x, theCenter.y, theCenter.z,
                     theLocalCriteria[ indexCenter ] );

          if ( theLocalCriteria[ indexCenter ] < theLocalCriteria[ indexBestCenter ] ) {
            flag = 1;
            theBestCenter = theCenter;
            indexBestCenter = indexCenter;
            if ( localverbose )
              fprintf( stdout, " => SUCCESS\n");
          }
          else {
            if ( localverbose )
              fprintf( stdout, " -> failure\n");
          }
          if ( localverbose )
            _print_neighorhood( theCenter.x, theCenter.y, theCenter.z, theLocalCriteria, DOUBLE, theLocalDim, NULL );

        }
      }
    }


  } while ( flag );

  
  center[0].x = theBestCenter.x + parcel[0].xmin;
  center[0].y = theBestCenter.y + parcel[0].ymin;
  center[0].z = theBestCenter.z + parcel[0].zmin;

  free( theLocalCriteria );
  free( theLocalLabels );
  return( 1 );
}


static int _get_centers_from_parcels( u16* theDist, u16* theLabels, int *theDim, 
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] )
{
  char *proc = "_get_centers_from_parcels";
  int i, x, y, z, l, v;
  typeParcel *parcel;

  parcel = (typeParcel *)malloc( nparcels*sizeof( typeParcel ) );
  if ( parcel == (typeParcel *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate parcel list\n", proc );
    return( -1 );
  }

  for ( i=0; i<nparcels; i++ ) {
    parcel[i].xmin = parcel[i].xmax = center[i].x;
    parcel[i].ymin = parcel[i].ymax = center[i].y;
    parcel[i].zmin = parcel[i].zmax = center[i].z;
    parcel[i].n = 0;
    parcel[i].sx = parcel[i].sy = parcel[i].sz = 0.0; 
  }

  for ( i=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ )
  for ( x=0; x<theDim[0]; x++, i++ ) {
    
    if ( theLabels[i] == 0 || theLabels[i] > nparcels ) continue;
    l = theLabels[i] - 1;

    if ( parcel[l].xmin > x ) parcel[l].xmin = x;
    if ( parcel[l].xmax < x ) parcel[l].xmax = x;
    if ( parcel[l].ymin > y ) parcel[l].ymin = y;
    if ( parcel[l].ymax < y ) parcel[l].ymax = y;
    if ( parcel[l].zmin > z ) parcel[l].zmin = z;
    if ( parcel[l].zmax < z ) parcel[l].zmax = z;

    parcel[l].n ++;

    parcel[l].sx += x;
    parcel[l].sy += y;
    parcel[l].sz += z;

  }

  
  for ( i=0; i<nparcels; i++ ) {
    if (_verbose_ >= 1 && (i%(nparcels/20)) == 0 )
      fprintf( stderr, "%s: parcel %d\ton\t%d\t(%d)\n", proc, i+1, nparcels, (100*(i+1))/nparcels);
    parcel[i].x = (int)( parcel[i].sx/(double)parcel[i].n + 0.5 );
    parcel[i].y = (int)( parcel[i].sy/(double)parcel[i].n + 0.5 );
    parcel[i].z = (int)( parcel[i].sz/(double)parcel[i].n + 0.5 );
    v = (parcel[i].z * theDim[1] + parcel[i].y) * theDim[0] + parcel[i].x;
    
    if ( !_force_exact_center_calculation_ && theLabels[v] == i+1 ) {

      center[i].x = parcel[i].x;
      center[i].y = parcel[i].y;
      center[i].z = parcel[i].z;
      
    }
    
    else {
      if ( _verbose_ >= 2 )
        fprintf( stderr, "%s: switch to iterative method for parcel %d\n", proc, i);
      
      if ( _get_parcel_center( &(parcel[i]), &(center[i]), i+1,
                               theDist, theLabels, theDim, inc ) <= 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: iterative calculation of center failed for parcel %d\n", proc, i+1 );
        return( -1 );
      }
    }

  }


  free( parcel );
  return( 1 );
}

static int _get_centers_from_parcels_u32( u32* theDist, u32* theLabels, int *theDim,
                                      typeParcelCenter *center, int nparcels, int inc[3][3][3] )
{
  char *proc = "_get_centers_from_parcels";
  int i, x, y, z, l, v;
  typeParcel *parcel;

  parcel = (typeParcel *)malloc( nparcels*sizeof( typeParcel ) );
  if ( parcel == (typeParcel *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate parcel list\n", proc );
    return( -1 );
  }

  for ( i=0; i<nparcels; i++ ) {
    parcel[i].xmin = parcel[i].xmax = center[i].x;
    parcel[i].ymin = parcel[i].ymax = center[i].y;
    parcel[i].zmin = parcel[i].zmax = center[i].z;
    parcel[i].n = 0;
    parcel[i].sx = parcel[i].sy = parcel[i].sz = 0.0;
  }

  for ( i=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ )
  for ( x=0; x<theDim[0]; x++, i++ ) {

    if ( theLabels[i] == 0 || theLabels[i] > nparcels ) continue;
    l = theLabels[i] - 1;

    if ( parcel[l].xmin > x ) parcel[l].xmin = x;
    if ( parcel[l].xmax < x ) parcel[l].xmax = x;
    if ( parcel[l].ymin > y ) parcel[l].ymin = y;
    if ( parcel[l].ymax < y ) parcel[l].ymax = y;
    if ( parcel[l].zmin > z ) parcel[l].zmin = z;
    if ( parcel[l].zmax < z ) parcel[l].zmax = z;

    parcel[l].n ++;

    parcel[l].sx += x;
    parcel[l].sy += y;
    parcel[l].sz += z;

  }


  for ( i=0; i<nparcels; i++ ) {
    if (_verbose_ >= 2 && (i%(nparcels/20)) == 0 )
      fprintf( stderr, "%s: parcel %d\ton\t%d\t(%d )\n", proc, i+1, nparcels, (100*(i+1))/nparcels);
    parcel[i].x = (int)( parcel[i].sx/(double)parcel[i].n + 0.5 );
    parcel[i].y = (int)( parcel[i].sy/(double)parcel[i].n + 0.5 );
    parcel[i].z = (int)( parcel[i].sz/(double)parcel[i].n + 0.5 );
    v = (parcel[i].z * theDim[1] + parcel[i].y) * theDim[0] + parcel[i].x;

    if ( !_force_exact_center_calculation_ && theLabels[v] == i+1 ) {

      center[i].x = parcel[i].x;
      center[i].y = parcel[i].y;
      center[i].z = parcel[i].z;

    }

    else {
      if ( _verbose_ >= 2 )
        fprintf( stderr, "%s: switch to iterative method for parcel %d\n", proc, i);

      if ( _get_parcel_center_u32( &(parcel[i]), &(center[i]), i+1,
                               theDist, theLabels, theDim, inc ) <= 0 ) {
        if ( _verbose_ )
          fprintf( stderr, "%s: iterative calculation of center failed for parcel %d\n", proc, i+1 );
        return( -1 );
      }
    }
  
  }


  free( parcel );
  return( 1 );
}





static int _get_parcels_from_centers( u16* theDist, u16* theLabels, int *theDim, 
                                      typeParcelCenter *center, int nparcels,
                                      int inc[3][3][3] )
{
  char *proc = "_get_parcels_from_centers";
  int MAXDIST = 65535;
  int NONINLIST, INLIST;
  
  typePointList thePointList;
  int i, j;
  int v = theDim[2]*theDim[1]*theDim[0];
  int nNewPoints, dmin;


  if ( 0 && _verbose_ >= 2 ) {
    fprintf( stderr, "%s:\n", proc );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][0][0], inc[0][0][1], inc[0][0][2],
             inc[1][0][0], inc[1][0][1], inc[1][0][2],
             inc[2][0][0], inc[2][0][1], inc[2][0][2] );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][1][0], inc[0][1][1], inc[0][1][2],
             inc[1][1][0], inc[1][1][1], inc[1][1][2],
             inc[2][1][0], inc[2][1][1], inc[2][1][2] );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][2][0], inc[0][2][1], inc[0][2][2],
             inc[1][2][0], inc[1][2][1], inc[1][2][2],
             inc[2][2][0], inc[2][2][1], inc[2][2][2] );
  }
  
  /* we need some additionnal values
   */
  if ( nparcels >= 65534 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: too many parcels\n", proc );
    return( -1 );
  }
  INLIST = nparcels+1;
  NONINLIST = nparcels+2;

  

  /* initialize labels image and distance map
     background : 0
     points other than seeds : MAXDIST
     seeds : 0
   */
  for ( i=0; i<v; i++ ) {
    if ( theLabels[ i ] ) {
      theLabels[ i ] = NONINLIST;
      theDist[ i ] = MAXDIST;
    }
    else {
      theDist[ i ] = 0;
    }
  }
  for ( i=0; i<nparcels; i++ ) {
    j = (center[i].z * theDim[1] + center[i].y) * theDim[0] + center[i].x;
    theLabels[ j ] = i+1;
    theDist[ j ] = 0;
  }
  
  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: initial centers\n", proc );
    for ( i=0; i<nparcels; i++ ) {
      fprintf( stderr, "center[#%3d] = %d %d %d\n", i, center[i].x, center[i].y, center[i].z );
    }
  }


  /* initialize point list
     get neighbors of parcels
     remove initial points
   */

  _init_point_list( &thePointList );
  for ( i=0; i<nparcels; i++ ) {
    j = (center[i].z * theDim[1] + center[i].y) * theDim[0] + center[i].x;
    if ( _add_point_to_list( &thePointList, center[i].x, center[i].y, center[i].z, 
                             j, i+1, 0, theDim ) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add point to list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }
  }
  nNewPoints = thePointList.nPoints;
  
  if ( _get_parcel_neighbors( &thePointList, nNewPoints, theLabels, theDim, INLIST, NONINLIST ) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add neighbors to list\n", proc );
    _free_point_list( &thePointList );
    return( -1 );
  }

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: %d/%d points at distance 0 in list (%d parcels)\n",
             proc, nNewPoints, thePointList.nPoints, nparcels );
  }

  if ( _remove_first_points_from_list( &thePointList, nNewPoints) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to remove points from list\n", proc );
    _free_point_list( &thePointList );
    return( -1 );
  }
  
  do {
    
    _update_distances_in_list( &thePointList, inc, MAXDIST, theDist, theLabels, theDim, nparcels );

    dmin = _get_min_distance( &thePointList, &nNewPoints );

    if ( 0 && _verbose_ >= 3 )
      fprintf( stderr, "%s: %d/%d points at distance %d in list (%d parcels)\n",
               proc, nNewPoints, thePointList.nPoints, dmin, nparcels );
    
    if ( _get_parcel_neighbors( &thePointList, nNewPoints, theLabels, theDim, INLIST, NONINLIST ) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add neighbors to list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }

    for ( i=0; i<nNewPoints; i++ ) {
      theDist[ thePointList.pt[i].i ] = thePointList.pt[i].d;
      theLabels[ thePointList.pt[i].i ] = thePointList.pt[i].l;
    }

    if ( _remove_first_points_from_list( &thePointList, nNewPoints) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to remove points from list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }

  } while ( thePointList.nPoints > 0 );



  _free_point_list( &thePointList );

  return( 1 );
}

static int _get_parcels_from_centers_u32( u32* theDist, u32* theLabels, int *theDim,
                                      typeParcelCenter *center, int nparcels,
                                      int inc[3][3][3] )
{
  char *proc = "_get_parcels_from_centers";
  int MAXDIST = 2147483647;
  int NONINLIST, INLIST;

  typePointList thePointList;
  int i, j;
  int v = theDim[2]*theDim[1]*theDim[0];
  int nNewPoints, dmin;


  if ( 0 && _verbose_ >= 2 ) {
    fprintf( stderr, "%s:\n", proc );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][0][0], inc[0][0][1], inc[0][0][2],
             inc[1][0][0], inc[1][0][1], inc[1][0][2],
             inc[2][0][0], inc[2][0][1], inc[2][0][2] );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][1][0], inc[0][1][1], inc[0][1][2],
             inc[1][1][0], inc[1][1][1], inc[1][1][2],
             inc[2][1][0], inc[2][1][1], inc[2][1][2] );
    fprintf( stderr, "%d %d %d  -  %d %d %d  -  %d %d %d\n",
             inc[0][2][0], inc[0][2][1], inc[0][2][2],
             inc[1][2][0], inc[1][2][1], inc[1][2][2],
             inc[2][2][0], inc[2][2][1], inc[2][2][2] );
  }

  /* we need some additionnal values
   */
  if ( nparcels >= 2147483646 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: too many parcels\n", proc );
    return( -1 );
  }
  INLIST = nparcels+1;
  NONINLIST = nparcels+2;



  /* initialize labels image and distance map
     background : 0
     points other than seeds : MAXDIST
     seeds : 0
   */
  for ( i=0; i<v; i++ ) {
    if ( theLabels[ i ] ) {
      theLabels[ i ] = NONINLIST;
      theDist[ i ] = MAXDIST;
    }
    else {
      theDist[ i ] = 0;
    }
  }
  for ( i=0; i<nparcels; i++ ) {
    j = (center[i].z * theDim[1] + center[i].y) * theDim[0] + center[i].x;
    theLabels[ j ] = i+1;
    theDist[ j ] = 0;
  }

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: initial centers\n", proc );
    for ( i=0; i<nparcels; i++ ) {
      fprintf( stderr, "center[#%3d] = %d %d %d\n", i, center[i].x, center[i].y, center[i].z );
    }
  }


  /* initialize point list
     get neighbors of parcels
     remove initial points
   */

  _init_point_list( &thePointList );
  for ( i=0; i<nparcels; i++ ) {
    j = (center[i].z * theDim[1] + center[i].y) * theDim[0] + center[i].x;
    if ( _add_point_to_list( &thePointList, center[i].x, center[i].y, center[i].z,
                             j, i+1, 0, theDim ) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add point to list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }
  }
  nNewPoints = thePointList.nPoints;

  if ( _get_parcel_neighbors_u32( &thePointList, nNewPoints, theLabels, theDim, INLIST, NONINLIST ) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to add neighbors to list\n", proc );
    _free_point_list( &thePointList );
    return( -1 );
  }

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: %d/%d points at distance 0 in list (%d parcels)\n",
             proc, nNewPoints, thePointList.nPoints, nparcels );
  }

  if ( _remove_first_points_from_list( &thePointList, nNewPoints) <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to remove points from list\n", proc );
    _free_point_list( &thePointList );
    return( -1 );
  }

  do {

    _update_distances_in_list_u32( &thePointList, inc, MAXDIST, theDist, theLabels, theDim, nparcels );

    dmin = _get_min_distance( &thePointList, &nNewPoints );

    if ( 0 && _verbose_ >= 3 )
      fprintf( stderr, "%s: %d/%d points at distance %d in list (%d parcels)\n",
               proc, nNewPoints, thePointList.nPoints, dmin, nparcels );

    if ( _get_parcel_neighbors_u32( &thePointList, nNewPoints, theLabels, theDim, INLIST, NONINLIST ) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to add neighbors to list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }

    for ( i=0; i<nNewPoints; i++ ) {
      theDist[ thePointList.pt[i].i ] = thePointList.pt[i].d;
      theLabels[ thePointList.pt[i].i ] = thePointList.pt[i].l;
    }
    
    if ( _remove_first_points_from_list( &thePointList, nNewPoints) <= 0 ) {
      if ( _verbose_ )
        fprintf( stderr, "%s: unable to remove points from list\n", proc );
      _free_point_list( &thePointList );
      return( -1 );
    }

  } while ( thePointList.nPoints > 0 );
  

    
  _free_point_list( &thePointList );
    
  return( 1 );
}





static int _remove_first_points_from_list( typePointList *thePointList, int npts )
{
  char *proc = "_remove_first_points_from_list";
  int i;
  int o = thePointList->nPoints;

  if ( npts > thePointList->nPoints ) {
    fprintf( stderr, "%s: too many points to be removed\n", proc );
    return( -1 );
  }
  
  for ( i=npts; i<thePointList->nPoints; i++ )
    thePointList->pt[i-npts] = thePointList->pt[i];
  thePointList->nPoints -= npts;

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: remove %d points from %d -> %d\n", proc, npts, o, thePointList->nPoints );
  }

  return( 1 );
}





static int _get_parcel_neighbors( typePointList *thePointList, int npts,
                                  u16* theLabels, int *theDim,
                                  int INLIST, int NONINLIST )
{
  char *proc = "_get_parcel_neighbors";
  int x, y, z;
  int i, j, k;
  int n, v;
  int o = thePointList->nPoints;
  int a = 0;


  for ( n=0; n<npts; n++ ) {

    x = thePointList->pt[n].x;
    y = thePointList->pt[n].y;
    z = thePointList->pt[n].z;

    if ( thePointList->pt[n].inside ) {
      if ( theDim[2] == 1 ) {
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (y+j)*theDim[0] + x+i;
          if ( theLabels[v] != NONINLIST ) continue;
          if ( _add_point_to_list( thePointList, x+i, y+j, z, v, -1, -1, theDim ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to add point to list\n", proc );
            return( -1 );
          }
          a ++;
          theLabels[v] = INLIST;
        }
      }
      else {
        for ( k=-1; k<=1; k++ )
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
          if ( theLabels[v] != NONINLIST ) continue;
          if ( _add_point_to_list( thePointList, x+i, y+j, z+k, v, -1, -1, theDim ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to add point to list\n", proc );
            return( -1 );
          }
          a ++;
          theLabels[v] = INLIST;
        }
      }
    }
    else {
      for ( k=-1; k<=1; k++ ) {
        if ( z+k < 0 || z+k >= theDim[2] ) continue;
        for ( j=-1; j<=1; j++ ) {
          if ( y+j < 0 || y+j >= theDim[1] ) continue;
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0] ) continue;
            v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
            if ( theLabels[v] != NONINLIST ) continue;
            if ( _add_point_to_list( thePointList, x+i, y+j, z+k, v, -1, -1, theDim ) <= 0 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to add point to list\n", proc );
              return( -1 );
            }
            a ++;
            theLabels[v] = INLIST;
          }
        }
      }
    }

  }

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: add %d points to %d -> %d\n", proc, a, o, thePointList->nPoints );
  }

  return( 1 );
}

static int _get_parcel_neighbors_u32( typePointList *thePointList, int npts,
                                  u32* theLabels, int *theDim,
                                  int INLIST, int NONINLIST )
{
  char *proc = "_get_parcel_neighbors";
  int x, y, z;
  int i, j, k;
  int n, v;
  int o = thePointList->nPoints;
  int a = 0;
  

  for ( n=0; n<npts; n++ ) {

    x = thePointList->pt[n].x;
    y = thePointList->pt[n].y;
    z = thePointList->pt[n].z;

    if ( thePointList->pt[n].inside ) {
      if ( theDim[2] == 1 ) {
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (y+j)*theDim[0] + x+i;
          if ( theLabels[v] != NONINLIST ) continue;
          if ( _add_point_to_list( thePointList, x+i, y+j, z, v, -1, -1, theDim ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to add point to list\n", proc );
            return( -1 );
          }
          a ++;
          theLabels[v] = INLIST;
        }
      }
      else {
        for ( k=-1; k<=1; k++ )
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
          if ( theLabels[v] != NONINLIST ) continue;
          if ( _add_point_to_list( thePointList, x+i, y+j, z+k, v, -1, -1, theDim ) <= 0 ) {
            if ( _verbose_ )
              fprintf( stderr, "%s: unable to add point to list\n", proc );
            return( -1 );
          }
          a ++;
          theLabels[v] = INLIST;
        }
      }
    }
    else {
      for ( k=-1; k<=1; k++ ) {
        if ( z+k < 0 || z+k >= theDim[2] ) continue;
        for ( j=-1; j<=1; j++ ) {
          if ( y+j < 0 || y+j >= theDim[1] ) continue;
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0] ) continue;
            v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
            if ( theLabels[v] != NONINLIST ) continue;
            if ( _add_point_to_list( thePointList, x+i, y+j, z+k, v, -1, -1, theDim ) <= 0 ) {
              if ( _verbose_ )
                fprintf( stderr, "%s: unable to add point to list\n", proc );
              return( -1 );
            }
            a ++;
            theLabels[v] = INLIST;
          }
        }
      }
    }

  }

  if ( 0 && _verbose_ >= 3 ) {
    fprintf( stderr, "%s: add %d points to %d -> %d\n", proc, a, o, thePointList->nPoints );
  }

  return( 1 );
}


static void _update_distances_in_list( typePointList *thePointList, int inc[3][3][3],
                                       int MAXDIST,
                                       u16* theDist, u16* theLabels, int *theDim,
                                       int nparcels )
{
  int n, l, d;
  typePoint *pt = thePointList->pt;
  int v, x, y, z, i, j, k;



  for ( n=0; n<thePointList->nPoints; n++ ) {

    if ( theDist[ pt[n].i ] == 0 ) {
      pt[n].d = 0;
      pt[n].l = theLabels[ pt[n].i ];
      continue;
    }

    x = pt[n].x;
    y = pt[n].y;
    z = pt[n].z;

    d = MAXDIST;
    l = 0;

    if ( pt[n].inside ) {
      if ( theDim[2] == 1 ) {
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (y+j)*theDim[0] + x+i;
          if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
          if ( theDist[v] + inc[1][j+1][i+1] < d ) {
            d = theDist[v] + inc[1][j+1][i+1];
            l = theLabels[v];
          }
        }
      }
      else {
        for ( k=-1; k<=1; k++ )
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
          if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
          if ( theDist[v] + inc[k+1][j+1][i+1] < d ) {
            d = theDist[v] + inc[k+1][j+1][i+1];
            l = theLabels[v];
          }
        }
      }
    }

    else {
      for ( k=-1; k<=1; k++ ) {
        if ( z+k < 0 || z+k >= theDim[2] ) continue;
        for ( j=-1; j<=1; j++ ) {
          if ( y+j < 0 || y+j >= theDim[1] ) continue;
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0] ) continue;
            v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
            if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
            if ( theDist[v] + inc[k+1][j+1][i+1] < d ) {
              d = theDist[v] + inc[k+1][j+1][i+1];
              l = theLabels[v];
            }
          }
        }
      }
    }

    pt[n].d = d;
    pt[n].l = l;

  }
}

static void _update_distances_in_list_u32( typePointList *thePointList, int inc[3][3][3],
                                       int MAXDIST,
                                       u32* theDist, u32* theLabels, int *theDim,
                                       int nparcels )
{
  int n, l, d;
  typePoint *pt = thePointList->pt;
  int v, x, y, z, i, j, k;


  
  for ( n=0; n<thePointList->nPoints; n++ ) {
    
    if ( theDist[ pt[n].i ] == 0 ) {
      pt[n].d = 0;
      pt[n].l = theLabels[ pt[n].i ];
      continue;
    }

    x = pt[n].x;
    y = pt[n].y;
    z = pt[n].z;

    d = MAXDIST;
    l = 0;
    
    if ( pt[n].inside ) {
      if ( theDim[2] == 1 ) {
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (y+j)*theDim[0] + x+i;
          if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
          if ( theDist[v] + inc[1][j+1][i+1] < d ) {
            d = theDist[v] + inc[1][j+1][i+1];
            l = theLabels[v];
          }
        }
      }
      else {
        for ( k=-1; k<=1; k++ )
        for ( j=-1; j<=1; j++ )
        for ( i=-1; i<=1; i++ ) {
          v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
          if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
          if ( theDist[v] + inc[k+1][j+1][i+1] < d ) {
            d = theDist[v] + inc[k+1][j+1][i+1];
            l = theLabels[v];
          }
        }
      }
    }

    else {
      for ( k=-1; k<=1; k++ ) {
        if ( z+k < 0 || z+k >= theDim[2] ) continue;
        for ( j=-1; j<=1; j++ ) {
          if ( y+j < 0 || y+j >= theDim[1] ) continue;
          for ( i=-1; i<=1; i++ ) {
            if ( x+i < 0 || x+i >= theDim[0] ) continue;
            v = (z+k)*theDim[1]*theDim[0] + (y+j)*theDim[0] + x+i;
            if ( theLabels[v] == 0 || theLabels[v] > nparcels ) continue;
            if ( theDist[v] + inc[k+1][j+1][i+1] < d ) {
              d = theDist[v] + inc[k+1][j+1][i+1];
              l = theLabels[v];
            }
          }
        }
      }
    }

    pt[n].d = d;
    pt[n].l = l;
  
  }
}





static int _get_min_distance( typePointList *thePointList, int *nb )
{
  int f, n, min;
  typePoint aux, *pt = thePointList->pt;

  min = pt[0].d;

  for ( n=1; n<thePointList->nPoints; n++ ) {
    if ( pt[n].d < min ) min = pt[n].d;
  }

  f = 0;
  while ( f < thePointList->nPoints && pt[f].d == min ) f++;
  
  for ( n=f; n<thePointList->nPoints; n++ ) {
    if ( pt[n].d == min ) {
      aux = pt[n];
      pt[n] = pt[f];
      pt[f] = aux;
      f ++;
    }
  }
  
  if ( 0 && _verbose_ >= 4 ) {
    fprintf( stderr, "min distance = %d, #pts = %d\n", min, f );
    for ( n=0; n<thePointList->nPoints; n++ ) {
      fprintf( stderr, "D[#%3d] = %9d (%d %d %d)\n", n, pt[n].d, pt[n].x, pt[n].y, pt[n].z );
    }
  }


  *nb = f;
  return( min );
}





static int _get_initial_parcel_centers( u16* theLabels, int *theDim, 
                                        typeParcelCenter *center, int nparcels )
{
  char *proc = "_get_initial_parcel_centers";
  int n, npoints;
  typePoint *pt, tmp;
  int i, j, x, y, z;

  npoints = _count_points( theLabels, USHORT, theDim );

  fprintf(stdout, "npoints = %d\n", npoints);
  fprintf(stdout, "nparcels = %d\n", nparcels);

  if ( npoints <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error or empty input image\n", proc );
    return( -1 );
  }

  if ( npoints <= nparcels ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: too many parcels or not enough points\n", proc );
    return( -1 );
  }
  
  pt = (typePoint *)malloc( npoints*sizeof( typePoint ) );
  if ( pt == (typePoint *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate point list\n", proc );
    return( -1 );
  }

  for ( i=0, j=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ ) 
  for ( x=0; x<theDim[0]; x++, i++ ) {
    if ( theLabels[ i ] ) {
      pt[j].x = x;
      pt[j].y = y;
      pt[j].z = z;
      j ++;
    }
  }

  if ( seedRandom >= 0 ) 
    srand( seedRandom );
  else 
    srand( time(0) );
  
  for ( n = npoints, j=0; j<nparcels; j++, n-- ) {
    i = j+( ((double)rand()/(double)RAND_MAX) * n + 0.5 ); 
    tmp = pt[j];   pt[j] = pt[i];   pt[i] = tmp;
  }
  
  for ( j=0; j<nparcels; j++ ) {
    center[j].x = pt[j].x;
    center[j].y = pt[j].y;
    center[j].z = pt[j].z;
  }

  free( pt );
  return( 1 );
}


static int _get_initial_parcel_centers_u32( u32* theLabels, int *theDim,
                                        typeParcelCenter *center, int nparcels )
{
  char *proc = "_get_initial_parcel_centers";
  int n, npoints;
  typePoint *pt, tmp;;
  int i, j, x, y, z;

  npoints = _count_points( theLabels, UINT, theDim );

  fprintf(stdout, "npoints = %d\n", npoints);
  fprintf(stdout, "nparcels = %d\n", nparcels);

  if ( npoints <= 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: error or empty input image\n", proc );
    return( -1 );
  }

  if ( npoints <= nparcels ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: too many parcels or not enough points\n", proc );
    return( -1 );
  }

  pt = (typePoint *)malloc( npoints*sizeof( typePoint ) );
  if ( pt == (typePoint *)NULL ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate point list\n", proc );
    return( -1 );
  }

  for ( i=0, j=0, z=0; z<theDim[2]; z++ )
  for ( y=0; y<theDim[1]; y++ )
  for ( x=0; x<theDim[0]; x++, i++ ) {
    if ( theLabels[ i ] ) {
      pt[j].x = x;
      pt[j].y = y;
      pt[j].z = z;
      j ++;
    }
  }

  if ( seedRandom >= 0 )
    srand( seedRandom );
  else
    srand( time(0) );

  for ( n = npoints, j=0; j<nparcels; j++, n-- ) {
    i = j+( ((double)rand()/(double)RAND_MAX) * n + 0.5 );
    tmp = pt[j];   pt[j] = pt[i];   pt[i] = tmp;
  }

  for ( j=0; j<nparcels; j++ ) {
    center[j].x = pt[j].x;
    center[j].y = pt[j].y;
    center[j].z = pt[j].z;
  }

  free( pt );
  return( 1 );
}





static int _count_points( void *theBuf, bufferType theBufType, int *theDim )
{
  char *proc = "_count_points";
  int i, v = theDim[0]*theDim[1]*theDim[2];
  int n = 0;

  switch ( theBufType ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such label type not handled yet\n", proc );
    }
    return( -1 );

  case UCHAR :
    {
      u8 *buf = (u8*)theBuf;
      for ( i=0; i<v; i++ )
	if ( buf[i] ) n++;
    }
    break;

  case USHORT :
    {
      u16 *buf = (u16*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  case UINT :
    {
      u32 *buf = (u32*)theBuf;
      for ( i=0; i<v; i++ )
        if ( buf[i] ) n++;
    }
    break;

  }

  return( n );
}





static int _check_free_borders( void *theBuf, bufferType theBufType, int *theDim )
{
  char *proc = "_check_free_borders";
  int x, y, z;

  switch ( theBufType ) {
    
  default :
    if ( _verbose_ ) {
      fprintf( stderr, "%s: such label type not handled yet\n", proc );
    }
    return( -1 );

  case UCHAR :
    {
      u8 *buf = (u8*)theBuf;

      if ( theDim[2] == 1 ) {

	for ( y=0; y<theDim[1]; y++ )
	  if ( buf[ y*theDim[0] ] || buf[ theDim[0]-1 + y*theDim[0] ] )
	    return( -1 );
	for ( x=0; x<theDim[0]; x++ )
	  if ( buf[ x ] || buf[ x + (theDim[1]-1)*theDim[0] ] ) 
	    return( -1 );

      }
      else {

	for ( y=0; y<theDim[1]; y++ ) 
        for ( x=0; x<theDim[0]; x++ )
	  if ( buf[ x + y*theDim[0] ] 
	       || buf[ x + y*theDim[0] + (theDim[2]-1)*theDim[0]*theDim[1] ] )
	    return( -1 );
	for ( z=0; z<theDim[2]; z++ )
        for ( x=0; x<theDim[0]; x++ )
	  if ( buf[ x + z*theDim[0]*theDim[1] ] 
	       || buf[ x + (theDim[1]-1)*theDim[0] + z*theDim[0]*theDim[1] ] )
	    return( -1 );
	for ( z=0; z<theDim[2]; z++ )
	for ( y=0; y<theDim[1]; y++ ) 
	  if ( buf[ y*theDim[0] + z*theDim[0]*theDim[1] ] 
	       || buf[ theDim[0]-1 + y*theDim[0] + z*theDim[0]*theDim[1] ] )
	    return( -1 );

      }

    }
    break;
  }

  return( 1 );
}










/************************************************************
 *
 * List management
 *
 ***********************************************************/

static void _init_point_list( typePointList *l )
{
  l->nPoints = 0;
  l->nAllocatedPoints = 0;
  l->pt = (typePoint *)NULL;
}

static void _free_point_list( typePointList *l )
{
  l->nPoints = 0;
  l->nAllocatedPoints = 0;
  if ( l->pt != (typePoint *)NULL )
    free( l->pt );
  l->pt = (typePoint *)NULL;
}

static int _add_point_to_list( typePointList *l, 
			       int x, int y, int z, int i, 
			       int label, int dist, int *theDim )
{
  char *proc = "_add_point_to_list";
  int n = l->nPoints;
  int newn;
  typePoint *pt = NULL;

  if ( n == l->nAllocatedPoints ) {

    newn = l->nAllocatedPoints + _nPointsToBeAllocated_;

    pt = (typePoint *)malloc( newn * sizeof( typePoint ) );
    if ( pt == NULL ) {
      if ( _verbose_ ) {
	fprintf( stderr, "%s: can not reallocate point list\n", proc );
	fprintf( stderr, "\t failed to add (%d,%d,%d) in list", x, y, z );
      }
      return( -1 );
    }

    if ( l->pt != NULL ) {
      (void)memcpy( (void*)pt, (void*)l->pt,  l->nAllocatedPoints * sizeof( typePoint ) );
      free( l->pt );
    }
    
    l->pt = pt;
    l->nAllocatedPoints = newn;
    
  }

  l->pt[ n ].x = x;
  l->pt[ n ].y = y;
  l->pt[ n ].z = z;
  l->pt[ n ].i = i;
  l->pt[ n ].l = label;
  l->pt[ n ].d = dist;

  l->pt[ n ].inside = 0;
  if ( x > 0 && x < theDim[0]-1 
       && y > 0 && y < theDim[1]-1 )
    l->pt[ n ].inside = 1;
  if ( theDim[2] > 1 && (z <= 0 || z >= theDim[2]-1) )
    l->pt[ n ].inside = 0;
  

  l->nPoints ++;

  return( l->nPoints );
}
