#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <YelTr3D.h>




static int _verbose_ = 1;

int main( int argc, char* argv[] )
{
  char nameTriangulation1[256];
  char nameTriangulation2[256];
  int nbNames = 0;
  int i;
  int closest = 0;

  typeTriangulation triangulation1;
  typeTriangulation triangulation2;
  int nv;
  double *distance = (double*)NULL;

  for ( i=1; i<argc; i++ ) {
    if ( argv[i][0] == '-' ) {
      
      if ( (strcmp ( argv[i], "-help" ) == 0) || 
	   (strcmp ( argv[i], "-h" ) == 0) ) {
	/*
	ErrorMessage( "help message\n", 1 );
	*/
      }
      
      else if ( (strcmp ( argv[i], "-verbose" ) == 0) || 
		(strcmp ( argv[i], "-v" ) == 0) ) {
	_verbose_ = 1;
      }

      else if ( (strcmp ( argv[i], "-no-verbose" ) == 0) || 
		(strcmp ( argv[i], "-nv" ) == 0) ) {
	_verbose_ = 0;
      }
      
      else if ( (strcmp ( argv[i], "-closest" ) == 0) || 
		(strcmp ( argv[i], "-c" ) == 0) ) {
	closest = 1;
      }
      
    }
    
    else if ( argv[i][0] != 0 ) {
      if ( nbNames == 0 ) {
	strcpy( nameTriangulation1, argv[i] );
      } 
      else if ( nbNames == 1 ) {
	strcpy( nameTriangulation2, argv[i] );
      } 
      else {
	fprintf( stderr, "too many names (%s)\n", argv[i] );
	/*
	ErrorMessage( nameImageIn, 0);
	*/
      }
      nbNames ++;
    }
  }

  initTriangulation( &triangulation1 );
  
  if ( readTriangulation( &triangulation1, nameTriangulation1 ) != 1 ) {
    fprintf( stderr, "error while reading %s\n", nameTriangulation1 );
    return( 0 );
  }

  initTriangulation( &triangulation2 );
  
  if ( readTriangulation( &triangulation2, nameTriangulation2 ) != 1 ) {
    freeTriangulation( &triangulation1 );
    fprintf( stderr, "error while reading %s\n", nameTriangulation2 );
    return( 0 );
  }

  nv = triangulation1.nVertex;
  if ( nv < triangulation2.nVertex ) nv = triangulation2.nVertex;

  distance = (double*)malloc( nv * sizeof(double) );
  if ( distance == (double*)NULL ) {
    freeTriangulation( &triangulation2 );
    freeTriangulation( &triangulation1 );
    fprintf( stderr, "error while allocating distances\n" );
    return( 0 );
  }

  if ( closest ) {
    
    if ( computeDistancesWithClosestVertex( &triangulation1, &triangulation2, distance ) != 1 ) {
      free( distance );
      freeTriangulation( &triangulation2 );
      freeTriangulation( &triangulation1 );    
      fprintf( stderr, "error while computing distances 1 -> 2 \n" );
      return( 0 );
    }
    
    printReport( distance, triangulation1.nVertex, 
		 nameTriangulation1, nameTriangulation2, 
		 "distance towards to the closest vertex" );

    if ( computeDistancesWithClosestVertex( &triangulation2, &triangulation1, distance ) != 1 ) {
      free( distance );
      freeTriangulation( &triangulation2 );
      freeTriangulation( &triangulation1 );    
      fprintf( stderr, "error while computing distances 2 -> 1 \n" );
      return( 0 );
    }
    
    printReport( distance, triangulation2.nVertex, 
		 nameTriangulation2, nameTriangulation1, 
		 "distance towards to the closest vertex" );


  }
  else {

    if ( computeDistancesWithIdenticalVertex( &triangulation1, &triangulation2, 
					      distance ) != 1 ) {
      free( distance );
      freeTriangulation( &triangulation2 );
      freeTriangulation( &triangulation1 );    
      fprintf( stderr, "error while computing distances 1 <-> 2 \n" );
      return( 0 );
    }

    printReport( distance, triangulation1.nVertex, 
		 nameTriangulation1, nameTriangulation2, 
		 "distance towards to the vertex of same index" );



  }

  freeTriangulation( &triangulation2 );
  freeTriangulation( &triangulation1 );






  free( distance );

  return(1);


}
