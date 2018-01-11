/*************************************************************************
 * test-inter-pnm.c - test program for image interpolation
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * July, 16 1999
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <iopnm.h>
#include <rl_interpol.h>
#include <typedefs.h>


static char program[256];
static char *usage = "image-1 image-2 generic-name\n\
[-nb %d][-noext] [-i %d] [-2D]"; 
static char *detail = "\n\
\t -nb %d : nombre d images crees\n\
\t -i %d  : indice de l image cree\n";


static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"\nErreur : %s",str);
  exit(0);
}


int main( int argc, char* argv[] )
{
  unsigned char *buffer1 = (unsigned char *)NULL;
  unsigned char *buffer2 = (unsigned char *)NULL;
  int bufferDims1[3] = {0,0,0};
  int bufferDims2[3] = {0,0,0};
  int nbytes1, nbytes2;

  short int *dist1 = (short int *)NULL;
  short int *dist2 = (short int *)NULL;
  int dimz = 256;
  typeDistanceMap theMap1;
  typeDistanceMap theMap2;
  double maxDistanceToBeUpdatedWith5x5x5Chamfer = -1.0;

  int nb=1;
  int indice=-1;
  int writeOriginaux=1;
  int nbNames = 0;
  char nameInput1[256];
  char nameInput2[256];
  char nameGeneric[256];
  


  int i, z;
  int ifirst, ilast;
  unsigned char *resBuf = (unsigned char *)NULL;
  double a, b;

  char outputName[256];
  

  strcpy( program, argv[0] );

  i = 1;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
	_ErrorParse( "message d'aide complet\n", 1);
      }
      
      else if ( strcmp ( argv[i], "-nb" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -nb...\n", 0 );
	if ( sscanf( argv[i],"%d",&nb ) <= 0 ) 
	  _ErrorParse( "parsing -nb...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-i" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -i...\n", 0 );
	if ( sscanf( argv[i],"%d",&indice ) <= 0 ) 
	  _ErrorParse( "parsing -i...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-noext" ) == 0 ) {
	writeOriginaux = 0;
      }
      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	_DistanceSetComputationTo2D();
      }
      else {
	sprintf( nameInput1, "unknown option %s\n", argv[i] );
	_ErrorParse( nameInput1, 0);
      }
    } 
    else {
      if (nbNames == 0) strcpy( nameInput1, argv[i] );
      else if (nbNames == 1) strcpy( nameInput2, argv[i] );
      else if (nbNames == 2) strcpy( nameGeneric, argv[i] );
      else _ErrorParse( "too much file names", 0);
      nbNames++;
    }
    i++;
  }
  
  if ( argc == 0 ) {
    _ErrorParse( "\n", 0 );
  }
  if ( nbNames != 3 ) {
    _ErrorParse( "bad number of file names\n", 0);
  }
  if ( nb <= 0 ) {
    fprintf( stderr, "Error: number of images to be created has to be positive (was equal to %d)\n", nb );
    exit(0);
  }
 

  /* lecture des images d'entree
   */
  buffer1 = _readPnmImage( nameInput1, &bufferDims1[0], &bufferDims1[1], &bufferDims1[2], &nbytes1 );
  buffer2 = _readPnmImage( nameInput2, &bufferDims2[0], &bufferDims2[1], &bufferDims2[2], &nbytes2 );

  if ( buffer1 == (unsigned char *)NULL || buffer2 == (unsigned char *)NULL ) {
    fprintf( stderr, "Error: while reading/allocating %s and %s\n", nameInput1, nameInput2 );
    exit( 0 );
  }

  if ( bufferDims1[0] != bufferDims2[0] || 
       bufferDims1[1] != bufferDims2[1] ||
       bufferDims1[2] != bufferDims2[2] ) {
    fprintf( stderr, "Error: %s and %s should have same dimensions\n", nameInput1, nameInput2 );
    free( buffer1 );
    free( buffer2 );
    exit( 0 );
  }

  if ( nbytes1 != 1 || nbytes2 != 1 ) {
    fprintf( stderr, "Error: unable to deal with such buffer type\n" );
    free( buffer1 );
    free( buffer2 );
    exit( 0 );
  }

  ifirst = 1;
  ilast = nb;
  
  if ( indice >= 1 && indice <= nb ) {
    ifirst = ilast = indice;
  }





  dist1 = (short int *)malloc( bufferDims1[0]*bufferDims1[1]*dimz * sizeof(short int) );
  dist2 = (short int *)malloc( bufferDims1[0]*bufferDims1[1]*dimz * sizeof(short int) ); 
  if ( dist1 == (short int *)NULL || dist2 == (short int *)NULL ) {
    fprintf( stderr, "Error: while allocating distances\n" );
    free( buffer1 );
    free( buffer2 );
    exit( 0 );
  }

  theMap1.buf = dist1;
  theMap1.dim[0] = bufferDims1[0];
  theMap1.dim[1] = bufferDims1[1];
  theMap1.dim[2] = dimz;
  theMap1.voxelSize[0] = 1.0;
  theMap1.voxelSize[1] = 1.0;
  theMap1.voxelSize[2] = 1.0;
  theMap1.multiplicativeCoefficient = 0.0;
  theMap2 = theMap1;
  theMap2.buf = dist2;
    

  resBuf = (unsigned char *)malloc( (ilast-ifirst+1)*bufferDims1[0]*bufferDims2[1]*bufferDims2[2] * sizeof(unsigned char) );
  if ( resBuf == (unsigned char *)NULL ) {
    fprintf( stderr, "Error: while allocating result\n" );
    free( buffer1 );
    free( buffer2 );
    free( dist1 );
    free( dist2 );
  }
  
  _DistanceSetNoVerbose();

  for ( z=0; z<bufferDims1[2]; z++ ) {
    fprintf( stderr, " processing slice #%d\n", z );

    theMap1.intensityMax =_InitDistanceImageFromSlice( (void*)&(buffer1[z*bufferDims1[0]*bufferDims1[1]]), 
						       UCHAR, bufferDims1[0], bufferDims1[1],
						       dist1, dimz );
    _ComputeSignedDistanceMap( &theMap1, maxDistanceToBeUpdatedWith5x5x5Chamfer );
    
    theMap2.intensityMax =_InitDistanceImageFromSlice( (void*)&(buffer2[z*bufferDims1[0]*bufferDims1[1]]), 
				 UCHAR, bufferDims1[0], bufferDims1[1],
				 dist2, dimz );
    _ComputeSignedDistanceMap( &theMap2, maxDistanceToBeUpdatedWith5x5x5Chamfer );

    _CombineTwoDistanceMaps2D( &theMap1, &theMap2 );


    for ( i=ifirst; i<=ilast; i++ ) {
      a = (double)(nb+1-i)/(double)(nb+1);
      b = (double)(i)/(double)(nb+1);
      fprintf( stderr, " creating image %3d of  slice #%d\n", i, z );
      if ( indice >= 1 && indice <= nb ) {
	_ComputeSliceFromDistances( (void*)&resBuf[bufferDims1[0]*bufferDims1[1]*z],
				    UCHAR, bufferDims1[0], bufferDims1[1],
				    dist1, dist2, dimz, a, b );
      } else {
	_ComputeSliceFromDistances( (void*)&(resBuf[ bufferDims1[0]*bufferDims1[1]*( (i-1)*bufferDims1[2] + z ) ]),
				    UCHAR, bufferDims1[0], bufferDims1[1],
				    dist1, dist2, dimz, a, b );
      }
    }
  }

  free( dist1 );
  free( dist2 );

  IoPnm_WriteGreyAsColor();
  IoPnm_SetMaxGreyValueTo255();

  if ( writeOriginaux != 0 ) {
    i = 0;
    sprintf( outputName, "%s.%d.ppm", nameGeneric, i );
    _writePnmImage( outputName, bufferDims1[0], bufferDims1[1], bufferDims1[2], 1, buffer1 );
  }

  for ( i=ifirst; i<=ilast; i++ ) {
    sprintf( outputName, "%s.%d.ppm", nameGeneric, i );
    if ( indice >= 1 && indice <= nb ) {
      _writePnmImage( outputName, bufferDims1[0], bufferDims1[1], bufferDims1[2], 1, resBuf );
    } else {
      _writePnmImage( outputName, bufferDims1[0], bufferDims1[1], bufferDims1[2], 1,
		      &(resBuf[ (i-1)*bufferDims1[0]*bufferDims1[1]*bufferDims1[2] ]) );
    }
  }
  
  if ( writeOriginaux != 0 ) {
    i = nb+1;
    sprintf( outputName, "%s.%d.ppm", nameGeneric, i );
    _writePnmImage( outputName, bufferDims1[0], bufferDims1[1], bufferDims1[2], 1, buffer2 );
  }



  free( buffer1 );
  free( buffer2 );
  free( resBuf );

  exit( 0 );
}
