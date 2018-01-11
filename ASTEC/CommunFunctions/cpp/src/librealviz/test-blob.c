#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <time.h>

#include <iopnm.h>
#include <connexe.h>
#include <convert.h>
#include <_predef_morpho.h>

#define PRINTTIMES 1


static char program[256];
static char *usage = "image-in image-out [-ms %d] [-mn %d] [-v] [-help]";
static char *details ="\n\
Connected components extraction\n\
-------------------------------\n\
The input image is suuposed to be binary:\n\
- object = {points with value > 0}\n\
- background = {points with value = 0}\n\
\n\
\t -ms | -minimal-size : specifies the minimal size of connected components\n\
\t -mn | -maximal-number : specifies the maximal number of ccs to be retained\n\
\t     if the specified number is positive, then the labels of the ccs are\n\
\t     sorted with respect to the cc's size (decreasing order). Thus the component\n\
\t     with label 1 is the largest one.\n\
\t -v | -verbose\n\
\t -nv | -no-verbose\n\
\t -h | -help : print this message";

static void ErrorMessage( char *str, int flag )
{
  (void)fprintf(stderr,"Usage: %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",details);
  (void)fprintf(stderr,"Error: %s",str);
  exit(0);
}







int _ComputeAbsDiff( void *buffer1, 
		     void *buffer2, 
		     void *bufferOut, bufferType type,
		     int *bufferDims );
int _ComputeAreasToBeCorrected( void *buffer1, bufferType type1,
				void *buffer2, bufferType type2,
				void *bufferOut, bufferType typeOut,
				int *bufferDims,
				int lowThreshold,
				int highThreshold,
				int iter_dilation_1,
				int iter_erosion_1,
				int iter_dilation_2 );











int main( int argc, char* argv[] )
{
  int i;
  int nbNames = 0;
  int status;
  char nameImageIn1[256];
  char nameImageIn2[256];
  char nameImageOut[256];

  int lowThreshold  = 100;
  int highThreshold = 100;
  int iter_dilation_1 = 2;
  int iter_erosion_1  = 4;
  int iter_dilation_2 = 8;

  void *bufferIn1 = (void*)NULL;
  void *bufferIn2 = (void*)NULL;
  void *bufferOut = (void*)NULL;
  int bufferIn1Dims[3] = {0,0,0};
  int bufferIn2Dims[3] = {0,0,0};

  int r;

  strcpy( program, argv[0] );

  if ( argc == 1 ) ErrorMessage( "\n", 1 );

  for ( i=1; i<argc; i++ ) {
    if ( argv[i][0] == '-' ) {

      if ( (strcmp ( argv[i], "-help" ) == 0) || 
	   (strcmp ( argv[i], "-h" ) == 0) ) {
	ErrorMessage( "help message\n", 1 );
      }
      
      else if ( (strcmp ( argv[i], "-verbose" ) == 0) || 
		(strcmp ( argv[i], "-v" ) == 0) ) {
	incrementVerboseInConnexe();
      }

      else if ( (strcmp ( argv[i], "-no-verbose" ) == 0) || 
		(strcmp ( argv[i], "-nv" ) == 0) ) {
	Connexe_noverbose();
      }
      
      /*
       * thresholds 
       */
      else if ( (strcmp ( argv[i], "-threshold" ) == 0) || 
		(strcmp ( argv[i], "-t" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -threshold...\n", 0 );
	status = sscanf( argv[i], "%d", &lowThreshold );
	if ( status <= 0 ) ErrorMessage( "parsing -threshold...\n", 0 );
	highThreshold = lowThreshold;
      }
      
      else if ( (strcmp ( argv[i], "-low-threshold" ) == 0) || 
		(strcmp ( argv[i], "-lt" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -low-threshold...\n", 0 );
	status = sscanf( argv[i], "%d", &lowThreshold );
	if ( status <= 0 ) ErrorMessage( "parsing -low-threshold...\n", 0 );
      }
      
      else if ( (strcmp ( argv[i], "-high-threshold" ) == 0) || 
		(strcmp ( argv[i], "-ht" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -high-threshold...\n", 0 );
	status = sscanf( argv[i],"%d", &highThreshold );
	if ( status <= 0 ) ErrorMessage( "parsing -high-threshold...\n", 0 );
      }

      
      /*
       * dilation #1
       */
      else if ( (strcmp ( argv[i], "-dilation-1" ) == 0) || 
		(strcmp ( argv[i], "-d1" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -dilation-1...\n", 0 );
	status = sscanf( argv[i], "%d", &iter_dilation_1 );
	if ( status <= 0 ) ErrorMessage( "parsing -dilation-1...\n", 0 );
      }


      
      /*
       * erosion #1
       */
      else if ( (strcmp ( argv[i], "-erosion-1" ) == 0) || 
		(strcmp ( argv[i], "-e1" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -erosion-1...\n", 0 );
	status = sscanf( argv[i], "%d", &iter_erosion_1 );
	if ( status <= 0 ) ErrorMessage( "parsing -erosion-1...\n", 0 );
      }


      
      /*
       * dilation #2
       */
      else if ( (strcmp ( argv[i], "-dilation-2" ) == 0) || 
		(strcmp ( argv[i], "-d2" ) == 0) ) {
	i += 1;
	if ( i >= argc)    ErrorMessage( "parsing -dilation-2...\n", 0 );
	status = sscanf( argv[i], "%d", &iter_dilation_2 );
	if ( status <= 0 ) ErrorMessage( "parsing -dilation-2...\n", 0 );
      }

      else {
	sprintf( nameImageIn1, "unknown option %s\n", argv[i] );
	ErrorMessage( nameImageIn1, 0);
      }
    }

    else if ( argv[i][0] != 0 ) {
      if ( nbNames == 0 ) {
	strcpy( nameImageIn1, argv[i] );
      } 
      else if ( nbNames == 1 ) {
	strcpy( nameImageIn2, argv[i] );
      } 
      else if ( nbNames == 2 ) {
	strcpy( nameImageOut, argv[i] );
      } 
      else {
	sprintf( nameImageIn1, "too many image name (%s)\n", argv[i] );
	ErrorMessage( nameImageIn1, 0);
      }
      nbNames ++;
    }
  }

  
  bufferIn1 = _readPnmImage( nameImageIn1, 
			    &bufferIn1Dims[0], &bufferIn1Dims[1], &bufferIn1Dims[2] );
  bufferIn2 = _readPnmImage( nameImageIn2, 
			    &bufferIn2Dims[0], &bufferIn2Dims[1], &bufferIn2Dims[2] );
  
  bufferOut = (void*)malloc( bufferIn1Dims[0] * bufferIn1Dims[1] * bufferIn1Dims[2] * sizeof(unsigned char) );


  r = _ComputeAreasToBeCorrected( bufferIn1, UCHAR,
				   bufferIn2, UCHAR,
				   bufferOut, UCHAR,
				   bufferIn1Dims,
				  lowThreshold, highThreshold,
				  iter_dilation_1,
				  iter_erosion_1,
				  iter_dilation_2 );
  switch ( r ) {
  case -1 :
    free( bufferOut );
    ErrorMessage( "error in computation", 0 );
  case 0 :
    free( bufferOut );
    fprintf( stderr, "no area to be corrected\n" );
    exit( 0 );
  }

  fprintf( stderr, "areas to be corrected = %d\n", r );

  _writePnmImage( nameImageOut, 
		  bufferIn1Dims[0], bufferIn1Dims[1], bufferIn1Dims[2], 
		  bufferOut );

  free( bufferOut );

}















int _ComputeAbsDiff( void *buffer1, 
		     void *buffer2, 
		     void *bufferOut, bufferType type,
		     int *bufferDims )
{
  int i;
  int v = bufferDims[0] * bufferDims[1] * bufferDims[2];

  
  
  switch ( type ) {
  default :
    return -1;
  case UCHAR :
    {
      unsigned char *theOut = (unsigned char *)bufferOut;
      unsigned char *the1   = (unsigned char *)buffer1;
      unsigned char *the2   = (unsigned char *)buffer2;
      for ( i=0; i<v; i++ ) {
	theOut[i] = (unsigned char)( ( the1[i] > the2[i] ) ? the1[i] - the2[i] : the2[i] - the1[i] );
      }
    }
    break;
  case USHORT :
    {
      unsigned short int *theOut = (unsigned short int *)bufferOut;
      unsigned short int *the1   = (unsigned short int *)buffer1;
      unsigned short int *the2   = (unsigned short int *)buffer2;
      for ( i=0; i<v; i++ ) {
	theOut[i] = (unsigned short int)( ( the1[i] > the2[i] ) ? the1[i] - the2[i] : the2[i] - the1[i] );
      }
    }
    break;
  }
  return 1;
}














int _ComputeAreasToBeCorrected( void *buffer1, bufferType type1,
				void *buffer2, bufferType type2,
				void *bufferOut, bufferType typeOut,
				int *bufferDims,
				int lowThreshold,
				int highThreshold,
				int iter_dilation_1,
				int iter_erosion_1,
				int iter_dilation_2 )
{

  void *bufferDiff = NULL;
  int n;

#ifdef PRINTTIMES
  int c0 = clock();
  int c1, c2;
#endif


  
  if ( type2 != type1 ) return( -1 );
  if ( typeOut != type1 ) {
    switch( type1 ) {
    default :
      return -1;
    case UCHAR :
      bufferDiff = malloc( bufferDims[0] * bufferDims[1] * bufferDims[2] 
			   * sizeof( unsigned char ) );
      break;
    case USHORT :
      bufferDiff = malloc( bufferDims[0] * bufferDims[1] * bufferDims[2] 
			   * sizeof( unsigned short int ) );
      break;
    }
    if ( bufferDiff == NULL ) {
      return -1;
    }
  } 
  else {
    bufferDiff = bufferOut;
  }

  
  /*
   * difference
   */
#ifdef PRINTTIMES
  c1 = clock();
#endif
  if ( _ComputeAbsDiff( buffer1, buffer2, bufferDiff, type1, bufferDims ) != 1 ) {
    if ( typeOut != type1 ) free( bufferDiff );
    return -1;
  }
#ifdef PRINTTIMES
  c2 = clock();
  fprintf( stderr, " Difference:              %f sec.\n", 
	   (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif


  /*
   * seuillage
   */
  
#ifdef PRINTTIMES
  c1 = clock();
#endif

  if ( lowThreshold <= 0 ||
       lowThreshold >= highThreshold ) {
    int i;
    int v = bufferDims[0] * bufferDims[1] * bufferDims[2];

    n = 0;

    switch ( type1 ) {
    default :
      if ( typeOut != type1 ) free( bufferDiff );
      return -1;
    case UCHAR :
      {
	unsigned char * theBuf = (unsigned char *)bufferDiff;
	for ( i=0; i<v; i++ ) {
	  if ( theBuf[i] >= highThreshold ) {
	    theBuf[i] = 255;
	    n ++;
	  }
	  else {
	    theBuf[i] = 0;
	  }
	}
      }
      break;
    case USHORT :
      {
	unsigned short int * theBuf = (unsigned short int  *)bufferDiff;
	for ( i=0; i<v; i++ ) {
	  if ( theBuf[i] >= highThreshold ) {
	    theBuf[i] = 65535;
	    n ++;
	  }
	  else {
	    theBuf[i] = 0;
	  }
	}
      }
      break;
    }
#ifdef PRINTTIMES
    c2 = clock();
    fprintf( stderr, " Thresholding:            %f sec.\n", 
	     (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif
  }
  else {
    
    n = HysteresisThresholdingWithAllParams( bufferDiff, type1,
					     bufferDiff, type1,
					     bufferDims,
					     (double)lowThreshold, 
					     (double)highThreshold,
					     4, 1, 1, 0, 1  );
    
    if ( n < 0 ) {
      if ( typeOut != type1 ) free( bufferDiff );
      return -1;
    }

#ifdef PRINTTIMES
    c2 = clock();
    fprintf( stderr, " Hysteresis Thresholding: %f sec.\n", 
	     (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif
  }

  if ( n == 0 ) {
    if ( typeOut != type1 ) free( bufferDiff );
    return 0;
  }




  /*
   * Dilatation #1
   */
#ifdef PRINTTIMES
  c1 = clock();
#endif

  switch( type1 ) {
  default :
    if ( typeOut != type1 ) free( bufferDiff );
    return -1;
  case UCHAR :
    _u8_predef_binary_Dilation( (u8*)bufferDiff, (u8*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_dilation_1 );
    break;
  case USHORT :
    _u16_predef_binary_Dilation( (u16*)bufferDiff, (u16*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_dilation_1 );
    break;
  }
#ifdef PRINTTIMES
    c2 = clock();
    fprintf( stderr, " Dilation #1:             %f sec.\n", 
	     (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif





  /*
   * Erosion #1
   */
#ifdef PRINTTIMES
  c1 = clock();
#endif

  switch( type1 ) {
  default :
    if ( typeOut != type1 ) free( bufferDiff );
    return -1;
  case UCHAR :
    _u8_predef_binary_Erosion( (u8*)bufferDiff, (u8*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_erosion_1 );
    break;
  case USHORT :
    _u16_predef_binary_Erosion( (u16*)bufferDiff, (u16*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_erosion_1 );
    break;
  }
#ifdef PRINTTIMES
    c2 = clock();
    fprintf( stderr, " Erosion #1:              %f sec.\n", 
	     (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif





  /*
   * Dilatation #2
   */
#ifdef PRINTTIMES
  c1 = clock();
#endif

  switch( type1 ) {
  default :
    if ( typeOut != type1 ) free( bufferDiff );
    return -1;
  case UCHAR :
    _u8_predef_binary_Dilation( (u8*)bufferDiff, (u8*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_dilation_2 );
    break;
  case USHORT :
    _u16_predef_binary_Dilation( (u16*)bufferDiff, (u16*)bufferDiff,
				bufferDims[0], bufferDims[1], bufferDims[2],
				8, iter_dilation_2 );
    break;
  }
#ifdef PRINTTIMES
  c2 = clock();
  fprintf( stderr, " Dilation #2:             %f sec.\n", 
	   (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif



  
  /*
   * connected components
   */
#ifdef PRINTTIMES
  c1 = clock();
#endif
  
  n = CountConnectedComponentsWithAllParams( bufferDiff, type1,
					     bufferDiff, type1,
					     bufferDims,
					     (double)1.0,
					     8, 1, 0, 0 );
#ifdef PRINTTIMES
  c2 = clock();
  fprintf( stderr, " Connected components:    %f sec.\n", 
	   (double)(c2-c1) / (double)CLOCKS_PER_SEC );
#endif

  if ( n < 0 ) {
    if ( typeOut != type1 ) free( bufferDiff );
    return -1;
  }
    
  if ( n == 0 ) {
    if ( typeOut != type1 ) free( bufferDiff );
    return 0;
  }





  if ( typeOut != type1 ) {
    ConvertBuffer( bufferDiff, type1, bufferOut, typeOut,
		   bufferDims[0] * bufferDims[1] * bufferDims[2] );
    free( bufferDiff );
  }


#ifdef PRINTTIMES
    c2 = clock();
    fprintf( stderr, " Total time:              %f sec.\n", 
	     (double)(c2-c0) / (double)CLOCKS_PER_SEC );
#endif
  return n;
}
