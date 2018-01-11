/*************************************************************************
 * blBuildTracks.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Ven 21 mar 2014 17:46:33 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>



static int _time_ = 1;
static int _clock_ = 1;

static int _verbose_ = 1;
static int _debug_ = 1;



#include <string-tools.h>

#include <bal-image.h>
#include <bal-biolib-tools.h>

#include <bal-biolib-tracker.h>





typedef struct local_parameter {

  /* file names
   */
  char *thedetection_format;
  char *thetransformation_format;
  char *thevesselness_format;

  char *restrack_name;

  int firstindex;
  int lastindex;

  enumUnitTransfo detection_unit;

  bal_blTrackerParameter tracking_parameters;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();




static char *program = NULL;

static char *usage = "[-detection-format | -detection %s]\n\
 [-detection-unit|-du voxel|real]\n\
 [-vesselness-format | -vesselness | -vessel %s]\n\
 [-transformation-format | -transformation | -trsf-format | -trsf %s]\n\
  -f[irst] %d -l[ast] %d\n\
 \n\
 [-half-neighborhood-unit | -hnu %lf %lf %lf]\n\
 [-half-neighborhood-voxel | -hnv %lf %lf %lf]\n\
 [-min-2D-distance-voxel %f]\n\
 [-min-2D-distance-unit %f]\n\
 [-cost-location future|past|both | -use-retraction]\n\
 [-cost-computation exact|approximate]\n\
 [-search-type force|buckets]\n\
 \n\
 [-min-vesselness %lf] [-max-vesselness %lf] [-min-track-length %d] [-max-tracks %d]\n\
 [-time|-notime] [-clock|-noclock]\n\
 [-v|-nv] [-help]";

static char *detail = "\
[-detection-format %s]  # format 'a la printf' of result detection files\n\
                  # must contain one '%d'\n\
[-vesselness-format %s] # format 'a la printf' of vesselness images\n\
                  # must contain one '%d'\n\
[-transformation-format %s] # format 'a la printf' of vesselness images\n\
                  # must contain one '%d'\n\
[-first %d]       # first value of the index in the format\n\
[-last %d]        # last value of the index in the format\n\
\n\
[-half-neighborhood-voxel | -hnv %lf %lf %lf] # half-size of the \n\
                  # neighborhood where correspondances are searched\n\
[-half-neighborhood-unit | -hnu %lf %lf %lf]  # same, in image units\n\
[-min-2D-distance-voxel %f] # below this distance, a maximal connection\n\
                  # is established, else the cost is computed\n\
[-min-2D-distance-unit %f]  #  same, in image units\n\
[-cost-location future|past|both | -use-retraction] # where to compute the cost\n\
                  # future = in next image (assumption is growth only)\n\
                  # both = in current and next images (allows retraction)\n\
[-cost-computation exact|approximate] # cost computation parameter\n\
[-search-type force|buckets]          # cost computation parameter\n\\n\
\n\
[-min-vesselness %lf]  # a connection below this threshold is discarded\n\
[-max-vesselness %lf]  # a connection above this threshold is discarded\n\
[-min-track-length %d] # paths shorter than this length are discarded\n\
[-max-tracks %d]       # maximal number of tracks\n\
[-v]              # increase verboseness\n\
\n";











int main( int argc, char *argv[] )
{
  local_parameter p;

  stringList theDetectionFileList;
  stringList theTransformationFileList;
  stringList theVesselnessFileList;
  int i;

  bal_blTrackerElementList theTrackerList;
  bal_blTrackList resList;
  
  double time_init, time_exit;
  double clock_init, clock_exit;

  time_init = _GetTime();
  clock_init = _GetClock();





  BAL_InitBlTrackList( &resList );


  /***************************************************
   *
   * parsing parameters
   *
   ***************************************************/
  program = argv[0];

  /* no arguments
   */
  if ( argc == 1 ) _ErrorParse( NULL, 0 );


  /* parsing parameters 
   */
  _InitParam( &p );
  _Parse( argc, argv, &p );
  

  

 

  /****************************************
   *
   * format case: building file name
   *
   ****************************************/

  initStringList( &theDetectionFileList );

  if ( buildStringListFromFormat( p.thedetection_format, 
				  p.firstindex, p.lastindex, 
				  &theDetectionFileList ) != 1 ) {
    _ErrorParse( "unable to build input detection file list from format\n", 0);
  }
  if ( 0 && _debug_ ) {
    printStringList( stderr,  &theDetectionFileList, (char*)NULL );
  }



  initStringList( &theTransformationFileList );

  if ( p.thetransformation_format != (char*)NULL ) {
    if ( buildStringListFromFormat( p.thetransformation_format, 
				    p.firstindex, p.lastindex, 
				    &theTransformationFileList ) != 1 ) {
      freeStringList( &theDetectionFileList );
      _ErrorParse( "unable to build input transformation file list from format\n", 0);
    }
    if ( 0 && _debug_ ) {
      printStringList( stderr,  &theTransformationFileList, (char*)NULL );
    }
  }



  initStringList( &theVesselnessFileList );
  
  if ( p.thevesselness_format != (char*)NULL ) {
    if ( buildStringListFromFormat( p.thevesselness_format, 
				    p.firstindex, p.lastindex, 
				    &theVesselnessFileList ) != 1 ) {
      if ( p.thetransformation_format != (char*)NULL ) freeStringList( &theTransformationFileList );
      freeStringList( &theDetectionFileList );
      _ErrorParse( "unable to build input vesselness file list from format\n", 0);
    }
    if ( 0 && _debug_ ) {
      printStringList( stderr,  &theVesselnessFileList, (char*)NULL );
    }
  }



  /****************************************
   *
   * filling the structure
   *
   ****************************************/

  BAL_InitBlTrackerElementList( &theTrackerList );

  if ( BAL_AllocBlTrackerElementList( &theTrackerList, 
				      p.lastindex - p.firstindex + 1 ) != 1 ) {
    if ( p.thevesselness_format != (char*)NULL ) freeStringList( &theVesselnessFileList );
    if ( p.thetransformation_format != (char*)NULL ) freeStringList( &theTransformationFileList );
    freeStringList( &theDetectionFileList );
    _ErrorParse( "unable to allocate tracker element list\n", 0);
  }

  theTrackerList.n = p.lastindex - p.firstindex + 1;

  for ( i=0; i<theTrackerList.n; i ++ ) {
    theTrackerList.data[i].index = p.firstindex + i;
    theTrackerList.data[i].detectionName = theDetectionFileList.data[i];
    theTrackerList.data[i].transformationName = theTransformationFileList.data[i];
    theTrackerList.data[i].vesselnessName = theVesselnessFileList.data[i];
    theTrackerList.data[i].readDetectionList.unit = p.detection_unit;
  }



  if ( BAL_ExtractBlTracks( &resList,
			    &theTrackerList,
			    &(p.tracking_parameters) ) != 1 ) {
    BAL_FreeBlTrackerElementList( &theTrackerList );
    if ( p.thevesselness_format != (char*)NULL ) freeStringList( &theVesselnessFileList );
    if ( p.thetransformation_format != (char*)NULL ) freeStringList( &theTransformationFileList );
    freeStringList( &theDetectionFileList );
    _ErrorParse( "error when processing\n", 0);
  }


  
  if ( BAL_WriteBlTrackList( &resList, p.restrack_name ) != 0 ) {
    BAL_FreeBlTrackList( &resList );
    _ErrorParse( "unable to write track file", 0 );
  }

  /*
  if ( BAL_ReadDetectionseBlTrackerElementList( &theTrackerList ) != 1 ) {
    BAL_FreeBlTrackerElementList( &theTrackerList );
    if ( p.thevesselness_format != (char*)NULL ) freeStringList( &theVesselnessFileList );
    if ( p.thetransformation_format != (char*)NULL ) freeStringList( &theTransformationFileList );
    freeStringList( &theDetectionFileList );
    _ErrorParse( "unable to read detections\n", 0);
  }
  */

  /* 1. build graph
     2. extract tracks
     1 et 2 a faire ensemble pour encapsuler la partie lemon ?
     3. select tracks 
   */

  
  BAL_FreeBlTrackerElementList( &theTrackerList );
  if ( p.thevesselness_format != (char*)NULL ) freeStringList( &theVesselnessFileList );
  if ( p.thetransformation_format != (char*)NULL ) freeStringList( &theTransformationFileList );
  freeStringList( &theDetectionFileList );



  /* ailleurs la convention etait que les tracks commencait a 0
   */
  {
    int t, i;
    for ( t=0; t<resList.n; t++ ) {
      for ( i=0; i<resList.data[t].detectionList.n; i++ )
	resList.data[t].detectionList.data[i].imageindex -= p.firstindex;
    }
  }


  if ( p.restrack_name != (char*)NULL ) {
    if ( BAL_WriteBlTrackList( &resList, p.restrack_name ) != 0 ) {
      BAL_FreeBlTrackList( &resList );
      _ErrorParse( "unable to write track file\n", 0 );
    }
  }
  
  BAL_FreeBlTrackList( &resList );



  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( _time_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, time_exit - time_init );

  if ( _clock_ ) 
    fprintf( stderr, "%s: elapsed time = %f\n", program, clock_exit - clock_init );


  return( 1 );
}











static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* general options 
       */
      if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
	   || (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( _verbose_ <= 0 )
	  _verbose_ = 1;
	else _verbose_ ++;
	BAL_IncrementVerboseInBalBiolibTracker();
      }
      
      else if ( (strcmp ( argv[i], "-no-verbose") == 0 && argv[i][11] == '\0')
	   || (strcmp ( argv[i], "-noverbose") == 0 && argv[i][11] == '\0')
	   || (strcmp ( argv[i], "-nv") == 0 && argv[i][3] == '\0') ) {
	_verbose_ = 0;
	BAL_SetVerboseInBalBiolibTracker( 0 );
      }
      else if ( strcmp ( argv[i], "--help" ) == 0 
		|| ( strcmp ( argv[i], "-help" ) == 0 && argv[i][5] == '\0' )
		|| ( strcmp ( argv[i], "--h" ) == 0 && argv[i][3] == '\0' )
		|| ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) ) {
	_ErrorParse( NULL, 1 );
      }
      
      else if ( (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') ) {
	_debug_ = 1;
      }
      
      else if ( strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0' ) {
	_time_ = 1;
      }
      else if ( strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0' ) {
	_time_ = 0;
      }
      else if ( strcmp ( argv[i], "-clock" ) == 0 && argv[i][6] == '\0' ) {
	_clock_ = 1;
      }
      else if ( strcmp ( argv[i], "-noclock" ) == 0 && argv[i][8] == '\0' ) {
	_clock_ = 0;
      }

      
      
      /* input related file names
       */
      else if ( strcmp ( argv[i], "-detection-format") == 0 
		|| strcmp ( argv[i], "-detection") == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -detection-format", 0 );
	p->thedetection_format = argv[i];
      }

      else if ( strcmp ( argv[i], "-vesselness-format") == 0 
		|| strcmp ( argv[i], "-vesselness") == 0  
		|| strcmp ( argv[i], "-vessel") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -vesselness-format", 0 );
	p->thevesselness_format = argv[i];
      }

      else if ( strcmp ( argv[i], "-transformation-format") == 0 
		|| strcmp ( argv[i], "-transformation") == 0  
		|| strcmp ( argv[i], "-trsf-format") == 0 
		|| strcmp ( argv[i], "-trsf") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -transformation-format", 0 );
	p->thetransformation_format = argv[i];
      }

      /* format related options
       */
      else if ( (strcmp ( argv[i], "-f" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-first" ) == 0 && argv[i][6] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -first ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->firstindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -first ...", 0 );
      }
      else if ( (strcmp ( argv[i], "-l" ) == 0 && argv[i][2] == '\0') 
		|| (strcmp ( argv[i], "-last" ) == 0 && argv[i][5] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -last ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->lastindex) );
	if ( status <= 0 ) _ErrorParse( "parsing -last ...", 0 );
      }

      
      /* search options
       */
      else if ( strcmp ( argv[i], "-detection-unit") == 0 
		|| (strcmp ( argv[i], "-dy") == 0 && argv[i][3] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "-detection-unit ...\n", 0 );
	if ( strcmp ( argv[i], "voxel") == 0 ) {
	  p->detection_unit = VOXEL_UNIT;
	}
	else if ( strcmp ( argv[i], "real") == 0 ) {
	  p->detection_unit = REAL_UNIT;
	}
	else {
	  _ErrorParse( "-detection-unit: unknown type ...\n", 0 );
	}
      }

     /* search options
       */
      else if ( strcmp ( argv[i], "-half-neighborhood-unit") == 0 
		|| (strcmp ( argv[i], "-hnu") == 0 && argv[i][4] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-unit ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_unit.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-unit ...", 0 );
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-unit ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_unit.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-unit ...", 0 );
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-unit ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_unit.z) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-unit ...", 0 );
      }

      else if ( strcmp ( argv[i], "-half-neighborhood-voxel") == 0 
		|| (strcmp ( argv[i], "-hnv") == 0 && argv[i][4] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-voxel ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_voxel.x) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-voxel ...", 0 );
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-voxel ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_voxel.y) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-voxel ...", 0 );
	i++;
	if ( i >= argc) _ErrorParse( "parsing -half-neighborhood-voxel ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.halfneighborhood_voxel.z) );
	if ( status <= 0 ) _ErrorParse( "parsing -half-neighborhood-voxel ...", 0 );
      }

      else if ( strcmp ( argv[i], "-search-type") == 0 ) {
	i++;
	if ( i >= argc) _ErrorParse( "-search-type ...\n", 0 );
	if ( strcmp ( argv[i], "force") == 0 ) {
	  p->tracking_parameters.searchType = _VOXEL_BRUTEFORCE_;
	}
	else if ( strcmp ( argv[i], "buckets") == 0 ) {
	  p->tracking_parameters.searchType = _VOXEL_BUCKETS_;
	}
	else {
	  _ErrorParse( "-search-type: unknown type ...\n", 0 );
	}
      }

      /* cost calculation
       */

      else if ( strcmp ( argv[i], "-min-2D-distance-voxel" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -min-2D-distance-voxel ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.costMin2Ddistance_voxel) );
	if ( status <= 0 ) _ErrorParse( "parsing -min-2D-distance-voxel ...", 0 );
      }

      else if ( strcmp ( argv[i], "-min-2D-distance-unit" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -min-2D-distance-unit ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.costMin2Ddistance_unit) );
	if ( status <= 0 ) _ErrorParse( "parsing -min-2D-distance-unit ...", 0 );
      }

      else if ( strcmp ( argv[i], "-cost-location") == 0 
		|| (strcmp ( argv[i], "-cl") == 0 && argv[i][3] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "-cost-location ...\n", 0 );
	if ( strcmp ( argv[i], "future") == 0 ) {
	  p->tracking_parameters.costLocation =  _FROM_NEXT_;
	}
	else if ( strcmp ( argv[i], "past") == 0 ) {
	  p->tracking_parameters.costLocation =  _FROM_PREV_;
	}
	else if ( strcmp ( argv[i], "both") == 0 ) {
	  p->tracking_parameters.costLocation =  _FROM_BOTH_;
	}
	else {
	  _ErrorParse( "-cost-location: unknown type ...\n", 0 );
	}
      }
      else if ( strcmp ( argv[i], "-use-retraction") == 0 
		|| strcmp ( argv[i], "-retraction") == 0 ) {
	p->tracking_parameters.costLocation =  _FROM_BOTH_;
      }

      else if ( strcmp ( argv[i], "-cost-computation") == 0 
		|| (strcmp ( argv[i], "-cc") == 0 && argv[i][3] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "-cost-computation ...\n", 0 );
	if ( strcmp ( argv[i], "exact") == 0 ) {
	  p->tracking_parameters.costCalculation =  _EXACT_;
	}
	else if ( strcmp ( argv[i], "approximate") == 0 
		  || strcmp ( argv[i], "approx") == 0) {
	  p->tracking_parameters.costCalculation =  _APPROXIMATE_;
	}
	else {
	  _ErrorParse( "-cost-computation: unknown type ...\n", 0 );
	}
      }

      

      /* tracks selection
       */
      else if ( strcmp ( argv[i], "-min-vesselness" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -min-vesselness ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.minVesselness) );
	if ( status <= 0 ) _ErrorParse( "parsing -min-vesselness ...", 0 );
      }
      else if ( strcmp ( argv[i], "-max-vesselness" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -max-vesselness ...\n", 0 );
	status = sscanf( argv[i], "%lf", &(p->tracking_parameters.maxVesselness) );
	if ( status <= 0 ) _ErrorParse( "parsing -max-vesselness ...", 0 );
      }
      else if ( strcmp ( argv[i], "-min-track-length" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -min-track-length ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->tracking_parameters.minTrackLength) );
	if ( status <= 0 ) _ErrorParse( "parsing -min-track-length ...", 0 );
      }
      else if ( strcmp ( argv[i], "-max-tracks" ) == 0  ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -max-tracks ...\n", 0 );
	status = sscanf( argv[i], "%d", &(p->tracking_parameters.maxTracks) );
	if ( status <= 0 ) _ErrorParse( "parsing -max-tracks ...", 0 );
      }

      /* unknown option
       */
      else {
	fprintf(stderr, "unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( p->restrack_name == (char*)NULL ) {
	p->restrack_name = argv[i];
      }
      else { 
	fprintf(stderr, "too many file names: '%s'\n",argv[i]);
      }
    }

  }
  
}





static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  if ( str != (char*)NULL )
    (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}





static void _InitParam( local_parameter *p )
{
  p->thedetection_format = (char*)NULL;
  p->thetransformation_format = (char*)NULL;
  p->thevesselness_format = (char*)NULL;

  p->restrack_name = (char*)NULL;

  p->firstindex = 0;
  p->lastindex = 0;

  p->detection_unit = VOXEL_UNIT;
 
  BAL_InitBlTrackerParameter( &(p->tracking_parameters) );

}





static double _GetTime() 
{
  struct timeval tv;
  gettimeofday(&tv, (void *)0);
  return ( (double) tv.tv_sec + tv.tv_usec*1e-6 );
}

static double _GetClock() 
{
  return ( (double) clock() / (double)CLOCKS_PER_SEC );
}
