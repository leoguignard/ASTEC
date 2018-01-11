/*************************************************************************
 *  -
 *
 * $Id: test-diameter.c,v 1.2 2004/06/10 09:23:33 greg Exp $
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Tue May 16 2000
 *
 *
 * ADDITIONS, CHANGES
 * 
 *
 */


#ifdef _INCLUDE_PELED_

#include  <stdlib.h>
#include  <stdio.h>
#include  <assert.h>
#include  <memory.h>
#include  <math.h>

#include  <gdiam.h>

#endif

#include <stdio.h>
#include <string.h>
#include <time.h>

#include <pick.h>
#include <rand.h>
#include <read.h>
#include <util.h>
#include <exact-diam.h>
#include <apprx-diam.h>

#include <peled-like.h>

#include <print.h>

#define _STRLENGTH_ 256


void _scan_coordinates( double** listOfPoints, int first, int last, int dim )
{
  int i,j;
  double c = 0;
  for (i=first; i<=last; i++ )
  for (j=0; j<dim; j++)  
    if ( c < listOfPoints[i][j] ) c = listOfPoints[i][j];
}


typedef enum {
  _TEST_,
  _EVAL_,
  _TIME_,
  _TIME_SCAN_,
  _COMP_
} enumMode;

typedef enum {
  _EXACT_,
  _QUADR_,
  _APPRX_,
  _PELED_
} enumMethod;


static char program[_STRLENGTH_];


static char *usage = "[-mode <test|eval|time|comp>] [-trials|-t %d] [-check|-nocheck]\n\
\n\
[-meth <exact|apprx|quadr|peled>] [-e|-epsilon %lf] [-er|-epsilon-ref %lf]\n\
[-tb | -tight-bounds | -ntb | -no-tight-bounds]\n\
[-remove-points-in-farthest | -rpif | -do-not-remove-points-in-farthest | -dnrpif]\n\
[-reduction-mode-diameter | -rmdiam %d]\n\
[-process-dbl-normals | -pdn | -do-not-process-dbl-normals | -dnpdn]\n\
[-reduction-mode-dblenorm | -rmdbnr %d]\n\
[-qforward | -qf | -qbackward | -qb]\n\
[-update-peled-diameter|-upd <quadratic | max-exact | max-exact-diam>]\n\
[-init-peled-diameter|-ipd <largest-dim | max-seg>]\n\
\n\
[-ply |-model %s]\n\
[-dist[rib] <c[ube] | e[llipse] | r[ellipse] | s[phere] | b[all] | d[iamcst]>\n\
[-dim %d] [-points|-p %d] [-vertices %d] [-random-rotation | -rr]\n\
\n\
[-init | -seed %ld]\n\
[-verbose | -v | -no-verbose | -nv]\n";

#ifdef _STATS_
static char *explain="\nCompute the diameter of sets of points and compute the complexity\n";
#else
static char *explain="\nCompute the diameter of sets of points\n";
#endif

static char *detail1 = "\
MODE\n\
 -mode:\n\
   test: infinite loop of tests (comparison with quadratic search)\n\
   eval: compute the diameter of a random distribution\n\
         or of points in a file\n\
         If -check is specified, compare with quadratic search\n\
   time: compute the average time to compute the diameter\n\
 -trials %d: number of tests to be done (not applicable for the TEST mode)\n\
 -check: relative to EVAL mode\n\
\n\
METHOD\n\
 -meth:\n\
    exact: our exact method\n\
    apprx: our approximate method\n\
    quadr: brute-force method\n\
 -epsilon %lf: precision for the approximate method.\n\
    The result consists in both lower and upper bound of the diameter\n\
    such that |upper-lower|<=epsilon\n\
 -tight-bounds: try to get tight bounds in approx method\n\
    else the upper bound may be often equal to lower+epsilon.\n\
\n\
 -remove-points-in-farthest: \n\
    dans boucle de la recherche iterative de la meilleure double normale\n\
    apres avoir calcule une double normale, on cherche le point le plus \n\
    de la sphere : lors de cette recherche on elimine les points a \n\
    l'interieur de la 'petite sphere'\n\
    while searching the farthest point of a sphere\n\
    (in the iterative serach of double normals) remove points from the\n\
    inner sphere\n\
\n\
 -reduction-mode-diameter %d: after the iterative search of double normals\n\
    try to remove points with the best found estimate\n\
    0 : no reduction\n\
    1 : reduction inside the inner sphere\n\
    2 : complete reduction\n\
    this is done just after the iterative search and before eventually\n\
    reducing the extremity candidates by looking at all found double normals\n\
\n";
static char *detail2 = "\
 -process-dbl-normals: reduce the number of extremities candidates\n\
\n\
 -reduction-mode-dblenorm  %d: within the reduction of extremities candidates\n\
    try to remove points with each double normal\n\
    0 : no reduction\n\
    1 : reduction inside the inner sphere\n\
    2 : complete reduction\n\
    ce n est fait que pour les points non candidat a l extremalite\n\
\n\
 -qbackward: process the found double normals in a backward order\n\
    ie from the largest (excepted the best one) to the smallest\n\
 -reduction-mode-dblenorm %d: reduction of points with each double normal\n\
\n\
 -update-peled-diameter: [quadratic | max-exact | max-exact-diam]\n\
\n\
DISTRIBUTION\n\
 -ply   %s: file containing 3D points (PLY format)\n\
 -model %s: id, other format\n\
 -distrib %s: random distribution\n\
   ball: in a ball\n\
   cube: in a cube\n\
   diamcst: in a 2D Reuleaux polygon\n\
   ellipse: on an ellipsoid\n\
   rellipse: on a regular ellipsoid\n\
   sphere: on a ball\n\
 -dim %d: dimension of the random set\n\
  -p %d: number of points of the random set\n\
 -vertices %d: gives the number of vertices of Reuleaux polygon\n\
   (only with '-dim 2')\n\
   '-vertices 1' -> triangle (2*1+1) = 3\n\
   '-vertices 2' -> pentagon (2*2+1) = 5\n\
\n\
PARAMETERS\n\
 -init |-seed %ld: seed for the sequence of pseudo-random integers\n";




static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s %s\n",program, usage, explain);
  if ( flag == 1 ) (void)fprintf(stderr,"%s%s",detail1,detail2);
  (void)fprintf(stderr,"Erreur : %s\n",str);
  exit(0);
}


static int _verbose_ = 0;






int main( int argc, char *argv[] )
{
  char modlname[_STRLENGTH_];
  
  enumDistribution typeDistribution = IN_CUBE;
#ifdef WIN32
  unsigned int seedRandom = time(0);
#else
  long int seedRandom = time(0);
#endif

  int _mode_;

  int _dim_ = 3;  
  double _max_diameter_ = 1.0;
  double _min_diameter_ = 0.2;
  int    _psommet_ = 1;
  int    _nbpoints_ = 1000;
  int    _trials_   = 1;
  double _epsilon_  = 0.01;
  double _epsilon_ref_  = 0.05;
  enumMode typeMode = _EVAL_;
  enumMethod typeMethod = _EXACT_;

  int i, status;

  double **listOfPoints = (double**)NULL;
  
  int c0, c1;
  int time0, time1, time2;

  typeCounter clock1, clock2;
  typeSegment pair1, pair2;
  double upper = 0.0;

  int checkTheResult = 0;
  int _nbpoints_max_ = 20000;
  int _nbpoints_min_ = 100;
  int _dim_max_      = 10;
  int _dim_min_      = 2;
  
  int do_compare_with_max_exact = 1;





  /* parse args
   */
  for ( i=1; i<argc; i++ ) {

    
    if ( argv[i][0] == '-' ) {
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
	_ErrorParse( "help message\n", 1 );
      }
      
      
      /* distributions
       */
      if ( (strcmp ( argv[i], "-distrib" ) == 0) ||
	   (strcmp ( argv[i], "-dist" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -distrib...\n", 0 );
	if ( (strcmp ( argv[i], "cube" ) == 0) ||
	     (strcmp ( argv[i], "c" ) == 0) ) {
	  typeDistribution = IN_CUBE;
	}
	else if ( (strcmp ( argv[i], "ellipse" ) == 0) ||
		  (strcmp ( argv[i], "e" ) == 0) ) {
	  typeDistribution = ON_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "rellipse" ) == 0) ||
		  (strcmp ( argv[i], "r" ) == 0) ) {
	  typeDistribution = ON_REG_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "gellipse" ) == 0) ||
		  (strcmp ( argv[i], "g" ) == 0) ) {
	  typeDistribution = ON_REG_ELLIPSOID;
	}
	else if ( (strcmp ( argv[i], "sphere" ) == 0) ||
		  (strcmp ( argv[i], "s" ) == 0) ) {
	  typeDistribution = ON_SPHERE;
	}
	 else if ( (strcmp ( argv[i], "ball" ) == 0) ||
		   (strcmp ( argv[i], "b" ) == 0) ) {
	   typeDistribution = IN_SPHERE;
	 } 
	else if ( (strcmp ( argv[i], "diamcst" ) == 0) ||
		  (strcmp ( argv[i], "d" ) == 0) ) {
	  typeDistribution = IN_CST_DIAMETER;
	} 
	else  {
	  sprintf( modlname, "unknown distribution = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      } 
    
      
      else if ( strcmp ( argv[i], "-model" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -model...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_MODEL;
      }
      else if ( strcmp ( argv[i], "-ply" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -ply...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_PLY_MODEL;
      }
      else if ( strcmp ( argv[i], "-pts" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -pts...\n", 0 );
	sprintf( modlname, "%s", argv[i] );
	typeDistribution = IN_PTS_MODEL;
      }
      

      /* number of points
       */
      else if ( (strcmp ( argv[i], "-p" ) == 0) ||
		(strcmp ( argv[i], "-points" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -points...\n", 0 );
	status = sscanf( argv[i],"%d",&_nbpoints_ );
	if ( status <= 0 ) _ErrorParse( "parsing -points...\n", 0 );
      }

      /* epsilon for approximation
       */
      else if ( (strcmp ( argv[i], "-e" ) == 0) ||
		(strcmp ( argv[i], "-epsilon" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -epsilon...\n", 0 );
	status = sscanf( argv[i],"%lf",&_epsilon_ );
	if ( status <= 0 ) _ErrorParse( "parsing -epsilon...\n", 0 );
      }

      else if ( (strcmp ( argv[i], "-er" ) == 0) ||
		(strcmp ( argv[i], "-epsilon-ref" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -epsilon-ref...\n", 0 );
	status = sscanf( argv[i],"%lf",&_epsilon_ref_ );
	if ( status <= 0 ) _ErrorParse( "parsing -epsilon-ref...\n", 0 );
      }

      else if ( (strcmp ( argv[i], "-tb" ) == 0) ||
		(strcmp ( argv[i], "-tight-bounds" ) == 0) ) {
	_DoTryToGetTightBounds();
      }
      else if ( (strcmp ( argv[i], "-ntb" ) == 0) ||
		(strcmp ( argv[i], "-no-tight-bounds" ) == 0) ) {
	_DoNotTryToGetTightBounds();
      }


      /* number of trials
       */
      else if ( (strcmp ( argv[i], "-t" ) == 0) ||
		(strcmp ( argv[i], "-trials" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -trials...\n", 0 );
	status = sscanf( argv[i],"%d",&_trials_ );
	if ( status <= 0 ) _ErrorParse( "parsing -trials...\n", 0 );
      }

      /* dimension
       */
      else if ( strcmp ( argv[i], "-dim" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -dim...\n", 0 );
	status = sscanf( argv[i],"%d",&_dim_ );
	if ( status <= 0 ) _ErrorParse( "parsing -dim...\n", 0 );
      }

      /* seed for random distribution
       */
      else if ( (strcmp ( argv[i], "-init" ) == 0) ||
		(strcmp ( argv[i], "-seed" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -seed...\n", 0 );
#ifdef WIN32
	status = sscanf( argv[i],"%u",&seedRandom );
#else
	status = sscanf( argv[i],"%ld",&seedRandom );
#endif
	if ( status <= 0 ) _ErrorParse( "parsing -seed...\n", 0 );
      }

      /* number of vertices for Reuleaux polygons = 2*p+1
	 _psommet_ = 1 => triangle.
       */
      else if ( strcmp ( argv[i], "-vertices" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -vertices...\n", 0 );
	status = sscanf( argv[i],"%d",&_psommet_ );
	if ( status <= 0 ) _ErrorParse( "parsing -vertices...\n", 0 );
      }
      
      else if ( (strcmp ( argv[i], "-random-rotation" ) == 0) ||
		(strcmp ( argv[i], "-rr" ) == 0) ) {
	_VerboseInPick(); 
	_DoApplyARandomRotation();
      }


      /* check the result with quadratic search
       */
      else if ( strcmp ( argv[i], "-check" ) == 0 ) {
	checkTheResult = 1;
      }
      else if ( strcmp ( argv[i], "-nocheck" ) == 0 ) {
	checkTheResult = 0;
      }
      
      
      /* try to remove points while processing 
       */
      else if ( (strcmp ( argv[i], "-remove-points-in-farthest" ) == 0) ||
		(strcmp ( argv[i], "-rpif" ) == 0) ) {
	_SetReductionModeInIterative( 1 );
      }
      else if ( (strcmp ( argv[i], "-do-not-remove-points-in-farthest" ) == 0) ||
		(strcmp ( argv[i], "-dnrpif" ) == 0) ) {
	_SetReductionModeInIterative( 0 );
      }


      /* try to remove points while processing 
       */
      else if ( (strcmp ( argv[i], "-process-dbl-normals" ) == 0) ||
		(strcmp ( argv[i], "-pdn" ) == 0) ) {
	_DoTryToReduceQ();
      }
      else if ( (strcmp ( argv[i], "-do-not-process-dbl-normals" ) == 0) ||
		(strcmp ( argv[i], "-dnpdn" ) == 0) ) {
	_DoNotTryToReduceQ();
      }

      else if ( (strcmp ( argv[i], "-reduction-mode-diameter" ) == 0) ||
		(strcmp ( argv[i], "-rmdiam" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -reduction-mode-diameter...\n", 0 );  
	status = sscanf( argv[i],"%d",&_mode_ );
	if ( status <= 0 ) _ErrorParse( "parsing -reduction-mode-diameter...\n", 0 );
	_SetReductionModeOfDiameter( _mode_ );
      }
      else if ( (strcmp ( argv[i], "-reduction-mode-dblenorm" ) == 0) ||
		(strcmp ( argv[i], "-rmdbnr" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -reduction-mode-dblenorm...\n", 0 );  
	status = sscanf( argv[i],"%d",&_mode_ );
	if ( status <= 0 ) _ErrorParse( "parsing -reduction-mode-dblenorm...\n", 0 );
	_SetReductionModeOfDbleNorm( _mode_ );
      }
      

      else if ( (strcmp ( argv[i], "-qforward" ) == 0) ||
		(strcmp ( argv[i], "-qf" ) == 0) ) {
	_SetQscanToForward();
      }
      else if ( (strcmp ( argv[i], "-qbackward" ) == 0) ||
		(strcmp ( argv[i], "-qb" ) == 0) ) {
	_SetQscanToBackward();
      }
      

      else if ( strcmp ( argv[i], "-init-peled-diameter" ) == 0 ||
		strcmp ( argv[i], "-ipd" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -init-peled-diameter...\n", 0 );  
	if ( strcmp ( argv[i], "largest-dim" ) == 0 ) {
	  set_init_peled_diameter_to_largest_dim();
	}
	else if ( strcmp ( argv[i], "max-seg" ) == 0 ) {
	  set_init_peled_diameter_to_maximal_seg();
	} else  {
	  sprintf( modlname, "unknown option for -init-peled-diameter = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      }
      
      else if ( strcmp ( argv[i], "-update-peled-diameter" ) == 0 ||
		strcmp ( argv[i], "-upd" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -update-peled-diameter...\n", 0 );  
	if ( strcmp ( argv[i], "quadratic" ) == 0 ) {
	  set_update_peled_diameter_to_quadratic();
	}
	else if ( strcmp ( argv[i], "max-exact" ) == 0 ) {
	  set_update_peled_diameter_to_max_exact();
	}
	else if ( strcmp ( argv[i], "max-exact-diam" ) == 0 ) {
	  set_update_peled_diameter_to_max_exact_with_diameter();
	} else  {
	  sprintf( modlname, "unknown option for -update-peled-diameter = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      }
      
    

      
      /* mode
       */
      else if ( strcmp ( argv[i], "-mode" ) == 0 ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -mode...\n", 0 );  
	if ( strcmp ( argv[i], "test" ) == 0 ) {
	  typeMode = _TEST_;
	}
	else if ( strcmp ( argv[i], "eval" ) == 0 ) {
	  typeMode = _EVAL_;
	}
	else if ( strcmp ( argv[i], "time" ) == 0 ) {
	  typeMode = _TIME_;
	}
	else if ( strcmp ( argv[i], "scan" ) == 0 ) {
	  typeMode = _TIME_SCAN_;
	}
	else if ( strcmp ( argv[i], "comp" ) == 0 ) {
	  typeMode = _COMP_;
	} else  {
	  sprintf( modlname, "unknown mode = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      }

#ifdef _INCLUDE_PELED_
      else if ( strcmp ( argv[i], "-no-max-exact" ) == 0 ) {
	do_compare_with_max_exact = 0;
      }
#endif

      else if ( (strcmp ( argv[i], "-method" ) == 0) ||
		(strcmp ( argv[i], "-meth" ) == 0) ) {
	i++;
	if ( i >= argc)    _ErrorParse( "parsing -method...\n", 0 );  
	if ( strcmp ( argv[i], "exact" ) == 0 ) {
	  typeMethod = _EXACT_;
	}
	else if ( strcmp ( argv[i], "quadr" ) == 0 ) {
	  typeMethod = _QUADR_;
	} 
	else if ( strcmp ( argv[i], "apprx" ) == 0 ) {
	  typeMethod = _APPRX_;
	}
	else if ( strcmp ( argv[i], "approx" ) == 0 ) {
	  typeMethod = _APPRX_;
	}
	else if ( strcmp ( argv[i], "peled" ) == 0 ) {
	  typeMethod = _PELED_;
	} else  {
	  sprintf( modlname, "unknown method = %s\n", argv[i] );
	  _ErrorParse( modlname, 0);
	}
      }

      else if ( (strcmp ( argv[i], "-local-no-verbose" ) == 0) ||
		(strcmp ( argv[i], "-lnv" ) == 0) ) {
	_verbose_ = 0;
      }
      else if ( (strcmp ( argv[i], "-local-verbose" ) == 0) ||
		(strcmp ( argv[i], "-lv" ) == 0) ) {
	if ( _verbose_ <= 0 )
	  _verbose_ = 1;
	else 
	  _verbose_ ++;
      }
      else if ( (strcmp ( argv[i], "-verbose" ) == 0) ||
		(strcmp ( argv[i], "-v" ) == 0) ) {
	_verbose_ = 1;
	_VerboseInApprxDiam();
	_VerboseWhenReducing();
	_VerboseInPeledLike();
      }
      else if ( (strcmp ( argv[i], "-no-verbose" ) == 0) ||
		(strcmp ( argv[i], "-nv" ) == 0) ) {
	_verbose_ = 0;
	_NoVerboseInApprxDiam();
	_NoVerboseWhenReducing();
	_NoVerboseInPeledLike();
      }


      /*--- option inconnue ---*/
      else {
	sprintf( modlname, "unknown option %s\n", argv[i] );
	_ErrorParse( modlname, 0);
      }
      
      
    } 

    /* arg without '-'
     */
    else if ( argv[i][0] != 0 ) {
      sprintf( modlname, "unknown option %s\n", argv[i] );
      _ErrorParse( modlname, 0);
    }
  }






#ifdef _STATS_
  typeMode = _TIME_;
  checkTheResult = 0;
  _InitScalarProductCounter();
#endif
  

  switch ( typeDistribution ) {
  case IN_PLY_MODEL :
  case IN_MODEL :
  case IN_PTS_MODEL :
    break;
  default :
    _SetRandomSeed( seedRandom );
  }

  _InitCounter( &clock1 );
  _InitCounter( &clock2 );



  switch ( typeMode ) {
  case _TEST_ :
    {
      int n, j, success, failure;
      double error=1e-12;

      success = failure = n = 0;

      fprintf( stderr, "%%   random seed      = %ld\n", _GetRandomSeed() );

      do {
	n ++;

	_nbpoints_ = _GetRandomIntNb( _nbpoints_min_, _nbpoints_max_ );
	_dim_      = _GetRandomIntNb( _dim_min_, _dim_max_ );
	if ( _dim_ == 2 ) 
	  typeDistribution = (enumDistribution)_GetRandomIntNb( IN_CUBE, IN_CST_DIAMETER );
	else
	  typeDistribution = (enumDistribution)_GetRandomIntNb( IN_CUBE, ON_REG_ELLIPSOID );

	_mode_ = _GetRandomIntNb( 0, 1 );
	_SetReductionModeInIterative( _mode_ );
	_mode_ = _GetRandomIntNb( 0, 2 );
	_SetReductionModeOfDiameter( _mode_ );
	_mode_ = _GetRandomIntNb( 0, 2 );
	_SetReductionModeOfDbleNorm( _mode_ );

	_mode_ = _GetRandomIntNb( 0, 1 );
	if ( _mode_ == 0 ) _DoNotTryToReduceQ();
	else               _DoTryToReduceQ();

	_mode_ = _GetRandomIntNb( 0, 1 );
	if ( _mode_ == 0 ) _SetQscanToForward();
	else               _SetQscanToBackward();

	_mode_ = _GetRandomIntNb( 0, 1 );
	if ( _mode_ == 0 ) _DoNotTryToGetTightBounds();
	else               _DoTryToGetTightBounds();
		     

	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	_ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );

	_epsilon_ = _epsilon_ref_ * _GetRandomDoubleNb();
	

	if ( _verbose_ ) {
	  fprintf( stderr, " %d points in dim=%d      \n", _nbpoints_, _dim_ );
	}

	upper = 0.0;
	switch ( typeMethod ) {
	default :
	case _EXACT_ :
	  (void)_ExactDiameterInOneList( &pair1, listOfPoints, 
					 0, _nbpoints_-1, _dim_ );
	  (void)_QuadraticDiameterInOneList( &pair2, listOfPoints, 
					     0, _nbpoints_-1, _dim_ );
	  break;
	case _PELED_ :
	  /* if ( n == 3 ) _VerboseInPeledLike(); */
	  (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					 0, _nbpoints_-1, _dim_ );
	  (void)_ExactDiameterInOneList( &pair2, listOfPoints, 
					     0, _nbpoints_-1, _dim_ );
	  /*
	  _NoVerboseInPeledLike();
	  _VerboseInPeledLike();
	  */
	  break;
	case _APPRX_ :
	  upper = _EstimeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_, _epsilon_ );
	  (void)_ExactDiameterInOneList( &pair2, listOfPoints, 
					 0, _nbpoints_-1, _dim_ );
	  break;
	}
	
	switch ( typeMethod ) {
	default :
	case _EXACT_ :
	  if ( ( pair1.extremity1 != pair2.extremity1 || pair1.extremity2 != pair2.extremity2 ) &&
	       ( pair1.extremity1 != pair2.extremity2 || pair1.extremity2 != pair2.extremity1 ) ) {
	    failure ++;
	    
	    _PrintEnv( stderr, _GetRandomSeed(), _nbpoints_, _max_diameter_,
		       _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
	    _PrintSegment( stderr, &pair2, "Points realizing the diameter (success)", _dim_ );
	    _PrintSegment( stderr, &pair1, "Points realizing the diameter (failure)", _dim_ );
	    
	  } else {
	    success ++;
	  }
	  break;
	case _APPRX_ :
	  printf( " dim=%2d pts=%5d E=%-9.7f | min=%-8.6f real=%-8.6f max=%-8.6f | max/real-1 - max/min-1 = %8.6f - %8.6f\n", 
		  _dim_, _nbpoints_, _epsilon_,
		  sqrt( pair1.squareDiameter), sqrt( pair2.squareDiameter), 
		  sqrt( upper ),
		  sqrt( upper ) / sqrt( pair2.squareDiameter) - 1.0,
		  sqrt( upper ) / sqrt( pair1.squareDiameter) - 1.0 );
	  
	  if ( sqrt( pair1.squareDiameter) > sqrt( pair2.squareDiameter) ||
	       (sqrt( upper ) / sqrt( pair1.squareDiameter) - 1.0 > _epsilon_ + error) ||
	       (sqrt( upper ) / sqrt( pair2.squareDiameter) - 1.0 > _epsilon_ + error) ) {
	    failure ++;
	    {
	      char name[50];
	      FILE *fp;

	      sprintf( name, "failure.%d.txt", failure );
	      fp = fopen(  name, "w" );	    
	      fprintf( fp, "\n" );
	      fprintf( fp, " #n = %6d - success %6d - failure %6d\n", n, success, failure );
	      fprintf( fp, " epsilon = %12.10f\n", _epsilon_ );
	      fprintf( fp, " dim=%2d pts=%5d E=%-9.7f | min=%-8.6f real=%-8.6f max=%-8.6f | max/real-1 - max/min-1 = %8.6f - %8.6f\n", 
		       _dim_, _nbpoints_, _epsilon_,
		       sqrt( pair1.squareDiameter), sqrt( pair2.squareDiameter), 
		       sqrt( upper ),
		       sqrt( upper ) / sqrt( pair2.squareDiameter) - 1.0,
		       sqrt( upper ) / sqrt( pair1.squareDiameter) - 1.0 );

	      _PrintEnv( fp, _GetRandomSeed(), _nbpoints_, _max_diameter_,
			 _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
	      
	      _PrintSegment( fp, &pair2, "Points realizing the diameter (success)", _dim_ );
	      _PrintSegment( fp, &pair1, "Points realizing the diameter (low bound)", _dim_ );
	      fprintf( fp, "\n" );
	      if ( (sqrt( upper ) > sqrt( pair1.squareDiameter) + _epsilon_ + error) ) {
		fprintf( fp, "ERROR: upper(=%12.10f) - lower(=%12.10f) > epsilon((=%12.10f)\n",
			 sqrt( upper ), sqrt( pair1.squareDiameter), _epsilon_ );
	      }
	      if ( (upper  < pair2.squareDiameter) ) {
		fprintf( fp, "ERROR: upper(=%12.10f) < real(=%12.10f)\n",
			 sqrt( upper ), sqrt( pair2.squareDiameter ) );
	      }
	      if ( (pair2.squareDiameter < pair1.squareDiameter) ) {
		fprintf( fp, "ERROR: real(=%12.10f) < lower(=%12.10f)\n",
			 sqrt( pair2.squareDiameter ), sqrt( pair1.squareDiameter ) );
	      }
	      fprintf( fp, "\n" );
	      fclose( fp );

	      sprintf( name, "failure.%d.pts",failure );
	      fp = fopen(  name, "w" );	  
	      fprintf( fp, "#\n" );
	      fprintf( fp, "#" );
	      for (i=0; i<argc; i++ ) fprintf( fp, " %s", argv[i] );
	      fprintf( fp, "\n" );
	      fprintf( fp, "#\n" );
	      fprintf( fp, "random seed  = %20ld\n", _GetRandomSeed() );
	      fprintf( fp, "random calls = %20ld\n", _GetRandomCalls() );
	      fprintf( fp, "points       = %d\n", _nbpoints_ );
	      fprintf( fp, "dim          = %d\n", _dim_ );
	      for (i=0; i<_nbpoints_; i++ ) {
		for (j=0; j<_dim_; j++ ) {
		  fprintf( fp, "%12.10f", listOfPoints[i][j] );
		  if ( j == _dim_-1) fprintf( fp, "\n" );
		  else               fprintf( fp, " " );
		}
	      }
	      fclose( fp );
	    }
	  } else {
	    success ++;
	  }
	}
	fprintf( stderr, " success %6d - failure %6d\r", success, failure );
	
	free( listOfPoints );
	
      } while( 1 );
      
    }
    break;













  case _COMP_ :

    {
#ifdef _INCLUDE_PELED_
    GPointPair   pair;
#endif
    double * points;

#ifdef WIN32
    fprintf( stderr, "Random seed = %u\n", seedRandom );
#else
    fprintf( stderr, "Random seed = %ld\n", seedRandom );
#endif
    printf( "\n" );
    _PrintEnv( stdout, _GetRandomSeed(), _nbpoints_, _max_diameter_,
	       _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
    printf( "\n" );

    

    for ( i = 0; i < _trials_; i++ ) {

      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	points = listOfPoints[0];
	break;
      case IN_PTS_MODEL:
	if ( i == 0 ) {
	  listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
	  points = listOfPoints[0];
	}
	break;
      case IN_PLY_MODEL:
	if ( i == 0 ) {
	  listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
	  points = listOfPoints[0];
	}
	break;
      case IN_MODEL:
	if ( i == 0 ) {
	  listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
	  points = listOfPoints[0];
	}
	break;
      }
      _VerboseInPick();
      _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );



      

      if ( do_compare_with_max_exact ) {
	time0 = time( NULL );   c0 = clock();
	(void)_ExactDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	c1 = clock();   time1 = time( NULL );
	_PrintSegment( stdout, &pair1, "Points realizing the diameter (Max Exact)", _dim_ );
	printf( "Time (Max Exact) = %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
	printf( "\n\n\n" );
      }




#ifdef _INCLUDE_PELED_
      if ( _dim_ == 3 ) {
	time0 = time( NULL );   c0 = clock();
	pair = gdiam_approx_diam_pair( (gdiam_real *)points, _nbpoints_, 0.0 );
	c1 = clock();   time1 = time( NULL );
	pair1.extremity1 = pair.p;
	pair1.extremity2 = pair.q;
	pair1.squareDiameter = pair.distance * pair.distance;
	_PrintSegment( stdout, &pair1, "Points realizing the diameter (Peled)", 3 );
	printf( "Time (Peled    )= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
	printf( "\n\n\n" );
      }
#endif





      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_largest_dim();
      set_update_peled_diameter_to_quadratic();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (000)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 000)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );


      if ( 0 ) {

      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_largest_dim();
      set_update_peled_diameter_to_quadratic();
      do_not_allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (001)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 001)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );

      }

      printf( "\n\n\n" );





      
      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_largest_dim();
      set_update_peled_diameter_to_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (010)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 010)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_largest_dim();
      set_update_peled_diameter_to_max_exact_with_diameter();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (020)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 020)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      
       _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_largest_dim();
      set_update_peled_diameter_to_smart_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (030)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 030)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      





      printf( "\n\n\n" );
      printf( "\n\n\n" );


      
      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_maximal_seg();
      set_update_peled_diameter_to_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (110)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 110)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_iterated_maximal_seg();
      set_update_peled_diameter_to_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (111)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 111)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      



      printf( "\n\n\n" );

      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_maximal_seg();
      set_update_peled_diameter_to_max_exact_with_diameter();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (120)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 120)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_iterated_maximal_seg();
      set_update_peled_diameter_to_max_exact_with_diameter();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (121)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 121)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      










      printf( "\n\n\n" );

      
      _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_maximal_seg();
      set_update_peled_diameter_to_smart_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (130)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 130)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

       _SetRandomSeed( seedRandom );
      set_evolution_peled_min_nb_points_to_grow();
      set_init_peled_diameter_to_iterated_maximal_seg();
      set_update_peled_diameter_to_smart_max_exact();
      allow_to_split_small_box();

      time0 = time( NULL );   c0 = clock();
      (void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					    0, _nbpoints_-1, _dim_ );
      c1 = clock();   time1 = time( NULL );
      _PrintSegment( stdout, &pair1, "Peled like (131)", _dim_ );
      print_peled_env();
      printf( "Time (Peled 131)= %f (%f) sec.\n\n", (c1-c0)/(double)CLOCKS_PER_SEC, (double)(time1-time0) );
      

      











    


      







      









      switch( typeDistribution ) {
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	if ( i != _trials_ - 1 )
	  break;
      default :
	free( listOfPoints );
      }      
    }

    }
    break;



















  default :
  case _EVAL_ :



    for ( i = 0; i < _trials_; i++ ) {

      if ( 0 && i == 8 ) {
	_verbose_ = 2;
	_VerboseInApprxDiam();
	_VerboseWhenReducing();
      }
      
      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	break;
      case IN_PTS_MODEL:
	if ( i == 0 )
	  listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
	break;
      case IN_PLY_MODEL:
	if ( i == 0 )
	  listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
	break;
      case IN_MODEL:
	if ( i == 0 )
	  listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
	break;
      }
      _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );
      
    

      c0 = clock();
      switch ( typeMethod ) {
      default :
      case _EXACT_ :
	(void)_ExactDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	break;
      case _PELED_ :
	/*
	(void)_Peled_1_LikeDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	_PrintSegment( stderr, &pair1, "Points realizing the diameter (Peled 1)", _dim_ );
	*/
	(void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
					      0, _nbpoints_-1, _dim_ );
	/*
	_PrintSegment( stderr, &pair1, "Points realizing the diameter (Peled 2)", _dim_ );
	*/
	break;
      case _APPRX_ :
	upper = _EstimeDiameterInOneList( &pair1, listOfPoints, 
					  0, _nbpoints_-1, _dim_, _epsilon_ );
	break;
      }
      c1 = clock();
      _AddToCounter( &clock1, c1-c0 );

      
      

      if ( checkTheResult  ) {
	/* if ( checkTheResult || _trials_ == 1 ) { */
	
	c0 = clock();
	(void)_QuadraticDiameterInOneList( &pair2, listOfPoints, 
					   0, _nbpoints_-1, _dim_ );
	c1 = clock();
	_AddToCounter( &clock2, c1-c0 );
      }
      

      
      if ( _verbose_ || _trials_ == 1 ) {
	if ( _verbose_ >= 2 ) fprintf( stdout, "\n" );
	switch ( typeMethod ) {
	default :
	case _EXACT_ :	
	case _PELED_ :	
	  _PrintSegment( stdout, &pair1, "Points realizing the diameter (1)", _dim_ );
	break;
	case _APPRX_ :
	  fprintf( stderr, "   |%-8.6f (lowr) - ", sqrt( pair1.squareDiameter) );
	  if ( checkTheResult ) 
	    fprintf( stderr, "%-8.6f (true) - ", sqrt( pair2.squareDiameter) );
	  
	  fprintf( stderr, "%-8.6f (uppr)| = ", sqrt( upper ) );
	  if ( checkTheResult ) 
	    fprintf( stderr, "%8.6f + %8.6f = ",
		     sqrt( pair2.squareDiameter) - sqrt( pair1.squareDiameter),
		     sqrt( upper ) - sqrt( pair2.squareDiameter) );
	  fprintf( stderr, "%f\n", sqrt( upper ) - sqrt( pair1.squareDiameter) );
	  if ( _verbose_ >= 2 || _trials_ == 1 ) 
	    _PrintSegment( stderr, &pair1, "Points realizing the diameter (low bound)", _dim_ );
	}
	if ( _verbose_ >= 2 || _trials_ == 1 ) {
	  if ( checkTheResult )  
	    _PrintSegment( stdout, &pair2, "Points realizing the diameter (2)", _dim_ );
	}
      }
    
      switch( typeDistribution ) {
      default :
	free( listOfPoints );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }
    }

  
    _PrintEnv( stdout, _GetRandomSeed(), _nbpoints_, _max_diameter_,
	       _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
    
    fprintf( stdout, "   computation time (1) = %g sec. on %d trials\n", 
	     _GetCounterAverage(&clock1, _trials_ )/(double)CLOCKS_PER_SEC, _trials_ );

    
    if ( checkTheResult ) {
      fprintf( stdout, "   computation time (2) = %g sec. on %d trials \n", 
	       _GetCounterAverage(&clock2, _trials_ )/(double)CLOCKS_PER_SEC, _trials_ );
    }
    
    fprintf( stdout, "\n" );

    break;








  case _TIME_ :
    

    switch( typeDistribution ) {
    default :
      break;
    case IN_PTS_MODEL:
      listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
      break;
    case IN_PLY_MODEL:
      listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
      break;
    case IN_MODEL:
      listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
      break;
    }


    time0 = time( NULL );
    c0 = clock();
    
    for ( i = 0; i < _trials_; i++ ) {

#ifdef _STATS_
      if ( _trials_ < 100 )
	fprintf( stderr, "1: %2d/%d \r", i, _trials_ );
      else if ( _trials_ < 1000 )
	fprintf( stderr, "1: %3d/%d \r", i, _trials_ );
      else if ( _trials_ < 10000 )
	fprintf( stderr, "1: %4d/%d \r", i, _trials_ );
      else
	fprintf( stderr, "1: %5d/%d \r", i, _trials_ );
#endif

      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	free( listOfPoints );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }
    }
    _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );

    c1 = clock();
    time1 = time( NULL );

    _AddToCounter( &clock1, c1-c0 );

    c0 = clock();

    for ( i = 0; i < _trials_; i++ ) {

#ifdef _STATS_
      if ( _trials_ < 100 )
	fprintf( stderr, "2: %2d/%d \r", i, _trials_ );
      else if ( _trials_ < 1000 )
	fprintf( stderr, "2: %3d/%d \r", i, _trials_ );
      else if ( _trials_ < 10000 )
	fprintf( stderr, "2: %4d/%d \r", i, _trials_ );
      else
	fprintf( stderr, "2: %5d/%d \r", i, _trials_ );
#endif

      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }

      switch ( typeMethod ) {
      default :
      case _EXACT_ :
	(void)_ExactDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	break;
      case _PELED_ :
	/*
	(void)_Peled_1_LikeDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	*/
	(void)_Peled_LikeDiameterInOneList( &pair1, listOfPoints, 
				       0, _nbpoints_-1, _dim_ );
	break;
      case _APPRX_ :
	upper = _EstimeDiameterInOneList( &pair1, listOfPoints, 
					  0, _nbpoints_-1, _dim_, _epsilon_ );
	break;
      case _QUADR_ :
	(void)_QuadraticDiameterInOneList( &pair1, listOfPoints, 
					   0, _nbpoints_-1, _dim_ );
      }
      
      switch( typeDistribution ) {
      default :
	free( listOfPoints );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }
      
    }

    c1 = clock();
    time2 = time( NULL );

    _AddToCounter( &clock2, c1-c0 );

    switch( typeDistribution ) {
    default :
      break;
    case IN_PTS_MODEL:
    case IN_PLY_MODEL:
    case IN_MODEL:
      free( listOfPoints );
      break;
    }

    _PrintEnv( stdout, _GetRandomSeed(), _nbpoints_, _max_diameter_,
	       _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
    
    fprintf( stdout, "   computation time = %g sec. on %d trials\n", 
	     (_GetCounterAverage(&clock2, _trials_ ) - 
	      _GetCounterAverage(&clock1, _trials_ ) )/(double)CLOCKS_PER_SEC, _trials_ );
    fprintf( stdout, "   computation time = %g sec. on %d trials\n",
	     (time2 - 2*time1 + time0)/(double)_trials_, _trials_ );
    fprintf( stdout, "\n" );

#ifdef _STATS_
    fprintf( stdout, "   average number of dot products per trial = %g\n", 
	     _GetScalarProductAverage( _trials_ ) );
    fprintf( stdout, "   average complexity = %g (dot product per point)\n", 
	     _GetScalarProductAverage( _trials_*_nbpoints_ ) );
    fprintf( stdout, "\n" );
#endif

    break;


  case _TIME_SCAN_ :
    

    switch( typeDistribution ) {
    default :
      break;
    case IN_PTS_MODEL:
      listOfPoints = (double**)_ReadPTSModel( modlname, &_nbpoints_, &_dim_ );
      break;
    case IN_PLY_MODEL:
      listOfPoints = (double**)_ReadPLYModel( modlname, &_nbpoints_, &_dim_ );
      break;
    case IN_MODEL:
      listOfPoints = (double**)_ReadModel( modlname, &_nbpoints_, &_dim_ );
      break;
    }


    time0 = time( NULL );
    c0 = clock();
    
    for ( i = 0; i < _trials_; i++ ) {
      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	free( listOfPoints );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }
    }
    
    _ApplyARandomRotation( listOfPoints, _nbpoints_, _dim_ );

    c1 = clock();
    time1 = time( NULL );

    _AddToCounter( &clock1, c1-c0 );

    c0 = clock();

    for ( i = 0; i < _trials_; i++ ) {

      switch( typeDistribution ) {
      default :
	listOfPoints = (double**)_PickPoints( _nbpoints_, _max_diameter_,
					      _min_diameter_, _psommet_, 
					      typeDistribution, _dim_ );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }

      _scan_coordinates( listOfPoints, 
			 0, _nbpoints_-1, _dim_ );
      
      switch( typeDistribution ) {
      default :
	free( listOfPoints );
	break;
      case IN_PTS_MODEL:
      case IN_PLY_MODEL:
      case IN_MODEL:
	break;
      }
      
    }

    c1 = clock();
    time2 = time( NULL );

    _AddToCounter( &clock2, c1-c0 );

    switch( typeDistribution ) {
    default :
      break;
    case IN_PTS_MODEL:
    case IN_PLY_MODEL:
    case IN_MODEL:
      free( listOfPoints );
      break;
    }

    _PrintEnv( stdout, _GetRandomSeed(), _nbpoints_, _max_diameter_,
	       _min_diameter_, _psommet_, typeDistribution, modlname, _dim_  );
    
    fprintf( stdout, "   computation time = %g sec. on %d trials\n", 
	     (_GetCounterAverage(&clock2, _trials_ ) - 
	      _GetCounterAverage(&clock1, _trials_ ) )/(double)CLOCKS_PER_SEC, _trials_ );
    fprintf( stdout, "   computation time = %g sec. on %d trials\n",
	     (time2 - 2*time1 + time0)/(double)_trials_, _trials_ );
    fprintf( stdout, "\n" );

#ifdef _STATS_
    fprintf( stdout, "   average number of dot products per trial = %g\n", 
	     _GetScalarProductAverage( _trials_ ) );
    fprintf( stdout, "   average complexity = %g (dot product per point)\n", 
	     _GetScalarProductAverage( _trials_*_nbpoints_ ) );
    fprintf( stdout, "\n" );
#endif

    break;
  }

  return( 0 );
}
