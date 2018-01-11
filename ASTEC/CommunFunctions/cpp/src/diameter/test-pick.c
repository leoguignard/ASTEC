#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <stdio.h>

#include <string.h>
#include <time.h>

#include <pick.h>
#include <rand.h>





#define _STRLENGTH_ 256

static char program[_STRLENGTH_];

static char *usage = "[-init | -seed %ld]\n\
[ -dist[rib] <c[ube] | e[llipse] | r[ellipse] | s[phere] | b[all] | d[iamcst]>]\n\
[ -points | -n %d] [-p %d]\n";

static char *detail = "";

static void _ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}


int main( int argc, char *argv[] )
{
  char filename[256];
  enumDistribution typeDistribution = IN_CUBE;
  long int seedRandom = time(0);

  int _dim_ = 2;  
  double _max_diameter_ = 1.0;
  double _min_diameter_ = 0.2;
  int    _psommet_ = 1;
  int    _nbpoints_ = 1000;
  
  int i, status;

  filename[0] = '\0';


  strcpy( program, argv[0] );
  for ( i=1; i<argc; i++ ) {
    if ( argv[i][0] == '-' ) {
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
	  fprintf( stderr, "unknown distribution = %s\n", argv[i] );
	}
      } 


      else if ( (strcmp ( argv[i], "-points" ) == 0) ||
		(strcmp ( argv[i], "-n" ) == 0)) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -points...\n", 0 );
	status = sscanf( argv[i],"%d",&_nbpoints_ );
	if ( status <= 0 ) _ErrorParse( "parsing -points...\n", 0 );
      }


      else if ( (strcmp ( argv[i], "-init" ) == 0) ||
		(strcmp ( argv[i], "-seed" ) == 0) ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -seed...\n", 0 );
	status = sscanf( argv[i],"%ld",&seedRandom );
	if ( status <= 0 ) _ErrorParse( "parsing -seed...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-p" ) == 0 ) {
	i += 1;
	if ( i >= argc)    _ErrorParse( "parsing -p...\n", 0 );
	status = sscanf( argv[i],"%d",&_psommet_ );
	if ( status <= 0 ) _ErrorParse( "parsing -p...\n", 0 );
      }
      
      /*--- option inconnue ---*/
      else {
	sprintf( filename, "unknown option %s\n", argv[i] );
	_ErrorParse( filename, 0);
      }
      
    } else if ( argv[i][0] != 0 ) {
      strcpy( filename, argv[i] );  
    }
  }

  _SetRandomSeed( seedRandom );

  if ( _dim_ == 2 && filename[0] != '\0' ) {
    unsigned char *buf = (unsigned char *)NULL;
    int d = 513;
    int m = 400;
    int c = 256;
    int x, y;
    
    int div;

    double ** listOfPoints;

    int fd;

    listOfPoints = _PickPoints( _nbpoints_, _max_diameter_, _min_diameter_,
				_psommet_, typeDistribution, _dim_ );

    div = _DivideListOfPoints( listOfPoints, _nbpoints_, _dim_ );
    
    buf = (unsigned char *)malloc( d*d * sizeof(unsigned char) );
    for ( i=0; i<d*d; i++ ) buf[i] = 0;

    for ( i=0; i<=div; i++ ) {
      x = (int)( listOfPoints[i][0] * m + c + 0.5 );
      y = (int)( listOfPoints[i][1] * m + c + 0.5 );
      if ( buf[x+y*d] == 0 ) buf[x+y*d] = 100;
      else                   buf[x+y*d]++;
    }
    
    for ( i=div+1; i<_nbpoints_; i++ ) {
      x = (int)( listOfPoints[i][0] * m + c + 0.5 );
      y = (int)( listOfPoints[i][1] * m + c + 0.5 );
      if ( buf[x+y*d] > 0 && buf[x+y*d] < 200 ) buf[x+y*d] = 255;
      if ( buf[x+y*d] == 255 ) continue;
      if ( buf[x+y*d] == 0 ) buf[x+y*d] = 200;
      else                   buf[x+y*d]++;
    }

    fd = creat( filename, 0644 );
    write( fd, buf, d*d * sizeof(unsigned char) );
    close( fd );
    
    free( buf );

  }
  
}
