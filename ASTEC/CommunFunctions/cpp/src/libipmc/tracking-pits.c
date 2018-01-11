/*************************************************************************
 * tracking-pits.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 19 jui 2013 22:21:01 CEST
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>

#include <tracking-tools.h>

#include <vt_names.h>
#include <vt_error.h>


static int _verbose_ = 1;
static int _debug_ = 0;


typedef struct local_par {
  
  char *theformat;
  char *resname;

  char *list;
  char *hist;

  char *description;

  int first;
  int last;

  int margin;
  int depth;

  int maxForwardNeighbors;
  int maxBackwardNeighbors;
  
  int removeBorderChains;

  char *matlabname;
  char *scilabname;
  char *textname;

  int print_time;

} local_par;


static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );

static char *usage = "[format-in] [-first %d] [-last %d]\n\
 [-margin %d] [-depth %d]\n\
 [-max-forward-neighbors | -maxfn %d] [-max-backward-neighbors | -maxbn %d]\n\
 [-removeBorderChains]\n\
 [-matlab %s] [-scilab %s] [-text %s] [-description|-desc %s]\n\
 [-list %s] [-histogram | -hist %s]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-time] [-notime]\n";

static char *detail = "\
 format-in    # format 'a la printf' of images to be processed, must contain a '%d'\n\
 -first %d    # first value of the index in the format\n\
 -last %d     # last value of the index in the format\n\
 -margin %d   # two circles are in the same chain if their center distance is less\n\
                than the sum of the radii + the margin\n\
 -depth  %d   # maximal increment for the next images\n\
                '2' means that one image can be skipped if no neighbors are found\n\
 -max-forward-neighbors  %d # selection of point chains\n\
                suppress pits with too many 'forward' neighbors\n\
                '1' keeps pits that have only one 'forward' neighbor\n\
                i.e. makes sure that there is no splitting\n\
 -max-backward-neighbors %d # selection of point chains\n\
                suppress pits with too many 'backward' neighbors\n\
                '1' keeps pits that have only one 'forward' neighbor\n\
                i.e. makes sure that there is no merging\n\
 -removeBorderChains #\n\
";




static char program[1024];

int main( int argc, char *argv[] )
{
  local_par par;
  circleListList readlist;
  chainList chainlist;
  statisticType stats;


  int nlist;
  char ch1desc[1024];
  FILE *f;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );


  if ( par.theformat == (char*)NULL ) {
    VT_ErrorParse( "no input format\n", 0 );
  }
  if ( par.last - par.first + 1 <= 0 ) {
    VT_ErrorParse( "last < first\n", 0 );
  }


  if ( par.description == (char*)NULL ) {
    sprintf( ch1desc, "ch1" );
  }
  else {
    sprintf( ch1desc, "ch1_%s", par.description );
  }



  /* allocation
   */
  nlist = par.last - par.first + 1;
  initCircleListList( &readlist );
  if ( allocCircleListList( &readlist, nlist ) != 1 ) {
    VT_ErrorParse( "allocation error\n", 0 );
  }


  /* lecture des points et construction des chaines
   */
  initChainList( &chainlist );

  if ( readCircleListList( &readlist, par.theformat, par.first, par.last ) != 1 ) {
    freeCircleListList( &readlist );
    VT_ErrorParse( "error when reading data set\n", 0 );
  }

  if ( chainListFromCircleListList( &readlist, &chainlist, 
				    par.depth, par.margin,
				    par.maxForwardNeighbors, par.maxBackwardNeighbors ) != 1 ) {
    freeChainList( &chainlist );
    freeCircleListList( &readlist );
    VT_ErrorParse( "error when building chain from data set\n", 0 );
  }


  /* selection des chaines
   */
  if ( par.removeBorderChains ) { 
    removeStartBoundChain( &chainlist, 0 );
    removeEndBoundChain( &chainlist, par.last-par.first );
  }


  /* outputs
   */

  if ( par.textname != NULL ) {
    f = fopen( par.textname, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.textname );
    }
    else { 
      printStatsChainList( f, &readlist, &chainlist, ch1desc );
      fclose( f );
    }
  }

  printStatsChainList( stderr, &readlist, &chainlist, ch1desc );

  initStatisticType( &stats );
  statsChainList( &stats,  &readlist, &chainlist );


  if ( par.scilabname != NULL ) {
    printChainResultXxxlab( &stats, ch1desc, par.scilabname, _SCILAB_ );
  }



  if ( par.hist != (char*)NULL ) {
    f = fopen( par.hist, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.hist );
    }
    else {
      printPlotHist( f, &(stats.length) );
      fclose( f );
    }
  }
  
  if ( par.list != (char*)NULL ) {
    f = fopen( par.list, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.list );
    }
    else {
      printPlotList( f, &(stats.length) );
      fclose( f );
    }
  }



  freeStatisticType( &stats );
  freeChainList( &chainlist );
  freeCircleListList( &readlist );




  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }

  return( 1 );
}

















static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, nb;
  int status;
  int maxchunks;
  char text[256];

  if ( VT_CopyName( _BaseName( program ), argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/

  for ( i=1, nb=0; i<argc; i++ ) {
    
    if ( argv[i][0] == '-' ) {

      /*--- arguments generaux ---*/
      if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) {
	VT_ErrorParse("\n", 0);
      }

      else if ( (strcmp ( argv[i], "-D" ) == 0 && argv[i][2] == '\0') 
		|| strcmp ( argv[i], "-debug" ) == 0 ) {
	_debug_ = 1;
      }


     else if ( strcmp ( argv[i], "-first" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -first...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->first) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -first...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-last" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -last...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->last) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -last...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-margin" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -margin...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->margin) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -margin...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-depth" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -depth...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->depth) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -depth...\n", 0 );
      }


      else if ( strcmp ( argv[i], "-histogram" ) == 0 
		|| strcmp ( argv[i], "-hist" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -histogram...\n", 0 );
	par->hist = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-list" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -list...\n", 0 );
	par->list = argv[i];
      }
  

      else if ( strcmp ( argv[i], "-max-forward-neighbors" ) == 0 
		|| (strcmp ( argv[i], "-maxfn" ) == 0 && argv[i][6] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -max-forward-neighbors...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->maxForwardNeighbors) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -max-forward-neighbors...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-max-backward-neighbors" ) == 0
		|| (strcmp ( argv[i], "-maxbn" ) == 0 && argv[i][6] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -max-backward-neighbors...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->maxBackwardNeighbors) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -max-backward-neighbors...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-removeBorderChains" ) == 0 ) {
	par->removeBorderChains = 1;
      }


      else if ( strcmp ( argv[i], "-description" ) == 0 || strcmp ( argv[i], "-desc" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -description...\n", 0 );
	par->description = argv[i];
      }


      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	par->matlabname = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-scilab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -scilab...\n", 0 );
	par->scilabname = argv[i];
      }
      else if ( strcmp ( argv[i], "-text" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -text...\n", 0 );
	par->textname = argv[i];
      }



      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
	setMaxChunks( 100 );
      }
      
      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
      setMaxChunks( 1 );
      }

      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-max-chunks", 0 );
	status = sscanf( argv[i], "%d", &maxchunks );
	if ( status <= 0 ) VT_ErrorParse( "-max-chunks", 0 );
	if ( maxchunks >= 1 ) setMaxChunks( maxchunks );
      }

      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	par->print_time = 0;
      }


      else {
	sprintf( text,"unknown option %s\n",argv[i] );
	VT_ErrorParse(text, 0);
      }

    }
    /*--- saisie des noms d'images ---*/
    else {
      if ( par->theformat == NULL ) {
        par->theformat = argv[i];
      }
      else if ( par->resname == NULL ) {
        par->resname = argv[i];
      }
      else {
        VT_ErrorParse("too much file names when parsing\n", 0 );
      }
    }
  }

}



static void VT_ErrorParse( char *str, int flag )
{
	(void)fprintf(stderr,"Usage : %s %s\n",_BaseName( program ), usage);
        if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
        (void)fprintf(stderr,"Erreur : %s",str);
        exit(0);
}



static void VT_InitParam( local_par *par )
{
  par->theformat = (char*)NULL;
  par->resname = (char*)NULL;

  par->list = (char*)NULL;
  par->hist = (char*)NULL;

  par->description = (char*)NULL;

  par->first = 1;
  par->last = 1;

  par->depth = 1;
  par->margin = 0;

  par->maxForwardNeighbors = -1;
  par->maxBackwardNeighbors = -1;

  par->removeBorderChains = 0;

  par->matlabname = (char*)NULL;
  par->scilabname = (char*)NULL;
  par->textname = (char*)NULL;

  par->print_time = 1;

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



static char *_BaseName( char *p )
{
  int l;
  if ( p == (char*)NULL ) return( (char*)NULL );
  l = strlen( p ) - 1;
  while ( l >= 0 && p[l] != '/' ) l--;
  if ( l < 0 ) l = 0;
  if ( p[l] == '/' ) l++;
  return( &(p[l]) );
}
