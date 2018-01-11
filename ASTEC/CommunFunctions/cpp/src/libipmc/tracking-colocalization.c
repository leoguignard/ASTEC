/*************************************************************************
 * tracking-colocalization.c -
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
  
  char *theformat1;
  char *theformat2;
  char *resname;

  char *colocalizationlist;
  char *colocalizationhist;

  char *description;

  int first;
  int last;

  int margin;
  int depth;

  int colocalization_margin;
  int colocalization_depth;

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



static char *usage = "[format-in-1] [format-in-2] [-first %d] [-last %d]\n\
 [-margin %d] [-depth %d]\n\
 [-colocalization-margin | -cmargin %d] [-colocalization-depth | -cdepth %d]\n\
 [-max-forward-neighbors | -maxfn %d] [-max-backward-neighbors | -maxbn %d]\n\
 [-removeBorderChains]\n\
 [-matlab %s] [-scilab %s] [-text %s] [-description|-desc %s]\n\
 [-colocalization-list | -clist %s] [-colocalization-histogram | -chist %s]\n\
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

  circleListList pitlist1;
  chainList pitchain1;
  statisticType stats1;

  circleListList pitlist2;
  chainList pitchain2;
  statisticType stats2;

  circleListList coloclist;
  chainList colocchain;
  statisticType colocstats;

  int nlist;
  char ch1desc[1024];
  char ch2desc[1024];
  char coldesc[1024];
  FILE *f;

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;



  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );


  if ( par.theformat1 == (char*)NULL || par.theformat2 == (char*)NULL) {
    VT_ErrorParse( "no input formats\n", 0 );
  }
  if ( par.last - par.first + 1 <= 0 ) {
    VT_ErrorParse( "last < first\n", 0 );
  }


  if ( par.description == (char*)NULL ) {
    sprintf( ch1desc, "ch1" );
    sprintf( ch2desc, "ch2" );
    sprintf( coldesc, "col" );
  }
  else {
    sprintf( ch1desc, "ch1_%s", par.description );
    sprintf( ch2desc, "ch2_%s", par.description );
    sprintf( coldesc, "col_%s", par.description );
  }



  /* allocation
   */
  nlist = par.last - par.first + 1;

  initCircleListList( &pitlist1 );
  if ( allocCircleListList( &pitlist1, nlist ) != 1 ) {
    VT_ErrorParse( "allocation error\n", 0 );
  }

  initCircleListList( &pitlist2 );
  if ( allocCircleListList( &pitlist2, nlist ) != 1 ) {
    VT_ErrorParse( "allocation error\n", 0 );
  }


  /* lecture des points et construction des chaines
   */
  fprintf( stderr, "... reading first data set\n" );

  initChainList( &pitchain1 );
  
  if ( readCircleListList( &pitlist1, par.theformat1, par.first, par.last ) != 1 ) {
    freeCircleListList( &pitlist1 );
    VT_ErrorParse( "error when reading first data set\n", 0 );
  }

  if ( chainListFromCircleListList( &pitlist1, &pitchain1, 
				    par.depth, par.margin,
				    par.maxForwardNeighbors, par.maxBackwardNeighbors ) != 1 ) {
    freeChainList( &pitchain1 );
    freeCircleListList( &pitlist1 );
    VT_ErrorParse( "error when building chain from first data set\n", 0 );
  }



  fprintf( stderr, "... reading second data set\n" );

  initChainList( &pitchain2 );

  if ( readCircleListList( &pitlist2, par.theformat2, par.first, par.last ) != 1 ) {
    freeChainList( &pitchain1 );
    freeCircleListList( &pitlist1 );
    freeCircleListList( &pitlist2 );
    VT_ErrorParse( "error when reading second data set\n", 0 );
  }

  if ( chainListFromCircleListList( &pitlist2, &pitchain2, 
				    par.depth, par.margin,
				    par.maxForwardNeighbors, par.maxBackwardNeighbors ) != 1 ) {
    freeChainList( &pitchain1 );
    freeCircleListList( &pitlist1 );
    freeChainList( &pitchain2 );
    freeCircleListList( &pitlist2 );
    VT_ErrorParse( "error when building chain from second data set\n", 0 );
  }



  /* selection des chaines
   */
  if ( par.removeBorderChains ) { 
    removeStartBoundChain( &pitchain1, 0 );
    removeEndBoundChain( &pitchain1, par.last-par.first );

    removeStartBoundChain( &pitchain2, 0 );
    removeEndBoundChain( &pitchain2, par.last-par.first );
  }





  /* colocalisation
   */
  fprintf( stderr, "... building colocalizations\n" );

  initCircleListList( &coloclist );
  initChainList( &colocchain );
  
  if ( colocalizationFromLists( &coloclist, &colocchain,
				&pitlist1, &pitlist2, 
				par.colocalization_depth, par.colocalization_margin,
				par.maxForwardNeighbors, par.maxBackwardNeighbors ) != 1 ) {
    freeChainList( &colocchain );
    freeCircleListList( &coloclist );
    freeChainList( &pitchain1 );
    freeCircleListList( &pitlist1 );
    freeChainList( &pitchain2 );
    freeCircleListList( &pitlist2 );
    VT_ErrorParse( "colocalization error\n", 0 );
  }

  
  /* selection des chaines
   */
  if ( par.removeBorderChains ) { 
    removeStartBoundChain( &colocchain, 0 );
    removeEndBoundChain( &colocchain, par.last-par.first );
  }



  /* outputs
   */

  if ( par.textname != (char*)NULL ) {
    f = fopen( par.textname, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.textname );
    }
    else { 
      printStatsChainList( f, &pitlist1, &pitchain1, ch1desc );
      printStatsChainList( f, &pitlist2, &pitchain2, ch2desc );
      printStatsChainList( f, &coloclist, &colocchain, coldesc );
      fclose( f );
    }
  }

  printStatsChainList( stderr, &pitlist1, &pitchain1, ch1desc );
  printStatsChainList( stderr, &pitlist2, &pitchain2, ch2desc );
  printStatsChainList( stderr, &coloclist, &colocchain, coldesc );
  





  initStatisticType( &stats1 );
  statsChainList( &stats1, &pitlist1, &pitchain1 );

  initStatisticType( &stats2 );
  statsChainList( &stats2, &pitlist2, &pitchain2 );

  initStatisticType( &colocstats );
  statsChainList( &colocstats, &coloclist, &colocchain );



  if ( par.colocalizationhist != (char*)NULL ) {
    f = fopen( par.colocalizationhist, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.colocalizationhist );
    }
    else {
      printPlotHist( f, &(colocstats.length) );
      fclose( f );
    }
  }
  
  if ( par.colocalizationlist != (char*)NULL ) {
    f = fopen( par.colocalizationlist, "w" );
    if ( f == (FILE*)NULL ) {
      if ( _verbose_ )
	fprintf( stderr, " ...unable to open '%s' for writing\n", par.colocalizationlist );
    }
    else {
      printPlotList( f, &(colocstats.length) );
      fclose( f );
    }
  }



  if ( par.scilabname != NULL ) {
    printColocalizationResultXxxlab( &stats1, ch1desc, 
				     &stats2, ch2desc,  
				     &colocstats, coldesc, 
				     par.scilabname, _SCILAB_ );
  }

  freeChainList( &colocchain );
  freeCircleListList( &coloclist );

  freeStatisticType( &stats1 );
  freeChainList( &pitchain1 );
  freeCircleListList( &pitlist1 );

  freeStatisticType( &stats2 );
  freeChainList( &pitchain2 );
  freeCircleListList( &pitlist2 );



  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName(program), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }
  
  return( 0 );
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
      
      else if ( strcmp ( argv[i], "-colocalization-histogram" ) == 0 
		|| strcmp ( argv[i], "-chist" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -colocalization-histogram...\n", 0 );
	par->colocalizationhist = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-colocalization-list" ) == 0 
		|| strcmp ( argv[i], "-clist" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -colocalization-list...\n", 0 );
	par->colocalizationlist = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-colocalization-margin" ) == 0 
		|| strcmp ( argv[i], "-cmargin" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -colocalization-margin...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->colocalization_margin) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -colocalization-margin...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-colocalization-depth" ) == 0 
		|| strcmp ( argv[i], "-cdepth" ) == 0 ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -colocalization-depth...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->colocalization_depth) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -colocalization-depth...\n", 0 );
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
      if ( par->theformat1 == NULL ) {
        par->theformat1 = argv[i];
      }
      else if ( par->theformat2 == NULL ) {
        par->theformat2 = argv[i];
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
  par->theformat1 = (char*)NULL;
  par->theformat2 = (char*)NULL;
  par->resname = (char*)NULL;

  par->colocalizationlist = (char*)NULL;
  par->colocalizationhist = (char*)NULL;

  par->description = (char*)NULL;

  par->first = 1;
  par->last = 1;

  par->depth = 1;
  par->margin = 0;

  par->colocalization_depth = 1;
  par->colocalization_margin = 0;

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
