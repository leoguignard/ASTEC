/*************************************************************************
 * tracking-profiles.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2013, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 25 sep 2013 18:00:00 CEST
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


/*
static int _verbose_ = 1;
*/
static int _debug_ = 0;


typedef struct local_par {
  
  char *thename;
  char *resname;
  char *theimageformat;

  int imageindexshift;

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

static char *usage = "%s %s [-image %s]\n\
 [-image-index-shift|-iis|-shift %d]\n\
 [-matlab %s] [-scilab %s] [-text %s] [-description|-desc %s]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-time] [-notime]\n";

static char *detail = "\
";




static char program[1024];

int main( int argc, char *argv[] )
{
  local_par par;
  circleListList readlist, filledlist;


  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );


  /* lecture de la liste
   */
  initCircleListList( &readlist );
  if ( readTrackedCircleList( &readlist, par.thename ) != 1 ) {
    freeCircleListList( &readlist );
    VT_ErrorParse( "error when reading data set\n", 0 );
  }


  /* gap filling 
   */
  initCircleListList( &filledlist );
  if ( fillTrackedCircleList( &filledlist, &readlist ) != 1 ) {
    freeCircleListList( &filledlist );
    freeCircleListList( &readlist );
    VT_ErrorParse( "error when filling gaps\n", 0 );
  }


  freeCircleListList( &readlist );


  /* had shift to image index
   */
  if ( par.imageindexshift != 0 ) {
    int i, j;
    for ( j=0; j<filledlist.n_data; j++ )
    for ( i=0; i<filledlist.data[j].n_data; i++ )
      filledlist.data[j].data[i].image += par.imageindexshift;
  }


  /* */
  if ( intensityTrackedCircleList( &filledlist, par.theimageformat ) != 1 ) {
    freeCircleListList( &filledlist );
    VT_ErrorParse( "error when computing intensity characteristics\n", 0 );
  }




  if ( par.resname != NULL ) {
    FILE *f = fopen( par.resname, "w" );
    if ( f == (FILE*)NULL ) {
      freeCircleListList( &filledlist );
      VT_ErrorParse( "error when opening output file\n", 0 );
    }
    printTrackedCircleList( f, &filledlist );
    fclose( f );
  }



  

  freeCircleListList( &filledlist );
  


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



      else if ( strcmp ( argv[i], "-image-index-shift" ) == 0 
		|| (strcmp ( argv[i], "-shift" ) == 0 && argv[i][6] == '\0') 
		|| (strcmp ( argv[i], "-iis" ) == 0 && argv[i][4] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image-index-shift...\n", 0 );
        status = sscanf( argv[i],"%d",&(par->imageindexshift) );
        if ( status <= 0 ) VT_ErrorParse( "parsing -image-index-shift...\n", 0 );
      }

      else if ( (strcmp ( argv[i], "-image" ) == 0 && argv[i][6] == '\0')
		|| (strcmp ( argv[i], "-format" ) == 0 && argv[i][6] == '\0') ) {
        i += 1;
        if ( i >= argc)    VT_ErrorParse( "parsing -image...\n", 0 );
	par->theimageformat = argv[i];
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
      if ( par->thename == NULL ) {
        par->thename = argv[i];
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
  par->thename = (char*)NULL;
  par->resname = (char*)NULL;
  par->theimageformat = (char*)NULL;

  par->imageindexshift = 0;

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
