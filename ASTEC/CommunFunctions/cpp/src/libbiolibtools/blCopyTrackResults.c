/*************************************************************************
 * blCopyTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Jeu 13 fév 2014 11:06:26 CET
 *
 * ADDITIONS, CHANGES
 *
 */

#include <stdlib.h>
#include <stdio.h>

#include <sys/time.h> /* gettimeofday() */
#include <time.h> /* clock() */
#include <string.h>

static int _verbose_ = 1;
static int _debug_ = 0;

#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thetrack_name;
  char *restrack_name;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );




static char *program = NULL;

static char *usage = "%s %s \n\
 [-v] [-help]";

static char *detail = "\
 -v : mode verbose\n\
\n";











int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_blTrackList list;


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
  

  
  BAL_InitBlTrackList( &list );
  if ( p.thetrack_name == (char*)NULL ) {
    _ErrorParse( "no detection file", 0 );
  }

  if ( BAL_ReadBlTrackList( &list, p.thetrack_name ) != 0 ) {
    _ErrorParse( "unable to read track file", 0 );
  }

  
  if ( BAL_WriteBlTrackList( &list, p.restrack_name ) != 0 ) {
    BAL_FreeBlTrackList( &list );
    _ErrorParse( "unable to write track file", 0 );
  }


  BAL_FreeBlTrackList( &list );
  return( 1 );
}













static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i, status;
  int inputisread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* transformation related file names
       */

      if ( (strcmp ( argv[i], "-verbose") == 0 && argv[i][8] == '\0')
		|| (strcmp ( argv[i], "-v") == 0 && argv[i][2] == '\0') ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
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

      /* unknown option
       */
      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( inputisread == 0 ) {
	p->thetrack_name = argv[i];
	inputisread = 1;
      }
      else if ( outputisread == 0 ) {
	p->restrack_name = argv[i];
	outputisread = 1;
      }
      else 
	fprintf(stderr,"too many file names: '%s'\n",argv[i]);
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
  p->thetrack_name = (char*)NULL;
  p->restrack_name = (char*)NULL;
}
