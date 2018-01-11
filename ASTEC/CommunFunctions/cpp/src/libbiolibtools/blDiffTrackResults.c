/*************************************************************************
 * blDiffTrackResults.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2014, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mer 19 mar 2014 22:09:36 CET
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
static int _debug_ = 1;

#include <bal-biolib-tools.h>






typedef struct local_parameter {

  /* file names
   */
  char *thetrack1_name;
  char *thetrack2_name;
  char *restrack_name;

} local_parameter;





/*------- Definition des fonctions statiques ----------*/
static void _ErrorParse( char *str, int flag );
static void _Parse( int argc, char *argv[], local_parameter *p );
static void _InitParam( local_parameter *par );




static char *program = NULL;

static char *usage = "%s %s -res %s\n\
 [-v] [-help]";

static char *detail = "\
 -v : mode verbose\n\
\n";











int main( int argc, char *argv[] )
{
  local_parameter p;
  bal_blTrackList thelist1, thelist2;
  int *equivlabel1, *equivlabel2;
  int i, j, equiv, n;

  bal_blDetectionList *declist1, *declist2;
  double epsilon = 0.001;

  int nequiv1, nnoequiv1;
  int nequiv2, nnoequiv2;

  bal_blTrackList reslist;
  char resname[1024];



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

  

  /***************************************************
   *
   * reading files
   *
   ***************************************************/

  BAL_InitBlTrackList( &thelist1 );
  BAL_InitBlTrackList( &thelist2 );


  if ( BAL_ReadBlTrackList( &thelist1, p.thetrack1_name ) != 0 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read track file '%s'\n", argv[0], p.thetrack1_name  );
    return( -1 );
  }

  if ( BAL_ReadBlTrackList( &thelist2, p.thetrack2_name ) != 0 ) {
    BAL_FreeBlTrackList( &thelist1 );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to read track file '%s'\n", argv[0], p.thetrack2_name  );
    return( -1 );
  }

  equivlabel1 = (int*)malloc( (thelist1.n + thelist2.n) * sizeof(int) );
  if ( equivlabel1 == (int*)NULL ) {
    BAL_FreeBlTrackList( &thelist2 );
    BAL_FreeBlTrackList( &thelist1 );
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to allocate label array\n", argv[0]  );
    return( -1 );
  }
  equivlabel2 = equivlabel1;
  equivlabel2 += thelist1.n;

  for ( i=0; i<thelist1.n; i++ ) equivlabel1[i] = -1;
  for ( i=0; i<thelist2.n; i++ ) equivlabel2[i] = -1;


  /***************************************************
   *
   * comparing tracks
   *
   ***************************************************/
  for ( i=0; i<thelist1.n; i++ ) {

    declist1 = &( thelist1.data[i].detectionList );

    for ( equiv=0, j=0; j<thelist2.n && equiv==0; j++){

      if ( equivlabel2[j] >= 0 ) continue;
      
      declist2 = &( thelist2.data[j].detectionList );
      
      if ( declist1->n != declist2->n ) continue;

      for ( n=0, equiv=1; n<declist1->n && equiv==1; n++ ) {

	if (  declist1->data[n].voxelcenter.x - epsilon > declist2->data[n].voxelcenter.x 
	      || declist1->data[n].voxelcenter.x + epsilon < declist2->data[n].voxelcenter.x 
	      || declist1->data[n].voxelcenter.y - epsilon > declist2->data[n].voxelcenter.y 
	      || declist1->data[n].voxelcenter.y + epsilon < declist2->data[n].voxelcenter.y 
	      || declist1->data[n].voxelcenter.z - epsilon > declist2->data[n].voxelcenter.z 
	      || declist1->data[n].voxelcenter.z + epsilon < declist2->data[n].voxelcenter.z ) 
	  equiv = 0;
      }

      if ( equiv == 0 ) continue;
      if ( equiv == 1 ) {
	equivlabel1[i] = j;
	equivlabel2[j] = i;
	if ( _debug_ ) {
	  fprintf( stderr, " track #%4d of '%s' is equal to track #%4d of '%s'\n",
		   i,  p.thetrack1_name, j,  p.thetrack2_name );
	}
      }
      
    }

  }



  if ( _debug_ ) {
    fprintf( stderr, "\n" );
    for ( i=0; i<thelist1.n; i++ ) {
      if ( equivlabel1[i] == -1 ) 
	fprintf( stderr, " track #%4d of '%s' has no equivalent\n",
		 i,  p.thetrack1_name );
    }
    fprintf( stderr, "\n" );
    for ( i=0; i<thelist2.n; i++ ) {
      if ( equivlabel2[i] == -1 ) 
	fprintf( stderr, " track #%4d of '%s' has no equivalent\n",
		 i,  p.thetrack2_name );
    }
    fprintf( stderr, "\n" );
  }




  for ( nequiv1=0, nnoequiv1=0, i=0; i<thelist1.n; i++ ) {
    if ( equivlabel1[i] >= 0 ) nequiv1 ++;
    else nnoequiv1 ++;
  }

  for ( nequiv2=0, nnoequiv2=0, i=0; i<thelist2.n; i++ ) {
    if ( equivlabel2[i] >= 0 ) nequiv2 ++;
    else nnoequiv2 ++;
  }

  if ( _verbose_ ) {
    fprintf( stderr, "Among %d tracks in '%s':\n", thelist1.n, p.thetrack1_name );
    fprintf( stderr, " - %4d tracks are also in '%s'\n", nequiv1, p.thetrack2_name );
    fprintf( stderr, " - %4d tracks are not  in '%s'\n", nnoequiv1, p.thetrack2_name );

    fprintf( stderr, "Among %d tracks in '%s':\n", thelist2.n, p.thetrack2_name );
    fprintf( stderr, " - %4d tracks are also in '%s'\n", nequiv2, p.thetrack1_name );
    fprintf( stderr, " - %4d tracks are not  in '%s'\n", nnoequiv2, p.thetrack1_name );
  }


 /***************************************************
   *
   * output
   *
   ***************************************************/

  if ( p.restrack_name != (char*) NULL ) {

    BAL_InitBlTrackList( &reslist );

    if ( BAL_AllocBlTrackList( &reslist, nequiv1 ) != 0 ) {
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate result track (common part)\n", argv[0] );
      return( -1 );
    }

    for ( reslist.n=0, i=0; i<thelist1.n; i++ ) {
      if ( equivlabel1[i] >= 0 ) {
	if ( BAL_CopyBlTrack( &(thelist1.data[i]),  &(reslist.data[reslist.n]) ) != 1 ) {
	  BAL_FreeBlTrackList( &reslist );
	  free( equivlabel1 );
	  BAL_FreeBlTrackList( &thelist2 );
	  BAL_FreeBlTrackList( &thelist1 );
	  if ( _verbose_ )
	     fprintf( stderr, "%s: unable to copy track #%d (common part)\n", argv[0], i );
	  return( -1 );
	}
	reslist.n ++;
      }
    }
    
    sprintf( resname, "%s-common.txt", p.restrack_name );
    if ( BAL_WriteBlTrackList( &reslist, resname ) != 0 ) {
      BAL_FreeBlTrackList( &reslist );
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write track file '%s'\n", argv[0], resname );
      return( -1 );
    }
    
    BAL_FreeBlTrackList( &reslist );    

    if ( BAL_AllocBlTrackList( &reslist, nnoequiv1 ) != 0 ) {
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate result track (part #1)\n", argv[0] );
      return( -1 );
    }

    for ( reslist.n=0, i=0; i<thelist1.n; i++ ) {
      if ( equivlabel1[i] < 0 ) {
	if ( BAL_CopyBlTrack( &(thelist1.data[i]),  &(reslist.data[reslist.n]) ) != 1 ) {
	  BAL_FreeBlTrackList( &reslist );
	  free( equivlabel1 );
	  BAL_FreeBlTrackList( &thelist2 );
	  BAL_FreeBlTrackList( &thelist1 );
	  if ( _verbose_ )
	     fprintf( stderr, "%s: unable to copy track #%d (part #1)\n", argv[0], i );
	  return( -1 );
	}
	reslist.n ++;
      }
    }
    
    sprintf( resname, "%s-part1.txt", p.restrack_name );
    if ( BAL_WriteBlTrackList( &reslist, resname ) != 0 ) {
      BAL_FreeBlTrackList( &reslist );
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write track file '%s'\n", argv[0], resname );
      return( -1 );
    }
    
    BAL_FreeBlTrackList( &reslist );    

    if ( BAL_AllocBlTrackList( &reslist, nnoequiv2 ) != 0 ) {
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to allocate result track (part #2)\n", argv[0] );
      return( -1 );
    }

    for ( reslist.n=0, i=0; i<thelist2.n; i++ ) {
      if ( equivlabel2[i] < 0 ) {
	if ( BAL_CopyBlTrack( &(thelist2.data[i]),  &(reslist.data[reslist.n]) ) != 1 ) {
	  BAL_FreeBlTrackList( &reslist );
	  free( equivlabel1 );
	  BAL_FreeBlTrackList( &thelist2 );
	  BAL_FreeBlTrackList( &thelist1 );
	  if ( _verbose_ )
	     fprintf( stderr, "%s: unable to copy track #%d (part #2)\n", argv[0], i );
	  return( -1 );
	}
	reslist.n ++;
      }
    }
    
    sprintf( resname, "%s-part2.txt", p.restrack_name );
    if ( BAL_WriteBlTrackList( &reslist, resname ) != 0 ) {
      BAL_FreeBlTrackList( &reslist );
      free( equivlabel1 );
      BAL_FreeBlTrackList( &thelist2 );
      BAL_FreeBlTrackList( &thelist1 );
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to write track file '%s'\n", argv[0], resname );
      return( -1 );
    }
    
    BAL_FreeBlTrackList( &reslist );    
  }






  free( equivlabel1 );
  BAL_FreeBlTrackList( &thelist2 );
  BAL_FreeBlTrackList( &thelist1 );



  return( 1 );
}














static void _Parse( int argc, char *argv[], local_parameter *p )
{
  int i;
  int input1isread = 0;
  int input2isread = 0;
  int outputisread = 0;

  program = argv[0];
	
  for ( i=1; i<argc; i++ ) {
  
    if ( argv[i][0] == '-' ) {

      /* general options
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

       else if ( strcmp ( argv[i], "-result") == 0
	   || (strcmp ( argv[i], "-res") == 0 && argv[i][4] == '\0') ) {
	i++;
	if ( i >= argc) _ErrorParse( "parsing -result", 0 );
	p->restrack_name = argv[i];
	outputisread = 1;
      }
      
      /* unknown option
       */
      else {
	fprintf(stderr,"unknown option: '%s'\n",argv[i]);
      }
    }
    
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( input1isread == 0 ) {
	p->thetrack1_name = argv[i];
	input1isread = 1;
      }
      else if ( input2isread == 0 ) {
	p->thetrack2_name = argv[i];
	input2isread = 1;
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
  p->thetrack1_name = (char*)NULL;
  p->thetrack2_name = (char*)NULL;
  p->restrack_name = (char*)NULL;

}
