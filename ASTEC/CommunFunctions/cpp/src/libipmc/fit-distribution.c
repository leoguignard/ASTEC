/*************************************************************************
 * fit-distribution.c -
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

#include <histogram.h>
#include <levenberg.h>

#include <vt_names.h>
#include <vt_error.h>


#include <fit-distribution-tools.h>


static int _verbose_ = 0;


typedef struct local_par {
  
  char *inputlist;
  char *inputhist;

  double initValues[9];

  int minlength;
  int maxlength;
  int removefromhistogram;

  char *matlabname;
  char *scilabname;

} local_par;



static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static void VT_PrintParam( FILE *f, local_par *par );
static char *_BaseName( char *p );

static char *usage = "[-list %s | [-histogram|-hist] %s] [-matlab %s] [-scilab %s]\n\
 [-init-distrib1 %lf %lf [%lf]] [-init-distrib2 %lf %lf [%lf]] [-init-distrib3 %lf %lf [%lf]]\n\
 [-min-length|-min %d] [-max-length|-max %d] [-remove-from-histogram|-rfh]\n\
 [-print-parameters]";

static char *detail = "\
";




static char program[1024];

int main( int argc, char *argv[] )
{
  local_par par;
  measureList theList;
  typeHistogram theHisto;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );

  /* read data 
   */

  if ( par.inputlist != (char*)NULL ) {

    initMeasureList( &theList );
    if ( readMeasureList( &theList, par.inputlist ) != 1 ) {
      freeMeasureList( &theList );
      VT_ErrorParse("unable to read input file\n", 0 );
    }
    
    /* build histogram
     */
    initHistogram( &theHisto );
    if ( build1DHistogramFromMeasureList( &theHisto, &theList ) != 1 ) {
      freeHistogram( &theHisto );
      freeMeasureList( &theList );
      VT_ErrorParse("unable to build histogram\n", 0 );
    }

  }


  if ( par.removefromhistogram && (par.minlength >= 0 || par.maxlength >=0) ) {
    switch ( theHisto.typeHisto ) {
    default :
      freeHistogram( &theHisto );
      freeMeasureList( &theList );
      VT_ErrorParse("such histogram type not handled yet\n", 0 );
    case SINT :
      {
	s32 *theBuf = (s32*)theHisto.data;
	int i;
	if ( par.minlength >= 0 && par.minlength < theHisto.xaxis.dim ) {
	  for ( i=0; i<par.minlength && i<theHisto.xaxis.dim; i++ ) 
	    theBuf[i] = 0;
	}
	if ( par.maxlength >= 0 && par.maxlength < theHisto.xaxis.dim ) {
	  for ( i=theHisto.xaxis.dim-1; i>par.maxlength && i>=0; i-- ) 
	    theBuf[i] = 0;
	}
      }
      break;
    case FLOAT :
      {
	r32 *theBuf = (r32*)theHisto.data;
	int i;
	if ( par.minlength >= 0 && par.minlength < theHisto.xaxis.dim ) {
	  for ( i=0; i<par.minlength && i<theHisto.xaxis.dim; i++ ) 
	    theBuf[i] = 0;
	}
	if ( par.maxlength >= 0 && par.maxlength < theHisto.xaxis.dim ) {
	  for ( i=theHisto.xaxis.dim-1; i>par.maxlength && i>=0; i-- ) 
	    theBuf[i] = 0;
	}
      }
      break;
    }
  }

  
  if ( par.matlabname != NULL ) {
    if ( fit3WeibullDistributionsOnCumulative( &theHisto, par.initValues, par.matlabname, _MATLAB_, par.minlength, par.maxlength ) != 1 ) {
      freeHistogram( &theHisto );
      freeMeasureList( &theList );
      VT_ErrorParse("unable to computing matlab result\n", 0 );
    }
  }

  if ( par.scilabname != NULL ) {
    if ( fit3WeibullDistributionsOnCumulative( &theHisto, par.initValues, par.scilabname, _SCILAB_, par.minlength, par.maxlength  ) != 1 ) {
      freeHistogram( &theHisto );
      freeMeasureList( &theList );
      VT_ErrorParse("unable to computing scilab result\n", 0 );
    }
  }

  freeHistogram( &theHisto );
  freeMeasureList( &theList );

  return( 1 );
}

















static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, nb;
  int printparams = 0;
  int status;
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
      
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	if ( _verbose_ <= 0 ) _verbose_ = 1;
	else _verbose_ ++;
	setVerboseInLevenberg( _verbose_ );
      }
      
      else if ( strcmp ( argv[i], "-nv" ) == 0 && argv[i][3] == '\0' ) {
	_verbose_ = 0;
	setVerboseInLevenberg( _verbose_ );
      }
      
      else if ( strcmp ( argv[i], "-list" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -list...\n", 0 );
	par->inputlist = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-histogram" ) == 0  || strcmp ( argv[i], "-hist" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -histogram...\n", 0 );
	par->inputhist = argv[i];
      }
      
      else if ( strcmp ( argv[i], "-min-length" ) == 0 
		|| (strcmp ( argv[i], "-min" ) == 0 && argv[i][4] == '\0') ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -min-length...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->minlength) );
	 if ( status <= 0 ) VT_ErrorParse( "parsing -min-length...\n", 0 );
      }
      
      else if ( strcmp ( argv[i], "-max-length" ) == 0 
		|| (strcmp ( argv[i], "-max" ) == 0 && argv[i][4] == '\0') ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -max-length...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->maxlength) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -max-length...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-remove-from-histogram" ) == 0 
		|| (strcmp ( argv[i], "-rfh" ) == 0 && argv[i][4] == '\0') ) {
	par->removefromhistogram = 1;
      }

      else if ( strcmp ( argv[i], "-init-distrib1" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib1...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[0]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib1...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib1...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[1]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib1...\n", 0 );
	i += 1;
	if ( i < argc ) {
	  status = sscanf( argv[i], "%lf", &(par->initValues[2]) );
	  if ( status <= 0 ) 
	    i -= 1;
	}
      }
	
      else if ( strcmp ( argv[i], "-init-distrib2" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib2...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[4]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib2...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib2...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[5]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib2...\n", 0 );
	i += 1;
	if ( i < argc ) {
	  status = sscanf( argv[i], "%lf", &(par->initValues[6]) );
	  if ( status <= 0 ) 
	    i -= 1;
	}
      }
      
      else if ( strcmp ( argv[i], "-init-distrib3" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib3...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[7]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib3...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init-distrib3...\n", 0 );
	status = sscanf( argv[i], "%lf", &(par->initValues[8]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -init-distrib3...\n", 0 );
	i += 1;
	if ( i < argc ) {
	  status = sscanf( argv[i], "%lf", &(par->initValues[9]) );
	  if ( status <= 0 ) 
	    i -= 1;
	}
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

      else if ( strcmp ( argv[i], "-print-parameters" ) == 0 ) {
	printparams = 1;
      }

      else {
	sprintf( text,"unknown option %s\n",argv[i] );
	VT_ErrorParse(text, 0);
      }

    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      VT_ErrorParse("too much file names when parsing\n", 0 );
    }
  }

  if ( printparams ) {
    VT_PrintParam( stderr, par );
    exit( 0 );
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
  par->inputlist = (char*)NULL;
  par->inputhist = (char*)NULL;

  par->initValues[0] = 0.6;
  par->initValues[1] = 2;
  par->initValues[2] = 2;

  par->initValues[3] = 0.2;
  par->initValues[4] = 4;
  par->initValues[5] = 2;

  par->initValues[6] = 0.2;
  par->initValues[7] = 10;
  par->initValues[8] = 1;

  par->removefromhistogram = 0;

  par->minlength = -1;
  par->maxlength = -1;

  par->matlabname = (char*)NULL;
  par->scilabname = (char*)NULL;
}



static void VT_PrintParam( FILE *f, local_par *par )
{
  fprintf( f, "input  list name = '%s'\n", par->inputlist );
  fprintf( f, "input  hist name = '%s'\n", par->inputlist );
  
  fprintf( f, "initial values distribution #1 = [%f %f %f]\n",
	   par->initValues[0], par->initValues[1], par->initValues[2] );
  fprintf( f, "initial values distribution #2 = [%f %f %f]\n",
	   par->initValues[3], par->initValues[4], par->initValues[5] );
  fprintf( f, "initial values distribution #3 = [%f %f %f]\n",
	   par->initValues[6], par->initValues[7], par->initValues[8] );

  fprintf( f, "min length = %d\n", par->minlength );
  fprintf( f, "max length = %d\n", par->maxlength );

  fprintf( f, "remove from histogram = %d\n", par->removefromhistogram );

  fprintf( f, "scilab name = '%s'\n", par->scilabname );
  fprintf( f, "matlab name = '%s'\n", par->matlabname );
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
