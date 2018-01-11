/*************************************************************************
 * test-max-diameter.c
 *
 * $Id: test-maxdiameter.c,v 1.1 1999/07/30 15:25:04 greg Exp $
 *
 * DESCRIPTION: 
 *
 *
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Jul 30 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 */

#include <vt_common.h>
#include <vt_elfparam.h>


typedef struct local_par {
  int testall;
  char filename[STRINGLENGTH];
  int nbpoints;
  int nbtests;
} local_par;


/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[nom-de-fichier] [-tests %d] [-points %d]";

static char *detail = "\
\t -v : mode verbose\n\
\t -D : mode debug\n";

static char program[STRINGLENGTH];

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  

  _TestsMaxDiameterInRandomList( par.nbpoints, par.nbtests, par.filename, par.testall );

}






#if defined(_ANSI_)
static void VT_Parse( int argc, char *argv[], local_par *par )
#else
  static void VT_Parse( argc, argv, par )
  int argc;
char *argv[];
local_par *par;
#endif
{
  int i, nb, status;
  int nb_estimates;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	VT_ErrorParse("do not deal with standard input", 0 );
      }
      /*--- arguments generaux ---*/
      else if ( strcmp ( argv[i], "-help" ) == 0 ) {
	VT_ErrorParse("\n", 1);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 ) {
	_VT_VERBOSE_ = 1;
      }
      else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	_VT_DEBUG_ = 1;
      }
      /*--- traitement eventuel de l'image d'entree ---*/

      /*---  ---*/
      else if ( strcmp ( argv[i], "-all" ) == 0 ) {
	par->testall = 1;
      }
      else if ( strcmp ( argv[i], "-tests" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -tests...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nbtests) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -tests...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-points" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -points...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nbpoints) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -points...\n", 0 );
      }


      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) { 
	strncpy( par->filename, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->filename,  ">" );  /* standart output */
  }


}

#if defined(_ANSI_)
static void VT_ErrorParse( char *str, int flag )
#else
  static void VT_ErrorParse( str, flag )
  char *str;
int flag;
#endif
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}

#if defined(_ANSI_)
static void VT_InitParam( local_par *par )
#else
  static void VT_InitParam( par )
  local_par *par;
#endif
{
  par->testall = 0;
  par->filename[0] = '\0';
  par->nbpoints = 1000;
  par->nbtests = 10;
}

