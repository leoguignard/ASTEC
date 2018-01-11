/*************************************************************************
 * generHisto.c - generation d'histogramme
 *
 * $Id: generHisto.c,v 1.6 2000/02/17 10:26:22 greg Exp $
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
 * May 1999
 *
 * Copyright Gregoire Malandain, INRIA
 *
 *
 * ADDITIONS, CHANGES:
 *
 *
 *
 * HELP MATLAB
 *
 * % matlab
 *
 * % pour changer de repertoire
 * >> cd /u/freud/0/
 *
 * % pour creer une figure 
 * >> figure
 *
 * % pour utiliser une figure existante (numero n)
 * >> figure(n)
 *
 * % pour ajouter une courbe dans une figure ou il y en a deja une
 * % (c'est une bascule)
 * >> hold
 *
 * % pour mettre une figure en couleur (ici en rouge)
 * >> plot (x,y, 'r')
 *
 * % pour dessiner avec des croix et en rouge
 * > plot (x,y, 'r+')
 *
 * % pour utiliser more pour l'aide
 * >> more on
 *
 * % pour avoir de l'aide 
 * >> help plot 
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 */

#include <vt_common.h>

static int _NBALLOCS_ = 500;

typedef struct {
  double vol;
  double brd;
  double radVol;
  double radSrf;
  double max;
  double med;
  double min;

  double value;
  double cumulatedValue;

  double xvalue;
  double yvalue;
} typeEntry;

typedef struct {
  double index;
  double value;
} typeElemHisto1D;

typedef struct {
  double xindex;
  double yindex;
  double value;
} typeElemHisto2D;



typedef enum {
  _HISTO_VOLUME_,
  _HISTO_RADIUS_,
  _HISTO_MAXFERET_,
  _HISTO_MEDFERET_,
  _HISTO_MINFERET_,
  _HISTO_MINSURMAX_,
  _HISTO_MINSMAX_MINSMED_
} enumOutput;

typedef enum {
  _SIMPLE_,
  _NORMALIZE_,
  _CUMULATIVE_,
  _CUMULATIVE_NORMALIZE_
} enumHisto;


typedef enum {
  _ALL_,
  _MATLAB_,
  _EXCEL_,
  _DEFAULT_
} enumFormatOutput;

typedef struct local_par {
  vt_names names;
  double pas;
  double max;

  double feretMin;
  double pixel;

  enumOutput typeOutput;
  enumHisto typeHisto;

  enumFormatOutput typeFormatOutput;
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



static void _PrintExcelPrn( typeEntry *theEntries,
			    int nbEntries,
			    char *name );

static void SortEntries( typeEntry *theEntries,
			 int left, 
			 int right );




static char *usage = "[file-in] [file-out]\n\
\t [-vol|-rad|-max|-med|-min | -minSmax] [-feretMin %lf]\n\
\t [-pixel %lf]\n\
\t [-pas %lf] [-vmax %lf] [-cum|[-ncum|-norm] [-matlab | -excel |-all]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'file-in' est '-', on prendra stdin\n\
\t si 'file-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n";

static char program[STRINGLENGTH];

#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  char string[256];
  local_par par;
  FILE *f, *fopen();
  int totalEntries = 0, nbEntries = 0;
  typeEntry *theEntries = (typeEntry *)NULL;
  int i, j;
  double theMax;
  int nbElems = 0;
  typeElemHisto1D *theHisto1D = (typeElemHisto1D *)NULL;

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  


  /*--- lecture du fichier d'entree ---*/
  if ( par.names.in[0] == '<' ) 
    VT_ErrorParse("not able to read standard input\n", 0 );


  {
    int ic;
    int readIsNotFinished = 1;
    typeEntry e;
    typeEntry *auxEntries = (typeEntry *)NULL;
    int allocatedEntries = 0;

    f = fopen( par.names.in, "r" );
    if ( f == (FILE*)NULL ) {
      sprintf ( string, "error in opening %s\n", par.names.in );
      VT_ErrorParse( string, 0 );
    }
    do {
      ic = fgetc( f );
      
      if ( ic == EOF ) {
	readIsNotFinished = 0;
	continue;
      }

      if ( (char)ic == '#' ) {
	do {
	  ic = fgetc( f );
	} while ( (char)ic != '\n' );
	continue;
      }
      
      if ( ungetc( ic, f ) != ic ) {
	fclose( f );
	sprintf ( string, "error in reading %s (ungetc)\n", par.names.in );
	VT_ErrorParse( string, 0 );
      }
      
      if ( fscanf( f, "%lf %lf %lf %lf %lf %lf %lf", 
		   &e.vol, &e.brd, &e.radVol, &e.radSrf, &e.max, &e.med, &e.min ) != 7 ) {
	fclose( f );
	sprintf ( string, "error in reading %s (entry #%d)\n", 
		  par.names.in, totalEntries-1 );
	VT_ErrorParse( string, 0 );
      }
      do {
	ic = fgetc( f );
      } while ( (char)ic != '\n' );
      totalEntries += 1;
      
      if ( totalEntries > allocatedEntries ) {
	allocatedEntries += _NBALLOCS_;
	auxEntries = (typeEntry *)realloc( theEntries, 
					   allocatedEntries*sizeof(typeEntry) );
	if ( auxEntries == (typeEntry *)NULL ) {
	  fclose( f );
	  if ( theEntries != (typeEntry *)NULL ) free( theEntries );
	  sprintf ( string, "error in allocating input array (size=%d)\n", 
		    allocatedEntries );
	  VT_ErrorParse( string, 0 );
	}
	theEntries = auxEntries;
      }
      theEntries[ totalEntries-1 ] = e;

    } while ( readIsNotFinished == 1 );
    
    fclose( f );
    fprintf( stderr, "%s: read %d entries in %s\n", 
	     argv[0], totalEntries, par.names.in );
    
    if ( totalEntries == 0 ) {
      sprintf ( string, "no entries found in %s\n", par.names.in );
      VT_ErrorParse( string, 0 );
    }
  }





  if ( par.typeFormatOutput == _EXCEL_ ) {
    _PrintExcelPrn( theEntries, totalEntries, par.names.out );
    free( theEntries );
    exit( 0 );
  }

  for ( i=0; i<totalEntries; i++ ) {
    theEntries[i].vol *= par.pixel * par.pixel * par.pixel;
    theEntries[i].radVol *= par.pixel;
    theEntries[i].max *= par.pixel;
    theEntries[i].med *= par.pixel;
    theEntries[i].min *= par.pixel;
  }

  /* supprime-t-on des valeurs ?
   */
  nbEntries = totalEntries;
  if ( par.feretMin > 0.0 ) {
    typeEntry e;
    for ( i=0; i<nbEntries; i++ ) {
      if ( theEntries[i].min < par.feretMin ) {
	e = theEntries[nbEntries-1];
	theEntries[nbEntries-1] = theEntries[i];
	theEntries[i] = e;
	nbEntries --;
	i--;
      }
    }
    fprintf( stderr, "%s: keep %d entries out of %d in %s (threshold was %f)\n", 
	     argv[0], nbEntries, totalEntries, par.names.in, par.feretMin );
  }





  /* construction des valeurs 
   */
  switch ( par.typeOutput ) {
  case _HISTO_MINSURMAX_ :
    for ( i=0; i<nbEntries; i++ ) {
      if ( theEntries[i].max > 0.0 )
	theEntries[i].value = theEntries[i].min / theEntries[i].max;
      else 
	theEntries[i].value = 1.0;
      if ( theEntries[i].value > 1.0 ) {
	fprintf( stderr, "Warning: min/max = %f for entry #%d\n", 
		 theEntries[i].value, i );
      }
    }
    break;
  case _HISTO_VOLUME_ :
    for ( i=0; i<nbEntries; i++ ) theEntries[i].value = theEntries[i].vol;
    break;
  case _HISTO_MAXFERET_ :
    for ( i=0; i<nbEntries; i++ ) theEntries[i].value = theEntries[i].max;
    break;
  case _HISTO_MEDFERET_ :
    for ( i=0; i<nbEntries; i++ ) theEntries[i].value = theEntries[i].med;
    break;
  case _HISTO_MINFERET_ :
    for ( i=0; i<nbEntries; i++ ) theEntries[i].value = theEntries[i].min;
    break;
  case _HISTO_RADIUS_ :
  default :
    for ( i=0; i<nbEntries; i++ ) theEntries[i].value = theEntries[i].radVol;
  }


  
  /* construction de l'histogramme
   */
  switch ( par.typeOutput ) {
  case _HISTO_VOLUME_ :
  case _HISTO_MAXFERET_ :
  case _HISTO_MEDFERET_ :
  case _HISTO_MINFERET_ :
  case _HISTO_RADIUS_ :
  default :
    /* histogramme 1D
     */
    theMax = theEntries[0].value;
    for ( i=1; i<nbEntries; i++ ) 
      if ( theMax < theEntries[i].value ) theMax = theEntries[i].value;
    fprintf( stderr, " ...   max value = %f\n", theMax );

    if ( theMax < par.max ) theMax = par.max;
    
    nbElems = (int)( theMax / par.pas + 0.5 );
    theHisto1D = (typeElemHisto1D *)malloc( (nbElems+1)*sizeof(typeElemHisto1D) );
    if ( theHisto1D == (typeElemHisto1D *)NULL ) {
      free( theEntries );
      VT_ErrorParse( "unable to allocate 1D histogram\n", 0 );
    }

    for ( i=0; i<=nbElems; i++ ) {
      /*
      theHisto1D[i].index = (double)i/(double)nbElems * theMax;
      */
      theHisto1D[i].index = (double)i * par.pas;
      theHisto1D[i].value = 0.0;
    }
    


    /* ici, on distribue chaque valeur de l'entree
       dans les cases
    */

    if ( 1 ) {
      for ( i=0; i<nbEntries; i++ ) {
	
	j = 0;
	do {
	  j++;
	} while ( j <= nbElems && theHisto1D[j].index < theEntries[i].value );
	/* on s'arrete quand 
	   j > nbElems || theHisto1D[j].index >= theEntries[i].value
	*/
	if ( j > nbElems ) {
	  theHisto1D[nbElems].value += 1.0;
	} else {
	  /* theHisto1D[j].index >= theEntries[i].value > theHisto1D[j-1].index 
	   */
	  theHisto1D[j-1].value += (theHisto1D[j].index - theEntries[i].value)
	    / (theHisto1D[j].index - theHisto1D[j-1].index);
	  theHisto1D[j].value += (theEntries[i].value - theHisto1D[j-1].index)
	    / (theHisto1D[j].index - theHisto1D[j-1].index);
	}
      }
    }



    if ( 0 ) {
      SortEntries( theEntries, 0, nbEntries-1 );
      theEntries[0].cumulatedValue = theEntries[0].value;
      for ( i=1; i<nbEntries; i++ ) {
	theEntries[i].cumulatedValue = theEntries[i-1].cumulatedValue + theEntries[i].value;
      }

      /* ... */

    }
    

    

    switch ( par.typeHisto ) {
    default :
      break;
    case _NORMALIZE_ :
      theMax = theHisto1D[0].value;
      for ( i=1; i<=nbElems; i++ )  
	if ( theHisto1D[i].value > theMax ) theMax = theHisto1D[i].value;
      for ( i=0; i<=nbElems; i++ )  
	 theHisto1D[i].value /= theMax;
      break;
    case _CUMULATIVE_ :
    case _CUMULATIVE_NORMALIZE_ :
      for ( i=1; i<=nbElems; i++ ) 
	theHisto1D[i].value += theHisto1D[i-1].value;
      if ( par.typeHisto == _CUMULATIVE_NORMALIZE_ ) 
	for ( i=0; i<=nbElems; i++ ) 
	  theHisto1D[i].value /= theHisto1D[nbElems].value;
      break;
    }


  } /* fin du switch */

  free( theEntries );


  if ( par.names.out[0] != '\0' && par.names.out[0] != '>') {
    switch ( par.typeFormatOutput ) {
    case _MATLAB_ :
      sprintf( string, "%s.m", par.names.out );
      break;
    default :
    case _DEFAULT_ :
      sprintf( string, "%s", par.names.out );
    }

    f = fopen( string, "w" );
    switch ( par.typeOutput ) {
    case _HISTO_VOLUME_ :
    case _HISTO_MAXFERET_ :
    case _HISTO_MEDFERET_ :
    case _HISTO_MINFERET_ :
    case _HISTO_RADIUS_ :
    default :

      switch ( par.typeFormatOutput ) {
      case _MATLAB_ :
	fprintf( f, "\n" );
	fprintf( f, "%% x = [ " );
	for ( i=0; i<=nbElems; i++ ) fprintf( f, "%g ", theHisto1D[i].index );
	fprintf( f, "];\n" );
	fprintf( f, "%s = [ ", par.names.out );
	for ( i=0; i<=nbElems; i++ ) fprintf( f, "%g ", theHisto1D[i].value );
	fprintf( f, "];\n" );
	fprintf( f, "\n" );
	fprintf( f, "%% figure;\n" );
	fprintf( f, "%% bar(x,%s);\n",par.names.out );
	fprintf( f, "\n" );
	break;

      case _ALL_ :
	for ( i=0; i<nbEntries; i++ )
	fprintf( f, "#%5d: %9.1f %8.3f %8.3f %8.3f %8.3f %8.3f : %f\n",
		 i+1, theEntries[i].vol, theEntries[i].radVol, theEntries[i].radSrf, 
		 theEntries[i].max, theEntries[i].med, theEntries[i].min,
		 theEntries[i].value );
	break;

      default :
      case _DEFAULT_ :
	for ( i=0; i<=nbElems; i++ ) {
	  fprintf( f, "%9f %f\n", theHisto1D[i].index, theHisto1D[i].value );
	}
      }
      
    }
  }



  /*--- operations eventuelles sur l'image d'entree ---*/

  
  /*--- initialisation de l'image resultat ---*/

  
  /*--- ecriture de l'image resultat ---*/

  
  /*--- liberations memoires ---*/
  switch ( par.typeOutput ) {
  case _HISTO_VOLUME_ :
  case _HISTO_MAXFERET_ :
  case _HISTO_MEDFERET_ :
  case _HISTO_MINFERET_ :
  case _HISTO_RADIUS_ :
  default :
    free( theHisto1D );
  }

  return( 1 );
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
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( program, argv[0] ) != 1 )
    VT_Error("Error while copying program name", (char*)NULL);
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  
  /*--- lecture des parametres ---*/
  i = 1; nb = 0;
  while ( i < argc ) {
    if ( argv[i][0] == '-' ) {
      if ( argv[i][1] == '\0' ) {
	if ( nb == 0 ) {
	  /*--- standart input ---*/
	  strcpy( par->names.in, "<" );
	  nb += 1;
	}
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

      /*---  ---*/
      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	par->typeFormatOutput = _MATLAB_;
      }
      else if ( strcmp ( argv[i], "-excel" ) == 0 ) {
	par->typeFormatOutput = _EXCEL_;
      }

      else if ( strcmp ( argv[i], "-all" ) == 0 ) {
	par->typeFormatOutput = _ALL_;
      }




      else if ( strcmp ( argv[i], "-vol" ) == 0 ) {
	par->typeOutput = _HISTO_VOLUME_;
      }
      else if ( strcmp ( argv[i], "-rad" ) == 0 ) {
	par->typeOutput = _HISTO_RADIUS_;
      }
      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
	par->typeOutput = _HISTO_MAXFERET_;
      }
      else if ( strcmp ( argv[i], "-med" ) == 0 ) {
	par->typeOutput = _HISTO_MEDFERET_;
      }
      else if ( strcmp ( argv[i], "-min" ) == 0 ) {
	par->typeOutput = _HISTO_MINFERET_;
      }
      
      else if ( strcmp ( argv[i], "-minSmax" ) == 0 || 
		strcmp ( argv[i], "-minsmax" ) == 0 ) {
	par->typeOutput = _HISTO_MINSURMAX_;
      }
      
      else if ( strcmp ( argv[i], "-minSmaxSmed" ) == 0 ) { 
	par->typeOutput = _HISTO_MINSMAX_MINSMED_;
      }
      
      else if ( strcmp ( argv[i], "-pixel" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -pixel...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pixel) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -pixel...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-feretMin" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -feretMin...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->feretMin) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -feretMin...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-pas" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -pas...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->pas) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -pas...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-vmax" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -vmax...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->max) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vmax...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-norm" ) == 0 ) {
	par->typeHisto = _NORMALIZE_;
      }
      else if ( strcmp ( argv[i], "-cum" ) == 0 ) {
	par->typeHisto = _CUMULATIVE_;
      }
       else if ( strcmp ( argv[i], "-ncum" ) == 0 ) {
	par->typeHisto = _CUMULATIVE_NORMALIZE_;
      }
     /*--- lecture du type de l'image de sortie ---*/

      /*--- option inconnue ---*/
      else {
	sprintf(text,"unknown option %s\n",argv[i]);
	VT_ErrorParse(text, 0);
      }
    }
    /*--- saisie des noms d'images ---*/
    else if ( argv[i][0] != 0 ) {
      if ( nb == 0 ) { 
	strncpy( par->names.in, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else if ( nb == 1 ) {
	strncpy( par->names.out, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else 
	VT_ErrorParse("too much file names when parsing\n", 0 );
    }
    i += 1;
  }
  
  /*--- s'il n'y a pas assez de noms ... ---*/
  if (nb == 0) {
    strcpy( par->names.in,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  if (nb == 1)
    strcpy( par->names.out, ">" );  /* standart output */
  

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
  VT_Names( &(par->names) );
  par->typeOutput = _HISTO_RADIUS_;
  par->pas = 0.1;
  par->max = 1.0;
  par->pixel = 1.0;
  par->feretMin = -1.0;
  par->typeHisto = _SIMPLE_;
  par->typeFormatOutput = _DEFAULT_;
}







static void _PrintExcelPrn( typeEntry *theEntries,
			    int nbEntries,
			    char *name )
{
  char string[256];
  FILE *f, *fopen();
  int i;

  sprintf( string, "%s.prn", name );
  f = fopen( string, "w" );
  /*
  fprintf ( f, "analyse   ELF       ATO       CRRA      vers      v_4.61\n" );
  fprintf ( f, "type de fi         0\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "echelle :1  0.065934mm\n" );
  fprintf ( f, "\n" );
  */
  /*
  fprintf ( f, "image_nb\t cell_nb\t volume\t proba\t Ft_max\t Ft_moy\t Ft_min\t Grossisse\t Ech_X\t Ech_Y\t nimage\t diam.moyen\n");
  
  for ( i=0; i<nbEntries; i++ ) {
    fprintf ( f, "1\t %9d\t %9d\t 1.0\t %9f\t %9f\t %9f\t 0\t 0.065934\t 0.065934\n",
	      i+1, (int)(theEntries[i].vol +0.5), 
	      theEntries[i].radVol,
	      theEntries[i].max,
	      theEntries[i].min );
	      
  }
  */
    fprintf ( f, "cell_nb\t volume\t Ft_max\t Ft_moy\t Ft_min\t diam.moyen\n");
  
  for ( i=0; i<nbEntries; i++ ) {
    fprintf ( f, "%9d\t %9f\t %9f\t %9f\t %9f\t %9f\n",
	      i+1, theEntries[i].vol, 
	      theEntries[i].max ,theEntries[i].med,
	      theEntries[i].min, 2.0*theEntries[i].radVol );
	      
  }


  /*
  fprintf ( f, "#fin_analyse\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, " minimas/maximas\n" );
  fprintf ( f, "surf min -  0.004347 max ->     6.899167\n" );
  fprintf ( f, "d_eq min -  0.074399 max ->     2.963831\n" );
  fprintf ( f, "\n" );
  fprintf ( f, "\n" );
  fprintf ( f, " pretraitement code pre_tre = 0\n" );
  fprintf ( f, "trame carree\n" );
  */
  fclose( f );
}








static void SortEntries( typeEntry *tab,
			 int left, 
			 int right )
{
  int i, last;
  typeEntry tmp;

  if ( left >= right ) return;

  tmp = tab[left];   tab[left] = tab[(left+right)/2];   tab[(left+right)/2] = tmp;
  
  last = left;
  for ( i = left+1; i <= right; i++ )       
    if ( tab[i].value < tab[left].value ) {
      tmp = tab[++last];   tab[last] = tab[i];   tab[i] = tmp;
    }

  tmp = tab[left];   tab[left] = tab[last];   tab[last] = tmp;
  
  SortEntries( tab, left, last-1 );
  SortEntries( tab, last+1, right );
}
