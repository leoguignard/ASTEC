#include <vt_common.h>
#include <vt_common.h>
#include <vt_contours.h>
#include <vt_recfilters4D.h>
#include <vt_image4D.h>
#include <vt_connexe4D.h>

typedef enum {
  EXTREMA = 1,
  NORME = 2
} typeOutput;


typedef struct {
  vt_names names;
  float sb;
  float sh;
  int type;
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

static char *usage = "[image-in] [image-out] [-base %s]\n\
\t [-sb %f] [-sh %f]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
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
  local_par par;

  FILE *f, *fopen();
  typeBoolean readingStandardInput = False;

  vt_name4D name, nameRes;
  char *base;
  vt_image4D image4D, imres4D;
  
  




  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  


  /*--- lecture des noms de l'image 4D d'entree ---*/
  VT_Name4D( &name );
  if ( (par.names.in[0] == '\0') || ((par.names.in[0] == '<') && (par.names.in[1] == '\0')) ) {
    f = stdin;
    readingStandardInput = True;
  } else {
    f = fopen( par.names.in , "r" );
  }

  if ( VT_ReadNameInr4D( f, &name ) != 1 ) {
    VT_ErrorParse("unable to read input names\n", 0);
  }

  if ( readingStandardInput == False ) fclose( f );





  /*--- lecture des images d'entree ---*/
  if ( VT_AllocAndReadImage4D( &image4D, &name ) != 1 ) {
    VT_FreeName4D( &name );
    VT_ErrorParse( "unable to read 4D image\n", 0 );
  }
  VT_FreeName4D( &name );
  if ( par.type == TYPE_UNKNOWN ) par.type = UCHAR;



  


  /*--- initialisation de l'image resultat ---*/
  VT_Name4D( &nameRes );
  
  if ( par.names.ext[0] == '\0' ) base = par.names.out;
  else                            base = par.names.ext;

  if ( VT_CreateName4D( &nameRes, par.names.out, ".inr", image4D.dimt ) != 1 ) {
    VT_FreeImage4D( &image4D );
    VT_ErrorParse("unable to create output names\n", 0);
  }
  if ( VT_AllocAndInitImage4D( &imres4D, &nameRes, 
			       image4D.dim.x, image4D.dim.y, image4D.dim.z,
			       image4D.dimt, par.type ) != 1 ) {
    
    VT_FreeImage4D( &image4D );
    VT_FreeName4D( &nameRes );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }


  if ( VT_Hysteresis4D( &image4D, &imres4D, par.sb, par.sh ) != 1 ) {
    VT_FreeImage4D( &image4D );
    VT_FreeImage4D( &imres4D );
    VT_FreeName4D( &nameRes );
    VT_ErrorParse("unable to compute hysteresis thresholding\n", 0);
  }

  VT_FreeImage4D( &image4D );



  if ( VT_WriteImage4D( &imres4D ) != 1 ) {
    VT_FreeImage4D( &imres4D );
    VT_FreeName4D( &nameRes );
    VT_ErrorParse("unable to write output image\n", 0);    
  }
  VT_FreeImage4D( &imres4D );


  VT_WriteNameInr4D( &nameRes, par.names.out );
  VT_FreeName4D( &nameRes );

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
    int o=0, s=0, r=0;
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



	else if ( strcmp ( argv[i], "-base" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -base...\n", 0 );
	  status = sscanf( argv[i],"%s", par->names.ext );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -base...\n", 0 );
	}



	/*--- traitement eventuel de l'image d'entree ---*/
	else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	  par->names.inv = 1;
	}
	else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	  par->names.swap = 1;
	}
	
	/*--- seuils ---*/
	else if ( strcmp ( argv[i], "-sb" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -sb...\n", 0 );
	  status = sscanf( argv[i],"%f",&(par->sb) );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -sb...\n", 0 );
	}
	else if ( strcmp ( argv[i], "-sh" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -sh...\n", 0 );
	  status = sscanf( argv[i],"%f",&(par->sh) );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -sh...\n", 0 );
	}



	
	/*--- lecture du type de l'image de sortie ---*/
	else if ( strcmp ( argv[i], "-r" ) == 0 ) {
	  r = 1;
	}
	else if ( strcmp ( argv[i], "-s" ) == 0 ) {
	  s = 1;
	}
	else if ( strcmp ( argv[i], "-o" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -o...\n", 0 );
	  status = sscanf( argv[i],"%d",&o );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -o...\n", 0 );
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

    /*--- type de l'image resultat ---*/
    if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
    if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
    if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
    if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
    if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
    if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
    if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program);
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
    par-> sb = 0;
    par-> sh = 0;
    par->type = TYPE_UNKNOWN;
}
