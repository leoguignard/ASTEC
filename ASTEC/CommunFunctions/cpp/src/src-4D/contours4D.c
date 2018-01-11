#include <vt_common.h>
#include <vt_common.h>
#include <vt_contours.h>
#include <vt_recfilters4D.h>
#include <vt_image4D.h>

typedef enum {
  EXTREMA = 1,
  NORME = 2
} typeOutput;


typedef struct {
  vt_names names;
  vt_contours rpar;
  vt_contours tpar;
  int type;
  typeOutput output;
  typeBoolean Tsmooth;
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
\t [-sigma %f ] [-st %f ] [-cont %d ]\n\
\t [-norme] [-tsmooth] [-3D]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -sigma       : sigma pour l'approximation de la gaussienne\n\
\t -st          : sigma selon T\n\
\t -cont        : points ajoutes aux bords\n\
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
  
  vt_image4D imtmp4D, *input4D;
  




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
  if ( par.type == TYPE_UNKNOWN ) par.type = image4D.type;



  




  /* lissage temporel */
  if ( ((par.tpar.dim != VT_4D) || (par.rpar.dim != VT_4D)) 
       && (par.Tsmooth == True) ) {

    vt_recfilters recPar;

    if ( VT_AllocAndInitImage4D( &imtmp4D, NULL,
			       image4D.dim.x, image4D.dim.y, image4D.dim.z,
			       image4D.dimt, FLOAT ) != 1 ) {
      VT_FreeImage4D( &image4D );
      VT_ErrorParse("unable to allocate auxiliary image\n", 0);
    }
    


    VT_RecFilters( &recPar );
    switch ( par.tpar.type_filter ) {
    case VT_RECFILTERS_DERICHE :
    case VT_RECGAUSSIAN_DERICHE :
      recPar.type_filter = par.tpar.type_filter;
      break;
    default :
      recPar.type_filter = VT_RECFILTERS_DERICHE;
    }
    recPar.value_coefficient = par.tpar.value_coefficient;
    recPar.length_continue = par.tpar.length_continue;
    recPar.derivative.x = VT_DERIVATIVE_0;
    recPar.derivative.y = VT_NODERIVATIVE;
    recPar.derivative.z = VT_NODERIVATIVE;


    if ( VT_RecFilterTOnImage4D( &image4D, &imtmp4D, &recPar ) != 1 ) {
      VT_FreeImage4D( &image4D );
      VT_FreeImage4D( &imtmp4D );
      VT_ErrorParse("unable to smooth input image\n", 0);
    }

    VT_FreeImage4D( &image4D );
    input4D = &imtmp4D;

    
  } else {

    input4D = &image4D;

  }




  /*--- initialisation de l'image resultat ---*/
  VT_Name4D( &nameRes );
  
  if ( par.names.ext[0] == '\0' ) base = par.names.out;
  else                            base = par.names.ext;

  if ( VT_CreateName4D( &nameRes, base, ".inr", input4D->dimt ) != 1 ) {
    VT_FreeImage4D( input4D );
    VT_ErrorParse("unable to create output names\n", 0);
  }
  if ( VT_AllocAndInitImage4D( &imres4D, &nameRes, 
			       input4D->dim.x, input4D->dim.y, input4D->dim.z,
			       input4D->dimt, par.type ) != 1 ) {
    
    VT_FreeImage4D( input4D );
    VT_FreeName4D( &nameRes );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }






  switch ( par.output ) {
  case EXTREMA :
  default :

    if ( VT_MaximaGradient4D( input4D, &imres4D, 
			      &(par.rpar), &(par.tpar) ) != 1 ) {

      VT_FreeImage4D( &imres4D );
      VT_FreeImage4D( input4D );
      VT_FreeName4D( &nameRes );
      VT_ErrorParse("unable to extract maxima from input image\n", 0);
    }
    break;

  case NORME :
    if ( VT_NormeGradient4DImage4D( input4D, &imres4D, 
				    &(par.rpar), &(par.tpar), 
				    VT_DERIVATIVE_1_CONTOURS ) != 1 ) {
    VT_FreeImage4D( &imres4D );
    VT_FreeImage4D( input4D );
    VT_FreeName4D( &nameRes );
    VT_ErrorParse("unable to compute norme from input image\n", 0);
    }

    break;
  }

  


  VT_FreeImage4D( input4D );




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
	  VT_Recfilters4DVerbose();
	  _VT_VERBOSE_ = 1;
	}
	else if ( strcmp ( argv[i], "-D" ) == 0 ) {
	  _VT_DEBUG_ = 1;
	}

	/*--- traitement eventuel de l'image d'entree ---*/
	else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	  par->names.inv = 1;
	}
	else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	  par->names.swap = 1;
	}
	

	else if ( strcmp ( argv[i], "-base" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -base...\n", 0 );
	  status = sscanf( argv[i],"%s", par->names.ext );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -base...\n", 0 );
	}
	


	/*--- sigma ---*/
	else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -sigma...\n", 0 );
	  status = sscanf( argv[i],"%f",&(par->rpar.value_coefficient.x) );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -sigma...\n", 0 );
	  par->rpar.value_coefficient.z = par->rpar.value_coefficient.y = par->rpar.value_coefficient.x;
	  par->tpar.value_coefficient.z = par->tpar.value_coefficient.y = par->tpar.value_coefficient.x = par->rpar.value_coefficient.x;

	  switch ( par->rpar.type_filter ) {
	  case TYPE_UNKNOWN :
	  case VT_RECGAUSSIAN_DERICHE :
	    par->rpar.type_filter = VT_RECGAUSSIAN_DERICHE;
	    break;
	  default :
	    VT_ErrorParse( "parsing -sigma...\n", 0 );
	  }
	  switch ( par->tpar.type_filter ) {
	  case TYPE_UNKNOWN :
	  case VT_RECGAUSSIAN_DERICHE :
	    par->tpar.type_filter = VT_RECGAUSSIAN_DERICHE;
	    break;
	  default :
	    VT_ErrorParse( "parsing -sigma...\n", 0 );
	  }
	}

	else if ( strcmp ( argv[i], "-st" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -st...\n", 0 );
	  status = sscanf( argv[i],"%f",&(par->tpar.value_coefficient.x) );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -st...\n", 0 );
	  par->tpar.value_coefficient.z = par->tpar.value_coefficient.y = par->tpar.value_coefficient.x;
	  switch ( par->tpar.type_filter ) {
	  case TYPE_UNKNOWN :
	  case VT_RECGAUSSIAN_DERICHE :
	    par->tpar.type_filter = VT_RECGAUSSIAN_DERICHE;
	    break;
	  default :
	    VT_ErrorParse( "parsing -st...\n", 0 );
	  }
	}

	/*--- bordure ---*/
	else if ( strcmp ( argv[i], "-cont" ) == 0 ) {
	  i += 1;
	  if ( i >= argc)    VT_ErrorParse( "parsing -cont...\n", 0 );
	  status = sscanf( argv[i],"%d",&(par->rpar.length_continue.x) );
	  if ( status <= 0 ) VT_ErrorParse( "parsing -cont...\n", 0 );
	  par->rpar.length_continue.z = par->rpar.length_continue.y = par->rpar.length_continue.x;
	  par->tpar.length_continue.z = par->tpar.length_continue.y = par->tpar.length_continue.x = par->rpar.length_continue.x;
	}
	


	
	else if ( strcmp ( argv[i], "-norme" ) == 0 ) {
	  par->output = NORME;
	}
	else if ( (strcmp ( argv[i], "-Tsmooth" ) == 0) || (strcmp ( argv[i], "-tsmooth" ) == 0) ) {
	  par->Tsmooth = True;
	}
	else if ( strcmp ( argv[i], "-3D" ) == 0 ) {
	  par->rpar.dim = par->tpar.dim = VT_3D;
	}
	else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	  par->rpar.dim = par->tpar.dim = VT_2D;
	}
	


	/*--- marta ? ---*/
	else if ( strcmp ( argv[i], "-marta" ) == 0 ) {
	  par->rpar.type_filter = VT_RECGAUSSIAN_MARTA;
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
    VT_Contours( &(par->rpar) );
    par->rpar.type_filter = VT_RECGAUSSIAN_DERICHE;
    VT_Contours( &(par->tpar) );
    par->tpar.type_filter = VT_RECGAUSSIAN_DERICHE;
    par->type = TYPE_UNKNOWN;
    par->output = EXTREMA;
    par->rpar.dim = par->tpar.dim = VT_4D;
    par->Tsmooth = False;
}
