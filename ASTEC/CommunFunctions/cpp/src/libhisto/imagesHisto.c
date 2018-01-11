/*************************************************************************
 * minimum.c -
 *
 * $Id: imagesHisto.c,v 1.2 2002/10/18 18:02:07 greg Exp $
 *
 * Copyright (c) INRIA 1999
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * 
 *
 * ADDITIONS, CHANGES
 *
 */

#include <vt_common.h>
#include <vt_sliceHisto.h>



typedef enum {
  _NONE_,
  _MOMENTS_,
  _DIRECTE_SSD_,
  _INVERSE_SSD_,
  _SYMETRIE_SSD_,
  _DIRECTE_Correlation_,
  _INVERSE_Correlation_,
  _SYMETRIE_Correlation_,
  _DIRECTE_Vraisemblance_,
  _INVERSE_Vraisemblance_,
  _SYMETRIE_Vraisemblance_,
  _MINIMUM_Vraisemblance_
} enumComputation;


typedef enum {
  _CONSTANT_,
  _LINEAR_
} enumCompensation;


typedef struct local_par {
  vt_names names;
  char auxname[ STRINGLENGTH ];
  char maskname[ STRINGLENGTH ];
  double psigma;
  int type;
  enumComputation typeComputation;
  enumComputation initComputation;
  enumCompensation typeCompensation;
  double a;
  double b;
  int lscales;
  int nscales;
} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-1] [image-2] [image-out]\n\
 'image-out' will be 'image-1' after intensity transformation\n\
    with an affine function: 'A' * intensity + 'B' \n\
    with A = a_arg * a_cmp and B = a_arg * b_cmp + b_arg\n\
    (a_cmp,b_cmp) are computed here\n\
    (a_arg,b_arg) are passed as arguments\n\
    typically they have been used to resample image-2\n\
\t [-matlab %s] [-sigma %lf]\n\
\t [-mode ssd|issd|sssd | corr|icorr|scorr | vrai|ivrai|svrai]\n\
\t [-function|-ftn  linear|cst] [-scales | -s %d] [-nscales | -ns %d]\n\
\t [-init moments] [-a %lf] [-b %lf]\n\
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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.2 $ $Date: 2002/10/18 18:02:07 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image1, *image2;
  vt_image *immask = NULL;
  int *histo1, max1;
  int *histo2, max2;
  int *theHisto1;
  int *theHisto2;
  int length;
  double theCoeff[2];

  double (*objective_func)(double *, void *) = NULL;
  typeIntensityTrsf dir_transfo_func = NULL;
  typeIntensityTrsf inv_transfo_func = NULL;
  int nparam = 0;

  int fd = 0;
  FILE *fp = NULL;
  char *longname;
  int id = 0;




  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #1\n", 0);
  
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #2\n", 0);
  


  if ( par.auxname[0] != '\0' ) {
    int i, k;
    
    longname = malloc( strlen(par.auxname)+ 10 );
    sprintf( longname, "%s.raw", par.auxname );
    fd = creat( longname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    sprintf( longname, "%s.m", par.auxname );
    fp = fopen( longname, "w" );
    free( longname );

    fprintf( fp, "\n" );
    fprintf( fp, "%% " );
    fprintf( fp, "\n" );
    fprintf( fp, "%%" );
    for ( i=0; i<argc; i++ ) fprintf( fp, " %s", argv[i] );
    fprintf( fp, "\n" );
    fprintf( fp, "%% " );
    fprintf( fp, "\n" );

    fprintf( fp, "\n" );
    k = strlen( par.auxname );
    for ( i = k-1; i >= 0 && par.auxname[i] != '/' ; i-- )
      ;
    fprintf( fp, "fid = fopen('%s.raw', 'r' );\n", &(par.auxname[i+1]) );
    fprintf( fp, "\n" );
  }





  /* calcul de l'histogramme #1
     puis #2
  */
  if ( par.maskname[0] != '\0' )
    immask = _VT_Inrimage( par.maskname );

  histo1 = _GetImageHisto( image1, immask, &max1 );
  histo2 = _GetImageHisto( image2, immask, &max2 );
  if ( max1 > max2 ) {
    length = max1+1;
    theHisto1 = histo1;
    theHisto2 = (int*)calloc(length, sizeof(int) );
    memcpy( theHisto2, histo2, (max2+1)*sizeof(int) );
    free( histo2 );
  }
  else if ( max2 > max1 ) {
    length = max2+1;
    theHisto2 = histo2;
    theHisto1 = (int*)calloc(length, sizeof(int) );
    memcpy( theHisto1, histo1, (max1+1)*sizeof(int) );
    free( histo1 );    
  }
  else {
    length = max1+1;
    theHisto1 = histo1;
    theHisto2 = histo2;
  }

  theHisto1[0] = 0;
  theHisto2[0] = 0;


  /* recalage des histogrammes
   */

  switch( par.typeCompensation ) {
  default :
  case _LINEAR_ :
    dir_transfo_func = _affine;
    inv_transfo_func = _inv_affine;
    nparam = 2;
    break;
  case _CONSTANT_ :
    dir_transfo_func = _constant; 
    inv_transfo_func = _inv_constant;
    nparam = 1;
  }


  switch ( par.initComputation ) {
  default :
  case _NONE_ :
    theCoeff[0] = 0.0;
    theCoeff[1] = 1.0;
    break;
 case _MOMENTS_ :
   theCoeff[0] = 0.0;
   theCoeff[1] = 1.0;
   switch( par.typeCompensation ) {
   default :
   case _LINEAR_ :
     _initAffineTrsfBetweenTwoHisto( theCoeff, theHisto2, length,
				     theHisto1, length );
     break;
   case _CONSTANT_ :
     _initConstantTrsfBetweenTwoHisto( theCoeff, theHisto2, length,
				       theHisto1, length );
     break;
   }
   printf( "Initial coefficients A=%lf B=%lf\n", theCoeff[1], theCoeff[0] );
  }

  switch ( par.typeComputation ) {
  default :
    objective_func = NULL;
    break;
  case _DIRECTE_SSD_ :    objective_func = _DirecteSSD;   break;
  case _INVERSE_SSD_ :    objective_func = _InverseSSD;   break;
  case _SYMETRIE_SSD_ :   objective_func = _SymetrieSSD;   break;
  case _DIRECTE_Correlation_ :    objective_func = _DirecteCorrelation;   break;
  case _INVERSE_Correlation_ :    objective_func = _InverseCorrelation;   break;
  case _SYMETRIE_Correlation_ :   objective_func = _SymetrieCorrelation;   break;
  case _DIRECTE_Vraisemblance_ :    objective_func = _DirecteVraisemblance;   break;
  case _INVERSE_Vraisemblance_ :    objective_func = _InverseVraisemblance;   break;
  case _SYMETRIE_Vraisemblance_ :   objective_func = _SymetrieVraisemblance;   break;
  case _MINIMUM_Vraisemblance_ :   objective_func = _MinimumVraisemblance;   break;
  }

  /* les coefficients s'appliquent a image 2
  */


  switch ( par.typeComputation ) {

  default :
  case _NONE_ :
    if ( par.initComputation == _NONE_ ) break;
    break;

  case _DIRECTE_SSD_ :
  case _INVERSE_SSD_ :
  case _SYMETRIE_SSD_ :

  case _DIRECTE_Correlation_ :
  case _INVERSE_Correlation_ :
  case _SYMETRIE_Correlation_ :

  case _DIRECTE_Vraisemblance_ :
  case _INVERSE_Vraisemblance_ :
  case _SYMETRIE_Vraisemblance_ :
  case _MINIMUM_Vraisemblance_ :
    
    if ( par.lscales <= 0 ) {
      _evalTrsfBetweenTwoHisto( theCoeff, theHisto2, length, theHisto1, length,
				dir_transfo_func, inv_transfo_func, nparam,
				objective_func, par.psigma );
    } 
    else {
      /* lscales : plus grande echelle
	 nscales : nombre d'echelle
      */
      _multiScaleEvalTrsfBetweenTwoHisto( theCoeff, theHisto2, length, theHisto1, length,
					    dir_transfo_func, inv_transfo_func, nparam,
					    objective_func, par.psigma, par.lscales, par.nscales );
    }
  }

  printf( "compensation de '%s' avec '%s' : j = %f + %f * i\n",
	  par.names.in, par.names.ext, theCoeff[0], theCoeff[1] );
  

  
  if ( par.auxname[0] != '\0' ) {
    _PrintOneHistoForMatlab( fd, fp, NULL, theHisto1, length, id++, _AFFINE_ );
    _PrintOneHistoForMatlab( fd, fp, NULL, theHisto2, length, id++, _AFFINE_ );
    _PrintOneHistoForMatlab( fd, fp, theCoeff, theHisto1, length, id++, _AFFINE_ );
    


    fprintf( fp, "\n" );
    fprintf( fp, "figure;\n" );
    fprintf( fp, "hold on;\n" );

    fprintf( fp, "plot( [0:%d], HISTO%d, 'b-' );\n", length-1, id-3 ); /* theHisto1 */
    fprintf( fp, "plot( [0:%d], HISTO%d, 'm-' );\n", length-1, id-2 ); /* theHisto2 */
    fprintf( fp, "plot( %f + %f*[0:%d], HISTO%d, 'r-' );\n", 
	     theCoeff[0], theCoeff[1], length-1, id-2 );
    fprintf( fp, "legend('image 1', 'image 2', 'image 2 transformed' );\n" );
    
    fprintf( fp, "hold off;\n" );
    fprintf( fp, "\n" );
  }
  


  /* les coefficients s'appliquent a image 2
     on les inverse pour les appliquer a image 1
     et etre coherent avec Robust_affine_equalization
   */
  
  theCoeff[0] = - theCoeff[0] / theCoeff[1];
  theCoeff[1] = 1.0 / theCoeff[1];
  printf( "Computed coefficients A=%lf B=%lf\n", theCoeff[1], theCoeff[0] );
  
  theCoeff[0] = par.a * theCoeff[0] + par.b;
  theCoeff[1] = par.a * theCoeff[1];
  
  printf( "Coefficients for interpolation A=%f B=%f\n", theCoeff[1], theCoeff[0] );
  



  if ( par.auxname[0] != '\0' ) {
    fprintf( fp, "\n" );
    fprintf( fp, "fclose(fid);\n" );
    fclose( fp );
    close( fd );
  }

  VT_FreeImage( image2 );
  VT_Free( (void**)&image2 );


  if ( par.names.out[0] != '\0' && par.names.out[0] != '>' ) {
    /*--- initialisation de l'image resultat ---*/

    printf( "transformation de '%s'\n", par.names.in );
     
    switch ( image1->type ) {
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
      break;
    case UCHAR :
      {
	u8 ***buf = (u8***)image1->array;
	int x,y,z,j;
	double v;
	for ( z=0; z<image1->dim.z; z++ ) {
	  for ( y=0; y<image1->dim.y; y++ )
	  for ( x=0; x<image1->dim.x; x++ ) {
	    v = theCoeff[0] + theCoeff[1] * buf[z][y][x];
	    if ( v < 0.0 ) {
	      buf[z][y][x] = 0.0;
	    } 
	    else {
	      j = (int)( v + 0.5 );
	      if ( j > 255 ) buf[z][y][x] = 255;
	      else           buf[z][y][x] = j;
	    }
	  }
	}
      }
      break;
    case USHORT :
      {
	u16 ***buf = (u16***)image1->array;
	int x,y,z,j;
	double v;
	for ( z=0; z<image1->dim.z; z++ ) {
	  for ( y=0; y<image1->dim.y; y++ )
	  for ( x=0; x<image1->dim.x; x++ ) {
	    v = theCoeff[0] + theCoeff[1] * buf[z][y][x];
	    if ( v < 0.0 ) {
	      buf[z][y][x] = 0.0;
	    } 
	    else {
	      j = (int)( v + 0.5 );
	      if ( j > 65535 ) buf[z][y][x] = 65535;
	      else           buf[z][y][x] = j;
	    }
	  }
	}
      }
      break;
     case SSHORT :
      {
	s16 ***buf = (s16***)image1->array;
	int x,y,z,j;
	double v;
	for ( z=0; z<image1->dim.z; z++ ) {
	  for ( y=0; y<image1->dim.y; y++ )
	  for ( x=0; x<image1->dim.x; x++ ) {
	    v = (theCoeff[0] + theCoeff[1] * (32768+(int)buf[z][y][x])) - 32768;
	    if ( v < 0.0 ) {
	      j = (int)( v - 0.5 );
	      if ( j < -32768 ) buf[z][y][x] = 32768;
	      else           buf[z][y][x] = j;
	    } 
	    else {
	      j = (int)( v + 0.5 );
	      if ( j > 32767 ) buf[z][y][x] = 32767;
	      else           buf[z][y][x] = j;
	    }
	  }
	}
      }
      break;
   }

    sprintf( image1->name, "%s", par.names.out );
    VT_WriteInrimage( image1 );

  }

  VT_FreeImage( image1 );
  VT_Free( (void**)&image1 );

  return( 1 );
}








static void VT_Parse( int argc, 
		      char *argv[], 
		      local_par *par )
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
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }





      else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	
	if ( strcmp ( argv[i], "moments" ) == 0 ||
	     strcmp ( argv[i], "mmnts" ) == 0 ) {
	  par->initComputation = _MOMENTS_;
	} 
      }


      else if ( strcmp ( argv[i], "-function" ) == 0 
		|| strcmp ( argv[i], "-ftn" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -function...\n", 0 );

	if ( strcmp ( argv[i], "linear" ) == 0 ) {
	  par->typeCompensation = _LINEAR_;
	}
	else if ( strcmp ( argv[i], "constant" ) == 0 
		  || strcmp ( argv[i], "cst" ) == 0 ) {
	  par->typeCompensation = _CONSTANT_;
	}
	else {
	  VT_ErrorParse( "unknown function...\n", 0 );
	}
      }

      else if ( strcmp ( argv[i], "-mode" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mode...\n", 0 );

	if ( strcmp ( argv[i], "ssd" ) == 0 ) {
	  par->typeComputation = _DIRECTE_SSD_;
	} 
	else if ( strcmp ( argv[i], "issd" ) == 0 ||
		  strcmp ( argv[i], "inv-ssd" ) == 0 ) {
	  par->typeComputation = _INVERSE_SSD_;
	} 
	else if ( strcmp ( argv[i], "sssd" ) == 0 ||
		  strcmp ( argv[i], "sym-ssd" ) == 0 ) {
	  par->typeComputation = _SYMETRIE_SSD_;
	} 
	
	else if ( strcmp ( argv[i], "corr" ) == 0 ) {
	  par->typeComputation = _DIRECTE_Correlation_;
	} 
	else if ( strcmp ( argv[i], "icorr" ) == 0 ||
		  strcmp ( argv[i], "inv-corr" ) == 0 ) {
	  par->typeComputation = _INVERSE_Correlation_;
	} 
	else if ( strcmp ( argv[i], "scorr" ) == 0 ||
		  strcmp ( argv[i], "sym-corr" ) == 0 ) {
	  par->typeComputation = _SYMETRIE_Correlation_;
	} 
	
	else if ( strcmp ( argv[i], "vrai" ) == 0 ) {
	  par->typeComputation = _DIRECTE_Vraisemblance_;
	} 
	else if ( strcmp ( argv[i], "ivrai" ) == 0 ||
		  strcmp ( argv[i], "inv-vrai" ) == 0 ) {
	  par->typeComputation = _INVERSE_Vraisemblance_;
	} 
	else if ( strcmp ( argv[i], "svrai" ) == 0 ||
		  strcmp ( argv[i], "sym-vrai" ) == 0 ) {
	  par->typeComputation = _SYMETRIE_Vraisemblance_;
	}
	else if ( strcmp ( argv[i], "mvrai" ) == 0 ||
		  strcmp ( argv[i], "min-vrai" ) == 0 ) {
	  par->typeComputation = _MINIMUM_Vraisemblance_;
	} 
	
	else {
	  VT_ErrorParse( "unknown mode...\n", 0 );
	}
	  
      }

      else if ( strcmp ( argv[i], "-matlab" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -matlab...\n", 0 );
	strncpy( par->auxname, argv[i], STRINGLENGTH );  
      }



      else if ( strcmp ( argv[i], "-mask" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mask...\n", 0 );
	strncpy( par->maskname, argv[i], STRINGLENGTH );  
      }



      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->psigma) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sigma...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-a" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -a...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->a) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -a...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-b" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -b...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->b) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -b...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-scales" ) == 0 
		|| strcmp ( argv[i], "-s" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -scales...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->lscales) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -scales...\n", 0 );
	if ( par->nscales == 0 ) par->nscales = par->lscales + 1;
      }

      else if ( strcmp ( argv[i], "-nscales" ) == 0 
		|| strcmp ( argv[i], "-ns" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -nscales...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->nscales) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -nscales...\n", 0 );
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
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
	nb += 1;
      }
      else if ( nb == 2 ) {
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
    VT_ErrorParse("not enough file names when parsing\n", 0 );
  }
  else if (nb == 1) {
    strcpy( par->names.ext,  "<" );  /* standart input */
    strcpy( par->names.out, ">" );  /* standart output */
  }
  else if (nb == 2)
    strcpy( par->names.out, ">" );  /* standart output */
  
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) )  par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) )  par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) )  par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) )  par->type = FLOAT;
  /* if ( par->type == TYPE_UNKNOWN ) VT_Warning("no specified type", program); */
}






static void VT_ErrorParse( char *str, int flag )
{
  (void)fprintf(stderr,"Usage : %s %s\n",program, usage);
  if ( flag == 1 ) (void)fprintf(stderr,"%s",detail);
  (void)fprintf(stderr,"Erreur : %s",str);
  exit(0);
}








static void VT_InitParam( local_par *par )
{
  VT_Names( &(par->names) );
  par->type = TYPE_UNKNOWN;
  par->psigma = 5.0;
  par->typeComputation = _NONE_;
  par->initComputation = _NONE_;
  par->typeCompensation = _LINEAR_;
  par->auxname[0] = '\0';
  par->maskname[0] = '\0';
  par->a = 1.0;
  par->b = 0.0;
  par->nscales = 0;
  par->lscales = 0;
}
