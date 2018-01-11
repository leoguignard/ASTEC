/*************************************************************************
 * slicesHisto-2.c -
 *
 * $Id: slicesHisto-2.c,v 1.5 2002/10/18 18:02:07 greg Exp $
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


typedef struct local_par {
  vt_names names;
  int threshold;
  int zref;
  int z;
  double psigma;
  int type;
  enumComputation typeComputation;
  enumComputation initComputation;
} local_par;







/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
\t [-matlab %s] [-sigma %lf] [-zref %d] [-z %d]\n\
\t [-init moments] [-threshold | -thres |-th %d] \n\
\t [-mode ssd|issd|sssd | corr|icorr|scorr | vrai|ivrai|svrai|mvrai]\n\
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
 $Revision: 1.5 $ $Date: 2002/10/18 18:02:07 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imcff;
  int **theHisto;
  int i, j, iref;

  int id = 0;


  int ncoeff = 2;
  double **theCoeff;
  int length, intensity_max = 0;

  int fd = 0;
  FILE *fp = NULL;
  char *longname;

  double (*func)(double *, void *) = NULL;

  double icoeff[2], coeff[2];


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  

  theHisto = _GetSlicesHisto( image, &intensity_max );
  for (i=0; i<image->dim.z; i++ ) {
    for ( j=0; j<par.threshold; j++ )
      theHisto[i][j] = 0;
  }
  length = intensity_max+1;


  if ( par.zref <= 0 || par.zref >= image->dim.z ) {
    iref = image->dim.z / 2; 
    /* iref =_GetReferenceSlice( theHisto, image->dim.z, intensity_max+1 ); */
  }
  else
    iref = par.zref;

  if ( par.names.ext[0] != '\0' ) {
    longname = malloc( strlen(par.names.ext)+ 10 );
    sprintf( longname, "%s.coeff.inr", par.names.ext );
    VT_InitImage( &imcff, longname, ncoeff, image->dim.z, 1, DOUBLE );
    free( longname );
  } else {
    VT_InitImage( &imcff, "coefficients.inr", ncoeff, image->dim.z, 1, DOUBLE );
  }

  VT_AllocImage( &imcff );
  theCoeff = ((double***)(imcff.array))[0];





  

  if ( par.names.ext[0] != '\0' ) {
    int i, k;
    
    longname = malloc( strlen(par.names.ext)+ 10 );
    sprintf( longname, "%s.raw", par.names.ext );
    fd = creat( longname, S_IWUSR|S_IRUSR|S_IRGRP|S_IROTH );
    sprintf( longname, "%s.m", par.names.ext );
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
    k = strlen( par.names.ext );
    for ( i = k-1; i >= 0 && par.names.ext[i] != '/' ; i-- )
      ;
    fprintf( fp, "fid = fopen('%s.raw', 'r' );\n", &(par.names.ext[i+1]) );
    fprintf( fp, "\n" );
  }



  _Print2DSlicesHistoForMatlab( fd, fp, theHisto, image->dim.z, length, id++ );





  if ( par.z < 0 ) {

    printf( "coupe de reference = %d\n", iref );
    printf( "initialisation = " );

    switch ( par.initComputation ) {
    default :
    case _NONE_ :
      printf( "aucune\n" );
      for (i=0; i<image->dim.z; i++ ) {
	theCoeff[i][0] = 0.0;
	theCoeff[i][1] = 1.0;
      }
      break;
    case _MOMENTS_ :
      printf( "moments\n" );
      theCoeff[iref][0] = 0.0;
      theCoeff[iref][1] = 1.0;
      for ( i=iref+1; i<image->dim.z; i++ ) {
	_initAffineTrsfBetweenTwoHisto( theCoeff[i], theHisto[i], length,
					theHisto[i-1], length );
	theCoeff[i][0] = theCoeff[i][0]*theCoeff[i-1][1] + theCoeff[i-1][0];
	theCoeff[i][1] = theCoeff[i][1]*theCoeff[i-1][1];
      }
      for ( i=iref-1; i>=0; i-- ) {
	_initAffineTrsfBetweenTwoHisto( theCoeff[i], theHisto[i], length,
					theHisto[i+1], length );
	theCoeff[i][0] = theCoeff[i][0]*theCoeff[i+1][1] + theCoeff[i+1][0];
	theCoeff[i][1] = theCoeff[i][1]*theCoeff[i+1][1];
      }
    }
    
  }











  switch ( par.typeComputation ) {
  default :
    func = NULL;
    if ( 0 ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      return( 1 );
    }
    break;
  case _DIRECTE_SSD_ :    func = _DirecteSSD;   break;
  case _INVERSE_SSD_ :    func = _InverseSSD;   break;
  case _SYMETRIE_SSD_ :   func = _SymetrieSSD;   break;
  case _DIRECTE_Correlation_ :    func = _DirecteCorrelation;   break;
  case _INVERSE_Correlation_ :    func = _InverseCorrelation;   break;
  case _SYMETRIE_Correlation_ :   func = _SymetrieCorrelation;   break;
  case _DIRECTE_Vraisemblance_ :    func = _DirecteVraisemblance;   break;
  case _INVERSE_Vraisemblance_ :    func = _InverseVraisemblance;   break;
  case _SYMETRIE_Vraisemblance_ :   func = _SymetrieVraisemblance;   break;
  case _MINIMUM_Vraisemblance_ :   func = _MinimumVraisemblance;   break;
  }


  switch ( par.typeComputation ) {

  default :
  case _NONE_ :
    if ( par.initComputation == _NONE_ ) break;


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

    if ( par.z >= 0 ) {

      /* comparaison de 2 coupes
       */

      _PrintOneHistoForMatlab( fd, fp, NULL, theHisto[iref], length, id++, _AFFINE_ );
      _PrintOneHistoForMatlab( fd, fp, NULL, theHisto[par.z], length, id++, _AFFINE_ );
      
      switch ( par.initComputation ) {
      default :
      case _NONE_ :
	icoeff[0] = 0.0; icoeff[1] = 1.0;
	break;
      case _MOMENTS_ :
	_initAffineTrsfBetweenTwoHisto( icoeff, theHisto[par.z], length,
					theHisto[iref], length );
	break;
      }

      _PrintOneHistoForMatlab( fd, fp, icoeff, theHisto[par.z], length, id++, _AFFINE_ );

      printf( "compensation de %d avec %d : j = %f + %f * i\n",
	      par.z, iref, icoeff[0], icoeff[1] );

      if ( func != NULL ) {
	coeff[0] = icoeff[0];
	coeff[1] = icoeff[1];	
	_evalTrsfBetweenTwoHisto( coeff, theHisto[par.z], length,
				  theHisto[iref], length,
				  _affine, _inv_affine, 2, func, par.psigma );
	printf( "compensation de %d avec %d : j = %f + %f * i\n",
		par.z, iref, coeff[0], coeff[1] );
      }
      _PrintOneHistoForMatlab( fd, fp, coeff, theHisto[par.z], length, id++, _AFFINE_ );
      

      fprintf( fp, "\n" );
      fprintf( fp, "figure;\n" );
      fprintf( fp, "hold on;\n" );

      fprintf( fp, "plot( [0:%d], HISTO%d, 'b-' );\n", length-1, id-4 );
      fprintf( fp, "plot( [0:%d], HISTO%d, 'r--' );\n", length-1, id-3 );
      if (  par.initComputation != _NONE_ )
	fprintf( fp, "plot( %f + %f*[0:%d], HISTO%d, 'r:' );\n", icoeff[0], icoeff[1], length-1, id-2 );
      fprintf( fp, "plot( %f + %f*[0:%d], HISTO%d, 'r-' );\n", coeff[0], coeff[1], length-1, id-1 );

      fprintf( fp, "hold off;\n" );
      fprintf( fp, "\n" );

    } 

    else {
      
      /* calcul sur toute l'image
       */

      if ( func != NULL ) {

	_evalTrsfsIn3DHistos( theCoeff, theHisto, length, image->dim.z, iref,
			      _affine, _inv_affine, func, par.psigma );

      }
    }

    break;

  }








  free ( theHisto );
  if ( par.names.ext[0] != '\0' ) VT_WriteInrimage( &imcff );



  


  if ( par.z < 0 && 
       par.names.out[0] != '\0' && par.names.out[0] != '>' &&
       ( par.initComputation != _NONE_ || par.typeComputation != _NONE_ ) ) {


    switch ( image->type ) {
    default :
      VT_ErrorParse("unable to deal with such image\n", 0);
      break;
    case UCHAR :
      {
	u8 ***buf = (u8***)image->array;
	int x,y,z,j;
	double v;
	for ( z=0; z<image->dim.z; z++ ) {
	  for ( y=0; y<image->dim.y; y++ )
	  for ( x=0; x<image->dim.x; x++ ) {
	    v = theCoeff[z][0] + theCoeff[z][1] * buf[z][y][x];
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
	u16 ***buf = (u16***)image->array;
	int x,y,z,j;
	double v;
	for ( z=0; z<image->dim.z; z++ ) {
	  for ( y=0; y<image->dim.y; y++ )
	  for ( x=0; x<image->dim.x; x++ ) {
	    v = theCoeff[z][0] + theCoeff[z][1] * buf[z][y][x];
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
    }
    sprintf( image->name, "%s", par.names.out );
    VT_WriteInrimage( image );
  




    theHisto = _GetSlicesHisto( image, &intensity_max );
    for (i=0; i<image->dim.z; i++ )
      theHisto[i][0] = 0;
    length = intensity_max+1;
    _Print2DSlicesHistoForMatlab( fd, fp, theHisto, image->dim.z, length, id++ );
    free ( theHisto );
  }





  if ( par.names.ext[0] != '\0' ) {
    fprintf( fp, "\n" );
    fprintf( fp, "fclose(fid);\n" );
    fclose( fp );
    close( fd );
  }







  VT_FreeImage( image );
  VT_Free( (void**)&image );
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
	     strcmp ( argv[i], "mmnts" ) == 0 ||
	     strcmp ( argv[i], "mmts" ) == 0 ) {
	  par->initComputation = _MOMENTS_;
	} 
	else {
	  VT_ErrorParse( "unknown init mode...\n", 0 );
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
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
      }



      else if ( strcmp ( argv[i], "-sigma" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sigma...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->psigma) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sigma...\n", 0 );
      }



      else if ( strcmp ( argv[i], "-threshold" ) == 0 ||
		strcmp ( argv[i], "-thres" ) == 0 ||
		strcmp ( argv[i], "-th" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -threshold...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->threshold) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -threshold...\n", 0 );
      }





      else if ( strcmp ( argv[i], "-zref" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -zref...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->zref) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -zref...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-z" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->z) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
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
  par->threshold = 0;
  par->type = TYPE_UNKNOWN;
  par->zref = -1;
  par->z = -1;
  par->psigma = 5.0;
  par->typeComputation = _NONE_;
  par->initComputation = _NONE_;
}
