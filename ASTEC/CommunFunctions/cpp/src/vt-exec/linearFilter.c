/*************************************************************************
 * linearFilter.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Fri Nov 30 21:58:59 CET 2012
 *
 * ADDITIONS, CHANGES
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <chunks.h>
#include <linearFiltering.h>
#include <linearFiltering-contours.h>

#include <vt_common.h>




/* #include <recbuffer.h> */
#include <zcross.h>

typedef enum {
  _FILTER_,
  _GRADIENT_MODULUS_,
  _HESSIAN_,
  _LAPLACIAN_,
  _ZCROSSINGS_HESSIAN_,
  _ZCROSSINGS_LAPLACIAN_,
  _GRADIENT_HESSIAN_,
  _GRADIENT_LAPLACIAN_,
  _EXTREMA_GRADIENT_
} enumOutput;

typedef struct local_par {
  vt_names names;
  
  int borderLengths[3];
  typeFilteringCoefficients filter[3];

  int sliceComputation;

  ImageType type;

  enumOutput typeOutput;

  int print_param;

  int print_time;
  
} local_par;





/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static void VT_PrintParam( FILE *theFile, local_par *par );
static double _GetTime();
static double _GetClock();
static char *_BaseName( char *p );




static char *usage = "[image-in] [image-out] [-x %d] [-y %d] [-z %d]\n\
 [-alpha %lf [%lf [%lf]]] [-sigma %lf [%lf [%lf]]] \n\
 [-border|cont %d [%d [%d]]] [-2D]\n\
 [-gaussian-type|-type [deriche|fidrich|young-1995|young-2002|gabor-young-2002|convolution]]\n\
 [-filter | -gradient | -gradient-modulus |-hessian | -laplacian | ...\n\
 ... | -zero-crossings-hessian | -zero-crossings-laplacian | ...\n\
 ... | -gradient-hessian | -gradient-laplacian | -gradient-extrema]\n\
 [-neg | -pos]\n\
 [-parallel|-no-parallel] [-max-chunks %d]\n\
 [-parallel-scheduling|-ps default|static|dynamic-one|dynamic|guided]\n\
 [-time] [-notime] [-param]\n\
 [-inv] [-swap] [-v] [-nv] [-D] [-help] [-h] [options-de-type]";

static char *detail = "\
 si 'image-in' est '-', on prendra stdin\n\
 si 'image-out' est absent, on prendra stdout\n\
 si les deux sont absents, on prendra stdin et stdout\n\
 -x | -y | -z : ordre de derivation selon X, Y ou Z\n\
 -alpha | -a  : alpha pour le filtre recursif de Deriche\n\
 -sigma       : sigma pour la gaussienne (approximation recursive ou convolution)\n\
 -border | -cont : points ajoutes aux bords\n\
-filter-type : \n\
 -inv : inverse 'image-in'\n\
 -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
 -v : mode verbose\n\
 -D : mode debug\n\
 options-de-type : -o 1    : unsigned char\n\
                   -o 2    : unsigned short int\n\
                   -o 2 -s : short int\n\
                   -o 4 -s : int\n\
                   -r      : float\n\
 si aucune de ces options n'est presente, on prend le type de 'image-in'\n";





static char program[STRINGLENGTH];

int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  int theDims[3];

  double time_init = _GetTime();
  double time_exit;
  double clock_init = _GetClock();
  double clock_exit;


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  if ( par.type == TYPE_UNKNOWN ) {
    switch ( par.typeOutput ) {
    default :
    case _FILTER_ :             par.type = FLOAT; break;
    case _GRADIENT_MODULUS_ :   par.type = FLOAT; break;
    case _HESSIAN_ :            par.type = FLOAT; break;
    case _LAPLACIAN_ :          par.type = FLOAT; break;
    case _ZCROSSINGS_HESSIAN_ :   par.type = UCHAR; break;
    case _ZCROSSINGS_LAPLACIAN_ : par.type = UCHAR; break;
    case _GRADIENT_HESSIAN_ :   par.type = FLOAT; break;
    case _GRADIENT_LAPLACIAN_ : par.type = FLOAT; break;
    case _EXTREMA_GRADIENT_ :   par.type = FLOAT; break;
    }
  }

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  if ( _VT_DEBUG_ || par.print_param )
    VT_PrintParam( stdout, &par );


  /*--- operations eventuelles sur l'image d'entree ---*/
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  VT_InitFromImage( &imres, image, par.names.out, image->type );
  if ( par.type != TYPE_UNKNOWN ) imres.type = par.type;
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  /*--- calcul ---*/
  theDims[0] = image->dim.x;
  theDims[1] = image->dim.y;
  theDims[2] = image->dim.z;

  switch ( par.typeOutput ) {

  default :
    VT_FreeImage( &imres ); 
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("such option not implemented yet\n", 0);
    exit( 1 );

  case _FILTER_ :
    if ( separableLinearFiltering( image->buf, image->type, 
				   imres.buf, imres.type, theDims,
				   par.borderLengths, par.filter ) != 1 ) {
      VT_FreeImage( &imres ); 
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to filter input image\n", 0);
      exit( 1 );
    }
    break;

  case _GRADIENT_MODULUS_ :
    if ( par.sliceComputation ) {
      if ( gradientModulus2D( image->buf, image->type, 
			      imres.buf, imres.type, theDims,
			      par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D gradient modulus \n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientModulus( image->buf, image->type, 
			    imres.buf, imres.type, theDims,
			    par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient modulus\n", 0);
	exit( 1 );
      }
    }
    break;

  case _HESSIAN_ :
    if ( par.sliceComputation ) {
      if ( gradientHessianGradient2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D gradient.Hessian.gradient\n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientHessianGradient( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient.Hessian.gradient\n", 0);
	exit( 1 );
      }
    }
    break;

  case _LAPLACIAN_ :
    if ( par.sliceComputation ) {
      if ( laplacian2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D laplacian \n", 0);
	exit( 1 );
      }
    }
    else {
      if ( laplacian( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute laplacian\n", 0);
	exit( 1 );
      }
    }
    break;

  case _ZCROSSINGS_HESSIAN_ :
    if ( par.sliceComputation ) {
      if ( gradientHessianGradientZeroCrossings2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D zero-crossings of gradient.Hessian.gradient \n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientHessianGradientZeroCrossings( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute zero-crossings of gradient.Hessian.gradient\n", 0);
	exit( 1 );
      }
    }
    break;

  case _ZCROSSINGS_LAPLACIAN_ :
    if ( par.sliceComputation ) {
      if ( laplacianZeroCrossings2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D zero-crossings of laplacian \n", 0);
	exit( 1 );
      }
    }
    else {
      if ( laplacianZeroCrossings( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute zero-crossings of laplacian\n", 0);
	exit( 1 );
      }
    }
    break;

  case _GRADIENT_HESSIAN_ :
    if ( par.sliceComputation ) {
      if ( gradientOnGradientHessianGradientZC2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient on 2D zero-crossings of gradient.Hessian.gradient \n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientOnGradientHessianGradientZC( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient on zero-crossings of gradient.Hessian.gradient\n", 0);
	exit( 1 );
      }
    }
    break;

  case _GRADIENT_LAPLACIAN_ :
    if ( par.sliceComputation ) {
      if ( gradientOnLaplacianZC2D( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient on 2D zero-crossings of laplacian\n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientOnLaplacianZC( image->buf, image->type, 
				     imres.buf, imres.type, theDims,
				     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute gradient on zero-crossings of laplacian\n", 0);
	exit( 1 );
      }
    }
    break;

  case _EXTREMA_GRADIENT_ :
    if ( par.sliceComputation ) {
      if ( gradientMaxima2D( image->buf, image->type, 
			     imres.buf, imres.type, theDims,
			     par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute 2D extrema of the gradient\n", 0);
	exit( 1 );
      }
    }
    else {
      if ( gradientMaxima( image->buf, image->type, 
			   imres.buf, imres.type, theDims,
			   par.borderLengths, par.filter ) != 1 ) {
	VT_FreeImage( &imres ); 
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to compute  extrema of the gradient\n", 0);
	exit( 1 );
      }
    }
    break;
  }

  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image );
  VT_FreeImage( &imres );
  VT_Free( (void**)&image );



  time_exit = _GetTime();
  clock_exit = _GetClock();

  if ( par.print_time ) { 
    fprintf( stderr, "%s: elapsed (real) time = %f\n", _BaseName( program ), time_exit - time_init );
    fprintf( stderr, "\t       elapsed (user) time = %f (processors)\n", clock_exit - clock_init );
    fprintf( stderr, "\t       ratio (user)/(real) = %f\n", (clock_exit - clock_init)/(time_exit - time_init) );
  }


  return( 1 );
}









static void VT_Parse( int argc, char *argv[], local_par *par )
{
  int i, j, nb, status;
  int o=0, s=0, r=0;
  int bool_contours = 0;
  int tmp;
  char text[STRINGLENGTH];
  
  if ( VT_CopyName( _BaseName( program ), argv[0] ) != 1 )
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
      else if ( strcmp ( argv[i], "-h" ) == 0 && argv[i][2] == '\0' ) {
	VT_ErrorParse("\n", 0);
      }
      else if ( strcmp ( argv[i], "-v" ) == 0 && argv[i][2] == '\0' ) {
	if ( _VT_VERBOSE_ <= 0 ) _VT_VERBOSE_ = 1;
	else                     _VT_VERBOSE_ ++;
      }
      else if ( strcmp ( argv[i], "-nv" ) == 0 ) {
	_VT_VERBOSE_ = 0;
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

      /*--- ordres de derivation ---*/
      else if ( strcmp ( argv[i], "-x" ) == 0 && argv[i][2] == '\0' ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -x...\n", 0 );
	status = sscanf( argv[i],"%d",&tmp );
	if ( status <= 0 ) VT_ErrorParse( "parsing -x...\n", 0 );
	switch ( tmp ) {
	default :
	  par->filter[0].derivative = NODERIVATIVE;   break;
	case 0 :
	  par->filter[0].derivative = DERIVATIVE_0;   break;
	case 1 :
	  par->filter[0].derivative = DERIVATIVE_1;   break;
	case 2 :
	  par->filter[0].derivative = DERIVATIVE_2;   break;
	case 3 :
	  par->filter[0].derivative = DERIVATIVE_3;   break;
	}
      }

      else if ( strcmp ( argv[i], "-y" ) == 0 && argv[i][2] == '\0' ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -y...\n", 0 );
	status = sscanf( argv[i],"%d",&tmp );
	if ( status <= 0 ) VT_ErrorParse( "parsing -y...\n", 0 );
	switch ( tmp ) {
	default :
	  par->filter[1].derivative = NODERIVATIVE;   break;
	case 0 :
	  par->filter[1].derivative = DERIVATIVE_0;   break;
	case 1 :
	  par->filter[1].derivative = DERIVATIVE_1;   break;
	case 2 :
	  par->filter[1].derivative = DERIVATIVE_2;   break;
	case 3 :
	  par->filter[1].derivative = DERIVATIVE_3;   break;
	}
      }

      else if ( strcmp ( argv[i], "-z" ) == 0 && argv[i][2] == '\0' ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -z...\n", 0 );
	status = sscanf( argv[i],"%d",&tmp );
	if ( status <= 0 ) VT_ErrorParse( "parsing -z...\n", 0 );
	switch ( tmp ) {
	default :
	  par->filter[2].derivative = NODERIVATIVE;   break;
	case 0 :
	  par->filter[2].derivative = DERIVATIVE_0;   break;
	case 1 :
	  par->filter[2].derivative = DERIVATIVE_1;   break;
	case 2 :
	  par->filter[2].derivative = DERIVATIVE_2;   break;
	case 3 :
	  par->filter[2].derivative = DERIVATIVE_3;   break;
	}
      }

      /*--- contours ? ---*/
      else if ( strcmp ( argv[i], "-edges" ) == 0 ) {
	bool_contours = 1;
      }

      /*--- alpha ---*/
      else if ( (strcmp ( argv[i], "-alpha" ) == 0) || 
		(strcmp ( argv[i], "-a" ) == 0 && argv[i][2] == '\0') ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing [-alpha|-a]...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->filter[0].coefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing [-alpha|-a]...\n", 0 );
	i ++;
	if ( i >= argc) {
	  par->filter[1].coefficient = par->filter[0].coefficient;
	  par->filter[2].coefficient = par->filter[0].coefficient;
	}
	else {
	  status = sscanf( argv[i], "%lf", &(par->filter[1].coefficient) );
	  if ( status <= 0 ) {
	    i--;
	    par->filter[1].coefficient = par->filter[0].coefficient;
	    par->filter[2].coefficient = par->filter[0].coefficient;
	  }
	  else {
	    i ++;
	    if ( i >= argc) par->filter[2].coefficient = 0.0;
	    else {
	      status = sscanf( argv[i], "%lf", &(par->filter[2].coefficient) );
	      if ( status <= 0 ) {
		i--;
		par->filter[2].coefficient = 0;
	      }
	    }
	  }
	}
	par->filter[0].type = ALPHA_DERICHE;
	par->filter[1].type = ALPHA_DERICHE;
	par->filter[2].type = ALPHA_DERICHE;
      }
      
      /*--- sigma ---*/
      else if ( (strcmp ( argv[i], "-sigma" ) == 0) ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -sigma ...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->filter[0].coefficient) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -sigma ...\n", 0 );
	i ++;
	if ( i >= argc) {
	  par->filter[1].coefficient = par->filter[0].coefficient;
	  par->filter[2].coefficient = par->filter[0].coefficient;
	}
	else {
	  status = sscanf( argv[i], "%lf", &(par->filter[1].coefficient) );
	  if ( status <= 0 ) {
	    i--;
	    par->filter[1].coefficient = par->filter[0].coefficient;
	    par->filter[2].coefficient = par->filter[0].coefficient;
	  }
	  else {
	    i ++;
	    if ( i >= argc) par->filter[2].coefficient = 0.0;
	    else {
	      status = sscanf( argv[i], "%lf", &(par->filter[2].coefficient) );
	      if ( status <= 0 ) {
		i--;
		par->filter[2].coefficient = 0;
	      }
	    }
	  }
	}
	for ( j=0; j<3; j++ ) {
	  switch ( par->filter[j].type ) {
	  case GAUSSIAN_DERICHE :
	  case GAUSSIAN_FIDRICH :
	  case GAUSSIAN_YOUNG_1995 :
	  case GAUSSIAN_YOUNG_2002 :
	  case GABOR_YOUNG_2002 :
	  case GAUSSIAN_CONVOLUTION :
	    break;
	  default :
	    par->filter[j].type = GAUSSIAN_CONVOLUTION;
	    break;
	  }
	}
      }

      
      /*--- bordure ---*/
      else if ( strcmp ( argv[i], "-cont" ) == 0 
		|| strcmp ( argv[i], "-border" ) == 0) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -cont...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->borderLengths[0]) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -cont...\n", 0 );
	i ++;
	if ( i >= argc) {
	  par->borderLengths[2] = par->borderLengths[1] = par->borderLengths[0];
	}
	else {
	  status = sscanf( argv[i], "%d", &(par->borderLengths[1]) );
	  if ( status <= 0 ) {
	    i--;
	    par->borderLengths[2] = par->borderLengths[1] = par->borderLengths[0];
	  }
	  else {
	    i ++;
	    if ( i >= argc) par->borderLengths[2] = 0;
	    else {
	      status = sscanf( argv[i], "%d", &(par->borderLengths[2]) );
	      if ( status <= 0 ) {
		i--;
		par->borderLengths[2] = 0;
	      }
	    }
	  }
	}
      }


      /*--- filter type ---*/
      else if ( strcmp ( argv[i], "-type" ) == 0 
		|| strcmp ( argv[i], "-gaussian-type" ) == 0 ) {
	i++;
	if ( strcmp ( argv[i], "deriche" ) == 0 ) {
	  par->filter[0].type = GAUSSIAN_DERICHE;
	  par->filter[1].type = GAUSSIAN_DERICHE;
	  par->filter[2].type = GAUSSIAN_DERICHE;
	}
	else if ( strcmp ( argv[i], "fidrich" ) == 0 ) {
	  par->filter[0].type = GAUSSIAN_FIDRICH;
	  par->filter[1].type = GAUSSIAN_FIDRICH;
	  par->filter[2].type = GAUSSIAN_FIDRICH;
	}
	else if ( strcmp ( argv[i], "young-1995" ) == 0 ) {
	  par->filter[0].type = GAUSSIAN_YOUNG_1995;
	  par->filter[1].type = GAUSSIAN_YOUNG_1995;
	  par->filter[2].type = GAUSSIAN_YOUNG_1995;
	}
	else if ( strcmp ( argv[i], "young-2002" ) == 0 ) {
	  par->filter[0].type = GAUSSIAN_YOUNG_2002;
	  par->filter[1].type = GAUSSIAN_YOUNG_2002;
	  par->filter[2].type = GAUSSIAN_YOUNG_2002;
	}
	else if ( strcmp ( argv[i], "gabor-young-2002" ) == 0 ) {
	  par->filter[0].type = GABOR_YOUNG_2002;
	  par->filter[1].type = GABOR_YOUNG_2002;
	  par->filter[2].type = GABOR_YOUNG_2002;
	}
	else if ( strcmp ( argv[i], "convolution" ) == 0  ) {
	  par->filter[0].type = GAUSSIAN_CONVOLUTION;
	  par->filter[1].type = GAUSSIAN_CONVOLUTION;
	  par->filter[2].type = GAUSSIAN_CONVOLUTION;
	}
	else {
	  VT_ErrorParse( "parsing -gaussian-type...\n", 0 );
	}
      }

      else if ( strcmp ( argv[i], "-filter" ) == 0 ) {
	par->typeOutput = _FILTER_;
      }
      else if ( (strcmp ( argv[i], "-gradient" ) == 0 && argv[i][9] == '\0')
		|| strcmp ( argv[i], "-gradient-modulus" ) == 0 ) {
	par->typeOutput = _GRADIENT_MODULUS_;
      }
      else if ( strcmp ( argv[i], "-hessian" ) == 0 ) {
	par->typeOutput = _HESSIAN_;
      }
      else if ( strcmp ( argv[i], "-laplacian" ) == 0 ) {
	par->typeOutput = _LAPLACIAN_;
      }
      else if ( strcmp ( argv[i], "-zero-crossings-hessian" ) == 0 
		|| strcmp ( argv[i], "-zcrossings-hessian" ) == 0 ) {
	par->typeOutput = _ZCROSSINGS_HESSIAN_;
      }
      else if ( strcmp ( argv[i], "-zero-crossings-laplacian" ) == 0 
		|| strcmp ( argv[i], "-zcrossings-laplacian" ) == 0 ) {
	par->typeOutput = _ZCROSSINGS_LAPLACIAN_;
      }
      else if ( strcmp ( argv[i], "-gradient-hessian" ) == 0 ) {
	par->typeOutput = _GRADIENT_HESSIAN_;
      }
      else if ( strcmp ( argv[i], "-gradient-laplacian" ) == 0 ) {
	par->typeOutput = _GRADIENT_LAPLACIAN_;
      }
      else if ( strcmp ( argv[i], "-gradient-extrema" ) == 0 ) {
	par->typeOutput = _EXTREMA_GRADIENT_;
      }


      else if ( strcmp ( argv[i], "-neg" ) == 0 ) {
	ZeroCrossings_Are_Negative();
      }
      else if ( strcmp ( argv[i], "-pos" ) == 0 ) {
	ZeroCrossings_Are_Positive();
      }

      /* 2D
       */
      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	par->sliceComputation = 1;
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

      /* parallelism
       */
      else if ( strcmp ( argv[i], "-parallel" ) == 0 ) {
	setMaxChunks( 100 );
      }
      
      else if ( strcmp ( argv[i], "-no-parallel" ) == 0 ) {
	setMaxChunks( 1 );
      }
      
      else if ( strcmp ( argv[i], "-max-chunks" ) == 0 ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-max-chunks", 0 );
	status = sscanf( argv[i], "%d", &tmp );
	if ( status <= 0 ) VT_ErrorParse( "-max-chunks", 0 );
	if ( tmp >= 1 ) setMaxChunks( tmp );
      }
      
      else if ( strcmp ( argv[i], "-parallel-scheduling" ) == 0 || 
		( strcmp ( argv[i], "-ps" ) == 0 && argv[i][3] == '\0') ) {
	i ++;
	if ( i >= argc)    VT_ErrorParse( "-parallel-scheduling", 0 );
	if ( strcmp ( argv[i], "default" ) == 0 ) {
	  setOpenMPScheduling( _DEFAULT_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "static" ) == 0 ) {
	  setOpenMPScheduling( _STATIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic-one" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_ONE_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "dynamic" ) == 0 ) {
	  setOpenMPScheduling( _DYNAMIC_SCHEDULING_ );
	}
	else if ( strcmp ( argv[i], "guided" ) == 0 ) {
	  setOpenMPScheduling( _GUIDED_SCHEDULING_ );
	}
	else {
	  fprintf( stderr, "unknown scheduling type: '%s'\n", argv[i] );
	  VT_ErrorParse( "-parallel-scheduling", 0 );
	}
      }

      else if ( (strcmp ( argv[i], "-param" ) == 0 && argv[i][6] == '\0') ) {
	par->print_param = 1;
      }

      else if ( (strcmp ( argv[i], "-time" ) == 0 && argv[i][5] == '\0') ) {
	par->print_time = 1;
      }
      else if ( (strcmp ( argv[i], "-notime" ) == 0 && argv[i][7] == '\0') 
		|| (strcmp ( argv[i], "-no-time" ) == 0 && argv[i][8] == '\0') ) {
	par->print_time = 0;
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
  
  /*--- ordres de derivation ---*/
  if ( bool_contours == 1 ) {
    if ( par->filter[0].derivative == DERIVATIVE_1 )
         par->filter[0].derivative =  DERIVATIVE_1_CONTOURS;
    if ( par->filter[1].derivative == DERIVATIVE_1 )
         par->filter[1].derivative =  DERIVATIVE_1_CONTOURS;
    if ( par->filter[2].derivative == DERIVATIVE_1 )
         par->filter[2].derivative =  DERIVATIVE_1_CONTOURS;
  }

    
  /*--- type de l'image resultat ---*/
  if ( (o == 1) && (s == 1) && (r == 0) ) par->type = SCHAR;
  if ( (o == 1) && (s == 0) && (r == 0) ) par->type = UCHAR;
  if ( (o == 2) && (s == 0) && (r == 0) ) par->type = USHORT;
  if ( (o == 2) && (s == 1) && (r == 0) ) par->type = SSHORT;
  if ( (o == 4) && (s == 1) && (r == 0) ) par->type = SINT;
  if ( (o == 0) && (s == 0) && (r == 1) ) par->type = FLOAT;
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
    VT_Names( &(par->names) );
    
    par->borderLengths[0] = 0;
    par->borderLengths[1] = 0;
    par->borderLengths[2] = 0;
    initFilteringCoefficients( &(par->filter[0]) );
    initFilteringCoefficients( &(par->filter[1]) );
    initFilteringCoefficients( &(par->filter[2]) );

    par->sliceComputation = 0;

    par->type = TYPE_UNKNOWN;

    par->typeOutput = _FILTER_;

    par->print_param = 0;

    par->print_time = 0;
}





static void VT_PrintParam( FILE *theFile, local_par *par )
{
  FILE *f = theFile;
  if ( theFile == (FILE*)NULL ) f = stderr;

  fprintf( f, "==================================================\n" );
  fprintf( f, "%s: parameters\n", _BaseName( program ) );

  VT_PrintNames( f, &(par->names) );
  
  fprintf( f, "- borders = [%d %d %d]\n",  par->borderLengths[0],
	   par->borderLengths[1], par->borderLengths[2] );
  
  printFilteringCoefficients( f, &(par->filter[0]), "filter along X" );
  printFilteringCoefficients( f, &(par->filter[1]), "filter along Y" );
  printFilteringCoefficients( f, &(par->filter[2]), "filter along Z" );

  fprintf( f, "- computation: " );
  switch ( par->typeOutput ) {
  default : fprintf( f, "unknown\n"); break;
  case _FILTER_ : fprintf( f, "_FILTER_\n"); break;
  case _GRADIENT_MODULUS_ : fprintf( f, "_GRADIENT_MODULUS_\n"); break;
  case _HESSIAN_ : fprintf( f, "_HESSIAN_\n"); break;
  case _LAPLACIAN_ : fprintf( f, "_LAPLACIAN_\n"); break;
  case _ZCROSSINGS_HESSIAN_ : fprintf( f, "_ZCROSSINGS_HESSIAN_\n"); break;
  case _ZCROSSINGS_LAPLACIAN_ : fprintf( f, "_ZCROSSINGS_LAPLACIAN_\n"); break;
  case _GRADIENT_HESSIAN_ : fprintf( f, "_GRADIENT_HESSIAN_\n"); break;
  case _GRADIENT_LAPLACIAN_ : fprintf( f, "_GRADIENT_LAPLACIAN_\n"); break;
  case _EXTREMA_GRADIENT_ : fprintf( f, "_EXTREMA_GRADIENT_\n"); break;
  }
  
  fprintf( f, "- output image type: " );
  switch ( par->type ) {
  default :     fprintf( f, "TYPE_UNKNOWN\n" ); break;
  case SCHAR :  fprintf( f, "SCHAR\n" ); break;
  case UCHAR :  fprintf( f, "UCHAR\n" ); break;
  case SSHORT : fprintf( f, "SSHORT\n" ); break;
  case USHORT : fprintf( f, "USHORT\n" ); break; 
  case UINT :   fprintf( f, "UINT\n" ); break; 
  case SINT :    fprintf( f, "INT\n" ); break;
  case ULINT :  fprintf( f, "ULINT\n" ); break;
  case FLOAT :  fprintf( f, "FLOAT\n" ); break;
  case DOUBLE : fprintf( f, "DOUBLE\n" ); break;
  }  

  fprintf( f, "- print time = %d\n",  par->print_time );

  fprintf( f, "==================================================\n" );
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
