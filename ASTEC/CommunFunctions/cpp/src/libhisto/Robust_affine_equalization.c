/*************************************************************************
 * minimum.c -
 *
 * $Id: Robust_affine_equalization.c,v 1.1 2003/09/19 16:57:18 greg Exp $
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

typedef struct local_par {
  vt_names names;
  int type;
  int threshold1;
  int threshold2;
  int max_iterations;
  double a_err;
  double fraction;
  double a;
  double b;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-1] [image-2] [image-out] [-a %lf] [-b %lf]\n\
\t [-t1 %d] [-t2 %d] [-a-err %lf] [-frac %lf] [-max %d] [-help] [-v]";

/* "\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]"; */

static char *detail = "\
 'image-out' will be 'image-1' after intensity transformation\n\
    with an affine function: 'A' * intensity + 'B' \n\
    with A = a_arg * a_cmp and B = a_arg * b_cmp + b_arg\n\
    (a_cmp,b_cmp) are computed here\n\
    (a_arg,b_arg) are passed as arguments\n\
    typically they have been used to resample image-2\n\
 -t1 %d      # threshold for 'image-1': defines admissble points\n\
    (only points with values in 'image-1' >= threshold are keeped)\n\
 -t2 %d      # threshold for 'image-2': defines admissble points\n\
    (only points with values in 'image-2' >= threshold are keeped)\n\
 -a-err %lf  # threshold on the difference between two\n\
    successive values of 'a'\n\
 -frac %lf   # fraction of admissible points that are kept for\n\
    the robust computation of 'a' and 'b'\n\
 -max %d     # maximum number of iterations\n\
\n\
 $Revision: 1.1 $ $Date: 2003/09/19 16:57:18 $ $Author: greg $\n";

static char program[STRINGLENGTH];


typedef struct {
  float val1;
  float val2;
  double error;
} typePoint;

 int errorcompare( const void *arg1, const void *arg2 )
{
  typePoint *pt1 = (typePoint *)arg1;
  typePoint *pt2 = (typePoint *)arg2;

  if ((*pt1).error > (*pt2).error)
    return (1);
  if ((*pt1).error < (*pt2).error)
    return (-1);
  return (0);
}



int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image1, *image2;
  int i;

  int n, npts, nkeep;
  typePoint *pts = NULL;

  double moy_1, moy_2, moy_12, moy_22, moy_11;
  double A, a, b, old_a, old_b;
  double a_err;
  
  int iteration;

  double val;

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
  

  if ( image1->dim.x != image2->dim.x
       || image1->dim.y != image2->dim.y
       || image1->dim.z != image2->dim.z ) {
    VT_ErrorParse( "images should have the same dimension\n", 0 );
  }

  pts = (typePoint *)malloc( (image1->dim.x*image1->dim.y*image1->dim.z)*sizeof(typePoint) );
  if ( pts == NULL ) 
    VT_ErrorParse( "unable to allocate points list\n", 0 );

  npts = 0;

  switch( image1->type ) {
  default :
    fprintf( stderr, "such image type ('%s') not handled yet\n", image1->name );
    VT_ErrorParse( "", 0 );

  case USHORT :
    {
      unsigned short int *buf1 = (unsigned short int *)image1->buf;

      switch( image2->type ) {
      default :
	fprintf( stderr, "such image type ('%s') not handled yet\n", image2->name );
	VT_ErrorParse( "", 0 );
      case UCHAR :
	{
	  unsigned char *buf2 = (unsigned char *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
      case USHORT :
	{
	  unsigned short int *buf2 = (unsigned short int *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
       case SSHORT :
	{
	  short int *buf2 = (short int *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
       }

    }
    break;

  case SSHORT :
    {
      short int *buf1 = (short int *)image1->buf;

      switch( image2->type ) {
      default :
	fprintf( stderr, "such image type ('%s') not handled yet\n", image2->name );
	VT_ErrorParse( "", 0 );
      case UCHAR :
	{
	  unsigned char *buf2 = (unsigned char *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
      case USHORT :
	{
	  unsigned short int *buf2 = (unsigned short int *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
       case SSHORT :
	{
	  short int *buf2 = (short int *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
       }

    }
    break;

  case FLOAT :
    {
      float *buf1 = (float *)image1->buf;

      switch( image2->type ) {
      default :
	fprintf( stderr, "such image type ('%s') not handled yet\n", image2->name );
	VT_ErrorParse( "", 0 );
      case FLOAT :
	{
	  float *buf2 = (float *)image2->buf;
	  for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	    if ( buf1[i] >= par.threshold1 && buf2[i] >= par.threshold2 ) {
	      pts[npts].val1 = buf1[i];
	      pts[npts].val2 = buf2[i];
	      npts ++;
	    }
	  }
	}
	break;
      }
    }
    break;
  }


  if ( _VT_VERBOSE_ ) {
    printf( "Keep %d points among %lu\n", npts, image1->dim.x*image1->dim.y*image1->dim.z );
    printf( "   thresholds were %d for '%s' and %d for '%s'\n", 
	    par.threshold1, image1->name, par.threshold2, image2->name );
    if ( npts <= 0 ) VT_ErrorParse( "no selected points\n", 0 );
  }


  moy_1  = 0.0;
  moy_2  = 0.0;
  moy_12 = 0.0;
  moy_11 = 0.0;
  moy_22 = 0.0;

  for ( n=0; n<npts; n++ ) {
    moy_1  += pts[n].val1;
    moy_11 += pts[n].val1 * pts[n].val1;
    moy_2  += pts[n].val2;
    moy_22 += pts[n].val2 * pts[n].val2;
    moy_12 += pts[n].val1 * pts[n].val2;
  }

  moy_1  /= (double)npts;
  moy_2  /= (double)npts;
  moy_12 /= (double)npts;
  moy_11 /= (double)npts;
  moy_22 /= (double)npts;


  /* Calcul de A */
  A = 0.5 * 
    ( (moy_22 - moy_2*moy_2) - (moy_11 - moy_1*moy_1) ) /
    ( moy_1*moy_2 - moy_12 );
  
  /* Calcul du coeff dir a */
  if (moy_12 - moy_1*moy_2 > 0)
    a = -A + sqrt(A*A+1);
  else
    a = -A - sqrt(A*A+1);
  
  /* Calcul de l'ordonnee a l'origine */
  b = moy_2 - a*moy_1;
  /* b = 0.0; */
  
  if ( _VT_VERBOSE_ ) {
    printf("Initial affine bias: y = %f * x + %f\n", a, b);
  }

  nkeep = (int)(par.fraction * npts + 0.5);
  
  if ( nkeep <= 0 )
    VT_ErrorParse( "no points are kept", 0 );
  if ( nkeep > npts ) nkeep = npts;

  if ( _VT_VERBOSE_ ) {
    printf( "   will further use %d points among %d\n", nkeep, npts );
  }


  iteration = 0;
  if ( iteration < par.max_iterations ) {
    do {
      
      iteration++;
      old_a = a;
      old_b = b;
      
      for ( n=0; n<npts;n++ ) 
	pts[n].error = fabs( pts[n].val2 - (a * pts[n].val1 + b) ) / ( 1 + b*b );
      
      qsort( pts, (size_t)npts, (size_t)sizeof(typePoint), errorcompare );
      
      moy_1  = 0.0;
      moy_2  = 0.0;
      moy_12 = 0.0;
      moy_11 = 0.0;
      moy_22 = 0.0;
      
      for ( n=0; n<nkeep; n++ ) {
	moy_1  += pts[n].val1;
	moy_11 += pts[n].val1 * pts[n].val1;
	moy_2  += pts[n].val2;
	moy_22 += pts[n].val2 * pts[n].val2;
	moy_12 += pts[n].val1 * pts[n].val2;
      }
      
      moy_1  /= (double)nkeep;
      moy_2  /= (double)nkeep;
      moy_12 /= (double)nkeep;
      moy_11 /= (double)nkeep;
      moy_22 /= (double)nkeep;
      
      /* Calcul de A */
      A = 0.5 * 
	( (moy_22 - moy_2*moy_2) - (moy_11 - moy_1*moy_1) ) /
	( moy_1*moy_2 - moy_12 );
      
      /* Calcul du coeff dir a */
      if (moy_12 - moy_1*moy_2 > 0)
	a = -A + sqrt(A*A+1);
      else
	a = -A - sqrt(A*A+1);
      
      /* Calcul de l'ordonnee a l'origine */
      b = moy_2 - a*moy_1;
      /* b = 0.0; */
      
      if ( old_a > a ) a_err = old_a - a;
      else             a_err = a - old_a;
      
      if ( _VT_VERBOSE_ ) {
	printf("Iteration #%d: y = %f * x + %f (error on 'a' is %g)\n", iteration, a, b, a_err );
      }
      
    } while( iteration < par.max_iterations && ( a_err > par.a_err ) );
  }





  if ( _VT_VERBOSE_ ) {
    printf( "Transform '%s' into '%s'\n", image1->name, par.names.out );
    printf( "Computed coefficients A=%lf B=%lf\n", a, b );
  }

  b = par.a * b + par.b;
  a = par.a *a;

  if ( _VT_VERBOSE_ ) {
    printf( "Coefficients for interpolation A=%lf B=%lf\n", a, b );
  }
  
    

  switch( image1->type ) {
  default :
    fprintf( stderr, "such image type ('%s') not handled yet\n", image1->name );
    VT_ErrorParse( "", 0 );

  case UCHAR :
    {
      unsigned char *buf1 = (unsigned char *)image1->buf;
      for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	val = a*buf1[i]+b;
	if ( val > 255.0 )  buf1[i] = 255;
	else if ( val > 0.0 ) buf1[i] = (unsigned char)(val+0.5);
	else                  buf1[i] = 0;
      }
    }
    break;
    
  case USHORT :
    {
      unsigned short int *buf1 = (unsigned short int *)image1->buf;
      for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	val = a*buf1[i]+b;
	if ( val > 65535.0 )  buf1[i] = 65535;
	else if ( val > 0.0 ) buf1[i] = (unsigned short int)(val+0.5);
	else                  buf1[i] = 0;
      }
    }
    break;
    
  case SSHORT :
    {
      short int *buf1 = (short int *)image1->buf;
      for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	val = a*buf1[i]+b;
	if ( val > 32767.0 )  buf1[i] = 32767;
	else if ( val > 0.0 ) buf1[i] = (unsigned short int)(val+0.5);
	else if ( val >  -32768.0 ) buf1[i] = (unsigned short int)(val-0.5);
	else                        buf1[i] = -32768;
      }
    }
    break;

   case FLOAT :
    {
      float *buf1 = (float *)image1->buf;
      for ( i=0; i<image1->dim.x*image1->dim.y*image1->dim.z; i++ ) {
	val = a*buf1[i]+b;
	buf1[i] = val;
      }
    }
    break;
 }


  sprintf( image1->name, "%s", par.names.out );
  VT_WriteInrimage( image1 );

  
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


      else if ( strcmp ( argv[i], "-t1" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -t1...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->threshold1) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -t1...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-t2" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -t2...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->threshold2) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -t2...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-a-err" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -a-err...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->a_err) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -a-err...\n", 0 );
      }
      else if ( strcmp ( argv[i], "-frac" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -frac...\n", 0 );
	status = sscanf( argv[i],"%lf",&(par->fraction) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -frac...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-max" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -max...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->max_iterations) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -max...\n", 0 );
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
  par->type = TYPE_UNKNOWN;
  par->threshold1 = 0;
  par->threshold2 = 0;

  par->max_iterations = 100;

  par->a_err = 0.0001;
  par->fraction = 0.75;

  par->a = 1.0;
  par->b = 0.0;

}
