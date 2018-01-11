/*************************************************************************
 * minimum.c -
 *
 * $Id: test-tensor-voting.c,v 1.1 2000/07/20 07:56:02 greg Exp $
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
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [image-out]\n\
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






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image imres;

  int i, j;
  int dimx = 257;
  int dimy = 257;

  /*
  double a =0.003;
  double b=2.85;
  */
  double a =0.003;
  double b=2.85;

  int x = dimx/2;
  int y = dimy/2;

  double vx = 1;
  double vy = 0;

  /*
    double sigma = sqrt( 1 / a );
  */
  double sigma = 5.0;
  double c = 0.01;

  double mn2, r2, theta, ct;

  float ***theBuf;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  if ( par.names.in[0] == '\0' || par.names.in[0] == '<' )
    sprintf( par.names.in, "strengthstick-1.inr" );
  VT_InitImage( &imres, par.names.in, dimx, dimy, 1, FLOAT );
  
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  theBuf = (float***)imres.array;

  for ( j=0; j<imres.dim.y; j++ )
  for ( i=0; i<imres.dim.x; i++ )
    theBuf[0][j][i] = 0.0;

  for ( j=0; j<imres.dim.y; j++ )
  for ( i=0; i<imres.dim.x; i++ ) {
    mn2 = (i-x)*(i-x) + (j-y)*(j-y);
    if ( mn2 == 0.0 ) {
      theBuf[0][j][i] = 1.0;
      continue;
    }
    ct = fabs( ((i-x)*vx + (j-y)*vy) / sqrt( mn2 ) );
    if ( ct < sqrt(2.0)/2.0 ) continue;

    /*
      theBuf[0][j][i] = 1.0;
    */
    theta = acos( ct );
    if ( ct == 1.0 || theta == 0.0 ) { 
      theBuf[0][j][i] = exp( - a * mn2 );
      continue;
    }
    r2 = mn2 / ( 4 * ( 1.0 - ct*ct ) );
    theBuf[0][j][i] = exp( - a * 4 * theta * theta * r2 ) *
                      exp( - b / r2 );
    /*
    if ( i == 120 ) {
      printf( "(%d %d) -> mn = %f R=%f theta = %f \n",
	      i, j, sqrt(mn2), sqrt(r2), theta );
      printf( "         -> L=2*R*theta=%f  C = %f * %f \n",
	      2*theta*sqrt(r2), exp( - a * 4 * theta * theta * r2 ),  
	      exp( - b / r2 ) );
    }
    */
  }
  
  
  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }

  printf( "c is %f should be %f\n", c, b * sigma * sigma );

  sprintf( imres.name, "strengthstick-2.inr" );

  for ( j=0; j<imres.dim.y; j++ )
  for ( i=0; i<imres.dim.x; i++ )
    theBuf[0][j][i] = 0.0;

  for ( j=0; j<imres.dim.y; j++ )
  for ( j=0; j<imres.dim.y; j++ )
  for ( i=0; i<imres.dim.x; i++ ) {
    mn2 = (i-x)*(i-x) + (j-y)*(j-y);
    if ( mn2 == 0.0 ) {
      theBuf[0][j][i] = 1.0;
      continue;
    }
    ct = fabs( ((i-x)*vx + (j-y)*vy) / sqrt( mn2 ) );
    if ( ct < sqrt(2.0)/2.0 ) continue;
    theta = acos( ct );

    if ( ct == 1.0 || theta == 0.0 ) { 
      theBuf[0][j][i] = exp( - mn2 / (sigma * sigma) );
      continue;
    }
    theBuf[0][j][i] = exp( - (1/(sigma*sigma)) *
			   ( theta * theta * mn2 / (1 -ct*ct) + c*
			     4 * (1 -ct*ct) / mn2 ) );

  }



  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }


  /*--- liberations memoires ---*/
  VT_FreeImage( &imres );
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
  /*
  if ( argc == 1 ) VT_ErrorParse("\n", 0 );
  */

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
  par->type = TYPE_UNKNOWN;
}
