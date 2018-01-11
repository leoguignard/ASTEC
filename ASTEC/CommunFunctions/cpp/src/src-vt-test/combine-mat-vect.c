/*************************************************************************
 * combine-mat-vect.c -
 *
 * $Id: combine-mat-vect.c,v 1.2 2000/08/16 15:22:19 greg Exp $
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
  char matrice[STRINGLENGTH];
  int type;
  int inverse;
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
static int  ReadMatrice( char *name, double *mat );
static int  InverseMat4x4( double *matrice, double *inv );




static char *usage = "[image-in] [image-out] [-mat %s]\n\
\t [-inv] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t -inv : inverse 'matrix'\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\t options-de-type : -o 1    : unsigned char\n\
\t                   -o 2    : unsigned short int\n\
\t                   -o 2 -s : short int\n\
\t                   -o 4 -s : int\n\
\t                   -r      : float\n\
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
\t $Revision: 1.2 $ $Date: 2000/08/16 15:22:19 $ $Author: greg $\n";


static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *imx = (vt_image*)NULL;
  vt_image *imy = (vt_image*)NULL;
  vt_image *imz = (vt_image*)NULL;
  char name[STRINGLENGTH];

  double mat[16];
  int i, j, k;
  double x, y, z;


  /*  0  1  2  3
      4  5  6  7
      8  9 10 11
     12 13 14 15 
     */
  for (i=0;i<16;i++) mat[i] = 0.0;
  mat[0] = mat[5] = mat[10] = mat[15] = 1.0;
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  
  /*--- lecture de l'image d'entree ---*/
  sprintf( name, "%s.x", par.names.in );
  imx = _VT_Inrimage( name );
  if ( imx == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read X input image\n", 0);

  sprintf( name, "%s.y", par.names.in );
  imy = _VT_Inrimage( name );
  if ( imy == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read Y input image\n", 0);
  
  sprintf( name, "%s.z", par.names.in );
  imz = _VT_Inrimage( name );
  if ( imz == (vt_image*)NULL ) 
    fprintf( stderr, " there is no Z image\n" );


  if ( (imx->type != imy->type) ||
       ( (imz != (vt_image*)NULL) && (imx->type != imz->type) ) ) {
    VT_ErrorParse(" images have different types\n", 0 );
  }




  

  if ( par.matrice[0] != '\0' ) {
    
    if ( ReadMatrice( par.matrice, mat ) != 1 ) {
      VT_FreeImage( imx );
      VT_Free( (void**)&imx );
      VT_FreeImage( imy );
      VT_Free( (void**)&imy );
      VT_FreeImage( imz );
      VT_Free( (void**)&imz );
      VT_ErrorParse("unable to read matrice\n", 0);
      }
    if ( par.inverse != 0 ) {
      double tmp[16];
      int rang;
      fprintf( stderr, " ... on inverse la matrice \n" );
      if ( (rang=InverseMat4x4( mat, tmp )) != 4 ) {
	fprintf( stderr, "Warning: la matrice %s est de rang %d.\n",
		 par.matrice, rang );
      }
      for (i=0;i<16;i++) mat[i] = tmp[i];
    } 

  } else {
    for (i=0;i<16;i++) mat[i] = 0.0;
    mat[0] = mat[5] = mat[10] = mat[15] = 1.0;
  }
  



  switch ( imx->type ) {
  default :
    VT_ErrorParse(" such type not handled yet\n", 0 );
    break;
  case FLOAT :
    {
      r32 ***theX = (r32***)imx->array;
      r32 ***theY = (r32***)imy->array;
      r32 ***theZ = (r32***)(imz != (vt_image*)NULL ? imz->array : NULL );
      
      if ( imz == (vt_image*)NULL ) {
	for ( k=0; k<imx->dim.z; k++ )
	for ( j=0; j<imx->dim.y; j++ )
	for ( i=0; i<imx->dim.x; i++ ) {
	  x = mat[0]*(i+theX[k][j][i]) + mat[1]*(j+theY[k][j][i]) + mat[3];
	  y = mat[4]*(i+theX[k][j][i]) + mat[5]*(j+theY[k][j][i]) + mat[7];
	  theX[k][j][i] = x-i;
	  theY[k][j][i] = y-j;
	}
      } else {
	for ( k=0; k<imx->dim.z; k++ )
	for ( j=0; j<imx->dim.y; j++ )
	for ( i=0; i<imx->dim.x; i++ ) {
	  x = mat[0]*(i+theX[k][j][i]) + mat[1]*(j+theY[k][j][i]) + mat[ 2]*(k+theZ[k][j][i]) + mat[ 3];
	  y = mat[4]*(i+theX[k][j][i]) + mat[5]*(j+theY[k][j][i]) + mat[ 6]*(k+theZ[k][j][i]) + mat[ 7];
	  z = mat[8]*(i+theX[k][j][i]) + mat[9]*(j+theY[k][j][i]) + mat[10]*(k+theZ[k][j][i]) + mat[11];
	  theX[k][j][i] = x-i;
	  theY[k][j][i] = y-j;
	  theZ[k][j][i] = z-k;	
	}
      }
    }
    break;
  }
  



  sprintf( imx->name, "%s.x", par.names.out);
  (void)VT_WriteInrimage( imx );
  VT_FreeImage( imx );
  VT_Free( (void**)&imx );
  
  sprintf( imy->name, "%s.y", par.names.out);
  (void)VT_WriteInrimage( imy );
  VT_FreeImage( imy );
  VT_Free( (void**)&imy );
  
  if ( imz != (vt_image*)NULL ) {
    sprintf( imz->name, "%s.z", par.names.out);
    (void)VT_WriteInrimage( imz );
    VT_FreeImage( imz );
    VT_Free( (void**)&imz );
  }


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

      /*--- matrice ---*/
      else if ( strcmp ( argv[i], "-mat" ) == 0 ) {
	i += 1;
	if ( i >= argc) VT_ErrorParse( "parsing -mat...\n", 0 );
	strncpy( par->matrice, argv[i], STRINGLENGTH );  
      }
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->inverse = 1;
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
  par->matrice[0] = '\0';
  par->inverse = 0;
}














#if defined(_ANSI_)
static int ReadMatrice( char *name, double *mat )
#else
static int ReadMatrice( name, mat )
char *name;
double *mat;
#endif
{
  FILE *fopen(), *fp;
  char text[STRINGLENGTH];
  int nbelts = 0;
  int status;
  
  /* lecture de 4 double par ligne
     On prevoit le cas ou la ligne commence par "O8 xxxxx ...
     */

  fp = fopen( name, "r" );
  if ( fp == NULL ) return( 0 );
  
  while ( (nbelts < 16) && (fgets( text, STRINGLENGTH, fp ) != NULL) ) {
    if ( (text[0] == 'O') && (text[1] == '8') ) {
      status = sscanf( &(text[2]), "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    } else {
      status = sscanf( text, "%lf %lf %lf %lf", 
		       &mat[nbelts+0], &mat[nbelts+1],
		       &mat[nbelts+2], &mat[nbelts+3] );
    }
    if ( status == 4 ) nbelts += 4;
  }
  fclose( fp );

  if ( _VT_DEBUG_ == 1 ) {
    fprintf( stderr, " lecture de la matrice %s\n", name );
    fprintf( stderr, " %d elements lus\n", nbelts );
    fprintf( stderr,"   %f %f %f %f\n", mat[0], mat[1], mat[2], mat[3] );
    fprintf( stderr,"   %f %f %f %f\n", mat[4], mat[5], mat[6], mat[7] );
    fprintf( stderr,"   %f %f %f %f\n", mat[8], mat[9], mat[10], mat[11] );
    fprintf( stderr,"   %f %f %f %f\n", mat[12], mat[13], mat[14], mat[15] );
  }
  if ( nbelts == 16 ) return ( 1 );
  return( 0 );
}









#define TINY 1e-12

#if defined(_ANSI_)
static int InverseMat4x4( double *matrice, double *inv )
#else
static int InverseMat4x4( matrice, inv )
double *matrice;
double *inv;
#endif
{
  register int i, j, k;
  int kmax, rang = 4;
  register double c, max;
  double mat [16];
  
  for (i=0; i<16; i++ ) {
    mat[i] = matrice[i] ;
    inv[i] = 0.0;
  }
  inv[0] = inv[5] = inv[10] = inv[15] = 1.0;
  
  for ( j=0; j<4; j++ ) {
    if ( (mat[j*4+j] > (-TINY)) && (mat[j*4+j] < TINY) ) {
      /* recherche du plus grand element non nul sur la colonne j */
      kmax = j;
      max = 0.0;
      for (k=j+1; k<4; k++ ) {
	c = ( mat[k*4+j] > 0.0 ) ? mat[k*4+j] : (-mat[k*4+j]) ;
	if ( (c > TINY) && (c > max) ) { max = c; kmax = k; }
      }
      if ( kmax == j ) {
	/* la ligne est nulle */
	rang --;
      } else {
	/* sinon, on additionne */
	for ( i=0; i<4; i++ ) {
	  mat[j*4+i] += mat[kmax*4+i];
	  inv[j*4+i] += inv[kmax*4+i];
	}
      }
    }
    if ( (mat[j*4+j] < (-TINY)) || (mat[j*4+j] > TINY) ) {
      /* les autres lignes */
      for (k=0; k<4; k++) {
	if ( k != j ) {
	  c = mat[k*4 + j] / mat[j*4 + j];
	  for ( i=0; i<4; i++ ) {
	    mat[k*4 + i] -= c * mat[j*4 + i];
	    inv[k*4 + i] -= c * inv[j*4 + i];
	  }
	}
      }
      /* la ligne */
      c = mat[j*4 + j];
      for ( i=0; i<4; i++ ) {
	mat[j*4 + i] /= c;
	inv[j*4 + i] /= c;
      }
    }
  }

  return( rang );
}
