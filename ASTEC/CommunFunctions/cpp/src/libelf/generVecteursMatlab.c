/*************************************************************************
 * removeCcOnBorder.c - met a zero les composantes connexes numerotees qui
 *                      "touchent" un bord (en X ou en Y)
 *
 * $Id: generVecteursMatlab.c,v 1.2 2000/03/01 17:38:24 greg Exp $
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
 */

#include <vt_common.h>

typedef struct local_par {
  vt_names names;
  int removeLabelsOnXBorder;
  int removeLabelsOnYBorder;
  int removeLabelsOnZBorder;
  int type;
} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifdef _UNUSED_
static void VT_Parse( int argc, char *argv[], local_par *par );
#endif
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );


static char *usage = "[image-in] [image-out]\n\
\t [-border | -xborder | -yborder | -zborder]\n\
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

#define NBMAXBORDERS 2
typedef struct typeComponent {
  int label;
  int nbBorders;
  int pt[NBMAXBORDERS][3];
  double dpt[NBMAXBORDERS][3];
  int nb;
  double norme;
} typeComponent;
  


#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *imlabels, *imcarac;
  int x, y, z, i, n;
  typeComponent *theCC = (typeComponent *)NULL;
  u8 *** theCrc = (u8 ***)NULL;
  double norm;
  double ctr[3];

  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  /*
    VT_Parse( argc, argv, &par );
  */
  
  /*--- lecture de l'image d'entree ---*/
  if ( 0 )
    printf( "labels = %s carac =%s\n", argv[1], argv[2] );


  imlabels = _VT_Inrimage( argv[1] );
  if ( imlabels == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read labels image\n", 0);

  imcarac = _VT_Inrimage( argv[2] );
  if ( imlabels == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read carac image\n", 0);



  /* get max
   */
  n = 0;
  switch ( imlabels->type ) {
  case UCHAR :
    {
      u8 *** theBuf = (u8 ***)imlabels->array;
      for ( z=0; z<imlabels->dim.z ; z++ )
      for ( y=0; y<imlabels->dim.y ; y++ )
      for ( x=0; x<imlabels->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];
    }
    break;
  case USHORT :
    {
      unsigned short int *** theBuf = (unsigned short int ***)imlabels->array;
      for ( z=0; z<imlabels->dim.z ; z++ )
      for ( y=0; y<imlabels->dim.y ; y++ )
      for ( x=0; x<imlabels->dim.x ; x++ ) 
	if ( theBuf[z][y][x] > n ) n = theBuf[z][y][x];
    }
    break;
  default :
    VT_ErrorParse( "unable to deal with such image type\n", 0);
  }


  /* allocation du tableau
   */
  theCC = (typeComponent *)malloc( (n+1) * sizeof(typeComponent) );
  if ( theCC == (typeComponent *)NULL ) {
    VT_ErrorParse( "unable to allocate auxiliary array\n", 0);
  }

  for (x=0; x<=n; x++ ) {
    theCC[x].label = 0;
    theCC[x].nbBorders = 0;
    theCC[x].nb = 0;
  }
  


  /* capture des frontieres
   */
  theCrc = (u8 ***)imcarac->array;

  
  switch ( imlabels->type ) {
  case UCHAR :
    {
      u8 *** theBuf = (u8 ***)imlabels->array;
      for ( z=0; z<imlabels->dim.z ; z++ )
      for ( y=0; y<imlabels->dim.y ; y++ )
      for ( x=0; x<imlabels->dim.x ; x++ ) {
	if ( theBuf[z][y][x] == 0 ) continue;
	if ( 0 )
	  printf( " processing ( %d %d )\n", theBuf[z][y][x], theCrc[z][y][x] );

	theCC[ (int)theBuf[z][y][x] ].nb ++;

	if ( theCrc[z][y][x] == 200 ) {
	  if ( theCC[ (int)theBuf[z][y][x] ].nbBorders >= NBMAXBORDERS ) {
	    fprintf( stderr, "found more than 2 borders for components #%d\n",
		     theBuf[z][y][x] );
	    continue;
	  }
	  i = theCC[ (int)theBuf[z][y][x] ].nbBorders;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][0] = x;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][1] = y;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][2] = z;
	  theCC[ (int)theBuf[z][y][x] ].nbBorders ++;
	  if ( 0 )
	    printf( " component %2d = (%3d %3d %3d)\n",
		    theBuf[z][y][x], x, y, z );
	}
      }
    }
    break;
  case USHORT :
    {
      unsigned short int *** theBuf = (unsigned short int ***)imlabels->array;
      for ( z=0; z<imlabels->dim.z ; z++ )
      for ( y=0; y<imlabels->dim.y ; y++ )
      for ( x=0; x<imlabels->dim.x ; x++ )  {
	if ( theBuf[z][y][x] == 0 ) continue;

	theCC[ (int)theBuf[z][y][x] ].nb ++;

	if ( theCrc[z][y][x] == 200 ) {
	  if ( theCC[ (int)theBuf[z][y][x] ].nbBorders >= NBMAXBORDERS ) {
	    fprintf( stderr, "found more than 2 borders for components #%d\n",
		     theBuf[z][y][x] );
	    continue;
	  }
	  i = theCC[ (int)theBuf[z][y][x] ].nbBorders;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][0] = x;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][1] = y;
	  theCC[ (int)theBuf[z][y][x] ].pt[i][2] = z;
	  theCC[ (int)theBuf[z][y][x] ].nbBorders ++;
	}
      }
    }
    break;
  default :
    VT_ErrorParse( "unable to deal with such image type\n", 0);
  }






  for (x=1; x<=n; x++ ) {
    if ( theCC[x].nbBorders != 2 ) {
       fprintf( stderr, "component #%d does not have 2 borders\n", x );
       continue;
    }

    
    norm =  (theCC[x].pt[0][0] - theCC[x].pt[1][0])*(theCC[x].pt[0][0] - theCC[x].pt[1][0]);
    norm += (theCC[x].pt[0][1] - theCC[x].pt[1][1])*(theCC[x].pt[0][1] - theCC[x].pt[1][1]);
    norm += (theCC[x].pt[0][2] - theCC[x].pt[1][2])*(theCC[x].pt[0][2] - theCC[x].pt[1][2]);
    theCC[x].norme = sqrt( (double)norm );

    /* 1. on centre le segment sur (0,0,0)
       2. on norme chauqe demi-segment
     */
    
    for ( y=0; y<3; y++ ) {
      ctr[y] =  (double)(theCC[x].pt[0][y] + theCC[x].pt[1][y]) / (double)2.0;
      theCC[x].dpt[0][y] = (double)theCC[x].pt[0][y] - ctr[y];
      theCC[x].dpt[1][y] = (double)theCC[x].pt[1][y] - ctr[y];
    }

    for ( z=0; z<2; z++ ) {
      norm = theCC[x].dpt[z][0] * theCC[x].dpt[z][0]
	+ theCC[x].dpt[z][1] * theCC[x].dpt[z][1]
	+ theCC[x].dpt[z][2] * theCC[x].dpt[z][2];

      norm = sqrt( norm );

      for ( y=0; y<3; y++ ) {
	theCC[x].dpt[z][y] /= norm;
      }
    }
    

  }







  {
    char name[256];
     FILE *f, *fopen();
    sprintf( name, "%s.m", argv[3] );
    f = fopen( name, "w" );
    
    /*
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    for (x=1; x<=n; x++ ) {
      if ( theCC[x].nbBorders != 2 ) continue;
      fprintf( f, "x = [ %f %f ];", theCC[x].pt[0][0], theCC[x].pt[1][0] );
      fprintf( f, "y = [ %f %f ];", theCC[x].pt[0][1], theCC[x].pt[1][1] );
      fprintf( f, "z = [ %f %f ];", theCC[x].pt[0][2], theCC[x].pt[1][2] );
      fprintf( f, "plot3(x,y,z);\n" );
    }
    fprintf( f, "grid on;\n" );
    fprintf( f, "xlabel('X');\n" );
    fprintf( f, "ylabel('Y');\n" );
    fprintf( f, "zlabel('Z');\n" );
    fprintf( f, "hold off;\n" );
    */
    fprintf( f, "\n" );
    fprintf( f, "%% image de labels          = %s\n", argv[1] );
    fprintf( f, "%% image de caracterisation = %s\n", argv[2] );
    fprintf( f, "\n\n" );
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "[u,v,w] = sphere(25);\n ");
    fprintf( f, "c=ones(size(w)); \n ");
    fprintf( f, "mesh (u,v,w, c);\n");
    fprintf( f, "hidden off;\n" );
    fprintf( f, "\n" );
    for (x=1; x<=n; x++ ) {
      if ( theCC[x].nbBorders != 2 ) continue;
      fprintf( f, "%% composante #%d\n", x );
      fprintf( f, "x = [ %f %f ];", theCC[x].dpt[0][0], theCC[x].dpt[1][0] );
      fprintf( f, "y = [ %f %f ];", theCC[x].dpt[0][1], theCC[x].dpt[1][1] );
      fprintf( f, "z = [ %f %f ];", theCC[x].dpt[0][2], theCC[x].dpt[1][2] );
      fprintf( f, "plot3(x,y,z);\n" );
    }
    fprintf( f, "grid on;\n" );
    fprintf( f, "axis square\n" );
    fprintf( f, "xlabel('X');\n" );
    fprintf( f, "ylabel('Y');\n" );
    fprintf( f, "zlabel('Z');\n" );
    fprintf( f, "hold off;\n" );

    fprintf( f, "\n\n" );
    fprintf( f, "i = [" );
    for (x=1; x<=n; x++ ) {
      fprintf( f, " %d", x );
    }
    fprintf( f, " ];\n" );
    fprintf( f, "n = [" );
    for (x=1; x<=n; x++ ) {
      fprintf( f, " %d", theCC[x].nb );
    }
    fprintf( f, " ];\n" );
    fprintf( f, "l = [" );
    for (x=1; x<=n; x++ ) {
      fprintf( f, " %f", theCC[x].norme );
    }
    fprintf( f, " ];\n" );
    fprintf( f, "figure;\n" );
    fprintf( f, "hold on;\n" );
    fprintf( f, "bar(i,n,'m');\n" );
    fprintf( f, "plot(i,l,'b+');\n" );
    fprintf( f, "legend('longueur','Nombre de points');\n");
    fprintf( f, "hold off;\n" );
  }
  return( 0 );
}




#ifdef _UNUSED_
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
      /*--- traitement eventuel de l'image d'entree ---*/
      else if ( strcmp ( argv[i], "-inv" ) == 0 ) {
	par->names.inv = 1;
      }
      else if ( strcmp ( argv[i], "-swap" ) == 0 ) {
	par->names.swap = 1;
      }
      


      else if ( strcmp ( argv[i], "-border" ) == 0 ) {
	par->removeLabelsOnXBorder = 1;
	par->removeLabelsOnYBorder = 1;
	par->removeLabelsOnZBorder = 1;
      }
      else if ( strcmp ( argv[i], "-xborder" ) == 0 ) {
	par->removeLabelsOnXBorder = 1;
      }
      else if ( strcmp ( argv[i], "-yborder" ) == 0 ) {
	par->removeLabelsOnYBorder = 1;
      }
      else if ( strcmp ( argv[i], "-zborder" ) == 0 ) {
	par->removeLabelsOnZBorder = 1;
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
#endif



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
  par->type = TYPE_UNKNOWN;
  par->removeLabelsOnXBorder = 0;
  par->removeLabelsOnYBorder = 0;
  par->removeLabelsOnZBorder = 0;
}
