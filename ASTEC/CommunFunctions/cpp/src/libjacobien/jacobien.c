/*************************************************************************
 * minimum.c -
 *
 * $Id: minimum.c,v 1.5 2000/08/16 16:31:56 greg Exp $
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
#include <vt_jacobien.h>

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
\t si aucune de ces options n'est presente, on prend le type de 'image-in'\n\
\n\
 $Revision: 1.5 $ $Date: 2000/08/16 16:31:56 $ $Author: greg $\n";

static char program[STRINGLENGTH];






int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image, imres;
  int x, y, z;
  int dimx, dimy, dimz, dimv;
  double n, sn;
  double nmin, nmax;

  typeListPoint3D theList;
  typeListPoint3D resList;
  float ***resBuf;

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
  
  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  
  VT_InitFromImage( &imres, image, par.names.out, FLOAT );
  imres.dim.v = 1;
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }

  if ( image->dim.v != 3 ) {
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("input image is not vectorial\n", 0);
  }

  dimx = image->dim.x;
  dimy = image->dim.y;
  dimz = image->dim.z;
  dimv = image->dim.v;


  if ( 0 ) {
    sn = nmin = nmax = 0;
    switch ( image->type ) {

    default :
      VT_FreeImage( image );
      VT_FreeImage( &imres );
      VT_Free( (void**)&image );
      VT_ErrorParse("such input image type is not handled\n", 0);

    case FLOAT :
      {
	float ***theBuf = (float***)image->array;
	for (z=0; z<dimz; z++)
	for (y=0; y<dimy; y++)
        for (x=0; x<dimx; x++) {
	  n = sqrt( theBuf[z][y][dimv*x]*theBuf[z][y][dimv*x] 
		    + theBuf[z][y][dimv*x+1]*theBuf[z][y][dimv*x+1] 
		    + theBuf[z][y][dimv*x+2]*theBuf[z][y][dimv*x+2] );
	  if ( x==0 && y==0 && z==0 )
	    nmin = nmax = n;
	  if ( nmin > n ) nmin = n;
	  if ( nmax < n ) nmax = n;
	  sn += n;
	}
	sn /= (dimx*dimy*dimz);
      }
      break;
      
    }
    
    fprintf( stdout, "normes: min = %g, max = %g, moy =%g\n", nmin, nmax, sn );
  }




  JCB_initListPoint3D( &theList );
  JCB_initListPoint3D( &resList );
  (void)JCB_get6Neighborhood( &theList );
  (void)JCB_allocListPoint3D( &resList, theList.n );
  resBuf = (float***)imres.array;

  
  switch ( image->type ) {

  default :
    VT_FreeImage( image );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image );
    VT_ErrorParse("such input image type is not handled\n", 0);
    
  case FLOAT :
    {
      float ***theBuf = (float***)image->array;
      for (z=0; z<dimz; z++)
      for (y=0; y<dimy; y++)
      for (x=0; x<dimx; x++) {

	if ( x==0 || x==dimx-1 || y==0 || y==dimy-1 || z==0 || z==dimz-1 ) {
	  resBuf[z][y][x] = 0;
	  continue;
	}

	resList.pts[0].x = theBuf[z][y][dimv*x];
	resList.pts[0].y = theBuf[z][y][dimv*x+1];
	resList.pts[0].z = theBuf[z][y][dimv*x+2];
	
	resList.pts[1].x = 1 + theBuf[z][y][dimv*(x+1)];
	resList.pts[1].y = theBuf[z][y][dimv*(x+1)+1];
	resList.pts[1].z = theBuf[z][y][dimv*(x+1)+2];
	
	resList.pts[2].x = -1 + theBuf[z][y][dimv*(x-1)];
	resList.pts[2].y = theBuf[z][y][dimv*(x-1)+1];
	resList.pts[2].z = theBuf[z][y][dimv*(x-1)+2];
	
	resList.pts[3].x = theBuf[z][(y+1)][dimv*x];
	resList.pts[3].y = 1 + theBuf[z][(y+1)][dimv*x+1];
	resList.pts[3].z = theBuf[z][(y+1)][dimv*x+2];

	resList.pts[4].x = theBuf[z][(y-1)][dimv*x];
	resList.pts[4].y = -1 + theBuf[z][(y-1)][dimv*x+1];
	resList.pts[4].z = theBuf[z][(y-1)][dimv*x+2];

	resList.pts[5].x = theBuf[(z+1)][y][dimv*x];
	resList.pts[5].y = theBuf[(z+1)][y][dimv*x+1];
	resList.pts[5].z = 1 + theBuf[(z+1)][y][dimv*x+2];

	resList.pts[6].x = theBuf[(z-1)][y][dimv*x];
	resList.pts[6].y = theBuf[(z-1)][y][dimv*x+1];
	resList.pts[6].z = -1 + theBuf[(z-1)][y][dimv*x+2];

	resBuf[z][y][x] = JCB_ComputeJacobien( &resList, &theList );

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
