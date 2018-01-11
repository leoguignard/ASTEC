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

typedef struct local_par {
  vt_names names;
  int type;
  int x, y, z;
  int label;
  char list[STRINGLENGTH];
} local_par;




/*------- Definition des fonctions statiques ----------*/
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );





static char *usage = "[image-in] [-pt %d %d %d] [-mask %s [-label %d]] [-list %s]\n\
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

#define LMAX 256




int main( int argc, char *argv[] )
{
  local_par par;
  vt_image *image;
  vt_image *immask;
  double sum = 0.0;
  int s, i, n = 0;
  int x, y, z;
  int lnb = 0;
  int lmin, lmax, labellist[LMAX];

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


  if ( par.list[0] != '\0' ) {
    FILE *flist = NULL, *fopen();
    int ret;

    flist = fopen( par.list, "r" );
    while ( (ret = fscanf( flist, "%d\n", &labellist[lnb])) != EOF ) {
      if ( ret == 1 ) lnb ++;
    }
    fclose( flist );
  }
  else {
    labellist[0] = par.label;
    lnb = 1;
  }

  lmin = lmax = labellist[0];
  for ( i=1; i<lnb; i++ ) {
    if ( lmin > labellist[i] ) lmin = labellist[i];
    if ( lmax < labellist[i] ) lmax = labellist[i];
  }

  if ( 0 ) {
    fprintf(stderr, " label(s) =" );
    for ( i=0; i<lnb; i++ )
      fprintf( stderr, " %d", labellist[i] );
    fprintf(stderr, "\n" );
  }

  /* pas de mask
   */

  if ( par.names.ext[0] == '\0' ) {

    switch ( image->type ) {

    default :
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("such input image type is not handled\n", 0);
      
    case FLOAT :
      {
	float ***theBuf = (float***)image->array;
	
	fprintf( stdout, "%g", theBuf[par.z][par.y][par.x] );
	
      }
      break;
      
    }

  }

  /* un mask
   */
  
  else {

    immask = _VT_Inrimage( par.names.ext );
    if ( immask == (vt_image*)NULL ) {
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to read mask image\n", 0);
    }

    switch ( immask->type ) {

    default :
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_FreeImage( immask );
      VT_Free( (void**)&immask );
      VT_ErrorParse("such mask image type is not handled\n", 0);

    case UCHAR :
      {
	u8 *** masBuf = (u8***)immask->array;
	
	switch ( image->type ) {

	default :
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_FreeImage( immask );
	  VT_Free( (void**)&immask );
	  VT_ErrorParse("such input image type is not handled\n", 0);
	  
	case FLOAT :
	  {
	    float ***theBuf = (float***)image->array;
	    
	    if ( lnb > 1 ) {
	      for ( z=0; z<image->dim.z; z++ )
	      for ( y=0; y<image->dim.y; y++ )
	      for ( x=0; x<image->dim.x; x++ ) {
		if ( masBuf[z][y][x] < lmin ) continue;
		if ( masBuf[z][y][x] > lmax ) continue;
		for ( i=0, s=0; i<lnb && s==0; i++ ) {
		  if ( masBuf[z][y][x] == labellist[i] ) {
		    sum += theBuf[z][y][x];
		    n ++;
		    s = 1;
		  }
		}
	      }
	    }
	    else if ( lnb == 1 ) {
	      for ( z=0; z<image->dim.z; z++ )
	      for ( y=0; y<image->dim.y; y++ )
	      for ( x=0; x<image->dim.x; x++ ) {
		if ( masBuf[z][y][x] == labellist[0] ) {
		  sum += theBuf[z][y][x];
		  n ++;
		}
	      }
	    }
	    else {
	      for ( z=0; z<image->dim.z; z++ )
	      for ( y=0; y<image->dim.y; y++ )
	      for ( x=0; x<image->dim.x; x++ ) {
		if ( masBuf[z][y][x] > 0 ) {
		  sum += theBuf[z][y][x];
		  n ++;
		}
	      }
	    }
	    
	  }
	  break;
	}

      } /* switch ( immask->type ), case UCHAR */
      break;

    } /* end of switch ( immask->type ) */

    fprintf( stdout, "sum= %g ", sum );
    fprintf( stdout, "n= %d ", n );
    fprintf( stdout, "nsum= %g ", sum/(double)n );
    fprintf( stdout, "\n" );

    VT_FreeImage( immask );
    VT_Free( (void**)&immask );
  }



  
  /*--- liberations memoires ---*/
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


      else if ( strcmp ( argv[i], "-p" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -p x...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->x) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -p x...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -p x y...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->y) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -p x y...\n", 0 );
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -p x y z...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->z) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -p x y z...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-mask" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -mask...\n", 0 );
	strncpy( par->names.ext, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-list" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -list...\n", 0 );
	strncpy( par->list, argv[i], STRINGLENGTH );  
      }

      else if ( strcmp ( argv[i], "-label" ) == 0 || strcmp ( argv[i], "-l" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -labe;...\n", 0 );
	status = sscanf( argv[i],"%d",&(par->label) );
	if ( status <= 0 ) VT_ErrorParse( "parsing -label...\n", 0 );
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
  par->x = 0;
  par->y = 0;
  par->z = 0;
  par->label = -1;
  par->list[0] = '\0';
}
