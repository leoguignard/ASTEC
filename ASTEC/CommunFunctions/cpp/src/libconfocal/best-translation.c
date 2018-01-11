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

double _average_ssd_u8( unsigned char *buf1, int dimx1,
			int x1_beg, int x1_end, int y1_beg, int y1_end,
			unsigned char *buf2, int dimx2,
			int x2_beg, int x2_end, int y2_beg, int y2_end,
			int xtrs_2_to_1, int ytrs_2_to_1 )
{
  int x1 ,y1, x2, y2;
  int bnx1 = x1_beg;
  int bpx1 = x1_end;
  int bny1 = y1_beg;
  int bpy1 = y1_end;
  int bnx2 = x2_beg;
  int bpx2 = x2_end;
  int bny2 = y2_beg;
  int bpy2 = y2_end;
  double v, d;
  double verr = 999999.0;

  /* x2 + xtrs_2_to_1 is paired with x1 
     x1 is between 0 and dimx1-1
     x2 is between 0 and dimx2-1
   */

  if ( bnx2 + xtrs_2_to_1 < bnx1 )       bnx2 = bnx1 - xtrs_2_to_1;
  else if ( bnx2 + xtrs_2_to_1 > bpx1 )  return( verr );
  else                                   bnx1 = bnx2 + xtrs_2_to_1;

  if ( bpx2 + xtrs_2_to_1 < bnx1 )       return( verr );
  else if ( bpx2 + xtrs_2_to_1 > bpx1 )  bpx2 = bpx1 - xtrs_2_to_1;
  else                                   bpx1 = bpx2 + xtrs_2_to_1;

  if ( bny2 + ytrs_2_to_1 < bny1 )       bny2 = bny1 - ytrs_2_to_1;
  else if ( bny2 + ytrs_2_to_1 > bpy1 )  return( verr );
  else                                   bny1 = bny2 + ytrs_2_to_1;

  if ( bpy2 + ytrs_2_to_1 < bny1 )       return( verr );
  else if ( bpy2 + ytrs_2_to_1 > bpy1 )  bpy2 = bpy1 - ytrs_2_to_1;
  else                                   bpy1 = bpy2 + ytrs_2_to_1;

  if ( 0 ) {
    printf( " X computation from [%d %d]=%d in #1 to [%d %d]=%d in #2\n",
	    bnx1, bpx1, bpx1-bnx1+1, bnx2, bpx2, bpx2-bnx2+1 );
    printf( " Y computation from [%d %d]=%d in #1 to [%d %d]=%d in #2\n",
	    bny1, bpy1, bpy1-bny1+1, bny2, bpy2, bpy2-bny2+1 );
  }

  v = 0.0;
  for ( y1 = bny1, y2 = bny2; y1 <= bpy1; y1 ++, y2 ++ ) 
  for ( x1 = bnx1, x2 = bnx2; x1 <= bpx1; x1 ++, x2 ++ ) {
    d = (double)buf1[y1 * dimx1 + x1] - (double)buf2[y2 * dimx2 + x2];
    v += d*d;
  }
  v /= ( bpy1 - bny1 + 1 ) * ( bpx1 - bnx1 + 1 );

  return( v );
}




typedef enum {
  _SSD_,
  _NONE_
} enumCriterium;


typedef struct local_par {
  vt_names names;
  int type;
  enumCriterium criterium;
  int xinit;
  int yinit;
  int x_neg_search, x_pos_search;
  int y_neg_search, y_pos_search;
  int color1;
  int color2;
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
  vt_image *image1, *image2, imres;
  int xmax, ymax;
  int x, y;
  double v, vmax;
  
  int ox1 = 0;
  int oy1 = 0;
  int ox2 = 0;
  int oy2 = 0;
  int dimx, dimy;
    
  
  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  xmax = par.xinit;
  ymax = par.yinit;

  /*--- lecture de l'image d'entree ---*/
  image1 = _VT_Inrimage( par.names.in );
  if ( image1 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #1\n", 0);
  image2 = _VT_Inrimage( par.names.ext );
  if ( image2 == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image #2\n", 0);

  if ( image1->type != image2->type )
    VT_ErrorParse("images must have the same type\n", 0);
  

  switch( image1->type ) {
  default :
    VT_ErrorParse("such type not handled yet\n", 0 );
    break;
  case UCHAR :
    {
      unsigned char *theBuf1 = (unsigned char *)image1->buf;
      unsigned char *theBuf2 = (unsigned char *)image2->buf;
      switch ( par.criterium ) {
      default : break;
      case _SSD_ :
	xmax = par.xinit;
	ymax = par.yinit;
	vmax = _average_ssd_u8( theBuf1, image1->dim.x,
				0, image1->dim.x-1, 0,image1->dim.y-1,
				theBuf2, image2->dim.x,
				0, image2->dim.x-1, 0,image2->dim.y-1,
				par.xinit, par.yinit );
	for ( y = par.y_neg_search; y <= par.y_pos_search; y++ )
        for ( x = par.x_neg_search; x <= par.x_pos_search; x++ ) {
	  v = _average_ssd_u8( theBuf1, image1->dim.x,
			       0, image1->dim.x-1, 0,image1->dim.y-1,
			       theBuf2, image2->dim.x,
			       0, image2->dim.x-1, 0,image2->dim.y-1,
			       par.xinit+x, par.yinit+y );
	  if ( vmax > v ) {
	    vmax = v;
	    xmax = par.xinit + x;
	    ymax = par.yinit + y;
	  }
	}
	break;
      } /* switch ( par.criterium ) */
    }
    break;
  }

  printf( "max: [%d %d] -> [%d %d]\n", par.xinit, par.yinit, xmax, ymax );

  
  if ( xmax < 0 ) ox1 = -xmax;
  else            ox2 = xmax;
  if ( image2->dim.x + xmax > image1->dim.x )
    dimx = image2->dim.x + xmax + ox1;
  else 
    dimx = image1->dim.x + ox2;

  if ( ymax < 0 ) oy1 = -ymax;
  else            oy2 = ymax;
  if ( image2->dim.y + ymax > image1->dim.y )
    dimy = image2->dim.y + ymax + oy1;
  else 
    dimy = image1->dim.y + oy2;
  
  if ( 1 ) {
    printf( "'%s' at (%d,%d) in [%d,%d]\n", image1->name, ox1, oy1, dimx, dimy );
    printf( "'%s' at (%d,%d) in [%d,%d]\n", image2->name, ox2, oy2, dimx, dimy );
  }



  /*--- initialisation de l'image resultat ---*/
  VT_Image( &imres );
  if ( par.names.out[0] == '\0' ) sprintf( par.names.out, "output.inr" );
  VT_InitFromImage( &imres, image1, par.names.out, image1->type );
  imres.dim.x = dimx;
  imres.dim.y = dimy;
  imres.dim.v = 3;

  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image1 );
    VT_Free( (void**)&image1 );
    VT_FreeImage( image2 );
    VT_Free( (void**)&image2 );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  

  switch( image1->type ) {
  default :
    VT_ErrorParse("such type not handled yet\n", 0 );
    break;
  case UCHAR :
    {
      unsigned char ***theBuf1 = (unsigned char ***)image1->array;
      unsigned char ***theBuf2 = (unsigned char ***)image2->array;
      unsigned char ***theBuf  = (unsigned char ***)imres.array;
      for ( y=0; y<dimy; y++ )
      for ( x=0; x<dimx; x++ )
	theBuf[0][y][ 3*x ] = theBuf[0][y][ 3*x+1 ] = theBuf[0][y][ 3*x+2 ] = 0;
      for ( y=0; y<image1->dim.y; y++ )
      for ( x=0; x<image1->dim.x; x++ )
	theBuf[0][y+oy1][3*(x+ox1) + par.color1] = theBuf1[0][y][x];
      for ( y=0; y<image2->dim.y; y++ )
      for ( x=0; x<image2->dim.x; x++ )
	theBuf[0][y+oy2][3*(x+ox2) + par.color2] = theBuf2[0][y][x];
    }
    
  }


  /*--- ecriture de l'image resultat ---*/
  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( image1 );
    VT_FreeImage( image2 );
    VT_FreeImage( &imres );
    VT_Free( (void**)&image1 );
    VT_Free( (void**)&image2 );
    VT_ErrorParse("unable to write output image\n", 0);
  }
  
  /*--- liberations memoires ---*/
  VT_FreeImage( image1 );
  VT_FreeImage( image2 );
  VT_FreeImage( &imres );
  VT_Free( (void**)&image1 );
  VT_Free( (void**)&image2 );
   
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

       else if ( strcmp ( argv[i], "-col1" ) == 0 ) {
	 i += 1;
	 if ( i >= argc)    VT_ErrorParse( "parsing -col1...\n", 0 );
	 if ( strcmp ( argv[i], "r" ) == 0 ) {
	   par->color1 = 0;
	 }
	 else if ( strcmp ( argv[i], "g" ) == 0 
		   || strcmp ( argv[i], "v" ) == 0 ) {
	   par->color1 = 1;
	 }
	 else if ( strcmp ( argv[i], "b" ) == 0 ) {
	   par->color1 = 2;
	 }
	 else {
	   fprintf( stderr, "unknown color '%s' for '-col1'\n", argv[i] );
	 }
       }
       else if ( strcmp ( argv[i], "-col2" ) == 0 ) {
	 i += 1;
	 if ( i >= argc)    VT_ErrorParse( "parsing -col2...\n", 0 );
	 if ( strcmp ( argv[i], "r" ) == 0 ) {
	   par->color2 = 0;
	 }
	 else if ( strcmp ( argv[i], "g" ) == 0 
		   || strcmp ( argv[i], "v" ) == 0 ) {
	   par->color2 = 1;
	 }
	 else if ( strcmp ( argv[i], "b" ) == 0 ) {
	   par->color2 = 2;
	 }
	 else {
	   fprintf( stderr, "unknown color '%s' for '-col2'\n", argv[i] );
	 }
       }
      
       else if ( strcmp ( argv[i], "-ssd" ) == 0 ) {
	 par->criterium = _SSD_;
       }
       else if ( strcmp ( argv[i], "-init" ) == 0 ) {
	 i += 1;
	 if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	 status = sscanf( argv[i],"%d",&(par->xinit) );
	 if ( status <= 0 ) VT_ErrorParse( "parsing -init...\n", 0 );
	 i += 1;
	 if ( i >= argc)    VT_ErrorParse( "parsing -init...\n", 0 );
	 status = sscanf( argv[i],"%d",&(par->yinit) );
	 if ( status <= 0 ) VT_ErrorParse( "parsing -init...\n", 0 );
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
  par->criterium = _NONE_;
  par->xinit = 0;
  par->yinit = 0;
  par->x_neg_search = -10;
  par->x_pos_search =  10;
  par->y_neg_search = -10;
  par->y_pos_search =  10;
  par->color1 = 0;
  par->color2 = 1;
}
