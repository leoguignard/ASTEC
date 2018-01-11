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


typedef struct {
  int x;
  int y;
  int z;
} typePoint;

typedef struct {
  int flag;
  typePoint s;
  typePoint e;
} typeSegment;



void _Draw( unsigned char * buf, int dimx, int dimy,
            typePoint *pt0, typePoint *pt1, unsigned char value )
{
  int dx, dy;
  int x, y, ix, iy;
  int i1, i2, g;

  if ( fabs( pt0->x - pt1->x ) >= fabs( pt0->y - pt1->y ) ) {

    if ( pt0->x >= pt1->x ) {
      x = (int)(pt1->x + 0.5);
      y = (int)(pt1->y + 0.5);
      dx = (int)(pt0->x + 0.5) - x;
      iy = ( pt0->y >= pt1->y ) ? 1 : -1 ;
    } 
    else {
      x = (int)(pt0->x + 0.5);
      y = (int)(pt0->y + 0.5);
      dx = (int)(pt1->x + 0.5) - x;
      iy = ( pt1->y >= pt0->y ) ? 1 : -1 ;
    }

    i1 = (int)(pt0->y + 0.5) - (int)(pt1->y + 0.5);
    if ( i1 < 0 ) i1 *= -1;
    i1 *= 2;
    
    g =  i1 - dx;
    i2 = i1 - 2*dx;

    do {
      buf[ x + y*dimx ] = value;
      x++;
      if ( g < 0 ) {
        g += i1;
      }
      else {
        y += iy;
        g += i2;
      }
    } while ( dx -- > 0 );
  }
  else {
    
    if ( pt0->y >= pt1->y ) {
      x = (int)(pt1->x + 0.5);
      y = (int)(pt1->y + 0.5);
      dy = (int)(pt0->y + 0.5) - y;
      ix = ( pt0->x >= pt1->x ) ? 1 : -1 ;
    } 
    else {
      x = (int)(pt0->x + 0.5);
      y = (int)(pt0->y + 0.5);
      dy = (int)(pt1->y + 0.5) - y;
      ix = ( pt1->x >= pt0->x ) ? 1 : -1 ;
    }

    i1 = (int)(pt0->x + 0.5) - (int)(pt1->x + 0.5);
    if ( i1 < 0 ) i1 *= -1;
    i1 *= 2;
    
    g =  i1 - dy;
    i2 = i1 - 2*dy;

    do {
      buf[ x + y*dimx ] = value;
      y++;
      if ( g < 0 ) {
        g += i1;
      }
      else {
        x += ix;
        g += i2;
      }
    } while ( dy -- > 0 );
  }

}












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
  vt_image *image;
  int max, i, j, n, x, y, z;
  typeSegment *seg;

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
  
  switch( image->type ) {
  default :
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("image type not handled yet\n", 0);
    break;
  case UCHAR :
    {
      unsigned char ***theBuf = (unsigned char ***)image->array;
      max = 0;
      for ( z=0; z<image->dim.z; z++ )
      for ( y=1; y<image->dim.y-1; y++ )
      for ( x=1; x<image->dim.x-1; x++ ) 
	if ( max < theBuf[z][y][x] ) max = theBuf[z][y][x];
      fprintf( stderr, "found %d components in %s\n", max, image->name );
      
      seg = (typeSegment *)malloc( (max+1)*sizeof( typeSegment ) );
      for (i=0; i<=max; i++ ) {
	seg[i].flag = 0;
	seg[i].s.x = seg[i].s.y = seg[i].s.z = -1;
	seg[i].e.x = seg[i].e.y = seg[i].e.z = -1;
      }

      for ( z=0; z<image->dim.z; z++ )
      for ( y=1; y<image->dim.y-1; y++ )
      for ( x=1; x<image->dim.x-1; x++ ) {
	if ( theBuf[z][y][x] == 0 ) continue;
	n = 0;
	for ( j = -1; j <= 1; j ++ )
	for ( i = -1; i <= 1; i ++ ) {
	  if ( theBuf[z][y+j][x+i] > 0 ) n++;
	}
	if ( n != 2 ) continue;
	if ( seg[ theBuf[z][y][x] ].s.x == -1 ) {
	  seg[ theBuf[z][y][x] ].s.x = x;
	  seg[ theBuf[z][y][x] ].s.y = y;
	  seg[ theBuf[z][y][x] ].s.z = z;
	}
	else if ( seg[ theBuf[z][y][x] ].e.x == -1 ) {
	  seg[ theBuf[z][y][x] ].e.x = x;
	  seg[ theBuf[z][y][x] ].e.y = y;
	  seg[ theBuf[z][y][x] ].e.z = z;
	  seg[ theBuf[z][y][x] ].flag = 1;
	}
	else {
	  seg[ theBuf[z][y][x] ].flag = -1;
	  fprintf( stderr, "component #%d has already 2 extremities:\n", 
		   theBuf[z][y][x] );
	  fprintf( stderr, "\t (%d %d %d) and (%d %d %d)\n",
		   seg[ theBuf[z][y][x] ].s.x = x,
		   seg[ theBuf[z][y][x] ].s.y = y,
		   seg[ theBuf[z][y][x] ].s.z = z,
		   seg[ theBuf[z][y][x] ].e.x = x,
		   seg[ theBuf[z][y][x] ].e.y = y,
		   seg[ theBuf[z][y][x] ].e.z = z );
	  fprintf( stderr, "\t extra extremity = (%d %d %d)\n", x, y, z );
	}
      }

      for (i=1; i<=max; i++ ) {
	if ( seg[i].flag == 1 ) {
	  fprintf( stderr, " connect #%d : (%d %d %d) and (%d %d %d)\n", i,
		   seg[i].s.x, seg[i].s.y, seg[i].s.z, 
		   seg[i].e.x, seg[i].e.y, seg[i].e.z );
	  _Draw( &(theBuf[ seg[i].s.z ][0][0] ),
		 image->dim.x, image->dim.y, 
		 &(seg[i].s), &(seg[i].e), 
		 (unsigned char)i );
	}
      }

    }
    break;
  }
  


  /*--- ecriture de l'image resultat ---*/
  strcpy( image->name, par.names.out );
  if ( VT_WriteInrimage( image ) == -1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to write output image\n", 0);
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
