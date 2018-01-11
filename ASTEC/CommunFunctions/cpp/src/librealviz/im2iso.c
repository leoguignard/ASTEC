/*************************************************************************
 * im2iso -
 *
 * $Id: im2iso.c,v 1.10 2000/08/16 16:31:51 greg Exp $
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
 *
 *
 */

#include <vt_common.h>
#include <vt_interpol.h>

#include <string.h>

typedef enum {
  _2D_ = 2,
  _3D_ = 3,
  _BINARY_
} enumComputation;

typedef enum {
  _FORWARD_,
  _BACKWARD_
} enumDirection;

typedef struct local_par {
  vt_names names;
  int type;

  float size[3];
  double refine;

  enumComputation typeComputation;
  enumDirection   theDirection;

} local_par;

/*------- Definition des fonctions statiques ----------*/
#ifndef NO_PROTO
static void VT_Parse( int argc, char *argv[], local_par *par );
static void VT_ErrorParse( char *str, int l );
static void VT_InitParam( local_par *par );
#else 
static void VT_Parse();
static void VT_ErrorParse();
static void VT_InitParam();
#endif

static char *usage = "[image-in] [image-out] [-vz %f %f %f] [-2D|-3D] \n\
\t [-forw|-back] [-refine %f]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\n\
\t Rend les images 3D isotropes par 'shape-based interpolation'\n\
\t avec une distance de chamfer 3x3x3\n\
\t\n\
\t si 'image-in' est '-', on prendra stdin\n\
\t si 'image-out' est absent, on prendra stdout\n\
\t si les deux sont absents, on prendra stdin et stdout\n\
\t\n\
\t -vz %f %f %f : taille du voxel de l'image d'entree\n\
\t                par defaut on utilise les tailles dans l'en-tete\n\
\t                de l'image d'entree\n\
\t -2D : pour une interpolation d'une image en niveau de gris\n\
\t       calcule des cartes de distance 2D (une par intensite)\n\
\t -3D : pour une interpolation d'une image en niveau de gris\n\
\t       calcule une carte de distance 3D (a partir de la surface d'elevation)\n\
\t -refine %f : en dessous de ce seuil, on utilise une distance de chamfer\n\
\t              5x5x5 pour affiner le calcul\n\
\t -forw : la premiere coupe de l'image resultat est celle de l'image d'entree\n\
\t -back : la derniere coupe de l'image resultat est celle de l'image d'entree\n\
\t\n\
\t -inv : inverse 'image-in'\n\
\t -swap : swap 'image-in' (si elle est codee sur 2 octets)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\n\
\t $Revision: 1.10 $ $Date: 2000/08/16 16:31:51 $ $Author: greg $\n";


static char program[STRINGLENGTH];














#if defined(_ANSI_)
int main( int argc, char *argv[] )
#else
  int main( argc, argv )
  int argc;
char *argv[];
#endif
{
  local_par par;
  vt_image *image, imres;

  int nbIntensites, minIntensite, maxIntensite;
  int *translationTable = (int*)NULL;


  vt_image *thePrev, *theNext, *theTmp, theImDist1, theImDist2;
  typeDistanceMap theDist;
  int newdimz;
  double z, zprevious, znext;
  int iz, izprevious, iznext;

  char distname[STRINGLENGTH];
  
  
  double a, b;


  /*--- initialisation des parametres ---*/
  VT_InitParam( &par );
  
  /*--- lecture des parametres ---*/
  VT_Parse( argc, argv, &par );
  

  /*--- lecture de l'image d'entree ---*/
  image = _VT_Inrimage( par.names.in );
  if ( image == (vt_image*)NULL ) 
    VT_ErrorParse("unable to read input image\n", 0);
  
  if ( (image->type != UCHAR) && (image->type != USHORT)  ) 
    VT_ErrorParse("unable to deal with such image\n", 0);

  /*--- operations eventuelles sur l'image d'entree ---*/
  /*
  if ( par.names.inv == 1 )  VT_InverseImage( image );
  if ( par.names.swap == 1 ) VT_SwapImage( image );
  */


  /*--- initialisation de l'image resultat ---*/
  if ( par.size[0] < 0 || par.size[1] < 0 || par.size[2] < 0 ) {
    par.size[0] = image->siz.x;
    par.size[1] = image->siz.y;
    par.size[2] = image->siz.z;
  }
  if ( par.size[2] <= par.size[0] ) 
    VT_ErrorParse("unable to deal with such voxel dimensions\n", 0);
  newdimz = 1 + (int)( ((double)(image->dim.z - 1)) * par.size[2] / par.size[0] );

  if ( 1 ) {
    fprintf(stderr, "%s: pixel size %f x %f, new slice thickness %f\n", 
	    par.names.out, par.size[0], par.size[1], par.size[0] );
    fprintf(stderr, "\t z dimension %d\n", newdimz );
  }
  



  VT_InitImage( &imres, par.names.out, image->dim.x, image->dim.y, newdimz, image->type );
  if ( VT_AllocImage( &imres ) != 1 ) {
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate output image\n", 0);
  }
  
  imres.siz.x = par.size[0];
  imres.siz.y = par.size[1];
  imres.siz.z = par.size[0];







  /* est-ce que l'image ne serait pas binaire par hasard ?
   */
  translationTable = VT_BuildTranslationTable( image, &minIntensite, &maxIntensite, &nbIntensites );
  if ( translationTable == (int*)NULL ) {
    VT_FreeImage( &imres );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("error in building translation table\n", 0);
  }

  if ( nbIntensites == 2 ) {
    par.typeComputation = _BINARY_;
    /* on s'assure que la plus petite intensite sera 0
     */
    if ( minIntensite != 0 ) {
      int i;
      int v = image->dim.x * image->dim.y * image->dim.z;
      switch ( image->type ) {
      case UCHAR :
	{
	  u8 *buf = (u8*)image->buf;
	  for (i=0;i<v;i++ ) buf[i] -= minIntensite;
	}
	break;
      case USHORT :
	{
	  u16 *buf = (u16*)image->buf;
	  for (i=0;i<v;i++ ) buf[i] -= minIntensite;
	}
	break;
      default :
	VT_FreeImage( &imres );
	VT_FreeImage( image );
	VT_Free( (void**)&image );
	VT_ErrorParse("unable to deal with such image type\n", 0);
      }
    }
  }


  if ( par.typeComputation == _2D_ ) {
    VT_ApplyTranslationTable( image, translationTable );
  } else {
    free( translationTable );
    translationTable = (int*)NULL;
  }









  /* initialisation des images de distances
   */
  switch ( par.typeComputation ) {
  case _BINARY_ :
    VT_InitImage( &theImDist1, (char*)NULL, image->dim.x, image->dim.y, 1, SSHORT );
    VT_InitImage( &theImDist2, (char*)NULL, image->dim.x, image->dim.y, 1, SSHORT );
    _DistanceSetComputationTo2D();
    break;
  case _2D_ :
    VT_InitImage( &theImDist1, (char*)NULL, image->dim.x, image->dim.y, nbIntensites-1, SSHORT );
    VT_InitImage( &theImDist2, (char*)NULL, image->dim.x, image->dim.y, nbIntensites-1, SSHORT );
    _DistanceSetComputationTo2D();
    break;
  case _3D_ :
    VT_InitImage( &theImDist1, (char*)NULL, image->dim.x, image->dim.y, 
		  maxIntensite, SSHORT );
    VT_InitImage( &theImDist2, (char*)NULL, image->dim.x, image->dim.y, 
		  maxIntensite, SSHORT );
  }
  
  if ( 1 ) {
    fprintf(stderr, " ... #intensities=%d min=%d max=%d\n",
	    nbIntensites, minIntensite, maxIntensite );
    fprintf(stderr, " ... mode " );
    switch ( par.typeComputation ) {
    case _BINARY_ :
      fprintf(stderr, "binary " ); break;
    case _2D_ :
      fprintf(stderr, "2D " ); break;
    case _3D_ :
      fprintf(stderr, "3D " ); break;
    }
    fprintf(stderr, " ---  distance map, dim along Z = %lu\n", theImDist1.dim.z );
  }













  if ( VT_AllocImage( &theImDist1 ) != 1 ) {
    VT_FreeImage( &imres );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate first distance image image\n", 0);
  }
  if ( VT_AllocImage( &theImDist2 ) != 1 ) {
    VT_FreeImage( &theImDist1 );
    VT_FreeImage( &imres );
    VT_FreeImage( image );
    VT_Free( (void**)&image );
    VT_ErrorParse("unable to allocate first distance image image\n", 0);
  }







  theDist.dim[0] = image->dim.x;
  theDist.dim[1] = image->dim.y;
  theDist.dim[2] = theImDist1.dim.z;

  theDist.voxelSize[0] = 1.0;
  theDist.voxelSize[1] = 1.0;
  theDist.voxelSize[2] = 1.0;
  theDist.multiplicativeCoefficient = 0.0;









  /* MODE FORWARD
     ============
  */
  if ( par.theDirection == _FORWARD_ ) {

    /* premier plan 
     */
    switch ( image->type ) {
    case UCHAR :
      (void)memcpy( imres.buf, image->buf, (image->dim.x * image->dim.y * sizeof(u8)) );
      break;
    case USHORT :
      (void)memcpy( imres.buf, image->buf, (image->dim.x * image->dim.y * sizeof(u16)) );
      break;
    default :
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to deal with such image type\n", 0);
    }
    

    _DistanceSetNoVerbose();



    /* distances du premier plan
     */

    sprintf( distname,"%s.dist.%d.inr", par.names.in, 0 );
    (void)VT_CopyName( theImDist2.name, distname );
    
    theDist.buf = theImDist2.buf; 
    theDist.intensityMax = VT_InitDistanceImageFromSlice( image, &theImDist2, 0 );
    _ComputeSignedDistanceMap( &theDist, par.refine );
    
    
    
    izprevious = -1;
    iznext = 0;
    thePrev = &theImDist1;
    theNext = &theImDist2;
    
    
    for ( iz = 1; iz<newdimz && iznext < image->dim.z; iz ++ ) {
      
      if ( _VT_VERBOSE_ )
	fprintf( stderr, "%s: processing slice #%3d",program, iz );
      
      z         = iz * par.size[0];
      
      zprevious = izprevious * par.size[2];
      znext     = iznext     * par.size[2];
      
      if ( z > znext ) {
	
	while ( (z > znext) && (iznext < image->dim.z) ) {
	  
	  theTmp = thePrev;   thePrev = theNext;    theNext = theTmp;
	  
	  izprevious ++;
	  iznext ++;
	  zprevious = izprevious * par.size[2];
	  znext     = iznext     * par.size[2];
	  
	  if ( iznext < image->dim.z ) {
	    
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, " distance...");
	    
	    sprintf( distname,"%s.dist.%d.inr", par.names.in, iznext );
	    (void)VT_CopyName( theNext->name, distname );
	    
	    
	    theDist.buf = theNext->buf;
	    theDist.intensityMax = VT_InitDistanceImageFromSlice( image, theNext, iznext );
	    _ComputeSignedDistanceMap( &theDist, par.refine );
	    
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, "done");
	    
	  }
	  else {
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, "                ");
	  }
	}
	
      }
      else {
	if ( _VT_VERBOSE_ )
	  fprintf( stderr, "                ");
      }
      
      if ( iznext < image->dim.z ) {
	
	a = znext - z;
	b = z - zprevious;
	if ( _VT_VERBOSE_ )
	  fprintf( stderr, " %f*[%3d] + %f*[%3d]", a, izprevious, b, iznext );
	VT_ComputeSliceFromDistances( thePrev, theNext, &imres, iz, a, b );
	
      }
      /*
      if ( _VT_VERBOSE_ )
	fprintf( stderr, "\n" );
      else
      */
      if ( _VT_VERBOSE_ )
	fprintf( stderr, "\r" );
      
    }
    
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "\n" );
  } 







  /* MODE BACKWARD
     =============
  */
  else {

    /* dernier plan 
     */
    switch ( image->type ) {
    case UCHAR :
      (void)memcpy( (void*)&(((u8***)(imres.array))[imres.dim.z-1][0][0]),
		    (void*)&(((u8***)(image->array))[image->dim.z-1][0][0]),
		    (image->dim.x * image->dim.y * sizeof(u8)) );
      break;
    case USHORT :
      (void)memcpy( (void*)&(((u16***)(imres.array))[imres.dim.z-1][0][0]),
		    (void*)&(((u16***)(image->array))[image->dim.z-1][0][0]),
		    (image->dim.x * image->dim.y * sizeof(u16)) );
      break;
    default :
      VT_FreeImage( &imres );
      VT_FreeImage( image );
      VT_Free( (void**)&image );
      VT_ErrorParse("unable to deal with such image type\n", 0);
    }
    

    _DistanceSetNoVerbose();



    /* distances du dernier plan
     */

    sprintf( distname,"%s.dist.%lu.inr", par.names.in, image->dim.z-1 );
    (void)VT_CopyName( theImDist2.name, distname );
    
    theDist.buf = theImDist2.buf; 
    theDist.intensityMax = VT_InitDistanceImageFromSlice( image, &theImDist2, image->dim.z-1 );
    _ComputeSignedDistanceMap( &theDist, par.refine );
    
    
    
    izprevious = -1;
    iznext = 0;
    thePrev = &theImDist1;
    theNext = &theImDist2;
    
    for ( iz = 1; iz<newdimz && iznext < image->dim.z; iz ++ ) {
      
      
      if ( _VT_VERBOSE_ )
	fprintf( stderr, "%s: processing slice #%3d",program, newdimz-1-iz );
      
      z         = iz * par.size[0];
      
      zprevious = izprevious * par.size[2];
      znext     = iznext     * par.size[2];
      
      if ( z > znext ) {
	
	while ( (z > znext) && (iznext < image->dim.z) ) {
	  
	  theTmp = thePrev;   thePrev = theNext;    theNext = theTmp;
	  
	  izprevious ++;
	  iznext ++;
	  zprevious = izprevious * par.size[2];
	  znext     = iznext     * par.size[2];
	  
	  if ( iznext < image->dim.z ) {
	    
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, " distance...");
	    
	    sprintf( distname,"%s.dist.%lu.inr", par.names.in, image->dim.z-1-iznext );
	    (void)VT_CopyName( theNext->name, distname );
	    
	    
	    theDist.buf = theNext->buf;
	    theDist.intensityMax = VT_InitDistanceImageFromSlice( image, theNext, image->dim.z-1-iznext );
	    _ComputeSignedDistanceMap( &theDist, par.refine );
	    
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, "done");
	    
	  }
	  else {
	    if ( _VT_VERBOSE_ )
	      fprintf( stderr, "                ");
	  }
	}
	
      }
      else {
	if ( _VT_VERBOSE_ )
	  fprintf( stderr, "                ");
      }
      
      if ( iznext < image->dim.z ) {
	
	a = znext - z;
	b = z - zprevious;
	if ( _VT_VERBOSE_ )
	  fprintf( stderr, " %f*[%3lu] + %f*[%3lu]", a, image->dim.z-1-izprevious, b, image->dim.z-1-iznext );
	VT_ComputeSliceFromDistances( thePrev, theNext, &imres, newdimz-1-iz, a, b );
	
      }
      /*
      if ( _VT_VERBOSE_ >= 2 )
	fprintf( stderr, "\n" );
      else
      */
      if ( _VT_VERBOSE_ )
	fprintf( stderr, "\r" );
      
    }
    
    if ( _VT_VERBOSE_ )
      fprintf( stderr, "\n" );

  }






  VT_FreeImage( &theImDist2 );
  VT_FreeImage( &theImDist1 );
  VT_FreeImage( image );
  VT_Free( (void**)&image );




  /* on met a jour les intensites
   */
  switch( par.typeComputation ) {

  case _BINARY_ :
    {
      int i;
      int v = imres.dim.x * imres.dim.y * imres.dim.z;
      switch ( imres.type ) {
      case UCHAR :
	{
	  u8 *buf = (u8*)imres.buf;
	  if ( minIntensite != 0 && maxIntensite != 255 ) {
	    for (i=0;i<v;i++ ) 
	      if ( buf[i] > 0 ) buf[i] = maxIntensite;
	      else              buf[i] = minIntensite;
	  }
	}
	break;
      case USHORT :
	{
	  u16 *buf = (u16*)imres.buf;
	  if ( minIntensite != 0 && maxIntensite != 65535 ) {
	    for (i=0;i<v;i++ ) 
	      if ( buf[i] > 0 ) buf[i] = maxIntensite;
	      else              buf[i] = minIntensite;
	  }
	}
	break;
      default :
	VT_FreeImage( &imres );
	VT_ErrorParse("unable to deal with such image type\n", 0);
      }
    }
    break;
  case _2D_ :
    {
      int i, j;
      for ( i=0; i<= maxIntensite; i++ ) {
	j = translationTable[ i ];
	if ( j >= 0 ) translationTable[ j ] = i;
      }
      VT_ApplyTranslationTable( &imres, translationTable );
      free( translationTable );
    }
    break;
  case _3D_ :
  default :
    break;
  }




  if ( VT_WriteInrimage( &imres ) == -1 ) {
    VT_FreeImage( &imres );
    VT_ErrorParse("unable to write output image\n", 0);
  }

  VT_FreeImage( &imres );

  return ( 0 );
}







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
	if ( _VT_VERBOSE_ <= 0 )
	  _VT_VERBOSE_ = 1;
	else
	  _VT_VERBOSE_  ++;
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

      else if ( strcmp ( argv[i], "-refine" ) == 0 ) {
	i += 1;
	if ( i >= argc)    VT_ErrorParse( "parsing -refine...\n", 0 );
	status = sscanf( argv[i],"%lf",&par->refine );
	if ( status <= 0 ) VT_ErrorParse( "parsing -refine...\n", 0 );
      }

      else if ( strcmp ( argv[i], "-vz" ) == 0 ) {
	if ( i+3 >= argc)    VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+1], "%f", &par->size[0] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+2], "%f", &par->size[1] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	status = sscanf( argv[i+3], "%f", &par->size[2] );
	if ( status <= 0 ) VT_ErrorParse( "parsing -vz...\n", 0 );
	i += 3;
      }

      else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	par->typeComputation = _2D_;
      }
      else if ( strcmp ( argv[i], "-3D" ) == 0 ) {
	par->typeComputation = _3D_;
      }

      else if ( strcmp ( argv[i], "-forw" ) == 0 ) {
	par->theDirection = _FORWARD_;
      }
      else if ( strcmp ( argv[i], "-back" ) == 0 ) {
	par->theDirection = _BACKWARD_;
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
  
}












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

  par->size[0] = par->size[1] = par->size[2] = -1.0;
  par->refine = -1.0;

  par->typeComputation = _3D_;
  par->theDirection = _FORWARD_;
}








