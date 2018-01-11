/*************************************************************************
 * cellfilter.c -
 *
 * $$
 *
 * Copyright (c) INRIA 2012
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Mar 12 21:57:10 CET 2012
 *
 * ADDITIONS, CHANGES
 *
 * 
 */

static int _verbose_ = 0;

#include <chamferdistance.h>

#include <vt_common.h>
#include <vt_morpho.h>
#include <morphotools.h>

#include <ccparameters.h>
#include <vt_cellfilter.h>




typedef enum {
  VT_DILATION = 1,
  VT_EROSION = 2,
  VT_CLOSING = 3,
  VT_OPENING  = 4,
} TypeOperation;





typedef struct local_par {
  vt_names names;
  TypeOperation type_operation;
  int nb_iterations; 
  Neighborhood connexite;
  DimType dim;
  int binary_mode;
  int radius;
  int sphere;

  int chamfer;

  int low_threshold;
  int high_threshold;
  int connectivity;
  int min_size;

} local_par;






static int Neighborhood2Int ( Neighborhood N )
{
  int connectivity = 26;
  switch ( N ) {
  case N04 :
    connectivity = 4; break;
  case N06 :
    connectivity = 6; break;
  case N08 :
    connectivity = 8; break;
  case N10 :
    connectivity = 10; break;
  case N18 :
    connectivity = 18; break;
  case N26 :
    connectivity = 26; break;
  }
  return( connectivity );
}



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

static char *usage = "[image-in] [image-out] [-i %d] [-con %d | -elt %s]\n\
\t [-dil | -ero | -fer | -clo | -ouv | -ope]\n\
\t [-2D] [-bin] [-sphere] [-radius | -R %d]\n\
\t [-lt %d] [-ht %d]  [-scc %d]\n\
\t [-inv] [-swap] [-v] [-D] [-help] [options-de-type]";

static char *detail = "\
\t if 'image-in' is equal to '-', we consider stdin\n\
\t if 'image-out' is not specified, we consider stdout\n\
\t if both are not specified, we consider stdin and stdout\n\
\t [-dil]      # dilation (default)\n\
\t [-ero]      # erosion\n\
\t [-fer,-clo] # closing (dilation then erosion)\n\
\t [-ouv,-ope] # opening (erosion then dilation)\n\
\t [-i %d]     # iterations number\n\
\t [-con %d]   # structuring element (connectivity): 4, 6, 8, 10, 26\n\
\t             # default = 26\n\
\t [-elt %s]   # user-defined structuring element\n\
\t       the file has to begin with the dimension of the box\n\
\t       which contains the structuring element, eg:\n\
\t       'XDIM=3'\n\
\t       'YDIM=3'\n\
\t       'ZDIM=1'\n\
\t       lines beginning with '#' are ignored\n\
\t       points of the structuring element are indicated with positive\n\
\t       numbers except for the center which should be indicated with '+'\n\
\t       if it belongs to the SE, else with '-'. Eg:\n\
\t       '1 1 1'\n\
\t       '1 + 1'\n\
\t       '1 1 1'\n\
\t       is the classical 3x3 dilation\n\
\t [-2D]       # slice by slice computation\n\
\t [-inv]      # inverse 'image-in'\n\
\t [-swap]     # swap bytes of 'image-in' (if encoded on 2 bytes)\n\
\t -v : mode verbose\n\
\t -D : mode debug\n\
\n\
 $Revision: 1.7 $ $Date: 2006/04/14 08:37:38 $ $Author: greg $\n";

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
	Neighborhood local_connexite;
	typeStructuringElement SE;
	int theDim[3];
	int dimension;

	vt_image subimage;
	vt_image subimres;
	char subimagename[256];

	typeParameter *theCC = (typeParameter *)NULL;
	int n, nCC = 0;

	int x, y, z;
	int xleftborder, xrightborder;
	int yleftborder, yrightborder;
	int zleftborder, zrightborder;
	int d;
	int dtotal = 0;

	/*--- initialisation des parametres ---*/
	VT_InitParam( &par );
	initStructuringElement( &SE );

	/*--- lecture des parametres ---*/
	VT_Parse( argc, argv, &par );

	/*--- lecture de l'image d'entree ---*/
	image = _VT_Inrimage( par.names.in );
	if ( image == (vt_image*)NULL ) 
		VT_ErrorParse("unable to read input image\n", 0);

	/*--- operations eventuelles sur l'image d'entree ---*/
	if ( par.names.inv == 1 )  VT_InverseImage( image );
	if ( par.names.swap == 1 ) VT_SwapImage( image );

	
	if ( par.names.ext[0] != '\0' ) {
	  if ( readStructuringElement( par.names.ext, &SE ) != 1 ) {
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to read structuring element\n", 0);
	  }
	}
	else if ( par.sphere ) {
	  if ( buildStructuringElementAs3DSphere( &SE, par.radius ) != 1 ) {
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to build structuring element as a sphere\n", 0);
	  }
	}

	dimension = image->dim.z == 1 ? 2 : 3;

	SE.nbIterations = par.nb_iterations;
	local_connexite = par.connexite;
	if ( par.dim == TwoD ) {
	  SE.dimension = 2;
	  dimension = 2;
	  switch ( par.connexite ) {	
	  default :
	    break;
	  case N06 :
	    local_connexite = N04;
	    break;
	  case N10 :
	  case N18 :
	  case N26 :
	    local_connexite = N08;
	  }
	}
	SE.connectivity = Neighborhood2Int( local_connexite );
	if ( par.radius > 0 ) {
	  SE.radius = par.radius;
	  if ( 0 ) 
	    printPseudoSphereDecomposition( stderr, par.radius, 3, (char*)NULL );
	}
	if ( par.binary_mode ) useBinaryMorphologicalOperations();



	/* component parameters
	 */
	theCC = ComputeParameterFromLabels( image, &nCC );
	fprintf(stderr,"found %d connect components in '%s'\n",
		nCC,  par.names.in );
	

	/*--- initialisation de l'image resultat ---*/
        VT_Image( &imres );
	VT_InitFromImage( &imres, image, par.names.out, image->type );
        if ( VT_AllocImage( &imres ) != 1 ) {
	    freeStructuringElement( &SE );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("unable to allocate output image\n", 0);
        }
	switch ( image->type ) {
	default :
	  VT_FreeImage( &imres );
	  freeStructuringElement( &SE );
	  VT_FreeImage( image );
	  VT_Free( (void**)&image );
	  VT_ErrorParse("such image type not handled in switch (init result image)\n", 0);
	  break;
	case UCHAR :
	  {
	     u8*** theBuf = (u8***)imres.array;
	     for ( z=0; z<imres.dim.z; z++ )
	     for ( y=0; y<imres.dim.y; y++ )
	     for ( x=0; x<imres.dim.x; x++ )
	       theBuf[z][y][x] = 0;
	  }
	  break;
	case USHORT :
	  {
	     u16*** theBuf = (u16***)imres.array;
	     for ( z=0; z<imres.dim.z; z++ )
	     for ( y=0; y<imres.dim.y; y++ )
	     for ( x=0; x<imres.dim.x; x++ )
	       theBuf[z][y][x] = 0;
	  }
	  break;
	case SSHORT :
	  {
	     s16*** theBuf = (s16***)imres.array;
	     for ( z=0; z<imres.dim.z; z++ )
	     for ( y=0; y<imres.dim.y; y++ )
	     for ( x=0; x<imres.dim.x; x++ )
	       theBuf[z][y][x] = 0;
	  }
	  break;
	}
	



	for ( n=1; n<=nCC; n++ ) {
	/* for ( n=1; n<=10; n++ ) { */


	  if ( theCC[n].volume <= 0 ) continue;

	  if ( _verbose_ )
	    fprintf( stderr, "processing component #%6d", n );

	  if ( 0 )
	    fprintf( stderr, "processing component #%5d [%3dx%3dx%3d]\n", n, 
		     theCC[n].ptmax[0]-theCC[n].ptmin[0]+1,
		     theCC[n].ptmax[1]-theCC[n].ptmin[1]+1,
		     theCC[n].ptmax[2]-theCC[n].ptmin[2]+1 );
	  

	  /* margin for morphological operations
	     we set the margin at 1, except at image border
	     (setting the border to 0 may cause problems when dimy or dimx = 1)
	     => basic morphological operations have to be corrected
	   */
	  xleftborder  = ( theCC[n].ptmin[0] > 0 ) ? 1 : 0;
	  xrightborder = ( theCC[n].ptmax[0] < image->dim.x-1 ) ? 1 : 0;
	  yleftborder  = ( theCC[n].ptmin[1] > 0 ) ? 1 : 0;
	  yrightborder = ( theCC[n].ptmax[1] < image->dim.y-1 ) ? 1 : 0;
	  zleftborder  = ( theCC[n].ptmin[2] > 0 ) ? 1 : 0;
	  zrightborder = ( theCC[n].ptmax[2] < image->dim.z-1 ) ? 1 : 0;
	  switch ( par.type_operation ) {
	  default :
	  case VT_EROSION :
	  case VT_OPENING :
	    break;
	  case VT_DILATION :
	  case VT_CLOSING :
	    if ( par.radius > 0 ) {
	      xleftborder  += par.radius;
	      xrightborder += par.radius;
	      yleftborder  += par.radius;
	      yrightborder += par.radius;
	      zleftborder  += par.radius;
	      zrightborder += par.radius;
	    }
	    else if ( par.nb_iterations > 0 ) {
	      xleftborder  += par.nb_iterations;
	      xrightborder += par.nb_iterations;
	      yleftborder  += par.nb_iterations;
	      yrightborder += par.nb_iterations;
	      zleftborder  += par.nb_iterations;
	      zrightborder += par.nb_iterations;
	    }
	  }

	  if ( 0 )
	    fprintf( stderr, "processing component #%5d [(%d)%3d(%d) x (%d)%3d(%d) x (%d)%3d(%d)]\n", n, 
		     xleftborder, theCC[n].ptmax[0]-theCC[n].ptmin[0]+1, xrightborder,
		     yleftborder, theCC[n].ptmax[1]-theCC[n].ptmin[1]+1, yrightborder,
		     zleftborder, theCC[n].ptmax[2]-theCC[n].ptmin[2]+1, zrightborder );
	  

	  /* auxiliary images
	   */

	  sprintf( subimagename, "component%d.hdr",n );
	  VT_InitImage( &subimage, subimagename, 
			theCC[n].ptmax[0]-theCC[n].ptmin[0]+1 + xleftborder + xrightborder,
			theCC[n].ptmax[1]-theCC[n].ptmin[1]+1 + yleftborder + yrightborder,
			theCC[n].ptmax[2]-theCC[n].ptmin[2]+1 + zleftborder + zrightborder, image->type );
	  if ( VT_AllocImage( &subimage ) != 1 ) {
	    VT_FreeImage( &imres );
	    freeStructuringElement( &SE );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    fprintf( stderr, "error for component #%d\n", n );
	    VT_ErrorParse("unable to allocate output image\n", 0);
	  }
	  
	  sprintf( subimagename, "res-component%d.hdr",n );
	  VT_InitImage( &subimres, subimagename, 
			theCC[n].ptmax[0]-theCC[n].ptmin[0]+1 + xleftborder + xrightborder,
			theCC[n].ptmax[1]-theCC[n].ptmin[1]+1 + yleftborder + yrightborder,
			theCC[n].ptmax[2]-theCC[n].ptmin[2]+1 + zleftborder + zrightborder, image->type );
	  if ( VT_AllocImage( &subimres ) != 1 ) {
	    VT_FreeImage( &subimres );
	    VT_FreeImage( &imres );
	    freeStructuringElement( &SE );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    fprintf( stderr, "error for component #%d\n", n );
	    VT_ErrorParse("unable to allocate output image\n", 0);
	  }



	  /* copy 
	   */

	  switch ( image->type ) {
	  default :
	    VT_FreeImage( &subimres );
	    VT_FreeImage( &subimage );
	    VT_FreeImage( &imres );
	    freeStructuringElement( &SE );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("such image type not handled in switch (copy subimage)\n", 0);
	    break;
	  case UCHAR :
	    {
	      u8*** theBuf = (u8***)image->array;
	      u8*** resBuf = (u8***)subimage.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ )
		resBuf[z][y][x] = 0;
	      for ( z=0; z<=theCC[n].ptmax[2]-theCC[n].ptmin[2]; z++ )
	      for ( y=0; y<=theCC[n].ptmax[1]-theCC[n].ptmin[1]; y++ )
	      for ( x=0; x<=theCC[n].ptmax[0]-theCC[n].ptmin[0]; x++ )
		if (  theBuf[theCC[n].ptmin[2]+z][theCC[n].ptmin[1]+y][theCC[n].ptmin[0]+x] == n )
		  resBuf[zleftborder+z][yleftborder+y][xleftborder+x] = n;

	    }
	    break;
	  case USHORT :
	    {
	      u16*** theBuf = (u16***)image->array;
	      u16*** resBuf = (u16***)subimage.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ )
		resBuf[z][y][x] = 0;
	      for ( z=0; z<=theCC[n].ptmax[2]-theCC[n].ptmin[2]; z++ )
	      for ( y=0; y<=theCC[n].ptmax[1]-theCC[n].ptmin[1]; y++ )
	      for ( x=0; x<=theCC[n].ptmax[0]-theCC[n].ptmin[0]; x++ )
		if (  theBuf[theCC[n].ptmin[2]+z][theCC[n].ptmin[1]+y][theCC[n].ptmin[0]+x] == n )
		  resBuf[zleftborder+z][yleftborder+y][xleftborder+x] = n;
	    }
	    break;
	  case SSHORT :
	    {
	      s16*** theBuf = (s16***)image->array;
	      s16*** resBuf = (s16***)subimage.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ )
		resBuf[z][y][x] = 0;
	      for ( z=0; z<=theCC[n].ptmax[2]-theCC[n].ptmin[2]; z++ )
	      for ( y=0; y<=theCC[n].ptmax[1]-theCC[n].ptmin[1]; y++ )
	      for ( x=0; x<=theCC[n].ptmax[0]-theCC[n].ptmin[0]; x++ )
		if (  theBuf[theCC[n].ptmin[2]+z][theCC[n].ptmin[1]+y][theCC[n].ptmin[0]+x] == n )
		  resBuf[zleftborder+z][yleftborder+y][xleftborder+x] = n;
	    }
	    break;
	  }



	  /* processing
	   */
	  theDim[0] = subimage.dim.x;
	  theDim[1] = subimage.dim.y;
	  theDim[2] = subimage.dim.z;

	  switch ( par.type_operation ) {
	  case VT_CLOSING :
	    if ( par.chamfer ) {
	      if ( morphologicalClosingWithDistance( subimage.buf, subimres.buf, subimage.type, theDim, par.radius, dimension ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in closing\n", 0);
	      }
	    }
	    else {
	      if ( morphologicalDilation( subimage.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in dilation (closing)\n", 0);
	      }
	      if ( morphologicalErosion( subimres.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in erosion (closing)\n", 0);
	      }
	    }
	    break;

	  case VT_OPENING :
	    if ( par.chamfer ) {
	      if ( morphologicalOpeningWithDistance( subimage.buf, subimres.buf, subimage.type, theDim, par.radius, dimension ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in opening\n", 0);
	      }
	    }
	    else {
	      if ( morphologicalErosion( subimage.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in erosion (opening)\n", 0);
	      }
	      if ( morphologicalDilation( subimres.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in dilation (opening)\n", 0);
	      }
	    }

	    /* post-processing 
	     */
	    if ( par.low_threshold > 0 ) {
	      d = removeExternalExtension( &subimage,  &subimres, par.low_threshold, par.high_threshold, 
					   par.connectivity, par.min_size );
	      if ( d > 0 ) {
		dtotal += d;
		if ( _verbose_ )
		  fprintf( stderr, " - remove %d components", d );
	      }
	    }
	    break;

	  case VT_EROSION :
	    if ( par.chamfer ) {
	      if ( morphologicalErosionWithDistance( subimage.buf, subimres.buf, subimage.type, theDim, par.radius, dimension ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in erosion\n", 0);
	      }
	    }
	    else {
	      if ( morphologicalErosion( subimage.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in erosion\n", 0);
	      }
	    }
	    break;

	  case VT_DILATION :
	  default :
	    if ( par.chamfer ) {
	      if ( morphologicalDilationWithDistance( subimage.buf, subimres.buf, subimage.type, theDim, par.radius, dimension ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in dilation\n", 0);
	      }
	    }
	    else {
	      if ( morphologicalDilation( subimage.buf, subimres.buf, subimage.type, theDim, &SE ) != 1 ) {
		VT_FreeImage( &subimres );
		VT_FreeImage( &subimage );
		VT_FreeImage( &imres );
		freeStructuringElement( &SE );
		VT_FreeImage( image );
		VT_Free( (void**)&image );
		VT_ErrorParse("error in dilation\n", 0);
	      }
	    }

	  }
	  
	  /* copy 
	   */

	  switch ( image->type ) {
	  default :
	    VT_FreeImage( &subimres );
	    VT_FreeImage( &subimage );
	    VT_FreeImage( &imres );
	    freeStructuringElement( &SE );
	    VT_FreeImage( image );
	    VT_Free( (void**)&image );
	    VT_ErrorParse("such image type not handled in switch (copy subimage)\n", 0);
	    break;
	  case UCHAR :
	    {
	      u8*** theBuf = (u8***)subimres.array;
	      u8*** resBuf = (u8***)imres.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ ) {
		if ( theBuf[z][y][x] == 0 ) continue;
		if ( x-xleftborder+theCC[n].ptmin[0] < 0 || x-xleftborder+theCC[n].ptmin[0] >= imres.dim.x ) continue;
		if ( y-yleftborder+theCC[n].ptmin[1] < 0 || y-yleftborder+theCC[n].ptmin[1] >= imres.dim.y ) continue;
		if ( z-zleftborder+theCC[n].ptmin[2] < 0 || z-zleftborder+theCC[n].ptmin[2] >= imres.dim.z ) continue;
		resBuf[z-zleftborder+theCC[n].ptmin[2]][y-yleftborder+theCC[n].ptmin[1]][x-xleftborder+theCC[n].ptmin[0]] = n;
	      }
	    }
	    break;
	  case USHORT :
	    {
	      u16*** theBuf = (u16***)subimres.array;
	      u16*** resBuf = (u16***)imres.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ ) {
		if ( theBuf[z][y][x] == 0 ) continue;
		if ( x-xleftborder+theCC[n].ptmin[0] < 0 || x-xleftborder+theCC[n].ptmin[0] >= imres.dim.x ) continue;
		if ( y-yleftborder+theCC[n].ptmin[1] < 0 || y-yleftborder+theCC[n].ptmin[1] >= imres.dim.y ) continue;
		if ( z-zleftborder+theCC[n].ptmin[2] < 0 || z-zleftborder+theCC[n].ptmin[2] >= imres.dim.z ) continue;
		resBuf[z-zleftborder+theCC[n].ptmin[2]][y-yleftborder+theCC[n].ptmin[1]][x-xleftborder+theCC[n].ptmin[0]] = n;
	      }
	    }
	    break;
	  case SSHORT :
	    {
	      s16*** theBuf = (s16***)subimres.array;
	      s16*** resBuf = (s16***)imres.array;
	      for ( z=0; z<subimage.dim.z; z++ )
	      for ( y=0; y<subimage.dim.y; y++ )
	      for ( x=0; x<subimage.dim.x; x++ ) {
		if ( theBuf[z][y][x] == 0 ) continue;
		if ( x-xleftborder+theCC[n].ptmin[0] < 0 || x-xleftborder+theCC[n].ptmin[0] >= imres.dim.x ) continue;
		if ( y-yleftborder+theCC[n].ptmin[1] < 0 || y-yleftborder+theCC[n].ptmin[1] >= imres.dim.y ) continue;
		if ( z-zleftborder+theCC[n].ptmin[2] < 0 || z-zleftborder+theCC[n].ptmin[2] >= imres.dim.z ) continue;
		resBuf[z-zleftborder+theCC[n].ptmin[2]][y-yleftborder+theCC[n].ptmin[1]][x-xleftborder+theCC[n].ptmin[0]] = n;
	      }
	    }
	    break;
	  }
	  


	  

	  VT_FreeImage( &subimres );
	  VT_FreeImage( &subimage );

	  if ( _verbose_ )
	    fprintf( stderr, "\n" );

	}

	if ( par.low_threshold > 0 ) {
	  fprintf( stderr, "%d components have been removed\n", dtotal );
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
    int connexite = 0;
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
	      _verbose_ = 1;
	      _VT_VERBOSE_ = 1;
	      MorphoTools_verbose();
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

	    /*--- traitement ---*/
	    else if ( strcmp ( argv[i], "-dil" ) == 0 ) {
		par->type_operation = VT_DILATION;
	    }
	    else if ( strcmp ( argv[i], "-ero" ) == 0 ) {
		par->type_operation = VT_EROSION;
	    }
	    else if ( (strcmp ( argv[i], "-fer" ) == 0) || (strcmp ( argv[i], "-clo" ) == 0) ) {
		par->type_operation = VT_CLOSING;
	    }
	    else if ( (strcmp ( argv[i], "-ouv" ) == 0) || (strcmp ( argv[i], "-ope" ) == 0) ) {
		par->type_operation = VT_OPENING;
	    }

	    else if ( strcmp ( argv[i], "-bin" ) == 0 ) {
	      par->binary_mode = 1;
	    }
	    else if ( strcmp ( argv[i], "-2D" ) == 0 ) {
	      par->dim = TwoD;
	    }

	    else if ( strcmp ( argv[i], "-elt" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -elt...\n", 0 );
		strncpy( par->names.ext, argv[i], STRINGLENGTH );  
	    }

	    else if ( strcmp ( argv[i], "-sphere" ) == 0 ) {
	      par->sphere = 1;
	    }
	    else if ( strcmp ( argv[i], "-radius" ) == 0 || strcmp ( argv[i], "-R" ) == 0) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -radius...\n", 0 );
		status = sscanf( argv[i],"%d",&par->radius );
		if ( status <= 0 ) VT_ErrorParse( "parsing -radius...\n", 0 );
	    }

	    else if ( strcmp ( argv[i], "-i" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -i...\n", 0 );
		status = sscanf( argv[i],"%d",&par->nb_iterations );
		if ( status <= 0 ) VT_ErrorParse( "parsing -i...\n", 0 );
	    }

	    else if ( strcmp ( argv[i], "-con" ) == 0 ) {
		i += 1;
		if ( i >= argc)    VT_ErrorParse( "parsing -con...\n", 0 );
		status = sscanf( argv[i],"%d",&connexite );
		if ( status <= 0 ) VT_ErrorParse( "parsing -con...\n", 0 );
	    }

	    else if ( strcmp ( argv[i], "-chamfer" ) == 0 ) {
	      par->chamfer = 1;
	    }

	    else if ( (strcmp ( argv[i], "-sb" ) == 0) || (strcmp ( argv[i], "-lt" ) == 0) ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse( "parsing -sb|lt...\n", 0 );
	      status = sscanf( argv[i],"%d",&(par->low_threshold) );
	      if ( status <= 0 ) VT_ErrorParse( "parsing -sb|lt...\n", 0 );
	      if ( par->high_threshold < par->low_threshold )
		par->high_threshold = par->low_threshold;
	    }
	    else if ( (strcmp ( argv[i], "-sh" ) == 0) || (strcmp ( argv[i], "-ht" ) == 0) ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse( "parsing -sh|ht...\n", 0 );
	      status = sscanf( argv[i],"%d",&(par->high_threshold) );
	      if ( status <= 0 ) VT_ErrorParse( "parsing -sh|ht...\n", 0 );
	    }	    
	    else if ( strcmp( argv[i], "-tcc" ) == 0 || strcmp( argv[i], "-scc" ) == 0 ) {
	      i += 1;
	      if ( i >= argc)    VT_ErrorParse(" parsing -tcc|-scc...\n", 0 );
	      status = sscanf( argv[i],"%d",&(par->min_size) );
	      if ( status <= 0 ) VT_ErrorParse(" parsing -tcc|-scc...\n", 0 );
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
    
    /*--- type de connexite ---*/
    switch ( connexite ) {
    case 4 :
	par->connexite = N04;   break;
    case 6 :
	par->connexite = N06;   break;
    case 8 :
	par->connexite = N08;   break;
    case 10 :
	par->connexite = N10;   break;
    case 18 :
	par->connexite = N18;   break;
    case 26 :
	par->connexite = N26;   break;
    }
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
  par->type_operation = VT_DILATION;
  par->nb_iterations = 1;
  par->connexite = N26;
  par->dim = ThreeD;
  par->binary_mode = 0;
  par->radius = 0;
  par->sphere = 0;

  par->chamfer = 0;

  par->low_threshold = 0;
  par->high_threshold = 0;
  par->connectivity = 6;
  par->min_size = 1;

}
