/*************************************************************************
 * rl_interpol.c - procedures for image interpolation
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * July, 6 1999
 *
 * ADDITIONS, CHANGES
 *
 *
 * - Fri Aug 11 15:49:52 MET DST 2000, G. Malandain
 *   ajout du cas USHORT dans _ComputeSliceFromDistances()
 *
 * - Thu Aug 10 20:51:22 MEST 2000, G. Malandain
 *   ajout du cas USHORT dans _InitDistanceImageFromSlice();
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <rl_interpol.h>






/* _InitDistanceImageFromSlice - passage d'un plan en gris a une image 3D binaire
 *
 *
 * cas 2D (dimz == 1)
 * ==================
 *
 * - les points de valeur positive (dans l'objet) sont mis a -infini
 * - les points de valeur 0        (dans le fond) sont mis a +infini
 *
 *
 *
 * cas 3D (dimz > 1)
 * =================
 *
 * L'image 3D binaire est une carte d'elevation. 
 * La regle est, pour une colonne "au-dessus" du point (x,y)
 * - les voxels d'indice z strictement inferieur a l'intensite 
 *   (le niveau de gris) sont mis a -infini
 * - les voxels d'indice z superieur ou egal a l'intensite 
 *   sont mis a +infini.
 *
 * l'intensite est donc definie par 
 *   z | ( v(x,y,z-1)<0 && v(x,y,z)>0 )
 *
 * En particulier 
 * - un point a 0 (le fond) donne une colonne de +infini
 * - un point de valeur >= dimz donne une colonne de -infini
 *
 *
 *
 */
int _InitDistanceImageFromSlice( const void *theSlice,
				 bufferType typeSlice,
				 const int dimx,
				 const int dimy,
				 short int *theDist,
				 const int dimz )
{
  char *proc = "_InitDistanceImageFromSlice";
  int z, i, v;
  int s = dimx*dimy;
  int imax = 0;
  

  /* cas 2D
   */
  if ( dimz == 1 ) {
    imax = 1;
    switch ( typeSlice ) {
    case UCHAR :
      {
	u8 *theBuf = (u8*)theSlice;
	for ( i=0; i<s; i++ ) {
	  if ( theBuf[i] > 0 ) theDist[i] = -32768;
	  else                 theDist[i] =  32767;
	}
      }
      break;
    case USHORT :
      {
	u16 *theBuf = (u16*)theSlice;
	for ( i=0; i<s; i++ ) {
	  if ( theBuf[i] > 0 ) theDist[i] = -32768;
	  else                 theDist[i] =  32767;
	}
      }
      break;
    default :
      fprintf( stderr, "%s: such type not handled in switch (2D)\n", proc );
      return( -1 );
    }

    return( imax );
  }
  
  

  /* cas 3D 
   */
  switch ( typeSlice ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theSlice;
      for ( i=0; i<s; i++ ) {
	v = theBuf[i];
	if ( imax < v ) imax = v;
	for ( z=0; z<v && z<dimz; z++ )
	  theDist[z*s + i] = -32768;
	for ( z=v; z<dimz; z++ )
	theDist[z*s + i] = 32767;
      }
    }
    break;
  case USHORT :
    {
      u16 *theBuf = (u16*)theSlice;
      for ( i=0; i<s; i++ ) {
	v = theBuf[i];
	if ( imax < v ) imax = v;
	for ( z=0; z<v && z<dimz; z++ )
	  theDist[z*s + i] = -32768;
	for ( z=v; z<dimz; z++ )
	theDist[z*s + i] = 32767;
      }
    }
    break;
  default :
    fprintf( stderr, "%s: such type not handled in switch (3D)\n", proc );
    return( -1 );
  }


  return( imax );
}

















/* _ComputeSliceFromDistances - passage d'une image 3D polaire en un plan en gris
 *
 * L'image 3D binaire est une interpolation entre 2 cartes d'elevation. 
 * La regle est, pour une colonne "au-dessus" du point (x,y) :
 * l'intensite i est caracterisee par 
 * - quel que z < i, v(x,y,z)<0
 * - quel que z >= i, v(x,y,z)>0
 *
 * On cherche donc la transition v(x,y,i-1)<0 && v(x,y,i)>0
 */
void _ComputeSliceFromDistances( void *theSlice,
				 bufferType typeSlice,
				 const int dimx,
				 const int dimy,
				 const short int *prevDist,
				 const short int *nextDist,
				 int dimz,
				 const double prevCoef,
				 const double nextCoef )
{
  char *proc = "_ComputeSliceFromDistances";
  int i;
  int s = dimx*dimy;
  int first, last, m;
  double vf, vl, vm;
  int max = 255;
  
  switch ( typeSlice ) {
  case UCHAR :
    max = 255;
    break;
  case USHORT :
    max = 65535;
    break;
  default :
    max = 255;
  }


  if ( dimz == 1 ) {
    switch ( typeSlice ) {
    case UCHAR :
      {
	u8 *theBuf = (u8*)theSlice;
	for ( i=0; i<s; i++ ) {
	  if ( prevCoef * (double)prevDist[ i ] + nextCoef * (double)nextDist[ i ]  > 0 )
	    theBuf[i] = 0;
	  else 
	    theBuf[i] = max;
	}
      }
      break;
    case USHORT :
      {
	u16 *theBuf = (u16*)theSlice;
	for ( i=0; i<s; i++ ) {
	  if ( prevCoef * (double)prevDist[ i ] + nextCoef * (double)nextDist[ i ]  > 0 )
	    theBuf[i] = 0;
	  else 
	    theBuf[i] = max;
	}
      }
      break;
    default :
      fprintf( stderr, "%s: such type not handled in switch (2D)\n", proc );
      return;
    }

    return;
  }




  /* cas 3D 
   */
  switch ( typeSlice ) {
  case UCHAR :
    {
      u8 *theBuf = (u8*)theSlice;
       for ( i=0; i<s; i++ ) {
	 first = 0;
	 last  = dimz-1;
	 vf = prevCoef*(double)prevDist[ i ] + nextCoef*(double)nextDist[ i ];
	 vl = prevCoef*(double)prevDist[ last*s + i ]  + nextCoef*(double)nextDist[ last*s + i ];

	 /* cas evident
	  */
	 if ( vf >= 0 ) {
	   theBuf[i] = 0; 
	   continue;
	 }
	 if ( vl <= 0.0 ) {
	   if ( dimz > max ) theBuf[i] = max;
	   else              theBuf[i] = dimz;
	   continue;
	 }

	 /* on a vf < 0 et vl > 0 
	    on cherche l'intensite par dichotomie
	 */
	 do {
	   m = (first+last)/2;
	   vm = prevCoef * (double)prevDist[ m*s + i ] + nextCoef * (double)nextDist[ m*s + i ];
	   if ( vm < 0.0 )    first = m;
	   else if (vm > 0.0) last = m; 
	 } while ( (last-first > 1) && (vm != 0.0) );

	 /* condition d'arret standard
	    valeur(i,z=first)<0 && valeur(i,z=last)>0 => intensite = last
	 */
	 if ( last-first == 1 ) {
	   if ( last > max ) theBuf[i] = max;
	   else              theBuf[i] = last;
	   continue;
	 }
	 
	 /* reste la condition d'arret vm == 0
	    on teste m-1 et m+1
	 */
	 vf = prevCoef*(double)prevDist[ (m-1)*s + i ] + nextCoef*(double)nextDist[ (m-1)*s + i ];
	 vl = prevCoef*(double)prevDist[ (m+1)*s + i ] + nextCoef*(double)nextDist[ (m+1)*s + i ];
	 if ( -vf > vl ) { 
	   if ( m > max ) theBuf[i] = max;
	   else           theBuf[i] = m;
	 } else if ( -vf < vl ) {
	   if ( m+1 > max ) theBuf[i] = max;
	   else             theBuf[i] = m+1;
	 } else {
	   if ( m > max ) theBuf[i] = max;
	   else           theBuf[i] = m;
	 }
       }
    }
    break;

  case USHORT :
    {
      u16 *theBuf = (u16*)theSlice;
       for ( i=0; i<s; i++ ) {
	 first = 0;
	 last  = dimz-1;
	 vf = prevCoef*(double)prevDist[ i ] + nextCoef*(double)nextDist[ i ];
	 vl = prevCoef*(double)prevDist[ last*s + i ]  + nextCoef*(double)nextDist[ last*s + i ];

	 /* cas evident
	  */
	 if ( vf >= 0 ) {
	   theBuf[i] = 0; 
	   continue;
	 }
	 if ( vl <= 0.0 ) {
	   if ( dimz > max ) theBuf[i] = max;
	   else              theBuf[i] = dimz;
	   continue;
	 }

	 /* on a vf < 0 et vl > 0 
	    on cherche l'intensite par dichotomie
	 */
	 do {
	   m = (first+last)/2;
	   vm = prevCoef * (double)prevDist[ m*s + i ] + nextCoef * (double)nextDist[ m*s + i ];
	   if ( vm < 0.0 )    first = m;
	   else if (vm > 0.0) last = m; 
	 } while ( (last-first > 1) && (vm != 0.0) );

	 /* condition d'arret standard
	    valeur(i,z=first)<0 && valeur(i,z=last)>0 => intensite = last
	 */
	 if ( last-first == 1 ) {
	   if ( last > max ) theBuf[i] = max;
	   else              theBuf[i] = last;
	   continue;
	 }
	 
	 /* reste la condition d'arret vm == 0
	    on teste m-1 et m+1
	 */
	 vf = prevCoef*(double)prevDist[ (m-1)*s + i ] + nextCoef*(double)nextDist[ (m-1)*s + i ];
	 vl = prevCoef*(double)prevDist[ (m+1)*s + i ] + nextCoef*(double)nextDist[ (m+1)*s + i ];
	 if ( -vf > vl ) { 
	   if ( m > max ) theBuf[i] = max;
	   else           theBuf[i] = m;
	 } else if ( -vf < vl ) {
	   if ( m+1 > max ) theBuf[i] = max;
	   else             theBuf[i] = m+1;
	 } else {
	   if ( m > max ) theBuf[i] = max;
	   else           theBuf[i] = m;
	 }
       }
    }
    break;

  default :
    fprintf( stderr, "%s: such type not handled in switch (3D)\n", proc );
    return;
  }
}
















/* _TestConversion - test unitaire
 *
 */
int _TestConversion( const unsigned char *theSlice,
		     const int dimx,
		     const int dimy ) 
{
  char *proc = "_TestConversion";
  short int *theDist = (short int *)NULL;
  unsigned char *resSlice = (unsigned char *)NULL;
  int dimz = 256;
  int i;
  int s = dimx*dimy;
  int initErreurs = 0;
  int distErreurs = 0;
  typeDistanceMap theMap;
  double maxDistanceToBeUpdatedWith5x5x5Chamfer = -1.0;


  theDist = (short int *)malloc( dimx*dimy*dimz * sizeof(short int) );
  if ( theDist == (short int *)NULL ) {
    fprintf( stderr, "%s: unable to allocate distance buffer\n", proc);
    return( 0 );
  }
  resSlice = (unsigned char *)malloc( dimx*dimy * sizeof(unsigned char) );
  if ( resSlice == (unsigned char *)NULL ) {
    fprintf( stderr, "%s: unable to allocate result buffer\n", proc);
    free( theDist );
    return( 0 );
  }



  _InitDistanceImageFromSlice( (void*)theSlice, UCHAR, dimx, dimy,
			       theDist, dimz );
  _ComputeSliceFromDistances( (void*)resSlice, UCHAR, dimx, dimy,
			      theDist, theDist, dimz,
			      (double)1.0, (double)0.0 );
  for ( i=0; i<s; i++ ) {
    if ( resSlice[i] != theSlice[i] ) {
      fprintf( stderr, "%s: erreur apres initialisation au point i=%7d (x=%4d,y=%4d)\n", 
	       proc, i, i-(i/dimx)*dimx, i/dimx );
      initErreurs ++;
    }
  }
  
  theMap.buf = theDist;
  theMap.dim[0] = dimx;
  theMap.dim[1] = dimy;
  theMap.dim[2] = dimz;
  theMap.voxelSize[0] = 1.0;
  theMap.voxelSize[1] = 1.0;
  theMap.voxelSize[2] = 1.0;
  theMap.multiplicativeCoefficient = 0.0;

  _DistanceSetNoVerbose();
  _ComputeSignedDistanceMap( &theMap, maxDistanceToBeUpdatedWith5x5x5Chamfer );

  _ComputeSliceFromDistances( (void*)resSlice, UCHAR, dimx, dimy,
			      theDist, theDist, dimz,
			      (double)1.0, (double)0.0 );
  for ( i=0; i<s; i++ ) {
    if ( resSlice[i] != theSlice[i] ) {
      fprintf( stderr, "%s: erreur apres calcul de distance au point i=%7d (x=%4d,y=%4d)\n", 
	       proc, i, i-(i/dimx)*dimx, i/dimx );
      distErreurs ++;
    }
  }
  
  if ( initErreurs == 0 && distErreurs == 0 ) {
    fprintf( stderr, "%s: tests were completed and successful.\n", proc );
  }


  free( theDist );
  free( resSlice );
  return( 1 );
}
