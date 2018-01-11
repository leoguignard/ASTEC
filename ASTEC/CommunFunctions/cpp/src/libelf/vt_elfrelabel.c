/*****************************************************************************
 * vt_elfrelabel.c - renumerote les composantes connexes par taille decroissante
 *
 * $Id: vt_elfrelabel.c,v 1.2 2000/03/09 16:20:55 greg Exp $
 *
 * Copyright (c) INRIA 2000
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 * 
 * CREATION DATE: 
 * Feb, 24 2000
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include <vt_elfrelabel.h>

static int _verbose_ = 1;

typedef struct {
  int label;
  int size;
} typeCC;

static void SortCCWithRespectToSize( typeCC *tab,
				     int left, 
				     int right );








int RelabelConnectedComponentsSortBySize( void *inputBuf,
					  bufferType typeIn,
					  int *theDim )
{
  char *proc = "RelabelConnectedComponentsSortBySize";
  int i, v;
  int lmax = 0;
  typeCC *theCC = (typeCC *)NULL;


  v = theDim[0]*theDim[1]*theDim[2];

  switch ( typeIn ) {
  case UCHAR : 
    {
      u8 *theBuf = (u8*)inputBuf;
      for ( i=0; i<v; i++ )
	if ( lmax < theBuf[i] ) lmax = theBuf[i];
    }
    break;
  case USHORT : 
    {
      u16 *theBuf = (u16*)inputBuf;
      for ( i=0; i<v; i++ )
	if ( lmax < theBuf[i] ) lmax = theBuf[i];
    }
    break;
  default :
    if ( _verbose_ ) {
      fprintf( stderr, " %s: can not deal with such image type (1).\n", proc );
    }
    return( -1 );
  }



  if ( lmax == 0 ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: null image.\n", proc );
    }
    return( -1 );
  }
  if ( lmax == 1 ) return( 1 );


  
  theCC = (typeCC*)malloc( (lmax+1)*sizeof(typeCC) );
  if ( theCC == (typeCC *)NULL ) {
    if ( _verbose_ ) {
      fprintf( stderr, " %s: can not allocate auxiliary array.\n", proc );
    }
  }
  
  for ( i=0; i<=lmax; i++ ) {
    theCC[i].label = i;
    theCC[i].size = 0;
  }



  switch ( typeIn ) {
  case UCHAR : 
    {
      u8 *theBuf = (u8*)inputBuf;
      
      for ( i=0; i<v; i++ ) {
	if ( theBuf[i] > 0 )
	  theCC[ (int)theBuf[i] ].size ++;
      }
    }
    break;
  case USHORT : 
    {
      u16 *theBuf = (u16*)inputBuf;
      
      for ( i=0; i<v; i++ ) {
	if ( theBuf[i] > 0 )
	  theCC[ (int)theBuf[i] ].size ++;
      }
    }
    break;
  default :
    if ( _verbose_ ) {
      fprintf( stderr, " %s: can not deal with such image type (2).\n", proc );
    }
    return( -1 );
  }




  SortCCWithRespectToSize( theCC, 1, lmax );
  /* ici, on a theCC[i] qui est la ieme composante
     selon la taille, et qui a le label theCC[i].label
     => on range les labels dans size
     theCC[ theCC[i].label ].size = i
     pour faire le changement
  */
  
  for ( i=1; i<=lmax; i++ ) 
    theCC[ theCC[i].label ].size = i;
  
  

  switch ( typeIn ) {
  case UCHAR : 
    {
      u8 *theBuf = (u8*)inputBuf;
      
      for ( i=0; i<v; i++ ) {
	if ( theBuf[i] > 0 )
	  theBuf[i] = theCC[ (int)theBuf[i] ].size;
      }
    }
    break;
  case USHORT : 
    {
      u16 *theBuf = (u16*)inputBuf;
      for ( i=0; i<v; i++ ) {
	if ( theBuf[i] > 0 )
	  theBuf[i] = theCC[ (int)theBuf[i] ].size;
      }
      free( theCC );
    }
    break;
  default :
    if ( _verbose_ ) {
      fprintf( stderr, " %s: can not deal with such image type (3).\n", proc );
    }
    return( -1 );
  }
  
  return( 1 );
}

























static void SortCCWithRespectToSize( typeCC *tab,
				     int left, 
				     int right )
{
  int i, last;
  typeCC tmp;

  if ( left >= right ) return;

  tmp = tab[left];   tab[left] = tab[(left+right)/2];   tab[(left+right)/2] = tmp;
  
  last = left;
  for ( i = left+1; i <= right; i++ )       
    if ( tab[i].size > tab[left].size ) {
      tmp = tab[++last];   tab[last] = tab[i];   tab[i] = tmp;
    }

  tmp = tab[left];   tab[left] = tab[last];   tab[last] = tmp;
  
  SortCCWithRespectToSize( tab, left, last-1 );
  SortCCWithRespectToSize( tab, last+1, right );
}



