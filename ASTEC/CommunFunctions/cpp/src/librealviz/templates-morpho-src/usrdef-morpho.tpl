

/***************************************************************************
 * usrdef-morpho.tpl - templates pour les procedures de morphologie mathematique
 *                     en niveau de gris dans des voisinages quelconques
 *
 * $Author: greg $
 * $Revision: 1.1 $
 * $Log: usrdef-morpho.tpl,v $
 * Revision 1.1  2001/08/10 10:31:28  greg
 * Ajout de la morphologie mathematique:
 * - fonctions de base
 * - class element structurant
 * Ajout du filtrage median
 *
 *
 *
 * $Id: usrdef-morpho.tpl,v 1.1 2001/08/10 10:31:28 greg Exp $
 ***************************************************************************/



template <class TYPE > static int _usrdef_valued_OPERATION( TYPE* T,
	TYPE* inputBuf, /* buffer to be resampled */
		      	 TYPE* resultBuf, /* result buffer */
		         unsigned int dx,
		         unsigned int dy,
		         unsigned int dz,
							  typeStructuringElementPoint *list,
							  int nb,
							  int iterations )
{
  int tmpMustBeAllocated = 0;
  TYPE *tmpBuf = NULL;
  TYPE *theBuf, *auxBuf, *resBuf, *fooBuf;
  int bufferSize;
  int centerBelongsToSE = 0;
  int i, j;
  int n;
  int x, y, z;
  int dimx = dx;
  int dimy = dy;
  int dimz = dz;

  int mdx = 0, pdx = 0;
  int mdy = 0, pdy = 0;
  int mdz = 0, pdz = 0;
  int born1x, born2x, born3x, born4x;
  int born1y, born2y, born3y, born4y;
  int born1z, born2z, born3z, born4z;

  int _INSIDE_Z_;
  int _INSIDE_Y_;

  TYPE v;
  int  vIsSet;



  if ( iterations <= 0 || nb <= 0 || list == NULL ) return( _morpho_ErrorInParameters_ );
  
  if ( iterations > 1 || inputBuf == resultBuf ) {
    tmpMustBeAllocated = 1;
  }

  bufferSize = dimx*dimy*dimz*sizeof( TYPE );

  if ( tmpMustBeAllocated ) {
    tmpBuf = (TYPE*)malloc( bufferSize );
    if ( tmpBuf == NULL ) return( _morpho_ErrMemFull_ );
  }
  //
  // iterations == 1 && inputBuf != resultBuf
  // 
  else {
    tmpBuf = resultBuf;
  }



  //
  // check whether (0,0,0) belongs to the structuring
  // element. If yes, put it at the beginning
  //
  for ( n=0; n<nb; n++ ) {
    if ( list[n].x == 0 && list[n].y == 0 && list[n].z == 0 ) {
      centerBelongsToSE = 1;
      if ( n != 0 ) {
	list[n].x = list[0].x;
	list[n].y = list[0].y;
	list[n].z = list[0].z;
	list[0].x = list[0].y = list[0].z = 0;
      }
    }
  }

  //
  // spatial extend of the structuring element
  // offset's computation
  //
  // dx is in the interval [mdx pdx]
  //
  mdx = pdx = list[0].x;
  mdy = pdy = list[0].y;
  mdz = pdz = list[0].z;
  list[0].o = list[0].x + list[0].y * dimx + list[0].z * dimx*dimy;
  for ( n=1; n<nb; n++ ) {
    if ( mdx > list[n].x ) mdx = list[n].x;
    if ( pdx < list[n].x ) pdx = list[n].x;
    if ( mdy > list[n].y ) mdy = list[n].y;
    if ( pdy < list[n].y ) pdy = list[n].y;
    if ( mdz > list[n].z ) mdz = list[n].z;
    if ( pdz < list[n].z ) pdz = list[n].z;
    list[n].o = list[n].x + list[n].y * dimx + list[n].z * dimx*dimy;
  }
  //
  // x >= born1x => x + mdx >= 0 
  // x <  born2x => x + pdx <  dimx
  // both => x coordinates of the SE belong to the image
  // 
  // x >= born3x => x + mdx >= dimx
  // x <  born4x => x + pdx <  0
  // One of these => x coordinates of the SE do not belong to the image
  born1x = (-mdx);   born2x = dimx-pdx;   born3x = dimx-mdx;   born4x = (-pdx);
  born1y = (-mdy);   born2y = dimy-pdy;   born3y = dimy-mdy;   born4y = (-pdy);
  born1z = (-mdz);   born2z = dimz-pdz;   born3z = dimz-mdz;   born4z = (-pdz);
  //
  //
  // 
  if ( born1x < 0    ) born1x = 0;
  if ( born2x > dimx ) born2x = dimx;



  /* iteration #1: (theBuf -> inputBuf) => (auxBuf -> tmpBuf)
     iteration #2: (theBuf -> tmpBuf)   => (auxBuf -> resultBuf)
     theBuf = auxBuf;
     auxBuf = resBuf;
     iteration #n: swap theBuf and auxBuf
     fooBuf = theBuf;
     theBuf = auxBuf;
     auxBuf = fooBuf;
     
     result in: theBuf
  */


  theBuf = inputBuf;
  auxBuf = tmpBuf;
  resBuf = resultBuf;


  for ( j=0; j<iterations; j++ ) {

    for (i=0, z=0; z<dimz; z++ ) {

      if ( z < born4z || z >= born3z ) {
	for ( y=0; y<dimy; y++ )
	for ( x=0; x<dimx; x++, i++ )
	  auxBuf[i] = 0;
	continue;
      }
      _INSIDE_Z_ = 0;
      if ( z >= born1z && z < born2z ) _INSIDE_Z_= 1;

      for ( y=0; y<dimy; y++ ) {
	_INSIDE_Y_ = 0;
	if ( _INSIDE_Z_== 1 && y >= born1y && y < born2y ) _INSIDE_Y_ = 1;

	// here list[n].x have to be tested
	// against mdx
	
	if ( centerBelongsToSE == 1 ) {
	  for ( x=0; x<born1x; x++, i++ ) {
	    v = theBuf[i];
	    if ( _INSIDE_Y_ == 1 ) {
	      for ( n=1; n<nb; n++ ) {
		if ( x+list[n].x < 0 ) continue;
		if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	      }
	    }
	    else {
	      for ( n=1; n<nb; n++ ) {
		if ( x+list[n].x < 0 ) continue;
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	      }
	    }
	    auxBuf[i] = v;
	  }
	}
	else {
	  for ( x=0; x<born1x; x++, i++ ) {
	    vIsSet = 0;
	    if ( _INSIDE_Y_ == 1 ) {
	      for ( n=0; n<nb; n++ ) {
		if ( x+list[n].x < 0 ) continue;
		if ( vIsSet == 0 ) {
		  v = theBuf[i+list[n].o];
		  vIsSet = 1;
		}
		else {
		  if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
		}
	      }
	    }
	    else {
	      for ( n=0; n<nb; n++ ) {
		if ( x+list[n].x < 0 ) continue;
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( vIsSet == 0 ) {
		  v = theBuf[i+list[n].o];
		  vIsSet = 1;
		}
		else {
		  if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
		}
	      }
	    }
	    auxBuf[i] = ( vIsSet == 1 ) ? v : 0;
	  }
	}

	// here no test on list[n].x
	// 
	if ( _INSIDE_Y_ == 1 ) {
	  for ( x=born1x; x<born2x; x++, i++ ) {
	    v = theBuf[i+list[0].o];
	    for ( n=1; n<nb; n++ ) {
	      if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	    }
	    auxBuf[i] = v;
	  }
	}
	else {
	  if ( centerBelongsToSE == 1 ) {
	    for ( x=born1x; x<born2x; x++, i++ ) {
	      v = theBuf[i];
	      for ( n=1; n<nb; n++ ) {
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	      }
	      auxBuf[i] = v;
	    }
	  }
	  else {
	    for ( x=born1x; x<born2x; x++, i++ ) {
	      vIsSet = 0;
	      for ( n=0; n<nb; n++ ) {
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( vIsSet == 0 ) {
		  v = theBuf[i+list[n].o];
		  vIsSet = 1;
		}
		else {
		  if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
		}
	      }
	      auxBuf[i] = ( vIsSet == 1 ) ? v : 0;
	    }
	  }
	}
      
	// here list[n].x have to be tested
	// against pdx
	
	if ( centerBelongsToSE == 1 ) {
	  for ( x=born2x; x<dimx; x++, i++ ) {
	    v = theBuf[i];
	    if ( _INSIDE_Y_ == 1 ) {
	      for ( n=1; n<nb; n++ ) {
		if ( x+list[n].x >= dimx ) continue;
		if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	      }
	    }
	    else {
	      for ( n=1; n<nb; n++ ) {
		if ( x+list[n].x >= dimx ) continue;
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
	      }
	    }
	    auxBuf[i] = v;
	  }
	}
	else {
	  for ( x=born2x; x<dimx; x++, i++ ) {
	    vIsSet = 0;
	    if ( _INSIDE_Y_ == 1 ) {
	      for ( n=0; n<nb; n++ ) {
		if ( x+list[n].x >= dimx ) continue;
		if ( vIsSet == 0 ) {
		  v = theBuf[i+list[n].o];
		  vIsSet = 1;
		}
		else {
		  if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
		}
	      }
	    }
	    else {
	      for ( n=0; n<nb; n++ ) {
		if ( x+list[n].x >= dimx ) continue;
		if ( y+list[n].y < 0 ) continue;
		if ( y+list[n].y >= dimy ) continue;
		if ( z+list[n].z < 0 ) continue;
		if ( z+list[n].z >= dimz ) continue;
		if ( vIsSet == 0 ) {
		  v = theBuf[i+list[n].o];
		  vIsSet = 1;
		}
		else {
		  if ( v _TEST_ theBuf[i+list[n].o] ) v = theBuf[i+list[n].o];
		}
	      }
	    }
	    auxBuf[i] = ( vIsSet == 1 ) ? v : 0;
	  }
	}


      } // loop on y 
    }   // loop on z

    // 
    // next iteration
    //
    if ( j == 0 ) {
      theBuf = auxBuf;
      auxBuf = resBuf;
    } else if ( j > 0 ) {
      fooBuf = theBuf;
      theBuf = auxBuf;
      auxBuf = fooBuf;
    }

  }

  if ( theBuf != resultBuf )
    memcpy( resultBuf, (void*)theBuf, bufferSize );

  if ( tmpMustBeAllocated ) free( tmpBuf );
  return( 1 );
}


