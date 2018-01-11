#include <vt_gerecfilters.h>


static void VT_Extrema3D( vt_image *resIm, 
			  float *gx, 
			  float *gy, 
			  float *gz, 
			  float *norme[3], 
			  int z );





int GE_RecFilterOnImage( vt_image *theIm, vt_image *resIm, vt_recfilters *par )
{ 
   
  vt_image auxIm;
  int auxIsAllocated = 1;
  register float *auxBuf;
  
  vt_ipt offset, local_dim;

  u16 ***u16array, *u16line;
  double *in, *out, *work;
  int initIsOk;
  int dim, dx, dxy;
  register int i, ind, x, y, z;
  int ifirst, ilast;
  char *proc="GE_RecFilterOnImage";
  


  if ( VT_Test2Image( resIm, theIm, proc ) == -1 ) return( -1 );
  if ( par->type_filter == TYPE_UNKNOWN ) {
    VT_Error("unknown type of recursive filter",proc);
    return( 0 );
  }

  if( theIm->type != USHORT ) {
    VT_Error( "unable to deal with such input image type", proc );
    return( - 1 );
  }

  u16array = (u16***)(theIm->array);
  u16line = (u16*)NULL;
  i = theIm->dim.x;
  if ( theIm->dim.y > i ) i = theIm->dim.y;
  if ( theIm->dim.z > i ) i = theIm->dim.z;
  u16line = (u16*)VT_Malloc( (unsigned int)(i * sizeof(u16)) );
  if ( u16line == (u16*)NULL ) {
    VT_Error( "unable to allocate auxiliary line", proc );
    return( -1 );
  }


  /*--- initialisation de l'image intermediaire ---*/
  VT_Image( &auxIm );
  VT_InitImage( &auxIm, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)FLOAT );
  if ( (resIm->type == FLOAT) && (resIm->buf != theIm->buf) ) {
    auxIm.buf = resIm->buf;
    auxIsAllocated = 0;
  }
  else {
    if ( VT_AllocImage( &auxIm ) != 1 ) {
      VT_Error( "unable to allocate auxiliary image", proc );
      VT_Free( (void**)&u16line );
	return( -1 );
    }
  }
  if ( VT_CopyImage( theIm, &auxIm ) != 1 ) {
    VT_Error( "unable to copy from input image", proc );
    if ( auxIsAllocated == 1 ) VT_FreeImage( &auxIm );
    VT_Free( (void**)&u16line );
    return( -1 );
  }
  auxBuf = (float*)(auxIm.buf);
  



  /*--- dimensions ---*/
  offset.x = offset.y = offset.z = 0;
  if ( par->length_continue.x > 0 ) offset.x = par->length_continue.x;
  if ( par->length_continue.y > 0 ) offset.y = par->length_continue.y;
  if ( par->length_continue.z > 0 ) offset.z = par->length_continue.z;
  local_dim.x = theIm->dim.x + 2 * offset.x;
  local_dim.y = theIm->dim.y + 2 * offset.y;
  local_dim.z = theIm->dim.z + 2 * offset.z;
  dx  = theIm->dim.x;
  dxy = theIm->dim.x * theIm->dim.y;
  



  /*--- allocations des tableaux de travail ---*/
  i = local_dim.x;
  if ( local_dim.y > i ) i = local_dim.y;
  if ( local_dim.z > i ) i = local_dim.z;

  in = out = work = (double*)NULL;
  in   = (double*)VT_Malloc( (unsigned int)(i * sizeof( double ) ) );
  out  = (double*)VT_Malloc( (unsigned int)(i * sizeof( double ) ) );
  work = (double*)VT_Malloc( (unsigned int)(i * sizeof( double ) ) );
  if ( (in == (double*)NULL) || (out == (double*)NULL) || (work == (double*)NULL) ) {
    VT_Error( "unable to allocate auxiliary arrays", proc );
    VT_Free( (void**)&in );
    VT_Free( (void**)&out );
    VT_Free( (void**)&work );
    if ( auxIsAllocated == 1 ) VT_FreeImage( &auxIm );
    VT_Free( (void**)&u16line );
    return( -1 );
  }



  


  /*--- filtrage recursif selon X ---*/
  if ( par->derivative.x != VT_NODERIVATIVE ) {
    initIsOk = VT_InitRecursiveCoefficients( (double)(par->value_coefficient.x), par->type_filter, par->derivative.x );
    if ( (initIsOk != 1) || (local_dim.x < 4) ) {
      if ( initIsOk != 1 )
	VT_Error( "unable to initialize coefficients to filter along X", proc );
      if ( local_dim.x < 4 )
	VT_Error( "too little dimension along X", proc );
    }
    else {
      dim = theIm->dim.x;
      for ( z = 0; z < theIm->dim.z; z ++ )
      for ( y = 0; y < theIm->dim.y; y ++ ) {

	  /*--- saisie d'une ligne ---*/
	  ind = z * dxy + y * dx;
	  for ( x = 0, i = offset.x; x < theIm->dim.x; x++, i++, ind++ ) {
	    in[i] = (double)auxBuf[ind];
	    u16line[x] = u16array[z][y][x];
	  }






	  /* modification de la ligne */
	  ifirst = ilast = -1;


	  for ( x = 1; x < dim; x ++ ) {
	    /* debut d'une zone de 0
	       on garde le dernier point non nul
	    */
	    if ( (u16line[x] == 0) && (u16line[x-1] > 0) )
	      ifirst = x - 1;
	    
	    /* fin d'une zone de zero 
	       on garde le premier point non nul
	    */
	    if ( (u16line[x] > 0) && (u16line[x-1] == 0) ) {
	      ilast = x;
	      
	      /* c'est au debut */
	      if ( ifirst < 0 ) {
		for ( i=offset.x; i < offset.x + ilast; i++ )
		  in[i] = in[ offset.x + ilast ];
	      } else {
		/* c'est au milieu */
		for ( i= offset.x+ifirst+1 ; i < offset.x+ilast; i++ ) {
		  in[i] = ( (double)(offset.x+ilast - i) * in[offset.x+ifirst]
			    + (double)(i-offset.x-ifirst) * in[offset.x+ilast] )
		    / (double)( ilast - ifirst );
		}
	      }
	    }
	  }
	  /* y'a des zeros a la fin */
	  if ( (u16line[dim-1] == 0) && (ifirst >= 0) ) {
	    for ( i= offset.x+ifirst+1; i < offset.x+dim; i++ )
	      in[i] = in[offset.x+ifirst];
	  }

	  

	  
	  /*--- ajout eventuel de points ---*/
	  x = theIm->dim.x + offset.x;
	  for ( i = 0; i < offset.x; i++, x++ ) {
	    in[i] = in[offset.x];
	    in[x] = in[theIm->dim.x + offset.x - 1];
	  }
	  /*--- traitement de la ligne ---*/
	  if ( VT_RecFilterOnLine( in, out, work, out, local_dim.x ) != 1 ) {
	    VT_Error( "unable to filter along X", proc );
	    VT_Free( (void**)&in );
	    VT_Free( (void**)&out );
	    VT_Free( (void**)&work );
	    if ( auxIsAllocated == 1 ) VT_FreeImage( &auxIm );
	    VT_Free( (void**)&u16line );
	    return( -1 );
	  }
	  /*--- copie de la ligne ---*/
	  ind = z * dxy + y * dx;
	  for ( x = 0, i = offset.x; x < theIm->dim.x; x++, i++, ind++ ) {
	    if ( u16line[x] == 0 ) {
	      auxBuf[ind] = 0.0;
	    } else {
	      auxBuf[ind] = (float)out[i];
	    }
	  }
	}
    }
  }




  /*--- filtrage recursif selon Y ---*/
  if ( par->derivative.y != VT_NODERIVATIVE ) {
	initIsOk = VT_InitRecursiveCoefficients( (double)(par->value_coefficient.y), par->type_filter, par->derivative.y );
	if ( (initIsOk != 1) || (local_dim.y < 4) ) {
	    if ( initIsOk != 1 )
		VT_Error( "unable to initialize coefficients to filter along Y", proc );
	    if ( local_dim.y < 4 )
		VT_Error( "too little dimension along Y", proc );
	}
	else {
	    dim = theIm->dim.y;
	    for ( z = 0; z < theIm->dim.z; z ++ )
	    for ( x = 0; x < theIm->dim.x; x ++ ) {
		/*--- saisie d'une ligne ---*/
		ind = z * dxy + x;
		for ( y = 0, i = offset.y; y < theIm->dim.y; y++, i++, ind += dx ) {
		    in[i] = (double)auxBuf[ind];
		    u16line[y] = u16array[z][y][x];
		}

		
	  /* modification de la ligne */
	  ifirst = ilast = -1;


	  for ( y = 1; y < dim; y ++ ) {
	    /* debut d'une zone de 0
	       on garde le dernier point non nul
	    */
	    if ( (u16line[y] == 0) && (u16line[y-1] > 0) )
	      ifirst = y - 1;
	    
	    /* fin d'une zone de zero 
	       on garde le premier point non nul
	    */
	    if ( (u16line[y] > 0) && (u16line[y-1] == 0) ) {
	      ilast = y;
	      
	      /* c'est au debut */
	      if ( ifirst < 0 ) {
		for ( i=offset.y; i < offset.y + ilast; i++ )
		  in[i] = in[ offset.y + ilast ];
	      } else {
		/* c'est au milieu */
		for ( i= offset.y+ifirst+1 ; i < offset.y+ilast; i++ ) {
		  in[i] = ( (double)(offset.y+ilast - i) * in[offset.y+ifirst]
			    + (double)(i-offset.y-ifirst) * in[offset.y+ilast] )
		    / (double)( ilast - ifirst );
		}
	      }
	    }
	  }
	  /* y'a des zeros a la fin */
	  if ( (u16line[dim-1] == 0) && (ifirst >= 0) ) {
	    for ( i= offset.y+ifirst+1; i < offset.y+dim; i++ )
	      in[i] = in[offset.y+ifirst];
	  }


		/*--- ajout eventuel de points ---*/
		y = theIm->dim.y + offset.y;
		for ( i = 0; i < offset.y; i++, y++ ) {
		    in[i] = in[offset.y];
		    in[y] = in[theIm->dim.y + offset.y - 1];
		}
		/*--- traitement de la ligne ---*/
		if ( VT_RecFilterOnLine( in, out, work, out, local_dim.y ) != 1 ) {
		    VT_Error( "unable to filter along Y", proc );
		    VT_Free( (void**)&in );
		    VT_Free( (void**)&out );
		    VT_Free( (void**)&work );
		    if ( auxIsAllocated == 1 ) VT_FreeImage( &auxIm );
		    VT_Free( (void**)&u16line );
		    return( -1 );
		}
		/*--- copie de la ligne ---*/
		ind = z * dxy + x;
		for ( y = 0, i = offset.y; y < theIm->dim.y; y++, i++, ind += dx ) {
		  if ( u16line[y] == 0 ) {
		    auxBuf[ind] = 0.0;
		  } else {
		    auxBuf[ind] = (float)out[i];
		  }
		}
	    }
	}
    }

    /*--- filtrage recursif selon Z ---*/
    if ( par->derivative.z != VT_NODERIVATIVE ) {
	initIsOk = VT_InitRecursiveCoefficients( (double)(par->value_coefficient.z), par->type_filter, par->derivative.z );
	if ( (initIsOk != 1) || (local_dim.z < 4) ) {
	    if ( initIsOk != 1 )
		VT_Error( "unable to initialize coefficients to filter along Z", proc );
	    if ( local_dim.z < 4 )
		VT_Error( "too little dimension along Z", proc );
	}
	else {
	    dim = theIm->dim.z;
	    for ( y = 0; y < theIm->dim.y; y ++ )
            for ( x = 0; x < theIm->dim.x; x ++ ) {
		/*--- saisie d'une ligne ---*/
		ind = y * dx + x;
		for ( z = 0, i = offset.z; z < theIm->dim.z; z++, i++, ind += dxy ) {
		    in[i] = (double)auxBuf[ind];
		    u16line[z] = u16array[z][y][x];
		}


	  /* modification de la ligne */
	  ifirst = ilast = -1;


	  for ( z = 1; z < dim; z ++ ) {
	    /* debut d'une zone de 0
	       on garde le dernier point non nul
	    */
	    if ( (u16line[z] == 0) && (u16line[z-1] > 0) )
	      ifirst = z - 1;
	    
	    /* fin d'une zone de zero 
	       on garde le premier point non nul
	    */
	    if ( (u16line[z] > 0) && (u16line[z-1] == 0) ) {
	      ilast = z;
	      
	      /* c'est au debut */
	      if ( ifirst < 0 ) {
		for ( i=offset.z; i < offset.z + ilast; i++ )
		  in[i] = in[ offset.z + ilast ];
	      } else {
		/* c'est au milieu */
		for ( i= offset.z+ifirst+1 ; i < offset.z+ilast; i++ ) {
		  in[i] = ( (double)(offset.z+ilast - i) * in[offset.z+ifirst]
			    + (double)(i-offset.z-ifirst) * in[offset.z+ilast] )
		    / (double)( ilast - ifirst );
		}
	      }
	    }
	  }
	  /* y'a des zeros a la fin */
	  if ( (u16line[dim-1] == 0) && (ifirst >= 0) ) {
	    for ( i= offset.z+ifirst+1; i < offset.z+dim; i++ )
	      in[i] = in[offset.z+ifirst];
	  }






		/*--- ajout eventuel de points ---*/
		z = theIm->dim.z + offset.z;
		for ( i = 0; i < offset.z; i++, z++ ) {
		    in[i] = in[offset.z];
		    in[z] = in[theIm->dim.z + offset.z - 1];
		}
		/*--- traitement de la ligne ---*/
		if ( VT_RecFilterOnLine( in, out, work, out, local_dim.z ) != 1 ) {
		    VT_Error( "unable to filter along Z", proc );
		    VT_Free( (void**)&in );
		    VT_Free( (void**)&out );
		    VT_Free( (void**)&work );
		    if ( auxIsAllocated == 1 ) VT_FreeImage( &auxIm );
		    VT_Free( (void**)&u16line );
		    return( -1 );
		}
		/*--- copie de la ligne ---*/
		ind = y * dx + x;
		for ( z = 0, i = offset.z; z < theIm->dim.z; z++, i++, ind += dxy ) {
		  if ( u16line[z] == 0 ) {
		    auxBuf[ind] = 0.0;
		  } else {
		    auxBuf[ind] = (float)out[i];
		  }
		}
		 
	    }
	}
    }
    
    /*--- liberations memoire ---*/
    VT_Free( (void**)&in );
    VT_Free( (void**)&out );
    VT_Free( (void**)&work );
    VT_Free( (void**)&u16line );

    /*--- liberations memoire et copie du resultat---*/
    if ( auxIsAllocated == 1 ) {
	if ( VT_CopyImage( &auxIm, resIm ) != 1 ) {
	    VT_Error( "unable to copy into output image", proc );
	    VT_FreeImage( &auxIm );
	    VT_Free( (void**)&u16line );
	    return( -1 );
	}
    }

    return( 1 );
}



















int GE_MaximaGradient( vt_image *theIm, vt_image *resIm, vt_contours *par )
{
    vt_recfilters rpar;
    vt_image theGx, theGy, theGz, theNo;
    float *lgx, *lgy, *lgz, *lno, *norme[3];
    register int i, z;
    int dxy;
    char *proc="GE_MaximaGradient";

    if ( VT_Test2Image( resIm, theIm, proc ) == -1 ) return( -1 );
    if( theIm->type != USHORT ) {
    VT_Error( "unable to deal with such input image type", proc );
    return( - 1 );
  }

    /*--- preparation des parametres de filtrage ---*/
    VT_RecFilters( &rpar );
    if ( par->value_coefficient.x > 0.0 ) rpar.value_coefficient.x = par->value_coefficient.x;
    if ( par->value_coefficient.y > 0.0 ) rpar.value_coefficient.y = par->value_coefficient.y;
    if ( par->value_coefficient.z > 0.0 ) rpar.value_coefficient.z = par->value_coefficient.z;
    if ( par->length_continue.x > 0 ) rpar.length_continue.x = par->length_continue.x;
    if ( par->length_continue.y > 0 ) rpar.length_continue.y = par->length_continue.y;
    if ( par->length_continue.z > 0 ) rpar.length_continue.z = par->length_continue.z;
    switch ( par->type_filter ) {
    case VT_RECFILTERS_DERICHE :
    case VT_RECGAUSSIAN_DERICHE :
	rpar.type_filter = par->type_filter;
	break;
    default :
	rpar.type_filter = VT_RECFILTERS_DERICHE;
    }




    /*--- initialisation ---*/
    dxy = theIm->dim.x * theIm->dim.y;
    norme[0] = norme[1] = norme[2] = (float*)NULL;
    VT_Image( &theGx );
    VT_Image( &theGy );
    VT_Image( &theGz );
    VT_Image( &theNo );


    /*--- extraction des maxima du gradient ---*/
    /*--- initialisation des images ---*/
    VT_InitImage( &theGz, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)FLOAT );
    VT_InitImage( &theGx, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)FLOAT );
    VT_InitImage( &theGy, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)FLOAT );
    VT_InitImage( &theNo, "", theIm->dim.x, theIm->dim.y, theIm->dim.z, (int)FLOAT );

    if ( VT_AllocImage( &theGx ) != 1 ) {
      VT_Error( "unable to allocate X gradient image", proc );
      return( -1 );
    }
    if ( VT_AllocImage( &theGy ) != 1 ) {
      VT_FreeImage( &theGx );
      VT_Error( "unable to allocate Y gradient image", proc );
      return( -1 );
    }
    if ( VT_AllocImage( &theNo ) != 1 ) {
      VT_FreeImage( &theGx );
      VT_FreeImage( &theGy );
      VT_Error( "unable to allocate norm image", proc );
      return( -1 );
    }
    if ( VT_AllocImage( &theGz ) != 1 ) {
      VT_FreeImage( &theNo );
      VT_FreeImage( &theGx );
      VT_FreeImage( &theGy );
      VT_Error( "unable to allocate Z gradient image", proc );
      return( -1 );
    }


    
    /*--- calcul des convolutions ---*/
    rpar.derivative.x = VT_DERIVATIVE_1_CONTOURS;   rpar.derivative.y = VT_DERIVATIVE_0;   rpar.derivative.z = VT_DERIVATIVE_0;
    if ( GE_RecFilterOnImage( theIm, &theGx, &rpar ) != 1 ) {
      VT_FreeImage( &theGx );
      VT_FreeImage( &theGy );
      VT_FreeImage( &theGz );
      VT_FreeImage( &theNo );
      VT_Error( "unable to compute X gradient image", proc );
      return( -1 );
    }
    rpar.derivative.x = VT_DERIVATIVE_0;   rpar.derivative.y = VT_DERIVATIVE_1_CONTOURS;   rpar.derivative.z = VT_DERIVATIVE_0;
    if ( GE_RecFilterOnImage( theIm, &theGy, &rpar ) != 1 ) {
      VT_FreeImage( &theGx );
      VT_FreeImage( &theGy );
      VT_FreeImage( &theGz );
      VT_FreeImage( &theNo );
      VT_Error( "unable to compute Y gradient image", proc );
      return( -1 );
    }
    rpar.derivative.x = VT_DERIVATIVE_0;   rpar.derivative.y = VT_DERIVATIVE_0;   rpar.derivative.z = VT_DERIVATIVE_1_CONTOURS;
    if ( GE_RecFilterOnImage( theIm, &theGz, &rpar ) != 1 ) {
      VT_FreeImage( &theGx );
      VT_FreeImage( &theGy );
      VT_FreeImage( &theGz );
      VT_FreeImage( &theNo );
      VT_Error( "unable to compute Z gradient image", proc );
      return( -1 );
    }
    




    /* calcul de la norme */
    lgx = (float*)theGx.buf;
    lgy = (float*)theGy.buf;
    lgz = (float*)theGz.buf;
    lno = (float*)theNo.buf;
    z = dxy * theIm->dim.z;
    for ( i = 0; i < z; i++, lno++, lgx++, lgy++, lgz++ )
      *lno = (float)VT_Sqrt( (double)((*lgx) * (*lgx) + (*lgy) * (*lgy) + (*lgz) * (*lgz)) );





    /*
    VT_CopyImage( &theNo, resIm );
    return( 1 );
    */



    /* premier et dernier plan : mise a zero */
    lgx = (float*)theGx.buf;
    for ( i = 0; i < dxy; i++ ) lgx[i] = 0.0;
    
    lgx = (float*)theGx.buf;
    lgx += (theIm->dim.z-1)*dxy;
    for ( i = 0; i < dxy; i++ ) lgx[i] = 0.0;
    



    /*--- autres plans ---*/
    /*--------------------*/
    lgx = (float*)theGx.buf;
    lgy = (float*)theGy.buf;
    lgz = (float*)theGz.buf;
    norme[1] = (float*)theNo.buf;
    

    for ( z = 1; z < theIm->dim.z-1; z++ ) {
      
      /*--- calcul des extrema du gradient du plan courant ---*/
      lgx += dxy;
      lgy += dxy;
      lgz += dxy;


      norme[0] = norme[1];
      norme[1] += dxy;
      norme[2] = norme[1] + dxy;

      VT_Extrema3D( &theGx, lgx, lgy, lgz, norme, z );
     
    }


    VT_FreeImage( &theGy );
    VT_FreeImage( &theGz );
    VT_FreeImage( &theNo );
    

    /*--- liberations memoire et copie du resultat---*/
    if ( VT_CopyImage( &theGx, resIm ) != 1 ) {
      VT_Error( "unable to copy into output image", proc );
      VT_FreeImage( &theGx );
      return( -1 );
    }

    VT_FreeImage( &theGx );

    return( 1 );
}











#if defined(_ANSI_)
static void VT_Extrema3D( vt_image *resIm, float *gx, float *gy, float *gz, float *norme[3], int z_offset )
#else
static void VT_Extrema3D( resIm, gx, gy, gz, norme, z_offset )
vt_image *resIm;
float *gx;
float *gy;
float *gz;
float *norme[3];
int z_offset;
#endif
{
    register double ngx, ngy, ngz, norm, norme_originale;
    register double xr, yr, zr, dx, dy, dz;
    register int x, y, z, xi, yi, zi, i, j, pos;
    int dimx, dimy, dimx1, dimy1, off;
    float *vl;
    double coef[2][2][2];
    double EPSILON = 0.00001; /* 1/2^16 = .00001525878906250000 */

    /*--- initialisation ---*/
    dimx = resIm->dim.x;
    dimy = resIm->dim.y;
    dimx1 = dimx - 1;
    dimy1 = dimy - 1;
    off = z_offset*dimx*dimy;
    vl = (float*)resIm->buf;

    /*--- on met les bords a zero ---*/
    for ( x=0, i=off, j=off+dimx*dimy1; x < dimx; x++, i++, j++ )
	*(vl+i) = *(vl+j) = 0.0;
    for ( y=1, i=off+dimx, j=off+2*dimx-1; y < dimy1; y++, i+=dimx, j+=dimx )
	*(vl+i) = *(vl+j) = 0.0;

    /*--- on se place correctement ---*/
    vl += off;
    z = 1;

    /*--- on parcourt le centre de l'image ---*/
    for ( y = 1; y < dimy1; y++ )
    for ( x=1, j=y*dimx+1; x < dimx1; x++, j++ ) {

	/*--- si la norme est trop petite, on ne considere pas le point ---*/
	norme_originale = *(norme[1]+j);
	if ( norme_originale <= EPSILON ) {
	    *(vl+j) = 0.0;
	    continue;
	}

	/*--- on normalise la norme du gradient ---*/
	ngx = *(gx+j) / norme_originale;
	ngy = *(gy+j) / norme_originale;
	ngz = *(gz+j) / norme_originale;

	/*--- ou trouve-t-on le premier point ? ---*/
	xr = (double)x + ngx;   yr = (double)y + ngy;   zr = (double)z + ngz;
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= 2) ) {
	   *(vl+j) = 0.0;
	  continue;
	}
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation lineaire : calcul des coefficients ---*/
	coef[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coef[1][0][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coef[0][1][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coef[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coef[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coef[1][0][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coef[0][1][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coef[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	pos = xi + yi * dimx;
	norm  = *(norme[zi]+pos)        * coef[0][0][0] + *(norme[zi]+pos+1)        * coef[1][0][0];
	norm += *(norme[zi]+pos+dimx)   * coef[0][1][0] + *(norme[zi]+pos+1+dimx)   * coef[1][1][0];
	norm += *(norme[zi+1]+pos)      * coef[0][0][1] + *(norme[zi+1]+pos+1)      * coef[1][0][1];
	norm += *(norme[zi+1]+pos+dimx) * coef[0][1][1] + *(norme[zi+1]+pos+1+dimx) * coef[1][1][1];

	/*--- comparaison par rapport au premier point : extrema ? ---*/
	if ( norme_originale <= norm ) {
	    *(vl+j) = 0.0;
	    continue;
	}

	/*--- second point ---*/
	xr = (double)x - ngx;   yr = (double)y - ngy;   zr = (double)z - ngz;
	if ( (xr < 0.0) || (xr >= dimx1) || (yr < 0.0) || (yr >= dimy1) || (zr < 0.0) || (zr >= 2) ) {
	   *(vl+j) = 0.0;
	  continue;
	}
	xi = (int)xr;           yi = (int)yr;           zi = (int)zr;
	dx = xr - (double)xi;   dy = yr - (double)yi;   dz = zr - (double)zi;
	/*--- interpolation lineaire : calcul des coefficients ---*/
	coef[0][0][0] = (1.0 - dx) * (1.0 - dy) * (1.0 - dz); /* xi,   yi,   zi   */
	coef[1][0][0] = dx         * (1.0 - dy) * (1.0 - dz); /* xi+1, yi,   zi   */
	coef[0][1][0] = (1.0 - dx) * dy         * (1.0 - dz); /* xi,   yi+1, zi   */
	coef[1][1][0] = dx         * dy         * (1.0 - dz); /* xi+1, yi+1, zi   */
	coef[0][0][1] = (1.0 - dx) * (1.0 - dy) * dz;         /* xi,   yi,   zi+1 */
	coef[1][0][1] = dx         * (1.0 - dy) * dz;         /* xi+1, yi,   zi+1 */
	coef[0][1][1] = (1.0 - dx) * dy         * dz;         /* xi,   yi+1, zi+1 */
	coef[1][1][1] = dx         * dy         * dz;         /* xi+1, yi+1, zi+1 */

	pos = xi + yi * dimx;
	norm  = *(norme[zi]+pos)        * coef[0][0][0] + *(norme[zi]+pos+1)        * coef[1][0][0];
	norm += *(norme[zi]+pos+dimx)   * coef[0][1][0] + *(norme[zi]+pos+1+dimx)   * coef[1][1][0];
	norm += *(norme[zi+1]+pos)      * coef[0][0][1] + *(norme[zi+1]+pos+1)      * coef[1][0][1];
	norm += *(norme[zi+1]+pos+dimx) * coef[0][1][1] + *(norme[zi+1]+pos+1+dimx) * coef[1][1][1];

	/*--- extrema ? ---*/
	if ( norme_originale < norm ) {
	    *(vl+j) = 0.0;
	    continue;
	}

	/*--- c'est un extrema ---*/
	*(vl+j) = (float)norme_originale;
    }
}

