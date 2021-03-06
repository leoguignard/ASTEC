




static int _Reech3DTriLinVectorField_DEFTYPE_TYPE ( void *parameter,
						    size_t first,
						    size_t last )
{
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef; 

  int i, j, k, l, ix, iy, iz;
  
  double xd, yd, zd, x, y, z;
  double dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  double v[3];
  double res;
  double v6, v5, v4;

  int rdimx=resDim[0], rdimy=resDim[1];

  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  int toffset1=tdimxy+tdimx+1, toffset2=tdimxy-tdimx-1;
  int t1dimx=tdimx-1, t1dimy=tdimy-1, t1dimz=tdimz-1;

  int ddimx=defDim[0], ddimy=defDim[1], ddimz=defDim[2];
  int ddimxy = ddimx * ddimy;
  int doffset1=ddimxy+ddimx+1, doffset2=ddimxy-ddimx-1;
  int d1dimx=ddimx-1, d1dimy=ddimy-1, d1dimz=ddimz-1;

  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  double borddimz = (double)tdimz-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;

  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* defz = theDef[2];
  DEFTYPE* dbuf;

  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;
  defz += first;
  
  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %d\r", k );

    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 )
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {
	
      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	zd = (double)k + *defz;
	defx ++;
	defy ++;
	defz ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[ 2] * k + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[ 6] * k + mat_bef[7];
	z = mat_bef[8] * i +  mat_bef[9] * j + mat_bef[10] * k + mat_bef[11];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	iz = (int)z;
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5)
	     || ( z <= -0.5 ) || ( z >= ddimz-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) &&
	     (z > 0.0) && (iz < d1dimz) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dz = z - iz;
	  dxdy = dx*dy;
	  dxdz = dx*dz;
	  dydz = dy*dz;
	  dxdydz = dxdy*dz;
	  
	  /* we have
	     v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
	     v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
	     v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
	     v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
	     v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v6 = dxdz-dxdydz;
	  v5 = dxdy-dxdydz;
	  v4 = dx-dxdy-v6;
	  
	  for ( l=0; l<3; l++ ) {
	    v[l] = 0;
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + iz * ddimxy + doffset1;
	    v[l] += dxdydz * (*dbuf);            /* tbuf(ix+1,iy+1,iz+1) */
	    dbuf --;
	    v[l] += (dydz-dxdydz) * (*dbuf);     /* tbuf(ix  ,iy+1,iz+1) */
	    dbuf -= d1dimx;
	    v[l] += v6 * (*dbuf);                /* tbuf(ix+1  ,iy,  iz+1) */
	    dbuf --;
	    v[l] += (dz-dydz-v6) * (*dbuf);      /* tbuf(ix  ,iy  ,iz+1) */
	    dbuf -= doffset2;
	    v[l] += v5 * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	    dbuf --;
	    v[l] += (dy-dydz-v5) * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf -= d1dimx;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf --;
	    v[l] += (1-dy-dz+dydz-v4) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	  }
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  
	  if ( (x < 0.0) || (ix == d1dimx) ) {
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = *dbuf;
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dz) * (*dbuf); /* (1-dz)* tbuf(ix,iy,iz) */
		  dbuf += ddimxy;
		  v[l] += dz * (*dbuf);     /* dz * tbuf(ix,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dy)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy*(1-dz) * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		  dbuf += doffset2+1;
		  v[l] += (1-dy)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf += ddimx;
		  v[l] += dy*dz * (*dbuf);        /* tbuf(ix,iy+1,iz+1) */
		}
	      }
	    }
	  }
	  else {
	    /* here we are sure that the border is either
	       along the Y or the Z axis */
	    dx = x - ix;
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx*(1-dz) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		  dbuf += ddimxy-1;
		  v[l] += (1-dx)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf ++;
		  v[l] += dx*dz * (*dbuf);        /* tbuf(ix+1,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<3; l++ ) {
		v[l] = 0;
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + iz * ddimxy;
		v[l] = (1-dx)*(1-dy) * (*dbuf); /* tbuf(ix,iy,iz) */
		dbuf ++;
		v[l] += dx*(1-dy) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		dbuf += d1dimx;
		v[l] += (1-dx)*dy * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		dbuf ++;
		v[l] += dx*dy * (*dbuf);        /* tbuf(ix+1,iy+1,iz) */
	      }
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	zd = z + v[2];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = zd;
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[2] * zd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[6] * zd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = mat_aft[8] * xd +  mat_aft[9] * yd + mat_aft[10] * zd + mat_aft[11];
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)x;
      iy = (int)y;
      iz = (int)z;
      tpt = (TYPE *)tbuf;
      
      /* are we on the border or not ? */
      if ( (x > 0.0) && (ix < t1dimx) &&
	   (y > 0.0) && (iy < t1dimy) &&
	   (z > 0.0) && (iz < t1dimz) ) {
	/* the corresponding point is in the box defined 
	   by (ix[+1],iy[+1],iz[+1]) */
	dx = x - ix;
	dy = y - iy;
	dz = z - iz;
	dxdy = dx*dy;
	dxdz = dx*dz;
	dydz = dy*dz;
	dxdydz = dxdy*dz;
	
	/* we have
	   v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
	   v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
	   v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
	   v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
	   v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
	   v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
	   v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
	   v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
	*/
	tpt += ix + iy * tdimx + iz * tdimxy + toffset1;
	v6 = dxdz-dxdydz;
	v5 = dxdy-dxdydz;
	v4 = dx-dxdy-v6;
	
	res = 0;
	res += dxdydz * (*tpt);            /* tbuf(ix+1,iy+1,iz+1) */
	tpt --;
	res += (dydz-dxdydz) * (*tpt);     /* tbuf(ix  ,iy+1,iz+1) */
	tpt -= t1dimx;
	res += v6 * (*tpt);                /* tbuf(ix+1  ,iy,  iz+1) */
	tpt --;
	res += (dz-dydz-v6) * (*tpt);      /* tbuf(ix  ,iy  ,iz+1) */
	tpt -= toffset2;
	res += v5 * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
	tpt --;
	res += (dy-dydz-v5) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
	tpt -= t1dimx;
	res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
	tpt --;
	res += (1-dy-dz+dydz-v4) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here, we are sure we are on some border */
      tpt += ix + iy * tdimx + iz * tdimxy;
      if ( (x < 0.0) || (ix == t1dimx) ) {
	if ( (y < 0.0) || (iy == t1dimy) ) {
	  if ( (z < 0.0) || (iz == t1dimz) ) {
	    *rbuf = *tpt;
	    continue;
	  }
	  dz = z - iz;
	  res  = (1-dz) * (*tpt); /* (1-dz)* tbuf(ix,iy,iz) */
	  tpt += tdimxy;
	  res += dz * (*tpt);     /* dz * tbuf(ix,iy,iz+1) */
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dy = y - iy;
	if ( (z < 0.0) || (iz == t1dimz) ) {
	  res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy,iz) */
	  tpt += tdimx;
	  res += dy * (*tpt);     /* dy * tbuf(ix,iy+1,iz) */
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dz = z - iz;
	res = (1-dy)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
	tpt += tdimx;
	res += dy*(1-dz) * (*tpt);    /* tbuf(ix,iy+1,iz) */
	tpt += toffset2+1;
	res += (1-dy)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
	tpt += tdimx;
	res += dy*dz * (*tpt);        /* tbuf(ix,iy+1,iz+1) */
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here we are sure that the border is either
	 along the Y or the Z axis */
      dx = x - ix;
      if ( (y < 0.0) || (iy == t1dimy) ) {
	if ( (z < 0.0) || (iz == t1dimz) ) {
	  res = (1-dx) * (*tpt); /* tbuf(ix,iy,iz) */
	  tpt ++;
	  res += dx * (*tpt);    /* tbuf(ix+1,iy,iz) */
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dz = z - iz;
	res = (1-dx)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
	tpt ++;
	res += dx*(1-dz) * (*tpt);    /* tbuf(ix+1,iy,iz) */
	tpt += tdimxy-1;
	res += (1-dx)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
	tpt ++;
	res += dx*dz * (*tpt);        /* tbuf(ix+1,iy,iz+1) */
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here we are sure that the border is along the Z axis */
      dy = y - iy;
      res = (1-dx)*(1-dy) * (*tpt); /* tbuf(ix,iy,iz) */
      tpt ++;
      res += dx*(1-dy) * (*tpt);    /* tbuf(ix+1,iy,iz) */
      tpt += t1dimx;
      res += (1-dx)*dy * (*tpt);    /* tbuf(ix,iy+1,iz) */
      tpt ++;
      res += dx*dy * (*tpt);        /* tbuf(ix+1,iy+1,iz) */
      *rbuf = (TYPE)_CONVERT_( res );
    }
  }
  return( 1 );
}





extern void Reech3DTriLinVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
			     int *theDim,  /* dimensions of this buffer */
			     void* resBuf, /* result buffer */
			     int *resDim,  /* dimensions of this buffer */
			     DEFTYPE** theDef, /* deformations */
                             int *defDim, /* dimensions of these buffers */
			     double* mat_aft,  /* transformation matrix */
			     double* mat_bef  /* transformation matrix */
			     )
{
  char *proc = "Reech3DTriLinVectorField_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = 1.0;
  p.bias = 0.0;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech3DTriLinVectorField_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}





/* Resampling procedure.

   Work for 3D images, not for vectorial ones.

   Cette procedure utilise un champ de deformations 
   backward, c'est-a-dire allant de l'image a reechantillonner
   a l'image reference. Une fois estime le point apres application
   du champ de deformation, on applique a celui-ci une matrice 
   de transformation.

   Entrees:
   - le buffer de l'image reference
   - une coupe (flottante) de l'image resultat
   - les 3 coupes de deformations (le vecteur) associees
   - la matrice

*/

static int _Reech3DTriLinVectorFieldgb_DEFTYPE_TYPE ( void *parameter,
						      size_t first,
						      size_t last )
{
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef;
  float gain = p->gain; 
  float bias = p->bias; 

  int i, j, k, l, ix, iy, iz;
  
  double xd, yd, zd, x, y, z;
  double dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  double v[3];
  double res;
  double v6, v5, v4;

  int rdimx=resDim[0], rdimy=resDim[1];

  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  int toffset1=tdimxy+tdimx+1, toffset2=tdimxy-tdimx-1;
  int t1dimx=tdimx-1, t1dimy=tdimy-1, t1dimz=tdimz-1;

  int ddimx=defDim[0], ddimy=defDim[1], ddimz=defDim[2];
  int ddimxy = ddimx * ddimy;
  int doffset1=ddimxy+ddimx+1, doffset2=ddimxy-ddimx-1;
  int d1dimx=ddimx-1, d1dimy=ddimy-1, d1dimz=ddimz-1;

  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  double borddimz = (double)tdimz-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;

  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* defz = theDef[2];
  DEFTYPE* dbuf;

  register double b=bias;
  register double g=gain;

  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;
  defz += first;

  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %d\r", k );

    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 )
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {
	
      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	zd = (double)k + *defz;
	defx ++;
	defy ++;
	defz ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[ 2] * k + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[ 6] * k + mat_bef[7];
	z = mat_bef[8] * i +  mat_bef[9] * j + mat_bef[10] * k + mat_bef[11];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	iz = (int)z;
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5)
	     || ( z <= -0.5 ) || ( z >= ddimz-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) &&
	     (z > 0.0) && (iz < d1dimz) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dz = z - iz;
	  dxdy = dx*dy;
	  dxdz = dx*dz;
	  dydz = dy*dz;
	  dxdydz = dxdy*dz;
	  
	  /* we have
	     v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
	     v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
	     v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
	     v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
	     v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v6 = dxdz-dxdydz;
	  v5 = dxdy-dxdydz;
	  v4 = dx-dxdy-v6;
	  
	  for ( l=0; l<3; l++ ) {
	    v[l] = 0;
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + iz * ddimxy + doffset1;
	    v[l] += dxdydz * (*dbuf);            /* tbuf(ix+1,iy+1,iz+1) */
	    dbuf --;
	    v[l] += (dydz-dxdydz) * (*dbuf);     /* tbuf(ix  ,iy+1,iz+1) */
	    dbuf -= d1dimx;
	    v[l] += v6 * (*dbuf);                /* tbuf(ix+1  ,iy,  iz+1) */
	    dbuf --;
	    v[l] += (dz-dydz-v6) * (*dbuf);      /* tbuf(ix  ,iy  ,iz+1) */
	    dbuf -= doffset2;
	    v[l] += v5 * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	    dbuf --;
	    v[l] += (dy-dydz-v5) * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf -= d1dimx;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf --;
	    v[l] += (1-dy-dz+dydz-v4) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	  }
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  
	  if ( (x < 0.0) || (ix == d1dimx) ) {
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = *dbuf;
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dz) * (*dbuf); /* (1-dz)* tbuf(ix,iy,iz) */
		  dbuf += ddimxy;
		  v[l] += dz * (*dbuf);     /* dz * tbuf(ix,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dy)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy*(1-dz) * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		  dbuf += doffset2+1;
		  v[l] += (1-dy)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf += ddimx;
		  v[l] += dy*dz * (*dbuf);        /* tbuf(ix,iy+1,iz+1) */
		}
	      }
	    }
	  }
	  else {
	    /* here we are sure that the border is either
	       along the Y or the Z axis */
	    dx = x - ix;
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx*(1-dz) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		  dbuf += ddimxy-1;
		  v[l] += (1-dx)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf ++;
		  v[l] += dx*dz * (*dbuf);        /* tbuf(ix+1,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<3; l++ ) {
		v[l] = 0;
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + iz * ddimxy;
		v[l] = (1-dx)*(1-dy) * (*dbuf); /* tbuf(ix,iy,iz) */
		dbuf ++;
		v[l] += dx*(1-dy) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		dbuf += d1dimx;
		v[l] += (1-dx)*dy * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		dbuf ++;
		v[l] += dx*dy * (*dbuf);        /* tbuf(ix+1,iy+1,iz) */
	      }
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	zd = z + v[2];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = zd;
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[2] * zd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[6] * zd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = mat_aft[8] * xd +  mat_aft[9] * yd + mat_aft[10] * zd + mat_aft[11];
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)x;
      iy = (int)y;
      iz = (int)z;
      tpt = (TYPE *)tbuf;
      
      /* are we on the border or not ? */
      if ( (x > 0.0) && (ix < t1dimx) &&
	   (y > 0.0) && (iy < t1dimy) &&
	   (z > 0.0) && (iz < t1dimz) ) {
	/* the corresponding point is in the box defined 
	   by (ix[+1],iy[+1],iz[+1]) */
	dx = x - ix;
	dy = y - iy;
	dz = z - iz;
	dxdy = dx*dy;
	dxdz = dx*dz;
	dydz = dy*dz;
	dxdydz = dxdy*dz;
	
	/* we have
	   v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
	   v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
	   v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
	   v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
	   v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
	   v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
	   v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
	   v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
	*/
	tpt += ix + iy * tdimx + iz * tdimxy + toffset1;
	v6 = dxdz-dxdydz;
	v5 = dxdy-dxdydz;
	v4 = dx-dxdy-v6;
	
	res = 0;
	res += dxdydz * (*tpt);            /* tbuf(ix+1,iy+1,iz+1) */
	tpt --;
	res += (dydz-dxdydz) * (*tpt);     /* tbuf(ix  ,iy+1,iz+1) */
	tpt -= t1dimx;
	res += v6 * (*tpt);                /* tbuf(ix+1  ,iy,  iz+1) */
	tpt --;
	res += (dz-dydz-v6) * (*tpt);      /* tbuf(ix  ,iy  ,iz+1) */
	tpt -= toffset2;
	res += v5 * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
	tpt --;
	res += (dy-dydz-v5) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
	tpt -= t1dimx;
	res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
	tpt --;
	res += (1-dy-dz+dydz-v4) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
	res = res * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here, we are sure we are on some border */
      tpt += ix + iy * tdimx + iz * tdimxy;
      if ( (x < 0.0) || (ix == t1dimx) ) {
	if ( (y < 0.0) || (iy == t1dimy) ) {
	  if ( (z < 0.0) || (iz == t1dimz) ) {
	    res = (double)(*tpt) * g + b;
	    *rbuf = (TYPE)_CONVERT_( res );
	    continue;
	  }
	  dz = z - iz;
	  res  = (1-dz) * (*tpt); /* (1-dz)* tbuf(ix,iy,iz) */
	  tpt += tdimxy;
	  res += dz * (*tpt);     /* dz * tbuf(ix,iy,iz+1) */
	  res = res * g + b;
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dy = y - iy;
	if ( (z < 0.0) || (iz == t1dimz) ) {
	  res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy,iz) */
	  tpt += tdimx;
	  res += dy * (*tpt);     /* dy * tbuf(ix,iy+1,iz) */
	  res = res * g + b;
	 *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dz = z - iz;
	res = (1-dy)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
	tpt += tdimx;
	res += dy*(1-dz) * (*tpt);    /* tbuf(ix,iy+1,iz) */
	tpt += toffset2+1;
	res += (1-dy)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
	tpt += tdimx;
	res += dy*dz * (*tpt);        /* tbuf(ix,iy+1,iz+1) */
	res = res * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here we are sure that the border is either
	 along the Y or the Z axis */
      dx = x - ix;
      if ( (y < 0.0) || (iy == t1dimy) ) {
	if ( (z < 0.0) || (iz == t1dimz) ) {
	  res = (1-dx) * (*tpt); /* tbuf(ix,iy,iz) */
	  tpt ++;
	  res += dx * (*tpt);    /* tbuf(ix+1,iy,iz) */
	  res = res * g + b;
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dz = z - iz;
	res = (1-dx)*(1-dz) * (*tpt); /* tbuf(ix,iy,iz) */
	tpt ++;
	res += dx*(1-dz) * (*tpt);    /* tbuf(ix+1,iy,iz) */
	tpt += tdimxy-1;
	res += (1-dx)*dz * (*tpt);    /* tbuf(ix,iy,iz+1) */
	tpt ++;
	res += dx*dz * (*tpt);        /* tbuf(ix+1,iy,iz+1) */
	res = res * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here we are sure that the border is along the Z axis */
      dy = y - iy;
      res = (1-dx)*(1-dy) * (*tpt); /* tbuf(ix,iy,iz) */
      tpt ++;
      res += dx*(1-dy) * (*tpt);    /* tbuf(ix+1,iy,iz) */
      tpt += t1dimx;
      res += (1-dx)*dy * (*tpt);    /* tbuf(ix,iy+1,iz) */
      tpt ++;
      res += dx*dy * (*tpt);        /* tbuf(ix+1,iy+1,iz) */
      res = res * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
    }
  }
  return( 1 );
}





extern void Reech3DTriLinVectorFieldgb_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
						      int *theDim,  /* dimensions of this buffer */
						      void* resBuf, /* result buffer */
						      int *resDim,  /* dimensions of this buffer */
						      DEFTYPE** theDef, /* deformations */
						      int *defDim, /* dimensions of these buffers */
						      double* mat_aft,  /* transformation matrix */
						      double* mat_bef,  /* transformation matrix */
						      float gain,
						      float bias )
{
  char *proc = "Reech3DTriLinVectorFieldgb_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = gain;
  p.bias = bias;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  fprintf( stderr, "%s: %d chunks\n", proc, chunks.n_allocated_chunks );

  /* processing
   */
  if ( processChunks( &_Reech3DTriLinVectorFieldgb_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}





static int _Reech3DNearestVectorField_DEFTYPE_TYPE ( void *parameter,
						    size_t first,
						    size_t last )
{
  char *proc = "_Reech3DNearestVectorField_DEFTYPE_TYPE";
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef; 

  int i, j, k, l, ix, iy, iz;
  
  double xd, yd, zd, x, y, z;
  double dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  double v[3];
  double v6, v5, v4;
  
  int rdimx=resDim[0], rdimy=resDim[1];
  
  int tdimx=theDim[0], tdimy=theDim[1], tdimz=theDim[2];
  int tdimxy=tdimx*tdimy;
  
  int ddimx=defDim[0], ddimy=defDim[1], ddimz=defDim[2];
  int ddimxy = ddimx * ddimy;
  int doffset1=ddimxy+ddimx+1, doffset2=ddimxy-ddimx-1;
  int d1dimx=ddimx-1, d1dimy=ddimy-1, d1dimz=ddimz-1;
  
  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  double borddimz = (double)tdimz-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *rbuf = (TYPE*)resBuf;
  
  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* defz = theDef[2];
  DEFTYPE* dbuf;
  
  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;
  defz += first;

  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  if ( _verbose_ > 2 ) {
    fprintf( stderr, "%s: Processing image\n", proc );
    fprintf( stderr, "   from point #%lu = (%d,%d,%d) to point #%lu = (%d,%d,%d)\n", 
	     first, ifirst, jfirst, kfirst,
	     last, ilast, jlast, klast );
  }
  if ( _verbose_ > 3 ) {
    if ( mat_bef == NULL )  fprintf( stderr, "%s: no 'before' matrix\n" ,proc );
    if ( mat_aft == NULL )  fprintf( stderr, "%s: no 'after' matrix\n" ,proc );
  }


  
  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "%s: Processing slice %d\n", proc, k );
    
    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 ) {
    if ( _verbose_ > 4 ) {
	if ( k == 512 ) {
	  fprintf( stderr, "%s: Processing line %d\n", proc, j );
	}
      }
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {

      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	zd = (double)k + *defz;
	defx ++;
	defy ++;
	defz ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[ 2] * k + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[ 6] * k + mat_bef[7];
	z = mat_bef[8] * i +  mat_bef[9] * j + mat_bef[10] * k + mat_bef[11];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	iz = (int)z;
	
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5)
	     || ( z <= -0.5 ) || ( z >= ddimz-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) &&
	     (z > 0.0) && (iz < d1dimz) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dz = z - iz;
	  dxdy = dx*dy;
	  dxdz = dx*dz;
	  dydz = dy*dz;
	  dxdydz = dxdy*dz;
	  
	  /* we have
	     v[7]=dxdydz;            coefficient of tbuf(ix+1,iy+1,iz+1)
	     v[6]=dxdz-dxdydz;       coefficient of tbuf(ix+1,iy,  iz+1)
	     v[5]=dxdy-dxdydz;       coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy-v[6];      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[3]=dydz-dxdydz;       coefficient of tbuf(ix  ,iy+1,iz+1)
	     v[2]=dz-dydz-v[6];      coefficient of tbuf(ix  ,iy  ,iz+1)
	     v[1]=dy-dydz-v[5];      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dy-dz+dydz-v[4]; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v6 = dxdz-dxdydz;
	  v5 = dxdy-dxdydz;
	  v4 = dx-dxdy-v6;
	  
	  for ( l=0; l<3; l++ ) {
	    v[l] = 0;
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + iz * ddimxy + doffset1;
	    v[l] += dxdydz * (*dbuf);            /* tbuf(ix+1,iy+1,iz+1) */
	    dbuf --;
	    v[l] += (dydz-dxdydz) * (*dbuf);     /* tbuf(ix  ,iy+1,iz+1) */
	    dbuf -= d1dimx;
	    v[l] += v6 * (*dbuf);                /* tbuf(ix+1  ,iy,  iz+1) */
	    dbuf --;
	    v[l] += (dz-dydz-v6) * (*dbuf);      /* tbuf(ix  ,iy  ,iz+1) */
	    dbuf -= doffset2;
	    v[l] += v5 * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	    dbuf --;
	    v[l] += (dy-dydz-v5) * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf -= d1dimx;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf --;
	    v[l] += (1-dy-dz+dydz-v4) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	  }
	  
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  
	  if ( (x < 0.0) || (ix == d1dimx) ) {
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = *dbuf;
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dz) * (*dbuf); /* (1-dz)* tbuf(ix,iy,iz) */
		  dbuf += ddimxy;
		  v[l] += dz * (*dbuf);     /* dz * tbuf(ix,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dy)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf += ddimx;
		  v[l] += dy*(1-dz) * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		  dbuf += doffset2+1;
		  v[l] += (1-dy)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf += ddimx;
		  v[l] += dy*dz * (*dbuf);        /* tbuf(ix,iy+1,iz+1) */
		}
	      }
	    }
	  }
	  else {
	    /* here we are sure that the border is either
	       along the Y or the Z axis */
	    dx = x - ix;
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      if ( (z < 0.0) || (iz == d1dimz) ) {
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		}
	      }
	      else {
		dz = z - iz;
		for ( l=0; l<3; l++ ) {
		  v[l] = 0;
		  dbuf = theDef[l];
		  dbuf += ix + iy * ddimx + iz * ddimxy;
		  v[l] = (1-dx)*(1-dz) * (*dbuf); /* tbuf(ix,iy,iz) */
		  dbuf ++;
		  v[l] += dx*(1-dz) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		  dbuf += ddimxy-1;
		  v[l] += (1-dx)*dz * (*dbuf);    /* tbuf(ix,iy,iz+1) */
		  dbuf ++;
		  v[l] += dx*dz * (*dbuf);        /* tbuf(ix+1,iy,iz+1) */
		}
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<3; l++ ) {
		v[l] = 0;
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + iz * ddimxy;
		v[l] = (1-dx)*(1-dy) * (*dbuf); /* tbuf(ix,iy,iz) */
		dbuf ++;
		v[l] += dx*(1-dy) * (*dbuf);    /* tbuf(ix+1,iy,iz) */
		dbuf += d1dimx;
		v[l] += (1-dx)*dy * (*dbuf);    /* tbuf(ix,iy+1,iz) */
		dbuf ++;
		v[l] += dx*dy * (*dbuf);        /* tbuf(ix+1,iy+1,iz) */
	      }
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	zd = z + v[2];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = zd;
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[2] * zd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[6] * zd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
	z = mat_aft[8] * xd +  mat_aft[9] * yd + mat_aft[10] * zd + mat_aft[11];
	if (( z <= -0.5 ) || ( z >= borddimz)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)(x+0.5);
      iy = (int)(y+0.5);
      iz = (int)(z+0.5);
      
      *rbuf = tbuf[ ix + iy * tdimx + iz * tdimxy ];
    }
    }
  }
  return( 1 );
}





extern void Reech3DNearestVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
						     int *theDim,  /* dimensions of this buffer */
						     void* resBuf, /* result buffer */
						     int *resDim,  /* dimensions of this buffer */
						     DEFTYPE** theDef, /* deformations */
						     int *defDim, /* dimensions of these buffers */
						     double* mat_aft,  /* transformation matrix */
						     double* mat_bef  /* transformation matrix */
						     )
{
  char *proc = "Reech3DNearestVectorField_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  if ( _verbose_ >= 2 ) {
    fprintf( stderr, "%s: has allocated %d chunks\n", proc, chunks.n_allocated_chunks );
    fprintf( stderr, "\t image indexes are in %lu - %lu\n", first, last ); 
  }


  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = 1.0;
  p.bias = 0.0;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech3DNearestVectorField_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}





static int _Reech2DTriLinVectorField_DEFTYPE_TYPE ( void *parameter,
						    size_t first,
						    size_t last )
{
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef; 

  int i, j, k, l, ix, iy;
  
  double xd, yd, x, y;
  double dx, dy, dxdy;
  double v[2];
  double res;
  double v1, v4;

  int rdimx=resDim[0], rdimy=resDim[1];

  int tdimx=theDim[0], tdimy=theDim[1];
  int tdimxy=tdimx*tdimy;
  int t1dimx=tdimx-1, t1dimy=tdimy-1;

  int ddimx=defDim[0], ddimy=defDim[1];
  int ddimxy = ddimx * ddimy;
  int d1dimx=ddimx-1, d1dimy=ddimy-1;

  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;

  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* dbuf;

  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;

  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %d\r", k );

    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 )
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {
	
      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	defx ++;
	defy ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[7];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dxdy = dx*dy;
	  
	  /* we have
	     v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v1 = dy-dxdy;
	  v4 = dx-dxdy;

	  for ( l=0; l<2; l++ ) {
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + k * ddimxy;
	    v[l] = 0;
	    v[l] += (1-dx-v1) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	    dbuf ++;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf += d1dimx;
	    v[l] += v1 * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf ++;
	    v[l] += dxdy * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	  }
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  if ( (x < 0.0) || (ix == d1dimx) ) {

	    /* we just look at y */
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l] = *dbuf;
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy) */
		dbuf += d1dimx;
		v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1) */
	      }
	    }
	  }
	  else {
	    dx = x - ix;
	    for ( l=0; l<2; l++ ) {
	      dbuf = theDef[l];
	      dbuf += ix + iy * ddimx + k * ddimxy;
	      v[l]  = (1-dx) * (*dbuf); /* (1-dx)* tbuf(ix,iy) */
	      dbuf ++;
	      v[l] += dx * (*dbuf);     /* dx * tbuf(ix+1,iy) */
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)x;
      iy = (int)y;
      tpt = (TYPE *)tbuf;
      
      /* are we on the border or not ? */
      if ( (x > 0.0) && (ix < t1dimx) &&
	   (y > 0.0) && (iy < t1dimy) ) {
	/* the corresponding point is in the box defined 
	   by (ix[+1],iy[+1],iz[+1]) */
	dx = x - ix;
	dy = y - iy;
	dxdy = dx*dy;
	
	/* we have
	   v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1,iz  )
	   v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  ,iz  )
	   v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1,iz  )
	   v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  ,iz  )
	*/
	tpt += ix + iy * tdimx + k * tdimxy;

	v1 = dy-dxdy;
	v4 = dx-dxdy;
	
	res = 0;
	res += (1-dx-v1) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
	tpt ++;
	res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
	tpt += t1dimx;
	res += (v1) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
	tpt ++;
	res += dxdy * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here, we are sure we are on some border */
      if ( (x < 0.0) || (ix == t1dimx) ) {
	/* we just look at y */
	if ( (y < 0.0) || (iy == t1dimy) ) {
	  *rbuf = *tpt;
	  continue;
	}
	dy = y - iy;
	res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy) */
	tpt += tdimx;
	res += dy * (*tpt);     /* dy * tbuf(ix,iy+1) */
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      dx = x - ix;
      res  = (1-dx) * (*tpt); /* (1-dx)* tbuf(ix,iy) */
      tpt ++;
      res += dx * (*tpt);     /* dx * tbuf(ix+1,iy) */
      *rbuf = (TYPE)_CONVERT_( res );
    }
  }
  return( 1 );
}





extern void Reech2DTriLinVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
						    int *theDim,  /* dimensions of this buffer */
						    void* resBuf, /* result buffer */
						    int *resDim,  /* dimensions of this buffer */
						    DEFTYPE** theDef, /* deformations */
						    int *defDim, /* dimensions of these buffers */
						    double* mat_aft,  /* transformation matrix */
						    double* mat_bef  /* transformation matrix */
						    )
{
  char *proc = "Reech2DTriLinVectorField_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = 1.0;
  p.bias = 0.0;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech2DTriLinVectorField_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}





static int _Reech2DTriLinVectorFieldgb_DEFTYPE_TYPE ( void *parameter,
						      size_t first,
						      size_t last )
{
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef;
  float gain = p->gain; 
  float bias = p->bias; 

  int i, j, k, l, ix, iy;
  
  double xd, yd, x, y;
  double dx, dy, dxdy;
  double v[2];
  double res;
  double v1, v4;

  int rdimx=resDim[0], rdimy=resDim[1];

  int tdimx=theDim[0], tdimy=theDim[1];
  int tdimxy=tdimx*tdimy;
  int t1dimx=tdimx-1, t1dimy=tdimy-1;

  int ddimx=defDim[0], ddimy=defDim[1];
  int ddimxy = ddimx * ddimy;
  int d1dimx=ddimx-1, d1dimy=ddimy-1;

  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  TYPE *tpt;
  register TYPE *rbuf = (TYPE*)resBuf;

  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* dbuf;

  register double b=bias;
  register double g=gain;

  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;

  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %d\r", k );

    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 )
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {
	
      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	defx ++;
	defy ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[7];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dxdy = dx*dy;
	  
	  /* we have
	     v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v1 = dy-dxdy;
	  v4 = dx-dxdy;

	  for ( l=0; l<2; l++ ) {
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + k * ddimxy;
	    v[l] = 0;
	    v[l] += (1-dx-v1) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	    dbuf ++;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf += d1dimx;
	    v[l] += v1 * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf ++;
	    v[l] += dxdy * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	  }
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  if ( (x < 0.0) || (ix == d1dimx) ) {

	    /* we just look at y */
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l] = *dbuf;
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy) */
		dbuf += d1dimx;
		v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1) */
	      }
	    }
	  }
	  else {
	    dx = x - ix;
	    for ( l=0; l<2; l++ ) {
	      dbuf = theDef[l];
	      dbuf += ix + iy * ddimx + k * ddimxy;
	      v[l]  = (1-dx) * (*dbuf); /* (1-dx)* tbuf(ix,iy) */
	      dbuf ++;
	      v[l] += dx * (*dbuf);     /* dx * tbuf(ix+1,iy) */
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)x;
      iy = (int)y;
      tpt = (TYPE *)tbuf;
      
      /* are we on the border or not ? */
      if ( (x > 0.0) && (ix < t1dimx) &&
	   (y > 0.0) && (iy < t1dimy) ) {
	/* the corresponding point is in the box defined 
	   by (ix[+1],iy[+1],iz[+1]) */
	dx = x - ix;
	dy = y - iy;
	dxdy = dx*dy;
	
	/* we have
	   v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1,iz  )
	   v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  ,iz  )
	   v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1,iz  )
	   v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  ,iz  )
	*/
	tpt += ix + iy * tdimx + k * tdimxy;

	v1 = dy-dxdy;
	v4 = dx-dxdy;
	
	res = 0;
	res += (1-dx-v1) * (*tpt); /* tbuf(ix  ,iy  ,iz  ) */
	tpt ++;
	res += v4 * (*tpt);                /* tbuf(ix+1,iy  ,iz  ) */
	tpt += t1dimx;
	res += (v1) * (*tpt);      /* tbuf(ix  ,iy+1,iz  ) */
	tpt ++;
	res += dxdy * (*tpt);                /* tbuf(ix+1,iy+1,iz  ) */
	res = res * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      /* here, we are sure we are on some border */
      if ( (x < 0.0) || (ix == t1dimx) ) {
	/* we just look at y */
	if ( (y < 0.0) || (iy == t1dimy) ) {
	  res = (double)(*tpt) * g + b;
	  *rbuf = (TYPE)_CONVERT_( res );
	  continue;
	}
	dy = y - iy;
	res  = (1-dy) * (*tpt); /* (1-dy)* tbuf(ix,iy) */
	tpt += tdimx;
	res += dy * (*tpt);     /* dy * tbuf(ix,iy+1) */
	res = (double)(*tpt) * g + b;
	*rbuf = (TYPE)_CONVERT_( res );
	continue;
      }
      dx = x - ix;
      res  = (1-dx) * (*tpt); /* (1-dx)* tbuf(ix,iy) */
      tpt ++;
      res += dx * (*tpt);     /* dx * tbuf(ix+1,iy) */
      res = res * g + b;
      *rbuf = (TYPE)_CONVERT_( res );
    }
  }
  return( 1 );
}





extern void Reech2DTriLinVectorFieldgb_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
						      int *theDim,  /* dimensions of this buffer */
						      void* resBuf, /* result buffer */
						      int *resDim,  /* dimensions of this buffer */
						      DEFTYPE** theDef, /* deformations */
						      int *defDim, /* dimensions of these buffers */
						      double* mat_aft,  /* transformation matrix */
						      double* mat_bef,  /* transformation matrix */
						      float gain,
						      float bias )
{
  char *proc = "Reech2DTriLinVectorFieldgb_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = gain;
  p.bias = bias;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  fprintf( stderr, "%s: %d chunks\n", proc, chunks.n_allocated_chunks );

  /* processing
   */
  if ( processChunks( &_Reech2DTriLinVectorFieldgb_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}





static int _Reech2DNearestVectorField_DEFTYPE_TYPE ( void *parameter,
						    size_t first,
						    size_t last )
{
  _VectorFieldResamplingParam *p = (_VectorFieldResamplingParam *)parameter;

  void* theBuf = p->theBuf;
  int *theDim = p->theDim;
  void* resBuf = p->resBuf;
  int *resDim = p->resDim;
  DEFTYPE** theDef = (DEFTYPE**)p->theDef;
  int *defDim = p->defDim;
  double* mat_aft = p->mat_aft; 
  double* mat_bef = p->mat_bef; 

  int i, j, k, l, ix, iy;
  
  double xd, yd, x, y;
  double dx, dy, dxdy;
  double v[2];
  double v1, v4;

  int rdimx=resDim[0], rdimy=resDim[1];

  int tdimx=theDim[0], tdimy=theDim[1];
  int tdimxy=tdimx*tdimy;

  int ddimx=defDim[0], ddimy=defDim[1];
  int ddimxy = ddimx * ddimy;
  int d1dimx=ddimx-1, d1dimy=ddimy-1;

  double borddimx = (double)tdimx-0.5, borddimy = (double)tdimy-0.5;
  TYPE *tbuf = (TYPE*)theBuf;
  register TYPE *rbuf = (TYPE*)resBuf;

  DEFTYPE* defx = theDef[0];
  DEFTYPE* defy = theDef[1];
  DEFTYPE* dbuf;

  int ifirst, jfirst, kfirst;
  int ilast, jlast, klast;

  k = kfirst = first / (rdimx*rdimy);
  j = jfirst = (first - kfirst*(rdimx*rdimy)) / rdimx;
  i = ifirst = (first - kfirst*(rdimx*rdimy) - jfirst*rdimx);

  klast = last / (rdimx*rdimy);
  jlast = (last - klast*(rdimx*rdimy)) / rdimx;
  ilast = (last - klast*(rdimx*rdimy) - jlast*rdimx);

  rbuf += first;
  defx += first;
  defy += first;

  if ( mat_bef == NULL ) {
    if ( resDim[0] != defDim[0] || resDim[1] != defDim[1] || resDim[2] != defDim[2] ) {
      fprintf( stderr, "deformation field should have the same dimension than result image\n" );
      return( -1 );
    }
  }

  for ( ; k<=klast; k++, j=0 ) {
    if ( _verbose_ > 1 )
      fprintf( stderr, "Processing slice %d\r", k );

    for ( ; (j<rdimy && k<klast) || (j<=jlast && k==klast); j++, i=0 )
    for ( ; (i<rdimx && (k<klast || (j<jlast && k==klast))) || (i<=ilast && j==jlast && k==klast); i++, rbuf++ ) {
	
      /* computation of the corresponding point after deformation */
      if ( mat_bef == NULL ) {
	xd = (double)i + *defx;
	yd = (double)j + *defy;
	defx ++;
	defy ++;
      }
      else {
	
	/* apply the first matrix */
	x = mat_bef[0] * i +  mat_bef[1] * j + mat_bef[3];
	y = mat_bef[4] * i +  mat_bef[5] * j + mat_bef[7];
	
	/* interpolate the vector deformation at (xd,yd,zd) */
	ix = (int)x;
	iy = (int)y;
	
	/* the point is outside the deformation field 
	 */
	if ( ( x <= -0.5 ) || ( x >= ddimx-0.5)
	     || ( y <= -0.5 ) || ( y >= ddimy-0.5) ) {
	  *rbuf = 0; 
	  continue; 
	}
	
	/* vector interpolation: are we on the border or not ? */
	if ( (x > 0.0) && (ix < d1dimx) &&
	     (y > 0.0) && (iy < d1dimy) ) {
	  /* the corresponding point is in the box defined 
	     by (ix[+1],iy[+1],iz[+1]) */
	  dx = x - ix;
	  dy = y - iy;
	  dxdy = dx*dy;
	  
	  /* we have
	     v[5]=dxdy;         coefficient of tbuf(ix+1,iy+1,iz  )
	     v[4]=dx-dxdy;      coefficient of tbuf(ix+1,iy  ,iz  )
	     v[1]=dy-dxdy;      coefficient of tbuf(ix  ,iy+1,iz  )
	     v[0]=1-dx-dy+dxdy; coefficient of tbuf(ix  ,iy  ,iz  )
	  */
	  
	  v1 = dy-dxdy;
	  v4 = dx-dxdy;

	  for ( l=0; l<2; l++ ) {
	    dbuf = theDef[l];
	    dbuf += ix + iy * ddimx + k * ddimxy;
	    v[l] = 0;
	    v[l] += (1-dx-v1) * (*dbuf); /* tbuf(ix  ,iy  ,iz  ) */
	    dbuf ++;
	    v[l] += v4 * (*dbuf);                /* tbuf(ix+1,iy  ,iz  ) */
	    dbuf += d1dimx;
	    v[l] += v1 * (*dbuf);      /* tbuf(ix  ,iy+1,iz  ) */
	    dbuf ++;
	    v[l] += dxdy * (*dbuf);                /* tbuf(ix+1,iy+1,iz  ) */
	  }
	}
	else {
	  
	  /* here, we are sure we are on some border */
	  if ( (x < 0.0) || (ix == d1dimx) ) {

	    /* we just look at y */
	    if ( (y < 0.0) || (iy == d1dimy) ) {
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l] = *dbuf;
	      }
	    }
	    else {
	      dy = y - iy;
	      for ( l=0; l<2; l++ ) {
		dbuf = theDef[l];
		dbuf += ix + iy * ddimx + k * ddimxy;
		v[l]  = (1-dy) * (*dbuf); /* (1-dy)* tbuf(ix,iy) */
		dbuf += d1dimx;
		v[l] += dy * (*dbuf);     /* dy * tbuf(ix,iy+1) */
	      }
	    }
	  }
	  else {
	    dx = x - ix;
	    for ( l=0; l<2; l++ ) {
	      dbuf = theDef[l];
	      dbuf += ix + iy * ddimx + k * ddimxy;
	      v[l]  = (1-dx) * (*dbuf); /* (1-dx)* tbuf(ix,iy) */
	      dbuf ++;
	      v[l] += dx * (*dbuf);     /* dx * tbuf(ix+1,iy) */
	    }
	  }
	}
	
	xd = x + v[0];
	yd = y + v[1];
	
      } /* ( mat_bef != NULL ) */
      
      /* computation of the corresponding point after matrix application */
      if ( mat_aft == NULL ) {
	x = xd;
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = yd;
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      else {
	x = mat_aft[0] * xd +  mat_aft[1] * yd + mat_aft[3];
	if (( x <= -0.5 ) || ( x >= borddimx)) { *rbuf = 0; continue; }
	y = mat_aft[4] * xd +  mat_aft[5] * yd + mat_aft[7];
	if (( y <= -0.5 ) || ( y >= borddimy)) { *rbuf = 0; continue; }
      }
      
      /* here, the point lies on the borders or completely inside
	 the image */
      ix = (int)(x+0.5);
      iy = (int)(y+0.5);
      
      *rbuf = tbuf[ ix + iy * tdimx + k * tdimxy ];
      
    }
  }
  return( 1 );
}





extern void Reech2DNearestVectorField_DEFTYPE_TYPE ( void* theBuf, /* buffer to be resampled */
						     int *theDim,  /* dimensions of this buffer */
						     void* resBuf, /* result buffer */
						     int *resDim,  /* dimensions of this buffer */
						     DEFTYPE** theDef, /* deformations */
						     int *defDim, /* dimensions of these buffers */
						     double* mat_aft,  /* transformation matrix */
						     double* mat_bef  /* transformation matrix */
						     )
{
  char *proc = "Reech2DNearestVectorField_DEFTYPE_TYPE";
  size_t first = 0;
  size_t last;
  int i;
  typeChunks chunks;
  _VectorFieldResamplingParam p;

  /* preparing parallelism
   */
  first = 0;
  last = (size_t)resDim[2] * (size_t)resDim[1] * (size_t)resDim[0] - 1;
  initChunks( &chunks );
  if ( buildChunks( &chunks, first, last, proc ) != 1 ) {
    if ( _verbose_ )
      fprintf( stderr, "%s: unable to compute chunks\n", proc );
    return;
  }

  p.theBuf = theBuf;
  p.theDim = theDim;
  p.resBuf = resBuf;
  p.resDim = resDim;
  p.theDef = (void**)theDef;
  p.defDim = defDim;
  p.mat_aft = mat_aft;
  p.mat_bef = mat_bef;
  p.gain = 1.0;
  p.bias = 0.0;

  for ( i=0; i<chunks.n_allocated_chunks; i++ ) 
    chunks.data[i].parameters = (void*)(&p);
  
  /* processing
   */
  if ( processChunks( &_Reech2DNearestVectorField_DEFTYPE_TYPE, &chunks, proc ) != 1 ) {
      if ( _verbose_ )
	fprintf( stderr, "%s: unable to resample image\n", proc );
      freeChunks( &chunks );
      return;
  }

   freeChunks( &chunks );
}
