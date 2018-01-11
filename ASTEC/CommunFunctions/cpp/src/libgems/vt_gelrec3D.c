
#include <vt_gelrec3D.h>

#define COMPUTE3D_NEW(x,y,z) prod = (x) * v1[0] + (y) * v1[1] + (z) * v1[2];  \
	                     if ( prod > 0.0 )      f = prod * prod * irp1;   \
                             else if ( prod < 0.0 ) f = prod * prod * irn1;   \
	                     else                   f  = 0.0;                 \
	                     prod = (x) * v2[0] + (y) * v2[1] + (z) * v2[2];  \
	                     if ( prod > 0.0 )      f += prod * prod * irp2;  \
                             else if ( prod < 0.0 ) f += prod * prod * irn2;  \
	                     prod = (x) * v3[0] + (y) * v3[1] + (z) * v3[2];  \
	                     f += prod * prod * imin;                         \
	                     r = (double)255.0 - (double)128.0 * f;           \
	                     if ( r <= 0.0 ) new = (u8)0;          \
	                     else if ( r >= 255.0 ) new = (u8)255; \
	                     else new = (u8)( (int)(r+0.5) )

#ifndef NO_PROTO
int VT_Reconstruct3D( vt_image *ext, vt_image *scl, vt_line *par )
#else
int VT_Reconstruct3D( ext, scl, par )
vt_image *ext;
vt_image *scl;
vt_line *par;
#endif
{
    char *proc="VT_Reconstruct3D";
    double *sigma=(double*)NULL;
    register int is, i, j, k, d;
    double step, s;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf;

    /*--- tests ---*/
    if ( VT_Test2Image( ext, scl, proc ) == -1 ) return( -1 );
    if ( (ext->type != UCHAR) || (scl->type != UCHAR) ) {
	VT_Error( "bad type for images", proc );
	return( -1 );
    }
    
    /*--- calcul des sigmas ---*/
    sigma = (double*)VT_Malloc( (unsigned int)( (par->nb_coeff + 1) * sizeof( double ) ) );
    sigma[0] = 0.0;
    sigma[1] = (double)(par->first_coeff);
    if ( par->nb_coeff > 1 ) {
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
	for (i = 1; i < par->nb_coeff; i ++ )
	    sigma[i + 1] = (double)(par->first_coeff) + (double)(i) * step;
    }

    /*--- buffers ---*/
    ebuf = (u8 ***)(ext->array);
    sbuf = (u8 ***)(scl->array);
    dimx = ext->dim.x;
    dimy = ext->dim.y;
    dimz = ext->dim.z;
    
    /*--- ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	s = sigma[ (int)sbuf[z][y][x] ] *  sigma[ (int)sbuf[z][y][x] ];
	is = (int)(s + 0.5);
	for ( k = 0; k <= is ; k ++ )
	for ( j = 0; j <= k ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	    d = i*i + j*j + k*k;
	    if ( d > s ) continue;
	    /* ( +/- i , . , . ) */
	    if ( x+i < dimx ) {
		if ( y+j < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y+j][x+i] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y+j][x+i] = (u8)255;
		}
		if ( y-j >=0 ) {
		    if ( z+k < dimz ) ebuf[z+k][y-j][x+i] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y-j][x+i] = (u8)255;
		}
		if ( y+k < dimy ) {
		    if ( z+j < dimz ) ebuf[z+j][y+k][x+i] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y+k][x+i] = (u8)255;
		}
		if ( y-k >= 0 ) {
		    if ( z+j < dimz ) ebuf[z+j][y-k][x+i] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y-k][x+i] = (u8)255;
		}
	    }
	    if ( x-i >= 0 ) {
		if ( y+j < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y+j][x-i] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y+j][x-i] = (u8)255;
		}
		if ( y-j >=0 ) {
		    if ( z+k < dimz ) ebuf[z+k][y-j][x-i] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y-j][x-i] = (u8)255;
		}
		if ( y+k < dimy ) {
		    if ( z+j < dimz ) ebuf[z+j][y+k][x-i] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y+k][x-i] = (u8)255;
		}
		if ( y-k >= 0 ) {
		    if ( z+j < dimz ) ebuf[z+j][y-k][x-i] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y-k][x-i] = (u8)255;
		}
	    }
	    /* ( +/- j , . , . ) */
	    if ( x+j < dimx ) {
		if ( y+i < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y+i][x+j] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y+i][x+j] = (u8)255;
		}
		if ( y-i < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y-i][x+j] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y-i][x+j] = (u8)255;
		}
		if ( y+k < dimy ) {
		    if ( z+i < dimz ) ebuf[z+i][y+k][x+j] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y+k][x+j] = (u8)255;
		}
		if ( y-k >= 0 ) {
		    if ( z+i < dimz ) ebuf[z+i][y-k][x+j] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y-k][x+j] = (u8)255;
		}
	    }
	    if ( x-j >= 0 ) {
		if ( y+i < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y+i][x-j] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y+i][x-j] = (u8)255;
		}
		if ( y-i < dimy ) {
		    if ( z+k < dimz ) ebuf[z+k][y-i][x-j] = (u8)255;
		    if ( z-k >= 0   ) ebuf[z-k][y-i][x-j] = (u8)255;
		}
		if ( y+k < dimy ) {
		    if ( z+i < dimz ) ebuf[z+i][y+k][x-j] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y+k][x-j] = (u8)255;
		}
		if ( y-k >= 0 ) {
		    if ( z+i < dimz ) ebuf[z+i][y-k][x-j] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y-k][x-j] = (u8)255;
		}
	    }
	    /* ( +/- k , . , . ) */
	    if ( x+k < dimx ) {
		if ( y+j < dimy ) {
		    if ( z+i < dimz ) ebuf[z+i][y+j][x+k] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y+j][x+k] = (u8)255;
		}
		if ( y-j >=0 ) {
		    if ( z+i < dimz ) ebuf[z+i][y-j][x+k] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y-j][x+k] = (u8)255;
		}
		if ( y+i < dimy ) {
		    if ( z+j < dimz ) ebuf[z+j][y+i][x+k] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y+i][x+k] = (u8)255;
		}
		if ( y-i >= 0 ) {
		    if ( z+j < dimz ) ebuf[z+j][y-i][x+k] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y-i][x+k] = (u8)255;
		}
	    }
	    if ( x-k >= 0 ) {
		if ( y+j < dimy ) {
		    if ( z+i < dimz ) ebuf[z+i][y+j][x-k] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y+j][x-k] = (u8)255;
		}
		if ( y-j >=0 ) {
		    if ( z+i < dimz ) ebuf[z+i][y-j][x-k] = (u8)255;
		    if ( z-i >= 0   ) ebuf[z-i][y-j][x-k] = (u8)255;
		}
		if ( y+i < dimy ) {
		    if ( z+j < dimz ) ebuf[z+j][y+i][x-k] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y+i][x-k] = (u8)255;
		}
		if ( y-i >= 0 ) {
		    if ( z+j < dimz ) ebuf[z+j][y-i][x-k] = (u8)255;
		    if ( z-j >= 0   ) ebuf[z-j][y-i][x-k] = (u8)255;
		}
	    }
	}
    }

    VT_Free( (void**)&sigma );
    return( 1 );
}

#ifndef NO_PROTO
int VT_ScaleReconstruct3D( vt_resline *res, vt_line *par )
#else
int VT_ScaleReconstruct3D( res, par )
vt_resline *res;
vt_line *par;
#endif
{
    double *sigma=(double*)NULL;
    register int borne, i, j, k;
    double v1[3], v2[3], v3[3], f, r, step;
    register double xd, yd, zd, irp1, irn1, irp2, irn2, imin, prod;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf, new;
    r32 ***x1buf, ***y1buf, ***z1buf;
    r32 ***x2buf, ***y2buf, ***z2buf;
    /*--- tests ---*/

    /*--- buffers ---*/
    ebuf = (u8 ***)(res->imext.array);
    sbuf = (u8 ***)(res->imscale.array);
    x1buf = (r32 ***)(res->imdirx.array);
    y1buf = (r32 ***)(res->imdiry.array);
    z1buf = (r32 ***)(res->imdirz.array);
    x2buf = (r32 ***)(res->imdir2x.array);
    y2buf = (r32 ***)(res->imdir2y.array);
    z2buf = (r32 ***)(res->imdir2z.array);
    dimx = res->imext.dim.x;
    dimy = res->imext.dim.y;
    dimz = res->imext.dim.z;
    
    /*--- calcul des sigmas ---*/
    sigma = (double*)VT_Malloc( (unsigned int)( (par->nb_coeff + 1) * sizeof( double ) ) );
    sigma[0] = 0.0;
    sigma[1] = (double)(par->first_coeff);
    if ( par->nb_coeff > 1 ) {
	step = (double)(par->last_coeff - par->first_coeff) / (double)(par->nb_coeff - 1);
	for (i = 1; i < par->nb_coeff; i ++ )
	    sigma[i + 1] = (double)(par->first_coeff) + (double)(i) * step;
    }

    /*--- initialisation ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) 
      ebuf[z][y][x] = (u8)0;

    /*--- iso-surface a 127.0 ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	/*--- vecteurs ---*/
	v1[0] = x1buf[z][y][x];   v1[1] = y1buf[z][y][x];   v1[2] = z1buf[z][y][x];
	v2[0] = x2buf[z][y][x];   v2[1] = y2buf[z][y][x];   v2[2] = z2buf[z][y][x];
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
	/*--- bornes : f vaut le sigma ---*/
	f = sigma[ (int)sbuf[z][y][x] ];
	borne = (int)(f + 0.5);
	if ( borne == 0 ) continue;

	/*--- 1 / rayons au carre ---*/
	imin = (double)1.0 / ( f * f );
	irp1 = irn1 = irp2 = irn2 = imin;

	for ( k = 0; k <= borne ; k ++ )
	for ( j = 0; j <= k ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	  /*--- coordonnees en double ---*/
	  xd = (double)i;
	  yd = (double)j;
	  zd = (double)k;
	  
	  /* ( + i , . , . ) */
	  if ( x+i < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, yd, zd );
		if (new > ebuf[z+k][y+j][x+i]) ebuf[z+k][y+j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, yd, (-zd) );
		if (new > ebuf[z-k][y+j][x+i]) ebuf[z-k][y+j][x+i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, (-yd), zd );
		if (new > ebuf[z+k][y-j][x+i]) ebuf[z+k][y-j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x+i]) ebuf[z-k][y-j][x+i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, zd, yd );
		if (new > ebuf[z+j][y+k][x+i]) ebuf[z+j][y+k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, zd, (-yd) );
		if (new > ebuf[z-j][y+k][x+i]) ebuf[z-j][y+k][x+i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, (-zd), yd );
		if (new > ebuf[z+j][y-k][x+i]) ebuf[z+j][y-k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x+i]) ebuf[z-j][y-k][x+i] = new;
	      }
	    }
	  }
	  /* ( - i , . , . ) */
	  if ( x-i >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), yd, zd );
		if (new > ebuf[z+k][y+j][x-i]) ebuf[z+k][y+j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), yd, (-zd) );
		if (new > ebuf[z-k][y+j][x-i]) ebuf[z-k][y+j][x-i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), (-yd), zd );
		if (new > ebuf[z+k][y-j][x-i]) ebuf[z+k][y-j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x-i]) ebuf[z-k][y-j][x-i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), zd, yd );
		if (new > ebuf[z+j][y+k][x-i]) ebuf[z+j][y+k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), zd, (-yd) );
		if (new > ebuf[z-j][y+k][x-i]) ebuf[z-j][y+k][x-i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), (-zd), yd );
		if (new > ebuf[z+j][y-k][x-i]) ebuf[z+j][y-k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x-i]) ebuf[z-j][y-k][x-i] = new;
	      }
	    }
	  }
	  /* ( + j , . , . ) */
	  if ( x+j < dimx ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, xd, zd );
		if (new > ebuf[z+k][y+i][x+j]) ebuf[z+k][y+i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, xd, (-zd) );
		if (new > ebuf[z-k][y+i][x+j]) ebuf[z-k][y+i][x+j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, (-xd), zd );
		if (new > ebuf[z+k][y-i][x+j]) ebuf[z+k][y-i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x+j]) ebuf[z-k][y-i][x+j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, zd, xd );
		if (new > ebuf[z+i][y+k][x+j]) ebuf[z+i][y+k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, zd, (-xd) );
		if (new > ebuf[z-i][y+k][x+j]) ebuf[z-i][y+k][x+j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, (-zd), xd );
		if (new > ebuf[z+i][y-k][x+j]) ebuf[z+i][y-k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x+j]) ebuf[z-i][y-k][x+j] = new;
	      }
	    }
	  }
	  /* ( - j , . , . ) */
	  if ( x-j >= 0 ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), xd, zd );
		if (new > ebuf[z+k][y+i][x-j]) ebuf[z+k][y+i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), xd, (-zd) );
		if (new > ebuf[z-k][y+i][x-j]) ebuf[z-k][y+i][x-j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), (-xd), zd );
		if (new > ebuf[z+k][y-i][x-j]) ebuf[z+k][y-i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x-j]) ebuf[z-k][y-i][x-j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), zd, xd );
		if (new > ebuf[z+i][y+k][x-j]) ebuf[z+i][y+k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), zd, (-xd) );
		if (new > ebuf[z-i][y+k][x-j]) ebuf[z-i][y+k][x-j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), (-zd), xd );
		if (new > ebuf[z+i][y-k][x-j]) ebuf[z+i][y-k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x-j]) ebuf[z-i][y-k][x-j] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x+k < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, yd, xd );
		if (new > ebuf[z+i][y+j][x+k]) ebuf[z+i][y+j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, yd, (-xd) );
		if (new > ebuf[z-i][y+j][x+k]) ebuf[z-i][y+j][x+k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, (-yd), xd );
		if (new > ebuf[z+i][y-j][x+k]) ebuf[z+i][y-j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x+k]) ebuf[z-i][y-j][x+k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, xd, yd );
		if (new > ebuf[z+j][y+i][x+k]) ebuf[z+j][y+i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, xd, (-yd) );
		if (new > ebuf[z-j][y+i][x+k]) ebuf[z-j][y+i][x+k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, (-xd), yd );
		if (new > ebuf[z+j][y-i][x+k]) ebuf[z+j][y-i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x+k]) ebuf[z-j][y-i][x+k] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x-k >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), yd, xd );
		if (new > ebuf[z+i][y+j][x-k]) ebuf[z+i][y+j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), yd, (-xd) );
		if (new > ebuf[z-i][y+j][x-k]) ebuf[z-i][y+j][x-k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), (-yd), xd );
		if (new > ebuf[z+i][y-j][x-k]) ebuf[z+i][y-j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x-k]) ebuf[z-i][y-j][x-k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), xd, yd );
		if (new > ebuf[z+j][y+i][x-k]) ebuf[z+j][y+i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), xd, (-yd) );
		if (new > ebuf[z-j][y+i][x-k]) ebuf[z-j][y+i][x-k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), (-xd), yd );
		if (new > ebuf[z+j][y-i][x-k]) ebuf[z+j][y-i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x-k]) ebuf[z-j][y-i][x-k] = new;
	      }
	    }
	  }
	}
    }
    VT_Free( (void**)&sigma );
    return( 1 );
}

#ifndef NO_PROTO
int VT_CsteReconstruct3D( vt_resline *res, double radius )
#else
int VT_CsteReconstruct3D( res, radius )
vt_resline *res;
double radius;
#endif
{
    register int borne, i, j, k;
    double v1[3], v2[3], v3[3], f, r;
    register double xd, yd, zd, irp1, irn1, irp2, irn2, imin, prod;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf, new;
    r32 ***x1buf, ***y1buf, ***z1buf;
    r32 ***x2buf, ***y2buf, ***z2buf;
    /*--- tests ---*/

    /*--- buffers ---*/
    ebuf = (u8 ***)(res->imext.array);
    sbuf = (u8 ***)(res->imscale.array);
    x1buf = (r32 ***)(res->imdirx.array);
    y1buf = (r32 ***)(res->imdiry.array);
    z1buf = (r32 ***)(res->imdirz.array);
    x2buf = (r32 ***)(res->imdir2x.array);
    y2buf = (r32 ***)(res->imdir2y.array);
    z2buf = (r32 ***)(res->imdir2z.array);
    dimx = res->imext.dim.x;
    dimy = res->imext.dim.y;
    dimz = res->imext.dim.z;
    
    /*--- bornes : f vaut le rayon ---*/
    f = radius;
    if ( f <= 0.0 ) f = 1.0;
    borne = (int)(f + 0.5);
    if ( borne == 0 ) borne = 1;

    /*--- 1 / rayons au carre ---*/
    imin = (double)1.0 / ( f * f );
    irp1 = irn1 = irp2 = irn2 = imin;

    /*--- initialisation ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) 
      ebuf[z][y][x] = (u8)0;

    /*--- iso-surface a 127.0 ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	/*--- vecteurs ---*/
	v1[0] = x1buf[z][y][x];   v1[1] = y1buf[z][y][x];   v1[2] = z1buf[z][y][x];
	v2[0] = x2buf[z][y][x];   v2[1] = y2buf[z][y][x];   v2[2] = z2buf[z][y][x];
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

	for ( k = 0; k <= borne ; k ++ )
	for ( j = 0; j <= k ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	  /*--- coordonnees en double ---*/
	  xd = (double)i;
	  yd = (double)j;
	  zd = (double)k;
	  
	  /* ( + i , . , . ) */
	  if ( x+i < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, yd, zd );
		if (new > ebuf[z+k][y+j][x+i]) ebuf[z+k][y+j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, yd, (-zd) );
		if (new > ebuf[z-k][y+j][x+i]) ebuf[z-k][y+j][x+i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, (-yd), zd );
		if (new > ebuf[z+k][y-j][x+i]) ebuf[z+k][y-j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x+i]) ebuf[z-k][y-j][x+i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, zd, yd );
		if (new > ebuf[z+j][y+k][x+i]) ebuf[z+j][y+k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, zd, (-yd) );
		if (new > ebuf[z-j][y+k][x+i]) ebuf[z-j][y+k][x+i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, (-zd), yd );
		if (new > ebuf[z+j][y-k][x+i]) ebuf[z+j][y-k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x+i]) ebuf[z-j][y-k][x+i] = new;
	      }
	    }
	  }
	  /* ( - i , . , . ) */
	  if ( x-i >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), yd, zd );
		if (new > ebuf[z+k][y+j][x-i]) ebuf[z+k][y+j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), yd, (-zd) );
		if (new > ebuf[z-k][y+j][x-i]) ebuf[z-k][y+j][x-i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), (-yd), zd );
		if (new > ebuf[z+k][y-j][x-i]) ebuf[z+k][y-j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x-i]) ebuf[z-k][y-j][x-i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), zd, yd );
		if (new > ebuf[z+j][y+k][x-i]) ebuf[z+j][y+k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), zd, (-yd) );
		if (new > ebuf[z-j][y+k][x-i]) ebuf[z-j][y+k][x-i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), (-zd), yd );
		if (new > ebuf[z+j][y-k][x-i]) ebuf[z+j][y-k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x-i]) ebuf[z-j][y-k][x-i] = new;
	      }
	    }
	  }
	  /* ( + j , . , . ) */
	  if ( x+j < dimx ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, xd, zd );
		if (new > ebuf[z+k][y+i][x+j]) ebuf[z+k][y+i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, xd, (-zd) );
		if (new > ebuf[z-k][y+i][x+j]) ebuf[z-k][y+i][x+j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, (-xd), zd );
		if (new > ebuf[z+k][y-i][x+j]) ebuf[z+k][y-i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x+j]) ebuf[z-k][y-i][x+j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, zd, xd );
		if (new > ebuf[z+i][y+k][x+j]) ebuf[z+i][y+k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, zd, (-xd) );
		if (new > ebuf[z-i][y+k][x+j]) ebuf[z-i][y+k][x+j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, (-zd), xd );
		if (new > ebuf[z+i][y-k][x+j]) ebuf[z+i][y-k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x+j]) ebuf[z-i][y-k][x+j] = new;
	      }
	    }
	  }
	  /* ( - j , . , . ) */
	  if ( x-j >= 0 ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), xd, zd );
		if (new > ebuf[z+k][y+i][x-j]) ebuf[z+k][y+i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), xd, (-zd) );
		if (new > ebuf[z-k][y+i][x-j]) ebuf[z-k][y+i][x-j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), (-xd), zd );
		if (new > ebuf[z+k][y-i][x-j]) ebuf[z+k][y-i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x-j]) ebuf[z-k][y-i][x-j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), zd, xd );
		if (new > ebuf[z+i][y+k][x-j]) ebuf[z+i][y+k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), zd, (-xd) );
		if (new > ebuf[z-i][y+k][x-j]) ebuf[z-i][y+k][x-j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), (-zd), xd );
		if (new > ebuf[z+i][y-k][x-j]) ebuf[z+i][y-k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x-j]) ebuf[z-i][y-k][x-j] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x+k < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, yd, xd );
		if (new > ebuf[z+i][y+j][x+k]) ebuf[z+i][y+j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, yd, (-xd) );
		if (new > ebuf[z-i][y+j][x+k]) ebuf[z-i][y+j][x+k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, (-yd), xd );
		if (new > ebuf[z+i][y-j][x+k]) ebuf[z+i][y-j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x+k]) ebuf[z-i][y-j][x+k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, xd, yd );
		if (new > ebuf[z+j][y+i][x+k]) ebuf[z+j][y+i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, xd, (-yd) );
		if (new > ebuf[z-j][y+i][x+k]) ebuf[z-j][y+i][x+k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, (-xd), yd );
		if (new > ebuf[z+j][y-i][x+k]) ebuf[z+j][y-i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x+k]) ebuf[z-j][y-i][x+k] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x-k >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), yd, xd );
		if (new > ebuf[z+i][y+j][x-k]) ebuf[z+i][y+j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), yd, (-xd) );
		if (new > ebuf[z-i][y+j][x-k]) ebuf[z-i][y+j][x-k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), (-yd), xd );
		if (new > ebuf[z+i][y-j][x-k]) ebuf[z+i][y-j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x-k]) ebuf[z-i][y-j][x-k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), xd, yd );
		if (new > ebuf[z+j][y+i][x-k]) ebuf[z+j][y+i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), xd, (-yd) );
		if (new > ebuf[z-j][y+i][x-k]) ebuf[z-j][y+i][x-k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), (-xd), yd );
		if (new > ebuf[z+j][y-i][x-k]) ebuf[z+j][y-i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x-k]) ebuf[z-j][y-i][x-k] = new;
	      }
	    }
	  }
	}
    }
    return( 1 );
}

#ifndef NO_PROTO
int VT_GreyReconstruct3D( vt_resline *res )
#else
int VT_GreyReconstruct3D( res )
vt_resline *res;
#endif
{
    register int borne, i, j, k;
    double v1[3], v2[3], v3[3], f, r;
    double rp1, rn1, rp2, rn2;
    register double xd, yd, zd, irp1, irn1, irp2, irn2, imin, prod;
    register int x, y, z;
    int dimx, dimy, dimz;
    u8 ***ebuf, ***sbuf, new;
    r32 ***x1buf, ***y1buf, ***z1buf, ***p1buf, ***n1buf;
    r32 ***x2buf, ***y2buf, ***z2buf, ***p2buf, ***n2buf;
    /*--- tests ---*/

    /*--- buffers ---*/
    ebuf = (u8 ***)(res->imext.array);
    sbuf = (u8 ***)(res->imscale.array);
    x1buf = (r32 ***)(res->imdirx.array);
    y1buf = (r32 ***)(res->imdiry.array);
    z1buf = (r32 ***)(res->imdirz.array);
    p1buf = (r32 ***)(res->imradp.array);
    n1buf = (r32 ***)(res->imradn.array);
    x2buf = (r32 ***)(res->imdir2x.array);
    y2buf = (r32 ***)(res->imdir2y.array);
    z2buf = (r32 ***)(res->imdir2z.array);
    p2buf = (r32 ***)(res->imrad2p.array);
    n2buf = (r32 ***)(res->imrad2n.array);
    dimx = res->imext.dim.x;
    dimy = res->imext.dim.y;
    dimz = res->imext.dim.z;
    
    /*--- initialisation ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) 
      ebuf[z][y][x] = (u8)0;

    /*--- iso-surface a 127.0 ----*/
    for ( z = 0; z < dimz; z ++ )
    for ( y = 0; y < dimy; y ++ )
    for ( x = 0; x < dimx; x ++ ) {
	if ( sbuf[z][y][x] == (u8)0 ) continue;
	/*--- vecteurs ---*/
	v1[0] = x1buf[z][y][x];   v1[1] = y1buf[z][y][x];   v1[2] = z1buf[z][y][x];
	v2[0] = x2buf[z][y][x];   v2[1] = y2buf[z][y][x];   v2[2] = z2buf[z][y][x];
	v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
	v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
	v3[2] = v1[0] * v2[1] - v1[1] * v2[0];
	/*--- bornes : f vaut le max ---*/
	f = p1buf[z][y][x];
	if ( f < n1buf[z][y][x] ) f =  n1buf[z][y][x];
	if ( f < p2buf[z][y][x] ) f =  p2buf[z][y][x];
	if ( f < n2buf[z][y][x] ) f =  n2buf[z][y][x];
	borne = (int)(f + 0.5);
	if ( borne == 0 ) continue;

	/*--- 1 / rayons au carre : f vaut min positif ---*/
	rp1 = (double)p1buf[z][y][x];   rn1 = (double)n1buf[z][y][x];
	rp2 = (double)p2buf[z][y][x];   rn2 = (double)n2buf[z][y][x];
	if ( (rp1 > 0.0) && (f > rp1) ) f = rp1;
	if ( (rn1 > 0.0) && (f > rn1) ) f = rn1;
	if ( (rp2 > 0.0) && (f > rp2) ) f = rp2;
	if ( (rn2 > 0.0) && (f > rn2) ) f = rn2;
	if ( rp1 <= 0.0 ) rp1 = f;
	if ( rn1 <= 0.0 ) rn1 = f;
	if ( rp2 <= 0.0 ) rp2 = f;
	if ( rn2 <= 0.0 ) rn2 = f;

	irp1 = (double)1.0 / ( rp1 * rp1 );
	irn1 = (double)1.0 / ( rn1 * rn1 );
	irp2 = (double)1.0 / ( rp2 * rp2 );
	irn2 = (double)1.0 / ( rn2 * rn2 );
	imin = (double)1.0 / ( f * f );

	for ( k = 0; k <= borne ; k ++ )
	for ( j = 0; j <= k ; j ++ )
	for ( i = 0; i <= j ; i ++ ) {
	  /*--- coordonnees en double ---*/
	  xd = (double)i;
	  yd = (double)j;
	  zd = (double)k;
	  
	  /* ( + i , . , . ) */
	  if ( x+i < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, yd, zd );
		if (new > ebuf[z+k][y+j][x+i]) ebuf[z+k][y+j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, yd, (-zd) );
		if (new > ebuf[z-k][y+j][x+i]) ebuf[z-k][y+j][x+i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( xd, (-yd), zd );
		if (new > ebuf[z+k][y-j][x+i]) ebuf[z+k][y-j][x+i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( xd, (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x+i]) ebuf[z-k][y-j][x+i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, zd, yd );
		if (new > ebuf[z+j][y+k][x+i]) ebuf[z+j][y+k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, zd, (-yd) );
		if (new > ebuf[z-j][y+k][x+i]) ebuf[z-j][y+k][x+i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( xd, (-zd), yd );
		if (new > ebuf[z+j][y-k][x+i]) ebuf[z+j][y-k][x+i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( xd, (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x+i]) ebuf[z-j][y-k][x+i] = new;
	      }
	    }
	  }
	  /* ( - i , . , . ) */
	  if ( x-i >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), yd, zd );
		if (new > ebuf[z+k][y+j][x-i]) ebuf[z+k][y+j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), yd, (-zd) );
		if (new > ebuf[z-k][y+j][x-i]) ebuf[z-k][y+j][x-i] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-xd), (-yd), zd );
		if (new > ebuf[z+k][y-j][x-i]) ebuf[z+k][y-j][x-i] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-yd), (-zd) );
		if (new > ebuf[z-k][y-j][x-i]) ebuf[z-k][y-j][x-i] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), zd, yd );
		if (new > ebuf[z+j][y+k][x-i]) ebuf[z+j][y+k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), zd, (-yd) );
		if (new > ebuf[z-j][y+k][x-i]) ebuf[z-j][y+k][x-i] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-xd), (-zd), yd );
		if (new > ebuf[z+j][y-k][x-i]) ebuf[z+j][y-k][x-i] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-xd), (-zd), (-yd) );
		if (new > ebuf[z-j][y-k][x-i]) ebuf[z-j][y-k][x-i] = new;
	      }
	    }
	  }
	  /* ( + j , . , . ) */
	  if ( x+j < dimx ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, xd, zd );
		if (new > ebuf[z+k][y+i][x+j]) ebuf[z+k][y+i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, xd, (-zd) );
		if (new > ebuf[z-k][y+i][x+j]) ebuf[z-k][y+i][x+j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( yd, (-xd), zd );
		if (new > ebuf[z+k][y-i][x+j]) ebuf[z+k][y-i][x+j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( yd, (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x+j]) ebuf[z-k][y-i][x+j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, zd, xd );
		if (new > ebuf[z+i][y+k][x+j]) ebuf[z+i][y+k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, zd, (-xd) );
		if (new > ebuf[z-i][y+k][x+j]) ebuf[z-i][y+k][x+j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( yd, (-zd), xd );
		if (new > ebuf[z+i][y-k][x+j]) ebuf[z+i][y-k][x+j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( yd, (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x+j]) ebuf[z-i][y-k][x+j] = new;
	      }
	    }
	  }
	  /* ( - j , . , . ) */
	  if ( x-j >= 0 ) {
	    if ( y+i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), xd, zd );
		if (new > ebuf[z+k][y+i][x-j]) ebuf[z+k][y+i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), xd, (-zd) );
		if (new > ebuf[z-k][y+i][x-j]) ebuf[z-k][y+i][x-j] = new;
	      }
	    }
	    if ( y-i < dimy ) {
	      if ( z+k < dimz ) {
		COMPUTE3D_NEW( (-yd), (-xd), zd );
		if (new > ebuf[z+k][y-i][x-j]) ebuf[z+k][y-i][x-j] = new;
	      }
	      if ( z-k >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-xd), (-zd) );
		if (new > ebuf[z-k][y-i][x-j]) ebuf[z-k][y-i][x-j] = new;
	      }
	    }
	    if ( y+k < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), zd, xd );
		if (new > ebuf[z+i][y+k][x-j]) ebuf[z+i][y+k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), zd, (-xd) );
		if (new > ebuf[z-i][y+k][x-j]) ebuf[z-i][y+k][x-j] = new;
	      }
	    }
	    if ( y-k >= 0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-yd), (-zd), xd );
		if (new > ebuf[z+i][y-k][x-j]) ebuf[z+i][y-k][x-j] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-yd), (-zd), (-xd) );
		if (new > ebuf[z-i][y-k][x-j]) ebuf[z-i][y-k][x-j] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x+k < dimx ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, yd, xd );
		if (new > ebuf[z+i][y+j][x+k]) ebuf[z+i][y+j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, yd, (-xd) );
		if (new > ebuf[z-i][y+j][x+k]) ebuf[z-i][y+j][x+k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( zd, (-yd), xd );
		if (new > ebuf[z+i][y-j][x+k]) ebuf[z+i][y-j][x+k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( zd, (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x+k]) ebuf[z-i][y-j][x+k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, xd, yd );
		if (new > ebuf[z+j][y+i][x+k]) ebuf[z+j][y+i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, xd, (-yd) );
		if (new > ebuf[z-j][y+i][x+k]) ebuf[z-j][y+i][x+k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( zd, (-xd), yd );
		if (new > ebuf[z+j][y-i][x+k]) ebuf[z+j][y-i][x+k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( zd, (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x+k]) ebuf[z-j][y-i][x+k] = new;
	      }
	    }
	  }
	  /* ( + k , . , . ) */
	  if ( x-k >= 0 ) {
	    if ( y+j < dimy ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), yd, xd );
		if (new > ebuf[z+i][y+j][x-k]) ebuf[z+i][y+j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), yd, (-xd) );
		if (new > ebuf[z-i][y+j][x-k]) ebuf[z-i][y+j][x-k] = new;
	      }
	    }
	    if ( y-j >=0 ) {
	      if ( z+i < dimz ) {
		COMPUTE3D_NEW( (-zd), (-yd), xd );
		if (new > ebuf[z+i][y-j][x-k]) ebuf[z+i][y-j][x-k] = new;
	      }
	      if ( z-i >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-yd), (-xd) );
		if (new > ebuf[z-i][y-j][x-k]) ebuf[z-i][y-j][x-k] = new;
	      }
	    }
	    if ( y+i < dimy ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), xd, yd );
		if (new > ebuf[z+j][y+i][x-k]) ebuf[z+j][y+i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), xd, (-yd) );
		if (new > ebuf[z-j][y+i][x-k]) ebuf[z-j][y+i][x-k] = new;
	      }
	    }
	    if ( y-i >= 0 ) {
	      if ( z+j < dimz ) {
		COMPUTE3D_NEW( (-zd), (-xd), yd );
		if (new > ebuf[z+j][y-i][x-k]) ebuf[z+j][y-i][x-k] = new;
	      }
	      if ( z-j >= 0   ) {
		COMPUTE3D_NEW( (-zd), (-xd), (-yd) );
		if (new > ebuf[z-j][y-i][x-k]) ebuf[z-j][y-i][x-k] = new;
	      }
	    }
	  }
	}
    }
    return( 1 );
}
