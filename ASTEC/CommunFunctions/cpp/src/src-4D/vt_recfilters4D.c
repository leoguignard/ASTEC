

#include <vt_recfilters4D.h>

static int _VERBOSE_ = 0;





int VT_MaximaGradient4D( vt_image4D *theIm, 
			 vt_image4D *resIm, 
			 vt_contours *spacePar,
			 vt_contours *timePar )
{
  char *proc = "VT_MaximaGradient4D";
  int t, z, y, x;
  vt_image4D theGx, theGy, theGz, theGt, theNorme;
  vt_recfilters space_par, time_par;
  int local_derivative;

  vt_name4D nameNorme;

  int traitement;
  double EPSILON = 0.00001; /* 1/2^16 = .00001525878906250000 */
  register double ngx, ngy, ngz, ngt, norm, norme_originale;
  r32 ****bufGx, ****bufGy, ****bufGz, ****bufGt, ****bufNorme;
  register double xr, yr, zr, tr, dx, dy, dz, dt;
  register int xi, yi, zi, ti;



  /*--- cas 3D ---*/
  /*--------------*/

  if ( (theIm->dimt < 4) || (spacePar->dim != VT_4D) || (timePar->dim != VT_4D) ) {

    for ( t=0; t<theIm->dimt; t++ ) {
      if ( _VERBOSE_ )
	fprintf( stderr, " %s: processing image #%3d\r", proc, t );
      if ( VT_MaximaGradient( &(theIm->images[t]),
			     &(resIm->images[t]),
			     spacePar ) != 1 ) {
	fprintf( stderr, "%s: 3D processing failed on image #%d.\n", proc, t );
	return( -1 );
      }
    }

    return( 1 );
  }




  /*--- cas 4D ---*/
  /*--------------*/


  /*--- allocation des images ---*/

  if ( VT_AllocAndInitImage4D( &theGx, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate X gradient image", proc );
    return( -1 );
  }
  if ( VT_AllocAndInitImage4D( &theGy, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate Y gradient image", proc );
    VT_FreeImage4D( &theGx );
    return( -1 );
  }
  if ( VT_AllocAndInitImage4D( &theGz, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate Z gradient image", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    return( -1 );
  }
  if ( VT_AllocAndInitImage4D( &theGt, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate T gradient image", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    return( -1 );
  }






  /*--- calcul des derivees et de la norme ---*/

  VT_RecFilters( &space_par );
  switch ( spacePar->type_filter ) {
  case VT_RECFILTERS_DERICHE :
  case VT_RECGAUSSIAN_DERICHE :
    space_par.type_filter = spacePar->type_filter;
    break;
  default :
    space_par.type_filter = VT_RECFILTERS_DERICHE;
  }
  space_par.value_coefficient = spacePar->value_coefficient;
  space_par.length_continue = spacePar->length_continue;

  

  VT_RecFilters( &time_par );
  switch ( timePar->type_filter ) {
  case VT_RECFILTERS_DERICHE :
  case VT_RECGAUSSIAN_DERICHE :
    time_par.type_filter = timePar->type_filter;
    break;
  default :
    time_par.type_filter = VT_RECFILTERS_DERICHE;
  }
  time_par.value_coefficient = timePar->value_coefficient;
  time_par.length_continue = timePar->length_continue;


  local_derivative = VT_DERIVATIVE_1_CONTOURS;





  space_par.derivative.x = VT_DERIVATIVE_0;
  space_par.derivative.y = VT_DERIVATIVE_0;
  space_par.derivative.z = VT_NODERIVATIVE;

  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along X and Y     \n", proc );
  if ( VT_RecFilter3DOnImage4D( theIm, &theGz, &space_par ) != 1 ) {
    VT_Error( "unable to smooth along X and Y", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }




  time_par.derivative.x = local_derivative;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_NODERIVATIVE;
  space_par.derivative.y = VT_NODERIVATIVE;
  space_par.derivative.z = VT_DERIVATIVE_0;

  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along Z     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &theGz, &theGt, &space_par ) != 1 ) {
    VT_Error( "unable to smooth along Z", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &theGt, &theGt, &time_par ) != 1 ) {
    VT_Error( "unable to derive along T", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }


  
  time_par.derivative.x = VT_DERIVATIVE_0;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_NODERIVATIVE;
  space_par.derivative.y = VT_NODERIVATIVE;
  space_par.derivative.z = local_derivative;


  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along Z     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &theGz, &theGz, &space_par ) != 1 ) {
    VT_Error( "unable to derive along Z", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &theGz, &theGz, &time_par ) != 1 ) {
    VT_Error( "unable to smooth along T", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }



  time_par.derivative.x = VT_DERIVATIVE_0;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_NODERIVATIVE;
  space_par.derivative.y = VT_NODERIVATIVE;
  space_par.derivative.z = VT_DERIVATIVE_0;

  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along Z     \n", proc );
  if ( VT_RecFilter3DOnImage4D( theIm, &theGx, &space_par ) != 1 ) {
    VT_Error( "unable to smooth along Z", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &theGx, &theGx, &time_par ) != 1 ) {
    VT_Error( "unable to smooth along T", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }



  
  space_par.derivative.x = VT_DERIVATIVE_0;
  space_par.derivative.y = local_derivative;
  space_par.derivative.z = VT_NODERIVATIVE;

  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along Y     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &theGx, &theGy, &space_par ) != 1 ) {
    VT_Error( "unable to derive along Y", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }



  space_par.derivative.x = local_derivative;
  space_par.derivative.y = VT_DERIVATIVE_0;
  space_par.derivative.z = VT_NODERIVATIVE;

  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along X     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &theGx, &theGx, &space_par ) != 1 ) {
    VT_Error( "unable to derive along X", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }


  VT_Name4D( &nameNorme );
  (void)VT_CreateName4D( &nameNorme, "norme.maxima4D", ".inr", theIm->dimt );

  if ( VT_AllocAndInitImage4D( &theNorme, &nameNorme,
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate modulus image", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    return( -1 );
  }
  if ( VT_NormeGradient4DImage4DWithDerivatives( &theGx, &theGy,
						 &theGz, &theGt, 
						 &theNorme ) != 1 ) {
    VT_Error( "unable to compute modulus image", proc );
    VT_FreeImage4D( &theGx );
    VT_FreeImage4D( &theGy );
    VT_FreeImage4D( &theGz );
    VT_FreeImage4D( &theGt );
    VT_FreeImage4D( &theNorme );
    return( -1 );
  }
    
  if ( 0 )
  (void)VT_WriteImage4D( &theNorme );
  VT_FreeName4D( &nameNorme );




  /*--- Extraction des maxima ---*/
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: maxima extraction    \n", proc );

  bufGx = (r32 ****)theGx.array;
  bufGy = (r32 ****)theGy.array;
  bufGz = (r32 ****)theGz.array;
  bufGt = (r32 ****)theGt.array;
  bufNorme = (r32 ****)theNorme.array;


  for ( t = 0; t < theIm->dimt; t ++ ) {
    if ( _VERBOSE_ )
      fprintf( stderr, " %s: processing image #%3d\r", proc, t );


    for ( z = 0; z < theIm->dim.z; z ++ )
    for ( y = 0; y < theIm->dim.y; y ++ )
    for ( x = 0; x < theIm->dim.x; x ++ ) {

      norme_originale = bufNorme[t][z][y][x];
      if ( norme_originale <= EPSILON ) {
	bufGx[t][z][y][x] = 0.0;
	continue;
      }


      /* traitement : X=1 Y=2 Z=4 T=8 */
      traitement = 15;

      if ( (x == 0) || (x == (theIm->dim.x - 1)) ) traitement -= 1;
      if ( (y == 0) || (y == (theIm->dim.y - 1)) ) traitement -= 2;
      if ( (z == 0) || (z == (theIm->dim.z - 1)) ) traitement -= 4;
      if ( (t == 0) || (t == (theIm->dimt - 1)) ) traitement -= 8;
      
      switch ( traitement ) {
	
      case 0 : /* coin */
	bufGx[t][z][y][x] = 0.0;
	break;
	
      case 1 : /* arete selon X */
	bufGx[t][z][y][x] = 0.0;
	break;
	
      case 2 : /* arete selon Y */
	bufGx[t][z][y][x] = 0.0;
	break;
	
      case 4 : /* arete selon Z */
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 8 : /* arete selon T */
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 3 : /* plan XY */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGy[t][z][y][x] * bufGy[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngy = bufGy[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 5 : /* plan XZ */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGz[t][z][y][x] * bufGz[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngz = bufGz[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 9 : /* plan XT */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 6 : /* plan YZ */
	
	norm = sqrt( bufGy[t][z][y][x] * bufGy[t][z][y][x] +
		     bufGz[t][z][y][x] * bufGz[t][z][y][x] );
	ngy = bufGy[t][z][y][x] / norm;
	ngz = bufGz[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 10 : /* plan YT */
	
	norm = sqrt( bufGy[t][z][y][x] * bufGy[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngy = bufGy[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 12 : /* plan ZT */
	
	norm = sqrt( bufGz[t][z][y][x] * bufGz[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngz = bufGz[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	bufGx[t][z][y][x] = 0.0;
	break;
	
	
      case 7 : /* volume XYZ */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGy[t][z][y][x] * bufGy[t][z][y][x] +
		     bufGz[t][z][y][x] * bufGz[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngy = bufGy[t][z][y][x] / norm;
	ngz = bufGz[t][z][y][x] / norm;
	
	xr = (double)x + ngx;
	yr = (double)y + ngy;
	zr = (double)z + ngz;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = t;         dt = 0.0;

	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi+1] * (1.0-dt) * dz * dy * dx;
	
	if ( norme_originale <= norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	

	xr = (double)x - ngx;
	yr = (double)y - ngy;
	zr = (double)z - ngz;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi+1] * (1.0-dt) * dz * dy * dx;
	
	if ( norme_originale < norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}

	bufGx[t][z][y][x] = norme_originale;

	break;
	
	
      case 11 : /* volume XYT */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGy[t][z][y][x] * bufGy[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngy = bufGy[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	xr = (double)x + ngx;
	yr = (double)y + ngy;
	tr = (double)t + ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = z;         dz = 0.0;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi+1] * dt * (1.0-dz) * dy * dx;
	
	if ( norme_originale <= norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	

	xr = (double)x - ngx;
	yr = (double)y - ngy;
	tr = (double)t - ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi+1] * dt * (1.0-dz) * dy * dx;
	
	if ( norme_originale < norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}

	bufGx[t][z][y][x] = norme_originale;

	break;
	
	
      case 13 : /* volume XZT */
	
	norm = sqrt( bufGx[t][z][y][x] * bufGx[t][z][y][x] +
		     bufGz[t][z][y][x] * bufGz[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngx = bufGx[t][z][y][x] / norm;
	ngz = bufGz[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	xr = (double)x + ngx;
	zr = (double)z + ngz;
	tr = (double)t + ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = y;         dy = 0.0;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi+1] * dt * dz * (1.0-dy) * dx;
	
	if ( norme_originale <= norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	

	xr = (double)x - ngx;
	zr = (double)z - ngz;
	tr = (double)t - ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi+1] * dt * dz * (1.0-dy) * dx;
	
	if ( norme_originale < norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}

	bufGx[t][z][y][x] = norme_originale;

	break;
	
	
      case 14 : /* volume YZT */
	
	norm = sqrt( bufGy[t][z][y][x] * bufGy[t][z][y][x] +
		     bufGz[t][z][y][x] * bufGz[t][z][y][x] +
		     bufGt[t][z][y][x] * bufGt[t][z][y][x] );
	ngy = bufGy[t][z][y][x] / norm;
	ngz = bufGz[t][z][y][x] / norm;
	ngt = bufGt[t][z][y][x] / norm;
	
	yr = (double)y + ngy;
	zr = (double)z + ngz;
	tr = (double)t + ngt;
	if ( (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = x;         dx = 0.0;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi+1][xi] * dt * dz * dy * (1.0-dx);
	
	if ( norme_originale <= norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	

	yr = (double)y - ngy;
	zr = (double)z - ngz;
	tr = (double)t - ngt;
	if ( (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi+1][xi] * dt * dz * dy * (1.0-dx);
	
	if ( norme_originale < norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}

	bufGx[t][z][y][x] = norme_originale;

	break;


	
	
      case 15 : /* volume 4D */
      default :
	ngx = bufGx[t][z][y][x] / norme_originale;
	ngy = bufGy[t][z][y][x] / norme_originale;
	ngz = bufGz[t][z][y][x] / norme_originale;
	ngt = bufGt[t][z][y][x] / norme_originale;
	
	xr = (double)x + ngx;
	yr = (double)y + ngy;
	zr = (double)z + ngz;
	tr = (double)t + ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi+1] * (1.0-dt) * dz * dy * dx
	  + bufNorme[ti+1][zi][yi+1][xi+1] * dt * (1.0-dz) * dy * dx
	  + bufNorme[ti+1][zi+1][yi][xi+1] * dt * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi+1][yi+1][xi] * dt * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi+1][xi+1] * dt * dz * dy * dx; 
	
	if ( norme_originale <= norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	

	xr = (double)x - ngx;
	yr = (double)y - ngy;
	zr = (double)z - ngz;
	tr = (double)t - ngt;
	if ( (xr < 0.0) || ( xr >= theIm->dim.x - 1) ||
	     (yr < 0.0) || ( yr >= theIm->dim.y - 1) ||
	     (zr < 0.0) || ( zr >= theIm->dim.z - 1) ||
	     (tr < 0.0) || ( tr >= theIm->dimt - 1) ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}
	
	xi = (int)xr;   dx = xr - (double)xi;
	yi = (int)yr;   dy = yr - (double)yi;
	zi = (int)zr;   dz = zr - (double)zi;
	ti = (int)tr;   dt = tr - (double)ti;
	
	norm = bufNorme[ti][zi][yi][xi] * (1.0-dt) * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi][xi+1] * (1.0-dt) * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi][yi+1][xi] * (1.0-dt) * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi][xi] * (1.0-dt) * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi][xi] * dt * (1.0-dz) * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi][yi+1][xi+1] * (1.0-dt) * (1.0-dz) * dy * dx
	  + bufNorme[ti][zi+1][yi][xi+1] * (1.0-dt) * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi][yi][xi+1] * dt * (1.0-dz) * (1.0-dy) * dx
	  + bufNorme[ti][zi+1][yi+1][xi] * (1.0-dt) * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi][yi+1][xi] * dt * (1.0-dz) * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi][xi] * dt * dz * (1.0-dy) * (1.0-dx) 
	  + bufNorme[ti][zi+1][yi+1][xi+1] * (1.0-dt) * dz * dy * dx
	  + bufNorme[ti+1][zi][yi+1][xi+1] * dt * (1.0-dz) * dy * dx
	  + bufNorme[ti+1][zi+1][yi][xi+1] * dt * dz * (1.0-dy) * dx
	  + bufNorme[ti+1][zi+1][yi+1][xi] * dt * dz * dy * (1.0-dx) 
	  + bufNorme[ti+1][zi+1][yi+1][xi+1] * dt * dz * dy * dx; 
	
	if ( norme_originale < norm ) {
	  bufGx[t][z][y][x] = 0.0;
	  break;
	}

	bufGx[t][z][y][x] = norme_originale;


      } /* fin du switch */
    } /* fin de la boucle sur x */

  
  } /* fin de la boucle sur t */


  VT_FreeImage4D( &theGy );
  VT_FreeImage4D( &theGz );
  VT_FreeImage4D( &theGt );
  VT_FreeImage4D( &theNorme );



  if ( VT_CopyImage4D( &theGx, resIm ) != 1 ) {
    VT_Error( "unable to copy into output image", proc );
    VT_FreeImage4D( &theGx );
    return( -1 );
  }




  VT_FreeImage4D( &theGx );

  return( 1 );
}











/* calcul d'un filtrage recursif
   selon T sur une image 4D 
*/
int VT_RecFilterTOnImage4D( vt_image4D *theIm, 
			    vt_image4D *resIm, 
			    vt_recfilters *par )
{
  char *proc = "VT_RecFilterTOnImage4D";
  int t, u, x, y, z;
  int dim, offset, local_dim;
  double *in, *out, *work;
  r32 ****resBuf = (r32****)NULL;
  

  
  if ( theIm->dimt < 4 ) {
    VT_Error( "too small T dimension", proc );
    return( 0 );
  }


  if ( theIm->dimt != resIm->dimt ) {
    VT_Error( "input and output have different T dimensions", proc );
    return( 0 );
  }
  
  for ( t=1; t<theIm->dimt; t++ ) {
    if ( theIm->images[0].type != theIm->images[t].type ) {
      VT_Error( "input images have different types", proc );
      return( 0 );
    }
    if ( resIm->images[0].type != resIm->images[t].type ) {
      VT_Error( "output images have different types", proc );
      return( 0 );
    }
  }

  if ( resIm->images[0].type != FLOAT ) {
    VT_Error( "unable to deal with such output image type", proc );
    return( 0 );
  }
  resBuf = (r32****)(resIm->array);




  dim = theIm->dimt;
  offset = 0;
  if ( par->length_continue.x > 0 ) offset = par->length_continue.x;
  local_dim = dim + 2 * offset;
  in = out = work = (double*)NULL;
  in   = (double*)VT_Malloc( (unsigned int)(local_dim * sizeof( double ) ) );
  out  = (double*)VT_Malloc( (unsigned int)(local_dim * sizeof( double ) ) );
  work = (double*)VT_Malloc( (unsigned int)(local_dim * sizeof( double ) ) );
  if ( (in == (double*)NULL) || (out == (double*)NULL) || (work == (double*)NULL) ) {
    VT_Error( "unable to allocate auxiliary arrays", proc );
    VT_Free( (void**)&in );
    VT_Free( (void**)&out );
    VT_Free( (void**)&work );
    return( -1 );
  }

  if ( par->derivative.x == VT_NODERIVATIVE ) {
    VT_Error( "nothing to compute", proc );
    VT_Free( (void**)&in );
    VT_Free( (void**)&out );
    VT_Free( (void**)&work );
    return( -1 );
  }

  if ( VT_InitRecursiveCoefficients( (double)(par->value_coefficient.x), 
				     par->type_filter, par->derivative.x ) != 1 ) {
    VT_Error( "unable to initialize coefficients to filter along T", proc );
    VT_Free( (void**)&in );
    VT_Free( (void**)&out );
    VT_Free( (void**)&work );
    return( -1 );
  }
    

  switch( theIm->images[0].type ) {
  case UCHAR :
    {
      u8 ****array = (u8****)(theIm->array);

      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	/* saisie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  in[u] = (double)(array[t][z][y][x]);
	/* ajout de points */
	t = theIm->dimt + offset;
	for ( u = 0; u < offset; u++, t++ ) {
	  in[u] = in[offset];
	  in[t] = in[theIm->dimt + offset - 1];
	}
	/*--- traitement de la ligne ---*/
	if ( VT_RecFilterOnLine( in, out, work, out, local_dim ) != 1 ) {
	  VT_Error( "unable to filter along T", proc );
	  VT_Free( (void**)&in );
	  VT_Free( (void**)&out );
	  VT_Free( (void**)&work );
	  return( -1 );
	}
	/* copie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  resBuf[t][z][y][x] = (r32)(out[u]);
      }
      
    }
    break;
    

  case USHORT :
    {
      u16 ****array = (u16****)(theIm->array);

      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	/* saisie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  in[u] = (double)(array[t][z][y][x]);
	/* ajout de points */
	t = theIm->dimt + offset;
	for ( u = 0; u < offset; u++, t++ ) {
	  in[u] = in[offset];
	  in[t] = in[theIm->dimt + offset - 1];
	}
	/*--- traitement de la ligne ---*/
	if ( VT_RecFilterOnLine( in, out, work, out, local_dim ) != 1 ) {
	  VT_Error( "unable to filter along T", proc );
	  VT_Free( (void**)&in );
	  VT_Free( (void**)&out );
	  VT_Free( (void**)&work );
	  return( -1 );
	}
	/* copie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  resBuf[t][z][y][x] = (r32)(out[u]);
      }
      
    }
    break;
    

  case FLOAT :
    {
      r32 ****array = (r32****)(theIm->array);

      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	/* saisie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  in[u] = (double)(array[t][z][y][x]);
	/* ajout de points */
	t = theIm->dimt + offset;
	for ( u = 0; u < offset; u++, t++ ) {
	  in[u] = in[offset];
	  in[t] = in[theIm->dimt + offset - 1];
	}
	/*--- traitement de la ligne ---*/
	if ( VT_RecFilterOnLine( in, out, work, out, local_dim ) != 1 ) {
	  VT_Error( "unable to filter along T", proc );
	  VT_Free( (void**)&in );
	  VT_Free( (void**)&out );
	  VT_Free( (void**)&work );
	  return( -1 );
	}
	/* copie de la ligne */
	for ( t=0, u=offset; t<theIm->dimt; t++, u++ )
	  resBuf[t][z][y][x] = (r32)(out[u]);
      }
      
    }
    break;


  default :
    VT_Error( "unable to deal with such input image type", proc );
    return( 0 );
    
  }
  

  VT_Free( (void**)&in );
  VT_Free( (void**)&out );
  VT_Free( (void**)&work );
  return( 1 );

}







/* calcul d'un filtrage recursif
   selon X, Y et Z sur toutes les
   images 3D d'une image 4D 
*/
int VT_RecFilter3DOnImage4D( vt_image4D *theIm, 
			     vt_image4D *resIm, 
			     vt_recfilters *par )
{
  char *proc = "VT_RecFilter3DOnImage4D";
  int t;

  for( t=0; t<theIm->dimt; t++ ) {

    if ( _VERBOSE_ )
      fprintf( stderr, " %s: processing image #%3d\r", proc, t );

    if ( VT_RecFilterOnImage( &(theIm->images[t]),
			      &(resIm->images[t]),
			      par ) != 1 ) {
      fprintf( stderr, "%s: processing failed on image #%d.\n", proc, t );
    }
  }
  return( 1 );
}






int VT_NormeGradient4DImage4D( vt_image4D *theIm,
			       vt_image4D *resIm,
			       vt_contours *spacePar,
			       vt_contours *timePar,
			       int derivative )
{
  char *proc="VT_NormeGradient4DImage4D";
  vt_image4D tmpIm, sumIm;
  vt_recfilters space_par, time_par;
  int local_derivative;
  int t, z, y, x;
  r32 ****tmp, ****sum;




  if ( (theIm->dimt < 4) || (spacePar->dim != VT_4D) || (timePar->dim != VT_4D) ) {

    for ( t=0; t<theIm->dimt; t++ ) {
      if ( _VERBOSE_ )
	fprintf( stderr, " %s: processing image #%3d\r", proc, t );
      if ( VT_NormeGradient( &(theIm->images[t]),
			     &(resIm->images[t]),
			     spacePar,
			     derivative ) != 1 ) {
	fprintf( stderr, "%s: 3D processing failed on image #%d.\n", proc, t );
	return( -1 );
      }
    }

    return( 1 );
  }








  if ( VT_AllocAndInitImage4D( &tmpIm, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate tmp image", proc );
    return( -1 );
  }
  if ( VT_AllocAndInitImage4D( &sumIm, NULL, 
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, FLOAT ) != 1 ) {
    VT_Error( "unable to allocate sum image", proc );
    VT_FreeImage4D( &tmpIm );
    return( -1 );
  }



  tmp = (r32****)tmpIm.array;
  sum = (r32****)sumIm.array;



  
  VT_RecFilters( &space_par );
  switch ( spacePar->type_filter ) {
  case VT_RECFILTERS_DERICHE :
  case VT_RECGAUSSIAN_DERICHE :
    space_par.type_filter = spacePar->type_filter;
    break;
  default :
    space_par.type_filter = VT_RECFILTERS_DERICHE;
  }
  space_par.value_coefficient = spacePar->value_coefficient;
  space_par.length_continue = spacePar->length_continue;

  

  VT_RecFilters( &time_par );
  switch ( timePar->type_filter ) {
  case VT_RECFILTERS_DERICHE :
  case VT_RECGAUSSIAN_DERICHE :
    time_par.type_filter = timePar->type_filter;
    break;
  default :
    time_par.type_filter = VT_RECFILTERS_DERICHE;
  }
  time_par.value_coefficient = timePar->value_coefficient;
  time_par.length_continue = timePar->length_continue;


  switch ( derivative ) {
  case VT_DERIVATIVE_1_CONTOURS :
    local_derivative = VT_DERIVATIVE_1_CONTOURS;
    break;
  default :
    local_derivative = VT_DERIVATIVE_1;
  }




  space_par.derivative.x = VT_DERIVATIVE_0;
  space_par.derivative.y = VT_DERIVATIVE_0;
  space_par.derivative.z = VT_NODERIVATIVE;
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along X and Y     \n", proc );
  if ( VT_RecFilter3DOnImage4D( theIm, &tmpIm, &space_par ) != 1 ) {
    VT_Error( "unable to smooth along X and Y", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }








  time_par.derivative.x = local_derivative;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_NODERIVATIVE;
  space_par.derivative.y = VT_NODERIVATIVE;
  space_par.derivative.z = VT_DERIVATIVE_0;


  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along Z     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &tmpIm, &sumIm, &space_par ) != 1 ) {
    VT_Error( "unable to smooth along Z", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &sumIm, &sumIm, &time_par ) != 1 ) {
    VT_Error( "unable to derive along T", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }









  time_par.derivative.x = VT_DERIVATIVE_0;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_NODERIVATIVE;
  space_par.derivative.y = VT_NODERIVATIVE;
  space_par.derivative.z = local_derivative;


  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along Z     \n", proc );
  if ( VT_RecFilter3DOnImage4D( &tmpIm, &tmpIm, &space_par ) != 1 ) {
    VT_Error( "unable to derive along Z", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &tmpIm, &tmpIm, &time_par ) != 1 ) {
    VT_Error( "unable to smooth along T", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }

  for ( t=0; t<theIm->dimt; t++ )
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ )
  for ( x=0; x<theIm->dim.x; x++ ) {
    sum[t][z][y][x] = sum[t][z][y][x] * sum[t][z][y][x] + 
      tmp[t][z][y][x] * tmp[t][z][y][x];
  }








  time_par.derivative.x = VT_DERIVATIVE_0;
  time_par.derivative.y = VT_NODERIVATIVE;
  time_par.derivative.z = VT_NODERIVATIVE;

  space_par.derivative.x = VT_DERIVATIVE_0;
  space_par.derivative.y = local_derivative;
  space_par.derivative.z = VT_DERIVATIVE_0;


  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along Y     \n", proc );
  if ( VT_RecFilter3DOnImage4D( theIm, &tmpIm, &space_par ) != 1 ) {
    VT_Error( "unable to derive along Y", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &tmpIm, &tmpIm, &time_par ) != 1 ) {
    VT_Error( "unable to smooth along T", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }

  for ( t=0; t<theIm->dimt; t++ )
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ )
  for ( x=0; x<theIm->dim.x; x++ ) {
    sum[t][z][y][x] += tmp[t][z][y][x] * tmp[t][z][y][x] ;
  }





  space_par.derivative.x = local_derivative;
  space_par.derivative.y = VT_DERIVATIVE_0;
  space_par.derivative.z = VT_DERIVATIVE_0;


  if ( _VERBOSE_ )
    fprintf( stderr, " %s: derivating along X     \n", proc );
  if ( VT_RecFilter3DOnImage4D( theIm, &tmpIm, &space_par ) != 1 ) {
    VT_Error( "unable to derive along X", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }
  if ( _VERBOSE_ )
    fprintf( stderr, " %s: smoothing along T     \n", proc );
  if ( VT_RecFilterTOnImage4D( &tmpIm, &tmpIm, &time_par ) != 1 ) {
    VT_Error( "unable to smooth along T", proc );
    VT_FreeImage4D( &tmpIm );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }

  for ( t=0; t<theIm->dimt; t++ )
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ )
  for ( x=0; x<theIm->dim.x; x++ ) {
    sum[t][z][y][x] = (r32)sqrt( (double)(sum[t][z][y][x] 
					  + tmp[t][z][y][x] * tmp[t][z][y][x]) );
  }



  VT_FreeImage4D( &tmpIm );



  if ( VT_CopyImage4D( &sumIm, resIm ) != 1 ) {
    VT_Error( "unable to copy into output image", proc );
    VT_FreeImage4D( &sumIm );
    return( -1 );
  }




  VT_FreeImage4D( &sumIm );
  return( 1 );
  
}






int VT_NormeGradient4DImage4DWithDerivatives( vt_image4D *theX, 
			       vt_image4D *theY,
			       vt_image4D *theZ,
			       vt_image4D *theT,
			       vt_image4D *theNorme )
{
  char *proc = "VT_NormeGradient4DImage4D";
  int t,z,y,x;

  if ( (theX->dimt != theY->dimt) || (theX->dimt != theZ->dimt) ||
       (theX->dimt != theT->dimt) || (theX->dimt != theNorme->dimt) ) {
    VT_Error( "Different t dimensions in images", proc );
    return( 0 );
  }



  for ( t=0; t < theX->dimt; t ++ ) {

    if ( VT_Test2Image( &(theX->images[t]), &(theY->images[t]), proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(theX->images[t]), &(theZ->images[t]), proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(theX->images[t]), &(theT->images[t]), proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(theX->images[t]), &(theNorme->images[t]), proc ) == -1 ) return( -1 );
    if ( ((theX->images[t]).type != (theY->images[t]).type) ||
	 ((theX->images[t]).type != (theZ->images[t]).type) ||
	 ((theX->images[t]).type != (theT->images[t]).type) ||
	 ((theX->images[t]).type != (theNorme->images[t]).type) ) {
      VT_Error( "Different types in images", proc );
      return( 0 );
    }



    switch ( (theX->images[t]).type ) {
    case FLOAT :
      {
	r32 ***bx = (r32***)(theX->images[t].array);
	r32 ***by = (r32***)(theY->images[t].array);
	r32 ***bz = (r32***)(theZ->images[t].array);
	r32 ***bt = (r32***)(theT->images[t].array);
	r32 ***bn = (r32***)(theNorme->images[t].array);

	for ( z=0; z < theX->images[t].dim.z; z++ )
        for ( y=0; y < theX->images[t].dim.y; y++ )
	for ( x=0; x < theX->images[t].dim.x; x++ ) {
	  bn[z][y][x] = (r32)sqrt( (double)bx[z][y][x] * (double)bx[z][y][x] +
				   (double)by[z][y][x] * (double)by[z][y][x] +
				   (double)bz[z][y][x] * (double)bz[z][y][x] +
				   (double)bt[z][y][x] * (double)bt[z][y][x] );
	}
      }
      break;
      
    default :
      VT_Error("such type not handled yet", proc );
      return( 0 );
    }

  }
 
  return( 1 );
}







int VT_NormeGradient3DImage4DWithDerivatives( vt_image4D *theX, 
			       vt_image4D *theY,
			       vt_image4D *theZ,
			       vt_image4D *theNorme )
{
  char *proc = "VT_NormeGradient3DImage4D";
  int t,z,y,x;

  if ( (theX->dimt != theY->dimt) || (theX->dimt != theZ->dimt) ||
       (theX->dimt != theNorme->dimt) ) {
    VT_Error( "Different t dimensions in images", proc );
    return( 0 );
  }



  for ( t=0; t < theX->dimt; t ++ ) {

    if ( VT_Test2Image( &(theX->images[t]), &(theY->images[t]), proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(theX->images[t]), &(theZ->images[t]), proc ) == -1 ) return( -1 );
    if ( VT_Test2Image( &(theX->images[t]), &(theNorme->images[t]), proc ) == -1 ) return( -1 );
    if ( ((theX->images[t]).type != (theY->images[t]).type) ||
	 ((theX->images[t]).type != (theZ->images[t]).type) ||
	 ((theX->images[t]).type != (theNorme->images[t]).type) ) {
      VT_Error( "Different types in images", proc );
      return( 0 );
    }



    switch ( (theX->images[t]).type ) {
    case FLOAT :
      {
	r32 ***bx = (r32***)(theX->images[t].array);
	r32 ***by = (r32***)(theY->images[t].array);
	r32 ***bz = (r32***)(theZ->images[t].array);
	r32 ***bn = (r32***)(theNorme->images[t].array);

	for ( z=0; z < theX->images[t].dim.z; z++ )
        for ( y=0; y < theX->images[t].dim.y; y++ )
	for ( x=0; x < theX->images[t].dim.x; x++ ) {
	  bn[z][y][x] = (r32)sqrt( (double)bx[z][y][x] * (double)bx[z][y][x] +
				   (double)by[z][y][x] * (double)by[z][y][x] +
				   (double)bz[z][y][x] * (double)bz[z][y][x] );
	}
      }
      break;
      
    default :
      VT_Error("such type not handled yet", proc );
      return( 0 );
    }

  }
 
  return( 1 );
}












void VT_Recfilters4DVerbose( )
{
  _VERBOSE_ = 1;
}

void VT_Recfilters4DNoVerbose( )
{
  _VERBOSE_ = 0;
}
