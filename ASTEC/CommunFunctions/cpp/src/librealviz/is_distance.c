/*************************************************************************
 * is_distance.c - calcul de distance signee
 *
 * $Id: is_distance.c,v 1.12 2000/08/11 16:34:57 greg Exp $
 *
 * DESCRIPTION: 
 *
 * Outils de calculs de distance
 *
 * Tous travaille sur des buffers images, supposes contigus
 * en memoire.
 *
 * pour calculer une carte de distance, on l'initialise d'abord
 * avec des grandes valeurs positives (_PLUS_INFINITY_) a l'interieur 
 * de l'objet et de grandes valeurs negatives (_MINUS_INFINITY_) 
 * a l'exterieur de l'objet.
 *
 * fonction -> _InitSignedDistanceMap()
 *
 * En fait, le signe importe peu, puisqu'on se servira de la carte
 * pour calculer un point le plus proche, ce qui est important
 * c'est qu'il est un changement de signe a l'interface et
 * que les distances y soient correctement definies.
 *
 * 
 * Ensuite, le calcul est entierement fait dans la fonction
 * _ComputeSignedDistanceMap()
 *
 * 
 * Les details du calcul sont les suivants:
 *
 *
 * On calcule ensuite les coefficients pour la carte de 
 * distance, qui dependent de la taille (en millimetres)
 * du voxel. 
 * Il y a des coefficients (masque 7x7x7) pour l'initialisation
 * des distances a l'interface, et des coefficients (masque 5x5x5)
 * pour la propagation en chamfrein 3x3x3 et 5x5x5. 
 * De plus, Pour garder une certain precision, les valeurs calculees
 * sont artificiellement multipliees par 100/min( taille du voxel ).
 * 
 * Il faudra donc multiplier les valeurs de la carte de distance obtenue
 * par min( taille du voxel )/100.0 (contenue dans la structure de 
 * parametres) pour avoir la distance en MILLIMETRES.
 *
 * fonction -> _InitDistanceCoefficients()
 *
 *
 *
 * Ensuite, il faut initialiser les distances a l'interface.
 *
 * fonction -> _InitBorderSignedDistanceMap()
 *
 * Cela est fait avec un masque 7x7x7, et des distances exactes
 * (modulo les arrondis de calcul) afin de garder une bonne 
 * precision de calcul du point le plus proche lorsque l'on
 * est effectivement proche de l'interface.
 * En outre, ce sont ces distances qui serviront de conditions
 * initiales pour la propagation
 *
 *
 *
 * On passe ensuite au calcul proprement dit, d'abord une
 * propagation avec un chamfrein 3x3x3
 *
 * fonction -> _ComputeSignedDistanceMapWithChamfer3x3x3()
 *
 *
 *
 * Cette premiere distance calculee, on peut affiner 
 * l'estimation avec une distance du chamfrein 5x5x5,
 * pour les distances inferieures a un certain seuil.
 *
 * fonction -> _UpdateSignedDistanceMapWithChamfer5x5x5()
 *
 *
 *
 * Il existe aussi une fonction d'interpolation de la 
 * distance 
 *
 * fonction -> _InterpoleSignedDistance()
 *
 * qui n'estime pas la distance a l'exterieur de l'image,
 * ce qu'il faudrait ajouter.
 *
 *
 *
 * AUTHOR:
 * Gregoire Malandain
 *
 * CREATION DATE:
 * Sat May 22 11:05:27 MET DST 1999
 *
 * Copyright Gregoire Malandain
 *
 *
 * ADDITIONS, CHANGES:
 *
 * - Mon Aug  7 10:35:21 MET DST 2000 (Gregoire Malandain)
 *   Ajout de _Update2DSignedDistanceMapWithChamfer5x5x5()
 *   et traitement des images binaires
 *
 * - Tue Jul  6 17:08:21 MET DST 1999 (Gregoire Malandain)
 *   Bug dans _InitBorderSignedDistanceMap(), 
 *   aux bords de l'image, mauvaise gestion de l'appartenance des
 *   voisins du point [tt] a l'image.
 *   
 *
 */



#include <is_distance.h>


static int _VERBOSE_ = 1;









void _DistanceSetVerbose ()
{
  _VERBOSE_ = 1;
}
void _DistanceSetNoVerbose ()
{
  _VERBOSE_ = 0;
}






typedef enum {
  _TWO_DIMENSIONAL_,
  _THREE_DIMENSIONAL_
} enumComputation;

static enumComputation typeComputation = _THREE_DIMENSIONAL_;

void _DistanceSetComputationTo2D ()
{
  typeComputation = _TWO_DIMENSIONAL_;
}

void _DistanceSetComputationTo3D ()
{
  typeComputation = _THREE_DIMENSIONAL_;
}










static void _InitBorder3DSignedDistanceMap( typeDistanceMap *theDist,
				   const typeDistanceCoefficients *par );
static void _InitBorder2DSignedDistanceMap( typeDistanceMap *theDist,
				   const typeDistanceCoefficients *par );

static void _Compute3DSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par );
static void _Compute2DSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par );


static void _Update3DSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate );
static void _Update2DSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate );








/* Fonction de calcul de la distance en un point reel

   Ce point est en coordonnees voxels, non millimetriques.
   
   En fait, on ne calcule pas la distance mais on
   evalue la valeur de la carte de distance
   par interpolation. 
   Cette valeur doit etre multipliee par une constante
   afin d'avoir la distance en MILLIMETRES.
   
   La valeur n'est pas estimee en dehors de l'image,
   ce qu'il faudrait faire de la facon suivante :
   pour les points en dehors de l'image,
   il faut estimer le plus proche point sur 
   le bord de l'image, estimer la valeur de la distance
   en ce plus proche point et y ajouter la distance
   pour l'atteindre. 
   Cela necessite en plus de passer la taille du voxel.
*/



int _InterpoleSignedDistance( double *d, 
			      const double x, const double y, const double z, 
			      const typeDistanceMap *theDist )
{
  int ix, iy, iz, t;
  double newx, newy, newz;
  double dx, dy, dz, dxdy,dxdz,dydz,dxdydz;
  double v6, v5, v4;
  double res = 0.0;

  /* on n'a pas traite le cas ou on sort de
     l'image
  */

  ix = (int)x;
  iy = (int)y;
  iz = (int)z;

  if ( (x >= 0.0) && (ix < theDist->dim[0]-1) &&
       (y >= 0.0) && (iy < theDist->dim[1]-1) &&
       (z >= 0.0) && (iz < theDist->dim[2]-1) ) {
    
    t = iz*theDist->dim[0]*theDist->dim[1] + iy*theDist->dim[0] + ix;

    dx = x - ix;
    dy = y - iy;
    dz = z - iz;
    dxdy = dx*dy;
    dxdz = dx*dz;
    dydz = dy*dz;
    dxdydz = dxdy*dz;
    
    v6 = dxdz-dxdydz;
    v5 = dxdy-dxdydz;
    v4 = dx-dxdy-v6;
    
    
    res += (1-dy-dz+dydz-v4) * theDist->buf[t];
    res += v4 *                theDist->buf[ t + 1 ];
    res += (dy-dydz-v5) *      theDist->buf[ t + theDist->dim[0] ];
    res += v5  *               theDist->buf[ t + theDist->dim[0] + 1 ];
    t += theDist->dim[0]*theDist->dim[1];
    res += (dz-dydz-v6) *      theDist->buf[t];
    res += v6 *                theDist->buf[ t + 1 ];
    res += (dydz-dxdydz) *     theDist->buf[ t + theDist->dim[0] ];
    res += dxdydz *            theDist->buf[ t + theDist->dim[0] + 1 ];
    
    *d = res;
    return( 1 );

  }

  /* on est en dehors de l'image
     
     On estime la distance par
     distance pour rejoindre le bord le plus proche de l'image
     + distance en ce bord
   */

  /* return( 0 ); */

  if ( x < 0.0 ) {
    res += x*theDist->voxelSize[0] * x*theDist->voxelSize[0];
    newx = 0.0;
    ix = 0;
  } else if ( ix >= theDist->dim[0]-1 ) {
    res += (x-(theDist->dim[0]-1))*theDist->voxelSize[0] *
      (x-(theDist->dim[0]-1))*theDist->voxelSize[0];
    newx = theDist->dim[0]-1;
    ix = theDist->dim[0]-1;
  } else {
    newx = x;
  }

  if ( y < 0.0 ) {
    res += y*theDist->voxelSize[1] * y*theDist->voxelSize[1];
    newy = 0.0;
    iy = 0;
  } else if ( iy >= theDist->dim[1]-1 ) {
    res += (y-(theDist->dim[1]-1))*theDist->voxelSize[1] *
      (y-(theDist->dim[1]-1))*theDist->voxelSize[1];
    newy = theDist->dim[1]-1;
    iy = theDist->dim[1]-1;
  } else {
    newy = y;
  }

  if ( z < 0.0 ) {
    res += z*theDist->voxelSize[2] * z*theDist->voxelSize[2];
    newz = 0.0;
    iz = 0;
  } else if ( iz >= theDist->dim[2]-1 ) {
    res += (z-(theDist->dim[2]-1))*theDist->voxelSize[2] *
      (z-(theDist->dim[2]-1))*theDist->voxelSize[2];
    newz = theDist->dim[2]-1;
    iz = theDist->dim[2]-1;
  } else {
    newz = z;
  }

  res = sqrt( res );
  t = iz*theDist->dim[0]*theDist->dim[1] + iy*theDist->dim[0] + ix;

  dx = newx - ix;
  dy = newy - iy;
  dz = newz - iz;


  if ( iz == 0 || iz == theDist->dim[2]-1 ) {
    
    if ( iy == 0 || iy == theDist->dim[1]-1 ) {

      if ( ix == 0 || ix == theDist->dim[0]-1 ) {
	res += theDist->buf[t];
      } else {
	res += (1.0 - dx) * theDist->buf[t] + dx * theDist->buf[t+1];
      }

    } else {
      
      if ( ix == 0 || ix == theDist->dim[0]-1 ) {
	res += (1.0 - dy) * theDist->buf[t] + dy * theDist->buf[t+theDist->dim[0]];
      } else {
	res += (1.0 - dx)*(1.0 - dy) * theDist->buf[t] + 
	  dx*(1.0 - dy) * theDist->buf[t+1] +
	  (1.0 - dx)*dy * theDist->buf[t+theDist->dim[0]] +
	  dx*dy * theDist->buf[t+theDist->dim[0]+1];
      }

    }

  } else {

    if ( iy == 0 || iy == theDist->dim[1]-1 ) {

      if ( ix == 0 || ix == theDist->dim[0]-1 ) {
	res += (1.0 - dz) * theDist->buf[t] + dz * theDist->buf[t+theDist->dim[0]*theDist->dim[1]];
      } else {
	res += (1.0 - dx)*(1.0 - dz) * theDist->buf[t] + 
	  dx*(1.0 - dz) * theDist->buf[t+1] + 
	  (1.0 - dx)*dz * theDist->buf[t+theDist->dim[0]*theDist->dim[1]] +
	  dx*dz * theDist->buf[t+theDist->dim[0]*theDist->dim[1]+1];
      }

    } else {
      
      if ( ix == 0 || ix == theDist->dim[0]-1 ) {
	res += (1.0 - dy)*(1.0 - dz) * theDist->buf[t] + 
	  dy*(1.0 - dz) * theDist->buf[t+theDist->dim[0]] + 
	  (1.0 - dy)*dz * theDist->buf[t+theDist->dim[0]*theDist->dim[1]] +
	  dy*dz * theDist->buf[t+theDist->dim[0]*theDist->dim[1]+theDist->dim[0]];
      } else {
	/* ce cas ne devrait pas survenir
	 */
	res += (1.0 - dz)*((1.0 - dx)*(1.0 - dy) * theDist->buf[t] + 
			   dx*(1.0 - dy) * theDist->buf[t+1] +
			   (1.0 - dx)*dy * theDist->buf[t+theDist->dim[0]] +
			   dx*dy * theDist->buf[t+theDist->dim[0]+1]);
	t += theDist->dim[0]*theDist->dim[1];
	res += dz*((1.0 - dx)*(1.0 - dy) * theDist->buf[t] + 
		   dx*(1.0 - dy) * theDist->buf[t+1] +
		   (1.0 - dx)*dy * theDist->buf[t+theDist->dim[0]] +
		   dx*dy * theDist->buf[t+theDist->dim[0]+1]);
      }
      
    }

  }
  *d = res;
  return( 2 );
}











static void _ComputeTimeFrom2Clocks( float c1, float c2,
				     int *hours, int *minutes, float *seconds )
{
  double d = ( (double)c2 / (double)CLOCKS_PER_SEC ) - 
    ( (double)c1 / (double)CLOCKS_PER_SEC );
  *hours = *minutes = 0;
  *seconds = 0.0;

  if ( d > 3600 ) {
    *hours = (int)(d / 3600);
    d -= *hours * 3600.0;
  }
  if ( d > 60 ) {
    *minutes = (int)(d / 60);
    d -= *minutes * 60.0;
  }
  *seconds = d;
}

static void _PrintTimeFrom2Clocks( float c1, float c2 )
{
  int h, m;
  float s;
  _ComputeTimeFrom2Clocks( c1, c2, &h, &m, &s );
  if ( h > 0 ) printf(" %d h", h);
  if ( m > 0 ) printf(" %d mn", m);
  printf(" %f s\n", s);
}








void _ComputeSignedDistanceMap( typeDistanceMap *theDist,
				const double maxDistUpdate )
{
  char *proc="_ComputeSignedDistanceMap";
  char *blank="                         ";
  typeDistanceCoefficients coeff;
  float exectime[2];

  if ( _VERBOSE_ ) {
    printf("...%s: initialisation\n", proc );
    exectime[0] = (float)clock();
  }

  _InitDistanceCoefficients( &coeff, theDist->voxelSize );
  theDist->multiplicativeCoefficient = coeff.multiplicativeCoefficient;
  _InitBorderSignedDistanceMap( theDist, &coeff );


  if ( _VERBOSE_ ) {
     exectime[1] = (float)clock();
     printf("...%s: ", blank );
     _PrintTimeFrom2Clocks( exectime[0], exectime[1] );
  }
  
  if ( _VERBOSE_ ) {
    printf("...%s: chamfer 3x3x3\n", proc );
    exectime[0] = (float)clock();
  }

  _ComputeSignedDistanceMapWithChamfer3x3x3( theDist, &coeff );

  if ( _VERBOSE_ ) {
     exectime[1] = (float)clock();
     printf("...%s: ", blank );
     _PrintTimeFrom2Clocks( exectime[0], exectime[1] );
  }
  

  if ( maxDistUpdate > 0.0 ) {

    if ( _VERBOSE_ ) {
      printf("...%s: chamfer 5x5x5 (distance <= %f mm)\n", proc, maxDistUpdate );
      exectime[0] = (float)clock();
    }

    _UpdateSignedDistanceMapWithChamfer5x5x5( theDist, &coeff, maxDistUpdate );

    if ( _VERBOSE_ ) {
      exectime[1] = (float)clock();
      printf("...%s: ", blank );
      _PrintTimeFrom2Clocks( exectime[0], exectime[1] );
    }
  }
}











void _InitBorderSignedDistanceMap( typeDistanceMap *theDist,
				   const typeDistanceCoefficients *par )
{
  if ( theDist->dim[2] > 1 ) {
    switch( typeComputation ) {
    default :
    case _THREE_DIMENSIONAL_ :
      _InitBorder3DSignedDistanceMap( theDist, par );
      break;
    case _TWO_DIMENSIONAL_ :
      _InitBorder2DSignedDistanceMap( theDist, par );
    }
  } else {
    _InitBorder2DSignedDistanceMap( theDist, par );
  }
}











/* Initialise les distances a l'interface
   
   L'interieur et l'exterieur sont deja initialises
   a + ou - l'infini : par exemple _MINUS_INFINITY_ et _PLUS_INFINITY_,
   comme cela est realise plus loin.
*/


static void _InitBorder3DSignedDistanceMap( typeDistanceMap *theDist,
				   const typeDistanceCoefficients *par )
{
  int t, x, y, z;
  int i, j, k;
  int tt, ii, jj, kk;

  v777 offset;
  int _insideZ_, _insideYZ_;

  short int *resBuf = theDist->buf;

  const int oy = theDist->dim[0];
  const int oz = theDist->dim[0]*theDist->dim[1];

  const int initMaxIndex = 7;
  const int initMiddleIndex = 3;

  const int initMIMinusOne = initMaxIndex-1;
  const int initMIPlusOne = initMiddleIndex+1;

  const int maxZinside = theDist->dim[2]-initMiddleIndex;
  const int maxYinside = theDist->dim[1]-initMiddleIndex;
  const int maxXinside = theDist->dim[0]-initMiddleIndex;

  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];
  const int dimZ = theDist->dim[2];


  /* initialisation des offsets par rapport au point central
     dans un voisinage 7x7x7
  */
  for ( k=0; k<initMaxIndex; k++ )
  for ( j=0; j<initMaxIndex; j++ )
  for ( i=0; i<initMaxIndex; i++ )
    offset[k][j][i] = (k-initMiddleIndex)*oz +
      (j-initMiddleIndex)*oy + (i-initMiddleIndex);
  

  /* on cherche le cas favorable ou il ne 
     faut pas faire de test
     
     comme on compare  (x,y,z) avec (x-1,y,z)
     il faut x >= initMiddleIndex+1 et x <theDist->dim[0]-initMiddleIndex
     pour que l'on puisse traiter
     (x,y,z) et (x-1,y,z) sans faire de tests
  */




  /* cas 3D
   */
  
  for ( t=0, z=0; z<dimZ; z++ ) {

    _insideZ_ = 0;
    if ( z>initMiddleIndex && z<maxZinside ) _insideZ_ = 1;

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( _insideZ_ == 1 && y>initMiddleIndex && y<maxYinside ) _insideYZ_ = 1;

      for ( x=0; x<dimX; x++, t++ ) {

	if ( _insideYZ_ == 1 && x>initMiddleIndex && x<maxXinside ) {



	  if ( resBuf[t] > 0 ) {

	    /* transition - | + selon X 
	     */
	    if ( resBuf[t-1] < 0 ) {
	      tt = t - 1;
	      for ( i=0, ii=initMIMinusOne; i<initMiddleIndex; i++, ii-- )
              for ( j=i; j<initMaxIndex-i; j++ )
	      for ( k=i; k<initMaxIndex-i; k++ ) {
		if ( resBuf[tt + offset[k][j][ii]] > par->init[k][j][ii] )
		  resBuf[tt + offset[k][j][ii]] = par->init[k][j][ii];
		if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
	      }
	    } /* fin de transition - | + selon X */


	    /* transition - | + selon Y 
	     */
	    if ( resBuf[t-oy] < 0 ) {
	      tt = t-oy;
	      for ( j=0, jj=initMIMinusOne; j<initMiddleIndex; j++, jj-- )
              for ( i=j; i<initMaxIndex-j; i++ )
	      for ( k=j; k<initMaxIndex-j; k++ ) {
		if ( resBuf[tt + offset[k][jj][i]] > par->init[k][jj][i] )
		  resBuf[tt + offset[k][jj][i]] = par->init[k][jj][i];
		if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
	      }
	    } /* fin de transition - | + selon Y */


	    /* transition - | + selon Z 
	     */
	    if ( resBuf[t-oz] < 0 ) {
	      tt = t-oz;
	      for ( k=0, kk=initMIMinusOne; k<initMiddleIndex; k++, kk-- )
              for ( j=k; j<initMaxIndex-k; j++ )
	      for ( i=k; i<initMaxIndex-k; i++ ) {
		if ( resBuf[tt + offset[kk][j][i]] > par->init[kk][j][i] )
		  resBuf[tt + offset[kk][j][i]] = par->init[kk][j][i];
		if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
	      }
	    }
	    

	  } else {
	    /* ici resBuf[t] est negatif 
	     */

	    /* transition + | - selon X 
	     */
	    if ( resBuf[t-1] > 0 ) {
	      tt = t-1;
	      for ( i=0, ii=initMIMinusOne; i<initMiddleIndex; i++, ii-- )
              for ( j=i; j<initMaxIndex-i; j++ )
	      for ( k=i; k<initMaxIndex-i; k++ ) {
		if ( resBuf[tt + offset[k][j][ii]] < -par->init[k][j][ii] )
		  resBuf[tt + offset[k][j][ii]] = -par->init[k][j][ii];
		if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = par->init[k][j][i];
	      }
	    } /* fin de transition + | - selon X */


	    /* transition + | - selon Y 
	     */
	    if ( resBuf[t-oy] > 0 ) {
	      tt = t-oy;
	      for ( j=0, jj=initMIMinusOne; j<initMiddleIndex; j++, jj-- )
              for ( i=j; i<initMaxIndex-j; i++ )
	      for ( k=j; k<initMaxIndex-j; k++ ) {
		if ( resBuf[tt + offset[k][jj][i]] < -par->init[k][jj][i] )
		  resBuf[tt + offset[k][jj][i]] = -par->init[k][jj][i];
		if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = par->init[k][j][i];
	      }
	    } /* fin de transition + | - selon Y */


	    /* transition + | - selon Z 
	     */
	    if ( resBuf[t-oz] > 0 ) {
	      tt = t-oz;
	      for ( k=0, kk=initMIMinusOne; k<initMiddleIndex; k++, kk-- )
              for ( j=k; j<initMaxIndex-k; j++ )
	      for ( i=k; i<initMaxIndex-k; i++ ) {
		if ( resBuf[tt + offset[kk][j][i]] < -par->init[kk][j][i] )
		  resBuf[tt + offset[kk][j][i]] = -par->init[kk][j][i];
		if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = par->init[k][j][i];
	      }
	      
	    } /* fin de transition + | - selon Z */

	  }




	} else {
	  /* ici, on sait qu'il faut tout tester
	   */


	  if ( x > 0 ) {
	    tt = t-1;

	    if ( resBuf[t] > 0 ) {
	      if ( resBuf[tt] < 0 ) {
		/* transition - | + selon X 
		 */
		/* Tue Jul  6 17:08:21 MET DST 1999
		   pour tt (ici [z,y,x-1]), il faut faire le test sur [z,y,x-1] et non
		   sur [z,y,x]
		   il faut changer x+i-initMiddleIndex 
		                en x-1+i-initMiddleIndex 
			      soit x+i-initMIPlusOne
		   id. pour y et z
		*/
		for ( i=initMIPlusOne; i<initMaxIndex; i++ ) {
		  if ( x+i-initMIPlusOne < 0 || x+i-initMIPlusOne >= dimX ) continue;
		  for ( j=initMIMinusOne-i; j<=i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( k=initMIMinusOne-i; k<=i; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[tt + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}

		for ( i=0; i<initMiddleIndex; i++ ) {
		  if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		  for ( j=i; j<initMaxIndex-i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( k=i; k<initMaxIndex-i; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}
		

	      } /* fin de transition - | + selon X */
	    } else {
	      if ( resBuf[tt] > 0 ) {
		/* transition + | - selon X 
		 */
		
		for ( i=initMIPlusOne; i<initMaxIndex; i++ ) {
		  if ( x+i-initMIPlusOne < 0 || x+i-initMIPlusOne >= dimX ) continue;
		  for ( j=initMIMinusOne-i; j<=i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( k=initMIMinusOne-i; k<=i; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[tt + offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}

		for ( i=0; i<initMiddleIndex; i++ ) {
		  if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		  for ( j=i; j<initMaxIndex-i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( k=i; k<initMaxIndex-i; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}

	      } /* fin de transition + | - selon X */
	    }
	  } /* fin des transitions selon X */
	  



	  if ( y > 0 ) {
	    tt = t-oy;
	    if ( resBuf[t] > 0 ) {
	      if ( resBuf[tt] < 0 ) {
		/* transition - | + selon Y 
		 */

		for ( j=initMIPlusOne; j<initMaxIndex; j++ ) {
		  if ( y+j-initMIPlusOne < 0 || y+j-initMIPlusOne >= dimY ) continue;
		  for ( i=initMIMinusOne-j; i<=j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    for ( k=initMIMinusOne-j; k<=j; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[tt + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}

		for ( j=0; j<initMiddleIndex; j++ ) {
		  if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		  for ( i=j; i<initMaxIndex-j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    for ( k=j; k<initMaxIndex-j; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}

	      } /* fin de transition - | + selon Y */
	    } else {
	      if ( resBuf[tt] > 0 ) {
		/* transition + | - selon Y 
		 */
		
		for ( j=initMIPlusOne; j<initMaxIndex; j++ ) {
		  if ( y+j-initMIPlusOne < 0 || y+j-initMIPlusOne >= dimY ) continue;
		  for ( i=initMIMinusOne-j; i<=j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    for ( k=initMIMinusOne-j; k<=j; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[tt+ offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}

		for ( j=0; j<initMiddleIndex; j++ ) {
		  if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		  for ( i=j; i<initMaxIndex-j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    for ( k=j; k<initMaxIndex-j; k++ ) {
		      if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		      if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}

	      } /* fin de transition + | - selon Y */
	    }
	  } /* fin des transitions selon Y */
	  


	  if ( z > 0 ) {
	    tt = t-oz;
	    if ( resBuf[t] > 0 ) {
	      if ( resBuf[tt] < 0 ) {
		/* transition - | + selon Z 
		 */


		for ( k=initMIPlusOne; k<initMaxIndex; k++ ) {
		  if ( z+k-initMIPlusOne < 0 || z+k-initMIPlusOne >= dimZ ) continue;
		  for ( j=initMIMinusOne-k; j<=k; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( i=initMIMinusOne-k; i<=k; i++ ) {
		      if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		      if ( resBuf[tt + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}
		  

		for ( k=0; k<initMiddleIndex; k++ ) {
		  if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		  for ( j=k; j<initMaxIndex-k; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( i=k; i<initMaxIndex-k; i++ ) {
		      if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		      if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}

	      } /* fin de transition - | + selon Z */
	    } else {
	      if ( resBuf[tt] > 0 ) {
		/* transition + | - selon Z 
		 */
		
		for ( k=initMIPlusOne; k<initMaxIndex; k++ ) {
		  if ( z+k-initMIPlusOne < 0 || z+k-initMIPlusOne >= dimZ ) continue;
		  for ( j=initMIMinusOne-k; j<=k; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( i=initMIMinusOne-k; i<=k; i++ ) {
		      if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		      if ( resBuf[tt + offset[k][j][i]] < -par->init[k][j][i] )
			resBuf[tt + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }
		}

		for ( k=0; k<initMiddleIndex; k++ ) {
		  if ( z+k-initMiddleIndex < 0 || z+k-initMiddleIndex >= dimZ ) continue;
		  for ( j=k; j<initMaxIndex-k; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    for ( i=k; i<initMaxIndex-k; i++ ) {
		      if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		      if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
			resBuf[t + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }
		}

	      } /* fin de transition + | - selon Z */
	    }
	  } /* fin des transitions selon Z */


	}


      } /* fin de la boucle sur les x */
    }
  }
}












/* Initialise les distances a l'interface
   
   L'interieur et l'exterieur sont deja initialises
   a + ou - l'infini : par exemple _MINUS_INFINITY_ et _PLUS_INFINITY_,
   comme cela est realise plus loin.
*/


static void _InitBorder2DSignedDistanceMap( typeDistanceMap *theDist,
				   const typeDistanceCoefficients *par )
{
  int t, x, y, z;
  int i, j, k;
  int tt, ii, jj;

  v777 offset;
  int _insideYZ_;

  short int *resBuf = theDist->buf;

  const int oy = theDist->dim[0];
  const int oz = theDist->dim[0]*theDist->dim[1];

  const int initMaxIndex = 7;
  const int initMiddleIndex = 3;

  const int initMIMinusOne = initMaxIndex-1;
  const int initMIPlusOne = initMiddleIndex+1;

  const int maxYinside = theDist->dim[1]-initMiddleIndex;
  const int maxXinside = theDist->dim[0]-initMiddleIndex;

  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];


  /* initialisation des offsets par rapport au point central
     dans un voisinage 7x7x7
  */
  for ( k=0; k<initMaxIndex; k++ )
  for ( j=0; j<initMaxIndex; j++ )
  for ( i=0; i<initMaxIndex; i++ )
    offset[k][j][i] = (k-initMiddleIndex)*oz +
      (j-initMiddleIndex)*oy + (i-initMiddleIndex);
  

  /* on cherche le cas favorable ou il ne 
     faut pas faire de test
     
     comme on compare  (x,y,z) avec (x-1,y,z)
     il faut x >= initMiddleIndex+1 et x <theDist->dim[0]-initMiddleIndex
     pour que l'on puisse traiter
     (x,y,z) et (x-1,y,z) sans faire de tests
  */



  k = initMiddleIndex;
  /* cas 2D
   */

  /* on sait qu'il n'y aura rien a faire
     pour les points d'intensite >= theDist.intensityMax
  */
  for ( t=0, z=0; z<theDist->intensityMax && z<theDist->dim[2]; z++ ) {

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( y>initMiddleIndex && y<maxYinside ) _insideYZ_ = 1;

      for ( x=0; x<dimX; x++, t++ ) {

	if ( _insideYZ_ == 1 && x>initMiddleIndex && x<maxXinside ) {



	  if ( resBuf[t] > 0 ) {

	    /* transition - | + selon X 
	     */
	    if ( resBuf[t-1] < 0 ) {
	      tt = t - 1;
	      for ( i=0, ii=initMIMinusOne; i<initMiddleIndex; i++, ii-- )
	      for ( j=i; j<initMaxIndex-i; j++ ) {
		if ( resBuf[tt + offset[k][j][ii]] > par->init[k][j][ii] )
		  resBuf[tt + offset[k][j][ii]] = par->init[k][j][ii];
		if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
	      }
	    } /* fin de transition - | + selon X */


	    /* transition - | + selon Y 
	     */
	    if ( resBuf[t-oy] < 0 ) {
	      tt = t-oy;
	      for ( j=0, jj=initMIMinusOne; j<initMiddleIndex; j++, jj-- )
	      for ( i=j; i<initMaxIndex-j; i++ ) {
		if ( resBuf[tt + offset[k][jj][i]] > par->init[k][jj][i] )
		  resBuf[tt + offset[k][jj][i]] = par->init[k][jj][i];
		if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
	      }
	    } /* fin de transition - | + selon Y */


	  } else {
	    /* ici resBuf[t] est negatif 
	     */

	    /* transition + | - selon X 
	     */
	    if ( resBuf[t-1] > 0 ) {
	      tt = t-1;
	      for ( i=0, ii=initMIMinusOne; i<initMiddleIndex; i++, ii-- )
	      for ( j=i; j<initMaxIndex-i; j++ ) {
		if ( resBuf[tt + offset[k][j][ii]] < -par->init[k][j][ii] )
		  resBuf[tt + offset[k][j][ii]] = -par->init[k][j][ii];
		if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = par->init[k][j][i];
	      }
	    } /* fin de transition + | - selon X */


	    /* transition + | - selon Y 
	     */
	    if ( resBuf[t-oy] > 0 ) {
	      tt = t-oy;
	      for ( j=0, jj=initMIMinusOne; j<initMiddleIndex; j++, jj-- )
              for ( i=j; i<initMaxIndex-j; i++ ) {
		if ( resBuf[tt + offset[k][jj][i]] < -par->init[k][jj][i] )
		  resBuf[tt + offset[k][jj][i]] = -par->init[k][jj][i];
		if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		  resBuf[t + offset[k][j][i]] = par->init[k][j][i];
	      }
	    } /* fin de transition + | - selon Y */


	  }




	} else {
	  /* ici, on sait qu'il faut tout tester
	   */


	  if ( x > 0 ) {
	    tt = t-1;

	    if ( resBuf[t] > 0 ) {
	      if ( resBuf[tt] < 0 ) {
		/* transition - | + selon X 
		 */
		/* Tue Jul  6 17:08:21 MET DST 1999
		   pour tt (ici [z,y,x-1]), il faut faire le test sur [z,y,x-1] et non
		   sur [z,y,x]
		   il faut changer x+i-initMiddleIndex 
		                en x-1+i-initMiddleIndex 
			      soit x+i-initMIPlusOne
		   id. pour y et z
		*/
		for ( i=initMIPlusOne; i<initMaxIndex; i++ ) {
		  if ( x+i-initMIPlusOne < 0 || x+i-initMIPlusOne >= dimX ) continue;
		  for ( j=initMIMinusOne-i; j<=i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    if ( resBuf[tt + offset[k][j][i]] > par->init[k][j][i] )
		      resBuf[tt + offset[k][j][i]] = par->init[k][j][i];
		  }
		}

		for ( i=0; i<initMiddleIndex; i++ ) {
		  if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		  for ( j=i; j<initMaxIndex-i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		      resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
		  }
		}
		

	      } /* fin de transition - | + selon X */
	    } else {
	      if ( resBuf[tt] > 0 ) {
		/* transition + | - selon X 
		 */
		
		for ( i=initMIPlusOne; i<initMaxIndex; i++ ) {
		  if ( x+i-initMIPlusOne < 0 || x+i-initMIPlusOne >= dimX ) continue;
		  for ( j=initMIMinusOne-i; j<=i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    if ( resBuf[tt + offset[k][j][i]] < -par->init[k][j][i] )
		      resBuf[tt + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }

		for ( i=0; i<initMiddleIndex; i++ ) {
		  if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		  for ( j=i; j<initMaxIndex-i; j++ ) {
		    if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		    if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		      resBuf[t + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }

	      } /* fin de transition + | - selon X */
	    }
	  } /* fin des transitions selon X */
	  



	  if ( y > 0 ) {
	    tt = t-oy;
	    if ( resBuf[t] > 0 ) {
	      if ( resBuf[tt] < 0 ) {
		/* transition - | + selon Y 
		 */

		for ( j=initMIPlusOne; j<initMaxIndex; j++ ) {
		  if ( y+j-initMIPlusOne < 0 || y+j-initMIPlusOne >= dimY ) continue;
		  for ( i=initMIMinusOne-j; i<=j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    if ( resBuf[tt + offset[k][j][i]] > par->init[k][j][i] )
		      resBuf[tt + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }

		for ( j=0; j<initMiddleIndex; j++ ) {
		  if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		  for ( i=j; i<initMaxIndex-j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    if ( resBuf[t + offset[k][j][i]] < -par->init[k][j][i] )
		      resBuf[t + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }

	      } /* fin de transition - | + selon Y */
	    } else {
	      if ( resBuf[tt] > 0 ) {
		/* transition + | - selon Y 
		 */
		
		for ( j=initMIPlusOne; j<initMaxIndex; j++ ) {
		  if ( y+j-initMIPlusOne < 0 || y+j-initMIPlusOne >= dimY ) continue;
		  for ( i=initMIMinusOne-j; i<=j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    if ( resBuf[tt+ offset[k][j][i]] < -par->init[k][j][i] )
		      resBuf[tt + offset[k][j][i]] = -par->init[k][j][i];
		    }
		  }

		for ( j=0; j<initMiddleIndex; j++ ) {
		  if ( y+j-initMiddleIndex < 0 || y+j-initMiddleIndex >= dimY ) continue;
		  for ( i=j; i<initMaxIndex-j; i++ ) {
		    if ( x+i-initMiddleIndex < 0 || x+i-initMiddleIndex >= dimX ) continue;
		    if ( resBuf[t + offset[k][j][i]] > par->init[k][j][i] )
		      resBuf[t + offset[k][j][i]] = par->init[k][j][i];
		    }
		  }

	      } /* fin de transition + | - selon Y */
	    }
	  } /* fin des transitions selon Y */
	  



	}


      } /* fin de la boucle sur les x */
    }
  }
}















void _ComputeSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par )
{
  if ( theDist->dim[2] > 1 ) {
    switch( typeComputation ) {
    default :
    case _THREE_DIMENSIONAL_ :
      _Compute3DSignedDistanceMapWithChamfer3x3x3( theDist, par );
      break;
    case _TWO_DIMENSIONAL_ :
      _Compute2DSignedDistanceMapWithChamfer3x3x3( theDist, par );
    }
  } else {
    _Compute2DSignedDistanceMapWithChamfer3x3x3( theDist, par );
  }
}








static void _Compute3DSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par )
{
  int t, x, y, z;
  int i, j, k;
  v333 offset;
  int _insideZ_, _insideYZ_;
  int currentEval, bestEval;

  short int *resBuf = theDist->buf;

  const int offsMaxIndex = 3;
  const int offsMiddleIndex = 1;
  
  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];
  const int dimZ = theDist->dim[2];
  
  const int dimXPlusOne = dimX + 1;
  const int dimYPlusOne = dimY + 1;

  const int dimXMinusOne = dimX - 1;
  const int dimYMinusOne = dimY - 1;
  const int dimZMinusOne = dimZ - 1;
  
 
  /* initialisation des offsets par rapport au point central
     dans un voisinage 3x3x3
  */
  for (k=0;k<offsMaxIndex;k++)
  for (j=0;j<offsMaxIndex;j++)
  for (i=0;i<offsMaxIndex;i++)
    offset[k][j][i] = (k-offsMiddleIndex)*dimX*dimY +
      (j-offsMiddleIndex)*dimX + (i-offsMiddleIndex);



  /* boucle forward
   */
  for ( t=0, z=0; z<dimZ; z++ ) {

    _insideZ_ = 0;
    if ( z>=offsMiddleIndex ) _insideZ_ = 1;

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( _insideZ_ == 1 && y>=offsMiddleIndex && y<dimYMinusOne ) _insideYZ_ = 1;
    
      for ( x=0; x<dimX; x++, t++ ) {

	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddleIndex && x<dimXMinusOne ) {



	  if ( resBuf[t] > par->maxInit ) {

	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[0][j][i]] + par->incr[1][1+j][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	    if ( currentEval < bestEval )
	      resBuf[t] = currentEval;

	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[0][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[0][j][i]] + par->incr[1][1+j][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][0][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][0]] > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[0][j][i]] - par->incr[1][1+j][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	    if ( currentEval > bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[0][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[0][j][i]] - par->incr[1][1+j][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][0][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][0]] < 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	  } 




	} else {




	  if ( resBuf[t] > par->maxInit ) {

	    if ( z > 0 ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  currentEval = resBuf[t + offset[0][j][i]] + par->incr[1][1+j][1+i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {
	    
	    if ( z > 0 ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  if ( resBuf[t + offset[0][j][i]] < 0 ) continue;
		  currentEval = resBuf[t + offset[0][j][i]] + par->incr[1][1+j][1+i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][0][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      if ( resBuf[t + offset[1][1][0]] > 0 ) {
		currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
		if ( currentEval < bestEval )
		  resBuf[t] = currentEval;
	      }
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    if ( z > 0 ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  currentEval = resBuf[t + offset[0][j][i]] - par->incr[1][1+j][1+i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    if ( z > 0 ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  if ( resBuf[t + offset[0][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[0][j][i]] - par->incr[1][1+j][1+i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( resBuf[t + offset[1][0][i]] > 0 ) continue;
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      if ( resBuf[t + offset[1][1][0]] < 0 ) {
		currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
		if ( currentEval > bestEval )
		  resBuf[t] = currentEval;
	      }
	    } 
	  }
	  
	}


      } /* fin de la boucle sur les x */
    }
  }




  /* boucle backward
   */
  for ( t=dimZ*dimY*dimX-1, 
	  z=dimZMinusOne; z>=0; z-- ) {

    _insideZ_ = 0;
    if ( z<dimZMinusOne ) _insideZ_ = 1;

    for ( y=dimYMinusOne; y>=0; y-- ) {

      _insideYZ_ = 0;
      if ( _insideZ_ == 1 && y>=1 && y<dimYMinusOne ) _insideYZ_ = 1;
    
      for ( x=dimXMinusOne; x>=0; x--, t-- ) {
	bestEval = resBuf[t];


	if ( _insideYZ_ == 1 && x>=1 && x<dimXMinusOne ) {

	  
	  
	  if ( resBuf[t] > par->maxInit ) {

	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[2][j][i]] + par->incr[3][1+j][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	    if ( currentEval < bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] > par->minInit ) {

	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[2][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[2][j][i]] + par->incr[3][1+j][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][2][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][2]] > 0 ) {
	      currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[2][j][i]] - par->incr[3][1+j][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	    if ( currentEval > bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( j=0; j<offsMaxIndex; j++ )
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[2][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[2][j][i]] - par->incr[3][1+j][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][2][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][2]] < 0 ) {
	      currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  }




	} else {





	  if ( resBuf[t] > par->maxInit ) {

	    if ( z<dimZMinusOne ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  currentEval = resBuf[t + offset[2][j][i]] + par->incr[3][1+j][1+i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    if ( z<dimZMinusOne ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  if ( resBuf[t + offset[2][j][i]] < 0 ) continue;
		  currentEval = resBuf[t + offset[2][j][i]] + par->incr[3][1+j][1+i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][2][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      if ( resBuf[t + offset[1][1][2]] > 0 ) {
		currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
		if ( currentEval < bestEval )
		  resBuf[t] = currentEval;
	      }
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    if ( z<dimZMinusOne ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  currentEval = resBuf[t + offset[2][j][i]] - par->incr[3][1+j][1+i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {

	    if ( z<dimZMinusOne ) {
	      for ( j=0; j<offsMaxIndex; j++ ) {
		if ( y+j < 1 || y+j >= dimYPlusOne ) continue;
		for ( i=0; i<offsMaxIndex; i++ ) {
		  if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		  if ( resBuf[t + offset[2][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[2][j][i]] - par->incr[3][1+j][1+i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][2][i]] > 0 ) continue;
		currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      if ( resBuf[t + offset[1][1][2]] < 0 ) {
		currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
		if ( currentEval > bestEval )
		  resBuf[t] = currentEval;
	      }
	    }
	  }




	}




      } /* fin de la boucle sur les x */
    }
  }
}
  













static void _Compute2DSignedDistanceMapWithChamfer3x3x3( typeDistanceMap *theDist,
					  const typeDistanceCoefficients *par )
{
  int t, x, y, z;
  int i, j, k;
  v333 offset;
  int _insideYZ_;
  int currentEval, bestEval;

  short int *resBuf = theDist->buf;

  const int offsMaxIndex = 3;
  const int offsMiddleIndex = 1;
  
  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];
  const int dimZ = theDist->dim[2];
  
  const int dimXPlusOne = dimX + 1;

  const int dimXMinusOne = dimX - 1;
  const int dimYMinusOne = dimY - 1;
  
 
  /* initialisation des offsets par rapport au point central
     dans un voisinage 3x3x3
  */
  for (k=0;k<offsMaxIndex;k++)
  for (j=0;j<offsMaxIndex;j++)
  for (i=0;i<offsMaxIndex;i++)
    offset[k][j][i] = (k-offsMiddleIndex)*dimX*dimY +
      (j-offsMiddleIndex)*dimX + (i-offsMiddleIndex);



  /* boucle forward
   */
  for ( t=0, z=0; z<theDist->intensityMax && z<dimZ; z++ ) {

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( y>=offsMiddleIndex && y<dimYMinusOne ) _insideYZ_ = 1;
    
      for ( x=0; x<dimX; x++, t++ ) {

	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddleIndex && x<dimXMinusOne ) {



	  if ( resBuf[t] > par->maxInit ) {

	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	    if ( currentEval < bestEval )
	      resBuf[t] = currentEval;

	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][0][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][0]] > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	    if ( currentEval > bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][0][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][0]] < 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	  } 




	} else {




	  if ( resBuf[t] > par->maxInit ) {

	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {
	    
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][0][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[1][0][i]] + par->incr[2][1][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      if ( resBuf[t + offset[1][1][0]] > 0 ) {
		currentEval = resBuf[t + offset[1][1][0]] + par->incr[2][2][1];
		if ( currentEval < bestEval )
		  resBuf[t] = currentEval;
	      }
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    if ( y > 0 ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( resBuf[t + offset[1][0][i]] > 0 ) continue;
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][0][i]] - par->incr[2][1][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x > 0 ) {
	      if ( resBuf[t + offset[1][1][0]] < 0 ) {
		currentEval = resBuf[t + offset[1][1][0]] - par->incr[2][2][1];
		if ( currentEval > bestEval )
		  resBuf[t] = currentEval;
	      }
	    } 
	  }
	  
	}


      } /* fin de la boucle sur les x */
    }
  }




  /* boucle backward
   */
  t=theDist->intensityMax*dimY*dimX-1;
  z=theDist->intensityMax-1;
  if ( theDist->intensityMax > dimZ ) {
    t=dimZ*dimY*dimX-1;
    z=dimZ-1;
  }
    
  for ( ; z>=0; z-- ) {

    for ( y=dimYMinusOne; y>=0; y-- ) {

      _insideYZ_ = 0;
      if ( y>=1 && y<dimYMinusOne ) _insideYZ_ = 1;
    
      for ( x=dimXMinusOne; x>=0; x--, t-- ) {
	bestEval = resBuf[t];


	if ( _insideYZ_ == 1 && x>=1 && x<dimXMinusOne ) {

	  
	  
	  if ( resBuf[t] > par->maxInit ) {

	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	    if ( currentEval < bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] > par->minInit ) {

	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][2][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][2]] > 0 ) {
	      currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( i=0; i<offsMaxIndex; i++ ) {
	      currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	    if ( currentEval > bestEval )
	      resBuf[t] = currentEval;
	    
	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( i=0; i<offsMaxIndex; i++ ) {
	      if ( resBuf[t + offset[1][2][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    if ( resBuf[t + offset[1][1][2]] < 0 ) {
	      currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  }




	} else {





	  if ( resBuf[t] > par->maxInit ) {

	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][2][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[1][2][i]] + par->incr[2][3][1+i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      if ( resBuf[t + offset[1][1][2]] > 0 ) {
		currentEval = resBuf[t + offset[1][1][2]] + par->incr[2][2][3];
		if ( currentEval < bestEval )
		  resBuf[t] = currentEval;
	      }
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {

	    if ( y<dimYMinusOne ) {
	      for ( i=0; i<offsMaxIndex; i++ ) {
		if ( x+i < 1 || x+i >= dimXPlusOne ) continue;
		if ( resBuf[t + offset[1][2][i]] > 0 ) continue;
		currentEval = resBuf[t + offset[1][2][i]] - par->incr[2][3][1+i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    if ( x<dimXMinusOne ) {
	      if ( resBuf[t + offset[1][1][2]] < 0 ) {
		currentEval = resBuf[t + offset[1][1][2]] - par->incr[2][2][3];
		if ( currentEval > bestEval )
		  resBuf[t] = currentEval;
	      }
	    }
	  }




	}




      } /* fin de la boucle sur les x */
    }
  }
}
  




















void _UpdateSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate )
{
  if ( theDist->dim[2] > 1 ) {
    switch( typeComputation ) {
    default :
    case _THREE_DIMENSIONAL_ :
      _Update3DSignedDistanceMapWithChamfer5x5x5( theDist, par, maxDistUpdate );
      break;
    case _TWO_DIMENSIONAL_ :
      _Update2DSignedDistanceMapWithChamfer5x5x5( theDist, par, maxDistUpdate );
    }
  } else {
    _Update2DSignedDistanceMapWithChamfer5x5x5( theDist, par, maxDistUpdate );
  }
}












static void _Update3DSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate )
{
  int t, x, y, z;
  int i, j, k;
  v555 offset;
  int _insideZ_, _insideYZ_;
  int currentEval, bestEval;

  short int *resBuf = theDist->buf;

  const int offsMax = 5;
  const int offsMiddle = 2;
  const int offsMiddlePlusOne = offsMiddle + 1;
  
  const int maxValUpdate = (int)( maxDistUpdate / par->multiplicativeCoefficient + 0.5 );
  
  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];
  const int dimZ = theDist->dim[2];

  const int dimXPlusTwo = dimX + 2;
  const int dimYPlusTwo = dimY + 2;
  
  const int dimXMinusTwo = dimX - 2;
  const int dimYMinusTwo = dimY - 2;
  const int dimZMinusTwo = dimZ - 2;
  
  /* initialisation des offsets par rapport au point central
     dans un voisinage 3x3x3
  */
  for (k=0;k<offsMax;k++)
  for (j=0;j<offsMax;j++)
  for (i=0;i<offsMax;i++)
    offset[k][j][i] = (k-offsMiddle)*dimX*dimY +
      (j-offsMiddle)*dimX + (i-offsMiddle);



  /* boucle forward
   */
  for ( t=0, z=0; z<dimZ; z++ ) {

    _insideZ_ = 0;
    if ( z>=offsMiddle ) _insideZ_ = 1;

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( _insideZ_ == 1 && y>=offsMiddle && y<dimYMinusTwo ) _insideYZ_ = 1;
    
      for ( x=0; x<dimX; x++, t++ ) {

	if ( resBuf[t] > maxValUpdate ) continue;
	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddle && x<dimXMinusTwo ) {



	  if ( resBuf[t] > par->maxInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[k][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[k][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }



	} else {




	  if ( resBuf[t] > par->maxInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ ) {
	      if ( z+k < offsMiddle ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ ) {
	      if ( z+k < offsMiddle ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[k][j][i]] < 0 ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ ) {
	      if ( z+k < offsMiddle ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( k=0; k<offsMiddle; k++ ) {
	      if ( z+k < offsMiddle ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[k][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }


	}


      } /* fin de la boucle sur les x */
    }
  }



  /* boucle backward
   */
  for ( t=dimZ*dimY*dimX-1, 
	  z=dimZ-1; z>=0; z-- ) {

    _insideZ_ = 0;
    if ( z<dimZMinusTwo ) _insideZ_ = 1;

    for ( y=dimY-1; y>=0; y-- ) {

      _insideYZ_ = 0;
      if ( _insideZ_ == 1 && y>=offsMiddle && y<dimYMinusTwo ) _insideYZ_ = 1;
    
      for ( x=dimX-1; x>=0; x--, t-- ) {

	if ( resBuf[t] > maxValUpdate ) continue;
	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddle && x<dimXMinusTwo ) {

	  
	  
	  if ( resBuf[t] > par->maxInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[k][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ )
	    for ( j=0; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[k][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }




	} else {




	  if ( resBuf[t] > par->maxInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ ) {
	      if ( z+k-offsMiddle >= dimZ ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ ) {
	      if ( z+k-offsMiddle >= dimZ ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[k][j][i]] < 0 ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] + par->incr[k][j][i];
		  if ( currentEval < bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ ) {
	      if ( z+k-offsMiddle >= dimZ ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( k=offsMiddlePlusOne; k<offsMax; k++ ) {
	      if ( z+k-offsMiddle >= dimZ ) continue;
	      for ( j=0; j<offsMax; j++ ) {
		if ( y+j < offsMiddle || y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[k][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[k][j][i]] - par->incr[k][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	      }
	    }
	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
		if ( y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }



	}



      } /* fin de la boucle sur les x */
    }
  }
}
  
























static void _Update2DSignedDistanceMapWithChamfer5x5x5( typeDistanceMap *theDist,
					       const typeDistanceCoefficients *par,
					       const double maxDistUpdate )
{
  int t, x, y, z;
  int i, j, k;
  v555 offset;
  int _insideYZ_;
  int currentEval, bestEval;

  short int *resBuf = theDist->buf;

  const int offsMax = 5;
  const int offsMiddle = 2;
  const int offsMiddlePlusOne = offsMiddle + 1;
  
  const int maxValUpdate = (int)( maxDistUpdate / par->multiplicativeCoefficient + 0.5 );
  
  const int dimX = theDist->dim[0];
  const int dimY = theDist->dim[1];
  const int dimZ = theDist->dim[2];

  const int dimXPlusTwo = dimX + 2;
  const int dimYPlusTwo = dimY + 2;
  
  const int dimXMinusTwo = dimX - 2;
  const int dimYMinusTwo = dimY - 2;
  
  /* initialisation des offsets par rapport au point central
     dans un voisinage 3x3x3
  */
  for (k=0;k<offsMax;k++)
  for (j=0;j<offsMax;j++)
  for (i=0;i<offsMax;i++)
    offset[k][j][i] = (k-offsMiddle)*dimX*dimY +
      (j-offsMiddle)*dimX + (i-offsMiddle);



  /* boucle forward
   */
  for ( t=0, z=0; z<theDist->intensityMax && z<dimZ; z++ ) {

    for ( y=0; y<dimY; y++ ) {

      _insideYZ_ = 0;
      if ( y>=offsMiddle && y<dimYMinusTwo ) _insideYZ_ = 1;
    
      for ( x=0; x<dimX; x++, t++ ) {

	if ( resBuf[t] > maxValUpdate ) continue;
	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddle && x<dimXMinusTwo ) {



	  if ( resBuf[t] > par->maxInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }



	} else {




	  if ( resBuf[t] > par->maxInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] > par->minInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }
	    
	  } else if ( resBuf[t] < -par->minInit ) {
	    
	    for ( j=0; j<offsMiddle; j++ ) {
	      if ( y+j < offsMiddle ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=0; i<offsMiddle; i++ ) {
	      if ( x+i < offsMiddle ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }


	}


      } /* fin de la boucle sur les x */
    }
  }



  /* boucle backward
   */
  for ( t=dimZ*dimY*dimX-1, 
	  z=dimZ-1; z>=0; z-- ) {

    for ( y=dimY-1; y>=0; y-- ) {

      _insideYZ_ = 0;
      if ( y>=offsMiddle && y<dimYMinusTwo ) _insideYZ_ = 1;
    
      for ( x=dimX-1; x>=0; x--, t-- ) {

	if ( resBuf[t] > maxValUpdate ) continue;
	bestEval = resBuf[t];

	if ( _insideYZ_ == 1 && x>=offsMiddle && x<dimXMinusTwo ) {

	  
	  
	  if ( resBuf[t] > par->maxInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
	      if ( currentEval < bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ )
	    for ( i=0; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
	      if ( currentEval > bestEval )
		bestEval = resBuf[t] = currentEval;
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }




	} else {




	  if ( resBuf[t] > par->maxInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] > par->minInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		if ( resBuf[t + offset[offsMiddle][j][i]] < 0 ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] + par->incr[offsMiddle][j][i];
		if ( currentEval < bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] < 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] + par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval < bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->maxInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
	      if ( y+j >= dimYPlusTwo ) continue;
	      for ( i=0; i<offsMax; i++ ) {
		if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		if ( currentEval > bestEval )
		  bestEval = resBuf[t] = currentEval;
	      }
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  } else if ( resBuf[t] < -par->minInit ) {

	    for ( j=offsMiddlePlusOne; j<offsMax; j++ ) {
		if ( y+j >= dimYPlusTwo ) continue;
		for ( i=0; i<offsMax; i++ ) {
		  if ( x+i < offsMiddle || x+i >= dimXPlusTwo ) continue;
		  if ( resBuf[t + offset[offsMiddle][j][i]] > 0 ) continue;
		  currentEval = resBuf[t + offset[offsMiddle][j][i]] - par->incr[offsMiddle][j][i];
		  if ( currentEval > bestEval )
		    bestEval = resBuf[t] = currentEval;
		}
	    }
	    for ( i=offsMiddlePlusOne; i<offsMax; i++ ) {
	      if ( x+i >= dimXPlusTwo ) continue;
	      if ( resBuf[t + offset[offsMiddle][offsMiddle][i]] > 0 ) continue;
	      currentEval = resBuf[t + offset[offsMiddle][offsMiddle][i]] - par->incr[offsMiddle][offsMiddle][i];
	      if ( currentEval > bestEval )
		resBuf[t] = currentEval;
	    }

	  }



	}



      } /* fin de la boucle sur les x */
    }
  }
}
  
























/* prepare la carte de distance a partir
   d'une image

   met les points a 0 a -infini : ce seront des 
   distances negatives

   met les points > 0 a +infini : ce seront des 
   distances positives
*/   

#define _MINUS_INFINITY_ -32768
#define _PLUS_INFINITY_   32767

void _InitSignedDistanceMap( const unsigned char *theBuf,
		       short int *resBuf,
		       const int *theDimension )
{
  int t, x, y, z;
  for ( t=0, z=0; z<theDimension[2]; z++ )
  for ( y=0; y<theDimension[1]; y++ )
  for ( x=0; x<theDimension[0]; x++, t++ ) {
    if ( theBuf[t] > 0 ) resBuf[t] = _PLUS_INFINITY_;
    else resBuf[t] = _MINUS_INFINITY_;
  }

}









/* Procedure d'initialisation des coefficients pour le calcul de distance
 

   Ce sont les coefficients pour initialiser
   les distances a l'interface entre le
   fond et l'objet. 
   Si la taille du voxel est 1, et si deux voxels
   6-connexes appartiennent l'un au fond,
   l'autre a l'objet, alors leur distance a 
   l'interface est de 0.5.
   
   on garde 2 chiffres de precision.
   la distance sera 
   la valeur calculee * MIN( voxel[i] ) / 100.0

   Pour initialiser, on suppose que la distance
   est estimee au centre de chaque voxel,
   les distances initiales se calculent 
   par rapport a la frontiere du voxel

*/

void _InitDistanceCoefficients( typeDistanceCoefficients *par,
				float *voxelSize /* float [3] */ 
				)
{
  double vx = voxelSize[0];
  double vy = voxelSize[1];
  double vz = voxelSize[2];
  double min;
  /* coordonnees du point le plus proche 
     sur le voxel
  */
  double x, y, z;
  int i,j,k;

  int initMaxIndex = 7;
  int initMiddleIndex = 3;

  int propMaxIndex = 5;
  int propMiddleIndex = 2;

  
  min = vx;
  if (min > vy) min = vy;
  if (min > vz) min = vz;

  par->multiplicativeCoefficient = min / 100.0;
  
  /* calcul des coefficients pour l'initialisation
     
     on calcule la distance a la surface du voxel
     (initMiddleIndex, initMiddleIndex, initMiddleIndex)
     
     pour cela, on calcule la distance du centre de chaque
     voxel a son plus proche voisin sur la surface du 
     voxel (initMiddleIndex, initMiddleIndex, initMiddleIndex).

     Le plus proche voisin (x,y,z) du voxel (i,j,k) a 
     des coordonnees qui peuvent s'exprimer independamment
     
  */

  for ( k=0; k<initMaxIndex; k++) 
  for ( j=0; j<initMaxIndex; j++)
  for ( i=0; i<initMaxIndex; i++) {
    
    if ( i < initMiddleIndex ) x = (double)initMiddleIndex - 0.5;
    else if ( i == initMiddleIndex ) x = (double)initMiddleIndex;
    else x = (double)initMiddleIndex + 0.5;

    if ( j < initMiddleIndex ) y = (double)initMiddleIndex - 0.5;
    else if ( j == initMiddleIndex ) y = (double)initMiddleIndex;
    else y = (double)initMiddleIndex + 0.5;

    if ( k < initMiddleIndex ) z = (double)initMiddleIndex - 0.5;
    else if ( k == initMiddleIndex ) z = (double)initMiddleIndex;
    else z = (double)initMiddleIndex + 0.5;

    par->init[k][j][i] = (int)( sqrt( (i-x)*(i-x)*vx*vx +
				     (j-y)*(j-y)*vy*vy +
				     (k-z)*(k-z)*vz*vz ) / 
			       par->multiplicativeCoefficient + 0.5 );
  }

  par->minInit = par->init[initMiddleIndex][initMiddleIndex][0];
  if ( par->minInit > par->init[initMiddleIndex][0][initMiddleIndex] )
    par->minInit = par->init[initMiddleIndex][0][initMiddleIndex];
  if ( par->minInit > par->init[0][initMiddleIndex][initMiddleIndex] )
    par->minInit = par->init[0][initMiddleIndex][initMiddleIndex];
  par->maxInit = par->init[0][0][0];


  /* calcul des coefficients pour la propagation de distance
       
     on calcule la distance au voxel central
     (initMiddleIndex, initMiddleIndex, initMiddleIndex)
     
  */

  for ( k=0; k<propMaxIndex; k++)
  for ( j=0; j<propMaxIndex; j++)
  for ( i=0; i<propMaxIndex; i++) {
    par->incr[k][j][i] = (int)( sqrt( (i-propMiddleIndex)*(i-propMiddleIndex)*vx*vx +
				      (j-propMiddleIndex)*(j-propMiddleIndex)*vy*vy +
				      (k-propMiddleIndex)*(k-propMiddleIndex)*vz*vz ) / 
				par->multiplicativeCoefficient + 0.5 );
  }

}













void _PrintDistanceCoefficients( typeDistanceCoefficients *par )
{
  int i,j,k;

  int initMaxIndex = 7;
  int initMiddleIndex = 3;

  int propMaxIndex = 5;
  int propMiddleIndex = 2;

  printf( "... coefficients d'initialisation : min=%d   max=%d\n", 
	  par->minInit, par->maxInit );

  for ( k=0; k<initMaxIndex; k++ ) {
    printf( "      plan %d\n",k-initMiddleIndex );
    for ( j=0; j<initMaxIndex; j++ ) {
      for ( i=0; i<initMaxIndex; i++ )
	printf("%5d ", par->init[k][j][i]);
      printf("   ---   ");
      for ( i=0; i<initMaxIndex; i++ )
	printf("%7f ", (double)(par->init[k][j][i])*par->multiplicativeCoefficient );
      printf("\n");
    }
    printf("\n");
  }

  printf( "... coefficients de propagation\n" );

  for ( k=0; k<propMaxIndex; k++ ) {
    printf( "      plan %d\n",k-propMiddleIndex );
    for ( j=0; j<propMaxIndex; j++ ) {
      for ( i=0; i<propMaxIndex; i++ )
	printf("%5d ", par->incr[k][j][i]);
      printf("   ---   ");
      for ( i=0; i<propMaxIndex; i++ )
	printf("%7f ", (double)(par->incr[k][j][i])*par->multiplicativeCoefficient );
      printf("\n");
    }
    printf("\n");
  }


}

























void _CombineTwoDistanceMaps2D( typeDistanceMap *theDist1,
				typeDistanceMap *theDist2 )
{
  typeDistanceMap *theMin, *theMax;
  typeDistanceCoefficients coeff;

  short int *minBuf, *maxBuf;
  int slice=theDist1->dim[0] * theDist1->dim[1];
  int z, i;
  int imin, imax;
  
  if ( typeComputation != _TWO_DIMENSIONAL_ )
    return;
  if ( theDist1->intensityMax == theDist2->intensityMax )
    return;

  _InitDistanceCoefficients( &coeff, theDist1->voxelSize );
  /* _PrintDistanceCoefficients( &coeff); */

  if ( theDist1->intensityMax > theDist2->intensityMax ) {
    theMax = theDist1;
    theMin = theDist2;
    imax = theDist1->intensityMax;
    imin = theDist2->intensityMax;
  } else {
    theMax = theDist2;
    theMin = theDist1;
    imax = theDist2->intensityMax;
    imin = theDist1->intensityMax;
  }
  minBuf = theMin->buf;
  maxBuf = theMax->buf;

  if ( imin == 0 ) {
    for (i=0;i<slice;i++)
      minBuf[i] = coeff.init[3][3][4];
    imin = 1;
  }


  for ( z=imin; z<imax; z++ )
  for (i=0;i<slice;i++)
    minBuf[z*slice+i] = minBuf[(z-1)*slice+i] + coeff.incr[2][2][3];
  
  theDist1->intensityMax = theDist2->intensityMax = imax;

}
