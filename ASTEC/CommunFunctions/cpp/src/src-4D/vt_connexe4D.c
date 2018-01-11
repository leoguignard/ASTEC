
#include <vt_connexe4D.h>

static int _VERBOSE_ = 1;







#define VT_F2I( F ) ( (F) >= 0.0 ? ((int)((F)+0.5)) : ((int)((F)-0.5)) )




static u16 _low_value_ = (u16)100;
static u16 _hig_value_ = (u16)200;

static int VT_LocalConnexe4D( vt_image4D *resIm,
			      int connexite );





int VT_Hysteresis4D( vt_image4D *theIm, 
		     vt_image4D *resIm,
		     float seuilBas,
		     float seuilHaut )
{
  char *proc = "VT_Hysteresis4D";
  vt_image4D auxIm, *outputIm;
  typeBoolean auxImIsAllocated = False;
  u16 ****resBuf;
  register int t, z, y, x;
  int nbcc, connexite = 2;
  

  
  if ( resIm->type != USHORT ) {
    
    if ( VT_AllocAndInitImage4D( &auxIm, NULL,
			       theIm->dim.x, theIm->dim.y, theIm->dim.z,
			       theIm->dimt, USHORT ) != 1 ) {
      VT_Error( "unable to allocate auxiliary image", proc );
      return( -1 );
    }
    auxImIsAllocated = True;
    outputIm = &auxIm;

  } else {
    outputIm = resIm;
  }
  resBuf = (u16****)outputIm->array;


  

  /* on transforme les entrees */
  switch ( theIm->type ) {
  case UCHAR :
    {
      u8**** theBuf = (u8****)theIm->array;
      int ish = VT_F2I( seuilHaut );
      int isb = VT_F2I( seuilBas );
      for ( t=0; t<theIm->dimt; t++ )
      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	if ( theBuf[t][z][y][x] >= ish )
	  resBuf[t][z][y][x] = _hig_value_;
	else { 
	  if ( theBuf[t][z][y][x] >= isb )
	    resBuf[t][z][y][x] = _low_value_;
	  else
	    resBuf[t][z][y][x] = 0;
	}
      }
    }
    break;
  case USHORT :
    {
      u16**** theBuf = (u16****)theIm->array;
      int ish = VT_F2I( seuilHaut );
      int isb = VT_F2I( seuilBas );
      for ( t=0; t<theIm->dimt; t++ )
      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	if ( theBuf[t][z][y][x] >= ish )
	  resBuf[t][z][y][x] = _hig_value_;
	else { 
	  if ( theBuf[t][z][y][x] >= isb )
	    resBuf[t][z][y][x] = _low_value_;
	  else
	    resBuf[t][z][y][x] = 0;
	}
      }
    }
    break;
  case FLOAT :
    {
      r32**** theBuf = (r32****)theIm->array;
      for ( t=0; t<theIm->dimt; t++ )
      for ( z=0; z<theIm->dim.z; z++ )
      for ( y=0; y<theIm->dim.y; y++ )
      for ( x=0; x<theIm->dim.x; x++ ) {
	if ( theBuf[t][z][y][x] >= seuilHaut )
	  resBuf[t][z][y][x] = _hig_value_;
	else { 
	  if ( theBuf[t][z][y][x] >= seuilBas )
	    resBuf[t][z][y][x] = _low_value_;
	  else
	    resBuf[t][z][y][x] = 0;
	}
      }
    }
    break;
  default :
    VT_Error( "unable to deal with such input type", proc );
    return( -1 );
  }



  if ( (nbcc = VT_LocalConnexe4D( outputIm,
			  connexite )) == -1 ) {
    VT_Error( "connected component extraction failed", proc );
    if ( auxImIsAllocated == True )
      VT_FreeImage4D( &auxIm );
    return( -1 );
   
  }
  
  if ( _VERBOSE_ ) {
    fprintf( stderr, "%s: found %3d %d-connected 4D components\n", 
	     proc, nbcc, connexite );
  }


  /*--- relabeling ---*/
  for ( t=0; t<theIm->dimt; t++ )
  for ( z=0; z<theIm->dim.z; z++ )
  for ( y=0; y<theIm->dim.y; y++ )
  for ( x=0; x<theIm->dim.x; x++ ) {
    if ( resBuf[t][z][y][x] > 0 )
      resBuf[t][z][y][x] = 65535;
  }


  


  if ( auxImIsAllocated == True ) {
    (void)VT_CopyImage4D( &auxIm, resIm );
    VT_FreeImage4D( &auxIm );
  }

  return( 1 );

}















#define _EQUIVALENCE_ARRAY_SIZE_ 65536

typedef struct _my_cc{
    int label;
    int size;
    int keep;
} _my_cc;

typedef struct _my_vois{
  int x;
  int y;
  int z;
  int t;
} _my_vois;



static int VT_LocalConnexe4D( vt_image4D *resIm,
			      int connexite )
{
  char *proc = "VT_LocalConnexe4D";
  _my_cc *cc = (_my_cc *)NULL;
  u16 ****resBuf;
  register int t, z, y, x, i, j, k;
  int nbVoisins, s;
  _my_vois v[40]; /* ( 3^4 - 1 ) / 2 */

  int nb_neighbors, label_neighbors[40];
  int valid_labels = 0, used_labels = 0;
  register int current_label;
  typeBoolean inside;

  

  if ( resIm->type != USHORT ) {
    VT_Error( "unable to deal with such input type", proc );
    return( -1 );
  }
  resBuf = (u16****)resIm->array;


    /*
   * Equivalence table initialization
   * 
   * We allocate the largest equivalence table which may be needed
   * for a labeling from 1 to 65535 (image of type unsigned short).
   * 
   * Dynamic allocation is possible if needed.
   */
  
  cc = (_my_cc *)VT_Malloc( (unsigned int)((int)_EQUIVALENCE_ARRAY_SIZE_ * sizeof(_my_cc)) );
  if ( cc == (_my_cc *)NULL ) {
    VT_Error( "allocation failed for equivalence array", proc );
    return( -1 );
  }
  cc[0].label = 0;
  




  /*--- initialisation des voisins ---*/
  nbVoisins = 0;
  connexite = 2; /* i <-> connexite par variete de codimension i */
  for ( t=-1; t<=0; t++ )
  for ( z=-1; z<=1; z++ )
  for ( y=-1; y<=1; y++ )
  for ( x=-1; x<=1; x++ ) {
    if ( (t == -1)
	 || ( (t == 0) && (z == -1) )
	 || ( (t == 0) && (z == 0) && (y == -1) )
	 || ( (t == 0) && (z == 0) && (y == 0) && (x == -1) ) ) {
      s = 0;
      if ( (x == 1) || (x == -1) ) s++;
      if ( (y == 1) || (y == -1) ) s++;
      if ( (z == 1) || (z == -1) ) s++;
      if ( (t == 1) || (t == -1) ) s++;
      if ( (s > 0) && (s <= connexite) ) {
	v[nbVoisins].x = x;
	v[nbVoisins].y = y;
	v[nbVoisins].z = z;
	v[nbVoisins].t = t;
	nbVoisins ++;
      }
    }
  }



		    
  /*--- traitement ---*/
  for ( t=0; t<resIm->dimt; t++ )
  for ( z=0; z<resIm->dim.z; z++ )
  for ( y=0; y<resIm->dim.y; y++ )
  for ( x=0; x<resIm->dim.x; x++ ) {
    
    if ( resBuf[t][z][y][x] < _low_value_ ) continue;
    
    inside = True;
    if ( (t == 0) ||
	 (z == 0) || (z == resIm->dim.z - 1) ||
	 (y == 0) || (y == resIm->dim.y - 1) ||
	 (x == 0) || (x == resIm->dim.x - 1) )
      inside = False;

    /*--- on cherche les voisins 
          deja parcourus et
	  deja numerotes         ---*/
    nb_neighbors = 0;
    
    if ( inside == True ) {

      for ( i = 0; i < nbVoisins; i ++ )
	if ( resBuf[t+v[i].t][z+v[i].z][y+v[i].y][x+v[i].x] > 0 )
	  label_neighbors[ nb_neighbors++ ] = cc[ (int)resBuf[t+v[i].t][z+v[i].z][y+v[i].y][x+v[i].x] ].label;

    } else {

      for ( i = 0; i < nbVoisins; i ++ ) {
	if ( (t+v[i].t >= 0) &&
	     (z+v[i].z >= 0) && (z+v[i].z < resIm->dim.z ) &&
	     (y+v[i].y >= 0) && (y+v[i].y < resIm->dim.y ) &&
	     (x+v[i].x >= 0) && (x+v[i].x < resIm->dim.x ) )
	  if ( resBuf[t+v[i].t][z+v[i].z][y+v[i].y][x+v[i].x] > 0 )
	    label_neighbors[ nb_neighbors++ ] = cc[ (int)resBuf[t+v[i].t][z+v[i].z][y+v[i].y][x+v[i].x] ].label;
      }

    }
    
    /*--- gestion des equivalences ---*/
    if ( nb_neighbors > 0 ) {
      /*
       * If there are some neighbors,
       * the equivalent label for all of them is the minimum
       * of all labels.
       */
      current_label = label_neighbors[0];
      for ( i = 1; i < nb_neighbors; i ++ )
	if ( label_neighbors[i] < current_label ) current_label = label_neighbors[i];
      /*
       * we have to report the equivalences inside the equivalence table
       */
      for ( i = 0; i < nb_neighbors; i ++ ) {
	if ( label_neighbors[i] != current_label ) {
	  k = label_neighbors[i];
	  /*
	   * si la classe n'a pas deja ete traitee
	   */
	  if ( cc[ k ].label != current_label ) {
	    cc[ current_label ].size += cc[ k ].size;;   cc[ k ].size = 0;
	    if ( cc[ k ].keep == 1 ) cc[ current_label ].keep = 1;
	    cc[ k ].label = current_label;
	    for ( j = k+1; j <= used_labels; j ++ )
	      if ( cc[ j ].label == k ) cc[ j ].label = current_label;
	  }
	}
      }
    } else {
      /*
       * No neighbors :
       * we use a new label for this point.
       * This is the only case of error in the whole procedure.
       * This could be a little bit improved as follows:
       * 1. we re-label the potential equivalence classes (without conditions of
       *    validity, i.e. without testing cc[ i ].keep and cc[ i ].size) to get
       *    the effective number of equivalence classes at this time.
       * 2. using this re-labeling, we change the image labels up to the current
       *    point.
       * 3. we re-build the equivalence table to have only one entry per equivalence
       *    classe.
       */
      current_label = ++used_labels;
      if ( used_labels > (_EQUIVALENCE_ARRAY_SIZE_ - 1) ) {
	VT_Free( (void**)(&cc) );
	VT_Error( "Too much used labels for connected components computation", proc );
	return( -1 );
      }
      cc[ current_label ].label = current_label;
      cc[ current_label ].size  = 0;
      cc[ current_label ].keep  = 0;
    }
    
    if ( resBuf[t][z][y][x] >= _hig_value_ ) cc[ current_label ].keep  = 1;
    resBuf[t][z][y][x] = (u16)current_label;
    cc[ current_label ].size ++;
    
  }


  /*--- pas de composantes ---*/
  if ( used_labels == 0 ) return( 0 );


  /*
   * counting and re-labeling all the valid connected components
   * this is not optimal, but more elegant because it may be used
   * in several further cases.
   *
   * The representative label of a equivalence's classe
   * (i.e. a connected component) is characterized by
   * ( cc[ i ].label == i ).
   * We have to check if this connected component is valid i.e.
   * - if it contains some points with (value >= _hig_value_),
   *   => (cc[ i ].keep == 1)
   * - if its size is large enough
   *   => (cc[ i ].size >= par->min_size) )
   * if yes, we increment the number of valid connected components
   * (++valid_labels) and give this value as a new label for the
   * connected component. 
   * if no, we give 0 as a new label and as a new size.
   *
   * If a label is not representative of its equivalence's classe
   * (i.e. cc[ i ].label < i), we give to it the new value
   * of its representative label (which may be 0 if the connected
   * component is not valid). Recall that the representative label
   * of an equivalence's classe is always parsed before the other
   * labels of the same equivalence's classe. Id. for the size.
   *
   * To keep recognizing the valid classes, the representative labels
   * have the flag keep equal to 1, while the others not.
   */

  for ( i = 1; i <= used_labels; i++ ) {
    if ( cc[ i ].label == i ) {
      if ( (cc[ i ].keep == 1) && (cc[ i ].size >= 1) )
	cc[ i ].label = ++valid_labels; 
      else {
	cc[ i ].label = cc[ i ].size = 0; 
	cc[ i ].keep  = 0;
      }
    }
    else {
      j = cc[ i ].label; 
      cc[ i ].label = cc[ j ].label;
      cc[ i ].size = cc[ j ].size;
      cc[ i ].keep  = 0;
    }
  }

  
  /*--- relabeling ---*/
  for ( t=0; t<resIm->dimt; t++ )
  for ( z=0; z<resIm->dim.z; z++ )
  for ( y=0; y<resIm->dim.y; y++ )
  for ( x=0; x<resIm->dim.x; x++ ) {
    resBuf[t][z][y][x] = (u16)cc[ (int)resBuf[t][z][y][x] ].label;
  }




  VT_Free( (void**)(&cc) );
  return( valid_labels );
}

