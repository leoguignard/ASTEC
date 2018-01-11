

#include <vt_image4D.h>

static int NAMESTOBEALLOCATED = 20;





void VT_Name4D( vt_name4D *name )
{
  name->nbAllocatedNames = 0;
  name->nbNames = 0;
  name->names = (char **)NULL;
}






int VT_AllocateName4D( vt_name4D *name )
{
  int i, l;
  char **names, *n;
  
  l = ( name->nbAllocatedNames + NAMESTOBEALLOCATED ) * STRINGLENGTH * sizeof(char)
    + ( name->nbAllocatedNames + NAMESTOBEALLOCATED ) * sizeof(char *);
  
  names = (char **)VT_Malloc( (unsigned int)l );
  if ( names == (char **)NULL ) {
    VT_Error( "Allocation failed", "VT_AllocateName4D" );
    return( -1 );
  }
  
  n = (char*)(names + name->nbAllocatedNames + NAMESTOBEALLOCATED);
  for ( i=0; i < name->nbAllocatedNames + NAMESTOBEALLOCATED; i++, n += STRINGLENGTH ) {
    names[i] = n;
  }

  for (i=0; i<name->nbAllocatedNames; i++) {
    VT_Strncpy( names[i], name->names[i], STRINGLENGTH );
  }

  VT_Free( (void**)&(name->names) );

  name->names = names;
  name->nbAllocatedNames += NAMESTOBEALLOCATED;


  return( 1 );
}







void VT_FreeName4D( vt_name4D *name )
{
  VT_Free( (void**)&(name->names) );
  name-> nbAllocatedNames = 0;
  name->nbNames = 0;
  name->names = (char **)NULL;
}









int VT_ReadNameInr4D( FILE *f,
		  vt_name4D *name4D )
{
  char *proc = "VT_ReadNameInr4D";
  char name[STRINGLENGTH];
  int status, nb=0;
  

/*
# 4D Inrimage File version 1.000000 
CONTENT 
*/

  do {
    status = fscanf( f, "%s", name );
  } while ( (status == 1) && (strcmp( name, "CONTENT" ) != 0) );
/*  
  Number of 4D INRIMAGE = 24
*/
  do {
    status = fscanf( f, "%s", name );
  } while ( (status == 1) && (strcmp( name, "=" ) != 0) );
  if ( fscanf( f, "%d", &nb ) != 1 ) {
    fprintf( stderr, "%s: error when reading\n", proc );
    return( 0 );
  }
  if ( nb > 0 ) NAMESTOBEALLOCATED = nb;
  if ( VT_AllocateName4D( name4D ) != 1 ) {
    VT_Error( "unable to reallocate name 4D", proc );
    return( 0 );
  }
 

/*
  ASCII
  Big endian
END
*/
  do {
    status = fscanf( f, "%s", name );
  } while ( (status == 1) && (strcmp( name, "END" ) != 0) );
/*
4D INRIMAGE kikinis.inr4D
  FileName
END
......blabla.inr
*/

  do {
    do {
      status = fscanf( f, "%s", name );
    } while ( (status == 1) && (strcmp( name, "END" ) != 0) );
    if ( (status = fscanf( f, "%s", name )) == 1 ) {
      VT_Strncpy( name4D->names[name4D->nbNames], name, STRINGLENGTH );
      name4D->nbNames ++;      
    }
  } while ( status == 1 );



  if ( _VT_DEBUG_  != 0 ) {
    int i;
    fprintf( stderr, " read %d names in %s\n", name4D->nbNames, proc );
    for (i=0; i<name4D->nbNames; i++ ) {
      fprintf( stderr, " #%3d: %s\n", i, name4D->names[i] );
    }
  }
  
  return( 1 );
}















int VT_ReadName4D( FILE *f,
		  vt_name4D *name4D )
{
  char *proc = "VT_ReadName4D";
  char name[STRINGLENGTH];
  
  
  while( fscanf( f, "%s", name ) == 1 ) {
    if ( name4D->nbAllocatedNames == name4D->nbNames ) {
      if ( VT_AllocateName4D( name4D ) != 1 ) {
	VT_Error( "unable to reallocate name 4D", proc );
	return( 0 );
      }
    }
    VT_Strncpy( name4D->names[name4D->nbNames], name, STRINGLENGTH );
    name4D->nbNames ++;
  }

  if ( _VT_DEBUG_  != 0 ) {
    int i;
    fprintf( stderr, " read %d names in %s\n", name4D->nbNames, proc );
    for (i=0; i<name4D->nbNames; i++ ) {
      fprintf( stderr, " #%3d: %s\n", i, name4D->names[i] );
    }
  }
  
  return( 1 );
}







int VT_CreateName4D( vt_name4D *name4D,
		     char *base,
		     char *suffixe,
		     int n )
{
  char *proc = "VT_CreateName4D";
  char *defaultBase = "tmp";
  char *mybase = base;
  int i;
  
  if ( base == (char*)NULL ) mybase = defaultBase;

  for ( i=0; i<n; i++ ) {
    if ( name4D->nbAllocatedNames == i ) {
      if ( VT_AllocateName4D( name4D ) != 1 ) {
	VT_Error( "unable to reallocate name 4D", proc );
      }
    }
    if ( suffixe != (char*)NULL ) {
      sprintf( name4D->names[name4D->nbNames],"%s.%03d%s", mybase,i,suffixe );
    } else {
      sprintf( name4D->names[name4D->nbNames],"%s.%03d", mybase,i );
    }
    name4D->nbNames ++;
  }
  return( 1 );
}



void VT_WriteNameInr4D( vt_name4D *name4D, char *nom )
{
  FILE *f, *fopen();
  int i;

  if ( (nom != (char*)NULL) && (nom[0] != '>') ) {
    f = fopen( nom, "w" );
  } else {
    f = stdout;
  }

  fprintf( f, "# 4D Inrimage File version 1.000000\n" );
  fprintf( f, "CONTENT\n" );
  fprintf( f, "  Number of 4D INRIMAGE = %d\n", name4D->nbNames );
  fprintf( f, "  ASCII\n" );
  fprintf( f, "  Big endian\n" );
  fprintf( f, "END\n" );

  for (i=0; i<name4D->nbNames; i++ ) {
    fprintf( f, "4D INRIMAGE kikinis.inr4D\n" );
    fprintf( f, "  FileName\n" );
    fprintf( f, "END\n" );
    fprintf( f, "%s\n", name4D->names[i] );
  }

  if ( (nom != (char*)NULL) && (nom[0] != '>') ) {
    fclose( f );
  }
}










void VT_WriteName4D( vt_name4D *name4D, char *nom )
{
  FILE *f, *fopen();
  int i;

  if ( (nom != (char*)NULL) && (nom[0] != '>') ) {
    f = fopen( nom, "w" );
  } else {
    f = stdout;
  }

  for (i=0; i<name4D->nbNames; i++ ) {
    fprintf( f, "%s\n", name4D->names[i] );
  }

  if ( (nom != (char*)NULL) && (nom[0] != '>') ) {
    fclose( f );
  }
}






int VT_AllocAndInitImage4D( vt_image4D *image,
			    vt_name4D *name4D,
			    int dimx,
			    int dimy,
			    int dimz,
			    int dimt,
			    int type )
{
  char *proc="VT_AllocAndInitImage4D";
  int i, j;
  int volumeArray = 0;

  if ( dimt <= 0 ) {
    VT_Error( "unable to allocate 4D image with dimt <= 0", proc );
    return( 0 );
  }
  

  image->images = (vt_image*)VT_Malloc( (unsigned int)(dimt * sizeof(vt_image)) );
  if ( image->images == (vt_image*)NULL ) {
    VT_Error( "images array allocation failed", proc );
    return( 0 );
  }



  for( i=0; i<dimt; i++ ) {

    if ( (name4D == (vt_name4D*)NULL) || (i >= name4D->nbNames) ) {
      VT_InitImage( &(image->images[i]), "",
		    dimx, dimy, dimz, type );
    } else {
      VT_InitImage( &(image->images[i]), name4D->names[i],
		    dimx, dimy, dimz, type );
    } 
    
    
    if ( VT_AllocImage( &(image->images[i]) ) != 1 ) {
      fprintf( stderr, " %s: allocation of image #%d failed\n.", proc, i );
      for ( j=0; j < i; j++ ) {
	VT_FreeImage( &(image->images[j]) );
      }
      VT_Free( (void**)&(image->images) );
      return( 0 );
    }

  }


  volumeArray = dimt;
  switch( type ) {
  case UCHAR :
    volumeArray *= sizeof( u8***);
    break;
  case USHORT :
    volumeArray *= sizeof( u16***);
    break;
  case FLOAT :
    volumeArray *= sizeof( r32***);
    break;
  case DOUBLE :
    volumeArray *= sizeof( r64***);
    break;
  default :
    VT_Error( "such 4D image type not handled yet", proc );
    for ( j=0; j < dimt; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }



  image->array = VT_Malloc( (unsigned int)(volumeArray) );
  if ( image->array == (void****)NULL ) {
    VT_Error( "allocation of temporal array failed", proc );
    for ( j=0; j < dimt; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }

  
  switch( type ) {
  case UCHAR :
    {
      u8 ****a = (u8****)image->array;
      for ( j=0; j < dimt; j++ )
	a[j] = (u8***)image->images[j].array;
    }
    break;
  case USHORT :
    {
      u16 ****a = (u16****)image->array;
      for ( j=0; j < dimt; j++ )
	a[j] = (u16***)image->images[j].array;
    }
    break;
  case FLOAT :
    {
      r32 ****a = (r32****)image->array;
      for ( j=0; j < dimt; j++ )
	a[j] = (r32***)image->images[j].array;
    }
    break;
  case DOUBLE :
    {
      r64 ****a = (r64****)image->array;
      for ( j=0; j < dimt; j++ )
	a[j] = (r64***)image->images[j].array;
    }
    break;
  default :
    VT_Error( "such 4D image type not handled yet", proc );
    for ( j=0; j < dimt; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }




  image->dim.x = dimx;
  image->dim.y = dimy;
  image->dim.z = dimz;
  image->dimt = dimt;
  image->type = type;


  return( 1 );

}









int VT_AllocAndReadImage4D( vt_image4D *image,
		     vt_name4D *name4D )
{
  char *proc="VT_AllocAndReadImage4D";
  int i, j;
  int volumeArray = 0;

  if ( name4D->nbNames <= 0 ) {
    VT_Error( "unable to allocate 4D image with dimt <= 0", proc );
    return( 0 );
  }
  

  image->images = (vt_image*)VT_Malloc( (unsigned int)(name4D->nbNames * sizeof(vt_image)) );
  if ( image->images == (vt_image*)NULL ) {
    VT_Error( "images array allocation failed", proc );
    return( 0 );
  }



  for( i=0; i<name4D->nbNames; i++ ) {

    if ( VT_ReadInrimage( &(image->images[i]), name4D->names[i] ) != 1 ) {
      fprintf( stderr, " %s: reading of image #%d failed\n.", proc, i );
      for ( j=0; j < i; j++ ) {
	VT_FreeImage( &(image->images[j]) );
      }
      VT_Free( (void**)&(image->images) );
      return( 0 );
    }

  }
  


  for( i=1; i<name4D->nbNames; i++ ) {
    if ( (image->images[i].dim.x != image->images[0].dim.x) ||
	 (image->images[i].dim.y != image->images[0].dim.y) ||
	 (image->images[i].dim.z != image->images[0].dim.z) ||
	 (image->images[i].dim.v != image->images[0].dim.v) ||
	 (image->images[0].dim.v != 1) ||
	 (image->images[i].type != image->images[0].type) ) {
      VT_Error( "Images have different dimensions or type", proc );
      for ( j=0; j < name4D->nbNames; j++ ) {
	VT_FreeImage( &(image->images[j]) );
      }
      VT_Free( (void**)&(image->images) );
      return( 0 );
    }
  }


  volumeArray = name4D->nbNames;
  switch( image->images[0].type ) {
  case UCHAR :
    volumeArray *= sizeof( u8***);
    break;
  case USHORT :
    volumeArray *= sizeof( u16***);
    break;
  case FLOAT :
    volumeArray *= sizeof( r32***);
    break;
  case DOUBLE :
    volumeArray *= sizeof( r64***);
    break;
  default :
    VT_Error( "such 4D image type not handled yet", proc );
    for ( j=0; j < name4D->nbNames; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }



  image->array = VT_Malloc( (unsigned int)(volumeArray) );
  if ( image->array == (void****)NULL ) {
    VT_Error( "allocation of temporal array failed", proc );
    for ( j=0; j < name4D->nbNames; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }

  
  switch( image->images[0].type ) {
  case UCHAR :
    {
      u8 ****a = (u8****)image->array;
      for ( j=0; j < name4D->nbNames; j++ )
	a[j] = (u8***)image->images[j].array;
    }
    break;
  case USHORT :
    {
      u16 ****a = (u16****)image->array;
      for ( j=0; j < name4D->nbNames; j++ )
	a[j] = (u16***)image->images[j].array;
    }
    break;
  case FLOAT :
    {
      r32 ****a = (r32****)image->array;
      for ( j=0; j < name4D->nbNames; j++ )
	a[j] = (r32***)image->images[j].array;
    }
    break;
  case DOUBLE :
    {
      r64 ****a = (r64****)image->array;
      for ( j=0; j < name4D->nbNames; j++ )
	a[j] = (r64***)image->images[j].array;
    }
    break;
  default :
    VT_Error( "such 4D image type not handled yet", proc );
    for ( j=0; j < name4D->nbNames; j++ ) VT_FreeImage( &(image->images[j]) );
    VT_Free( (void**)&(image->images) );
    return( 0 );
  }





  image->dim.x = image->images[0].dim.x;
  image->dim.y = image->images[0].dim.y;
  image->dim.z = image->images[0].dim.z;
  image->dimt = name4D->nbNames;
  image->type = image->images[0].type;


  return( 1 );

}





void VT_FreeImage4D( vt_image4D *image )
{
  int i;

  for (i=0; i<image->dimt; i++ ) 
    VT_FreeImage( &(image->images[i]) );
  VT_Free( (void**)&(image->images) );
  image->dimt = 0;
  image->dim.x = 0;
  image->dim.y = 0;
  image->dim.z = 0;
  
}






int VT_CopyImage4D( vt_image4D *theIm,
		    vt_image4D *resIm )
{
  char *proc = "VT_CopyImage4D";
  int t;

  for( t=0; t<theIm->dimt; t++ ) {

    if ( VT_CopyImage( &(theIm->images[t]),
		       &(resIm->images[t]) ) != 1 ) {
      fprintf( stderr, "%s: copy failed on image #%d.\n", proc, t );
      return( -1 );
    }
  }
  return( 1 );
}










int VT_WriteImage4D( vt_image4D *image ) 
{
  char *proc= "VT_WriteImage4D";
  int i;

  for (i=0; i<image->dimt; i++ ) {
    if ( VT_WriteInrimage( &(image->images[i]) ) != 1 ) {
      fprintf( stderr, "%s: writing failed on image #%d.\n", proc, i );
       return( -1 );
    }
  }
  return( 1 );

}
