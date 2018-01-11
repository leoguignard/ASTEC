/*
 * mt_mergeSegmentations.c -
 *
 * $Id: mt_mergeSegmentations.c,v 1.0 2014/04/02 11:55:34 gael Exp $
 *
 *
 * AUTHOR:
 * Gael Michelin (gael.michelin@inria.fr)
 *
 * CREATION DATE:
 * 2014/04/02
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#include <mt_mergeSegmentations.h>

//static int _verbose_ = 0;


int MT_BackgroundMask( vt_image imageIn, 
					   vt_image *imageMask, int value)
{
  char *proc = "MT_Compute3DMultiScale";
  int dim;
  bufferType type;
  unsigned char *theBuf;
  unsigned char *inU8;
  unsigned short int *inU16;
  int i;
  
  theBuf=(unsigned char *)imageMask->buf;
  dim=imageMask->dim.x*imageMask->dim.y*imageMask->dim.z;
  
  // Mask computation
   
  type=imageIn.type;
  switch (type) {
  case UCHAR:
	fprintf(stdout, "UCHAR\n");
    inU8=(unsigned char *)imageIn.buf;
	for (i=0; i<dim; i++)
	  theBuf[i]=((int)inU8[i] == value) ? (unsigned char)255 : (unsigned char)0;
	break;
  case USHORT:
  case SSHORT:
	fprintf(stdout, "USHORT\n");
    inU16=(unsigned short int *)imageIn.buf;
	for (i=0; i<dim; i++)
	  theBuf[i]=((int)inU16[i] == value) ? (unsigned char)255 : (unsigned char)0;
	fprintf(stdout, "i=%d\n", i);
	break;
  default:
	fprintf(stderr, "%s: input image type not handled yet\n", proc);
	return(0);
  }
  return(1);
}


int MT_ExtOnBackground(vt_image imageIn,
							  vt_image imageExt, 
							  vt_image imageMask,
							  int dilatation,			//argument de dilatation pr chacun des labels de imageExt
							  int *labelsExt,			//liste des labels de imageExt qui intersectent la region non nulle d'imageMask
							  int *volumes,				//volumes des labels de labelsExt
							  int *volumesIntersections,//volume des intersections entre les labelsExt et le background
							  int *nbVoisins,			//nombre de labels de imageIn en contact avec les labelsExt
							  int *n					//nombre de labelsExt comptes... (size des tableaux)
							  )
{
  //TODO
  return(1);
}

