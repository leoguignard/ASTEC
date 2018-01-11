/* */


#ifndef _YelTrtroisD_h_
#define _YelTrtroisD_h_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  double x;
  double y;
  double z;
} typeVertex;

typedef struct {
  int vertex[3];
  int neighbor[3];
} typeTriangle;

typedef struct {
  int nVertex;
  typeVertex *vertex;
  int nTriangle;
  typeTriangle *triangle;
} typeTriangulation;




extern void initTriangulation( typeTriangulation *t );
extern void freeTriangulation( typeTriangulation *t );
extern int readTriangulation( typeTriangulation *t,
			      char *filename );
extern int printTriangulation( typeTriangulation *triangulation,
			       char *filename );

extern int computeDistancesWithIdenticalVertex( typeTriangulation *t1, 
					 typeTriangulation *t2,
					 double *d );

extern int computeDistancesWithClosestVertex( typeTriangulation *t1, 
					      typeTriangulation *t2,
					      double *d );

extern void printReport( double *d, int nd, char *t1, char *t2, char *desc );

#ifdef __cplusplus
}
#endif


#endif
