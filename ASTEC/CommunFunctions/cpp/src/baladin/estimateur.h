#ifndef ESTIMATEUR_H
#define ESTIMATEUR_H

#ifdef __cplusplus
extern "C" {
#endif

#include <behavior.h>

#include <baladin.h>
#include <matrix.h>



typedef struct {
  typeField x;     /* position    */
  typeField y;
  typeField u;     /* déplacement */
  typeField v;
  typeField rho;   /* similarite  */
} typeScalarWeighted2DDisplacement;

typedef struct {
  typeField x;     /* position    */
  typeField y;
  typeField z;
  typeField u;     /* déplacement */
  typeField v;
  typeField w;
  typeField rho;   /* similarite  */
} typeScalarWeighted3DDisplacement;

typedef enum {
  _UnknownField_,
  _ScalarWeighted2DDisplacement_,
  _ScalarWeighted3DDisplacement_
} enumFieldType;


typedef struct {
  int n_pairs;       /* nombre de vecteurs */
  int n_allocated_pairs;
  enumFieldType type;
  void *pairs;
} FIELD;






double LS_Trsf_Estimation ( _MATRIX *T, FIELD *field, 
			    enumTypeTransfo transfo,
			    PARAM *param );
double LSSW_Trsf_Estimation ( _MATRIX *T, FIELD *field, 
			      enumTypeTransfo transfo,
			      PARAM *param );
double Trsf_Estimation ( _MATRIX *T, FIELD *field, 
			 enumTypeTransfo transfo, 
			 enumTypeEstimator estimator,
			 PARAM *param );
double Trsf_Trimmed_Estimation ( _MATRIX *T, FIELD * field, 
				 enumTypeTransfo transfo, 
				 PARAM *param, 
				 int do_some_writing );


double Estimate_Transformation( _MATRIX *T, FIELD * field, 
				PARAM *param, 
				int do_some_writing );

/* FIELD management */

int Allocate_Field ( FIELD * field, enumFieldType type, int npoints );

void Free_Field ( FIELD * field );

/* MISC */

void CreateFileDef(FIELD *field, 
		   char *nom_image, char *nom_champ );

void PrintField( FILE *f, FIELD *field );

#ifdef __cplusplus
}
#endif

#endif
