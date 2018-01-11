#ifndef VECTEUR_H
#define VECTEUR_H

#ifdef __cplusplus
extern "C" {
#endif


#include <baladin.h>
#include <balimage.h>
#include <estimateur.h>

extern void set_at_field_computation_print_field_( int i );

void CalculChampVecteurs (FIELD *field, 
			  bal_image *inrimage_flo, BLOCS *blocs_flo,
			  bal_image *inrimage_ref, BLOCS *blocs_ref,
			  PARAM *param );

int CalculAttributsBlocs (bal_image *inrimage, BLOCS *blocs,
			  int seuil_bas, int seuil_haut, float seuil_pourcent, 
			  PARAM *param );

#ifdef __cplusplus
}
#endif

#endif
