#ifndef MISC_H
#define MISC_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>
#include <stdio.h>

#include <behavior.h>

#include <baladin.h>

void PrintListOfBlocks( FILE *f, BLOCS* b );
void PrintListOfValidBlocks( FILE *f, BLOCS* b );
void PrintPointersOfBlocks( FILE *f, BLOCS* b );
void PrintPointersOfValidBlocks( FILE *f, BLOCS* b );
void PrintParametersOfBlocks( FILE *f, BLOCS* b );
void PrintParametersForBlocks( FILE *f, PARAM *p );

#ifdef __cplusplus
}
#endif

#endif
