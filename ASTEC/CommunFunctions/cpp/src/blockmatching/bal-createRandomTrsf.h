/*************************************************************************
 * applyTrsf.c -
 *
 * $Id$
 *
 * Copyright (c) INRIA 2012, all rights reserved
 *
 * AUTHOR:
 * Gregoire Malandain (gregoire.malandain@inria.fr)
 * 
 * CREATION DATE: 
 * Mon Nov 19 17:45:00 CET 2012
 *
 *
 * ADDITIONS, CHANGES
 *
 *
 *
 *
 */

#ifndef BAL_CREATERANDOMTRSF_H
#define BAL_CREATERANDOMTRSF_H

#include <bal-transformation-tools.h>

#ifdef __cplusplus
extern "C" {
#endif

int createRandomTrsf(
        char *restrsf_name,
        char *template_name,
        bal_doublePoint fixedpoint,
        enumTypeTransfo transformation_type,
        int print);

#ifdef __cplusplus
}
#endif

#endif //BAL_CREATERANDOMTRSF_H
