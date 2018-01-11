/*************************************************************************
 * invTrsf.c -
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



#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <bal-blockmatching.h>
#include <bal-cropImage.h>
#include <bal-invTrsf.h>




int main(int argc, char *argv[])
{

    char* workingDir = "";

    char* imgIn = "";
    char* imgFirst = "";
    char* imgLast = "";
    char* ellipseFormatedLsit = "";
    int numHead=0;
    int xOffset=0;
    int yOffset=0;

    // 1) Setup Dirs
    //  - datas = x4_dil200_10microns_100Hz_images_analyze/

    //  - resultats = TRACKING
    //  - centers = TRACKING/HEADnum/CENTERS
    //  - transfo = TRACKING/HEADnum/TRANSFO
    //  - tmp-transfo = TRACKING/HEADnum/TMP-TRANSFO
    //  - logs = TRACKING/HEADnum/LOGS
    //  - images = TRACKING/HEADnum/IMAGES
    //  - tmp-images = TRACKING/HEADnum/TMP-IMAGES

//    // 2) Dimensions des imagettes
//    bal_integerPoint dim;
//    dim.x = 2 * xOffset;
//    dim.y = 2 * yOffset;

//    // 3) foreach images
//    // - recupere xpos et ypos tronqués
//    // - x et y tronqués
//    bal_integerPoint origin;
//    origin.x = xPos-xOffset;
//    origin.y = yPos-yOffset;

//    bal_integerPoint slice;
//    slice.x = 0;
//    slice.y = 0;

//    // Image I
//    if(!cropImage(
//            "4_dil200_10microns_100Hz_images_firstImage", // HDR
//            "imageNumerotee dans rep de sortie", // -flo.HDR
//            "path vers le fichier de transformations", // -flo.TRSF
//            NULL,
//            NULL,
//            origin,
//            dim,
//            NULL
//            ))
//    {
//            fprint("Whoahowohwohwohoa");
//    }

//    // Image I+1
//    if(!cropImage(
//            "4_dil200_10microns_100Hz_images_firstImage+1", // HDR
//            "imageNumerotee dans rep de sortie", // -ref.HDR
//            "path vers le fichier de transformations", // -ref.TRSF
//            NULL,
//            NULL,
//            origin,
//            dim,
//            slice
//            ))
//    {
//            fprint("Whoahowohwohwohoa");
//    }

//    if(!invTrsf(
//            "path vers le fichier de transformations", // -ref.TRSF
//            "path vers le fichier de transformations", // -ref-inv.TRSF
//            NULL,
//            NULL,
//            NULL
//            ))
//    {
//        fprint("Whoahowohwohwohoa");
//    }


    exit( 0 );
}




