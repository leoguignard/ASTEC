###########################################################################
###########################################################################
## Copyright (C) 2018  Guignard Leo <guingardl__@__janelia.hhmi.org>     ##
##                                                                       ##
## This program is free software: you can redistribute it and/or modify  ##
## it under the terms of the GNU General Public License as published by  ##
## the Free Software Foundation, either version 3 of the License, or     ##
## (at your option) any later version.                                   ##
##                                                                       ##
## This program is distributed in the hope that it will be useful,       ##
## but WITHOUT ANY WARRANTY; without even the implied warranty of        ##
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         ##
## GNU General Public License for more details.                          ##
##                                                                       ##
## You should have received a copy of the GNU General Public License     ##
## along with this program.  If not, see <http://www.gnu.org/licenses/>. ##
###########################################################################
###########################################################################

import os, sys
sys.path.append('CommunFunctions') 
from ImageHandling import SpatialImage, imread, imsave
from scipy import ndimage as nd
import numpy as np
from cpp_wrapping import reech3d

def compute_volumes(im):
    labels = np.unique(im)
    volume = nd.sum(np.ones_like(im), im, index=np.int16(labels))
    return dict(zip(labels, volume))


def croping(image_input, impage_output,downsize):
    ''' Automatically crop and resample an images
    image_input : path to the input image
    image_output : path to the output image
    downsize : voxel size of the resampled image (in \mu m)
    '''
    shape_begin=imread(image_input).shape #Image Size
    image_main = reech3d(image_input, (shape_begin[0]/np.float(downsize), shape_begin[1]/np.float(downsize), shape_begin[2])) #Downsampling
    vxsize=image_main.resolution
    im_maxed = image_main.max(axis=2)
    thr = np.mean(im_maxed)
    im_th = np.zeros((im_maxed.shape[0], im_maxed.shape[1], 1), dtype=np.uint8)
    im_th[im_maxed > thr] = 1
    comp_co = nd.label(im_th)[0]
    volumes = compute_volumes(comp_co)
    volumes.pop(0)
    label = volumes.keys()[np.argmax(volumes.values())]
    bb = nd.find_objects(comp_co)[label-1]
    bb2 = (slice(max(bb[0].start - 40, 0), min(image_main.shape[0], bb[0].stop + 40), None), slice(max(bb[1].start - 40, 0), min(image_main.shape[1], bb[1].stop + 40), None), slice(0, image_main.shape[2]))
    out=SpatialImage(image_main[bb2])
    out.voxelsize=(float("{0:.1f}".format(image_main.resolution[0])), float("{0:.1f}".format(image_main.resolution[1])), float("{0:.1f}".format(image_main.resolution[2])))
    imsave(impage_output, out.astype(np.uint16))
