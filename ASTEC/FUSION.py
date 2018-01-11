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

import os
import os.path
import sys
sys.path.append('../CommunFunctions') #TO ADD /Functions/
import numpy as np
from ImageHandling import SpatialImage, imread, imsave
from crop import compute_volumes,croping
import math
from cpp_wrapping import reech, linear_registration, apply_trsf
from scipy import ndimage as nd
from time import time as time_estimation



# -- replace this directory to match your file layout --


def axis_rotation_matrix(axis, angle, min_space=None, max_space=None):
    """ Return the transformation matrix from the axis and angle necessary
    axis : axis of rotation ("X", "Y" or "Z")
    angle : angle of rotation (in degree)
    min_space : coordinates of the bottom point (usually (0, 0, 0))
    max_space : coordinates of the top point (usually im shape)
    """
    I = np.linalg.inv
    D = np.dot
    if axis not in ["X", "Y", "Z"]:
        raise Exception("Unknown axis : "+ str(axis))
    rads = math.radians(angle)
    s = math.sin(rads)
    c = math.cos(rads)

    centering = np.identity(4)
    if min_space is None and max_space is not None:
        min_space = np.array([0.,0.,0.])

    if max_space is not None:
        space_center = (max_space-min_space)/2.
        offset = -1.*space_center
        centering[:3,3] = offset

    rot = np.identity(4)
    if axis=="X":
        rot = np.array([ [1., 0., 0., 0.],
                            [0., c, -s,  0.],
                            [0., s,  c,  0.],
                            [0., 0., 0., 1.] ])
    elif axis=="Y":
        rot = np.array([ [c,   0., s,  0.],
                            [0.,  1., 0., 0.],
                            [-s,  0., c,  0.],
                            [0.,  0., 0., 1.] ])

    elif axis=="Z":
        rot = np.array([ [c, -s,  0., 0.],
                            [s,  c,  0., 0.],
                            [0., 0., 1., 0.],
                            [0., 0., 0., 1.] ])

    return D(I(centering), D(rot, centering))
     

def histogram(image, nbins=256):
    """Return histogram of image.
        
        Unlike `np.histogram`, this function returns the centers of bins and
        does not rebin integer arrays. For integer arrays, each integer value has
        its own bin, which improves speed and intensity-resolution.
        
        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.
        
        Returns
        -------
        hist : array
        The values of the histogram.
        bin_centers : array
        The values at the center of the bins.
        """
    
    # For integer types, histogramming with bincount is more efficient.
    if np.issubdtype(image.dtype, np.integer):
        offset = 0
        if np.min(image) < 0:
            offset = np.min(image)
        hist = np.bincount(image.ravel() - offset)
        bin_centers = np.arange(len(hist)) + offset
        
        # clip histogram to start with a non-zero bin
        idx = np.nonzero(hist)[0][0]
        return hist[idx:], bin_centers[idx:]
    else:
        hist, bin_edges = np.histogram(image.flat, nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
        return hist, bin_centers


def threshold_otsu(image, nbins=256):
    """Return threshold value based on Otsu's method.
        
        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.
        
        Returns
        -------
        threshold : float
        Threshold value.
        
        References
        ----------
        .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method
        
        Examples
        --------
        >>> from skimage.data import camera
        >>> image = camera()
        >>> thresh = threshold_otsu(image)
        >>> binary = image > thresh
        """
    hist, bin_centers = histogram(image, nbins)
    hist = hist.astype(float)
    
    # class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]
    
    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:])**2
    
    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]
    return threshold

def exp_func(x, length=500, speed=5):
    """ Decay function used to take into account the remotness to the camera
    x : value to compute
    length : lenght of the function
    speed : speed of the function
    """

    return .1+np.exp(-((np.float32(x)*speed/length)))

def build_mask(im, direction):
    """Return the mask on a given image from the decay function
    im : intensity image (SpatialImage)
    direction : if True the camera is in the side of the first slices in Z
    """
    th = threshold_otsu(im)
    im_th=np.zeros_like(im)
    im_th[im>th]=1
    if direction==False:
        im_th=im_th[:,:,-1::-1]
    im_th_sum=np.cumsum(im_th, axis=2)
    if direction==False:
        im_th_sum=im_th_sum[:,:,-1::-1]
    mask = exp_func(im_th_sum, np.max(im_th_sum))
    return mask


def fusion(images_input, image_output,temporary_path,ori, mirrors=False):
    """Compute the fusion of a set of raw data at a given time point
    images_input : list of raw images
    image_output : name of the fused image
    temporary_path : path where the temporary images will be written
    ori : orientation of the rotation between the consecutive raw images
    mirrors : if True a mirror was used to correct the X-Y transformation"""

    # References Images
    references_files=[image_input.replace('.inr','_ref.inr') for image_input in images_input]
    im_flos_t0=[imread(image_input) for image_input in images_input]
    im_ref_t0=im_flos_t0[0]; #The First angle is the starting reference point
    voxelsize = im_ref_t0.voxelsize;
    print voxelsize
    im_ref_t0.resolution = voxelsize
    imsave(references_files[0], im_ref_t0) #Save the first angle as a reference 
    im_flos_t0=im_flos_t0[1:len(im_flos_t0)];
    
    for i in range(len(im_flos_t0)):
        if i%2==0:
            if mirrors:
                im=SpatialImage(im_flos_t0[i])
            else:
                im=SpatialImage(im_flos_t0[i].copy())[-1::-1,:,:]
        else :
            im=SpatialImage(im_flos_t0[i])
        im.voxelsize=voxelsize
        imsave(references_files[i+1], im)
 

    #Rotation Matrix
    rotation_files=[image_input.replace('.inr','_rot.txt') for image_input in images_input]
    if ori=='left':
        a=270.
    else:
        a=90.
    for i, im in enumerate(im_flos_t0):
        if i==0:
            angle=0.
        else:
            angle=a
        print 'Angle used :'+str(angle)+' for ' + rotation_files[i+1]
        rot = axis_rotation_matrix(axis="Y", angle=angle, min_space=(0,0,0),
                                   max_space=np.multiply(im.shape[:3], im.resolution))
        np.savetxt(rotation_files[i+1],rot)


    registration_files=[image_input.replace('.inr','_reg.inr') for image_input in images_input];
    print 'Linear Reech in ' +registration_files[0]
    reech(references_files[0], registration_files[0], voxelsize[0])
    #imsave((registration_files[0]).replace('.inr','.tiff'),imread(registration_files[0])) #TEMPORARY 
    

    trsf_files=[image_input.replace('.inr','.trsf') for image_input in images_input];
    for i in range(1, len(images_input)):
        print 'Linear Registration in ' +registration_files[i]
        linear_registration(registration_files[0], references_files[i], rotation_files[i], registration_files[i], trsf_files[i])
        #imsave((registration_files[i]).replace('.inr','.tiff'),imread(registration_files[i])) #TEMPORARY
        
    #Mask    
    mask_files=[image_input.replace('.inr','_mask.inr') for image_input in images_input];
    mask=build_mask(im_ref_t0, True)
    mask.resolution=voxelsize
    temporary_mask=mask_files[0].replace('.inr','_temp.inr')
    imsave(temporary_mask, mask)
    #imsave(temporary_mask.replace('.inr','.tiff'), mask) #TEMPORARY
    reech(temporary_mask, mask_files[0], voxelsize[0])
    full_mask=imread(mask_files[0])
    #imsave((mask_files[0]).replace('.inr','.tiff'),full_mask) #TEMPORARY
    for i in range(1, len(images_input)):  
        print ' Calcul Mask in '+ mask_files[i]
        im=imread(references_files[i])
        if i%2==1:
            direction=False
        else:
            direction=True
        mask=build_mask(im, direction)
        mask.resolution=voxelsize
        temporary_mask=mask_files[i].replace('.inr','_temp.inr')
        imsave(temporary_mask, mask)
        del mask
        apply_trsf(temporary_mask, trsf_files[i], mask_files[i], registration_files[0]) 
        full_mask+=imread(mask_files[i])
        #imsave((mask_files[i]).replace('.inr','.tiff'),imread(mask_files[i])) #TEMPORARY 

    final=np.zeros_like(full_mask)
    final=final.astype(np.uint16)
    for i in range(len(images_input)):
        final+=(imread(registration_files[i])*imread(mask_files[i]))/full_mask

    im_th=np.zeros((final.shape[0], final.shape[1], 1), dtype=np.uint16)
    im_max=np.max(final, axis=2)
    th=np.mean(im_max)
    im_th[im_max>th]=1
    comp_co = nd.label(im_th)[0]
    volumes = compute_volumes(comp_co)
    volumes.pop(0)        
    label = volumes.keys()[np.argmax(volumes.values())]
    bb = nd.find_objects(comp_co)[label-1]
    bb2 = (slice(max(bb[0].start - 40, 0), min(final.shape[0], bb[0].stop + 40), None), slice(max(bb[1].start - 40, 0), min(final.shape[1], bb[1].stop + 40), None), slice(0, final.shape[2]))
    imsave(image_output, SpatialImage(final[bb2]))
    #imsave(image_output.replace('.inr','.tiff'), SpatialImage(final[bb2])) #TEMPORARY

    
def fusion_process(time_angles_files,output_file,temporary_path,
                   ori, resolution, target_resolution, delay,
                   ext_im, mirrors=False):
    """Compute the fusions of a given time-series of raw data
    time_angles_files : ??
    output_file : ??
    image_output : name of the fused image
    temporary_path : path where the temporary images will be written
    ori : orientation of the rotation between the consecutive raw images
    resolution : resolution of the raw images in $\mu m$
    target_resolution : resolution of the final fused image in $\mu m$
    delay : time to add in the name of the fused files if the movie have been splited
    ext_im : extension of the image ('.inr', '.tiff', '.tif', '.h5')
    mirrors : if True a mirror was used to correct the X-Y transformation"""
    
    tozip=False #Zip files ?
    if ext_im.lower()=='.zip':
        tozip=True
    
    start_process=time_estimation()
    output_path=output_file[:output_file.rfind('/')]
    os.system('mkdir -p ' + output_path) #Create output folder
    
    os.system("rm -rf "+temporary_path) #Delte the temporary folder if previous created
    os.system('mkdir -p ' + temporary_path) #Create temporary folder
    downsize=target_resolution/resolution[0] #Downsampling
     

    ### Pre-treatment (Unzip if necessary and conver in inr format )
    angle_paths=[temporary_path+'ANGLE_'+str(a)+'/' for a in range(len(time_angles_files))]
    [os.system('mkdir -p ' + angle_path) for angle_path in angle_paths] #Create Angle Temporary Directory
        
    inr_files=[]
    for time_angle_file,angle_path in zip(time_angles_files,angle_paths):
        if tozip: #Unzip if necessary
            image_file=time_angle_file[time_angle_file.rfind('/')+1:]
            print 'Unzip '+time_angle_file
            f=os.listdir(angle_path)
            os.system('unzip ' + time_angle_file  + ' -d ' + angle_path)
            [image_file for image_file in os.listdir(angle_path) if image_file not in f]
            ext_im=image_file[image_file.find('.'):]
            image_file=angle_path+image_file    
        else:
        	image_file=time_angle_file
        print 'Convert in inr '+image_file
        im=imread(image_file)
        im.resolution=resolution
        image_file=image_file[image_file.rfind('/')+1:]
        inr_file=angle_path+image_file.replace(ext_im, '.inr').replace('\\', '')
        imsave(inr_file, im)
        inr_files.append(inr_file)

    
    ### Croping process
    cropped_files=[inr_file.replace('.inr','_cropp.inr') for inr_file in inr_files]
    for inr_file,cropped_file in zip(inr_files,cropped_files):
        print 'Crop ' + inr_file
        croping(inr_file, cropped_file,downsize)
  

    ### Fusion process 
    print ' Fusion  on '+str(cropped_files)
    fusion(cropped_files,output_file, temporary_path , ori, mirrors)
    print 'Fusion done in '+output_file

    return time_estimation()-start_process
    
        
#List angles images files in a folder 
def read_raw_data(datapath):
    print ' Read rawdata path ' + datapath
    Extensions=[]; FileNames=[]; ErrorReturn=-1,0,0,'',''
    for f in os.listdir(datapath):
        if not f.startswith('.'):
            extension=f[f.find('.'):]
            if extension not in Extensions:
                Extensions.append(f[f.rfind('.'):])
            FileNames.append(f[:f.rfind('.')])
 
    if len(FileNames)==0:
        print 'Error this folder '+datapath+' contains no files'
        return ErrorReturn
    if len(Extensions)!=1:
        print ' Error this folder '+datapath+' contains several different extensions '+str(Extensions)
        return ErrorReturn
    if extension.lower()!='.zip' and extension.lower()!='.tif' and extension.lower()!='.zip' and extension.lower()!='.hdr' and extension.lower()!='.f5':
        print ' Error this folder '+datapath+' does not have valid angles images (.tif,.tiff.hdr.f5,.zip) '+str(Extensions)
        return ErrorReturn
    lens=len(FileNames[0])
    for i in range(len(FileNames)):
        if len(FileNames[i])!=lens:
            print ' Error file names in this folder '+datapath+' does not have the same length'
            return ErrorReturn
        
    #Only one files
    if len(FileNames)==1:
        #Look for 3 consecutives digit starting from the end
        filename=FileNames[0]; idx=len(filename)-3; startName=filename[0:idx]; endName='';
        if len(filename)<3:
            print ' Error the file in this folder '+datapath+' does not have 3 digits to define time steps'
            return ErrorReturn
        while len(startName)>=3:
            try : 
                timestep=int(filename[len(startName):len(startName)+3])
                return 1,timestep,timestep,extension,datapath+'/'+startName+'$TIME'+filename[len(startName)+3:len(filename)]
            except ValueError:
                startName=startName[0:len(startName)-1];
        print ' Error the file in this folder '+datapath+' does not have 3 digits to define time steps'
        return ErrorReturn
        
    FileNames.sort()
    #Look for the first common name
    startName=''; common=1; idx=0
    while common==1:
        common=1; checkLetter=FileNames[0][idx]
        for i in range(len(FileNames)):
            if FileNames[i][idx]!=checkLetter:
                common=0
        if common==1:
            startName+=checkLetter
        idx+=1
    
    #Look for the last common name
    endName=''; common=1; idx=lens-1
    while common==1:
        common=1; checkLetter=FileNames[0][idx]
        for i in range(len(FileNames)):
            if FileNames[i][idx]!=checkLetter:
                common=0
        if common==1:
            endName=checkLetter+endName
        idx-=1
        
    #Look for the definition of times (should be on 3 digits)
    digittime=lens-len(startName)-len(endName)
    if digittime>3:
        print ' Error file names in this folder '+datapath+' does not have 3 digits to define time steps'
        return ErrorReturn
    if digittime<3:
        if len(startName)+digittime<3:
            print ' Error file names in this folder '+datapath+' does not have 3 digits to define time steps'
            return ErrorReturn
        startName=startName[0:len(startName)-3+digittime]
    
    begin=None; idx=0
    for i in range(len(FileNames)):
        try : 
            timestep=int(FileNames[i][len(startName):lens-len(endName)])
            if begin is None:
                begin=timestep
            else:
                if timestep!=begin+idx:
                    print ' Error file names in this folder '+datapath+' does not have consecutive 3 digits to define time steps'
                    return ErrorReturn
            idx+=1
        except ValueError:
            print ' Error file names in this folder '+datapath+' does not always have 3 digits to define time steps'
            return ErrorReturn
    end=timestep
    return 1,begin,end,extension,datapath+'/'+startName+'$TIME'+endName
