path_to_bins = 'ASTEC/CommunFunctions/cpp/build/bin/'

path_filters = path_to_bins + 'recfilters'
path_linearfilters = path_to_bins + 'linearFilter'
path_reech3d = path_to_bins + 'reech3d'
path_apply_trsf = path_to_bins + "applyTrsf"
path_block = path_to_bins + "blockmatching"
path_morpho = path_to_bins + "morpho"
path_regional_max = path_to_bins + 'regionalmax'
path_connexe = path_to_bins + 'connexe'
path_watershed = path_to_bins + 'watershed'
path_gradient_norm = path_to_bins + 'norme_gradient'



import os
from ImageHandling import imread, imsave, SpatialImage

def recfilter(path_input, path_output='tmp.inr', filter_value=2, lazy=False):
    ''' Perform a gaussian filtering on an intensity image
    path_input : path to the image to filter
    path_output : path to the temporary output image
    filter_value : sigma of the gaussian filter
    lazy : do not return the output image if True
    '''
    os.system(path_filters + ' ' + path_input +\
              ' ' + path_output +\
              ' -cont 10 -sigma ' + str(filter_value) +\
              ' -x 0 -y 0 -z 0 -o 2')
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out  


def linearfilter(path_input, path_output='tmp.inr', filter_value=2, rad_min=1, lazy=False):
    ''' Perform a gaussian filtering on an intensity image
    path_input : path to the image to filter
    path_output : path to the temporary output image
    filter_value : sigma of the gaussian filter
    rad_min : TO REMOVE, NOT USED
    lazy : do not return the output image if True
    '''
    os.system(path_linearfilters + ' ' + path_input +\
              ' ' + path_output +\
              ' -cont 10 -sigma ' + str(filter_value) +\
              ' -x 0 -y 0 -z 0 -o 2')
    if not lazy:
        out = imread(path_output)
        os.system('rm ' + path_output)
        return out        

def regionalmax(path_input, path_output, h_min):
    ''' Perform the h-minima operation on a given image
    path_input : path to the input image
    path_output : path to the output image
    h_min : h-minima parameter value
    '''
    temp_out=path_output.replace('.inr','_out_regionalmax.inr')
    os.system(path_regional_max + ' ' + path_input +\
              ' -diff ' + path_output +' '+temp_out+\
              ' -h ' + str(h_min) +\
              ' -inv')
    
    os.system('rm ' + temp_out)

def connexe(path_input, path_output, high_th):
    ''' Perform the connected componant operation
    path_input : path to the input image
    path_output : path to the output image
    high_th : high threshold for an hysteresis filtering on the componants
    '''
    os.system(path_connexe + ' ' + path_input +\
              ' ' + path_output +\
              ' -lt 1 -ht ' + str(high_th) +\
              ' -labels -o 2')


def watershed(path_seeds, path_int, path_output=None, lazy=True):
    ''' Perform the watershed operation
    path_seeds : path to the seeds image
    path_int : path to the intensity image
    path_output : path to the output image
    lazy : do not return the output image if True
    '''    
    if type(path_seeds)!=str:
        imsave("seeds.inr", path_seeds)
        path_seeds = "seeds.inr"
    if type(path_int)!=str:
        imsave("intensity.inr", path_int)
        path_int = "intensity.inr"

    if path_output is None:
        lazy = False
        path_output = 'seg.inr'
 
    os.system(path_watershed + ' ' + path_seeds +\
              ' ' + path_int +\
              ' ' + path_output \
              )
    if not lazy:
        out=imread(path_output)
        os.system('rm seeds.inr intensity.inr seg.inr')
        return out

def reech3d(im_path, output_shape):  
    ''' Perform a resampling operation
    im_path : path to the image to resample
    output_shape : desired output shape
    '''
    tmp_file=im_path.replace('.inr','_temp.inr')
    os.system(path_reech3d + ' '+im_path+'  '+tmp_file+
              ' -x ' + str(int(output_shape[0])) +
              ' -y ' + str(int(output_shape[1])) +
              ' -z ' + str(int(output_shape[2])))
    out = imread(tmp_file)
    os.system('rm '+tmp_file)
    return  out 

def reech(path_flo, path_output, voxelsize):
    ''' Perform a resampling operation
    path_flo : path to the image to resample
    path_output : path to the output image
    voxelsize : size of the voxels in $\mu m**3$ (x, y, z)
    '''
    os.system(path_reech3d +
              " " + path_flo + 
              " " + path_output +
              " -linear"
              " -iso " + str(voxelsize))

def non_linear_registration(image_file_flo,image_file_ref, affine_image, affine_trsf,vectorfield_image,vectorfield_trsf):
    ''' Compute the non-linear transformation that register the floating image onto the reference image
    image_file_flo : path to the floating image
    image_file_ref : path to the reference image
    affine_image : path to the floating image after affine registration
    affine_trsf : path to the affine transformation
    vectorfield_image : path to the floating image after affine o non-linear registration
    vectorfield_trsf : path to the non-linear registration (affine o non-linear)
    '''
    os.system(path_block +
              " -ref " + image_file_ref+
              " -flo " + image_file_flo+
              " -res " + affine_image +
              " -res-trsf " +affine_trsf+
              " -trsf-type affine" +
              " -estimator wlts" +
              " -py-gf" +
              " -pyramid-highest-level 5" +
              " -pyramid-lowest-level 3" +
              " -lts-fraction 0.55")
    os.system(path_block +
              " -ref " + image_file_ref+
              " -flo " + image_file_flo+
              " -res " + vectorfield_image+
              " -res-trsf " + vectorfield_trsf+
              " -init-trsf " +affine_trsf+
              " -trsf-type vectorfield" +
              " -estimator wlts" +
              " -py-gf" +
              " -pyramid-highest-level 5" +
              " -pyramid-lowest-level 3" +
              " -elastic-sigma 2.0 2.0 2.0" +
              " -fluid-sigma 2.0 2.0 2.0")
    

def linear_registration(path_ref, path_flo, path_trsf, path_output, path_output_trsf):
    ''' Compute the linear transformation that register the floating image onto the reference image
    path_flo : path to the floating image
    path_ref : path to the reference image
    path_output : path to the floating image after affine registration
    path_trsf : path to the initial registration
    path_output_trsf : path to the affine transformation
    '''
    os.system(path_block +
              " -ref " + path_ref +
              " -flo " + path_flo +
              " -res " + path_output +
              " -res-trsf " + path_output_trsf +
              " -init-trsf " + path_trsf +
              " -trsf-type affine" +
              " -estimator wlts" +
              " -pyramid-highest-level 6" +
              " -pyramid-lowest-level 3" +
              " -lts-fraction 0.55")


def apply_trsf(path_flo, path_trsf, path_output="tmp_seeds.inr", 
               template=None, nearest=True, lazy=True):
    ''' Apply a transformation to a given image
    path_flo : path to the floating image
    path_trsf : path to the transformation
    path_output : path to the output image
    template : path to the template image
    nearest : do not interpolate (take the nearest value) if True, to use when applying on label images
    '''
    command_line = path_apply_trsf + " " + path_flo + " " + path_output 
    command_line += " -trsf " + path_trsf
    if not template is None:
        command_line += " -template " + template
    if not nearest is None:
        command_line += " -nearest"
    os.system(command_line)
    if not lazy:
        out=imread(path_output)
        if path_output=='tmp_seeds.inr':
            os.system('rm -f tmp_seeds.inr')
        return out


def find_local_minima(path_out, path_ref, h_min, mask=None, sigma=2):
    ''' Find local minima in an intensity image
    path_out : path to the output seeds image
    path_ref : path to the reference intensity image
    h_min : value of the h-minima operator value
    mask : mask on the intensity image
    sigma : value of the gaussian filter in voxels
    '''
    from os import path
    path_mask_out=path_out.replace('.inr','_mask_'+str(h_min)+'.inr')
    tmp_min=path_out.replace('.inr','_local_minima_out.inr')
    tmp_filt=path_out.replace('.inr','_local_minima_filter'+str(sigma)+'.inr') 
    if not path.exists(tmp_filt) and mask==None:
        recfilter(path_ref, tmp_filt, filter_value=sigma, lazy=True)
    if mask==None:
        os.system(path_regional_max + ' ' + tmp_filt + ' ' +\
                  ' -diff ' + path_mask_out + ' ' +\
                  tmp_min + ' ' +\
                  '-h ' + str(h_min) + ' ' +\
                  '-inv')
    else:
        os.system(path_regional_max + ' ' + mask + ' ' +
                  '-diff ' + path_mask_out + ' ' +
                  tmp_min + ' ' +
                  '-h ' + str(h_min))
    os.system(path_connexe + ' ' + path_mask_out + ' ' +
              path_out + ' ' +
              '-sb 1 -sh ' + str(h_min) +
              ' -labels -o 2')
    try:
        im=imread(path_out.replace('\\', ''))
    except:
        im=None
    os.system('rm -f '+tmp_filt+' '+tmp_min);
    return im, path_mask_out

def morpho(image_input,image_output,paramstre):
	''' Morphological operation
	'''
	os.system(path_morpho+' '+image_input+' '+' '+image_output+' '+ paramstre)

	
	
	

def outer_detection(im_ref_tmp, radius, seg_ref_tmp):
    ''' Compute the detection of the outer of the embryo
    im_ref_tmp : intensity image for the outer detection (SpatialImage)
    radius : radius of the grey closing to perform
    seg_ref_tmp : segmented reference image (SpatialImage)
    '''
    from copy import deepcopy
    if radius!='0':
        imsave("tmp_bounds.inr", im_ref_tmp)
        os.system(path_morpho + " tmp_bounds.inr closed.inr -clo -R " + radius)
        im=imread("closed.inr")
    else:
        im=deepcopy(im_ref_tmp)
    imax=np.max(im)
    h=np.histogram(im, range=(0, imax), bins=imax)
    cumhist=np.cumsum(h[0][::-1]),h[1][::-1]
    vol=np.sum(seg_ref_tmp!=1)#*1.10
    low=np.max(cumhist[1][cumhist[0]>vol])
    im_th=np.zeros_like(im)
    #im=imread("closed.inr")
    im_th[im>=low]=1 # Cytoplasm
    if radius!='0':
        imsave("tmp.inr", SpatialImage(im_th))
        os.system((path_morpho + " tmp.inr closing.inr -clo -R " + radius))
    else:
        imsave("closing.inr", SpatialImage(im_th))
    os.system((path_morpho + " closing.inr erode.inr -ero -R 5"))
    imE=imread("closing.inr")
    imE=nd.binary_fill_holes(imE)
    mask=np.uint8(imE)
    bounds=nd.binary_dilation(mask, structure=nd.generate_binary_structure(3, 1))-mask
    im_refB=im_ref_tmp.copy()
    im_refB[bounds.astype(np.bool)]=np.max(im_ref_tmp)
    imsave('tmp.inr', SpatialImage(im_refB))
    os.system(path_filters + " tmp.inr out_bounds.inr -x 0 -y 0 -z 0 -sigma 1 -o 2")
    return imread('out_bounds.inr'), bounds.astype(np.bool)

def gradient_norm(image_input,gradient_output):
    ''' Perform the gradient norm of an intensity image
    im_input : input image (SpatialImage)
    path_input : path to the input image
    path_output : path to the output image
    '''
    os.system(path_gradient_norm + ' ' + image_input + ' ' + gradient_output + ' -sigma 1')