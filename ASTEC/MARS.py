import os, sys
#sys.path.append('CommunFunctions') 
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import imread, imsave,SpatialImage
import numpy as np
from scipy import ndimage as nd
from cpp_wrapping import recfilter, regionalmax, connexe, watershed


def mars_segmentation(image_input, segmentation_output, sigma, h_min, sigma_ws, th=0):
    """ Perform the watershed segmentation of a given intensity image
    image_input : path to the input image
    segmentation_output : path to the final segmetented image
    sigma : value of the gaussian filter on the intensity image used for the h-minima operation (in voxels)
    h_min : value of the h-minima operation parameter (in bit)
    cleaning_needed : True if a removing of the cells touching the border of the image is needed
    sigma_ws : value of the gaussian filter on the intensity image used for the watershed (in voxel)
    th:Threshold for output cleaning
    """
    os.system('mkdir -p '+segmentation_output[:segmentation_output.rfind('/')])
    print 'Process Segmentation of '+image_input
    if image_input.split('.')[-1]=='tif' or image_input.split('.')[-1]=='tiff':
        print ' Convert in inr ' + image_input
        image_input_inr=''
        for n in image_input.split('.')[:-1]:
            image_input_inr+=n
        image_input_inr+='.inr'
        imsave(image_input_inr, imread(image_input))
        image_input=image_input_inr

    #### Definition of paths to the different outputs ####
    path_gSigma=segmentation_output.replace('.inr','g' + str(sigma) + '.inr') # Path to the smoothed image for the local minima detection
    path_g_5=segmentation_output.replace('.inr','g'+str(sigma_ws)+'.inr') # Path to the smoothed image for the watershed
    path_rm=segmentation_output.replace('.inr','rm_s'+str(sigma)+'_h'+str(h_min)+'.inr') # Path to the regionalmax image
    path_cc=segmentation_output.replace('.inr','c_s'+str(sigma)+'_h'+str(h_min)+'.inr') # Path to the seeds image
    
    print 'Filter with sigma='+str(sigma)+' in ' + path_gSigma
    recfilter(image_input, path_gSigma, sigma, lazy=True)
    
    print 'Filter with sigma='+str(sigma_ws)+' in ' + path_g_5
    recfilter(image_input, path_g_5, sigma_ws, lazy=True)
    
    print 'Find Local minnima with h_min='+str(h_min)+' in ' + path_rm
    regionalmax(path_gSigma, path_rm, h_min)
    
    print 'Find connex composant in ' + path_cc
    connexe(path_rm, path_cc, h_min)
    
    print 'Process watershed segmentation in ' + segmentation_output
    watershed(path_cc, path_g_5, segmentation_output)
    
    #Delete Temporary files
    os.system('rm ' + path_gSigma + ' ' + path_g_5 + ' ' + path_rm + ' ' + path_cc)
    
    print 'Segmentation done'

