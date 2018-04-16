#!/usr/bin/python2.7


#Mars first time point segmentation
from definitions import *


from ImageHandling import imread, imsave, SpatialImage
from MARS import mars_segmentation 
from lineage import timeNamed,timesNamed

import os
if not os.path.isdir(segmented_Path):
    os.mkdir(segmented_Path)  

os.system("cp -f "+astec_Path+"definitions.py "+segmented_Path )
#os.system("cp -f "+astec_Path+"2-mars.py "+segmented_Path )




fused_file =timeNamed(fused_files,begin)  #First Fused time step
segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment
reconstruct_file=None

if mars_method == 1:
    print "Starting with fused files..."
if mars_method == 2:
    print "Starting with GACE files..."
    reconstruct_file =timeNamed(reconstruct_files,begin)  


#Apply Automatic MARS Segmentation
mars_segmentation(fused_file, mars_file, mars_sigma1, mars_h_min, mars_sigma2, 
				  method=mars_method, reconstructed_image=reconstruct_file,
				  sigma_membrane=mars_sigma_membrane, sensitivity=mars_sensitivity, manual=mars_manual, manual_sigma=mars_manual_sigma, 
				  hard_thresholding=mars_hard_thresholding, hard_threshold=mars_hard_threshold, sigma_TV=mars_sigma_TV, sigma_LF=mars_sigma_LF, sample=mars_sample)

