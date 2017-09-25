#Mars first time point segmentation
from definitions import *


from ImageHandling import imread, imsave, SpatialImage
from MARS import mars_segmentation 
from lineage import timeNamed,timesNamed

import os
if not os.path.isdir(segmented_Path):
    os.mkdir(segmented_Path)  


#Parameters for MARS segmentation
sigma1 = 0.6 / target_resolution   #sigma 1 (0.6um)
sigma2 = 0.15 / target_resolution #sigma 2 (0.15um) 
h_min = 4   # H min initialisation to ease correction

fused_file =timeNamed(fused_files,begin)  #First Fused time step
segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment

mars_segmentation(fused_file, segmentation_file, sigma1, h_min, sigma2) #Apply Automatic MARS Segmentation

imsave(segmentation_file.replace('.inr','_mars.tiff'),imread(segmentation_file)) #To have it in TIFF
