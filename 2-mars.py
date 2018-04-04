#Mars first time point segmentation
from definitions import *


from ImageHandling import imread, imsave, SpatialImage
from MARS import mars_segmentation 
from lineage import timeNamed,timesNamed

import os
if not os.path.isdir(segmented_Path):
    os.mkdir(segmented_Path)  

os.system("cp -f "+astec_Path+"definitions.py "+segmented_Path )
os.system("cp -f "+astec_Path+"2-mars.py "+segmented_Path )




fused_file =timeNamed(fused_files,begin)  #First Fused time step
segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment
reconstruct_file=None

if Mars_method == 1:
    print "Starting with fused files..."
if Mars_method == 2:
    print "Starting with GACE files..."
    reconstruct_file =timeNamed(reconstruct_files,begin)  


#Apply Automatic MARS Segmentation
mars_segmentation(fused_file, mars_file, sigma1_mars, h_min_mars, sigma2_mars, 
				  method=Mars_method, reconstructed_image=reconstruct_file,
				  sigma_membrane=sigma_membrane, sensitivity=sensitivity, manual=manual, manual_sigma=manual_sigma, 
				  hard_thresholding=hard_thresholding, hard_threshold=hard_threshold, sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample)

# MERGE Sophia - Montpellier 
# COMMENTARY (GAEL) : POURQUOI ON A BESOIN DE SAUVER EN TIFF ICI ? C'EST LE INR QUI EST LU EN ENTREE DE 3-manualcorrection.py.....
#imsave(mars_file.replace('.inr','_mars.tiff'),imread(mars_file)) #To have it in TIFF
