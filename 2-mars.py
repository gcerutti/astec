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



### Parameters for MARS segmentation
sigma1 = 0.6 # sigma 1 (0.6um) in real coordinates
sigma2 = 0.15 # sigma 2 (0.15um) in real coordinates
h_min = 4   # H min initialisation to ease correction


### Gace Parameters (if Mars_method is set to 2 or 3):

# membrane_renforcement
sigma_membrane=0.9 # membrane enhancement parameter (in real units, a priori 0.9 um is a good choice for data like Patrick/Ralph/Aquila)

# anisotropicHist /!\ critical step
sensitivity=0.99 # membrane binarization parameter, /!\ if failure, one should enter in "manual" mode of the function anisotropicHist via activation of 'manual' option

manual=False     # By default, this parameter is set to False. If failure, (meaning that thresholds are very bad, meaning that the binarized image is very bad),
				 # set this parameter to True and relaunch the computation on the test image. If the method fails again, "play" with the value of manual_sigma... and good luck.
manual_sigma=15  # Axial histograms fitting initialization parameter for the computation of membrane image binarization axial thresholds (this parameter is used iif manual = True).
				 # One may need to test different values of manual_sigma. We suggest to test values between 5 and 25 in case of initial failure. Good luck.

hard_thresholding=False  # If the previous membrane threshold method failed, one can force the thresholding with a "hard" threshold applied on the whole image. To do so, this option must be set to True.
hard_threshold=1.0       # If hard_thresholding = True, the enhanced membranes image is thresholded using this parameter (value 1 seems to be ok for time-point t001 of Aquila embryo for example).

# TVmembrane
sigma_TV=3.6     # parameter which defines the voting scale for membrane structures propagation by tensor voting method (real coordinates). 
				 # This parameter shoud be set between 3 um (little cells) and 4.5 um(big gaps in the binarized membrane image)
sigma_LF=0.9     # Smoothing parameter for reconstructed image (in real coordinates). It seems that the default value = 0.9 um is ok for classic use.
sample=0.2       # Parameter for tensor voting computation speed optimisation (do not touch if not bewared)





fused_file =timeNamed(fused_files,begin)  #First Fused time step
segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment
reconstruct_file=None

if Mars_method == 1:
    print "Starting with fused files..."
if Mars_method == 2:
    print "Starting with GACE files..."
    reconstruct_file =timeNamed(reconstruct_files,begin)  


#Apply Automatic MARS Segmentation
mars_segmentation(fused_file, segmentation_file, sigma1, h_min, sigma2, 
				  method=Mars_method, reconstructed_image=reconstruct_file,
				  sigma_membrane=sigma_membrane, sensitivity=sensitivity, manual=manual, manual_sigma=manual_sigma, 
				  hard_thresholding=hard_thresholding, hard_threshold=hard_threshold, sigma_TV=sigma_TV, sigma_LF=sigma_LF, sample=sample)

# MERGE Sophia - Montpellier 
# COMMENTARY (GAEL) : POURQUOI ON A BESOIN DE SAUVER EN TIFF ICI ? C'EST LE INR QUI EST LU EN ENTREE DE 3-manualcorrection.py.....
imsave(segmentation_file.replace('.inr','_mars.tiff'),imread(segmentation_file)) #To have it in TIFF
