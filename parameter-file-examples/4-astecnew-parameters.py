############################################################
##
## useful parameters for segmentation propagation
##
## to be used with 4-astec-new.py
##
############################################################



# PATH_EMBRYO = ''

## ##### explanation #####
##
## path to the embryo data
## e.g. '/media/DATA/171107-Karine-St8'
## if not present, the current directory is used
##



# EN = ''

## ##### explanation #####
##
## Embryo Name
## CRBM naming format is YYMMDD-SaintOfTheDays-Stage
## eg: '171107-Karine-St8'
## (automatically extracted from PATH_EMBRYO if not provided)
##



begin = 0
end = 0
delta = 1
raw_delay = 0

## ##### explanation #####
##
## time points to be processed
##
## begin: first time point
## end: last time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
##
## For the fusion step, not giving the 'begin' and 'end' values
## causes the fusion of all found data in raw data directories
## according these 4 directories are different
##
## When testing or tuning parameters, it is advised to set 'begin' and 'end'
## at the same value 
##



# EXP_FUSE = ''

## ##### explanation #####
##
## fusion directory name
##
## This fusion directory contains the images to be segmented
## default value is
## EXP_FUSE = 'RELEASE'
##



# EXP_SEG = ''

## ##### explanation #####
##
## Segmentation results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
## default value is
## EXP_SEG = 'RELEASE'
##



# result_image_suffix = 'mha'
# default_image_suffix = 'mha'

## ##### explanation #####
##
## defines the output image format
##
## result_image_suffix is for all output images
## default_image_suffix is for all images (including auxiliary ones)
## default are 'inr'
##
## 'mha' should be prefered (readable by fiji)
##



######################################################################
##
## reconstruction parameters
##
######################################################################



# astec_intensity_transformation = 'Normalization_to_u8'
# astec_intensity_enhancement = None


## ##### explanation #####
##
## ASTEC method (nothing but a cell-based seeded watershed) may be applied on
## a transformed input image or on the original image (eg the result of
## the fusion step). This transformed imaged is made of a combination
## (by the maximum operator) of the transformed intensity image and a
## membrane-enhanced image.
## The transformed intensity image can be
## - None
## - Identity (the input image)
## - Normalization_to_u8 (a normalized version of the input image on 1 byte)
## - Cell_Normalization_to_u8 (a cell-based normalized version of the input image on 1 byte)
## The membrane-enhanced image can be
## - None
## - GACE: 'Global Automated Cell Extractor'
## - GLACE:
##
## This choice can be tuned with the two parameters
## - astec_intensity_transformation
## - astec_intensity_enhancement
##
## astec_intensity_transformation can be chosen in [None, 'Identity', 'Normalization_to_u8', 'Cell_Normalization_to_u8']
## astec_intensity_enhancement can be chosen in [None, 'GACE', 'GLACE']
##
## The 'Cell_Normalization_to_u8' intensity transformation method can be tuned with
## - astec_cell_normalization_min_method
## - astec_cell_normalization_min_method
## both variable are to be chosen in ['global', 'cell', 'cellborder', 'cellinterior', 'voxel']
## Choosing both of them as 'global' comes to choose 'Normalization_to_u8'
##



# astec_keep_reconstruction = True

## ##### explanation #####
##
## Previous tuning enables to perform the watershed segmentation on an image that is not
## the original image. If one does not want to keep such images, please turn the next variable to False


# astec_cell_normalization_min_method = 'cellinterior'
# astec_cell_normalization_max_method = 'cellborder'
# normalization_min_percentile = 0.01
# normalization_max_percentile = 0.99
# cell_normalization_sigma = 5.0

## ##### explanation #####
##
## to be written



# astec_sigma_membrane = 0.9
# astec_hard_thresholding = 1.0
# astec_hard_threshold = 1.0
# astec_sensitivity = 0.99
# astec_manual = False
# astec_manual_sigma = 15
# astec_sigma_TV = 3.6
# astec_sigma_LF = 0.9
# astec_sample = 0.2

## ##### explanation #####
##
## GACE/GLACE parameters
## GACE/GLACE methods are made of three steps. GLACE is nothing but a cell-based GACE.
## 1. membrane detection
## 2. membrane segmentation
## 3. tensor voting
##
## - membrane detection:
## this is the gaussian sigma that is used to compute image derivatives
## (in real units, a priori 0.9 um is a good choice for data like Patrick/Ralph/Aquila)
## astec_sigma_membrane=0.9
##
## - membrane segmentation:
## segmentation of the computed extrema is done by thresholding.
## Either astec_hard_thresholding is set to True, and a global hard threshold (astec_hard_threshold)
## is applied to the whole extrema image. This is not advised, and should be used only when the
## adaptative thresholding has failed.
## Or astec_hard_thresholding is set to False, meaning that an adaptative thresholding is computed
## by histogram fitting. This adaptative thresholding is governed by three parameters
## - astec_sensitivity:
##   membrane binarization parameter, use larger values (smaller than or equal to 1.0) to increase the quantity
##   of binarized membranes to be used for tensor voting. Default value is 0.99.
##   # /!\ if failure, one should enter in "manual" mode of the function
##   anisotropicHist via activation of 'manual' option
## - astec_manual:
##   By default, this parameter is set to False. If failure, (meaning that thresholds are very bad,
##   meaning that the binarized image is very bad), set this parameter to True and relaunch the
##   computation on the test image. If the method fails again, "play" with the value of manual_sigma...
##   and good luck.
## - astec_manual_sigma:
##   Axial histograms fitting initialization parameter for the computation of membrane image binarization
##   axial thresholds (this parameter is used if manual = True). One may need to test different values of
##   manual_sigma. We suggest to test values between 5 and 25 in case of initial failure. Good luck.
##
## - tensor voting:
## tensor voting is governed by 3 parameters
## - astec_sigma_TV:
##   parameter which defines the voting scale for membrane structures propagation by tensor voting method (real
##   coordinates). This parameter shoud be set between 3 um (little cells) and 4.5 um
##   (big gaps in the binarized membrane image)
## - astec_sample:
##   Parameter for tensor voting computation speed optimisation (do not touch if not bewared)
##   set the fraction of binarized membrane image used for tensor voting (default is 0.2).
## - astec_sigma_LF:
##   Additional smoothing parameter for reconstructed image (in real coordinates).
##   It seems that the default value = 0.9 um is ok for standard use.
##



######################################################################
##
## SEGMENTATION PROPAGATION PARAMETERS
##
######################################################################


# propagation_strategy = ""

## ##### explanation #####
##
## Different segmentation propagation strategies are implemented
## - propagation_strategy = 'seeds_from_previous_segmentation'
##   In this strategy, a seeded watershed is used to segment the image at timepoint t+1
##   Seeds are built from the cells of the segmentation at timepoint t
##   The segmentation at timepoint t is then deformed onto the image at timepoint t+1
##   and cells are eroded

# previous_seg_erosion_cell_iterations = 10
# previous_seg_erosion_background_iterations = 25
# previous_seg_erosion_cell_min_size = 1000

## ##### explanation #####
##
## Parameters that control the seed construction from the deformed
## segmentation at timepoint t
## - previous_seg_erosion_cell_iterations: erosion size for the cells
## - previous_seg_erosion_background_iterations: erosion size for the background
## - previous_seg_erosion_cell_min_size: cells above this minimum size are discarded
##   for seed building
##



