############################################################
##
## useful parameters for segmentation propagation
##
## to be used with 5-postcorrection-new.py
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

## ##### explanation #####
##
## time points to be processed
##
## begin: first time point
## end: last time point
## delta: delta/time interval between two time points (if one does not want
##        to deal with every single time point) (default = 1)
##



# EXP_SEG = ''

## ##### explanation #####
##
## Segmentation results will be stored in 'PATH_EMBRYO'/SEG/SEG_'EXP_SEG'
## default value is
## EXP_SEG = 'RELEASE'
##
# EXP_SEG = ''



## ##### explanation #####
##
## Post-correction results will be stored in 'PATH_EMBRYO'/POST/SEG_'EXP_POST'
## default value is
## EXP_POST = 'RELEASE'
##



# result_image_suffix = 'mha'
# default_image_suffix = 'mha'

## ##### explanation #####
##
## defines the output image format
##
## result_image_suffix is for all output images
## default_image_suffix is for all images (including auxiliary ones)
## default is 'mha'
##
## 'mha' should be prefered (readable by fiji)
##



# result_lineage_suffix = 'xml'

## ##### explanation #####
##
## defines the output lineage format
##
## default is 'xml'
##
## 'xml' should be prefered (there are portability issues between python 2.x and 3.x)
##



######################################################################
##
## post-correction parameters
##
######################################################################


# volume_minimal_value = 2000

## ##### explanation #####
##
## volume (in voxels) threshold
## end branches candidate for fusion either end before the end of the sequence
## or have a leaf (terminal) cell whose volume is below this threshold
##



# lifespan_minimal_value = 25

## ##### explanation #####
##
## end branches whose length is (stricly) below this threshold
## will be deleted
##



# test_early_division = False

## ##### explanation #####
##
## if test_early_division is set to True
## it will be tested whether the division that gives rise to the end branch
## is too close to the first division of its sister branch, or along the mother trunk
## (too close means less than lifespan_minimal_value)
##



# test_volume_correlation = True
# correlation_threshold = 0.9

## ##### explanation #####
##
## if test_volume_correlation is set to True
## it will be tested whether there is an anti-correlation between the cell volumes
## of the end branch and the ones of its sister branch
## the branch is deleted is the Pearson correlation coefficient is less than - correlation_threshold
##
