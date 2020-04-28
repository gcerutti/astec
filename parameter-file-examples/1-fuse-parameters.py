############################################################
##
## useful parameters for fusion
##
## to be used with 1-fuse.py
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
# delta = 1

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



acquisition_orientation = 'left'
acquisition_mirrors = False
acquisition_leftcamera_stacking_direction = 'increasing'
acquisition_resolution = (.17, .17, 1.)
# raw_delay = 0

## ##### explanation #####
##
## raw data parameters
##
## acquisition_orientation: image orientation ('right' or 'left')
##   gives the rotation (wrt to the Y axis) of the left camera frame of stack #0
##   to be aligned with the the left camera frame of stack #1
##   - 'right': +90 degrees
##   - 'left': -90 degrees
## acquisition_mirrors: mirroring of the right camera image along the X-axis
## 	 Right camera images may have to be mirrored along the X-axis to be aligned
##   with the left camera images
##   - 'True': the mirroring is done
##   - 'False': the mirroring is still to be done
##   Since it should depend on the apparatus, all acquisitions performed by the same
##   microscope should be the same tuning
## acquisition_leftcamera_stacking_direction
##   gives the order of stacking of in the Z direction
##   'increasing': from the high-contrasted images to the fuzzy ones
##   'decreasing': from the fuzzy images to the high-contrasted ones
## acquisition_resolution: acquisition voxel size
##   e.g. raw_resolution = (.21, .21, 1.)
## raw_delay: increment to to be added to the time values
##   (values in range [begin,end]) when generating the
##   fused image names (default is 0)
##   eg: acquisition at time point 't' results in the fused image at time 't+raw_delay'
## acquisition_leftcamera_stacking_direction
##
## To determine the configuration (raw_ori,raw_resolution) (ie ('left',False),
## ('left',True), ('right', False), or ('right', True)), it is advised to perform the fusion
## fusion for only one time point (by setting 'begin' and 'end' at the same value) with 
## a large 'target_resolution'
##



# DIR_RAWDATA = ''
# DIR_LEFTCAM_STACKZERO = ''
# DIR_RIGHTCAM_STACKZERO = ''
# DIR_LEFTCAM_STACKONE = ''
# DIR_RIGHTCAM_STACKONE = ''

# DIR_RAWDATA_CHANNEL_2 = ''
# DIR_LEFTCAM_STACKZERO_CHANNEL_2 = ''
# DIR_RIGHTCAM_STACKZERO_CHANNEL_2 = ''
# DIR_LEFTCAM_STACKONE_CHANNEL_2 = ''
# DIR_RIGHTCAM_STACKONE_CHANNEL_2 = ''

# DIR_RAWDATA_CHANNEL_3 = ''
# DIR_LEFTCAM_STACKZERO_CHANNEL_3 = ''
# DIR_RIGHTCAM_STACKZERO_CHANNEL_3 = ''
# DIR_LEFTCAM_STACKONE_CHANNEL_3 = ''
# DIR_RIGHTCAM_STACKONE_CHANNEL_3 = ''

## ##### explanation #####
##
## Paths to the raw data
##
## allows to specify paths to the raw (eq acquisition) data ## the 4 paths are built as follows
## <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_LEFTCAM_STACKZERO>
## <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_RIGHTCAM_STACKZERO>
## <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_LEFTCAM_STACKONE>
## <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_RIGHTCAM_STACKONE>
##
## default values are
##
## DIR_RAWDATA = 'RAWDATA'
## DIR_LEFTCAM_STACKZERO = 'LC/Stack0000'
## DIR_RIGHTCAM_STACKZERO = 'RC/Stack0000'
## DIR_LEFTCAM_STACKONE = 'LC/Stack0001'
## DIR_RIGHTCAM_STACKONE = 'RC/Stack0001'
##
## When temporary files are kept (option -k), LEFTCAM_STACKZERO, RIGHTCAM_STACKZERO
## LEFTCAM_STACKONE, and RIGHTCAM_STACKONE related files are respectively stored in
## directories FUSE/FUSE_<EXP_FUSE>/TEMP_<XXX>/ANGLE_[0,1,2,3]
## 
## multi-channel fusion
## --------------------
## Fusion transformations computed for the first channel are used for the fusion of the
## other channels
##
## paths to the raw data (of the other channels) have also to be specified
## with the variables (X = 2 or 3) and the same path construction
## eg: <PATH_EMBRYO>/<DIR_RAWDATA_CHANNEL_X>/<DIR_LEFTCAM_STACKZERO_CHANNEL_X>
## DIR_RAWDATA_CHANNEL_X
## DIR_LEFTCAM_STACKZERO_CHANNEL_X
## DIR_RIGHTCAM_STACKZERO_CHANNEL_X
## DIR_LEFTCAM_STACKONE_CHANNEL_X
## DIR_RIGHTCAM_STACKONE_CHANNEL_X
##
## if DIR_RAWDATA_CHANNEL_X is not given, it is replaced by DIR_RAWDATA or its default value
## if any of the four path DIR_dirCAM_STACKx_CHANNEL_X is not given, it is also replaced
## by the main channel variable (DIR_dirCAM_STACKx) or its default value
##



# acquisition_leftcam_image_prefix = ''
# acquisition_rightcam_image_prefix = ''

## ##### explanation #####
##
## raw data file names
##
## allows to specify prefix name of the raw data
## may be useful when several camera acquisitions are stored into
## the same directory. Eg, if one has set
##   DIR_RAWDATA='raw_data'
##   DIR_LEFTCAM_STACKZERO = 'Stack_0_Channel_0'
##   DIR_RIGHTCAM_STACKZERO = 'Stack_0_Channel_0'
##   DIR_LEFTCAM_STACKONE = 'Stack_1_Channel_0'
##   DIR_RIGHTCAM_STACKONE = 'Stack_1_Channel_0'
## meaning that stacks #0 of both left and right cameras are stored into the
## same directory 'raw_data/Stack_0_Channel_0'
##
## File names for the left camera acquisition are then built with
##   <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_LEFTCAM_STACKZERO>/<acquisition_leftcam_image_prefix>+<TIME>
##   <PATH_EMBRYO>/<DIR_RAWDATA>/<DIR_LEFTCAM_STACKONE>/<acquisition_leftcam_image_prefix>+<TIME>
## where <TIME> is a 3-digit number (with 0 padding). The same for the right camera
## acquisition files. Default values are
##
## acquisition_leftcam_image_prefix = 'Cam_Left_00'
## acquisition_rightcam_image_prefix = 'Cam_Right_00'
##



# EXP_FUSE = ''
# EXP_FUSE_CHANNEL_2 = ''
# EXP_FUSE_CHANNEL_3 = ''

## ##### explanation #####
##
## fusion directory names
##
## Fusion results will be stored in <PATH_EMBRYO>/FUSE/FUSE_<EXP_FUSE>
## default value is
## EXP_FUSE = 'RELEASE'
##
## for multi-channel fusion, a fusion result per additional channel has to be given
## EXP_FUSE_CHANNEL_2 = ''
## EXP_FUSE_CHANNEL_3 = ''
##



target_resolution = .3

## ##### explanation #####
##
## isotropic voxel size of the fusion result (fused image)
##



# RESULT_IMAGE_SUFFIX_FUSE = 'mha'
# result_image_suffix = 'mha'
# default_image_suffix = 'mha'

## ##### explanation #####
##
## defines the output image format
##
## RESULT_IMAGE_SUFFIX_FUSE is for the fusion result image(s)
## result_image_suffix is for all output images
## default_image_suffix is for all images (including auxiliary ones)
## default are 'inr'
##
## 'mha' should be prefered (readable by fiji)
##



#####################################################################
##
## fusion parameters (advanced)
##
######################################################################

## ##### explanation #####
##
## the fusion of the 4 acquisitions follows a number of steps
##
## 1. Optionally, a slit line correction.
##    Some Y lines may appear brighter or darker in the acquisition, which
##    may cause artifacts in the reconstructed (ie fused) image, which, in
##    turn, may impair further segmentation steps.
## 2. a change of resolution in the X and Y directions only (Z remains unchanged)
##    it allows to decrease the data volume if the new pixel size is larger
##    than the acquisition one
## 3. Optionally, a crop of the resampled acquisitions
##    it allows to decrease the volume of data
##    the crop is based on the analysis of a MIP view (in the Z direction) of
##    the volume
## 4. Optionally, a mirroring of the 'right' image
##    depends of the <raw_mirrors> value (see supra)
## 5. Linear registration of the 3 last images on the first one (considered as the reference)
##    The reference image is resampled again, to get an isotropic voxel
##    (same voxel size in the 3 directions: X, Y, Z)
## 6. Linear combination of images, weighted by an ad-hoc function
## 7. Crop of the fused image
##    still based on the analysis of a MIP view (in the Z direction)
##



# acquisition_slit_line_correction = True

##
## step 1. Slit line corrrection
## default is False
##



##
## step 2. resampling in X and Y with 'target_resolution'
## Z remains unchanged
##



# raw_crop = True
# raw_margin_x_0 = 40
# raw_margin_x_1 = 40
# raw_margin_y_0 = 40
# raw_margin_y_1 = 40

##
## step 3. raw data cropping parameters
##
## raw_crop: if False, then the acquisition images are not cropped
##   (default is True)
## raw_margin_x_0: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'left' X direction
##   (default is 40)
## raw_margin_x_1: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'right' X direction
##   (default is 40)
## raw_margin_y_0: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'left' Y direction
##   (default is 40)
## raw_margin_y_1: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'right' Y direction
##   (default is 40)
##



##
## step 4. mirroring
## depends on the <raw_mirrors> value
##



# fusion_strategy = 'direct-fusion'

# fusion_preregistration_compute_registration = False
# fusion_preregistration_transformation_type = 'translation'
# fusion_preregistration_transformation_estimation_type = 'wlts'
# fusion_preregistration_lts_fraction = 0.55
# fusion_preregistration_pyramid_highest_level = 6
# fusion_preregistration_pyramid_lowest_level = 3
# fusion_preregistration_normalization = True

# fusion_registration_compute_registration = True
# fusion_registration_transformation_type = 'affine'
# fusion_registration_transformation_estimation_type = 'wlts'
# fusion_registration_lts_fraction = 0.55
# fusion_registration_pyramid_highest_level = 6
# fusion_registration_pyramid_lowest_level = 3
# fusion_registration_normalization = True

# fusion_stack_preregistration_compute_registration = True
# fusion_stack_preregistration_transformation_type = 'affine'
# fusion_stack_preregistration_transformation_estimation_type = 'wlts'
# fusion_stack_preregistration_lts_fraction = 0.55
# fusion_stack_preregistration_pyramid_highest_level = 6
# fusion_stack_preregistration_pyramid_lowest_level = 3
# fusion_stack_preregistration_normalization = True

# fusion_stack_registration_compute_registration = True
# fusion_stack_registration_transformation_type = 'vectorfield'
# fusion_stack_registration_transformation_estimation_type = 'wlts'
# fusion_stack_registration_lts_fraction = 1.0
# fusion_stack_registration_pyramid_highest_level = 6
# fusion_stack_registration_pyramid_lowest_level = 3
# fusion_stack_registration_normalization = True

##
## step 5. co-registration of acquisitions ('direct-fusion' and 'hierarchical-fusion')
##

## ##### explanation #####
##
## there are two ways to perform the fusion of the 4 acquisitions:
##
## 1. 'direct-fusion'
##    each acquisition is linearly co-registered with the first acquisition (stack #0, left camera).
##    Then weights and images are transformed thanks to the computed transformations.
##    Finally a weighted linear combination gives the result.
## 2. 'hierarchical-fusion'
##    from the couple (left camera, right camera), each stack is reconstructed, following the same scheme than
##    the direct fusion but with only 2 images. Then stack#1 is (non-)linearly co-registered with stack #0.
##    Images and weights associated with stack#1 are then (non-)linearly transformed.
##    Finally a weighted linear combination gives the result.
##

##
## fusion_preregistration_* and fusion_registration_* control the co-registration of two acquisitions
## it is either used in the 'direct-fusion' method (to co-register each acquisition onto the first one)
## or in the 'hierarchical-fusion' method (to co-register couple of opposite acquisitions to reconstruct stacks)
##
## fusion_registration_pyramid_highest_level: highest level of the pyramid image for registration
##   registration is done hierarchically with a pyramid of images at each level, image dimensions are divided by 2.
##   'fusion_registration_pyramid_highest_level = 6' means that registration starts with images whose dimensions
##   are 1/64th of the original image
## fusion_registration_pyramid_lowest_level: lowest level of the pyramid image for registration
##

##
## fusion_stack_preregistration_* and fusion_stack_registration_* control the co-registration of two stacks
## it is only used in the 'hierarchical-fusion' method (to co-register the reconstructed stacks)
##



# fusion_weighting = 'guignard-weighting'
# fusion_weighting_channel_1 = 'guignard-weighting'
# fusion_weighting_channel_2 = 'guignard-weighting'
# fusion_weighting_channel_3 = 'guignard-weighting'

##
## step 6. weighted linear combination of the resampled co-registered acquisitions
##

## ##### explanation #####
##
## the weighted linear combination of the 4 co-registered stacks can be tuned by the 'fusion_weighting'
## variable
##
## - 'uniform': uniform (or constant) weighting, it comes to the average of the resampled co-registered stacks
## - 'ramp': the weights are linearly increasing or decreasing along the Z axis
## - 'corner': the weights are constant in a corner portion of the stack, defined by two diagonals in the XZ-section
## - 'guignard': original historical weighting function, described in Leo Guignard's Phd thesis, that puts more weight
##    to sections close to the camera and take also account the traversed material
##
## The variable 'fusion_weighting' allows to set the fusion weighting for all the channels to be processed. Using
## the variables 'fusion_weighting_channel_X' allows to set different weighting schemes for each channel.
##



# fusion_crop = True
# fusion_margin_x_0 = 40
# fusion_margin_x_1 = 40
# fusion_margin_y_0 = 40
# fusion_margin_y_1 = 40

##
## step 7. cropping of the fused image
##
## fusion_crop: if False, then the acquisition images are not cropped
##   (default is True)
## fusion_margin_x_0: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'left' X direction
##   (default is 40)
## fusion_margin_x_1: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'right' X direction
##   (default is 40)
## fusion_margin_y_0: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'left' Y direction
##   (default is 40)
## fusion_margin_y_1: parameter for margin of the bounding box computed
##   for the cropping of the  raw acquisition image in 'right' Y direction
##   (default is 40)
##
