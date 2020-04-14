#
# A filter-based membrane enhancement procedure is applied on the fused image
# the fused image is also casted (converted) into a 1 byte image
# the reconstructed image  (2019-Tutorial100_fuse_t000_membrane.inr)
# is built from those two image and is located in SEG/SEG_SEG02/RECONSTRUCTION/
#
# Fused images are in FUSE/FUSE_RELEASE
#

PATH_EMBRYO = '.'

EN = '2019-Tutorial100'

begin = 0

EXP_SEG = 'MARS_TEST04'
mars_intensity_transformation = 'Normalization_to_u8'
mars_intensity_enhancement = 'GACE'
default_image_suffix = 'mha'
