PATH_EMBRYO = '.'

EN = '2019-Tutorial100'



begin = 0

#
# A filter-based membrane enhancement procedure is applied on the fused image
# the reconstructed image (2019-Tutorial100_fuse_t000_membrane.inr)
# is located in SEG/SEG_SEG02/RECONSTRUCTION/
#

EXP_SEG = 'PARAM04'
mars_intensity_transformation = 'Normalization_to_u8'
mars_intensity_enhancement = 'GACE'
default_image_suffix = 'mha'
