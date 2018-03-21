# ASTEC  intra-sequence registration
from definitions import *
if not os.path.isdir(intrareg_Path):
    os.mkdir(intrareg_Path)  
os.system("cp -f "+astec_Path+"definitions.py "+intrareg_Path )
os.system("cp -f "+astec_Path+"1.5-intraregistration.py "+intrareg_Path ) 

#Import Functions
from ImageHandling import imread, imsave, SpatialImage
from lineage import timeNamed,timesNamed #, read_lineage_tree,write_lineage_tree,
from cpp_wrapping import rigid_registration, multiple_trsfs, change_multiple_trsfs
#from ASTEC import segmentation_propagation
from REGISTRATION import compute_intra_sequence_two_by_two_registration, compute_intra_sequence_respatialization_trsfs, compute_optimized_intra_sequence_respatialization_trsfs, intra_sequence_realignment, compose_transformation_stack_with_a_transformation

iso=1.0

### Parameters:



### Stuff:

# Step-by-step


compute_intra_sequence_two_by_two_registration(fused_files, intrareg_step_files, begin=begin, end=end, delta=delta, verbose=True)

# Recomputing the trsfs with a unique reference

time_point_ref=compute_intra_sequence_respatialization_trsfs(intrareg_step_files, intrareg_multiple_files, begin=begin, end=end, verbose=True)

# Next step : usage of function "changeMultipleTrsfs" for 
#	- computing the minimal window for building a 3D+t subsequence of images of the full sequence 
#	- building a template image into which the original images will be transformed


if intrareg_multiple_files.count('$TIME')==2:
	intrareg_comp_format=intrareg_multiple_files.split('$TIME')
	intrareg_comp_format=intrareg_comp_format[0]+'$TIME'+intrareg_comp_format[1]+('%03d'%time_point_ref)+intrareg_comp_format[2]
else:
	intrareg_comp_format=intrareg_multiple_files

compute_optimized_intra_sequence_respatialization_trsfs(postsegment_files, intrareg_comp_format, intrareg_change_template, intrareg_change_files, begin, end, threshold=2, iso=1.0, margin=10, verbose=True)

# Resampling a sequence of segmented images in a unique referential


intra_sequence_realignment(postsegment_files, "TMP_COMPOSE/foo_compose_t$TIME.mha", intrareg_change_files, 
						   template_image=intrareg_change_template, begin=begin, end=end, delta=delta, nearest=True, visu=True, verbose=True)

################################################################################################
#### or for example, if one wants to apply transformation on the sequence of fused images : ####
################################################################################################
intra_sequence_realignment(fused_files, "TMP_COMPOSE/foo_fused_compose_t$TIME.mha", intrareg_change_files, 
						   template_image=intrareg_change_template, begin=begin, end=end, delta=delta, nearest=False, visu=False, verbose=True)


#############################################################################################################################
#### or for example, if an optimal template has not been computed before or if one does not need to keep this template : ####
#############################################################################################################################
#intra_sequence_realignment(postsegment_files, "TMP/foo_t$TIME.mha", intrareg_comp_format, 
#						   template_image=None, begin=begin, end=end, delta=delta, nearest=True, iso=1.0, threshold=2, margin=10, visu=True, verbose=True)


# If needed to compose with any transformation to rotate the movie onto a specific orientation

if os.path.exists(intrareg_germinal_file):
	# The following code will enable to compute a resampling of the sequence following a specific orientation 
	# (usually, this orientation is given by registering the reference time-point 'time_point_ref' of the present
	# sequence onto the reference embryo 140317-Patrick-St8 and by composing this transformation with the 
	# "germinal transformation" of Patrick)
	if not os.path.isdir(intrareg_germinal_Path):
		os.mkdir(intrareg_germinal_Path) 
	compose_transformation_stack_with_a_transformation(intrareg_comp_format, intrareg_germinal_file, intrareg_germinal_files, begin, end, verbose=True)
	intra_sequence_realignment(postsegment_files, "TMP/foo_germinal_t$TIME.mha", intrareg_germinal_files, 
						   template_image=None, begin=begin, end=end, delta=delta, nearest=True, iso=1.0, threshold=2, margin=10, visu=True, verbose=True)
