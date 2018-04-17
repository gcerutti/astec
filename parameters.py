

# File for parameters and nomenclature settings

PATH_EMBRYO=''	# Must be the path to the embryo data 
				# eg: '/media/DATA/171107-Karine-St8'
				# Can also be set by providing option '-e' in the command line
EN=''			# Embryo Name (format YYMMDD-SaintOfTheDays-Stage)
				# eg: '171107-Karine-St8'
				# (automatically extracted from PATH_EMBRYO if not provided)


EXP_FUSE=''	# Workspace for the considered step; this workspace
			# is included in the repository corresponding to the
			# ASTEC step ([FUSE|SEG|POST]) and prefixed by this
			# repository name (eg. FUSE_EXP_NO_CROP or SEG_RELEASE)
			# ----> use replaceDIR_EXP_FUSE(path_fuse_exp[...], EXP_FUSE) 
			# method from <astec-package>/nomenclature.py to get the right
			# path for the given fusion experience
EXP_REG=''	#
EXP_MARS=''	#
EXP_SEG=''	#
EXP_POST=''	#


#TIME=''		# Time-point of an embryo snapshot 
#TIMEREF=''		# For registration, time-point of reference snapshot
#TIMEFLO=''		# For registration, time-point of floating snapshot

##########################
### GENERAL PARAMETERS ###
##########################

begin=5 				 # First Point
end=5   				 # Last Point
delta=1     			 # Delta between two time points (if one does not want
						 # to deal with every single time point) (default = 1)
target_resolution = .3   # Isotropic resolution of the final fused and 
						 # segmented images (default = 0.3)


##########################
### RAWDATA DEFINITION ###
##########################

raw_ori = 'left' 				# if im2 angle - im1 angle < 0 => right
raw_resolution = (.17, .17, 1.) # Resolution of the raw images (here are the 
								# known values for 140317-Patrick-St8)
#raw_resolution = (.21, .21, 1.) # Resolution of the raw images for Karine
raw_delay = 0 					# If the time stamps in the folder are not the
								# actual time stamps in the global movie
raw_mirrors = False  			# Depends on the acquisition protocol, value 
								# can be set to True or False
								#  - the standard value of this parameter is 
								#    False
								#  - in case of axial symmetry between left 
								#    and right cameras, then set to True




#########################
### FUSION PARAMETERS ###
#########################

fusion_margin_x_0 = 40 # margin_x_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'left' x direction 
fusion_margin_x_1 = 40 # margin_x_1 [default=40]: parameter for margin of the
				       # bounding box computed for the cropping of the 
				       # resampled image in 'right' x direction
fusion_margin_y_0 = 40 # margin_y_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'top' y direction
fusion_margin_y_1 = 40 # margin_y_1 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'bottom' y direction
fusion_crop = True     # crop [default=True]: if False, then the resampled 
					   # image is not cropped ; if True, then image is cropped




#######################
### MARS PARAMETERS ###
#######################

# Modules choice
mars_method=1 			# 1 for 'Classic' method
			  			# 2 for 'Gace' method
# General parameters for MARS segmentation
mars_sigma1 = 0.6  		# sigma 1 (0.6um) in real coordinates
mars_sigma2 = 0.15 		# sigma 2 (0.15um) in real coordinates
mars_h_min = 4     		# H min initialisation to ease correction
# Gace Parameters (if mars_method is set to 2):
# membrane_renforcement
mars_sigma_membrane=0.9 # membrane enhancement parameter (in real units, a 
						# priori 0.9 um is a good choice for data like 
						# Patrick/Ralph/Aquila)
# anisotropicHist /!\ critical step
mars_sensitivity=0.99   # membrane binarization parameter, /!\ if failure,
						# one should enter in "manual" mode of the function
						# anisotropicHist via activation of 'manual' option

mars_manual=False     	# By default, this parameter is set to False. If 
						# failure, (meaning that thresholds are very bad, 
						# meaning that the binarized image is very bad),
				 		# set this parameter to True and relaunch the 
				 		# computation on the test image. If the method fails
				 		# again, "play" with the value of manual_sigma... 
				 		# and good luck.
mars_manual_sigma=15    # Axial histograms fitting initialization parameter 
						# for the computation of membrane image binarization
						# axial thresholds (this parameter is used iif 
						# manual = True).
						# One may need to test different values of 
						# manual_sigma. We suggest to test values between 5 and
						# 25 in case of initial failure. Good luck.

mars_hard_thresholding=False  # If the previous membrane threshold method 
							  # failed, one can force the thresholding with a
							  # "hard" threshold applied on the whole image. 
							  # To do so, this option must be set to True.
mars_hard_threshold=1.0       # If hard_thresholding = True, the enhanced 
							  # membranes image is thresholded using this 
							  # parameter (value 1 seems to be ok for 
							  # time-point t001 of Aquila embryo for example).

# Tensor voting framework
mars_sigma_TV=3.6     # parameter which defines the voting scale for membrane
					  # structures propagation by tensor voting method (real
					  # coordinates). 
				 	  # This parameter shoud be set between 3 um (little cells)
				 	  # and 4.5 um(big gaps in the binarized membrane image)
mars_sigma_LF=0.9     # Smoothing parameter for reconstructed image (in real
					  # coordinates). It seems that the default value = 0.9 um
					  # is ok for classic use.
mars_sample=0.2       # Parameter for tensor voting computation speed 
					  # optimisation (do not touch if not bewared)



####################################
### MANUAL CORRECTION PARAMETERS ###
####################################

mancor_mapping_file='' # path to mapping file for manual correction of the 
					   # mars segmentation. See above the syntax of this file.
					   # - 1 line per label association
					   # - background label has value 1
					   # - the character '#' denotes commented lines 
					   # See file "mapping.txt" in the astec project to get an
					   # example
'''
# EXAMPLE OF mancor_mapping_file CONTENT:
# here the input label 8 will be mapped with new value 7, etc...
8 7
9 2  
4 64 
29 23
# ... etc ...
# background labels
30 1 
89 1 
'''

###########################################
### SEGMENTATION PROPAGATION PARAMETERS ###
###########################################

# Modules choice
astec_method=1
# Glace parameters (if astec_method is set to 2 or 3):
# membrane_renforcement
sigma_membrane=0.9 # membrane enhancement parameter (in real units, a priori
				   # 0.9 um is a good choice for data like Patrick/Ralph/...)
# anisotropicHist /!\ critical step
sensitivity=0.99 # membrane binarization parameter, /!\ if failure, one should
				 # enter in "manual" mode of the function anisotropicHist via 
				 # activation of 'manual' option

manual=False     # By default, this parameter is set to False. If failure, 
				 # (meaning that thresholds are very bad, meaning that the 
				 # binarized image is very bad),
				 # set this parameter to True and relaunch the computation on
				 # the test image. If the method fails again, "play" with the
				 # value of manual_sigma... and good luck.
manual_sigma=15  # Axial histograms fitting initialization parameter for the
				 # computation of membrane image binarization axial thresholds
				 # (this parameter is used iif manual = True).
				 # One may need to test different values of manual_sigma. We
				 # suggest to test values between 5 and 25 in case of initial
				 # failure. Good luck.

hard_thresholding=False  # If the previous membrane threshold method failed,
						 # one can force the thresholding with a "hard" 
						 # threshold applied on the whole image. To do so, this
						 # option must be set to True.
hard_threshold=1.0       # If hard_thresholding = True, the enhanced membranes
						 # image is thresholded using this parameter (value 1
						 # seems to be ok for time-point t001 of Aquila embryo
						 # for example).

# TVmembrane
sigma_TV=3.6     # parameter which defines the voting scale for membrane 
				 # structures propagation by tensor voting method (real 
				 # coordinates). 
				 # This parameter shoud be set between 3 um (little cells) and
				 # 4.5 um(big gaps in the binarized membrane image)
sigma_LF=0.9     # Smoothing parameter for reconstructed image (in real
				 # coordinates). It seems that the default value = 0.9 um is ok
				 # for classic use.
sample=0.2       # Parameter for tensor voting computation speed optimisation
				 # (do not touch if not bewared)




