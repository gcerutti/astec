

# File for parameters and nomenclature settings

PATH_EMBRYO=''	
                # Must be the path to the embryo data 
				# eg: '/media/DATA/171107-Karine-St8'
				# Can also be set by providing option '-e' in the command line
EN=''			
                # Embryo Name 
                # CRBM naming format is YYMMDD-SaintOfTheDays-Stage
				# eg: '171107-Karine-St8'
				# (automatically extracted from PATH_EMBRYO if not provided)


DIR_TO_LEFTCAM_STACKZERO=''
DIR_TO_LEFTCAM_STACKONE=''
DIR_TO_RIGHTCAM_STACKZERO=''
DIR_TO_RIGHTCAM_STACKONE=''

EXP_FUSE=''	
            # Workspace for the considered step; this workspace
			# is included in the repository corresponding to the
			# ASTEC step ([FUSE|SEG|POST]) and prefixed by this
			# repository name (eg. FUSE_EXP_NO_CROP or SEG_RELEASE)
			# ----> use replaceDIR_EXP_FUSE(path_fuse_exp[...], EXP_FUSE) 
			# method from <astec-package>/nomenclature.py to get the right
			# path for the given fusion experience
EXP_REG=''	
            #
EXP_MARS=''	
            #
EXP_SEG=''	
            #
EXP_POST=''	
            #


#TIME=''		# Time-point of an embryo snapshot 
#TIMEREF=''		# For registration, time-point of reference snapshot
#TIMEFLO=''		# For registration, time-point of floating snapshot

##########################
### GENERAL PARAMETERS ###
##########################

begin=1 				 
                         # First time point
end=5   				 
                         # Last time point
delta=1     			 
                         # Delta between two time points (if one does not want
						 # to deal with every single time point) (default = 1)
target_resolution = .3   
                         # Isotropic resolution of the final fused and 
						 # segmented images (default = 0.3)


##########################
### RAWDATA DEFINITION ###
##########################

raw_ori = 'left' 				
                                # if im2 angle - im1 angle < 0 => right
raw_resolution = (.17, .17, 1.) 
                                # Resolution of the raw images (here are the 
								# known values for 140317-Patrick-St8)
#raw_resolution = (.21, .21, 1.) 
                                # Resolution of the raw images for Karine
raw_delay = 0 					
                                # If the time stamps in the folder are not the
								# actual time stamps in the global movie
raw_mirrors = False  			
                                # Depends on the acquisition protocol, value 
								# can be set to True or False
								#  - the standard value of this parameter is 
								#    False
								#  - in case of axial symmetry between left 
								#    and right cameras, then set to True




#########################
### FUSION PARAMETERS ###
#########################


raw_margin_x_0 = 40 
                       # margin_x_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # raw acquisition image in 'left' x direction 
raw_margin_x_1 = 40 
                       # margin_x_1 [default=40]: parameter for margin of the
				       # bounding box computed for the cropping of the 
				       # raw acquisition image in 'right' x direction
raw_margin_y_0 = 40 
                       # margin_y_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # raw acquisition image in 'top' y direction
raw_margin_y_1 = 40 
                       # margin_y_1 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # raw acquisition image in 'bottom' y direction
raw_crop = True     
                       # crop [default=True]: 
                       # if False, then the resampled image is not cropped ; 
					   # if raw acquisition, then image is cropped


fusion_margin_x_0 = 40 
                       # margin_x_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'left' x direction 
fusion_margin_x_1 = 40 
                       # margin_x_1 [default=40]: parameter for margin of the
				       # bounding box computed for the cropping of the 
				       # resampled image in 'right' x direction
fusion_margin_y_0 = 40 
                       # margin_y_0 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'top' y direction
fusion_margin_y_1 = 40 
                       # margin_y_1 [default=40]: parameter for margin of the
					   # bounding box computed for the cropping of the 
					   # resampled image in 'bottom' y direction
fusion_crop = True     
                       # crop [default=True]: 
                       # if False, then the resampled image is not cropped ; 
					   # if True, then image is cropped




#######################
### MARS PARAMETERS ###
#######################

# Modules choice
mars_method=1 			
                        # 1 for 'Classic' method
			  			# 2 for 'Gace' method
# General parameters for MARS segmentation
mars_sigma1 = 0.6  		
                        # sigma 1 (0.6um) in real coordinates
mars_sigma2 = 0.15 		
                        # sigma 2 (0.15um) in real coordinates
mars_h_min = 4     		
                        # H min initialisation to ease correction
# Gace Parameters (if mars_method is set to 2):
# membrane_renforcement
mars_sigma_membrane=0.9 
                        # membrane enhancement parameter (in real units, a 
						# priori 0.9 um is a good choice for data like 
						# Patrick/Ralph/Aquila)
# anisotropicHist /!\ critical step
mars_sensitivity=0.99   
                        # membrane binarization parameter, /!\ if failure,
						# one should enter in "manual" mode of the function
						# anisotropicHist via activation of 'manual' option

mars_manual=False     	
                        # By default, this parameter is set to False. If 
						# failure, (meaning that thresholds are very bad, 
						# meaning that the binarized image is very bad),
				 		# set this parameter to True and relaunch the 
				 		# computation on the test image. If the method fails
				 		# again, "play" with the value of manual_sigma... 
				 		# and good luck.
mars_manual_sigma=15    
                        # Axial histograms fitting initialization parameter 
						# for the computation of membrane image binarization
						# axial thresholds (this parameter is used iif 
						# manual = True).
						# One may need to test different values of 
						# manual_sigma. We suggest to test values between 5 and
						# 25 in case of initial failure. Good luck.

mars_hard_thresholding=False  
                              # If the previous membrane threshold method 
							  # failed, one can force the thresholding with a
							  # "hard" threshold applied on the whole image. 
							  # To do so, this option must be set to True.
mars_hard_threshold=1.0       
                              # If hard_thresholding = True, the enhanced 
							  # membranes image is thresholded using this 
							  # parameter (value 1 seems to be ok for 
							  # time-point t001 of Aquila embryo for example).

# Tensor voting framework
mars_sigma_TV=3.6     
                      # parameter which defines the voting scale for membrane
					  # structures propagation by tensor voting method (real
					  # coordinates). 
				 	  # This parameter shoud be set between 3 um (little cells)
				 	  # and 4.5 um(big gaps in the binarized membrane image)
mars_sigma_LF=0.9     
                      # Smoothing parameter for reconstructed image (in real
					  # coordinates). It seems that the default value = 0.9 um
					  # is ok for classic use.
mars_sample=0.2       
                      # Parameter for tensor voting computation speed 
					  # optimisation (do not touch if not bewared)



####################################
### MANUAL CORRECTION PARAMETERS ###
####################################

mancor_seg_file='' 	   
                       # segmentation file to be manually corrected
					   # If not provided, then looking for the output of MARS
mancor_mapping_file='' 
                       # path to mapping file for manual correction of the 
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

astec_membrane_reconstruction_method=0
# Membrane reconstruction module choice
# 0 for 'Classic' method
# 1 for 'Glace' method
# 2 for 'Gace' method
# If not set or set to 0, the input fused image is not processed for 
# membrane structures enhancement.
# If set to 1, the GLACE reconstruction method is going to be called
# If set to 2, the GACE reconstruction method is going to be called

astec_fusion_u8_method=0
# Selection of the method which converts the fused image into a 8 bits
# images for the segmentation propagation. 
# If set to 0 (default), calling the historical "to_u8" method
# If set to 1, calling the mc-adhocFuse program which enhances the fused 
# image while converting it to u8 knowing the semgnetation propagation from
# previous time point

astec_flag_hybridation=False
# If set to True and if the membrane_reconstruction_method parameter is 
# provided and not equal to 0, then the reconstructed gray level image
# used for the segmentation propagation framework is goind to be an hybridation 
# between the original fused image and the result of image reconstruction by
# the specified method.

astec_keep_reconstruct_files=False 
# Set it to True in order to keep a copy
# Flag enabling to keep a copy of graylevel files provided to the watershed



# General parameters for segmentation propagation
astec_sigma1 = 0.6  		
                            # sigma 1 (0.6um) in real coordinates
astec_sigma2 = 0.15 		
                            # sigma 2 (0.15um) in real coordinates
astec_h_min_min = 4
                            # H min initialisation to ease correction
astec_h_min_max = 18   		
                            # H min initialisation to ease correction

# Glace Parameters (if astec_membrane_reconstruction_method is set to 1 or 2):
# membrane_renforcement
astec_sigma_membrane=0.9
                        # membrane enhancement parameter (in real units, a
						# priori 0.9 um is a good choice for data like 
						# Patrick/Ralph/Aquila)
# anisotropicHist /!\ critical step
astec_sensitivity=0.99  
                        # membrane binarization parameter, /!\ if failure,
						# one should enter in "manual" mode of the function
						# anisotropicHist via activation of 'manual' option

astec_manual=False     	
                        # By default, this parameter is set to False. If 
						# failure, (meaning that thresholds are very bad, 
						# meaning that the binarized image is very bad),
				 		# set this parameter to True and relaunch the 
				 		# computation on the test image. If the method fails
				 		# again, "play" with the value of manual_sigma... 
				 		# and good luck.
astec_manual_sigma=15   
                        # Axial histograms fitting initialization parameter 
						# for the computation of membrane image binarization
						# axial thresholds (this parameter is used iif 
						# manual = True).
						# One may need to test different values of 
						# manual_sigma. We suggest to test values between 5 and
						# 25 in case of initial failure. Good luck.

astec_hard_thresholding=False 
                              # If the previous membrane threshold method 
							  # failed, one can force the thresholding with a
							  # "hard" threshold applied on the whole image. 
							  # To do so, this option must be set to True.
astec_hard_threshold=1.0      
                              # If hard_thresholding = True, the enhanced 
							  # membranes image is thresholded using this 
							  # parameter (value 1 seems to be ok for 
							  # time-point t001 of Aquila embryo for example).

# Tensor voting framework
astec_sigma_TV=3.6    
                      # parameter which defines the voting scale for membrane
					  # structures propagation by tensor voting method (real
					  # coordinates). 
				 	  # This parameter shoud be set between 3 um (little cells)
				 	  # and 4.5 um(big gaps in the binarized membrane image)
astec_sigma_LF=0.9    
                      # Smoothing parameter for reconstructed image (in real
					  # coordinates). It seems that the default value = 0.9 um
					  # is ok for classic use.
astec_sample=0.2      
                      # Parameter for tensor voting computation speed 
					  # optimisation (do not touch if not bewared)
astec_rayon_dil=3.6   
                      # dilatation ray for propagated ROI from time t to t+1
					  # (default: 3.6, in real coordinates) 

# Fused image conversion parameters (if astec_fusion_u8_method is set to 1)
astec_min_percentile=0.01   
                      # mc-adhocFuse parameter of type %f (default: 0.01)
astec_max_percentile=0.99   
                      # mc-adhocFuse parameter of type %f (default: 0.99)
astec_min_method='cellinterior'
                      # mc-adhocFuse param. (default: 'cellinterior')
				      # taken in global|cell|cellborder|cellinterior|voxel
astec_max_method='cellborder'
                      # mc-adhocFuse parameter (default: 'cellborder')
					  # taken in global|cell|cellborder|cellinterior|voxel
astec_sigma_hybridation=5.0 
                      # mc-adhocFuse parameter of type %f (default: 5.0)

# Default parameters (for classical use, default values should not be changed)
astec_RadiusOpening=20 		
                      # (using the approximation of a sphere of radius 20
					  # voxels as a structuring element)
astec_Thau= 25 				
                      # s(c)=h2+(c).N2(c) >t identical
astec_MinVolume=1000 		
                      # Miss Suppressing cells with to small volumes (not
					  # yet in supdata)
astec_VolumeRatioBigger=0.5 
                      # If a cell in St+1 is at least 50% bigger than its
					  # progeny in St+1, 
astec_VolumeRatioSmaller=0.1
                      # Cells in St+1 that are 10% or smaller than their
					  # equivalent in St+1 are tagged for correction
astec_MorphosnakeIterations=10 
                      # Then, an active contour algorithm is applied
				      # using the dilated shape of c, obtained by 
					  # iterating 10 times
astec_NIterations=200 		
                      # The algorithm is applied up to stability 
					  # (at th voxels) or after n iterations 
					  # (in our case th = 103 and n = 200). 
astec_DeltaVoxels=10**3  	
                      # y (at th voxels)
astec_Volum_Min_No_Seed=100 
                      # Then, if the volume of c is greater than 100 
					  # voxels (2.7 um3)
astec_nb_proc=10 			
                      # Number of processor ...
astec_nb_proc_ace=7   		
                      # number of processors for ACE (7 is recommanded)



##################################
### POST-CORRECTION PARAMETERS ###
##################################

postcor_Volume_Threshold=10000 	
                      # volume low threshold for final cells 
					  # (in voxels) 
postcor_Soon=True 				
                      # True if the cell life span has to be 
					  # taken into account
postcor_ShortLifespan=25 		
                      # (length < SL time points, in our case 
					  # SL = 10)
postcor_PearsonThreshold=0.9; 	
                      # If they are anticorrelated (Pearson 
					  # correlation under -0.9)



