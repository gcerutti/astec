# File for parameters and nomenclature settings

begin=5 #First Point
end=5 #Last Point

# RAWDATA DEFINITION

delta = 1 	 				# Delta between two time points (if one does not want to fuse every single time point)
ori = 'left' 				# if im2 angle - im1 angle < 0 => right
resolution = (.17, .17, 1.) # Resolution of the raw images (here are the known values for 140317-Patrick-St8)
resolution = (.21, .21, 1.) # Resolution of the raw images for Karine
delay = 0 					# If the time stamps in the folder are not the actual time stamps in the global movie
mirrors = False  			# Depends on the acquisition protocol, value can be set to True or False
							#  - the standard value of this parameter is False
							#  - in case of axial symmetry between left and right cameras, then set to True


# FUSION PARAMETERS

target_resolution = .3   # Isotropic resolution of the final fused image
fusion_dilation_x_0 = 120 #  dilation_x_0 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'left' x direction 
fusion_dilation_x_1 = 40 #  dilation_x_1 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'right' x direction
fusion_dilation_y_0 = 80 #  dilation_y_0 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'top' y direction
fusion_dilation_y_1 = 40 #  dilation_y_1 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'bottom' y direction
fusion_no_crop = False   #  no_crop [default=False]: if True, then the resampled image is not cropped


# MARS PARAMETERS

# modules choice
Mars_methods=['Classic','Gace','Hybridation']
Mars_method=1 # 1 for 'Classic' method, 2 for 'Gace' method, 3 for 'Hybridation' method 

### Parameters for MARS segmentation
sigma1_mars = 0.6  # sigma 1 (0.6um) in real coordinates
sigma2_mars = 0.15 # sigma 2 (0.15um) in real coordinates
h_min_mars = 4     # H min initialisation to ease correction


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







#FIND PATH AND EMBRYO NAME
import sys,os

astec_Path=os.getcwd()
#astec_Path="/home/gmicheli/TEST_ASTEC/ASTEC/160708-Aquila-St8/ASTEC-180326"
astec_Path="/home/gmicheli/MyData/DIGEM/TEST_ASTEC/RAW/171107-Karine-St8/astec-package"

tab_Path=astec_Path.split('/')
EN=tab_Path[len(tab_Path)-2]  #Embryo Name
if  EN.count('-')!=2:
	print 'ERROR in astec path '+astec_Path+ ' or embryo name '+EN
	print  '-> path should be /media/DATA/<embryoname>/ASTEC-YYMMDD'
	print '    where embryoname=YYMMDD-SaintOfTheDays-Stage'
	quit()

datapath=astec_Path[:astec_Path.rfind('/')+1]
astec_Path+='/'

print 'Embryo '+EN+' from '+str(begin)+' to '+str(end)

# ASTEC LIBARRIES
sys.path.append(astec_Path+"ASTEC") #Add the ASTEC Function
sys.path.append(astec_Path+'ASTEC/CommunFunctions')



#Image Path definition
rawdata_Path=datapath+"RAWDATA/"

path_angle1=rawdata_Path+"LC/Stack0000" # 1st image from the left camera, good quality image at the beginning
path_angle2=rawdata_Path+"RC/Stack0000" # 1st image from the right camera
path_angle3=rawdata_Path+"LC/Stack0001" # 2nd image from the left camera
path_angle4=rawdata_Path+"RC/Stack0001" # 2nd from the right camera


#FUSION DATA 
fuse_Path=datapath+"FUSE/" #path fused images
fused_files=fuse_Path+EN+'_fuse_t$TIME.inr' #  fused images names

#INTRA REGISTRATION DATA 
intrareg_Path=datapath+"FUSE/REG/" #path intra registration data
intrareg_step_files=intrareg_Path+EN+'_reg_t$TIMEFLO_t$TIMEREF.trsf' #  intra registration step-by-step trsf file names
intrareg_multiple_files=intrareg_Path+EN+'_reg_compose_t$TIME_t$TIME.trsf' #  intra registration composed trsf file names
intrareg_change_files=intrareg_Path+EN+'_reg_compose_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_change_template=intrareg_Path+EN+'_reg_compose_template.inr.gz' #  intra registration template file name for recentered trsfs
iso_intra_registration=1.0

#postsegment_files=datapath+"GLACE/SEG/POST/"+EN+'_glas_seg_post_t$TIME.inr' #Segmentation output files

#INTRA REGISTRATION COMPOSED WITH ROTATION SO THAT GERMINAL CELLS ARE AT THE DOWN OF THE IMAGE # NOT USED NOW
intrareg_germinal_Path=intrareg_Path+"COMPOSE_GERMINAL/" # path intra registration data under germinal cells reorientation constraint
intrareg_germinal_file=intrareg_germinal_Path+EN+'_germinal.trsf' # transformation (rotation) to compose with all the intra-registration trsf files 
intrareg_germinal_files=intrareg_germinal_Path+EN+'_germinal_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_germinal_template=intrareg_germinal_Path+EN+'_germinal_template.inr.gz' #  intra registration template file name for recentered trsfs

#RECONSTRUCTION DATA 
reconstruct_Path=fuse_Path+"RECONSTRUCTION/" #path reconstructed images
reconstruct_files=reconstruct_Path+EN+'_rec_t$TIME.inr' #  reconstructed images names

#SEGMENTED DATA 
segmented_Path=fuse_Path+"SEG/" #segmented images
mars_file=segmented_Path+EN+'_fuse_mars_t$TIME.inr' #Segmentation output files
segmentation_files=segmented_Path+EN+'_fuse_seg_t$TIME.inr' #Segmentation output files
lineage_tree_filename=segmented_Path+EN+'_fuse_seg_lin_tree.pkl' #The main lineage tree file output 
lineage_tree_test_filename=segmented_Path+EN+'_fuse_seg_lin_tree.test' #The main lineage tree test file output 

#MAPPING FILE
mapping_path=segmented_Path+EN+"_fuse_seg_t$TIME.map"

#POST SEGMENTED DATA 
postsegment_Path=segmented_Path+"POST/" #post segmentation images
postsegment_files=postsegment_Path+EN+'_fuse_seg_post_t$TIME.inr' #Segmentation output files
post_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.pkl' #The main lineage tree file output 
post_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.test' #The main lineage tree test file output


#CELL NAMES
name_file=datapath+EN+"-names.txt"
named_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.pkl' #The main lineage tree file output 
named_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.test' #The main lineage tree test file output 


###################################################################
####### VIRTUAL EMBRYO SHOULD BE TAKEN OUT OF ASTEC PACKAGE #######
###################################################################

#VIRTUAL EMBRYO
vo_login="emmanuel.faure"
vo_passwd="ascidie"
vo_lib_Path=astec_Path+"3DCloudEmbryo/" #Path for scripts library
vo_Path=postsegment_Path+"VO/" #Output path for mesh data
vo_files=vo_Path+EN+'_fuse_seg_post_vo_t$TIME.obj' # output files

