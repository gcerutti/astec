
# RAWDATA DEFINITION
begin=0
end=3 #Last Point

delta = 1 # Delta between two time points (if one does not want to fuse every single time point)
ori = 'right' # if im2 angle - im1 angle < 0 => right
resolution = (.17, .17, 1.) # Resolution of the raw images
delay = 0 # If the time stamps in the folder are not the actual time stamps in the global movie
mirrors = False  #TO COMMENT
target_resolution = .3 # Isotropic resolution of the final fused image

#FIND PATH AND EMBRYO NAME
import sys,os

astec_Path=os.getcwd()
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

#RECONSTRUCTION DATA 
reconstruct_Path=fuse_Path+"RECONSTRUCTION/" #path reconstructed images
reconstruct_files=reconstruct_Path+EN+'_rec_t$TIME.inr' #  reconstructed images names

#SEGMENTED DATA 
segmented_Path=fuse_Path+"SEG/" #segmented images
segmentation_files=segmented_Path+EN+'_fuse_seg_t$TIME.inr' #Segmentation output files
lineage_tree_filename=segmented_Path+EN+'_fuse_seg_lin_tree.pkl' #The main lineage tree file output 
lineage_tree_test_filename=segmented_Path+EN+'_fuse_seg_lin_tree.test' #The main lineage tree test file output 


#POST SEGMENTED DATA 
postsegment_Path=segmented_Path+"POST/" #post segmentation images
postsegment_files=postsegment_Path+EN+'_fuse_seg_post_t$TIME.inr' #Segmentation output files
post_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.pkl' #The main lineage tree file output 
post_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree.test' #The main lineage tree test file output


#CELL NAMES
name_file=datapath+EN+"-names.txt"
named_lineage_tree_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.pkl' #The main lineage tree file output 
named_lineage_tree_test_filename=postsegment_Path+EN+'_fuse_seg_post_lin_tree_named.test' #The main lineage tree test file output 


#VIRTUAL EMBRYO
vo_login="emmanuel.faure"
vo_passwd="ascidie"
vo_lib_Path=astec_Path+"3DCloudEmbryo/" #Path for scripts library
vo_Path=postsegment_Path+"VO/" #Output path for mesh data
vo_files=vo_Path+EN+'_fuse_seg_post_vo_t$TIME.obj' # output files

