# File defining the nomenclature for ASTEC experiments

import sys,os


# Definition of main stage directories of ASTEC package:
DIR_STAGE_RAW 	='RAWDATA'# RAW DATA stored in PATH_EMBRYO/DIR_STAGE_RAW
DIR_STAGE_FUSE	='FUSE'	# Fusion experiments stored in 
						# PATH_EMBRYO/DIR_STAGE_FUSE/FLAG_EXP_FUSE
DIR_STAGE_REG	='REG'	# Registration experiments stored in 
						# PATH_EMBRYO/DIR_STAGE_REG/FLAG_EXP_REG
DIR_STAGE_MARS	='SEG'	# Mars experiments stored in
						# PATH_EMBRYO/DIR_STAGE_MARS/FLAG_EXP_MARS
DIR_STAGE_SEG	='SEG'	# Segmentation propagation experiments stored in
						# PATH_EMBRYO/DIR_STAGE_SEG/FLAG_EXP_SEG
DIR_STAGE_POST	='POST'	# Post correction experiments stored in
						# PATH_EMBRYO/DIR_STAGE_SEG/FLAG_EXP_POST



# Flags to be replaced in paths
FLAG_PATH_EMBRYO='$PATHEMBRYO'# Must be the path to the embryo data 
FLAG_EN='$EN'				# Embryo Name (format YYMMDD-SaintOfTheDays-Stage)

FLAG_EXP_FUSE='$FUSE'		# defining fusion experiments subdirectory
FLAG_EXP_REG='$REG'			# defining registration experiments subdirectory
FLAG_EXP_MARS='$MARS'		# defining mars experiments subdirectory
FLAG_EXP_SEG='$SEG'			# defining seg propagation experiments subdirectory
FLAG_EXP_POST='$POST'		# defining post correction experiments subdirectory

FLAG_TIME='$TIME'			# Time-point of an embryo snapshot 
FLAG_TIMEREF='$TIMEREF'		# For registration,time-point of reference snapshot
FLAG_TIMEFLO='$TIMEFLO'		# For registration,time-point of floating snapshot

#FLAG_WORKSPACE='$WORKSPACE'# Workspace for the considered step; this workspace
						# is included in the repository corresponding to the
						# ASTEC step ([FUSE|SEG|POST]) and prefixed by this
						# repository name (eg. FUSE_EXP_NO_CROP or SEG_RELEASE)


### RAW DATA 
path_rawdata=os.path.join(FLAG_PATH_EMBRYO, DIR_STAGE_RAW)
# 1st image from the left camera, good quality image at the beginning
path_rawdata_angle1=os.path.join(path_rawdata,"LC"+os.path.sep+"Stack0000") 
# 1st image from the right camera
path_rawdata_angle2=os.path.join(path_rawdata,"RC"+os.path.sep+"Stack0000") 
# 2nd image from the left camera
path_rawdata_angle3=os.path.join(path_rawdata,"LC"+os.path.sep+"Stack0001") 
# 2nd from the right camera
path_rawdata_angle4=os.path.join(path_rawdata,"RC"+os.path.sep+"Stack0001") 


### FUSION DATA 
# path fused images
path_fuse          =os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_FUSE) 
# path for fusion workspace
path_fuse_exp      =os.path.join(path_fuse,DIR_STAGE_FUSE+'_'+FLAG_EXP_FUSE)
# fused images names
path_fuse_exp_files=os.path.join(path_fuse_exp,FLAG_EN+'_fuse_t$TIME.inr') 
# logfile
path_fuse_logfile  =os.path.join(path_fuse_exp,'1-fuse.log')	




### MARS DATA 
# path mars images
path_mars          =os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_MARS) 
# path for mars workspace
path_mars_exp      =os.path.join(path_mars,DIR_STAGE_MARS+'_'+FLAG_EXP_MARS)
# mars images names
path_mars_exp_files=os.path.join(path_mars_exp,FLAG_EN+'_mars_t$TIME.inr') 
# mars temporary images path (for reconstruction)
path_mars_exp_reconstruct=os.path.join(path_mars_exp,"RECONSTRUCTION") 
# mars temporary images file names (for reconstruction)
path_mars_exp_reconstruct_files=\
os.path.join(path_mars_exp_reconstruct,FLAG_EN+'_rec_t$TIME.inr') 
# logfile
path_mars_logfile  =os.path.join(path_mars_exp,'2-mars.log')	




### SEGMENTATION DATA 
# path segmented images
path_seg          =os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_SEG) 
# path for segmentation workspace
path_seg_exp      =os.path.join(path_seg,DIR_STAGE_SEG+'_'+FLAG_EXP_SEG)
# segmentated images names
path_seg_exp_files=os.path.join(path_seg_exp,FLAG_EN+'_seg_t$TIME.inr') 
# cells lineage file
path_seg_exp_lineage=os.path.join(path_seg_exp,FLAG_EN+'_seg_lineage.pkl') 
# segmentation temporary images path (for reconstruction)
path_seg_exp_reconstruct=os.path.join(path_seg_exp,"RECONSTRUCTION") 
# segmentation temporary images file name (for reconstruction)
path_seg_exp_reconstruct_files=\
os.path.join(path_seg_exp_reconstruct,FLAG_EN+'_rec_t$TIME.inr') 
# cells lineage test file
path_seg_exp_lineage_test=os.path.join(path_seg_exp,FLAG_EN+'_seg_lineage.test') 
# logfiles
path_seg_logfile  =os.path.join(path_seg_exp,'4-astec.log')	

### MANUAL CORRECTION DATA 
# logfiles
path_mancor_logfile  =os.path.join(path_seg_exp,'3-manualcorrection.log')	




### POST CORRECTION DATA 
# path post corrected images
path_post          =os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_POST) 
# path for post correction workspace
path_post_exp      =os.path.join(path_post, DIR_STAGE_POST+'_'+FLAG_EXP_POST)
# post corrected images names
path_post_exp_files=os.path.join(path_post_exp,FLAG_EN+'_post_t$TIME.inr') 
# post corrected cells lineage file
path_post_exp_lineage=os.path.join(path_post_exp,FLAG_EN+'_post_lineage.pkl') 
# post corrected cells lineage test file
path_post_exp_lineage_test=os.path.join(\
									path_post_exp,FLAG_EN+'_post_lineage.test')
# logfile
path_post_logfile  =os.path.join(path_post_exp,'5-postcorrection.log')	





### INTRA REGISTRATION DATA 
path_intrareg=os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_REG) # Path intra registration data
path_intrareg_exp=os.path.join(path_intrareg,DIR_STAGE_REG+'_'+FLAG_EXP_REG) # Path intra registration data
path_intrareg_step_files=os.path.join(path_intrareg_exp, \
		FLAG_EN+'_reg_t$TIMEFLO_t$TIMEREF.trsf') # Intra registration step-by-step
											# trsf file names
path_intrareg_multiple_files=os.path.join(path_intrareg, \
		FLAG_EN+'_reg_compose_t$TIME_t$TIME.trsf') # Intra registration composed 
											  #trsf file names
path_intrareg_change_files=os.path.join(path_intrareg, \
		FLAG_EN+'_reg_compose_t$TIME.trsf') # Intra registration recentered trsf 
									   # file names
path_intrareg_change_template=os.path.join(path_intrareg, \
		FLAG_EN+'_reg_compose_template.inr.gz') # Intra registration template file
										   # name for recentered trsfs
iso_intra_registration=1.0		# Parameter for intra registration resampled 
								# images resolution
















#postsegment_files=datapath+"GLACE/SEG/POST/"+EN+'_glas_seg_post_t$TIME.inr' #Segmentation output files

#INTRA REGISTRATION COMPOSED WITH ROTATION SO THAT GERMINAL CELLS ARE AT THE DOWN OF THE IMAGE # NOT USED NOW
intrareg_germinal_Path=path_intrareg_exp+"COMPOSE_GERMINAL/" # path intra registration data under germinal cells reorientation constraint
intrareg_germinal_file=intrareg_germinal_Path+FLAG_EN+'_germinal.trsf' # transformation (rotation) to compose with all the intra-registration trsf files 
intrareg_germinal_files=intrareg_germinal_Path+FLAG_EN+'_germinal_t$TIME.trsf' #  intra registration recentered trsf file names
intrareg_germinal_template=intrareg_germinal_Path+FLAG_EN+'_germinal_template.inr.gz' #  intra registration template file name for recentered trsfs




#RECONSTRUCTION DATA 
reconstruct_Path=os.path.join(path_fuse_exp,"RECONSTRUCTION") #path reconstructed images
reconstruct_files=os.path.join(reconstruct_Path,FLAG_EN+'_rec_t$TIME.inr') #  reconstructed images names



#SEGMENTED DATA 
segmented_Path=os.path.join(FLAG_PATH_EMBRYO,DIR_STAGE_SEG) #segmented images
mars_file=segmented_Path+FLAG_EN+'_fuse_mars_t$TIME.inr' #Segmentation output files
segmentation_files=segmented_Path+FLAG_EN+'_fuse_seg_t$TIME.inr' #Segmentation output files
lineage_tree_filename=segmented_Path+FLAG_EN+'_fuse_seg_lin_tree.pkl' #The main lineage tree file output 
lineage_tree_test_filename=segmented_Path+FLAG_EN+'_fuse_seg_lin_tree.test' #The main lineage tree test file output 

#MAPPING FILE
mapping_path=segmented_Path+FLAG_EN+"_fuse_seg_t$TIME.map"

#POST SEGMENTED DATA 
postsegment_Path=segmented_Path+"POST/" #post segmentation images
postsegment_files=postsegment_Path+FLAG_EN+'_fuse_seg_post_t$TIME.inr' #Segmentation output files
post_lineage_tree_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree.pkl' #The main lineage tree file output 
post_lineage_tree_test_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree.test' #The main lineage tree test file output


#CELL NAMES
name_file=os.path.join(FLAG_PATH_EMBRYO,FLAG_EN+"-names.txt")
named_lineage_tree_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree_named.pkl' #The main lineage tree file output 
named_lineage_tree_test_filename=postsegment_Path+FLAG_EN+'_fuse_seg_post_lin_tree_named.test' #The main lineage tree test file output 


###################################################################
####### VIRTUAL EMBRYO SHOULD BE TAKEN OUT OF ASTEC PACKAGE #######
###################################################################

#VIRTUAL EMBRYO
#vo_login="emmanuel.faure"
#vo_passwd="ascidie"
#vo_lib_Path=astec_Path+"3DCloudEmbryo/" 		   # Path for scripts library
#vo_Path=postsegment_Path+"VO/" 					   # Output path for mesh data
#vo_files=vo_Path+EN+'_fuse_seg_post_vo_t$TIME.obj' # Output files





##############################################################
################ FUNCTIONS FOR FLAGS REPLACEMENT #############
##############################################################

def replaceTIME(filename,time,_format="%03d"):
    """
    Replaces all the occurences of "$TIME" by its value given by time of 
    type int at the format specified by _format argument (default="%03d")
    """
    time_point=_format%time #Format time on 3 digit
    return filename.replace(FLAG_TIME, time_point)

def replaceTIMEFLO(filename,time,_format="%03d"):
    """
    Replaces all the occurences of "$TIMEFLO" by its value given by time of 
    type int at the format specified by _format argument (default="%03d")
    """
    time_point=_format%time #Format time on 3 digit
    return filename.replace(FLAG_TIMEFLO, time_point)

def replaceTIMEREF(filename,time,_format="%03d"):
    """
    Replaces all the occurences of "$TIMEREF" by its value given by time of 
    type int at the format specified by _format argument (default="%03d")
    """
    time_point=_format%time #Format time on 3 digit
    return filename.replace(FLAG_TIMEREF, time_point)

def replaceTimes(filename,d,_format="%03d"):
    """
    Replaces all the occurences of each key from the dictionnary d given in 
    parameter with its corresponding value of type int at the format 
    specified by _format argument (default="%03d")
    """
    for k,time in d.iteritems():
        assert type(k)==str, "Unexpected non str dictionnary key '%s'"%str(k)
        if type(time)!=int:
            print "Non int val '%s' specified for key '%s'. Trying a cast."%(str(time),k)
            time=int(time)
        assert filename.find(k)>=0
        time_point=_format%time
        filename=filename.replace(k, time_point)
    return filename

def replaceEN(filename,embryoname, check_name=False):
    """
    Replaces all the occurences of "$EN" by its value given by EN of type str
    check_name: parameter enabling to check the EN name composition wrt the
    			nomenclature (default=False)
    """
    if check_name:
        assert embryoname.count('-')!=2, 'ERROR in embryo name %s\n\
            -> EN should follow the template YYMMDD-SaintOfTheDays-Stage'\
            %embryoname
    return filename.replace(FLAG_EN, embryoname)

def replacePATH_EMBRYO(filename,path):
    """
    Replaces all the occurences of "$EMBRYOPATH" by its value given by EN of type str
    """
    #assert path, "Specified embryo path should not be empty."
    return filename.replace(FLAG_PATH_EMBRYO, path.rstrip(os.path.sep))

#def replaceWORKSPACE(filename,path):
#    """
#    Replaces all the occurences of "$WORKSPACE" by its value given by EN of type str
#    """
#    assert path, "Specified workspace should not be empty."
#    return filename.replace(FLAG_WORKSPACE, path)

def replaceEXP(filename,flag,path):
    """
    Replaces all the occurences of "$FUSE" by its value given by path of type str
	If path is empty, "$FUSE" is replaced by "RELEASE" and a message is printed
    """
    if not path:
        print ""
        path='RELEASE'
    return filename.replace(flag, path)

def replaceFlags(filename, parameters, check_name=False):
    """
    Function that replaces in the specified str filename the flags found
    (which verify the regular expression r'\$[A-Z]*', e.g. '$EN' or $EXP_SEG)
    with the corresponding str found in the parameters (of type 'module'):
    	FLAG_PATH_EMBRYO -> parameters.PATH_EMBRYO
    	FLAG_EN -> parameters.EN
    	FLAG_EXP_FUSE -> parameters.EXP_FUSE
    	FLAG_EXP_REG -> parameters.EXP_REG
    	FLAG_EXP_MARS -> parameters.EXP_MARS
    	FLAG_EXP_SEG -> parameters.EXP_SEG
    	FLAG_EXP_POST -> parameters.EXP_POST
    """
    import re
    found=re.findall(r'\$[A-Z]*',filename)
    for flag in found:
        if flag==FLAG_PATH_EMBRYO:
            filename=replacePATH_EMBRYO(filename, parameters.PATH_EMBRYO)
        elif flag==FLAG_EN:
            if not parameters.EN:
                embryoname=\
              parameters.PATH_EMBRYO.rstrip(os.path.sep).split(os.path.sep)[-1]
                filename=replaceEN(filename, embryoname, check_name=check_name)
            else:
                filename=replaceEN(filename, parameters.EN, \
                                   check_name=check_name)
        elif flag==FLAG_EXP_FUSE:
            filename=replaceEXP(filename, flag, parameters.EXP_FUSE)
        elif flag==FLAG_EXP_REG:
            filename=replaceEXP(filename, flag, parameters.EXP_REG)
        elif flag==FLAG_EXP_MARS:
            filename=replaceEXP(filename, flag, parameters.EXP_MARS)
        elif flag==FLAG_EXP_SEG:
            filename=replaceEXP(filename, flag, parameters.EXP_SEG)
        elif flag==FLAG_EXP_POST:
            filename=replaceEXP(filename, flag, parameters.EXP_POST)
        else:
        	print "replaceFlags: flag '%s' was not replaced"%flag
    return filename