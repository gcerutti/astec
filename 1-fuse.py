#Fusion Process
from definitions import *

print "begin="+str(begin)


from FUSION import read_raw_data,fusion_process
from lineage import timeNamed,timesNamed
import os

if not os.path.isdir(fuse_Path):
    os.mkdir(fuse_Path)  


#Search for image format in different angle folders
success1,begin1,end1,ext_im1,path_im1=read_raw_data(path_angle1) 
success2,begin2,end2,ext_im2,path_im2=read_raw_data(path_angle2)
success3,begin3,end3,ext_im3,path_im3=read_raw_data(path_angle3) 
success4,begin4,end4,ext_im4,path_im4=read_raw_data(path_angle4) 

if not success1==success2==success3==success4==1 :
     print 'Error in your files, please double check your path files '
elif not begin1==begin2==begin3==begin4:
     print 'Error in your angles file do not start at the same time point'
elif not end1==end2==end3==end4:
    print 'Error in your angles file do not end at the same time point'
elif not ext_im1==ext_im2==ext_im3==ext_im4 :
    print 'Error in your angles file do not have the same extension'
else:
    begin=begin1;end=end1;ext_im=ext_im1;
    print 'Process Fusion from ' + str(begin)+ ' to ' + str(end)
    angles_files=[path_im1, path_im2, path_im3, path_im4] #Combine Angle Path
    temporary_path=fuse_Path+"TEMP_$TIME/" #Temporary Path
    #PROCESS THE FUSION
    for time in range(begin, end+1, delta): # Interation on time steps
        fused_file=timeNamed(fused_files,time+delay)
        if not os.path.isfile(fused_file):
            time_angles_files=[timeNamed(angle_file + ext_im,time) for angle_file in angles_files]
            temporary_time_path=timeNamed(temporary_path,time) #Temporary Path for this time point
            print temporary_time_path
            time_process=fusion_process(time_angles_files,
                       fused_file,  
                       temporary_time_path,
                       ori, resolution,target_resolution, delay,
                       ext_im1, mirrors = mirrors, targetResolution=target_resolution)
            print "Time point " + str(time) + " takes " + str(time_process) + " to compute\n\n\n"
            os.system("rm -rf "+temporary_time_path) ### Cleaning temporary files
