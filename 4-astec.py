# ASTEC  segmentation propagation
from definitions import *
if not os.path.isdir(segmented_Path):
    os.mkdir(segmented_Path)  
os.system("cp -f "+astec_Path+"definitions.py "+segmented_Path )
os.system("cp -f "+astec_Path+"4-astec.py "+segmented_Path ) 

#Import Functions
from ImageHandling import imread, imsave, SpatialImage
from ASTEC import segmentation_propagation
from lineage import read_lineage_tree,write_lineage_tree,timeNamed,timesNamed
from lineage_test import pkl_lineage_test, imageDict

### Parameters:
sigma1 = 0.6 #/ target_resolution   #sigma 1 (0.6um)

h_min_min=2 #(in our case hmin = 2 and hmax = 18).
h_min_max=18 #(in our case hmin = 2 and hmax = 18).
#h_min_max=50 #(in our case hmin = 2 and hmax = 50).
RadiusOpening=20 #(using the approximation of a sphere of radius 20 voxels as a structuring element)
Thau= 25 # s(c)=h2+(c).N2(c) >t IDENTIQUE 
MinVolume=1000 #Miss Suppressing cells with to small volumes (not yet in supdata)
VolumeRatioBigger=0.5#If a cell in St+1 is at least 50% bigger than its progeny in St+1, 
VolumeRatioSmaller=0.1 #Cells in St+1 that are 10% or more smaller than their equivalent in St+1 are tagged for correction. 
MorphosnakeIterations=10 #Then, an active contour algorithm is applied using the dilated shape of c, obtained by iterating 10 times
NIterations=200 # The algorithm is applied up to stability (at th voxels) or after n iterations (in our case th = 103 and n = 200). 
DeltaVoxels=10**3  #y (at th voxels)
Volum_Min_No_Seed=100 #Then, if the volume of c is greater than 100 voxels (2.7 um3)
    
nb_proc=10 # Number of processor ...


lin_tree_information=read_lineage_tree(lineage_tree_filename) # Read the lineage tree (in case it was previously created)

####SAFETY CHECK AFTER RELAUNCH
if 'lin_tree' in lin_tree_information:
    import numpy as np
    cellat={}
    for y in lin_tree_information['lin_tree']:
        t=y/10**4
        if t not in cellat:
            cellat[t]=1
        else:
            cellat[t]+=1

    restart=-1
    t=begin
    while restart==-1 and t<=end:
        time_segment=t+delta #Time point of Segmentation 
        segmentation_file=timeNamed(segmentation_files,time_segment) #Output Segmentation file
        if not os.path.isfile(segmentation_file): 
            print 'Miss segmentation file at ' +str(t) + ' -> '+segmentation_file
            restart=t
        else:
            if cellat[t]==0:
                    print 'Miss lineage at '+str(t)
                    restart=t
            else:
                try :
                    seg=imread(segmentation_file)
                except IOError:
                    print ' Error in '+segmentation_file
                    restart=t
        if restart==-1:
            print ' Safe segmentation at '+str(t)
        t+=1
    print ' --> Restart at '+str(restart)
    begin=restart
   



### PROCESS PROPAGATION SEGMENTATION 
for t in range(begin, end):
    time_segment=t+delta #Time point of Segmentation 
    print 'Starting the segmentation at ' + str(time_segment)
    fused_file_ref=timeNamed(fused_files,t) #Previous image file
    fused_file=timeNamed(fused_files,time_segment) #Actual image file to segment
    segmentation_file_ref=timeNamed(segmentation_files,t) #Previous Segmentation file
    segmentation_file=timeNamed(segmentation_files,time_segment) #Output Segmentation file
    temporary_folder=timeNamed(segmented_Path+'TEMP_$TIME/',t) #  TEMPORARY FOLDER
    os.system("mkdir -p " + temporary_folder )#  TEMPORARY FOLDER

    vf_file=timesNamed(temporary_folder+'VF_t$TIME1_on_t$TIME2.inr','$TIME1',t,'$TIME2',time_segment) #VECTOR FIELDS FILE
    h_min_files=timeNamed(temporary_folder+'h_min_t$TIME_h$HMIN_s$SIGMA.inr.gz',time_segment)  #HMIN FILES
    seed_file=timeNamed(temporary_folder+'Seed_t$TIME.inr',t) #SEEDS FILE

    #PROCESS PROGATION SEGMENTATION
    seg_from_opt_h, lin_tree_information=segmentation_propagation(t,fused_file_ref,segmentation_file_ref, fused_file, seed_file,vf_file , h_min_files, h_min_min,h_min_max, sigma1, lin_tree_information, delta, nb_proc,
        RadiusOpening=RadiusOpening,Thau=Thau,MinVolume=MinVolume,VolumeRatioBigger=VolumeRatioBigger,VolumeRatioSmaller=VolumeRatioSmaller,MorphosnakeIterations=MorphosnakeIterations,NIterations=NIterations,DeltaVoxels=DeltaVoxels)
    
    #SAVE OUTPUT
    print 'Write the segmentation in ' + segmentation_file
    imsave(segmentation_file, seg_from_opt_h)
    write_lineage_tree(lineage_tree_filename,lin_tree_information) #Save the current lineage tree
    os.system("rm -rf  " + temporary_folder ) #DELETE TEMPORATY FILES

print 'ASTEC SEGMENTATION DONE'

### PROCESS LINEAGE TREE FILE VERIFICATION
print 'PROCESS LINEAGE TREE VERIFICATION'
image_dict_seg=imageDict(segmentation_files.replace("$TIME","*"))
rapport=pkl_lineage_test(lin_tree_information, image_dict_seg, file_out=lineage_tree_test_filename)

print 'LINEAGE TREE FILE VERIFICATION DONE'
