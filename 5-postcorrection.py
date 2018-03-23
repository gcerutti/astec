# Post correction of the segmentation
from definitions import *
if not os.path.isdir(postsegment_Path):
    os.mkdir(postsegment_Path) 
os.system("cp -f "+astec_Path+"definitions.py "+postsegment_Path )
os.system("cp -f "+astec_Path+"5-postcorrection.py "+postsegment_Path ) 


from lineage import write_tlp_from_lin_tree,read_lineage_tree,write_lineage_tree,timeNamed,timesNamed
from post_correction import apply_cell_fusion,remove_too_little_branches
from lineage_test import pkl_lineage_test, imageDict

lin_tree_information=read_lineage_tree(lineage_tree_filename) # Read the lineage tree (in case it was previously created)



# Lineage Tree Post-correction
### PARAMETERS
Volume_Threshold=10000 # volume low threshold for final cells (in voxels) 
Soon=True #True if the cell life span has to be taken into account
ShortLifespan=25 # (length < SL time points, in our case SL = 10)
PearsonThreshold=0.9; #If they are anticorrelated (Pearson correlation under -0.9)


### CORRECT THE LINEAGE TREE WITH THE CELLS VOLUMES
lin_tree_cor, new_volumes, to_fuse, been_fused=remove_too_little_branches(lin_tree_information['lin_tree'], lin_tree_information['volumes_information'], Volume_Threshold, soon=Soon)


### APPLYING THE CORRECTION ON THE IMAGES
apply_cell_fusion(lin_tree_information['lin_tree'],lin_tree_information['volumes_information'], to_fuse,segmentation_files,postsegment_files,begin, end, delta)

#SAVE THE NEW LINEAGE TREE
lin_tree_information['lin_tree']=lin_tree_cor


### SAVE THE FINAL LINEAGE TREE
write_lineage_tree(post_lineage_tree_filename,lin_tree_information)

### PROCESS LINEAGE TREE FILE VERIFICATION
print 'PROCESS LINEAGE TREE VERIFICATION'
image_dict_seg_post=imageDict(postsegment_files.replace("$TIME","*"))
rapport=pkl_lineage_test(lin_tree_information, image_dict_seg_post, file_out=post_lineage_tree_test_filename)

print 'LINEAGE TREE VERIFICATION DONE'

