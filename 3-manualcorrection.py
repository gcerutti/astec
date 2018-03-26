#Manual correction of the first time point
from definitions import *

from ImageHandling import imread, imsave, SpatialImage
from lineage import timeNamed,timesNamed
import numpy as np

segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment


seg = imread(segmentation_file)

mapping = np.arange(np.max(seg)+1)

'''
# Example of mapping setting:
mapping[8] = 7	  # This mapping is going to set voxels of label 8 with the new value 7
mapping[9] = 4
# ... etc ...
mapping[29] = 23
mapping[30] = 1   # 1 corresponds to the background label
mapping[89] = 1	  # 1 corresponds to the background label
'''

seg_corrected = mapping[seg]

cells_list = np.unique(seg_corrected)
nb_cells = len(cells_list)

print "list of cells by ids: " + str(cells_list)
print "total cells: " + str(nb_cells)

imsave(segmentation_file, SpatialImage(seg_corrected, voxelsize=seg.voxelsize).astype(np.uint16)) # Save into segmentation file as inr /!\ HERE WE OVERWRITE ON THE ORIGINAL SEGMENTATION FILE (GAEL)

# MERGE Sophia - Montpellier (UTILITE ?)
# <<<<
#imsave(segmentation_file.replace('.inr','_mars.tiff'), imread(segmentation_file))
# <<<<


