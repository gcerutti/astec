#Manual correction of the first time point
from definitions import *

from ImageHandling import imread, imsave, SpatialImage
from lineage import timeNamed,timesNamed
import numpy as np

segmentation_file = timeNamed(segmentation_files,begin) #First time step to segment


seg = imread(segmentation_file)

mapping = np.arange(np.max(seg)+1)
'''
mapping[8] = 7
mapping[9] = 4
mapping[11] = 6
mapping[21] = 14
mapping[22] = 16
mapping[19] = 1
mapping[24] = 13
mapping[29] = 23
mapping[30] = 1
mapping[33] = 27
mapping[43] = 34
mapping[47] = 42
mapping[45] = 42
mapping[54] = 48
mapping[56] = 42
mapping[67] = 62
mapping[65] = 61
mapping[69] = 66
mapping[82] = 73
mapping[78] = 72
mapping[88] = 80
mapping[86] = 79
mapping[89] = 1
mapping[87] = 1
'''
seg_corrected = mapping[seg]


imsave(segmentation_file, SpatialImage(seg_corrected, voxelsize=seg.voxelsize).astype(np.uint16)) # Save into segmentation file as inr
