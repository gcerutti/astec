#!/usr/bin/python2.7


import os, sys, imp

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__), \
		"ASTEC","CommunFunctions"))
sys.path.append(os.path.join(os.path.dirname(__file__), \
		"ASTEC","CommunFunctions"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"), \
    "CommunFunctions","cpp"), "Unable to find the 'cpp' library link in %s,\
    please install properly the library and make a logical link to its bin \
    repository at the path %s."%\
    (os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"), \
    "CommunFunctions"), \
    os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"), \
    "CommunFunctions","cpp"))

from optparse import OptionParser

from nomenclature import *
from REGISTRATION import readLUT
from ImageHandling import imread, imsave, SpatialImage
import numpy as np


### Options parsing

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="p",
                  help="python file containing parameters definition",\
                  metavar="FILE")
parser.add_option("-e", "--embryo-rep", dest="e",
                  help="path to the embryo data", metavar="PATH")
parser.add_option("-m", "--mapping-file", dest="m",
                  help="path to the manual cell to cell mapping file", metavar="MAPPING")
parser.add_option("-q", "--quiet",
                  action="store_false", dest="verbose", default=True,
                  help="don't print status messages to stdout")

(options, args) = parser.parse_args()

parameters_file=options.p

### Parameters file
if not parameters_file:
    parameters_file=raw_input('Provide the parameters file: ')
try:
    assert os.path.isfile(parameters_file)
except ValueError:
    print "Provided path '%s' was not found. Exiting."%parameters_file
    sys.exit(1)

### Reading parameters from parameters_file, put in the instance p
try:
    p = imp.load_source('*', parameters_file) # p is a structure which contains
                                              # all the parameters defined in
                                              # parameters_file
except ValueError:
    print "Unable to load the source '%s'. Exiting."%parameters_file
    sys.exit(1)


### Particular options parsing: redefining p.PATH_EMBRYO

if options.e:
    if p.PATH_EMBRYO:
        test=raw_input("Caution: parameter PATH_EMBRYO '%s' defined in '%s'\
              will be overwritten by '%s'. Are you sure? (Y/n)"\
              %(p.PATH_EMBRYO,parameters_file,options.e))
        if test != 'Y' and test != 'y':
            print "Usage: option -e should be used when PATH_EMBRYO option is\
             kept empty in parameters file. Exiting."
            sys.exit(1) 
    p.PATH_EMBRYO = options.e
    # Embryo name from path: defining p.EN (also set if options.e is provided)
    if p.EN:
        test=raw_input("Caution: parameter EN '%s' defined in '%s'\
              will be overwritten. Are you sure? (Y/n)"\
              %(p.PATH_EMBRYO,parameters_file))
        if test != 'Y' and test != 'y':
            print "Usage: option '-e' should be used when PATH_EMBRYO and EN \
             options are kept empty in parameters file. Exiting."
            sys.exit(1) 
    p.EN='' # replaceFlags function will deal with empty p.EN field

### Particular options parsing: redefining p.mancor_mapping_file

if options.m:
    assert not p.mancor_mapping_file, "Usage: option -m should be used when \
    mancor_mapping_file option is kept empty in parameters file. Exiting."
    p.mancor_mapping_file=options.m

### Building paths from nomenclature.py and parameters file


path_mars_exp = replaceFlags(path_mars_exp, p)
print "Mars data will be searched in directory %s"%replaceFlags(path_mars_exp,
																 p)
assert os.path.isdir(path_mars_exp), "Provided fuse directory '%s' not found"\
									 %path_mars_exp

path_mars_exp_files = replaceFlags(path_mars_exp_files, p)

path_seg = replaceFlags(path_seg, p)
path_seg_exp = replaceFlags(path_seg_exp, p)
path_seg_exp_files = replaceFlags(path_seg_exp_files, p)
path_log_file = replaceFlags(path_mancor_logfile, p)


### Segmentation directory and subdirectory construction

if not os.path.isdir(path_seg):
    os.mkdir(path_seg)  
if not os.path.isdir(path_seg_exp):
    os.mkdir(path_seg_exp)  


### Log file

if os.path.exists(os.getcwd()+os.path.sep+'.git'):
    os.system('echo "# astec-package version: `git describe`" >> %s'\
        %path_log_file)
else:
    os.system('echo "# astec-package version was not found, meaning it was not\
     cloned from the inria forge project" >> %s'%path_log_file)
with open(path_log_file, 'a') as log_file:
    log_file.write('# Python executable: '+sys.executable+'\n')
    log_file.write('# Working directory: %s\n'%os.getcwd())
    log_file.write('# Parameters file: %s\n'% parameters_file)
    log_file.write('# Mapping file: %s\n'% p.mancor_mapping_file)
    log_file.write('# Embryo path: %s\n'% p.PATH_EMBRYO)
    log_file.write('# Embryo name: %s\n'% p.EN)
    log_file.write('# Command line:\n'+(' '.join(sys.argv))+'\n\n\n')

### Copy of parameters file

os.system("cp -f "+parameters_file+" "+path_seg_exp )
os.system("cp -f "+p.mancor_mapping_file+" "+path_seg_exp )





#######################################
### Manual correction Process Stuff ###
#######################################

#Manual correction of the first time point

mars_file = replaceTIME(path_mars_exp_files,p.begin+p.raw_delay) #Input seg
segmentation_file = replaceTIME(path_seg_exp_files,p.begin+p.raw_delay) #Output

seg = imread(mars_file)

mapping = np.arange(np.max(seg)+1)

'''
# Example of mapping setting:
mapping[8] = 7	  # This mapping is going to set voxels of label 8 with the new value 7
mapping[9] = 4
# ... etc ...
mapping[29] = 23
mapping[30] = 1   # 1 corresponds to the background label
mapping[89] = 1	  # 1 corresponds to the background label

---> mapping text file:
# Random commentary
8 7
9
4
# ... etc ...
29 23
30 1 
89 1
'''

theMap={}
if p.mancor_mapping_file:
    assert os.path.exists(p.mancor_mapping_file), "File '%s' not found."%p.mancor_mapping_file
    theMap=readLUT(p.mancor_mapping_file)
for k,v in theMap.iteritems():
    mapping[k] = v

seg_corrected = mapping[seg]

cells_list = np.unique(seg_corrected)
nb_cells = len(cells_list)

print "list of cells by ids: " + str(cells_list)
print "total cells: " + str(nb_cells)

imsave(segmentation_file, SpatialImage(seg_corrected, voxelsize=seg.voxelsize).astype(np.uint16)) # Save into segmentation file as inr /!\ HERE WE OVERWRITE ON THE ORIGINAL SEGMENTATION FILE (GAEL)

