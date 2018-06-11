#!/usr/bin/python2.7


import os, sys, imp

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions"))
assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions","cpp")), "Unable to find the 'cpp' library link in %s,\
    please install properly the library and make a logical link to its bin \
    repository at the path %s."%\
    (os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions"), \
    os.path.join(os.path.dirname(__file__),"ASTEC", \
    "CommunFunctions","cpp"))

from optparse import OptionParser

from nomenclature import *
from MARS import mars_segmentation 


### Options parsing

parser = OptionParser()
parser.add_option("-p", "--parameters", dest="p",
                  help="python file containing parameters definition",\
                  metavar="FILE")
parser.add_option("-e", "--embryo-rep", dest="e",
                  help="path to the embryo data", metavar="PATH")
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




### Building paths from nomenclature.py and parameters file


path_fuse_exp = replaceFlags(path_fuse_exp, p)
print "Fused data will be searched in directory %s"%replaceFlags(path_fuse_exp,
																 p)
assert os.path.isdir(path_fuse_exp), "Provided fuse directory '%s' not found"\
									 %path_fuse_exp
path_fuse_exp_files = replaceFlags(path_fuse_exp_files, p)

path_mars = replaceFlags(path_mars, p)
path_mars_exp = replaceFlags(path_mars_exp, p)
path_mars_exp_files = replaceFlags(path_mars_exp_files, p)
path_mars_exp_reconstruct_files=replaceFlags(path_mars_exp_reconstruct_files,p)
path_log_file = replaceFlags(path_mars_logfile, p)


### Mars directory and subdirectory construction

if not os.path.isdir(path_mars):
    os.mkdir(path_mars)  
if not os.path.isdir(path_mars_exp):
    os.mkdir(path_mars_exp)  


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
    log_file.write('# Embryo path: %s\n'% p.PATH_EMBRYO)
    log_file.write('# Embryo name: %s\n'% p.EN)
    log_file.write('# Command line:\n'+(' '.join(sys.argv))+'\n\n\n')

### Copy of parameters file

os.system("cp -f "+parameters_file+" "+path_mars_exp )





##########################
### Mars Process Stuff ###
##########################

#Mars first time point segmentation
#
# GM: attribute raw_delay to be tested
#
fused_file = replaceTIME(path_fuse_exp_files,p.begin+p.raw_delay) # fused file
mars_file = replaceTIME(path_mars_exp_files,p.begin+p.raw_delay) # result
reconstruct_file = None

if p.mars_method == 1:
    print "Starting with fused files..."
if p.mars_method == 2:
    print "Starting with GACE files..."
    reconstruct_file = \
      replaceTIME(path_mars_exp_reconstruct_files,p.begin+p.raw_delay) # Gace


#Apply Automatic MARS Segmentation

mars_segmentation(fused_file, mars_file, 
				  p.mars_sigma1, p.mars_h_min, p.mars_sigma2, 
				  method=p.mars_method, reconstructed_image=reconstruct_file,
				  sigma_membrane=p.mars_sigma_membrane, 
				  sensitivity=p.mars_sensitivity, manual=p.mars_manual, 
				  manual_sigma=p.mars_manual_sigma, 
				  hard_thresholding=p.mars_hard_thresholding, 
				  hard_threshold=p.mars_hard_threshold, 
				  sigma_TV=p.mars_sigma_TV, sigma_LF=p.mars_sigma_LF, 
				  sample=p.mars_sample)

