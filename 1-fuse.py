#!/usr/bin/python2.7

import os, sys, imp

assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC"))
sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))

from optparse import OptionParser

from nomenclature import *
from FUSION import read_raw_data,fusion_process
from lineage import timeNamed,timesNamed


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

print "Raw data will be searched in directory %s"%replaceFlags(path_rawdata,p)

path_angle1 = replaceFlags(path_rawdata_angle1, p)
path_angle2 = replaceFlags(path_rawdata_angle2, p)
path_angle3 = replaceFlags(path_rawdata_angle3, p)
path_angle4 = replaceFlags(path_rawdata_angle4, p)

path_fuse = replaceFlags(path_fuse, p)
path_fuse_exp = replaceFlags(path_fuse_exp, p)

path_fuse_exp_files = replaceFlags(path_fuse_exp_files, p)


path_log_file = replaceFlags(path_fuse_logfile, p)

### Fusion directory and subdirectory construction

if not os.path.isdir(path_fuse):
    os.mkdir(path_fuse)  
if not os.path.isdir(path_fuse_exp):
    os.mkdir(path_fuse_exp)  


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

os.system("cp -f "+parameters_file+" "+path_fuse_exp )

############################
### Fusion Process Stuff ###
############################

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
    temporary_path=os.path.join(path_fuse_exp,"TEMP_$TIME") #Temporary Path
    #PROCESS THE FUSION
    for time in range(begin, end+1, p.delta): # Interation on time steps

        fused_file=replaceTIME(path_fuse_exp_files,time+p.raw_delay)

        if not os.path.isfile(fused_file):
            time_angles_files=[timeNamed(angle_file + ext_im,time) \
                               for angle_file in angles_files]
            temporary_time_path=timeNamed(temporary_path,time) # Temporary Path
                                                               # for this t-p
            print temporary_time_path
            time_process=fusion_process(time_angles_files,
                       fused_file,  
                       temporary_time_path,
                       p.raw_ori, p.raw_resolution, p.target_resolution, 
                       p.raw_delay, ext_im1, 
                       mirrors = p.raw_mirrors, 
                       margin_x_0=p.fusion_margin_x_0, 
                       margin_x_1=p.fusion_margin_x_1, 
                       margin_y_0=p.fusion_margin_y_0, 
                       margin_y_1=p.fusion_margin_y_1, 
                       crop=p.fusion_crop )
            print "Time point " + str(time) + " takes " + str(time_process) +\
                  " to compute\n\n\n"
            os.system("rm -rf "+temporary_time_path) # Cleaning temporary files
