
import os
import sys
import imp
import re
import time
import subprocess
import getpass
import shutil

#
#
#
#
#


class Monitoring(object):

    def __init__(self):
        self.verbose = 1
        self.debug = 0
        self.logfile = None
        self.keepTemporaryFiles = False
        self.forceResultsToBeBuilt = False

    def update_from_args(self, args):
        self.verbose = args.verbose
        self.debug = args.debug
        self.keepTemporaryFiles = args.keepTemporaryFiles
        self.forceResultsToBeBuilt = args.forceResultsToBeBuilt

    def write_parameters(self, logfilename):
        with open(logfilename, 'a') as logfile:
            logfile.write("\n")
            logfile.write('Monitoring parameters\n')
            logfile.write('- verbose is ' + str(self.verbose)+'\n')
            logfile.write('- debug is ' + str(self.debug)+'\n')
            logfile.write('- logfile is ' + str(self.logfile)+'\n')
            logfile.write('- keepTemporaryFiles is ' + str(self.keepTemporaryFiles)+'\n')
            logfile.write('- forceResultsToBeBuilt is ' + str(self.forceResultsToBeBuilt)+'\n')
            logfile.write("\n")
            return

    def print_parameters(self):
        print("")
        print('Monitoring parameters')
        print('- verbose is ' + str(self.verbose))
        print('- debug is ' + str(self.debug))
        print('- logfile is ' + str(self.logfile))
        print('- keepTemporaryFiles is ' + str(self.keepTemporaryFiles))
        print('- forceResultsToBeBuilt is ' + str(self.forceResultsToBeBuilt))
        print("")

    def copy(self, m):
        self.verbose = m.verbose
        self.debug = m.debug
        self.logfile = m.logfile
        self.keepTemporaryFiles = m.keepTemporaryFiles
        self.forceResultsToBeBuilt = m.forceResultsToBeBuilt

    @staticmethod
    def to_console(text):
        print(text)

    def to_log(self, text):
        with open(self.logfile, 'a') as logfile:
            logfile.write(text+'\n')

    def to_log_and_console(self, text, verboseness=0):
        if self.verbose >= verboseness or self.debug > 0:
            self.to_console(text)
        self.to_log(text)


#
#
#
#
#


class Experiment(object):

    def __init__(self):
        self.embryoPath = None
        self.embryoName = None
        self.firstTimePoint = -1
        self.lastTimePoint = -1
        self.deltaTimePoint = 1

    #
    #
    #
    def update_from_args(self, args):
        if args.embryoPath is None:
            return
        if not os.path.isdir(args.embryoPath):
            print ("Experiment.updateFromArgs: '" + args.embryoPath + "' is not a valid directory. Exiting.")
            sys.exit(1)
        self.embryoPath = args.embryoPath
        return

    #
    #
    #
    def update_from_file(self, parameter_file):

        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Experiment.updateFromFile: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        if hasattr(parameters, 'PATH_EMBRYO'):
            if parameters.PATH_EMBRYO is not None:
                if not os.path.isdir(parameters.PATH_EMBRYO):
                    print ("Experiment.updateFromFile: '" + parameters.PATH_EMBRYO
                           + "' is not a valid directory. Exiting.")
                    sys.exit(1)
                self.embryoPath = parameters.PATH_EMBRYO

        if hasattr(parameters, 'EN'):
            if parameters.EN is not None:
                self.embryoName = parameters.EN

        if hasattr(parameters, 'begin'):
            if parameters.begin is not None:
                self.firstTimePoint = parameters.begin

        if hasattr(parameters, 'end'):
            if parameters.end is not None:
                self.lastTimePoint = parameters.end

        if hasattr(parameters, 'delta'):
            if parameters.delta is not None:
                self.deltaTimePoint = parameters.delta
        return

    #
    #
    #
    def write_parameters(self, logfilename):
        with open(logfilename, 'a') as logfile:
            logfile.write("\n")
            logfile.write('Experiment parameters\n')
            logfile.write('- embryoPath is ' + str(self.embryoPath)+'\n')
            logfile.write('- embryoName is ' + str(self.embryoName)+'\n')
            logfile.write('- firstTimePoint is ' + str(self.firstTimePoint)+'\n')
            logfile.write('- lastTimePoint is ' + str(self.lastTimePoint)+'\n')
            logfile.write('- deltaTimePoint is ' + str(self.deltaTimePoint)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('Experiment parameters')
        print('- embryoPath is ' + str(self.embryoPath))
        print('- embryoName is ' + str(self.embryoName))
        print('- firstTimePoint is ' + str(self.firstTimePoint))
        print('- lastTimePoint is ' + str(self.lastTimePoint))
        print('- deltaTimePoint is ' + str(self.deltaTimePoint))
        print("")

#
#
#
#
#


def get_parameter_file(parameter_file):
    """
    check if the given parameter file is valid, otherwise ask for a file name
    :param parameter_file:
    :return:
    """
    if parameter_file is not None and os.path.isfile(parameter_file):
        return parameter_file
    new_parameter_file = raw_input('   Provide the parameter file: ')
    if os.path.isfile(new_parameter_file) is not False:
        print ("getParameterFile: '"+new_parameter_file+"' is not a valid file. Exiting.")
        sys.exit(1)
    return new_parameter_file

#
#
#
#
#


def write_history_information(logfile_name, experiment, parameter_file, start_time, path_to_exe):
    with open(logfile_name, 'a') as logfile:
        logfile.write("\n")
        logfile.write("# "+time.strftime("%a, %d %b %Y %H:%M:%S", start_time)+"\n")
        logfile.write("# Embryo path: '"+str(experiment.embryoPath)+"'\n")
        logfile.write("# Embryo name: '"+str(experiment.embryoName)+"'\n")
        logfile.write("# Parameter file: '" + str(parameter_file) + "'\n")
        logfile.write("# Command line: '"+" ".join(sys.argv) + "'\n")
        logfile.write("# Working directory: '"+str(os.getcwd())+"'\n")
        logfile.write("# User: '" + str(getpass.getuser()) + "'\n")
        logfile.write("# Python executable: " + sys.executable + "\n")
        logfile.write("# ASTEC version: ")
        if not os.path.exists(path_to_exe+os.path.sep+'.git'):
            logfile.write("not found\n")
        else:
            pipe = subprocess.Popen("cd "+path_to_exe+"; git describe; cd "+str(os.getcwd()),
                                    shell=True, stdout=subprocess.PIPE).stdout
            o = pipe.next()
            v = o.split('\n')
            logfile.write(str(v[0]+"\n"))
        logfile.write("# \n")
    return

#
#
#
#
#


def copy_date_stamped_file(thefile, directory, timestamp):
    d = time.strftime("%Y-%m-%d-%H:%M:%S", timestamp)
    resfile = directory+os.path.sep+re.sub(r'(\.*).py', r'\1', thefile.split(os.path.sep)[-1])+'-'+d+'.py'
    shutil.copy2(thefile, resfile)
