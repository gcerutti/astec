
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

    def update_from_args(self, args):
        self.verbose = args.verbose
        self.debug = args.debug
        self.keepTemporaryFiles = args.keepTemporaryFiles
        self.forceResultsToBeBuilt = args.forceResultsToBeBuilt

    def copy(self, m):
        if m is None:
            return
        self.verbose = m.verbose
        self.debug = m.debug
        self.logfile = m.logfile
        self.keepTemporaryFiles = m.keepTemporaryFiles
        self.forceResultsToBeBuilt = m.forceResultsToBeBuilt

    @staticmethod
    def to_console(text):
        print(text)

    def to_log(self, text):
        if self.logfile is not None:
            with open(self.logfile, 'a') as logfile:
                logfile.write(text+'\n')

    def to_log_and_console(self, text, verboseness=0):
        if self.verbose >= verboseness or self.debug > 0:
            self.to_console(text)
        self.to_log(text)


monitoring = Monitoring()


#
#
#
#
#


class ExperimentGenericSubdirectory(object):

    #
    # how to build sub-directories of /path/to/experiment
    # <main_directory>/<sub_directory_prefix><sub_directory_suffix>
    #
    #
    def __init__(self):
        self._main_directory = None
        self._sub_directory_prefix = None
        self._sub_directory_suffix = None
        self._directory = None

    #
    #
    #

    def _write_directory_list(self, logfile):
        self._set_directory()
        if self._directory is None:
            logfile.write("             None\n")
        else:
            if len(self._directory) <= 0:
                logfile.write("             Empty list\n")
            else:
                for i in range(len(self._directory)):
                    logfile.write("             #" + str(i) + ": '" + str(self.get_directory(i)) + "'\n")

    def write_parameters_in_file(self, directory_type, logfile):
        self._set_directory()
        if type(directory_type) is not str:
            logfile.write('ExperimentGenericSubdirectory.set_default: unknown arg type.\n')
        if directory_type.lower() == 'fuse':
            logfile.write("- subpath/to/fusion is \n")
            self._write_directory_list(logfile)
        elif directory_type.lower() == 'intrareg':
            logfile.write("- subpath/to/intraregistration is \n")
            self._write_directory_list(logfile)
        elif directory_type.lower() == 'mars':
            logfile.write("- subpath/to/mars is \n")
            self._write_directory_list(logfile)
        elif directory_type.lower() == 'post':
            logfile.write("- subpath/to/postcorrection is '" + str(self.get_directory()) + "'\n")
            self._write_directory_list(logfile)
        elif directory_type.lower() == 'seg':
            logfile.write("- subpath/to/segmentation is \n")
            self._write_directory_list(logfile)
        else:
            logfile.write('ExperimentGenericSubdirectory.write_parameters_in_file: unknown arg.\n')
        return

    def write_parameters(self, directory_type, logfilename):
        self._set_directory()
        with open(logfilename, 'a') as logfile:
            self.write_parameters_in_file(directory_type, logfile)
        return

    def _print_directory_list(self):
        self._set_directory()
        if self._directory is None:
            print("             None")
        else:
            if len(self._directory) <= 0:
                print("             Empty list")
            else:
                for i in range(len(self._directory)):
                    print("             #" + str(i) + ": '" + str(self.get_directory(i)) + "'")

    def print_parameters(self, directory_type):
        self._set_directory()
        if type(directory_type) is not str:
            print('ExperimentGenericSubdirectory.set_default: unknown arg type.')
        if directory_type.lower() == 'fuse':
            print("- subpath/to/fusion is")
            self._print_directory_list()
        elif directory_type.lower() == 'intrareg':
            print("- subpath/to/intraregistration is")
            self._print_directory_list()
        elif directory_type.lower() == 'mars':
            print("- subpath/to/mars is")
            self._print_directory_list()
        elif directory_type.lower() == 'post':
            print("- subpath/to/postcorrection is")
            self._print_directory_list()
        elif directory_type.lower() == 'seg':
            print("- subpath/to/segmentation is")
            self._print_directory_list()
        else:
            print('ExperimentGenericSubdirectory.print_parameters: unknown arg.\n')
        return

    #
    # set default value of sub-directory template according to type
    #

    def set_default(self, directory_type):

        if type(directory_type) is not str:
            print("ExperimentGenericSubdirectory.set_default: unknown arg type. Exiting.")
            sys.exit(1)
        if directory_type.lower() == 'fuse':
            if self._main_directory is None:
                self._main_directory = 'FUSE'
            if self._sub_directory_prefix is None:
                self._sub_directory_prefix = 'FUSE_'
            if self._sub_directory_suffix is None:
                self._sub_directory_suffix = 'RELEASE'
        elif directory_type.lower() == 'intrareg':
            if self._main_directory is None:
                self._main_directory = 'INTRAREG'
            if self._sub_directory_prefix is None:
                self._sub_directory_prefix = 'INTRAREG_'
            if self._sub_directory_suffix is None:
                self._sub_directory_suffix = 'RELEASE'
        elif directory_type.lower() == 'mars':
            if self._main_directory is None:
                self._main_directory = 'SEG'
            if self._sub_directory_prefix is None:
                self._sub_directory_prefix = 'SEG_'
            if self._sub_directory_suffix is None:
                self._sub_directory_suffix = 'RELEASE'
        elif directory_type.lower() == 'post':
            if self._main_directory is None:
                self._main_directory = 'POST'
            if self._sub_directory_prefix is None:
                self._sub_directory_prefix = 'POST_'
            if self._sub_directory_suffix is None:
                self._sub_directory_suffix = 'RELEASE'
        elif directory_type.lower() == 'seg':
            if self._main_directory is None:
                self._main_directory = 'SEG'
            if self._sub_directory_prefix is None:
                self._sub_directory_prefix = 'SEG_'
            if self._sub_directory_suffix is None:
                self._sub_directory_suffix = 'RELEASE'
        else:
            print("ExperimentGenericSubdirectory.set_default: unknown arg.")
        return

    #
    #
    #

    def _set_directory(self):
        """

        :return:
        """
        proc = "ExperimentGenericSubdirectory.set_directory"
        if self._directory is not None:
            return
        self._directory = []
        if type(self._sub_directory_suffix) == str:
            subdir = str(self._sub_directory_prefix) + str(self._sub_directory_suffix)
            self._directory.append(os.path.join(str(self._main_directory), subdir))
        elif type(self._sub_directory_suffix) == list or type(self._sub_directory_suffix) == tuple:
            for s in self._sub_directory_suffix:
                subdir = str(self._sub_directory_prefix) + str(s)
                self._directory.append(os.path.join(str(self._main_directory), subdir))
        else:
            monitoring.to_log_and_console("Warning: " + proc + ", unhandled _sub_directory_suffix type")

    def get_number_directories(self):
        #
        # self._directory is None
        # build directory list
        #
        if self._directory is None:
            self._set_directory()

        if len(self._directory) <= 0:
            return 0

        if type(self._directory) == str:
            return 1

        if type(self._directory) == list or type(self._directory) == tuple:
            return len(self._directory)

        return 0

    def get_directory(self, i=0):
        """
        return the ith directory
        :param i:
        :return:
        """
        proc = "ExperimentGenericSubdirectory.get_directory"

        #
        # self._directory is None
        # build directory list
        #
        if self._directory is None:
            self._set_directory()

        if len(self._directory) <= 0:
            return None

        #
        # self._directory is not None and not empty
        #

        if type(self._directory) == str:
            return self._directory
        # list = ['a', 'b']
        elif type(self._directory) == list:
            if i < len(self._directory):
                return self._directory[i]
            else:
                monitoring.to_log_and_console("Warning: " + proc + ", index out of range")
                return None
        # tuple = ('a', 'b')
        elif type(self._directory) == tuple:
            if i < len(self._directory):
                return self._directory[i]
            else:
                monitoring.to_log_and_console("Warning: " + proc + ", index out of range")
                return None
        else:
            monitoring.to_log_and_console("Warning: " + proc + ", unhandled _directory type")
            return None

        monitoring.to_log_and_console("Warning: " + proc + ", this should not be reached")
        return None

    #
    # read values of sub-directory template according to type
    #

    def update_from_parameters(self, directory_type, parameters):

        if type(directory_type) is not str:
            print("ExperimentGenericSubdirectory.set_default: unknown arg type. Exiting.")
            sys.exit(1)

        if directory_type.lower() == 'fuse':
            if hasattr(parameters, 'EXP_FUSE'):
                if parameters.EXP_FUSE is not None:
                    self._sub_directory_suffix = parameters.EXP_FUSE
        elif directory_type.lower() == 'intrareg':
            if hasattr(parameters, 'EXP_INTRAREG'):
                if parameters.EXP_INTRAREG is not None:
                    self._sub_directory_suffix = parameters.EXP_INTRAREG
        elif directory_type.lower() == 'mars':
            if hasattr(parameters, 'EXP_SEG'):
                if parameters.EXP_SEG is not None:
                    self._sub_directory_suffix = parameters.EXP_SEG
            if hasattr(parameters, 'EXP_MARS'):
                if parameters.EXP_MARS is not None:
                    self._sub_directory_suffix = parameters.EXP_MARS
        elif directory_type.lower() == 'post':
            if hasattr(parameters, 'EXP_POST'):
                if parameters.EXP_POST is not None:
                    self._sub_directory_suffix = parameters.EXP_POST
        elif directory_type.lower() == 'seg':
            if hasattr(parameters, 'EXP_SEG'):
                if parameters.EXP_SEG is not None:
                    self._sub_directory_suffix = parameters.EXP_SEG
        else:
            print("ExperimentGenericSubdirectory.set_default: unknown arg.")

        return

#
#
#
#
#


class Experiment(object):

    def __init__(self):

        self.embryo_path = None
        self.embryoName = None

        self.first_time_point = -1
        self.last_time_point = -1
        self.delta_time_point = 1
        self.delay_time_point = 0

        self.time_digits = 3

        #
        # stage
        # eg FUSE, MARS, MAN-CORRECTION, SEG, POST-CORRECTION, INTRAREG, PROPERTIES
        #
        self.stage = None

        #
        # sub-directories
        #
        self.fusion = ExperimentGenericSubdirectory()
        self.mars = ExperimentGenericSubdirectory()
        self.seg = ExperimentGenericSubdirectory()
        self.post = ExperimentGenericSubdirectory()
        self.intrareg = ExperimentGenericSubdirectory()

        self.fusion.set_default('fuse')
        self.mars.set_default('mars')
        self.seg.set_default('seg')
        self.post.set_default('post')
        self.intrareg.set_default('intrareg')

        #
        #
        #
        self.path_logdir = None
        self.path_history_file = None
        self.path_log_file = None

    #
    #
    #
    def write_parameters(self, logfilename):
        with open(logfilename, 'a') as logfile:
            logfile.write("\n")
            logfile.write('Experiment parameters\n')
            logfile.write('- embryo_path is ' + str(self.embryo_path)+'\n')
            logfile.write('- embryoName is ' + str(self.embryoName)+'\n')
            logfile.write('- first_time_point is ' + str(self.first_time_point)+'\n')
            logfile.write('- last_time_point is ' + str(self.last_time_point)+'\n')
            logfile.write('- delta_time_point is ' + str(self.delta_time_point)+'\n')
            logfile.write('- delay_time_point is ' + str(self.delay_time_point)+'\n')
            logfile.write('- time_digits is ' + str(self.time_digits) + '\n')
            self.fusion.write_parameters_in_file('FUSE', logfile)
            self.mars.write_parameters_in_file('MARS', logfile)
            self.seg.write_parameters_in_file('SEG', logfile)
            self.post.write_parameters_in_file('POST', logfile)
            self.intrareg.write_parameters_in_file('INTRAREG', logfile)
            logfile.write('- path_logdir is ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file is ' + str(self.path_history_file) + '\n')
            logfile.write('- path_log_file is ' + str(self.path_log_file) + '\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('Experiment parameters')
        print('- embryo_path is ' + str(self.embryo_path))
        print('- embryoName is ' + str(self.embryoName))
        print('- first_time_point is ' + str(self.first_time_point))
        print('- last_time_point is ' + str(self.last_time_point))
        print('- delta_time_point is ' + str(self.delta_time_point))
        print('- delay_time_point is ' + str(self.delay_time_point))
        print('- time_digits is ' + str(self.time_digits))
        self.fusion.print_parameters('FUSE')
        self.mars.print_parameters('MARS')
        self.seg.print_parameters('SEG')
        self.post.print_parameters('POST')
        self.intrareg.print_parameters('INTRAREG')
        print('- path_logdir is ' + str(self.path_logdir))
        print('- path_history_file is ' + str(self.path_history_file))
        print('- path_log_file is ' + str(self.path_log_file))
        print("")

    #
    #
    #
    def get_time_format(self, time_digits=3):
        """

        :param time_digits:
        :return:
        """
        form = "%0" + str(time_digits) + "d"
        return form

    def get_time_index(self, index, time_digits=3):
        """

        :param index:
        :param time_digits:
        :return:
        """
        ind = '{:0{width}d}'.format(index, width=time_digits)
        return ind

    def get_image_suffix(self, directory_type, subdirectory_type='fuse'):
        """

        :param directory_type:
        :param subdirectory_type:
        :return:
        """
        if directory_type.lower() == 'fuse':
            return self.embryoName + '_fuse'
        elif directory_type.lower() == 'intrareg':
            if subdirectory_type.lower() == 'fuse':
                return self.embryoName + '_intrareg_fuse'
            elif subdirectory_type.lower() == 'post':
                return self.embryoName + '_intrareg_post'
            elif subdirectory_type.lower() == 'seg':
                return self.embryoName + '_intrareg_seg'
            else:
                print("Experiment.get_image_name: unknown arg #2.")
                return None
        elif directory_type.lower() == 'post':
            return self.embryoName + '_post'
        elif directory_type.lower() == 'seg':
            return self.embryoName + '_seg'
        else:
            print("Experiment.get_image_name: unknown arg.")
            return None

    def get_image_name(self, index, directory_type, subdirectory_type='fuse'):
        """

        :param index:
        :param directory_type:
        :param subdirectory_type:
        :return:
        """
        ind = self.get_time_index(index)
        suffix = self.get_image_suffix(directory_type, subdirectory_type)
        if suffix is None:
            return None
        return suffix + '_t' + str(ind)

    def get_image_format(self, directory_type, subdirectory_type='fuse'):
        """

        :param directory_type:
        :param subdirectory_type:
        :return:
        """
        suffix = self.get_image_suffix(directory_type, subdirectory_type)
        if suffix is None:
            return None
        return suffix + '_t' + self.get_time_format()

    def get_movie_name(self, first, last, directory_type, subdirectory_type='fuse'):
        """

        :param first:
        :param last:
        :param directory_type:
        :param subdirectory_type:
        :return:
        """
        f = self.get_time_index(first)
        l = self.get_time_index(last)
        suffix = self.get_image_suffix(directory_type, subdirectory_type)
        if suffix is None:
            return None
        return suffix + '_t' + str(f) + '-' + str(l)

    #
    #
    #
    def update_from_args(self, args):
        """

        :param args:
        :return:
        """
        if args.embryo_path is None:
            return
        if not os.path.isdir(args.embryo_path):
            print ("Experiment.updateFromArgs: '" + args.embryo_path + "' is not a valid directory. Exiting.")
            sys.exit(1)
        self.embryo_path = args.embryo_path
        return

    #
    #
    #
    def update_from_file(self, parameter_file):
        """

        :param parameter_file:
        :return:
        """

        proc = "Experiment.updateFromFile"
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print proc + ": '" + parameter_file + "' is not a valid file."
            print "\t Exiting."
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        if hasattr(parameters, 'PATH_EMBRYO'):
            if parameters.PATH_EMBRYO is not None:
                if not os.path.isdir(parameters.PATH_EMBRYO):
                    print ("Experiment.updateFromFile: '" + parameters.PATH_EMBRYO
                           + "' is not a valid directory. Exiting.")
                    sys.exit(1)
                self.embryo_path = parameters.PATH_EMBRYO
            else:
                self.embryo_path = os.getcwd()
        else:
            self.embryo_path = os.getcwd()

        if hasattr(parameters, 'EN'):
            if parameters.EN is not None:
                self.embryoName = parameters.EN
            else:
                embryo_path = self.embryo_path
                self.embryoName = embryo_path.split(os.path.sep)[-1]
        else:
            embryo_path = self.embryo_path
            self.embryoName = embryo_path.split(os.path.sep)[-1]

        if hasattr(parameters, 'begin'):
            if parameters.begin is not None:
                self.first_time_point = parameters.begin
            else:
                print proc + ": it is mandatory to specify the first time point"
                print "\t Exiting."
                sys.exit(1)

        if hasattr(parameters, 'end'):
            if parameters.end is not None:
                self.last_time_point = parameters.end
            else:
                print proc + ": it is mandatory to specify the last time point"
                print "\t Exiting."
                sys.exit(1)

        if hasattr(parameters, 'delta'):
            if parameters.delta is not None:
                self.delta_time_point = parameters.delta

        if hasattr(parameters, 'raw_delay'):
            if parameters.raw_delay is not None:
                self.delay_time_point = parameters.raw_delay

        self.fusion.update_from_parameters('FUSE', parameters)
        self.mars.update_from_parameters('MARS', parameters)
        self.seg.update_from_parameters('SEG', parameters)
        self.post.update_from_parameters('POST', parameters)
        self.intrareg.update_from_parameters('INTRAREG', parameters)

        return

    #
    #
    #
    def update_from_stage(self, stage, executable, timestamp):
        """

        :param stage:
        :param executable:
        :param timestamp:
        :return:
        """
        if stage.lower() == 'fuse':
            self.path_logdir = os.path.join(self.fusion.get_directory(), "LOGS")
        elif stage.lower() == 'mars':
            self.path_logdir = os.path.join(self.mars.get_directory(), "LOGS")
        elif stage.lower() == 'man-correction' or stage.lower() == 'manual-correction':
            self.path_logdir = os.path.join(self.seg.get_directory(), "LOGS")
        elif stage.lower() == 'seg':
            self.path_logdir = os.path.join(self.seg.get_directory(), "LOGS")
        elif stage.lower() == 'post' or stage.lower() == 'post-correction':
            self.path_logdir = os.path.join(self.post.get_directory(), "LOGS")
        elif stage.lower() == 'intrareg' or stage.lower() == 'intraregistration':
            self.path_logdir = os.path.join(self.intrareg.get_directory(), "LOGS")
        elif stage.lower() == 'properties':
            self.path_logdir = os.path.join(self.intrareg.get_directory(), "LOGS")
        else:
            print("Experiment.update_from_stage: unknown stage.")

        if executable is None:
            local_executable = 'unknown'
        else:
            local_executable = os.path.basename(executable)
            if local_executable[-3:] == '.py':
                local_executable = local_executable[:-3]

        self.path_history_file = os.path.join(self.path_logdir, str(local_executable)+'-history.log')
        if timestamp is None:
            d = time.strftime("%Y-%m-%d-%H:%M:%S", time.localtime())
        else:
            d = time.strftime("%Y-%m-%d-%H:%M:%S", timestamp)
        self.path_log_file = os.path.join(self.path_logdir, str(local_executable) + '-' + str(d) + '.log')

        return


########################################################################################
#
#
#
########################################################################################


def _fullname(prefix, desc):
    if prefix is not None:
        return prefix + desc
    else:
        return desc


def _fulldesc(prefix, desc):
    return '- ' + _fullname(prefix, desc) + ' = '


class RegistrationParameters(object):

    def __init__(self):
        #
        # prefix is for naming the parameters
        #
        self.prefix = None
        #
        #
        #
        self.compute_registration = True

        #
        # parameters
        #
        self.transformation_type = 'affine'
        self.transformation_estimation_type = 'wlts'
        self.lts_fraction = 0.55
        self.pyramid_highest_level = 6
        self.pyramid_lowest_level = 3
        self.normalization = True

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('RegistrationParameters\n')
            logfile.write(_fulldesc(None, 'prefix')+str(self.prefix)+'\n')
            logfile.write(_fulldesc(self.prefix, 'compute_registration') + str(self.compute_registration) + '\n')

            logfile.write(_fulldesc(self.prefix, 'transformation_type')+str(self.transformation_type)+'\n')
            logfile.write(_fulldesc(self.prefix, 'transformation_estimation_type')
                          + str(self.transformation_estimation_type)+'\n')
            logfile.write(_fulldesc(self.prefix, 'lts_fraction')+str(self.lts_fraction)+'\n')
            logfile.write(_fulldesc(self.prefix, 'pyramid_highest_level')+str(self.pyramid_highest_level)+'\n')
            logfile.write(_fulldesc(self.prefix, 'pyramid_lowest_level')+str(self.pyramid_lowest_level)+'\n')
            logfile.write(_fulldesc(self.prefix, 'normalization')+str(self.normalization)+'\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('RegistrationParameters')

        print(_fulldesc(None, 'prefix') + str(self.prefix))
        print(_fulldesc(self.prefix, 'compute_registration') + str(self.compute_registration))

        print(_fulldesc(self.prefix, 'transformation_type') + str(self.transformation_type))
        print(_fulldesc(self.prefix, 'transformation_estimation_type') + str(self.transformation_estimation_type))
        print(_fulldesc(self.prefix, 'lts_fraction') + str(self.lts_fraction))
        print(_fulldesc(self.prefix, 'pyramid_highest_level') + str(self.pyramid_highest_level))
        print(_fulldesc(self.prefix, 'pyramid_lowest_level') + str(self.pyramid_lowest_level))
        print(_fulldesc(self.prefix, 'normalization') + str(self.normalization))

        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        if hasattr(parameters, _fullname(self.prefix, 'compute_registration')):
            if getattr(parameters, _fullname(self.prefix, 'compute_registration'), 'None') is not None:
                self.compute_registration = getattr(parameters, _fullname(self.prefix, 'compute_registration'))

        if hasattr(parameters, _fullname(self.prefix, 'transformation_type')):
            if getattr(parameters, _fullname(self.prefix, 'transformation_type'), 'None') is not None:
                self.transformation_type = getattr(parameters, _fullname(self.prefix, 'transformation_type'))

        if hasattr(parameters, _fullname(self.prefix, 'transformation_estimation_type')):
            if getattr(parameters, _fullname(self.prefix, 'transformation_estimation_type'), 'None') is not None:
                self.transformation_estimation_type = getattr(parameters,
                                                              _fullname(self.prefix, 'transformation_estimation_type'))

        if hasattr(parameters, _fullname(self.prefix, 'lts_fraction')):
            if getattr(parameters, _fullname(self.prefix, 'lts_fraction'), 'None') is not None:
                self.lts_fraction = getattr(parameters, _fullname(self.prefix, 'lts_fraction'))

        if hasattr(parameters, _fullname(self.prefix, 'pyramid_highest_level')):
            if getattr(parameters, _fullname(self.prefix, 'pyramid_highest_level'), 'None') is not None:
                self.pyramid_highest_level = getattr(parameters, _fullname(self.prefix, 'pyramid_highest_level'))
        if hasattr(parameters, _fullname(self.prefix, 'pyramid_lowest_level')):
            if getattr(parameters, _fullname(self.prefix, 'pyramid_lowest_level'), 'None') is not None:
                self.pyramid_lowest_level = getattr(parameters, _fullname(self.prefix, 'pyramid_lowest_level'))

        if hasattr(parameters, _fullname(self.prefix, 'normalization')):
            if getattr(parameters, _fullname(self.prefix, 'normalization'), 'None') is not None:
                self.normalization = getattr(parameters, _fullname(self.prefix, 'normalization'))


########################################################################################
#
#
#
########################################################################################

def get_parameter_file(parameter_file):
    """
    check if the given parameter file is valid, otherwise ask for a file name
    :param parameter_file: the parameter file name to be tested
    :return: the parameter file name
    """
    if parameter_file is not None and os.path.isfile(parameter_file):
        return parameter_file
    new_parameter_file = raw_input('   Provide the parameter file: ')
    if os.path.isfile(new_parameter_file) is not False:
        print ("getParameterFile: '"+new_parameter_file+"' is not a valid file. Exiting.")
        sys.exit(1)
    return new_parameter_file


########################################################################################
#
#
#
########################################################################################

def _write_git_information(path, logfile, desc):
    """

    :param path:
    :param logfile:
    :param desc:
    :return:
    """
    gitremote = 'git remote get-url origin'
    gitdescribe = 'git describe'
    gitlog = 'git log -n 1 --format="Commit (tag, date, ref): %H -- %cD -- %D"'

    logfile.write(str(desc) + " path: ")

    pipe = subprocess.Popen("cd " + path + "; " + "pwd" + "; cd " + str(os.getcwd()),
                            shell=True, stdout=subprocess.PIPE).stdout
    o = pipe.next()
    v = o.split('\n')
    logfile.write(str(v[0] + "\n"))

    logfile.write(str(desc) + " repository: ")

    # if not os.path.exists(path + os.path.sep + '.git'):
    #     logfile.write("not found\n")
    #     return

    pipe = subprocess.Popen("cd " + path + "; " + gitremote + "; cd " + str(os.getcwd()),
                            shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (stdoutData, stderrData) = pipe.communicate()

    if len(stderrData) > 6:
        if stderrData[0:6] == "fatal:":
            logfile.write("no a git repository\n")
        else:
            logfile.write("no a git repository?\n")
    elif len(stderrData) > 0:
        logfile.write("no a git repository?!\n")
    else:
        v = stdoutData.split('\n')
        logfile.write(str(v[0] + "\n"))

        logfile.write(str(desc) + " version: ")
        pipe = subprocess.Popen("cd " + path + "; " + gitdescribe + "; cd " + str(os.getcwd()),
                                shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (stdoutData, stderrData) = pipe.communicate()
        if len(stderrData) > 0:
            logfile.write("no version found\n")
        else:
            v = stdoutData.split('\n')
            logfile.write(str(v[0] + "\n"))

        pipe = subprocess.Popen("cd " + path + "; " + gitlog + "; cd " + str(os.getcwd()),
                                shell=True, stdout=subprocess.PIPE).stdout
        o = pipe.next()
        v = o.split('\n')
        logfile.write("#\t" + str(v[0] + "\n"))

    return


def write_history_information(logfile_name,
                              experiment=None,
                              parameter_file=None,
                              start_time=None,
                              path_to_exe=None,
                              path_to_vt=None):
    """
    Write history information
    :param logfile_name: file where to write information
    :param experiment: experiment values (class)
    :param parameter_file: name of the parameter file
    :param start_time: start time
    :param path_to_exe: path to Astec repository, used to get git version
    :param path_to_vt: path to VT build repository, used to get git version
    :return:
    """
    with open(logfile_name, 'a') as logfile:
        logfile.write("\n")
        if start_time is not None:
            logfile.write("# "+time.strftime("%a, %d %b %Y %H:%M:%S", start_time)+"\n")
        if experiment is not None:
            logfile.write("# Embryo path: '"+str(experiment.embryo_path)+"'\n")
            logfile.write("# Embryo name: '"+str(experiment.embryoName)+"'\n")
        if parameter_file is not None:
            logfile.write("# Parameter file: '" + str(parameter_file) + "'\n")
        logfile.write("# Command line: '"+" ".join(sys.argv) + "'\n")
        logfile.write("# Working directory: '"+str(os.getcwd())+"'\n")
        logfile.write("# User: '" + str(getpass.getuser()) + "'\n")
        logfile.write("# Python executable: " + sys.executable + "\n")
        if path_to_exe is not None:
            _write_git_information(path_to_exe, logfile, "# ASTEC")
        if path_to_vt is not None:
            _write_git_information(path_to_vt, logfile, "# VT")
        logfile.write("# \n")
    return

########################################################################################
#
#
#
########################################################################################


def copy_date_stamped_file(thefile, directory, timestamp):
    """
    Copy a file to the designated directory while adding a time stamp to its name
    :param thefile:
    :param directory:
    :param timestamp:
    :return:
    """
    d = time.strftime("%Y-%m-%d-%H:%M:%S", timestamp)
    if len(thefile.split('.')) > 1:
        ext = thefile.split('.')[-1]
        filename = re.sub(r'(\.*).' + ext, r'\1', thefile.split(os.path.sep)[-1]) + '-' + d + '.' + ext
    else:
        filename = thefile.split(os.path.sep)[-1] + '-' + d
    resfile = os.path.join(directory, filename)
    shutil.copy2(thefile, resfile)


########################################################################################
#
#
#
########################################################################################


def read_lut(filename):
    """
    Return a dictionnary of integer key-to-key correspondances
    :param filename:
    :return:
    """
    # proc = 'read_lut'
    lut = {}

    if not os.path.isfile(filename):
        # monitoring.to_log_and_console(proc + ": file '" + str(filename) + "' does not exists", 0)
        return lut

    f = open(filename)
    for line in f:
        li = line.strip()
        if li.startswith('#'):
            continue
        info = li.split()
        if len(info) == 2:
            # if not lut.has_key(int(info[0])):
            #   lut[int(info[0])] = None
            lut[int(info[0])] = int(info[1])
    f.close()

    return lut


########################################################################################
#
# image file utilities
#
########################################################################################


recognized_extensions = ['.zip', '.h5', '.tif', '.tiff', '.TIF', '.TIFF', '.inr', '.inr.gz', '.mha', '.mha.gz']


def get_extension(filename):
    """ Return the file extension. Must be in the set of recognized extensions.
    :param filename:
    :return: None in case of unrecognized extension,
             else the recognized extension (begins with '.')
    """
    for e in recognized_extensions:
        if len(filename) < len(e):
            continue
        if filename[len(filename)-len(e):len(filename)] == e:
            return e
    return None

#
#
#
#
#


def add_suffix(filename, suffix, new_dirname=None, new_extension=None):
    """
    Add a suffix to a filename (ie before the extension)
    :param filename:
    :param suffix: suffix to be added
    :param new_dirname: change the directory name of the file
    :param new_extension: change the extension of the file
    :return: the transformed file name
    """
    proc = 'add_suffix'
    if filename is None:
        print(proc + ": was called with '" + str(filename) + "'")
        return
    b = os.path.basename(filename)
    d = os.path.dirname(filename)
    e = get_extension(b)
    if e is None:
        print(proc + ": file extension of '"+str(filename)+"' was not recognized")
        print("\t Exiting")
#        monitoring.to_log_and_console(proc + ": file extension of '"+str(filename)+"' was not recognized", 0)
#        monitoring.to_log_and_console("\t Exiting", 0)
        sys.exit(1)
    new_basename = b[0:len(b)-len(e)]
    new_basename += suffix
    if new_extension is None:
        new_basename += e
    else:
        if new_extension[0] == '.':
            new_basename += new_extension
        else:
            new_basename += '.' + new_extension
    if new_dirname is None:
        res_name = os.path.join(d, new_basename)
    else:
        res_name = os.path.join(new_dirname, new_basename)
    return res_name


#
#
#
#
#

def find_file(data_path, file_prefix, local_monitoring=None, verbose=True):
    """
    find a file in a directory with a given prefix. The suffix is unknown

    :param data_path:
    :param file_prefix:
    :param local_monitoring:
    :param verbose:
    :return:
    """
    proc = "find_file"

    if not os.path.isdir(data_path):
        if local_monitoring is not None:
            local_monitoring.to_log_and_console("Error:")
            local_monitoring.to_log_and_console(proc + ": '" + str(data_path) + "' is not a valid directory ?!")
            local_monitoring.to_log_and_console("\t Exiting.")
        else:
            print(proc + ": '" + str(data_path) + "' is not a valid directory ?!")
            print("\t Exiting.")
        sys.exit(1)

    if file_prefix is None:
        print(proc + ": file prefix was 'None'?!")
        print("\t Exiting.")
        sys.exit(1)

    #
    # if there is any extension, remove if from the file_prefix length
    # recall that the '.' is part of the extension
    #
    extension = get_extension(file_prefix)
    if extension is not None:
        length_file_prefix = len(file_prefix) - len(extension)
    else:
        length_file_prefix = len(file_prefix)

    #
    # get all file names beginning by the given prefix followed by '.'
    #
    file_names = []
    for f in os.listdir(data_path):
        if len(f) <= length_file_prefix:
            pass
        if f[0:length_file_prefix] == file_prefix[0:length_file_prefix] and f[length_file_prefix] == '.':
            file_names.append(f)

    if len(file_names) == 0:
        if local_monitoring is not None:
            local_monitoring.to_log_and_console(proc + ": no image with name '" + str(file_prefix)
                                                + "' was found in '" + str(data_path) + "'", 4)
        elif verbose is True:
            print(proc + ": no image with name '" + str(file_prefix) + "' was found in '" + str(data_path) + "'")
        return None

    if len(file_names) > 1:
        if local_monitoring is not None:
            local_monitoring.to_log_and_console("\t " + proc + ": warning")
            local_monitoring.to_log_and_console("\t several images with name '" + str(file_prefix) + "' were found in")
            local_monitoring.to_log_and_console("\t    '" + str(data_path) + "'")
            local_monitoring.to_log_and_console("\t    -> "+str(file_names))
            local_monitoring.to_log_and_console("\t returned file is '" + str(file_names[0]) + "'")
        else:
            print(proc + ": several images with name '"
                  + str(file_prefix) + "' were found in '" + str(data_path) + "'")
            print("\t "+str(file_names))
        # return None

    return file_names[0]


#
#
#
#
#

def get_file_suffix(experiment, data_path, file_format, flag_time=None):
    """

    :param experiment:
    :param data_path:
    :param file_format:
    :param flag_time:
    :return:
    """

    proc = "get_file_suffix"

    if not os.path.isdir(data_path):
        monitoring.to_log_and_console(proc + ": weird, data path '" + str(data_path) + "' is not a valid directory", 0)
        return None

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    suffixes = {}
    nimages = 0
    nfiles = 0

    if flag_time is not None:
        flag = flag_time
    else:
        flag = "$TIME"

    #
    # get and count suffixes for images
    #
    for current_time in range(first_time_point + experiment.delay_time_point + experiment.delta_time_point,
                              last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        time_point = '{:0{width}d}'.format(current_time, width=experiment.time_digits)
        file_prefix = file_format.replace(flag, time_point)

        for f in os.listdir(data_path):
            if len(f) <= len(file_prefix):
                pass
            if f[0:len(file_prefix)] == file_prefix and f[len(file_prefix)] == '.':
                suffix = f[len(file_prefix) + 1:len(f)]
                suffixes[suffix] = suffixes.get(suffix, 0) + 1
                nfiles += 1

        nimages += 1

    for s, n in suffixes.items():
        if n == nimages:
            return s

    if nfiles < nimages:
        monitoring.to_log_and_console(proc + ": weird, not enough images '" + str(file_format)
                                      + "' were found in '" + str(data_path) + "'", 0)
        monitoring.to_log_and_console("\t Exiting.", 0)
        exit(1)

    monitoring.to_log_and_console(proc + ": no common suffix for '" + str(file_format)
                                  + "' was found in '" + str(data_path) + "'", 2)
    monitoring.to_log_and_console("\t time point range was ["+str(first_time_point)+", "+str(last_time_point)+"]")
    return None
