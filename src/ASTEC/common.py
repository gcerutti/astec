
import os
import sys
import imp
import re
import time
import subprocess
import getpass
import shutil

import CommunFunctions.cpp_wrapping as cpp_wrapping

#
#
#
#
#


def _timestamp_to_str(timestamp=None):
    """
    build a string from a time stamp, e.g. time.localtime()
    :param timestamp:
    :return:
    """
    if timestamp is None:
        timestamp = time.localtime()
    d = time.strftime("%Y-%m-%d-%H-%M-%S", timestamp)
    return d

##################################################
#
# Monitoring processing
#
##################################################


class Monitoring(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):
        self.verbose = 1
        self.debug = 0
        self.log_filename = None
        self.keepTemporaryFiles = False
        self.forceResultsToBeBuilt = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('Monitoring parameters')
        print('- verbose is ' + str(self.verbose))
        print('- debug is ' + str(self.debug))
        print('- log_filename is ' + str(self.log_filename))
        print('- keepTemporaryFiles is ' + str(self.keepTemporaryFiles))
        print('- forceResultsToBeBuilt is ' + str(self.forceResultsToBeBuilt))
        print("")

    def write_parameters(self):
        if self.log_filename is not None:
            with open(self.log_filename, 'a') as logfile:
                logfile.write("\n")
                logfile.write('Monitoring parameters\n')
                logfile.write('- verbose is ' + str(self.verbose)+'\n')
                logfile.write('- debug is ' + str(self.debug)+'\n')
                logfile.write('- log_filename is ' + str(self.log_filename)+'\n')
                logfile.write('- keepTemporaryFiles is ' + str(self.keepTemporaryFiles)+'\n')
                logfile.write('- forceResultsToBeBuilt is ' + str(self.forceResultsToBeBuilt)+'\n')
                logfile.write("\n")
        return

    def update_execution_time(self, start_time=None, end_time=None):
        if self.log_filename is not None:
            with open(self.log_filename, 'a') as logfile:
                logfile.write('# Total execution time = ' + str(time.mktime(end_time) - time.mktime(start_time))
                              + ' sec\n')
                logfile.write("\n\n")
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_args(self, args):
        self.verbose = args.verbose
        self.debug = args.debug
        self.keepTemporaryFiles = args.keepTemporaryFiles
        self.forceResultsToBeBuilt = args.forceResultsToBeBuilt

    ############################################################
    #
    # setters
    #
    ############################################################

    def set_log_filename(self, experiment, cli_name=None, timestamp=None):
        """
        set the full logfile name. The log directory name is built
        :param experiment:
        :param cli_name:
        :param timestamp:
        :return:
        """

        if not isinstance(experiment, Experiment):
            print("Monitoring.set_log_filename" + ": weird type ('" + str(type(experiment)) + "') for 'experiment'.")
            print("\t Exiting.")
            sys.exit(1)

        if cli_name is None:
            local_executable = 'unknown'
        else:
            local_executable = os.path.basename(cli_name)
            if local_executable[-3:] == '.py':
                local_executable = local_executable[:-3]
        log_filename = local_executable + '-' + _timestamp_to_str(timestamp) + '.log'
        log_dirname = experiment.get_log_dirname()
        if not os.path.isdir(log_dirname):
            os.makedirs(log_dirname)
        self.log_filename = os.path.join(log_dirname, log_filename)
        return

    ############################################################
    #
    # misc
    #
    ############################################################

    def copy(self, m):
        """
        make a copy from an other object
        :param m:
        :return:
        """
        if m is None:
            return
        self.verbose = m.verbose
        self.debug = m.debug
        self.log_filename = m.log_filename
        self.keepTemporaryFiles = m.keepTemporaryFiles
        self.forceResultsToBeBuilt = m.forceResultsToBeBuilt

    @staticmethod
    def to_console(text):
        print(text)

    def to_log(self, text):
        if self.log_filename is not None:
            with open(self.log_filename, 'a') as logfile:
                logfile.write(text+'\n')

    def to_log_and_console(self, text, verboseness=0):
        if self.verbose >= verboseness or self.debug > 0:
            self.to_console(text)
        self.to_log(text)


monitoring = Monitoring()


##################################################
#
# sub-directories management utilities
#
##################################################


#
# utilities
#

def _get_directory(directories, i=0):
    """
    return the ith directory
    :param directories:
    :param i:
    :return:
    """
    proc = "_get_directory"

    if len(directories) <= 0:
        return None

    #
    # directories is not None and not empty
    #

    if type(directories) == str:
        return directories
    # list = ['a', 'b']
    elif type(directories) == list:
        if i < len(directories):
            return directories[i]
        else:
            monitoring.to_log_and_console("Warning: " + proc + ", index out of range")
            return None
    # tuple = ('a', 'b')
    elif type(directories) == tuple:
        if i < len(directories):
            return directories[i]
        else:
            monitoring.to_log_and_console("Warning: " + proc + ", index out of range")
            return None

    monitoring.to_log_and_console("Warning: " + proc + ", unhandled _directory type")
    return None

#
#
#


def _print_directory_list(directories):
    if directories is None:
        print("       None")
    else:
        if len(directories) <= 0:
            print("       Empty list")
        else:
            for i in range(len(directories)):
                print("       #" + str(i) + ": '" + str(_get_directory(directories, i)) + "'")
    return

#
#
#


def _write_directory_list_in_file(directories, logfile):
    if directories is None:
        logfile.write("       None\n")
    else:
        if len(directories) <= 0:
            logfile.write("       Empty list\n")
        else:
            for i in range(len(directories)):
                logfile.write("      #" + str(i) + ": '" + str(_get_directory(directories, i)) + "'\n")
    return

#
#
#


def _make_directory_list(directories, min_index=-1):
    if directories is None:
        return
    if len(directories) <= 0:
        return
    if min_index <= 0:
        length = len(directories)
    else:
        length = min(min_index, len(directories))
    for i in range(length):
        d = _get_directory(directories, i)
        if not os.path.isdir(d):
            os.makedirs(d)
    return


def _rmtree_directory_list(directories, min_index=-1):
    if directories is None:
        return
    if len(directories) <= 0:
        return
    if min_index <= 0:
        length = len(directories)
    else:
        length = min(min_index, len(directories))
    for i in range(length):
        d = _get_directory(directories, i)
        if os.path.isdir(d):
            shutil.rmtree(d)
    return

#
#
#


##################################################
#
# rawdata directory management
#
##################################################

class RawdataChannel(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, c=None):

        self._parent_directory = None
        self._main_directory = 'RAWDATA'

        #
        # raw data directories
        #
        if c == 0:
            self.angle0_sub_directory = os.path.join('LC', 'Stack0000')
            self.angle1_sub_directory = os.path.join('RC', 'Stack0000')
            self.angle2_sub_directory = os.path.join('LC', 'Stack0001')
            self.angle3_sub_directory = os.path.join('RC', 'Stack0001')
        else:
            self.angle0_sub_directory = None
            self.angle1_sub_directory = None
            self.angle2_sub_directory = None
            self.angle3_sub_directory = None

        #
        # raw data file names
        # assumed to be the same for all channels
        #
        self._angle0_file_prefix = "Cam_Left_00"
        self._angle1_file_prefix = "Cam_Right_00"
        self._angle2_file_prefix = "Cam_Left_00"
        self._angle3_file_prefix = "Cam_Right_00"

        #
        # temporary_paths
        #
        self.tmp_directory = list()

        #
        #
        #
        self.fusion_weighting = 'guignard-weighting'

        #
        #
        #
        self._time_digits_for_acquisition = 3

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self, desc=None):
        if desc is not None:
            print("  - " + str(desc) + " =")
        if self.is_empty():
            print('    RawdataChannel empty')
            return
        else:
            print('    RawdataChannel')
        print("    - _parent_directory = " + str(self._parent_directory))
        print("    - _main_directory = " + str(self._main_directory))

        print('    - angle0_sub_directory = ' + str(self.angle0_sub_directory))
        print('    - angle1_sub_directory = ' + str(self.angle1_sub_directory))
        print('    - angle2_sub_directory = ' + str(self.angle2_sub_directory))
        print('    - angle3_sub_directory = ' + str(self.angle3_sub_directory))

        print('    - _angle0_file_prefix = ' + str(self._angle0_file_prefix))
        print('    - _angle1_file_prefix = ' + str(self._angle1_file_prefix))
        print('    - _angle2_file_prefix = ' + str(self._angle2_file_prefix))
        print('    - _angle3_file_prefix = ' + str(self._angle3_file_prefix))

        for j in range(0, len(self.tmp_directory)):
            print('    - tmp_directory #' + str(j) + ' = ' + str(self.tmp_directory[j]))

        print('    - fusion_weighting = ' + str(self.fusion_weighting))
        print('    - _time_digits_for_acquisition = ' + str(self._time_digits_for_acquisition))
        return

    def write_parameters_in_file(self, logfile, desc=None):
        if desc is not None:
            logfile.write("  - " + str(desc) + " =\n")
        logfile.write('    RawdataChannel')
        if self.is_empty():
            logfile.write(' empty\n')
            return
        else:
            logfile.write('\n')
        logfile.write("    - _parent_directory = " + str(self._parent_directory) + "\n")
        logfile.write("    - _main_directory = " + str(self._main_directory) + "\n")

        logfile.write('    - angle0_sub_directory = ' + str(self.angle0_sub_directory)+'\n')
        logfile.write('    - angle1_sub_directory = ' + str(self.angle1_sub_directory)+'\n')
        logfile.write('    - angle2_sub_directory = ' + str(self.angle2_sub_directory)+'\n')
        logfile.write('    - angle3_sub_directory = ' + str(self.angle3_sub_directory)+'\n')

        logfile.write('    - _angle0_file_prefix = ' + str(self._angle0_file_prefix) + '\n')
        logfile.write('    - _angle1_file_prefix = ' + str(self._angle1_file_prefix) + '\n')
        logfile.write('    - _angle2_file_prefix = ' + str(self._angle2_file_prefix) + '\n')
        logfile.write('    - _angle3_file_prefix = ' + str(self._angle3_file_prefix) + '\n')

        for j in range(0, len(self.tmp_directory)):
            logfile.write('    - tmp_directory #' + str(j) + ' = ' + str(self.tmp_directory[j]) + '\n')

        logfile.write('    - fusion_weighting = ' + str(self.fusion_weighting) + '\n')
        logfile.write('    - _time_digits_for_acquisition = ' + str(self._time_digits_for_acquisition) + '\n')
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, channel_id, parameters):

        if hasattr(parameters, 'DIR_RAWDATA'):
            self._main_directory = getattr(parameters, 'DIR_RAWDATA')
        elif hasattr(parameters, 'DIR_RAWDATA_CHANNEL_' + str(channel_id)):
            self._main_directory = getattr(parameters, 'DIR_RAWDATA_CHANNEL_' + str(channel_id))

        if hasattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL_' + str(channel_id)):
            self.angle0_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL_' + str(channel_id))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL' + str(channel_id)):
            self.angle0_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKZERO_CHANNEL' + str(channel_id))
        elif channel_id == 1 and hasattr(parameters, 'DIR_LEFTCAM_STACKZERO'):
            self.angle0_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKZERO')

        if hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL_' + str(channel_id)):
            self.angle1_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL_' + str(channel_id))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL' + str(channel_id)):
            self.angle1_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKZERO_CHANNEL' + str(channel_id))
        elif channel_id == 1 and hasattr(parameters, 'DIR_RIGHTCAM_STACKZERO'):
            self.angle1_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKZERO')

        if hasattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL_' + str(channel_id)):
            self.angle2_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL_' + str(channel_id))
        elif hasattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL' + str(channel_id)):
            self.angle2_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKONE_CHANNEL' + str(channel_id))
        elif channel_id == 1 and hasattr(parameters, 'DIR_LEFTCAM_STACKONE'):
            self.angle2_sub_directory = getattr(parameters, 'DIR_LEFTCAM_STACKONE')

        if hasattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL_' + str(channel_id)):
            self.angle3_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL_' + str(channel_id))
        elif hasattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL' + str(channel_id)):
            self.angle3_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKONE_CHANNEL' + str(channel_id))
        elif channel_id == 1 and hasattr(parameters, 'DIR_RIGHTCAM_STACKONE'):
            self.angle3_sub_directory = getattr(parameters, 'DIR_RIGHTCAM_STACKONE')

        if hasattr(parameters, 'fusion_weighting'):
            self.fusion_weighting = getattr(parameters, 'fusion_weighting')
        elif hasattr(parameters, 'fusion_weighting_channel_' + str(channel_id)):
            self.fusion_weighting = getattr(parameters, 'fusion_weighting_channel_' + str(channel_id))

        return

    ############################################################
    #
    # getters
    #
    ############################################################

    def get_angle_path(self, angle_id):
        if angle_id is 0:
            return os.path.join(self._parent_directory, self._main_directory, self.angle0_sub_directory)
        elif angle_id is 1:
            return os.path.join(self._parent_directory, self._main_directory, self.angle1_sub_directory)
        elif angle_id is 2:
            return os.path.join(self._parent_directory, self._main_directory, self.angle2_sub_directory)
        elif angle_id is 3:
            return os.path.join(self._parent_directory, self._main_directory, self.angle3_sub_directory)
        return None

    def get_image_name(self, angle_id, time_value):
        if angle_id is 0:
            return self._angle0_file_prefix + self.timepoint_to_str(time_value)
        elif angle_id is 1:
            return self._angle1_file_prefix + self.timepoint_to_str(time_value)
        elif angle_id is 2:
            return self._angle2_file_prefix + self.timepoint_to_str(time_value)
        elif angle_id is 3:
            return self._angle3_file_prefix + self.timepoint_to_str(time_value)
        return None

    ############################################################
    #
    # setters
    #
    ############################################################

    def set_parent_directory(self, parentdir):
        if not os.path.isdir(parentdir):
            print("RawdataChannel.set_parent_directory: '" + str(parentdir) + "' is not a valid directory.")
            print("\t Exiting.")
            sys.exit(1)
        self._parent_directory = parentdir
        return

    def set_time_digits_for_acquisition(self, time_digits=3):
        self._time_digits_for_acquisition = time_digits
        return

    ############################################################
    #
    # misc
    #
    ############################################################

    def identical_sub_directories(self, c):
        if self.angle0_sub_directory != c.angle0_sub_directory \
                and self.angle1_sub_directory != c.angle1_sub_directory \
                and self.angle2_sub_directory != c.angle2_sub_directory \
                and self.angle3_sub_directory != c.angle3_sub_directory:
            return False
        return True

    def is_empty(self):
        if self.angle0_sub_directory is None or self.angle1_sub_directory is None or self.angle2_sub_directory is None \
                or self.angle3_sub_directory is None:
            return True
        return False

    def sub_directories_are_different(self):
        if self.angle0_sub_directory != self.angle1_sub_directory \
                and self.angle0_sub_directory != self.angle2_sub_directory \
                and self.angle0_sub_directory != self.angle3_sub_directory \
                and self.angle1_sub_directory != self.angle2_sub_directory \
                and self.angle1_sub_directory != self.angle3_sub_directory \
                and self.angle2_sub_directory != self.angle3_sub_directory:
            return True
        return False

    def timepoint_to_str(self, i):
        return '{:0{width}d}'.format(i, width=self._time_digits_for_acquisition)

#
#
#


class RawdataSubdirectory(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):
        #
        # Channels
        #
        self._n_max_channels = 3
        self.channel = list()
        for i in range(self._n_max_channels):
            c = RawdataChannel(i)
            self.channel.append(c)
        self._time_digits_for_acquisition = 3
        self.set_time_digits_for_acquisition(self._time_digits_for_acquisition)
        return

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        for i in range(self._n_max_channels):
            self.channel[i].print_parameters('channel #' + str(i))
        print('  - _time_digits_for_acquisition = ' + str(self._time_digits_for_acquisition))
        return

    def write_parameters_in_file(self, logfile):
        for i in range(self._n_max_channels):
            self.channel[i].write_parameters_in_file(logfile, 'channel #' + str(i))
        logfile.write('  - _time_digits_for_acquisition = ' + str(self._time_digits_for_acquisition) + '\n')
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameters):
        for i in range(self._n_max_channels):
            self.channel[i].update_from_parameters(i+1, parameters)
        return

    ############################################################
    #
    # setters
    #
    ############################################################

    def get_number_channels(self):
        for i in range(self._n_max_channels):
            if self.channel[i].is_empty():
                return i
        return None

    def get_tmp_directory(self, i, channel_id=0):
        """
        :return:
        """
        if self.channel[channel_id].tmp_directory is None:
            return None
        return _get_directory(self.channel[channel_id].tmp_directory, i)

    def get_time_digits_for_acquisition(self):
        return self._time_digits_for_acquisition

    ############################################################
    #
    # setters
    #
    ############################################################

    def set_parent_directory(self, parentdir):
        if not os.path.isdir(parentdir):
            print("RawdataSubdirectory.set_parent_directory: '" + str(parentdir) + "' is not a valid directory.")
            print("\t Exiting.")
            sys.exit(1)
        for i in range(self._n_max_channels):
            self.channel[i].set_parent_directory(parentdir)
        return

    def set_time_digits_for_acquisition(self, time_digits=3):
        self._time_digits_for_acquisition = time_digits
        for i in range(self._n_max_channels):
            self.channel[i].set_time_digits_for_acquisition(time_digits)
        return


##################################################
#
# sub-directories management
#
##################################################

class GenericSubdirectory(object):
    """
    Class for defining sub-directories of the form
    <main_directory>/<sub_directory_prefix><sub_directory_suffix>
    e.g. <FUSE>/<FUSE_><RELEASE>
    """

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):

        #
        # how to build sub-directory names
        # _parent_directory is /path/to/experiment
        # ./<_main_directory>/<_sub_directory_prefix><_sub_directory_suffix>
        #
        self._parent_directory = None
        self._main_directory = None
        self._sub_directory_prefix = None
        self._sub_directory_suffix = None

        #
        # sub-directories
        # 1. main result directory(ies)
        # 2. log directory (subdirectory of #1)
        # 3. temporary directory (subdirectory of #1)
        # 4. reconstruction directory (subdirectory of #1)
        #    used for segmentation, when image have been pre-processed before segmentation
        #
        self._directory = None
        self._sub_directory = None
        self._log_directory = None
        self._tmp_directory = None
        self._rec_directory = None

        #
        # file related args
        # <_file_prefix><_file_suffix>_t<time_point>
        #
        self._file_prefix = None
        self._file_suffix = None
        self._time_prefix = '_t'
        self._time_digits_for_filename = 3

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        self._set_directory()
        self._set_log_directory()
        self._set_rec_directory()
        print("    - _parent_directory = " + str(self._parent_directory))
        print("    - _main_directory = " + str(self._main_directory))
        print("    - _sub_directory_prefix = " + str(self._sub_directory_prefix))
        print("    - _sub_directory_suffix = " + str(self._sub_directory_suffix))
        #
        print("    - _directory = ")
        _print_directory_list(self._directory)
        print("    - _log_directory = ")
        _print_directory_list(self._log_directory)
        print("    - _tmp_directory = ")
        _print_directory_list(self._tmp_directory)
        print("    - _rec_directory = ")
        _print_directory_list(self._rec_directory)
        #
        print("    - _file_prefix = " + str(self._file_prefix))
        print("    - _file_suffix = " + str(self._file_suffix))
        print("    - _time_digits_for_filename = " + str(self._time_digits_for_filename))
        return

    def _write_directory_list_in_file(self, logfile):
        self._set_directory()
        if self._directory is None:
            logfile.write("    None\n")
        else:
            if len(self._directory) <= 0:
                logfile.write("    Empty list\n")
            else:
                for i in range(len(self._directory)):
                    logfile.write("      #" + str(i) + ": '" + str(self.get_directory(i)) + "'\n")

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        self._set_log_directory()
        self._set_rec_directory()
        #
        logfile.write("    - _parent_directory = " + str(self._parent_directory) + "\n")
        logfile.write("    - _main_directory = " + str(self._main_directory) + "\n")
        logfile.write("    - _sub_directory_prefix = " + str(self._sub_directory_prefix) + "\n")
        logfile.write("    - _sub_directory_suffix = " + str(self._sub_directory_suffix) + "\n")
        #
        logfile.write("    - _directory = " + "\n")
        _write_directory_list_in_file(self._directory, logfile)
        logfile.write("    - _log_directory = " + "\n")
        _write_directory_list_in_file(self._log_directory, logfile)
        logfile.write("    - _tmp_directory = " + "\n")
        _write_directory_list_in_file(self._tmp_directory, logfile)
        logfile.write("    - _rec_directory = " + "\n")
        _write_directory_list_in_file(self._rec_directory, logfile)
        #
        logfile.write("    - _file_prefix = " + str(self._file_prefix) + "\n")
        logfile.write("    - _file_suffix = " + str(self._file_suffix) + "\n")
        logfile.write("    - _time_digits_for_filename = " + str(self._time_digits_for_filename) + "\n")
        return

    ############################################################
    #
    # update
    #
    ############################################################

    ############################################################
    #
    # getters
    #
    ############################################################

    def get_directory(self, i=0):
        """
        return the ith directory
        :param i:
        :return:
        """
        if self._directory is None:
            self._set_directory()
        return _get_directory(self._directory, i)

    def get_directory_suffix(self):
        return self._sub_directory_suffix

    def get_file_prefix(self):
        return self._file_prefix

    def get_file_suffix(self):
        return self._file_suffix

    def get_image_format(self):
        return self._file_prefix + self._file_suffix + self._time_prefix + self._get_time_format()

    def get_image_name(self, time_value):
        if time_value is None:
            return None
        return self._file_prefix + self._file_suffix + self._time_prefix + self.timepoint_to_str(time_value)

    def get_log_directory(self, i=0):
        """
        :return:
        """
        if self._log_directory is None:
            self._set_log_directory()
        return _get_directory(self._log_directory, i)

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

    def get_number_sub_directories(self):
        #
        # self._directory is None
        # build directory list
        #
        if self._sub_directory is None:
            self._set_sub_directory()

        if len(self._sub_directory) <= 0:
            return 0

        if type(self._sub_directory) == str:
            return 1

        if type(self._sub_directory) == list or type(self._sub_directory) == tuple:
            return len(self._sub_directory)

        return 0

    def get_sub_directory(self, i=0):
        """
        return the ith directory
        :param i:
        :return:
        """
        if self._sub_directory is None:
            self._set_directory()
        return _get_directory(self._sub_directory, i)

    def _get_time_format(self):
        form = "%0" + str(self._time_digits_for_filename) + "d"
        return form

    def get_time_prefix(self):
        return self._time_prefix

    def get_tmp_directory(self, i=0):
        """
        :return:
        """
        #
        # temporary directories are built with a timepoint variable
        #
        if self._tmp_directory is None:
            return None
        return _get_directory(self._tmp_directory, i)

    def get_rec_directory(self, i=0):
        """
        return the 'reconstruction' sub-directory.
        When segmenting an image, a reconstructed image can be used as support for segmentation.
        :return:
        """
        if self._rec_directory is None:
            self._set_rec_directory()
        return _get_directory(self._rec_directory, i)

    ############################################################
    #
    # setters
    #
    ############################################################

    def _set_directory(self, force=False):
        """
        build the list of directories from the directory suffixes
        :return:
        """
        if self._directory is not None and force is False:
            return

        self._set_sub_directory(force)
        #
        # empty list
        #
        self._directory = []
        for i in range(self.get_number_sub_directories()):
            self._directory.append(os.path.join(str(self._parent_directory), self.get_sub_directory(i)))
        return

    def set_directory_suffix(self, directory_suffix):
        self._sub_directory_suffix = directory_suffix
        self._set_directory(force=True)

    def set_file_prefix(self, file_prefix):
        self._file_prefix = file_prefix
        return

    def _set_log_directory(self):
        """
        :return:
        """
        if self._log_directory is not None:
            return
        #
        # empty list
        #
        self._set_directory()
        self._log_directory = []
        for i in range(len(self._directory)):
            self._log_directory.append(os.path.join(self.get_directory(i), "LOGS"))
        return

    def set_parent_directory(self, parentdir):
        if not os.path.isdir(parentdir):
            print("GenericSubdirectory.set_parent_directory: '" + str(parentdir) + "' is not a valid directory.")
            print("\t Exiting.")
            sys.exit(1)
        self._parent_directory = parentdir
        return

    def _set_rec_directory(self):
        """
        :return:
        """
        self._set_directory()
        self._rec_directory = []
        for i in range(len(self._directory)):
            self._rec_directory.append(os.path.join(self.get_directory(i), "RECONSTRUCTION"))
        return

    def set_rec_directory_to_tmp(self):
        """
        :return:
        """
        self._set_directory()
        self._rec_directory = self._tmp_directory
        return

    def _set_sub_directory(self, force=False):
        """
        build the list of directories from the directory suffixes
        :return:
        """
        proc = "GenericSubdirectory._set_sub_directory"
        if self._sub_directory is not None and force is False:
            return
        #
        # empty list
        #
        self._sub_directory = []
        if self._sub_directory_suffix is None:
            monitoring.to_log_and_console("Warning: " + proc + ", _sub_directory_suffix is None")
        elif type(self._sub_directory_suffix) == str:
            subdir = str(self._sub_directory_prefix) + str(self._sub_directory_suffix)
            self._sub_directory.append(os.path.join(str(self._main_directory), subdir))
        elif type(self._sub_directory_suffix) == list or type(self._sub_directory_suffix) == tuple:
            for s in self._sub_directory_suffix:
                subdir = str(self._sub_directory_prefix) + str(s)
                self._sub_directory.append(os.path.join(str(self._main_directory), subdir))
        else:
            monitoring.to_log_and_console("Warning: " + proc + ", unhandled _sub_directory_suffix type")
        return

    def set_tmp_directory(self, timepoint):
        """
        :return:
        """
        self._set_directory()
        self._tmp_directory = []
        for i in range(len(self._directory)):
            self._tmp_directory.append(os.path.join(self.get_directory(i), "TEMP_" + self.timepoint_to_str(timepoint)))
        return
    
    def set_time_digits_for_filename(self, time_digits=3):
        self._time_digits_for_filename = time_digits
        return

    ############################################################
    #
    # misc
    #
    ############################################################

    def get_history_filename(self, cli_name):
        if cli_name is None:
            local_cli_name = 'unknown'
        else:
            local_cli_name = os.path.basename(cli_name)
            if local_cli_name[-3:] == '.py':
                local_cli_name = local_cli_name[:-3]
        return os.path.join(self.get_log_directory(), local_cli_name + '-history.log')

    #
    # manage directories
    #

    def make_directory(self):
        self._set_directory()
        _make_directory_list(self._directory)

    def make_tmp_directory(self):
        _make_directory_list(self._tmp_directory)

    def make_rec_directory(self):
        self._set_rec_directory()
        _make_directory_list(self._rec_directory)

    def rmtree_tmp_directory(self):
        _rmtree_directory_list(self._tmp_directory)

    #
    # misc
    #

    def timepoint_to_str(self, i):
        if type(i) is int:
            return '{:0{width}d}'.format(i, width=self._time_digits_for_filename)
        elif type(i) is str:
            return '{:0{width}d}'.format(int(i), width=self._time_digits_for_filename)
        else:
            print("GenericSubdirectory.timepoint_to_str: type '" + str(type(i)) + "' not handled yet.")
            return None

#
#
# FUSION sub-directory
#
#


class FuseSubdirectory(GenericSubdirectory):

    def __init__(self):
        GenericSubdirectory.__init__(self)
        self._main_directory = 'FUSE'
        self._sub_directory_prefix = 'FUSE_'
        self._sub_directory_suffix = 'RELEASE'
        self._file_suffix = "_fuse"

        self._xzsection_directory = list()

    def print_parameters(self):
        self._set_directory()
        print("  - subpath/to/fusion is")
        GenericSubdirectory.print_parameters(self)
        return

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        logfile.write("  - subpath/to/fusion is \n")
        GenericSubdirectory.write_parameters_in_file(self, logfile)

    def update_from_parameters(self, parameters):
        if hasattr(parameters, 'EXP_FUSE'):
            if parameters.EXP_FUSE is not None:
                self._sub_directory_suffix = parameters.EXP_FUSE
        return

    def set_xzsection_directory(self, time_value):
        t = "XZSECTION_" + self.timepoint_to_str(time_value)
        for c in range(self.get_number_directories()):
            d = os.path.join(self.get_directory(c), t)
            self._xzsection_directory.append(d)
        return

    def get_xzsection_directory(self, channel_id=0):
        d = _get_directory(self._xzsection_directory, channel_id)
        if not os.path.isdir(d):
            os.makedirs(d)
        return d

#
#
# INTRAREG sub-directory
#
#


class IntraregSubdirectory(GenericSubdirectory):

    def __init__(self):
        GenericSubdirectory.__init__(self)
        self._main_directory = 'INTRAREG'
        self._sub_directory_prefix = 'INTRAREG_'
        self._sub_directory_suffix = 'RELEASE'
        self._file_suffix = "_intrareg"

    def print_parameters(self):
        self._set_directory()
        print("  - subpath/to/intraregistration is")
        GenericSubdirectory.print_parameters(self)
        return

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        logfile.write("  - subpath/to/intraregistration is \n")
        GenericSubdirectory.write_parameters_in_file(self, logfile)

    def update_from_parameters(self, parameters):
        if hasattr(parameters, 'EXP_INTRAREG'):
            if parameters.EXP_INTRAREG is not None:
                self._sub_directory_suffix = parameters.EXP_INTRAREG
        return

#
#
# MARS sub-directory
#
#


class MarsSubdirectory(GenericSubdirectory):

    def __init__(self):
        GenericSubdirectory.__init__(self)
        self._main_directory = 'SEG'
        self._sub_directory_prefix = 'SEG_'
        self._sub_directory_suffix = 'RELEASE'
        self._file_suffix = "_mars"

    def print_parameters(self):
        self._set_directory()
        print("  - subpath/to/mars is")
        GenericSubdirectory.print_parameters(self)
        return

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        logfile.write("  - subpath/to/mars is \n")
        GenericSubdirectory.write_parameters_in_file(self, logfile)

    def update_from_parameters(self, parameters):
        if hasattr(parameters, 'EXP_SEG'):
            if parameters.EXP_SEG is not None:
                self._sub_directory_suffix = parameters.EXP_SEG
        if hasattr(parameters, 'EXP_MARS'):
            if parameters.EXP_MARS is not None:
                self._sub_directory_suffix = parameters.EXP_MARS
        return

#
#
# Astec sub-directory
#
#


class AstecSubdirectory(GenericSubdirectory):

    def __init__(self):
        GenericSubdirectory.__init__(self)
        self._main_directory = 'SEG'
        self._sub_directory_prefix = 'SEG_'
        self._sub_directory_suffix = 'RELEASE'
        self._file_suffix = "_seg"

    def print_parameters(self):
        self._set_directory()
        print("  - subpath/to/seg is")
        GenericSubdirectory.print_parameters(self)
        return

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        logfile.write("  - subpath/to/seg is \n")
        GenericSubdirectory.write_parameters_in_file(self, logfile)

    def update_from_parameters(self, parameters):
        if hasattr(parameters, 'EXP_SEG'):
            if parameters.EXP_SEG is not None:
                self._sub_directory_suffix = parameters.EXP_SEG
        return

#
#
# POST sub-directory
#
#


class PostSubdirectory(GenericSubdirectory):

    def __init__(self):
        GenericSubdirectory.__init__(self)
        self._main_directory = 'POST'
        self._sub_directory_prefix = 'POST_'
        self._sub_directory_suffix = 'RELEASE'
        self._file_suffix = "_post"

    def print_parameters(self):
        self._set_directory()
        print("  - subpath/to/postcorrection is")
        GenericSubdirectory.print_parameters(self)
        return

    def write_parameters_in_file(self, logfile):
        self._set_directory()
        logfile.write("  - subpath/to/postcorrection is \n")
        GenericSubdirectory.write_parameters_in_file(self, logfile)

    def update_from_parameters(self, parameters):
        if hasattr(parameters, 'EXP_POST'):
            if parameters.EXP_POST is not None:
                self._sub_directory_suffix = parameters.EXP_POST
        return


##################################################
#
# full experiment management
#
##################################################


class Experiment(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):

        self._embryo_path = None
        self._embryo_name = None

        self.first_time_point = -1
        self.last_time_point = -1
        self.delta_time_point = 1
        self.delay_time_point = 0

        self._time_digits_for_filename = 3
        self._time_digits_for_cell_id = 4

        #
        # sub-directories
        # 1. specialized ones
        #    for automated initialisation
        # 2. generic one
        #    copy from the targeted specialized one
        #
        self.rawdata_dir = RawdataSubdirectory()
        self.fusion_dir = FuseSubdirectory()
        self.mars_dir = MarsSubdirectory()
        self.astec_dir = AstecSubdirectory()
        self.post_dir = PostSubdirectory()
        self.intrareg_dir = IntraregSubdirectory()

        self.working_dir = GenericSubdirectory()

        #
        #
        #
        self.set_time_digits_for_filename(self._time_digits_for_filename)

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'mha'
        self.default_image_suffix = 'mha'

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('Experiment parameters')

        print('- _embryo_path is ' + str(self._embryo_path))
        print('- _embryo_name is ' + str(self._embryo_name))

        print('- first_time_point is ' + str(self.first_time_point))
        print('- last_time_point is ' + str(self.last_time_point))
        print('- delta_time_point is ' + str(self.delta_time_point))
        print('- delay_time_point is ' + str(self.delay_time_point))

        print('- _time_digits_for_filename is ' + str(self._time_digits_for_filename))
        print('- _time_digits_for_cell_id is ' + str(self._time_digits_for_cell_id))

        print('- raw data directory is')
        self.rawdata_dir.print_parameters()
        print('- fusion directory is')
        self.fusion_dir.print_parameters()
        print('- mars directory is')
        self.mars_dir.print_parameters()
        print('- segmentation directory is')
        self.astec_dir.print_parameters()
        print('- post-correction directory is')
        self.post_dir.print_parameters()
        print('- intra-registration directory is')
        self.intrareg_dir.print_parameters()
        print('- working directory is')
        self.working_dir.print_parameters()

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))
        print("")
        return

    def write_parameters(self, log_filename=None):
        if log_filename is not None:
            local_log_filename = log_filename
        else:
            local_log_filename = monitoring.log_filename
        if local_log_filename is not None:
            with open(local_log_filename, 'a') as logfile:
                logfile.write("\n")
                logfile.write('Experiment parameters\n')

                logfile.write('- _embryo_path is ' + str(self.get_embryo_path())+'\n')
                logfile.write('- _embryo_name is ' + str(self.get_embryo_name())+'\n')

                logfile.write('- first_time_point is ' + str(self.first_time_point)+'\n')
                logfile.write('- last_time_point is ' + str(self.last_time_point)+'\n')
                logfile.write('- delta_time_point is ' + str(self.delta_time_point)+'\n')
                logfile.write('- delay_time_point is ' + str(self.delay_time_point)+'\n')

                logfile.write('- _time_digits_for_filename is ' + str(self._time_digits_for_filename) + '\n')
                logfile.write('- _time_digits_for_cell_id is ' + str(self._time_digits_for_cell_id) + '\n')

                logfile.write('- raw data directory is \n')
                self.rawdata_dir.write_parameters_in_file(logfile)
                logfile.write('- fusion directory is \n')
                self.fusion_dir.write_parameters_in_file(logfile)
                logfile.write('- mars directory is \n')
                self.mars_dir.write_parameters_in_file(logfile)
                logfile.write('- segmentation directory is \n')
                self.astec_dir.write_parameters_in_file(logfile)
                logfile.write('- post-correction directory is \n')
                self.post_dir.write_parameters_in_file(logfile)
                logfile.write('- intra-registration directory is \n')
                self.intrareg_dir.write_parameters_in_file(logfile)
                logfile.write('- working directory is \n')
                self.working_dir.write_parameters_in_file(logfile)

                logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
                logfile.write('- default_image_suffix = ' + str(self.default_image_suffix) + '\n')

                logfile.write("\n")
        return

    def update_history_at_start(self, cli_name=None, start_time=None, parameter_file=None, path_to_vt=None):
        history_filename = self.working_dir.get_history_filename(cli_name)
        with open(history_filename, 'a') as logfile:
            logfile.write("\n")
            if start_time is not None:
                logfile.write("# " + time.strftime("%a, %d %b %Y %H:%M:%S", start_time) + "\n")
            logfile.write("# Embryo path: '" + str(self.get_embryo_path()) + "'\n")
            logfile.write("# Embryo name: '" + str(self.get_embryo_name()) + "'\n")
            if parameter_file is not None:
                logfile.write("# Parameter file: '" + str(parameter_file) + "'\n")
            logfile.write("# Command line: '" + " ".join(sys.argv) + "'\n")
            logfile.write("# Working directory: '" + str(os.getcwd()) + "'\n")
            logfile.write("# User: '" + str(getpass.getuser()) + "'\n")
            logfile.write("# Python executable: " + sys.executable + "\n")
            if cli_name is not None:
                _write_git_information(os.path.dirname(cli_name), logfile, "# ASTEC")
            if path_to_vt is not None:
                _write_git_information(path_to_vt, logfile, "# VT")
            logfile.write("# \n")
        return

    def update_history_execution_time(self, cli_name=None, start_time=None, end_time=None):
        history_filename = self.working_dir.get_history_filename(cli_name)
        with open(history_filename, 'a') as logfile:
            logfile.write('# Total execution time = ' + str(time.mktime(end_time) - time.mktime(start_time)) + ' sec\n')
            logfile.write("\n\n")
        return

    def copy_stamped_file(self, timestamp=None, thefile=None, directory=None):
        if thefile is None:
            return
        if directory is not None:
            local_directory = directory
        else:
            local_directory = self.working_dir.get_log_directory()

        d = _timestamp_to_str(timestamp)
        if len(thefile.split('.')) > 1:
            ext = thefile.split('.')[-1]
            filename = re.sub(r'(\.*).' + ext, r'\1', thefile.split(os.path.sep)[-1]) + '-' + d + '.' + ext
        else:
            filename = thefile.split(os.path.sep)[-1] + '-' + d

        resfile = os.path.join(local_directory, filename)
        shutil.copy2(thefile, resfile)
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_args(self, args):
        """

        :param args:
        :return:
        """
        if args.embryo_path is None:
            return
        self._set_embryo_path(args.embryo_path)
        return

    #
    #
    #

    def _update_embryo_path_from_parameters(self, parameters):
        proc = '_update_embryo_path_from_parameters'
        if hasattr(parameters, 'PATH_EMBRYO'):
            if parameters.PATH_EMBRYO is not None:
                if not os.path.isdir(parameters.PATH_EMBRYO):
                    print(proc + ": '" + parameters.PATH_EMBRYO + "' is not a valid directory. Exiting.")
                    sys.exit(1)
                self._set_embryo_path(parameters.PATH_EMBRYO)
            else:
                self._set_embryo_path(os.getcwd())
        else:
            self._set_embryo_path(os.getcwd())
        return

    def _embryo_name_from_embryo_path(self):
        sp = self._embryo_path.split(os.path.sep)
        if len(sp) is 0:
            return
        if sp[-1] is not '':
            self._set_embryo_name(sp[-1])
            return
        if len(sp) >= 2 and sp[-2] is not '':
            self._set_embryo_name(sp[-2])
            return
        return

    def _update_embryo_name_from_parameters(self, parameters):
        self._update_embryo_path_from_parameters(parameters)
        if hasattr(parameters, 'EN'):
            if parameters.EN is not None:
                self._set_embryo_name(parameters.EN)
            else:
                self._embryo_name_from_embryo_path()
        else:
            self._embryo_name_from_embryo_path()
        return

    def update_from_parameters(self, parameter_file):
        """

        :param parameter_file:
        :return:
        """
        proc = "Experiment.update_from_parameters"

        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print(proc + ": '" + parameter_file + "' is not a valid file.")
            print("\t Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self._update_embryo_path_from_parameters(parameters)
        self._update_embryo_name_from_parameters(parameters)

        if hasattr(parameters, 'begin'):
            if parameters.begin is not None:
                self.first_time_point = parameters.begin
            else:
                print(proc + ": it is mandatory to specify the first time point")
                print("\t Exiting.")
                sys.exit(1)

        if hasattr(parameters, 'end'):
            if parameters.end is not None:
                self.last_time_point = parameters.end
            else:
                print(proc + ": it is mandatory to specify the last time point")
                print("\t Exiting.")
                sys.exit(1)

        if hasattr(parameters, 'delta'):
            if parameters.delta is not None:
                self.delta_time_point = parameters.delta

        if hasattr(parameters, 'raw_delay'):
            if parameters.raw_delay is not None:
                self.delay_time_point = parameters.raw_delay

        if hasattr(parameters, 'time_digits_for_filename'):
            if parameters.time_digits_for_filename is not None:
                self.set_time_digits_for_filename(parameters.time_digits_for_filename)
                
        if hasattr(parameters, 'time_digits_for_cell_id'):
            if parameters.time_digits_for_cell_id is not None:
                self._time_digits_for_cell_id = parameters.time_digits_for_cell_id

        self.rawdata_dir.update_from_parameters(parameters)
        self.fusion_dir.update_from_parameters(parameters)
        self.mars_dir.update_from_parameters(parameters)
        self.astec_dir.update_from_parameters(parameters)
        self.post_dir.update_from_parameters(parameters)
        self.intrareg_dir.update_from_parameters(parameters)

        self._fix_fusion_directories()

        #
        # set result_image_suffix
        #
        if hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.RESULT_IMAGE_SUFFIX_FUSE
        if hasattr(parameters, 'result_image_suffix'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.result_image_suffix

        #
        # set result_image_suffix to default_image_suffix
        # if default_image_suffix is given and not result_image_suffix
        #
        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix
                if not hasattr(parameters, 'result_image_suffix') \
                        and not hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
                    self.result_image_suffix = parameters.default_image_suffix
        return

    ############################################################
    #
    # getters
    #
    ############################################################

    def get_embryo_name(self):
        if self._embryo_name is None:
            self._set_embryo_name(self._embryo_path.split(os.path.sep)[-1])
        return self._embryo_name

    def get_embryo_path(self):
        return self._embryo_path

    def get_log_dirname(self):
        return os.path.join(self.working_dir.get_log_directory())

    def get_segmentation_image(self, time_value):
        proc = 'Experiment.get_segmentation_image'
        #
        # try to find segmentation image
        #
        if time_value is None:
            return None
        seg_name = self.astec_dir.get_image_name(time_value)
        seg_image = find_file(self.astec_dir.get_directory(), seg_name, callfrom=proc, local_monitoring=None,
                              verbose=False)
        if seg_image is not None:
            return os.path.join(self.astec_dir.get_directory(), seg_image)
        #
        # try to find mars image
        #
        mars_name = self.mars_dir.get_image_name(time_value)
        mars_image = find_file(self.mars_dir.get_directory(), mars_name, callfrom=proc, local_monitoring=None,
                               verbose=False)
        if mars_image is not None:
            #
            # copy mars image with the same suffix
            #
            seg_name += "." + mars_image[len(mars_name)+1:]
            monitoring.to_log_and_console("    .. " + proc + ": copy '" + str(mars_image) + "' into '"
                                          + str(seg_name) + "'", 2)
            mars_image = os.path.join(self.mars_dir.get_directory(), mars_image)
            seg_image = os.path.join(self.astec_dir.get_directory(), seg_name)
            shutil.copy2(mars_image, seg_image)
            return seg_image
        #
        #
        #
        monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                      + str(time_value), 2)
        monitoring.to_log_and_console("       \t" + "1. was looking for file '" + str(seg_name) + "'", 2)
        monitoring.to_log_and_console("       \t" + "   in directory '" + str(self.astec_dir.get_directory()) + "'", 2)
        monitoring.to_log_and_console("       \t" + "2. was looking for file '" + str(mars_name) + "'", 2)
        monitoring.to_log_and_console("       \t" + "   in directory '" + str(self.mars_dir.get_directory()) + "'", 2)

        return None

    def get_time_digits_for_filename(self):
        return self._time_digits_for_filename

    def get_time_digits_for_cell_id(self):
        return self._time_digits_for_cell_id

    def get_time_format(self):
        form = "%0" + str(self._time_digits_for_filename) + "d"
        return form

    def get_time_index(self, index):
        ind = '{:0{width}d}'.format(index, width=self._time_digits_for_filename)
        return ind

    ############################################################
    #
    # setters
    #
    ############################################################

    def _set_embryo_name(self, embryo_name):
        self._embryo_name = embryo_name
        #
        # set file prefix
        #
        self.fusion_dir.set_file_prefix(embryo_name)
        self.mars_dir.set_file_prefix(embryo_name)
        self.astec_dir.set_file_prefix(embryo_name)
        self.post_dir.set_file_prefix(embryo_name)
        self.intrareg_dir.set_file_prefix(embryo_name)
        return

    def _set_embryo_path(self, embryo_path):
        if not os.path.isdir(embryo_path):
            print("Experiment._set_embryo_path: '" + str(embryo_path) + "' is not a valid directory.")
            print("\t Exiting.")
            sys.exit(1)
        self._embryo_path = embryo_path
        #
        # set parent directories
        #
        self.rawdata_dir.set_parent_directory(embryo_path)
        self.fusion_dir.set_parent_directory(embryo_path)
        self.mars_dir.set_parent_directory(embryo_path)
        self.astec_dir.set_parent_directory(embryo_path)
        self.post_dir.set_parent_directory(embryo_path)
        self.intrareg_dir.set_parent_directory(embryo_path)
        return

    def set_fusion_tmp_directory(self, time_value):
        if self.rawdata_dir.get_number_channels() is not self.fusion_dir.get_number_directories():
            print("Experiment.set_fusion_tmp_directory: number of channels is different from fusion directories.")
            print("\t Exiting.")
            sys.exit(1)
        t = "TEMP_" + self.fusion_dir.timepoint_to_str(time_value)
        for c in range(self.fusion_dir.get_number_directories()):
            self.rawdata_dir.channel[c].tmp_directory = []
            d = os.path.join(self.fusion_dir.get_directory(c), t, "ANGLE_0")
            self.rawdata_dir.channel[c].tmp_directory.append(d)
            d = os.path.join(self.fusion_dir.get_directory(c), t, "ANGLE_1")
            self.rawdata_dir.channel[c].tmp_directory.append(d)
            d = os.path.join(self.fusion_dir.get_directory(c), t, "ANGLE_2")
            self.rawdata_dir.channel[c].tmp_directory.append(d)
            d = os.path.join(self.fusion_dir.get_directory(c), t, "ANGLE_3")
            self.rawdata_dir.channel[c].tmp_directory.append(d)
            d = os.path.join(self.fusion_dir.get_directory(c), t)
            self.rawdata_dir.channel[c].tmp_directory.append(d)
            for d in self.rawdata_dir.channel[c].tmp_directory:
                if not os.path.isdir(d):
                    os.makedirs(d)

    def set_time_digits_for_filename(self, time_digits=3):
        self._time_digits_for_filename = time_digits
        self.fusion_dir.set_time_digits_for_filename(time_digits)
        self.mars_dir.set_time_digits_for_filename(time_digits)
        self.astec_dir.set_time_digits_for_filename(time_digits)
        self.post_dir.set_time_digits_for_filename(time_digits)
        self.intrareg_dir.set_time_digits_for_filename(time_digits)

    ############################################################
    #
    # misc
    #
    ############################################################

    def _fix_fusion_directories(self):
        n_channels = self.rawdata_dir.get_number_channels()
        if n_channels is 1:
            return
        n_dirs = self.fusion_dir.get_number_directories()
        if n_dirs >= n_channels:
            return

        dir_suffix = self.fusion_dir.get_directory_suffix()
        new_dir_suffix = []
        if type(dir_suffix) is str:
            for c in range(n_channels):
                if c is 0:
                    new_dir_suffix.append(dir_suffix)
                else:
                    new_dir_suffix.append(dir_suffix + "_CHANNEL_" + str(c+1))
        elif type(dir_suffix) is list:
            for c in range(n_channels):
                if c < n_dirs:
                    new_dir_suffix.append(dir_suffix[c])
                else:
                    new_dir_suffix.append(dir_suffix[0] + "_CHANNEL_" + str(c+1))
        else:
            print("Experiment._fix_fusion_directories: type '" + str(type(dir_suffix)) + "' not handled yet.")
            print("\t Exiting.")
            sys.exit(1)
        self.fusion_dir.set_directory_suffix(new_dir_suffix)

    def remove_fusion_tmp_directory(self):
        self.rawdata_dir.tmp_directory = []
        for c in range(self.fusion_dir.get_number_directories()):
            shutil.rmtree(self.rawdata_dir.channel[c].tmp_directory[4])


########################################################################################
#
# Registration
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

    ############################################################
    #
    # initialisation
    #
    ############################################################

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
        self.pyramid_highest_level = 6
        self.pyramid_lowest_level = 3
        self.gaussian_pyramid = False
        self.transformation_type = 'affine'
        self.elastic_sigma = 4.0
        self.transformation_estimation_type = 'wlts'
        self.lts_fraction = 0.55
        self.fluid_sigma = 4.0

        self.normalization = True

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('RegistrationParameters')

        print(_fulldesc(None, 'prefix') + str(self.prefix))
        print(_fulldesc(self.prefix, 'compute_registration') + str(self.compute_registration))

        print(_fulldesc(self.prefix, 'pyramid_highest_level') + str(self.pyramid_highest_level))
        print(_fulldesc(self.prefix, 'pyramid_lowest_level') + str(self.pyramid_lowest_level))
        print(_fulldesc(self.prefix, 'gaussian_pyramid') + str(self.gaussian_pyramid))

        print(_fulldesc(self.prefix, 'transformation_type') + str(self.transformation_type))

        print(_fulldesc(self.prefix, 'elastic_sigma') + str(self.elastic_sigma))

        print(_fulldesc(self.prefix, 'transformation_estimation_type') + str(self.transformation_estimation_type))
        print(_fulldesc(self.prefix, 'lts_fraction') + str(self.lts_fraction))
        print(_fulldesc(self.prefix, 'fluid_sigma') + str(self.fluid_sigma))

        print(_fulldesc(self.prefix, 'normalization') + str(self.normalization))

        print("")

    def write_parameters(self, log_filename=None):
        if log_filename is not None:
            local_log_filename = log_filename
        else:
            local_log_filename = monitoring.log_filename
        if local_log_filename is not None:
            with open(local_log_filename, 'a') as logfile:
                logfile.write("\n")
                logfile.write('RegistrationParameters\n')
                logfile.write(_fulldesc(None, 'prefix')+str(self.prefix)+'\n')
                logfile.write(_fulldesc(self.prefix, 'compute_registration') + str(self.compute_registration) + '\n')

                logfile.write(_fulldesc(self.prefix, 'pyramid_highest_level') + str(self.pyramid_highest_level) + '\n')
                logfile.write(_fulldesc(self.prefix, 'pyramid_lowest_level') + str(self.pyramid_lowest_level) + '\n')
                logfile.write(_fulldesc(self.prefix, 'gaussian_pyramid') + str(self.gaussian_pyramid) + '\n')

                logfile.write(_fulldesc(self.prefix, 'transformation_type') + str(self.transformation_type) + '\n')

                logfile.write(_fulldesc(self.prefix, 'elastic_sigma') + str(self.elastic_sigma) + '\n')

                logfile.write(_fulldesc(self.prefix, 'transformation_estimation_type')
                              + str(self.transformation_estimation_type)+'\n')
                logfile.write(_fulldesc(self.prefix, 'lts_fraction')+str(self.lts_fraction)+'\n')
                logfile.write(_fulldesc(self.prefix, 'fluid_sigma') + str(self.fluid_sigma) + '\n')

                logfile.write(_fulldesc(self.prefix, 'normalization')+str(self.normalization)+'\n')

                logfile.write("\n")
        return

    ############################################################
    #
    # update
    #
    ############################################################

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


def blockmatching(image_ref, image_flo, image_output, trsf_output, trsf_init, parameters, other_options=None):
    """

    :param image_ref:
    :param image_flo:
    :param image_output:
    :param trsf_output:
    :param trsf_init:
    :param parameters:
    :param other_options:
    :return:
    """

    proc = "common.blockmatching"

    #
    # parameter type checking
    #

    if not isinstance(parameters, RegistrationParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    cpp_wrapping.blockmatching(image_ref, image_flo, image_output, trsf_output, trsf_init,
                               py_hl=parameters.pyramid_highest_level, py_ll=parameters.pyramid_lowest_level,
                               gaussian_pyramid=parameters.gaussian_pyramid,
                               transformation_type=parameters.transformation_type,
                               elastic_sigma=parameters.elastic_sigma,
                               transformation_estimator=parameters.transformation_estimation_type,
                               lts_fraction=parameters.lts_fraction, fluid_sigma=parameters.fluid_sigma,
                               normalization=parameters.normalization, other_options=other_options,
                               monitoring=monitoring)
    return


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
    if new_parameter_file is None or new_parameter_file is '':
        print("getParameterFile: no parameter file. Exiting.")
        sys.exit(1)
    if os.path.isfile(new_parameter_file) is not False:
        print("getParameterFile: '"+new_parameter_file+"' is not a valid file. Exiting.")
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

def find_file(data_path, file_prefix, callfrom=None, local_monitoring=None, verbose=True):
    """
    find a file in a directory with a given prefix. The suffix is unknown

    :param data_path:
    :param file_prefix:
    :param local_monitoring:
    :param verbose:
    :param callfrom:
    :return:
    """
    proc = "find_file"

    if not os.path.isdir(data_path):
        if local_monitoring is not None:
            local_monitoring.to_log_and_console("Error:")
            local_monitoring.to_log_and_console(proc + ": '" + str(data_path) + "' is not a valid directory ?!")
            if callfrom is not None:
                local_monitoring.to_log_and_console("\t call from '" + str(callfrom) + "'")
            local_monitoring.to_log_and_console("\t Exiting.")
        else:
            print(proc + ": '" + str(data_path) + "' is not a valid directory ?!")
            if callfrom is not None:
                print("\t call from '" + str(callfrom) + "'")
            print("\t Exiting.")
        sys.exit(1)

    if file_prefix is None:
        print(proc + ": file prefix was 'None'?!")
        if callfrom is not None:
            print("\t call from '" + str(callfrom) + "'")
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
    filenames = []
    for f in os.listdir(data_path):
        if len(f) <= length_file_prefix:
            pass
        if f[0:length_file_prefix] == file_prefix[0:length_file_prefix] and f[length_file_prefix] == '.':
            filenames.append(f)

    if len(filenames) == 0:
        if local_monitoring is not None:
            local_monitoring.to_log_and_console(proc + ": no image with name '" + str(file_prefix)
                                                + "' was found in '" + str(data_path) + "'", 4)
        elif verbose is True:
            print(proc + ": no image with name '" + str(file_prefix) + "' was found in '" + str(data_path) + "'")
        return None

    if len(filenames) > 1:
        if local_monitoring is not None:
            local_monitoring.to_log_and_console("\t " + proc + ": warning")
            local_monitoring.to_log_and_console("\t several images with name '" + str(file_prefix) + "' were found in")
            local_monitoring.to_log_and_console("\t    '" + str(data_path) + "'")
            local_monitoring.to_log_and_console("\t    -> "+str(filenames))
            local_monitoring.to_log_and_console("\t returned file is '" + str(filenames[0]) + "'")
        elif verbose is True:
            print(proc + ": several images with name '"
                  + str(file_prefix) + "' were found in '" + str(data_path) + "'")
            print("\t "+str(filenames))
            print("\t returned file is '" + str(filenames[0]) + "'")
        # return None

    return filenames[0]


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

        time_point = experiment.get_time_index(current_time)
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
        monitoring.to_log_and_console("\t found "+str(nfiles)+" images instead of "+str(nimages))
        monitoring.to_log_and_console("\t Exiting.", 0)
        exit(1)

    monitoring.to_log_and_console(proc + ": no common suffix for '" + str(file_format)
                                  + "' was found in '" + str(data_path) + "'", 2)
    monitoring.to_log_and_console("\t time point range was ["+str(first_time_point)+", "+str(last_time_point)+"]")
    return None
