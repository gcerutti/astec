
import os
import imp
import sys

import commonTools
import nomenclature
import lineage

from CommunFunctions.ImageHandling import imread

#
#
#
#
#

monitoring = commonTools.Monitoring()


########################################################################################
#
# classes
# - computation environment
# - computation parameters
#
########################################################################################


class AstecEnvironment(object):

    def __init__(self):

        #
        # segmentation data paths
        #
        self.path_seg_exp = None
        self.path_seg_exp_files = None
        self.path_seg_exp_lineage = None

        #
        #
        #
        self.path_reconstruction = None
        self.temporary_path = None

        #
        #
        #
        self.path_logdir = None
        self.path_history_file = None
        self.path_log_file = None

    def update_from_file(self, parameter_file, start_time):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        self.path_seg_exp = nomenclature.replaceFlags(nomenclature.path_seg_exp, parameters)
        self.path_seg_exp_files = nomenclature.replaceFlags(nomenclature.path_seg_exp_files, parameters)
        self.path_seg_exp_lineage = nomenclature.replaceFlags(nomenclature.path_seg_exp_lineage, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_seg_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_seg_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_seg_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecEnvironment\n')

            logfile.write('- path_seg_exp = ' + str(self.path_seg_exp) + '\n')
            logfile.write('- path_seg_exp_files = ' + str(self.path_seg_exp_files) + '\n')
            logfile.write('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage) + '\n')

            logfile.write('- path_reconstruction = ' + str(self.path_reconstruction) + '\n')
            logfile.write('- temporary_path = ' + str(self.temporary_path) + '\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecEnvironment')

        print('- path_seg_exp = ' + str(self.path_seg_exp))
        print('- path_seg_exp_files = ' + str(self.path_seg_exp_files))
        print('- path_seg_exp_lineage = ' + str(self.path_seg_exp_lineage))

        print('- path_reconstruction = ' + str(self.path_reconstruction))
        print('- temporary_path = ' + str(self.temporary_path))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class AstecParameters(object):

    def __init__(self):
        #
        #
        #
        self.acquisition_delay = 0

        #
        # images suffixes/formats
        #
        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('AstecParameters\n')

            logfile.write('- acquisition_delay = ' + str(self.acquisition_delay) + '\n')

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = '+str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('AstecParameters')

        print('- acquisition_delay = ' + str(self.acquisition_delay))

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # acquisition parameters
        #

        #
        # images suffixes/formats
        #
        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix


########################################################################################
#
# some internal procedures
#
########################################################################################

########################################################################################
#
# some internal procedures
#
########################################################################################

#
#
#
#
#


def astec_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "astec_control"
    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_seg_exp):
        monitoring.to_log_and_console(proc + ": weird, '" + str(environment.path_seg_exp) + "' does not exists", 1)
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    # re-read the lineage file, if any
    # and check whether any time point should be re-computed
    #

    firstTimePoint = experiment.firstTimePoint + experiment.delayTimePoint
    lastTimePoint = experiment.lastTimePoint + experiment.delayTimePoint

    lineage_tree_information = lineage.read_lineage_tree(environment.path_seg_exp_lineage)

    if len(lineage_tree_information) > 0 and  'lin_tree' in lineage_tree_information:
        monitoring.to_log_and_console("    .. test '" + str(environment.path_seg_exp_lineage) + "'", 1)
        cellat = {}
        for y in lineage_tree_information['lin_tree']:
            t=y/10**4
            if t not in cellat:
                cellat[t] = 1
            else:
                cellat[t] += 1

        restart=-1
        t=firstTimePoint
        while restart==-1 and t<=lastTimePoint:
            #
            # possible time point of segmentation, test if ok
            #
            time_value = t + experiment.deltaTimePoint
            acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))
            segmentation_file = environment.path_seg_exp_files.replace(nomenclature.FLAG_TIME, acquisition_time)
            if not os.path.isfile(os.path.join(environment.path_seg_exp, segmentation_file)):
                monitoring.to_log_and_console("       image '" + segmentation_file + "' not found", 1)
                restart=t
            else:
                if cellat[t] == 0:
                    monitoring.to_log_and_console("       lineage of image '" + segmentation_file + "' not found", 1)
                    restart=t
                else:
                    try :
                        segmentation_image = imread(segmentation_file)
                    except IOError:
                        monitoring.to_log_and_console("       error in image '" + segmentation_file + "'", 1)
                        restart=t
            #
            #
            #

            if restart == -1:
                monitoring.to_log_and_console("       time '" + str(t) + "' seems ok", 1)
            t+=1
        firstTimePoint=restart
        monitoring.to_log_and_console("    " + proc + ": restart computation at time '" + str(firstTimePoint)
                                      + "'", 1)
    else:
        monitoring.to_log_and_console("    " + proc + ": start computation at time '" + str(firstTimePoint)
                                      + "'", 1)
