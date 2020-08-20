
import os
import imp
import sys
import time
import numpy as np
from scipy import ndimage as nd

import common
import reconstruction

import CommunFunctions.cpp_wrapping as cpp_wrapping
from CommunFunctions.ImageHandling import imread, imsave

#
#
#
#
#

monitoring = common.Monitoring()


########################################################################################
#
# classes
# - computation environment
# - computation parameters
#
########################################################################################


#
#
#
#
#

class WatershedParameters(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self, obj=None):
        #
        #
        # gaussian standard deviation: smoothing input image for seed extraction
        # h-value for regional minima extraction
        # threshold for regional minima labeling (has to be in [1, h] range
        # gaussian standard deviation: smoothing input image for watershed
        #
        #
        if obj is not None and isinstance(obj, WatershedParameters):
            self.watershed_seed_sigma = obj.watershed_seed_sigma
            self.watershed_seed_hmin = obj.watershed_seed_hmin
            self.watershed_seed_high_threshold = obj.watershed_seed_high_threshold
            self.watershed_membrane_sigma = obj.watershed_membrane_sigma
        else:
            self.watershed_seed_sigma = 0.6
            self.watershed_seed_hmin = 4
            self.watershed_seed_high_threshold = None
            self.watershed_membrane_sigma = 0.15

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('WatershedParameters')

        print('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma))
        print('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin))
        print('- watershed_seed_high_threshold = ' + str(self.watershed_seed_high_threshold))
        print('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma))

        print("")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('WatershedParameters\n')

            logfile.write('- watershed_seed_sigma = ' + str(self.watershed_seed_sigma) + '\n')
            logfile.write('- watershed_seed_hmin = ' + str(self.watershed_seed_hmin) + '\n')
            logfile.write('- watershed_seed_high_threshold = ' + str(self.watershed_seed_high_threshold) + '\n')
            logfile.write('- watershed_membrane_sigma = ' + str(self.watershed_membrane_sigma) + '\n')

            logfile.write("\n")
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # watershed parameters
        #
        if hasattr(parameters, 'mars_sigma1'):
            if parameters.mars_sigma1 is not None:
                self.watershed_seed_sigma = parameters.mars_sigma1
        if hasattr(parameters, 'mars_watershed_seed_sigma'):
            if parameters.mars_watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.mars_watershed_seed_sigma
        if hasattr(parameters, 'astec_sigma1'):
            if parameters.astec_sigma1 is not None:
                self.watershed_seed_sigma = parameters.astec_sigma1
        if hasattr(parameters, 'astec_watershed_seed_sigma'):
            if parameters.astec_watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.astec_watershed_seed_sigma
        if hasattr(parameters, 'watershed_seed_sigma'):
            if parameters.watershed_seed_sigma is not None:
                self.watershed_seed_sigma = parameters.watershed_seed_sigma

        if hasattr(parameters, 'mars_h_min'):
            if parameters.mars_h_min is not None:
                self.watershed_seed_hmin = parameters.mars_h_min
        if hasattr(parameters, 'mars_watershed_seed_hmin'):
            if parameters.mars_watershed_seed_hmin is not None:
                self.watershed_seed_hmin = parameters.mars_watershed_seed_hmin
        if hasattr(parameters, 'watershed_seed_hmin'):
            if parameters.watershed_seed_hmin is not None:
                self.watershed_seed_hmin = parameters.watershed_seed_hmin

        if hasattr(parameters, 'watershed_seed_high_threshold'):
            if parameters.watershed_seed_high_threshold is not None:
                self.watershed_seed_high_threshold = parameters.watershed_seed_high_threshold

        if hasattr(parameters, 'mars_sigma2'):
            if parameters.mars_sigma2 is not None:
                self.watershed_membrane_sigma = parameters.mars_sigma2
        if hasattr(parameters, 'mars_watershed_membrane_sigma'):
            if parameters.mars_watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.mars_watershed_membrane_sigma
        if hasattr(parameters, 'astec_sigma2'):
            if parameters.astec_sigma2 is not None:
                self.watershed_membrane_sigma = parameters.astec_sigma2
        if hasattr(parameters, 'astec_watershed_membrane_sigma'):
            if parameters.astec_watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.astec_watershed_membrane_sigma
        if hasattr(parameters, 'watershed_membrane_sigma'):
            if parameters.watershed_membrane_sigma is not None:
                self.watershed_membrane_sigma = parameters.watershed_membrane_sigma


#
#
#
#
#

class SeedEditionParameters(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):
        #
        #
        # gaussian standard deviation: smoothing input image for seed extraction
        # h-value for regional minima extraction
        # threshold for regional minima labeling (has to be in [1, h] range
        # gaussian standard deviation: smoothing input image for watershed
        #
        #
        self.seed_edition_dir = None
        self.seed_edition_file = None

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('SeedEditionParameters')

        print('- seed_edition_dir = ' + str(self.seed_edition_dir))
        print('- seed_edition_file = ' + str(self.seed_edition_file))

        print("")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('SeedEditionParameters\n')

            logfile.write('- seed_edition_dir = ' + str(self.seed_edition_dir) + '\n')
            logfile.write('- seed_edition_file = ' + str(self.seed_edition_file) + '\n')

            logfile.write("\n")
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # seed edition parameters
        #
        if hasattr(parameters, 'seed_edition_dir'):
            if parameters.seed_edition_dir is not None:
                self.seed_edition_dir = parameters.seed_edition_dir

        if hasattr(parameters, 'seed_edition_file'):
            if parameters.seed_edition_file is not None:
                self.seed_edition_file = parameters.seed_edition_file
        if hasattr(parameters, 'seed_edition_files'):
            if parameters.seed_edition_files is not None:
                self.seed_edition_file = parameters.seed_edition_files

    ############################################################
    #
    # misc
    #
    ############################################################

    def n_seed_editions(self):
        #
        # number of elements
        #
        if self.seed_edition_file is None:
            return 0
        if type(self.seed_edition_file) is str:
            return 1
        if type(self.seed_edition_file) is list:
            if len(self.seed_edition_file) == 2 and type(self.seed_edition_file[0]) is str \
                    and type(self.seed_edition_file[1]) is str:
                return 1
            for e in self.seed_edition_file:
                if type(e) is not list or len(e) != 2 or type(e[0]) is not str or type(e[1]) is not str:
                    return 0
            return len(self.seed_edition_file)
        return 0

    def seed_edition(self, i=0):
        if i >= self.n_seed_editions():
            return None, None
        if type(self.seed_edition_file) is str:
            if self.seed_edition_dir is None:
                return self.seed_edition_file, None
            else:
                return os.path.join(self.seed_edition_dir, self.seed_edition_file), None
        if type(self.seed_edition_file) is list:
            if len(self.seed_edition_file) == 2 and type(self.seed_edition_file[0]) is str \
                    and type(self.seed_edition_file[1]) is str:
                if self.seed_edition_dir is None:
                    return self.seed_edition_file[0], self.seed_edition_file[1]
                else:
                    return os.path.join(self.seed_edition_dir, self.seed_edition_file[0]), \
                           os.path.join(self.seed_edition_dir, self.seed_edition_file[1])
            if self.seed_edition_dir is None:
                return self.seed_edition_file[i][0], self.seed_edition_file[i][1]
            else:
                return os.path.join(self.seed_edition_dir, self.seed_edition_file[i][0]), \
                       os.path.join(self.seed_edition_dir, self.seed_edition_file[i][1])


#
#
#
#
#

class MarsParameters(WatershedParameters, SeedEditionParameters, reconstruction.ReconstructionParameters):

    def __init__(self):

        ############################################################
        #
        # initialisation
        #
        ############################################################

        WatershedParameters.__init__(self)
        SeedEditionParameters.__init__(self)
        reconstruction.ReconstructionParameters.__init__(self)
        #
        #
        #
        self.first_time_point = -1
        self.last_time_point = -1

        return

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('MarsParameters')

        print('- first_time_point = ' + str(self.first_time_point))
        print('- last_time_point = ' + str(self.last_time_point))

        print("")

        WatershedParameters.print_parameters(self)
        SeedEditionParameters.print_parameters(self)
        reconstruction.ReconstructionParameters.print_parameters(self)
        return

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('MarsParameters\n')

            logfile.write('- first_time_point = ' + str(self.first_time_point) + '\n')
            logfile.write('- last_time_point = ' + str(self.last_time_point) + '\n')

        WatershedParameters.write_parameters(self, log_file_name)
        SeedEditionParameters.write_parameters(self, log_file_name)
        reconstruction.ReconstructionParameters.write_parameters(self, log_file_name)

        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        #
        #
        if hasattr(parameters, 'mars_begin'):
            self.first_time_point = parameters.mars_begin
        if hasattr(parameters, 'mars_end'):
            self.last_time_point = parameters.mars_end

        #
        # reconstruction methods
        # backward compatibility
        #

        if hasattr(parameters, 'mars_method'):
            if parameters.mars_method == 1:
                self.intensity_transformation = 'Identity'
                self.intensity_enhancement = None
            elif parameters.mars_method == 2:
                self.intensity_transformation = None
                self.intensity_enhancement = 'GACE'

        #
        # watershed parameters
        # seed edition parameters
        # reconstruction parameters
        #
        WatershedParameters.update_from_parameters(self, parameter_file)
        SeedEditionParameters.update_from_parameters(self, parameter_file)
        reconstruction.ReconstructionParameters.update_from_parameters(self, parameter_file)

        return


########################################################################################
#
# some internal procedures
#
########################################################################################

def build_seeds(input_image, difference_image, output_seed_image, experiment, parameters,
                operation_type='min', check_background_label=False):
    """
    Extract regional minima or maxima from an image and label them.

    :param input_image:
    :param difference_image:
    :param output_seed_image:
    :param experiment:
    :param parameters:
    :param operation_type: 'min' for regional h-minima, 'max' for regional h-maxima
    :param check_background_label: check whether the largest seed (assumed to be the background one) is labeled 1
    :return:
    """

    proc = "build_seeds"

    #
    # variable checking
    #
    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if isinstance(parameters, WatershedParameters) is False:
        monitoring.to_log_and_console(proc + ": bad type for 'parameters' parameter", 1)
        sys.exit(1)

    if operation_type.lower() is 'min' and input_image is None:
        monitoring.to_log_and_console(proc + ": null input image for h-min computation", 1)
        sys.exit(1)

    if operation_type.lower() is 'max' and difference_image is None:
        monitoring.to_log_and_console(proc + ": null difference image for h-max computation", 1)
        sys.exit(1)

    if operation_type.lower() != 'min' and operation_type.lower() != 'max':
        monitoring.to_log_and_console(proc + ": operation type '" + str(operation_type) + "' not handled", 1)
        sys.exit(1)

    #
    # difference_image is the output of 'regional_minima'. It is the difference between
    # the reconstructed image after addition of height h and the original image. This way, valley are transformed
    # into hills. Note that the resulting image has its values in [0, h].
    # Such an image can be used as an input for further regional minima computation with *smaller* h
    #
    # Used by both mars and astec
    #

    if os.path.isfile(output_seed_image) and monitoring.forceResultsToBeBuilt is False:
        return

    if input_image is not None:
        monitoring.to_log_and_console("    .. seed extraction '" + str(input_image).split(os.path.sep)[-1]
                                      + "' with h = " + str(parameters.watershed_seed_hmin), 2)
    else:
        monitoring.to_log_and_console("    .. seed extraction '" + str(difference_image).split(os.path.sep)[-1]
                                      + "' with h = " + str(parameters.watershed_seed_hmin), 2)

    #
    # regional extrema computation
    #
    # the result (of the regional minima) is a difference between the original image and the reconstruction
    # by the same image lowered by hmin. Then the result image is in the range [0, hmin]
    #

    #
    # seed extraction from a membrane image:
    # 1. linear smoothing of input image before regional extrema computation
    # 2. h-minima computation
    #
    # seed extraction from a difference image:
    # h-maxima computation
    #

    if parameters.watershed_seed_sigma > 0.0:
        monitoring.to_log_and_console("       smoothing '" + str(input_image).split(os.path.sep)[-1]
                                      + "' with sigma = " + str(parameters.watershed_seed_sigma), 2)
        seed_preimage = common.add_suffix(input_image, "_smoothing_for_seeds",
                                          new_dirname=experiment.working_dir.get_tmp_directory(0),
                                          new_extension=experiment.default_image_suffix)
        if not os.path.isfile(seed_preimage) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(input_image, seed_preimage, parameters.watershed_seed_sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
    else:
        seed_preimage = input_image

    if not os.path.isfile(seed_preimage):
        monitoring.to_log_and_console("       '" + str(seed_preimage).split(os.path.sep)[-1] + "' does not exist",
                                      2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # name the difference image if required
    #
    hmin = parameters.watershed_seed_hmin
    if difference_image is None:
        local_difference_image = common.add_suffix(seed_preimage, "_seed_diff_h" + str('{:03d}'.format(hmin)),
                                                   new_dirname=experiment.working_dir.get_tmp_directory(0),
                                                   new_extension=experiment.default_image_suffix)
    else:
        local_difference_image = difference_image

    #
    #
    #
    if operation_type.lower() == 'min':
        monitoring.to_log_and_console("       extract regional minima '"
                                      + str(seed_preimage).split(os.path.sep)[-1] + "' with h = " + str(hmin), 2)
        if not os.path.isfile(local_difference_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.regional_minima(seed_preimage, local_difference_image, h=hmin,
                                         monitoring=monitoring)
    else:
        monitoring.to_log_and_console("       extract regional maxima '"
                                      + str(seed_preimage).split(os.path.sep)[-1] + "' with h = " + str(hmin), 2)
        if not os.path.isfile(local_difference_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.regional_maxima(seed_preimage, local_difference_image, h=hmin,
                                         monitoring=monitoring)

    #
    # check whether the computation succeeds
    #

    if not os.path.isfile(local_difference_image):
        monitoring.to_log_and_console("       '" + str(local_difference_image).split(os.path.sep)[-1]
                                      + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # hysteresis thresholding
    # to extract only the extrema that have a height of hmin, the high threshold should be hmin
    # to get more extrema, one can use a smaller hight threshold
    #

    high_threshold = parameters.watershed_seed_hmin
    if parameters.watershed_seed_high_threshold is not None and 0 < parameters.watershed_seed_high_threshold \
            < parameters.watershed_seed_hmin:
        high_threshold = parameters.watershed_seed_high_threshold

    monitoring.to_log_and_console("       label regional extrema '"
                                  + str(local_difference_image).split(os.path.sep)[-1] + "' with threshold = "
                                  + str(high_threshold), 2)

    if not os.path.isfile(output_seed_image) or monitoring.forceResultsToBeBuilt is True:
        cpp_wrapping.connected_components(local_difference_image, output_seed_image, high_threshold=high_threshold,
                                          monitoring=monitoring)

    #
    # check whether the largest connected component has been labeled '1' (background)
    #
    # get the list of labels (remove 0)
    # get the list of volumes
    # get the list of indexes associated with the maximal volume
    #
    if check_background_label:
        seeds = imread(output_seed_image)
        labels = list(np.unique(seeds))
        labels.pop(0)
        volumes = nd.sum(np.ones_like(seeds), seeds, index=np.int16(labels))
        indexmax = [i for i, j in enumerate(volumes) if j == max(volumes)]
        if len(indexmax) > 1:
            monitoring.to_log_and_console("       several regional extrema have a maximal count", 2)
        if int(labels[indexmax[0]]) is not 1:
            monitoring.to_log_and_console("       relabel seed #" + str(labels[indexmax[0]]) + " into 1 (background)", 2)
            newlabel = max(labels) + 1
            seeds[seeds == 1] = newlabel
            seeds[seeds == labels[indexmax[0]]] = 1
            imsave(output_seed_image, seeds)
        del seeds

    if difference_image is None:
        os.remove(local_difference_image)

    return


def watershed(seed_image, membrane_image, result_image, experiment, parameters):
    """

    :param seed_image:
    :param membrane_image:
    :param result_image:
    :param experiment:
    :param parameters:
    :return:
    """

    proc = "watershed"

    #
    # variable checking
    #
    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, WatershedParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    if parameters.watershed_membrane_sigma > 0.0:
        monitoring.to_log_and_console("    .. smoothing '" + str(membrane_image).split(os.path.sep)[-1] +
                                      "' with sigma = " + str(parameters.watershed_membrane_sigma), 2)
        height_image = common.add_suffix(membrane_image, "_smoothing_for_watershed",
                                         new_dirname=experiment.working_dir.get_tmp_directory(0),
                                         new_extension=experiment.default_image_suffix)
        if not os.path.isfile(height_image) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.linear_smoothing(membrane_image, height_image, parameters.watershed_membrane_sigma,
                                          real_scale=True, other_options="-o 2", monitoring=monitoring)
        if not os.path.isfile(height_image):
            monitoring.to_log_and_console("       smoothing '" + str(membrane_image).split(os.path.sep)[-1] +
                                          "' seems to have failed, use non-smoothed image", 2)
            height_image = membrane_image
    else:
        height_image = membrane_image

    #
    # watershed
    #

    if not os.path.isfile(seed_image):
        monitoring.to_log_and_console("       '" + str(seed_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)
    if not os.path.isfile(height_image):
        monitoring.to_log_and_console("       '" + str(height_image).split(os.path.sep)[-1] + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    monitoring.to_log_and_console("    .. watershed '" + str(height_image).split(os.path.sep)[-1] + "'", 2)

    if not os.path.isfile(result_image) or monitoring.forceResultsToBeBuilt is True \
            or (isinstance(parameters, MarsParameters) and parameters.n_seed_editions()) > 0:
        cpp_wrapping.watershed(seed_image, height_image, result_image, monitoring=monitoring)

    return


def _seed_correction(seed_image, corrected_seed_image, parameters):
    """

    :param seed_image:
    :param corrected_seed_image:
    :param parameters:
    :return:
    """

    proc = "_seed_correction"

    if isinstance(parameters, SeedEditionParameters) is False:
        monitoring.to_log_and_console(proc + ": bad type for 'parameters' parameter", 1)
        sys.exit(1)

    if parameters.n_seed_editions() == 0:
        return seed_image

    ifile = seed_image
    ofile = corrected_seed_image
    for i in range(parameters.n_seed_editions()):
        fusion, seeds = parameters.seed_edition(i)
        monitoring.to_log_and_console("       correction ['" + str(fusion).split(os.path.sep)[-1] + "', '"
                                      + str(seeds).split(os.path.sep)[-1] + "']", 2)
        cpp_wrapping.mc_seed_edit(ifile, ofile, fusion, seeds, monitoring=monitoring)
        ifile = ofile

    return corrected_seed_image


def _mars_watershed(template_image, membrane_image, mars_image, experiment, parameters):
    """

    :param template_image: fused image name, to name the other images after it
    :param membrane_image: input image for processing
    :param mars_image:
    :param experiment:
    :param parameters:
    :return:
    """

    proc = "_mars_watershed"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    #
    # computation of seed image
    # - [smoothing]
    # - regional minima
    # - hysteresis thresholding
    #

    seed_image = common.add_suffix(template_image, "_seed_h" + str('{:03d}'.format(parameters.watershed_seed_hmin)),
                                   new_dirname=experiment.working_dir.get_tmp_directory(0),
                                   new_extension=experiment.default_image_suffix)

    if not os.path.isfile(seed_image) or monitoring.forceResultsToBeBuilt is True:
        build_seeds(membrane_image, None, seed_image, experiment, parameters, check_background_label=True)

    #
    # seed correction (if any)
    #

    corrected_seed_image = common.add_suffix(template_image, "_seed_h"
                                             + str('{:03d}'.format(parameters.watershed_seed_hmin)) + "_corrected",
                                             new_dirname=experiment.working_dir.get_tmp_directory(0),
                                             new_extension=experiment.default_image_suffix)
    result_seed_image = _seed_correction(seed_image, corrected_seed_image, parameters)

    #
    #
    #
    watershed(result_seed_image, membrane_image, mars_image, experiment, parameters)

    return


def _volume_diagnosis(mars_image, ncells=10):
    """

    :param mars_image:
    :param ncells:
    :return:
    """

    proc = "_volume_diagnosis"

    #
    # variable checking
    #
    if not os.path.isfile(mars_image):
        monitoring.to_log_and_console("    "+proc+": error, '"+str(mars_image)+"' was not found", 2)
        return

    #
    #
    #
    image = imread(mars_image)
    labels = np.unique(image)
    volumes = nd.sum(np.ones_like(image), image, index=np.int16(labels))
    list_for_sort = list()
    for i in range(len(labels)):
        list_for_sort.append([volumes[i], labels[i]])

    #
    # statistics without the background
    #
    m = np.mean(volumes[1:])
    s = np.std(volumes[1:])

    list_for_sort.sort()

    monitoring.to_log_and_console("    .. diagnosis on cell volumes, smallest cells to be looked at", 1)
    monitoring.to_log_and_console("       mean cell volume = " + str(m) + ", standard deviation = " + str(s), 1)
    for i in range(len(labels)):
        if i <= ncells or list_for_sort[i][0] <= m - 2*s:
            monitoring.to_log_and_console('       cell #'+'{:3d}'.format(list_for_sort[i][1])+' volume='
                                          + '{:10.1f}'.format(list_for_sort[i][0]), 1)

    del image
    return

########################################################################################
#
#
#
########################################################################################


#
#
#
#
#

def mars_process(current_time, experiment, parameters):
    """
    MARS segmentation of one timepoint image.
    :param current_time:
    :param experiment:
    :param parameters:
    :return:
    """

    proc = "mars_process"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, MarsParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # nothing to do if the segmentation image exists
    #
    mars_dir = experiment.mars_dir.get_directory(0)
    mars_name = experiment.mars_dir.get_image_name(current_time)
    mars_image = common.find_file(mars_dir, mars_name, file_type='image', callfrom=proc, local_monitoring=None,
                                  verbose=False)

    if mars_image is not None:
        mars_image = os.path.join(mars_dir, mars_image)
        if monitoring.forceResultsToBeBuilt is False and parameters.n_seed_editions() == 0:
            monitoring.to_log_and_console('    mars image already existing', 2)
            #
            # compute diagnosis anyway
            #
            _volume_diagnosis(mars_image)
            return
        else:
            monitoring.to_log_and_console('    mars image already existing, but forced', 2)
    else:
        mars_image = os.path.join(mars_dir, mars_name + '.' + experiment.result_image_suffix)

    #
    #
    #

    input_dir = experiment.fusion_dir.get_directory(0)
    input_name = experiment.fusion_dir.get_image_name(current_time)

    input_image = common.find_file(input_dir, input_name, file_type='image', callfrom=proc, local_monitoring=monitoring)

    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '" + str(input_dir) + "'", 2)
        monitoring.to_log_and_console("       skip time " + str(current_time), 2)
        return

    #
    # common.find_file() return the file name in the directory
    #
    input_image = os.path.join(input_dir, input_image)

    #
    # build the 'membrane' image to be segmented
    # this 'membrane' image is computed from the input image:
    # - it may be done by intensity normalisation and/or
    # - membrane extraction
    #

    reconstruction.monitoring.copy(monitoring)
    membrane_image = reconstruction.build_membrane_image(current_time, experiment, parameters)

    if membrane_image is None or not os.path.isfile(membrane_image):
        monitoring.to_log_and_console("       '" + str(membrane_image).split(os.path.sep)[-1]
                                      + "' does not exist", 2)
        monitoring.to_log_and_console("\t Exiting.")
        sys.exit(1)

    #
    # compute the seeded watershed
    #

    _mars_watershed(input_image, membrane_image, mars_image, experiment, parameters)

    #
    #
    #
    _volume_diagnosis(mars_image)

    return


#
#
#
#
#


def mars_control(experiment, parameters):
    """

    :param experiment: class describing the experiment
    :param parameters: class describing the mars parameters
    :return:
    """

    proc = "mars_control"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, MarsParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # make sure that the result directory exists
    #

    experiment.mars_dir.make_directory()

    monitoring.to_log_and_console('', 1)

    #
    #
    #

    if parameters.first_time_point < 0 or parameters.last_time_point < 0:
        monitoring.to_log_and_console("... time interval does not seem to be defined in the parameter file")
        monitoring.to_log_and_console("    set parameters 'begin' and 'end'")
        monitoring.to_log_and_console("\t Exiting")
        sys.exit(1)

    if parameters.first_time_point > parameters.last_time_point:
        monitoring.to_log_and_console("... weird time interval: 'begin' = " + str(parameters.first_time_point)
                                      + ", 'end' = " + str(parameters.last_time_point))

    for time_value in range(parameters.first_time_point + experiment.delay_time_point,
                            parameters.last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

        acquisition_time = experiment.working_dir.timepoint_to_str(time_value)

        #
        # start processing
        #
        monitoring.to_log_and_console('... mars processing of time ' + acquisition_time, 1)
        start_time = time.time()

        #
        # set and make temporary directory
        # - there is one temporary directory per timepoint
        #

        experiment.mars_dir.set_tmp_directory(time_value)
        experiment.mars_dir.make_tmp_directory()

        if parameters.keep_reconstruction is False:
            experiment.mars_dir.set_rec_directory_to_tmp()

        #
        # processing
        #

        mars_process(time_value, experiment, parameters)

        #
        # cleaning
        #

        if monitoring.keepTemporaryFiles is False:
            experiment.mars_dir.rmtree_tmp_directory()

        #
        # end processing for a time point
        #
        end_time = time.time()
        monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
        monitoring.to_log_and_console('', 1)

    return
