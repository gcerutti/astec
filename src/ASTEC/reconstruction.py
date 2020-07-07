
import os
import imp
import sys

import ace
import common
import CommunFunctions.cpp_wrapping as cpp_wrapping

monitoring = common.Monitoring()


########################################################################################
#
#
#
########################################################################################

class ReconstructionParameters(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):
        #
        #
        #
        self.intensity_transformation = 'Identity'
        self.intensity_enhancement = None
        self.cell_normalization_min_method = 'cellinterior'
        self.cell_normalization_max_method = 'cellborder'

        self.normalization_min_percentile = 0.01
        self.normalization_max_percentile = 0.99
        self.cell_normalization_sigma = 5.0

        #
        # membrane enhancement parameters
        #
        self.ace = ace.AceParameters()

        #
        # registration parameters
        #
        self.registration = []

        self.registration.append(common.RegistrationParameters())
        self.registration[0].prefix = 'linear_registration'
        self.registration[0].pyramid_highest_level = 5
        self.registration[0].pyramid_lowest_level = 3
        self.registration[0].gaussian_pyramid = True
        self.registration[0].transformation_type = 'affine'
        self.registration[0].transformation_estimation_type = 'wlts'
        self.registration[0].lts_fraction = 1.0

        self.registration.append(common.RegistrationParameters())
        self.registration[1].prefix = 'nonlinear_registration_'
        self.registration[1].pyramid_highest_level = 5
        self.registration[1].pyramid_lowest_level = 3
        self.registration[1].gaussian_pyramid = True
        self.registration[1].transformation_type = 'vectorfield'
        self.registration[1].elastic_sigma = 2.0
        self.registration[1].transformation_estimation_type = 'wlts'
        self.registration[1].lts_fraction = 1.0
        self.registration[1].fluid_sigma = 2.0
        #
        #
        #
        self.keep_reconstruction = True

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('ReconstructionParameters')

        print('- intensity_transformation = ' + str(self.intensity_transformation))
        print('- intensity_enhancement = ' + str(self.intensity_enhancement))
        print('- cell_normalization_min_method = ' + str(self.cell_normalization_min_method))
        print('- cell_normalization_max_method = ' + str(self.cell_normalization_max_method))

        print('- normalization_min_percentile = ' + str(self.normalization_min_percentile))
        print('- normalization_max_percentile = ' + str(self.normalization_max_percentile))
        print('- cell_normalization_sigma = ' + str(self.cell_normalization_sigma))

        self.ace.print_parameters()

        for p in self.registration:
            p.print_parameters()

        print('- keep_reconstruction = ' + str(self.keep_reconstruction))
        print("")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('ReconstructionParameters\n')

            logfile.write('- intensity_transformation = ' + str(self.intensity_transformation) + '\n')
            logfile.write('- intensity_enhancement = ' + str(self.intensity_enhancement) + '\n')
            logfile.write('- cell_normalization_min_method = ' + str(self.cell_normalization_min_method) + '\n')
            logfile.write('- cell_normalization_max_method = ' + str(self.cell_normalization_max_method) + '\n')

            logfile.write('- normalization_min_percentile = ' + str(self.normalization_min_percentile) + '\n')
            logfile.write('- normalization_max_percentile = ' + str(self.normalization_max_percentile) + '\n')
            logfile.write('- cell_normalization_sigma = ' + str(self.cell_normalization_sigma) + '\n')

            self.ace.write_parameters(log_file_name)

            for p in self.registration:
                p.write_parameters(log_file_name)

            logfile.write('- keep_reconstruction = ' + str(self.keep_reconstruction) + '\n')

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
        # reconstruction method
        #

        if hasattr(parameters, 'intensity_transformation'):
            self.intensity_transformation = parameters.intensity_transformation
        if hasattr(parameters, 'mars_intensity_transformation'):
            self.intensity_transformation = parameters.mars_intensity_transformation
        if hasattr(parameters, 'astec_intensity_transformation'):
            self.intensity_transformation = parameters.astec_intensity_transformation

        if hasattr(parameters, 'intensity_enhancement'):
            self.intensity_enhancement = parameters.intensity_enhancement
        if hasattr(parameters, 'mars_intensity_enhancement'):
            self.intensity_enhancement = parameters.mars_intensity_enhancement
        if hasattr(parameters, 'astec_intensity_enhancement'):
            self.intensity_enhancement = parameters.astec_intensity_enhancement

        if hasattr(parameters, 'cell_normalization_min_method'):
            self.cell_normalization_min_method = parameters.cell_normalization_min_method
        if hasattr(parameters, 'astec_cell_normalization_min_method'):
            self.cell_normalization_min_method = parameters.astec_cell_normalization_min_method

        if hasattr(parameters, 'cell_normalization_max_method'):
            self.cell_normalization_max_method = parameters.cell_normalization_max_method
        if hasattr(parameters, 'astec_cell_normalization_max_method'):
            self.cell_normalization_max_method = parameters.astec_cell_normalization_max_method

        if hasattr(parameters, 'normalization_min_percentile'):
            self.normalization_min_percentile = parameters.normalization_min_percentile
        if hasattr(parameters, 'mars_normalization_min_percentile'):
            self.normalization_min_percentile = parameters.mars_normalization_min_percentile
        if hasattr(parameters, 'astec_normalization_min_percentile'):
            self.normalization_min_percentile = parameters.astec_normalization_min_percentile

        if hasattr(parameters, 'normalization_max_percentile'):
            self.normalization_max_percentile = parameters.normalization_max_percentile
        if hasattr(parameters, 'mars_normalization_max_percentile'):
            self.normalization_max_percentile = parameters.mars_normalization_max_percentile
        if hasattr(parameters, 'astec_normalization_max_percentile'):
            self.normalization_max_percentile = parameters.astec_normalization_max_percentile

        if hasattr(parameters, 'cell_normalization_sigma'):
            self.cell_normalization_sigma = parameters.cell_normalization_sigma
        if hasattr(parameters, 'astec_cell_normalization_sigma'):
            self.cell_normalization_sigma = parameters.astec_cell_normalization_sigma

        #
        #
        #
        self.ace.update_from_parameters(parameter_file)

        #
        #
        #
        if hasattr(parameters, 'keep_reconstruction'):
            if parameters.keep_reconstruction is not None:
                self.keep_reconstruction = parameters.keep_reconstruction
        if hasattr(parameters, 'mars_keep_reconstruction'):
            if parameters.mars_keep_reconstruction is not None:
                self.keep_reconstruction = parameters.mars_keep_reconstruction
        if hasattr(parameters, 'astec_keep_reconstruction'):
            if parameters.astec_keep_reconstruction is not None:
                self.keep_reconstruction = parameters.astec_keep_reconstruction


#
#
#
#
#

def get_deformation_from_current_to_previous(current_time, experiment, parameters, previous_time):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
    :return:
    """

    proc = 'get_deformation_from_current_to_previous'

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ReconstructionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # image to be registered
    # floating image = image at previous time
    # reference image = image at current time
    #

    current_name = experiment.fusion_dir.get_image_name(current_time)
    current_image = common.find_file(experiment.fusion_dir.get_directory(), current_name, file_type='image',
                                     callfrom=proc, local_monitoring=None, verbose=False)
    if current_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time " + str(current_time), 2)
        return None
    current_image = os.path.join(experiment.fusion_dir.get_directory(), current_image)

    previous_name = experiment.fusion_dir.get_image_name(previous_time)
    previous_image = common.find_file(experiment.fusion_dir.get_directory(), previous_name, file_type='image',
                                      callfrom=proc, local_monitoring=None, verbose=False)
    if previous_image is None:
        monitoring.to_log_and_console("    .. " + proc + " no fused image was found for time " + str(previous_time), 2)
        return None
    previous_image = os.path.join(experiment.fusion_dir.get_directory(), previous_image)

    #
    # auxiliary file names
    #
    affine_image = common.add_suffix(previous_image, "_affine", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                     new_extension=experiment.default_image_suffix)
    affine_trsf = common.add_suffix(previous_image, "_affine", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                    new_extension="trsf")
    vector_image = common.add_suffix(previous_image, "_vector", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                     new_extension=experiment.default_image_suffix)
    vector_trsf = common.add_suffix(previous_image, "_vector", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                    new_extension="trsf")

    if os.path.isfile(vector_trsf) and not monitoring.forceResultsToBeBuilt:
        return vector_trsf

    #
    # linear registration
    #
    common.blockmatching(current_image, previous_image, affine_image, affine_trsf, None, parameters.registration[0])

    if not os.path.isfile(affine_image) or not os.path.isfile(affine_trsf):
        monitoring.to_log_and_console("    .. " + proc + " linear registration failed for time " + str(current_time), 2)
        return None

    #
    # non-linear registration
    #
    common.blockmatching(current_image, previous_image, vector_image, vector_trsf, affine_trsf,
                         parameters.registration[1])
    if not os.path.isfile(vector_image) or not os.path.isfile(vector_trsf):
        monitoring.to_log_and_console("    .. " + proc + " non-linear registration failed for time " +
                                      str(current_time), 2)
        return None

    return vector_trsf

#
#
#
#
#


def get_previous_deformed_segmentation(current_time, experiment, parameters, previous_time=None):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
    :return:
    """

    #
    #
    #

    proc = 'get_previous_deformed_segmentation'

    #
    # it will generate (in the temporary_path directory)
    #
    # files related to the deformation computation
    # - $EN_fuse_t$TIME_affine.inr
    # - $EN_fuse_t$TIME_affine.trsf
    # - $EN_fuse_t$TIME_vector.inr
    # - $EN_fuse_t$TIME_vector.trsf
    #
    # the deformed segmentation
    # - $EN_[mars,seg]_t$TIME_deformed.inr
    #

    prev_segimage = experiment.get_segmentation_image(previous_time)
    if prev_segimage is None:
        if previous_time is None:
            monitoring.to_log_and_console("    .. " + proc + ": was called with 'previous_time = None'", 2)
        else:
            monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                          + str(previous_time), 2)
        return None

    prev_def_segimage = common.add_suffix(prev_segimage, "_deformed", new_dirname=experiment.working_dir.get_tmp_directory(0),
                                               new_extension=experiment.default_image_suffix)

    if os.path.isfile(prev_def_segimage):
        return prev_def_segimage

    deformation = get_deformation_from_current_to_previous(current_time, experiment, parameters, previous_time)

    if deformation is None:
        monitoring.to_log_and_console("    .. " + proc + ": no deformation was found for time "
                                      + str(current_time) + " towards " + str(previous_time), 2)
        return None

    monitoring.to_log_and_console("    .. resampling of '" + str(prev_segimage).split(os.path.sep)[-1] + "'", 2)

    cpp_wrapping.apply_transformation(prev_segimage, prev_def_segimage, deformation,
                                      interpolation_mode='nearest', monitoring=monitoring)

    return prev_def_segimage


########################################################################################
#
#
#
########################################################################################


def build_membrane_image(current_time, experiment, parameters, previous_time=None):
    """

    :param current_time:
    :param experiment:
    :param parameters:
    :param previous_time:
    :return:
    """
    #
    # build input image(s) for segmentation
    # 0. build names
    # 1. do image enhancement if required
    # 2. transform input image if required
    # 3. mix the two results
    #
    # If any transformation is performed on the input image, the final result image (input of the watershed)
    # will be named (in the 'temporary_path' directory)
    # - 'input_image'_membrane
    # if fusion is required
    # - 'input_image'_enhanced
    #   is the membrance-enhanced image (by tensor voting method)
    # - 'input_image'_intensity
    #   is the intensity-transformed image (by normalization)
    #

    proc = "build_membrane_image"

    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, ReconstructionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)
        
    #
    # get fused image
    #
    input_dir = experiment.fusion_dir.get_directory(0)
    input_name = experiment.fusion_dir.get_image_name(current_time)

    input_image = common.find_file(input_dir, input_name, file_type='image', callfrom=proc, local_monitoring=monitoring)

    if input_image is None:
        monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '" + str(input_dir) + "'", 2)
        return None
    input_image = os.path.join(input_dir, input_image)

    #
    # set image names
    #
    # there can be intensity enhancement (ie filtering based membrane enhancement)
    #
    # if no intensity (i.e. membrane) enhancement
    #    if no intensity
    #
    #
    enhanced_image = None
    intensity_image = None
    membrane_image = None

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement.lower() == 'none':

        #
        # no membrane enhancement image
        # only set intensity_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none' \
                or parameters.intensity_transformation.lower() == 'identity':
            #
            # nothing to do
            #
            intensity_image = input_image
        elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
            intensity_image = common.add_suffix(input_image, "_membrane",
                                                new_dirname=experiment.working_dir.get_rec_directory(0),
                                                new_extension=experiment.default_image_suffix)
            #
            # 'reconstruction' directory has to be made
            #
            experiment.working_dir.make_rec_directory()
        else:
            monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)

    elif parameters.intensity_enhancement.lower() == 'gace' or parameters.intensity_enhancement.lower() == 'glace':
        #
        # only set enhanced_image
        # or set the 3 names: intensity_image, enhanced_image, and membrane_image
        #
        if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none':
            enhanced_image = common.add_suffix(input_image, "_membrane",
                                               new_dirname=experiment.working_dir.get_rec_directory(0),
                                               new_extension=experiment.default_image_suffix)
        elif parameters.intensity_transformation.lower() == 'identity':
            intensity_image = input_image
            enhanced_image = common.add_suffix(input_image, "_enhanced",
                                               new_dirname=experiment.working_dir.get_tmp_directory(0),
                                               new_extension=experiment.default_image_suffix)
            membrane_image = common.add_suffix(input_image, "_membrane",
                                               new_dirname=experiment.working_dir.get_rec_directory(0),
                                               new_extension=experiment.default_image_suffix)
        elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'global_normalization_to_u8' \
                or parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
            intensity_image = common.add_suffix(input_image, "_intensity",
                                                new_dirname=experiment.working_dir.get_tmp_directory(0),
                                                new_extension=experiment.default_image_suffix)
            enhanced_image = common.add_suffix(input_image, "_enhanced",
                                               new_dirname=experiment.working_dir.get_tmp_directory(0),
                                               new_extension=experiment.default_image_suffix)
            membrane_image = common.add_suffix(input_image, "_membrane",
                                               new_dirname=experiment.working_dir.get_rec_directory(0),
                                               new_extension=experiment.default_image_suffix)
        else:
            monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                          + str(parameters.intensity_transformation) + "'", 2)
            return None
        #
        # 'reconstruction' directory has to be made
        #
        experiment.working_dir.make_rec_directory()
    else:
        monitoring.to_log_and_console("    unknown membrane enhancement method: '"
                                      + str(parameters.intensity_enhancement) + "'", 2)
        return None

    #
    # computation of membrane enhancement image, if required
    #

    if parameters.intensity_enhancement is None or parameters.intensity_enhancement.lower() == 'none':
        pass
    elif parameters.intensity_enhancement.lower() == 'gace':
        if (not os.path.isfile(enhanced_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. membrane global enhancement of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            ace.monitoring.copy(monitoring)
            ace.global_membrane_enhancement(input_image, enhanced_image, experiment,
                                            temporary_path=experiment.working_dir.get_tmp_directory(0),
                                            parameters=parameters.ace)
    elif parameters.intensity_enhancement.lower() == 'glace':
        if (not os.path.isfile(enhanced_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. membrane cell enhancement of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            ace.monitoring.copy(monitoring)
            previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, experiment, parameters,
                                                                                previous_time)
            if previous_deformed_segmentation is None:
                monitoring.to_log_and_console("    .. " + proc + ": switch to 'ace' ", 1)
                ace.global_membrane_enhancement(input_image, enhanced_image, experiment,
                                                temporary_path=experiment.working_dir.get_tmp_directory(0),
                                                parameters=parameters.ace)
            else:
                ace.cell_membrane_enhancement(input_image, previous_deformed_segmentation, enhanced_image, experiment,
                                              temporary_path=experiment.working_dir.get_tmp_directory(0),
                                              parameters=parameters.ace)
    else:
        monitoring.to_log_and_console("    unknown membrane enhancement method: '"
                                      + str(parameters.intensity_enhancement) + "'", 2)
        return None

    #
    # computation of transformed intensity image if required
    #

    arit_options = "-o 2"

    if parameters.intensity_transformation is None or parameters.intensity_transformation.lower() == 'none':
        pass
    elif parameters.intensity_transformation.lower() == 'identity':
        pass
    elif parameters.intensity_transformation.lower() == 'normalization_to_u8' \
            or parameters.intensity_transformation.lower() == 'global_normalization_to_u8':
        if (not os.path.isfile(intensity_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. intensity global normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            cpp_wrapping.global_normalization_to_u8(input_image, intensity_image,
                                                    min_percentile=parameters.normalization_min_percentile,
                                                    max_percentile=parameters.normalization_max_percentile,
                                                    other_options=None, monitoring=monitoring)
        arit_options = "-o 1"
    elif parameters.intensity_transformation.lower() == 'cell_normalization_to_u8':
        if (not os.path.isfile(intensity_image) and (membrane_image is None or not os.path.isfile(membrane_image))) \
                or monitoring.forceResultsToBeBuilt is True:
            monitoring.to_log_and_console("    .. intensity cell-based normalization of '"
                                          + str(input_image).split(os.path.sep)[-1] + "'", 2)
            if previous_time is None:
                monitoring.to_log_and_console("       previous time point was not given", 2)
                monitoring.to_log_and_console("    .. " + proc + ": switch to 'global normalization' ", 1)
                cpp_wrapping.global_normalization_to_u8(input_image, intensity_image,
                                                        min_percentile=parameters.normalization_min_percentile,
                                                        max_percentile=parameters.normalization_max_percentile,
                                                        other_options=None, monitoring=monitoring)
            else:
                previous_deformed_segmentation = get_previous_deformed_segmentation(current_time, experiment, parameters,
                                                                                    previous_time)
                cpp_wrapping.cell_normalization_to_u8(input_image, previous_deformed_segmentation, intensity_image,
                                                      min_percentile=parameters.normalization_min_percentile,
                                                      max_percentile=parameters.normalization_max_percentile,
                                                      cell_normalization_min_method=parameters.cell_normalization_min_method,
                                                      cell_normalization_max_method=parameters.cell_normalization_max_method,
                                                      sigma=parameters.cell_normalization_sigma,
                                                      other_options=None, monitoring=monitoring)
        arit_options = "-o 1"
    else:
        monitoring.to_log_and_console("    unknown intensity transformation method: '"
                                      + str(parameters.intensity_transformation) + "'", 2)
        return None

    #
    # fusion of the two images (if fusion is required)
    #

    if enhanced_image is None:
        return intensity_image
    else:
        if intensity_image is None:
            return enhanced_image
        else:
            #
            # mix images
            #
            arit_options += " -max"
            if not os.path.isfile(membrane_image) or monitoring.forceResultsToBeBuilt is True:
                monitoring.to_log_and_console("       fusion of intensity and enhancement", 2)
                cpp_wrapping.arithmetic_operation(intensity_image, enhanced_image, membrane_image,
                                                  other_options=arit_options)
            return membrane_image

    #
    # should not reach this point
    #
