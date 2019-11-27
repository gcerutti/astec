import os
import imp
import shutil
import sys
import time

import commonTools
from CommunFunctions.ImageHandling import imread
import CommunFunctions.cpp_wrapping as cpp_wrapping

#
#
#
#
#

monitoring = commonTools.Monitoring()


########################################################################################
#
# classes
# - computation parameters
#
########################################################################################


class IntraRegParameters(object):

    def __init__(self):

        #
        # Co-registration parameters
        #

        self.registration = commonTools.RegistrationParameters()
        self.registration.transformation_type = 'rigid'
        self.registration.prefix = 'intra_registration_'

        #
        # intra-sequence transformation parameters
        # reference, ie image to be still (up to a translation)
        # while composing transformation
        # 'reference_transformation_file' is the transformation
        # that will be resampled the whole series afterwards
        # (defined for the reference image)
        #

        self.reference_index = None
        self.reference_transformation_file = None
        self.reference_transformation_angles = None

        #
        # intra-sequence transformation parameters
        # input template, ie how to define the useful information to be kept
        #
        self.template_type = "FUSION"
        self.template_threshold = None
        self.margin = None

        #
        # output template
        #
        self.resolution = 0.6

        #
        # force rebuilding of template and of transformations versus a reference
        # useful when transformations have already been computed for fused image as template
        # they can re-used for segmentation images as template
        #
        self.rebuild_template = False

        #
        # resampling parameters
        #
        self.sigma_segmentation_images = 1.0

        self.resample_fusion_images = True
        self.resample_segmentation_images = False
        self.resample_post_segmentation_images = False

        self.movie_fusion_images = True
        self.movie_segmentation_images = False
        self.movie_post_segmentation_images = False

        self.xy_movie_fusion_images = []
        self.xz_movie_fusion_images = []
        self.yz_movie_fusion_images = []

        self.xy_movie_segmentation_images = []
        self.xz_movie_segmentation_images = []
        self.yz_movie_segmentation_images = []

        self.xy_movie_post_segmentation_images = []
        self.xz_movie_post_segmentation_images = []
        self.yz_movie_post_segmentation_images = []

        #
        # images suffixes/formats
        #

        self.result_image_suffix = 'inr'
        self.default_image_suffix = 'inr'

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('IntraRegParameters\n')

            #
            # Co-registration parameters
            #

            self.registration.write_parameters(log_file_name)

            #
            # intra-sequence transformation parameters
            #

            logfile.write('- reference_index = ' + str(self.reference_index) + '\n')
            logfile.write('- reference_transformation_file = ' + str(self.reference_transformation_file) + '\n')
            logfile.write('- reference_transformation_angles = ' + str(self.reference_transformation_angles) + '\n')

            logfile.write('- template_type = ' + str(self.template_type) + '\n')
            logfile.write('- template_threshold = ' + str(self.template_threshold) + '\n')
            logfile.write('- margin = ' + str(self.margin) + '\n')

            logfile.write('- resolution = ' + str(self.resolution) + '\n')

            logfile.write('- rebuild_template = ' + str(self.rebuild_template) + '\n')

            #
            # resampling parameters
            #

            logfile.write('- sigma_segmentation_images = ' + str(self.sigma_segmentation_images) + '\n')
            logfile.write('- resample_fusion_images = ' + str(self.resample_fusion_images) + '\n')
            logfile.write('- resample_segmentation_images = ' + str(self.resample_segmentation_images) + '\n')
            logfile.write('- resample_post_segmentation_images = ' + str(self.resample_post_segmentation_images) + '\n')

            #
            # movie parameters
            #

            logfile.write('- movie_fusion_images = ' + str(self.movie_fusion_images) + '\n')
            logfile.write('- movie_segmentation_images = ' + str(self.movie_segmentation_images) + '\n')
            logfile.write('- movie_post_segmentation_images = ' + str(self.movie_post_segmentation_images) + '\n')

            logfile.write('- xy_movie_fusion_images = ' + str(self.xy_movie_fusion_images) + '\n')
            logfile.write('- xz_movie_fusion_images = ' + str(self.xz_movie_fusion_images) + '\n')
            logfile.write('- yz_movie_fusion_images = ' + str(self.yz_movie_fusion_images) + '\n')

            logfile.write('- xy_movie_segmentation_images = ' + str(self.xy_movie_segmentation_images) + '\n')
            logfile.write('- xz_movie_segmentation_images = ' + str(self.xz_movie_segmentation_images) + '\n')
            logfile.write('- yz_movie_segmentation_images = ' + str(self.yz_movie_segmentation_images) + '\n')

            logfile.write('- xy_movie_post_segmentation_images = ' + str(self.xy_movie_post_segmentation_images) + '\n')
            logfile.write('- xz_movie_post_segmentation_images = ' + str(self.xz_movie_post_segmentation_images) + '\n')
            logfile.write('- yz_movie_post_segmentation_images = ' + str(self.yz_movie_post_segmentation_images) + '\n')

            #
            # images suffixes/formats
            #

            logfile.write('- result_image_suffix = ' + str(self.result_image_suffix) + '\n')
            logfile.write('- default_image_suffix = ' + str(self.default_image_suffix) + '\n')

            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('IntraRegParameters')

        #
        # Co-registration parameters
        #

        self.registration.print_parameters()

        #
        # intra-sequence transformation parameters
        #

        print('- reference_index = ' + str(self.reference_index))
        print('- reference_transformation_file = ' + str(self.reference_transformation_file))
        print('- reference_transformation_angles = ' + str(self.reference_transformation_angles))

        print('- template_type = ' + str(self.template_type))
        print('- template_threshold = ' + str(self.template_threshold))
        print('- margin = ' + str(self.margin))

        print('- resolution = ' + str(self.resolution))

        print('- rebuild_template = ' + str(self.rebuild_template))

        #
        # resampling parameters
        #

        print('- sigma_segmentation_images = ' + str(self.sigma_segmentation_images))
        print('- resample_fusion_images = ' + str(self.resample_fusion_images))
        print('- resample_segmentation_images = ' + str(self.resample_segmentation_images))
        print('- resample_post_segmentation_images = ' + str(self.resample_post_segmentation_images))

        #
        # movie parameters
        #

        print('- movie_fusion_images = ' + str(self.movie_fusion_images))
        print('- movie_segmentation_images = ' + str(self.movie_segmentation_images))
        print('- movie_post_segmentation_images = ' + str(self.movie_post_segmentation_images))

        print('- xy_movie_fusion_images = ' + str(self.xy_movie_fusion_images))
        print('- xz_movie_fusion_images = ' + str(self.xz_movie_fusion_images))
        print('- yz_movie_fusion_images = ' + str(self.yz_movie_fusion_images))

        print('- xy_movie_segmentation_images = ' + str(self.xy_movie_segmentation_images))
        print('- xz_movie_segmentation_images = ' + str(self.xz_movie_segmentation_images))
        print('- yz_movie_segmentation_images = ' + str(self.yz_movie_segmentation_images))

        print('- xy_movie_post_segmentation_images = ' + str(self.xy_movie_post_segmentation_images))
        print('- xz_movie_post_segmentation_images = ' + str(self.xz_movie_post_segmentation_images))
        print('- yz_movie_post_segmentation_images = ' + str(self.yz_movie_post_segmentation_images))

        #
        # images suffixes/formats
        #

        print('- result_image_suffix = ' + str(self.result_image_suffix))
        print('- default_image_suffix = ' + str(self.default_image_suffix))

        print("")

    def update_from_args(self, args):
        self.reference_transformation_file = args.reference_transformation_file
        self.reference_transformation_angles = args.reference_transformation_angles

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        margin_is_updated = False
        parameters = imp.load_source('*', parameter_file)

        #
        # co-registration parameters
        #

        self.registration.update_from_file(parameter_file)

        #
        # intra-sequence transformation parameters
        #

        if hasattr(parameters, 'intra_registration_reference_index'):
            if parameters.intra_registration_reference_index is not None:
                self.reference_index = parameters.intra_registration_reference_index
        if hasattr(parameters, 'intra_registration_reference_resampling_transformation_file'):
            if parameters.intra_registration_reference_resampling_transformation_file is not None:
                self.reference_transformation_file = parameters.intra_registration_reference_resampling_transformation_file
        if hasattr(parameters, 'intra_registration_reference_resampling_transformation_angles'):
            if parameters.intra_registration_reference_resampling_transformation_angles is not None:
                self.reference_transformation_angles = parameters.intra_registration_reference_resampling_transformation_angles

        if hasattr(parameters, 'intra_registration_template_type'):
            if parameters.intra_registration_template_type is not None:
                self.template_type = parameters.intra_registration_template_type
        if hasattr(parameters, 'intra_registration_template_threshold'):
            if parameters.intra_registration_template_threshold is not None:
                self.template_threshold = parameters.intra_registration_template_threshold
        if hasattr(parameters, 'intra_registration_margin'):
            if parameters.intra_registration_margin is not None:
                self.margin = parameters.intra_registration_margin
                margin_is_updated = True

        if hasattr(parameters, 'intra_registration_resolution'):
            if parameters.intra_registration_resolution is not None:
                self.resolution = parameters.intra_registration_resolution

        if hasattr(parameters, 'intra_registration_rebuild_template'):
            if parameters.intra_registration_rebuild_template is not None:
                self.rebuild_template = parameters.intra_registration_rebuild_template

        #
        # resampling parameters
        #

        if hasattr(parameters, 'intra_registration_sigma_segmentation_images'):
            if parameters.intra_registration_sigma_segmentation_images is not None:
                self.sigma_segmentation_images = parameters.intra_registration_sigma_segmentation_images

        if hasattr(parameters, 'intra_registration_resample_fusion_images'):
            if parameters.intra_registration_resample_fusion_images is not None:
                self.resample_fusion_images = parameters.intra_registration_resample_fusion_images
        if hasattr(parameters, 'intra_registration_resample_segmentation_images'):
            if parameters.intra_registration_resample_segmentation_images is not None:
                self.resample_segmentation_images = parameters.intra_registration_resample_segmentation_images
        if hasattr(parameters, 'intra_registration_resample_post_segmentation_images'):
            if parameters.intra_registration_resample_post_segmentation_images is not None:
                self.resample_post_segmentation_images = parameters.intra_registration_resample_post_segmentation_images

        if hasattr(parameters, 'intra_registration_movie_fusion_images'):
            if parameters.intra_registration_movie_fusion_images is not None:
                self.movie_fusion_images = parameters.intra_registration_movie_fusion_images
        if hasattr(parameters, 'intra_registration_movie_segmentation_images'):
            if parameters.intra_registration_movie_segmentation_images is not None:
                self.movie_segmentation_images = parameters.intra_registration_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_movie_post_segmentation_images'):
            if parameters.intra_registration_movie_post_segmentation_images is not None:
                self.movie_post_segmentation_images = parameters.intra_registration_movie_post_segmentation_images

        if hasattr(parameters, 'intra_registration_xy_movie_fusion_images'):
            if parameters.intra_registration_xy_movie_fusion_images is not None:
                self.xy_movie_fusion_images = parameters.intra_registration_xy_movie_fusion_images
        if hasattr(parameters, 'intra_registration_xz_movie_fusion_images'):
            if parameters.intra_registration_xz_movie_fusion_images is not None:
                self.xz_movie_fusion_images = parameters.intra_registration_xz_movie_fusion_images
        if hasattr(parameters, 'intra_registration_yz_movie_fusion_images'):
            if parameters.intra_registration_yz_movie_fusion_images is not None:
                self.yz_movie_fusion_images = parameters.intra_registration_yz_movie_fusion_images

        if hasattr(parameters, 'intra_registration_xy_movie_segmentation_images'):
            if parameters.intra_registration_xy_movie_segmentation_images is not None:
                self.xy_movie_segmentation_images = parameters.intra_registration_xy_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_xz_movie_segmentation_images'):
            if parameters.intra_registration_xz_movie_segmentation_images is not None:
                self.xz_movie_segmentation_images = parameters.intra_registration_xz_movie_segmentation_images
        if hasattr(parameters, 'intra_registration_yz_movie_segmentation_images'):
            if parameters.intra_registration_yz_movie_segmentation_images is not None:
                self.yz_movie_segmentation_images = parameters.intra_registration_yz_movie_segmentation_images

        if hasattr(parameters, 'intra_registration_xy_movie_post_segmentation_images'):
            if parameters.intra_registration_xy_movie_post_segmentation_images is not None:
                self.xy_movie_post_segmentation_images = parameters.intra_registration_xy_movie_post_segmentation_images
        if hasattr(parameters, 'intra_registration_xz_movie_post_segmentation_images'):
            if parameters.intra_registration_xz_movie_post_segmentation_images is not None:
                self.xz_movie_post_segmentation_images = parameters.intra_registration_xz_movie_post_segmentation_images
        if hasattr(parameters, 'intra_registration_yz_movie_post_segmentation_images'):
            if parameters.intra_registration_yz_movie_post_segmentation_images is not None:
                self.yz_movie_post_segmentation_images = parameters.intra_registration_yz_movie_post_segmentation_images

        #
        # images suffixes/formats
        #

        if hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.RESULT_IMAGE_SUFFIX_FUSE
        if hasattr(parameters, 'result_image_suffix'):
            if parameters.result_image_suffix is not None:
                self.result_image_suffix = parameters.result_image_suffix

        if hasattr(parameters, 'default_image_suffix'):
            if parameters.default_image_suffix is not None:
                self.default_image_suffix = parameters.default_image_suffix
                if not hasattr(parameters, 'result_image_suffix') \
                        and not hasattr(parameters, 'RESULT_IMAGE_SUFFIX_FUSE'):
                    self.result_image_suffix = parameters.default_image_suffix

        #
        # default for margin is none (ie no margin)
        # which is convenient for fusion images
        # however, when dealing with segmentation images as template, a little margin will be great
        # thus, define the default as 10 (if no specification from the user)
        #
        if hasattr(parameters, 'intra_registration_template_type'):
            if parameters.intra_registration_template_type.lower() == 'segmentation' \
                    or parameters.intra_registration_template_type.lower() == 'seg' \
                    or parameters.intra_registration_template_type.lower() == 'post-segmentation' \
                    or parameters.intra_registration_template_type.lower() == 'post_segmentation' \
                    or parameters.intra_registration_template_type.lower() == 'post':
                if margin_is_updated is False:
                    parameters.intra_registration_margin = 10


########################################################################################
#
# define paths and file names
#
########################################################################################


def _get_cotrsf_path(experiment):
    """

    :param experiment:
    :return:
    """
    cotrsf_path = os.path.join(experiment.embryo_path, experiment.intrareg.get_directory(), 'CO-TRSFS')
    return cotrsf_path


def _get_cotrsf_path_name(experiment, floating_index, reference_index):
    """

    :param experiment:
    :param floating_index:
    :param reference_index:
    :return:
    """
    flo = experiment.get_time_index(floating_index)
    ref = experiment.get_time_index(reference_index)
    co_trsf_name = experiment.embryoName + '_intrareg_flo' + str(flo) + '_ref' + str(ref) + '.trsf'
    return os.path.join(_get_cotrsf_path(experiment), co_trsf_name)


def _get_cotrsf_path_format(experiment):
    """

    :param experiment:
    :return:
    """
    form = experiment.get_time_format()
    co_trsf_format = experiment.embryoName + '_intrareg_flo' + form + '_ref' + form + '.trsf'
    return os.path.join(_get_cotrsf_path(experiment), co_trsf_format)


def _get_trsf_path(experiment):
    """

    :param experiment:
    :return:
    """
    firstindex = experiment.first_time_point + experiment.delay_time_point
    lastindex = experiment.last_time_point + experiment.delay_time_point
    trsf_dir = 'TRSFS' + "_t" + str(firstindex) + "-" + str(lastindex)
    trsf_path = os.path.join(experiment.embryo_path, experiment.intrareg.get_directory(), trsf_dir)
    return trsf_path


def _get_trsf_name(experiment, index):
    """

    :param experiment:
    :param index:
    :return:
    """
    ind = experiment.get_time_index(index)
    trsf_name = experiment.embryoName + '_intrareg_t' + str(ind) + '.trsf'
    return trsf_name


def _get_trsf_format(experiment):
    """

    :param experiment:
    :return:
    """
    form = experiment.get_time_format()
    trsf_format = experiment.embryoName + '_intrareg_t' + form + '.trsf'
    return trsf_format


def _get_template_path_name(experiment, parameters):
    """

    :param experiment:
    :return:
    """
    firstindex = experiment.first_time_point + experiment.delay_time_point
    lastindex = experiment.last_time_point + experiment.delay_time_point
    result_template = "template" + "_t" + str(firstindex) + "-" + str(lastindex) + "." + parameters.result_image_suffix
    return os.path.join(_get_trsf_path(experiment), result_template)


########################################################################################
#
# some internal procedures
#
########################################################################################

def _check_data(experiment, suffix=None):
    """
    Check whether all the images (from the first time point to the last one) exist
    :param experiment:
    :param suffix:
    :return:
    """

    proc = "_check_data"

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    path_fusion = os.path.join(experiment.embryo_path, experiment.fusion.get_directory())

    for current_time in range(first_time_point + experiment.delta_time_point,
                              last_time_point + 1, experiment.delta_time_point):

        input_name = experiment.get_image_name(current_time, 'fuse')

        if suffix is None:
            input_image = commonTools.find_file(path_fusion, input_name, local_monitoring=monitoring, callfrom=proc)

            if input_image is None:
                monitoring.to_log_and_console("    .. image '" + input_name + "' not found in '"
                                              + str(path_fusion) + "'", 2)
                return False

        else:
            input_image = input_name + '.' + suffix
            if not os.path.isfile(os.path.join(path_fusion, input_image)):
                monitoring.to_log_and_console("    .. image '" + input_image + "' not found in '"
                                              + str(path_fusion) + "'", 2)
                return False

    return True


########################################################################################
#
#
#
########################################################################################


def _coregistration_control(experiment, parameters):
    """
    Perform the co-registration of any couple of two successive images
    Resulting transformations are computed with fused image at t as floating image
    and used image at t+delta_t as reference image

    :param experiment:
    :param parameters:
    :return:
    """

    proc = "_coregistration_control"

    #    if _check_data(experiment) is False:
    #        monitoring.to_log_and_console(proc + ": error, some fused data are missing", 1)
    #        monitoring.to_log_and_console("\t Exiting")
    #        sys.exit(1)

    path_intrareg_cotrsf = _get_cotrsf_path(experiment)
    if not os.path.isdir(path_intrareg_cotrsf):
        os.makedirs(path_intrareg_cotrsf)

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    path_fusion = os.path.join(experiment.embryo_path, experiment.fusion.get_directory())

    #
    # loop on time
    #

    for reference_time in range(first_time_point + experiment.delta_time_point,
                                last_time_point + 1, experiment.delta_time_point):

        floating_time = reference_time - experiment.delta_time_point
        trsf_name = _get_cotrsf_path_name(experiment, floating_time, reference_time)

        if not os.path.isfile(trsf_name) or monitoring.forceResultsToBeBuilt is True:

            floating_prefix = experiment.get_image_name(floating_time, 'fuse')
            floating_name = commonTools.find_file(path_fusion, floating_prefix, local_monitoring=monitoring,
                                                  callfrom=proc)
            if floating_name is None:
                monitoring.to_log_and_console(proc + ": error, image '" + str(floating_prefix) + "' was not found", 1)
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)

            reference_prefix = experiment.get_image_name(reference_time, 'fuse')
            reference_name = commonTools.find_file(path_fusion, reference_prefix, local_monitoring=monitoring,
                                                   callfrom=proc)
            if reference_name is None:
                monitoring.to_log_and_console(proc + ": error, image '" + str(reference_prefix) + "' was not found", 1)
                monitoring.to_log_and_console("\t Exiting")
                sys.exit(1)

            floating_image = os.path.join(path_fusion, floating_name)
            reference_image = os.path.join(path_fusion, reference_name)
            monitoring.to_log_and_console("       co-registering '" + floating_name + "'", 2)
            cpp_wrapping.linear_registration(reference_image, floating_image, None, trsf_name, None,
                                             py_hl=parameters.pyramid_highest_level,
                                             py_ll=parameters.pyramid_lowest_level,
                                             transformation_type=parameters.transformation_type,
                                             transformation_estimator=parameters.transformation_estimation_type,
                                             lts_fraction=parameters.lts_fraction,
                                             normalization=parameters.normalization,
                                             monitoring=monitoring)

    return


def _transformations_from_reference(experiment, parameters, temporary_dir):
    """
    Combine the transformations issued from the co-registration of pairs of successive images
    to get transformations from one given image
    :param experiment:
    :param parameters:
    :return:
    """

    if not os.path.isdir(temporary_dir):
        os.makedirs(temporary_dir)

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    #
    #
    #
    build_trsf = False
    if monitoring.forceResultsToBeBuilt is True:
        build_trsf = True
    else:
        for i in range(first_time_point + experiment.delta_time_point, last_time_point + 1,
                       experiment.delta_time_point):
            trsf_name = _get_trsf_name(experiment, i)
            if not os.path.isfile(os.path.join(temporary_dir, trsf_name)):
                build_trsf = True
                break

    #
    #
    #
    if build_trsf is True:
        if parameters.reference_index is None:
            reference_index = first_time_point
        else:
            reference_index = parameters.reference_index

        format_input = _get_cotrsf_path_format(experiment)

        trsf_format = _get_trsf_format(experiment)
        format_output = os.path.join(temporary_dir, trsf_format)

        cpp_wrapping.multiple_trsfs(format_input, format_output, first_time_point, last_time_point,
                                    reference_index, trsf_type=parameters.registration.transformation_type,
                                    monitoring=monitoring)

    return


def _is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def _get_reference_transformation(parameters, temporary_dir):
    """

    :param parameters:
    :return:
    """

    proc = "_get_reference_transformation"

    if parameters.reference_transformation_file is not None:
        if os.path.isfile(parameters.reference_transformation_file):
            return parameters.reference_transformation_file

    if parameters.reference_transformation_angles is not None:
        if type(parameters.reference_transformation_angles) == str:
            s = parameters.reference_transformation_angles.split(' ')
            create_trsf_option = "-angle-unit degree"
            l = len(create_trsf_option)
            i = 0
            while i < len(s):
                if s[i].lower() == 'x':
                    if i+1 < len(s) and _is_float(s[i+1]):
                        create_trsf_option += ' -xrotation ' + s[i + 1]
                    else:
                        monitoring.to_log_and_console(proc + ": weird value for 'X rotation' -> " + str(s[i + 1]), 1)
                elif s[i].lower() == 'y':
                    if i + 1 < len(s) and _is_float(s[i + 1]):
                        create_trsf_option += ' -yrotation ' + s[i + 1]
                    else:
                        monitoring.to_log_and_console(proc + ": weird value for 'Y rotation' -> " + str(s[i + 1]), 1)
                elif s[i].lower() == 'z':
                    if i + 1 < len(s) and _is_float(s[i + 1]):
                        create_trsf_option += ' -zrotation ' + s[i + 1]
                    else:
                        monitoring.to_log_and_console(proc + ": weird value for 'Z rotation' -> " + str(s[i + 1]), 1)
                i += 2

            if len(create_trsf_option) > l:
                trsf_name = os.path.join(temporary_dir, "reference_transformation.trsf")
                cpp_wrapping.create_trsf(trsf_name, other_options=create_trsf_option, monitoring=monitoring)
                if os.path.isfile(trsf_name):
                    return trsf_name
                else:
                    monitoring.to_log_and_console(proc + ": unable to compute reference transformation", 1)
            else:
                monitoring.to_log_and_console(proc + ": unable to translate reference transformation angles", 1)

    return None


def _transformations_and_template(experiment, parameters, temporary_dir):
    """
    From transformations from one given image, compute the template to resample all images.
    :param experiment:
    :param parameters:
    :param temporary_dir:
    :return:
    """

    proc = "_transformations_and_template"

    if isinstance(parameters, IntraRegParameters) is False:
        monitoring.to_log_and_console(proc + ": bad type for 'parameters' parameter", 1)
        sys.exit(1)

    result_dir = _get_trsf_path(experiment)
    result_template = _get_template_path_name(experiment, parameters)

    if os.path.isfile(result_template) and monitoring.forceResultsToBeBuilt is not True \
            and parameters.rebuild_template is not True:
        return

    monitoring.to_log_and_console("       warning: this stage may be long", 2)

    if not os.path.isdir(result_dir):
        os.makedirs(result_dir)

    #
    # which format to use?
    # if requested, try to use segmentation image (threshold should be set to 2)
    # else use fusion images
    #
    template_format = None

    if parameters.template_type.lower() == 'segmentation' or parameters.template_type.lower() == 'seg':
        #
        # check whether segmentation image share a common suffix
        #
        path_template_format = os.path.join(experiment.embryo_path, experiment.seg.get_directory())
        suffix = commonTools.get_file_suffix(experiment, path_template_format, experiment.get_image_format('seg'),
                                             flag_time=experiment.get_time_format())

        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + experiment.seg.get_directory() + "'", 1)
            monitoring.to_log_and_console("\t switch to fused images as templates", 1)

        else:
            monitoring.to_log_and_console("       ... build template from segmentation images of '"
                                          + experiment.seg.get_directory() + "'", 2)
            template_name = experiment.get_image_format('seg') + '.' + suffix
            template_format = os.path.join(path_template_format, template_name)

    elif parameters.template_type.lower() == 'post-segmentation' \
            or parameters.template_type.lower() == 'post_segmentation' or parameters.template_type.lower() == 'post':
        #
        # check whether post-corrected segmentation image share a common suffix
        #
        path_template_format = os.path.join(experiment.embryo_path, experiment.post.get_directory())
        suffix = commonTools.get_file_suffix(experiment, path_template_format, experiment.get_image_format('post'),
                                             flag_time=experiment.get_time_format())

        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + experiment.post.get_directory() + "'", 1)
            monitoring.to_log_and_console("\t switch to fused images as templates", 1)

        else:
            monitoring.to_log_and_console("       ... build template from post-segmentation images of '"
                                          + experiment.post.get_directory() + "'", 2)
            template_name = experiment.get_image_format('post') + '.' + suffix
            template_format = os.path.join(path_template_format, template_name)

    #
    # use fusion images to build the template
    #
    if template_format is None:
        #
        # check whether fusion image share a common suffix
        #
        path_template_format = os.path.join(experiment.embryo_path, experiment.fusion.get_directory())
        suffix = commonTools.get_file_suffix(experiment, path_template_format, experiment.get_image_format('fuse'),
                                             flag_time=experiment.get_time_format())
        if suffix is None:
            monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                          + experiment.fusion.get_directory() + "'", 1)
            monitoring.to_log_and_console("\t Exiting", 1)
            sys.exit(1)
        else:
            monitoring.to_log_and_console("       ... build template from fusion images of '"
                                          + experiment.fusion.get_directory() + "'", 2)
            template_name = experiment.get_image_format('fuse') + '.' + suffix
            template_format = os.path.join(path_template_format, template_name)

    #
    # other parameters
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    trsf_format = _get_trsf_format(experiment)
    format_input = os.path.join(temporary_dir, trsf_format)
    format_output = os.path.join(result_dir, trsf_format)

    if parameters.reference_index is None:
        reference_index = first_time_point
    else:
        reference_index = parameters.reference_index

    #
    #
    #

    cpp_wrapping.change_multiple_trsfs(format_input, format_output, first_time_point, last_time_point, reference_index,
                                       result_template, trsf_type=parameters.registration.transformation_type,
                                       resolution=parameters.resolution, threshold=parameters.template_threshold,
                                       margin=parameters.margin, format_template=template_format,
                                       reference_transformation=_get_reference_transformation(parameters, temporary_dir),
                                       monitoring=monitoring)

    return


def _resample_images(experiment, parameters, template_image, directory_type, interpolation_mode='linear'):
    """
    resample all images given a set of transformations and a template
    :param experiment:
    :param parameters:
    :param template_image:
    :param directory_type:
    :param interpolation_mode:
    :return:
    """

    proc = "_resample_images"

    #
    # in case the template has been gziped, or copied into an other format
    #

    b = os.path.basename(template_image)
    d = os.path.dirname(template_image)
    local_template_name = commonTools.find_file(d, b, local_monitoring=monitoring, callfrom=proc)
    if local_template_name is None:
        monitoring.to_log_and_console(proc + ": template '" + str(b) + "' was not found in '" + str(d) + "'", 1)
        monitoring.to_log_and_console("\t resampling will not be done")
        return
    local_template_image = os.path.join(d, local_template_name)

    #
    #
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    #
    #
    #

    if directory_type.lower() == 'fuse':
        subdir = experiment.fusion
    elif directory_type.lower() == 'post':
        subdir = experiment.post
    elif directory_type.lower() == 'seg':
        subdir = experiment.seg
    else:
        monitoring.to_log_and_console(proc + ": unknown directory type '" + str(directory_type) + "'", 1)
        monitoring.to_log_and_console("\t resampling will not be done")
        return

    #
    # loop on directories
    #

    trsf_dir = _get_trsf_path(experiment)

    for idir in range(subdir.get_number_directories()):

        dir_input = subdir.get_directory(idir)
        monitoring.to_log_and_console("     . resampling '" + str(dir_input) + "'", 2)
        dir_input = os.path.join(experiment.embryo_path, subdir.get_directory(idir))
        dir_output = os.path.join(experiment.embryo_path, experiment.intrareg.get_directory(),
                                  subdir.get_directory(idir))
        if not os.path.isdir(dir_output):
            os.makedirs(dir_output)

        #
        # loop on images
        #

        for t in range(first_time_point + experiment.delay_time_point,
                       last_time_point + experiment.delay_time_point + 1, experiment.delta_time_point):

            input_name = experiment.get_image_name(t, directory_type)
            output_name = experiment.get_image_name(t, 'intrareg', directory_type)

            output_image = commonTools.find_file(dir_output, output_name, local_monitoring=None, verbose=False,
                                                 callfrom=proc)

            if output_image is None or monitoring.forceResultsToBeBuilt is True or parameters.rebuild_template is True:

                input_image = commonTools.find_file(dir_input, input_name, local_monitoring=monitoring, callfrom=proc)
                output_image = os.path.join(dir_output, output_name + '.' + str(parameters.result_image_suffix))
                trsf_name = os.path.join(trsf_dir, _get_trsf_name(experiment, t))

                if input_image is None:
                    monitoring.to_log_and_console(proc + ": image '" + str(input_name) + "' was not found in '"
                                                  + str(dir_input) + "'", 1)
                elif not os.path.isfile(trsf_name):
                    monitoring.to_log_and_console(proc + ": transformation '" + str(_get_trsf_name(experiment, t))
                                                  + "' was not found in '" + str(trsf_dir) + "'", 1)
                else:
                    monitoring.to_log_and_console("       resampling '" + str(input_image) + "'", 2)
                    input_image = os.path.join(dir_input, input_image)
                    cpp_wrapping.apply_transformation(input_image, output_image, trsf_name, local_template_image,
                                                      interpolation_mode=interpolation_mode,
                                                      cell_based_sigma=parameters.sigma_segmentation_images,
                                                      monitoring=monitoring)

    return


def _make_movies(experiment, parameters, directory_type, xylist, xzlist, yzlist):
    """

    :param experiment:
    :param parameters:
    :param directory_type:
    :param xylist:
    :param xzlist:
    :param yzlist:
    :return:
    """

    proc = "_make_movies"

    #
    #
    #

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    #
    #
    #

    if directory_type.lower() == 'fuse':
        subdir = experiment.fusion
    elif directory_type.lower() == 'post':
        subdir = experiment.post
    elif directory_type.lower() == 'seg':
        subdir = experiment.seg
    else:
        monitoring.to_log_and_console(proc + ": unknown directory type '" + str(directory_type) + "'", 1)
        monitoring.to_log_and_console("\t resampling will not be done")
        return

    #
    # should we switch to default behavior?
    #

    xy = []
    xz = []
    yz = []

    if len(xylist) > 0 or len(xzlist) > 0 or len(yzlist) > 0:
        xy = xylist
        xz = xzlist
        yz = yzlist

    #
    # loop on directories
    #

    for idir in range(subdir.get_number_directories()):

        dir_input = os.path.join(experiment.intrareg.get_directory(), subdir.get_directory(idir))
        monitoring.to_log_and_console("     . movies from '" + str(dir_input) + "'", 2)
        dir_input = os.path.join(experiment.embryo_path, experiment.intrareg.get_directory(),
                                 subdir.get_directory(idir))
        dir_output = os.path.join(experiment.embryo_path, experiment.intrareg.get_directory(), 'MOVIES',
                                  subdir.get_directory(idir))
        if not os.path.isdir(dir_output):
            os.makedirs(dir_output)

        #
        # default behavior
        #
        if len(xy) == 0 and len(xz) == 0 and len(yz) == 0:
            #
            # read the first image and set XY slice in the middle
            #
            first_prefix = experiment.get_image_name(first_time_point, 'intrareg', directory_type)
            first_name = commonTools.find_file(dir_input, first_prefix, local_monitoring=None, verbose=False,
                                               callfrom=proc)
            if first_name is None:
                monitoring.to_log_and_console(proc + ": no file '" + str(first_prefix) + "' in '" + str(dir_input) + "'", 1)
                monitoring.to_log_and_console("\t movies will not be done")
                return
            first_image = imread(os.path.join(dir_input, first_name))
            xy.append(int(first_image.shape[2] / 2))
            del first_image

        input_format = experiment.get_image_format('intrareg', directory_type)
        input_format = os.path.join(dir_input, input_format + '.' + str(parameters.result_image_suffix))

        #
        # processing
        #

        if len(xy) > 0:
            for s in xy:
                name_output = experiment.get_movie_name(first_time_point, last_time_point, 'intrareg', directory_type)
                name_output += '_xy' + '{:0{width}d}'.format(s, width=4)
                name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
                if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True \
                        or parameters.rebuild_template is True:
                    monitoring.to_log_and_console("       process xy=" + str(s), 2)
                    cpp_wrapping.crop_sequence(input_format, name_output, first_time_point, last_time_point, 'xy', s,
                                               monitoring=monitoring)

        if len(xz) > 0:
            for s in xz:
                name_output = experiment.get_movie_name(first_time_point, last_time_point, 'intrareg', directory_type)
                name_output += '_xz' + '{:0{width}d}'.format(s, width=4)
                name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
                if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True \
                        or parameters.rebuild_template is True:
                    monitoring.to_log_and_console("       process xz=" + str(s), 2)
                    cpp_wrapping.crop_sequence(input_format, name_output, first_time_point, last_time_point, 'xz', s,
                                               monitoring=monitoring)

        if len(yz) > 0:
            for s in yz:
                name_output = experiment.get_movie_name(first_time_point, last_time_point, 'intrareg', directory_type)
                name_output += '_yz' + '{:0{width}d}'.format(s, width=4)
                name_output = os.path.join(dir_output, name_output + '.' + str(parameters.result_image_suffix))
                if os.path.isfile(name_output) is False or monitoring.forceResultsToBeBuilt is True \
                        or parameters.rebuild_template is True:
                    monitoring.to_log_and_console("       process yz=" + str(s), 2)
                    cpp_wrapping.crop_sequence(input_format, name_output, first_time_point, last_time_point, 'yz', s,
                                               monitoring=monitoring)

    return


########################################################################################
#
#
#
########################################################################################


def intraregistration_control(experiment, parameters):
    """

    :param experiment:
    :param parameters:
    :return:
    """

    proc = "intraregistration_control"

    if isinstance(experiment, commonTools.Experiment) is False:
        monitoring.to_log_and_console(proc + ": bad type for 'experiment' parameter", 1)
        sys.exit(1)
    if isinstance(parameters, IntraRegParameters) is False:
        monitoring.to_log_and_console(proc + ": bad type for 'parameters' parameter", 1)
        sys.exit(1)

    #
    # start processing
    #
    start_time = time.time()

    #
    # if template does not exists,
    # 1. compute transformations between successive images
    #    INTRAREG/INTRAREG_<EXP_INTRAREG>/CO-TRSFS directory
    # 2. compose transformations wrt a reference (default = first time point)
    #    INTRAREG/INTRAREG_<EXP_INTRAREG>/CO-TRSFS/TEMP directory
    # 3. re-compute transformations (ie translations) and template that includes all transformed images
    #    INTRAREG/INTRAREG_<EXP_INTRAREG>/TRSFS_t<first>-<last> directory
    #

    result_dir = _get_trsf_path(experiment)
    result_template = _get_template_path_name(experiment, parameters)

    if not os.path.isdir(result_dir) or os.path.isfile(result_template) \
            or parameters.rebuild_template is True or monitoring.forceResultsToBeBuilt is True:

        if experiment.delta_time_point > 1:
            monitoring.to_log_and_console(proc + ": warning, delta_time=" + str(experiment.delta_time_point)
                                          + ", this step may be fragile", 1)
        #
        # co-registration of any 2 successive images
        # will fill the INTRAREG/INTRAREG_<EXP_INTRAREG>/CO-TRSFS directory
        #

        monitoring.to_log_and_console("    .. co-registrations", 2)
        _coregistration_control(experiment, parameters.registration)

        #
        # composition of transformations by propagation
        # results in a set of transformation from the reference image (given by the reference_index)
        # towards each image
        # will fill the INTRAREG/INTRAREG_<EXP_INTRAREG>/TEMP directory
        #

        monitoring.to_log_and_console("    .. transformation composition", 2)

        temporary_dir = os.path.join(_get_cotrsf_path(experiment), 'TEMP')
        _transformations_from_reference(experiment, parameters, temporary_dir)

        #
        # re-composition of transformations
        # template creation for further resampling operations
        # will fill the INTRAREG/INTRAREG_<EXP_INTRAREG>/TRSFS_t<first>-<last> directory
        #

        monitoring.to_log_and_console("    .. transformation recomposition and template generation", 2)

        _transformations_and_template(experiment, parameters, temporary_dir)

        if monitoring.keepTemporaryFiles is False:
            shutil.rmtree(temporary_dir)

    #
    # template image and resampling transformations have been computed
    # resample images if required
    #
    # NOTE: il y a un probleme de coherence, puisque la premiere image de segmentation peut etre issue
    # de mars et donc etre nommee differemment. Pour le reechantillonage, on pourrait utiliser
    # reconstruction.get_segmentation_image().
    #

    if parameters.resample_fusion_images is True or parameters.movie_fusion_images is True \
            or len(parameters.xy_movie_fusion_images) > 0 or len(parameters.xz_movie_fusion_images) > 0 \
            or len(parameters.yz_movie_fusion_images) > 0:
        monitoring.to_log_and_console("    .. resampling fusion images", 2)
        _resample_images(experiment, parameters, result_template, 'fuse')

    if parameters.resample_segmentation_images is True or parameters.movie_segmentation_images is True \
            or len(parameters.xy_movie_segmentation_images) > 0 or len(parameters.xz_movie_segmentation_images) > 0 \
            or len(parameters.yz_movie_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. resampling segmentation images", 2)
        _resample_images(experiment, parameters, result_template, 'seg', interpolation_mode='nearest')

    if parameters.resample_post_segmentation_images is True or parameters.movie_post_segmentation_images is True \
            or len(parameters.xy_movie_post_segmentation_images) > 0 \
            or len(parameters.xz_movie_post_segmentation_images) > 0 \
            or len(parameters.yz_movie_post_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. resampling post-segmentation images", 2)
        _resample_images(experiment, parameters, result_template, 'post', interpolation_mode='nearest')

    #
    # make 3D=2D+t images = movie of evolving slices with respect to time
    #

    if parameters.movie_fusion_images is True or len(parameters.xy_movie_fusion_images) > 0 \
            or len(parameters.xz_movie_fusion_images) > 0 or len(parameters.yz_movie_fusion_images) > 0:
        monitoring.to_log_and_console("    .. movies from fusion images", 2)
        _make_movies(experiment, parameters, 'fuse', parameters.xy_movie_fusion_images,
                     parameters.xz_movie_fusion_images, parameters.yz_movie_fusion_images)

    if parameters.movie_segmentation_images is True or len(parameters.xy_movie_segmentation_images) > 0 \
            or len(parameters.xz_movie_segmentation_images) > 0 or len(parameters.yz_movie_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. movies from segmentation images", 2)
        _make_movies(experiment, parameters, 'seg', parameters.xy_movie_segmentation_images,
                     parameters.xz_movie_segmentation_images, parameters.yz_movie_segmentation_images)

    if parameters.movie_post_segmentation_images is True or len(parameters.xy_movie_post_segmentation_images) > 0 \
            or len(parameters.xz_movie_post_segmentation_images) > 0 \
            or len(parameters.yz_movie_post_segmentation_images) > 0:
        monitoring.to_log_and_console("    .. movies from post-segmentation images", 2)
        _make_movies(experiment, parameters, 'post', parameters.xy_movie_post_segmentation_images,
                     parameters.xz_movie_post_segmentation_images, parameters.yz_movie_post_segmentation_images)

    #
    # end processing
    #

    end_time = time.time()

    monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
    monitoring.to_log_and_console('', 1)

    return
