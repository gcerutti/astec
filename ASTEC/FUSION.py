
import os
import imp
import sys
import time
import math
import shutil
import subprocess
import numpy as np
from scipy import ndimage as nd

import commonTools
import nomenclature
from CommunFunctions.ImageHandling import SpatialImage, imread, imsave
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
# - computation environment
# - computation parameters
#
########################################################################################


class FusionEnvironment(object):

    def __init__(self):
        #
        # raw data
        #
        self.path_angle1 = None
        self.path_angle2 = None
        self.path_angle3 = None
        self.path_angle4 = None

        self.path_angle1_files = None
        self.path_angle2_files = None
        self.path_angle3_files = None
        self.path_angle4_files = None

        #
        # fused data
        #
        self.path_fuse = None
        self.path_fuse_exp = None
        self.path_fuse_exp_files = None

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

        self.path_angle1 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle1, parameters)
        self.path_angle2 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle2, parameters)
        self.path_angle3 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle3, parameters)
        self.path_angle4 = nomenclature.replaceFlags(nomenclature.path_rawdata_angle4, parameters)

        self.path_angle1_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle1_files, parameters)
        self.path_angle2_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle2_files, parameters)
        self.path_angle3_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle3_files, parameters)
        self.path_angle4_files = nomenclature.replaceFlags(nomenclature.path_rawdata_angle4_files, parameters)

        self.path_fuse = nomenclature.replaceFlags(nomenclature.path_fuse, parameters)
        self.path_fuse_exp = nomenclature.replaceFlags(nomenclature.path_fuse_exp, parameters)

        self.path_fuse_exp_files = nomenclature.replaceFlags(nomenclature.path_fuse_exp_files, parameters)

        self.path_logdir = nomenclature.replaceFlags(nomenclature.path_fuse_logdir, parameters)
        self.path_history_file = nomenclature.replaceFlags(nomenclature.path_fuse_historyfile, parameters)
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_fuse_logfile, parameters, start_time)

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('FusionEnvironment\n')
            logfile.write('- path_angle1 = ' + str(self.path_angle1)+'\n')
            logfile.write('- path_angle2 = ' + str(self.path_angle2)+'\n')
            logfile.write('- path_angle3 = ' + str(self.path_angle3)+'\n')
            logfile.write('- path_angle4 = ' + str(self.path_angle4)+'\n')

            logfile.write('- path_angle1_files = ' + str(self.path_angle1_files)+'\n')
            logfile.write('- path_angle2_files = ' + str(self.path_angle2_files)+'\n')
            logfile.write('- path_angle3_files = ' + str(self.path_angle3_files)+'\n')
            logfile.write('- path_angle4_files = ' + str(self.path_angle4_files)+'\n')

            logfile.write('- path_fuse = ' + str(self.path_fuse)+'\n')
            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp)+'\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files)+'\n')

            logfile.write('- path_logdir = ' + str(self.path_logdir) + '\n')
            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('FusionEnvironment')
        print('- path_angle1 = ' + str(self.path_angle1))
        print('- path_angle2 = ' + str(self.path_angle2))
        print('- path_angle3 = ' + str(self.path_angle3))
        print('- path_angle4 = ' + str(self.path_angle4))

        print('- path_angle1_files = ' + str(self.path_angle1_files))
        print('- path_angle2_files = ' + str(self.path_angle2_files))
        print('- path_angle3_files = ' + str(self.path_angle3_files))
        print('- path_angle4_files = ' + str(self.path_angle4_files))

        print('- path_fuse = ' + str(self.path_fuse))
        print('- path_fuse_exp = ' + str(self.path_fuse_exp))
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files))

        print('- path_logdir = ' + str(self.path_logdir))
        print('- path_history_file = ' + str(self.path_history_file))
        print('- path_log_file = ' + str(self.path_log_file))
        print("")


#
#
#
#
#


class FusionParameters(object):

    def __init__(self):
        #
        # acquisition parameters
        #
        self.acquisition_orientation = 'left'
        self.acquisition_mirrors = False
        self.acquisition_resolution = (0.17, 0.17, 1.0)
        self.acquisition_delay = 0

        #
        # Correction of slit lines
        #
        self.acquisition_slit_line_correction = False

        #
        # fused image parameters
        #
        self.target_resolution = (0.3, 0.3, 0.3)

        #
        # Cropping of acquisition images (before fusion)
        #
        self.acquisition_cropping = True
        self.acquisition_cropping_margin_x_0 = 40
        self.acquisition_cropping_margin_x_1 = 40
        self.acquisition_cropping_margin_y_0 = 40
        self.acquisition_cropping_margin_y_1 = 40

        #
        # Registration parameters
        #
        self.registration_transformation_type = 'affine'
        self.registration_transformation_estimation_type = 'wlts'
        self.registration_lts_fraction = 0.55
        self.registration_pyramid_highest_level = 6
        self.registration_pyramid_lowest_level = 3

        #
        # Cropping of fused image (after fusion)
        #
        self.fusion_cropping = True
        self.fusion_cropping_margin_x_0 = 40
        self.fusion_cropping_margin_x_1 = 40
        self.fusion_cropping_margin_y_0 = 40
        self.fusion_cropping_margin_y_1 = 40

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('FusionParameters\n')

            logfile.write('- acquisition_orientation = '+str(self.acquisition_orientation)+'\n')
            logfile.write('- acquisition_mirrors     = '+str(self.acquisition_mirrors)+'\n')
            logfile.write('- acquisition_resolution  = '+str(self.acquisition_resolution)+'\n')
            logfile.write('- acquisition_delay       = ' + str(self.acquisition_delay)+'\n')
            logfile.write('- target_resolution  = '+str(self.target_resolution)+'\n')

            logfile.write('- acquisition_cropping = '+str(self.acquisition_cropping)+'\n')
            logfile.write('- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0)+'\n')
            logfile.write('- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1)+'\n')
            logfile.write('- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0)+'\n')
            logfile.write('- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1)+'\n')

            logfile.write('- acquisition_slit_line_correction = '+str(self.acquisition_slit_line_correction)+'\n')

            logfile.write('- registration_transformation_type = ' + str(self.registration_transformation_type) + '\n')
            logfile.write('- registration_transformation_estimation_type = '
                          + str(self.registration_transformation_estimation_type) + '\n')
            logfile.write('- registration_lts_fraction = ' + str(self.registration_lts_fraction) + '\n')
            logfile.write('- registration_pyramid_highest_level = '
                          + str(self.registration_pyramid_highest_level) + '\n')
            logfile.write('- registration_pyramid_lowest_level = '
                          + str(self.registration_pyramid_lowest_level) + '\n')

            logfile.write('- fusion_cropping = '+str(self.fusion_cropping)+'\n')
            logfile.write('- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0)+'\n')
            logfile.write('- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1)+'\n')
            logfile.write('- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0)+'\n')
            logfile.write('- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1)+'\n')
            logfile.write("\n")
        return

    def print_parameters(self):
        print("")
        print('FusionParameters')

        print('- acquisition_orientation = '+str(self.acquisition_orientation))
        print('- acquisition_mirrors     = '+str(self.acquisition_mirrors))
        print('- acquisition_resolution  = '+str(self.acquisition_resolution))
        print('- acquisition_delay       = ' + str(self.acquisition_delay))
        print('- target_resolution  = '+str(self.target_resolution))

        print('- acquisition_cropping = '+str(self.acquisition_cropping))
        print('- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0))
        print('- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1))
        print('- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0))
        print('- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1))

        print('- acquisition_slit_line_correction = '+str(self.acquisition_slit_line_correction))

        print('- registration_transformation_type = ' + str(self.registration_transformation_type))
        print('- registration_transformation_estimation_type = '
              + str(self.registration_transformation_estimation_type))
        print('- registration_lts_fraction = ' + str(self.registration_lts_fraction))
        print('- registration_pyramid_highest_level = ' + str(self.registration_pyramid_highest_level))
        print('- registration_pyramid_lowest_level = ' + str(self.registration_pyramid_lowest_level))

        print('- fusion_cropping = '+str(self.fusion_cropping))
        print('- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0))
        print('- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1))
        print('- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0))
        print('- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1))
        print("")

    def update_from_file(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print ("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        #
        # acquisition parameters
        #
        if hasattr(parameters, 'raw_ori'):
            if parameters.raw_ori is not None:
                self.acquisition_orientation = parameters.raw_ori

        if hasattr(parameters, 'raw_mirrors'):
            if parameters.raw_mirrors is not None:
                self.acquisition_mirrors = parameters.raw_mirrors

        if hasattr(parameters, 'raw_resolution'):
            if parameters.raw_resolution is not None:
                self.acquisition_resolution = parameters.raw_resolution

        if hasattr(parameters, 'raw_delay'):
            if parameters.raw_delay is not None:
                self.acquisition_delay = parameters.raw_delay

        #
        # correction of slit lines
        #
        if hasattr(parameters, 'acquisition_slit_line_correction'):
            if parameters.acquisition_slit_line_correction is not None:
                self.acquisition_slit_line_correction = parameters.acquisition_slit_line_correction

        #
        # fused image parameters
        #
        if hasattr(parameters, 'target_resolution'):
            if parameters.target_resolution is not None:
                self.target_resolution = parameters.target_resolution

        #
        # Cropping of acquisition images (before fusion)
        #
        if hasattr(parameters, 'raw_crop'):
            if parameters.raw_crop is not None:
                self.acquisition_cropping = parameters.raw_crop
        if hasattr(parameters, 'raw_margin_x_0'):
            if parameters.raw_margin_x_0 is not None:
                self.acquisition_cropping_margin_x_0 = parameters.raw_margin_x_0
        if hasattr(parameters, 'raw_margin_x_1'):
            if parameters.raw_margin_x_1 is not None:
                self.acquisition_cropping_margin_x_1 = parameters.raw_margin_x_1
        if hasattr(parameters, 'raw_margin_y_0'):
            if parameters.raw_margin_y_0 is not None:
                self.acquisition_cropping_margin_y_0 = parameters.raw_margin_y_0
        if hasattr(parameters, 'raw_margin_y_1'):
            if parameters.raw_margin_y_1 is not None:
                self.acquisition_cropping_margin_y_1 = parameters.raw_margin_y_1

        #
        # Cropping of fused image (after fusion)
        #
        if hasattr(parameters, 'fusion_crop'):
            if parameters.fusion_crop is not None:
                self.fusion_cropping = parameters.fusion_crop
        if hasattr(parameters, 'fusion_margin_x_0'):
            if parameters.fusion_margin_x_0 is not None:
                self.fusion_cropping_margin_x_0 = parameters.fusion_margin_x_0
        if hasattr(parameters, 'fusion_margin_x_1'):
            if parameters.fusion_margin_x_1 is not None:
                self.fusion_cropping_margin_x_1 = parameters.fusion_margin_x_1
        if hasattr(parameters, 'fusion_margin_y_0'):
            if parameters.fusion_margin_y_0 is not None:
                self.fusion_cropping_margin_y_0 = parameters.fusion_margin_y_0
        if hasattr(parameters, 'fusion_margin_y_1'):
            if parameters.fusion_margin_y_1 is not None:
                self.fusion_cropping_margin_y_1 = parameters.fusion_margin_y_1




########################################################################################
#
# some internal procedures
#
########################################################################################


__recognized_extensions__ = ['.zip', '.h5', '.tif', '.tiff', '.TIF', '.TIFF', '.inr', '.inr.gz', '.mha', '.mha.gz']
__extension_to_be_converted__ = ['.h5', '.tif', '.tiff', '.TIF', '.TIFF']


def _get_extension(filename):
    """ Return the file extension. Must be in the set of recognized extensions.
    :param filename:
    :return: None in case of unrecognized extension,
             else the recognized extension (begins with '.')
    """
    for e in __recognized_extensions__:
        if len(filename) < len(e):
            continue
        if filename[len(filename)-len(e):len(filename)] == e:
            return e
    return None


def _add_suffix(filename, suffix, new_dirname=None, new_extension=None):
    """
    Add a suffix to a filenename (ie before the extension)
    :param filename:
    :param suffix: suffix to be added
    :param new_dirname: change the directory name of the file
    :param new_extension: change the extension of the file
    :return: the transformed file name
    """
    b = os.path.basename(filename)
    d = os.path.dirname(filename)
    e = _get_extension(b)
    if e is None:
        monitoring.to_log_and_console("_add_suffix: file extension of '"+str(filename)+"' was not recognized", 0)
        monitoring.to_log_and_console("\t Exiting", 0)
        sys.exit(1)
    new_basename = b[0:len(b)-len(e)]
    new_basename += suffix
    if new_extension is None:
        new_basename += e
    else:
        new_basename += new_extension
    if new_dirname is None:
        res_name = os.path.join(d, new_basename)
    else:
        res_name = os.path.join(new_dirname, new_basename)
    return res_name


def _read_image_name(data_path, temporary_path, file_name, resolution):
    """
    Read an image. Eventually, unzip a compressed file, and convert the image
    to a format known by executables
    :param data_path: path to data directory
    :param temporary_path: directory for temporary file
    :param file_name: image file
    :param resolution: resolution of the result image
            required to write the output image with the right resolution
    :return:
    """

    proc = "_read_image_name"
    default_extension = '.inr'

    #
    # test whether the extension is zip
    #
    f = file_name
    full_name = os.path.join(data_path, f)

    if f[len(f)-4:len(f)] == '.zip':

        prefix = f[0:len(f)-4]
        #
        # unzipping
        #
        monitoring.to_log_and_console("    .. unzipping '" + str(f) + "'", 2)
        cmd = 'unzip '+os.path.join(data_path, f) + ' -d ' + str(temporary_path)

        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        file_names = []
        for f in os.listdir(temporary_path):
            if len(f) <= len(prefix):
                pass
            if f[0:len(prefix)] == prefix:
                file_names.append(f)
        if len(file_names) == 0:
            monitoring.to_log_and_console(proc + ": no image with name '" + str(prefix)
                                          + "' was found in '" + str(temporary_path) + "'")
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)

        if len(file_names) > 1:
            monitoring.to_log_and_console(proc + ": several images with name '"
                                          + str(prefix) + "' were found in '" + str(temporary_path) + "'")
            monitoring.to_log_and_console("\t " + str(file_names))
            monitoring.to_log_and_console("\t Exiting")
            sys.exit(1)
        #
        #
        #
        f = file_names[0]
        full_name = os.path.join(temporary_path, f)

    #
    # test whether the file has to be converted into a more 'readable' format
    #
    for extension in __extension_to_be_converted__:
        if f[len(f)-len(extension):len(f)] == extension:
            prefix = f[0:len(f) - len(extension)]
            monitoring.to_log_and_console("    .. converting '" + str(f) + "'", 2)
            image = imread(full_name)
            image.resolution = resolution
            full_name = os.path.join(temporary_path, prefix) + default_extension
            imsave(full_name, image)
            break

    return full_name


def _find_image_name(data_path, image_name):
    """

    :param data_path:
    :param image_name:
    :return:
    """
    proc = "_find_image_name"

    file_names = []
    for f in os.listdir(data_path):
        if len(f) <= len(image_name):
            pass
        if f[0:len(image_name)] == image_name:
            file_names.append(f)

    if len(file_names) == 0:
        monitoring.to_log_and_console(proc + ": no image with name '" + str(image_name)
                                      + "' was found in '" + str(data_path) + "'", 4)
        return None

    if len(file_names) > 1:
        monitoring.to_log_and_console(proc + ": several images with name '"
                                      + str(image_name) + "' were found in '" + str(data_path) + "'")
        monitoring.to_log_and_console("\t "+str(file_names))
        return None

    return file_names[0]


def _analyze_data_directory(data_dir):
    """
    Parse a directory containing images
    :param data_dir:
    :return:
    1. the common prefix of image file names
    2. the number of characters used to encoded the variable part
       (time points)
    3. the list of the variable parts
    4. the common suffix of image file names
       may be longer than just the file extension
    """

    proc = "_analyze_data_directory"
    images = []
    extensions = []
    #
    # recognize images and extensions
    #
    for f in os.listdir(data_dir):
        for e in __recognized_extensions__:
            if f[len(f)-len(e):len(f)] == e:
                if e not in extensions:
                    extensions.append(e)
                    if len(extensions) > 1:
                        print proc + ": several image extensions were found in '" + data_dir + "'"
                        print "\t -> " + str(extensions)
                        print "\t Exiting."
                        sys.exit(1)
                images.append(f)

    #
    # one image case
    #

    if len(images) == 1:
        suffix = extensions[0]
        time_length = 0
        im = images[0]
        length = len(im) - 1 - len(suffix)
        for i in range(0, 3):
            if '0' <= im[length-i] <= '9':
                time_length += 1
            else:
                break
        prefix = im[0:len(im) - time_length - len(suffix)]
        time_points = im[len(im) - time_length - len(suffix):len(im) - len(suffix)]
        return prefix, time_length, time_points, suffix

    #
    # several images
    # 1. check that image names are of the same length
    # 2. get prefix = common part at beginning
    # 3. get suffix = common part at end
    # 4. get length for time point encoding
    # 5. get list of time points
    #

    for i in range(1, len(images)):
        if len(images[0]) != len(images[i]):
            print proc + ": image names are of different lengths in '" + data_dir + "'"
            print "\t -> " + images[0] + ", " + images[i]
            print "\t Exiting."
            sys.exit(1)

    prefix = ''
    for j in range(0, len(images[0])):
        ok = True
        for i in range(1, len(images)):
            if images[0][j] != images[i][j]:
                ok = False
                break
        if ok is True:
            prefix += images[0][j]
        else:
            break

    suffix = ''
    for j in range(len(images[0]) - 1, -1, -1):
        ok = True
        for i in range(1, len(images)):
            if images[0][j] != images[i][j]:
                ok = False
                break
        if ok is True:
            suffix += images[0][j]
        else:
            break
    suffix = suffix[::-1]

    time_length = len(images[0]) - len(prefix) - len(suffix)

    time_points = []
    for i in range(0, len(images)):
        time_points.append(images[i][len(images[i]) - time_length - len(suffix):len(images[i]) - len(suffix)])

    return prefix, time_length, time_points, suffix


########################################################################################
#
# cropping
#
########################################################################################


def _crop_spatial_image(image, margin_x_0=40, margin_x_1=40, margin_y_0=40, margin_y_1=40):
    """
    Crop a spatial image in XY plane
    :param image:
    :param margin_x_0:
    :param margin_x_1:
    :param margin_y_0:
    :param margin_y_1:
    :return:
    """

    #
    # MIP projection
    #
    mip_image = image.max(axis=2)

    #
    # get a threshold = image mean of the MIP image
    #
    threshold = np.mean(mip_image)

    #
    # create a binary image and apply the threshold
    # the get the connected component (4-connectivity)
    #
    bin_image = np.zeros((mip_image.shape[0], mip_image.shape[1], 1), dtype=np.uint8)
    bin_image[mip_image > threshold] = 1
    del mip_image
    cc_image, cc_n = nd.label(bin_image)

    #
    # compute the volumes of each connected component
    # and create a dictionary of tuples (label, volume)
    #
    labels = np.unique(cc_image)
    volumes = nd.sum(np.ones_like(cc_image), cc_image, index=np.int16(labels))
    dict_volumes = dict(zip(labels, volumes))

    #
    # remove the background
    # then get the label associated to the largest connected component
    dict_volumes.pop(0)
    max_label = dict_volumes.keys()[np.argmax(dict_volumes.values())]

    #
    # get the bounding boxes for all objects
    # it is not necessary to searched for all labels
    # seems that there is no bounding box computed for label #0
    #
    # boundingBoxes = nd.find_objects(ccImage, max_label=maxLabel)
    # maxBox = boundingBoxes[int(maxLabel)-1]
    #
    max_box = nd.find_objects(cc_image, max_label=max_label)[int(max_label)-1]
    del cc_image
    xmin = max(max_box[0].start - margin_x_0, 0)
    xmax = min(image.shape[0], max_box[0].stop + margin_x_1)
    ymin = max(max_box[1].start - margin_y_0, 0)
    ymax = min(image.shape[1], max_box[1].stop + margin_y_1)
    new_box = (slice(xmin, xmax, None),
               slice(ymin, ymax, None),
               slice(0, image.shape[2]))

    new_image = SpatialImage(image[new_box])
    new_image._set_resolution(image._get_resolution())

    monitoring.to_log_and_console("       crop from [0," + str(image.shape[0]) + "]x[0,"
                                  + str(image.shape[1]) + "] to [" + str(xmin) + ","
                                  + str(xmax) + "]x[" + str(ymin) + "," + str(ymax) + "]", 2)

    return new_image


def _crop_disk_image(the_image, res_image, margin_x_0=40, margin_x_1=40, margin_y_0=40, margin_y_1=40):
    """
    Crop an image on disk
    :param the_image:
    :param res_image:
    :param margin_x_0:
    :param margin_x_1:
    :param margin_y_0:
    :param margin_y_1:
    :return:
    """
    #
    # read input image
    #
    image = imread(the_image)

    new_image = _crop_spatial_image(image, margin_x_0, margin_x_1, margin_y_0, margin_y_1)
    del image

    # imsave(resImage, newImage.astype(np.uint16))
    imsave(res_image, new_image)
    del new_image
    return


########################################################################################
#
# computation of a rotation matrix
#
########################################################################################


def _axis_rotation_matrix(axis, angle, min_space=None, max_space=None):
    """ Return the transformation matrix from the axis and angle necessary
    axis : axis of rotation ("X", "Y" or "Z")
    angle : angle of rotation (in degree)
    min_space : coordinates of the bottom point (usually (0, 0, 0))
    max_space : coordinates of the top point (usually im shape)
    """
    i = np.linalg.inv
    d = np.dot
    if axis not in ["X", "Y", "Z"]:
        raise Exception("Unknown axis : " + str(axis))
    rads = math.radians(angle)
    s = math.sin(rads)
    c = math.cos(rads)

    centering = np.identity(4)
    if min_space is None and max_space is not None:
        min_space = np.array([0., 0., 0.])

    if max_space is not None:
        space_center = (max_space-min_space)/2.
        offset = -1.*space_center
        centering[:3, 3] = offset

    rot = np.identity(4)
    if axis == "X":
        rot = np.array([[1., 0., 0., 0.],
                        [0., c, -s, 0.],
                        [0., s, c, 0.],
                        [0., 0., 0., 1.]])
    elif axis == "Y":
        rot = np.array([[c,   0., s,  0.],
                        [0., 1., 0., 0.],
                        [-s, 0., c, 0.],
                        [0., 0., 0., 1.]])

    elif axis == "Z":
        rot = np.array([[c, -s,  0., 0.],
                        [s, c, 0., 0.],
                        [0., 0., 1., 0.],
                        [0., 0., 0., 1.]])

    return d(i(centering), d(rot, centering))


########################################################################################
#
# function for the ad hoc computation of weights
# for the linear combination of images of the 4 cameras
#
########################################################################################


def _histogram(image, nbins=256):
    """Return histogram of image.

        Unlike `np.histogram`, this function returns the centers of bins and
        does not rebin integer arrays. For integer arrays, each integer value has
        its own bin, which improves speed and intensity-resolution.

        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

        Returns
        -------
        hist : array
        The values of the histogram.
        bin_centers : array
        The values at the center of the bins.
        """

    proc = "_histogram"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    # For integer types, histogramming with bincount is more efficient.
    if np.issubdtype(image.dtype, np.integer):
        offset = 0
        if np.min(image) < 0:
            offset = np.min(image)
        hist = np.bincount(image.ravel() - offset)
        bin_centers = np.arange(len(hist)) + offset

        # clip histogram to start with a non-zero bin
        idx = np.nonzero(hist)[0][0]
        return hist[idx:], bin_centers[idx:]
    else:
        hist, bin_edges = np.histogram(image.flat, nbins)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.
        return hist, bin_centers


def _threshold_otsu(image, nbins=256):
    """Return threshold value based on Otsu's method.

        Parameters
        ----------
        image : array
        Input image.
        nbins : int
        Number of bins used to calculate histogram. This value is ignored for
        integer arrays.

        Returns
        -------
        threshold : float
        Threshold value.

        References
        ----------
        .. [1] Wikipedia, http://en.wikipedia.org/wiki/Otsu's_Method
        """

    proc = "_threshold_otsu"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    hist, bin_centers = _histogram(image, nbins)
    hist = hist.astype(float)

    # class probabilities for all possible thresholds
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]
    # class means for all possible thresholds
    mean1 = np.cumsum(hist * bin_centers) / weight1
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / weight2[::-1])[::-1]

    # Clip ends to align class 1 and class 2 variables:
    # The last value of `weight1`/`mean1` should pair with zero values in
    # `weight2`/`mean2`, which do not exist.
    variance12 = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance12)
    threshold = bin_centers[:-1][idx]
    return threshold


def _exp_func(x, length=500, speed=5):
    """ Decay function used to take into account the remotness to the camera
    x : value to compute
    length : lenght of the function
    speed : speed of the function
    """

    return .1+np.exp(-((np.float32(x)*speed)/length))


def _build_mask(image, direction):
    """Return the mask on a given image from the decay function
    im : intensity image (SpatialImage)
    direction : if True the camera is in the side of the first slices in Z
    """

    proc = "_build_mask"
    if not isinstance(image, SpatialImage):
        print proc + ": argument image is not an ndarray"
        return

    th = _threshold_otsu(image)
    im_th = np.zeros_like(image)
    im_th[image > th] = 1
    if direction is False:
        im_th = im_th[:, :, -1::-1]
    im_th_sum = np.cumsum(im_th, axis=2)
    if direction is False:
        im_th_sum = im_th_sum[:, :, -1::-1]
    mask = _exp_func(im_th_sum, np.max(im_th_sum))
    return mask


########################################################################################
#
#
#
########################################################################################


def fusion_process(input_images, fused_image, temporary_paths, parameters):
    """
    
    :param input_images:
    :param fused_image:
    :param temporary_paths:
    :param parameters:
    :return:
    """
    proc = 'fusion_process'

    if monitoring.debug > 1:
        print ""
        print proc + " was called with:"
        print "- input_images = " + str(input_images)
        print "- fused_image = " + str(fused_image)
        print "- temporary_paths = " + str(temporary_paths)
        print ""

    #
    # nothing to do if the fused image exists
    #
    if os.path.isfile(fused_image) and monitoring.forceResultsToBeBuilt is False:
        return

    #
    # how to copy a list:
    # NB: 'res_images = inputImages' acts like pointers
    #
    res_images = input_images[:]

    #
    # slit line correction
    # this correction has to be done on original data (without resampling)
    # Crop could be done beforehand to reduce the computational burden
    #
    if parameters.acquisition_slit_line_correction is True:
        the_images = res_images[:]
        res_images = []
        for i in range(0, len(the_images)):
            res_images.append(_add_suffix(input_images[i], "_line_corrected", new_dirname=temporary_paths[i]))

        for i in range(0, len(the_images)):
            monitoring.to_log_and_console("    .. correcting slit lines of '"
                                          + the_images[i].split(os.path.sep)[-1] + "'", 2)
            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                cpp_wrapping.slitline_correction(the_images[i], res_images[i])
            else:
                monitoring.to_log_and_console("       already existing", 2)

    #
    # to do: linear filtering to compensate for resolution change
    # for a change of voxel size from x0 to x1
    # smooth with a Gaussian of sigma = \sqrt(2)^(ln(x0/x1) / ln(2))
    #

    #
    # first change of resolution
    # - for X and Y: target resolution (supposed to be larger than original)
    # - for Z: original resolution (supposed to be larger than target)
    #

    the_images = res_images[:]
    res_images = []
    for i in range(0, len(the_images)):
        res_images.append(_add_suffix(input_images[i], "_resample", new_dirname=temporary_paths[i]))

    for i in range(0, len(the_images)):

        im = imread(the_images[i])
        if type(parameters.target_resolution) == int or type(parameters.target_resolution) == float:
            resampling_resolution = [parameters.target_resolution, parameters.target_resolution, im.voxelsize[2]]
        elif (type(parameters.target_resolution) == list or type(parameters.target_resolution) == tuple) \
                and len(parameters.target_resolution) == 3:
            resampling_resolution = [parameters.target_resolution[0], parameters.target_resolution[1], im.voxelsize[2]]
        else:
            monitoring.to_log_and_console(proc+': unable to set target resolution for first resampling', 0)
            monitoring.to_log_and_console("\t target resolution was '"+str(parameters.target_resolution)+"'", 0)
            monitoring.to_log_and_console("\t image resolution was '" + str(im.voxelsize) + "'", 0)
            monitoring.to_log_and_console("Exiting.", 0)
            sys.exit(1)
        del im

        monitoring.to_log_and_console("    .. resampling '" + the_images[i].split(os.path.sep)[-1]
                                      + "' at " + str(resampling_resolution), 2)
        if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
            cpp_wrapping.apply_transformation(the_images[i], res_images[i], the_transformation=None,
                                              template_image=None,
                                              voxel_size=resampling_resolution, nearest=False, monitoring=monitoring)
        else:
            monitoring.to_log_and_console("       already existing", 2)

    #
    # 2D crop of resampled acquisition images
    #

    if parameters.acquisition_cropping is True:
        the_images = res_images[:]
        res_images = []
        for i in range(0, len(the_images)):
            res_images.append(_add_suffix(input_images[i], "_crop", new_dirname=temporary_paths[i]))

        for i in range(0, len(the_images)):
            monitoring.to_log_and_console("    .. cropping '" + the_images[i].split(os.path.sep)[-1], 2)
            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                _crop_disk_image(the_images[i], res_images[i],
                                 parameters.acquisition_cropping_margin_x_0,
                                 parameters.acquisition_cropping_margin_x_1,
                                 parameters.acquisition_cropping_margin_y_0,
                                 parameters.acquisition_cropping_margin_y_1)
            else:
                monitoring.to_log_and_console("       already existing", 2)


    #
    # Mirroring of 'right' images if required
    #

    if parameters.acquisition_mirrors is False:
        the_images = res_images[:]
        res_images = []
        for i in range(0, len(the_images)):
            if i == 0 or i == 2:
                res_images.append(the_images[i])
            else:
                res_images.append(_add_suffix(input_images[i], "_mirror", new_dirname=temporary_paths[i]))

        for i in range(0, len(the_images)):
            if i == 0 or i == 2:
                continue
            monitoring.to_log_and_console("    .. mirroring  #" + str(i) + " '" + the_images[i].split(os.path.sep)[-1], 2)
            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                the_im = imread(the_images[i])
                res_im = SpatialImage(the_im.copy())[-1::-1, :, :]
                res_im._set_resolution(the_im._get_resolution())
                imsave(res_images[i], res_im)
                del the_im
                del res_im
            else:
                monitoring.to_log_and_console("       already existing", 2)

    #
    # 1. Putting all images in a common reference
    # - resampling of first image in an isotropic grid = reference image
    # - co-registration of other images
    # 2. Compute weights with an ad-hoc method
    #
    the_images = res_images[:]
    res_images = []
    init_trsfs = []
    res_trsfs = []
    unreg_weight_images = []
    weight_images = []

    for i in range(0, len(the_images)):
        res_images.append(_add_suffix(input_images[i], "_reg", new_dirname=temporary_paths[i]))
        init_trsfs.append(_add_suffix(input_images[i], "_init", new_dirname=temporary_paths[i], new_extension=".trsf"))
        res_trsfs.append(_add_suffix(input_images[i], "_reg", new_dirname=temporary_paths[i], new_extension=".trsf"))
        unreg_weight_images.append(_add_suffix(input_images[i], "_init_weight", new_dirname=temporary_paths[i]))
        weight_images.append(_add_suffix(input_images[i], "_weight", new_dirname=temporary_paths[i]))

    if parameters.acquisition_orientation == 'left':
        default_angle = 270.0
    else:
        default_angle = 90.0

    for i in range(0, len(the_images)):

        if i == 0:
            #
            # resampling first image
            #
            monitoring.to_log_and_console("    .. resampling '" + the_images[i].split(os.path.sep)[-1]
                                          + "' at " + str(parameters.target_resolution), 2)
            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                cpp_wrapping.apply_transformation(the_images[i], res_images[i], the_transformation=None,
                                                  template_image=None,
                                                  voxel_size=parameters.target_resolution, nearest=False,
                                                  monitoring=monitoring)
            else:
                monitoring.to_log_and_console("       already existing", 2)
        else:
            #
            # other images:
            # - set initial rotation
            # - register images
            #
            monitoring.to_log_and_console("    .. co-registering '" + the_images[i].split(os.path.sep)[-1], 2)

            if i == 1:
                angle = 0.0
            else:
                angle = default_angle
            monitoring.to_log_and_console("       angle used for '" + init_trsfs[i].split(os.path.sep)[-1]
                                          + "' is " + str(angle), 2)

            im = imread(the_images[i])
            rotation_matrix = _axis_rotation_matrix(axis="Y", angle=angle, min_space=(0, 0, 0),
                                                    max_space=np.multiply(im.shape[:3], im.resolution))
            del im

            np.savetxt(init_trsfs[i], rotation_matrix)
            del rotation_matrix

            if not os.path.isfile(res_images[i]) or monitoring.forceResultsToBeBuilt is True:
                #
                # a tow-fold registration, translation then affine, could be preferable
                #
                cpp_wrapping.linear_registration(res_images[0], the_images[i], res_images[i],
                                                 res_trsfs[i], init_trsfs[i],
                                                 py_hl=parameters.registration_pyramid_highest_level,
                                                 py_ll=parameters.registration_pyramid_lowest_level,
                                                 transformation_type=parameters.registration_transformation_type,
                                                 transformation_estimator=parameters.registration_transformation_estimation_type,
                                                 lts_fraction=parameters.registration_lts_fraction,
                                                 monitoring=monitoring)

            else:
                monitoring.to_log_and_console("       already existing", 2)

        #
        # compute weighting masks
        # - mask is computed on an untransformed image
        #   however, resolution may have changes, or it can be cropped
        #   or it can be mirrored
        # - mask are then transformed with the computed transformation
        #
        if i % 2 == 1:
            direction = False
        else:
            direction = True

        im = imread(the_images[i])
        unreg_weight = _build_mask(im, direction)
        unreg_weight._set_resolution(im._get_resolution())
        imsave(unreg_weight_images[i], unreg_weight)
        del im
        del unreg_weight

        monitoring.to_log_and_console("    .. resampling '" + unreg_weight_images[i].split(os.path.sep)[-1], 2)
        if i == 0:
            if not os.path.isfile(weight_images[i]) or monitoring.forceResultsToBeBuilt is True:
                cpp_wrapping.apply_transformation(unreg_weight_images[i], weight_images[i],
                                                  the_transformation=None, template_image=None,
                                                  voxel_size=parameters.target_resolution,
                                                  nearest=False, monitoring=monitoring)
            else:
                monitoring.to_log_and_console("       already existing", 2)
        else:
            if not os.path.isfile(weight_images[i]) or monitoring.forceResultsToBeBuilt is True:
                cpp_wrapping.apply_transformation(unreg_weight_images[i], weight_images[i],
                                                  the_transformation=res_trsfs[i], template_image=res_images[0],
                                                  voxel_size=None, nearest=False, monitoring=monitoring)
            else:
                monitoring.to_log_and_console("       already existing", 2)

        if i == 0:
            full_mask = imread(weight_images[i])
        else:
            full_mask += imread(weight_images[i])

    #
    # compute fused image as a linear combination of co-registered images
    # the sun of weights have been precomputed to mimic historical behavior
    #
    # do not forget to cast the result on 16 bits
    #
    if monitoring.debug > 0:
        tmp_mask_image = _add_suffix(fused_image, "_mask_sum", new_dirname=temporary_paths[4])
        imsave(tmp_mask_image, full_mask)

    monitoring.to_log_and_console("    .. combining images", 2)
    for i in range(0, len(the_images)):
        if i == 0:
            full_image = (imread(weight_images[i]) * imread(res_images[i])) / full_mask
        else:
            full_image += (imread(weight_images[i]) * imread(res_images[i])) / full_mask

    del full_mask

    full_image = full_image.astype(np.uint16)

    if monitoring.debug > 0:
        tmp_fused_image = _add_suffix(fused_image, "_uncropped_fusion", new_dirname=temporary_paths[4])
        imsave(tmp_fused_image, full_image)

    #
    # fused image can be cropped as well
    #
    if parameters.fusion_cropping is True:
        monitoring.to_log_and_console("    .. cropping '" + fused_image.split(os.path.sep)[-1], 2)
        full_image = _crop_spatial_image(full_image,
                                         parameters.fusion_cropping_margin_x_0,
                                         parameters.fusion_cropping_margin_x_1,
                                         parameters.fusion_cropping_margin_y_0,
                                         parameters.fusion_cropping_margin_y_1)

    imsave(fused_image, full_image)
    del full_image

    return


#
#
#
#
#


def fusion_preprocess(input_images, fused_image, time_point, environment, parameters):
    """

    :param input_images:
    :param fused_image:
    :param time_point:
    :param environment:
    :param parameters:
    :return:
    """

    proc = "fusion_preprocess"

    if monitoring.debug > 1:
        print ""
        print proc + " was called with:"
        print "- input_images = " + str(input_images)
        print "- fused_image = " + str(fused_image)
        print "- time_point = " + str(time_point)
        print ""

    monitoring.to_log_and_console('... fusion of time ' + time_point, 1)

    if os.path.isfile(fused_image):
        if not monitoring.forceResultsToBeBuilt:
            monitoring.to_log_and_console('    already existing', 2)
            return
        else:
            monitoring.to_log_and_console('    already existing, but forced', 2)

    #
    # start processing
    #

    start_time = time.time()

    #
    # directory for auxiliary files
    #
    # ANGLE_0: LC/Stack0000
    # ANGLE_1: RC/Stack0000
    # ANGLE_2: LC/Stack0001
    # ANGLE_3: LR/Stack0001
    #

    temporary_paths = list()

    temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_0"))
    temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_1"))
    temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_2"))
    temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_3"))
    temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME"))

    #
    # recall that time_point is a string here
    # nomenclature.replaceTIME() can not be used
    #
    for i in range(0, len(temporary_paths)):
        temporary_paths[i] = temporary_paths[i].replace(nomenclature.FLAG_TIME, time_point)
        if not os.path.isdir(temporary_paths[i]):
            os.makedirs(temporary_paths[i])

    #
    # get image file names
    # - may involve unzipping and conversion
    #
    monitoring.to_log_and_console('    get original images', 2)

    images = list()

    images.append(_read_image_name(environment.path_angle1, temporary_paths[0],
                                   input_images[0],
                                   parameters.acquisition_resolution))
    images.append(_read_image_name(environment.path_angle2, temporary_paths[1],
                                   input_images[1],
                                   parameters.acquisition_resolution))
    images.append(_read_image_name(environment.path_angle3, temporary_paths[2],
                                   input_images[2],
                                   parameters.acquisition_resolution))
    images.append(_read_image_name(environment.path_angle4, temporary_paths[3],
                                   input_images[3],
                                   parameters.acquisition_resolution))

    #
    #
    #
    monitoring.to_log_and_console('    fuse images', 2)

    fusion_process(images, fused_image, temporary_paths, parameters)

    #
    # remove temporary files if required
    #

    if monitoring.keepTemporaryFiles is False:
        shutil.rmtree(temporary_paths[4])

    #
    # end processing
    #

    end_time = time.time()
    monitoring.to_log_and_console('    computation time = ' + str(end_time - start_time) + ' s', 1)
    monitoring.to_log_and_console('', 1)

    return


#
#
#
#
#


def fusion_control(experiment, environment, parameters):
    """

    :param experiment:
    :param environment:
    :param parameters:
    :return:
    """

    proc = 'fusion_control'
    default_width = 3

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir(environment.path_fuse_exp):
        os.makedirs(environment.path_fuse_exp)

    if not os.path.isdir(environment.path_logdir):
        os.makedirs(environment.path_logdir)

    monitoring.to_log_and_console('', 1)

    #
    # if data directories are different, parse them
    #

    if environment.path_angle1 != environment.path_angle2 \
            and environment.path_angle1 != environment.path_angle3 \
            and environment.path_angle1 != environment.path_angle4 \
            and environment.path_angle2 != environment.path_angle3 \
            and environment.path_angle2 != environment.path_angle4 \
            and environment.path_angle3 != environment.path_angle4:

        prefix1, time_length1, time_points1, suffix1 = _analyze_data_directory(environment.path_angle1)
        prefix2, time_length2, time_points2, suffix2 = _analyze_data_directory(environment.path_angle2)
        prefix3, time_length3, time_points3, suffix3 = _analyze_data_directory(environment.path_angle3)
        prefix4, time_length4, time_points4, suffix4 = _analyze_data_directory(environment.path_angle4)

        if monitoring.debug > 0:
            print ""
            print "analysis of '" + str(environment.path_angle1) + "'"
            print "   -> " + prefix1
            print "   -> " + str(time_length1)
            print "   -> " + str(time_points1)
            print "   -> " + suffix1
            print "analysis of '" + str(environment.path_angle2) + "'"
            print "   -> " + prefix2
            print "   -> " + str(time_length2)
            print "   -> " + str(time_points2)
            print "   -> " + suffix2
            print "analysis of '" + str(environment.path_angle3) + "'"
            print "   -> " + prefix3
            print "   -> " + str(time_length3)
            print "   -> " + str(time_points3)
            print "   -> " + suffix3
            print "analysis of '" + str(environment.path_angle4) + "'"
            print "   -> " + prefix4
            print "   -> " + str(time_length4)
            print "   -> " + str(time_points4)
            print "   -> " + suffix4
            print ""

        #
        # loop over acquisitions
        # 1. case where all acquisition have to be processed
        #    begin < 0 or end < 0 or begin > end or delta < 0
        # 2. only a few acquisitions have to be processed
        #

        extra_zeros = ''
        if time_length1 < default_width:
            extra_zeros = (default_width - time_length1) * '0'

        if experiment.firstTimePoint < 0 or experiment.lastTimePoint < 0 or experiment.deltaTimePoint < 0 \
                or experiment.firstTimePoint > experiment.lastTimePoint:

            for time_point in time_points1:

                #
                # fused image name
                #

                fused_image = environment.path_fuse_exp_files.replace(nomenclature.FLAG_TIME, extra_zeros + time_point)

                #
                # input image names
                #

                images = list()

                images.append(prefix1 + time_point + suffix1)
                im = prefix2 + time_point + suffix2
                if time_point not in time_points2:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle2 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)
                im = prefix3 + time_point + suffix3
                if time_point not in time_points3:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle3 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)
                im = prefix4 + time_point + suffix4
                if time_point not in time_points4:
                    print proc + ": image '" + im + "' not found in '" + environment.path_angle4 + "'"
                    print "\t Exiting."
                    sys.exit(1)
                else:
                    images.append(im)

                #
                # process
                #

                fusion_preprocess(images, fused_image, extra_zeros + time_point, environment, parameters)

        else:

            for time_value in range(experiment.firstTimePoint, experiment.lastTimePoint + 1, experiment.deltaTimePoint):

                acquisition_time = str('{:0{width}d}'.format(time_value, width=time_length1))
                fused_time = str('{:0{width}d}'.format(time_value + parameters.acquisition_delay, width=time_length1))

                #
                # fused image name
                #

                fused_image = environment.path_fuse_exp_files.replace(nomenclature.FLAG_TIME, extra_zeros + fused_time)

                #
                # input image names
                #

                images = list()

                im = prefix1 + acquisition_time + suffix1
                if acquisition_time not in time_points1:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.path_angle1 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix2 + acquisition_time + suffix2
                if acquisition_time not in time_points2:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.path_angle2 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix3 + acquisition_time + suffix3
                if acquisition_time not in time_points3:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.path_angle3 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)
                im = prefix4 + acquisition_time + suffix4
                if acquisition_time not in time_points4:
                    monitoring.to_log_and_console("    .. image '" + im + "' not found in '"
                                                  + environment.path_angle4 + "'", 2)
                    monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                    continue
                else:
                    images.append(im)

                #
                # process
                #

                fusion_preprocess(images, fused_image, extra_zeros + acquisition_time, environment, parameters)

    #
    # here data directories are not different, we have to rely on built names
    #

    else:

        for time_value in range(experiment.firstTimePoint, experiment.lastTimePoint+1, experiment.deltaTimePoint):

            acquisition_time = str('{:0{width}d}'.format(time_value, width=default_width))

            #
            # fused image name
            #

            fused_image = nomenclature.replaceTIME(environment.path_fuse_exp_files,
                                                   time_value+parameters.acquisition_delay)

            #
            # input image names
            #

            images = list()

            name = nomenclature.replaceTIME(environment.path_angle1_files, time_value)
            im = _find_image_name(environment.path_angle1, name)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.path_angle1 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle2_files, time_value)
            im = _find_image_name(environment.path_angle2, name)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.path_angle2 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle3_files, time_value)
            im = _find_image_name(environment.path_angle3, name)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.path_angle3 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            name = nomenclature.replaceTIME(environment.path_angle4_files, time_value)
            im = _find_image_name(environment.path_angle4, name)
            if im is None:
                monitoring.to_log_and_console("    .. image '" + name + "' not found in '"
                                              + environment.path_angle4 + "'", 2)
                monitoring.to_log_and_console("       skip time " + str(acquisition_time), 2)
                continue
            else:
                images.append(im)

            #
            # process
            #

            fusion_preprocess(images, fused_image, acquisition_time, environment, parameters)

    return
