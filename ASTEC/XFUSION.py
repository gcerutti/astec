
import os, imp, sys
import time
import math
import copy
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


class FusionEnvironment( object ):

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

    def updateFromFile(self, parameterFile, starttime ):
        if ( parameterFile == None ):
            return
        if ( not os.path.isfile( parameterFile ) ):
            print ("Error: '" + parameterFile + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameterFile )

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
        self.path_log_file = nomenclature.replaceFlags(nomenclature.path_fuse_logfile, parameters, starttime )

    def writeParameters( self, logfileName ):
        with open(logfileName, 'a') as logfile:
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

    def printParameters( self ):
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

class FusionParameters( object ):

    def __init__(self):
        #
        # acquisition parameters
        #
        self.acquisition_orientation = 'left'
        self.acquisition_mirrors = False
        self.acquisition_resolution = ( 0.17, 0.17, 1.0 )
        self.acquisition_delay = 0

        #
        # fused image parameters
        #
        self.target_resolution = ( 0.3, 0.3, 0.3 )

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

    def writeParameters( self, logfileName ):
        with open(logfileName, 'a') as logfile:
            logfile.write("\n")
            logfile.write( 'FusionParameters\n')

            logfile.write( '- acquisition_orientation = '+str(self.acquisition_orientation)+'\n' )
            logfile.write( '- acquisition_mirrors     = '+str(self.acquisition_mirrors)+'\n' )
            logfile.write( '- acquisition_resolution  = '+str(self.acquisition_resolution)+'\n' )
            logfile.write( '- acquisition_delay       = ' + str(self.acquisition_delay)+'\n' )
            logfile.write( '- target_resolution  = '+str(self.target_resolution)+'\n' )

            logfile.write( '- acquisition_cropping = '+str(self.acquisition_cropping)+'\n' )
            logfile.write( '- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0)+'\n' )
            logfile.write( '- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1)+'\n' )
            logfile.write( '- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0)+'\n' )
            logfile.write( '- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1)+'\n' )

            logfile.write( '- registration_transformation_type = ' + str(self.registration_transformation_type) + '\n')
            logfile.write( '- registration_transformation_estimation_type = ' + str(self.registration_transformation_estimation_type) + '\n')
            logfile.write( '- registration_lts_fraction = ' + str(self.registration_lts_fraction) + '\n')
            logfile.write( '- registration_pyramid_highest_level = ' + str(self.registration_pyramid_highest_level) + '\n')
            logfile.write( '- registration_pyramid_lowest_level = ' + str(self.registration_pyramid_lowest_level) + '\n')

            logfile.write( '- fusion_cropping = '+str(self.fusion_cropping)+'\n' )
            logfile.write( '- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0)+'\n' )
            logfile.write( '- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1)+'\n' )
            logfile.write( '- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0)+'\n' )
            logfile.write( '- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1)+'\n' )
            logfile.write("\n")
        return

    def printParameters( self ):
        print("")
        print( 'FusionParameters')

        print( '- acquisition_orientation = '+str(self.acquisition_orientation) )
        print( '- acquisition_mirrors     = '+str(self.acquisition_mirrors) )
        print( '- acquisition_resolution  = '+str(self.acquisition_resolution) )
        print( '- acquisition_delay       = ' + str(self.acquisition_delay) )
        print( '- target_resolution  = '+str(self.target_resolution) )

        print( '- acquisition_cropping = '+str(self.acquisition_cropping) )
        print( '- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0) )
        print( '- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1) )
        print( '- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0) )
        print( '- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1) )

        print( '- registration_transformation_type = ' + str(self.registration_transformation_type))
        print( '- registration_transformation_estimation_type = ' + str(self.registration_transformation_estimation_type))
        print( '- registration_lts_fraction = ' + str(self.registration_lts_fraction))
        print( '- registration_pyramid_highest_level = ' + str(self.registration_pyramid_highest_level))
        print( '- registration_pyramid_lowest_level = ' + str(self.registration_pyramid_lowest_level))

        print( '- fusion_cropping = '+str(self.fusion_cropping) )
        print( '- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0) )
        print( '- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1) )
        print( '- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0) )
        print( '- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1) )
        print("")

    def updateFromFile( self, parameterFile ):
        if ( parameterFile == None ):
            return
        if ( not os.path.isfile( parameterFile ) ):
            print ("Error: '" + parameterFile + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameterFile )

        #
        # acquisition parameters
        #
        if hasattr(parameters, 'raw_ori'):
            if ( parameters.raw_ori != None ):
                self.acquisition_orientation = parameters.raw_ori

        if hasattr(parameters, 'raw_mirrors'):
            if ( parameters.raw_mirrors != None ):
                self.acquisition_mirrors = parameters.raw_mirrors

        if hasattr(parameters, 'raw_resolution'):
            if ( parameters.raw_resolution != None ):
                self.acquisition_resolution = parameters.raw_resolution

        if hasattr(parameters, 'raw_delay'):
            if ( parameters.raw_delay != None ):
                self.acquisition_delay = parameters.raw_delay

        #
        # fused image parameters
        #
        if hasattr(parameters, 'target_resolution'):
            if ( parameters.target_resolution != None ):
                self.target_resolution = parameters.target_resolution

        #
        # Cropping of acquisition images (before fusion)
        #
        if hasattr(parameters, 'raw_crop'):
            if ( parameters.raw_crop != None ):
                self.acquisition_cropping = parameters.raw_crop
        if hasattr(parameters, 'raw_margin_x_0'):
            if ( parameters.raw_margin_x_0 != None ):
                self.acquisition_cropping_margin_x_0 = parameters.raw_margin_x_0
        if hasattr(parameters, 'raw_margin_x_1'):
            if ( parameters.raw_margin_x_1 != None ):
                self.acquisition_cropping_margin_x_1 = parameters.raw_margin_x_1
        if hasattr(parameters, 'raw_margin_y_0'):
            if ( parameters.raw_margin_y_0 != None ):
                self.acquisition_cropping_margin_y_0 = parameters.raw_margin_y_0
        if hasattr(parameters, 'raw_margin_y_1'):
            if ( parameters.raw_margin_y_1 != None ):
                self.acquisition_cropping_margin_y_1 = parameters.raw_margin_y_1

        #
        # Cropping of fused image (after fusion)
        #
        if hasattr(parameters, 'fusion_crop'):
            if ( parameters.fusion_crop != None ):
                self.fusion_cropping = parameters.fusion_crop
        if hasattr(parameters, 'fusion_margin_x_0'):
            if ( parameters.fusion_margin_x_0 != None ):
                self.fusion_cropping_margin_x_0 = parameters.fusion_margin_x_0
        if hasattr(parameters, 'fusion_margin_x_1'):
            if ( parameters.fusion_margin_x_1 != None ):
                self.fusion_cropping_margin_x_1 = parameters.fusion_margin_x_1
        if hasattr(parameters, 'fusion_margin_y_0'):
            if ( parameters.fusion_margin_y_0 != None ):
                self.fusion_cropping_margin_y_0 = parameters.fusion_margin_y_0
        if hasattr(parameters, 'fusion_margin_y_1'):
            if ( parameters.fusion_margin_y_1 != None ):
                self.fusion_cropping_margin_y_1 = parameters.fusion_margin_y_1










########################################################################################
#
# some internal procedures
#
########################################################################################


recognizedExtensions = [ '.zip', '.h5','.tif','.tiff','.TIF','.TIFF', '.inr', '.inr.gz', '.mha', '.mha.gz' ]

def _getExtension( filename ):
    """ Return the file extension. Must be in the set of recognized extensions.
    :param filename:
    :return: None in case of unrecognized extension,
             else the recognized extension (begins with '.')
    """
    for e in recognizedExtensions:
        if len( filename ) < len( e ):
            continue
        if filename[len(filename)-len(e):len(filename)] == e:
            return e
    return None





def _addSuffix( filename, suffix, newDirname=None, newExtension=None ):
    """ Add a suffix to a filenename (ie before the extension)
    :param filename:
    :param suffix: suffix to be added
    :param newDirname: change the dirname of the fil
    :param newExtension: change the extension of the file
    :return: the transformed filename
    """
    b=os.path.basename(filename)
    d=os.path.dirname(filename)
    e = _getExtension( b )
    if e == None:
        monitoring.toLogAndConsole( "_addSuffix: file extension of '"+str(filename)+"' was not recognized", 0 )
        monitoring.toLogAndConsole( "\t Exiting", 0 )
        sys.exit( 1 )
    newBasename=b[0:len(b)-len(e)]
    newBasename += suffix
    if newExtension==None:
        newBasename += e
    else:
        newBasename += newExtension
    if newDirname == None:
        resName = os.path.join(d, newBasename)
    else:
        resName = os.path.join(newDirname, newBasename)
    return resName





extensionToBeConverted=['.h5','.tif','.tiff','.TIF','.TIFF']

def _readImageName( datapath, temporary_path, prefix, resolution ):

    proc=_readImageName
    defaultExtension = '.inr'

    fileNames=[]
    for f in os.listdir(datapath):
        if len( f ) <= len( prefix ):
            pass
        if f[0:len(prefix)] == prefix:
            fileNames.append(f)

    if len(fileNames)==0:
        monitoring.toLogAndConsole( proc+": no image with name '"+str(prefix)+"' was found in '"+str(datapath)+"'", 0)
        monitoring.toLogAndConsole( "\t Exiting", 0)
        sys.exit(1)

    if len(fileNames) > 1:
        monitoring.toLogAndConsole(proc + ": several images with name '" + str(prefix) + "' were found in '" + str(datapath) + "'")
        monitoring.toLogAndConsole( "\t "+str(fileNames) )
        monitoring.toLogAndConsole("\t Exiting")
        sys.exit(1)

    #
    # test whether the extension is zip
    #
    f=fileNames[0]
    extension = f[len(prefix):len(f)]
    fullName = os.path.join( datapath, f )

    if extension == '.zip':

        #
        # unzipping
        #
        monitoring.toLogAndConsole( "    .. unzipping '"+str(f)+ "'", 2 )
        cmd='unzip '+os.path.join( datapath, f )+' -d '+str(temporary_path)

        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        fileNames = []
        for f in os.listdir(temporary_path):
            if len(f) <= len(prefix):
                pass
            if f[0:len(prefix)] == prefix:
                fileNames.append(f)
        if len(fileNames) == 0:
            monitoring.toLogAndConsole(proc + ": no image with name '" + str(prefix) + "' was found in '" + str(temporary_path) + "'")
            monitoring.toLogAndConsole("\t Exiting")
            sys.exit(1)

        if len(fileNames) > 1:
            monitoring.toLogAndConsole(proc + ": several images with name '" + str(prefix) + "' were found in '" + str(temporary_path) + "'")
            monitoring.toLogAndConsole("\t " + str(fileNames))
            monitoring.toLogAndConsole("\t Exiting")
            sys.exit(1)
        #
        #
        #
        f = fileNames[0]
        fullName = os.path.join(temporary_path, f)

    #
    #
    #
    extension = f[len(prefix):len(f)]

    #
    # test whether the file has to be converted into a more 'readable' format
    #
    if extension in extensionToBeConverted:
        monitoring.toLogAndConsole( "    .. converting '"+str(f)+ "'", 2 )
        image = imread( os.path.join( temporary_path, f) )
        image.resolution = resolution
        fullName = os.path.join( temporary_path, prefix)+defaultExtension
        imsave(fullName, image)

    return fullName




def _cropNdImage( image, margin_x_0=40, margin_x_1=40, margin_y_0=40, margin_y_1=40 ):

    #
    # MIP projection
    #
    mipImage = image.max(axis=2)

    #
    # get a threshold = image mean of the MIP image
    #
    threshold = np.mean(mipImage)

    #
    # create a binary image and apply the threshold
    # the get the connected component (4-connectivity)
    #
    binImage = np.zeros((mipImage.shape[0], mipImage.shape[1], 1), dtype=np.uint8)
    binImage[mipImage > threshold] = 1
    ccImage, nCC = nd.label(binImage)

    #
    # compute the volumes of each connected component
    # and create a dictionary of tuples (label, volume)
    #
    labels = np.unique( ccImage )
    volumes = nd.sum(np.ones_like(ccImage), ccImage, index=np.int16(labels))
    dictVolumes = dict(zip(labels, volumes))

    #
    # remove the background
    # then get the label associated to the largest connected component
    dictVolumes.pop(0)
    maxLabel = dictVolumes.keys()[np.argmax(dictVolumes.values())]

    #
    # get the bounding boxes for all objects
    # it is not necessary to searched for all labels
    # seems that there is no bounding box computed for label #0
    #
    # boundingBoxes = nd.find_objects(ccImage, max_label=maxLabel)
    # maxBox = boundingBoxes[int(maxLabel)-1]
    #
    maxBox = nd.find_objects(ccImage, max_label=maxLabel)[int(maxLabel)-1]
    xmin = max(maxBox[0].start - margin_x_0, 0)
    xmax = min(image.shape[0], maxBox[0].stop + margin_x_1)
    ymin = max(maxBox[1].start - margin_y_0, 0)
    ymax = min(image.shape[1], maxBox[1].stop + margin_y_1)
    newBox = (slice(xmin, xmax, None),
               slice(ymin, ymax, None),
               slice(0, image.shape[2]))

    newImage = SpatialImage(image[newBox])
    newImage._set_resolution( image._get_resolution() )

    monitoring.toLogAndConsole("       crop from [0,"+str(image.shape[0])+"]x[0,"+str(image.shape[1])+"] to ["+
                               str(xmin)+","+str(xmax)+"]x["+str(ymin)+","+str(ymax)+"]", 2)


    return newImage



def _cropDiskImage( theImage, resImage, margin_x_0=40, margin_x_1=40, margin_y_0=40, margin_y_1=40 ):

    #
    # read input image
    #
    image = imread( theImage )

    newImage = _cropNdImage( image, margin_x_0, margin_x_1, margin_y_0, margin_y_1 )

    # imsave(resImage, newImage.astype(np.uint16))
    imsave(resImage, newImage )





def _axis_rotation_matrix(axis, angle, min_space=None, max_space=None):
    """ Return the transformation matrix from the axis and angle necessary
    axis : axis of rotation ("X", "Y" or "Z")
    angle : angle of rotation (in degree)
    min_space : coordinates of the bottom point (usually (0, 0, 0))
    max_space : coordinates of the top point (usually im shape)
    """
    I = np.linalg.inv
    D = np.dot
    if axis not in ["X", "Y", "Z"]:
        raise Exception("Unknown axis : "+ str(axis))
    rads = math.radians(angle)
    s = math.sin(rads)
    c = math.cos(rads)

    centering = np.identity(4)
    if min_space is None and max_space is not None:
        min_space = np.array([0.,0.,0.])

    if max_space is not None:
        space_center = (max_space-min_space)/2.
        offset = -1.*space_center
        centering[:3,3] = offset

    rot = np.identity(4)
    if axis=="X":
        rot = np.array([ [1., 0., 0., 0.],
                            [0., c, -s,  0.],
                            [0., s,  c,  0.],
                            [0., 0., 0., 1.] ])
    elif axis=="Y":
        rot = np.array([ [c,   0., s,  0.],
                            [0.,  1., 0., 0.],
                            [-s,  0., c,  0.],
                            [0.,  0., 0., 1.] ])

    elif axis=="Z":
        rot = np.array([ [c, -s,  0., 0.],
                            [s,  c,  0., 0.],
                            [0., 0., 1., 0.],
                            [0., 0., 0., 1.] ])

    return D(I(centering), D(rot, centering))





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

        Examples
        --------
        >>> from skimage.data import camera
        >>> image = camera()
        >>> thresh = threshold_otsu(image)
        >>> binary = image > thresh
        """
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

    return .1+np.exp(-((np.float32(x)*speed/length)))




def _build_mask(im, direction):
    """Return the mask on a given image from the decay function
    im : intensity image (SpatialImage)
    direction : if True the camera is in the side of the first slices in Z
    """
    th = _threshold_otsu(im)
    im_th=np.zeros_like(im)
    im_th[im>th]=1
    if direction==False:
        im_th=im_th[:,:,-1::-1]
    im_th_sum=np.cumsum(im_th, axis=2)
    if direction==False:
        im_th_sum=im_th_sum[:,:,-1::-1]
    mask = _exp_func(im_th_sum, np.max(im_th_sum))
    return mask










########################################################################################
#
#
#
########################################################################################


def fusionImages( inputImages, fusedImage, temporary_paths, environment, parameters ):

    proc = 'fusionImages';

    #
    # nothing to do if the fused image exists
    #
    if os.path.isfile(fusedImage) and monitoring.forceResultsToBeBuilt == False:
        return

    #
    # copie de liste
    # NB: 'theImages = inputImages' acts like pointers
    #
    theImages = inputImages[:]
    resImages = []

    #
    # to do: correct for planes
    #

    #
    # to do: linear filtering to compensate for resolution change
    # for a change of voxel size from x0 to x1
    # smooth with a Gaussian of sigma = \sqrt(2)^( ln(x0/x1) / ln(2) )
    #

    #
    # first change of resolution
    # - for X and Y: target resolution (supposed to be larger than original)
    # - for Z: original resolution (supposed to be larger than target)
    #

    for i in range(0, len(theImages)):
        resImages.append( _addSuffix( inputImages[i], "_resample", newDirname=temporary_paths[i] ) )

    for i in range(0, len(theImages)):

        im = imread( theImages[i] )
        if type(parameters.target_resolution)==int or type(parameters.target_resolution)==float:
            resampling_resolution = [parameters.target_resolution, parameters.target_resolution, im.voxelsize[2]]
        elif (type(parameters.target_resolution)==list or type(parameters.target_resolution)==tuple) and len(parameters.target_resolution) == 3:
            resampling_resolution = [parameters.target_resolution[0], parameters.target_resolution[1], im.voxelsize[2]]
        else:
            monitoring.toLogAndConsole( proc+': unable to set target resolution for first resampling', 0)
            monitoring.toLogAndConsole( "\t target resolution was '"+str(parameters.target_resolution)+"'", 0)
            monitoring.toLogAndConsole( "\t image resolution was '" + str(im.voxelsize) + "'", 0)
            monitoring.toLogAndConsole( "Exiting.", 0)
            sys.exit( 1 )
        del im

        monitoring.toLogAndConsole( "    .. resampling '"+theImages[i].split(os.path.sep)[-1]+ "' at "+str(resampling_resolution), 2)
        if not os.path.isfile( resImages[i] ) or monitoring.forceResultsToBeBuilt == True:
            cpp_wrapping.applyTrsfCLI(theImages[i], resImages[i], theTrsf=None, templateImage=None,
                                      voxelsize=resampling_resolution, nearest=False, monitoring=monitoring )
        else:
            monitoring.toLogAndConsole("       already existing", 2 )

    #
    # 2D crop of resampled acquisition images
    #


    if parameters.acquisition_cropping == True:
        theImages = resImages[:]
        resImages = []
        for i in range(0, len(theImages)):
            resImages.append(_addSuffix(inputImages[i], "_crop", newDirname=temporary_paths[i]))

        for i in range(0, len(theImages)):
            monitoring.toLogAndConsole("    .. cropping '" + theImages[i].split(os.path.sep)[-1], 2)
            if not os.path.isfile(resImages[i]) or monitoring.forceResultsToBeBuilt == True:
                _cropDiskImage( theImages[i], resImages[i],
                                parameters.acquisition_cropping_margin_x_0,
                                parameters.acquisition_cropping_margin_x_1,
                                parameters.acquisition_cropping_margin_y_0,
                                parameters.acquisition_cropping_margin_y_1)
            else:
                monitoring.toLogAndConsole("       already existing", 2)

    #
    # Mirroring of 'right' images if required
    #

    if parameters.acquisition_mirrors == False:
        theImages = resImages[:]
        resImages = []
        for i in range(0, len(theImages)):
            if i == 0 or i == 2:
                resImages.append( theImages[i] )
            else:
                resImages.append(_addSuffix(inputImages[i], "_mirror", newDirname=temporary_paths[i]))

        for i in range(0, len(theImages)):
            if i == 0 or i == 2:
                continue
            monitoring.toLogAndConsole("    .. mirroring '" + theImages[i].split(os.path.sep)[-1], 2)
            if not os.path.isfile(resImages[i]) or monitoring.forceResultsToBeBuilt == True:
                theIm = imread( theImages[i] )
                resIm = SpatialImage(theIm.copy())[-1::-1, :, :]
                resIm._set_resolution(theIm._get_resolution())
                imsave( resImages[i], resIm )
            else:
                monitoring.toLogAndConsole("       already existing", 2)

    #
    # 1. Putting all images in a common reference
    # - resampling of first image in an isotropic grid = reference image
    # - co-registration of other images
    # 2. Compute weights with an ad-hoc method
    #
    theImages = resImages[:]
    resImages = []
    initTrsfs = []
    resTrsfs = []
    unregWeightImages = []
    weightImages = []

    for i in range(0, len(theImages)):
        resImages.append(_addSuffix(inputImages[i], "_reg", newDirname=temporary_paths[i]))
        initTrsfs.append(_addSuffix(inputImages[i], "_init", newDirname=temporary_paths[i], newExtension=".trsf"))
        resTrsfs.append(_addSuffix(inputImages[i], "_reg", newDirname=temporary_paths[i], newExtension=".trsf"))
        unregWeightImages.append(_addSuffix(inputImages[i], "_initweight", newDirname=temporary_paths[i]))
        weightImages.append(_addSuffix(inputImages[i], "_weight", newDirname=temporary_paths[i]))

    if parameters.acquisition_orientation == 'left':
        defaultAngle = 270.0
    else:
        defaultAngle = 90.0

    for i in range(0, len(theImages)):

        if i==0:
            #
            # resampling first image
            #
            monitoring.toLogAndConsole( "    .. resampling '" + theImages[i].split(os.path.sep)[-1] + "' at " + str(parameters.target_resolution), 2)
            if not os.path.isfile(resImages[i]) or monitoring.forceResultsToBeBuilt == True:
                cpp_wrapping.applyTrsfCLI(theImages[i], resImages[i], theTrsf=None, templateImage=None,
                                          voxelsize=parameters.target_resolution, nearest=False, monitoring=monitoring)
            else:
                monitoring.toLogAndConsole("       already existing", 2)
        else:
            #
            # other images:
            # - set initial rotation
            # - register images
            #
            monitoring.toLogAndConsole("    .. co-registering '" + theImages[i].split(os.path.sep)[-1], 2 )

            if i==1:
                angle = 0.0
            else:
                angle =defaultAngle
            monitoring.toLogAndConsole("       angle used for '" + initTrsfs[i].split(os.path.sep)[-1]+"' is "+str(angle), 2)

            im = imread(theImages[i])
            rotationMatrix = _axis_rotation_matrix(axis="Y", angle=angle, min_space=(0, 0, 0),
                                       max_space=np.multiply(im.shape[:3], im.resolution))
            del im

            np.savetxt(initTrsfs[i], rotationMatrix)

            if not os.path.isfile(resImages[i]) or monitoring.forceResultsToBeBuilt == True:
                #
                # a tow-fold registration, translation then affine, could be preferable
                #
                cpp_wrapping.singleRegistrationCLI( resImages[0], theImages[i], resImages[i], resTrsfs[i], initTrsfs[i],
                                                    py_hl=parameters.registration_pyramid_highest_level,
                                                    py_ll=parameters.registration_pyramid_lowest_level,
                                                    trsf_type=parameters.registration_transformation_type,
                                                    trsf_estimator=parameters.registration_transformation_estimation_type,
                                                    lts_fraction=parameters.registration_lts_fraction,
                                                    monitoring=monitoring)
            else:
                monitoring.toLogAndConsole("       already existing", 2)

        #
        # compute weighting masks
        # - mask is computed on an untransformed image
        #   however, resolution may have changes, or it can be cropped
        #   or it can be mirrored
        # - mask are then transformed with the computed transformation
        #
        if i%2==1:
            direction=False
        else:
            direction=True

        im = imread(theImages[i])
        unRegMask = _build_mask(im, direction)
        unRegMask._set_resolution(im._get_resolution())
        imsave( unregWeightImages[i], unRegMask )
        del im

        monitoring.toLogAndConsole("    .. resampling '" + unregWeightImages[i].split(os.path.sep)[-1], 2)
        if i==0:
            if not os.path.isfile(weightImages[i]) or monitoring.forceResultsToBeBuilt == True:
                cpp_wrapping.applyTrsfCLI(unregWeightImages[i], weightImages[i], theTrsf=None, templateImage=None,
                                          voxelsize=parameters.target_resolution, nearest=False, monitoring=monitoring)
            else:
                monitoring.toLogAndConsole("       already existing", 2)
        else:
            if not os.path.isfile(weightImages[i]) or monitoring.forceResultsToBeBuilt == True:
                cpp_wrapping.applyTrsfCLI(unregWeightImages[i], weightImages[i], theTrsf=resTrsfs[i], templateImage=resImages[0],
                                          voxelsize=None, nearest=False, monitoring=monitoring)
            else:
                monitoring.toLogAndConsole("       already existing", 2)

        if i == 0:
            fullMask = imread( weightImages[i] )
        else:
            fullMask += imread( weightImages[i] )

    #
    # compute fused image as a linear combination of co-registered images
    # the sun of weights have been precomputed to mimic historical behavior
    #
    # do not forget to cast the result on 16 bits
    #
    if monitoring.debug > 0:
        tmpMaskImage = _addSuffix(fusedImage, "_mask_sum", newDirname=temporary_paths[4])
        imsave(tmpMaskImage, fullMask)

    monitoring.toLogAndConsole("    .. combining images", 2)
    for i in range(0, len(theImages)):
        if i==0:
            fullImage = ( imread(weightImages[i]) * imread(resImages[i]) ) / fullMask
        else:
            fullImage += ( imread(weightImages[i]) * imread(resImages[i]) ) / fullMask

    del fullMask

    fullImage = fullImage.astype(np.uint16)

    if monitoring.debug > 0:
        tmpFusedImage = _addSuffix(fusedImage, "_uncropped_fusion", newDirname=temporary_paths[4])
        imsave(tmpFusedImage, fullImage)

    #
    # fused image can be cropped as well
    #
    if parameters.fusion_cropping == True:
        monitoring.toLogAndConsole("    .. cropping '" + fusedImage.split(os.path.sep)[-1], 2)
        fullImage= _cropNdImage( fullImage,
                                  parameters.fusion_cropping_margin_x_0,
                                  parameters.fusion_cropping_margin_x_1,
                                  parameters.fusion_cropping_margin_y_0,
                                  parameters.fusion_cropping_margin_y_1)

    imsave(fusedImage, fullImage)
    del fullImage

    return









def fusionProcess( experiment, environment, parameters ):

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir( environment.path_fuse_exp ):
        os.makedirs( environment.path_fuse_exp )

    if not os.path.isdir( environment.path_logdir ):
        os.makedirs( environment.path_logdir )

    monitoring.toLogAndConsole('', 1)

    #
    # loop over acquisitions
    #

    for timePoint in range( experiment.firstTimePoint, experiment.lastTimePoint+1, experiment.deltaTimePoint):

        #
        # result image
        #

        fusedImage = nomenclature.replaceTIME(environment.path_fuse_exp_files, timePoint)

        monitoring.toLogAndConsole('... fusion of time ' + str(timePoint), 1)

        if os.path.isfile(fusedImage):
            if not monitoring.forceResultsToBeBuilt:
                monitoring.toLogAndConsole('    already existing', 2 )
                continue
            else:
                monitoring.toLogAndConsole('    already existing, but forced', 2)

        #
        # start processing
        #

        starttime = time.time()


        #
        # directory for auxiliary files
        #
        temporary_paths = []


        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_0" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_1" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_2" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_3" ) )
        temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME"))

        for i in range( 0, len(temporary_paths) ):
            temporary_paths[i] = nomenclature.replaceTIME( temporary_paths[i], timePoint )
            if not os.path.isdir(temporary_paths[i]):
                os.makedirs(temporary_paths[i])


        #
        # get image file names
        # - may involve unzipping and conversion
        #
        monitoring.toLogAndConsole('    get original images', 2 )

        images=[]

        images.append( _readImageName( environment.path_angle1, temporary_paths[0],
                        nomenclature.replaceTIME(environment.path_angle1_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle2, temporary_paths[1],
                        nomenclature.replaceTIME(environment.path_angle2_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle3, temporary_paths[2],
                        nomenclature.replaceTIME(environment.path_angle3_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle4, temporary_paths[3],
                        nomenclature.replaceTIME(environment.path_angle4_files,timePoint), parameters.acquisition_resolution ) )


        #
        #
        #
        monitoring.toLogAndConsole('    fuse images', 2 )

        fusionImages( images, fusedImage, temporary_paths, environment, parameters )










        if monitoring.keepTemporaryFiles == False:
            os.rmdir(temporary_path)
        #
        # end processing
        #

        endtime = time.time()
        monitoring.toLogAndConsole('    computation time = '+str(endtime-starttime)+ ' s', 1)
        monitoring.toLogAndConsole('', 1)


    return











