
import os, imp, sys
import time
import subprocess

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





#
#
#
#
#

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
#
#
########################################################################################


recognizedExtensions = [ '.zip', '.h5','.tif','.tiff','.TIF','.TIFF', '.inr', '.inr.gz', '.mha', '.mha.gz' ]

def _getExtension( filename ):
    for e in recognizedExtensions:
        if len( filename ) < len( e ):
            continue
        if filename[len(filename)-len(e):len(filename)] == e:
            return e
    return None





def _addSuffix( filename, suffix, extension=None ):
    e = _getExtension( filename )
    if e == None:
        print( "_addSuffix: file extension of '"+str(filename)+"' was not recognized")
        print( "Exiting" )
        sys.exit( 1 )
    newname=filename[0:len(filename)-len(e)]
    newname += suffix
    if extension==None:
        newname += e
    else:
        newname += extension
    return( newname )





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
        print( proc+": no image with name '"+str(prefix)+"' was found in '"+str(datapath)+"'")
        print( "\t Exiting")
        sys.exit(1)

    if len(fileNames) > 1:
        print(proc + ": several images with name '" + str(prefix) + "' were found in '" + str(datapath) + "'")
        print( "\t "+str(fileNames) )
        print("\t Exiting")
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
        if monitoring.verbose >= 2:
            print( "    .. unzipping '"+str(f)+ "'")
        cmd='unzip '+os.path.join( datapath, f )+' -d '+str(temporary_path)

        subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)

        fileNames = []
        for f in os.listdir(temporary_path):
            if len(f) <= len(prefix):
                pass
            if f[0:len(prefix)] == prefix:
                fileNames.append(f)
        if len(fileNames) == 0:
            print(proc + ": no image with name '" + str(prefix) + "' was found in '" + str(temporary_path) + "'")
            print("\t Exiting")
            sys.exit(1)

        if len(fileNames) > 1:
            print(proc + ": several images with name '" + str(prefix) + "' were found in '" + str(temporary_path) + "'")
            print("\t " + str(fileNames))
            print("\t Exiting")
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
        if monitoring.verbose >= 2:
            print( "    .. converting '"+str(f)+ "'")
        image = imread( os.path.join( temporary_path, f) )
        image.resolution = resolution
        fullName = os.path.join( temporary_path, prefix)+defaultExtension
        imsave(fullName, image)

    return fullName





def fusionImages( inputImages, fusedImage, temporary_paths, environment, parameters ):

    proc = 'fusionImages';
    print( 'fusionImages was called with')
    print( '- inputImages='+str(inputImages) )
    print( '- fusedImage=' + str(fusedImage))
    print( '- temporary_paths=' + str(temporary_paths))

    theImages = inputImages
    print( 'theImages ='+str(theImages))
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

    resImages = []
    for i in range(0, len(inputImages)):
        print( "_addSuffix to '"+inputImages[i]+"'")
        resImages.append( _addSuffix( inputImages[i], "_resample1" ) )

    for i in range(0, len(inputImages)):
        im = imread( theImages[i] )

        if type(parameters.target_resolution)==int or type(parameters.target_resolution)==float:
            resampling_resolution = [parameters.target_resolution, parameters.target_resolution, im.voxelsize[2]]
        elif (type(parameters.target_resolution)==list or type(parameters.target_resolution)==tuple) and len(parameters.target_resolution) == 3:
            resampling_resolution = [parameters.target_resolution[0], parameters.target_resolution[1], im.voxelsize[2]]
        else:
            print( proc+': unable to set target resolution for first resampling')
            print( "\t target resolution was '"+str(parameters.target_resolution)+"'")
            print( "\t image resolution was '" + str(im.voxelsize) + "'")
            print( "Exiting.")
            sys.exit( 1 )

        if monitoring.verbose >= 2:
            print( "    .. resampling '"+theImages[i].split(os.path.sep)[-1]+ "' at "+str(resampling_resolution))
        cpp_wrapping.applyTrsfCLI(theImages[i], resImages[i], theTrsf=None, templateImage=None,
                                  voxelsize=resampling_resolution, nearest=False, monitoring=monitoring )


    pass





def fusionProcess( experiment, environment, parameters ):

    #
    # make sure that the result directory exists
    #

    if not os.path.isdir( environment.path_fuse_exp ):
        os.makedirs( environment.path_fuse_exp )

    if (monitoring.verbose > 1):
        print('')

    #
    # loop over acquisitions
    #

    for timePoint in range( experiment.firstTimePoint, experiment.lastTimePoint+1, experiment.deltaTimePoint):

        #
        # result image
        #

        fusedImage = nomenclature.replaceTIME(environment.path_fuse_exp_files, timePoint)

        if (monitoring.verbose > 1):
            print('... fusion of time ' + str(timePoint))

        if os.path.isfile(fusedImage):
            if not monitoring.forceResultsToBeBuilt:
                if (monitoring.verbose > 1):
                    print('    already existing')
                continue
            else:
                if (monitoring.verbose > 1):
                    print('    already existing, but forced')

        #
        # start processing
        #

        starttime = time.time()


        #
        # directory for auxiliary files
        #
        temporary_paths = []


        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_1" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_2" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_3" ) )
        temporary_paths.append( os.path.join( environment.path_fuse_exp, "TEMP_$TIME", "ANGLE_4" ) )
        temporary_paths.append(os.path.join(environment.path_fuse_exp, "TEMP_$TIME"))

        for i in range( 0, len(temporary_paths) ):
            temporary_paths[i] = nomenclature.replaceTIME( temporary_paths[i], timePoint )
            if not os.path.isdir(temporary_paths[i]):
                os.makedirs(temporary_paths[i])


        #
        # get image file names
        # - may involve unzipping and conversion
        #
        if (monitoring.verbose > 1):
            print('    get original images' )
        images=[]

        images.append( _readImageName( environment.path_angle1, temporary_paths[0],
                        nomenclature.replaceTIME(environment.path_angle1_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle2, temporary_paths[1],
                        nomenclature.replaceTIME(environment.path_angle2_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle3, temporary_paths[2],
                        nomenclature.replaceTIME(environment.path_angle3_files,timePoint), parameters.acquisition_resolution ) )
        images.append( _readImageName( environment.path_angle4, temporary_paths[3],
                        nomenclature.replaceTIME(environment.path_angle4_files,timePoint), parameters.acquisition_resolution ) )

        print images


        #
        #
        #
        if (monitoring.verbose > 1):
            print('    fuse images' )

        fusionImages( images, fusedImage, temporary_paths, environment, parameters )










        if monitoring.keepTemporaryFiles == False:
            os.rmdir(temporary_path)
        #
        # end processing
        #

        endtime = time.time()
        if monitoring.verbose > 1:
            print( '    computation time = '+str(endtime-starttime)+ 'sec')
            print( '' )


    return











