
import os, imp

import commonTools
import nomenclature

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

            logfile.write('- path_fuse = ' + str(self.path_fuse)+'\n')
            logfile.write('- path_fuse_exp = ' + str(self.path_fuse_exp)+'\n')
            logfile.write('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files)+'\n')

            logfile.write('- path_history_file = ' + str(self.path_history_file)+'\n')
            logfile.write('- path_log_file = ' + str(self.path_log_file)+'\n')
            logfile.write("\n")
        return

    def printParameters( self ):
        print("\n")
        print('FusionEnvironment\n')
        print('- path_angle1 = ' + str(self.path_angle1)+'\n')
        print('- path_angle2 = ' + str(self.path_angle2)+'\n')
        print('- path_angle3 = ' + str(self.path_angle3)+'\n')
        print('- path_angle4 = ' + str(self.path_angle4)+'\n')

        print('- path_fuse = ' + str(self.path_fuse)+'\n')
        print('- path_fuse_exp = ' + str(self.path_fuse_exp)+'\n')
        print('- path_fuse_exp_files = ' + str(self.path_fuse_exp_files)+'\n')

        print('- path_history_file = ' + str(self.path_history_file)+'\n')
        print('- path_log_file = ' + str(self.path_log_file)+'\n')
        print("\n")





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
        print("\n")
        print( 'FusionParameters\n')
        print( '- acquisition_orientation = '+str(self.acquisition_orientation)+'\n' )
        print( '- acquisition_mirrors     = '+str(self.acquisition_mirrors)+'\n' )
        print( '- acquisition_resolution  = '+str(self.acquisition_resolution)+'\n' )
        print( '- acquisition_delay       = ' + str(self.acquisition_delay)+'\n' )
        print( '- target_resolution  = '+str(self.target_resolution)+'\n' )
        print( '- acquisition_cropping = '+str(self.acquisition_cropping)+'\n' )
        print( '- acquisition_cropping_margin_x_0 = '+str(self.acquisition_cropping_margin_x_0)+'\n' )
        print( '- acquisition_cropping_margin_x_1 = '+str(self.acquisition_cropping_margin_x_1)+'\n' )
        print( '- acquisition_cropping_margin_y_0 = '+str(self.acquisition_cropping_margin_y_0)+'\n' )
        print( '- acquisition_cropping_margin_y_1 = '+str(self.acquisition_cropping_margin_y_1)+'\n' )
        print( '- fusion_cropping = '+str(self.fusion_cropping)+'\n' )
        print( '- fusion_cropping_margin_x_0 = '+str(self.fusion_cropping_margin_x_0)+'\n' )
        print( '- fusion_cropping_margin_x_1 = '+str(self.fusion_cropping_margin_x_1)+'\n' )
        print( '- fusion_cropping_margin_y_0 = '+str(self.fusion_cropping_margin_y_0)+'\n' )
        print( '- fusion_cropping_margin_y_1 = '+str(self.fusion_cropping_margin_y_1)+'\n' )
        print("\n")

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

