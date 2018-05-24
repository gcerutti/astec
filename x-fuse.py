#!/usr/bin/python2.7

import os, sys
import time
from argparse import ArgumentParser


#
# local imports
# add ASTEC subdirectory
#

# assert os.path.isdir(os.path.join(os.path.dirname(__file__),"ASTEC")), \
#    'ASTEC directory not found'
# sys.path.append(os.path.join(os.path.dirname(__file__),"ASTEC"))


import ASTEC.commonTools as commonTools
import ASTEC.XFUSION as XFUSION




#
#
#
#
#

def _setOptions( parser ):
    if not isinstance( parser, ArgumentParser ):
        print ( 'not an ArgumentParser' )
        return
    #
    # common parameters
    #

    parser.add_argument( '-p', '--parameters',
                         action='store', dest='parameterFile', const=None,
                         help='python file containing parameters definition' )
    parser.add_argument('-e', '--embryo-rep',
                        action='store', dest='embryoPath', const=None,
                        help='path to the embryo data' )

    #
    # control parameters
    #

    parser.add_argument('-k', '--keep-temporary-files',
                        action='store_const', dest='keepTemporaryFiles',
                        default=False, const=True,
                        help='keep temporary files' )

    parser.add_argument( '-v', '--verbose',
                         action='count', dest='verbose', default=1,
                         help='incrementation of verboseness' )
    parser.add_argument( '-nv', '--no-verbose',
                         action='store_const', dest='verbose', const=0,
                         help='no verbose at all')
    parser.add_argument( '-d', '--debug',
                         action='count', dest='debug', default=0,
                         help='incrementation of debug level' )
    parser.add_argument( '-nd', '--no-debug',
                         action='store_const', dest='debug', const=0,
                         help='no debug information' )

    return





#
#
# main 
#
#
if __name__ == '__main__':


    #
    # initialization
    #
    starttime = time.localtime()
    monitoring = commonTools.Monitoring()
    experiment = commonTools.Experiment()
    parameters = XFUSION.FusionParameters()
    environment = XFUSION.FusionEnvironment()

    #
    # reading command line arguments
    #
    parser = ArgumentParser( description = 'Fusion of multiple acquisitions' )
    _setOptions( parser )
    args = parser.parse_args()

    monitoring.updateFromArgs( args )
    experiment.updateFromArgs( args )

    #
    # reading parameter files
    # and updating parameters
    #
    parameterFile= commonTools.getParameterFile( args.parameterFile )
    experiment.updateFromFile( parameterFile )
    parameters.updateFromFile( parameterFile )
    environment.updateFromFile( parameterFile, starttime )

    #
    # make fusion directory and subdirectory if required
    #
    if not os.path.isdir( environment.path_fuse_exp ):
        os.makedirs( environment.path_fuse_exp )



    #
    # write information in history file
    #
    commonTools.writeHistoryInformation( environment.path_history_file,
                                         experiment,
                                         parameterFile,
                                         starttime,
                                         os.path.dirname(__file__) )

    #
    # write log file
    #
    experiment.writeParameters( environment.path_log_file )
    environment.writeParameters( environment.path_log_file )
    parameters.writeParameters( environment.path_log_file )

    #
    # copy parameter file
    #
    commonTools.copyDateStampedFile( parameterFile, environment.path_fuse_exp, starttime )

    #
    # prepare processing
    #
    XFUSION.monitoring.copy( monitoring )
    XFUSION.monitoring.logfile = environment.path_log_file


    #
    # processing
    #
    XFUSION.fusionProcess(experiment, environment, parameters)


    #
    # end of execution
    # write execution time in both log and history file
    #
    endtime = time.localtime()
    with open(environment.path_log_file, 'a') as logfile:
        logfile.write("\n")
        logfile.write('Total execution time = '+str(time.mktime(endtime)-time.mktime(starttime))+' sec\n')
        logfile.write("\n")

    with open(environment.path_history_file, 'a') as logfile:
        logfile.write('# Total execution time = '+str(time.mktime(endtime)-time.mktime(starttime))+' sec\n')
        logfile.write("\n\n")
