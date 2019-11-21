#!/usr/bin/python2.7

import os
import time
from argparse import ArgumentParser

#
# local imports
# add ASTEC subdirectory
#


import ASTEC.commonTools as commonTools
import ASTEC.INTRAREGISTRATION as INTRAREG
from ASTEC.CommunFunctions.cpp_wrapping import path_to_vt


#
#
#
#
#


def _set_options(my_parser):
    """


#
#
# main call
#
#


if __name__ == '__main__':
    main()

#
#
# main call
#
#


if __name__ == '__main__':
    main()
    :param my_parser:
    :return:
    """
    proc = "_set_options"
    if not isinstance(my_parser, ArgumentParser):
        print proc + ": argument is not of type ArgumentParser"
        return
    #
    # common parameters
    #

    my_parser.add_argument('-p', '--parameters',
                           action='store', dest='parameterFile', const=None,
                           help='python file containing parameters definition')
    my_parser.add_argument('-e', '--embryo-rep',
                           action='store', dest='embryo_path', const=None,
                           help='path to the embryo data')

    #
    # control parameters
    #

    my_parser.add_argument('-k', '--keep-temporary-files',
                           action='store_const', dest='keepTemporaryFiles',
                           default=False, const=True,
                           help='keep temporary files')

    my_parser.add_argument('-f', '--force',
                           action='store_const', dest='forceResultsToBeBuilt',
                           default=False, const=True,
                           help='force building of results')

    my_parser.add_argument('-v', '--verbose',
                           action='count', dest='verbose', default=2,
                           help='incrementation of verboseness')
    my_parser.add_argument('-nv', '--no-verbose',
                           action='store_const', dest='verbose', const=0,
                           help='no verbose at all')
    my_parser.add_argument('-d', '--debug',
                           action='count', dest='debug', default=0,
                           help='incrementation of debug level')
    my_parser.add_argument('-nd', '--no-debug',
                           action='store_const', dest='debug', const=0,
                           help='no debug information')

    #
    # specific args
    #
    my_parser.add_argument('-t', '--reference-transformation',
                           action='store', dest='reference_transformation_file', const=None,
                           help='resampling transformation to be applied to the reference image')
    my_parser.add_argument('-a', '--reference-angles',
                           action='store', dest='reference_transformation_angles', const=None,
                           help='angles wrt to X, Y and Z axis to build the reference resampling transformation,' +
                                'it is a string formed by the axis name then the angles, eg "X 70 Z -120"')
    return


#
#
# main function
#
#


def main():

    #
    # initialization
    #
    start_time = time.localtime()
    monitoring = commonTools.Monitoring()
    experiment = commonTools.Experiment()
    parameters = INTRAREG.IntraRegParameters()

    #
    # reading command line arguments
    #
    parser = ArgumentParser(description='Fused sequence intra-registration')
    _set_options(parser)
    args = parser.parse_args()

    monitoring.update_from_args(args)
    experiment.update_from_args(args)
    parameters.update_from_args(args)

    #
    # reading parameter files
    # and updating parameters
    #
    parameterFile = commonTools.get_parameter_file(args.parameterFile)
    experiment.update_from_file(parameterFile)
    experiment.update_from_stage("intrareg", __file__, start_time)

    if not os.path.isdir(experiment.path_logdir):
        os.makedirs(experiment.path_logdir)

    parameters.update_from_file(parameterFile)

    #
    # write history information in history file
    #
    commonTools.write_history_information(experiment.path_history_file,
                                          experiment,
                                          parameterFile,
                                          start_time,
                                          os.path.dirname(__file__),
                                          path_to_vt())

    #
    # define log file
    # and write some information
    #
    monitoring.logfile = experiment.path_log_file
    INTRAREG.monitoring.copy(monitoring)

    monitoring.write_parameters(monitoring.logfile)
    experiment.write_parameters(monitoring.logfile)
    parameters.write_parameters(monitoring.logfile)

    #
    # copy parameter file
    #
    commonTools.copy_date_stamped_file(parameterFile, experiment.path_logdir, start_time)

    #
    # processing
    #
    INTRAREG.intraregistration_control(experiment, parameters)

    #
    # end of execution
    # write execution time in both log and history file
    #
    endtime = time.localtime()
    with open(experiment.path_log_file, 'a') as logfile:
        logfile.write("\n")
        logfile.write('Total execution time = '+str(time.mktime(endtime)-time.mktime(start_time))+' sec\n')
        logfile.write("\n")

    with open(experiment.path_history_file, 'a') as logfile:
        logfile.write('# Total execution time = '+str(time.mktime(endtime)-time.mktime(start_time))+' sec\n')
        logfile.write("\n\n")


#
#
# main call
#
#


if __name__ == '__main__':
    main()