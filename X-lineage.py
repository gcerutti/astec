#!/usr/bin/python2.7

# import os
import cPickle as pkl
import time
import sys

from argparse import ArgumentParser
from xml.etree import ElementTree

#
# local imports
# add ASTEC subdirectory
#


# import ASTEC.commonTools as commonTools
# import ASTEC.TEMPLATE as TEMPLATE
# import ASTEC.nomenclature as nomenclature
# from ASTEC.CommunFunctions.cpp_wrapping import path_to_vt

import ASTEC.xmlTools as xmlTools

#
#
#
#
#


def _set_options(my_parser):
    proc = "_set_options"
    if not isinstance(my_parser, ArgumentParser):
        print proc + ": argument is not of type ArgumentParser"
        return
    #
    # common parameters
    #

    my_parser.add_argument('-ipl', '--input-pkl-lineage',
                           action='store', dest='inputPklLineage', const=None,
                           help='pkl file containing the lineage')

    my_parser.add_argument('-ixl', '--input-xml-lineage',
                           action='store', dest='inputXmlLineage', const=None,
                           help='xml file containing the lineage')

    my_parser.add_argument('-opl', '--output-pkl-lineage',
                           action='store', dest='outputPklLineage', const=None,
                           help='pkl file containing the lineage')

    my_parser.add_argument('-oxl', '--output-xml-lineage',
                           action='store', dest='outputXmlLineage', const=None,
                           help='xml file containing the lineage')

    my_parser.add_argument('-ext', '--extract-feature',
                           action='store', nargs='*', dest='outputFeatures', const=None,
                           help="features to be extracted from the lineage: 'lineage', 'h_min', 'volume', 'sigma'" +
                                ",'label_in_time', 'barycenter', 'principal-value', 'principal-vector', 'fate'" +
                                ", 'all-cells', 'name', 'history', 'contact'")

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

    #
    #
    #

    start_time = time.localtime()

    #
    # reading command line arguments
    #

    parser = ArgumentParser(description='X-lineage')
    _set_options(parser)
    args = parser.parse_args()


    #
    # read input file
    #

    lineagedict = {}

    if args.inputXmlLineage is not None:
        xmltree = ElementTree.parse(args.inputXmlLineage)
        lineagedict = xmlTools.xml2dict(xmltree)

        print str(lineagedict)
        # sys.exit()

    elif args.inputPklLineage is not None:
        lineagefile = open(args.inputPklLineage, 'r')
        lineagedict = pkl.load(lineagefile)
        lineagefile.close()
        # print str(lineagedict)

    else:
        print "error: no input file"

    #
    # select features if requiread
    #

    outputdict = {}

    if args.outputFeatures is not None:

        #
        # search for required features
        #

        for feature in args.outputFeatures:

            # print "search feature '" + str(feature) + "'"
            target_key = xmlTools.keydictionary[feature]

            for searchedkey in target_key['input_keys']:
                if searchedkey in lineagedict:
                    # print "found feature '" + str(ok) + "'"
                    outputdict[target_key['output_key']] = lineagedict[searchedkey]
                    break
            else:
                print "error: feature '" + str(feature) + "' not found in dictionary"

    else:

        #
        # copy dictionary
        #

        # print "copy dictionary"
        outputdict = lineagedict

    if outputdict == {}:
        print "error: empty input dictionary ?! ... exiting"
        sys.exit()

    #
    # produces outputs
    #

    if args.outputPklLineage is not None:
        lineagefile = open(args.outputPklLineage, 'w')
        pkl.dump(outputdict, lineagefile)
        lineagefile.close()

    if args.outputXmlLineage is not None:
        xmltree = xmlTools.dict2xml(outputdict)
        xmltree.write(args.outputXmlLineage)


    endtime = time.localtime()

    print '# Total execution time = '+str(time.mktime(endtime)-time.mktime(start_time))+' sec\n'

