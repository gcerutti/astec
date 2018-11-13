
import os

import cPickle as pkl
import xml.etree.ElementTree as ElementTree
import numpy as np

import commonTools

#
#
#
#
#

monitoring = commonTools.Monitoring()


########################################################################################
#
# key correspondences
#
# Examples
# - from 'full_properties_Samson_MN20.pkl', keys are
#   ['volumes information',
#    'Cells labels in time',
#    'Barycenters',
#    'Fate',
#    'All Cells',
#    'Principal values',
#    'Names',
#    'cell_cell_contact_information',
#    'Lineage tree',
#    'Cells history',
#    'Principal vectors']
# - from 'new_lineage_tree_MN20.pkl', keys are cell labels
# - from '140317-Patrick-St8_seg_lineage.pkl', keys are
#   ['h_mins_information', 'lin_tree', 'volumes_information', 'sigmas_information']
#
#
########################################################################################


keydictionary = {'lineage': {'output_key': 'lin_tree',
                             'input_keys': ['lineage_tree', 'lin_tree', 'Lineage tree']},
                 'h_min': {'output_key': 'cell_h_min',
                           'input_keys': ['cell_h_min', 'h_mins_information']},
                 'volume': {'output_key': 'cell_volume',
                            'input_keys': ['cell_volume', 'volumes_information', 'volumes information']},
                 'sigma': {'output_key': 'cell_sigma',
                           'input_keys': ['cell_sigma', 'sigmas_information']},
                 'label_in_time': {'output_key': 'cell_labels_in_time',
                                   'input_keys': ['cell_labels_in_time', 'Cells labels in time']},
                 'barycenter': {'output_key': 'cell_barycenter',
                                'input_keys': ['cell_barycenter', 'Barycenters']},
                 'fate': {'output_key': 'cell_fate',
                          'input_keys': ['cell_fate', 'Fate']},
                 'all-cells': {'output_key': 'all_cells',
                               'input_keys': ['all_cells', 'All Cells']},
                 'principal-value': {'output_key': 'cell_principal_values',
                                     'input_keys': ['cell_principal_values', 'Principal values']},
                 'name': {'output_key': 'cell_name',
                          'input_keys': ['cell_name', 'Names']},
                 'contact': {'output_key': 'cell_contact_surface',
                             'input_keys': ['cell_contact_surface', 'cell_cell_contact_information']},
                 'history': {'output_key': 'cell_history',
                             'input_keys': ['cell_history', 'Cells history']},
                 'principal-vector': {'output_key': 'cell_principal_vectors',
                                      'input_keys': ['cell_principal_vectors', 'Principal vectors']}}


########################################################################################
#
# to translate a dictionary into XML
#
########################################################################################

#
# from stackoverflow.com
# questions/3095434/inserting-newlines-in-xml-file-generated-via-xml-etree-elementtree-in-python
#
#
# used for pretty printing
#


def _indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            _indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i


########################################################################################
#
# types
# 'lineage':  'lineage_tree' :
#     dictionary de liste de int
#     lineage_tree.590002 = <type 'list'>
# 'h_min' : 'cell_h_min' :
# 'volume' : 'cell_volume' :
#     dictionary de int
#     cell_volume.590002 = <type 'int'>
#     590002: 236936
# 'sigma': 'cell_sigma':
# 'label_in_time': 'cell_labels_in_time'
#     dictionary de liste de numpy.int64
#     cell_labels_in_time.1 = <type 'list'>
#     1: [10002, 10003, 10004, ..., 10082]
# 'barycenter': 'cell_barycenter'
#     dictionary de numpy.ndarray de numpy.float64
#     cell_barycenter.590002 = <type 'numpy.ndarray'>
#     590002: array([ 258.41037242,  226.74975943,  303.67167927])
# 'fate': 'cell_fate'
#     dictionary de str
#     cell_fate.590002 = <type 'str'>
#     590002: 'Mesoderm Notochord 1'
# 'all-cells': 'all_cells'  # liste de toutes les cellules ?
#     liste de numpy.int64
#     all_cells = <type 'list'>
# 'principal-value': 'cell_principal_values'
#     dictionary de liste de numpy.float64
#     cell_principal_values.590002 = <type 'list'>
#     590002: [1526.0489371146978, 230.60881177650205, 91.063513300019849]
# 'name': 'cell_name'
#     dictionary de str
#     cell_name.590002 = <type 'str'>
#     590002: 'a9.0012_'
# 'contact': 'cell_contact_surface',
#     dictionary de dictionary de int
#     cell_contact_surface.590002.590019 = <type 'int'>
#     590002: {590001: 1808, 590003: 1436, 590004: 5012, ..., 590225: 2579}
# 'history': 'cell_history'
#     dictionary de numpy.ndarray de numpy.int64
#     cell_history.590002 = <type 'numpy.ndarray'>
#     590002: array([510002, 520002, 530002, 540002, 550002, 560002, 570002, 580002,
#         590002, 600002, 610002, 620002, 630002, 640002, 650002, 660002,
#         670002, 680002, 690002, 700002, 710002, 720002, 730002, 740002,
#         750002, 760002, 770002, 780002, 790002, 800002, 810002, 820002,
#         830002, 840002, 850002, 860002, 870002, 880002])
# 'principal-vector': 'cell_principal_vectors'    # liste de numpy.ndarray
#     dictionary de liste de numpy.ndarray de numpy.float64
#     cell_principal_vectors.590002 = <type 'list'>
#     590002: [array([ 0.17420991, -0.74923203,  0.63898534]),
#         array([-0.24877611,  0.59437038,  0.7647446 ]),
#         array([ 0.95276511,  0.29219037,  0.08284582])]
#
########################################################################################

def _set_xml_element_text(element, value):
    """

    :param element:
    :param value:
    :return:
    """
    proc = "_set_xml_element_text"

    #
    # dictionary : recursive call
    #   dictionary element may be list, int, numpy.ndarray, str
    # list : may be int, numpy.int64, numpy.float64, numpy.ndarray
    #

    print "value is '" + str(value) + "'"

    if type(value) == dict:
        # print proc + ": type is dict"
        for k, v in value.iteritems():
            _dict2xml(element, k, v)

    elif type(value) == list:

        #
        # empty list
        #

        if len(value) == 0:
            element.text = repr(value)
        #
        # 'lineage', 'label_in_time', 'all-cells', 'principal-value'
        #

        elif type(value[0]) in (int, float, np.int64, np.float64):
            # element.text = str(value)
            element.text = repr(value)

        #
        # 'principal-vector' case
        #  liste de numpy.ndarray de numpy.float64
        #
        elif type(value[0]) == np.ndarray:
            text = "["
            for i in range(len(value)):
                # text += str(list(value[i]))
                text += repr(list(value[i]))
                if i < len(value)-1:
                    text += ", "
                    if i > 0 and i % 10 == 0:
                        text += "\n  "
            text += "]"
            element.text = text
            del text

        else:
            monitoring.to_log_and_console(proc + ": error, element list type ('" + str(type(value[0]))
                                          + "') not handled yet")

    #
    # 'barycenter', 'cell_history'
    #
    elif type(value) == np.ndarray:
        # element.text = str(list(value))
        element.text = repr(list(value))

    #
    # 'volume', 'contact'
    #
    elif type(value) == int:
        # element.text = str(value)
        element.text = repr(value)

    #
    # 'fate', 'name'
    #
    elif type(value) == str:
        element.text = repr(value)

    else:
        monitoring.to_log_and_console(proc + ": element type '" + str(type(value))
                                      + "' not handled yet, uncomplete translation")


#
#
#


def _dict2xml(parent, tag, value):
    """

    :param parent:
    :param tag:
    :param value:
    :return:
    """

    #
    # integers can not be XML tags
    #
    if type(tag) in (int, np.int64):
        child = ElementTree.Element('cell', attrib={'cell-id': str(tag)})
    else:
        child = ElementTree.Element(str(tag))

    _set_xml_element_text(child, value)

    parent.append(child)
    return parent


#
# procedure d'appel et creation de la racine
#


def dict2xml(dictionary, defaultroottag='data'):
    """

    :param dictionary:
    :param defaultroottag:
    :return:
    """

    proc = "dict2xml"

    if type(dictionary) is not dict:
        monitoring.to_log_and_console(proc + ": error, input is of type '" + str(type(dictionary)) + "'")
        return None

    #
    # s'il n'y a qu'un seul element dans le dictionnaire, on appelle la racine
    # d'apres celui-ci (eg lineage_tree), sinon on cree une racine qui contient
    # tous les elements
    #

    if len(dictionary) == 1:

        roottag = dictionary.keys()[0]
        root = ElementTree.Element(roottag)
        _set_xml_element_text(root, dictionary[roottag])

    elif len(dictionary) > 1:

        root = ElementTree.Element(defaultroottag)
        for k, v in dictionary.iteritems():
            _dict2xml(root, k, v)

    else:
        monitoring.to_log_and_console(proc + ": error, empty dictionary ?!")
        return None

    _indent(root)
    tree = ElementTree.ElementTree(root)

    return tree


########################################################################################
#
# to translate a XML tree into dictionary
#
########################################################################################


def _set_dictionary_value(root):
    """

    :param root:
    :return:
    """
    if len(root) == 0:

        #
        # pas de branche, on renvoie la valeur
        #

        # return ast.literal_eval(root.text)
        return eval(root.text)

    else:

        dictionary = {}

        for child in root:

            # print "child.tag=" + str(child.tag)
            # print "len(child)=" + str(len(child))
            # print "child.text=" + str(child.text)

            key = child.tag
            if child.tag == 'cell':
                key = np.int64(child.attrib['cell-id'])
            dictionary[key] = _set_dictionary_value(child)

    return dictionary


def xml2dict(tree):
    """

    :param tree:
    :return:
    """

    proc = "xml2dict"

    root = tree.getroot()

    dictionary = {}

    for k, v in keydictionary.iteritems():

        if root.tag == v['output_key']:
            monitoring.to_log_and_console(proc + ": process root.tag = '" + str(root.tag) + "'", 3)
            dictionary[str(root.tag)] = _set_dictionary_value(root)
            break
    else:
        for child in root:
            monitoring.to_log_and_console(proc + ": process child.tag = '" + str(child.tag) + "'", 3)
            dictionary[str(child.tag)] = _set_dictionary_value(child)

    return dictionary


########################################################################################
#
# to read a set of files into a dictionary
#
########################################################################################


#
# update dictionary from what has been read
#

def _update_read_dictionary(propertiesdict, tmpdict, filename):
    """

    :param propertiesdict:
    :param tmpdict:
    :return:
    """
    proc = "_update_read_dictionary"
    unknownkeys = []

    for tmpkey in tmpdict.keys():
        foundkey = False

        for k in keydictionary.keys():
            # print "       compare '" + str(tmpkey) + "' with '" + str(k) + "'"
            if tmpkey in keydictionary[k]['input_keys']:
                outputkey = keydictionary[k]['output_key']
                monitoring.to_log_and_console("   ... recognized key '" + str(outputkey) + "'", 2)
                #
                # update if key already exists, else just create the dictionary entry
                #
                if outputkey in propertiesdict.keys():
                    if type(propertiesdict[outputkey]) is dict and type(tmpdict[tmpkey]) is dict:
                        propertiesdict[outputkey].update(tmpdict[tmpkey])
                    elif type(propertiesdict[outputkey]) is list and type(tmpdict[tmpkey]) is list:
                        propertiesdict[outputkey].append(tmpdict[tmpkey])
                    else:
                        monitoring.to_log_and_console(proc + ": error, can not update property '" + str(outputkey)
                                                      + "'")
                else:
                    propertiesdict[outputkey] = tmpdict[tmpkey]
                foundkey = True
                break

        if foundkey is False:
            unknownkeys.append(tmpkey)

    if len(unknownkeys) > 0 and len(unknownkeys) == len(tmpdict.keys()):
        #
        # no key was found
        # it is assumed it's a lineage tree: add some test here ?
        #
        monitoring.to_log_and_console("   ... assume '" + str(filename) + "' is a lineage", 1)
        outputkey = keydictionary['lineage']['output_key']
        if outputkey in propertiesdict.keys():
            if type(propertiesdict[outputkey]) is dict and type(tmpdict) is dict:
                propertiesdict[outputkey].update(tmpdict)
            else:
                monitoring.to_log_and_console(proc + ": error, can not update property '" + str(outputkey) + "'")
        else:
            propertiesdict[outputkey] = tmpdict

    elif len(unknownkeys) > 0:
        #
        # some keys were found
        #
        monitoring.to_log_and_console("   ... unrecognized key(s) are '" + str(unknownkeys) + "'", 1)

    return propertiesdict


#
# types issued from the reading of xml files may be erroneous
# fix it
#

def _set_types_from_xml(propertiesdict):
    """

    :param propertiesdict:
    :return:
    """

    if propertiesdict == {}:
        return {}

    if 'cell_barycenter' in propertiesdict.keys():
        monitoring.to_log_and_console("   ... translate types of 'cell_barycenter'", 1)
        for c in propertiesdict['cell_barycenter'].keys():
            propertiesdict['cell_barycenter'][c] = np.array(propertiesdict['cell_barycenter'][c])

    if 'cell_history' in propertiesdict.keys():
        monitoring.to_log_and_console("   ... translate types of 'cell_history'", 1)
        for c in propertiesdict['cell_history'].keys():
            propertiesdict['cell_history'][c] = np.array(propertiesdict['cell_history'][c])

    if 'cell_principal_vectors' in propertiesdict.keys():
        monitoring.to_log_and_console("   ... translate types of 'cell_principal_vectors'", 1)
        for c in propertiesdict['cell_principal_vectors'].keys():
            for v in range(len(propertiesdict['cell_principal_vectors'][c])):
                propertiesdict['cell_principal_vectors'][c][v] \
                    = np.array(propertiesdict['cell_principal_vectors'][c][v])

    return propertiesdict


#
#
#

def read_dictionary(inputfilenames):
    """

    :param inputfilenames:
    :return:
    """
    proc = 'read_dictionary'

    if inputfilenames is None:
        monitoring.to_log_and_console(proc + ": error, no input files")
        return {}

    propertiesdict = {}

    #
    # read xml files
    #

    for filename in inputfilenames:

        if not os.path.isfile(filename):
            monitoring.to_log_and_console(proc + ": error, file '" + str(filename) + "' does not exist")
            continue

        if filename[len(filename) - 3:len(filename)] == "xml":
            monitoring.to_log_and_console("... reading '" + str(filename) + "'", 1)
            inputxmltree = ElementTree.parse(filename)
            tmpdict = xml2dict(inputxmltree)
            propertiesdict = _update_read_dictionary(propertiesdict, tmpdict, filename)
            del tmpdict

    #
    # translation of xml may take place here
    #
    propertiesdict = _set_types_from_xml(propertiesdict)

    for filename in inputfilenames:

        if filename[len(filename) - 3:len(filename)] == "pkl":
            monitoring.to_log_and_console("... reading '" + str(filename) + "'", 1)
            inputfile = open(filename, 'r')
            tmpdict = pkl.load(inputfile)
            inputfile.close()
            propertiesdict = _update_read_dictionary(propertiesdict, tmpdict, filename)
            del tmpdict

    #
    #
    #

    for filename in inputfilenames:
        if filename[len(filename) - 3:len(filename)] == "xml":
            continue
        elif filename[len(filename) - 3:len(filename)] == "pkl":
            continue
        else:
            monitoring.to_log_and_console(proc + ": error: extension not recognized for '" + str(filename) + "'")

    return propertiesdict


########################################################################################
#
# comparison of two dictionaries
#
########################################################################################

def compare_dictionaries(d1, d2, features, name1, name2):
    """

    :param d1:
    :param d2:
    :param features:
    :param name1:
    :param name2:
    :return:
    """

    monitoring.to_log_and_console("\n", 1)
    monitoring.to_log_and_console("... comparison between '" + str(name1) + "' and '" + str(name2) + "'", 1)

    #
    # 1. find common keys
    #

    unpairedkeys1 = []
    unpairedkeys2 = []
    pairedkeys = []
    unrecognizedkeys1 = []
    unrecognizedkeys2 = []

    for k1 in d1.keys():

        #
        # loop on known dictionary
        #
        recognizedkey = False
        for k in keydictionary.keys():

            if k1 in keydictionary[k]['input_keys']:
                recognizedkey = True
                pairedkey = False
                #
                # got it, try to pair it
                #
                for k2 in d2.keys():
                    if k2 in keydictionary[k]['input_keys']:
                        pairedkey = True
                        pairedkeys.append([k1,k2])
                        break
                if pairedkey is False:
                    unpairedkeys1.append(k1)
                break

        if recognizedkey is False:
            unrecognizedkeys1.append(k1)

    #
    #
    #

    for k2 in d2.keys():

        #
        # loop on known dictionary
        #
        recognizedkey = False
        for k in keydictionary.keys():
            if k2 in keydictionary[k]['input_keys']:
                recognizedkey = True
                pairedkey = False
                #
                # got it, try to pair it
                #
                for k1 in d1.keys():
                    if k1 in keydictionary[k]['input_keys']:
                        pairedkey = True
                        # pairedkeys.append([k1,k2])
                        break
                if pairedkey is False:
                    unpairedkeys2.append(k2)
                break

        if recognizedkey is False:
            unrecognizedkeys2.append(k2)

    #
    #
    #

    print_summary = True

    #
    #
    #

    if print_summary is True:
        monitoring.to_log_and_console("    found keys", 1)
        if len(pairedkeys) > 0:
            monitoring.to_log_and_console("    ... common keys to '" + str(name1) + "' and '" + str(name2) + "'", 1)
            for p in pairedkeys:
                monitoring.to_log_and_console("        " + str(p[0]) + " <-> " + str(p[1]), 1)
        if len(unpairedkeys1) > 0:
            monitoring.to_log_and_console("    ... keys in '" + str(name1) + "' and not in '" + str(name2) + "'", 1)
            for k in unpairedkeys1:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unpairedkeys2) > 0:
            monitoring.to_log_and_console("    ... keys not in '" + str(name1) + "' and in '" + str(name2) + "'", 1)
            for k in unpairedkeys2:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unrecognizedkeys1) > 0:
            monitoring.to_log_and_console("    ... keys in '" + str(name1) + "' not recognized", 1)
            for k in unrecognizedkeys1:
                monitoring.to_log_and_console("        " + str(k), 1)
        if len(unrecognizedkeys2) > 0:
            monitoring.to_log_and_console("    ... keys in '" + str(name2) + "' not recognized", 1)
            for k in unrecognizedkeys2:
                monitoring.to_log_and_console("        " + str(k), 1)


    if features is None or len(features) == 0:
        pass

    monitoring.to_log_and_console("\n", 1)
    return


########################################################################################
#
#
#
########################################################################################





########################################################################################
#
# utilities for debugging, etc.
#
########################################################################################

def print_type(d, t=None, desc=None):

    if desc is None:
        desc = ""
    if t is None:
        t = ""

    if type(d) is dict:

        print "type of " + desc + " is " + str(t) + str(type(d))
        for k in d.keys():
            print_type(d[k], t + str(type(d)) + ".", desc + "." + str(k))

    elif type(d) in (list, np.array, np.ndarray):
        print "type of " + desc + " is " + str(t) + str(type(d))
        print_type(d[0], t + str(type(d)) + ".", desc + "[0]")

    else:
        print "type of " + desc + " is " + str(t) + str(type(d))

    return




