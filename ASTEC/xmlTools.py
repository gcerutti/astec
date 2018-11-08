import xml.etree.ElementTree as ElementTree
import numpy as np
import ast

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


keydictionary = {'lineage': {'output_key': 'lineage_tree',
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


#
#
#


def _set_xml_element_text(element, value):
    """

    :param element:
    :param value:
    :return:
    """
    proc = "_set_xml_element_text"

    if type(value) in (int, float, str, np.float64, np.ndarray):
        # print "   - value is '" + str(value) + "'"
        element.text = str(value)

    elif type(value) == dict:
        # print proc + ": type is dict"
        for k, v in value.iteritems():
            _dict2xml(element, k, v)
    elif type(value) == list:
        if type(value[0]) in (int, float, str, np.float64, np.ndarray):
            element.text = str(value)
        else:
            pass
    else:
        print proc + ": type '" + str(type(value)) + "' not handled yet, uncomplete translation"


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
        print proc + ": error, input is of type '" + str(type(dictionary)) + "'"
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
        print proc + ": error, empty dictionary ?!"
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

    dictionary = {}

    for child in root:

        # print child.tag
        # print "length = '" + str(len(list(child))) + "'"

        key = child.tag

        if child.tag == 'cell':
            key = np.int64(child.attrib['cell-id'])

        if len(list(child)) == 0:
            # print "child.text = '" + str(child.text) + "'"
            dictionary[key] = ast.literal_eval(child.text)
        else:
            dictionary[key] = _set_dictionary_value(child)

        # print "'" + str(key) + "' = '" + str(dictionary[key]) + "'"

    return dictionary



def xml2dict(tree):
    """

    :param tree:
    :return:
    """

    root = tree.getroot()

    dictionary = {}

    for k, v in keydictionary.iteritems():
        if root.tag == v['output_key']:
            dictionary[str(root.tag)] = _set_dictionary_value(root)
            break
    else:
        for child in root:
            dictionary[str(child.tag)] = _set_dictionary_value(child)

    return dictionary











