
import os
import imp
import sys

from scipy.stats.stats import pearsonr

import common
import properties as properties


#
#
#
#
#

monitoring = common.Monitoring()


########################################################################################
#
# classes
# - computation parameters
#
########################################################################################


#
#
#

class PostCorrectionParameters(object):

    ############################################################
    #
    # initialisation
    #
    ############################################################

    def __init__(self):

        ############################################################
        #
        # initialisation
        #
        ############################################################

        self.soon = False
        self.volume_minimal_value = 2000
        self.lifespan_minimal_value = 25
        self.pearson_threshold = 0.9

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('PostCorrectionParameters')

        print('- soon = ' + str(self.soon))
        print('- volume_minimal_value = ' + str(self.volume_minimal_value))
        print('- lifespan_minimal_value = ' + str(self.lifespan_minimal_value))
        print('- pearson_threshold = ' + str(self.pearson_threshold))

        print("")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('PostCorrectionParameters\n')

            logfile.write('- soon = ' + str(self.soon) + '\n')
            logfile.write('- volume_minimal_value = ' + str(self.volume_minimal_value) + '\n')
            logfile.write('- lifespan_minimal_value = ' + str(self.lifespan_minimal_value) + '\n')
            logfile.write('- pearson_threshold = ' + str(self.pearson_threshold) + '\n')

            print("")
        return

    ############################################################
    #
    # update
    #
    ############################################################

    def update_from_parameters(self, parameter_file):
        if parameter_file is None:
            return
        if not os.path.isfile(parameter_file):
            print("Error: '" + parameter_file + "' is not a valid file. Exiting.")
            sys.exit(1)

        parameters = imp.load_source('*', parameter_file)

        if hasattr(parameters, 'postcorrection_soon'):
            if parameters.postcorrection_soon is not None:
                self.soon = parameters.postcorrection_soon
        if hasattr(parameters, 'postcor_Soon'):
            if parameters.postcor_Soon is not None:
                self.soon = parameters.postcor_Soon

        if hasattr(parameters, 'volume_minimal_value'):
            if parameters.volume_minimal_value is not None:
                self.volume_minimal_value = parameters.volume_minimal_value
        if hasattr(parameters, 'postcorrection_volume_minimal_value'):
            if parameters.postcorrection_volume_minimal_value is not None:
                self.volume_minimal_value = parameters.postcorrection_volume_minimal_value
        if hasattr(parameters, 'postcor_Volume_Threshold'):
            if parameters.postcor_Volume_Threshold is not None:
                self.volume_minimal_value = parameters.postcor_Volume_Threshold

        if hasattr(parameters, 'lifespan_minimal_value'):
            if parameters.lifespan_minimal_value is not None:
                self.lifespan_minimal_value = parameters.lifespan_minimal_value
        if hasattr(parameters, 'postcorrection_lifespan_minimal_value'):
            if parameters.postcorrection_lifespan_minimal_value is not None:
                self.lifespan_minimal_value = parameters.postcorrection_lifespan_minimal_value
        if hasattr(parameters, 'postcor_lifespan_minimal_value'):
            if parameters.postcor_lifespan_minimal_value is not None:
                self.lifespan_minimal_value = parameters.postcor_lifespan_minimal_value


########################################################################################
#
# some internal procedures
#
########################################################################################

def _get_volumes(cell, direct_lineage, volume):
    """
    Return the volumes of cell n and its progeny up to division
    n : starting cell
    vol : dictionary of volumes
    lin_tree : lineage tree
    """
    tmp = []
    while len(direct_lineage.get(cell, '')) == 1:
        tmp.append(volume[cell])
        cell = direct_lineage[cell][0]
    return tmp


def _suspected_early_division(direct_lineage, reverse_lineage, parameters):
    proc = "_suspected_early_division"


def _get_branch_to_be_removed(direct_lineage, reverse_lineage, volume, experiment, parameters):

    proc = "_get_branch_to_be_removed"
    #
    # leaves are nodes without lineage
    # get the largest time point value (last leaves)
    #
    nodes = list(set(direct_lineage.keys()).union(set([v for values in direct_lineage.values() for v in values])))
    leaves = set(nodes) - set(direct_lineage.keys())
    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
    last_time = max(leaves) / 10**time_digits_for_cell_id

    #
    # a branch (ended) by a leave is candidate for deletion if
    # - it ends before the last time point (so before the end of the series), or
    # - the volume of the leave cell is too small
    #
    candidates_for_deletion = [leave for leave in leaves if (((leave/10**time_digits_for_cell_id) < last_time)
                                                             or (volume[leave] < parameters.volume_minimal_value))]
    #
    # find the branch to be deleted
    # - the last bifurcation should be as near as possible to the last time
    #
    first_branch_time_point = 0
    branch_to_be_deleted = []
    for leaf in candidates_for_deletion:
        #
        # get a branch from leaf 'leaf'
        # get back from the leaf until its parents has more than 1 daughter
        #
        le = leaf
        branch = [le]
        while len(direct_lineage.get(reverse_lineage.get(le, ''), '')) == 1:
            le = reverse_lineage[le]
            branch.append(le)
        branch.reverse()
        #
        # we already got a branch that begins later
        # note that it also compare the cell label
        # we could compare the time of the branch start and have an other criteria
        # to distinguish between two branches beginning at the same time
        #
        if min(branch) / 10**time_digits_for_cell_id < first_branch_time_point:
            continue
        if 0 < len(branch_to_be_deleted) < len(branch):
            continue
        #
        # it's a small branch, choose it if it is better
        #
        if len(branch) < parameters.lifespan_minimal_value:
            branch_to_be_deleted = branch
            first_branch_time_point = min(branch) / 10**time_digits_for_cell_id
            print(" (1) ... branche_max: " + str(min(branch_to_be_deleted)) + " / " + str(branch_to_be_deleted))
            continue
        #
        # test whether the branch begins at the very beginning (it has no parent)
        #
        if reverse_lineage.get(min(branch), '') == '':
            continue
        #
        # not a small branch, test a Pearson correlation coefficient
        #
        sisters = direct_lineage[reverse_lineage[min(branch)]]
        print("min = " + str(min(branch)) + ", sisters = " + str(sisters))
        if len(sisters) > 2:
            monitoring.to_log_and_console(str(proc) + ": cell #" + str(reverse_lineage[min(branch)]) +
                                          "divides in more than 2 branches, skip it")
            continue
        vol0 = _get_volumes(sisters[0], direct_lineage, volume)
        vol1 = _get_volumes(sisters[1], direct_lineage, volume)
        min_length = min(len(vol0), len(vol1))
        #
        # pearsonr() returns
        # - Pearson's correlation coefficient
        # - Two-tailed p-value.
        # a correlation coefficient closed to -1 means there is a negative correlation between
        # the two volume array (one is increasing, the other is decreasing)
        #
        pearson_correlation = pearsonr(vol0[:min_length-1], vol1[:min_length-1])
        print(str(min_length) + " = " + str(vol0) + " - " + str(vol1) + "=" + str(pearson_correlation))
        if pearson_correlation[0] < -parameters.pearson_threshold:
            branch_to_be_deleted = branch
            first_branch_time_point = min(branch) / 10 ** time_digits_for_cell_id
            continue
        #
        # it is not a small branch,
        # its volume evolution is not negatively correlated with its sister branch
        #
        if not parameters.soon:
            continue

    print("... branche_max: " + str(min(branch_to_be_deleted)) + " / " + str(branch_to_be_deleted))

    return []


def _remove_small_branches(direct_lineage, volume, experiment, parameters):
    proc = "_remove_small_branches"
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, PostCorrectionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    #
    #
    reverse_lineage = {v: k for k, values in direct_lineage.iteritems() for v in values}
    while True:
        branch = _get_branch_to_be_removed(direct_lineage, reverse_lineage, volume, experiment, parameters)
        if len(branch) == 0:
            break
        monitoring.to_log_and_console(str(proc) + "    ... check branch starting at " + str(min(branch)))

########################################################################################
#
#
#
########################################################################################


def postcorrection_process(experiment, parameters):

    proc = "postcorrection_process"
    #
    # parameter type checking
    #

    if not isinstance(experiment, common.Experiment):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'experiment' variable: "
                                      + str(type(experiment)))
        sys.exit(1)

    if not isinstance(parameters, PostCorrectionParameters):
        monitoring.to_log_and_console(str(proc) + ": unexpected type for 'parameters' variable: "
                                      + str(type(parameters)))
        sys.exit(1)

    #
    # read the lineage file, if any
    #
    segmentation_dir = experiment.astec_dir.get_directory()
    lineage_tree_file = common.find_file(segmentation_dir, experiment.astec_dir.get_file_name("_lineage"),
                                         file_type='lineage', callfrom=proc, verbose=False)

    if lineage_tree_file is None:
        monitoring.to_log_and_console(str(proc) + ": unable to find lineage file in " + str(segmentation_dir))
        sys.exit(1)
    elif os.path.isfile(os.path.join(segmentation_dir, lineage_tree_file)):
        lineage_tree_path = os.path.join(segmentation_dir, lineage_tree_file)
        lineage_tree_information = properties.read_dictionary(lineage_tree_path)
    else:
        monitoring.to_log_and_console(str(proc) + ": " + str(lineage_tree_file) + " is not a file?")
        sys.exit(1)

    #
    # test lineage
    #
    monitoring.to_log_and_console("    .. test lineage", 1)
    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
    properties.check_volume_lineage(lineage_tree_information, time_digits_for_cell_id=time_digits_for_cell_id)

    sys.exit(1)
    #
    # get lineage and volume dictionaries
    #
    dict_lineage = {}
    dict_volume = {}
    for key in lineage_tree_information.keys():
        if key == properties.keydictionary['lineage']['output_key']:
            dict_lineage = lineage_tree_information[properties.keydictionary['lineage']['output_key']]
    for key in lineage_tree_information.keys():
        if key == properties.keydictionary['volume']['output_key']:
            dict_volume = lineage_tree_information[properties.keydictionary['volume']['output_key']]

    if dict_lineage == {}:
        monitoring.to_log_and_console(str(proc) + ": empty lineage information in " + str(lineage_tree_file))
        sys.exit(1)
    if dict_volume == {}:
        monitoring.to_log_and_console(str(proc) + ": empty volume information in " + str(lineage_tree_file))
        sys.exit(1)
    #
    #
    #
    _remove_small_branches(dict_lineage, dict_volume, experiment, parameters)