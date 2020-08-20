
import os
import imp
import sys
import copy
import time
import operator
import shutil
import numpy as np

from scipy.stats.stats import pearsonr

import common
import properties as properties
import CommunFunctions.cpp_wrapping as cpp_wrapping
from CommunFunctions.ImageHandling import imread, imsave, SpatialImage


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

        self.test_early_division = False
        self.volume_minimal_value = 2000
        self.lifespan_minimal_value = 25
        self.pearson_threshold = 0.9

        # diagnosis
        self.test_lineage = False

    ############################################################
    #
    # print / write
    #
    ############################################################

    def print_parameters(self):
        print("")
        print('PostCorrectionParameters')

        print('- soon = ' + str(self.test_early_division))
        print('- volume_minimal_value = ' + str(self.volume_minimal_value))
        print('- lifespan_minimal_value = ' + str(self.lifespan_minimal_value))
        print('- pearson_threshold = ' + str(self.pearson_threshold))

        print('- test_lineage = ' + str(self.test_lineage))

        print("")

    def write_parameters(self, log_file_name):
        with open(log_file_name, 'a') as logfile:
            logfile.write("\n")
            logfile.write('PostCorrectionParameters\n')

            logfile.write('- test_early_division = ' + str(self.test_early_division) + '\n')
            logfile.write('- volume_minimal_value = ' + str(self.volume_minimal_value) + '\n')
            logfile.write('- lifespan_minimal_value = ' + str(self.lifespan_minimal_value) + '\n')
            logfile.write('- pearson_threshold = ' + str(self.pearson_threshold) + '\n')

            logfile.write('- test_lineage = ' + str(self.test_lineage) + '\n')

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

        if hasattr(parameters, 'test_early_division'):
            if parameters.test_early_division is not None:
                self.test_early_division = parameters.test_early_division
        if hasattr(parameters, 'postcorrection_test_early_division'):
            if parameters.postcorrection_test_early_division is not None:
                self.test_early_division = parameters.postcorrection_test_early_division
        if hasattr(parameters, 'postcorrection_soon'):
            if parameters.postcorrection_soon is not None:
                self.test_early_division = parameters.postcorrection_soon
        if hasattr(parameters, 'postcor_Soon'):
            if parameters.postcor_Soon is not None:
                self.test_early_division = parameters.postcor_Soon

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

        if hasattr(parameters, 'test_lineage'):
            if parameters.test_lineage is not None:
                self.test_lineage = parameters.test_lineage
        if hasattr(parameters, 'postcorrection_test_lineage'):
            if parameters.postcorrection_test_lineage is not None:
                self.test_lineage = parameters.postcorrection_test_lineage
        if hasattr(parameters, 'postcor_test_lineage'):
            if parameters.postcor_test_lineage is not None:
                self.test_lineage = parameters.postcor_test_lineage


########################################################################################
#
# fusion procedure
#
########################################################################################

def _map_cell_fusion(astec_image, post_image, current_time, cells_to_be_fused, time_digits_for_cell_id=4):
    #
    #
    #
    labels_to_be_fused = cells_to_be_fused.get(current_time, '')
    #
    # no fusion to be done, just copy the file
    #
    if labels_to_be_fused == '':
        monitoring.to_log_and_console('    no cell fusion to be done for time #' + str(current_time), 2)
        astec_ext = common.get_image_extension(astec_image)
        post_ext = common.get_image_extension(post_image)
        if astec_ext == post_ext:
            shutil.copy2(astec_image, post_image)
        else:
            im = imread(astec_image)
            imsave(post_image, im)
        return

    #
    # build a mapping
    #
    im = imread(astec_image)
    mapping = range(np.max(im) + 1)
    for new, old in labels_to_be_fused:
        newlabel = new % 10 ** time_digits_for_cell_id
        oldlabel = old % 10 ** time_digits_for_cell_id
        mapping = [newlabel if i == oldlabel else i for i in mapping]
    #
    # after list comprehensions, mapping is of type list
    # needs to be casted into np.array to be applied on image
    #
    mapping = np.array(mapping, dtype=np.uint16)

    im = SpatialImage(mapping[im], voxelsize=im.voxelsize)
    imsave(post_image, im)
    return


########################################################################################
#
# some internal procedures
#
########################################################################################


def _test_early_division(direct_lineage, reverse_lineage, cell, lifespan_minimal_value, first_time_point, last_time_point, time_digits_for_cell_id):
    proc = "_test_early_division"

    #
    # get all the progeny of the mother cell of 'cell'
    #
    mother_cell = reverse_lineage[cell]
    sisters = copy.deepcopy(direct_lineage[mother_cell])
    if len(sisters) > 2:
        monitoring.to_log_and_console(str(proc) + ": cell #" + str(reverse_lineage[cell]) +
                                      "divides in more than 2 branches, skip it")
        return False
    sisters.remove(cell)
    #
    # extract the sister branch from the sister cell if 'cell'
    # until a leaf or a bifurcation
    #
    c = sisters[0]
    sister_branch = [c]
    while len(direct_lineage.get(c, [])) == 1:
        c = direct_lineage[c][0]
        sister_branch.append(c)

    #
    # the sister branch is too short (and ends before the end of the sequence)
    #
    if len(sister_branch) < lifespan_minimal_value \
            and sister_branch[-1]/10**time_digits_for_cell_id != last_time_point:
        return True
    #
    # extract the mother branch in the reverse order until the first bifurcation
    #
    mother_branch = [mother_cell]
    c = mother_cell
    while c in reverse_lineage.keys():
        c = reverse_lineage[c]
        #
        # TODO: test to be added
        #
        #if len(direct_lineage.get(c, [])) > 1:
        #    break
        mother_branch.append(c)

    #
    # the mother branch is too short (and begins after the beginning of the sequence)
    #
    if len(mother_branch) < lifespan_minimal_value \
            and mother_branch[-1]/10**time_digits_for_cell_id != first_time_point:
        return True

    return False


########################################################################################
#
#
#
########################################################################################


def _get_leaves(direct_lineage, reverse_lineage, volume, experiment, parameters):
    """
    Return a zip of tuple (leaf, branch_length, division_cell, division_time)
    each leaf defines a simple branch (issued from the division_cell) that is candidate
    for deletion, ie
    - either it finishes before the end of the sequence
    - or the volume of the leaf cell is too small
    :param direct_lineage:
    :param reverse_lineage:
    :param volume:
    :param experiment:
    :param parameters:
    :return:
    """

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
    #
    # leaves are nodes without lineage
    # get the largest time point value (last leaves)
    #
    nodes = list(set(direct_lineage.keys()).union(set([v for values in direct_lineage.values() for v in values])))
    leaves = set(nodes) - set(direct_lineage.keys())
    last_time = max(leaves) / 10 ** time_digits_for_cell_id
    first_time = min(set(direct_lineage.keys())) / 10 ** time_digits_for_cell_id

    #
    # a branch (ended) by a leave is candidate for deletion if
    # - it ends before the last time point (so before the end of the series), or
    # - the volume of the leave cell is too small
    #
    candidates_for_deletion = [leaf for leaf in leaves if (((leaf / 10 ** time_digits_for_cell_id) < last_time)
                                                             or (volume[leaf] < parameters.volume_minimal_value))]

    lengths = []
    division_cells = []
    division_times = []
    for leaf in candidates_for_deletion:
        #
        # get a branch from leaf 'leaf'
        # get back from the leaf until its parents has more than 1 daughter
        # (so the parent of two branches does not belong to the extracted branch)
        # reverse it
        # the branch begins then with the first point after the bifurcation
        #
        cell = leaf
        branch = [cell]
        while len(direct_lineage.get(reverse_lineage.get(cell, ''), '')) == 1:
            cell = reverse_lineage[cell]
            branch.append(cell)
        lengths.append(len(branch))
        division_cell = reverse_lineage.get(branch[-1], '')
        division_cells.append(division_cell)
        if division_cell is not '':
            division_times.append(division_cell / 10 ** time_digits_for_cell_id)
        else:
            branch.reverse()
            if len(branch) <= 4:
                monitoring.to_log_and_console("        ... found an orphan branch: " + str(branch), 3)
            else:
                monitoring.to_log_and_console("        ... found an orphan branch: " + str(branch[0:2]) + "..."
                                              + str(branch[-2:]), 3)
            division_times.append(first_time - 1)

    #
    # branches is an array of tuples (#leaf_id, branch_length, #division_id)
    #
    branches = zip(candidates_for_deletion, lengths, division_cells, division_times)
    return branches


def _fuse_branch(lineage, volume, surfaces, labels_to_be_fused, division_cell, branch, experiment):

    proc = "_fuse_branch"
    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()

    #
    # check division cell progeny
    #
    progeny = copy.deepcopy(lineage[division_cell])

    if len(progeny) > 2:
        monitoring.to_log_and_console(str(proc) + ": cell #" + str(division_cell) +
                                      " divides in more than 2 branches, skip it")
        return
    elif len(progeny) == 1:
        monitoring.to_log_and_console(str(proc) + ": cell #" + str(division_cell)
                                      + " has only one daughter?! this is weird")
        return

    progeny.remove(branch[0])

    #
    # follow the progeny
    #
    ibranch = 0
    iprogeny = 0
    sister_branch = []

    while True:
        cell_time = branch[ibranch] / 10 ** time_digits_for_cell_id
        labels_to_be_fused.setdefault(cell_time, []).append((progeny[iprogeny], branch[ibranch]))
        sister_branch.append(progeny[iprogeny])

        #
        # update volumes
        #
        volume[progeny[iprogeny]] += volume[branch[ibranch]]
        volume.pop(branch[ibranch], None)
        #
        # update surfaces
        # 1. a. update surfaces for common neighbors
        #    b. add new neighbors
        # 2. suppress branch[ibranch] from progeny[iprogeny] neighbors
        # 3. update neighbors's contact surface
        #
        for cell in surfaces[branch[ibranch]]:
            if cell in surfaces[progeny[iprogeny]]:
                surfaces[progeny[iprogeny]][cell] += surfaces[branch[ibranch]][cell]
            else:
                surfaces[progeny[iprogeny]][cell] = surfaces[branch[ibranch]][cell]
        if branch[ibranch] in surfaces[progeny[iprogeny]].keys():
            del(surfaces[progeny[iprogeny]][branch[ibranch]])
        for cell in surfaces[progeny[iprogeny]].keys():
            if cell in surfaces.keys():
                surfaces[cell][progeny[iprogeny]] = surfaces[progeny[iprogeny]][cell]
                if branch[ibranch] in surfaces[cell].keys():
                    del(surfaces[cell][branch[ibranch]])
        surfaces.pop(branch[ibranch], None)
        #
        # update lineage
        #
        if ibranch == 0:
            lineage[division_cell].remove(branch[ibranch])
            lineage.pop(branch[ibranch], None)
        else:
            lineage.pop(branch[ibranch], None)

        #
        # next cell, are we done?
        #
        ibranch += 1
        if ibranch >= len(branch):
            break

        progeny = lineage.get(progeny[iprogeny], '')
        if progeny == '':
            break

        iprogeny = -1
        if len(progeny) >= 2:
            for i in range(len(progeny)):
                if branch[ibranch] in surfaces[progeny[i]].keys():
                    if iprogeny == -1:
                        iprogeny = i
                    elif surfaces[progeny[i]][branch[ibranch]] > surfaces[progeny[iprogeny]][branch[ibranch]]:
                        iprogeny = i
            if iprogeny == -1:
                monitoring.to_log_and_console(str(proc) + ": cell #" + str(branch[ibranch]) + " has no neighbors in "
                                              + str(progeny) + " ?! this is weird")
                break
        else:
            iprogeny = 0

    #
    #
    #
    if len(branch) <= 4:
        strbranch = str(branch)
    else:
        strbranch = str(branch[0:2]) + "..." + str(branch[-2:])
    if len(sister_branch) <= 4:
        strsister = str(sister_branch)
    else:
        strsister = str(sister_branch[0:2]) + "..." + str(sister_branch[-2:])

    monitoring.to_log_and_console("            ... " + strbranch + " was fused with " + strsister, 3)

    return


def _prune_lineage_tree(lineage, volume, surfaces, experiment, parameters):
    proc = "_prune_lineage_tree"
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

    time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
    first_time = min(set(lineage.keys())) / 10 ** time_digits_for_cell_id
    last_time = max(set([v for values in lineage.values() for v in values])) / 10 ** time_digits_for_cell_id

    labels_to_be_fused = {}
    undeletable_leaves = set()

    iteration = 1
    previous_iteration_division_time = last_time
    while True:

        reverse_lineage = {v: k for k, values in lineage.iteritems() for v in values}

        #
        # get branches that end before the last time point or
        # with a leaf cell that has a small volume
        #
        branches = _get_leaves(lineage, reverse_lineage, volume, experiment, parameters)

        #
        # branches that occurs after iteration_division_time
        # are no more candidates for deletion
        #
        for leaf, branch_length, division_cell, division_time in branches:
            if division_time > previous_iteration_division_time:
                undeletable_leaves.add(leaf)

        candidates = [(le, b, dc, dt) for le, b, dc, dt in branches if le not in undeletable_leaves]
        candidates = sorted(candidates, key=operator.itemgetter(1))
        candidates = sorted(candidates, key=operator.itemgetter(3), reverse=True)

        monitoring.to_log_and_console("       iteration #" + str(iteration) + ", found " + str(len(candidates))
                                      + " candidates among " + str(len(branches)) + " end branches", 2)

        # candidates is an array of tuples (#leaf_id, branch_length, #division_cell, division_time)
        # sorted by decreasing division cell id (and then by increasing length)
        n_fusions = 0
        iteration_division_time = first_time
        for leaf, branch_length, division_cell, division_time in candidates:

            if leaf in undeletable_leaves:
                continue

            #
            # do nothing if
            # - mother_cell is the root (branch begins from the very beginning)
            # - mother_cell has less than 2 progeny (one branch has been fused with the other)
            #   => end branches have to be recomputed
            #
            # break (and recompute branches) when
            # - mother has no more lineage (belongs to a fused branch)
            #

            if division_time < iteration_division_time:
                break
            if division_cell == '':
                monitoring.to_log_and_console("         branch ending at " + str(leaf) + " begins at root. Skip it", 4)
                undeletable_leaves.add(leaf)
                continue
            if len(lineage[division_cell]) < 2:
                # this is probably the sister branch of a branch that has been fused with
                # monitoring.to_log_and_console("            mother cell ("+ str(division_cell) + ") of branch ending at "
                #                              + str(leaf) + " has " + str(len(lineage[division_cell]))
                #                              + " progeny. Skip it", 4)
                continue

            if branch_length < parameters.lifespan_minimal_value:
                pass
            else:
                continue

            # we found a branch to be removed
            if division_time > iteration_division_time:
                iteration_division_time = division_time

            #
            # get the whole branch from the leaf
            #
            le = leaf
            branch = [le]
            while reverse_lineage.get(le, '') != division_cell:
                le = reverse_lineage[le]
                branch.append(le)
            #
            # reverse it, so it begins from the first point after the division cell
            #
            branch.reverse()

            #
            # fuse the branch
            # lineage, volume and surfaces are passed by their ids, so the modifications done
            # by _fuse_branch() are kept
            #
            _fuse_branch(lineage, volume, surfaces, labels_to_be_fused, division_cell, branch, experiment)
            n_fusions += 1

        #
        # end of an iteration
        # no deletion? end of the loop
        #
        if n_fusions == 0:
            break

        monitoring.to_log_and_console("         " + str(n_fusions) + " branches have been fused at division time "
                                      + str(iteration_division_time), 1)
        previous_iteration_division_time = iteration_division_time
        iteration += 1

    return lineage, volume, surfaces, labels_to_be_fused


########################################################################################
#
#
#
########################################################################################

def contact_surface_computation(experiment):
    """

    :param experiment: commonTools.Experiment
    :return:
    """

    proc = 'contact_surface_computation'

    #
    # get directory name where to find co-registered images of the sequence
    # as well as the common image suffix
    #
    astec_dir = experiment.astec_dir.get_directory()

    #
    # is there a post-segmentation directory in the intra-registration directory ?
    #
    if not os.path.isdir(astec_dir):
        monitoring.to_log(proc + ": '" + str(astec_dir) + "' does not exist")
        return None

    monitoring.to_log_and_console("       will compute contact surface properties from '" + str(astec_dir) + "'", 0)

    #
    # build name format for (post-corrected) segmentation images
    #
    name_format = experiment.astec_dir.get_file_prefix() + experiment.astec_dir.get_file_suffix() \
                  + experiment.astec_dir.get_time_prefix() + experiment.get_time_format()

    suffix = common.get_file_suffix(experiment, astec_dir, name_format, flag_time=experiment.get_time_format())
    if suffix is None:
        monitoring.to_log_and_console(proc + ": no consistent naming was found in '"
                                      + str(astec_dir) + "'", 1)
        monitoring.to_log_and_console("Exiting.", 0)
        sys.exit(1)
    name_format += "." + str(suffix)
    template_format = os.path.join(astec_dir, name_format)

    #
    #
    #
    output_name = experiment.astec_dir.get_file_prefix() + experiment.astec_dir.get_file_suffix() + "_surfaces"
    output_name = os.path.join(astec_dir, output_name)

    if os.path.isfile(output_name + ".xml"):
        if not monitoring.forceResultsToBeBuilt:
            monitoring.to_log_and_console('       xml file already existing', 2)
            return output_name + ".xml"
        else:
            monitoring.to_log_and_console('       xml file already existing, but forced', 2)

    first_time_point = experiment.first_time_point + experiment.delay_time_point
    last_time_point = experiment.last_time_point + experiment.delay_time_point

    options = "-feature contact -surface-estimation 6-neighbors"
    cpp_wrapping.cell_properties(template_format, output_name + ".xml", first_time_point, last_time_point,
                                 diagnosis_file=output_name + ".txt", other_options=options, monitoring=monitoring)

    return output_name + ".xml"


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
    if parameters.test_lineage:
        monitoring.to_log_and_console("   ... test lineage (before)", 1)
        time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
        properties.check_volume_lineage(lineage_tree_information, time_digits_for_cell_id=time_digits_for_cell_id)

    #
    # get lineage and volume dictionaries
    #
    dict_lineage = properties.get_dictionary_entry(lineage_tree_information, 'lineage')
    if dict_lineage == {}:
        monitoring.to_log_and_console(str(proc) + ": empty lineage information in " + str(lineage_tree_file))
        sys.exit(1)

    dict_volume = properties.get_dictionary_entry(lineage_tree_information, 'volume')
    if dict_volume == {}:
        monitoring.to_log_and_console(str(proc) + ": empty volume information in " + str(lineage_tree_file))
        sys.exit(1)

    #
    #
    #
    monitoring.to_log_and_console("   ... compute contact surfaces", 1)
    surface_file = contact_surface_computation(experiment)
    surface_information = properties.read_dictionary(surface_file)
    dict_surface = properties.get_dictionary_entry(surface_information, 'contact')
    if dict_surface == {}:
        monitoring.to_log_and_console(str(proc) + ": empty surface information in " + str(surface_file))
        sys.exit(1)

    #
    #
    #
    monitoring.to_log_and_console("   ... lineage pruning", 1)
    lineage, volume, surfaces, cells_to_be_fused = _prune_lineage_tree(dict_lineage, dict_volume, dict_surface,
                                                                       experiment, parameters)
    # sys.exit(1)
    # _remove_small_branches(dict_lineage, dict_volume, experiment, parameters)

    #
    # save lineage tree
    #
    new_lineage_tree = {properties.keydictionary['lineage']['output_key']: lineage,
                        properties.keydictionary['volume']['output_key']: volume}
    segmentation_dir = experiment.post_dir.get_directory()
    lineage_tree_path = os.path.join(segmentation_dir, experiment.post_dir.get_file_name("_lineage") + "."
                                     + experiment.result_lineage_suffix)
    properties.write_dictionary(lineage_tree_path, new_lineage_tree)

    #
    # test lineage
    #
    if parameters.test_lineage:
        monitoring.to_log_and_console("   ... test lineage (after)", 1)
        time_digits_for_cell_id = experiment.get_time_digits_for_cell_id()
        properties.check_volume_lineage(new_lineage_tree, time_digits_for_cell_id=time_digits_for_cell_id)

    sys.exit(1)

    #
    # apply cell fusion
    #
    first_time = experiment.first_time_point + experiment.delay_time_point
    last_time = experiment.last_time_point + experiment.delay_time_point
    for current_time in range(first_time, last_time + 1, experiment.delta_time_point):

        acquisition_time = experiment.get_time_index(current_time)
        monitoring.to_log_and_console('... cell fusion of time #' + acquisition_time, 1)

        #
        # nothing to do if the post-segmentation image exists
        #
        post_dir = experiment.post_dir.get_directory()
        post_name = experiment.post_dir.get_image_name(current_time)
        post_image = common.find_file(post_dir, post_name, file_type='image', callfrom=proc, local_monitoring=None,
                                      verbose=False)

        if post_image is not None:
            if monitoring.forceResultsToBeBuilt is False:
                monitoring.to_log_and_console('    post-segmentation image already existing', 2)
                continue
            else:
                monitoring.to_log_and_console('    post-segmentation image already existing, but forced', 2)
        post_image = os.path.join(post_dir, post_name + "." + experiment.result_image_suffix)

        astec_dir = experiment.astec_dir.get_directory()
        astec_name = experiment.astec_dir.get_image_name(current_time)
        astec_image = common.find_file(astec_dir, astec_name, file_type='image', callfrom=proc, local_monitoring=None,
                                       verbose=False)
        if astec_image is None:
            monitoring.to_log_and_console("    .. " + proc + ": no segmentation image was found for time "
                                          + str(current_time), 2)
            return
        astec_image = os.path.join(astec_dir, astec_image)
        #
        #
        #

        # print(str(cells_to_be_fused))
        # print("")
        # print("")
        _map_cell_fusion(astec_image, post_image, current_time, cells_to_be_fused,
                         time_digits_for_cell_id=time_digits_for_cell_id)

    return
