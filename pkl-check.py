#!/usr/bin/env/python

import argparse
import os
import numpy as np
import warnings
try:
    from TGMMlibraries import lineageTree
except Exception as e:
    warnings.warn('TGMMlibraries required and not loaded properly.')
    print('error raised: {:s}\n'.format(e))
    exit()

def get_daughters(n):
    if len(n) == 0:
        return n, n
    else:
        l = n[0]
        cycle = int(n.split('.')[0][1:])
        if not n.split('.')[1][-1].isdigit():
            end = n.split('.')[1][-1]
            id_ = int(n.split('.')[1][:-1])
        else:
            end = ''
            id_ = int(n.split('.')[1])
        return '%s%d.%04d%s'%(l, cycle+1, id_*2-1, end), '%s%d.%04d%s'%(l, cycle+1, id_*2, end)

def write_wrong_daughters(p, wrong_daugthers):
    with open(os.path.join(p, 'wrong_daugthers.csv'), 'w') as f:
        f.write('embryo name, id mother, name mother, in problematic cells list?\n'
                ', id daughter 1, name daughter 1, in problematic cells list?\n'
                ', id daughter 2, name daughter 2, in problematic cells list?\n\n\n')
        for em , vals in wrong_daugthers.iteritems():
            f.write(em + ',')
            for i, lame in enumerate(vals):
                f.write('{:d}, {:s}, {:s}\n'.format(lame[0][0][0], lame[1][0], 'yes' if lame[0][0][1] else 'no'))
                f.write(',{:d}, {:s}, {:s}\n'.format(lame[0][1][0], lame[1][1], 'yes' if lame[0][1][1] else 'no'))
                f.write(',{:d}, {:s}, {:s}\n,\n'.format(lame[0][2][0], lame[1][2], 'yes' if lame[0][2][1] else 'no'))
                if i < len(vals)-1:
                    f.write(',')
            f.write('\n')
        f.close()

def write_wrong_linkages(p, incons_within_cycle):
    with open(os.path.join(p, 'wrong_linkages.csv'), 'w') as f:
        f.write('embryo name,\n,cell id, [...]\n')
        f.write(',cell name, [...]\n')
        f.write(',in problematic cells list?, [...]\n')
        f.write(',is one of the sister somewhere else in the lineage tree?\n\n\n')
        for em, vals in incons_within_cycle.iteritems():
            if 0<len(vals):
                f.write(em + ',\n')
                for lame in vals:
                    # f.write(',')
                    tmp = np.array(lame[:-1])
                    f.write((',{:s}'*len(tmp)).format(*tmp[:, 0]))
                    f.write('\n')
                    f.write((',{:s}'*len(tmp)).format(*tmp[:, 1]))
                    f.write('\n')
                    f.write((',{:s}'*len(tmp)).format(*tmp[:, 2]))
                    f.write('\n')
                    f.write(','+lame[-1]+'\n,\n')
                f.write('\n')
        f.close()

def find_inconsistancies(lT, em):
    incons_within_cycle = {}
    wrong_daugthers = {}
    incons_within_cycle[em] = []
    wrong_daugthers[em] = []
    first_cells = list(lT.time_nodes[lT.t_b])
    to_treat = first_cells
    full_name_set = set(lT.name.values())
    while to_treat!=[]:
        c = to_treat.pop()
        cycle = lT.get_cycle(c)
        names_in_cycle = [lT.name[ci] for ci in cycle]
        name_set = set(names_in_cycle)
        if not len(name_set) == 1:
            if names_in_cycle[-1] != '':
                incons_within_cycle[em] += [[(lT.lT2pkl[ci], lT.name[ci],ci in lT.prob_cells) for ci in cycle]]
                if 2<=len(name_set):
                    for ni in sorted(name_set)[:-1]:
                        d1, d2 = get_daughters(ni)
                        if not d1 in name_set and d1 in full_name_set:
                            incons_within_cycle[em][-1] += [d1 + ' found']
                        elif not d2 in name_set and d2 in full_name_set:
                            incons_within_cycle[em][-1] += [d2 + ' found']
                        else:
                            incons_within_cycle[em][-1] += ['sister NOT found']
        ds = lT.successor.get(cycle[-1], [])
        if 2==len(ds) and not tuple(sorted([lT.name[ds[0]], lT.name[ds[1]]])) == get_daughters(lT.name[cycle[-1]]):
            if not(lT.name[ds[1]] == lT.name[ds[0]] == ''):
                prob_cells = [lT.lT2pkl[cycle[-1]]]+[lT.lT2pkl[di] for di in ds]
                prob_cells = [(ci, (lT.pkl2lT[ci] in lT.prob_cells)) for ci in prob_cells]
                wrong_daugthers[em] += [(prob_cells, [lT.name[lT.pkl2lT[ci[0]]] for ci in prob_cells])]
        to_treat += ds
    return incons_within_cycle, wrong_daugthers


def main():
    parser = argparse.ArgumentParser(description='Read a pkl ASTEC file and check for inconsistencies.' +
                                                 'To work, the lineage tree library is necessary')
    parser.add_argument('-i', '--input', help='input ASTEC .pkl file', required=True)
    parser.add_argument('-o', '--output', help='path to the output folder that will contain the csv files)', 
                            required=False, default = '.')
    
    args = parser.parse_args()
    
    if os.path.exists(args.output) and not os.path.isdir(args.output):
        print('Output has to be a folder ({:s})'.format(args.output))
    if not os.path.exists(args.output):
        print('Creating the ouput folder: {:s}'.format(args.output))
        os.makedirs(args.output)
    
    em = os.path.basename(args.input).split('.')[0]
    lT = lineageTree(args.input, ASTEC=True)

    incons_within_cycle, wrong_daugthers = find_inconsistancies(lT, em)
    write_wrong_daughters(args.output, wrong_daugthers)
    write_wrong_linkages(args.output, incons_within_cycle)


if __name__ == '__main__':
    main()