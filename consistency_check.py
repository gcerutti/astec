import argparse
import cPickle as pkl
import numpy as np

def main(embryo_name, start_time, output):
    with open(embryo_name) as f:
        datas = pkl.load(f)
    lin_tree = datas['cell_lineage']
    fates3 = datas["cell_fate_3"]
    names = datas['cell_name']
    vol = datas['cell_volume']
    surf_ex = datas['cell_contact_surface']

    cell_per_time = {}
    for c, d in lin_tree.iteritems():
        cell_per_time.setdefault(c//10**4, []).append(c)
        for ci in d:
            cell_per_time.setdefault(ci//10**4, []).append(ci)
    cell_per_time = {k: len(set(v)) for k, v in cell_per_time.iteritems()}
    first_time = min([k for k, v in cell_per_time.iteritems() if v == start_time])

    surfaces = {c:np.sum(v.values()) for c, v in surf_ex.iteritems()}
    inv_lin_tree = {vi:k for k, v in lin_tree.iteritems() for vi in v}
    # to_remove = []
    to_remove = [c for c in lin_tree if c//10**4 < first_time]
    no_name = []
    same_name = []
    no_vol = []
    for c in lin_tree:
        if not c in names:
            to_remove += [c]
            if not c in datas.get('unknown_key', []):
                no_name += [c]
        if not c in vol:
            to_remove += [c]
            if not c in datas.get('unknown_key', []):
                no_vol += [c]
        if len(lin_tree[c]) == 2:
            if names.get(lin_tree[c][0], 'a') == names.get(lin_tree[c][1], 'b'):
                tmp = []
                same_name += lin_tree[c]
                to_treat = list(lin_tree[c])
                while to_treat != []:
                    tmp += [to_treat.pop()]
                    to_treat += lin_tree.get(tmp[-1], [])
                # print tmp
                to_remove += tmp

    to_remove += datas.get('unknown_key', [])
    to_remove = set(to_remove)

    for c in to_remove:
        surf_ex.pop(c, -1)
        vol.pop(c, -1)
        names.pop(c, -1)
        fates3.pop(c, -1)
        surf_ex.pop(c, -1)
        surfaces.pop(c, -1)
        if c in inv_lin_tree and inv_lin_tree[c] in lin_tree:
            tmp = lin_tree[inv_lin_tree[c]]
            if len(tmp) == 1:
                lin_tree.pop(inv_lin_tree[c])
            else:
                tmp.remove(c)
        for v in surf_ex.itervalues():
            v.pop(c, -1)
        lin_tree.pop(c, -1)

    nodes = set(lin_tree).union(inv_lin_tree)
    to_remove_surf = []
    for c in surf_ex:
        if not c in nodes:
            # tmp += [c]
            to_remove_surf += [c]

    for c in to_remove_surf:
        surf_ex.pop(c)

    been_removed = []
    for c, v in surf_ex.iteritems():
        for ci in to_remove_surf + list(to_remove):
            if ci in v:
                v.pop(ci)
                been_removed += [ci]

    been_removed = set(been_removed)

    to_remove_surf = set(to_remove_surf).difference(nodes)

    branches_to_correct=[k for k in lin_tree.keys() if len(lin_tree[k])>1 and len(list(set(lin_tree[k])))<2]
    for bc in branches_to_correct:
        dau=lin_tree[bc][0]
        lin_tree.update({bc:[dau]})

    with open(output, 'w') as f:
        f.write('Cells with no name')
        print 'Cells with no name'
        for c in no_name:
            f.write(',{:d}'.format(c))
            print '{:d}'.format(c),
        f.write('\nSister cells with same name')
        print '\nSister cells with same name'
        for c in same_name:
            f.write(',{:d}'.format(c))
            print '{:d}'.format(c),
        f.write('\nCells with no volume')
        print '\nCells with no volume'
        for c in no_vol:
            f.write(',{:d}'.format(c))
            print '{:d}'.format(c),
        f.write('\nCells that exists as neighbors but not as cells')
        print '\nCells that exists as neighbors but not as cells'
        for c in to_remove_surf:
            f.write(',{:d}'.format(c))
            print '{:d}'.format(c),
        f.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Output the potential errors in a manually corrected lineage tree')
    parser.add_argument('-i', '--input', help='path to the pkl file', required=True)
    parser.add_argument('-s', '--start-time', help='number of cells the lineage tree should start to be checked at (optional, default value is 64)',
                            required=False, default=64, type=int)
    parser.add_argument('-o', '--output', help='output csv file', required=True)
    
    args = parser.parse_args()
    
    main(args.input, args.start_time, args.output)