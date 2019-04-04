from ASTEC.CommunFunctions.ImageHandling import SpatialImage, imread, imsave
import numpy as np
import argparse
from scipy import ndimage as nd

def main():
    parser = argparse.ArgumentParser(description='Extract KTR from intensity and segmentation.')
    parser.add_argument('-i', '--intensity', help='input KTR intensity file', required=True)
    parser.add_argument('-s', '--segmentation', help='input segmented file', required=True)
    parser.add_argument('-o', '--output', help='output csv file', required=True)
    
    args = parser.parse_args()

    #read images
    intensity = imread(args.intensity)
    seg = imread(args.segmentation)

    # compute bounding boxes to speed up the process
    bboxes = nd.find_objects(seg)

    # list the cell ids 
    c_ids = np.unique(seg)

    # build the avarage for each cell id
    avgs = {}
    for c in c_ids:
        # we assume that 1 is the background
        if c != 1:
            bb = bboxes[c-1]
            avgs[c] = np.mean(intensity[bb][seg[bb] == c])

    # write the averages in the output file
    with open(args.output, 'w') as f:
        f.write('c_id, avg\n')
        for c, avg in avgs.iteritems():
            f.write('{c:d}, {avg:.4f}\n'.format(c=c, avg=float(avg)))
        f.close()

if __name__ == '__main__':
    main()
