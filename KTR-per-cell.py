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

    intensity = imread(args.intensity)
    seg = imread(args.segmentation)

    bboxes = nd.find_objects(seg)
    c_ids = np.unique(seg)

    avgs = {}
    for c in c_ids:
        if c != 1:
            bb = bboxes[c-1]
            avgs[c] = np.mean(intensity[bb][intensity[bb] == c])

    with open(args.output, 'w') as f:
        f.write('c_id, avg\n')
        for c, avg in avgs.iteritems():
            f.write('{c:d}, {avg:.4f}\n'.format(c=c, avg=float(avg)))
        f.close()

if __name__ == '__main__':
    main()
