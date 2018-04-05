import os, sys
#sys.path.append('CommunFunctions') 
sys.path.append(os.path.join(os.path.dirname(__file__),"CommunFunctions"))
from ImageHandling import SpatialImage, imread, imsave
from scipy import ndimage as nd
import numpy as np
from cpp_wrapping import reech3d

def compute_volumes(im):
    labels = np.unique(im)
    volume = nd.sum(np.ones_like(im), im, index=np.int16(labels))
    return dict(zip(labels, volume))


def croping(image_input, impage_output,downsize, dilation_x_0=40, dilation_x_1=40, dilation_y_0=40,dilation_y_1=40, no_crop=False):
    ''' Automatically crop and resample an image
    image_input : path to the input image
    image_output : path to the output image
    downsize : voxel size of the resampled image (in \mu m)
    # Crop options:
      no_crop [default=False]: if True, then the resampled image is not cropped
    # Dilation parameters: 
      dilation_x_0 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'left' x direction 
      dilation_x_1 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'right' x direction
      dilation_y_0 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'top' y direction
      dilation_y_1 [default=40]: parameter for dilation of the bounding box computed for the cropping of the resampled image in 'bottom' y direction
    '''
    shape_begin=imread(image_input).shape #Image Size
    image_main = reech3d(image_input, (shape_begin[0]/np.float(downsize), shape_begin[1]/np.float(downsize), shape_begin[2])) #Downsampling
    if no_crop:
        out=SpatialImage(image_main)
    else:
        vxsize=image_main.resolution
        im_maxed = image_main.max(axis=2)
        thr = np.mean(im_maxed)
        im_th = np.zeros((im_maxed.shape[0], im_maxed.shape[1], 1), dtype=np.uint8)
        im_th[im_maxed > thr] = 1
        comp_co = nd.label(im_th)[0]
        volumes = compute_volumes(comp_co)
        volumes.pop(0)
        label = volumes.keys()[np.argmax(volumes.values())]
        bb = nd.find_objects(comp_co)[label-1]
        #change 40-->120 et 40-->80 nouvelle box
        # dilation parameter:
        # 
        bb2 = (slice(max(bb[0].start - dilation_x_0, 0), min(image_main.shape[0], bb[0].stop + dilation_x_1), None), 
               slice(max(bb[1].start - dilation_y_0, 0), min(image_main.shape[1], bb[1].stop + dilation_y_1), None), 
               slice(0, image_main.shape[2]))
        out=SpatialImage(image_main[bb2])
    out.voxelsize=(float("{0:.1f}".format(image_main.resolution[0])), float("{0:.1f}".format(image_main.resolution[1])), float("{0:.1f}".format(image_main.resolution[2])))
    imsave(impage_output, out.astype(np.uint16))
