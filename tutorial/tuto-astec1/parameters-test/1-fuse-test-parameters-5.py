#
# left and right camera images are in the same directory
# there are two channels
#
# only two time points are reconstructed
#

PATH_EMBRYO = '.'

EN = '2019-Tutorial100'

DIR_RAWDATA = 'RAWDATA'
DIR_LEFTCAM_STACKZERO = 'stack_0_channel_0'
DIR_RIGHTCAM_STACKZERO = 'stack_0_channel_0'
DIR_LEFTCAM_STACKONE = 'stack_1_channel_0'
DIR_RIGHTCAM_STACKONE = 'stack_1_channel_0'

DIR_LEFTCAM_STACKZERO_CHANNEL2 = 'stack_0_channel_1'
DIR_RIGHTCAM_STACKZERO_CHANNEL2 = 'stack_0_channel_1'
DIR_LEFTCAM_STACKONE_CHANNEL2 = 'stack_1_channel_1'
DIR_RIGHTCAM_STACKONE_CHANNEL2 = 'stack_1_channel_1'

begin = 0
end = 1

acquisition_orientation = 'right'
acquisition_mirrors = False
acquisition_resolution = (1., 1., 1.)

target_resolution = 1.0

EXP_FUSE = 'TEST05'
