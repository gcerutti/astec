#
# identical to tutorial parameters
# but fusion is made in a hierarchical way
# 1. stacks are reconstructed independently
# 2. they are co-registered
# 3. a global reconstruction is performed
#
# only two time points are reconstructed
#

PATH_EMBRYO = '.'

EN = '2019-Tutorial100'

begin = 0
end = 1

acquisition_orientation = 'right'
acquisition_mirrors = False
acquisition_resolution = (1., 1., 1.)

target_resolution = 1.0

fusion_strategy = 'hierarchical-fusion'

EXP_FUSE = 'TEST06'
