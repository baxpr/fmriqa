#!/bin/bash

sh run_spm12.sh /usr/local/MATLAB/MATLAB_Runtime/v92 \
out_path /OUTPUTS \
t1_file ../INPUTS/t1.nii.gz \
seg_file ../INPUTS/seg.nii.gz \
fmri_file ../INPUTS/fmri.nii.gz \
project UNK_PROJ \
subject 600000000001 \
session 600000000001 \
scan 50362
