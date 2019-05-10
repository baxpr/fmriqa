
singularity \
run \
--cleanenv \
--bind INPUTS:/INPUTS \
--bind OUTPUTS:/OUTPUTS \
baxpr-fmriqa-master-v4.2.0.simg \
out_path /OUTPUTS \
t1_file INPUTS/t1.nii.gz \
seg_file INPUTS/seg.nii.gz \
fmri_file INPUTS/fmri.nii.gz \
project UNK_PROJ \
subject UNK_SUBJ \
session UNK_SESS \
scan UNK_SCAN
