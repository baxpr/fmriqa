# Functional MRI QA pipeline

## Building

1. Test the matlab code before compiling: `src/testmatlab.m`
1. Compile: `compile_matlab.sh`
1. Test the compiled runtime: `bin/test_compiled_matlab.sh`
1. Build the Singularity container: `Singularity.v4.2.0`, <https://www.singularity-hub.org/collections/2945>

## Usage

1. See `test_sing_container.sh`.

## Inputs

The inputs must all be provided, in the correct order. Paths are with respect to the container root.

1. Name of the output directory 
1. Filename of the T1 structural image (.nii.gz)
1. Filename of the segmented T1 image (.nii.gz), typically the SEG output of a MultiAtlas or SLANT pipeline
1. Filename of the 4D fMRI (.nii.gz)
1. XNAT project label
1. XNAT subject label
1. XNAT session label
1. XNAT scan label (of the fMRI)

## Processing

1. Motion realignment and creation of mean fMRI
1. Coregister T1 to mean fMRI
1. Compute SNR and quality metrics
1. Carpet plots, graphical report

## Outputs

```
fmriqa.pdf                               PDF report
rp_fmri.txt                              Realignment parameters (SPM12 style)
fmriqa_stats.csv                         Summary stats
fmriqa_stats_wide.csv                    Summary stats in wide format (XNAT/REDCap compatible)
FD.txt                                   Framewise displacement time series
DVARS.txt                                DVARS time series
global.txt                               Global mean time series
meanfmri.nii.gz                          Mean fMRI image after realignment
median_voxel_displacement_mm.txt         Framewise displacement, median over voxels
temporal_snr.nii.gz                      Temporal signal-to-noise ratio image
voxel_displacement_mm_95prctile.nii.gz   Framewise displacement image (95th percentile over time)
```