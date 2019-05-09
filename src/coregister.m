function rcseg_file = coregister( ...
	out_path,t1_file,seg_file,meanfmri_file)


%% Apply seg file as mask to get brain-only image
Vt1 = spm_vol(t1_file);
Vseg = spm_vol(seg_file);
spm_check_orientations([Vt1; Vseg]);

Yt1 = spm_read_vols(Vt1);
Yseg = spm_read_vols(Vseg);

% Brain values are >=4
Yt1(Yseg(:)<4) = 0;

Vout = rmfield(Vt1,'pinfo');
Vout.fname = fullfile(out_path,'brain.nii');
spm_write_vol(Vout,Yt1);


%% Align brain center of mass to fmri

% Image centers of mass
fmri_com = find_center_of_mass(meanfmri_file);
t1_com = find_center_of_mass(t1_file);

% How to move the source image
source_shift = fmri_com - t1_com;
source_shift_matrix = spm_matrix(source_shift);

% Move via batch reorient
matlabbatch = [];
tag = 0;

tag = tag + 1;
matlabbatch{tag}.spm.util.reorient.srcfiles = {seg_file};
matlabbatch{tag}.spm.util.reorient.transform.transM = source_shift_matrix;
matlabbatch{tag}.spm.util.reorient.prefix = 'c';

tag = tag + 1;
matlabbatch{tag}.spm.util.reorient.srcfiles = {t1_file};
matlabbatch{tag}.spm.util.reorient.transform.transM = source_shift_matrix;
matlabbatch{tag}.spm.util.reorient.prefix = 'c';

spm_jobman('run',matlabbatch);

% Get filenames for reoriented images
[t1_p,t1_n,t1_e] = fileparts(t1_file);
[seg_p,seg_n,seg_e] = fileparts(seg_file);
ct1_file = fullfile(t1_p,['c' t1_n t1_e]);
cseg_file = fullfile(seg_p,['c' seg_n seg_e]);


%% Coregister and resample
matlabbatch = [];
tag = 0;

% Coregister
tag = tag + 1;
matlabbatch{tag}.spm.spatial.coreg.estimate.ref = {meanfmri_file};
matlabbatch{tag}.spm.spatial.coreg.estimate.source = {ct1_file};
matlabbatch{tag}.spm.spatial.coreg.estimate.other = {cseg_file};
matlabbatch{tag}.spm.spatial.coreg.estimate.eoptions = struct( ...
	'cost_fun', 'nmi', ...
	'sep', [4 2], ...
	'tol', [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001], ...
	'fwhm', [7 7] );

% Reslice
tag = tag + 1;
matlabbatch{tag}.spm.spatial.coreg.write.ref = {meanfmri_file};
matlabbatch{tag}.spm.spatial.coreg.write.source = {cseg_file};
matlabbatch{tag}.spm.spatial.coreg.write.roptions = struct( ...
	'interp', 0, ...
	'wrap', [0 0 0], ...
	'mask', 0, ...
	'prefix', 'r' );

spm_jobman('run',matlabbatch)

% Filename for coregistered image
rcseg_file = fullfile(seg_p,['rc' seg_n seg_e]);


