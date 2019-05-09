function [t1_file,seg_file,fmri_file,origfmri_file] = setup( ...
	out_path, ...
	t1_file, ...
	seg_file, ...
	fmri_file ...
	)

% Copy files to output directory, unzip, and give consistent names
[~,n,e] = fileparts(t1_file);
copyfile(t1_file,out_path);
t1_file = fullfile(out_path,[n e]);
if strcmp('.gz',t1_file(end-2:end))
	system(['gunzip -f ' t1_file]);
	t1_file = t1_file(1:end-3);
end
movefile(t1_file,fullfile(out_path,'t1.nii'));
t1_file = fullfile(out_path,'t1.nii');

[~,n,e] = fileparts(seg_file);
copyfile(seg_file,out_path);
seg_file = fullfile(out_path,[n e]);
if strcmp('.gz',seg_file(end-2:end))
	system(['gunzip -f ' seg_file]);
	seg_file = seg_file(1:end-3);
end
movefile(seg_file,fullfile(out_path,'seg.nii'));
seg_file = fullfile(out_path,'seg.nii');

[~,n,e] = fileparts(fmri_file);
copyfile(fmri_file,out_path);
fmri_file = fullfile(out_path,[n e]);
if strcmp('.gz',fmri_file(end-2:end))
	system(['gunzip -f ' fmri_file]);
	fmri_file = fmri_file(1:end-3);
end
movefile(fmri_file,fullfile(out_path,'fmri.nii'));
fmri_file = fullfile(out_path,'fmri.nii');

% Make an extra copy of the fmri so we can get the original headers even
% after SPM realignment has been at work
copyfile(fmri_file,fullfile(out_path,'origfmri.nii'));
origfmri_file = fullfile(out_path,'origfmri.nii');

