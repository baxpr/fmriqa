function [fmri_file,rfmri_file,meanfmri_file,rp_file] = realignment(fmri_file)

[fmri_p,fmri_n,func_e] = fileparts(fmri_file);

% SPM job
matlabbatch = [];
tag = 1;
matlabbatch{tag}.spm.spatial.realign.estwrite.data = {{fmri_file}};
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{tag}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{tag}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{tag}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{tag}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{tag}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{tag}.spm.spatial.realign.estwrite.roptions.prefix = 'r';
spm_jobman('run',matlabbatch)

% Filename of realignment params
rp_file = fullfile(fmri_p,['rp_' fmri_n '.txt']);

% Filenames for realigned images
meanfmri_file = fullfile(fmri_p,['mean' fmri_n func_e]);
rfmri_file = fullfile(fmri_p,['r' fmri_n func_e]);

