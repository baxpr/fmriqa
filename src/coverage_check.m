function coverage_check( ...
	out_path, ...
	seg_file, ...
	meanfmri_file, ...
	rfmri_file ...
	)

% First coreg reslice the mean fmri and the realigned fmri to the
% coregistered seg file. Use first volume of realigned fmri to save memory
clear matlabbatch
matlabbatch{1}.spm.spatial.coreg.write.ref = {seg_file};
matlabbatch{1}.spm.spatial.coreg.write.source = {meanfmri_file};
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'r';

matlabbatch{2}.spm.spatial.coreg.write.ref = {seg_file};
matlabbatch{2}.spm.spatial.coreg.write.source = {[rfmri_file ',1']};
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'r';

spm_jobman('run',matlabbatch)

[p,n,e] = fileparts(meanfmri_file);
rmeanfmri_file = fullfile(p,['r' n e]);

[p,n,e] = fileparts(rfmri_file);
rrfmri_file = fullfile(p,['r' n e]);

% Load the resliced images
Vseg = spm_vol(seg_file);
[Yseg,XYZ] = spm_read_vols(Vseg);
Yrmean = spm_read_vols(spm_vol(rmeanfmri_file));
Yrrfmri = spm_read_vols(spm_vol(rrfmri_file));

% Find the gray matter
Ygm = zeros(size(Yseg));
Ygm(:) = (Yseg(:)>0) & ~ismember(Yseg(:), ...
	[4 11 40 41 44 45 49 50 51 52]);
Vgm = rmfield(Vseg,'pinfo');
Vgm.fname = fullfile(out_path,'graytemp.nii');
spm_write_vol(Vgm,Ygm);

% Then, any values that are 0 or nan in the resliced fmri will be marked
Ylost1 = Ygm & (Yrmean==0);
Ylost2 = ~Ylost1 & Ygm & (Yrrfmri==0);

V1 = rmfield(Vseg,'pinfo');
V1.fname = fullfile(out_path,'lost_coverage.nii');
spm_write_vol(V1,Ylost1);

V2 = rmfield(Vseg,'pinfo');
V2.fname = fullfile(out_path,'lost_realign.nii');
spm_write_vol(V2,Ylost2);

% Show the gray matter image
spm_check_registration(Vgm.fname);

% Overlay the lost voxels images
spm_orthviews('Xhairs','off');
spm_orthviews('addcolouredimage',1,V1.fname,[0 1 0]);
spm_orthviews('addcolouredimage',1,V2.fname,[1 0 0]);
title({'Green = lack of coverage','Red = lost due to movement'})

% Move to center of mass (of gray matter image). We'll move over a bit in x
% to get away from the exact centerline
com = (XYZ * Ygm(:) / sum(Ygm(:)));
spm_orthviews('Reposition',com+[8 0 0]');


% Capture the graphical output
set(spm_figure('FindWin','Graphics'), ...
	'PaperOrientation','portrait', ...
	'PaperUnits','inches', ...
	'PaperSize',[8.5 11], ...
	'PaperType','usletter', ...
	'PaperPositionMode','auto' ...
	)
print( ...
	spm_figure('FindWin','Graphics'), ...
	'-dpng', ...
	'-r600', ...
	fullfile(out_path,'output4.png') ...
	);

