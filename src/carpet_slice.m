function [carpetsliceplot,sliceimg] = carpet_slice(fmri_file)

% Carpet plot, slice vs volume

% Load images
Vfmri = spm_vol(fmri_file);
spm_check_orientations(Vfmri);
Yfmri = spm_read_vols(Vfmri);

% Threshold fmri - make non-brain voxels nan
Ym = mean(Yfmri,4);
threshold = spm_antimode(Ym(:));
Ymask = nan(size(Ym));
Ymask(Ym>threshold) = 1;
Yfmri = 100 / mean(Ym(Ymask(:)==1)) * Yfmri;
for v = 1:size(Yfmri,4)
	Yfmri(:,:,:,v) = Yfmri(:,:,:,v) .* Ymask;
end

carpetsliceplot = squeeze(nanmean(nanmean(Yfmri,1),2));

sliceimg = squeeze(Ym(round(size(Ym,1)/2),:,:))';
