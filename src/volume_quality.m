function [med_displc,glob,FD,DVARS,mvt_file,tsnr_file,stats] = volume_quality( ...
	out_path, ...
	meanfmri_file, ...
	fmri_file, ...
	rp_file ...
	)


% Load motion params
rp = load(rp_file);
rpdeg = rp;
rpdeg(:,4:6) = rpdeg(:,4:6) * 180 / pi;
nvol = size(rp,1);

% Use mean fMRI to get brain voxel mask
Vmean = spm_vol(meanfmri_file);
[Ymean,XYZ] = spm_read_vols(Vmean);
threshold = spm_antimode(Ymean(:));
mask = Ymean > threshold;


%% Voxel displacements

% Keep just in-mask coords
XYZmask = XYZ(:,mask(:)>0);

% Add 1 to make homogeneous coords
XYZmask = [XYZmask; ones(1,size(XYZmask,2))];

% Compute voxel displacements at each in-mask voxel
displ = zeros(size(XYZmask,2),nvol);
for t = 2:nvol
	m = spm_matrix(rp(t,:));
	mprev = spm_matrix(rp(t-1,:));
	pXYZmask = (mprev\m) * XYZmask;
	displ(:,t) = sqrt(sum((pXYZmask-XYZmask).^2));
end

% Median displacement at each volume
med_displc = prctile(displ,50)';
save(fullfile(out_path,'median_voxel_displacement_mm.txt'),'med_displc','-ascii')
med_displc(1) = nan;

% Overall 95th percentile of voxel displacement
displc95 = prctile(displ(:),95);

% Movement summary at each voxel
Ymvt = zeros(size(Ymean));
Ymvt(mask(:)>0) = prctile(displ',95);
Vmvt = rmfield(Vmean,'pinfo');
mvt_file = fullfile(out_path,'voxel_displacement_mm_95prctile.nii');
Vmvt.fname = mvt_file;
spm_write_vol(Vmvt,Ymvt);


%% Framewise displacement

% We'll calculate the frame-to-frame displacement of points on a sphere of
% radius 50. The matrix at a point gives its displacement relative to zero
% such that
%     Xt = Mt * X0
% We also know
%     Xtm1 = Mtm1 * X0
% We want M such that
%     Xt  = M * Xtm1
% So we get
%     Mt * X0 = M * Mtm1 * X0
%     (Mt*X0) * inv(Mtm1*X0) = M = Mt * X0 * inv(X0) * inv(Mtm1) = Mt/Mtm1

% Initialize some variables used for the Power method
amm = nan(nvol,1);
bmm = nan(nvol,1);
gmm = nan(nvol,1);
FD = zeros(nvol,1);

% Examine each time point separately
for t = 2:nvol
	
	% Displacement from previous vol, Power et al 2012 method. This is a
	% sum of displacements over the six degrees of freedom, therefore more
	% approximate than the above.
	radius = 50;
	fdM = spm_matrix(rp(t,:)) / spm_matrix(rp(t-1,:));
	param = spm_imatrix(fdM);
	qa = spm_matrix([0 0 0 param(4) 0 0]) * [0 radius 0 0]';
	amm(t) = sqrt(sum( (qa-[0 radius 0 0]').^2 ));
	qb = spm_matrix([0 0 0 0 param(5) 0]) * [0 0 radius 0]';
	bmm(t) = sqrt(sum( (qb-[0 0 radius 0]').^2 ));
	qg = spm_matrix([0 0 0 0 0 param(6)]) * [radius 0 0 0]';
	gmm(t) = sqrt(sum( (qg-[radius 0 0 0]').^2 ));
	
	FD(t) = sum(abs( ...
		[rp(t,1)-rp(t-1,1) rp(t,2)-rp(t-1,2) rp(t,3)-rp(t-1,3) ...
		amm(t) bmm(t) gmm(t)] ));
	
end

% Save to file
save(fullfile(out_path,'FD.txt'),'FD','-ascii');
FD(1) = nan;


%% DVARS

% Frame by frame RMS difference in voxel intensities, percent change units
% relative to the mean intensity of the brain. Compute for each volume
% separately to save loading the whole fmri into memory
V1 = spm_vol([fmri_file ',1']);
Y1 = spm_read_vols(V1);
DVARS = zeros(nvol,1);
for t = 2:nvol
	Vt = spm_vol([fmri_file ',' num2str(t)]);
	spm_check_orientations([V1;Vt]);
	Yt = spm_read_vols(Vt);
	DVARS(t) = 100  ./ mean(Ymean(mask)) .* sqrt(mean( (Yt(mask)-Y1(mask)).^2 ));
	V1 = Vt;
	Y1 = Yt;
end
save(fullfile(out_path,'DVARS.txt'),'DVARS','-ascii')
DVARS(1) = nan;


%% Global time series and temporal SNR
Y = spm_read_vols(spm_vol(fmri_file));
osize = size(Y);
Y = reshape(Y,[],osize(4))';
Y(:,~mask) = nan;
gm = nanmean(Y(:));
Y = 100 / gm * Y;
glob = nanmean(Y,2);
save(fullfile(out_path,'global.txt'),'glob','-ascii')

% Linear detrend before computing SNR
%tsnr = nanmean(Y) ./ nanstd(detrend(Y));

% Robust TSNR
tsnr = nanmean(Y) ./ (1.253*mad(detrend(Y)));

% Save TSNR image
Ytsnr = reshape(tsnr,osize(1:3));
Vtsnr = rmfield(Vmean,'pinfo');
tsnr_file = fullfile(out_path,'temporal_snr.nii');
Vtsnr.fname = tsnr_file;
spm_write_vol(Vtsnr,Ytsnr);


%% Fill stats structure and write to file
stats = struct( ...
	'fd_mean',nanmean(FD), ...
	'dvars_mean',nanmean(DVARS), ...
	'tsnr_robust_median',nanmedian(tsnr), ...
	'global_temporal_stddev',std(glob), ...
	'voxel_displacement_mm_95prctile',displc95, ...
	'maxtrans_firstvol_mm_deprecated',max(max(abs(rpdeg(:,1:3)))), ...
	'maxrot_firstvol_deg_deprecated',max(max(abs(rpdeg(:,4:6)))) ...
	);

fid = fopen(fullfile(out_path,'fmriqa_v4_stats.csv'),'wt');
fn = fieldnames(stats);
for s = 1:length(fn)
	fprintf(fid,'%s,%0.3f\n',fn{s},stats.(fn{s}));
end
fclose(fid);
