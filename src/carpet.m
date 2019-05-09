function [carpetplot,labels,mask] = carpet(seg_file,fmri_file)

% Carpet plots of fMRI data as described by
%    Power JD. A simple but useful way to assess fMRI scan qualities.
%    Neuroimage 2017 154:150-158
%    doi: 10.1016/j.neuroimage.2016.08.009
%    PubMed PMID: 27510328
%    PubMed Central PMCID: PMC5296400
%    https://www.ncbi.nlm.nih.gov/pubmed/27510328


% Load images
Vseg = spm_vol(seg_file);
Vfmri = spm_vol(fmri_file);
spm_check_orientations([Vseg; Vfmri]);
Yseg = spm_read_vols(Vseg);
Yfmri = spm_read_vols(Vfmri);

% Threshold fmri and apply to masks
Ym = mean(Yfmri,4);
threshold = spm_antimode(Ym(:));
Yseg(Ym<threshold) = 0;

% Define compartments
clear cpts
cpts.Lfrontal = [101 105 113 119 121 125 137 139 141 143 147 153 163 165 ...
	179 187 191 205];
cpts.Rfrontal = [100 104 112 118 120 124 136 138 140 142 146 152 162 164 ...
	178 186 190 204];
cpts.Lparietal = [107 167 169 175 195 199];
cpts.Rparietal = [106 166 168 174 194 198];
cpts.Lsensorimotor = [151 183 193 149 177];
cpts.Rsensorimotor = [150 182 192 148 176];
cpts.Loccipital = [109 115 129 135 145 157 161 197];
cpts.Roccipital = [108 114 128 134 144 156 160 196];
cpts.Ltemporal = [117 123 133 155 171 181 185 201 203 207];
cpts.Rtemporal = [116 122 132 154 170 180 184 200 202 206];
cpts.subcortical = [31 32 35 47 48 57 58 59 60 102 103 172 173 75 76 36 ...
	37 23 30 55 56 61 62];
cpts.cerebellum = [38 39 71 72 73];
cpts.white = [40 41 44 45];
cpts.csf = [4 11 49 50 51 52];

% Get index images for each compartment
clear masks
cptnames = fieldnames(cpts);
for c = 1:length(cptnames)
	masks.(cptnames{c}) = ismember( Yseg, cpts.(cptnames{c}) );
end

% Split white matter compartment into inner and outer by eroding
nhood = nan(3,3,3);
nhood(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
nhood(:,:,2) = [0 1 0; 1 1 1; 0 1 0];
nhood(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
masks.white_inner = imerode(masks.white,nhood);
masks.white_outer = masks.white & ~masks.white_inner;
masks = rmfield(masks,'white');
masknames = fieldnames(masks);

% Smooth and extract for each compartment
clear cptplots
for c = 1:length(masknames)
	cptplots.(masknames{c}) = smooth_and_extract(Yfmri,masks.(masknames{c}));
end

% Save memory
clear Yfmri

% Final plot
cp1 = [ ...
	cptplots.Lfrontal; ...
	cptplots.Lsensorimotor; ...
	cptplots.Lparietal; ...
	cptplots.Loccipital; ...
	cptplots.Ltemporal ...
	];
cp2 = [ ...
	cptplots.Rfrontal; ...
	cptplots.Rsensorimotor; ...
	cptplots.Rparietal; ...
	cptplots.Roccipital; ...
	cptplots.Rtemporal ...
	];
cp3 = cptplots.subcortical;
cp4 = cptplots.cerebellum;
cp5 = cptplots.white_outer;
cp6 = cptplots.white_inner;
cp7 = cptplots.csf;

% Save memory
clear cptplots

carpetplot = [cp1; cp2; cp3; cp4; cp5; cp6; cp7];

labels = [ ...
	1 * ones(size(cp1,1),1); ...
	2 * ones(size(cp2,1),1); ...
	3 * ones(size(cp3,1),1); ...
	4 * ones(size(cp4,1),1); ...
	5 * ones(size(cp5,1),1); ...
	6 * ones(size(cp6,1),1); ...
	7 * ones(size(cp7,1),1) ...
	];

mask = zeros(size(Yseg));
for c = 1:length(masknames)
	switch masknames{c}
		case {'Lfrontal','Lsensorimotor','Lparietal','Loccipital','Ltemporal'}
			mask(masks.(masknames{c})) = 1;
		case {'Rfrontal','Rsensorimotor','Rparietal','Roccipital','Rtemporal'}
			mask(masks.(masknames{c})) = 2;
		case 'subcortical'
			mask(masks.(masknames{c})) = 3;
		case 'cerebellum'
			mask(masks.(masknames{c})) = 4;
		case 'white_outer'
			mask(masks.(masknames{c})) = 5;
		case 'white_inner'
			mask(masks.(masknames{c})) = 6;
		case 'csf'
			mask(masks.(masknames{c})) = 7;
		otherwise
	end
end

