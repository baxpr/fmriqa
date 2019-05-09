function vercheck(matlab_ver,spm_path,spm_ver,vbm_ver)

% Check versions and initialize
%
% Examples:
%    matlab_ver = 'R2015a'
%    spm_ver = 'SPM12 (6225)'
%    vbm_ver = 'vbm8'
%
% spm_ver and vbm_ver can be [] if not needed

if verLessThan('matlab',matlab_ver)
	error('Want Matlab %s or newer, but got %s',matlab_ver,version('-release'));
end

if ~isempty(spm_ver)

	addpath(spm_path);
	
	spm_got = spm('version');
	if ~strcmp(spm_ver,spm_got)
		error('Want %s, but got %s',spm_ver,spm_got);
	end

	if ~isempty(vbm_ver)
		if ~exist([spm_path '/toolbox/' vbm_ver],'dir')
			error(['VBM8 toolbox not found in ' spm_path '/toolbox/' vbm_ver])
		end
	end

end

