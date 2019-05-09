%MATLAB
addpath(genpath('${code_path}'))
fmriqa_v4_1( ...
	'${magick_path}', ...
	'${spm12_path}', ...
	'${temp_dir}', ...
	'${t1_file}', ...
	'${seg_file}', ...
	'${fmri_file}', ...
	'${project}', ...
	'${subject}', ...
	'${session}', ...
	'${scan}' ...
	);
