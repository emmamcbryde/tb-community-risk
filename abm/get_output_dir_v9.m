function outDir = get_output_dir_v9()
%GET_OUTPUT_DIR_V9 Return the ABM output folder and create it if needed.

thisFile = mfilename('fullpath');
[thisDir, ~, ~] = fileparts(thisFile);

outDir = fullfile(thisDir, 'output');

if ~exist(outDir, 'dir')
    mkdir(outDir);
end
end