function [cached_path, is_cached] = cached_filepath(atlas, tasktype, subject, scan_dir)
%locate where the cached file is stored

%% ensure inputs are correct form. Constructing path to cached data
if contains(atlas, 'desikan', 'IgnoreCase', true)
    %cached_data_folder = '~/Documents/MATLAB/brain_data_preprocess/data/cached_desikan';    
    cached_data_folder = 'data/cached_desikan';
    % if local path not working, be sure to run startup.m
else
    error("Atlas %s not supported yet", atlas)
end

if contains(tasktype, 'REST1', 'IgnoreCase', true)
    tasktype = 'rfMRI_REST1';
elseif contains(tasktype, 'REST2', 'IgnoreCase', true)
    tasktype = 'rfMRI_REST2';
else
    error("Task %s not supported yet (or typo)", tasktype)
end

if isnumeric(subject)
    subject = num2str(subject);
end

if contains(scan_dir, 'LR', 'IgnoreCase', true)
    scan_dir = "LR";
elseif contains(scan_dir, 'RL', 'IgnoreCase', true)
    scan_dir = "RL";
end

filename    = strcat(subject, '_', scan_dir, '.mat');
cached_path = fullfile(cached_data_folder, tasktype, filename);
is_cached   = isfile(cached_path);

end

