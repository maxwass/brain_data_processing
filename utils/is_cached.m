function [cached, cached_filename] = is_cached(subject, atlas, chosen_roi, path2fmri)
%subject = string/char array of subject id
%atlas :: in {'desikan', 'destrieux'}
%path2fmri :: path to raw fmri file
% chosen_roi is a struct of ROI's from cortical and subcortical regions
%   chosen_roi.cortical:    [int]
%   chosen_roi.subcortical: [int]
% include_subcortical :: logical
% subcortical_first   :: logical

load_and_check = false;

% tasktype: str 'rfMRI_REST1'. Possibly include REST2 in future
if contains(path2fmri, "REST1")
    tasktype = "rfMRI_REST1";
elseif contains(path2fmri, "REST2")
    tasktype = "rfMRI_REST2";
else
    error('Have not yet handled non-rest taskes: %s', path2fmri);
end

% scan_dir = scan direction: str in {LR, RL}
if contains(path2fmri, "LR")
    scan_dir = "LR"; 
elseif contains(path2fmri, "RL")
    scan_dir = "RL";
else
    error('Could not detect scan dir from fmri file name');
end


cached_filename = sprintf('cached_%s/%s/%s_%s.mat', atlas, tasktype, subject, scan_dir);

if isfile(cached_filename)
    cached = true;
    
    %should we load it and check params?
    if load_and_check
        cached_fields = load(cached_filename, 'chosen_roi', 'scan', 'tasktype', 'atlas');
        if contains(cached_fields.scan, "RL")
            cached_scan_dir = "RL";
        else
            cached_scan_dir = "LR";
        end
    
        if (isequal(chosen_roi.cortical, cached_fields.chosen_roi.cortical) && ...
        isequal(chosen_roi.subcortical, cached_fields.chosen_roi.subcortical) && ...
        isequal(scan_dir, cached_scan_dir) && ...
        isequal(tasktype, cached_fields.tasktype) && ...
        isequal(atlas, cached_fields.atlas) )
            correct = true;
        else
            correct = false;
            error('Incorrect parameters in cached file...');
        end
    end
  
else
    cached = false;
end



end

