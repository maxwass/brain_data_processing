function [roi_idxs] = get_roi_idxs(atlas, include_subcortical)
%% Find the subset of indices of interest from the full fMRI signals and SCs
% full fMRI signals and SCs are assumed to be ordered in the same way.

%% subcortical nodes assumed to be the first 19 nodes. Then cortical nodes.

if contains(atlas, "desikan", 'IgnoreCase', true)
    if include_subcortical
        roi_idxs = (1:87);
    else
        roi_idxs = (20:87);
    end
else
    error("Atlas %s not supported yet", atlas);

end

