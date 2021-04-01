function [ave_node_val] = average_node_value(atlas, roi_idxs)
% Compute the average node value over all mean vectors of all scans.
% To be used mean centering.

if contains(atlas, "desikan", 'IgnoreCase', true)
    sum_stat_file = load('summary_stats/fmri_desikan_summary_stats.mat');
    average_vector_ = mean([sum_stat_file.lrs.mean_vectors, sum_stat_file.rls.mean_vectors], 2);
    average_vector  = average_vector_(roi_idxs,:);
    ave_node_val = mean(average_vector);
else
    error("Atlas %s not supported yet", atlas)
end
end

