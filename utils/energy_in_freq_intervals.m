function [energy_contributions] = energy_in_freq_intervals(x_freq, intervals, y)
%find energy in the intervals (idxs of eigenvalues) given

%% Inputs
%x_freq :: [float] NxM - each column is a signal. Treated as frequency
%           signals
% y :: vector of values to which the intervals refer
%intervals is a [cell] of 1x2 [int]. Assumed to be non-overlapping valid
%   intervals. Interpreted as the indices of the eigenvalues to keep.

%% Outputs
% energy_contributions : [float] length(intervals)xM energy in
% each interval for each signal (column vector)

    [~, num_obsvs] = size(x_freq);
    
    energy_contributions = zeros(length(intervals),num_obsvs);
    for i = 1:length(intervals)

        %filter each signal (column) for components in this frequency
        %interval
        [idxs_to_keep] = ~intervals_to_logical_vec({intervals{i}}, y);
        x_freq_band = x_freq(idxs_to_keep,:);
        
        %x_freq_band = freq_filtering_idx(x_freq, {intervals{i}});
        
        %compute energy in this frequency interval for each signal
        energy_contributions(i,:) = vecnorm(x_freq_band,2).^2;
    end
    
end

