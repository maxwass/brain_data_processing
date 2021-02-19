function [which_idxs_remove, raw_threshold] = preprocess_signals(x, GFT, filter_params)
%which signals to keep or remove

which_metric   = filter_params.name;
% which_metric determines which metric to use to determine whether to keep
% or not
%  'energy'
%       => filter_params = {threshold :: float}
%  'freq_distribution' 
%       => filter_params = {ranges :: {[int]}, threshold :: float, patient_info...}
%  'active_frontal_lobe'
%       =>'active_frontal_lobe' => {threshold :: float}

%filter_params fields
%   name = which_metric filter to use
%   threshold = threshold value
%   use_percentile = use percentile or raw threshold?
%   ranges = for freq filtering

filter_params.use_percentile;
% use_percentile := use threshold as percentile or raw threshold?

%%Outputs
%which_idxs_remove - [int] indices removed by threshold
%raw_threshold = threshold if ~use_per, prctile(_,threshold) if use_per


% filter_params.threshold is expected to be a number in [0,100]. If
% normalized thresholding is used (currently everything except use_per
% cases) ==> we divide it by 100 bc we are thresholding fractions
if filter_params.use_percentile && ~(0<filter_params.threshold && filter_params.threshold<100)
	error('percentile must be in [0,100]: %f', filter_params.threshold);
end
if (filter_params.threshold<1)
	error('normalization used in threshold...expected number in [0,100]: %f', filter_params.threshold);
end


raw_threshold = filter_params.threshold;  %change if per used 


if ismember("None", which_metric)
    which_idxs_remove = [];
    
elseif ismember("freq_distribution", which_metric)
    
    %take GFT of signals, break up into freq ranges. Pick metric to map
    %from energy distribution to decision
    if length(filter_params.ranges)~=3
        error('freq_distribution required 3 ranges');
    end
    
    energy = vecnorm(x,2).^2;
    x_freq = GFT*x;
    energy_contributions = energy_in_freq_intervals(x_freq, filter_params.ranges);
    lpf_energy = energy_contributions(1,:);
    lpf_energy_frac = lpf_energy./energy;
    
    if filter_params.use_percentile
        raw_threshold = prctile(lpf_energy_frac, filter_params.threshold); %must be in [0,100]
        which_idxs_remove = find(lpf_energy_frac > raw_threshold);
    else
        filter_params.threshold=filter_params.threshold/100; %place in [0,1]
        which_idxs_remove = find(lpf_energy_frac > filter_params.threshold);
    end


elseif ismember("energy", which_metric)
    %remove observations with total energy above ___
    
    energy = vecnorm(x,2).^2;
    if filter_params.use_percentile
        raw_threshold = prctile(energy, filter_params.threshold);
        which_idxs_remove = find(energy > raw_threshold);
    else
        which_idxs_remove = find(energy > filter_params.threshold);
    end


elseif ismemmber("active_frontal_lobe", which_metric)
    % remove observations where clear thinking is occuring
    error('active frontal lobe not yet implemented');
else
    error('filter %s incorrect or not implimented',which_metric);
end


end

