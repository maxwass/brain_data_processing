function [which_idxs_remove, raw_cutoff] = preprocess_signals(x, GFT, which_metric, filter_params, use_per)
%which signals to keep or remove

% which_metric determines which metric to use to determine whether to keep
% or not
%  'energy'
%       => filter_params = {cutoff :: float}
%  'freq_distribution' 
%       => filter_params = {ranges :: {[int]}, cutoff :: float, patient_info...}
%  'active_frontal_lobe'
%       =>'active_frontal_lobe' => {cutoff :: float}

% use_perc := use cutoff as percentile or raw threshold?

%%Outputs
%which_idxs_remove - [int] indices removed by cutoff
%raw_cutoff = cutoff if ~use_per, prctile(_,cutoff) if use_per


% filter_params.cutoff is expected to be a number in [0,100]. If
% normalized thresholding is used (currently everything except use_per
% cases) ==> we divide it by 100 bc we are thresholding fractions
if use_per && ~(0<filter_params.cutoff && filter_params.cutoff<100)
	error('percentile must be in [0,100]: %f', filter_params.cutoff);
end
if (filter_params.cutoff<1)
	error('normalization used in cutoff...expected number in [0,100]: %f', filter_params.cutoff);
end


raw_cutoff = filter_params.cutoff;  %change if per used 


if ismember("None", which_metric)
    which_idxs_remove = [];
    
elseif ismember("freq_distribution", which_metric)
    
    %take GFT of signals, break up into freq ranges. Pick metric to map
    %from energy distribution to decision
    if length(filter_params.ranges)~=3
        error('freq_distribution required 3 ranges');
    end

    %fp = filter_params;
    %[x_freq, ~, ~] = apply_GFT(x, fp.subject, fp.atlas, fp.include_subcortical, fp.GSO);
    
    energy = vecnorm(x,2).^2;
    x_freq = GFT*x;
    energy_contributions = energy_in_freq_intervals(x_freq, filter_params.ranges);
    lpf_energy = energy_contributions(1,:);
    lpf_energy_frac = lpf_energy./energy;
    
    if use_per
        raw_cutoff = prctile(lpf_energy_frac, filter_params.cutoff); %must be in [0,100]
        which_idxs_remove = find(lpf_energy_frac > raw_cutoff);
    else
        filter_params.cutoff=filter_params.cutoff/100; %place in [0,1]
        which_idxs_remove = find(lpf_energy_frac > filter_params.cutoff);
    end


elseif ismember("energy", which_metric)
    %remove observations with total energy above ___
    
    energy = vecnorm(x,2).^2;
    if use_per
        raw_cutoff = prctile(energy, filter_params.cutoff);
        which_idxs_remove = find(energy > raw_cutoff);
    else
        which_idxs_remove = find(energy > filter_params.cutoff);
    end


elseif ismemmber("active_frontal_lobe", which_metric)
    % remove observations where clear thinking is occuring
    error('active frontal lobe not yet implemented');
end


end

