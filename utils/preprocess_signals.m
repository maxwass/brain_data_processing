function [processed_x] = preprocess_signals(x, which_metric, filter_params, use_per)
%which signals to keep or remove

% which_metric determines which metric to use to determine whether to keep
% or not
%  'energy'
%       => filter_params = {cutoff :: float}
%  'freq_distribution' 
%       => filter_params = {ranges :: {[int]}, cutoff :: float}
%  'active_frontal_lobe'
%       =>'active_frontal_lobe' => {cutoff :: float}

% use_perc := use cutoff as percentile or raw threshold?

if use_per && ~(0<filter_params.cutoff<1)
	error('percentile must be in [0,1]: %f', filter_params.cutoff);
end


if ismember("freq_distribution", which_metric)
    
    %take GFT of signals, break up into freq ranges. Pick metric to map
    %from energy distribution to decision
    if length(filter_params.range)~=3
        error('freq_distribution required 3 ranges');
    end
    
    [low_freq_interval, med_freq_interval, high_freq_interval} = ...
        deal(filter_params.ranges{:}}
    
    [x_freq, ~, ~] = apply_GFT(x, subject, atlas, include_subcortical, GSO);
    [lpf_signal] = freq_filtering(x_freq, {low_freq_interval});
    [mpf_signal] = freq_filtering(x_freq, {med_freq_interval});
    [hpf_signal] = freq_filtering(x_freq, {high_freq_interval});
    
    energy = vecnorm(x_freq,2).^2;
    lpf_energy = vecnorm(lpf_signal,2).^2;
    mpf_energy = vecnorm(mpf_signal,2).^2;
    hpf_energy = vecnorm(hpf_signal,2).^2;
    
    lpf_energy_contrib = lpf_energy./energy;
    mpf_energy_contrib = mpf_energy./energy;
    hpf_energy_contrib = hpf_energy./energy;
    
    if use_per
        p = prctile(lpf_energy_contrib, filter_params.cutoff);
        which_idxs = lpf_energy_contrib > p;
    else
        which_idxs = (lpf_energy_contrib > filter_params.cutoff);
    end


elseif ismember("energy", which_metric)
    %remove observations with total energy above ___
    
    energy = vecnorm(x,2).^2;
    if use_per
        p = prctile(energy, filter_params.cutoff);
        which_idxs = energy > p;
    else
        which_idxs = (energy > filter_params.cutoff);
    end


elseif ismemmber("active_frontal_lobe", which_metric)
    % remove observations where clear thinking is occuring
    error('active frontal lobe not yet implemented');
end

processed_x = x;
processed_x(which_idxs) = [];


end

