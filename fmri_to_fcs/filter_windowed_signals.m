function [filtered_signals] = filter_windowed_signals(signals, subject, atlas, include_subcortical, which_filters)
%Each column in signals is a average windowed signal over dtseries.



[N, num_signals] = size(signals);

filter_idxs = true(size(1,num_signals));

if ismember("energy", which_filters)
    signal_norms   = vecnorm(signals,2);
    low_energy_idxs    = (signal_norms < mean(signal_norms) );
    filter_idxs = filter_idxs && low_energy_idxs;
end
if ismember("high_freq_energy", which_filters)
    %extract_Sc
    [A] = extract_sc(subject, atlas, include_subcortical);
    [A_norm, L, L_norm] = compute_GSOs(A);
    
    high_freq_energy = 
    
    filter_idxs = filter_idxs && low_energy_idxs;
    
end
outputArg1 = inputArg1;
filtered_signals = inputArg2;
end


filtered_signals = signals(low_energy_idxs);
end

