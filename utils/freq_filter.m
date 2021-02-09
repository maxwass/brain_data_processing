function [filtered_signal] = freq_filter(signals_freq, freq_interval_idxs)
% signals :: N x num_obsvs [float] of signals in GFT domain
% freq_interval_idxs :: 1 x 2 [int] of frequencies to include


num_freq_idxs = length(freq_interval_idxs);
if length(freq_interval_idxs) > 2
    error('freq_filter: must only give one interval range - and thus 2 idxs, given %d',num_freq_idxs);
end

low = freq_interval_idxs(1);
high= freq_interval_idxs(2);
if low > high %note: low can equal high -> only include one freq
    error('freq_filter: interval idxs must be montonically increasing');
end
    
%only rows (freq components) corresponding to interval are one
mask = zeros(size(signals_freq));
mask( low:high , :) = 1;

filtered_signal = mask.*signals_freq;

end

