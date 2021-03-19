function [filter_signals_freq] = freq_filtering_top(signals_freq, t)

[N,~] = size(signals_freq);

if t>=N
    error('When choosing top %d freq components, %d must be less than # eigenvalues %d', t,t,N);
end

%first find which are the top frequencies by energy


signal_freq_energy = signals_freq.^2;
median_energies = median(signal_freq_energy,2);
[~, top_freqs] = maxk(median_energies,t);


%{
[B, I] = sort(signal_freq_energy,1); %sort each column in descending

%find t most common frequencies in top t+5(?) freqs
search_range = min(t+5,N);
I_top = I(1:search_range,:);
I_top_vec = I_top(:);

top = zeros(t,1);
for j = 1:t
    most_common = mode(I_top_vec);
    top(j)      = most_common;
    I_top_vec(I_top_vec == most_common) = [];
end
%}

%% construct filtered signal
filter_signals_freq = zeros(size(signals_freq));
filter_signals_freq(top_freqs,:) = signals_freq(top_freqs,:); 

end
