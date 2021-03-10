function [bins, energy_in_bin] = binned_energy_in_freq_intervals(freqs, mean_freq_signal, n)
%freqs :: columns are eigenvalues
% mean_freq_signal :: columns are freq components of mean signals
% bin_width :: float dictating how large the intervals should be

    min_freq = min(freqs, [], 'all');
    max_freq = max(freqs, [], 'all');
    
    %bins = (min_freq-.0001):bin_width:(max_freq+.0001);
    bins = linspace(min_freq,max_freq,n);
    energy_in_bin = zeros(length(bins)-1,1);

    for l = 1:(length(bins)-1)
        low  = bins(l);
        high = bins(l+1);
    
        % find energy within this bin (e.g. freq range)
        idxs_in_bin  = (freqs>=low) & (freqs<high);
        energy_in_bin(l) = sum(mean_freq_signal(idxs_in_bin).^2, 'all');
    end  

end


