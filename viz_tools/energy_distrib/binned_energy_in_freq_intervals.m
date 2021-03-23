function [energy_in_bin] = binned_energy_in_freq_intervals(bins, freqs, mean_freq_signal)
%freqs :: columns are eigenvalues
% mean_freq_signal :: columns are freq components of mean signals
% bin_width :: float dictating how large the intervals should be

    num_bins = length(bins)-1;
    energy_in_bin = zeros(num_bins,1);

    for l = 1:num_bins
        low  = bins(l);
        high = bins(l+1);
    
        % find energy within this bin (e.g. freq range)
        idxs_in_bin  = (low<=freqs) & (freqs<high);
        energy_in_bin(l) = sum(mean_freq_signal(idxs_in_bin).^2, 'all');
    end  

end


