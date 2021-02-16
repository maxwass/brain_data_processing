function [which_idxs_remove, raw_cutoff] = preprocess_fcs(x, GFT, fc, which_metric, filter_params, use_per)
% which_metric determines which metric to use to determine whether to keep
% or not
%  'eig'
%       => filter_params = {which_eig :: int, cutoff :: float}
%  'psd' <-- positive semi-definite metric
%       =>'filter_params' => {psd_metric :: str, cutoff :: float}

% use_perc := use cutoff as percentile or raw threshold?

if use_per && ~(0<filter_params.cutoff && filter_params.cutoff<100)
	error('percentile must be in [0,1]: %f', filter_params.cutoff);
end

raw_cutoff = filter_params.cutoff; %change later if per used


switch which_metric
    case 'None'
        which_idxs_remove = [];
        
    case 'eig'
        k = filter_params.which_eig;
        max_eigs = apply_to_tensor_slices(@(c) eigs(c,k), fc);
        kth_eig = max_eigs(:,end); %column vector: entry i <=>  kth largest eigenvalue by magnitude for ith fc
        
        if use_per
            raw_cutoff = prctile(kth_eig, filter_params.cutoff);
            which_idxs_remove = find(kth_eig > raw_cutoff);
        else
            which_idxs_remove = find(kth_eig > filter_params.cutoff);
        end
        
        %ensure this is an array of logicals
    case 'psd'
        [~, ~, num_windows] = size(fc);
        
        dists = zeros(1,num_windows-1);
        
        %measures change/distance between fc
        for i = 1:(num_windows-1)
            dists(i) = psd_dist(fc(:,:,i), fc(:,:,i+1), filter_params.psd_metric);
        end
        
    case 'freq_distribution'
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
        
        
        if use_per
            raw_cutoff = prctile(lpf_energy_frac, filter_params.cutoff); %must be in [0,100]
            which_idxs_remove = find(lpf_energy_frac > raw_cutoff);
        else
            filter_params.cutoff=filter_params.cutoff/100; %place in [0,1]
            which_idxs_remove = find(lpf_energy_frac > filter_params.cutoff);
        end

        
end

end

