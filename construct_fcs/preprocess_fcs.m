function [which_idxs_remove, raw_threshold] = preprocess_fcs(subsets, GFT, fc, filter_params)

which_metric =  filter_params.name;

% which_metric determines which metric to use to determine whether to keep
% or not
%  'eig'
%       => filter_params = {which_eig :: int, threshold :: float}
%  'psd' <-- positive semi-definite metric
%       =>'filter_params' => {psd_metric :: str, threshold :: float}

% use_perc := use threshold as percentile or raw threshold?


% filter_params.threshold is expected to be a number in [0,100]. If
% normalized thresholding is used (currently everything except use_per
% cases) ==> we divide it by 100 bc we are thresholding fractions
if filter_params.use_percentile && ~(0<filter_params.threshold && filter_params.threshold<100)
	error('percentile must be in [0,100]: %f', filter_params.threshold);
end
if (filter_params.threshold<1)
	error('normalization used in threshold...expected number in [0,100]: %f', filter_params.threshold);
end

raw_threshold = filter_params.threshold; %change later if per used


switch which_metric
    case 'None'
        which_idxs_remove = [];
        
    case 'eig'
        k = filter_params.which_eig;
        max_eigs = apply_to_tensor_slices(@(c) eigs(c,k), fc);
        kth_eig = max_eigs(:,end); %column vector: entry i <=>  kth largest eigenvalue by magnitude for ith fc
        
        if filter_params.use_percentile
            raw_threshold = prctile(kth_eig, filter_params.threshold);
            which_idxs_remove = find(kth_eig > raw_threshold);
        else
            which_idxs_remove = find(kth_eig > filter_params.threshold);
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
        
        ave_signal_subsets = apply_to_tensor_slices(@(x) mean(x,2), subsets);
        energy = vecnorm(ave_signal_subsets,2).^2;
        x_freq = GFT*ave_signal_subsets;
        [num_eigs, ~] = size(GFT);
        y = 1:length(num_eigs);
        energy_contributions = energy_in_freq_intervals(x_freq, filter_params.ranges, y);
        lpf_energy = energy_contributions(1,:);
        lpf_energy_frac = lpf_energy./energy;
        
        if filter_params.use_percentile
            raw_threshold = prctile(lpf_energy_frac, filter_params.threshold); %must be in [0,100]
            which_idxs_remove = find(lpf_energy_frac > raw_threshold);
        else
            filter_params.threshold=filter_params.threshold/100; %place in [0,1]
            which_idxs_remove = find(lpf_energy_frac > filter_params.threshold);
        end
    otherwise
        error('filter %s incorrect or not implimented',which_metric);
end

        
end


