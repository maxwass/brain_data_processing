function [processed_fcs, which_idxs_remove] = preprocess_fcs(fcs, which_metric, metric_params)
% which_metric determines which metric to use to determine whether to keep
% or not
%  'eig'
%       => filter_params = {which_eig :: int, cutoff :: float}
%  'psd' <-- positive semi-definite metric
%       =>'filter_params' => {psd_metric :: str, cutoff :: float}

% use_perc := use cutoff as percentile or raw threshold?

if use_per && ~(0<filter_params.cutoff<1)
	error('percentile must be in [0,1]: %f', filter_params.cutoff);
end

switch which_metric
    case 'eig'
        k = filter_params.which_eig;
        max_eigs = apply_to_tensor_slices(@(c) eigs(c,k), fcs);
        kth_eig = max_eigs(:,end); %column vector: entry i <=>  kth largest eigenvalue by magnitude for ith fc
        
        if use_per
            p = prctile(kth_eig, metric_params.cutoff);
            which_idxs_remove = kth_eig > p;
        else
            which_idxs_remove = (kth_eig > filter_params.cutoff);
        end
        
        %ensure this is an array of logicals
    case 'psd'
        [~, ~, num_windows] = size(fcs);
        
        dists = zeros(1,num_windows-1);
        
        %measures change/distance between fcs
        for i = 1:(num_windows-1)
            dists(i) = psd_dist(fcs(:,:,i), fcs(:,:,i+1), filter_params.psd_metric);
        end
end

processed_fcs = fcs;
processed_fcs(which_idxs_remove) = [];

end

