function [energy_in_bins] = energy_in_variation_bins(bins, scan_info_obj, GSO)
    
	[GFT, eigenvals, S] = extract_GFT(subject_id, atlas, include_subcortical, GSO);
	A = extract_sc(subject_id, atlas, include_subcortical);

	V = GFT';
	% ith position is # 0-crossings of ith freq component
	zcs = zero_crossings(V, A);
	tvs = total_variation(V, A);
    
    
    roi_idxs = get_roi_idxs(atlas, include_subcortical);

    %[dtseries] = load_functional_dtseries(atlas, path2fmri, subject, rawdatafolder, chosen_roi);
    x  = load_functional_dtseries(subject_id, atlas, tasktype, scan_dir, raw_hcp_datafolder);
    x = x(roi_idxs,:);
    x_mean = mean(x,2);
    mean_freq_signal = GFT*x_mean;
    
    
    energy_in_bins = zeros(length(bins),1);
    
    num_bins = length(bins)-1;
    for l = 1:num_bins
        low  = bins(l);
        high = bins(l+1);
    
        % find energy within this bin (e.g. freq range)
        idxs_in_bin  = (low<=tvs) & (tvs<high);
        energy_in_bins(l) = sum(mean_freq_signal(idxs_in_bin).^2, 'all');
    end  



end

