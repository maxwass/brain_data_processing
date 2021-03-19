function [energy_in_bins] = energy_in_variation_bins(bins, subject_id, atlas, include_subcortical, GSO, chosen_roi, path2fmri)
    
	[GFT, eigenvals] = extract_GFT(subject_id, atlas, include_subcortical, GSO);
	A = extract_sc(subject_id, atlas, include_subcortical);
	S = which_GSO(GSO, A);

	V = GFT';
	L_sign = diag(sum(sign(A),2))-sign(A);
	% ith position is # 0-crossings of ith freq component
	zero_crossings = (1/4)*total_variation(sign(V), L_sign);

	L = diag(sum(A,2)) - A;
	total_variations = total_variation(V, L);
    
    
    [dtseries] = load_functional_dtseries(atlas, path2fmri, subject, rawdatafolder, chosen_roi);
    
    
    
    num_bins = length(bins)-1;
    for l = 1:num_bins
        low  = bins(l);
        high = bins(l+1);
    
        % find energy within this bin (e.g. freq range)
        idxs_in_bin  = (low<=total_variations) & (total_variations<high);
        energy_in_bins(l) = sum(mean_freq_signal(idxs_in_bin).^2, 'all');
    end  



end

