function [tvs_bin_summary, zcs_bin_summary] = binned_variation_statistics(bins, GSO)
%% 1) Collect tvs and zcs of freq componenets of all SCs (in GSO form)
%  2) Split these into respective bins
%  3) Compute statistics (min/max/median/ave of tvs/zcs in each bin


scs_file = load('scs_desikan.mat');
subject_list = scs_file.subject_list;
num_bins = length(bins)-1;

all_tvs = cell(num_bins,1);
all_zcs = cell(num_bins,1);

for idx = 1:length(subject_list)
    subject_id = subject_list(idx);
    filename = sprintf("%d.mat", subject_id);
    filepath = fullfile("viz_scs", "scs_spectral_info", filename);
    gso_struct = load(filepath).gso_struct;
    s = gso_struct.(GSO);
    [eigvals, tvs, zcs] = deal(s.eigenvals, s.total_variations, s.zero_crossings);
   
    
    for bin_idx = 1:num_bins
        [low, high] = deal(bins(bin_idx), bins(bin_idx+1));
        evals_idxs = (low <= eigvals) & (eigvals < high);
        
        all_tvs{bin_idx} = [all_tvs{bin_idx}, tvs(evals_idxs)'];
        all_zcs{bin_idx} = [all_zcs{bin_idx}, zcs(evals_idxs)'];
    end
    
    
end

% take mean/meadian/min/max of tvs of each bin
tvs_bin_summary.medians = cellfun(@median, all_tvs);
tvs_bin_summary.means = cellfun(@mean, all_tvs);
tvs_bin_summary.stds = cellfun(@std, all_tvs);
tvs_bin_summary.mins = min_max_cellfun(@min, all_tvs);
tvs_bin_summary.maxs = min_max_cellfun(@max, all_tvs);


zcs_bin_summary.medians = cellfun(@median, all_zcs);
zcs_bin_summary.means = cellfun(@mean, all_zcs);
tvs_bin_summary.stds = cellfun(@std, all_zcs);
zcs_bin_summary.mins = min_max_cellfun(@min, all_zcs);
zcs_bin_summary.maxs = min_max_cellfun(@max, all_zcs);


end

function out = min_max_cellfun(anon_f, cell_arr)

    out = zeros(1,length(cell_arr));
    for l = 1:length(cell_arr)
        o = anon_f(cell_arr{l});
        if isempty(o)
            out(l) = nan;
        else  
            out(l) = o;
        end
        
    end
    
end
