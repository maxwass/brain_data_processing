%% Plot energy in bin of eigenvalues. Overlay range of some variation 
%   metric inside this interval.
close; clear;

bin_width = .05;
b = 200;
include_subcortical = false;
if include_subcortical
	title_txt = sprintf('Cortical & Subcortical');
else
	title_txt = sprintf('Cortical');
    %tv_ranges.L = [0,405];
    %tv_ranges.L_norm = 
end



cols=1;
t = tiledlayout(4,cols);
fs = 20;
title(t,sprintf('Summed Energy Distribution Over All patients: %s Nodes', title_txt), 'FontSize', fs+3);


axes = gobjects(4,cols);

GSO = 'L';
axes(1,1) = plot_histogram(GSO, include_subcortical, b, false);
if cols==2
    axes(1,2) = plot_histogram(GSO, ~include_subcortical, b, false);
    linkaxes(axes(1,:),'xy')
end


GSO = 'L_norm';
axes(2,1) = plot_histogram(GSO, include_subcortical, b, false);
if cols==2
    axes(2,2) = plot_histogram(GSO, ~include_subcortical, b, false);
    linkaxes(axes(2,:),'xy')
end


GSO = 'A';
axes(3,1) = plot_histogram(GSO, include_subcortical, b, false);
if cols==2
    axes(3,2) = plot_histogram(GSO, ~include_subcortical, b, false);
    linkaxes(axes(3,:),'xy')
end


GSO = 'A_norm';
axes(4,1) = plot_histogram(GSO, include_subcortical, b, false);
if cols==2
    axes(4,2) = plot_histogram(GSO, ~include_subcortical, b, false);
    linkaxes(axes(4,:),'xy')
    xlabel(axes(4,2), 'Eigenvalues');
end


xlabel(axes(4,1), 'Eigenvalues', 'FontSize',fs);




%set(axes, 'YScale', 'log');
%ylabel(axes, 'Log_10( Energy )');
%set(axes, 'YScale', 'linear');

grid(axes, 'MINOR')
%ylim(axes, [0,10^11]);


function ax = plot_histogram(GSO, include_subcortical, num_bins, remove_first)
    filename = sprintf('GSO-%s_desikan_REST1--freq_distrib', GSO);
    if include_subcortical
        filepath = fullfile("energy_distrib", "ed_data", "GFT_using_subcortical", filename);
    else
         filepath = fullfile("energy_distrib", "ed_data", "GFT_only_cortical", filename);
    end
    GSO_data = load(filepath);
    
    min_freq = min(GSO_data.freqs, [], 'all');
    max_freq = max(GSO_data.freqs, [], 'all');
    bins = linspace(min_freq,max_freq, num_bins);

    [energy_in_bin] = binned_energy_in_freq_intervals(bins, GSO_data.freqs, GSO_data.mean_freq_signal);
    [tvs_bin_summary, zcs_bin_summary] = binned_variation_statistics(bins, GSO);
    
    
    if remove_first
        energy_in_bin(1) = 0;
    end
    
    
    ax = nexttile();
    yyaxis left;
    histogram(ax, 'BinEdges', bins, 'BinCounts', energy_in_bin, 'DisplayName', 'Summed Energy');
    
    yyaxis right;
    x = 0.5 * (bins(1:end-1) + bins(2:end));
    min_tv = min(tvs_bin_summary.mins); %(~isnan(tvs_bin_summary.mins)));
    max_tv = max(tvs_bin_summary.maxs); %(~maxisnan(tvs_bin_summary.maxs)));
    tvs_label = sprintf('medians/mins(%.1f)/maxs(%.1f) TV', min_tv, max_tv);
    e1 = errorbar(ax, x, tvs_bin_summary.medians, tvs_bin_summary.mins, tvs_bin_summary.maxs, 'DisplayName', tvs_label, 'Color', 'g');
    %hold on;
    %e2 = errorbar(ax, x, zcs_bin_summary.medians, zcs_bin_summary.mins, zcs_bin_summary.maxs, 'DisplayName', 'Min/Median/Max ZCS', 'Color', 'r');
    ylim([min_tv, max_tv]);
    
    ylabel('Variation');
    hold off;
    legend;
    
    text(ax, .02, .88, GSO, 'Units', "Normalized", 'FontSize',25);
  
end

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




