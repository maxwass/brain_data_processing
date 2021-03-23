%% plot energy distribution

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


