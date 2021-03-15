%% plot energy distribution

close;

bin_width = .05;
b = 200;
include_subcortical = false;

t = tiledlayout(4,2);
title(t,'Energy Distribution over all patients')


axes = gobjects(4,2);

GSO = 'L_norm';
axes(1,1) = plot_histogram(GSO, include_subcortical, b, false);
axes(1,2) = plot_histogram(GSO, ~include_subcortical, b, false);
% overlay total variations/zero crossings?
linkaxes(axes(1,:),'xy')

GSO = 'L';
axes(2,1) = plot_histogram(GSO, include_subcortical, b, false);
axes(2,2) = plot_histogram(GSO, ~include_subcortical, b, false);
linkaxes(axes(2,:),'xy')

GSO = 'A';
axes(3,1) = plot_histogram(GSO, include_subcortical, b, false);
axes(3,2) = plot_histogram(GSO, ~include_subcortical, b, false);
linkaxes(axes(3,:),'xy')

GSO = 'A_norm';
axes(4,1) = plot_histogram(GSO, include_subcortical, b, false);
axes(4,2) = plot_histogram(GSO, ~include_subcortical, b, false);
linkaxes(axes(4,:),'xy')
%set(axes, 'YScale', 'log');
%ylabel(axes, 'Log_10( Energy )');
%set(axes, 'YScale', 'linear');

grid(axes, 'MINOR')
%ylim(axes, [0,10^11]);


function ax = plot_histogram(GSO, include_subcortical, bin_width, remove_first)
    filename = sprintf('GSO-%s_desikan_REST1--freq_distrib', GSO);
    if include_subcortical
        filepath = fullfile("energy_distrib", "ed_data", "GFT_using_subcortical", filename);
    else
         filepath = fullfile("energy_distrib", "ed_data", "GFT_only_cortical", filename);
    end
    GSO_data = load(filepath);
    
    
    [bins, energy_in_bin] = binned_energy_in_freq_intervals(GSO_data.freqs, GSO_data.mean_freq_signal, bin_width);
    if remove_first
        energy_in_bin(1) = 0;
    end
    
    ax = nexttile();
    histogram(ax, 'BinEdges', bins, 'BinCounts', energy_in_bin);
    
    title_txt = sprintf('%s', GSO);
    if remove_first
        title_txt = sprintf('%s. 0 Freq Removed.',title_txt);
    end
    if include_subcortical
        title_txt = sprintf('%s w/ subcortical',title_txt);
    else
        title_txt = sprintf('%s only cortical ',title_txt);
    end
    title(ax, title_txt);
    
end


