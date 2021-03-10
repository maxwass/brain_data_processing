%% plot energy distribution

bin_width = .1;
b = 200;
t = tiledlayout(2,2);

axes = gobjects(4,1);

GSO = 'L_norm';
axes(1) = plot_histogram(GSO, b, false);

GSO = 'L';
axes(2) = plot_histogram(GSO, b, false);

GSO = 'A';
axes(3) = plot_histogram(GSO, b, false);

GSO = 'A_norm';
axes(4) = plot_histogram(GSO, b, false);

set(axes, 'YScale', 'log');
%set(axes, 'YScale', 'linear');


function ax = plot_histogram(GSO, bin_width, remove_first)
    filename = sprintf('GSO-%s_desikan_REST1--freq_distrib', GSO);
    GSO_data = load(filename);
    [bins, energy_in_bin] = binned_energy_in_freq_intervals(GSO_data.freqs, GSO_data.mean_freq_signal, bin_width);
    if remove_first
        energy_in_bin(1) = 0;
    end
    ax = nexttile();
    histogram(ax, 'BinEdges', bins, 'BinCounts', energy_in_bin);
    xlabel(ax, 'freq')
    title_txt = sprintf('%s', GSO);
    if remove_first
        title_txt = sprintf('%s. 0 Freq Removed',title_txt);
    end
    title(ax, title_txt);
    ylabel(ax, 'Energy');
end


