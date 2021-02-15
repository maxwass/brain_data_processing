function plot_scalars(ax, signals, signals_freq, fcs, which_scalars, extras)
%% Compute and plot scalars to plot on single axes
% for (possibly windowed) processed signals from (single) fmri scan,
%  plot scalar metrics to ax

%inputs: 
% signals: time signals
% signals_freq: same signals^ in freq space
% fcs: 

%% 
[N, ~, num_windows] = size(fcs);
[num_evals, num_obsvs] = size(signals_freq);


low_freq_cutoffs  = floor(num_evals/3);
med_freq_cutoffs  = 2*low_freq_cutoffs;

low_freq_interval = [1, low_freq_cutoffs];
med_freq_interval = [low_freq_cutoffs+1,  med_freq_cutoffs];
high_freq_interval= [med_freq_cutoffs+1, num_evals];
intervals = {low_freq_interval, med_freq_interval, high_freq_interval};
interval_labels = ['lpf', 'mpf', 'hpf'];
line_colors = ['r','m','c'];




%% scalars directly from signals/frequency signals
if isequal(which_scalars, "dtseries_energy")
    energy = vecnorm(signals_freq,2).^2;
    [energy_in_range] = energy_in_freq_intervals(signals_freq, intervals);
    
    plot(ax, energy, 'k', 'LineWidth', 1, 'DisplayName', 'signal energy');
    hold(ax, 'on');
    
    for j = 1:length(intervals)
        interval = intervals{j};
        line_color = line_colors(j);
        interval_label = sprintf('%s [%d,%d]',interval_labels(j), interval(1), interval(2));
        plot(ax, energy_in_range(j,:), line_color, 'LineWidth', 1, 'DisplayName', interval_label);
    end
    
    
    
    ylabel(ax, 'Energy', 'FontSize', 15);
    set(ax, 'YColor', 'k')
    hold(ax, 'off');

elseif isequal(which_scalars, "wdw_ave_energy")

    %{
    [lpf_signal] = freq_filtering(signals_freq, {low_freq_interval});
    [mpf_signal] = freq_filtering(signals_freq, {med_freq_interval});
    [hpf_signal] = freq_filtering(signals_freq, {high_freq_interval});

    energy = vecnorm(signals_freq,2).^2;
    lpf_energy = vecnorm(lpf_signal,2).^2;
    mpf_energy = vecnorm(mpf_signal,2).^2;
    hpf_energy = vecnorm(hpf_signal,2).^2;
    %}
    energy = vecnorm(signals_freq,2).^2;
    [energy_in_range] = energy_in_freq_intervals(signals_freq, intervals);
    %yyaxis(ax,'left'); %main metric on left
    plot(ax, energy, 'k', 'LineWidth', 1, 'DisplayName', 'signal energy');
    hold(ax, 'on');
    ylabel(ax, 'Energy', 'FontSize', 15);
    set(ax, 'YColor', 'k')
    
    %yyaxis(ax,'right'); %main metric on left
    for j = 1:length(intervals)
        interval = intervals{j};
        line_color = line_colors(j);
        interval_label = sprintf('%s [%d,%d]',interval_labels(j), interval(1), interval(2));
        %plot(ax, lpf_energy./energy, 'r', 'LineWidth', 1, 'DisplayName', lpr_int);
        plot(ax, energy_in_range(j,:), line_color, 'LineWidth', 1, 'DisplayName', interval_label);
    end
    
    %{
    interval_label = sprintf('%s [%d,%d]',interval_labels{i}, low_freq_interval(1), low_freq_interval(2));
    %plot(ax, lpf_energy./energy, 'r', 'LineWidth', 1, 'DisplayName', lpr_int);
    plot(ax, lpf_energy, 'r', 'LineWidth', 1, 'DisplayName', interval_label);
    
    mpr_int = sprintf('mpf [%d,%d]', med_freq_interval(1), med_freq_interval(2));
    %plot(ax, mpf_energy./energy, 'm', 'LineWidth', 1, 'DisplayName', mpr_int);
    plot(ax, mpf_energy, 'm', 'LineWidth', 1, 'DisplayName', mpr_int);
    
    hpr_int = sprintf('hpf [%d,%d]', high_freq_interval(1), high_freq_interval(2));
    %plot(ax, hpf_energy./energy, 'c', 'LineWidth', 1, 'DisplayName', hpr_int);
    plot(ax, hpf_energy, 'c', 'LineWidth', 1, 'DisplayName', hpr_int);
    %}

    %ylabel(ax, '% energy in freq ranges', 'FontSize', 15);
    %set(ax, 'YColor', 'k') % set y axis color to same as line color
    hold(ax, 'off');

    
elseif isequal(which_scalars, "fc_traj")
    %max_eigs is the maximum magnitude eigenvector for each covarianace matrix
    %max_eig  = @(c) eigs(c,1);
    max_two_eigs  = @(c) eigs(c,2);
    max_eigs = apply_to_tensor_slices(max_two_eigs, fcs);

    % *_diffs is a vector of distance between adjacent matrices
    airm_diffs   = zeros(num_windows-1, 1);
    stein_diffs  = zeros(num_windows-1, 1);
    herdin_diffs = zeros(num_windows-1, 1);
    inner_product    = zeros(num_windows, 1);


    for j = 1:(num_windows-1)
        %airm_diffs(j)   = psd_dist( covs(:,:,j), covs(:,:,j+1), 'airm');
        stein_diffs(j)  = psd_dist( fcs(:,:,j), fcs(:,:,j+1), 'stein');
        herdin_diffs(j) = psd_dist( fcs(:,:,j), fcs(:,:,j+1), 'herdin');
    end
    yyaxis(ax,'left'); %main metric on left
    plot(ax, max_eigs(1,:), 'k', 'LineWidth', 1, 'DisplayName', '1 max |eigs|');
    hold(ax,'on');
    plot(ax, max_eigs(2,:), 'k:','LineWidth', 1, 'DisplayName', '2 max |eigs|');
    ylabel(ax, 'Max Eigs', 'FontSize', 15);
    set(ax, 'YColor', 'k')
    
    yyaxis(ax,'right');
    %plot(ax, [stein_diffs;mean(stein_diffs)], 'LineWidth', 1, 'DisplayName', 'stein');
    herdin_diffs_rescale = (herdin_diffs/max(herdin_diffs))*mean(max_eigs(1,:));
    x_diffs = 1.5:1:(num_windows-.5);
    plot(ax, x_diffs, herdin_diffs, 'm', 'LineWidth', 1, 'DisplayName', 'herdin');
    ylabel(ax, 'Covariance Measures/Distances', 'FontSize', 15);
    set(ax, 'YColor', 'm') % set y axis color to same as line color
    
    hold(ax, 'off');
end

xlabel(ax, 'Windows')
%ylim(ax, [0,100])
legend(ax)


%% fc scalars



%% perform plotting on given ax
%{
yyaxis(ax,'left'); %main metric on left
plot(ax, energy, 'k', 'LineWidth', 1, 'DisplayName', 'signal energy');
hold(ax, 'on');
lpr_int = sprintf('lpf [%d,%d]',low_freq_interval(1), low_freq_interval(2));
plot(ax, lpf_energy, 'k-', 'LineWidth', 1, 'DisplayName', lpr_int);
mpr_int = sprintf('mpf [%d,%d]', med_freq_interval(1), med_freq_interval(2));
plot(ax, mpf_energy, 'k--', 'LineWidth', 1, 'DisplayName', mpr_int);
hpr_int = sprintf('hpf [%d,%d]', high_freq_interval(1), high_freq_interval(2));
plot(ax, hpf_energy, 'k:', 'LineWidth', 1, 'DisplayName', hpr_int);

ylabel(ax, 'energy in freq range', 'FontSize', 15);
set(ax, 'YColor', 'k') % set y axis color to same as line color
%}


%{
yyaxis(ax,'right');
plot(ax, max_eigs(1,:), 'c', 'LineWidth', 1, 'DisplayName', '1 max |eigs|');
plot(ax, max_eigs(2,:), 'c--o','LineWidth', 1, 'DisplayName', '2 max |eigs|');
%plot(ax, apply_to_tensor_slices(@norm, max_eigs), 'LineWidth', 1, 'DisplayName', '2 max |eigs|');
%plot(ax, abs(max_eigs), '--y', 'LineWidth', 1, 'DisplayName', 'max |eigs|');
%plot(ax, [stein_diffs;mean(stein_diffs)], 'LineWidth', 1, 'DisplayName', 'stein');
herdin_diffs_rescale = (herdin_diffs/max(herdin_diffs))*mean(max_eigs(1,:));
plot(ax, [herdin_diffs_rescale;mean(herdin_diffs_rescale)], 'c:', 'LineWidth', 1, 'DisplayName', 'herdin');

ylabel(ax, 'Covariance Measures/Distances', 'FontSize', 15);
set(ax, 'YColor', 'c') % set y axis color to same as line color

hold(ax, 'off');
xlabel(ax, 'Window')
%ylim(ax, [0,100])
legend(ax)
%}

end
