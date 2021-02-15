function plot_signal_energy(ax, x, x_freq, intervals, interval_labels, line_colors, normalize, optional)
%plot energies of signals broken down by frequency ranges in intervals
%(indices of eigenvalues)

%%Inputs
% ax - axis to plot on
% x  - time signals
% x_freq -
% intervals = [cell] of 1x2 [int] defining the freq ranges to use
% interval_labels - [cell] of strings label for legend
% optional = struct with optional parameters such as
%            -cutoff
%            -indices of removed samples

%% compute total energy at each time point, as well as energy in each freq range
energy          = vecnorm(x_freq,2).^2;
energy_in_range = energy_in_freq_intervals(x_freq, intervals);


%use these handles to add to legend. Don't want to add all data
%points added on x axis (to denote removed data) to on legend. Only one.
if isempty(optional.remove_idxs)
    legend_handles = gobjects(1,length(intervals) + 1); %total energy + energy in each range
else
    legend_handles = gobjects(1,length(intervals) + 2); %^+data points
end

if normalize
    yyaxis(ax,'left');
end
legend_handles(1) = plot(ax, energy, 'k', 'LineWidth', 1, 'DisplayName', 'signal energy');
ylabel(ax, 'Energy', 'FontSize', 12);
hold(ax, 'on');

%show normalized distribution of energies on right y axis  
if normalize
    yyaxis(ax,'right');
    energy_in_range = energy_in_range./energy;
end

%% plot energy in each frequency range as separate lines
for j = 1:length(intervals)
	interval = intervals{j};
	line_color = line_colors{j};
	interval_label = sprintf('%s [%d,%d]',interval_labels{j}, interval(1), interval(2));
    legend_handles(j+1) = plot(ax, energy_in_range(j,:), line_color, 'LineWidth', 1, 'DisplayName', interval_label);
end


%% if a raw threshold was used, place line showing the threshold value
if normalize
    optional.raw_cutoff = optional.raw_cutoff/100;
end
if optional.plot_cutoff_line || normalize
    yline(ax, optional.raw_cutoff, 'g', 'LineWidth', 2);
end


%find where values are the smallest so as to place text in convient
%location to not block data
filt = [1,1,1,1,1,1,1];
smoothed_energies = conv(energy, filt./length(filt), 'same'); %convolve with ave filter
[~, x_text] = min(smoothed_energies);
[y_text, ~] = max(smoothed_energies);

%place message
text(ax, x_text, .9*y_text, optional.message, 'FontSize', 12, 'Color', 'g');
%text(ax, .3, .9, optional.message, 'FontSize', 12, 'Color', 'g');


%% highlight xticks that could be filtered out
for j = 1:length(optional.remove_idxs)
    removed_idx = optional.remove_idxs(j);
    h = plot(ax, removed_idx,0,'bx','MarkerSize',12, 'DisplayName', ...
        'filtered point');
end

%if no data points were removed then h will be undefined
if ~isempty(optional.remove_idxs)
    legend_handles(end) = h;
end

%% misc plotting
%place xticks on outside
%xlabel(ax, 'time','FontSize', 15);
if normalize
    ylabel(ax, 'Fraction Energy in Range', 'FontSize', 12);
    set(ax, 'YColor', 'r')
end
legend(ax, legend_handles);
hold(ax, 'off');

end

