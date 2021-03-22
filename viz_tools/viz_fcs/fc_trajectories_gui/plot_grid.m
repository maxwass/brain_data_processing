function plot_grid(app)

%{
dtseries_summary_fc = [];
if use_all_signals
    dtseries_summary_fc = app.dtseries_mean_center_all;
elseif app.raw_filter_rm_frcom_fc_comp.Value
    dtseries_summary_fc = app.dtseries_mean_center_all';
    dtseries_summary_fc(:, app.removed_dtseries_idxs) =[];
    summary_fc = cov(z');
elseif use_in_range
end
%}
[~,num_obsvs] = size(app.dtseries_frequency_filtered);
if app.raw_filter_rm_from_fc_comp.Value
    kept_dtseries_idxs = setdiff(1:num_obsvs, app.removed_dtseries_idxs);
    dtseries_summary_fc = app.dtseries_frequency_filtered(:, kept_dtseries_idxs);
else
    dtseries_summary_fc = app.dtseries_frequency_filtered;
end

if strcmp(app.PlotFCs.Value, 'Corr')
	fcs = app.corrs;
    summary_fc = corr(dtseries_summary_fc');
elseif strcmp(app.PlotFCs.Value, 'Cov')
	fcs = app.covs;
    summary_fc = cov(dtseries_summary_fc');
else
	error('Unrecognized FC selection %s\n', app.PlotScalars.Value)
end

fc_col = 1;
freq_col = 2;


which_scalars    = app.PlotScalars.Value;
which_raw_filter = app.raw_filter.Value;
which_wdw_filter = app.wdw_filter.Value;


%% We reuse same axes (tiled layout) for faster updating. Only create new axes
% if change_flag is set externally (new patient is loaded, new
% windowsize/movesize, first time running, etc)

if app.change_axes==true
    
    app.tile_plot = tiledlayout(app.Panel, app.NumWindowPlotSpinner.Value + 1, app.cols,'TileSpacing','compact','Padding','none');
    app.axes.wdw_features = gobjects(app.NumWindowPlotSpinner.Value, app.num_features_plot); %doesn't include scalars on top
    
    %first row is scalars
    scalar_span = [1,app.cols-1];
    app.axes.scalars = nexttile(app.tile_plot, 1, scalar_span);
    app.axes.summary_fc  = nexttile(app.tile_plot);
    %% create axes and store 
    % all the rest of the rows are summaries of windows (fcs, freq repr of
    % mean vectors...)
    feature_span = [1,app.cols/2];
    for ftr_row = 1:(app.NumWindowPlotSpinner.Value)
        for feature = 1:app.num_features_plot
            %compute which tile we're on
            tile_num = ftr_row*app.cols + (app.cols/2)*(feature-1) + 1;
            
            %place a long axis here
            app.axes.wdw_features(ftr_row, feature)= nexttile(app.tile_plot,tile_num,feature_span);
        end
    end

    %% which scalar plots should we do?
    if isequal(which_scalars, "signal_freq_distrib")
        x_mean = mean(app.dtseries(app.roi_idxs,:),2);
        x      = app.dtseries(app.roi_idxs,:) - app.ave_node_val;%-x_mean; %show *all* signals, even those filtered out
        x_freq = app.GFT*x;
        optional.cutoff         = app.raw_filter_cutoff.Value;
        optional.raw_cutoff     = app.filter_dtseries_raw_cutoff; %convert from percentile to value in computation
        optional.use_percentile = app.raw_filter_use_percentile.Value;
        optional.percentile_val = app.raw_filter_cutoff.Value;
        optional.remove_idxs    = app.removed_dtseries_idxs;
        optional.cutoff_line_y_axis = 'left';
        if app.normalize_scalar_distrib_Button.Value
                optional.cutoff_line_y_axis = 'right';
        end
 
        %if raw energy threshold is used, plot line which denotes cutoff
        if isequal(which_raw_filter, 'energy')
            optional.plot_cutoff_line = true;
            optional.message = sprintf('energy > %.0f \n use prcntile?  %s', ...
                optional.raw_cutoff, string(optional.use_percentile));
        %must plot normalized by total energy for threshold to be straight line    
        elseif isequal(which_raw_filter, 'freq_distribution')
            optional.plot_cutoff_line = true;
            optional.message = sprintf('frac energy in low freq > %.2f\n use prcntile? %s', ...
                optional.raw_cutoff, string(optional.use_percentile));
        elseif isequal(which_raw_filter, 'None')
            optional.plot_cutoff_line = false;
            optional.message = sprintf('No filtering done.');
        else
            error('unrecognized filter used: %s', which_raw_filter);
        end
        [~,num_obsvs] = size(x);
        xticks_ = 1:50:num_obsvs;
        
	elseif isequal(which_scalars, "wdw_ave_energy")
        x = app.ave_signal_windows;
        x_freq = app.GFT*x;
        optional.cutoff         = app.wdw_filter.Value;
        optional.raw_cutoff     = app.filter_wdw_raw_cutoff;  %convert from percentile to value
        optional.use_percentile = app.wdw_filter_use_percentile.Value;
        optional.percentile_val = app.wdw_filter_cutoff.Value;
        optional.remove_idxs    = app.removed_window_idxs;
        optional.cutoff_line_y_axis = 'left';
        if app.normalize_scalar_distrib_Button.Value
                optional.cutoff_line_y_axis = 'right';
        end
        
        %if raw energy threshold is used, plot line which denotes cutoff
        threshold_filter = ["energy", "eig_1", "eig_2"];
        if ismember(which_wdw_filter, threshold_filter)
            optional.plot_cutoff_line = true;
            optional.message = sprintf('%s > %.0f || prcntile? %s', ...
                which_raw_filter, optional.raw_cutoff, ...
                string(optional.use_percentile));
            
        elseif isequal(which_wdw_filter, 'freq_distribution')
            optional.plot_cutoff_line = false;
            optional.message = sprintf('frac energy in low freq > %.2f %%\n use prcntile? %s', ...
                optional.raw_cutoff, string(optional.use_percentile));
        else
            error('unrecognized filter used: %s', which_wdw_filter);
        end  
        xticks_ = 1:app.num_windows;
        
    elseif isequal(which_scalars, "fc_traj")
        fc_traj = apply_to_tensor_slices(@(x) cov(x'), app.signal_windows);
        xticks_ = 1:app.num_windows;
        
    end
    
    
    if ismember(which_scalars, ["signal_freq_distrib","wdw_ave_energy"]) 
        plot_signal_energy(app.axes.scalars, x, x_freq, ...
            app.energy_freq_plot.intervals, ...
            app.energy_freq_plot.interval_labels, ...
            app.energy_freq_plot.line_colors, ...
            app.normalize_scalar_distrib_Button.Value,...
            optional);
    else
        plot_fc_metrics(app.axes.scalars, fc_traj);
    end
        
    text(app.axes.scalars,   double(app.low_index),  -3, 'L', 'Color', 'magenta', 'FontSize', 18);
    text(app.axes.scalars,   double(app.high_index), -3, 'H', 'Color', 'blue',   'FontSize', 18);
    xticks(app.axes.scalars, xticks_);
    
    app.change_axes=false;
end

%summary overall FC
ax = app.axes.summary_fc;
imagesc(ax, summary_fc);
set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
xlim(ax,[1,length(app.eigvals)]); xlim(ax,'manual');
ylim(ax,[1,length(app.eigvals)]); ylim(ax, 'manual');
daspect(ax,[1 1 1]);
txt = sprintf("%s: %s", app.PlotFCs.Value, "all");
%title(ax, txt,'FontSize', 12); %make string to indicate what used: all, filterd, specified range




%% add marker to x axis to denote low and high idxs
scalar_plot_ch = get(app.axes.scalars,'Children');
for idx_ch = 1:length(scalar_plot_ch)
    ch = scalar_plot_ch(idx_ch);
    if isa(ch, 'matlab.graphics.primitive.Text')
        if strcmp(ch.String, 'L')
            ch.Position = [app.low_index, 0, 0];
        elseif strcmp(ch.String, 'H')
            ch.Position = [app.high_index, 0, 0];
        end
    end
end



%% FCs: Column 1
row=0;
for i = app.low_index:app.high_index
	row = row + 1;
	ax = app.axes.wdw_features(row,fc_col);
    fc = fcs(:,:,i);
	imagesc(ax, fc);
    
	
    set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
    yl = sprintf('%d',i);
    ylh = ylabel(ax, yl, 'FontSize', 20);
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation',0, 'Position',ylp-[2*ext(3) 0 0])
    if i==app.low_index
        set(ylh, 'Color', 'magenta');
    elseif i==app.high_index
        set(ylh, 'Color', 'blue');
    end
	
    
    xlim(ax,[1,length(app.eigvals)]); xlim(ax,'manual');
    ylim(ax,[1,length(app.eigvals)]); ylim(ax, 'manual');
    daspect(ax,[1 1 1]);
    %pbaspect(ax,[1 1 1]);
    
    
end
fc_axes = [app.axes.summary_fc; app.axes.wdw_features(1:end,fc_col)];
percentile = app.fcscolorsSpinner.Value; 
[cl] = find_cov_clim(fcs, percentile); %ciel needed?
set(fc_axes, 'ColorMap', app.colormap_file.colorMap);
set(fc_axes,'CLim',[-cl cl]);
cb = colorbar(fc_axes(end),'West');
%cb.Layout.Tile = 'west';
%title(fc_axes(2), app.PlotFCs.Value,'FontSize', 12);
linkaxes(fc_axes,'xy');
           


%% Freq Representations: Column 2
%% eigenvalues can  be very spread over interval, making viz hard. 
% Find largest num_jumps_delete jumps between eigenvalues and shorten them
% to be space length. New shifted eigenvalues(plot_eigs) are for plotting
% only.
y_values_to_plot = 'energy'; % 'magnitude', 'raw';

centered_energies = vecnorm(app.ave_signal_windows,2).^2;

% L2 norming
%ave_signal_windows_freq = app.GFT*(app.ave_signal_windows./vecnorm(app.ave_signal_windows));
ave_signal_windows_freq = app.GFT*(app.ave_signal_windows);

if isequal(y_values_to_plot, 'raw')
    y_values = ave_signal_windows_freq;
    y_max = max(max(y_values)); %use percentile instead?
    y_min = min(min(y_values));
elseif isequal(y_values_to_plot, 'magnitude')
    y_values = abs(ave_signal_windows_freq);
    y_max = max(max(y_values));
    y_min = 0;
elseif isequal(y_values_to_plot, 'energy')
    y_values = ave_signal_windows_freq.^2;
    y_max = max(max(y_values));
    y_min = 0;
end

%% remove jumps between large gaps in x-axis (i.e. eigenvalues/total variations)
x_values_to_plot = 'Total Variations'; %'eigvals'; 'zcs';
use_eig_order = contains(app.VariationDropDown.Value, "eig", "IgnoreCase", true);
use_tv_l_order = contains(app.VariationDropDown.Value, "Lz", "IgnoreCase", true);
use_tv_a_order = contains(app.VariationDropDown.Value, "Az", "IgnoreCase", true);
use_zero_cross_order = contains(app.VariationDropDown.Value, "Zero-Crossings", "IgnoreCase", true);

if use_tv_l_order
    x_values_to_plot = 'Total Variations';
    x_values = app.total_variations;
elseif use_eig_order
    x_values_to_plot = 'Eigenvalues';
    x_values = app.eigvals;
elseif use_tv_a_order
    error('not implimented');
elseif use_zero_cross_order
    x_values_to_plot = 'Zero-Crossings';
    x_values = app.zero_crossings;
end
[x_labels, ~] = create_labels(x_values, app.sig_digs, app.eig_label_sparse_middle_jumps);

num_jumps_delete = app.xbreaksSpinner.Value;
space = 0.03*(max(x_values)-min(x_values));
x_diffs = diff(x_values);
[largest_diffs,left_idx_largest_diffs] = maxk(x_diffs,num_jumps_delete);

%sort so that we can do order sensitive operations next
[left_idx_largest_diffs, B] = sort(left_idx_largest_diffs, 'ascend');
largest_diffs = largest_diffs(B);

%shift eigvals over from largest differnce to smallest difference
% must do it this way to keep eigenvalues consistant (?)
for i = num_jumps_delete:-1:1
    x_values(left_idx_largest_diffs(i)+1:end) = x_values(left_idx_largest_diffs(i)+1:end) - (largest_diffs(i) - space);
end

%% place plots
row=0;
for i = app.low_index:app.high_index
	row = row + 1;
	ax = app.axes.wdw_features(row,freq_col);
	%freq_plot = stem(ax, plot_eigs, signal_f_i, 'MarkerEdgeColor','green', 'color', 'k');
	tv_plot = stem(ax, x_values, y_values(:,i), 'MarkerEdgeColor','green', 'color', 'k');
    %set(ax,'xtick',[]);
    %set(ax,'xticklabel',[]);
    if i ~=app.high_index
        set(ax,'yticklabel',[]);
    end
    
    ylim(ax,[y_min,y_max])
    xlim(ax,[min(x_values)-.1, max(x_values)]);
    str = sprintf('Filtered Energy (mean centered): %.0f', centered_energies(i) );
    text(ax, 'String', str, 'Units', 'normalized', 'Position', [.6, .1])
	%xlabel(ax,'freq');
    ylabel(ax, y_values_to_plot);
end
freq_axes = app.axes.wdw_features(1:end,freq_col);
linkaxes(freq_axes,'xy');

% place vertical lines denoting where jumps were made
ax_offset = 0.1*(app.y_max-app.y_min); %how far from botton of plot should line axtend?
text_offset = .1; %how far offr horrizontal line should text be?
for j = num_jumps_delete:-1:1
 
    % new x positions of removed larger interval
	x0 = x_values(left_idx_largest_diffs(j))+space/4;
	x1 = x_values(left_idx_largest_diffs(j))+3*space/4;
	
    % draw vertical lines to denote where jump occured
    line(freq_axes(end), [x0, x0], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
	line(freq_axes(end), [x1, x1], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
	
    %  text indication how large jump was
    gap = num2str(round(largest_diffs(j),5));
	%gap = gap(2:end); 
	line(freq_axes(end), [x0, x1], [app.y_min + ax_offset/2, app.y_min + ax_offset/2], 'Color','blue');
	text(freq_axes(end), (x0+x1)/2,app.y_min+(ax_offset/2)+text_offset,gap, 'FontSize',10);
end


%%IF ONLY 1 PLOT THEN NO EVENS AND ODDS!

%make odd plots have xticks as the idxs

[x_values_sort, x_values_sort_idx] = sort(x_values); %must be in increasing sorted order
[x_labels_sort] = x_labels(x_values_sort_idx);
if length(unique(x_values_sort)) < length(x_values)
    %dont adjust x-axis...  repeated elements -> Zero-Crossings most likely
else
    xticks(freq_axes, x_values_sort); 
    xticklabels(freq_axes, x_labels_sort);
end
set(freq_axes,'fontsize',7);
%xlabel(freq_axes(1:2:end), 'Eigval IDX');
xtickangle(freq_axes, 45)
xlabel(freq_axes(end), x_values_to_plot);
ylabel(freq_axes(end), y_values_to_plot);
    
    
%{
xticks(freq_axes(1:2:end), plot_eigs);
xticklabels(freq_axes(1:2:end), app.eig_idxs);
set(freq_axes(1:2:end),'fontsize',7);
%xlabel(freq_axes(1:2:end), 'Eigval IDX');
xtickangle(freq_axes(1:2:end), 45)


%if app.NumWindowPlotSpinner.Value>1
xticks(freq_axes, plot_eigs)
xticklabels(freq_axes, app.eig_labels);
set(freq_axes,'fontsize',7);
xtickangle(freq_axes, 45)
%xlabel(freq_axes(2:2:end), 'Eigvalue (freq)');
%}

%{
%make even plots have the xticks as the eigenvalues
xticks(freq_axes(2:2:end), plot_eigs)
xticklabels(freq_axes(2:2:end), app.eig_labels);
set(freq_axes(2:2:end),'fontsize',7);
xtickangle(freq_axes(2:2:end), 45)
%xlabel(freq_axes(2:2:end), 'Eigvalue (freq)');
%}



if app.raw_filter_rm_from_fc_comp.Value
    mean_over = 'filtered';
else
    mean_over = 'all'; 
end
freq_txt = sprintf('$GFT_{%s}( \\bar{x}_{window(%d)} - \\bar{x}_{%s} )$', app.GSODropDown.Value, app.low_index, mean_over);
text(freq_axes(1), 'String', freq_txt, 'Units', 'normalized', 'Position', [.55, .9], 'FontSize', 15,'Interpreter','Latex');

%title(freq_axes(1), freq_txt,'FontSize', 15,'Interpreter','Latex');

end

