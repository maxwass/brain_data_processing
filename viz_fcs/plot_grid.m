function plot_grid(app)

if strcmp(app.NormFreqDropDown.Value, 'L2 Norm')
    %each signal is now L2 normed
	ave_signal_windows_freq_reprs = app.ave_signal_windows_freq./vecnorm(app.ave_signal_windows_freq);
else
	ave_signal_windows_freq_reprs = app.ave_signal_windows_freq;
end

%{
dtseries_summary_fc = [];
if use_all_signals
    dtseries_summary_fc = app.plot_dtseries;
elseif app.raw_filter_rm_from_fc_comp.Value
    dtseries_summary_fc = app.plot_dtseries';
    dtseries_summary_fc(:, app.removed_dtseries_idxs) =[];
    summary_fc = cov(z');
elseif use_in_range
end
%}
dtseries_summary_fc = app.dtseries_frequency_filtered;
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
    %%create axes and store 
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

    if isequal(which_scalars, "signal_freq_distrib")
        x = app.plot_dtseries; %show *all* signals, even those filtered out
        x_freq = app.GFT*x;
        optional.cutoff         = app.raw_filter_cutoff.Value;
        optional.raw_cutoff     = app.filter_dtseries_raw_cutoff; %convert from percentile to value in computation
        optional.use_percentile = app.raw_filter_use_percentile.Value;
        optional.percentile_val = app.raw_filter_cutoff.Value;
        optional.remove_idxs    = app.removed_dtseries_idxs;

        %if raw energy threshold is used, plot line which denotes cutoff
        if isequal(which_raw_filter, 'energy')
            optional.plot_cutoff_line = true;
            optional.message = sprintf('energy > %.0f || prcntile? %s', ...
                optional.raw_cutoff, string(optional.use_percentile));
        %must plot normalized by total energy for threshold to be straight line    
        elseif isequal(which_raw_filter, 'freq_distribution')
            optional.plot_cutoff_line = false;
            optional.message = sprintf('frac energy in low freq > %.2f\n use prcntile? %s', ...
                optional.raw_cutoff, string(optional.use_percentile));
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
    end
    plot_signal_energy(app.axes.scalars, x, x_freq, ...
            app.energy_freq_plot.intervals, ...
            app.energy_freq_plot.interval_labels, ...
            app.energy_freq_plot.line_colors, ...
            app.normalize_scalar_distrib_Button.Value,...
            optional);
        
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
title(fc_axes(2), app.PlotFCs.Value,'FontSize', 12);
linkaxes(fc_axes,'xy');
           


%% Freq Representations: Column 2
%% eigenvalues can  be very spread over interval, making viz hard. 
% Find largest num_jumps_delete jumps between eigenvalues and shorten them
% to be space length. New shifted eigenvalues(plot_eigs) are for plotting
% only.
y_max_freq = max(max(ave_signal_windows_freq_reprs)); %use percentile instead?
y_min_freq = min(min(ave_signal_windows_freq_reprs));

num_jumps_delete = app.xbreaksSpinner.Value;
space = 0.03*(max(app.eigvals)-min(app.eigvals));
eig_diffs = diff(app.eigvals);
[largest_diffs,left_idx_largest_diffs] = maxk(eig_diffs,num_jumps_delete);

%sort so that we can do order sensitive operations next
[left_idx_largest_diffs, B] = sort(left_idx_largest_diffs, 'ascend');
largest_diffs = largest_diffs(B);

%shift eigvals over from largest differnce to smallest difference
% must do it this way to keep eigenvalues consistant (?)
plot_eigs = app.eigvals;
for i = num_jumps_delete:-1:1
    plot_eigs(left_idx_largest_diffs(i)+1:end) = plot_eigs(left_idx_largest_diffs(i)+1:end) - (largest_diffs(i) - space);
end

row=0;
% place plots
for i = app.low_index:app.high_index
	row = row + 1;
	ax = app.axes.wdw_features(row,freq_col);
    %yyaxis(ax,'right');
	freq_plot = stem(ax, plot_eigs, ave_signal_windows_freq_reprs(:,i), 'MarkerEdgeColor','green', 'color', 'k');
	
    set(ax,'xtick',[]);
    set(ax,'xticklabel',[]);
    if i ~=app.high_index
        set(ax,'yticklabel',[]);
    end
    
    ylim(ax,[y_min_freq,y_max_freq])
    xlim(ax,[min(plot_eigs), max(plot_eigs)]);
    str = sprintf('L2 Norm (mean centered): %.0f', norm(ave_signal_windows_freq_reprs(:,i), 2) );
    text(ax, 'String', str, 'Units', 'normalized', 'Position', [.6, .95])
	%xlabel(ax,'freq');
end
freq_axes = app.axes.wdw_features(1:end,freq_col);
linkaxes(freq_axes,'xy');

% place vertical lines denoting where jumps were made
ax_offset = 0.1*(app.y_max-app.y_min); %how far from botton of plot should line axtend?
text_offset = .1; %how far offr horrizontal line should text be?
for j = num_jumps_delete:-1:1
 
    % new x positions of removed larger interval
	x0 = plot_eigs(left_idx_largest_diffs(j))+space/4;
	x1 = plot_eigs(left_idx_largest_diffs(j))+3*space/4;
	
    % draw vertical lines to denote where jump occured
    line(freq_axes(end), [x0, x0], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
	line(freq_axes(end), [x1, x1], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
	
    %  text indication how large jump was
    gap = num2str(round(largest_diffs(j),5));
	%gap = gap(2:end); 
	line(freq_axes(end), [x0, x1], [app.y_min + ax_offset/2, app.y_min + ax_offset/2], 'Color','blue');
	text(freq_axes(end), (x0+x1)/2,app.y_min+(ax_offset/2)+text_offset,gap, 'FontSize',10);
end

%make odd plots have xticks as the idxs
xticks(freq_axes(1:2:end), plot_eigs);
xticklabels(freq_axes(1:2:end), app.eig_idxs);
set(freq_axes(1:2:end),'fontsize',7);
%xlabel(freq_axes(1:2:end), 'Eigval IDX');
xtickangle(freq_axes(1:2:end), 45)


%make even plots have the xticks as the eigenvalues
xticks(freq_axes(2:2:end), plot_eigs)
xticklabels(freq_axes(2:2:end), app.eig_labels);
set(freq_axes(1:2:end),'fontsize',7);
xtickangle(freq_axes(2:2:end), 45)

%xlabel(freq_axes(2:2:end), 'Eigvalue (freq)');



freq_title = sprintf("GFT of %s (x wind - total mean) on %s", app.NormFreqDropDown.Value,  app.GSODropDown.Value);
title(freq_axes(1), freq_title,'FontSize', 15);

end

