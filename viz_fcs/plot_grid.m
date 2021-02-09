function plot_grid(app)


%% pull out (most) inputs (stored as private variables in app)
low_index = app.low_index;
high_index = app.high_index;
num_windows_plot = int64(app.NumWindowPlotSpinner.Value);
num_features_plot = app.cols; % fcs | freq repr | __

signal_windows     = app.signal_windows;
ave_signal_windows = app.ave_signal_windows;
ave_signal_windows_freq_reprs_raw    = app.ave_signal_windows_freq_reprs_raw;
ave_signal_windows_freq_reprs_normed = app.ave_signal_windows_freq_reprs_raw./vecnorm(ave_signal_windows_freq_reprs_raw);

if strcmp(app.NormFreqReprButtonGroup.SelectedObject.Text, 'L2 Norm')
	ave_signal_windows_freq_reprs = ave_signal_windows_freq_reprs_normed;
else
	ave_signal_windows_freq_reprs = ave_signal_windows_freq_reprs_raw;
end

which_fc = app.WhichFCButtonGroup.SelectedObject.Text;
if strcmp(which_fc, 'Corr')
	fcs = app.corrs;
elseif strcmp(which_fc, 'Cov')
	fcs = app.covs;
else
	error('Unrecognized FC selection %s\n', which_fc)
end

[N, ~, ~] = size(fcs);
fc_col = 1;
freq_col = 2;


which_scalars = app.ScalarsDropDown.Value;



%% We reuse same axes (tiled layout) for faster updating. Only create new axes
% if change_flag is set externally (new patient is loaded, new
% windowsize/movesize, first time running, etc)
sz = [num_windows_plot, num_features_plot];
if app.change_axes==true
    
    app.tile_plot = tiledlayout(app.Panel, num_windows_plot + 1, num_features_plot,'TileSpacing','compact','Padding','none');
    app.tiled_axes = gobjects(num_windows_plot, num_features_plot);
    for row = 1:(num_windows_plot+1)
        for col = 1:num_features_plot
            %lin_idx = (row-1)*num_features_plot + col;
            app.tiled_axes(row,col)= nexttile(app.tile_plot);
        end
    end
    
    scalar_plot_ax      = nexttile(app.tile_plot, 1, [1 2]);
    app.tiled_axes(1,1) = scalar_plot_ax;
    
    % for scalars, always use raw (each vector not normed to unit energy)
    plot_scalars(scalar_plot_ax, ave_signal_windows, ave_signal_windows_freq_reprs_raw, fcs, which_scalars)
    text(scalar_plot_ax,  double(low_index),  -3, 'L', 'Color', 'magenta', 'FontSize', 18);
    text(scalar_plot_ax,  double(high_index), -3, 'H', 'Color', 'blue',   'FontSize', 18);
    xticks(scalar_plot_ax, 1:2:app.num_windows)
    %set(ax, 'XTick', 1:20:num_windows)
    app.change_axes=false;
end




%% place scalar summary plots in first row all they way across columns
%add marker to x axis to denote low and high idxs
scalar_plot_ax = app.tiled_axes(1,1);
scalar_plot_ch = get(scalar_plot_ax,'Children');
for idx_ch = 1:length(scalar_plot_ch)
    ch = scalar_plot_ch(idx_ch);
    if isa(ch, 'matlab.graphics.primitive.Text')
        if strcmp(ch.String, 'L')
            ch.Position = [low_index, -3, 0];
        elseif strcmp(ch.String, 'H')
            ch.Position = [high_index, -3, 0];
        end
    end
end



%% FCs: Column 1
row=1;
for i = low_index:high_index
	row = row + 1;
	ax = app.tiled_axes(row,fc_col);%nexttile(app.tile_plot, idx);
    fc = fcs(:,:,i);
	imagesc(ax, fc);
    
	
    set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
    yl = sprintf('%d',i);
    ylh = ylabel(ax, yl, 'FontSize', 20);
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation',0, 'Position',ylp-[2*ext(3) 0 0])
    if i==low_index
        set(ylh, 'Color', 'magenta');
    elseif i==high_index
        set(ylh, 'Color', 'blue');
    end
	
    
    xlim(ax,[1,length(app.eigvals)]); xlim(ax,'manual');
    ylim(ax,[1,length(app.eigvals)]); ylim(ax, 'manual');
    daspect(ax,[1 1 1]);
    %pbaspect(ax,[1 1 1]);
    
    
end
fc_axes = app.tiled_axes(2:end,fc_col);
percentile = app.fcscolorsSpinner.Value; 
[cl] = find_cov_clim(fcs, percentile); %ciel needed?
set(fc_axes, 'ColorMap', app.colormap_file.colorMap);
set(fc_axes,'CLim',[-cl cl]);
cb = colorbar(fc_axes(end),'West');
%cb.Layout.Tile = 'west';
title(fc_axes(1), which_fc,'FontSize', 15);
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

row=1;
% place plots
for i = low_index:high_index
	row = row + 1;
	ax = app.tiled_axes(row,freq_col);
    %yyaxis(ax,'right');
	freq_plot = stem(ax, plot_eigs, ave_signal_windows_freq_reprs(:,i), 'MarkerEdgeColor','green', 'color', 'k');
	
    set(ax,'xtick',[]);
    set(ax,'xticklabel',[]);
    if i ~=high_index
        set(ax,'yticklabel',[]);
    end
    
    ylim(ax,[y_min_freq,y_max_freq])
    xlim(ax,[min(plot_eigs), max(plot_eigs)]);
    str = sprintf('L2 Norm (mean centered): %.0f', norm(ave_signal_windows_freq_reprs_raw(:,i), 2) );
    text(ax, 'String', str, 'Units', 'normalized', 'Position', [.6, .95])
	%xlabel(ax,'freq');
end
freq_axes = app.tiled_axes(2:end,freq_col);
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

xticks(freq_axes(end), plot_eigs)
xticklabels(freq_axes(end), app.eig_labels)
xtickangle(freq_axes(end), 45)
freq_title = sprintf("GFT of %s (x wind - total mean) on %s", app.NormFreqReprButtonGroup.SelectedObject.Text,  app.GSOButtonGroup.SelectedObject.Text);
title(freq_axes(1), freq_title,'FontSize', 15);

end

