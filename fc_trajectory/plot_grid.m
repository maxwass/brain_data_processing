function plot_grid(app)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sz = [int64(app.NumWindowPlotSpinner.Value), app.cols];
if app.tiled_axes==-1
    app.tile_plot = tiledlayout(app.Panel, int64(app.NumWindowPlotSpinner.Value), app.cols,'TileSpacing','compact','Padding','none');
    app.tiled_axes = gobjects(int64(app.NumWindowPlotSpinner.Value), app.cols);
    for row = 1:int64(app.NumWindowPlotSpinner.Value)
        for col = 1:app.cols
            lin_idx = (row-1)*app.cols + col;
            app.tiled_axes(row,col)= nexttile(app.tile_plot);
        end
    end
end



corr_axes = zeros(int64(app.NumWindowPlotSpinner.Value),1);
row = 0;
col = 1;
for i = app.low_index:app.high_index
	row = row + 1;
	ax = app.tiled_axes(row,col);%nexttile(app.tile_plot, idx);
	imagesc(ax, app.corrs(:,:,i));
	set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
    yl = sprintf('%d',i);
    ylh = ylabel(ax, yl, 'FontSize', 20);
    ylp = get(ylh, 'Position');
    ext = get(ylh,'Extent');
    set(ylh, 'Rotation',0, 'Position',ylp-[ext(3) 0 0])
    %set(ylh,'rotation',0,'VerticalAlignment','middle')
	colormap(ax, app.colormap_file.colorMap);
	%drawnow;
	corr_axes(row)=ax;
end
title(corr_axes(1), "Correlations",'FontSize', 15);
linkaxes(corr_axes,'xy');
           

k = 3;
[diffs,left_index] = maxk(diff(app.eigvals),k);
row = 0;
col = 2;

diff_ = diff(app.eigvals)<.005;
start_chunk_label = find(diff_,1,'first');
end_chunk_label   = find(diff_,1,'last');
middle_chunk_label = ceil((end_chunk_label-start_chunk_label)/2);

eig_labels = {'0'};
eig_labels(1)=[];
sig_digs = 4;
%s = num2str(app.eigvals(1));
%last_idx = min(length(s),sig_digs);
%eig_labels = [s(1:last_idx)];
for i = 1:length(app.eigvals)
    e = app.eigvals(i);
    s = num2str(e); % '%0.3940'
    last_idx = min(length(s),sig_digs);
    s = s(1:last_idx);
    %{
    if s(1)=='1'
        s = s(1:last_idx); % '1.45'
    else
        s = s(2:last_idx+1); % '.39'
    end
    if (i > start_chunk_label) && (i < end_chunk_label) && (i~=middle_chunk_label)
        s = '';
    end
    ls = length(s);
    for miss_chars = sig_digs:-1:(ls+1)
        s = strcat(s,'0');
    end
    %}
    eig_labels{i}=s;
end
    
%shift eigvals over
plot_eigs = app.eigvals;
space = .03;
for i = k:-1:1
    plot_eigs(left_index(i)+1:end) = plot_eigs(left_index(i)+1:end) - (diffs(i) - space);
end


freq_axes = zeros(int64(app.NumWindowPlotSpinner.Value),1);
for i = app.low_index:app.high_index
	row = row + 1;
	ax = app.tiled_axes(row,col);
    %yyaxis(ax,'right');
	freq_plot = stem(ax, plot_eigs, app.aw_freq_normed(:,i), 'MarkerEdgeColor','green', 'color', 'k');
	
    set(ax,'xtick',[]);
    set(ax,'xticklabel',[]);
    set(ax,'yticklabel',[]);
    %xticklabels(ax, eig_labels)
    %xtickangle(ax, 45)
    %place vertical line to denote skip
    if row == int64(app.NumWindowPlotSpinner.Value)
        for j = k:-1:1
            %xline(ax, plot_eigs(left_index(j))+space/2);
            ax_offset = .15;
        
            x0 = plot_eigs(left_index(j))+space/4;
            x1 = plot_eigs(left_index(j))+3*space/4;
            line(ax, [x0, x0], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
            line(ax, [x1, x1], [app.y_min-ax_offset, app.y_min+ax_offset], 'Color','red');
            gap = num2str(round(diffs(j),2));
            gap = gap(2:end); 
            line(ax, [x0, x1], [app.y_min + ax_offset/2, app.y_min + ax_offset/2], 'Color','blue');
            text(ax, (x0+x1)/2,app.y_min+(ax_offset/2)+0.03,gap, 'FontSize',10);
        end
    end
    
    
    %yyaxis(ax,'right');
    %ylabel(ax,'|freq signal|', 'FontSize',15, 'Color', 'k');
    ylim(ax,[app.y_min,app.y_max])
    xlim(ax,[-.01, max(plot_eigs)]);
	%xlabel(ax,'freq');
    
	freq_axes(row)=ax;
end

freq_title = sprintf("GFT of NORM'ed (x_wind-total_mean) on %s", app.GSOButtonGroup.SelectedObject.Text);
title(freq_axes(1), freq_title,'FontSize', 15);
xticks(freq_axes(end), plot_eigs)
xticklabels(freq_axes(end), eig_labels)
xtickangle(freq_axes(end), 45)

linkaxes(freq_axes,'xy');

for i = 1:int64(app.NumWindowPlotSpinner.Value)
    c_ax = corr_axes(i);
    daspect(c_ax,[1 1 1]);
    pbaspect(c_ax,[1 1 1]); 
    
    f_ax = freq_axes(i);
    %daspect(f_ax,[1 1 1]);
    %pbaspect(f_ax,[6 2 1]); 
end

cbar = colorbar(corr_axes(int64(app.NumWindowPlotSpinner.Value)),'SouthOutside');
cbar.Limits = [-1,1];

drawnow;
end

