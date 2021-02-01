function fc_traj_plot_callback(obj, event, app)
%Callback for fc_trajectory app 
%update plots

if strcmp(app.PlaySwitch.Value,'Play')
    %update indices for next time
    jumpsize = int64(app.JumpSpinner.Value);
    % num_windows is the total # windows in fc trajectory, nothing to do
    % with plotting
    app.low_index     = mod(app.low_index + jumpsize, app.num_windows);
    app.high_index    = mod(app.low_index + app.num_windows_plot-1, app.num_windows);
    plot_grid(app);
    
end

%{
rows=app.rows;
subject = app.SubjectIDsListBox.Value;
atlas   = app.AtlasButtonGroup.SelectedObject.Text;
scan_dir    = app.SelectScanButtonGroup.SelectedObject.Text;
include_subcortical_text = app.BrainAreasButtonGroup.SelectedObject.Text;
windowsize = int64(app.windowsizeSpinner.Value);
movesize = int64(app.movesizeSpinner.Value);
jumpsize = int64(app.JumpSpinner.Value);
num_windows_total = app.num_windows;


%app.tile_plot = tiledlayout(app.Panel, rows,app.num_windows_plot,'TileSpacing','none','Padding','none');

if include_subcortical_text
    txt = sprintf('patient %s %s1. Windowsize %d, Movesize %d || %s || Cortical+Subcortical', subject, scan_dir, windowsize, movesize, atlas);
else
    txt = sprintf('patient %s %s1. Windowsize %d, Movesize %d || %s || Only Cortical', subject, scan_dir, windowsize, movesize, atlas);
end

%title(app.tile_plot,txt, 'FontSize', 15);
            
app.low_index = mod(app.low_index+jumpsize, num_windows_total);
high_index    = mod(app.low_index+app.num_windows_plot-1,num_windows_total);
            
corr_axes = zeros(app.num_windows_plot);
idx=-1;
for i = 1:app.num_windows_plot
    idx = app.low_index+i;
	%ax = nexttile(app.tile_plot);
	ax = app.tiled_axes(1,i);
    imagesc(ax, app.corrs(:,:,idx));
	set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
    colormap(ax, app.colormap_file.colorMap);
    ttl = sprintf('Window %d',idx);
    title(ax, ttl);
	corr_axes(i)=ax;
end
linkaxes(corr_axes,'xy');
cbar = colorbar(corr_axes(app.num_windows_plot));
cbar.Limits = [-1,1];            
freq_axes = zeros(app.num_windows_plot);
drawnow;

idx=-1;
for i = 1:app.num_windows_plot
    idx = app.low_index+i;
    %ax = nexttile(app.tile_plot);
    ax = app.tiled_axes(2,i);
	freq_plot = stem(ax, app.eigvals, app.aw_freq_normed(:,idx), 'MarkerEdgeColor','green', 'color', 'k');
	ylim(ax,[app.y_min,app.y_max])
	xlabel(ax,'freq');
	if idx==app.low_index
        ylabel(ax,'|freq signal|');
    else
        set(ax,'yticklabel',[]);
    end
    freq_axes(i)=ax;
end
linkaxes(freq_axes,'xy');
            
axis equal;

for i = 1:app.num_windows_plot
    c_ax = corr_axes(i);
    daspect(c_ax,[1 1 1]);
    pbaspect(c_ax,[1 1 1]); 
    
    f_ax = freq_axes(i);
    daspect(f_ax,[1 1 1]);
    pbaspect(f_ax,[1 1 1]); 
end

drawnow;
end
%}

