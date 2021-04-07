%% distribution of functions of FCs
clear; close all;
load("correlations_colormap.mat");
subject_list = str2num(load('data_accounting/hcp_1200_subject_list.mat').hcp1200_subject_list); % 1113x1 int array
atlas = 'desikan';
task = "REST1";
include_subcortical = false;
if ~include_subcortical
    which_nodes = "Cortical";
else
    which_nodes = "Cortical+Sub-Cortical";
end

%% Real FCs
%[fcs, subject_ids, scan_dirs] = f(subject_list, atlas, task, include_subcortical);
%save('viz_tools/viz_fcs/fc_metric_distrib_data/fc_metrics.mat', 'fcs', 'subject_ids', 'scan_dirs');
load('viz_tools/viz_fcs/fc_metric_distrib_data/fc_metrics.mat');
[N, ~] = size(fcs{1});
fcs_tensor = zeros(N, N, length(fcs));
for idx = 1:length(fcs)
    fcs_tensor(:, :, idx) = fcs{idx};
end
[frob, max_svs, max_abs_val] = compute_fc_metrics(fcs_tensor);
[real_ms.name, real_ms.frobs, real_ms.max_svs, real_ms.max_abs_vals] = ...
    deal('Real', frob, max_svs, max_abs_val);


%% PS FCs
scs_upper_tri = load_processed_scs(atlas, include_subcortical, 'log');
scs = apply_to_tensor_slices(@(x) x+x', scs_upper_tri);
c1 = [.5, .5, .2];
c2 = [.01, .2, -.1];
c3 = [.49, .0163, -1.3061e-04]; %mid pass
coeffs_list = {c1, c2, c3};

[metrics_list, fields, names] = compute_all_ps_metrics(coeffs_list, scs);

%% add on Real fcs to end
metrics_list{end+1} = real_ms;

%% Plot on same axis on log scale
fig = figure();
t = tiledlayout(1, length(fields));
title(t, 'Distributions of FC metrics. Real & PS (H^2). Desikan w/ only Cortical Nodes.', 'FontSize', 30);
field_axes = gobjects(length(fields),1);

for field_idx = 1:length(fields)
    ax = nexttile();
    hold on;
    for m_idx = 1:length(metrics_list)
        metrics = metrics_list{m_idx};
        x = log10(metrics.(fields(field_idx)));
        histogram(ax, x, 'DisplayName', metrics.name);
    end
    legend;
    field_axes(field_idx) = ax;
    xlabel(ax, 'Log 10 Scale', 'FontSize', 25);
    title(names(field_idx), 'FontSize', 25);
    hold off;
end
ylabel(field_axes(1), 'Counts', 'FontSize', 25);

%% Plots for FCs ordered according to some metric
%{
name = names(2);
val = vals{2};

fig = figure();
num_per_group = 5;
t = tiledlayout(3,num_per_group);
title(t, sprintf('FCs across Range of %s', name),'FontSize', 30);

[sorted_vals, sort_idxs] = sort(val);
subject_ids = sort_cell_by(subject_ids, sort_idxs);
scan_dirs = sort_cell_by(scan_dirs, sort_idxs);
fcs = sort_cell_by(fcs, sort_idxs);

percentile = 95;
fs = 25;
% low 
axes_list = place_fcs(t, fcs, sorted_vals, subject_ids, scan_dirs, 1, num_per_group, colorMap, percentile, name);
ylabel(axes_list(1), 'Lowest', 'FontSize', fs);

% median
start = idivide(int16(length(fcs)), 2)-idivide(int16(num_per_group), 2);
axes_list = place_fcs(t, fcs, sorted_vals, subject_ids, scan_dirs, start, start+num_per_group-1, colorMap, percentile, name);
ylabel(axes_list(1), 'Median', 'FontSize', fs);

% high
axes_list = place_fcs(t, fcs, sorted_vals, subject_ids, scan_dirs, length(fcs)-num_per_group+1, length(fcs), colorMap, percentile, name);
ylabel(axes_list(1), 'Highest', 'FontSize', fs);
%}


function [metrics_list, fields, names] = compute_all_ps_metrics(coeffs_list, scs)

    metrics_list = cell(length(coeffs_list), 1);
    for c_idx = 1:length(coeffs_list)
        coeffs = coeffs_list{c_idx};
        [ms.frobs, ms.max_svs, ms.max_abs_vals] = compute_ps_metrics(coeffs, scs);
        ms.name = join(string(coeffs));
        metrics_list{c_idx} = ms;
    end
    
    fields = ["frobs", "max_svs", "max_abs_vals"];
    names = ["Frob Norm", "Max Singular Value", "Max Absolute Value Entrywise"];

end

function [frob, max_svs, max_abs_val] = compute_ps_metrics(coeffs, scs)
    % H^2 = ensamble cov of diffused white signals
    ps_c1_fcs_tensor = apply_to_tensor_slices(@(x) diffusion_filter(coeffs, x)^2, scs);
    [frob, max_svs, max_abs_val] = compute_fc_metrics(ps_c1_fcs_tensor);

end

function [H] = diffusion_filter(cs, S)
    H = zeros(size(S));
    for l = 1:length(cs)
        H = H + cs(l)*S^l; 
    end

end

function [frob, max_svs, max_abs_val] = compute_fc_metrics(fcs_tensor)

    frob = apply_to_tensor_slices(@(x) norm(x,'fro'), fcs_tensor);
    max_svs = apply_to_tensor_slices(@(x) norm(x), fcs_tensor);
    max_abs_val = apply_to_tensor_slices(@(x) max(abs(x), [], 'all'), fcs_tensor);
    %max_val = apply_to_tensor_slices(@(x) max(x, [], 'all'), fcs_tensor);
    %min_val = apply_to_tensor_slices(@(x) min(x, [], 'all'), fcs_tensor);

end

% assume fcs already in sorted order
function fc_axes = place_fcs(t, fcs, fc_metrics, subject_ids, scan_dirs, start_idx, end_idx, colormap, percentile, name)

fc_axes = gobjects(end_idx-start_idx+1,1);

for idx = start_idx:end_idx
    ax = nexttile(t);
    imagesc(ax, fcs{idx});
    daspect(ax,[1 1 1]);
    xlabel(sprintf('%s: %.2f', name, fc_metrics(idx)), 'FontSize', 15);
    title(sprintf('%d-%s', subject_ids{idx}, scan_dirs{idx}), 'FontSize', 15);
    fc_axes(idx-start_idx+1) = ax;
end

clim = get_cov_clim(fcs, start_idx:end_idx, percentile);
set(fc_axes,'CLim',[-clim clim]);
set(fc_axes, 'ColorMap', colormap);
cb = colorbar(fc_axes(end), 'EastOutside');
set(cb, 'FontSize', 15);
set(fc_axes, 'YtickLabel', []);
set(fc_axes, 'XtickLabel', []);




end


function lim = get_cov_clim(fcs, index_range, percentile)

[N, ~] = size(fcs{1});
tensor = zeros(N, N, length(index_range));
for idx = 1:length(index_range)
    tensor(:, :, idx) = fcs{index_range(idx)};
end

[lim] = find_cov_clim(tensor, percentile);

end

% include diagonal elements
function v = get_upper_tri_as_vec(x)

    m  = (1:size(x,1)) >= (1:size(x,2)).';
    v  = At(m);

end

function c_sorted = sort_cell_by(c, sorted_indices)
    c_sorted = cell(size(c));
    for l = 1:length(sorted_indices)
        c_sorted{l} = c{sorted_indices(l)};
    end
end

function place_summary_points(ax, x, markertype)
    y = -15;
    fs = 15;
    scatter(ax, [max(x), min(x), mean(x), median(x)], [0,0,0,0], markertype, 'DisplayName', '');
    text(max(x), y, 'max', 'FontSize', fs);
    text(min(x), y, 'min', 'FontSize', fs);
    text(mean(x), y, 'mean', 'FontSize', fs);
    text(median(x),y, 'median', 'FontSize', fs);

end


function [fcs, subject_ids, scan_dirs] = f(subject_list, atlas, task, include_subcortical)
%% go through all scans and find the variations and energies wrt to GSO
    scan_directions = ["LR", "RL"];

    fcs = {};
    subject_ids = {};
    scan_dirs = {};
    [total_scans_processed, scans_rejected] = deal(0,0);
    [missing_scs, missing_fmris, missing_atlas] = deal(0,0,0);

    
    DEBUG = false;

    for idx = 1:length(subject_list)
        
        fprintf('%d/%d\n', idx, length(subject_list));
        
        for scan_dir_idx = 1:length(scan_directions)
            scan_direction = scan_directions(scan_dir_idx);
            scan_info = ScanInfo(subject_list(idx), atlas, task, scan_direction, include_subcortical);
            try 
                fcs{end+1} = scan_info.compute_fc();
                subject_ids{end+1} = subject_list(idx);
                scan_dirs{end+1} = scan_directions(scan_dir_idx);
                total_scans_processed = total_scans_processed + 1;
            catch ME
                if ~contains(ME.identifier,'DoesNotExist')
                    rethrow(ME);
                elseif contains(ME.identifier,'DoesNotExist:Atlas')
                    missing_atlas = missing_atlas+1;
                elseif contains(ME.identifier,'DoesNotExist:fMRI')
                    missing_fmris = missing_fmris+1;
                elseif contains(ME.identifier,'DoesNotExist:SC')
                    missing_scs = missing_scs + 1;
                end
                if DEBUG
                    fprintf('%s: %s\n', ME.message, ME.identifier);
                end
                scans_rejected = scans_rejected + 1;
                continue
            end
        end

    end
    
	fprintf('Summary: %d total_scans_processed\n', total_scans_processed);
	fprintf('%d scans_rejected: %d atlas | %d fmris | %d scs\n', scans_rejected, missing_atlas, missing_fmris, missing_scs/2);

end

function scs_upper_tri = load_processed_scs(atlas, include_subcortical, edge_transform)
    cortical_idxs = get_roi_idxs(atlas, include_subcortical);

    if isequal(atlas, 'desikan')
        scs_file = load('scs_desikan.mat');
    else
        error('atlas %s unsupported', atlas)
    end
    
    scs_upper_tri = scs_file.scs(cortical_idxs, cortical_idxs, :);

    if isequal(edge_transform, 'log')
        scs_upper_tri = log(scs_upper_tri+1);
    end


end


