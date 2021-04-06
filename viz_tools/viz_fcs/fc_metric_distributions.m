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

%[fcs, subject_ids, scan_dirs] = f(subject_list, atlas, task, include_subcortical);
%save('viz_tools/viz_fcs/fc_metric_distrib_data/fc_metrics.mat', 'fcs', 'subject_ids', 'scan_dirs');
load('viz_tools/viz_fcs/fc_metric_distrib_data/fc_metrics.mat');
[N, ~] = size(fcs{1});
fcs_tensor = zeros(N, N, length(fcs));
for idx = 1:length(fcs)
    fcs_tensor(:, :, idx) = fcs{idx};
end

%fcs = cellfun(@corrcov, fcs, 'UniformOutput', false);

frob = apply_to_tensor_slices(@(x) norm(x,'fro'), fcs_tensor);
max_sv = apply_to_tensor_slices(@(x) norm(x), fcs_tensor);
max_abs_val = apply_to_tensor_slices(@(x) max(abs(x), [], 'all'), fcs_tensor);
max_val = apply_to_tensor_slices(@(x) max(x, [], 'all'), fcs_tensor);
min_val = apply_to_tensor_slices(@(x) min(x, [], 'all'), fcs_tensor);

vals = {frob, max_abs_val};%, max_sv, max_val, min_val};
names = ["Frob Norm", "Max Abs Value Entrywise"];%, "Max Singular Value", "Max Value Entrywise", "Min Abs Value Entrywise"];

%{
fig = figure();
t = tiledlayout(length(vals),1);
title(t, 'Distributions of FC metrics. Desikan w/ only Cortical Nodes.', 'FontSize', 30);
for idx = 1:length(vals)
    ax = nexttile();
    x = vals{idx};
    name = names(idx);
    histogram(ax, x);
    hold on;
    place_summary_points(ax, x)
    yyaxis right;
    [~, stats] = cdfplot(x);
    %txt = sprintf("%s. Max %.2f, Min %.2f, Mean %.2f, Median %.2f", ...
    %    name, max(x), min(x), mean(x), median(x));
    
    title(name, 'FontSize', 25);
end
%}

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

function place_summary_points(ax, x)
    y = -15;
    fs = 20;
    scatter(ax, [max(x), min(x), mean(x), median(x)], [0,0,0,0], '*');
    text(max(x), y,'max', 'FontSize', fs);
    text(min(x), y,'min', 'FontSize', fs);
    text(mean(x), y,'mean', 'FontSize', fs);
    text(median(x),y,'median', 'FontSize', fs);

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

