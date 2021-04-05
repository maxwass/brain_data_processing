%% distribution of functions of FCs
clear; close;
subject_list = str2num(load('data_accounting/hcp_1200_subject_list.mat').hcp1200_subject_list); % 1113x1 int array
atlas = 'desikan';
task = "REST1";
include_subcortical = false;
if ~include_subcortical
    which_nodes = "Cortical";
else
    which_nodes = "Cortical+Sub-Cortical";
end

%fcs = f(subject_list, atlas, task, include_subcortical);
load('viz_tools/viz_fcs/fc_metric_distrib_data/fcs.mat');
[N, ~] = size(fcs{1});
fcs_tensor = zeros(N, N, length(fcs));
for idx = 1:length(fcs)
    fcs_tensor(:, :, idx) = fcs{idx};
end


frob = apply_to_tensor_slices(@(x) norm(x,'fro'), fcs_tensor);
max_sv = apply_to_tensor_slices(@(x) norm(x), fcs_tensor);
max_abs_val = apply_to_tensor_slices(@(x) max(abs(x), [], 'all'), fcs_tensor);
max_val = apply_to_tensor_slices(@(x) max(x, [], 'all'), fcs_tensor);
min_val = apply_to_tensor_slices(@(x) min(x, [], 'all'), fcs_tensor);

vals = {frob, max_abs_val};%, max_sv, max_val, min_val};
names = ["Frob Norm", "Max Abs Value Entrywise"];%, "Max Singular Value", "Max Value Entrywise", "Min Abs Value Entrywise"];


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

function place_summary_points(ax, x)
    y = -15;
    fs = 20;
    scatter(ax, [max(x), min(x), mean(x), median(x)], [0,0,0,0], '*');
    text(max(x), y,'max', 'FontSize', fs);
    text(min(x), y,'min', 'FontSize', fs);
    text(mean(x), y,'mean', 'FontSize', fs);
    text(median(x),y,'median', 'FontSize', fs);

end


function [fcs] = f(subject_list, atlas, task, include_subcortical)
%% go through all scans and find the variations and energies wrt to GSO
    scan_directions = ["LR", "RL"];

    fcs = {};
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

