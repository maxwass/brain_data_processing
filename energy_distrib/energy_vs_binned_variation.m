%% Energy vs Binned Variations (TV/ZCs/Eigvals)

subject_list = str2num(load('data/hcp_1200_subject_list.mat').hcp1200_subject_list); % 1113x1 int array
atlas = 'desikan';
task = "REST1";
include_subcortical = false;


GSOs = ["L", "L_norm", "A", "A_norm"];
variation_metrics = ["total_variation_L", "zero_crossings", "eigenvalues"];

for row = 1:length(GSOs)
    GSO = GSOs(row);
    
    for var_metric_idx = 1:length(variation_metrics)
        var_metric = variation_metrics(var_metric_idx);
        filename = strcat('GSO-', GSO, '-VariatonMetric-', var_metric, '.mat');
        filepath = fullfile('energy_distrib/ed_data_new', filename);
        if isfile(filepath)
            variations_energies = load(filepath);
            [variations, energies] = deal(variations_energies.variations, variations_energies.energies);
        else
            [variations, energies] = f(subject_list, atlas, task, include_subcortical, GSO, var_metric);
            % save
            save(filepath,'variations', 'energies');
        end
        fprintf('GSO %s, variation metric %s\n', GSO, var_metric);
        min_var    = min(min(variations));
        max_var    = max(max(variations));
        min_energy = min(min(energies));
        max_energy = max(max(energies));
        fprintf('      min var: %f, max var: %f, min energy %f, max energy: %f\n\n', min_var, max_var, min_energy, max_energy);
        %min_max = bin_struct.(GSO).(var_metric);
        %bins = linspace(min_max.min, min_max.max, num_bins);
        %place histogram
        
    end
end


% Each pair of (GSO, var_metric) also has a range for the bins
bin_struct = struct("L", [], "L_norm", [], "A", [], "A_norm", []);
bin_struct.L      = struct("total_variation_L", mms(0, 421), "zero_crossings", mms(0, 760),   "eigenvalues", mms(0,421));
bin_struct.L_norm = struct("total_variation_L", mms(0, 270), "zero_crossings", mms(0, 794.1), "eigenvalues", mms(0, 1.26));
bin_struct.A      = struct("total_variation_L", mms(0, 256), "zero_crossings", mms(0, 794.1), "eigenvalues", mms(-34,204));
bin_struct.A_norm = struct("total_variation_L", mms(0, 270), "zero_crossings", mms(0, 794),   "eigenvalues", mms(-.7,1.0));

num_bins = 200;
t = tiledlayout(length(GSOs), length(variation_metrics));
axes_grid = gobjects(length(GSOs), length(variation_metrics));

for row = 1:length(GSOs)
    GSO = GSOs(row);
    tilenum_start = (row-1)*length(variation_metrics);
    axes_grid(row,:) = create_GSO_hist(tilenum_start, GSO, variation_metrics);
    ylabel(axes_grid(row,1), GSO, 'FontSize', 25);
end

for col = 1:length(variation_metrics)
    title(axes_grid(1,col),    variation_metrics(col), 'FontSize', 25);
    %xlabel(axes_grid(end,col), 
end


function [axes_row] =  create_GSO_hist(tilenum_start, GSO, variation_metrics)
    
    axes_row = gobjects(length(variation_metrics),1);
    for var_metric_idx = 1:length(variation_metrics)
        var_metric = variation_metrics(var_metric_idx);
        % load data for particular GSO and var metric
        filename = strcat('GSO-', GSO, '-VariatonMetric-', var_metric, '.mat');
        filepath = fullfile('energy_distrib/ed_data_new', filename);
        variations_energies = load(filepath);
        [variations, energies] = deal(variations_energies.variations, variations_energies.energies);
        
        % binnerize data
        [bin_edges, bin_counts] = binnerize(variations, energies, bin_struct.(GSO).(var_metric));
        
        % perform plotting
        ax = nexttile(tilenum_start+var_metric_idx);
        histogram(ax, 'BinEdges',bin_edges,'BinCounts', bin_counts);
        axes_row(variation_metrix) = ax;

    end

end

function [bin_counts] = binnerize(variations, energies, bins_min_max)

    bin_edges = linspace(bins_min_max.min, bins_min_max.max);
    thing_to_count = energies(:); %need to go to 1D?
    thing_to_bin = variations(:);

    which_bins = discretize(thing_to_bin, bin_edges);
    %B = accumarray(which_bins,data,sz,fun,fillval);
    bin_counts = groupsummary(thing_to_count, which_bins, 'sum');

end


function [variations, energies] = f(subject_list, atlas, task, include_subcortical, GSO, variation_metric)
%% go through all scans and find the variations and energies wrt to GSO

    roi_idxs = get_roi_idxs(atlas, include_subcortical);
    num_rois = length(roi_idxs);
    ave_node_val = average_node_value(atlas, get_roi_idxs(atlas, include_subcortical));
    scan_directions = ["LR", "RL"];

    variations = zeros(num_rois, 2*length(subject_list));
    energies   = zeros(num_rois, 2*length(subject_list));
    total_scans_processed = 0;
    scans_rejected = 0;
    missing_scs = 0;
    missing_fmris = 0;
    missing_atlas = 0;
    
    DEBUG = false;

    for idx = 1:length(subject_list)
        
        if DEBUG
            fprintf('%d/%d\n', idx, length(subject_list));
        end
        
        for scan_dir_idx = 1:length(scan_directions)
            scan_direction = scan_directions(scan_dir_idx);
            scan_info = ScanInfo(subject_list(idx), atlas, task, scan_direction, include_subcortical);
            try 
                [variations(:, idx), energies(:, idx)] = compute_variations_energies(scan_info, variation_metric, GSO, ave_node_val);
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
    
    if DEBUG
        fprintf('Summary: %d total_scans_processed\n', total_scans_processed);
        fprintf('%d scans_rejected: %d atlas | %d fmris | %d scs\n', scans_rejected, missing_atlas, missing_fmris, missing_scs/2);
    end
    
    variations = variations(:, 1:total_scans_processed);
    energies   = energies(:, 1:total_scans_processed);

    
end

function [variations, energies] = compute_variations_energies(scan_info, variation_metric, GSO, ave_node_val)
    
    A = scan_info.extract_sc();
    [GFT, evals_vec, ~] = scan_info.extract_GFT(GSO);
    V = GFT';
    
    if isequal(variation_metric, "total_variation_L")
        variations = total_variation(V,A);
    elseif isequal(variation_metric, "total_variation_A")
        error("total variation A not yet implimented");
    elseif isequal(variation_metric, "zero_crossings")
        variations = zero_crossings(V,A);
    elseif isequal(variation_metric, "eigenvalues")
        variations = evals_vec;
    else
        error("unrecognized variation metric '%s'.", variation_metric);
    end
    
    %GFT(x1+...+xN) = GFT(x1)+...+GFT(xN)
    x_mean      = mean(scan_info.load_functional_dtseries(), 2) - ave_node_val;
    x_mean_hat  = GFT*x_mean;
    energies    = x_mean_hat.^2;

end

function s = mms(min, max)
    s = struct("min", min, "max", max);
end
        