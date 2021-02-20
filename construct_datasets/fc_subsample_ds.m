function [fcs_tensor, scs_tensor] = fc_subsample_ds(rand_seed)

% for further time optimization, make all cached data be loaded in one
% file. Batch to limit memory overflow (all files would be ~7 gb).


%% construct fcs for patient by 
% 1) removing some raw signals
% 2) perform frequency filtering
% 3) computing fcs via windowing or resampling
% 4) removing some fc's

%% Define high level hp's
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
include_subcortical = false;
GSO = "L";

%% which rois to consider
if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
    sum_stat_file = load('data/fmri_desikan_summary_stats.mat');
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;
    error("Atlas " + atlas + " has not been implimented yet...")
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end

if include_subcortical
    num_rois = length(chosen_roi.cortical) + length(chosen_roi.subcortical);
    roi_idxs = (1:num_rois);
else
    num_rois = length(chosen_roi.cortical);
    roi_idxs = (20:num_rois);
end


%% Mean center: compute scalar s to subtract from all signals for mean centering: x-s*1's
average_vector_ = mean([sum_stat_file.lrs.mean_vectors, sum_stat_file.rls.mean_vectors], 2);
average_vector  = average_vector_(roi_idxs,:);
average_node_value = mean(average_vector); %use this to mean center!


%% Define inputs needed for each step in pipeline

%% Preprocess inputs
%create ranges to be used for frequency energy thresholding in ranges
low_freq_cutoffs  = floor(num_rois/3);
med_freq_cutoffs  = 2*low_freq_cutoffs;
low_freq_interval = [1, low_freq_cutoffs];
med_freq_interval = [low_freq_cutoffs+1,  med_freq_cutoffs];
high_freq_interval= [med_freq_cutoffs+1, num_rois];
ranges = {low_freq_interval, med_freq_interval, high_freq_interval};

%preprocess_filter_params = struct('name',"freq_distribution", 'threshold', 90, 'use_percentile', true, 'ranges', ranges);
preprocess_filter_params.name            = 'freq_distribution';
preprocess_filter_params.threshold       = 90;
preprocess_filter_params.use_percentile  = true;
preprocess_filter_params.ranges          = ranges;


%% Frequency filtering inputs
%which to *keep*
freq_filt_flag = false;
freq_filtering_intervals = {[1,num_rois]};


%Subset construction inputs
%subset_contstruction = 'windowing'; %
%subset_construction_options = struct('windowsize', 400, 'movesize', 300);

subset_contstruction = 'sampling';
subset_construction_options = struct('num_subsets', 5, 'sps', 400, 'with_replacement', false);
rng(rand_seed); %for random sampling


%% Post processing inputs
fc_filter_params = struct('name', 'eig', 'which_eig', 1, 'threshold', 90, 'use_percentile', true);


%% determine which patients to do this for
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

% preallocate if slow
fcs = {};
scs = {};
subject_ids = {};

num_iters = 0;
total_time = 0;
num_scans = 0;
for subject_idx = 1:length(subject_list)
    subject     = char(num2str(subject_list(subject_idx,:)));
    path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    %path2fmris = {path_to_LR1, path_to_RL1};
    [cached_LR, ~] = is_cached(subject, atlas, chosen_roi, path_to_LR1);
    [cached_RL, ~] = is_cached(subject, atlas, chosen_roi, path_to_RL1);
    path2fmris = {};
    if ~cached_LR && ~cached_RL %if neither have, go to next
        continue
    end
    if cached_LR
        path2fmris{end+1} = path_to_LR1;
    end
    if cached_RL
        path2fmris{end+1} = path_to_RL1;
    end
   
    
    if freq_filt_flag
        [GFT, evals_vec] = extract_GFT(subject, atlas, include_subcortical, GSO);
        iGFT = GFT';
    end
    
    start = tic();
    fprintf('starting %d/%d: (%d scans)...', subject_idx, length(subject_list), length(path2fmris));
    
    for path2fmris_idx = 1:length(path2fmris)
        path2fmri = path2fmris{path2fmris_idx};
        
    
        %% load and center 'raw' fmri data. 
        x        = process_fmri(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);
        x        = x - average_node_value; %

        if ~include_subcortical
            x = x(20:end, :);
        end

        %{
        %% perform pre-preprocessing of 'raw' signals. This removes some subset of raw signals.
        [which_raw_idxs_remove, ~] = preprocess_signals(x, GFT, preprocess_filter_params);
        x(:, which_raw_idxs_remove) = [];
        %}
        
        %% frequency filter signals
        if freq_filt_flag
            [x] = iGFT*freq_filtering(GFT*x, freq_filtering_intervals);
        end
        
        
        %% compute (sampled/windowed) subsets of signals and their respective fcs
        if isequal(subset_contstruction, 'windowing')
            x_subsets = windowed_signals(x, subset_construction_options.windowsize, subset_construction_options.movesize);
        elseif isequal(subset_contstruction, 'sampling')
            x_subsets = random_subsets(x,  subset_construction_options.num_subsets, subset_construction_options.sps, subset_construction_options.with_replacement);
        else
            error('incorrect option: %s', subset_contstruction);
        end
        
        covs = apply_to_tensor_slices(@(z) cov(z'), x_subsets);
        [~,~,num_windows] = size(covs);

        %{
        %% post-processing on fcs
        [which_fc_idxs_remove, ~] = preprocess_fcs(x_subsets, GFT, covs, fc_filter_params);
        covs(:,:, which_fc_idxs_remove) = [];
        %}
        
        %% place into tensors
        [A] = extract_sc(subject, atlas, include_subcortical);
        
        num_scans = num_scans + 1;
        %place covs into cov tensor
        fcs{end+1} = covs;
        %place A's into sc tensor
        scs{end+1} = A; %repmat(A, 1,1, num_windows); %to save memory, we can repmat in python later
        subject_ids{end+1} = subject;
        
        %optional viz
        %corrs = apply_to_tensor_slices(@corrcov, covs);
        %inspect_fcs_sc(corrs, 99, A);
        
       
    end
    
    stop = toc(start);
    elapsed_time = stop;
    num_iters = num_iters + 1;
    num_iters_left = length(subject_list) - subject_idx;
    total_time = total_time + elapsed_time;
    ave_patient_time = total_time/num_iters;
    expected_time_left = ave_patient_time*num_iters_left/3600; %in hrs
    txt = sprintf('elapsed time: %.2f | ave time: %.1f | expected time remain %.2f (hrs)=== %%\n', elapsed_time, ave_patient_time, expected_time_left);
    disp(txt);
end

[~,~, fcs_per_scan] = size(covs);
num_fcs_total = length(fcs) * fcs_per_scan;
fcs_tensor = zeros(num_rois, num_rois, num_fcs_total);
scs_tensor = zeros(num_rois, num_rois, num_fcs_total);
    
for l = 1:length(fcs)
	global_start_idx = (l-1)*fcs_per_scan + 1;
    global_end_idx   = global_start_idx + fcs_per_scan - 1;
    fcs_tensor(:,:, global_start_idx:global_end_idx) = fcs{l};
	scs_tensor(:,:, global_start_idx:global_end_idx) = repmat(scs{l},1,1,fcs_per_scan);
end


% atlas/tasktype/GSO/include_subcortical/sampling_method

%
%filename = sprintf('


end




function inspect_fcs_sc(fcs, percentile, sc)
    [N,~,num_fcs] = size(fcs);

    f = figure();
    t = tiledlayout('flow');
    axes_list = gobjects(num_fcs,1);


    for i = 1:num_fcs
        ax = nexttile(t);
        axes_list(i) = ax;
    
        fc = fcs(:,:,i);
        imagesc(ax, fc);
    
	
        set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
        yl = sprintf('%d',i);
        ylh = ylabel(ax, yl, 'FontSize', 20);
        ylp = get(ylh, 'Position');
        ext = get(ylh,'Extent');
        set(ylh, 'Rotation',0, 'Position',ylp-[2*ext(3) 0 0])
	
    
        xlim(ax,[1,N]); xlim(ax,'manual');
        ylim(ax,[1,N]); ylim(ax, 'manual');
        daspect(ax,[1 1 1]);
        %pbaspect(ax,[1 1 1]);
    
    
    end
    [cl] = find_cov_clim(fcs, percentile);
    colormap_file = load("viz_fcs/correlations_colormap.mat");
    set(axes_list, 'ColorMap', colormap_file.colorMap);
    set(axes_list,'CLim',[-cl cl]);
    cb = colorbar(axes_list(end),'West', 'AxisLocation','out');
    linkaxes(axes_list,'xy');
    
    ax = nexttile(t);
    title('raw sc');
    imagesc(ax, sc);
    xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    colormap(ax, 'gray');
    set(ax,'CLim',[0 max(max(sc))]);
    %cb = colorbar(ax,'East');
    colorbar(ax, 'AxisLocation','out')
    
    ax = nexttile(t);
    title('thresholded sc');
    imagesc(ax, sc>0);
    xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    colormap(ax, 'gray');
    set(ax,'CLim',[0 1]);
    %cb = colorbar(ax,'East');
    colorbar(ax, 'AxisLocation','out')
    
end