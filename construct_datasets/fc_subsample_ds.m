function fc_subsample_ds(intervals_to_keep, include_subcortical)
%% Define inputs needed for each step in pipeline

%% patient and scan info
info = struct("atlas", 'desikan', "tasktype", 'rfMRI_REST1', "include_subcortical", include_subcortical, "GSO", 'L', "rand_seed", 10);
rng(info.rand_seed);


%% Preprocess inputs
preprocess_filter.filter          = false;
preprocess_filter.name            = 'freq_distribution';
preprocess_filter.threshold       = 90;
preprocess_filter.use_percentile  = true;
preprocess_filter.splits          = 3;

%% Frequency filtering inputs
freq_filter.filter = true;
freq_filters.intervals_to_keep = intervals_to_keep;
if info.include_subcortical
    freq_filter.intervals_to_keep = intervals_to_keep;% {[2,29],[59,87]}; %
else
    freq_filter.intervals_to_keep = intervals_to_keep; %{[2,68]};
end
freq_filter.intervals_txt = intervals_to_string(freq_filters.intervals_to_keep);


%% Subset construction inputs
subset_construction = struct("name", 'full');
%subset_construction = struct("name", 'windowing', 'windowsize', 400, 'movesize', 400);
%subset_construction = struct("name", 'sampling', "num_subsets", 3, "sps", 700, "with_replacement", false);

%% Post processing inputs
fc_filter = struct('filter', false, 'name', 'eig', 'which_eig', 1, 'threshold', 90, 'use_percentile', true);

%% create unique filename based on filters and save
filepath = unique_filename(preprocess_filter, freq_filter, subset_construction, fc_filter, info);

%% Create actual dataset
[data, chosen_roi] = create_dataset(preprocess_filter, freq_filter, subset_construction, fc_filter, info);

save(filepath, "data", "subset_construction", "preprocess_filter", "freq_filter", "fc_filter", "info", "chosen_roi", '-v7');

end

function filepath = unique_filename(preprocess_filter, freq_filter, subset_construction, fc_filter, info)

%% combine info and filters for unique filename

if isequal(subset_construction.name, "windowing")
    subset_construction_txt = sprintf("windowsize%d_movesize%d", ...
        subset_construction.windowsize,...
        subset_construction.movesize);
elseif isequal(subset_construction.name, "sampling")
    subset_construction_txt = sprintf("numsubsets%d_sps%d_withreplacement%d", ...
        subset_construction.num_subsets,...
        subset_construction.sps,...
        subset_construction.with_replacement);
else
    %using full to make fc
    subset_construction_txt = "";
end

if preprocess_filter.filter
    error("not implimented");
end

if freq_filter.filter
    freq_intervals_kept_txt = sprintf("FreqIntervalsKept%s", freq_filter.intervals_txt);
else
    freq_intervals_kept_txt = sprintf("FreqIntervalsKept1-87");
end

if fc_filter.filter
    error("not implimented");
end

filename = sprintf("%s%s", subset_construction_txt, freq_intervals_kept_txt);

%% each subset technique has own folder
if isequal(subset_construction.name, "windowing")
    filepath = fullfile("data", "cached_subset_datasets_REST1", "windowing", filename);
elseif isequal(subset_construction.name, "sampling")
    filepath = fullfile("data", "cached_subset_datasets_REST1", "sampling", filename);
elseif isequal(subset_construction.name, "full")
    filepath = fullfile("data", "cached_subset_datasets_REST1", "full", filename);
else
    error('Dont know where to save data to');
end

end



function [data, chosen_roi] = create_dataset(preprocess_filter, freq_filter, subset_construction, fc_filter, info)
% subset_construction ::  struct. 
%  If name == 'windowing'
%       fields "windowsize" :: int, "movesize" :: int
%  If name == 'sampling'
%       fields "num_subsets" :: int, "sps" :: int, "with_replacement" :: logical
% in {'windowing', 'sampling'}
% freq_filter :: struct
%   fields: filter :: logical, intervals_to_keep :: [cell] of integer
%   ranges to keep
% preprocess_filter :: struct
%   fields: filter, name, threshold, use_percentile, ranges
% fc_filter :: struct



% for further time optimization, make all cached data be loaded in one
% file. Batch to limit memory overflow (all files would be ~7 gb).


%% construct fcs for patient by 
% 1) removing some raw signals
% 2) perform frequency filtering
% 3) computing fcs via windowing or resampling
% 4) removing some fc's


%% Define high level hp's
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

[atlas, tasktype, include_subcortical, GSO] = ...
    deal(info.atlas, info.tasktype, info.include_subcortical, info.GSO);

%% which rois to consider
if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
    sum_stat_file      = load('data/fmri_desikan_summary_stats.mat');
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;
    error("Atlas " + atlas + " has not been implimented yet...")
    elsea
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end

if include_subcortical
    num_rois = length(chosen_roi.cortical) + length(chosen_roi.subcortical);
    roi_idxs = (1:num_rois);
else
    num_rois = length(chosen_roi.cortical);
    roi_idxs = (20:num_rois);
end

if preprocess_filter.filter
    preprocess_filter.ranges = construct_freq_ranges(num_rois, preprocess_filter.splits);
end


%% Mean center: compute scalar s to subtract from all signals for mean centering: x-s*1's
average_vector_ = mean([sum_stat_file.lrs.mean_vectors, sum_stat_file.rls.mean_vectors], 2);
average_vector  = average_vector_(roi_idxs,:);
average_node_value = mean(average_vector); %use this to mean center!


%% determine which patients to do this for
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list   = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

%% loop through all scans
% preallocate if slow
data = repmat(struct('subject_id',-1,'lr_fcs',[], 'rl_fcs', [],'sc',zeros(num_rois,num_rois)),length(subject_list),1);
fcs = {};
scs = {};
fcs_per_scan = {};
subject_ids = zeros(3000,1,'uint32');



num_iters = 0;
total_time = 0;
num_scans = 0; %total scans
total_fcs = 0;
for subject_idx = 1:length(subject_list)
    subject_int = subject_list(subject_idx,:);
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
   
    % if performing freq filtering, extract once for all scans (all share
    % same sc)
    if freq_filter.filter
        [GFT, evals_vec] = extract_GFT(subject, atlas, include_subcortical, GSO);
        iGFT = GFT';
    end
    
    %% extract once for both lr/rl
    [A] = extract_sc(subject, atlas, include_subcortical);
    
    %% populate non-fc components of struct
    begin_flag = true;
    data(subject_idx).subject_id = subject_int;
    data(subject_idx).sc = A;
    
    start = tic();
    fprintf('starting %d/%d: (%d scans)...', subject_idx, length(subject_list), length(path2fmris));
    
    for path2fmris_idx = 1:length(path2fmris)
        path2fmri = path2fmris{path2fmris_idx};
        
    
        %% load and center 'raw' fmri data. 
        x        = load_functional_dtseries(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);
        x        = x - average_node_value; %

        if ~include_subcortical
            x = x(20:end, :);
        end

        %% perform pre-preprocessing of 'raw' signals. This removes some subset of raw signals.
        if preprocess_filter.filter
            [which_raw_idxs_remove, ~] = samples_to_remove(x, GFT, preprocess_filter);
            x(:, which_raw_idxs_remove) = [];
        end
        
        %% frequency filter signals
        if freq_filter.filter
            %[x] = iGFT*freq_filtering_idx(GFT*x, freq_filter.intervals_to_keep);
            [idxs_to_remove] = ~intervals_to_logical_vec(freq_filter.intervals_to_keep, eigvals);
            x_hat = GFT*x;
            x_hat(idxs_to_remove,:) = 0;
            x = iGFT*x_hat;
        end
                
        %% compute (sampled/windowed) subsets of signals and their respective fcs
        if isequal(subset_construction.name, 'windowing')
            x_subsets = windowed_signals(x, subset_construction.windowsize, subset_construction.movesize);
            covs = apply_to_tensor_slices(@(z) cov(z'), x_subsets);
            [~,~,num_fcs] = size(covs);
        elseif isequal(subset_construction.name, 'sampling')
            x_subsets = random_subsets(x,  subset_construction.num_subsets, subset_construction.sps, subset_construction.with_replacement);
            covs = apply_to_tensor_slices(@(z) cov(z'), x_subsets);
            [~,~,num_fcs] = size(covs);
        elseif isequal(subset_construction.name, 'full')
            covs = cov(x');
            num_fcs = 1;
        else
            error('incorrect option: %s', subset_construction);
        end
        
        %% post-processing on fcs
        if fc_filter.filter
            [which_fc_idxs_remove, ~] = preprocess_fcs(x_subsets, GFT, covs, fc_filter);
            covs(:,:, which_fc_idxs_remove) = [];
            error('correct tensor construction...relies on same number of fcs per patient');
        end
        
        
        % initialize data struct now that we know the sizes
        if (subject_idx == 1) && begin_flag
            s = struct('subject_id',-1,'contains_lr',false, 'contains_rl', false, 'lr_fcs', zeros(num_rois,num_rois,num_fcs), 'rl_fcs', zeros(num_rois,num_rois,num_fcs), 'sc', zeros(num_rois,num_rois));
            data = repmat(s,length(subject_list),1);
            data(1).sc = A;
            data(1).subject_id = subject_int;
            begin_flag = false;
        end
        
        %populate covs portion of struct
        if contains(path2fmri, 'LR')       
            data(subject_idx).contains_lr = true;
            data(subject_idx).lr_fcs = covs;
        elseif contains(path2fmri, 'RL')
            data(subject_idx).contains_rl = true;
            data(subject_idx).rl_fcs = covs;
        else
            error('unrecgonized scan dir');
        end
        
        
        num_scans = num_scans + 1;
        fcs{end+1} = covs;
        scs{end+1} = A; %repmat(A, 1,1, num_windows); %to save memory, we can repmat in python later
        subject_ids(num_scans)  = int32(str2double(subject)); %correct index?

        total_fcs = total_fcs + num_fcs;
        if num_fcs < 1
            error('%s: filtered out ALL fcs.', subject);
        end
        
        
        
        
        
        
        %optional viz
        %corrs = apply_to_tensor_slices(@corrcov, covs);
        %inspect_fcs_sc(corrs, 99, A);
        %inspect_fcs_sc(corrcov(covs), 99, A);
        
       
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

%{
fcs_tensor = zeros(num_rois, num_rois, total_fcs);
scs_tensor = zeros(num_rois, num_rois, total_fcs);
matlab_indices_tensor = zeros(total_fcs,2, 'uint16');

global_start_idx = 1;
for l = 1:length(fcs)
    global_end_idx = global_start_idx + fcs_per_scan{l} - 1; % may be different number of scans if filtered
    fcs_tensor(:,:, global_start_idx:global_end_idx) = fcs{l};
	scs_tensor(:,:, global_start_idx:global_end_idx) = repmat(scs{l},1,1,fcs_per_scan{l});
    matlab_indices_tensor(l,:) = [global_start_idx, global_end_idx];
    global_start_idx =  global_end_idx+1;
    
end

python_indices_tensor = matlab_indices_tensor-1;

subject_ids = subject_ids(1:num_scans); %check that num_scans correct to use
%}

end

function  txt = intervals_to_string(intervals)
    txt = '';
    for j = 1:length(intervals)
        interval_j = intervals{j};
        new_txt = sprintf('%d-%d_', interval_j(1), interval_j(2));
        txt = strcat(txt, new_txt);
    end
    
    txt = txt(1:end-1);
end

function [ranges] = construct_freq_ranges(num_rois, splits)
%num_rois must be divisible by splits

%create ranges to be used for frequency energy thresholding in ranges
interval_length = floor(num_rois/splits);
ranges = cell(splits,1);

first = 1;
last  = interval_length;
for l = 1:(splits-1)
    ranges(l) = [first, last];
    first = first + interval_length;
    last  = min(last + interval_length, num_rois);
end

ranges{end} = [first, num_rois]; %last one always has equal to/more than rest

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
