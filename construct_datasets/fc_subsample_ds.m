%% construct fcs for patient by 
% 1) removing some raw signals
% 2) perform frequency filtering
% 3) computing fcs via windowing or resampling
% 4) removing some fc's

clear; clc; close all;

%% Define high level hp's
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
include_subcortical = false;
GSO = "L";

if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;h
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end

if include_subcortical
    num_rois = length(chosen_roi.cortical) + length(chosen_roi.subcortical);
else
    num_rois = length(chosen_roi.cortical);
end


%% Define inputs needed for each step in pipeline
%Preprocess inputs
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


%Frequency filtering inputs
%which to *keep*
freq_filtering_intervals = {[2,num_rois]};


%Subset construction inputs
subset_contstruction = 'windowing'; %'sampling'
subset_construction_options = struct('windowsize', 150, 'movesize', 100);
%subset_construction_options = struct('num_subsets', 5, 'sps', 400, 'with_replacement', false);


%Post processing inputs
fc_filter_params = struct('name', 'eig', 'which_eig', 1, 'threshold', 90, 'use_percentile', true);


%% determine which patients to do this for
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

for subject_idx = 1:length(subject_list)
    subject     = char(num2str(subject_list(subject_idx,:)));
    path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    path2fmris = {path_to_LR1, path_to_RL1};
    [GFT, evals_vec] = extract_GFT(subject, atlas, include_subcortical, GSO);
    iGFT = GFT';
    
    for path2fmris_idx = 1:length(path2fmris)
        path2fmri = path2fmris{path2fmris_idx};
        
    
        %% load and center 'raw' fmri data. 
        x        = process_fmri(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);
        x_center = x - mean(x,2);

        if ~include_subcortical
            x_center = x_center(20:end, :);
        end

        %% perform pre-preprocessing of 'raw' signals. This removes some subset of raw signals.
        [which_raw_idxs_remove, ~] = preprocess_signals(x_center, GFT, preprocess_filter_params);
        x_center(:, which_raw_idxs_remove) = [];
        processed_x = x_center;
        
        %% frequency filter signals
        [freq_filtered_x] = iGFT*freq_filtering(GFT*processed_x, freq_filtering_intervals);
        
        
        %% compute (sampled/windowed) subsets of signal and their respective fcs
        if isequal(subset_contstruction, 'windowing')
            x_subsets = windowed_signals(freq_filtered_x, subset_construction_options.windowsize, subset_construction_options.movesize);
        elseif isequal(subset_contstruction, 'sampling')
            x_subsets = random_subsets(freq_filtered_x,  subset_construction_options.num_subsets, subset_construction_options.sps, subset_construction_options.with_replacement);
        else
            error('incorrect option: %s', subset_contstruction);
        end
        
        covs = apply_to_tensor_slices(@(z) cov(z'), x_subsets);
        [~,~,num_windows] = size(covs);

        %% post-processing on fcs
        [which_fc_idxs_remove, ~] = preprocess_fcs(x_subsets, GFT, covs, fc_filter_params);
        covs(:,:, which_fc_idxs_remove) = [];
        
    end
end