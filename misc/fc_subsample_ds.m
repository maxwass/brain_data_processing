%% construct fcs for patient by 
% 1) removing some raw signals
% 2) computing fcs via windowing or resampling
% 3) removing some fc's

clear; clc; close all;

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

brain_dataset = load('data/brain_dataset_sc_fc_pairs.mat');
subject_list = int2str(brain_dataset.final_subject_list); %all subjects with scs

subject = subject_list(1,:);
path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];

chosen_roi = load('data/desikan_roi', 'roi').roi;
%chosen_roi.subcortical = [];
include_subcortical = false;
GSO = "L";

%% load and center 'raw' fmri data
path2fmri = path_to_LR1;

x        = process_fmri(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);
x_center = x - mean(x,2);

if ~include_subcortical
    x_center = x_center(20:end, :);
end

%% perform pre-preprocessing of 'raw' signals
which_raw_metric = 'freq_distribution';
keep_under    = 0.9;
raw_filter_params = struct('cutoff', keep_under);
use_per       = true;

processed_x   = preprocess_signals(x_center, which_raw_metric, raw_filter_params, use_per);

%% compute (sampled/windowed) subsets of signal and their respective fcs
[windowsize, movesize] = deal(30, 20);
[num_subsets, sps] = deal(100, 30); %sps: samples per subset
window_or_sample = 'window';

if isequal(window_or_sample, 'window')
    x_subsets = windowed_signals(x_center, windowsize, movesize);
elseif isequal(window_or_sample, 'sample')
    x_subsets = random_subsets(signals, num_subsets, sps);
else
    error('incorrect option: %s', window_or_sample);
end
fc_covs   = apply_to_tensor_slices(@(z) cov(z'), x_subsets);

%% pre-process fcs
which_fc_metric  = 'eig';
fc_filter_params = struct('which_eig', 1);

[processed_covs, which_idxs_remove] = ...
    preprocess_fcs(fc_covs, which_fc_metric, fc_filter_params);

%% return processed_covs
