%% create fc window datasets

clear; clc;

%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
include_subcortical = true;
subcortical_first = true;


%% It is a directory of processed fmri signals according to the desikan atlas. 
%One .mat file per scan (one for LR, one for RL, if they exist).
%This directory is NOT included in repo: ~ 1GB
cached_data_folder_desikan   = '~/Documents/MATLAB/brain_data_preprocess/cached_brain_data_desikan';
cached_data_folder_destrieux = '~/Documents/MATLAB/brain_data_preprocess/cached_brain_data_destrieux';

if atlas =="desikan"
    cached_data_folder = cached_data_folder_desikan;
else
    cached_data_folder = cached_data_folder_destrieux;
end

%% determine which patients to do this for
% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
%load('data/hcp_1200_subject_list.mat');
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

windowed_covs_lr = cell(size(subject_list));
windowed_covs_rl = cell(size(subject_list));

for i = 1:length(subject_list)
    subject = char(string(subject_list(i)));
    cached_lr_filename = [cached_data_folder '/' subject '_LR.mat'];
    cached_rl_filename = [cached_data_folder '/' subject '_RL.mat'];
    
    if isfile(cached_lr_filename)
        file     = load(cached_lr_filename);
        dtseries = file.dtseries;
        mean_centered_dtseries = dtseries - mean(dtseries,2);
        
        %pre-filtering on raw fmri signals
        filtered_dtseries   = filter_dtseries(mean_centered_dtseries, subject, atlas, include_subcortical, which_filters);
        
        
        [signal_windows] = windowed_signals(filtered_dtseries, windowsize, movesize);
        [covs]            = construct_fcs(signal_windows);
        
        
        %post-filtering on averaged window'ed vectors
        filter_idxs = filter_windowed_signals(signal_windows, subject, atlas, include_subcortical, which_filters);
        
        
        %filter_idxs is a row vector: ith entry true = keep column i
        filter_covs = covs(:, :, boolean(filter_idxs)); %apply columnwise!
        
        
        windowed_covs_lr{i} = filter_covs;
        
    end
    
    if isfile(cached_rl_filename)
        rl_file = load(cached_rl_filename);
    end
end
    
%% Pipeline

% load data

%create fc's w/ optional filter
% create average signals in window
[signal_windows] = windowed_signals(dtseries, windowsize, movesize);

% NEED
% filter
%[metrics] = some metric of ave_signals

% pre ave filters


% post ave filters
% ex) only take ave signals with energy less than average
%        energies = norm of cols 
%        above_ave_energies = energies >= mean(energies)
%        energies(above_ave_energies) = []
% ex) low pass filter ave signals
%        freq ave signals = GFT(ave_signals)
%        lpf_freq_ave_signals = set high freqs to zero
    

% compute covariances


%save to individual mat file