%% average energy distribution in freq domain with different GSOs

% construct (freq, energy) pairs for each patient
close all;
clear; clc;

%% directory to read fmri from, and write processed data to
% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

% It is a directory of processed fmri signals according to the desikan atlas. 
%One .mat file per scan (one for LR, one for RL, if they exist).
%This directory is NOT included in repo: ~ 1GB
cached_desikan   = '~/Documents/MATLAB/brain_data_preprocess/data/cached_desikan';
cached_destrieux = '~/Documents/MATLAB/brain_data_preprocess/data/cached_destrieux';


%% determine which patients to do this for
% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
%load('data/hcp_1200_subject_list.mat');
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

%% determine which nodes to do this for
include_subcortical = false;
roi_idxs = get_roi_idxs(atlas, include_subcortical);
num_rois = length(roi_idxs);
[ave_node_val] = average_node_value(atlas, roi_idxs);

%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
sub_tasktype = 'REST1';
GSO = 'L';
%GSO = 'L_norm';
%GSO = 'A';
%GSO = 'A_norm';


if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_desikan;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_destrieux;
else
	error("Atlas %s not supported yet", atlas)
end


%% perform loading and saving
% load each fmri file, save processed dtseries to .mat file, then update
%  completed array for time estimate
% save full cov (87x87) and remove subcortical later

num_subjects  = length(subject_list);
mean_freq_signal = zeros(num_rois, 2*num_subjects);
freqs            = zeros(num_rois, 2*num_subjects);

counter = 0;
for i_index = 1:length(subject_list)
    subject = char(string(subject_list(i_index)));
    fprintf('%d: subject %s\n', i_index, subject);
       
    [cached_path_lr, is_cached_lr] = cached_filepath(atlas, tasktype, subject, "LR");
    [cached_path_rl, is_cached_rl] = cached_filepath(atlas, tasktype, subject, "RL");

    A = extract_sc(subject, atlas, include_subcortical);
	[GFT, eigvals, ~] = extract_GFT(subject, atlas, include_subcortical, GSO);
        
        
    % attempt load lr
    if has_cached_lr
        dtseries_lr = load(cached_path_lr).dtseries;
        x_mean_lr = mean(dtseries_lr, 2) - ave_node_val ; %across rows
        x_mean_lr = x_mean_lr(roi_idxs, :);
        x_mean_freq_lr = GFT*x_mean_lr;  
        counter = counter + 1;
        mean_freq_signal(:, counter)  = x_mean_freq_lr;
        freqs(:, counter) = eigvals;
    end
    
    % attempt load rl
    if has_cached_rl
        dtseries_rl = load(cached_path_rl).dtseries;
        x_mean_rl = mean(dtseries_rl, 2) - ave_node_val ; %across rows
        x_mean_rl = x_mean_rl(roi_idxs, :);
        x_mean_freq_rl = GFT*x_mean_rl; 

        counter = counter + 1;
        mean_freq_signal(:, counter) = x_mean_freq_rl;
        freqs(:, counter) = eigvals;
    end
    
end

% allocated larger array, remove unused part
mean_freq_signal = mean_freq_signal(:, 1:counter);
freqs = freqs(:, 1:counter);

% save this data
filename = sprintf('GSO-%s_%s_%s--freq_distrib', GSO, atlas, sub_tasktype);
if include_subcortical
    filepath = fullfile("energy_distrib", "ed_data", "GFT_using_subcortical", filename);
else
    filepath = fullfile("energy_distrib", "ed_data", "GFT_only_cortical", filename);
end
save(filepath, 'mean_freq_signal', 'freqs');
   
