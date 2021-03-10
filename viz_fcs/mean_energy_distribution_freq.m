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
cached_desikan   = '~/Documents/MATLAB/brain_data_preprocess/cached_desikan';
cached_destrieux = '~/Documents/MATLAB/brain_data_preprocess/cached_destrieux';


%% determine which patients to do this for
% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
%load('data/hcp_1200_subject_list.mat');
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

%% determine which nodes to do this for
include_subcortical = true;
subcortical_first = true;
if include_subcortical==0
    roi_idxs = (20:87); %subcortical are first. See README
else
	roi_idxs = (1:87);
end
num_rois = length(roi_idxs);

%% values to compute
%what is the average value in *each* brain region (roi) in each scan
%  average across rows
%mean_rois = zeros(num_rois, num_subjects);

%what is the average value of each *vector observation* in each scan
%  average across columns
%num_obsvs = 1200;
%mean_obsvs = zeros(1, num_subjects * num_obsvs);
sum_stat_file = load('data/fmri_desikan_summary_stats.mat');
   
average_vector_ = mean([sum_stat_file.lrs.mean_vectors, sum_stat_file.rls.mean_vectors], 2);
average_vector  = average_vector_(roi_idxs,:);
average_node_value = mean(average_vector);

%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
sub_tasktype = 'REST1';
GSO = 'A_norm';

if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_desikan;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_destrieux;
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end


%% perform loading and saving
% load each fmri file, save processed dtseries to .mat file, then update
%  completed array for time estimate
% save full cov (87x87) and remove subcortical later

num_subjects  = length(subject_list);
ave_node_val = 70000;


mean_freq_signal = zeros(num_rois, 3000);
freqs = zeros(num_rois, 3000);

counter = 0;
num_patients_with_disc_graph = 0;

for i_index = 1:length(subject_list)
    subject = char(string(subject_list(i_index)));
    fprintf('%d: subject %s\n', i_index, subject);
    
    
    %subject = subject_list(i_index,:); %must be char array for [...] to work later
       
	cached_filename_lr = [cached_data_folder,'/', tasktype, '/', subject,'_LR.mat'];
	cached_filename_rl = [cached_data_folder,'/', tasktype, '/', subject,'_RL.mat'];
    
    has_cached_lr = isfile(cached_filename_lr);
    has_cached_rl = isfile(cached_filename_rl);
    
    A = extract_sc(subject, atlas, include_subcortical);
    D_vec = sum(A,2);
    if any(~sum(A,2))
        d = sum(~sum(A,2));
        warning('   subject %s has %d nodes with no edges\n', subject, d);
        num_patients_with_disc_graph = num_patients_with_disc_graph + 1;
        fprintf('   so far %d patients have disconnected graphs\n', num_patients_with_disc_graph);
    end
	[GFT, eigvals] = extract_GFT(subject, atlas, include_subcortical, GSO);
        
        
    % attempt load lr
    if has_cached_lr
        dtseries_lr = load(cached_filename_lr).dtseries;
        x_mean_lr = mean(dtseries_lr, 2) - ave_node_val ; %across rows
        x_mean_freq_lr = GFT*x_mean_lr;  
        counter = counter + 1;
        mean_freq_signal(:, counter)  = x_mean_freq_lr;
        freqs(:, counter) = eigvals;
    end
    
    % attempt load rl
    if has_cached_rl
        dtseries_rl = load(cached_filename_rl).dtseries;
        x_mean_rl = mean(dtseries_rl, 2) - ave_node_val ; %across rows
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
save(filename, 'mean_freq_signal', 'freqs');
   