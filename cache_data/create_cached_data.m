%% load raw fmri files, process with desikan atlas. Save to mat file for 
% fast access (don't need to read fmri file now)

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
%load('data_accounting/hcp_1200_subject_list.mat');
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
sub_tasktype = 'REST1';
include_subcortical = true;
subcortical_first = true;


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
num_completed = 0;
num_processed = 0;
total_time    = 0.0;

already_completed = false(size(subject_list));
processed = false(size(subject_list));
times = zeros(size(subject_list));


counter = 0;
for i_index = 1:length(subject_list)
    subject = char(string(subject_list(i_index)));
    
    %subject = subject_list(i_index,:); %must be char array for [...] to work later
       
	cached_filename_lr = [cached_data_folder,'/', tasktype, '/', subject,'_LR.mat'];
	cached_filename_rl = [cached_data_folder,'/', tasktype, '/', subject,'_RL.mat'];
    
    
    has_cached_lr = isfile(cached_filename_lr);
    has_cached_rl = isfile(cached_filename_rl);
    
    path_to_fmri_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_fmri_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    %check which fmri scans exist - can have neither, one, or both
    has_fmri_LR1 = isfile(path_to_fmri_LR1);
    has_fmri_RL1 = isfile(path_to_fmri_RL1);
    
    % attempt lr
    time_lr = 0.0;
    load_and_save_lr = has_fmri_LR1 && ~has_cached_lr;
    if load_and_save_lr
        start = tic;
        fmri_path = path_to_fmri_LR1;
        scan = ['LR_' sub_tasktype];
        dtseries = process_fmri(atlas, fmri_path, subject, raw_hcp_datafolder, chosen_roi);
        par_save_dtseries(cached_filename_lr, atlas, scan, tasktype, chosen_roi, dtseries)
        time_lr = toc(start);   
    end
    
    % attempt rl
    time_rl = 0.0;
    load_and_save_rl = has_fmri_RL1 && ~has_cached_rl;
    if load_and_save_rl
        start = tic;
        fmri_path = path_to_fmri_RL1;
        scan = ['RL_' sub_tasktype];
        dtseries = process_fmri(atlas, fmri_path, subject, raw_hcp_datafolder, chosen_roi);
        par_save_dtseries(cached_filename_rl, atlas, scan, tasktype, chosen_roi, dtseries)
        time_rl = toc(start);
    end
    runtime = time_lr + time_rl;
    
    
    % we didnt do any loading, don't include in time averaging
    if ~load_and_save_lr && ~load_and_save_rl
        num_completed = num_completed + 1; % number *already* completed  
        fprintf('%d: %s already processed\n', i_index, subject);
        
    else

        num_processed = num_processed + 1; % number we processed in this execution of program
        total_time  = total_time + runtime;
        ave_runtime = total_time/num_processed;
        num_remain  = num_subjects - num_processed - num_completed;
        time_remain = num_remain*ave_runtime/3600;
        txt = sprintf('%d: %s - runtime: %.2f | ave runtime: %.2f | ~time remain %.2f (hrs)\n',...
            i_index, subject, runtime, ave_runtime, time_remain);
        disp(txt)

    end

    %fprintf('%d: %s - processed in %.2f\n', i_index, subject, runtime);

end


function par_save_dtseries(filename, atlas, scan, tasktype, chosen_roi, dtseries)
    save(filename, 'atlas', 'scan', 'tasktype', 'chosen_roi', 'dtseries');
end
