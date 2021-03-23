% Author: Max Wasserman, maxw14k@gmail.com, 1/15/21
%   Modified script from Zhenghwu Zhang, now a professor at UNC
%   With help from Marty Cole, current PhD student at UofR


%% Load raw fmri files, compute, and save FULL cov for each patients LR and/or RL files.


%% setup
clear all;
clc;

% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

%Where to save matlab outputs. We save one .mat file per patient.

%UPDATE THIS TO BE LOCAL FOLDER
desikan_write_folder   = '~/Desktop/geom_dl/data/brain_data/fcs_desikan_subcortical_cortical'; %subcortical first
subcortical_first = true;
destrieux_write_folder = '~/Desktop/geom_dl/data/brain_data/fcs_destrieux';

% cifti_matlab/ includes external libraries for reading .nii data. 
%addpath('fc_processing/cifti_matlab')


%load subject list (list of strings of subject ids)
load('data_accounting/hcp_1200_subject_list.mat')
subject_list = hcp1200_subject_list; %loaded from hcp_1200_subject_list.mat


atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';


[total_time, num_iters, num_iters_left] = deal(0.0, 0, 1114);

%% computation
for i_index=1:length(subject_list)
    subject = subject_list(i_index,:); %must be char array for [...] to work later

    if(strcmp(atlas,"desikan"))
        chosen_roi     = load('data/desikan_roi_zhengwu', 'roi').roi;
        write_filename = [desikan_write_folder,'/',subject,'.mat'];
    elseif(strcmp(atlas,"destrieux"))
        chosen_roi     = load('data/destrieux_roi_zhengwu', 'roi').roi;
        write_filename = [destrieux_write_folder,'/',subject,'.mat'];
    else
        error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
    end
    

    files = dir(write_filename);
    [file_exist, ~] = size(files); %if any files exist in patients dir
    if file_exist
        continue %skip if already exists
    end
   
    
    path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    %check which fmri scans exist - can have neither, one, or both
    has_LR1 = isfile(path_to_LR1);
    has_RL1 = isfile(path_to_RL1);
    
    disp([num2str(i_index) ': patient ' subject ' - has LR? ' num2str(has_LR1) ' | has RL? ' num2str(has_RL1)]) 
    
    % attempt lr
    [fc_cov_lr, time_lr] = deal(nan, 0.0);
    if has_LR1
        start = tic;
        name = 'LR';
        
        dtseries_lr = load_functional_dtseries(atlas, path_to_LR1, subject, rawdatafolder, chosen_roi);
        fc_cov_lr = cov(dtseries_lr');
        fc_corr_lr= corrcov(fc_cov_lr);
        
        time_lr = toc(start);
        disp(['   time for LR: ' num2str(time_lr)]);   
    end
    
    % attempt rl
    [fc_cov_rl, time_rl] = deal(nan, 0.0);
    if has_RL1
        start = tic;
        name = 'RL';
        
        dtseries_rl = load_functional_dtseries(atlas, path_to_RL1, subject, raw_hcp_datafolder, chosen_roi);
        fc_cov_rl = cov(dtseries_rl');
        fc_corr_rl= corrcov(fc_cov_rl);
        
        time_rl = toc(start);
        disp(['   time for RL: ' num2str(time_rl)]);
    end
    
    
    % print intermediate status update for user:
    %  elapsed_time = time spent in this patient
    %  num_iters = # iterations performed (we may have skipped already
    %   computed patients)
    %  total_time = time spent since beginning of program
    %  ave_patient_time = average time spent on each patient so far
    %  expected_time_left = approx how much time (in hours)
    elapsed_time = time_lr + time_rl; 
    num_iters = num_iters + 1;
    num_iters_left = length(subject_list) - i_index;
    total_time = total_time + elapsed_time;
    ave_patient_time = total_time/num_iters;
    expected_time_left = ave_patient_time*num_iters_left/3600; %in hrs
    
    txt = sprintf('===elapsed time: %.2f | ave time: %.1f | expected time remain %.2f (hrs)=== %%\n', elapsed_time, ave_patient_time, expected_time_left);
    disp(txt)
    %disp(['   total time: ' num2str(elapsed_time) ', ave time: ' num2str(ave_patient_time)]);
    
    % for each patient save a file containing its functional connnectivies
    %  (covariance, but also correlation for convienece) and other metadata
    subject = string(subject);
    tasktype = string(tasktype);
    save(write_filename, 'subject', 'fc_cov_lr', 'fc_corr_lr', 'fc_cov_rl', 'fc_corr_rl','tasktype', 'ChosenROI_cortical', 'ChosenROI_subcortical', 'atlas', 'subcortical_first')

end