%%%%%%%%%%%%%%%%%%%%%%5
% This script only supports two parcellations: (1) desikan, (2) destrieux;

clear;

%% setup
%datafolder = '/Users/maxwasserman/Documents/MATLAB/brain_data_preprocess/brain_data';
datafolder = '/Volumes/Elements/brain_data';
%where to save matlab outputs
%desikan_write_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_subset';
%desikan_write_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_cortical_subcortical'; %cortical first
desikan_write_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_subcortical_cortical'; %subcortical first
subcortical_first = true;
destrieux_write_folder = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_destrieux';
%external library for reading .nii data
addpath('/Users/maxwasserman/Documents/MATLAB/brain_data_preprocess/functional_connectivity_processing/cifti_matlab')


tasktype='rfMRI_REST1';
%tasktype = 'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING';
%'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';



%load subject list (list of strings of subject ids)
load Allprocessedid.mat
nSubject = size(subjList); %loaded from Allprocessedid.mat
nSubject = nSubject(1);

[total_time, num_iters, num_iters_left] = deal(0.0, 0, 1114);

%% computation
for i_index=1:nSubject
    subject = subjList(i_index,:);
    
    % Decide on atlas
    atlas = "desikan"; %"destrieux"
    if(strcmp(atlas,"desikan"))                 %missing 4 and 39 are corpus collusum
        ChosenROI_cortical    =   setxor(0:70,[0, 4, 39]);
        ChosenROI_subcortical =   setxor(1:21,[1, 2]);
        write_filename = [desikan_write_folder,'/',subject,'.mat'];
    elseif(strcmp(atlas,"destrieux"))           %missing 42 and 117 are corpus collusum
        ChosenROI_cortical    =   setxor(0:150,[0, 42, 117]); %may have to reorder if you want all left hemisphere regions together for viz
        ChosenROI_subcortical = []; %setxor(1:21,[1, 2]) % same as desikan
        write_filename = [destrieux_write_folder,'/',subject,'.mat'];
    else
        error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
    end
    
    files = dir(write_filename);
    [file_exist, ~] = size(files);
    if file_exist
        continue %skip if already done
    end
    
    % finding which files this patient has
    path_to_LR1 = [datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];

    has_LR1 = isfile(path_to_LR1);
    has_RL1 = isfile(path_to_RL1);

    
    disp([num2str(i_index) ': patient ' subject ' - has LR? ' num2str(has_LR1) ' | has RL? ' num2str(has_RL1)]) 
    
    % attempt lr
    [fc_cov_lr, time_lr] = deal(nan, 0.0);
    if has_LR1
        start = tic;
        name = 'LR';
        [fc_cov_lr, fc_corr_lr] = compute_fmri_cov(atlas, path_to_LR1, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical);
        time_lr = toc(start);
        disp(['   time for LR: ' num2str(time_lr)]);   
    end
    
    % attempt rl
    [fc_cov_rl, time_rl] = deal(nan, 0.0);
    if has_RL1
        start = tic;
        name = 'RL';
        [fc_cov_rl, fc_corr_rl] = compute_fmri_cov(atlas, path_to_RL1, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical);
        time_rl = toc(start);
        disp(['   time for RL: ' num2str(time_rl)]);
    end
    
    
    % intermediate status update
    elapsed_time = time_lr + time_rl;
    num_iters = num_iters + 1;
    num_iters_left = nSubject - i_index;
    total_time = total_time + elapsed_time;
    ave_patient_time = total_time/num_iters;
    expected_time_left = ave_patient_time*num_iters_left/3600; %in hrs
    
    txt = sprintf('===elapsed time: %.2f | ave time: %.1f | expected time remain %.2f (hrs)=== %%\n', elapsed_time, ave_patient_time, expected_time_left);
    disp(txt)
    %disp(['   total time: ' num2str(elapsed_time) ', ave time: ' num2str(ave_patient_time)]);

    save(write_filename, 'subject', 'fc_cov_lr', 'fc_corr_lr', 'fc_cov_rl', 'fc_corr_rl','tasktype', 'ChosenROI_cortical', 'ChosenROI_subcortical', 'atlas', 'subcortical_first')

end