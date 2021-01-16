%% test whether properly placed fcs and scs into dataset tensors

clear all;

%dataset variables
% fcs_* - tensor of fc matrices
% *_scs - tensor of sc matrices
%
load '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/brain_dataset_cov.mat'
clear('atlas', 'ChosenROI_cortical', 'ChosenROI_subcortical', 'tasktype')

%use this to load each individuals .mat file with fcs
fc_data_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_subcortical_cortical';

%scs
% subject_list (1065x1 int64)
% scs (87x87x1065 double)
sc_file = load('~/Documents/MATLAB/brain_data_preprocess/data/scs_desikan.mat');
SCs = sc_file.scs; 
subject_list_sc = sc_file.subject_list;

%loop over every subject in final dataset
%*Invariant*the ith entry (subject id) in the final_subject_list 
%  corresponds to the ith matrix in the other tensors
for idx =1:length(final_subject_list)
    
    subject = final_subject_list(idx);
    fprintf("%d: subject %d...testing...\n", idx, subject)
    
    %load patient mat file with fc data: fcs_*. This is the data straight 
    % from fmri processing.
    fc_data_mat_file = sprintf('%s/%s.mat',fc_data_folder, num2str(subject));
    fcs_and_metadata =  matfile(fc_data_mat_file);
    
 
    s = ~(string(subject)==fcs_and_metadata.subject);

    % if any lr is not equal, lr will be true
    lr = ~all(all(fcs_cov_lr(:,:,idx)==fcs_and_metadata.fc_cov_lr));
    lr = lr || ~all(all(fcs_corr_lr(:,:,idx)==fcs_and_metadata.fc_corr_lr));
    
    %account for possible nan. 
    %fc_cov_lr will be single NaN. 
    %fcs_cov_lr will be a matrix of NaN's.
    if all(all(isnan(fcs_and_metadata.fc_cov_lr))) && all(any(isnan(fcs_cov_lr(:,:,idx))))
        lr = false;
    end
    
    % if any rl is not equal, rl will be true
    rl = ~all(all(fcs_cov_rl(:,:,idx)==fcs_and_metadata.fc_cov_rl));
    rl = rl || ~all(all(fcs_cov_rl(:,:,idx)==fcs_and_metadata.fc_cov_rl));
    
    %account for possible nan. 
    %fc_cov_rl will be single NaN. 
    %fcs_cov_rl will be a matrix of NaN's.
    if all(all(isnan(fcs_and_metadata.fc_cov_rl))) && all(any(isnan(fcs_cov_rl(:,:,idx))))
        rl = false;
    end

    %SCs is the data given by Zhengwu. The scs that have corresponding fcs
    %were placed in the raw_scs tensor. Use subject to find indices in
    %respective subject_lists to compare (what should be equal) scs.
    index_in_subject_list_sc = int64(find(subject_list_sc == subject));
    sc     = ~all(all(SCs(:,:,index_in_subject_list_sc)==raw_scs(:,:,idx)));
    
    %check that idx'th matrix is equal to idx'th subjects
    any_not_equal = s || lr || rl || sc;
    if any_not_equal 
        fprintf("subject: %d, lr: %d, rl: %d, sc: %d\n", s, lr, rl, sc);
        error('SOMETHING IS WRONG AT IDX %d SUBJECT %d\n', idx, subject);
    end
end

fprintf("All subject seem to match their original data files. Success.")
