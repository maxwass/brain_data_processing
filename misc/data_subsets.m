%% go from full data tensors (with nans for missing data) to subsets

clear; clc;
load brain_dataset_sc_fc_pairs.mat
load hcp_1200_subjects_missing_fc_data.mat

%% to lr only

have_lr_load = setdiff(final_subject_list, missing_LR);

have_lr_ids = [];
have_lr_idxs = [];
missing_lr_ids = [];
missing_lr_idxs = [];


% which subjects have lr?
have_lr_ids  = setdiff(final_subject_list, missing_LR); %ids in list but not missing_lr
[log, have_LR_idxs] = ismember(have_lr_ids, final_subject_list);
%which subject which are in this subset are missing lr?
missing_lr = intersect(final_subject_list, missing_LR); 


% which subjects have rl?
have_rl_ids  = setdiff(final_subject_list, missing_RL); %ids in list but not missing_lr
[log, have_RL_idxs] = ismember(have_rl_ids, final_subject_list);
%which subject which are in this subset are missing lr?
missing_rl = intersect(final_subject_list, missing_RL); 



%{
for idx = 1:length(final_subject_list)
    cov_lr = fcs_cov_lr(:,:,idx);
    if all(all(isnan(cov_lr)))
        missing_lr_ids  = union(missing_lr_ids, final_subject_list(idx));
        missing_lr_idxs = union(missing_lr_idxs, idx);
    else        
        have_lr_ids  = union(have_lr_ids, final_subject_list(idx));
        have_lr_idxs = union(have_lr_idxs, idx);
    end
end
%}


%% create dataset of lrs

%fcs_cov_lr_full = (fcs_cov_lr(:,:,missing_lr_idx) = []);



a = all(intersect(have_lr_load,have_lr_ids) == have_lr_load);

%given list of indices to remove, create full dataset
fcs_cov_lr_full  = fcs_cov_lr(:,:,have_lr_idxs);
fcs_corr_lr_full = fcs_corr_lr(:,:,have_lr_idxs);
raw_scs_lr = raw_scs(have_lr_idxs);
transform_scs_lr = transform_scs(have_lr_idxs);
subject_list_lr  = final_subject_list(have_lr_idxs);



 