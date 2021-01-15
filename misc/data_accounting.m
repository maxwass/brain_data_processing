%% Accounting on which files we have (or need) for which patients
clear; clc;

%% SCS. Given by Zhengwu -> Yang -> Max.
% all_id (1065x1 double)
% loaded_tensor_sub (4-D double)
load('HCP_subcortical_CMData_desikan.mat')
SCs = squeeze(loaded_tensor_sub(:,:,1,:)); % (87x87x1065 double)
subject_list_sc = int64(all_id);
clear('mode3', 'loaded_bin_network_sub','loaded_tensor_sub', 'all_id');


%% OLD correlation fcs. Given by Zhengwu->Yang->Max. Patients w/ fc are subset of those w/ sc
% variables in 'desikan_fc_all.mat'
%  all_fc (1x1065 cell): 
%    ->missing_data (empty cells) @ [239,297,351,387,639,870,1064]
%  ChosenROI_cortical (1x68 double)
%  ChosenROI_subcortical (1x68 double)
%  subjList (1065x1 double)
load('desikan_fc_all.mat')
subject_list_fc = int64(subjList); %subject_list_fc == subject_list_sc
FCs = reshape(cell2mat(all_fc),87,87,[]); %(87x87x1058 double)
clear('all_fc', 'subjList')

%which patients are missing in FCs relative to SCs. If we can find FC's for
%these patients...more data samples
missing_fcs_index = [239,297,351,387,639,870,1064]; %indexes of patients WITH scs and WITHOUT (old) fcs
missing_fcs_ids   = subject_list_sc(missing_fcs_index);


%% new FCs (downloaded from HCP_1200 server and did local computation)
% variables in 'hcp_1200_subject_list.mat'
%  hcp1200_subject_list (1113x1 double)
load hcp_1200_subject_list.mat %subject list (patient ids)
subject_list_hcp1200 = int64(double(string(hcp1200_subject_list)));
clear('hcp1200_subject_list')

%check that provided scs patients are all in the HCP_1200 dataset
sc_subset_hcp = all(ismember(subject_list_sc, subject_list_hcp1200)); %True


% variables in 'subjects_missing_data.mat'
%  missing_LR (1x31 double)
%  missing_RL (1x19 couble)
%  missing_LR_and_RL (1x17 double)
%  type -  (char array) which fmri files {'msmall', 'no_msmall'}
load subjects_missing_data.mat


%% Creating non-intersecting sets:
% missing_only_LR   = patients in HCP_1200 only missing 1 LR scan (not both)
% missing_only_RL   = patients in HCP_1200 only missing 1 RL scan (not both)
% missing_LR_and_RL = patients in HCP_1200 only missing both LR and RL scans (not only 1)

missing_only_LR = setdiff(missing_LR,missing_LR_and_RL); % 14
missing_only_RL = setdiff(missing_RL,missing_LR_and_RL); % 2
missing_LR_and_RL = intersect(missing_LR, missing_RL); %17

%% which patients are missing fcs?
missing_one_fc = union(missing_only_LR, missing_only_RL);  %patients missing exactly one scan
missing_any_fc = union(missing_one_fc, missing_LR_and_RL); %patients missing one or both scans

% exist_fc = exist_both_fc U exist_exactly_one_fc
exist_both_fc = setdiff(subject_list_hcp1200, missing_any_fc); %patients with both scans
exist_exactly_one_fc = union(missing_only_LR, missing_only_RL); %patients with exactly 1 scan
exist_any_fc = setdiff(subject_list_hcp1200, missing_LR_and_RL); %patients with at least one fc


%% What data do we need for more training samples?
exist_sc = subject_list_sc;
missing_sc = setdiff(subject_list_hcp1200, exist_sc); %no fc anyway, so low priority

% create disjoint sets of patients
% exist_any_fc_and_sc = {exist_exactly_one_fc_and_sc}U{exist_both_fc_and_sc}
exist_exactly_one_fc_and_sc = intersect(exist_exactly_one_fc, exist_sc);
exist_both_fc_and_sc = intersect(exist_both_fc, exist_sc); 
exist_any_fc_and_sc = intersect(exist_any_fc, exist_sc);

% For training set, we need an sc and (at least one) fc
path2folder = '/Users/maxwasserman/Documents/MATLAB/brain_data_preprocess/data';
%save([path2folder '/fc_and_sc_sets.mat'], 'exist_any_fc_and_sc', 'exist_both_fc_and_sc');
fprintf("%d patients HAVE FC and HAVE SC (= Training Set size)\n", length(exist_any_fc_and_sc));

%high priority patients - already have fc, just need sc
exist_any_fc_missing_sc = intersect(exist_any_fc, missing_sc);
fprintf("%d patients HAVE FC but DONT HAVE SC\n", length(exist_any_fc_missing_sc));

%high priority patients - already have sc, just need fc
% only one patient that has an sc does NOT have an fc (179548)
missing_fc_exist_sc = intersect(missing_LR_and_RL, exist_sc);
fprintf("%d patients DONT HAVE FC but HAVE SC\n", length(missing_fc_exist_sc));

%lowest priority patients - dont have sc OR fc
missing_fc_missing_sc = intersect(missing_LR_and_RL, missing_sc);
fprintf("%d patients DONT HAVE FC and DONT HAVE SC\n", length(missing_fc_missing_sc));


