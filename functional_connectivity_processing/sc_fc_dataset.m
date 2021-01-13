%% zip together sc's and fc's
% 1) check that patient id are equal
% 2) check that (sub)cortical indices are identical
% 3) resave with combined


clear

%fc_data_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_cortical_subcortical';
fc_data_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_subcortical_cortical';

% load in relevant files

% all_id (1065x1 double)
% loaded_tensor_sub (4-D double)
load('HCP_subcortical_CMData_desikan.mat')
clear('mode3', 'loaded_bin_network_sub'); %unused variables
SCs = squeeze(loaded_tensor_sub(:,:,1,:)); % (87x87x1065 double)
subject_list_sc = int64(all_id);

% all_fc (1x1065 cell)
% ChosenROI_cortical (1x68 double)
% ChosenROI_subcortical (1x68 double)
% subjList (1065x1 double)
load('desikan_fc_all.mat')
subject_list_fc = int64(subjList);
FCs = reshape(cell2mat(all_fc),87,87,[]); %(87x87x1058 double)

%which patients are missing in FCs relative to SCs
% 150019 - missing LR
% 160931 - missing LR
% 173233 - missing LR
% 179548 - missing both
% 351938 - missing LR
% 693461 - missing LR
% 995174 - missing LR
missing_fcs_index = [239,297,351,387,639,870,1064]; %indexes of patients with scs missing fcs
missing_fcs_ids   = subject_list_sc(missing_fcs_index);


desikan_data_folder = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan';

load hcp_1200_subject_list.mat %subjList
subject_list_hcp1200 = int64(double(string(hcp1200_subject_list))); % (1113x1 double)
clear(hcp1200_subject_list)


%fcs here are strict subset of scs - outside knowledge - given missing ones

%are scs subset of hcp1200? If so, ven diagrams more easy - all within big
%hcp circle
sc_subset_hcp = all(ismember(subject_list_sc, subject_list_hcp1200)); %True

%are the same ones missing from scs and missing LR/RL/both?
load subjects_missing_data.mat %currently for no_msmaii


missing_only_LR = setdiff(missing_LR,missing_LR_and_RL); % 14
missing_only_RL = setdiff(missing_RL,missing_LR_and_RL); % 2
missing_LR_and_RL = intersect(missing_LR, missing_RL); %17

missing_one_fc = union(missing_only_LR, missing_only_RL); %patients only missing exactly one scan
missing_any_fc = union(missing_one_fc, missing_LR_and_RL); %patients missing any scans

exist_exactly_one_fc = union(missing_only_LR, missing_only_RL);
exist_fc = setdiff(subject_list_hcp1200, missing_LR_and_RL); %patients with at least one fc
exist_both_fc = setdiff(subject_list_hcp1200, missing_any_fc); %patients with both scans


%samples we CAN use - sc and at least one fc
% all but one (179548) that have scs, have at least one fc
exist_exactly_one_fc_and_sc = intersect(exist_exactly_one_fc, subject_list_sc); %sc and exactly one fc
exist_both_fc_and_sc = intersect(exist_both_fc, subject_list_sc); %sc and both fc
exist_fc_and_sc = intersect(exist_fc, subject_list_sc); %sc and at least one fc

%in order to get more data, which scs do we need?
no_sc = setdiff(subject_list_hcp1200, subject_list_sc); %no fc anyway, so low priority
exist_fc_no_sc = intersect(exist_fc, no_sc); %highest priority
no_fc_exist_sc = intersect(missing_LR_and_RL, subject_list_sc); %high priority
no_fc_no_sc = intersect(missing_LR_and_RL, no_sc);



%save full cov (87x87) and remove subcortical later
N = 87;
N_cortical = 68;
num_samples = length(exist_fc_and_sc);
final_subject_list = exist_fc_and_sc;

fcs_cov_lr = zeros(N,N,num_samples);
fcs_cov_rl = zeros(N,N,num_samples);
fcs_cov_mean = zeros(N,N,num_samples);

fcs_corr_lr = zeros(N,N,num_samples);
fcs_corr_rl = zeros(N,N,num_samples);
fcs_corr_mean = zeros(N,N,num_samples);

raw_scs = zeros(N,N,num_samples);
transform_scs = zeros(N,N,num_samples);

%combine sc's with their corresponding fcs
for i_index = 1:num_samples
    subject = exist_fc_and_sc(i_index);
    
    %load relevent fcs and metadata from .mat file in data folder
    fprintf('loading %d: %d fcs...\n', i_index, subject);
    subject_file = sprintf('%s/%s.mat',fc_data_folder, num2str(subject));
    
    %struct with variables from .mat file
    fcs_and_metadata =  matfile(subject_file);
    

    fc_cov_lr = fcs_and_metadata.fc_cov_lr;
    fc_cov_rl = fcs_and_metadata.fc_cov_rl;
    
    fc_corr_lr = fcs_and_metadata.fc_corr_lr;
    fc_corr_rl = fcs_and_metadata.fc_corr_rl;
    
    has_lr = ~any(any(isnan(fc_cov_lr)));
    has_rl = ~any(any(isnan(fc_cov_rl)));

    if has_lr && has_rl
        %cov's
        fcs_cov_lr(:,:,i_index)   = fc_cov_lr;
        fcs_cov_rl(:,:,i_index)   = fc_cov_rl;
        fcs_cov_mean(:,:,i_index) = (fc_cov_lr + fc_cov_rl)/2;
        %corr's
        fcs_corr_lr(:,:,i_index)   = fc_corr_lr;
        fcs_corr_rl(:,:,i_index)   = fc_corr_rl;
        fcs_corr_mean(:,:,i_index) = (fc_corr_lr + fc_corr_rl)/2;
        
    elseif has_lr
        %cov's
        fcs_cov_lr(:,:,i_index)   = fc_cov_lr;
        fcs_cov_rl(:,:,i_index)   = nan(size(fc_cov_lr)); %use existing  lr for size
        fcs_cov_mean(:,:,i_index) = fc_cov_lr;
        %corrs
        fcs_corr_lr(:,:,i_index)   = fc_corr_lr;
        fcs_corr_rl(:,:,i_index)   = nan(size(fc_corr_lr)); %use existing  lr for size
        fcs_corr_mean(:,:,i_index) = fc_corr_lr;
    elseif has_rl
        %cov
        fcs_cov_rl(:,:,i_index)   = fc_cov_rl;
        fcs_cov_lr(:,:,i_index)   = nan(size(fc_cov_rl)); %use existing rl for size
        fcs_cov_mean(:,:,i_index) = fc_cov_rl;
        %corr
        fcs_corr_rl(:,:,i_index)   = fc_corr_rl;
        fcs_corr_lr(:,:,i_index)   = nan(size(fc_corr_rl)); %use existing rl for size
        fcs_corr_mean(:,:,i_index) = fc_corr_rl;
    else
        error('constructing dataset but %d subject (%s) has neither LR or RL fc',i_index, subject)
    end

    
    %pick out relevant sc in SCs
    if ~ismember(subject, subject_list_sc)
        error('Subject %d NOT in subject_list_sc',subject);
    end
    index_in_subject_list_sc = int64(find(subject_list_sc == subject));
    raw_sc = SCs(:,:,index_in_subject_list_sc); %between 0 and 19843.00
    raw_scs(:,:,i_index) = raw_sc;
    
    %see Yang code/transform_sc.m... sc seems to have subcrotical for first 19
    %[sc, sc_cortical] = transform_sc(raw_sc);
    transform_scs(:,:,i_index) = log(raw_sc+raw_sc'+1);

    
    clear fcs_and_metadata
    
end

%metadata is same for all subject, use last patient
fcs_and_metadata =  matfile(subject_file);
ChosenROI_cortical = fcs_and_metadata.ChosenROI_cortical;
ChosenROI_subcortical = fcs_and_metadata.ChosenROI_subcortical;
subcortical_first = fcs_and_metadata.subcortical_first;
atlas = fcs_and_metadata.atlas;
tasktype = string(fcs_and_metadata.tasktype);
clear fcs_and_metadata


%save data and metadata

%metadata:ChosenROI_cortical, Chosen_RIO_subcortical, subcortical_first, 
% atlas, tasktype, final_subject_list

%data
% final_subject_list - int array of subject id's
% fcs_cov_lr/_rl/mean - double tensor 87x87xnum_samples (1064)
% fcs_corr_lr/rl/mean - double tensor 87x87xnum_samples (1064)
% raw_scs, transform_scs - int tensor/double tensor 87x87xnum_samples (1064)


save('brain_dataset_cov.mat',... 
    'fcs_cov_lr', 'fcs_cov_rl',  'fcs_cov_mean',...
    'fcs_corr_lr', 'fcs_corr_rl', 'fcs_corr_mean',...
    'raw_scs', 'transform_scs',...
    'final_subject_list',...
    'ChosenROI_cortical', 'ChosenROI_subcortical', 'subcortical_first',...
    'atlas','tasktype');

