%% zip together sc's and fc's
% 1) check that patient id are equal
% 2) check that (sub)cortical indices are identical
% 3) resave with combined


clear

fc_data_folder   = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_desikan_cortical_subcortical';


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

fcs_cov_lr = zeros(N,N,num_smaples);
fcs_cov_rl = zeros(N,N,num_smaples);

scs = zeros(N,N,num_smaples);
scs_cortical = zeros(N,N,num_smaples);

%%PROBLEM: SCS and FCS LABELED DIFFERENTLY - change this



%combine sc's with their corresponding fcs
for i_index = 1:num_samples
    subject = exist_fc_and_sc(i_index);
    
    %load relevent fcs from .mat file in data folder
    subject_file = sprintf('%s/%s.mat',fc_data_folder, num2str(subject));
    fcs_and_metadata =  matfile(subject_file);
    
    %TODO: check for existance, and pick out (and potentially cast) lr, rl,
    %...
    
    %TODO
    %check to see if fc_corr = mean(fc_corr_lr, fc_corr_rl) is close to fc
    %given by Yang
    
    
    
    %pick out relevant sc in SCs
    if ~ismember(subject, subject_list_sc)
        error('Subject %d NOT in subject_list_sc',subject);
    end
    index_in_subject_list_sc = int64(find(subject_list_sc == subject));
    raw_sc = SCs(:,:,index_in_subject_list_sc);
    [sc, sc_cortical] = transform_sc(raw_sc);
    
    scs(:,:,i_index) = sc;
    scs_cortical(:,:,i_index) = sc_cortical;
    
    
    
    
end



%TODO: 1114-1065 = 49 BUT only 48 in_hcp_not_sc?
%figure out how many have at least one LR/RL and sc

%now lets figure out how many have at least one fMRI scan (LR,RL) and a sc
no_fc = union(missing_LR, missing_RL);
no_sc = setdiff(subject_list_hcp1200, subject_list_sc); %sc IS subset of hcp
no_fc_no_sc = intersect(no_fc, no_sc);

has_fc = setdiff(subject_list_hcp1200, no_fc);
has_sc = subject_list_sc;
has_fc_and_sc = intersect(has_fc, has_sc);


%for i_index = 1:length(hc1200_subject_list)
    %pull out first subject
    %does this subject have an sc?
    %does this subject have an LR?
    %does this subject have an RL?
    
    %If no sc OR neither LR/RL -> exclude

%    disp('NOT YET IMPL')
%end



