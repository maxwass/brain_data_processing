%% zip together sc's and fc's for dataset of (sc, fc) pairs

clear; clc;

%'/fcs_desikan_subcortical_cortical' is the output of compute_fcs.m.
% It is a directory of processed fcs. One .mat file per patient.
%This directory is NOT included in repo: 245 MB. To download this directory:
%https://drive.google.com/drive/folders/1LQocAEeEI02704KrmIior4HEcXx8KxSk?usp=sharing
fc_data_folder   = '~/Desktop/geom_dl/data/brain_data/fcs_desikan_subcortical_cortical';


%% load in relevant files
%scs
% subject_list (1065x1 int64)
% scs (87x87x1065 double)
sc_file = load('data/scs_desikan.mat');
SCs = sc_file.scs; 
subject_list_sc = sc_file.subject_list;

% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
load('data/hcp_1200_subject_list.mat'); %subjList
subject_list_hcp1200 = int64(double(string(hcp1200_subject_list))); % (1113x1 double)
clear('hcp1200_subject_list')

% exist_any_fc_and_sc (1064x1 int64)  subjects with fc(s) and sc
% exist_both_fc_and_sc (1051x1 int64)  subjects with both fcs and sc
load('data/fc_and_sc_sets.mat');
missing_one_scan = setdiff(exist_any_fc_and_sc, exist_both_fc_and_sc);

%switch this if only like to consider patients with both functional scans
dataset = exist_any_fc_and_sc; %exist_both_fc_and_sc


%save full cov (87x87) and remove subcortical later
N = 87;
N_cortical = 68;
num_samples = length(dataset);
final_subject_list = dataset;

tensor_shape = [N, N, num_samples];

fcs_cov_lr = zeros(tensor_shape);
fcs_cov_rl = zeros(tensor_shape);
fcs_cov_mean = zeros(tensor_shape);

fcs_corr_lr = zeros(tensor_shape);
fcs_corr_rl = zeros(tensor_shape);
fcs_corr_mean = zeros(tensor_shape);

raw_scs = zeros(tensor_shape);
transform_scs = zeros(tensor_shape);

%combine sc's with their corresponding fcs
for i_index = 1:num_samples
    subject = dataset(i_index);
    
    %load relevent fcs and metadata from .mat file in data folder
    fprintf('loading %d: %d fcs...\n', i_index, subject);
    subject_file = sprintf('%s/%s.mat',fc_data_folder, num2str(subject));
    
    %struct with variables from .mat file
    fcs_and_metadata =  matfile(subject_file);

    fc_cov_lr = fcs_and_metadata.fc_cov_lr;
    fc_cov_rl = fcs_and_metadata.fc_cov_rl;
    
    fc_corr_lr = fcs_and_metadata.fc_corr_lr;
    fc_corr_rl = fcs_and_metadata.fc_corr_rl;
    
    %fc_cov_* saved as nan matrix if doesnt exist
    has_lr = ~any(any(isnan(fc_cov_lr)));
    has_rl = ~any(any(isnan(fc_cov_rl)));
    
    if has_lr && has_rl
        if ismember(subject,missing_one_scan)
            error('Subject %d is missing a scan but in wrong part!',subject)
        end
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
    
    %index order between final_subject_list and subject_list_sc can be
    %different
    %find index of subject in the subject list (which indexes the SCs)
    index_in_subject_list_sc = int64(find(subject_list_sc == subject));
    raw_sc = SCs(:,:,index_in_subject_list_sc); %between 0 and 19843.00
    raw_scs(:,:,i_index) = raw_sc;
    
    %see Yang code/transform_sc.m... sc seems to have subcrotical for first 19
    %[sc, sc_cortical] = transform_sc(raw_sc);
    transform_scs(:,:,i_index) = log(raw_sc+raw_sc'+1);

    
    clear fcs_and_metadata
    
end

%% prepare metadata
%metadata is same for all subject, use last patient
fcs_and_metadata      =  matfile(subject_file);
ChosenROI_cortical    = fcs_and_metadata.ChosenROI_cortical;
ChosenROI_subcortical = fcs_and_metadata.ChosenROI_subcortical;
subcortical_first     = fcs_and_metadata.subcortical_first;
atlas    = char(fcs_and_metadata.atlas);
tasktype = char(fcs_and_metadata.tasktype); %must use char for python loading
clear fcs_and_metadata


%% save data and metadata


%OLD correlation fcs file 
% ChosenROI_cortical (1x68 double)
% ChosenROI_subcortical (1x68 double)
old_fc_file = load('data/correlations_desikan_old.mat');
ChosenROI_cortical = old_fc_file.ChosenROI_cortical;
ChosenROI_subcortical = old_fc_file.ChosenROI_subcortical;
clear('old_fc_file')


%metadata:ChosenROI_cortical, Chosen_RIO_subcortical, subcortical_first, 
% atlas, tasktype, final_subject_list

%data
% final_subject_list - int array of subject id's
% fcs_cov_lr/_rl/mean - double tensor 87x87xnum_samples (1064)
% fcs_corr_lr/rl/mean - double tensor 87x87xnum_samples (1064)
% raw_scs, transform_scs - int tensor/double tensor 87x87xnum_samples (1064)


%note the invariant enforced:
% INVARIANT: ith entry in final_subject_list corresponds to the ith entry in
%             all data tensors fcs_* and *_scs

save('data/brain_dataset_sc_fc_pairs.mat',... 
    'fcs_cov_lr', 'fcs_cov_rl',  'fcs_cov_mean',...
    'fcs_corr_lr', 'fcs_corr_rl', 'fcs_corr_mean',...
    'raw_scs', 'transform_scs',...
    'final_subject_list',...
    'ChosenROI_cortical', 'ChosenROI_subcortical', 'subcortical_first',...
    'atlas','tasktype');
