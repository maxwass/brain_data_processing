clear
close all

%move to matlab folder, currently in geom_dl folder
path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo)); %recursively adds all repo files/folders

atlas = char('desikan');
tasktype = char('rfMRI_REST1');
subcortical_first_fc = false;
subcortical_first_sc = true;


load('data/desikan_fc_all.mat') %f1058 unctional graphs
load('data/scs.mat') %1065 structural graphs
subject_list = int64(all_id); %1065x1
SCs = squeeze(loaded_tensor_sub(:,:,1,:));

missing_data = [239,297,351,387,639,870,1064];

num_patients = 1065;%-length(missing_data);
N_full = 87;  N_cortical = 68;
tensor_full_shape = [N_full, N_full, num_patients];
tensor_cortical_shape = [N_cortical, N_cortical, num_patients];

%initialize matrice which are to be saved in .mat file
% each row is a resize (1,N^2) matrix
%func_corr_full = zeros(num_patients,N_full^2);
%func_corr_cortical  = zeros(num_patients,N_cortical^2);
fcs_corr_full = zeros(tensor_full_shape);
fcs_corr_cortical  = zeros(tensor_cortical_shape);


%struct_full     = zeros(num_patients,N_full^2);
%struct_cortical      = zeros(num_patients,N_cortical^2);
%struct_sparsity_full = zeros(1, num_patients);
%struct_sparsity_cortical = zeros(1, num_patients);
scs_full     = zeros(tensor_full_shape);
scs_cortical = zeros(tensor_cortical_shape);

scs_sparsity_full     = zeros(1,num_patients);
scs_sparsity_cortical = zeros(1,num_patients);



for i = 1 : 1065
    
    if ismember(i,missing_data)
        continue
    end
    
    %func_corr is a correlation matrix of functional brain activity
    fc = all_fc{i};
    fcs_corr_full(:,:,i)     = fc;
    fcs_corr_cortical(:,:,i) = fc(1:N_cortical,1:N_cortical);
    
    
    %struct_conn is an upper triangular matrix
    % apply transformation to make symmetric and decrease variation between
    % small and large entries.
    % Calculate sparsities of each structural network
    sc = SCs(:,:,i);
    transformed_sc           = log(sc+sc'+1);
    transformed_sc_cortical  = transformed_sc(20:end,20:end);
    scs_full(:,:,i)          = transformed_sc;
    scs_cortical(:,:,i)      = transformed_sc_cortical;
    scs_sparsity_full(i)     = (1/2)*length(find(transformed_sc>0))/(N_full*(N_full-1)/2);
    scs_sparsity_cortical(i) = (1/2)*length(find(transformed_sc_cortical>0))/(N_cortical*(N_cortical-1)/2);
    
    
    % ----save all 4 data structs to .mat file. ---
    % -read in with python script
    % -make a func to do reading in. Call func in hp_search_brain_data.py
    %   file and train.
    
end

%remove rows where patients have missing data
fcs_corr_full(:,:,missing_data) = [];
fcs_corr_cortical(:,:,missing_data) = [];

scs_full(:,:,missing_data) = [];
scs_cortical(:,:,missing_data) = [];
scs_sparsity_full(missing_data) = [];
scs_sparsity_cortical(missing_data) = [];
subject_list(missing_data) = [];

save('old_brain_datataset_corr.mat',... 
    'fcs_corr_full', 'fcs_corr_cortical',...
    'scs_full', 'scs_cortical',...
    'scs_sparsity_full', 'scs_sparsity_cortical',...
    'subject_list',...
    'ChosenROI_cortical', 'ChosenROI_subcortical',...
    'subcortical_first_fc','subcortical_first_sc',...
    'atlas','tasktype');


%save('old_brain_data','func_corr_full','func_corr_sub', 'struct_full', 'struct_sub', 'struct_sparsity_full', 'struct_sparsity_sub');
