function [fc_cov] = load_clean_compute_cov(atlas, fmri_path, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical)%,datafolder, data_write_folder, ntotal, ncortex)
    
    ncortex = 64984;
    ntotal = 96854; %total number of voxels?

    nROIChosenROI = length(ChosenROI_cortical) + length(ChosenROI_subcortical); %68 in desikan cortical +  19 subcortical GM (gray matter) structures

    
    %% load the subject specific label
    if(strcmp(atlas,"desikan"))
        subjlab = ft_read_cifti([datafolder '/' subject '/' subject '.aparc.32k_fs_LR.dlabel.nii']);
        segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc']);
    elseif(strcmp(atlas,"destrieux"))
        subjlab = ft_read_cifti([datafolder '/' subject '/' subject '.aparc.a2009s.32k_fs_LR.dlabel.nii']);
        segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc_a2009s']);
    else
        error("Atlas " + atlas + " not found. Use Desikan or destrieux.")
    end
    
    
    %% load data
    if(strcmp(tasktype,'rfMRI_REST1'))
        %resting fMRI data
        subj_tdata = ft_read_cifti(fmri_path); %LR scans left to right
    else
        error('only REST1 handled at the moment')
    end

    %% initialize data structure dtseries ('data time series'?) for storing observations
    nTR = size(subj_tdata.dtseries,2);
    dtseries = zeros(nROIChosenROI,nTR);
    disp(['Subject ' subject ' ' tasktype ' ' name ' has ' num2str(nTR) ' timepoints']);
    
    %% loop over each brain region and insert all observations for region 
    %   into row of dtseries
    
    %cortical structure
    for j_index=1:length(ChosenROI_cortical)
        r = ChosenROI_cortical(j_index);
        indices = logical([(segmentation_atlas==r); zeros(ntotal-ncortex,1)]);
        dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1); % averaging activity value in particular region
    end

    % subcortical structure:
    for j=1:length(ChosenROI_subcortical)
        r = ChosenROI_subcortical(j);
        j_index = j_index + 1; % use old int index from above, and increment
        indices = logical(subj_tdata_LR1.brainstructure==r);
        dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
    end
    
    % covariance computed using *all* observations (no windowing!)
    fc_cov = cov(dtseries');
end

%%%%%%%%%%%%%%%%%%%%%%5
% Zhengwu Zhang
% 03/22/2017
% Calculate The Covariance Trajectory
% This script only supports two parcellations: (1) desikan, (2) destrieux;

clear all;

%% setup
%datafolder = '/Users/maxwasserman/Documents/MATLAB/brain_data_preprocess/brain_data';
datafolder = '/Volumes/Elements/brain_data';
tasktype='rfMRI_REST1';
%tasktype = 'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING';
%'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';

data_write_folder = '/Users/maxwasserman/Desktop/geom_dl/data/brain_data/fcs_matlabs';

%external library for reading .nii data
addpath('/Users/maxwasserman/Documents/MATLAB/brain_data_preprocess/Correlation_Calculation/cifti_matlab')

%TODO: ensure subject id's for structural graphs MATCH these subject id's

%load subject list (list of strings of subject ids)
load Allprocessedid.mat


%% Decide on atlas
atlas = "desikan"; %"destrieux"
if(strcmp(atlas,"desikan"))                 %missing 4 and 39 are corpus collusum
    ChosenROI_cortical    =   setxor(0:70,[0, 4, 39]); % og task ones [28,29,45,68,69,102,103,119,142,143];
    ChosenROI_subcortical =   []; %setxor(1:21,[1, 2]);
elseif(strcmp(atlas,"destrieux"))           %missing 42 and 117 are corpus collusum
    ChosenROI_cortical    =   setxor(0:150,[0, 42, 117]);
    ChosenROI_subcortical = []; %setxor(1:21,[1, 2])
else
    error("Atlas " + atlas + " not found. Use Desikan or destrieux.")
end

%  Deskian: cortical:    everything (up to 70) but zero and missing (4,39) 
%           subcortical: all but 1,2 (may have to reorder the left, right, left, right sequence
%                         if want all left first, then right)


%useful code:
% unique(subjlab.x100206_aparc_a2009s(~isnan(subjlab.x100206_aparc_a2009s)))
%  - run code, break point, open subjlab and subj_tdata_somethnig and look
%  at fields. brain_strcture and brain_structure_label??
%  - all labels for Destrieux
%  - use setdiff to find missing labels (corpus collosum)
% 
% Destreuix: cortical:  everything but 0 and missing (42, 117)
%            subcortical: -all but 1,2 (these correspond to Left and Right hemisphere) 
%                         -(may have to reorder the left, right, left, right sequence
%                            if want all left first, then right)
%                         - same as desikan

nROIChosenROI = length(ChosenROI_cortical) + length(ChosenROI_subcortical); %68 in desikan cortical +  19 subcortical GM (gray matter) structures

ncortex = 64984;
ntotal = 96854; %total number of voxels?


%% computation
nSubject = size(subjList); %loaded from Allprocessedid.mat
nSubject = nSubject(1);

tic;
% memory issues with parfor
for i_index=1:nSubject

subject = subjList(i_index,:);


%% finding which files this patient has
path_to_LR1 = [datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
path_to_RL1 = [datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];

has_LR1 = isfile(path_to_LR1);
has_RL1 = isfile(path_to_RL1);

%%%%%%% start func


%% load the subject specific label
if(strcmp(atlas,"desikan"))
    subjlab = ft_read_cifti([datafolder '/' subject '/' subject '.aparc.32k_fs_LR.dlabel.nii']);
    segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc']);
elseif(strcmp(atlas,"destrieux"))
    subjlab = ft_read_cifti([datafolder '/' subject '/' subject '.aparc.a2009s.32k_fs_LR.dlabel.nii']);
    segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc_a2009s']);
else
    error("Atlas " + atlas + " not found. Use Desikan or destrieux.")
end


%% reading the data
[subj_tdata_LR1, subj_tdata_RL1] = deal(NaN, NaN); %one line initialization

if(strcmp(tasktype,'rfMRI_REST1'))
    %resting fMRI data
    if has_LR1
        subj_tdata_LR1 = ft_read_cifti(path_to_LR1); %LR scans left to right
    end
    if has_LR2
        subj_tdata_RL1 = ft_read_cifti(path_to_RL1);
    end
%else
%    %task fMRI data
%    subj_tdata_LR1 = ft_read_cifti([datafolder '/' subject '/' tasktype '_LR_Atlas_MSMAll.dtseries.nii']);
%    subj_tdata_RL1 = ft_read_cifti([datafolder '/' subject '/' tasktype '_RL_Atlas_MSMAll.dtseries.nii']);
end

%% use atlas to compute average functional activity in each defined region

% Zhengwu notes
% Calculate correlations for desikan + subcortical structures
% (note: here LR stands for left and right hemispheres, not phase
% encoding.)

% number of time point
if has_LR1
    nTR_LR = size(subj_tdata_LR1.dtseries,2);
    disp(['Subject ' subject ' ' tasktype ' LR1 has ' num2str(nTR_LR) ' timepoints']);
    dtseriesLR = zeros(nROIChosenROI,nTR_LR);
end
if has_RL1:
    nTR_RL = size(subj_tdata_RL1.dtseries,2);
    disp(['Subject ' subject ' ' tasktype ' RL1 has ' num2str(nTR_RL) ' timepoints']);
    dtseriesRL = zeros(nROIChosenROI,nTR_RL);
end
%display(['Subject ' subject ' ' tasktype ' LR2 has ' num2str(nTR_LR) ' timepoints']);


%cortical structure
for j_index=1:length(ChosenROI_cortical)
    r = ChosenROI_cortical(j_index);
    indices = logical([(segmentation_atlas==r); zeros(ntotal-ncortex,1)]);
    if has_LR1
        dtseriesLR(j_index,:) = mean(subj_tdata_LR1.dtseries(indices,:),1); % averaging activity value in particular region
    end
    if has_RL1
        dtseriesRL(j_index,:) = mean(subj_tdata_RL1.dtseries(indices,:),1);
    end
end

% subcortical structure:
for j=1:length(ChosenROI_subcortical)
    r = ChosenROI_subcortical(j);
    j_index = j_index + 1; % use old int index from above, and increment
    indices = logical(subj_tdata_LR1.brainstructure==r);
    if has_LR1
        dtseriesLR(j_index,:) = mean(subj_tdata_LR1.dtseries(indices,:),1);
    end
    if has_RL1
        dtseriesRL(j_index,:) = mean(subj_tdata_RL1.dtseries(indices,:),1);
    end
end

%% compute covariance matrix
if has_LR1
    cov_LR1 = cov(dtseriesLR');
end
if has_RL1
    cov_RL1 = cov(dtseriesRL');
end

%%%% end func

%% save as we go
disp("Saving data...");
filename = [data_write_folder,'/',subject,'.mat'];
save filename cov_LR1 cov_RL1 tasktype ChosenROI_cortical ChosenROI_subcortical


%if only one, use one
%if both, average or fancy resampling



%% Covariance trajectory is separate problem. Standard solution to find
%   functional connectivity is simply compute covariance simultaneously for
%   all observtions.



%calculate time serior covariance trajectory
wind_start_idx = 1;
wind_end_idx = wind_start_idx + windowsize;

idx = 0; 
while(wind_end_idx < nTR_LR)
    idx = idx + 1;
    %TODO: How are these combined to create final covariance matrix? Whats the
    %temporal averaging alg? Simple average? Verify that is will provide
    %covariance matrix? Where is correlation computed?
    cov_traj_LR(idx,:,:) = cov(dtseriesLR(:,wind_start_idx:wind_end_idx)'); %TODO: check if this is covariance***
    cov_traj_RL(idx,:,:) = cov(dtseriesRL(:,wind_start_idx:wind_end_idx)');
    wind_start_idx = wind_start_idx + movesize;
    wind_end_idx = wind_end_idx + movesize;
end

allcov_traj_LR{i_index} = cov_traj_LR;
allcov_traj_RL{i_index} = cov_traj_RL;

end
toc;

%save the data
eval(['save allcov_traj_' tasktype winds ' allcov_traj_LR allcov_traj_RL ChosenROI_cortical ChosenROI_subcortical windowsize movesize'])