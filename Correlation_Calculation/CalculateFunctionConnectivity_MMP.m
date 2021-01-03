%%%%%%%%%%%%%%%%%%%%%%5
% Zhengwu Zhang
% 03/22/2017
% Calculate The Covariance Trajectory
% This script only supports two parcellations: (1) desikan, (2) destrieux;

addpath('./cifti_matlab')

%load subject
load Allprocessedid.mat
subjList = allprocessid(1:500);
nSubject = size(subjList);

%ROI index for caculating the covariance trajectory

%those rois are for the motion tasks
ChosenROI_cortical = [28,29,45,68,69,102,103,119,142,143];
ChosenROI_subcortical = [];


nROIChosenROI = length(ChosenROI_cortical) + length(ChosenROI_subcortical); %68 in desikan plus and 19 subcortical GM structures

ncortex = 64984;
ntotal = 96854;

%parameters for calculating the covariance trajectory
windowsize = 30;
movesize = 10;

winds = sprintf('_wind%d_move%d_',windowsize,movesize);

%datafolder = '/media/samsi/Extend2/HCP_fMRI/';
%datafolder = 'F:\HCP_fMRI';
datafolder = './';
tasktype='rfMRI_REST1';
%tasktype = 'tfMRI_GAMBLING';
%tasktype = 'tfMRI_MOTOR';
%tasktype = 'tfMRI_GAMBLING';
%tasktype = 'tfMRI_SOCIAL';
%tasktype = 'tfMRI_LANGUAGE';

tic;
% memory issues with parfor
for i_index=1:nSubject

subject = num2str(subjList(i_index));

%load the subject specific label
%desikan atlas, values 4 and 39 are not in atlas -- correspond to corpus collosum
%subjlab = ft_read_cifti([datafolder subject '/' subject '.aparc.32k_fs_LR.dlabel.nii']);
%segmentation_atlas = eval(['subjlab.x' num2s?tr(subject) '_aparc']);

%destrieux atlas, values 4 and 39 are not in atlas -- correspond to corpus collosum
subjlab = ft_read_cifti([datafolder subject '/' subject '.aparc.a2009s.32k_fs_LR.dlabel.nii']);
segmentation_atlas = eval(['subjlab.x' num2str(subject) '_aparc_a2009s']);

%reading the data
if(strcmp(tasktype,'rfMRI_REST1'))
    %resting fMRI data
    subj_tdata_LR1 = ft_read_cifti([datafolder subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii']);
    subj_tdata_RL1 = ft_read_cifti([datafolder subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii']);
else
    %task fMRI data
    subj_tdata_LR1 = ft_read_cifti([datafolder subject '/' tasktype '_LR_Atlas_MSMAll.dtseries.nii']);
    subj_tdata_RL1 = ft_read_cifti([datafolder subject '/' tasktype '_RL_Atlas_MSMAll.dtseries.nii']);
end



% Calculate correlations for desikan + subcortical structures
% (note: here LR stands for left and right hemispheres, not phase
% encoding.)

% number of time point
nTR_LR = size(subj_tdata_LR1.dtseries,2);
nTR_RL = size(subj_tdata_RL1.dtseries,2);

display(['Subject ' subject ' ' tasktype ' LR1 has ' num2str(nTR_LR) ' timepoints']);
display(['Subject ' subject ' ' tasktype ' LR2 has ' num2str(nTR_LR) ' timepoints']);

dtseriesLR = zeros(nROIChosenROI,nTR_LR);
dtseriesRL = zeros(nROIChosenROI,nTR_RL);

%cortical structure
% values 4 and 39 are not in atlas -- correspond to corpus collosum
for j_index=1:length(ChosenROI_cortical)
    r = ChosenROI_cortical(j_index);
    indices = logical([(segmentation_atlas==r); zeros(ntotal-ncortex,1)]);
    dtseriesLR(j_index,:) = mean(subj_tdata_LR1.dtseries(indices,:),1);
    dtseriesRL(j_index,:) = mean(subj_tdata_RL1.dtseries(indices,:),1);
end

% subcortical structure:
for j=1:length(ChosenROI_subcortical)
    r = ChosenROI_subcortical(j);
    j_index = j_index + 1;
    indices = logical(subj_tdata_LR1.brainstructure==r);
    dtseriesLR(j_index,:) = mean(subj_tdata_LR1.dtseries(indices,:),1);
    dtseriesRL(j_index,:) = mean(subj_tdata_RL1.dtseries(indices,:),1);    
end

%calculate time serior covariance trajectory
wind_start_idx = 1;
wind_end_idx = wind_start_idx + windowsize;

idx = 0; 
while(wind_end_idx < nTR_LR)
    idx = idx + 1;
    cov_traj_LR(idx,:,:) = cov(dtseriesLR(:,wind_start_idx:wind_end_idx)');
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