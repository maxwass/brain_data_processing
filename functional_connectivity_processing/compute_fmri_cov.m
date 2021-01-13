function [fc_cov, fc_corr] = compute_fmri_cov(atlas, fmri_path, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical)
%% inputs:
% atlas: str in {'desikan', 'destrieux'}
% fmri_path: str to fMRI .nii file
% datafolder: str path to all brain data (atlas', fmri's, etc). Currently
%               on external hd
% subject: str subject id
% tasktype: str 'rfMRI_REST1'. Possibly include REST2 in future
% name: str in {LR1, RL1}
% ChosenROI_cortical: [int]
% ChosenROI_subcortical: [int]
%% outputs:
% fc_cov = covariance matrix computed over *all* time points


%% constants
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
    
    
%% load fMRI data
if(strcmp(tasktype,'rfMRI_REST1'))
    %resting fMRI data
    subj_tdata = ft_read_cifti(fmri_path); %LR scans left to right...~1 GB memory
else
    error('only REST1 handled at the moment')
end

%% initialize data structure dtseries ('data time series'?) for storing observations
nTR = size(subj_tdata.dtseries,2);
dtseries = zeros(nROIChosenROI,nTR);
%disp(['Subject ' subject ' ' tasktype ' ' name ' has ' num2str(nTR) ' timepoints']);
disp(['   ' name ' has ' num2str(nTR) ' timepoints']);

%% loop over each brain region and insert all observations for region into row of dtseries
% Zhengwu notes
% Calculate correlations for desikan + subcortical structures
% (note: here LR stands for left and right hemispheres, not phase
% encoding.) ???????

%Note that cortical indices will be FIRST 68 in matrix. Subcortical will be
%LAST 19 - CHANGING

total_ROIs = length(ChosenROI_subcortical) + length(ChosenROI_cortical);


for j_index = 1:total_ROIs
    if j_index <= 19 %subcortical
        r = ChosenROI_subcortical(j_index);
        indices = logical(subj_tdata.brainstructure==r);
        dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
        
    elseif (20 <= j_index) && (j_index <= 87) %cortical
        j = j_index-length(ChosenROI_subcortical);
        r = ChosenROI_cortical(j);
        indices = logical([(segmentation_atlas==r); zeros(ntotal-ncortex,1)]);
        dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
        
    else
        error('more indices %d than ROI %d',j_index ,total_ROIs);
    end
end
    

%{
% subcortical structure:
for j=1:length(ChosenROI_subcortical)
    r = ChosenROI_subcortical(j);
    j_index = j_index + 1; % use old int index from above, and increment
    indices = logical(subj_tdata.brainstructure==r);
    dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
end

%cortical structure
for j_index=1:length(ChosenROI_cortical)
    r = ChosenROI_cortical(j_index);
    indices = logical([(segmentation_atlas==r); zeros(ntotal-ncortex,1)]);
    %place average activity in region r into j_index row
    dtseries(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
end
%}

    
%% covariance computed using *all* observations (no windowing!)
fc_cov = cov(dtseries');
fc_corr = corrcov(fc_cov);

end