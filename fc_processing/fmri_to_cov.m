function [fc_cov, fc_corr] = fmri_to_cov(atlas, fmri_path, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical)
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

    
%% load the subject specific label (atlas)
% (note: here LR stands for left and right hemispheres, NOT phase
% encoding - aka directionality of scanning left to right or right to
% left.)
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
    %resting fMRI data - ~0.9 GB memory
    subj_tdata = ft_read_cifti(fmri_path); %LR scans left to right and vica versa
else
    error('only REST1 handled at the moment')
end

%% initialize data structure dtseries ('data time series'?) for storing observations
nTR = size(subj_tdata.dtseries,2); %number of observations
dtseries = zeros(nROIChosenROI,nTR);
disp(['   ' name ' has ' num2str(nTR) ' timepoints']);

%% loop over each brain region and insert all observations for region into row of dtseries
total_ROIs = length(ChosenROI_subcortical) + length(ChosenROI_cortical);

for roi_index = 1:total_ROIs
    if roi_index <= 19 %subcortical
        roi = ChosenROI_subcortical(roi_index);
        %create mask for brain region r
        indices = logical(subj_tdata.brainstructure==roi);
        %average activity value in masked region => functional signal
        pre = mean(subj_tdata.dtseries(indices,:),1);
        
    elseif (20 <= roi_index) && (roi_index <= 87) %cortical
        %roi_index should be 1 when we first enter this loop -> subtract
        %off previous increments from above subcortical case.
        roi = ChosenROI_cortical( roi_index - length(ChosenROI_subcortical) );
        %create mask for brain region r
        indices = logical([(segmentation_atlas==roi); zeros(ntotal-ncortex,1)]);
        %average activity value in masked region => functional signal
        pre = mean(subj_tdata.dtseries(indices,:),1);
        
    else
        error('more indices %d than ROI %d',roi_index ,total_ROIs);
    end
    
    %average activity value in masked region => functional signal
    dtseries(roi_index,:) = mean(subj_tdata.dtseries(indices,:),1);
    %check that pre == 
    if ~ all(pre == dtseries(roi_index,:))
        fprintf('change is invalid!')
    end
    
end
    
dtseries_check = zeros(nROIChosenROI,nTR);

% subcortical structure:
for j_index=1:length(ChosenROI_subcortical)
    roi = ChosenROI_subcortical(j_index);
    indices = logical(subj_tdata.brainstructure==roi);
    dtseries_check(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
end

%cortical structure
for j=1:length(ChosenROI_cortical)
    j_index = j_index + 1; % use old int index from above, and increment
    roi = ChosenROI_cortical(j);
    indices = logical([(segmentation_atlas==roi); zeros(ntotal-ncortex,1)]);
    %place average activity in region r into j_index row
    dtseries_check(j_index,:) = mean(subj_tdata.dtseries(indices,:),1);
end

if ~all(all(dtseries==dtseries_check))
    fprintf('Zhengwu code different from yours!')
end


    
%% covariance computed using *all* observations (no windowing!)
fc_cov = cov(dtseries');
fc_corr = corrcov(fc_cov);

end