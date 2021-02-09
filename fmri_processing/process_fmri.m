function [dtseries] = process_fmri(atlas, path2fmri, subject, rawdatafolder, chosen_roi)
%% inputs:
% atlas: str in {'desikan', 'destrieux'}
% fmri_path: str to fMRI .nii file
% rawdatafolder: str path to all brain data (atlas', fmri's, etc). Currently
%               on external hd
% subject: str subject id
% chosen_roi is a struct of ROI's from cortical and subcortical regions
%   chosen_roi.cortical:    [int]
%   chosen_roi.subcortical: [int]
%   chosen_roi.subcortical_first: logical
%% outputs:
% dtseries: num_roi x num_obvs = vector observations at each time point

%% first check if we have cached this processing
[cached, cached_filename] = is_cached(subject, atlas, chosen_roi, path2fmri);

if cached
    dtseries = load(cached_filename, 'dtseries').dtseries;
    return
end


%% constants
ncortex = 64984;
ntotal = 96854; %total number of voxels?


num_chosen_roi = length(chosen_roi.cortical) + length(chosen_roi.subcortical); %68 in desikan cortical +  19 subcortical GM (gray matter) structures

    
%% load the subject specific label (atlas)
[segmentation_atlas] = load_atlas(atlas, subject, rawdatafolder);
    
%% load fMRI data
[fmri_data_struct] = load_fmri(path2fmri);

%% initialize data structure dtseries ('data time series'?) for storing observations
nTR = size(fmri_data_struct.dtseries,2); %number of observations
dtseries = double(zeros(num_chosen_roi,nTR));
%disp(['   ' name ' has ' num2str(nTR) ' timepoints']);


%% loop over each brain region and insert all observations for region into row of dtseries
total_ROIs = length(chosen_roi.subcortical) + length(chosen_roi.cortical);

startprocess = tic;
for roi_index = 1:total_ROIs
    
    if roi_index <= 19 %subcortical
        roi = chosen_roi.subcortical(roi_index);
        
        %create mask for brain region r
        indices = logical(fmri_data_struct.brainstructure==roi);
    
    elseif (20 <= roi_index) && (roi_index <= 87) %cortical
        %roi_index should be 1 when we first enter this loop -> subtract
        %off previous increments from above subcortical case.
        roi = chosen_roi.cortical( roi_index - length(chosen_roi.subcortical) );
        
        %create mask for brain region r
        indices = logical([(segmentation_atlas==roi); zeros(ntotal-ncortex,1)]);
    
    else
        error('more indices %d than ROI %d',roi_index ,total_ROIs);
    end
    
    %average activity value in masked region => functional signal
    dtseries(roi_index,:) = mean(fmri_data_struct.dtseries(indices,:),1);
end

stopprocess = toc(startprocess);
%fprintf('time for processing data: %.2f\n', stopprocess);

    
end

