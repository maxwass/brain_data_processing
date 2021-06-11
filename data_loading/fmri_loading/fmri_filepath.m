function [path2fmri] = fmri_filepath(rawdatafolder, atlas, tasktype, subject, scan_dir)

%saved as rfMRI_REST* vs REST*
if ~contains(tasktype,'rfMRI_')
    tasktype=strcat("rfMRI_", tasktype);
end
if isnumeric(subject)
    subject = num2str(subject);
end

%% which ROI's are we analyzing? and with which atlas?
if contains(atlas, 'desikan', 'IgnoreCase', true)
    chosen_roi = load('data/desikan_roi_zhengwu.mat');
    %num_chosen_roi = 87;
    num_chosen_roi = length(chosen_roi.cortical) + length(chosen_roi.subcortical); %68 in desikan cortical +  19 subcortical GM (gray matter) structures
else
    error('Atlas %s not supported yet', atlas);
end
    

%% load fMRI data
%path2fmri = [rawdatafolder '/' subject '/' tasktype '_' scan_dir '_Atlas_hp2000_clean.dtseries.nii'];
filename  = strcat(tasktype, '_', scan_dir, '_Atlas_hp2000_clean.dtseries.nii');
path2fmri = char(fullfile(rawdatafolder,subject, filename));

end

