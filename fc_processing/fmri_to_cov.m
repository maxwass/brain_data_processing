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

%loads fmri file and averages activity values in atlas regions
dtseries = load_fmri(atlas, fmri_path, subject, tasktype, datafolder, name, ChosenROI_cortical, ChosenROI_subcortical);
    
%% covariance computed using *all* observations (no windowing!)
fc_cov = cov(dtseries');
fc_corr = corrcov(fc_cov);

end