% scalar plot testing
close all; clear;

path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo));

%load in first patient which is saved on machine for faster testing
subject = 100206;
load('dtseries_testing.mat','current_dtseries');
atlas = "desikan";
include_subcortical = false;
GSO = 'L';

[num_rois, ~] = size(current_dtseries);
roi_idxs = (20:num_rois);
if include_subcortical
    roi_idxs = 1:num_rois;
end

%% mean center raw signals
dtseries            = current_dtseries(roi_idxs,:);
mean_signal         = mean(dtseries,2);
dtseries_center     = dtseries - mean_signal;

%% compute frequency representation of windowed signals
[signals_freq, GFT, evals] = apply_GFT(dtseries_center,subject, atlas, include_subcortical, GSO);

%% window signal
windowsize = 20;
movesize = 5;

[covs, aw]   = windowed_fcs(dtseries_center, windowsize, movesize);
corrs = corrcov_tensor(covs);
[N, ~, num_windows] = size(covs);

%% load scs from saved tensor file
A = extract_sc(subject, atlas, include_subcortical);
A_upp_tri = triu(A,0);
A_upp_tri_norm = A_upp_tri./norm(A_upp_tri,2); %only consider upper triangular part

ax = axes;
plot_scalars(ax, signals_freq, covs, corrs, A_upp_tri_norm)
set(ax, 'XTick', 1:20:num_windows)

text(10,0,'low','Color','red','FontSize',8)