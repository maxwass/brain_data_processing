%freq plots
clear; clc; close all;

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

brain_dataset = load('data/brain_dataset_sc_fc_pairs.mat');
subject_list = int2str(brain_dataset.final_subject_list); %all subjects with scs

subject = subject_list(1,:);
path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];

chosen_roi = load('data/desikan_roi', 'roi').roi;
%chosen_roi.subcortical = [];
include_subcortical = false;
GSO = "L";

[windowsize, movesize] = deal(30, 20);

path2fmri = path_to_LR1;


[dtseries]    = process_fmri(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);

%mean center
dtseries = dtseries - mean(dtseries,2);

[signal_windows] = windowed_signals(dtseries, windowsize, movesize);
[fc_covs]     = construct_fcs(signal_windows);
[fc_corr]     = apply_to_tensor_slices(@corrcov, fc_covs);
col_mean      = @(x) mean(x,2);
[ave_signals] = apply_to_tensor_slices(col_mean, signal_windows);

if ~include_subcortical
    cortical_indices = {'20:end','20:end',':'};
    dtseries       = dtseries(20:end, :);
    signal_windows = signal_windows(20:end, :, :);
    fc_covs        = fc_covs(20:end,20:end, :);
    fc_corr        = fc_corr(20:end,20:end, :);
    ave_signals    = ave_signals(20:end, :);
end

%% GFT computation
%[all_signals_freq, ~, ~]       = apply_GFT(dtseries,    subject, atlas, include_subcortical, GSO);
[ave_signals_freq, GFT, evals] = apply_GFT(ave_signals, subject, atlas, include_subcortical, GSO);

%plot scalars!
ax = axes;
%plot_scalars(ax, dtseries, all_signals_freq, covs, corrs, S); %adjust this for all signals!
S = extract_sc(subject, atlas, include_subcortical);
plot_scalars(ax, ave_signals, ave_signals_freq, fc_covs, fc_corr, S);


%{
%% Which signal windows are the most similar? EXCLUDE ZERO FREQ COMPONENT
D   = pdist(signals_freq(2:end,:)', 'euclidean');
D_s = squareform(D);
%h   = heatmap(D_s);

% K-NN similarity graph: find smallest values in each row. Ignore self
% similarity (make it very large)
k=10;
[B,I]=mink(D_s + eye(size(D_s))*100000,k,2);

%Do some clustering to find REST vs THINK windows
Y = tsne(ave_signals_freq_lr(2:end,:)','Algorithm','barneshut','NumPCAComponents',15, 'NumDimensions',3);

%kmeans: https://www.mathworks.com/help/stats/kmeans.html
%add image to each datapoint: https://www.mathworks.com/matlabcentral/answers/311205-how-to-add-images-to-data-points
[idx,C] = kmeans(X,2);

%% Plot freq representation
norms = sqrt(sum(ave_signals_freq.^2, 1)); %should i disregard 0 freq component?

%divide each column by respective norm
signals_freq_norm = bsxfun(@rdivide, ave_signals_freq, norms);

%find max for fixed y axis limits (excluding 0 component)
sig_concat = abs([signals_freq_norm, ave_signals_freq_rl_norm]);
y_max = max(max(sig_concat(2:end,:)));

evals = diag(evals);

num_windows = 10;
t = tiledlayout(2,num_windows);

for idx = 1:num_windows
    nexttile(t)
    f_lr = abs(signals_freq_norm(:,idx));
    s_lr = stem(evals(2:end), f_lr(2:end), 'MarkerEdgeColor','green','color', 'k');
    ylim([0, y_max])
    %freq_mag_lr = ave_signals_freq_lr(:,idx)/norm(ave_signals_freq_lr(2:end,idx));
    %s_lr = stem(evals(2:end), abs(freq_mag_lr(2:end)), 'MarkerEdgeColor','green');
    if idx==1
        ylabel('|freq signal|');
    end
    title(sprintf('LR1: Window %d', idx));
    
    nexttile(num_windows+idx)
    %normalize WIHTOUT 0 freq component?
    f_rl = abs(ave_signals_freq_rl_norm(:,idx));
    s_rl = stem(evals(2:end), f_rl(2:end), 'MarkerEdgeColor','green', 'color', 'k');
    ylim([0, y_max])
    %freq_mag_rl = ave_signals_freq_rl(:,idx)/norm(ave_signals_freq_rl(2:end,idx));
    %s_rl = stem(evals(2:end), abs(freq_mag_rl(2:end)), 'MarkerEdgeColor','green', 'color', 'k');
    xlabel('freq');
    if idx==1
        ylabel('|freq signal|');
    end
    title(sprintf('RL1: Window %d', idx));
end
%}