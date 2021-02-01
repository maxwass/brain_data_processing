%freq plots
path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo));

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

brain_dataset = load('data/brain_dataset_sc_fc_pairs.mat');
subject_list = int2str(brain_dataset.final_subject_list); %all subjects with scs

subject = subject_list(1,:);
path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
name = 'LR';
path2fmri = path_to_LR1;
if exist('dtseries_lr','var') ==0
    dtseries_lr = load_fmri(atlas, path2fmri, subject, tasktype, raw_hcp_datafolder, name, chosen_roi.cortical, chosen_roi.subcortical);
    [~, fc_corr_lr, ave_signals_lr] = windowed_fcs(dtseries_lr, fc_traj_params.windowsize, fc_traj_params.movesize);
end
name = 'RL';
path2fmri = path_to_RL1;
if exist('dtseries_rl','var') ==0
    dtseries_rl = load_fmri(atlas, path2fmri, subject, tasktype, raw_hcp_datafolder, name, chosen_roi.cortical, chosen_roi.subcortical);
    [~, fc_corr_rl, ave_signals_rl] = windowed_fcs(dtseries_rl, fc_traj_params.windowsize, fc_traj_params.movesize);
end


include_subcortical = false;
A_full = brain_dataset.transform_scs(:,:,i_index);
A_cortical = A_full(20:end,20:end);
if include_subcortical ==1
    A = A_full;
else
    A = A_cortical;
end


%% GFT computation
D_vec  = sum(A,2);
D = diag(D_vec);
D_norm = diag(D_vec.^(-.5));
L = D-A;
L_norm = D_norm*L*D_norm;
[evecs, evals] = eig(L);
GFT = transpose(evecs);

%% transform into GF space with GFT

if include_subcortical==1
    ave_signals_freq_lr = GFT*ave_signals_lr; %cols are signals
    ave_signals_freq_rl = GFT*ave_signals_rl; %cols are signals
else
    ave_signals_freq_lr = GFT*ave_signals_lr(20:end,:); %cols are signals
    ave_signals_freq_rl = GFT*ave_signals_rl(20:end,:); %cols are signals
end

%% Which signal windows are the most similar? EXCLUDE ZERO FREQ COMPONENT
D   = pdist(ave_signals_freq_lr(2:end,:)', 'euclidean');
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
lr_norms = sqrt(sum(ave_signals_freq_lr.^2, 1)); %should i disregard 0 freq component?
rl_norms = sqrt(sum(ave_signals_freq_rl.^2, 1));

%divide each column by respective norm
ave_signals_freq_lr_norm = bsxfun(@rdivide, ave_signals_freq_lr, lr_norms);
ave_signals_freq_rl_norm = bsxfun(@rdivide, ave_signals_freq_rl, rl_norms);

%find max for fixed y axis limits (excluding 0 component)
sig_concat = abs([ave_signals_freq_lr_norm, ave_signals_freq_rl_norm]);
y_max = max(max(sig_concat(2:end,:)));

evals = diag(evals);

num_windows = 10;
t = tiledlayout(2,num_windows);

for idx = 1:num_windows
    nexttile(t)
    f_lr = abs(ave_signals_freq_lr_norm(:,idx));
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