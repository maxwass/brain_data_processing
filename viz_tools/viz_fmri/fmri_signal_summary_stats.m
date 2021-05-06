%% Find summary statistics of fmri data
% For all patient with any fc scans (even without sc scan) in REST1,
% view mean distribution
% -view stdv distribution
% -view median distribution
% -view min/max distribution


close all;
clear; clc;
%% load raw fmri files, process with desikan atlas. Save to mat file for 
% fast access (don't need to read fmri file now)


%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
sub_tasktype = 'REST1';
include_subcortical = true;
subcortical_first = true;
GSO = 'L';

summary_stats_file = 'fmri_desikan_summary_stats.mat';

    
if isfile(summary_stats_file)
    load(summary_stats_file);
else


%% directory to read fmri from, and write processed data to
% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

% It is a directory of processed fmri signals according to the desikan atlas. 
%One .mat file per scan (one for LR, one for RL, if they exist).
%This directory is NOT included in repo: ~ 1GB
cached_desikan   = '~/Documents/MATLAB/brain_data_preprocess/data/cached_desikan';
cached_destrieux = '~/Documents/MATLAB/brain_data_preprocess/data/cached_destrieux';


%% determine which patients to do this for
% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
subject_list_char = load('data_accounting/hcp_1200_subject_list.mat').hcp1200_subject_list; % 1113x1 char array
subject_list = zeros(size(subject_list_char));
for i = 1:length(subject_list_char)
    subject_list(i) = str2double(subject_list_char(i,:));
end

%fc_sc_set_file = load('fc_and_sc_sets.mat');
%subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)



if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu');
	cached_data_folder = cached_desikan;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu');
	cached_data_folder = cached_destrieux;
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end


%% perform loading and saving
% load each fmri file, save processed dtseries to .mat file, then update
%  completed array for time estimate
% save full cov (87x87) and remove subcortical later

num_subjects  = length(subject_list);

filenames = {};
filenames_lr = {};
filenames_rl = {};
for i_index = 1:length(subject_list)
    subject = char(string(subject_list(i_index)));
    cached_filename_lr = [cached_data_folder,'/', tasktype, '/', subject,'_LR.mat'];
	cached_filename_rl = [cached_data_folder,'/', tasktype, '/', subject,'_RL.mat'];

    has_cached_lr = isfile(cached_filename_lr);
    has_cached_rl = isfile(cached_filename_rl);
    
    % attempt load lr
    if has_cached_lr
        filenames_lr{end+1} = cached_filename_lr;
        filenames{end+1}    = cached_filename_lr;
    end
    
    % attempt load rl
    if has_cached_rl
        filenames_rl{end+1} = cached_filename_rl;
        filenames{end+1}    = cached_filename_rl;
    end
end

lrs = compute_summary_stats(filenames_lr, include_subcortical);
rls = compute_summary_stats(filenames_rl, include_subcortical);
save(summary_stats_file, 'lrs', 'rls', 'include_subcortical', 'filenames_lr', 'filenames_rl');

end




fs = 8;
f = figure();
t = tiledlayout(6,1);
nbins = 30;
txt = sprintf('# lrs: %d, # rls: %d', length(filenames_lr), length(filenames_rl));
title(t,txt);

ax = nexttile();
plot_hist(ax, lrs.mean_distrib, lrs.mean_distrib, 'mean distribution', 'energy', fs, nbins)

ax = nexttile();
plot_hist(ax, lrs.median_distrib, rls.median_distrib, 'median distribution', 'energy', fs, nbins)

ax = nexttile();
plot_hist(ax, lrs.stdv_distrib, rls.stdv_distrib, 'stdv distribution', 'stdv', fs, nbins)

ax = nexttile();
plot_hist(ax, lrs.max_distrib, rls.max_distrib, 'max distribution', 'energy', fs, nbins)

ax = nexttile();
plot_hist(ax, lrs.min_distrib, rls.min_distrib, 'min distribution', 'energy', fs, nbins)

ax = nexttile();
plot_hist(ax, lrs.range_distrib, rls.range_distrib, 'range distribution', 'energy', fs, nbins)
    



fs = 8;
f1 = figure();
t1 = tiledlayout(4,1);
txt = sprintf('# lrs: %d, # rls: %d', length(filenames_lr), length(filenames_rl));
title(t1,txt);

ax = nexttile();
y = [mean(lrs.mean_vectors,2), mean(rls.mean_vectors,2)];
bar(ax, y);
xlabel(ax, 'ROI');
ylabel(ax, 'value');
title('Node values in Mean of MEAN vector over (LR or RL) observations', 'FontSize', 18);

ax = nexttile();
y = [mean(lrs.median_vectors,2), mean(rls.median_vectors,2)];
bar(ax, y);
xlabel(ax, 'ROI');
ylabel(ax, 'value');
title('Node values in Mean vector of MEDIAN vector over (LR or RL) observations', 'FontSize', 18);


ax = nexttile();
means = [lrs.mean_vectors, rls.mean_vectors];
y = mean(means,2);
bar(ax, y);
hold on;
overall_ave_node_val = mean(y);
yline(ax, overall_ave_node_val, 'Label', sprintf('ave node val %f', overall_ave_node_val));
xlabel(ax, 'ROI');
ylabel(ax, 'value');
hold off;
title('Node values in Mean vector of MEAN vector over *all* observations', 'FontSize', 18);


ax = nexttile();
means = [lrs.median_vectors, rls.median_vectors];
y = mean(means,2);
bar(ax, y);
hold on;
overall_ave_node_val = mean(y);
yline(ax, overall_ave_node_val, 'Label', sprintf('ave node val %.0f', overall_ave_node_val));
xlabel(ax, 'ROI');
ylabel(ax, 'value');
hold off;
title('Node values in Mean vector of MEDIAN vector over *all* observations', 'FontSize', 18);



function s = compute_summary_stats(filenames, include_subcortical)
    
    dtseries = load(filenames{1}).dtseries;
    if ~include_subcortical
        dtseries = dtseries(20:end,:);
    end
    [num_roi, num_obsvs] = size(dtseries);

    nf = length(filenames);
    s.mean_distrib   = zeros(nf,1);
    s.stdv_distrib   = zeros(nf,1);
    s.max_distrib    = zeros(nf,1);
    s.min_distrib    = zeros(nf,1);
    s.median_distrib = zeros(nf,1);
    s.range_distrib  = zeros(nf,1);
    s.mean_vectors   = zeros(num_roi, nf);
    s.median_vectors = zeros(num_roi, nf);

    for f = 1:length(filenames)
        filename = filenames{f};
        dtseries = load(filename).dtseries;
        if ~include_subcortical
            dtseries = dtseries(20:end,:);
        end
    
        energies          = vecnorm(dtseries,2).^2;
        s.mean_distrib(f)   = mean(energies);
        s.stdv_distrib(f)   = std(energies);
        s.max_distrib(f)    = max(energies);
        s.min_distrib(f)    = min(energies);
        s.median_distrib(f) = median(energies);
        s.range_distrib(f)  = max(energies)-min(energies);
        s.mean_vectors(:,f)     = mean(dtseries,2);
        s.median_vectors(:,f)   = median(dtseries,2);
    end
    
end

function plot_hist(ax, x_lr, x_rl, display_name, xlabel, fs, nbins)
    histogram(ax, x_lr, nbins, 'DisplayName', [display_name '- LR'], 'Normalization', 'probability','FaceAlpha',.4);
    hold on;
    m_lr = mean(x_lr);
    s_lr = std(x_lr);
    xline(ax, m_lr,'-','ave', 'LineWidth',1, 'Color', 'r', 'FontSize', fs, 'DisplayName', sprintf('mean = %.0f, stdv = %.0f',m_lr,s_lr));
    
    histogram(ax, x_rl, nbins, 'DisplayName', [display_name '- RL'], 'Normalization', 'probability','FaceAlpha',.4);
    m_rl = mean(x_rl);
    s_rl = std(x_rl);
    xline(ax, m_rl,'-','ave', 'LineWidth',1, 'Color', 'b', 'FontSize', fs, 'DisplayName', sprintf('mean = %.0f, stdv = %.0f',m_rl,s_rl));
    
    hold off;
    %xlabel(ax, [xlabel]);
    ylabel(ax, 'frequency');
    legend;

end
