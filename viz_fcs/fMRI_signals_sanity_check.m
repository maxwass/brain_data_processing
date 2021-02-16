%% sanity checks on fmri time series data
close all;
%% load raw fmri files, process with desikan atlas. Save to mat file for 
% fast access (don't need to read fmri file now)

clear; clc;

%% directory to read fmri from, and write processed data to
% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

% It is a directory of processed fmri signals according to the desikan atlas. 
%One .mat file per scan (one for LR, one for RL, if they exist).
%This directory is NOT included in repo: ~ 1GB
cached_desikan   = '~/Documents/MATLAB/brain_data_preprocess/cached_desikan';
cached_destrieux = '~/Documents/MATLAB/brain_data_preprocess/cached_destrieux';


%% determine which patients to do this for
% new FCs (downloaded from HCP_1200 server and did local computation)
% hcp1200_subject_list (1x1113 double): 
%load('data/hcp_1200_subject_list.mat');
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

%% values to compute
%what is the average value in *each* brain region (roi) in each scan
%  average across rows
%mean_rois = zeros(num_rois, num_subjects);

%what is the average value of each *vector observation* in each scan
%  average across columns
%num_obsvs = 1200;
%mean_obsvs = zeros(1, num_subjects * num_obsvs);




%% parameters to change for different tasks, atlas', etc
atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';
sub_tasktype = 'REST1';
include_subcortical = true;
subcortical_first = true;
GSO = 'L_norm';


if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_desikan;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;
	cached_data_folder = cached_destrieux;
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end


%% perform loading and saving
% load each fmri file, save processed dtseries to .mat file, then update
%  completed array for time estimate
% save full cov (87x87) and remove subcortical later

num_subjects  = length(subject_list);

counter = 0;
for i_index = 1:length(subject_list)
    subject = char(string(subject_list(i_index)));
    
    %subject = subject_list(i_index,:); %must be char array for [...] to work later
       
	cached_filename_lr = [cached_data_folder,'/', tasktype, '/', subject,'_LR.mat'];
	cached_filename_rl = [cached_data_folder,'/', tasktype, '/', subject,'_RL.mat'];
    
    
    has_cached_lr = isfile(cached_filename_lr);
    has_cached_rl = isfile(cached_filename_rl);
    
    % attempt load lr
    if has_cached_lr
        dtseries_lr = load(cached_filename_lr).dtseries;
    end
    
    % attempt load rl
    if has_cached_rl
        dtseries_rl = load(cached_filename_rl).dtseries;
    end
    
    
    %plotting
    [N, num_obsvs] = size(dtseries_lr);
    mean_roi_lr   = mean(dtseries_lr, 2); %across rows
    stdv_roi_lr   = std(dtseries_lr, 0, 2);
    mean_obsvs_lr = mean(dtseries_lr, 1); %down cols
    stdv_obsvs_lr = std(dtseries_lr, 0, 1);
    
    [N, num_obsvs] = size(dtseries_rl);
    mean_roi_rl    = mean(dtseries_rl, 2); %across rows
    stdv_roi_rl    = std(dtseries_rl, 0, 2);
    mean_obsvs_rl  = mean(dtseries_rl, 1); %down cols
    stdv_obsvs_rl  = std(dtseries_rl, 0, 1);
    
    
    
    if has_cached_lr && has_cached_rl
        scan_dir = 'LR';
        f = figure();
        t = tiledlayout(f, 4,1);
        title_txt = sprintf('Summaries of Node values in fmri scans\nSubject %s | Task: %s | left LR, right RL', subject, tasktype);
        title(t, title_txt, 'FontSize', 25);
        
        ax_energies = nexttile(t);
        plot_energies(ax_energies, dtseries_lr, dtseries_rl);
        
        ax_obsvs = nexttile(t);
        plot_ave_val(ax_obsvs, mean_obsvs_lr, mean_obsvs_rl);

        ax_roi = nexttile(t);
        plot_mean_vec(ax_roi,mean_roi_lr, stdv_roi_lr, mean_roi_rl, stdv_roi_rl);

        ax_freq = nexttile(t);
        plot_mean_vec_freq(ax_freq, mean_roi_lr, mean_roi_rl, subject, atlas, include_subcortical, GSO);
    end

end

function plot_ave_val(ax, data_lr, data_rl)
        histogram(ax, data_lr, 'DisplayName', 'lr');
        hold on;
        histogram(ax,data_rl, 'DisplayName', 'rl');
        hold off;
        xlabel(ax, 'ave value over *all nodes* for *each* observation');
        ylabel(ax, 'frequency');
        legend;
end

function plot_mean_vec(ax, mean_data_lr, stdv_data_lr, mean_data_rl, stdv_data_rl);
    transp = .2;
    y = [mean_data_lr, mean_data_rl];
    
    x = 1:length(mean_data_lr);
    b_lr = bar(ax, y);
    %b_lr.FaceAlpha = transp;
	xlabel(ax, 'ROI');
	ylabel(ax, 'value');
	title('Node values in mean vector over *all* observations', 'FontSize', 18);
    y_lr = yline(ax, mean(mean_data_lr),'-',sprintf('lr ave = %.1f',mean(mean_data_lr)), 'LineWidth',1, 'Color', 'r', 'FontSize', 12);
    y_lr.LabelHorizontalAlignment = 'left';
    y_rl = yline(ax, mean(mean_data_rl),'-',sprintf('lr ave = %.1f',mean(mean_data_rl)), 'LineWidth',1, 'Color', 'b', 'FontSize', 12);
    %hold on;
    %b_rl = bar(ax, mean_data_rl);
    %b_rl.FaceAlpha = transp;
    
    %hold off;
	
    %hold on;
	%er = errorbar(ax, x, mean_data, -stdv_data ,stdv_data); 
	%er.Color = [0 0 0];                            
	%er.LineStyle = 'none'; 
	%yline(ax, mean(mean_data),'-',sprintf('average value = %.1f',mean(mean_data)), 'LineWidth',3, 'Color', 'r', 'FontSize', 15);
	%hold off;
end

function plot_mean_vec_freq(ax, mean_vec_lr, mean_vec_rl, subject, atlas, include_subcortical, GSO)
    
    [GFT, evals] = extract_GFT(subject, atlas, include_subcortical, GSO);
    mean_vec_freq_lr = GFT*mean_vec_lr;
    mean_vec_freq_rl = GFT*mean_vec_rl;
    plot_eigs = diag(evals);
	freq_plot_lr = stem(ax, plot_eigs, mean_vec_freq_lr, 'MarkerEdgeColor','green', 'color', 'k', 'DisplayName', 'lr');
	hold on;
    freq_plot_rl = stem(ax, plot_eigs, mean_vec_freq_rl, 'MarkerEdgeColor','red', 'color', 'k', 'DisplayName', 'rl');
    hold off;
    legend;
    
    xlabel(ax, sprintf('graph frequency of GSO %s', GSO));
	%ylabel(ax_freq, 'freq signal');
	txt = sprintf('Freq repr of *mean* vector over *all observations*');
	if isequal(GSO,'L')
        N = length(mean_vec_lr);
        m_lr = mean_vec_freq_lr(1)*sqrt(N)/N;
        m_rl = mean_vec_freq_rl(1)*sqrt(N)/N;
        txt = sprintf('%s: 0-freq component lr/rl %10.1f./%10.1f. x_f(0)*sqrt(N)/N =%10.1f/%10.1f.', txt, mean_vec_freq_lr(1),mean_vec_freq_rl(1), m_lr, m_rl);
    end
	title(ax, txt, 'FontSize', 18);
end

function plot_energies(ax, time_series_1, time_series_2)
    energies_1 = vecnorm(time_series_1).^2;
    energies_2 = vecnorm(time_series_2).^2;
    
    histogram(ax, time_series_1, 'DisplayName', 'lr');
    hold on;
    histogram(ax, time_series_2, 'DisplayName', 'rl');
    hold off;

	xlabel(ax, 'Energy');
	ylabel(ax, 'frequency');
	legend;


end
