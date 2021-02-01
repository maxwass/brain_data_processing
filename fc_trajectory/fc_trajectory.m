%% setup
clear;
clc;

% Create colormap that is green for negative, red for positive,
% and a chunk inthe middle that is black.
%https://www.mathworks.com/matlabcentral/answers/81352-colormap-with-both-positive-and-negative-values
greenColorMap = [zeros(1, 132), linspace(0, 1, 124)];
redColorMap = [linspace(1, 0, 124), zeros(1, 132)];
colorMap = [redColorMap; greenColorMap; zeros(1, 256)]';


path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo)); %recursively adds all repo files/folders

% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';


path2repo = '~/Documents/MATLAB/brain_data_preprocess'; %CHANGE THIS
addpath(genpath(path2repo)); %recursively adds all repo files/folders

% The directory brain_data/ is assumed to have a folder for each patient. The
%  name of each directory should be exactly equal to the patient id. All
%  files (atlas, fmri scans, etc) should reside in this top level folder 
%  (e.g. dont make subfolder inside the subject's folder)
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

%Where to save matlab outputs. We save one .mat file per patient.
plots_dir   = '~/Documents/MATLAB/brain_data_preprocess/plots/fc_trajs';
subcortical_first = true;

%load subject list (list of strings of subject ids)
%load('data/hcp_1200_subject_list.mat')

%Where to save matlab outputs. We save one .mat file per patient.
plots_dir   = '~/Documents/MATLAB/brain_data_preprocess/plots/fc_trajs';
subcortical_first = true;

%load subject list (list of strings of subject ids)
%load('data/hcp_1200_subject_list.mat')
%subject_list = hcp1200_subject_list; %loaded from hcp_1200_subject_list.mat
brain_dataset = load('data/brain_dataset_sc_fc_pairs.mat');
subject_list = int2str(brain_dataset.final_subject_list); %all subjects with scs


atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';


fc_traj_params = struct('colorMap', colorMap, 'windowsize', 30, 'movesize', 10);


[total_time, num_iters, num_iters_left] = deal(0.0, 0, 1114);
%%
for i_index=1:length(subject_list)
    
    subject = subject_list(i_index,:); %must be char array for [...] to work later
       
    if(strcmp(atlas,"desikan"))                 %missing 4 and 39 are corpus collusum. WHAT IS 0??
        ChosenROI_cortical    =   setxor(0:70,[0, 4, 39]);
        ChosenROI_subcortical =   setxor(1:21,[1, 2]); %1,2 are general Left/Right hemisphere
        plot_write_path = plots_dir + "/desikan_cortical_" + string(subject);
    elseif(strcmp(atlas,"destrieux"))           %missing 42 and 117 are corpus collusum
        ChosenROI_cortical    =   setxor(0:150,[0, 42, 117]); %may have to reorder if you want all left hemisphere regions together for viz
        ChosenROI_subcortical = []; %setxor(1:21,[1, 2]) % same as desikan
        filename = plots_dir + "/destrieux_cortical_" + string(subject);
    else
        error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
    end
    
    chosen_roi = struct('cortical', ChosenROI_cortical, 'subcortical', ChosenROI_subcortical);

    
    path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    %check which fmri scans exist - can have neither, one, or both
    has_LR1_scan = isfile(path_to_LR1);
    has_RL1_scan = isfile(path_to_RL1);
    
    disp([num2str(i_index) ': patient ' subject ' - has LR? ' num2str(has_LR1_scan) ' | has RL? ' num2str(has_RL1_scan)]) 
    
    %if we've already created plot(s), skip
    lr_filename_fc  = plot_write_path + "_LR1.jpg";
    missing_LR1_fc_plot = ~isfile(lr_filename_fc);
    
    lr_filename_freq = plot_write_path + "_LR1_freq.jpg";
    missing_LR1_freq_plot = ~isfile(lr_filename_freq);

    rl_filename_fc = plot_write_path + "_RL1.jpg";
    missing_RL1_fc_plot = ~isfile(rl_filename_fc);
    
    rl_filename_freq = plot_write_path + "_RL1_freq.jpg";
    missing_RL1_freq_plot = ~isfile(rl_filename_freq);

  
    %% LR: load, create windowed fcs and/or freqs, save image
    if missing_LR1_fc_plot || missing_LR1_freq_plot
        name = 'LR';
        path2fmri = path_to_LR1;
        %load fmri data
        start = tic;
        dtseries_lr = load_fmri(atlas, path2fmri, subject, tasktype, raw_hcp_datafolder, name, chosen_roi.cortical, chosen_roi.subcortical);
        %created windowed fcs and ave signals
        [~, fc_corr_w, ave_signals_w] = windowed_fcs(dtseries_lr, fc_traj_params.windowsize, fc_traj_params.movesize);
        size_fc = size(fc_corr_w);
        time = toc(start);
        disp(['   time for ' name ' load + windowing: ' num2str(time)]);
        
        if missing_LR1_fc_plot
            fprintfp('   creating windowed fc plot...');
            start = tic;
            plot_cov_traj(fc_corr_w, fc_traj_params, lr_filename_fc, subject, name);
            stop = toc(start);
            fprintf('%s s\n', num2str(stop));
        end
        if missing_LR1_freq_plot
            fprintf('   creating freq plot...');
            start = tic;
            A = brain_dataset.transform_scs(:,:,i_index);
            plot_freq_traj(A, ave_signals_w, fc_traj_params, lr_filename_freq, subject, name);
            stop = toc(start);
            fprintf('%s s\n', num2str(stop));
        end
        
        clear dtseries_lr
    end
    
    %% RL: load, create windowed fcs, save image
    if missing_RL1_fc_plot
        name = 'RL';
        plot_filename_cov = rl_filename_fc;
        path = path_to_RL1;
        plot_cov_traj(fc_traj_params, plot_filename_cov, atlas, path, subject, tasktype, raw_hcp_datafolder, name, chosen_roi)
    end
    
    clear dtseries
end


function plot_cov_traj(windowed_fcs, fc_traj_params,  plot_filename, subject, name) 

	%% plot grid of covariances
    rows = 10; %CHANGE TO FLOW in tiledleyout so dont need to set this
    cols = ceil((1190/fc_traj_data.movesize)/rows); 
	txt = sprintf('patient %s %s1. Windowsize %d, Movesize %d', subject, name, fc_traj_params.windowsize, fc_traj_data.movesize);
        
    %ensures images are full sized
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    t = tiledlayout(rows,cols,'TileSpacing','compact','Padding','none');
    title(t,txt, 'FontSize', 15)
        
	start_plot = tic;
	num_covs = size_fc(3);
    for c = 1:num_covs
        nexttile;
        corr = windowed_fcs(:,:,c);
        corr_cortical = corr(20:end, 20:end);
        imagesc(corr_cortical);
        set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
        if c<cols
        	title(sprintf('%d: lr', c), 'FontSize', 8)
        end
        axis image
    end
	% Apply the colormap.
	colormap(fc_traj_data.colorMap);
	colorbar
	caxis([-1,1])
   
	exportgraphics(t,plot_filename,'BackgroundColor','white','Resolution',600)
	time_plot = toc(start_plot);
	disp(['   time for ' name ' plotting: ' num2str(time_plot)]);
end

function plot_freq_traj(A, ave_signals, fc_traj_params,  plot_filename, subject, name) 

    D = diag(sum(A,2));
    D_= D.^(-1/2);
    L = D-A;
    L_norm = D_*L*D;
    [evals, evecs] = eig(L_norm);
    GFT = transpose(evecs);
    
    
    

	%% plot grid of freq signals
    rows = 10; %CHANGE TO FLOW in tiledleyout so dont need to set this
    cols = ceil((1190/fc_traj_data.movesize)/rows); 
	txt = sprintf('patient %s %s1. Windowsize %d, Movesize %d', subject, name, fc_traj_params.windowsize, fc_traj_data.movesize);
        
    %ensures images are full sized
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    t = tiledlayout(rows,cols,'TileSpacing','compact','Padding','none');
    title(t,txt, 'FontSize', 15)
        
	start_plot = tic;
	num_covs = size_fc(3);
    for c = 1:num_covs
        nexttile;
        corr = ave_signals(:,:,c);
        corr_cortical = corr(20:end, 20:end);
        imagesc(corr_cortical);
        set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
        if c<cols
        	title(sprintf('%d: lr', c), 'FontSize', 8)
        end
        axis image
    end
	% Apply the colormap.
	colormap(fc_traj_data.colorMap);
	colorbar
	caxis([-1,1])
   
	exportgraphics(t,plot_filename,'BackgroundColor','white','Resolution',600)
	time_plot = toc(start_plot);
	disp(['   time for ' name ' plotting: ' num2str(time_plot)]);

end


