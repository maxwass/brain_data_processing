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

%Where to save matlab outputs. We save one .mat file per patient.
plots_dir   = '~/Documents/MATLAB/brain_data_preprocess/plots/fc_trajs';
subcortical_first = true;

%load subject list (list of strings of subject ids)
load('data/hcp_1200_subject_list.mat')
subject_list = hcp1200_subject_list; %loaded from hcp_1200_subject_list.mat

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1'; %'tfMRI_GAMBLING'; 'tfMRI_MOTOR'; 'tfMRI_GAMBLING'; 'tfMRI_SOCIAL'; 'tfMRI_LANGUAGE';


cov_traj_data = struct('colorMap', colorMap, 'windowsize', 30, 'movesize', 10);


[total_time, num_iters, num_iters_left] = deal(0.0, 0, 1114);

%% computation
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
    has_LR1 = isfile(path_to_LR1);
    has_RL1 = isfile(path_to_RL1);
    
    disp([num2str(i_index) ': patient ' subject ' - has LR? ' num2str(has_LR1) ' | has RL? ' num2str(has_RL1)]) 
    
    %if we've already created LR image, skip it
    lr_filename = plot_write_path + "_LR1.jpg";
    if isfile(lr_filename)
        has_LR1 = false;
    end
    
    %if we've already created RL image, skip it
    rl_filename = plot_write_path + "_RL1.jpg";
    if isfile(rl_filename)
        has_RL1 = false;
    end
     
    
    %% LR: load, create windowed fcs, save image
    if has_LR1
        name = 'LR';
        plot_filename = lr_filename;
        path = path_to_LR1;
        plot_cov_traj(cov_traj_data, plot_filename, atlas, path, subject, tasktype, raw_hcp_datafolder, name, chosen_roi)
    end
    
    %% RL: load, create windowed fcs, save image
    if has_RL1
        name = 'RL';
        plot_filename = rl_filename;
        path = path_to_RL1;
        plot_cov_traj(cov_traj_data, plot_filename, atlas, path, subject, tasktype, raw_hcp_datafolder, name, chosen_roi)
    end
    
end


function plot_cov_traj(cov_traj_data,  plot_filename, atlas, path2fmri, subject, tasktype, raw_hcp_datafolder, name, chosen_roi)
%% attempt: load, create windowed fcs
	start = tic;
	dtseries = load_fmri(atlas, path2fmri, subject, tasktype, raw_hcp_datafolder, name, chosen_roi.cortical, chosen_roi.subcortical);
	[fc_cov, fc_corr] = windowed_fcs(dtseries, cov_traj_data.windowsize, cov_traj_data.movesize);
	size_fc = size(fc_cov);
	time = toc(start);
	disp(['   time for ' name ' loading: ' num2str(time)]);   

	%% plot grid of covariances
    rows = 10; %CHANGE TO FLOW in tiledleyout so dont need to set this
    cols = ceil((1190/cov_traj_data.movesize)/rows); 
	txt = sprintf('patient %s %s1. Windowsize %d, Movesize %d', subject, name, cov_traj_data.windowsize, cov_traj_data.movesize);
        
    %ensures images are full sized
    figure('units','normalized','outerposition',[0 0 1 1], 'visible', 'off');
    t = tiledlayout(rows,cols,'TileSpacing','compact','Padding','none');
    title(t,txt, 'FontSize', 15)
        
	start_plot = tic;
	num_covs = size_fc(3);
    for c = 1:num_covs
        nexttile;
        corr = fc_corr(:,:,c);
        corr_cortical = corr(20:end, 20:end);
        imagesc(corr_cortical);
        set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
        if c<cols
        	title(sprintf('%d: lr', c), 'FontSize', 8)
        end
        axis image
    end
	% Apply the colormap.
	colormap(cov_traj_data.colorMap);
	colorbar
	caxis([-1,1])
   
	exportgraphics(t,plot_filename,'BackgroundColor','white','Resolution',600)
	time_plot = toc(start_plot);
	disp(['   time for ' name ' plotting: ' num2str(time_plot)]);
    
    clear dtseries fc_cov fc_corr t

end

