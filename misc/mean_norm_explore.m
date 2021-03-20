%% Effect of Mean Normalization/Removing DC component/etc on cov

% https://stats.stackexchange.com/questions/391838/does-mean-centering-reduce-covariance
% any 'global' scalar shift of data shouldnt change covariance. Cov(x-a,y-b) = Cov(x,y)  
% Removing DC component IS NOT a scalar shift. When you remove the DC
% component of *each* signal, this is like scalar shifting each signal by a different amount. 
%  Thus not same as 'global' scalar shift'. 


%load in some data
clear; clc; close all;

%% Define high level hp's
raw_hcp_datafolder = '/Volumes/Elements/brain_data';

sum_stat_file = load('fmri_desikan_summary_stats.mat');
colormap_file = load("viz_fcs/correlations_colormap.mat");

atlas = "desikan"; %"destrieux"
tasktype='rfMRI_REST1';
include_subcortical = false;
GSO = "L";

if(strcmp(atlas,"desikan"))
    chosen_roi         = load('data/desikan_roi_zhengwu', 'roi').roi;
elseif(strcmp(atlas,"destrieux"))
    chosen_roi         = load('data/destrieux_roi_zhengwu', 'roi').roi;h
else
	error("Atlas " + atlas + " not found. Use Desikan or Destrieux.")
end

if include_subcortical
    num_rois = length(chosen_roi.cortical) + length(chosen_roi.subcortical);
else
    num_rois = length(chosen_roi.cortical);
end


%% Define inputs needed for each step in pipeline

%Frequency filtering inputs
%which to *keep*
freq_filtering_intervals = {[2,num_rois]};

%% determine which patients to do this for
fc_sc_set_file = load('fc_and_sc_sets.mat');
subject_list = fc_sc_set_file.exist_any_fc_and_sc; % (1064x1 int64)

for subject_idx = 1:length(subject_list)
    subject     = char(num2str(subject_list(subject_idx,:)));
    path_to_LR1 = [raw_hcp_datafolder '/' subject '/' tasktype '_LR_Atlas_hp2000_clean.dtseries.nii'];
    path_to_RL1 = [raw_hcp_datafolder '/' subject '/' tasktype '_RL_Atlas_hp2000_clean.dtseries.nii'];
    
    path2fmris = {path_to_LR1, path_to_RL1};
    [GFT, evals_vec] = extract_GFT(subject, atlas, include_subcortical, GSO);
    iGFT = GFT';
    
    for path2fmris_idx = 1:length(path2fmris)
        path2fmri = path2fmris{path2fmris_idx};
        
    
        %% load and center 'raw' fmri data. 
        x        = load_functional_dtseries(atlas, path2fmri, subject, raw_hcp_datafolder, chosen_roi);
        %x_center = x - mean(x,2);

        %remove subcortical if needed
        [num_rois,~] = size(x);
        if include_subcortical==0
            roi_idxs = (20:num_rois); %subcortical are first. See README
        else
            roi_idxs = (1:num_rois);
        end
        
        average_vector_ = mean([sum_stat_file.lrs.mean_vectors, sum_stat_file.rls.mean_vectors], 2);
        average_vector  = average_vector_(roi_idxs,:);
        average_node_value = mean(average_vector);
        
        
        
 
        f = figure();
        t = tiledlayout('flow');
        axes_list = gobjects(6,1);
        
        axes_list(1) = nexttile();
        raw_fc = corr(x(roi_idxs,:)'); 
        imagesc(raw_fc); %raw cov
        title('raw corr');
        
        axes_list(2) = nexttile();
        all_scalar_center = x(roi_idxs,:) - average_node_value;
        global_scalar_center_fc = corr(all_scalar_center');
        imagesc(global_scalar_center_fc);
        title('AA GLOBAL Scalar : subtract ave *node* val over ALL mean vectors for ALL subject');
        
        axes_list(3) = nexttile();
        local_scalar_center = x(roi_idxs,:) - mean(mean( x(roi_idxs,:),2) );
        local_scalar_center_fc = corr(local_scalar_center');
        imagesc(local_scalar_center_fc);
        title('AA LOCAL Scalar: subtract ave *node* val over mean vector of LOCAL subject');
        
        axes_list(4) = nexttile();
        all_vector_center = x(roi_idxs,:) - average_vector;
        global_vector_center_fc = corr(all_vector_center');
        imagesc(global_vector_center_fc);
        title('BB GLOBAL Vector: subtract ave mean vector *vector* val over ALL subject');
        
        axes_list(5) = nexttile();
        local_vector_center = x(roi_idxs, :) - mean(x(roi_idxs, :),2);
        local_vector_center_fc = corr(local_vector_center');
        imagesc( local_vector_center_fc );
        title('BB LOCAL Vector: subtract ave mean vector *vector* val over LOCAL subject');
        
        axes_list(6) = nexttile();
        x_cort = x(roi_idxs,:);
        x_cort_freq = GFT*x_cort;
        x_cort_freq(1,:) = 0; %remove DC component
        x_cort_filt = iGFT*x_cort_freq;
        dc_remove_fc = corr(x_cort_filt');
        imagesc( dc_remove_fc );
        title('Remove DC component - SHOULD EQUAL SCALAR TRANSFORM');
            
        set(axes_list, 'ColorMap', colormap_file.colorMap);
        set(axes_list,'CLim',[-1 1]);
        cb = colorbar(axes_list(end),'WestOutside');
        linkaxes(axes_list,'xy');
        
        %scalar shifts should NOT affect cov (and thus corr) matrix
        local_scalar_center_fc_same  = all(all( abs(raw_fc-local_scalar_center_fc) <.001));
        global_scalar_center_fc_same = all(all( abs(raw_fc-global_scalar_center_fc)<.001));
        fc_remove_fc_same = all(all(abs(raw_fc - dc_remove_fc)<.001));
        
        
        
        


        
    end
end


