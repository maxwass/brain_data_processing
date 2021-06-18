%% Define inputs needed for each step in pipeline
clear; close;

%% patient and scan info
atlas = 'desikan';
task = 'REST';
include_subcortical = false;
rand_seed = 10;
rng(rand_seed);

%% Subset construction inputs
subset_construction = struct("name", 'full'); %'sep');
%subset_construction = struct("name", 'windowing', 'windowsize', 400, 'movesize', 400);
%subset_construction = struct("name", 'sampling', "num_subsets", 3, "sps", 700, "with_replacement", false);

%% create unique filename based on filters and save
%filepath = unique_filename(subset_construction, "all_fc_sc");
filepath = 'data/created_datasets/all_fc_sc.mat';

%% Create actual dataset
[data] = create_dataset_full(atlas, include_subcortical);
save(filepath, "data", "atlas", "include_subcortical", "task", '-v7');


function filepath = unique_filename(subset_construction, filename)

%% combine info and filters for unique filename

if isequal(subset_construction.name, "windowing")
    subset_construction_txt = sprintf("windowsize%d_movesize%d", ...
        subset_construction.windowsize,...
        subset_construction.movesize);
elseif isequal(subset_construction.name, "sampling")
    subset_construction_txt = sprintf("numsubsets%d_sps%d_withreplacement%d", ...
        subset_construction.num_subsets,...
        subset_construction.sps,...
        subset_construction.with_replacement);
elseif isequal(subset_construction.name, "full")
    %using full to make fc
    subset_construction_txt = "full";
elseif isequal(subset_construction.name, "sep")
    subset_construction_txt = "sep";
end


%% each subset technique has own folder
dataset_folder = "subsample_datasets";
if isequal(subset_construction.name, "windowing")
    filepath = fullfile("data", dataset_folder, "windowing", filename);
elseif isequal(subset_construction.name, "sampling")
    filepath = fullfile("data", dataset_folder, "sampling", filename);
elseif isequal(subset_construction.name, "full")
    filepath = fullfile("data", dataset_folder, "full", filename);
elseif isequal(subset_construction.name, "sep")
    filepath = fullfile("data", dataset_folder, "sep", filename);
else
    error('Dont know where to save data to');
end

end

function [data] = create_dataset_full(atlas, include_subcortical)

%% which rois to consider
roi_idxs = get_roi_idxs(atlas, include_subcortical);
num_rois = length(roi_idxs);

%% determine which patients to do this for
subject_list = load('hcp_1200_subject_list.mat').hcp1200_subject_list;

%% loop through all scans
% preallocate if slow
which_scans_exist = struct('REST1_LR', false, 'REST2_LR', false, 'REST1_RL', false, 'REST2_RL', false);
fcs_individual = struct('REST1_LR', [], 'REST2_LR', [], 'REST1_RL', [], 'REST2_RL', []);
fcs_task_grouped = struct('REST1', [], 'REST2', []);
fcs_scandir_grouped = struct('LR', [], 'RL', []);
fcs_all = struct('all', []);
fcs = struct('individual', fcs_individual, 'task_grouped', fcs_task_grouped, 'scandir_grouped', fcs_scandir_grouped, 'fcs_all', fcs_all);
data = repmat(struct('subject_id',-1,'fcs', fcs, 'sc', zeros(num_rois,num_rois), 'which_scans_exist', which_scans_exist), length(subject_list),1);
subject_ids = zeros(3000,1,'uint32');

% num_scans :: number of subjects we include in dataset
[num_iters, num_scans, total_time] = deal(0, 0, 0);
for subject_idx = 1:length(subject_list)
    subject_int = str2double(subject_list(subject_idx,:));
    p = patient(subject_int);
    [all, at_least_one_of_each, at_least_one_per_task, none] = p.which_fcs_exist(atlas, include_subcortical);
    
    start = tic();
    if p.exist_sc(atlas) && ~none
        
        num_scans = num_scans + 1;
        fprintf('idx/subjects included/total %d/%d/%d: ', subject_idx, num_scans, length(subject_list));

    
        %% populate struct
        subject_ids(num_scans)     = int32(subject_int); %correct index?
        data(num_scans).subject_id = subject_int;
        data(num_scans).sc         = p.sc(atlas, include_subcortical);
        
        data(num_scans).fcs = populate_fcs(p, atlas, include_subcortical);
        
        %data(num_scans).fc         = p.full_rest_fc(atlas, include_subcortical, mean_norm);

        [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = p.rest_scans(atlas, include_subcortical);
        data(num_scans).which_scans_exist = struct('REST1_LR', lr_rest1.exist(), 'REST2_LR', lr_rest2.exist(), 'REST1_RL', rl_rest1.exist(), 'REST2_RL', rl_rest2.exist());
  
        %optional viz
        %p.viz_fcs(atlas, include_subcortical, 99.5, mean_norm);
    end
            
	elapsed_time = toc(start);
 
    %elapsed_time = stop;

    num_iters_left = length(subject_list) - subject_idx;
    total_time = total_time + elapsed_time;
    ave_patient_time = total_time/num_scans;
    expected_time_left = ave_patient_time*num_iters_left/60; %in hrs
    fprintf('elapsed time: %.2f | ave time: %.1f | expected time remain %.2f (mins)=== %%\n', elapsed_time, ave_patient_time, expected_time_left);
end


data = data(1:num_scans);

end


function [fcs_struct] = populate_fcs(p, atlas, include_subcortical)
    [rest1_lr_fc, rest1_rl_fc, rest2_lr_fc, rest2_rl_fc] = p.individual_fc(atlas, include_subcortical, 'cov');
    fcs_individual = struct('REST1_LR', rest1_lr_fc, 'REST2_LR', rest2_lr_fc, 'REST1_RL', rest1_rl_fc, 'REST2_RL', rest2_rl_fc);
    
    [rest1_fc, rest2_fc] = p.task_grouped_fc(atlas, include_subcortical, 'cov');
    fcs_task_grouped = struct('REST1', rest1_fc, 'REST2', rest2_fc);
    
    [lr_fc, rl_fc] = p.direction_grouped_fc(atlas, include_subcortical, 'cov');
    fcs_scandir_grouped = struct('LR', lr_fc, 'RL', rl_fc);
    
    [total_fc] = p.concat_all_fc(atlas, include_subcortical, 'cov');
    fcs_all = struct('all', total_fc);
    
    fcs_struct = struct('individual', fcs_individual, 'task_grouped', fcs_task_grouped, 'scandir_grouped', fcs_scandir_grouped, 'fcs_all', fcs_all);
end



function [data] = create_dataset_sep(atlas, include_subcortical)

%% which rois to consider
roi_idxs = get_roi_idxs(atlas, include_subcortical);
num_rois = length(roi_idxs);

%% determine which patients to do this for
subject_list = load('hcp_1200_subject_list.mat').hcp1200_subject_list;

%% loop through all scans
% preallocate if slow
which_fcs = struct('lr1', false, 'lr2', false, 'rl1', false, 'rl2', false);
data = repmat(struct('subject_id',-1,'fc_REST1',[], 'fc_REST2', [], 'sc', zeros(num_rois,num_rois), 'which_fcs', which_fcs),length(subject_list),1);
subject_ids = zeros(3000,1,'uint32');

[num_iters, num_scans, total_time] = deal(0, 0, 0);
for subject_idx = 1:length(subject_list)
    subject_int = str2double(subject_list(subject_idx,:));
    p = patient(subject_int);
    [all, at_least_one_of_each, at_least_one_per_task, none] = p.which_fcs_exist(atlas, include_subcortical);
    
    start = tic();
    if p.exist_sc(atlas) &&  ~none %at_least_one_per_task
        
        fprintf('%d/%d: ...', subject_idx, length(subject_list));
        num_scans = num_scans + 1;
    
        %% populate struct
        subject_ids(num_scans)       = int32(subject_int); %correct index?
        data(num_scans).subject_id = subject_int;
        data(num_scans).sc         = p.sc(atlas, include_subcortical);
        
        data(num_scans).fc_REST1   = p.rest_fc(atlas, include_subcortical, 'REST1');
        data(num_scans).fc_REST2   = p.rest_fc(atlas, include_subcortical, 'REST2');
        
        [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = p.rest_scans(atlas, include_subcortical);
        data(num_scans).which_fcs = struct('lr1', lr_rest1.exist(), 'lr2', lr_rest2.exist(), 'rl1', rl_rest1.exist(), 'rl2', rl_rest2.exist());
  
        %optional viz
        %p.viz_fcs(atlas, include_subcortical, 99.5);
    end
            
	elapsed_time = toc(start);
 
    %elapsed_time = stop;

    num_iters_left = length(subject_list) - subject_idx;
    total_time = total_time + elapsed_time;
    ave_patient_time = total_time/num_scans;
    expected_time_left = ave_patient_time*num_iters_left/3600; %in hrs
    fprintf('elapsed time: %.2f | ave time: %.1f | expected time remain %.2f (hrs)=== %%\n', elapsed_time, ave_patient_time, expected_time_left);
end


data = data(1:num_scans);

end


function inspect_fcs_sc(p, atlas, include_subcortical, percentile)
    % display all fcs as well as overall average

    lr_rest1 = ScanInfo(p.subject_id, atlas, 'REST1', 'LR', include_subcortical);
    lr_rest2 = ScanInfo(p.subject_id, atlas, 'REST2', 'LR', include_subcortical);
    rl_rest1 = ScanInfo(p.subject_id, atlas, 'REST1', 'RL', include_subcortical);
    rl_rest2 = ScanInfo(p.subject_id, atlas, 'REST2', 'RL', include_subcortical);
    
    possible_scans = [lr_rest1, lr_rest2, rl_rest1, rl_rest2];
    actual_scans = [lr_rest1];
    actual_scans(end) = [];
    for scan_idx = 1:length(possible_scans)
        scan = possible_scans(scan_idx);
        if scan.exist()
            actual_scans(end+1) = scan;
        end
    end
    
    N = length(get_roi_idxs(actual_scans(1).atlas, actual_scans(1).include_subcortical));
    fcs_tensor = zeros(N, N, length(actual_scans));
    
    figure;
    t = tiledlayout('flow');
    axes_list = gobjects(length(actual_scans),1);
    
 
    for scan_idx = 1:length(actual_scans)
        scan = actual_scans(scan_idx);
        fc   = scan.compute_fc();
        fcs_tensor(:,:,scan_idx) = fc;

        ax = nexttile(t);
        axes_list(scan_idx) = ax;
        
        imagesc(ax, fc);
        set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
        yl = sprintf('%s || %s',scan.scan_direction, scan.task);
        ylh = ylabel(ax, yl, 'FontSize', 20);
        %ylp = get(ylh, 'Position');
        %ext = get(ylh,'Extent');
        %set(ylh, 'Rotation',0, 'Position',ylp-[2*ext(3) 0 0])
	
    
        xlim(ax,[1,N]); xlim(ax,'manual');
        ylim(ax,[1,N]); ylim(ax, 'manual');
        daspect(ax,[1 1 1]);
        %pbaspect(ax,[1 1 1]);
    end
    [cl] = find_cov_clim(fcs_tensor, percentile);
    colormap_file = load("viz_fcs/correlations_colormap.mat");
    set(axes_list, 'ColorMap', colormap_file.colorMap);
    set(axes_list,'CLim',[-cl cl]);
    cb = colorbar(axes_list(end),'East', 'AxisLocation','out');
    linkaxes(axes_list,'xy');
    
    % |/ overall fc ax
    overall_fc = p.full_rest_fc(atlas, include_subcortical);
    full_ax = nexttile(t);
    ax = full_ax;
    %axes_list(end) = ax;
    imagesc(ax, overall_fc);
	set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);
	ylh = ylabel(ax, 'Concat/Ave FC', 'FontSize', 20);

	xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    
    overall_fc_expand = zeros(N, N, 2);
    overall_fc_expand(:,:,1) = overall_fc;
    overall_fc_expand(:,:,2) = overall_fc;
    [cl_ave] = find_cov_clim(overall_fc_expand, percentile);
    colormap_file = load("viz_fcs/correlations_colormap.mat");
    set(full_ax, 'ColorMap', colormap_file.colorMap);
    set(full_ax,'CLim',[-cl_ave cl_ave]);
    cb = colorbar(full_ax,'East', 'AxisLocation','out');
    % ^^ overall fc ax
    
    ax = nexttile(t);
    title('raw sc');
    imagesc(ax, p.sc());
    xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    colormap(ax, 'gray');
    set(ax,'CLim',[0 max(max(sc))]);
    %cb = colorbar(ax,'East');
    colorbar(ax, 'AxisLocation','out')
    
    ax = nexttile(t);
    title('thresholded sc');
    imagesc(ax, p.sc()>0);
    xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    colormap(ax, 'gray');
    set(ax,'CLim',[0 1]);
    %cb = colorbar(ax,'East');
    colorbar(ax, 'AxisLocation','out')
    
end
