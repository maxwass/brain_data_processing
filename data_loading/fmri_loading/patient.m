classdef patient
    properties
        subject_id int64 {mustBePositive}
        
    end
    methods
        function obj = patient(subject_id)

            %filter subject
            if ~isnumeric(subject_id)
                obj.subject_id = str2double(subject_id);
            else
                obj.subject_id = subject_id;
            end

        end
        
        function dtseries = load_functional_dtseries(obj, atlas, task, scan_direction, include_subcortical)
            dtseries = ScanInfo(obj.subject_id, atlas, task, scan_direction, include_subcortical).load_functional_dtseries();
        end
        
        %% all different ways of computing FCs
        function [rest1_lr_fc, rest1_rl_fc, rest2_lr_fc, rest2_rl_fc] = individual_fc(obj, atlas, include_subcortical, which_fc)
            rest1_lr_fc = ScanInfo(obj.subject_id, atlas, 'REST1', 'LR', include_subcortical).compute_fc(which_fc);
            rest1_rl_fc = ScanInfo(obj.subject_id, atlas, 'REST1', 'RL', include_subcortical).compute_fc(which_fc);
            rest2_rl_fc = ScanInfo(obj.subject_id, atlas, 'REST2', 'RL', include_subcortical).compute_fc(which_fc);
            rest2_lr_fc = ScanInfo(obj.subject_id, atlas, 'REST2', 'LR', include_subcortical).compute_fc(which_fc);
        end
        
        function [total_fc] = concat_all_fc(obj, atlas, include_subcortical, which_fc)
            [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);
            total_dtseries = concat_dtseries({lr_rest1, rl_rest1, lr_rest2, rl_rest2});
            
            total_fc = compute_fc_if_exist(total_dtseries, which_fc);
        end
        
        function [rest1_fc, rest2_fc] = task_grouped_fc(obj, atlas, include_subcortical, which_fc)
            [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);

            rest1_dtseries = concat_dtseries({lr_rest1, rl_rest1});
            rest2_dtseries = concat_dtseries({lr_rest2, rl_rest2});
            
            rest1_fc = compute_fc_if_exist(rest1_dtseries, which_fc);
            rest2_fc = compute_fc_if_exist(rest2_dtseries, which_fc);
            
        end
        
        function [lr_fc, rl_fc] = direction_grouped_fc(obj, atlas, include_subcortical, which_fc)
            [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);

            lr_dtseries = concat_dtseries({lr_rest1, lr_rest2});
            rl_dtseries = concat_dtseries({rl_rest1, rl_rest2});
            
            lr_fc = compute_fc_if_exist(lr_dtseries, which_fc);
            rl_fc = compute_fc_if_exist(rl_dtseries, which_fc);
            
        end
        %{
        function [fc] = full_rest_fc(obj, atlas, include_subcortical, mean_norm)
            [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);
            
            full_lr_dtseries = concat_dtseries(lr_rest1, lr_rest2, mean_norm);
            full_rl_dtseries = concat_dtseries(rl_rest1, rl_rest2, mean_norm);
            
            if any(any(isnan(full_lr_dtseries))) && any(any(isnan(full_rl_dtseries)))
                fc = NaN;
            elseif any(any(isnan(full_lr_dtseries)))
                fc = cov(full_rl_dtseries');
            elseif any(any(isnan(full_rl_dtseries)))
                fc = cov(full_lr_dtseries');
            else
                fc_lr = cov(full_lr_dtseries');
                fc_rl = cov(full_rl_dtseries');
                fc = (fc_lr + fc_rl)/2;
                
                
                % do plotting here
                %{
                if lr_rest1.exist() && lr_rest2.exist() && rl_rest1.exist() && rl_rest2.exist()
                    
                    all_fc_plot(lr_rest1, lr_rest2, rl_rest1, rl_rest2, 99.5)
 
                end
                %}
                
            end
        end
        %}
        %
        
        function viz_fcs(obj, atlas, include_subcortical, percentile, which_fc)
            %{
            [all, ~, ~, ~] = obj.which_fcs_exist(atlas, include_subcortical);
            
            if all
                [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);
                all_fc_plot(lr_rest1, lr_rest2, rl_rest1, rl_rest2, mean_norm, percentile, which_fc);
            end
            %}
            f = figure();
            t = tiledlayout(2,6);
            title(t, sprintf('%d - %s', obj.subject_id, which_fc), 'FontSize', 25);
            fs_title = 25;
            
            [rest1_lr_fc, rest1_rl_fc, rest2_lr_fc, rest2_rl_fc] = obj.individual_fc(atlas, include_subcortical, which_fc);
            [total_fc] = obj.concat_all_fc(atlas, include_subcortical, which_fc);
            [rest1_fc, rest2_fc] = obj.task_grouped_fc(atlas, include_subcortical, which_fc);
            [lr_fc, rl_fc] = obj.direction_grouped_fc(atlas, include_subcortical, which_fc);
            
            all_fcs = {rest1_lr_fc, rest1_rl_fc, rest2_lr_fc, rest2_rl_fc, total_fc, rest1_fc, rest2_fc, lr_fc, rl_fc};
            non_nan_fcs = {};
            for l = 1:length(all_fcs)
               if sum(size(all_fcs{l}))>2
                   non_nan_fcs{end+1} = all_fcs{l};
               end
            end
            all_real_fcs = cat(3, non_nan_fcs{:});
            [cl] = find_cov_clim(all_real_fcs, percentile); % if corr, set to 100
            [N, ~, ~] = size(all_real_fcs);
            
            
            % plot all individuals
            fcs = {rest1_lr_fc, rest1_rl_fc, rest2_lr_fc, rest2_rl_fc};
            names = {'REST1-LR', 'REST1-RL', 'REST2-LR', 'REST2-RL'};
            for idx = 1:length(fcs)
                ax = nexttile(idx);
                fc_to_plot = fcs{idx};
                if sum(size(fc_to_plot))<=2
                    fc_to_plot = zeros(N, N);
                end
                imagesc(ax, fc_to_plot);
                adjust_colors_shape(ax, cl, N);
                title(names{idx}, 'FontSize', fs_title)
            end
            
            % plot total
            ax = nexttile([2, 2]);
            imagesc(ax, total_fc);
            adjust_colors_shape(ax, cl, N);
            title('total timeseries', 'FontSize', fs_title);
            cb_large = colorbar(ax, 'EastOutside', 'FontSize',25); 
            
            % plot all groups
            fcs = {rest1_fc, rest2_fc};
            names = {'REST1', 'REST2'};
            ax_offset = 6;
            for idx = 1:length(fcs)
                ax = nexttile(idx+ax_offset);
                fc_to_plot = fcs{idx};
                if sum(size(fc_to_plot))<=2
                    fc_to_plot = zeros(N, N);
                end
                imagesc(ax, fc_to_plot);
                adjust_colors_shape(ax, cl, N);
                title(names{idx}, 'FontSize', fs_title)
            end
            
            fcs = {lr_fc, rl_fc};
            names = {'LR', 'RL'};
            ax_offset = 8;
            for idx = 1:length(fcs)
                ax = nexttile(idx+ax_offset);
                fc_to_plot = fcs{idx};
                if sum(size(fc_to_plot))<=2
                    fc_to_plot = zeros(N, N);
                end
                imagesc(ax, fc_to_plot);
                adjust_colors_shape(ax, cl, N);
                title(names{idx}, 'FontSize', fs_title);
            end
            
            
        end
        
        function [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = rest_scans(obj, atlas, include_subcortical)
            lr_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'LR', include_subcortical);
            lr_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'LR', include_subcortical);
            rl_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'RL', include_subcortical);
            rl_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'RL', include_subcortical);
        end
        
        function [all, at_least_one_of_each, at_least_one_per_task, none] = which_fcs_exist(obj, atlas, include_subcortical)
            [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = obj.rest_scans(atlas, include_subcortical);
            
            all = false;
            if lr_rest1.exist() && lr_rest2.exist() && rl_rest1.exist() && rl_rest2.exist()
                all = true;
            end
            
            at_least_one_of_each = false;
            if (lr_rest1.exist() || lr_rest2.exist()) && (rl_rest1.exist() || rl_rest2.exist())
                at_least_one_of_each = true;
            end
            
            none = false;
            if ~lr_rest1.exist() && ~lr_rest2.exist() && ~rl_rest1.exist() && ~rl_rest2.exist()
                none = true;
            end
            
            at_least_one_per_task = false;
            if (lr_rest1.exist() || rl_rest1.exist()) && (lr_rest2.exist() || rl_rest2.exist())
               at_least_one_per_task = true; 
            end
        end
        
        function f = exist_sc(obj, atlas)
            try
                extract_sc(obj.subject_id, atlas, false);
                f = true;
            catch ME
                f = false;
            end
            
        end
        
        function A = sc(obj, atlas, include_subcortical)
            [A] = extract_sc(obj.subject_id, atlas, include_subcortical);
        end
         
        function [is_eq] = eq(obj, other)
            is_eq = isequal(obj.subject_id, other.subject_id);
            
            
        end
    end
end


function [fc_or_nan] = compute_fc_if_exist(dtseries, which_fc)
   
    if sum(size(dtseries))>2
        % time series exists. Which fc to compute?
        if strcmp(which_fc, 'cov')
            fc_or_nan = cov(dtseries');
        elseif strcmp(which_fc, 'corr')
            fc_or_nan = corr(dtseries');

        else
            error('unrecognized which_fc %s', which_fc);
        end
        
    else
        % time series is NaN
        fc_or_nan = NaN;
     end

end
% place all dtseries into a list.
% then concatentate them all
function [total_dtseries] = concat_dtseries(scan_info_list)

    dtseries = {};
    for l = 1:length(scan_info_list)
        scan_info = scan_info_list{l};
        
        if scan_info.exist() 
           scan_dtseries = scan_info.load_functional_dtseries();
           dtseries{l} = scan_dtseries - mean(scan_dtseries,2);
        end

    end
    
    if ~isempty(dtseries)
        total_dtseries = horzcat(dtseries{:});
    else 
        total_dtseries = NaN;
    end
    
end

function [full_dtseries] = concat_dtseries_2(scan_info_1, scan_info_2, mean_norm)

           
	if ~scan_info_1.exist() && ~scan_info_2.exist()
        full_dtseries = NaN;
	elseif ~scan_info_1.exist()
        full_dtseries = scan_info_2.load_functional_dtseries();
	elseif ~scan_info_2.exist()
        full_dtseries = scan_info_1.load_functional_dtseries();
    else
        dt1 = scan_info_1.load_functional_dtseries();
        dt2 = scan_info_2.load_functional_dtseries();
        if mean_norm
           full_dtseries = cat(2, dt1-mean(dt1,2), dt2-mean(dt2,2));
        else
            full_dtseries = cat(2, dt1, dt2);
        end
        [rows, cols] = size(full_dtseries);
        if rows > cols
            disp('ERROR: Concatenate incorrectly')
        end
	end
    
end

function all_fc_plot(lr_rest1, lr_rest2, rl_rest1, rl_rest2, mean_norm, percentile, which_fc)
f = figure();
t = tiledlayout(2,4);


if strcmp(which_fc,'cov')
    fc_lr1 = lr_rest1.compute_fc();
    fc_lr2 = lr_rest2.compute_fc();
    fc_lr_full = cov(concat_dtseries(lr_rest1, lr_rest2, mean_norm)');

    fc_rl1 = rl_rest1.compute_fc();
    fc_rl2 = rl_rest2.compute_fc();
    fc_rl_full = cov(concat_dtseries(rl_rest1, rl_rest2, mean_norm)');

    fc = (fc_lr_full + fc_rl_full)/2;
elseif strcmp(which_fc, 'corr')
    fc_lr1 = corrcov(lr_rest1.compute_fc());
    fc_lr2 = corrcov(lr_rest2.compute_fc());
    fc_lr_full = corrcov(cov(concat_dtseries(lr_rest1, lr_rest2, mean_norm)'));

    fc_rl1 = corrcov(rl_rest1.compute_fc());
    fc_rl2 = corrcov(rl_rest2.compute_fc());
    fc_rl_full = corrcov(cov(concat_dtseries(rl_rest1, rl_rest2, mean_norm)'));

    fc = (fc_lr_full + fc_rl_full)/2;
    
    percentile = 100; % now normalized between 0 and 1
else
    error('unrecognized which_fc argument %s', which_fc);
end


[cl] = find_cov_clim(cat(3, fc_lr1, fc_lr2, fc_rl1, fc_rl2), percentile);
[cl_large] = find_cov_clim(cat(3, fc_lr_full, fc_rl_full, fc), percentile);
[N, ~] = size(fc);

fs_x = 15;
fs_title = 30;

% lrs
ax = nexttile(1);
imagesc(ax, fc_lr1);
adjust_colors_shape(ax, cl, N);
title('REST1', 'FontSize', fs_title)
ylabel(ax, 'LR', 'FontSize', fs_title);
xlabel(sum_stat_cov_txt(fc_lr1), 'FontSize', fs_x);

ax = nexttile(2);
imagesc(ax, fc_lr2);
title('REST2','FontSize', fs_title);
adjust_colors_shape(ax, cl, N);
cb = colorbar(ax, 'SouthOutside', 'FontSize', 18);
xlabel(sum_stat_cov_txt(fc_lr2), 'FontSize', fs_x);


ax = nexttile(3);
imagesc(ax, fc_lr_full);
title('[REST1, REST2]', 'FontSize', fs_title);
adjust_colors_shape(ax, cl_large, N);
xlabel(sum_stat_cov_txt(fc_lr_full), 'FontSize', fs_x);

ax = nexttile(4, [2,1]);
imagesc(ax, fc);
title('Concat and Ave', 'FontSize', fs_title);
adjust_colors_shape(ax, cl_large, N);
cb_large = colorbar(ax, 'EastOutside', 'FontSize',25);
xlabel(sum_stat_cov_txt(fc), 'FontSize', fs_x);

% rls
ax = nexttile(5);
imagesc(ax, fc_rl1);
adjust_colors_shape(ax, cl, N);
ylabel(ax, 'RL', 'FontSize', fs_title);
xlabel(sum_stat_cov_txt(fc_rl1), 'FontSize', fs_x);

ax = nexttile(6);
imagesc(ax, fc_rl2);
adjust_colors_shape(ax, cl, N);
xlabel(sum_stat_cov_txt(fc_rl2), 'FontSize', fs_x);

ax = nexttile(7);
imagesc(ax, fc_rl_full);
adjust_colors_shape(ax, cl_large, N);
xlabel(sum_stat_cov_txt(fc_rl_full), 'FontSize', fs_x);

end



function adjust_colors_shape(ax, cl, N)
	set(ax,'xticklabel',[]); set(ax,'yticklabel',[]);

	xlim(ax,[1,N]); xlim(ax,'manual');
	ylim(ax,[1,N]); ylim(ax, 'manual');
	daspect(ax,[1 1 1]);
    
    colormap_file = load("viz_fcs/correlations_colormap.mat");
    set(ax, 'ColorMap', colormap_file.colorMap);
    set(ax,'CLim',[-cl cl]);
end
 

function [txt] = sum_stat_cov_txt(fc)
    min_ = min(fc, [], 'all');
    max_ = max(fc, [], 'all');
    median_ = median(fc, 'all');
    mean_ = mean(fc, 'all');
    
    txt = sprintf('min %.1e | median %.1e | max %.1e', min_, median_, max_); 
end