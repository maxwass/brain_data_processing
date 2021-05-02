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
        
        function [fc] = full_rest_fc(obj, atlas, include_subcortical)
            lr_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'LR', include_subcortical);
            lr_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'LR', include_subcortical);
            rl_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'RL', include_subcortical);
            rl_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'RL', include_subcortical);
            
            full_lr_dtseries = concat_dtseries(lr_rest1, lr_rest2);
            full_rl_dtseries = concat_dtseries(rl_rest1, rl_rest2);
            
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
        
        function viz_fcs(obj, atlas, include_subcortical, percentile)
            [all, ~, ~] = obj.which_fcs_exist(atlas, include_subcortical);
            
            if all
                all_fc_plot(lr_rest1, lr_rest2, rl_rest1, rl_rest2, percentile);
            end
            
        end
        
        function [lr_rest1, lr_rest2, rl_rest1, rl_rest2] = rest_scans(obj, atlas, include_subcortical)
            lr_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'LR', include_subcortical);
            lr_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'LR', include_subcortical);
            rl_rest1 = ScanInfo(obj.subject_id, atlas, 'REST1', 'RL', include_subcortical);
            rl_rest2 = ScanInfo(obj.subject_id, atlas, 'REST2', 'RL', include_subcortical);
        end
        
        function [all, at_least_one_of_each, none] = which_fcs_exist(obj, atlas, include_subcortical)
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


function [full_dtseries] = concat_dtseries(scan_info_1, scan_info_2)

           
	if ~scan_info_1.exist() && ~scan_info_2.exist()
        full_dtseries = NaN;
	elseif ~scan_info_1.exist()
        full_dtseries = scan_info_2.load_functional_dtseries();
	elseif ~scan_info_2.exist()
        full_dtseries = scan_info_1.load_functional_dtseries();
    else
        dt1 = scan_info_1.load_functional_dtseries();
        dt2 = scan_info_2.load_functional_dtseries();
        full_dtseries = cat(2, dt1, dt2);
        [rows, cols] = size(full_dtseries);
        if rows > cols
            disp('ERROR: Concatenate incorrectly')
        end
	end
    
end

function all_fc_plot(lr_rest1, lr_rest2, rl_rest1, rl_rest2, percentile)
f = figure();
t = tiledlayout(2,4);


fc_lr1 = lr_rest1.compute_fc();
fc_lr2 = lr_rest2.compute_fc();
fc_lr_full = cov(concat_dtseries(lr_rest1, lr_rest2)');

fc_rl1 = rl_rest1.compute_fc();
fc_rl2 = rl_rest2.compute_fc();
fc_rl_full = cov(concat_dtseries(rl_rest1, rl_rest2)');

fc = (fc_lr_full + fc_rl_full)/2;


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