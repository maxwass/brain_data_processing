%% Inspect if means are different out of REST1/REST2
close all;

atlas = 'desikan';
include_subcortical = false;

p = patient(993675);

lr1_dt = p.load_functional_dtseries(atlas, 'REST1', 'LR', include_subcortical);
lr2_dt = p.load_functional_dtseries(atlas, 'REST2', 'LR', include_subcortical);

[mean_lr1, stdv_lr1] = mean_stdv_signal(lr1_dt);
[mean_lr2, stdv_lr2] = mean_stdv_signal(lr2_dt);

ax = axis();
plot_mean_vec(ax, mean_lr1, stdv_lr1, mean_lr2, stdv_lr2);


function [mean_roi, stdv_roi] = mean_stdv_signal(dtseries)
    mean_roi    = mean(dtseries, 2); %across rows
    stdv_roi    = std(dtseries, 0, 2);
end

function plot_mean_vec(ax, mean_data_1, stdv_data_1, mean_data_2, stdv_data_2);
    transp = .2;
    y = [mean_data_1, mean_data_2];
    
    x = 1:length(mean_data_1);
    b_lr = bar(x, y);
    %b_lr.FaceAlpha = transp;
	xlabel('ROI', 'FontSize', 20);
	ylabel('value', 'FontSize', 20);
	title('Ave Node Value: REST1 vs REST2', 'FontSize', 18);
    y_lr = yline(mean(mean_data_1),'-',sprintf('REST1 ave = %.1f',mean(mean_data_1)), 'LineWidth',1, 'Color', 'r', 'FontSize', 12);
    y_lr.LabelHorizontalAlignment = 'left';
    hold on;
    y_rl = yline(mean(mean_data_2),'-',sprintf('REST2 ave = %.1f',mean(mean_data_2)), 'LineWidth',1, 'Color', 'b', 'FontSize', 12);
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